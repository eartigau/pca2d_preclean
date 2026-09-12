"""A window to run the pipeline, for people who should not have to read the
command line first.

    pca2d-gui

It shows the objects that are in the data root, lets one of them be run on its
own or several of them together against one observer basis, exposes the
handful of settings that change a result, writes the command it is about to
run so that it can be copied into a terminal, and streams the run's own log
into the window with the colours it would have in a terminal.

Nothing here decides anything: every option maps to one key of config.yaml or
one flag of `pca2d-preclean`, the command is shown before it runs, and the
settings can be exported as a variant file, which is how a run is made
reproducible (variants/README.md).
"""
from __future__ import annotations

import glob
import json
import os
import queue
import re
import subprocess
import sys
import threading

HOME_STATE = os.path.expanduser("~/.pca2d_gui.json")
ANSI = re.compile(r"\033\[(\d+)m")
#: the logger's colours, and how this window paints them
LEVELS = {"32": ("info", "#1a7f37"), "34": ("value", "#0a58ca"),
          "33": ("warn", "#b26a00"), "31": ("error", "#c1121f")}
STAGES = ("cube", "fit", "figures", "correct", "lbl")
#: what the window can change, and the config key each one is
OPTIONS = [
    ("n_star", "star components", "twoframe.n_star", "int",
     "0 keeps the star fixed at its per-parity template, which is the nominal"),
    ("n_earth", "observer components", "twoframe.n_earth", "int", ""),
    ("mean", "static part", "twoframe.mean", ("star", "offset", "full", "iterate"),
     "star: one star spectrum per order parity, taken out before the fit"),
    ("star_basis", "star basis", "twoframe.star_basis", ("spline", "grid"), ""),
    ("velocity_term", "fit a velocity per exposure", "twoframe.velocity_term",
     "bool", "fitted, never divided out"),
    ("iters", "sweeps at most", "twoframe.iters", "int", ""),
    ("shrink", "divide only what is significant", "correct.shrink", "bool",
     "each observer component, column by column, all exposures together"),
    ("mask", "samples a corrected file blanks", "correct.mask",
     ("common", "exposure", "none"),
     "common: one set of lines for the whole campaign"),
    ("width_kms", "high pass (km/s)", "highpass.width_kms", "float", ""),
    ("dv", "grid step (km/s)", "domain.dv", "float", ""),
    ("nightly_stack", "coadd each night", "input.nightly_stack",
     ("auto", "true", "false"), "a memory decision, never a modelling one"),
    ("run", "run LBL (hours)", "lbl.run", "bool", ""),
]


def objects_in(root, pattern="*t.fits"):
    """[(name, number of files)] for every object folder under the data root."""
    if not root or not os.path.isdir(root):
        return []
    out = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            continue
        n = len(glob.glob(os.path.join(path, pattern)))
        if n:
            out.append((name, n))
    return out


def instrument_of(root, name, pattern="*t.fits"):
    """INSTRUME of the first file that opens, or '?'."""
    from astropy.io import fits
    for path in sorted(glob.glob(os.path.join(root, name, pattern)))[:3]:
        try:
            with fits.open(path) as hdulist:
                value = str(hdulist[0].header.get("INSTRUME", "")).strip()
        except Exception:                                     # noqa: BLE001
            continue
        if value:
            return value
    return "?"


def build_command(state):
    """The `pca2d-preclean` command a set of choices means, as a list.

    One object or several, the stages asked for, the variant if one is chosen,
    and the two counts, which are flags rather than config keys so that a run
    says in its own name what it fitted.
    """
    argv = ["pca2d-preclean"]
    objects = list(state.get("objects") or [])
    if len(objects) > 1:
        argv += ["--objects", ",".join(objects)]
    elif objects:
        argv += ["--object", objects[0]]
    if state.get("config"):
        argv += ["--config", state["config"]]
    if state.get("data_dir"):
        argv += ["--data-dir", state["data_dir"]]
    if state.get("out_dir"):
        argv += ["--out-dir", state["out_dir"]]
    if state.get("variant") and state["variant"] != "(none)":
        argv += ["--variant", state["variant"]]
    if state.get("n_star") not in (None, ""):
        argv += ["--n-star", str(state["n_star"])]
    if state.get("n_earth") not in (None, ""):
        argv += ["--n-earth", str(state["n_earth"])]
    stages = [s for s in STAGES if state.get("stage_" + s)]
    if stages and len(stages) != len(STAGES):
        argv += ["--stages", ",".join(stages)]
    if state.get("dry_run"):
        argv.append("--dry-run")
    return argv


def variant_yaml(state, defaults=None):
    """The settings that differ from the configuration, as a variant file.

    Only what was changed: a variant is the nominal plus its own lines, and a
    file that repeats the nominal says nothing (variants/README.md).
    """
    defaults = defaults or {}
    out = {}
    for key, _label, path, kind, _help in OPTIONS:
        if key not in state or state[key] in (None, ""):
            continue
        value = state[key]
        if kind == "int":
            value = int(value)
        elif kind == "float":
            value = float(value)
        elif kind == "bool":
            value = bool(value)
        elif value in ("true", "false"):
            value = value == "true"
        section, name = path.split(".")
        if defaults.get(section, {}).get(name, object()) == value:
            continue
        out.setdefault(section, {})[name] = value
    return out


def _read_state():
    try:
        with open(HOME_STATE) as handle:
            return json.load(handle)
    except Exception:                                         # noqa: BLE001
        return {}


def _write_state(state):
    try:
        with open(HOME_STATE, "w") as handle:
            json.dump(state, handle, indent=1)
    except OSError:
        pass


class App:
    """The window itself."""

    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk

        self.tk, self.ttk = tk, ttk
        self.root = root
        self.proc = None
        self.lines = queue.Queue()
        self.saved = _read_state()
        root.title("pca2d-preclean")
        root.geometry("1180x760")
        style = ttk.Style()
        for theme in ("aqua", "clam", "default"):
            if theme in style.theme_names():
                style.theme_use(theme)
                break
        style.configure("Head.TLabel", font=("Helvetica", 13, "bold"))
        style.configure("Hint.TLabel", foreground="#666")

        self.vars = {}
        self._build_top(root)
        panes = ttk.Panedwindow(root, orient="horizontal")
        panes.pack(fill="both", expand=False, padx=10)
        left, right = ttk.Frame(panes), ttk.Frame(panes)
        panes.add(left, weight=1)
        panes.add(right, weight=2)
        self._build_objects(left)
        self._build_options(right)
        self._build_command(root)
        self._build_log(root)
        self.refresh_objects()
        self.root.after(80, self._drain)
        self.root.protocol("WM_DELETE_WINDOW", self._close)

    # ---- widgets ------------------------------------------------------
    def _build_top(self, parent):
        ttk, tk = self.ttk, self.tk
        frame = ttk.Frame(parent)
        frame.pack(fill="x", padx=10, pady=(10, 6))
        ttk.Label(frame, text="pca2d-preclean", style="Head.TLabel").grid(
            row=0, column=0, sticky="w")
        ttk.Label(frame, style="Hint.TLabel",
                  text="two-frame precleaning, then LBL. Pick a data root, pick"
                       " objects, look at the command, run it.").grid(
            row=0, column=1, columnspan=3, sticky="w", padx=8)
        for i, (key, label, default) in enumerate((
                ("data_dir", "data root", self.saved.get("data_dir", "data")),
                ("config", "config", self.saved.get("config", "config.yaml")),
                ("out_dir", "output root (optional)",
                 self.saved.get("out_dir", "")))):
            ttk.Label(frame, text=label).grid(row=i + 1, column=0, sticky="w",
                                              pady=2)
            var = tk.StringVar(value=default)
            self.vars[key] = var
            entry = ttk.Entry(frame, textvariable=var, width=70)
            entry.grid(row=i + 1, column=1, sticky="we", padx=6)
            var.trace_add("write", lambda *_: self._sync())
            ttk.Button(frame, text="Browse",
                       command=lambda k=key: self._browse(k)).grid(row=i + 1,
                                                                   column=2)
        ttk.Button(frame, text="Rescan", command=self.refresh_objects).grid(
            row=1, column=3, padx=4)
        frame.columnconfigure(1, weight=1)

    def _build_objects(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text="objects")
        box.pack(fill="both", expand=True, pady=4)
        self.tree = ttk.Treeview(box, columns=("files", "instrument"),
                                 show="tree headings", selectmode="extended",
                                 height=11)
        self.tree.heading("#0", text="object")
        self.tree.heading("files", text="files")
        self.tree.heading("instrument", text="instrument")
        self.tree.column("#0", width=170)
        self.tree.column("files", width=60, anchor="e")
        self.tree.column("instrument", width=90)
        self.tree.pack(fill="both", expand=True, padx=6, pady=6)
        self.tree.bind("<<TreeviewSelect>>", lambda _e: self._sync())
        ttk.Label(box, style="Hint.TLabel",
                  text="select several to fit them together against ONE observer"
                       " basis; each keeps its own star spectrum").pack(
            anchor="w", padx=6, pady=(0, 6))

    def _build_options(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text="settings")
        box.pack(fill="both", expand=True, pady=4)
        grid = ttk.Frame(box)
        grid.pack(fill="both", expand=True, padx=8, pady=6)
        for i, (key, label, path, kind, hint) in enumerate(OPTIONS):
            row, col = i % 6, (i // 6) * 3
            ttk.Label(grid, text=label).grid(row=row, column=col, sticky="w",
                                             pady=2)
            default = self.saved.get(key, self._config_default(path, kind))
            if kind == "bool":
                var = tk.BooleanVar(value=bool(default))
                widget = ttk.Checkbutton(grid, variable=var)
            elif isinstance(kind, tuple):
                var = tk.StringVar(value=str(default))
                widget = ttk.Combobox(grid, textvariable=var, values=list(kind),
                                      width=10, state="readonly")
            else:
                var = tk.StringVar(value="" if default is None else str(default))
                widget = ttk.Entry(grid, textvariable=var, width=12)
            widget.grid(row=row, column=col + 1, sticky="w", padx=(6, 18))
            var.trace_add("write", lambda *_: self._sync())
            self.vars[key] = var
            if hint:
                self._tooltip(widget, hint)
        run = ttk.Frame(box)
        run.pack(fill="x", padx=8, pady=(0, 8))
        ttk.Label(run, text="stages").grid(row=0, column=0, sticky="w")
        for i, stage in enumerate(STAGES):
            var = tk.BooleanVar(value=self.saved.get("stage_" + stage, True))
            self.vars["stage_" + stage] = var
            ttk.Checkbutton(run, text=stage, variable=var,
                            command=self._sync).grid(row=0, column=i + 1, padx=3)
        ttk.Label(run, text="variant").grid(row=1, column=0, sticky="w", pady=4)
        self.vars["variant"] = tk.StringVar(value=self.saved.get("variant",
                                                                 "(none)"))
        self.variant_box = ttk.Combobox(run, textvariable=self.vars["variant"],
                                        values=self._variants(), width=18,
                                        state="readonly")
        self.variant_box.grid(row=1, column=1, columnspan=3, sticky="w")
        self.vars["variant"].trace_add("write", lambda *_: self._sync())

    def _build_command(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text="the command this runs")
        box.pack(fill="x", padx=10, pady=6)
        self.command = tk.Text(box, height=2, wrap="word",
                               font=("Menlo", 11), background="#f6f6f6",
                               relief="flat")
        self.command.pack(fill="x", padx=6, pady=6)
        bar = ttk.Frame(parent)
        bar.pack(fill="x", padx=10)
        self.run_button = ttk.Button(bar, text="Run", command=self.start)
        self.run_button.pack(side="left")
        self.stop_button = ttk.Button(bar, text="Stop", command=self.stop,
                                      state="disabled")
        self.stop_button.pack(side="left", padx=4)
        ttk.Button(bar, text="Dry run", command=lambda: self.start(dry=True)).pack(
            side="left", padx=4)
        ttk.Button(bar, text="Export YAML...", command=self.export).pack(
            side="left", padx=12)
        ttk.Button(bar, text="Save log...", command=self.save_log).pack(side="left")
        ttk.Button(bar, text="Open outputs", command=self.open_outputs).pack(
            side="left", padx=4)
        self.status = ttk.Label(bar, text="idle", style="Hint.TLabel")
        self.status.pack(side="right")

    def _build_log(self, parent):
        ttk, tk = self.ttk, self.tk
        box = ttk.Labelframe(parent, text="output")
        box.pack(fill="both", expand=True, padx=10, pady=8)
        self.log = tk.Text(box, wrap="none", font=("Menlo", 11),
                           background="white", relief="flat")
        bar = ttk.Scrollbar(box, command=self.log.yview)
        self.log.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        self.log.pack(fill="both", expand=True, padx=6, pady=6)
        for _code, (name, colour) in LEVELS.items():
            self.log.tag_configure(name, foreground=colour)
        self.log.tag_configure("plain", foreground="#222")

    def _tooltip(self, widget, text):
        tk = self.tk

        def enter(_event):
            self.status.configure(text=text)

        def leave(_event):
            self.status.configure(text="idle" if not self.proc else "running")

        widget.bind("<Enter>", enter)
        widget.bind("<Leave>", leave)

    # ---- state --------------------------------------------------------
    def _config_default(self, path, kind):
        section, name = path.split(".")
        cfg = self._config()
        value = (cfg.get(section) or {}).get(name)
        if value is None and kind == "bool":
            return False
        return value

    def _config(self):
        if getattr(self, "_cfg", None) is None:
            try:
                from .config import load_config
                self._cfg = load_config(self.vars["config"].get()
                                        if "config" in self.vars else "config.yaml")
            except Exception:                                 # noqa: BLE001
                self._cfg = {}
        return self._cfg

    def _variants(self):
        folder = os.path.join(os.path.dirname(
            os.path.abspath(self.vars["config"].get() if "config" in self.vars
                            else "config.yaml")), "variants")
        names = sorted(os.path.splitext(os.path.basename(p))[0]
                       for p in glob.glob(os.path.join(folder, "*.yaml")))
        return ["(none)"] + names

    def state(self):
        out = {"objects": [self.tree.item(i, "text")
                           for i in self.tree.selection()]}
        for key, var in self.vars.items():
            out[key] = var.get()
        return out

    def _sync(self, *_args):
        state = self.state()
        self.command.delete("1.0", "end")
        self.command.insert("1.0", " ".join(build_command(state)))
        _write_state({k: v for k, v in state.items() if k != "objects"})

    def refresh_objects(self):
        root = self.vars["data_dir"].get()
        self.tree.delete(*self.tree.get_children())
        found = objects_in(root)
        for name, n in found:
            self.tree.insert("", "end", text=name, values=(n, "..."))
        self._cfg = None
        self.variant_box.configure(values=self._variants())
        self._sync()
        if found:
            threading.Thread(target=self._instruments, args=(root, found),
                             daemon=True).start()

    def _instruments(self, root, found):
        for item, (name, _n) in zip(self.tree.get_children(), found):
            value = instrument_of(root, name)
            self.root.after(0, lambda i=item, v=value:
                            self.tree.set(i, "instrument", v))

    # ---- running ------------------------------------------------------
    def start(self, dry=False):
        if self.proc is not None:
            return
        state = self.state()
        if not state["objects"]:
            self._write("no object selected\n", "error")
            return
        state["dry_run"] = dry
        instruments = {self.tree.set(i, "instrument")
                       for i in self.tree.selection()} - {"...", "?", ""}
        if len(instruments) > 1:
            self._write("these objects are not from one instrument (%s): a joint"
                        " fit needs one domain and one grid, and the run will"
                        " refuse them\n" % ", ".join(sorted(instruments)), "warn")
        argv = build_command(state)
        self._write("\n%s\n" % (" ".join(argv)), "value")
        env = dict(os.environ, PCA2D_COLOUR="1", PYTHONUNBUFFERED="1")
        try:
            self.proc = subprocess.Popen(
                argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=env,
                cwd=os.path.dirname(os.path.abspath(self.vars["config"].get()))
                or None)
        except OSError as exc:
            self._write("cannot start it: %s\n" % exc, "error")
            self.proc = None
            return
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.status.configure(text="running")
        threading.Thread(target=self._reader, daemon=True).start()

    def _reader(self):
        for line in self.proc.stdout:
            self.lines.put(line)
        code = self.proc.wait()
        self.lines.put("\n[finished with code %d]\n" % code)
        self.proc = None
        self.root.after(0, self._finished)

    def _finished(self):
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        self.status.configure(text="idle")

    def stop(self):
        if self.proc is not None:
            self.proc.terminate()
            self._write("\n[stopped]\n", "warn")

    def _drain(self):
        while True:
            try:
                line = self.lines.get_nowait()
            except queue.Empty:
                break
            self._write(line)
        self.root.after(80, self._drain)

    def _write(self, line, forced=None):
        tag = forced or "plain"
        if forced is None:
            codes = ANSI.findall(line)
            for code in codes:
                if code in LEVELS:
                    tag = LEVELS[code][0]
                    break
        text = ANSI.sub("", line)
        if text.startswith("\r"):
            self.log.delete("end-2l", "end-1l")
            text = text.lstrip("\r")
        self.log.insert("end", text, tag)
        self.log.see("end")
        if "| stage " in text:
            self.status.configure(text=text.split("| ", 1)[-1].strip()[:60])

    # ---- the two buttons that write things ----------------------------
    def export(self):
        from tkinter import filedialog, messagebox

        import yaml
        state = self.state()
        body = variant_yaml(state, self._config())
        if not body:
            messagebox.showinfo("nothing to export",
                                "every setting is the configuration's own, so a"
                                " variant file would say nothing.")
            return
        path = filedialog.asksaveasfilename(
            title="save these settings as a variant",
            initialdir=os.path.join(os.path.dirname(os.path.abspath(
                self.vars["config"].get())), "variants"),
            defaultextension=".yaml", initialfile="mine.yaml")
        if not path:
            return
        head = ("# written by pca2d-gui: the nominal configuration plus these\n"
                "# lines. Run it with: pca2d-preclean --object NAME --variant %s\n"
                % os.path.splitext(os.path.basename(path))[0])
        with open(path, "w") as handle:
            handle.write(head)
            yaml.safe_dump(body, handle, sort_keys=False, default_flow_style=False)
        self._write("wrote %s\n" % path, "value")
        self.variant_box.configure(values=self._variants())

    def open_outputs(self):
        """Show the output root in the file browser: a run leaves one PDF and a
        folder of corrected spectra, and they are easier to find by looking."""
        root = self.vars["out_dir"].get() or (
            (self._config().get("output") or {}).get("directory") or "outputs")
        path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(
            self.vars["config"].get())), root))
        if not os.path.isdir(path):
            self._write("no %s yet\n" % path, "warn")
            return
        opener = ("open" if sys.platform == "darwin"
                  else "explorer" if os.name == "nt" else "xdg-open")
        subprocess.Popen([opener, path])

    def save_log(self):
        from tkinter import filedialog
        path = filedialog.asksaveasfilename(title="save the output",
                                            defaultextension=".log")
        if path:
            with open(path, "w") as handle:
                handle.write(self.log.get("1.0", "end"))
            self._write("wrote %s\n" % path, "value")

    def _browse(self, key):
        from tkinter import filedialog
        if key == "config":
            path = filedialog.askopenfilename(title="configuration",
                                              filetypes=[("YAML", "*.yaml")])
        else:
            path = filedialog.askdirectory(title=key.replace("_", " "))
        if path:
            self.vars[key].set(path)
            if key == "data_dir":
                self.refresh_objects()

    def _close(self):
        if self.proc is not None:
            self.proc.terminate()
        self.root.destroy()


def main(argv=None):
    import tkinter as tk

    root = tk.Tk()
    App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
