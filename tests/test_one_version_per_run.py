"""A run is one version of the code, whatever is edited while it goes.

On 2026-09-16 a joint run started at 10:06, the code was edited at 10:46, and
at 10:52 its LBL stage imported lbl.py for the first time: the new file, which
asked the config.py loaded at 10:06 for a function that did not exist yet.
"""

import subprocess
import sys
import time

from pca2d import provenance


def test_every_stage_module_is_loaded_before_the_first_stage():
    code = ("import sys; from pca2d import provenance; provenance.freeze();"
            "need = ['pca2d.lbl', 'pca2d.lbltemplate', 'pca2d.rvpages',"
            " 'pca2d.twoframe', 'pca2d.build', 'pca2d.reconstruct',"
            " 'pca2d.figures.bundle', 'pca2d.storage', 'pca2d.joint'];"
            "missing = [m for m in need if m not in sys.modules];"
            "print('missing', missing);"
            "print('tk', 'tkinter' in sys.modules);"
            "sys.exit(1 if missing else 0)")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, timeout=120)
    assert done.returncode == 0, done.stdout + done.stderr
    assert "tk False" in done.stdout, "a run has no screen and needs no tkinter"


def test_a_source_edited_mid_run_is_noticed(tmp_path):
    (tmp_path / "a.py").write_text("x = 1\n")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b.py").write_text("y = 2\n")
    (tmp_path / "__pycache__").mkdir()
    (tmp_path / "__pycache__" / "a.cpython.pyc.py").write_text("")
    snapshot = provenance.fingerprint(str(tmp_path))
    assert "__pycache__" not in " ".join(snapshot)
    assert provenance.drift(snapshot, str(tmp_path)) == []

    time.sleep(0.01)
    (tmp_path / "a.py").write_text("x = 10\n")
    (tmp_path / "c.py").write_text("z = 3\n")
    assert provenance.drift(snapshot, str(tmp_path)) == ["a.py", "c.py"]


def test_the_warning_is_said_once(monkeypatch):
    from pca2d import cli

    said, lines = [], []
    monkeypatch.setattr(provenance, "drift", lambda snapshot: ["lbl.py"])
    monkeypatch.setattr(cli, "log", lambda text, level="info": lines.append(text))
    cli.say_if_code_moved({}, "stage lbl", said)
    cli.say_if_code_moved({}, "stage lbl", said)
    assert len(lines) == 1 and "lbl.py" in lines[0]
    assert "stage lbl" in lines[0]
