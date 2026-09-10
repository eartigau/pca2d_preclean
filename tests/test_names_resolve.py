"""Every name a module calls at run time is one it can actually reach.

Written the day a run spent fourteen minutes on the fit and then died on
`NameError: name 'log' is not defined`, in four figure scripts at once: the
commit that put a timestamp on every printed line replaced their `print` with
`log` and imported it in none of them. Nothing caught it because the figure
scripts are launched as subprocesses at the end of a long run, which is the
worst possible place to learn that a module does not import.

This is a static check and costs milliseconds. It does not replace running the
thing; it removes the class of failure where running the thing is the only way
to find out.
"""

import ast
import importlib
import pathlib

import pytest

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "pca2d"
MODULES = sorted(PACKAGE.rglob("*.py"))


def toplevel_names(tree):
    """What a module's own body binds: imports, defs, classes, assignments."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            names.update(a.asname or a.name for a in node.names)
        elif isinstance(node, ast.Import):
            names.update((a.asname or a.name).split(".")[0] for a in node.names)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    names.add(target.id)
        elif isinstance(node, (ast.arg, ast.Name)) and isinstance(
                getattr(node, "ctx", None), ast.Store):
            names.add(getattr(node, "id", getattr(node, "arg", "")))
    return names


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_the_logger_is_imported_wherever_it_is_called(path):
    tree = ast.parse(path.read_text())
    calls = {node.func.id for node in ast.walk(tree)
             if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)}
    if "log" not in calls:
        return
    assert "log" in toplevel_names(tree), (
        "%s calls log() and imports it from nowhere; it will raise NameError"
        " the first time that line is reached" % path.name)


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_nothing_local_is_called_log(path):
    """The other half of the same commit's damage.

    bundle.run(cmd, log) took the list of failures in a parameter called `log`,
    which shadowed the logger every line of that function then called: a
    failing figure raised TypeError instead of landing on the bundle's last
    page. A name that is a function at module level and a list inside one
    function is a trap whether or not it is being sprung today.
    """
    tree = ast.parse(path.read_text())
    imports_log = any(
        isinstance(node, ast.ImportFrom)
        and any(alias.name == "log" for alias in node.names)
        for node in ast.walk(tree))
    if not imports_log:
        return
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        arguments = [a.arg for a in node.args.args + node.args.kwonlyargs]
        assigned = [t.id for inner in ast.walk(node)
                    if isinstance(inner, ast.Assign)
                    for t in inner.targets if isinstance(t, ast.Name)]
        assert "log" not in arguments + assigned, (
            "%s:%d %s() rebinds `log`, which the module imports as the logger"
            % (path.name, node.lineno, node.name))


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_every_module_imports(path):
    """Import each one for real, by its dotted name.

    By name and not by file path: inside a package, `from .logger import log`
    only resolves when the module knows which package it is in, and loading it
    from a path does not tell it.
    """
    if path.name == "__init__.py":
        return
    dotted = ".".join(path.relative_to(PACKAGE.parent).with_suffix("").parts)
    importlib.import_module(dotted)
