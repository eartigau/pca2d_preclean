"""A console_script entry point returns an exit status, not a result.

Both of these are wired into pyproject as console_scripts, so whatever main()
returns is handed to sys.exit. twoframe.main returned the two bases, which
printed a pair of arrays and exited 1 after a perfectly good fit; and
reconstruct.main returned the number of files written, so correcting 316
spectra exited 316. Nothing could be chained after either of them, and a
`set -e` script died at the first.
"""

import inspect

from pca2d import reconstruct, twoframe


def returns_only_none(func):
    """Every `return` in the function body is bare or returns None."""
    import ast

    tree = ast.parse(inspect.getsource(func))
    for node in ast.walk(tree):
        if isinstance(node, ast.Return):
            if node.value is None:
                continue
            if isinstance(node.value, ast.Constant) and node.value.value is None:
                continue
            return False
    return True


def test_the_fit_entry_point_returns_nothing():
    assert returns_only_none(twoframe.main)


def test_the_apply_entry_point_returns_nothing():
    assert returns_only_none(reconstruct.main)
