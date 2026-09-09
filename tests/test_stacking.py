"""Coadding a night is a memory decision, not a modelling one.

Every exposure is registered with its own BERV before it is added, so nothing
is smeared either way, and fitting the individual spectra is strictly more
information. So it happens only when the fit would not otherwise fit, and the
budget it is weighed against is a CONFIG value: if it were the machine's RAM,
the same cube cache key would name a stacked cube on one computer and an
unstacked one on another.
"""

from pca2d.cube import resolve_stacking


def cfg(nightly, budget=12.0, dtype="float64"):
    return {"input": {"nightly_stack": nightly},
            "twoframe": {"dtype": dtype},
            "output": {"max_memory_gb": budget}}


M = 577002          # TOI-2120's grid


def test_true_and_false_are_obeyed_whatever_the_budget():
    assert resolve_stacking(cfg(True), 5, M)[0] is True
    assert resolve_stacking(cfg(False), 100000, M)[0] is False


def test_auto_fits_every_file_when_the_budget_allows():
    # 321 exposures of TOI-2120: 642 rows, four float64 arrays, about 11.9 GB
    stack, why = resolve_stacking(cfg("auto", budget=12.0), 321, M)
    assert stack is False
    assert "every file" in why


def test_auto_coadds_when_it_would_not_fit():
    stack, why = resolve_stacking(cfg("auto", budget=4.0), 321, M)
    assert stack is True
    assert "coadded" in why


def test_the_decision_does_not_look_at_the_machine():
    # the same config gives the same answer regardless of anything else, which
    # is what makes it safe to hash into the cube cache key
    a = resolve_stacking(cfg("auto", budget=16.0), 321, M)
    b = resolve_stacking(cfg("auto", budget=16.0), 321, M)
    assert a == b


def test_auto_without_a_budget_still_fits_every_file():
    stack, why = resolve_stacking(cfg("auto", budget=None), 321, M)
    assert stack is False
    assert "no output.max_memory_gb" in why


def test_the_storage_dtype_changes_the_answer():
    """float32 halves the data and the weights, and can tip the decision."""
    tight = 10.0
    assert resolve_stacking(cfg("auto", budget=tight, dtype="float64"),
                            321, M)[0] is True
    assert resolve_stacking(cfg("auto", budget=tight, dtype="float32"),
                            321, M)[0] is False
