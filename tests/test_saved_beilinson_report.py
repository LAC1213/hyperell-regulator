"""Portable exports retain precision, scope, conventions and failed attempts."""

from copy import deepcopy
from decimal import Decimal
import hashlib
import json

import pytest

from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.experiments import beilinson_report as report


def saved_row(label="249.a.249.1", cond=249):
    return dict(label=label, cond=cond, status="ok", regulator_digits=100,
                L_digits=100, L_second_derivative="5.067627040347058713642075941",
                det="33.48773151266455986457282728", c_plus="6.445925843562459201",
                c_minus="9.137006586242914998", model="u^6 + 4*u^5 + 4*u^4 + 2*u^3 + 1",
                triple=dict(P=["-1", "0"], Q=["infinity", "-1"], R=["0", "-1"],
                            order_P=7, order_Q=14, verdict_P=True, verdict_Q=True),
                comparison=dict(normalisation=report.NORMALISATION, ratio="28.0000000000000000000000000000000000000001",
                                rational="28", relative_agreement="3.57142857e-42"),
                settings=dict(report_version=b.REPORT_VERSION, regulator_digits=100,
                              lvalue_digits=100, with_lvalue=True, allow_undecided=True))


def test_scope_precision_policy_and_failures():
    current = saved_row()
    high_conductor = saved_row("1300.a.1.1", 1300)
    low_precision = saved_row("277.a.277.1", 277)
    low_precision["settings"]["regulator_digits"] = 15
    obsolete = saved_row("295.a.295.1", 295)
    obsolete["settings"]["report_version"] = b.REPORT_VERSION - 1
    failure = dict(label="997.a.997.1", cond=997, status="numerical_failure",
                   error="quadrature limit", settings=deepcopy(current["settings"]))
    selected = report.select_results([high_conductor, low_precision, failure, obsolete, current])
    assert [r["label"] for r in selected] == [current["label"], failure["label"]]
    assert selected[0]["comparison"] == current["comparison"]
    assert Decimal(selected[0]["absolute_ratio_error"]) == Decimal("1e-40")
    assert "absolute_ratio_error" not in current
    assert selected[1] == failure


def test_rejects_stale_conventions_and_insufficient_lvalue_precision():
    row = saved_row()
    row["comparison"]["normalisation"] = "2pi c+ c- L''(std,1)"
    with pytest.raises(ValueError, match="convention"):
        report.select_results([row])
    row = saved_row()
    row["L_digits"] = 15
    with pytest.raises(ValueError, match="100-digit L-value"):
        report.select_results([row])


def test_exports_exact_current_model_and_functions_without_integrating(tmp_path):
    row = saved_row()
    source = tmp_path/"results.json"
    source.write_text(json.dumps([row]))
    original = source.read_bytes()
    cache = tmp_path/"curves.json"
    # Original LMFDB model: y^2+(x^3+1)y=x^2+x.
    cache.write_text(json.dumps([dict(label=row["label"], cond=249,
                                    eqn=[[0, 1, 1], [1, 0, 0, 1]])]))
    out = tmp_path/"reports"
    report.export(source, out, cache)
    assert source.read_bytes() == original
    rows = json.loads((out/"beilinson_100_digits.json").read_text())
    exported = rows[0]
    assert exported["det"] == row["det"]
    assert exported["algebraic_data"]["image_points"] == [["-1", "0"], ["0", "-1"], "infinity"]
    assert exported["algebraic_data"]["shift"] == "0"
    assert exported["algebraic_data"]["scale"] == "1"
    assert set(exported["algebraic_data"]["functions"]) == {"P", "Q"}
    tex = (out/"beilinson_100_digits.tex").read_text()
    assert r"q=\pi c^+c^-L^{\prime\prime}(\mathrm{std},1)/D" in tex
    assert r"2\operatorname{Re}(z)" in tex
    assert hashlib.sha256(original).hexdigest() in tex
    inputs = report.explorer_inputs(out/"beilinson_100_digits.json", tmp_path/"inputs.json")
    entry = inputs["curves"][0]
    assert entry["numerical"]["det"] == row["det"]
    # A rational guess with only 42 matching digits is retained as a guess.
    assert entry["q"] is None and not entry["comparison_verified"]
    assert entry["rational_guess"] == "28"
    assert entry["quotient"] == row["comparison"]["ratio"]
    assert entry["pole"] == "infinity" and len(entry["exact_functions"]) == 2
    assert "u" not in entry["F"]
    rows[0]["comparison"]["relative_agreement"] = "1e-101"
    (out/"beilinson_100_digits.json").write_text(json.dumps(rows))
    verified = report.explorer_inputs(out/"beilinson_100_digits.json", tmp_path/"inputs.json")
    assert verified["curves"][0]["q"] == "28"
    assert verified["curves"][0]["comparison_verified"]


def test_tex_preserves_undecided_verdict_and_reports_failure():
    row = saved_row()
    row["triple"]["verdict_P"] = None
    row["algebraic_data"] = dict(model_latex="u^5+1", shift=None, scale="1", pole="infinity",
                                 functions={k:dict(latex="u+v") for k in ("P", "Q")})
    failure = dict(label="997.a.997.1", status="numerical_failure", error="cut_system_failed & timeout")
    tex = report.render([row, failure], "abc123")
    assert row["label"] + r"$^{\dagger}$" in tex
    assert r"cut\_system\_failed \& timeout" in tex
    assert "997.a.997.1" in tex


def test_unshifted_model_keeps_its_original_variable_and_infinity_pole():
    row = saved_row("277.a.277.1", 277)
    row["model"] = "x^6 + 2*x^5 + 3*x^4 + 4*x^3 - x^2 - 2*x + 1"
    row["triple"] = dict(P=["-1", "0"], Q=["0", "-1"], R=["infinity", "-1"],
                         order_P=15, order_Q=3)
    curve = dict(label=row["label"], cond=277, eqn=[[0, -1, -1], [1, 1, 1, 1]])
    data = report.algebraic_data(row, curve)
    assert data["shift"] is None and data["scale"] == "1"
    assert data["pole"] == "infinity"
    assert data["model"].startswith("u^6 + 2*u^5")
    assert data["image_points"] == [["-1", "0"], ["0", "-1"], "infinity"]
