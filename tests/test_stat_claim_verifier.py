"""Direct-mode tests for StatClaimVerifier.

The statistics are checked against values computed by hand, because the whole
point of the contract is that its numbers are reproducible.
"""

import json

DATASET_URL = "https://data.example.org/ab_test.csv"
REPORT_URL = "https://research.example.org/report-42"


def make_csv(n_a, x_a, n_b, x_b, header="user_id,cohort,converted"):
    rows = [header]
    for i in range(n_a):
        rows.append(f"t{i},treatment,{'yes' if i < x_a else 'no'}")
    for i in range(n_b):
        rows.append(f"c{i},control,{'yes' if i < x_b else 'no'}")
    return "\n".join(rows)


def deploy_helpers_removed():
    pass


def register(contract, **over):
    kwargs = dict(
        report_url=REPORT_URL,
        dataset_url=DATASET_URL,
        claim_text="treatment converts 12% better than control, significant at alpha 0.05",
        group_field="cohort",
        group_a="treatment",
        group_b="control",
        outcome_field="converted",
        success_value="yes",
        alpha_bp=500,
        claimed_direction=1,
        min_effect_bp=100,
        subgroups=1,
    )
    kwargs.update(over)
    return contract.register_claim(**kwargs)


def test_register_validates_inputs(direct_deploy, direct_vm):
    c = direct_deploy("contracts/stat_claim_verifier.py")
    assert register(c) == "1"

    import pytest
    with direct_vm.expect_revert("group_a and group_b must differ"):
        register(c, group_b="treatment")
    with direct_vm.expect_revert("alpha_bp must be within"):
        register(c, alpha_bp=0)
    with direct_vm.expect_revert("claim_text too short"):
        register(c, claim_text="short")


def test_attest_counts_exactly(direct_deploy, direct_vm):
    csv_text = make_csv(1000, 150, 1000, 100)
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": csv_text})
    direct_vm.mock_llm(
        r"Decide whether this dataset can bear on this claim",
        json.dumps({"applicable": "yes", "reason": "cohort column matches"}),
    )
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    out = json.loads(c.attest(cid))
    assert out["n_a"] == 1000 and out["x_a"] == 150
    assert out["n_b"] == 1000 and out["x_b"] == 100
    assert out["mode"] == "deterministic"
    assert out["applicable"] == "yes"


def test_statistics_match_hand_calculation(direct_deploy, direct_vm):
    # 150/1000 vs 100/1000
    #   diff = +500 bp
    #   z^2 = (50000^2 * 2000) / (1e6 * 250 * 1750) = 11428.6 -> 11428 milli
    #   se  = sqrt(250*1750 / (2000*1e6)) = 0.014790 -> CI half = 290 bp
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)
    stats = json.loads(c.evaluate(cid))
    assert stats["diff_bp"] == 500
    assert stats["z_sq_milli"] == 11428
    assert stats["z_crit_sq_milli"] == 3841
    assert stats["ci_bp"] == [211, 789]
    assert stats["mdes_bp"] == 414


def test_verdict_supported(direct_deploy, direct_vm):
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "SUPPORTED"
    assert c.is_supported(cid) is True


def test_verdict_significant_but_small(direct_deploy, direct_vm):
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c, min_effect_bp=600)
    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "SIGNIFICANT_BUT_SMALL"


def test_verdict_direction_mismatch(direct_deploy, direct_vm):
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c, claimed_direction=-1)
    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "DIRECTION_MISMATCH"


def test_verdict_refuted_when_powered(direct_deploy, direct_vm):
    # 10000/100000 vs 10000/100000: gap is exactly 0, the 95% interval is
    # [-26, +26], which rules out the 100 bp effect the claim asserts.
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(100000, 10000, 100000, 10000)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c, min_effect_bp=100)
    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "REFUTED"


def test_verdict_underpowered_when_small(direct_deploy, direct_vm):
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(40, 12, 40, 10)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "UNDERPOWERED"


def test_verdict_insufficient_when_dataset_does_not_fit(direct_deploy, direct_vm):
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset",
                       json.dumps({"applicable": "partial", "reason": "outcome is a proxy"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "INSUFFICIENT"


def test_appeal_reopens_and_reruns(direct_deploy, direct_vm):
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)
    c.evaluate(cid)
    c.finalize(cid)

    msg = c.appeal(cid, "control group includes returning users")
    assert "1/2" in msg
    state = json.loads(c.get_claim(cid))
    assert state["state"] == "ATTESTED"
    assert state["verdict"] == ""

    c.attest(cid)
    c.evaluate(cid)
    assert c.finalize(cid) == "SUPPORTED"


def test_semantic_column_resolution_when_header_is_odd(direct_deploy, direct_vm):
    csv_text = make_csv(200, 40, 200, 20, header="uid,arm,outcome_flag")
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": csv_text})
    direct_vm.mock_llm(r"A dataset has this header row", json.dumps({"group_idx": 1, "outcome_idx": 2}))
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c, group_field="cohort", outcome_field="converted")
    out = json.loads(c.attest(cid))
    assert out["mode"] == "semantic"
    assert out["n_a"] == 200 and out["x_a"] == 40


def test_rejects_unknown_group(direct_deploy, direct_vm):
    csv_text = "user_id,cohort,converted\nu1,treatment,yes\nu2,treatment,no"
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": csv_text})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    import pytest
    with direct_vm.expect_revert("one of the declared groups has no rows"):
        c.attest(cid)


def test_critical_values_and_bonferroni(direct_deploy, direct_vm):
    # One deployment, several claims: the harness only allows one contract
    # instance per file, so the levels are exercised through distinct claims.
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    c = direct_deploy("contracts/stat_claim_verifier.py")

    cases = [
        (500, 1, 3841),   # alpha 0.05
        (100, 1, 6634),   # alpha 0.01
        (20, 1, 9549),    # alpha 0.002
        (500, 3, 6634),   # 500/3 = 166 bp -> stricter neighbour 0.01
        (10, 10, 10827),  # 1 bp effective -> clamp to the strictest row
    ]
    got = []
    for alpha_bp, subgroups, _expected in cases:
        cid = register(c, alpha_bp=alpha_bp, subgroups=subgroups)
        c.attest(cid)
        c.evaluate(cid)
        got.append(json.loads(c.get_claim(cid))["stats"]["z_crit_sq_milli"])

    assert got == [expected for _a, _s, expected in cases]
