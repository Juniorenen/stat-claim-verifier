"""Consensus behaviour: does the validator actually reject a bad leader?"""

import json

from test_stat_claim_verifier import DATASET_URL, REPORT_URL, make_csv, register


def _mock_common(vm, applicable="yes"):
    vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 150, 1000, 100)})
    vm.mock_llm(r"Decide whether this dataset",
                json.dumps({"applicable": applicable, "reason": "ok"}))


def test_validator_agrees_with_faithful_leader(direct_deploy, direct_vm):
    _mock_common(direct_vm, "yes")
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)
    assert direct_vm.run_validator() is True


def test_validator_rejects_different_applicability(direct_deploy, direct_vm):
    _mock_common(direct_vm, "yes")
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)

    # A validator that independently concludes the dataset does not fit the
    # claim must disagree, even though every count matches.
    direct_vm.clear_mocks()
    _mock_common(direct_vm, "no")
    assert direct_vm.run_validator() is False


def test_validator_rejects_different_counts(direct_deploy, direct_vm):
    _mock_common(direct_vm, "yes")
    c = direct_deploy("contracts/stat_claim_verifier.py")
    cid = register(c)
    c.attest(cid)

    # Same claim, different dataset behind the URL: counts diverge.
    direct_vm.clear_mocks()
    direct_vm.mock_web(DATASET_URL, {"status": 200, "body": make_csv(1000, 149, 1000, 100)})
    direct_vm.mock_llm(r"Decide whether this dataset", json.dumps({"applicable": "yes"}))
    assert direct_vm.run_validator() is False
