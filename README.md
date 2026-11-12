# StatClaimVerifier

An Intelligent Contract that adjudicates **statistical claims** — the kind that
research reports, A/B write-ups, grant applications and token dashboards make
all day — by recomputing the test itself instead of trusting the prose.

A claimant points at a public dataset and at a claim made about it:

> "treatment converts 12% better than control, significant at alpha = 0.05"

The contract extracts the group/outcome counts, runs a two-proportion test in
exact integer arithmetic, and issues a **graded verdict**, not a boolean.

---

## Why this is a primitive and not a one-off demo

Every escrow, bounty, grant and insurance contract on GenLayer eventually needs
the answer to "is this number real?". Today each of them either trusts a
report's prose or bolts on a bespoke "ask an LLM if this looks right" check,
which is exactly the thin wrapper the review guidelines exclude.

This contract answers a narrower and harder question than "is this true?":

**Given this dataset, does it support this specific statistical claim, at this
significance level, with this effect size?**

That is a computation with a right answer, and the contract states what it did
with enough precision that anyone can re-derive it.

---

## Verdicts

The output is one of six states. The distinctions are the point:

| Verdict | Meaning |
|---|---|
| `SUPPORTED` | Significant, in the claimed direction, and the effect clears the materiality floor |
| `SIGNIFICANT_BUT_SMALL` | Real effect, but the claim oversells it — below the stated minimum effect |
| `DIRECTION_MISMATCH` | Significant, but in the opposite direction to what is claimed |
| `REFUTED` | Not significant **and** the confidence interval is tight enough to rule out the claimed effect |
| `UNDERPOWERED` | Not significant, but the interval is too wide to conclude anything — the study could not have seen the effect |
| `INSUFFICIENT` | The dataset does not measure what the claim asserts (judged, not assumed) |

The `REFUTED` / `UNDERPOWERED` split is the one most human reviews get wrong:
absence of significance is not evidence of absence, and only the confidence
interval can tell you which one you are looking at. Both are computed, not
guessed.

---

## Where the non-determinism lives

This is the part that matters for consensus, so it is deliberately small.

Deterministic (every validator must produce the identical value):

- CSV/TSV parsing, including quoted fields and `""` escapes
- Column resolution by exact, then substring, then case-insensitive match
- Group/success counting
- The whole hypothesis test, the interval, the power floor
- Bonferroni correction and the decision tree

Non-deterministic (two places, both re-derived by validators and compared):

1. **Column resolution fallback.** Only when the header cannot be matched
   mechanically (`arm` for cohort, `outcome_flag` for converted). The LLM
   returns two column indices; the counts are then computed deterministically
   from those indices.
2. **Applicability.** Does this dataset bear on this claim at all? This is
   judgement — a proxy metric, a missing control arm, a cohort definition that
   drifts from what the claim says. Validators independently reach their own
   `yes` / `partial` / `no` and must agree exactly on that field.

In both cases the validator re-runs the whole leader function and compares
`group_idx`, `outcome_idx`, `n_a`, `x_a`, `n_b`, `x_b` and `applicable`
field by field. A leader whose counts cannot be reproduced is rejected and the
network rotates.

---

## The statistics, in integers

No floats anywhere: float summation is not reproducible across nodes and this
contract's whole value is that its numbers are.

With `N = n_a + n_b`, `X = x_a + x_b`:

```
z² = (x_a·n_b − x_b·n_a)² · N / (n_a·n_b·X·(N−X))
```

Computed as a ratio of exact integers, stored as `z²·1000`, and compared
against a table of critical values squared — so the significance decision needs
no division and no square root.

The interval half-width needs one square root, done with integer Newton
(`_isqrt`) on the scaled variance, and the power floor uses
`z_alpha + z_beta = 2.801585` for 80% power.

Critical values (`z²·1000`, two-sided):

| alpha | z | z² |
|---|---|---|
| 0.001 | 3.290527 | 10827 |
| 0.005 | 2.807034 | 7879 |
| 0.01 | 2.575829 | 6634 |
| 0.02 | 2.326348 | 5411 |
| 0.05 | 1.959964 | 3841 |
| 0.10 | 1.644854 | 2705 |

Multiple subgroups are handled with a Bonferroni correction: the effective
alpha is `alpha / subgroups`, and the table lookup rounds to the stricter
neighbouring level rather than interpolating.

---

## Contract API

```python
register_claim(
    report_url: str,          # where the claim is published
    dataset_url: str,         # the public CSV/TSV the claim rests on
    claim_text: str,          # the assertion, quoted
    group_field: str,         # e.g. "cohort"
    group_a: str,             # e.g. "treatment"
    group_b: str,             # e.g. "control"
    outcome_field: str,       # e.g. "converted"
    success_value: str,       # e.g. "yes"
    alpha_bp: u32,            # 500 = 0.05
    claimed_direction: i32,   # 1 = A better, -1 = B better, 0 = unspecified
    min_effect_bp: i32,       # smallest effect that matters, in bp
    subgroups: u32,           # how many comparisons the claim makes at once
) -> claim_id

attest(claim_id)     # fetch + extract (non-deterministic, validator-checked)
evaluate(claim_id)   # the test itself (pure, reproducible)
finalize(claim_id)   # apply the decision tree, write the verdict
appeal(claim_id, note)  # reopen a claim, max 2, requires a stated reason

get_claim(claim_id) -> json   # full record, counts, stats, verdict
count_claims() -> u32
is_supported(claim_id) -> bool  # for other contracts composing on this one
```

Lifecycle: `OPEN → ATTESTED → EVALUATED → FINAL`, with `appeal` returning a
claim from `FINAL`/`EVALUATED` to `ATTESTED` and clearing the computed fields.

---

## Example

Dataset at `https://data.example.org/ab_test.csv`:

```
user_id,cohort,converted
t0,treatment,yes
...
```

```python
cid = contract.register_claim(
    report_url="https://research.example.org/report-42",
    dataset_url="https://data.example.org/ab_test.csv",
    claim_text="treatment converts 12% better than control, significant at alpha 0.05",
    group_field="cohort", group_a="treatment", group_b="control",
    outcome_field="converted", success_value="yes",
    alpha_bp=500, claimed_direction=1, min_effect_bp=100, subgroups=1,
)
contract.attest(cid)    # {"n_a": 1000, "x_a": 150, ...}
contract.evaluate(cid)  # {"diff_bp": 500, "z_sq_milli": 11428, "ci_bp": [211, 789], ...}
contract.finalize(cid)  # "SUPPORTED"
```

150/1000 vs 100/1000 → +500 bp, z² = 11.43 against a 3.84 threshold, interval
[211, 789] bp. Hand-checked in the test suite.

---

## Tests

`16 passing`, Direct Mode (in-process, no Docker):

```
tests/test_stat_claim_verifier.py   13 tests
  - input validation on registration
  - exact counts on a clean dataset
  - statistics compared against hand calculation (diff, z², CI, MDES)
  - every verdict branch reachable and asserted
  - critical values per alpha, and Bonferroni rounding to the stricter level
  - semantic column resolution when headers are non-standard
  - appeal reopens and re-runs the pipeline
  - rejection when a declared group has no rows

tests/test_consensus.py              3 tests
  - validator agrees with a faithful leader
  - validator rejects a leader whose applicability judgement differs
  - validator rejects a leader whose counts cannot be reproduced
```

Run:

```bash
python -m venv .venv
.venv/Scripts/pip install genlayer-test    # or .venv/bin/pip on POSIX
.venv/Scripts/python -m pytest tests/ -q
```

`genvm-lint check contracts/stat_claim_verifier.py` → lint passed, validation
passed, 8 methods (3 view, 5 write).

### Note on Windows

`genlayer-test` 0.29.2 calls `os.unlink()` on the temp file it has just pointed
fd 0 at. POSIX allows that; Windows raises `WinError 32`. `tests/conftest.py`
replaces that injector with an equivalent that defers deletion to interpreter
exit. It is test-harness only and does not touch the contract.

### Note on the runner pin

The header pins `py-genlayer:1jb45aa8…`, the hash used by the official
Prediction Market example. `genvm-lint` reports a newer runner available
(`1zr6nqk5…`); the pin was left alone so behaviour matches the documented
example. Bumping it is a one-line change if the submission pipeline prefers
latest.

---

## Limits, stated plainly

- **CSV/TSV only.** The contract refuses HTML and PDF rather than guessing at
  tables. A PDF path would need layout-aware extraction, and pretending a CSV
  parser handles it would be worse than declining.
- **Two-group, binary-outcome tests.** That is a large fraction of real A/B and
  cohort claims, but it is not every statistical claim.
- **Wald interval, normal approximation.** Fine at these sample sizes, weak for
  tiny counts or rates near 0/1; `_MIN_GROUP_N = 30` blocks the worst of it.
- **The dataset URL is trusted to be the right dataset.** Applicability
  judgement catches a mismatched dataset often, but a claimant who links a
  genuine dataset that measures something adjacent will only be caught by that
  judgement step, not by a proof.
- **`alpha_bp` outside the table gets the conservative neighbour**, not an
  interpolated critical value. Documented, deterministic, slightly stricter
  than requested.
