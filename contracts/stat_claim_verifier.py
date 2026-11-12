# { "Depends": "py-genlayer:1jb45aa8ynh2a9c9xn3b7qqh8sm5q93hwfp7jqmwsfhh8jpz09h6" }


from genlayer import *

import json
import typing
from dataclasses import dataclass


# ---------------------------------------------------------------- constants

# Two-sided critical values of the standard normal, squared, x1000.
# (alpha in basis points -> z^2 * 1000), ordered by increasing alpha.
_Z_CRIT_TABLE_SQ_MILLI = [
    (10, 10827),    # 0.001 -> 3.290527^2
    (20, 9549),     # 0.002 -> 3.090232^2
    (50, 7879),     # 0.005 -> 2.807034^2
    (100, 6634),    # 0.01  -> 2.575829^2
    (200, 5411),    # 0.02  -> 2.326348^2
    (500, 3841),    # 0.05  -> 1.959964^2
    (1000, 2705),   # 0.10  -> 1.644854^2
]

_Z_ALPHA_MICRO = 1959964      # 1.959964 * 1e6, used for the 95% interval
_Z_ALPHA_BETA_MICRO = 841621  # 0.841621 * 1e6, 80% power
_MIN_GROUP_N = 30
_MAX_SUBGROUPS = 50

STATE_OPEN = "OPEN"
STATE_ATTESTED = "ATTESTED"
STATE_EVALUATED = "EVALUATED"
STATE_FINAL = "FINAL"


# ------------------------------------------------------------- pure helpers

def _isqrt(n: int) -> int:
    """Integer square root (Newton). Deterministic across nodes."""
    if n <= 0:
        return 0
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x


def _z_crit_sq_milli(alpha_bp: int, subgroups: int) -> int:
    """Critical z^2 x1000 for a Bonferroni-corrected two-sided alpha.

    The table is ordered by increasing alpha; we want the largest alpha we are
    allowed to use, i.e. the last row whose alpha does not exceed the effective
    level. That keeps the test conservative under multiple comparisons.
    """
    effective = alpha_bp // max(1, subgroups)
    chosen = _Z_CRIT_TABLE_SQ_MILLI[0][1]
    for bound, value in _Z_CRIT_TABLE_SQ_MILLI:
        if bound <= effective:
            chosen = value
    return chosen


def _split_csv_line(line: str):
    """Minimal RFC4180-ish splitter: honours double quotes and "" escapes."""
    out = []
    cur = []
    in_quotes = False
    i = 0
    n = len(line)
    while i < n:
        ch = line[i]
        if in_quotes:
            if ch == '"':
                if i + 1 < n and line[i + 1] == '"':
                    cur.append('"')
                    i += 2
                    continue
                in_quotes = False
                i += 1
                continue
            cur.append(ch)
            i += 1
            continue
        if ch == '"':
            in_quotes = True
            i += 1
            continue
        if ch == ',':
            out.append("".join(cur))
            cur = []
            i += 1
            continue
        cur.append(ch)
        i += 1
    out.append("".join(cur))
    cells = []
    for cell in out:
        cells.append(cell.strip())
    return cells


def _norm(value: str) -> str:
    return value.strip().strip('"').strip("'").lower()


def _resolve_column(header, wanted: str) -> int:
    """Exact match, then containment, then case-insensitive. -1 when absent."""
    target = _norm(wanted)
    if target == "":
        return -1
    for idx, name in enumerate(header):
        if _norm(name) == target:
            return idx
    for idx, name in enumerate(header):
        if target in _norm(name):
            return idx
    return -1


def _has_usable_header(row) -> bool:
    if len(row) < 2:
        return False
    for cell in row:
        if _norm(cell) == "":
            return False
    return True


def _count_groups(
    rows,
    group_idx: int,
    outcome_idx: int,
    group_a: str,
    group_b: str,
    success_value: str,
):
    """Deterministic tally over parsed rows -> (n_a, x_a, n_b, x_b, skipped)."""
    ga = _norm(group_a)
    gb = _norm(group_b)
    sv = _norm(success_value)
    n_a = 0
    x_a = 0
    n_b = 0
    x_b = 0
    skipped = 0
    width = max(group_idx, outcome_idx) + 1
    for row in rows:
        if len(row) < width:
            skipped += 1
            continue
        g = _norm(row[group_idx])
        o = _norm(row[outcome_idx])
        if g == ga:
            n_a += 1
            if o == sv:
                x_a += 1
        elif g == gb:
            n_b += 1
            if o == sv:
                x_b += 1
        else:
            skipped += 1
    return n_a, x_a, n_b, x_b, skipped


# ------------------------------------------------------------------ storage

@allow_storage
@dataclass
class Claim:
    claim_id: str
    submitter: Address
    report_url: str
    dataset_url: str
    claim_text: str
    group_field: str
    group_a: str
    group_b: str
    outcome_field: str
    success_value: str
    alpha_bp: u32
    claimed_direction: i32          # 1 = A better, -1 = B better, 0 = unspecified
    min_effect_bp: i32
    subgroups: u32
    state: str
    # extraction results
    group_idx: i32
    outcome_idx: i32
    n_a: u32
    x_a: u32
    n_b: u32
    x_b: u32
    skipped_rows: u32
    applicable: str
    extraction_mode: str
    # deterministic statistics
    diff_bp: i32
    z_sq_milli: bigint
    z_crit_sq_milli: bigint
    ci_lo_bp: i32
    ci_hi_bp: i32
    mdes_bp: i32
    # verdict + lifecycle
    verdict: str
    appeals: u32
    attested_at: str
    evaluated_at: str


class StatClaimVerifier(gl.Contract):
    claims: TreeMap[str, Claim]
    next_id: u32

    def __init__(self):
        self.next_id = u32(0)

    # ----------------------------------------------------------- registration

    @gl.public.write
    def register_claim(
        self,
        report_url: str,
        dataset_url: str,
        claim_text: str,
        group_field: str,
        group_a: str,
        group_b: str,
        outcome_field: str,
        success_value: str,
        alpha_bp: u32,
        claimed_direction: i32,
        min_effect_bp: i32,
        subgroups: u32,
    ) -> str:
        a_bp = int(alpha_bp)
        direction = int(claimed_direction)
        min_eff = int(min_effect_bp)
        sg = int(subgroups)

        if not report_url.startswith("http"):
            raise gl.vm.UserError("report_url must be an http(s) URL")
        if not dataset_url.startswith("http"):
            raise gl.vm.UserError("dataset_url must be an http(s) URL")
        if len(claim_text.strip()) < 12:
            raise gl.vm.UserError("claim_text too short to be a testable claim")
        if _norm(group_a) == _norm(group_b):
            raise gl.vm.UserError("group_a and group_b must differ")
        if _norm(success_value) == "":
            raise gl.vm.UserError("success_value must not be empty")
        if a_bp <= 0 or a_bp > 2000:
            raise gl.vm.UserError("alpha_bp must be within (0, 2000]")
        if direction not in (-1, 0, 1):
            raise gl.vm.UserError("claimed_direction must be -1, 0 or 1")
        if min_eff < 0:
            raise gl.vm.UserError("min_effect_bp must be >= 0")
        if sg < 1 or sg > _MAX_SUBGROUPS:
            raise gl.vm.UserError("subgroups must be within [1, 50]")

        cid = str(int(self.next_id) + 1)
        self.claims[cid] = Claim(
            claim_id=cid,
            submitter=gl.message.sender_address,
            report_url=report_url,
            dataset_url=dataset_url,
            claim_text=claim_text,
            group_field=group_field,
            group_a=group_a,
            group_b=group_b,
            outcome_field=outcome_field,
            success_value=success_value,
            alpha_bp=u32(a_bp),
            claimed_direction=i32(direction),
            min_effect_bp=i32(min_eff),
            subgroups=u32(sg),
            state=STATE_OPEN,
            group_idx=i32(-1),
            outcome_idx=i32(-1),
            n_a=u32(0),
            x_a=u32(0),
            n_b=u32(0),
            x_b=u32(0),
            skipped_rows=u32(0),
            applicable="",
            extraction_mode="",
            diff_bp=i32(0),
            z_sq_milli=bigint(0),
            z_crit_sq_milli=bigint(0),
            ci_lo_bp=i32(0),
            ci_hi_bp=i32(0),
            mdes_bp=i32(0),
            verdict="",
            appeals=u32(0),
            attested_at="",
            evaluated_at="",
        )
        self.next_id = u32(int(self.next_id) + 1)
        return cid

    # -------------------------------------------------------------- extraction

    @gl.public.write
    def attest(self, claim_id: str) -> str:
        claim = self.claims.get(claim_id)
        if claim is None:
            raise gl.vm.UserError("unknown claim")
        if claim.state not in (STATE_OPEN, STATE_ATTESTED):
            raise gl.vm.UserError("claim already evaluated; only OPEN/ATTESTED can be attested")

        dataset_url = claim.dataset_url
        group_field = claim.group_field
        outcome_field = claim.outcome_field
        group_a = claim.group_a
        group_b = claim.group_b
        success_value = claim.success_value
        claim_text = claim.claim_text

        def leader_fn():
            response = gl.nondet.web.get(dataset_url)
            status = getattr(response, "status_code", 200)
            if status is not None and int(status) >= 400:
                raise gl.vm.UserError(f"[EXTERNAL] dataset fetch failed: {status}")
            body = response.body
            text = body.decode("utf-8") if isinstance(body, (bytes, bytearray)) else str(body)
            if len(text.strip()) == 0:
                raise gl.vm.UserError("[EXPECTED] dataset is empty")

            lines = [ln for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n") if ln.strip() != ""]
            if len(lines) < 3:
                raise gl.vm.UserError("[EXPECTED] dataset needs a header row plus at least two data rows")

            header = _split_csv_line(lines[0])
            if not _has_usable_header(header):
                raise gl.vm.UserError("[EXPECTED] dataset header row is not usable")

            g_idx = _resolve_column(header, group_field)
            o_idx = _resolve_column(header, outcome_field)
            mode = "deterministic"

            if g_idx < 0 or o_idx < 0:
                # Language is genuinely required here: the headers do not match
                # the field names by any mechanical rule.
                listing = " | ".join(f"{i}:{h}" for i, h in enumerate(header))
                prompt = f"""A dataset has this header row (index:name):

{listing}

We need two columns:
- the column holding the experimental group / cohort label (we look for "{group_field}")
- the column holding the binary outcome we count successes in (we look for "{outcome_field}")

Return JSON only: {{"group_idx": <int or -1>, "outcome_idx": <int or -1>}}"""
                picked = gl.nondet.exec_prompt(prompt, response_format="json")
                if not isinstance(picked, dict):
                    raise gl.vm.UserError("[LLM_ERROR] column resolution did not return JSON")
                g_idx = int(picked.get("group_idx", -1))
                o_idx = int(picked.get("outcome_idx", -1))
                mode = "semantic"

            if g_idx < 0 or o_idx < 0:
                raise gl.vm.UserError("[EXPECTED] declared columns are not present in the dataset")
            if g_idx >= len(header) or o_idx >= len(header):
                raise gl.vm.UserError("[EXPECTED] resolved column index out of range")
            if g_idx == o_idx:
                raise gl.vm.UserError("[EXPECTED] group and outcome resolve to the same column")

            rows = [_split_csv_line(ln) for ln in lines[1:]]
            n_a, x_a, n_b, x_b, skipped = _count_groups(
                rows, g_idx, o_idx, group_a, group_b, success_value
            )
            if n_a == 0 or n_b == 0:
                raise gl.vm.UserError("[EXPECTED] one of the declared groups has no rows")

            # Applicability is judgement, not lookup: does this dataset actually
            # answer this claim, on its own terms?
            sample = "\n".join(lines[:6])[:1500]
            applicability_prompt = f"""Claim under audit:
"{claim_text}"

Dataset columns: {", ".join(header)}
Rows in group "{group_a}": {n_a}
Rows in group "{group_b}": {n_b}
Rows matching neither group (excluded): {skipped}
Success value counted in "{outcome_field}": "{success_value}"
First rows:
{sample}

Decide whether this dataset can bear on this claim as written.
"yes"      = rows, groups and outcome line up with what the claim asserts
"partial"  = it can be computed, but a caveat materially limits what it proves
"no"       = the dataset does not measure what the claim asserts

Return JSON only: {{"applicable": "yes"|"partial"|"no", "reason": "<one sentence>"}}"""
            verdict_obj = gl.nondet.exec_prompt(applicability_prompt, response_format="json")
            if not isinstance(verdict_obj, dict):
                raise gl.vm.UserError("[LLM_ERROR] applicability check did not return JSON")
            applicable = _norm(str(verdict_obj.get("applicable", "")))
            if applicable not in ("yes", "partial", "no"):
                raise gl.vm.UserError("[LLM_ERROR] applicability verdict out of domain")

            return {
                "group_idx": g_idx,
                "outcome_idx": o_idx,
                "n_a": n_a,
                "x_a": x_a,
                "n_b": n_b,
                "x_b": x_b,
                "skipped_rows": skipped,
                "applicable": applicable,
                "extraction_mode": mode,
            }

        def validator_fn(leaders_res: gl.vm.Result) -> bool:
            if not isinstance(leaders_res, gl.vm.Return):
                return False
            leader = leaders_res.calldata
            if not isinstance(leader, dict):
                return False
            try:
                mine = leader_fn()
            except Exception:
                return False
            # Counts and column bindings must agree exactly: they are the
            # deterministic part, computed from the same source.
            for field in ("group_idx", "outcome_idx", "n_a", "x_a", "n_b", "x_b", "applicable"):
                if leader.get(field) != mine.get(field):
                    return False
            return True

        result = gl.vm.run_nondet_unsafe(leader_fn, validator_fn)

        claim.group_idx = i32(int(result["group_idx"]))
        claim.outcome_idx = i32(int(result["outcome_idx"]))
        claim.n_a = u32(int(result["n_a"]))
        claim.x_a = u32(int(result["x_a"]))
        claim.n_b = u32(int(result["n_b"]))
        claim.x_b = u32(int(result["x_b"]))
        claim.skipped_rows = u32(int(result["skipped_rows"]))
        claim.applicable = str(result["applicable"])
        claim.extraction_mode = str(result["extraction_mode"])
        claim.state = STATE_ATTESTED
        import datetime
        claim.attested_at = datetime.datetime.now().isoformat()
        self.claims[claim_id] = claim
        return json.dumps(
            {
                "n_a": int(claim.n_a),
                "x_a": int(claim.x_a),
                "n_b": int(claim.n_b),
                "x_b": int(claim.x_b),
                "applicable": claim.applicable,
                "mode": claim.extraction_mode,
            },
            sort_keys=True,
        )

    # ------------------------------------------------------------ statistics

    @gl.public.write
    def evaluate(self, claim_id: str) -> str:
        claim = self.claims.get(claim_id)
        if claim is None:
            raise gl.vm.UserError("unknown claim")
        if claim.state != STATE_ATTESTED:
            raise gl.vm.UserError("claim must be ATTESTED before it can be evaluated")

        n_a = int(claim.n_a)
        x_a = int(claim.x_a)
        n_b = int(claim.n_b)
        x_b = int(claim.x_b)
        if n_a <= 0 or n_b <= 0:
            raise gl.vm.UserError("degenerate group sizes")
        if x_a > n_a or x_b > n_b:
            raise gl.vm.UserError("success count exceeds group size")

        big_n = n_a + n_b
        big_x = x_a + x_b

        # Difference in proportion, in basis points, exact.
        raw = (x_a * n_b - x_b * n_a)          # sign and magnitude of the gap
        denom = n_a * n_b
        diff_bp = (raw * 10000) // denom

        # z^2 = (x_a*n_b - x_b*n_a)^2 * N / (n_a*n_b*X*(N-X)), saturated cases guard.
        if big_x == 0 or big_x == big_n:
            z_sq_milli = bigint(0)
            se_micro = 0
        else:
            num = raw * raw * big_n
            den = denom * big_x * (big_n - big_x)
            z_sq_milli = bigint((num * 1000) // den)
            # se^2 = X(N-X) / (N*n_a*n_b)
            se_sq_micro = (big_x * (big_n - big_x) * 1000000000000) // (big_n * denom)
            se_micro = _isqrt(se_sq_micro)

        crit = int(_z_crit_sq_milli(int(claim.alpha_bp), int(claim.subgroups)))
        half_bp = (_Z_ALPHA_MICRO * se_micro) // 100000000
        mdes_bp = ((_Z_ALPHA_MICRO + _Z_ALPHA_BETA_MICRO) * se_micro) // 100000000

        claim.diff_bp = i32(diff_bp)
        claim.z_sq_milli = z_sq_milli
        claim.z_crit_sq_milli = bigint(crit)
        claim.ci_lo_bp = i32(diff_bp - half_bp)
        claim.ci_hi_bp = i32(diff_bp + half_bp)
        claim.mdes_bp = i32(mdes_bp)
        claim.state = STATE_EVALUATED
        import datetime
        claim.evaluated_at = datetime.datetime.now().isoformat()
        self.claims[claim_id] = claim

        return json.dumps(
            {
                "diff_bp": int(claim.diff_bp),
                "ci_bp": [int(claim.ci_lo_bp), int(claim.ci_hi_bp)],
                "z_sq_milli": int(claim.z_sq_milli),
                "z_crit_sq_milli": int(claim.z_crit_sq_milli),
                "mdes_bp": int(claim.mdes_bp),
                "alpha_effective_bp": int(claim.alpha_bp) // int(claim.subgroups),
            },
            sort_keys=True,
        )

    @gl.public.write
    def finalize(self, claim_id: str) -> str:
        claim = self.claims.get(claim_id)
        if claim is None:
            raise gl.vm.UserError("unknown claim")
        if claim.state != STATE_EVALUATED:
            raise gl.vm.UserError("claim must be EVALUATED before it can be finalized")

        n_a = int(claim.n_a)
        n_b = int(claim.n_b)
        diff = int(claim.diff_bp)
        z_sq = int(claim.z_sq_milli)
        crit = int(claim.z_crit_sq_milli)
        mdes = int(claim.mdes_bp)
        min_eff = int(claim.min_effect_bp)
        claimed_dir = int(claim.claimed_direction)

        if claim.applicable != "yes":
            verdict = "INSUFFICIENT"
        elif n_a < _MIN_GROUP_N or n_b < _MIN_GROUP_N:
            verdict = "UNDERPOWERED"
        elif z_sq >= crit:
            direction = 1 if diff > 0 else (-1 if diff < 0 else 0)
            if claimed_dir != 0 and direction != claimed_dir:
                verdict = "DIRECTION_MISMATCH"
            elif min_eff > 0 and (diff if diff >= 0 else -diff) < min_eff:
                verdict = "SIGNIFICANT_BUT_SMALL"
            else:
                verdict = "SUPPORTED"
        else:
            # Not significant. The question is whether the interval is tight
            # enough to rule out an effect the claim would care about. If the
            # whole interval sits below the materiality threshold, the study
            # had the precision to look and found nothing: REFUTED. If the
            # interval straddles the threshold, the study simply could not tell.
            threshold = min_eff if min_eff > 0 else 100
            if claimed_dir >= 0:
                if int(claim.ci_hi_bp) < threshold:
                    verdict = "REFUTED"
                else:
                    verdict = "UNDERPOWERED"
            else:
                if int(claim.ci_lo_bp) > -threshold:
                    verdict = "REFUTED"
                else:
                    verdict = "UNDERPOWERED"

        claim.verdict = verdict
        claim.state = STATE_FINAL
        self.claims[claim_id] = claim
        return verdict

    # --------------------------------------------------------------- appeals

    @gl.public.write
    def appeal(self, claim_id: str, note: str) -> str:
        claim = self.claims.get(claim_id)
        if claim is None:
            raise gl.vm.UserError("unknown claim")
        if int(claim.appeals) >= 2:
            raise gl.vm.UserError("appeal budget exhausted")
        if claim.state == STATE_OPEN:
            raise gl.vm.UserError("nothing to appeal before attestation")
        if len(note.strip()) < 8:
            raise gl.vm.UserError("appeal note must state what is wrong with the extraction")

        claim.appeals = u32(int(claim.appeals) + 1)
        claim.state = STATE_ATTESTED
        claim.verdict = ""
        claim.z_sq_milli = bigint(0)
        claim.z_crit_sq_milli = bigint(0)
        claim.evaluated_at = ""
        self.claims[claim_id] = claim
        return f"appeal {int(claim.appeals)}/2 accepted; re-attest and re-evaluate"

    # ------------------------------------------------------------------ views

    @gl.public.view
    def get_claim(self, claim_id: str) -> str:
        claim = self.claims.get(claim_id)
        if claim is None:
            return "{}"
        return json.dumps(
            {
                "claim_id": claim.claim_id,
                "submitter": str(claim.submitter),
                "dataset_url": claim.dataset_url,
                "claim_text": claim.claim_text,
                "state": claim.state,
                "counts": {
                    "n_a": int(claim.n_a),
                    "x_a": int(claim.x_a),
                    "n_b": int(claim.n_b),
                    "x_b": int(claim.x_b),
                    "skipped_rows": int(claim.skipped_rows),
                },
                "applicable": claim.applicable,
                "extraction_mode": claim.extraction_mode,
                "stats": {
                    "diff_bp": int(claim.diff_bp),
                    "ci_bp": [int(claim.ci_lo_bp), int(claim.ci_hi_bp)],
                    "z_sq_milli": int(claim.z_sq_milli),
                    "z_crit_sq_milli": int(claim.z_crit_sq_milli),
                    "mdes_bp": int(claim.mdes_bp),
                },
                "verdict": claim.verdict,
                "appeals": int(claim.appeals),
                "attested_at": claim.attested_at,
                "evaluated_at": claim.evaluated_at,
            },
            sort_keys=True,
        )

    @gl.public.view
    def count_claims(self) -> u32:
        return self.next_id

    @gl.public.view
    def is_supported(self, claim_id: str) -> bool:
        """Convenience for other contracts composing on this one."""
        claim = self.claims.get(claim_id)
        if claim is None:
            return False
        return claim.verdict == "SUPPORTED"


# Notes (kept out of the runner header on purpose:
# GenVM concatenates every leading comment line into the runner
# descriptor, so a multi-line header is not valid JSON).
# StatClaimVerifier — an Intelligent Contract that adjudicates statistical claims.
#
# A claimant points at a public dataset and at a claim made about it
# ("treatment converts 12% better than control, significant at alpha=0.05").
# The contract extracts the counts, runs the two-proportion test itself in
# exact integer arithmetic, and issues a graded verdict instead of a boolean.
#
# Design split that matters for consensus:
#   - Everything that CAN be deterministic IS deterministic: CSV parsing,
#     column resolution when headers are well-formed, the hypothesis test,
#     the confidence interval, the power floor, the decision tree.
#   - Non-determinism is confined to two places where language is genuinely
#     involved: resolving non-standard headers, and judging whether the
#     dataset is applicable to the claim at all. Both are re-derived
#     independently by validators and compared field by field.
