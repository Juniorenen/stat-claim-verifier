#!/bin/bash
GLJS="C:/Users/heade/AppData/Local/Author Software/nvm/installs/v24.21.0/node_modules/genlayer/dist/index.js"
C="0xFA70F90e35Ad790f88f8c331fd8d6A45A922f449"
DATASET="https://raw.githubusercontent.com/Juniorenen/stat-claim-verifier/main/datasets/example_ab.csv"
REPORT="https://raw.githubusercontent.com/Juniorenen/stat-claim-verifier/main/README.md"
CLAIM="Treatment converts at least 10 percentage points better than control, significant at alpha 0.05"
cd /c/Users/heade/genlayer-statclaim || exit 1

step () {
  local label="$1"; shift
  node "$GLJS" "$@" > "e2e_$label.log" 2>&1
  local tx
  tx=$(grep -o "0x[0-9a-f]\{64\}" "e2e_$label.log" | head -1)
  echo "--- $label tx=$tx"
  sleep 9
  node "$GLJS" receipt "$tx" 2>&1 | grep -E "execution_result|payload" | head -2
}

step register write "$C" register_claim --args "$REPORT" "$DATASET" "$CLAIM" cohort treatment control converted yes 500 1 100 1
step attest write "$C" attest --args 1
step evaluate write "$C" evaluate --args 1
step finalize write "$C" finalize --args 1
