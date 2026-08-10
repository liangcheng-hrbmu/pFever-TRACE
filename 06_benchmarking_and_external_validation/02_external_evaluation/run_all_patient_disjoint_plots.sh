#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "$SCRIPT_DIR/plot_patient_disjoint_broad5_withLinearSVC_main.R"
Rscript "$SCRIPT_DIR/plot_patient_disjoint_broad5_withLinearSVC_supp.R"
Rscript "$SCRIPT_DIR/plot_patient_disjoint_222_ours_vs_145_withLinearSVC_final.R"

echo "[DONE] All patient-disjoint figures generated."
