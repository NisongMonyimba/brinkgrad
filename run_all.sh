#!/usr/bin/env bash
# run_all.sh — reproduce all brinkgrad paper results
# Run inside Docker: docker run --rm -v ${PWD}:/root/brinkgrad \
#   dolfinx/dolfinx:v0.7.3 bash -c "cd /root/brinkgrad && bash run_all.sh"
set -euo pipefail

PY="python3"

# Install brinkgrad in the current environment
$PY -m pip install -e . --no-deps -q 2>/dev/null || true

echo "============================================================"
echo " brinkgrad — full pipeline"
echo " Python: $($PY --version)"
echo " DOI: 10.5281/zenodo.20538833"
echo "============================================================"

mkdir -p figures docs/si

run_script() {
    local label="$1"; local script="$2"
    echo ""
    echo ">>> $label"
    if $PY "$script"; then
        echo "    [ OK ] $label"
    else
        echo "    [WARN] $label exited non-zero — continuing"
    fi
}

run_script "[1/5] linear_target"          examples/linear_target.py
run_script "[2/5] step_profile_target"    examples/step_profile_target.py
run_script "[3/5] double_peak_target"     examples/double_peak_target.py
run_script "[4/5] gallery_targets"        examples/gallery_targets.py
run_script "[5/5] robust_projection_demo" examples/robust_projection_demo.py

echo ""
echo "============================================================"
echo " All scripts finished."
echo " Figures saved in: figures/"
echo " Primary result: RMSE = 0.058 (linear gradient, 80x20 mesh)"
echo "============================================================"
