#!/usr/bin/env bash
set -euo pipefail
DRY_RUN=0; NO_COMMIT=0
for arg in "$@"; do
  [[ $arg == "--dry-run"   ]] && DRY_RUN=1
  [[ $arg == "--no-commit" ]] && NO_COMMIT=1
done

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

echo "════════════════════════════════════════════════════════"
echo "  micrograd fix-and-rerun  $(date '+%Y-%m-%d %H:%M:%S')"
echo "════════════════════════════════════════════════════════"

if ! python3 -c "import dolfinx" 2>/dev/null; then
  echo "[ERROR] dolfinx not importable. Run: conda activate fenicsx"; exit 1
fi
echo "[OK] dolfinx importable"

F1="micrograd/gradient_optimizer.py"
F2="micrograd/utilities.py"
F3="micrograd/binary_validation.py"

echo ""; echo "── Patch 1/3: $F1 ──"
cp "$F1" "${F1}.bak"
python3 - "$F1" << 'PY'
import sys, pathlib
path = pathlib.Path(sys.argv[1])
src  = path.read_text()
patches = [
    ("w_f=1e-7, w_c=5e4",
     "w_f=1e-3, w_c=5e1"),
    ("def run(self, max_iter=80, beta_continuation=(1,2,4,8,16,32,64), move=0.2,",
     "def run(self, max_iter=400, beta_continuation=(1,2,4,8,16), move=0.05,"),
    ("n_betas = len(beta_continuation); iters_per_beta = max_iter // n_betas",
     "n_betas = len(beta_continuation); iters_per_beta = min(80, max_iter // n_betas)"),
]
changed = []
for old, new in patches:
    if old in src:
        src = src.replace(old, new, 1); changed.append(old[:60])
    else:
        print(f"  [WARN] not found (already patched?): {old[:70]}")
path.write_text(src)
for c in changed: print(f"  patched: {c!r}")
print(f"  {len(changed)}/{len(patches)} patches applied")
PY

echo ""; echo "── Patch 2/3: $F2 ──"
cp "$F2" "${F2}.bak"
python3 - "$F2" << 'PY'
import sys, pathlib
path = pathlib.Path(sys.argv[1])
src  = path.read_text()
old = "def alpha(r): return alpha_min + (alpha_max-alpha_min)*(1.0-r)/(1.0+r)"
new = (
    "def alpha(r):\n"
    "    # floor prevents RMSE=inf when binary_validation sets alpha_min=0\n"
    "    _raw = alpha_min + (alpha_max - alpha_min) * (1.0 - r) / (1.0 + r)\n"
    "    return ufl.max_value(_raw, 1e-4)"
)
if old in src:
    src = src.replace(old, new, 1); path.write_text(src)
    print("  patched: alpha() floor added")
else:
    print("  [WARN] alpha() line not found — already patched or source changed.")
PY

echo ""; echo "── Patch 3/3: $F3 ──"
cp "$F3" "${F3}.bak"
python3 - "$F3" << 'PY'
import sys, pathlib
path = pathlib.Path(sys.argv[1])
src  = path.read_text()
old = (
    "    opt = GradientGeneratorOptimizer(\n"
    "        Lx=2000e-6, Ly=500e-6, nx=80, ny=20,\n"
    "        target_expr=lambda x: x[1]/500e-6,\n"
    "        w_f=1e-7, w_c=5e4, V_star=0.5\n"
    "    )\n"
    "    opt.run(max_iter=400,\n"
    "            beta_continuation=[1,2,4,8,16,32,64],\n"
    "            move=0.2)"
)
new = (
    "    opt = GradientGeneratorOptimizer(\n"
    "        Lx=2000e-6, Ly=500e-6, nx=80, ny=20,\n"
    "        target_expr=lambda x: x[1]/500e-6,\n"
    "        w_f=1e-3, w_c=5e1, V_star=0.5\n"
    "    )\n"
    "    opt.run(max_iter=400,\n"
    "            beta_continuation=[1,2,4,8,16],\n"
    "            move=0.05)"
)
if old in src:
    src = src.replace(old, new, 1); path.write_text(src)
    print("  patched: __main__ hyperparams corrected")
else:
    print("  [WARN] __main__ block not matched exactly — check binary_validation.py manually.")
PY

echo ""; echo "── Diffs ──"
for f in "$F1" "$F2" "$F3"; do
  echo "--- $f ---"; diff "${f}.bak" "$f" || true; echo
done

[[ $DRY_RUN -eq 1 ]] && {
  echo "[DRY RUN] restoring originals."
  cp "${F1}.bak" "$F1"; cp "${F2}.bak" "$F2"; cp "${F3}.bak" "$F3"; exit 0
}

echo "── Running optimisation + binary validation (400 iters, ~5-15 min) ──"
python3 - << 'PY'
import sys, json, pathlib, numpy as np
sys.path.insert(0, '.')
from micrograd.binary_validation import run_binary_validation
from micrograd import GradientGeneratorOptimizer
from micrograd.utilities import helmholtz_filter, heaviside_projection

opt = GradientGeneratorOptimizer(
    Lx=2000e-6, Ly=500e-6, nx=80, ny=20,
    target_expr=lambda x: x[1]/500e-6,
    w_f=1e-3, w_c=5e1, V_star=0.5,
)
opt.run(max_iter=400, beta_continuation=[1,2,4,8,16], move=0.05)

helmholtz_filter(opt.rho, opt.rho_filt, opt.V_rho, opt.r_filter)
heaviside_projection(opt.rho_filt, opt.rho_phys, beta=16)

rho = opt.rho_phys.x.array.copy()
gray  = float(np.mean((rho > 0.05) & (rho < 0.95)))
fvol  = float(np.mean(rho > 0.5))
print(f"rho_phys: gray={gray:.3f}  fluid-vol={fvol:.3f}  range=[{rho.min():.3f},{rho.max():.3f}]")
if gray > 0.10:
    print(f"[WARN] gray={gray:.1%} still high — consider more iterations.")

results = run_binary_validation(opt, verbose=True)
results.update(gray_frac=gray, fluid_vol=fvol)
pathlib.Path("/tmp/binary_validation.json").write_text(json.dumps(results, indent=2))

import math
g  = results["rmse_gray"]
bs = results["rmse_binary_stokes"]
dr = (bs-g)/g*100 if math.isfinite(bs) and g>0 else float("nan")
print(f"\nVERDICT: {'VALIDATED <15%' if math.isfinite(dr) and abs(dr)<15 else f'GAP={dr:+.1f}%' if math.isfinite(dr) else 'Binary Stokes still inf'}")
PY

echo ""; echo "── Updating manuscript/macros.tex ──"
python3 - << 'PY'
import json, re, math, pathlib
p = pathlib.Path("/tmp/binary_validation.json")
if not p.exists(): print("  [SKIP] no JSON"); exit()
r = json.loads(p.read_text())
def fmt(v, d=4): return r"\infty" if not math.isfinite(v) else f"{v:.{d}f}"
g=r["rmse_gray"]; bb=r["rmse_binary_brinkman"]; bs=r["rmse_binary_stokes"]
updates = [
    ("rmseTopOpt",         fmt(g,3)),
    ("rmseTopOptPP",       fmt(g*100,1)),
    ("rmseBinaryBrinkman", fmt(bb,3)),
    ("rmseBinaryStokes",   fmt(bs,3)),
    ("grayFrac",           f"{r.get('gray_frac',0):.3f}"),
    ("fluidVol",           f"{r.get('fluid_vol',0):.3f}"),
]
mac = pathlib.Path("manuscript/macros.tex")
if not mac.exists(): print(f"  [SKIP] {mac} not found"); exit()
src = mac.read_text(encoding="utf-8")
for name, val in updates:
    pat = r"\\newcommand\{\\"+name+r"\}\{[^}]+\}"
    rep = f"\\\\newcommand{{\\\\{name}}}{{{val}}}"
    src = re.sub(pat, rep, src) if re.search(pat, src) \
          else src.rstrip()+f"\n\\newcommand{{\\{name}}}{{{val}}}\n"
mac.write_text(src, encoding="utf-8")
print(f"  macros.tex updated: gray={fmt(g,4)}  bb={fmt(bb,4)}  bs={fmt(bs,4)}")
PY

[[ $NO_COMMIT -eq 1 ]] && { echo "[--no-commit] done."; exit 0; }

echo ""; echo "── git commit + push ──"
VERDICT=$(python3 -c "
import json,math
r=json.load(open('/tmp/binary_validation.json'))
g=r.get('rmse_gray',0); bs=r.get('rmse_binary_stokes',float('nan'))
dr=(bs-g)/g*100 if math.isfinite(bs) and g>0 else float('nan')
print('VALIDATED' if math.isfinite(dr) and abs(dr)<15 else f'GAP={dr:+.1f}%' if math.isfinite(dr) else 'inf-pending')
" 2>/dev/null || echo "unknown")

git add micrograd/gradient_optimizer.py micrograd/utilities.py \
        micrograd/binary_validation.py manuscript/macros.tex 2>/dev/null || true

git diff --cached --quiet && { echo "  Nothing to commit."; exit 0; }

git commit -m "fix: rescale weights, alpha floor, beta schedule; verdict=${VERDICT}

- w_f: 1e-7→1e-3, w_c: 5e4→5e1 (ratio 5e11→5e4)
- move: 0.2→0.05
- iters_per_beta: capped at 80 (was 28 at 200 iters/7 levels)
- alpha(): ufl.max_value floor=1e-4 (fixes RMSE=inf in Binary Stokes)
- binary_validation __main__: corrected hyperparams
"
git push origin main && echo "  Pushed."

echo ""
echo "════════════════════════════════════════════════════════"
echo "  Done.  /tmp/binary_validation.json"
echo "════════════════════════════════════════════════════════"
