

import numpy as np
import itertools
import pandas as pd
from scipy.optimize import root, brentq

#Parameter grids
N_vals   = list(range(2, 16))                                  # {2, 3, …, 15}
lam_vals = np.logspace(np.log10(0.001), np.log10(100), 20)    # 20 log-spaced
R_vals   = np.logspace(np.log10(10),    np.log10(10000), 20)  # 20 log-spaced

rho_grid = np.linspace(0.05, 0.95, 200)  # belief-space sweep

#Self-contained P2 evaluator
def compute_P2(rho_star, N_val, lam_val, R_val, q0_guess=0.5):
    p = rho_star

    # Guard: denominators that appear in A_X, B_X require p ∉ {0, 1}
    if p <= 0.0 or p >= 1.0:
        return np.nan, np.nan

    # Inner structural term (used by D_Y)
    def _inner(q):
        # Mathematica: Y[q1, r2, p_eq, l, c2, n] has inner term depending on
        return (
            p * (N_val * p**2
                 + N_val * (-1 + (1-p)**N_val) * (-1+q) * p   # r_pub = p
                 - p * (-1+q) * (1 + (-1+(1-p)**N_val) * p)   # r_pub = p
                 + N_val * p * (-1 + (-1+q) * p))              # r_pub = p
            - p**N_val * (q - N_val*q + p*(p + (-2+N_val)*q + p*(-1+q)*p))
        )

    def DY(q):
        """D_Y = inner^2 / ((N-1)^2*(1-p)^4*p^4)"""
        t = _inner(q)
        denom = (N_val - 1)**2 * (1-p)**4 * p**4
        if denom == 0:
            return np.nan
        return t**2 / denom

    def AX(q):
        """A_X coefficient at symmetric benchmark rho = p."""
        rho = p
        denom1 = (N_val - 1) * (p - 1) * p
        denom2 = (N_val - 1) * p
        denom3 = (N_val - 1) * (p - 1)**2
        if denom1 == 0 or denom2 == 0 or denom3 == 0:
            return np.nan
        return (
            q + (-1 + p**(-2+N_val)) * q
            + N_val * (p - p**N_val) * q / denom1
            + (1-q) * (-1 + N_val * (1 - (1-p)**(N_val-1)) * rho / denom2)
            - ((-1+q) * (
                (N_val-1) * p**2
                + p**N_val * (-1 + rho)
                - (N_val-1) * (-1 + (1-p)**N_val) * rho
                - p * (-2 + N_val + N_val*rho)
            )) / denom3
        )

    def BX(q):
        """B_X coefficient."""
        denom = (N_val - 1) * (p - 1) * p**2
        if denom == 0:
            return np.nan
        return (
            (1 - (1-p)**N_val + (-1 - N_val + (1-p)**N_val)*p + p**N_val*(-1+N_val+p)) * q
            / denom
        )

    def DX(q):
        """D_X coefficient at symmetric benchmark."""
        rho = p
        inner_a = (
            p**N_val * (q - N_val*q + p*(p + (-2+N_val)*q + p*(-1+q)*rho))
            + p * (-N_val*p**2 - N_val*(-1+(1-p)**N_val)*(-1+q)*rho
                   + p*(-1+q)*(1+(-1+(1-p)**N_val)*rho)
                   + N_val*p*(1 + rho - q*rho))
        )
        inner_b = (
            p**N_val * (p - N_val*q - p*rho + p*q*(-1+N_val+rho))
            + p * (-(N_val-1)*p**2 + q + rho
                   + (-N_val - (N_val-1)*(1-p)**N_val*(-1+q) - q + N_val*q)*rho
                   + p*(-2 + N_val*(1 + rho - q*rho)))
        )
        denom = (N_val-1)**2 * (p-1)**4 * p**3
        if denom == 0:
            return np.nan
        return inner_a * inner_b / denom

    def X_full(q):
        sj = p * R_val
        ax = AX(q);  bx = BX(q);  dx = DX(q)
        if any(v is None or np.isnan(v) for v in [ax, bx, dx]):
            return np.nan
        return ax + bx*p - 2*lam_val*(1-p)*p/(sj+3)*dx

    def Y_full(q):
        sj = p * R_val
        dy = DY(q)
        if dy is None or np.isnan(dy):
            return np.nan
        return -lam_val*(1-p)*p/(sj+3) * dy

    # U(q) = Y*q^2 + X*q 
    def U(q):
        y = Y_full(q);  x = X_full(q)
        if y is None or x is None or np.isnan(y) or np.isnan(x):
            return np.nan
        return y*q**2 + x*q

    #F(q) = q + X/(2Y)
    def F_scalar(q):
        y = Y_full(q)
        if y is None or y == 0 or np.isnan(y):
            return np.inf
        x = X_full(q)
        if x is None or np.isnan(x):
            return np.inf
        return q + x / (2*y)

    def F_vec(q_arr):
        return [F_scalar(q_arr[0])]

    #Strategy 1: bracket scan + brentq
    q_star   = np.nan
    scan_pts = np.linspace(0.02, 0.98, 50)
    f_vals   = np.array([F_scalar(q) for q in scan_pts])
    finite   = np.isfinite(f_vals)

    bracket = None
    for i in range(len(scan_pts) - 1):
        if finite[i] and finite[i+1] and f_vals[i] * f_vals[i+1] < 0:
            bracket = (scan_pts[i], scan_pts[i+1])
            break

    if bracket is not None:
        try:
            q_star = brentq(F_scalar, bracket[0], bracket[1], xtol=1e-12, rtol=1e-10)
        except Exception:
            pass

    #Strategy 2: hybr fallback (warm-start then fixed guesses)
    if np.isnan(q_star):
        guesses = [q0_guess] + [g for g in [0.5, 0.3, 0.7, 0.1, 0.9] if abs(g - q0_guess) > 0.05]
        for q0 in guesses:
            try:
                sol = root(F_vec, [q0], method='hybr', tol=1e-10)
                if sol.success and 0.0 < sol.x[0] < 1.0:
                    if abs(F_scalar(sol.x[0])) < 1e-7:
                        q_star = sol.x[0]
                        break
            except Exception:
                continue

    # Strategy 3: boundary equilibrium (q*=0 or q*=1)
    if np.isnan(q_star):
        u0 = U(1e-9)    # proxy for q=0 (avoid exact 0 division issues)
        u1 = U(1-1e-9)  # proxy for q=1
        if np.isnan(u0) and np.isnan(u1):
            return np.nan, np.nan
        elif np.isnan(u0):
            q_star = 1.0
        elif np.isnan(u1):
            q_star = 0.0
        else:
            q_star = 1.0 if u1 >= u0 else 0.0

    #Evaluate P2 at q*
    ax = AX(q_star)
    bx = BX(q_star)
    if ax is None or bx is None or np.isnan(ax) or np.isnan(bx):
        return np.nan, q_star

    N0 = (ax + bx*p)**2
    N1 = 2*bx*(ax + bx*p)
    N2 = bx**2

    D0 = p*(1-p)
    D1 = 1 - 2*p

    if D0 == 0:
        return np.nan, q_star

    return N2/D0 - N1*D1/D0**2 + N0*(1/D0**2 + D1**2/D0**3), q_star

#Main sweep 
results        = []   # list of dicts, one per (N, lam, R) configuration
nan_failures   = 0    # points where all three strategies failed
total_attempts = 0    # total (rho*, N, lam, R) evaluations

total = len(N_vals) * len(lam_vals) * len(R_vals)
done  = 0

print(f"Sweeping {total} parameter configurations …")
print(f"  N ∈ {{{N_vals[0]}, …, {N_vals[-1]}}},  "
      f"λ ∈ [{lam_vals[0]:.4f}, …, {lam_vals[-1]:.1f}],  "
      f"R ∈ [{R_vals[0]:.1f}, …, {R_vals[-1]:.0f}]")
print()

q_prev = 0.5   

for N_val, lam_val, R_val in itertools.product(N_vals, lam_vals, R_vals):
    p2_grid = []
    q_rho   = q_prev   

    for rho in rho_grid:
        total_attempts += 1
        val, q_sol = compute_P2(rho, N_val, lam_val, R_val, q0_guess=q_rho)
        if not np.isnan(q_sol):
            q_rho  = q_sol
            q_prev = q_sol
        if np.isnan(val):
            nan_failures += 1
        p2_grid.append(val)

    p2_arr  = np.array(p2_grid, dtype=float)
    valid   = p2_arr[~np.isnan(p2_arr)]
    n_valid = len(valid)

    if n_valid == 0:
        min_p2     = np.nan
        argmin_rho = np.nan
        any_neg    = False
    else:
        idx_min    = np.nanargmin(p2_arr)
        min_p2     = p2_arr[idx_min]
        argmin_rho = rho_grid[idx_min]
        any_neg    = bool(np.any(valid < -1e-10))

    results.append({
        'N'          : N_val,
        'lambda'     : lam_val,
        'R'          : R_val,
        'min_P2'     : min_p2,
        'rho_at_min' : argmin_rho,
        'P2_negative': any_neg,
        'n_valid'    : n_valid,
    })

    done += 1
    if done % 500 == 0 and done != total:
        print(f"  {done}/{total} done …")

print(f"  {done}/{total} done.")
print()

# Summary
neg_cases = [r for r in results if r['P2_negative']]
print("=" * 60)
print(f"Total configurations evaluated : {total}")

#  Solver reliability report 
nan_rate = 100 * nan_failures / total_attempts if total_attempts > 0 else 0.0
print(f"Truly unresolvable points      : {nan_failures} / {total_attempts}  ({nan_rate:.3f}%)")
print(f"  (boundary equilibria q*∈{{0,1}} are valid — not counted as failures)")
if nan_failures == 0:
    print("    ✓ Every point resolved: interior root or boundary equilibrium.")
elif nan_rate < 1.0:
    print(f"    ⚠ {nan_rate:.3f}% of points unresolvable — results are reliable.")
else:
    print(f"    ✗ {nan_rate:.1f}% of points unresolvable — interpret with caution.")
print()

# Sign of P2
if not neg_cases:
    print(">>> P2 is NON-NEGATIVE for ALL configurations tested (threshold: < -1e-10).")
    print("    Informational dominance FAILS everywhere in the parameter space explored.")
else:
    print(f">>> P2 IS NEGATIVE for some configurations (threshold: < -1e-10).")
    print(f"    {len(neg_cases)} / {total} configurations show informational dominance.")
    print("    Informational dominance HOLDS in those regions.")

print()

# Global minimum across the entire experiment 
all_min_p2 = [r['min_P2'] for r in results if not np.isnan(r['min_P2'])]
if all_min_p2:
    global_min_rec = min(
        (r for r in results if not np.isnan(r['min_P2'])),
        key=lambda r: r['min_P2']
    )
    global_min = global_min_rec['min_P2']
    print(f"Global minimum of P₂ over all (N, λ, R, ρ*) : {global_min:.6e}")
    print(f"  achieved at  N={global_min_rec['N']},  "
          f"λ={global_min_rec['lambda']:.5g},  "
          f"R={global_min_rec['R']:.5g},  "
          f"ρ*≈{global_min_rec['rho_at_min']:.4f}")

df = pd.DataFrame(results)
df.to_csv("P2_parameter_sweep_results.csv", index=False)
print("Results saved to P2_parameter_sweep_results.csv")
