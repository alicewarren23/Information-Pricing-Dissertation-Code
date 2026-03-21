

import numpy as np
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from itertools import groupby
from matplotlib.lines import Line2D
from scipy.optimize import brentq


#Regimes

REGIME_LABELS = {1: '1', 2: '2', 4: '4', 5: '5', 6: '6'}

REGIME_FULLNAMES = {
    1: '(1) Cooperative',
    2: '(2) Competitive',
    4: '(4) Symbiotic',
    5: '(5) Exploitative',
    6: '(6) No-deal',
}

FEASIBLE_REGIMES = {1, 4}


#Belief estimation

def compute_rho_hat(sequence, w_R):
    n_11 = n_10 = n_01 = n_00 = 0
    for prev, curr in zip(sequence[:-1], sequence[1:]):
        if   prev == 1 and curr == 1: n_11 += 1
        elif prev == 1 and curr == 0: n_10 += 1
        elif prev == 0 and curr == 1: n_01 += 1
        elif prev == 0 and curr == 0: n_00 += 1

    if w_R == 1:
        return (n_11 + 1) / (n_11 + n_10 + 2)
    elif w_R == 0:
        return (n_01 + 1) / (n_01 + n_00 + 2)
    else:
        raise ValueError("w_R must be 0 or 1")


#Equilibrium strategies 

def compute_p_eq(rho, N):
    """Compute the symmetric pre-transaction equilibrium mixing probability."""
    if rho <= 1.0 / N:
        return 0.0
    elif rho >= 1.0 - 1.0 / N:
        return 1.0
    else:
        def f(p):
            if p <= 0 or p >= 1:
                return float('inf')
            num   = p - p ** N
            denom = 1 - p ** N - (1 - p) ** N
            if abs(denom) < 1e-14:
                return float('inf')
            return num / denom - rho
        try:
            return brentq(f, 1e-10, 1 - 1e-10)
        except Exception:
            return 0.5


def compute_q1_pre(rho1, rho, N):
    """Compute the seller's pre-transaction best-response."""
    if rho1 < 1.0 / N - 1e-12:
        return 0.0
    elif rho1 > 1.0 - 1.0 / N + 1e-12:
        return 1.0
    else:
        if   rho1 > rho + 1e-12: return 1.0
        elif rho1 < rho - 1e-12: return 0.0
        else:                     return 1.0


# Payoff coefficients

def compute_A_N(rho, p, N):
    """Coefficient A_N(rho, p) in the quadratic payoff expansion."""
    if abs(p) < 1e-14:
        return -(N * (N - 3) / (2 * (N - 1))) * rho - 1.0 / (N - 1)
    elif abs(1 - p) < 1e-14:
        return ((N - 3) * N / (2 * (N - 1))) * rho - N / 2.0 + 1
    else:
        p_N   = p ** N
        omp_N = (1 - p) ** N
        numerator = (
            rho * (p_N * (N * (p - 1) - 2 * p + 1) + omp_N +
                   p * ((N - 2) * omp_N + 2) - 1) +
            p_N * (N * (-p) + N + 2 * p - 1) - p ** 2
        )
        denominator = (N - 1) * (p - 1) ** 2 * p ** 2
        return numerator / denominator


def compute_B_N(rho, p, N):
    """Coefficient B_N(rho, p) in the quadratic payoff expansion."""
    if abs(p) < 1e-14:
        return N * rho - 1
    elif abs(1 - p) < 1e-14:
        return (N * (N + 1) / (2 * (N - 1))) * rho - N / 2.0
    else:
        p_N   = p ** N
        omp_N = (1 - p) ** N
        term1 = (rho * (p * (p_N + omp_N - 1) - N * (omp_N + p - 1))
                 / ((N - 1) * (p - 1) ** 2 * p))
        term2 = (-p_N + N * (p - 1) + 1) / ((N - 1) * (p - 1) ** 2)
        return term1 + term2


def compute_C_N(rho, p, N):
    """Coefficient C_N(rho, p) in the quadratic payoff expansion."""
    if abs(p) < 1e-14 or abs(1 - p) < 1e-14:
        return 0.0
    p_N   = p ** N
    omp_N = (1 - p) ** N
    term1 = (rho * ((N * (-p) + N + p) * p_N + p * (-((N - 1) * omp_N) - 1))
             / ((N - 1) * (p - 1) ** 2 * p))
    term2 = ((N * (p - 1) - p) * p_N + p) / ((N - 1) * (p - 1) ** 2 * p)
    return term1 + term2


def compute_D_N(rho, p, N):
    """Coefficient D_N(rho, p) in the quadratic payoff expansion."""
    if abs(p) < 1e-14 or abs(1 - p) < 1e-14:
        return 0.0
    p_N   = p ** N
    omp_N = (1 - p) ** N
    term1 = (p_N - (N - 1) * p ** 2 + (N - 2) * p) / ((N - 1) * (p - 1) ** 2)
    term2 = (rho * (p_N - (N - 1) * omp_N - N * p + N - 1)
             / ((N - 1) * (p - 1) ** 2))
    return term1 - term2


#Post-transaction strategies and payoffs

def compute_q_eq(rho2, p_eq, N):
    """Compute the symmetric post-transaction equilibrium strategy."""
    A_N = compute_A_N(rho2, p_eq, N)
    B_N = compute_B_N(rho2, p_eq, N)
    if abs(A_N) < 1e-14:
        return 0.5
    return max(0.0, min(1.0, -B_N / A_N))


def compute_q1_post(rho1, rho2, p_eq, N):
    """Compute the seller's post-transaction best-response."""
    q_eq   = compute_q_eq(rho2, p_eq, N)
    A_N    = compute_A_N(rho1, p_eq, N)
    B_N    = compute_B_N(rho1, p_eq, N)
    alpha1 = A_N * q_eq + B_N
    return 1.0 if alpha1 > 1e-12 else 0.0


def compute_U1_post(q1, q2, p_eq, rho1, N):
    """Seller payoff after the transaction."""
    A_N = compute_A_N(rho1, p_eq, N)
    B_N = compute_B_N(rho1, p_eq, N)
    C_N = compute_C_N(rho1, p_eq, N)
    D_N = compute_D_N(rho1, p_eq, N)
    return q1 * (A_N * q2 + B_N) + C_N * q2 + D_N


def compute_U_i_general(q_i, q_j, p_eq, rho_i, N):
    """General payoff for player i given own and opponent's strategies."""
    A_N = compute_A_N(rho_i, p_eq, N)
    B_N = compute_B_N(rho_i, p_eq, N)
    C_N = compute_C_N(rho_i, p_eq, N)
    D_N = compute_D_N(rho_i, p_eq, N)
    return q_i * (A_N * q_j + B_N) + C_N * q_j + D_N


#Price-of-information bounds

def compute_slice(rho1_grid, rho2_fixed, rho_public, N):
    psi_min_vals, psi_max_vals = [], []
    p_eq = compute_p_eq(rho_public, N)

    for rho1 in rho1_grid:
        q1_pre  = compute_q1_pre(rho1, rho_public, N)
        q_eq    = compute_q_eq(rho2_fixed, p_eq, N)
        q1_post = compute_q1_post(rho1, rho2_fixed, p_eq, N)

        U1_pre  = compute_U_i_general(q1_pre, p_eq,    p_eq, rho1, N)
        U1_post = compute_U1_post(q1_post, q_eq, p_eq, rho1, N)

        U2_pre  = compute_U_i_general(p_eq,    q1_pre,  p_eq, rho1, N)
        U2_post = compute_U_i_general(q_eq,    q1_post, p_eq, rho1, N)

        psi_min_vals.append(U1_pre  - U1_post)
        psi_max_vals.append(U2_post - U2_pre)

    return np.array(psi_min_vals), np.array(psi_max_vals)


# Regime classification and plotting

def classify(psi_min, psi_max):
    """Classify a (Psi_min, Psi_max) pair into one of the six regimes."""
    if   psi_max > 0 and psi_min > 0 and psi_max > psi_min: return 1
    elif psi_max > 0 and psi_min > 0 and psi_min > psi_max: return 2
    elif psi_min < 0 < psi_max:                              return 4
    elif psi_max < 0 and psi_min < 0 and psi_max > psi_min: return 5
    elif psi_max < 0 and psi_min < 0 and psi_max < psi_min: return 6
    else:                                                    return 0


def plot_state_slice(ax, rho1_grid, psi_min_vals, psi_max_vals,
                     rho2_fixed, state_label, N):
    res     = len(rho1_grid)
    regions = np.array([classify(psi_min_vals[k], psi_max_vals[k])
                        for k in range(res)])

    prev_reg  = regions[0]
    seg_start = rho1_grid[0]
    for k in range(1, res):
        if regions[k] != prev_reg or k == res - 1:
            end  = rho1_grid[k] if k < res - 1 else rho1_grid[-1]
            feas = prev_reg in FEASIBLE_REGIMES
            ax.axvspan(seg_start, end,
                       color='#ffe600' if feas else '#aaaaaa',
                       alpha=0.25, linewidth=0)
            seg_start = rho1_grid[k]
            prev_reg  = regions[k]

    ax.plot(rho1_grid, psi_min_vals, color='#2E86AB', lw=2.2, zorder=3,
            label=r'$\Psi_{\min}$')
    ax.plot(rho1_grid, psi_max_vals, color='#E84855', lw=2.2, zorder=3,
            label=r'$\Psi_{\max}$')

    ax.axhline(0, color='#6B6B6B', lw=0.9, ls='--', alpha=0.6, zorder=2)

    y_lo = min(psi_min_vals.min(), psi_max_vals.min())
    y_hi = max(psi_min_vals.max(), psi_max_vals.max())
    pad  = (y_hi - y_lo) * 0.08
    ax.set_ylim(y_lo - pad, y_hi + pad * 3)
    y_lo_ax, y_hi_ax = ax.get_ylim()

    placed = set()
    for reg_id, grp in groupby(enumerate(regions), key=lambda x: x[1]):
        idxs = [i for i, _ in grp]
        if reg_id not in REGIME_LABELS:
            continue
        if idxs[0] > 0:
            ax.axvline(rho1_grid[idxs[0]], color='#555555', lw=1.0,
                       ls='--', alpha=0.7, zorder=4)
        if reg_id not in placed:
            mid_x = rho1_grid[idxs[len(idxs) // 2]]
            ax.text(mid_x, y_hi_ax - (y_hi_ax - y_lo_ax) * 0.05,
                    REGIME_LABELS[reg_id],
                    ha='center', va='top', fontsize=10, fontweight='bold',
                    color='#333333', zorder=5)
            placed.add(reg_id)

    for thr in [1.0 / N, 1.0 - 1.0 / N]:
        ax.axvline(thr, color='#6B6B6B', lw=0.9, ls=':', alpha=0.8, zorder=2)

    ax.set_xlim(0, 1)
    ax.set_xlabel(r'$\hat{\rho}_{\mathrm{seller}}$', fontsize=12)
    ax.set_title(
        rf'State $w = {state_label}$: fixed '
        rf'$\hat{{\rho}}_{{\mathrm{{buyer}}}} = {rho2_fixed:.3g}$',
        fontsize=13, pad=10,
    )

    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(True, which='major', linestyle='-',  linewidth=0.6, alpha=0.35)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.45)
    ax.tick_params(axis='both', labelsize=10)

    return regions


#Markov sequence simulation and belief estimation ──────────────

def simulate_and_estimate(T=500, p11_true=0.75, p10_true=0.30, seed=42):
   
    np.random.seed(seed)
    seq2    = np.zeros(T, dtype=int)
    seq2[0] = np.random.randint(0, 2)
    for t in range(1, T):
        if seq2[t - 1] == 1:
            seq2[t] = int(np.random.rand() < p11_true)
        else:
            seq2[t] = int(np.random.rand() < p10_true)

    rho2_w1 = compute_rho_hat(seq2, w_R=1)
    rho2_w0 = compute_rho_hat(seq2, w_R=0)

    return rho2_w1, rho2_w0, seq2


# Main

def main():
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams.update({
        'axes.spines.top':   False,
        'axes.spines.right': False,
    })

    N          = 3
    rho_public = 0.5
    T          = 500
    seed       = 42
    n_grid     = 500

    #Simulate buyer's history and form posteriors
    rho2_w1, rho2_w0, _ = simulate_and_estimate(
        T=T, p11_true=0.75, p10_true=0.30, seed=seed
    )

    rho1_grid = np.linspace(0.0, 1.0, n_grid)

    #Compute bounds for each state
    psi_min_1, psi_max_1 = compute_slice(rho1_grid, rho2_w1, rho_public, N)
    psi_min_0, psi_max_0 = compute_slice(rho1_grid, rho2_w0, rho_public, N)

    #Plot 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2), sharey=False)

    regs1 = plot_state_slice(ax1, rho1_grid, psi_min_1, psi_max_1,
                             rho2_w1, state_label=1, N=N)
    regs2 = plot_state_slice(ax2, rho1_grid, psi_min_0, psi_max_0,
                             rho2_w0, state_label=0, N=N)

    ax1.set_ylabel(r'$\Psi$', fontsize=12)

    all_present = sorted((set(regs1) | set(regs2)) & set(REGIME_FULLNAMES))

    leg_handles = [
        Line2D([0], [0], color='#2E86AB', lw=2.2, label=r'$\Psi_{\min}$'),
        Line2D([0], [0], color='#E84855', lw=2.2, label=r'$\Psi_{\max}$'),
        mpatches.Patch(color='#ffe600', alpha=0.4,
                       label='Feasible (trade possible)'),
        mpatches.Patch(color='#aaaaaa', alpha=0.4, label='Infeasible'),
    ]
    for reg_id in all_present:
        leg_handles.append(
            mpatches.Patch(color='none', label=REGIME_FULLNAMES[reg_id],
                           linewidth=0)
        )

    fig.legend(handles=leg_handles[:4], loc='lower center', ncol=4,
               fontsize=10, frameon=True, framealpha=0.95,
               bbox_to_anchor=(0.5, -0.04))

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.17)

    plt.savefig("markov_price_bounds.png", dpi=300, bbox_inches="tight")
    print("Figure saved to: markov_price_bounds.png")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
