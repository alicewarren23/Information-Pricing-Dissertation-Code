

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl


def simulate_markov_chain(R_max, p11, p01, seed=123):
    """Simulate a two-state Markov chain of length R_max."""
    np.random.seed(seed)
    w = np.zeros(R_max + 1)
    w[0] = 1
    for t in range(1, R_max + 1):
        if w[t - 1] == 1:
            w[t] = np.random.rand() < p11
        else:
            w[t] = np.random.rand() < p01
    return w


def compute_decomposition(w, R_max, rho_star):
    S_I = 0
    S_J = 0
    M_R = 0

    errors = []
    martingale_terms = []
    bias_terms = []

    for t in range(1, R_max + 1):
        J_t = 1 if w[t - 1] == 1 else 0
        I_t = 1 if (w[t - 1] == 1 and w[t] == 1) else 0

        S_I += I_t
        S_J += J_t

        D_t = I_t - J_t * rho_star
        M_R += D_t

        rho_hat = (S_I + 1) / (S_J + 2)
        noise = M_R / (S_J + 2)
        bias = (1 - 2 * rho_star) / (S_J + 2)
        total_error = rho_hat - rho_star

        errors.append(total_error)
        martingale_terms.append(noise)
        bias_terms.append(bias)

    return (
        np.array(errors),
        np.array(martingale_terms),
        np.array(bias_terms),
    )


def plot_decomposition(R_vals, errors, martingale_terms, bias_terms,
                       save_path="decomposition_plot.png"):
    COLOR_ERROR = '#2E86AB'
    COLOR_MART  = '#E84855'
    COLOR_BIAS  = '#F4A261'
    COLOR_ZERO  = '#6B6B6B'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)

    for ax, xlim in zip((ax1, ax2), (1000, 200)):
        ax.plot(R_vals, errors,
                color=COLOR_ERROR, linewidth=2.0, zorder=3,
                label=r'Total error $\hat{\rho}_R - \rho^*$')
        ax.plot(R_vals, martingale_terms,
                linestyle='--', color=COLOR_MART, linewidth=1.8, zorder=3,
                label=r'Martingale term $M_R / (S_J+2)$')
        ax.plot(R_vals, bias_terms,
                linestyle=':', color=COLOR_BIAS, linewidth=1.8, zorder=3,
                label=r'Bias term $(1-2\rho^*)/(S_J+2)$')
        ax.axhline(0, color=COLOR_ZERO, linewidth=0.8, linestyle='-',
                   alpha=0.5, zorder=1)

        ax.set_xlim(0, xlim)
        ax.set_xlabel(r'Sample size $R$', fontsize=12)
        ax.set_title(rf'Martingale decomposition ($R \leq {xlim}$)',
                     fontsize=13, pad=10)

        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.grid(True, which='major', linestyle='-',  linewidth=0.6, alpha=0.35)
        ax.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.45)
        ax.tick_params(axis='both', labelsize=10)
        ax.legend(frameon=True, fontsize=10)

    ax1.set_ylabel(r'Error components', fontsize=12)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def main():
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams.update({
        'axes.spines.top':   False,
        'axes.spines.right': False,
    })

    # Parameters 
    R_max    = 1000
    p11      = 0.7
    p01      = 0.3
    rho_star = p11 

    # Simulate 
    w = simulate_markov_chain(R_max, p11, p01, seed=123)

    #Decompose 
    errors, martingale_terms, bias_terms = compute_decomposition(
        w, R_max, rho_star
    )
    R_vals = np.arange(1, R_max + 1)

    #Plot
    plot_decomposition(R_vals, errors, martingale_terms, bias_terms,
                       save_path="decomposition_plot.png")


if __name__ == "__main__":
    main()
