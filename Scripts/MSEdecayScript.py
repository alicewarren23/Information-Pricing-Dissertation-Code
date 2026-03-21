

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl


def compute_mse_curves(p11, p01, R_max=100, n_points=1000):
    pi1 = p01 / (p01 + (1 - p11))
    C_var = p11 * (1 - p11) / pi1

    R = np.linspace(1, R_max, n_points)
    mse_second_order = C_var / R + 0.1 / R ** 2
    mse_asymptotic   = C_var / R

    return R, mse_second_order, mse_asymptotic


def plot_mse_decay(R, mse_second_order, mse_asymptotic,
                   save_path="mse_decay_plot.png"):
    
    fig, ax = plt.subplots(figsize=(7, 4.8))

    ax.plot(
        R, mse_second_order,
        linewidth=2.2,
        color='#2E86AB',
        label=r'$\delta_R \approx C_{\mathrm{Var}}/R + O(R^{-2})$',
        zorder=3,
    )
    ax.plot(
        R, mse_asymptotic,
        linestyle='--',
        linewidth=1.8,
        color='#C91919',
        alpha=0.85,
        label=r'$C_{\mathrm{Var}}/R$ asymptote',
        zorder=3,
    )

    ax.axhline(0, color='#6B6B6B', lw=0.8, ls='-', alpha=0.4, zorder=1)

    ax.set_xlabel(r'History length $R$', fontsize=12)
    ax.set_ylabel(r'Mean squared error $\delta_R$', fontsize=12)
    ax.set_title(r'Theoretical decay of mean squared error', fontsize=13, pad=10)

    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(True, which='major', linestyle='-',  linewidth=0.6, alpha=0.35)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.45)
    ax.tick_params(axis='both', labelsize=10)

    ax.legend(frameon=True, fontsize=10)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def main():
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams.update({
        'axes.spines.top':   False,
        'axes.spines.right': False,
    })

    p11 = 0.7
    p01 = 0.3

    R, mse_second_order, mse_asymptotic = compute_mse_curves(
        p11=p11, p01=p01, R_max=100, n_points=1000
    )

    plot_mse_decay(R, mse_second_order, mse_asymptotic,
                   save_path="mse_decay_plot.png")


if __name__ == "__main__":
    main()
