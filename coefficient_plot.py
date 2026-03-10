import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.rcParams.update(mpl.rcParamsDefault)

#Parameters
N   = 3
lam = 1
R   = 100.0


#Structural coefficient helpers

def _inner(q, p, rho):
    return (
        p*(N*p**2
           + N*(-1+(1-p)**N)*(-1+q)*rho
           - p*(-1+q)*(1+(-1+(1-p)**N)*rho)
           + N*p*(-1+(-1+q)*rho))
        - p**N*(q - N*q + p*(p + (-2+N)*q + p*(-1+q)*rho))
    )


def DY_val(q, p, rho=None):
    if rho is None:
        rho = p
    t = _inner(q, p, rho)
    return t**2 / ((N-1)**2*(1-p)**4*p**4)


def DX_val(q, p, rho=None):
    if rho is None:
        rho = p

    inner_a = (p**N*(q - N*q + p*(p + (-2+N)*q + p*(-1+q)*rho))
               + p*(-N*p**2 - N*(-1+(1-p)**N)*(-1+q)*rho
                   + p*(-1+q)*(1+(-1+(1-p)**N)*rho)
                   + N*p*(1 + rho - q*rho)))

    inner_b = (p**N*(p - N*q - p*rho + p*q*(-1+N+rho))
               + p*(-(N-1)*p**2 + q + rho
                    + (-N - (N-1)*(1-p)**N*(-1+q) - q + N*q)*rho
                    + p*(-2 + N*(1 + rho - q*rho))))

    return inner_a * inner_b / ((N-1)**2*(p-1)**4*p**3)


def AX_val(q, p, rho=None):
    if rho is None:
        rho = p
    return (
        q + (-1 + p**(-2+N))*q
        + N*(p - p**N)*q / ((N-1)*(p-1)*p)
        + (1-q)*(-1 + N*(1-(1-p)**(N-1))*rho / ((N-1)*p))
        - ((-1+q)*(
            (N-1)*p**2
            + p**N*(-1+rho)
            - (N-1)*(-1+(1-p)**N)*rho
            - p*(-2+N+N*rho)
        )) / ((N-1)*(p-1)**2)
    )


def BX_val(q, p):
    return (
        (1 - (1-p)**N + (-1-N+(1-p)**N)*p + p**N*(-1+N+p)) * q
        / ((N-1)*(p-1)*p**2)
    )


#Full X and Y

def X_full(q, p):
    rho = p
    sj  = rho * R
    AX  = AX_val(q, p)
    BX  = BX_val(q, p)
    DX  = DX_val(q, p)

    return AX + BX*rho - 2*lam*(1-rho)*rho/(sj+3)*DX


def Y_full(q, p):
    rho = p
    sj  = rho * R
    return -lam*(1-rho)*rho/(sj+3) * DY_val(q, p)


#Equilibrium strategy

def q_star(p):

    q_guess = 0.5

    X = X_full(q_guess, p)
    Y = Y_full(q_guess, p)

    if abs(Y) < 1e-12:
        return np.nan

    q_int = -X / (2*Y)

    return np.clip(q_int, 0, 1)


# Analytical curvature coefficient P2

def P2_val(p):

    qs = q_star(p)

    if np.isnan(qs):
        return np.nan

    AX = AX_val(qs, p)
    BX = BX_val(qs, p)

    N0 = (AX + BX*p)**2
    N1 = 2*BX*(AX + BX*p)
    N2 = BX**2

    D0 = p*(1-p)
    D1 = 1 - 2*p

    return N2/D0 - N1*D1/D0**2 + N0*(1/D0**2 + D1**2/D0**3)


def full_coeff(p):

    qs = q_star(p)

    if np.isnan(qs):
        return np.nan, np.nan

    p2 = P2_val(p)
    dy = DY_val(qs, p)

    coeff = 3*p2 / (4*lam*dy)

    return p2, coeff


#Sweep and plot

if __name__ == "__main__":

    rho_vals = np.linspace(0.1, 0.9, 400)

    P2_vals = []
    coeff_vals = []

    for r in rho_vals:

        p2, coeff = full_coeff(r)

        P2_vals.append(p2)
        coeff_vals.append(coeff)

    P2_vals = np.array(P2_vals)
    coeff_vals = np.array(coeff_vals)

    print("Minimum P2 across belief grid:", np.nanmin(P2_vals))


    #Plot

    COLOR_P2    = '#2E86AB' 
    COLOR_COEFF = '#E84855' 
    COLOR_ZERO  = '#6B6B6B'   # neutral grey for reference line
    FILL_ALPHA  = 0.12

    plt.rcParams.update({
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))

    #Left Plot: P2
    ax0 = axes[0]
    ax0.plot(rho_vals, P2_vals, color=COLOR_P2, linewidth=2.2, zorder=3)
    ax0.fill_between(rho_vals, P2_vals, 0,
                     where=(P2_vals >= 0), color=COLOR_P2, alpha=FILL_ALPHA, zorder=2)
    ax0.fill_between(rho_vals, P2_vals, 0,
                     where=(P2_vals < 0),  color=COLOR_COEFF, alpha=FILL_ALPHA, zorder=2)
    ax0.axhline(0, color=COLOR_ZERO, linestyle='--', linewidth=1, zorder=1)
    ax0.set_xlabel(r'True parameter $\rho^*$', fontsize=12)
    ax0.set_ylabel(r'$P_2(\rho^*)$', fontsize=12)
    ax0.set_title(r'Curvature coefficient $P_2(\rho^*)$', fontsize=13, pad=10)
    ax0.set_xlim(0.1, 0.9)
    ax0.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax0.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax0.grid(True, which='major', linestyle='-', linewidth=0.6, alpha=0.35)
    ax0.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.45)
    ax0.tick_params(axis='both', labelsize=10)

    #Right Plot: full coefficient
    ax1 = axes[1]
    ax1.plot(rho_vals, coeff_vals, color=COLOR_COEFF, linewidth=2.2, zorder=3)
    ax1.fill_between(rho_vals, coeff_vals, 0,
                     where=(np.array(coeff_vals) >= 0), color=COLOR_COEFF, alpha=FILL_ALPHA, zorder=2)
    ax1.fill_between(rho_vals, coeff_vals, 0,
                     where=(np.array(coeff_vals) < 0),  color=COLOR_P2, alpha=FILL_ALPHA, zorder=2)
    ax1.axhline(0, color=COLOR_ZERO, linestyle='--', linewidth=1, zorder=1)
    ax1.set_xlabel(r'True parameter $\rho^*$', fontsize=12)
    ax1.set_ylabel(r'$\frac{3P_2}{4\lambda D_Y}$', fontsize=13)
    ax1.set_title(r'Payoff coefficient $\frac{3P_2}{4\lambda D_Y}$', fontsize=13, pad=10)
    ax1.set_xlim(0.1, 0.9)
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.grid(True, which='major', linestyle='-', linewidth=0.6, alpha=0.35)
    ax1.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.45)
    ax1.tick_params(axis='both', labelsize=10)

    fig.suptitle(
        rf'Symmetric belief benchmark: $\hat{{\rho}}_{{pub}}=\rho^*$, $N={N}$, $\lambda={lam}$, $R={int(R)}$',
        fontsize=13, y=1.02
    )

    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig("equilibrium_curvature.pdf", dpi=300, bbox_inches="tight")
    plt.show()
