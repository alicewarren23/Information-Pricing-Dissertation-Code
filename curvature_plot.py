
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

# Parameters
N   = 3
lam = 1
R   = 100.0

#Structural coefficients
def _inner(q, p, rho):
    """Shared inner term used in both D_Y and D_X."""
    return (
        p*(N*p**2
           + N*(-1+(1-p)**N)*(-1+q)*rho
           - p*(-1+q)*(1+(-1+(1-p)**N)*rho)
           + N*p*(-1+(-1+q)*rho))
        - p**N*(q - N*q + p*(p + (-2+N)*q + p*(-1+q)*rho))
    )

def DY_val(q, p, rho=None):
    """D_Y = inner^2 / ((N-1)^2*(1-p)^4*p^4)"""
    if rho is None:
        rho = p   # symmetric benchmark
    t = _inner(q, p, rho)
    return t**2 / ((N-1)**2*(1-p)**4*p**4)


def DX_val(q, p, rho=None):
    """D_X = inner_a * inner_b / ((N-1)^2*(1-p)^4*p^3)"""
    if rho is None:
        rho = p
    # inner_a = -inner (sign-flipped version used in Mathematica X)
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
    """A_X coefficient (rho_pub = rho; for symmetric benchmark rho = p_eq = p)."""
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
    """B_X coefficient: q multiplies the entire numerator"""
    return (
        (1 - (1-p)**N + (-1-N+(1-p)**N)*p + p**N*(-1+N+p)) * q
        / ((N-1)*(p-1)*p**2)
    )


def X_full(q, p):
    """Full X(q, rho=p, p_eq=p)"""
    rho = p
    sj  = rho * R
    AX  = AX_val(q, p)
    BX  = BX_val(q, p)
    DX  = DX_val(q, p)
    return AX + BX*rho - 2*lam*(1-rho)*rho/(sj+3)*DX


def Y_full(q, p):
    """Full Y(q, rho=p, p_eq=p)"""
    rho = p
    sj  = rho * R
    return -lam*(1-rho)*rho/(sj+3) * DY_val(q, p)


def find_q_star(p, q_guess=0.5):
    """Interior equilibrium q* at rho=p_eq=p."""
    def F(q):
        y = Y_full(q, p)
        if abs(y) < 1e-12:
            return np.inf
        return q + X_full(q, p) / (2*Y_full(q, p))
    sol = root(F, q_guess, method='hybr')
    if not sol.success:
        return np.nan
    return sol.x[0]


#Analytical P2 

def P2_val(rho_star):
    p  = rho_star
    qs = find_q_star(p)
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


def full_coeff(rho_star):
    """ Full curvature coefficient  3*P2 / (4*lambda*D_Y)"""
    p  = rho_star
    qs = find_q_star(p)
    if np.isnan(qs):
        return np.nan, np.nan

    p2  = P2_val(p)
    dy  = DY_val(qs, p)

    coeff = 3*p2 / (4*lam*dy)
    return p2, coeff


# Sweep over rho^*
rho_vals = np.linspace(0.05, 0.95, 300)

P2_vals    = []
coeff_vals = []
q_guess    = 0.5

for r in rho_vals:
    qs = find_q_star(r, q_guess)
    if np.isnan(qs):
        P2_vals.append(np.nan)
        coeff_vals.append(np.nan)
        continue
    q_guess = qs

    p2, coeff = full_coeff(r)
    P2_vals.append(p2)
    coeff_vals.append(coeff)

P2_vals    = np.array(P2_vals)
coeff_vals = np.array(coeff_vals)

#Plot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Left
axes[0].plot(rho_vals, P2_vals, color='steelblue')
axes[0].axhline(0, color='black', linestyle='--', linewidth=1)
axes[0].set_xlabel(r'True parameter $\rho^*$')
axes[0].set_ylabel(r'$P_2(\rho^*)$')
axes[0].set_title(r'Curvature coefficient $P_2(\rho^*)$')
axes[0].minorticks_on()
axes[0].grid(True, which='minor', linestyle=':', linewidth=0.8, alpha=0.7)

# Right
axes[1].plot(rho_vals, coeff_vals, color='darkorange')
axes[1].axhline(0, color='black', linestyle='--', linewidth=1)
axes[1].set_xlabel(r'True parameter $\rho^*$')
axes[1].set_ylabel(r'$\frac{3P_2}{4\lambda D_Y}$')
axes[1].set_title(r'Informational dominance coefficient $\frac{3P_2}{4\lambda D_Y}$')
axes[1].minorticks_on()
axes[1].grid(True, which='minor', linestyle=':', linewidth=0.8, alpha=0.7)

fig.suptitle(
    rf'Symmetric belief benchmark: $\hat{{\rho}}_{{pub}}=\rho^*$, $N={N}$, $\lambda={lam}$, $R={int(R)}$',
    fontsize=12
)
plt.tight_layout()
plt.show()
