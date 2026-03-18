
import numpy as np
import matplotlib.pyplot as plt

# SIMULATION FUNCTION

def simulate_markov_chain(R, p11, p01, seed=None):
  
    if seed is not None:
        np.random.seed(seed)

    w = np.zeros(R + 1)
    w[0] = 1  

    for t in range(1, R + 1):
        if w[t - 1] == 1:
            w[t] = np.random.rand() < p11
        else:
            w[t] = np.random.rand() < p01

    return w

# MARTINGALE DECOMPOSITION

def martingale_decomposition(w, rho_star):

    S_I = 0  
    S_J = 0  
    M_R = 0

    errors = []
    martingale_terms = []
    bias_terms = []

    for t in range(1, len(w)):

        J_t = 1 if w[t - 1] == 1 else 0
        I_t = 1 if (w[t - 1] == 1 and w[t] == 1) else 0

        S_I += I_t
        S_J += J_t

        # martingale increment
        D_t = I_t - J_t * rho_star
        M_R += D_t

        # estimator
        rho_hat = (S_I + 1) / (S_J + 2)

        # decomposition
        noise = M_R / (S_J + 2)
        bias = (1 - 2 * rho_star) / (S_J + 2)
        total_error = rho_hat - rho_star

        errors.append(total_error)
        martingale_terms.append(noise)
        bias_terms.append(bias)

    return errors, martingale_terms, bias_terms

# PLOTTING

def plot_results(R_vals, errors, martingale_terms, bias_terms):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, xlim in zip((ax1, ax2), (len(R_vals), 200)):
        ax.plot(R_vals, errors, label='Total Error')
        ax.plot(R_vals, martingale_terms, linestyle='--', label='Martingale Term')
        ax.plot(R_vals, bias_terms, linestyle=':', label='Bias Term')

        ax.axhline(0, linewidth=0.8)
        ax.set_xlim(0, xlim)
        ax.set_xlabel('Sample Size R')
        ax.set_ylabel('Error')
        ax.set_title(f'Martingale Decomposition (R ≤ {xlim})')
        ax.grid()
        ax.legend()

    plt.tight_layout()
    plt.show()

# MAIN SCRIPT

def main():
    # PARAMETERS
    R_max = 1000
    p11 = 0.7
    p01 = 0.3
    rho_star = p11

    # SIMULATION
    w = simulate_markov_chain(R_max, p11, p01, seed=123)

    # DECOMPOSITION
    errors, martingale_terms, bias_terms = martingale_decomposition(w, rho_star)

    # PLOT
    R_vals = np.arange(1, R_max + 1)
    plot_results(R_vals, errors, martingale_terms, bias_terms)


if __name__ == "__main__":
    main()
