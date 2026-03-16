import numpy as np
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-whitegrid")

p11 = 0.7
p01 = 0.3

pi1 = p01 / (p01 + (1 - p11))
Cvar = p11 * (1 - p11) / pi1


R = np.linspace(1, 100, 1000)

mse_second_order = Cvar / R + 0.1 / R**2
mse_asymptotic = Cvar / R


plt.figure(figsize=(10,6))

plt.plot(R, mse_second_order,
         linewidth=2.5,
         color="#2E86AB",
         label=r"$\delta_R \approx C_{\mathrm{Var}}/R$")

plt.plot(R, mse_asymptotic,
         linestyle="--",
         linewidth=2,
         color="black",
         alpha=0.8,
         label=r"$1/R$ decay rate")

plt.xlabel(r"History length $R$", fontsize=13)
plt.ylabel(r"Mean squared error $\delta_R$", fontsize=13)
plt.title("Theoretical decay of mean squared error", fontsize=14)
plt.grid(True, linestyle="--", linewidth=0.6, alpha=0.7)
plt.legend(frameon=True)
plt.tight_layout()

plt.savefig("mse_decay_plot.png", dpi=300)
plt.show()
