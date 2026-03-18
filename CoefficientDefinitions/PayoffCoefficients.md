## Payoff Coefficients 

The document contains the payoff coefficients used player $i$'s in the volatility-neutral payoff (equation 4.8) when facing player $j$ and N-1 equally informed players. 

$$
\Theta(\hat\rho_i, \vec{p}_{-i}) = q_j A_N(\hat\rho_i, p_{eq}) + B_N(\hat\rho_i, p_{eq})
$$


$$
\Omega(\hat\rho_i, \vec{p}_{-i}) = C_{N}(\hat\rho_i, p_{eq}) + D_N(\hat\rho_i, p_{eq})
$$


where 

$$
A_N(\hat\rho_i, p)  = \frac{\hat\rho_i(p^N(N(p-1) - 2p + 1) + (1-p)^N + p((N-2)(1-p)^N + 2) -1)}{(N-1)(p-1)^2p^2} + \frac{p^N(N(-p) + N + 2p -1) -p^2}{(N-1)(p-1)^2p^2}
$$

$$
B_N(\hat\rho_i, p) = \frac{\hat\rho_i(p(p^N + (1-p)^N -1) - N((1-p)^N + p - 1))}{(N-1)(p-1)^2p} + \frac{-p^N + N(p-1)+1}{(N-1)(p-1)^2}
$$

$$
C_N(\hat\rho_i, p) = \frac{\hat\rho_i((N(-p) + N + p) - N((1-p)^N + p-1))}{(N-1)(p-1)^2p} + \frac{(N(p-1) -p)p^N + p}{(N-1)(p-1)^2}
$$

$$
D_N(\hat\rho_i, p) = \frac{p^N - (N-1)p^2 + (N-2)p}{(N-1)(p-1)^2} + \frac{\hat\rho_i(p^N - ((N-1)(1-p)^N) - Np + N -1)}{(N-1)(p-1)^2}
$$
