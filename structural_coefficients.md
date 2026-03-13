## Structural Coefficient Notebook 

This documents contains the coefficients that appears in the volatility-averse payoff for a player $i$ when facing a player $j$ and $N-2$ equally informed players.(equation 4.5). 

The following expressions the variables are given, 

- $N$ = the number of players in the game 
- $\lambda$ = the voltility aversion coefficient
- $S_{J,i}$ = the number of transitions in player $i$s string from state 1 to any other state (0 or 1).
- $\rho_{pub}$ = the estimator of the true parameter held by the N-2 equally informed players
- $p_{eq}$ = the strategy played by the N-2 equally informed players
- $q_j$ = the strategy played my player $j$
--------



$$
\boxed{Y^{(i)} =\frac{\lambda(\hat\rho_i - 1)\hat\rho_i}{S_{J,i} + 3}D_Y}
$$

where 

$$
D_Y=  \frac{1}{(N-1)^2(1-p_{eq})^4 p_{eq}^4} (p_{eq}(Np_{eq}^2 + N(-1+(1-p_{eq})^ N)(-1+q_j)\hat\rho_{pub}-p_{eq}(-1+q_i)(1+(-1+(1-p_{eq})^N\hat\rho_{pub})+ Np_{eq}(-1+(-1+q_j)\hat\rho_{pub})) - p_{eq}^N(q_j - Nq_j+p_{eq}(p_{eq} +(-2+N)q_j + p_{eq}(-1+q_j)\hat\rho_{pub})))^2
$$

$$
\boxed{X^{(i)} = A_X + B_X \hat\rho_i -  \frac{2\lambda(1 - \hat\rho_i)\hat\rho_i}{S_{J,i} + 3} D_X}
$$

where 

$$
A_X  = q_j + (-1+p_{eq}^{-2+N})q_j + \frac{N(p_{eq}-p_{eq}^N)q_j}{(-1+N)(-1+p_{eq})p_{eq}}  + (1-q_j)(-1 + \frac{n(1-(1-p_{eq})^{N-1})\hat\rho_{pub}}{(-1+N)p_{eq}})- \frac{(-1+q_j)((-1+N)p_{eq}^2 + p_{eq}^{N}(-1+\hat\rho_{pub})- (-1+N)(-1 + (1-p_{eq})^N)\hat\rho_{pub}-p_{eq}(-2+N+N\hat\rho_{pub}))}{(-1+N)(-1+p_{eq})^2}
$$

$$
B_X = \frac{(1-(1-p_{eq})^N+(-1-N + (1-p_{eq})^N)p_{eq}+p_{eq}^N(-1+
N+p_{eq}))q_j}{(-1+N)(-1+p_{eq})p_{eq}^2}
$$

and 

$$
D_X = \frac{1}{(-1+N)^2(-1+p_{eq})^4p_{eq}^3}  (p_{eq}^N(q_j - Nq_j + p_{eq}(p_{eq}+(-2+N)q_j + p_{eq}(-1+q_j)\hat\rho_{pub})) + p_{eq}(-Np_{eq}^2 - N(-1 + (1-p_{eq})^N)(-1 + q_j)\hat\rho_{pub} + p_{eq}(-1+q_j)(1+(-1+(1-p_{eq})^N)\hat\rho_{pub})+Np_{eq}(1+\hat\rho_{pub} - q_j\hat\rho_{pub})))(p_{eq}(p_{eq}- Nq_j - p_{eq}\hat\rho_{pub}+ p_{eq}q_j(-1+N + \hat\rho_{pub})) + p_{eq}(-((-1+N)p_{eq}^2) + q_j + \hat\rho_{pub}+ (-N-(-1+N))(1-p_{eq})^N (-1+q_j) - q_j + Nq_j)\hat\rho_{pub} + p_{eq}(-2+N(1+ \hat\rho_{pub} - q_j \hat\rho_{pub}))))
$$

$$
\boxed{R^{(i)} = A_R + B_R \hat\rho_i - \frac{\lambda(\hat\rho_i - 1)\hat\rho_i}{S_{J,i} + 3}D_R}
$$

where

$$
A_R = (-1 + \frac{N(-p_{eq}+p_{eq}^N)}{(-1+N)(-1+ p_{eq})p_{eq}})q_j + \frac{(-1+q_j)((-1+N)p_{eq}^2 + p_{eq}^N(-1+ \hat\rho_{pub})-(-1+N)(-1+(1-p_{eq})^N)\hat\rho_{pub}(-2+N+N\hat\rho_{pub}))}{(-1+N)(-1+p_{eq})^2}
$$

$$
B_R = \frac{N(p_{eq} - p_{eq}^N)q_j}{(-1+N)(-1+p_{eq})p_{eq}}
$$

and 

$$
D_R = \frac{(-1+q_j)((-1+N)p_{eq}^2+ p_{eq}^N (-1+ \hat\rho_{pub})- (-1+N)(-1 + (1- p_{eq})^N)\hat\rho_{pub}- p_{eq}(-2+ N + N\hat\rho_{pub}))}{(-1+N)(-1+p_{eq})^2})-1 + \frac{n(-p_{eq}+p_{eq}^N)}{(-1+N)(-1+p_{eq})p_{eq}}q_j 
$$
