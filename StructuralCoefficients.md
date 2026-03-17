## Structural Coefficient Notebook

This document contains the structural coefficients that appear in the volatility-averse payoff for player $i$ when facing one opponent $j$
and $N-2$ equally informed players (Equation 4.5 of the dissertation). These coefficients are derived symbolically from the payoff representation
in Data as Commodity (Appendix B) and expressed in closed form for the
asymmetric belief configuration used in the dissertation.

The following expressions use the variables:

- $N$ : total number of players in the game
- $\lambda$ : volatility-aversion coefficient
- $S_{J,i}$ : number of transitions from state 1 in player $i$'s observed string
- $\rho_{pub}$ : estimator of the true parameter held by the $N-2$ equally informed players
- $p_{eq}$ : strategy played by the $N-2$ equally informed players
- $q_j$ : strategy played by player $j$
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
\boxed{X^{(i)} = A_R + B_R \hat\rho_i - \frac{\lambda(\hat\rho_i - 1)\hat\rho_i}{S_{J,i} + 3}D_R}
$$

where

$$
A_K = (-1 + \frac{N(-p_{eq}+p_{eq}^N)}{(-1+N)(-1+ p_{eq})p_{eq}})q_j + \frac{(-1+q_j)((-1+N)p_{eq}^2 + p_{eq}^N(-1+ \hat\rho_{pub})-(-1+N)(-1+(1-p_{eq})^N)\hat\rho_{pub}(-2+N+N\hat\rho_{pub}))}{(-1+N)(-1+p_{eq})^2}
$$

$$
B_K = \frac{N(p_{eq} - p_{eq}^N)q_j}{(-1+N)(-1+p_{eq})p_{eq}}
$$

and 

$$
D_K = \frac{(-1+q_j)((-1+N)p_{eq}^2+ p_{eq}^N (-1+ \hat\rho_{pub})- (-1+N)(-1 + (1- p_{eq})^N)\hat\rho_{pub}- p_{eq}(-2+ N + N\hat\rho_{pub}))}{(-1+N)(-1+p_{eq})^2})-1 + \frac{n(-p_{eq}+p_{eq}^N)}{(-1+N)(-1+p_{eq})p_{eq}}q_j 
$$
