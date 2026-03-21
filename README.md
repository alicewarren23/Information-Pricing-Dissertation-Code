# Numerical Illustrations for Dissertation

This repository contains the Python code used to generate the numerical figures appearing in my MSci Mathematics dissertation.

Four figures are reproduced:

1. Price-of-information bounds across equilibrium regimes analysed in Chapter 2.
2. Martingale decomposition of the Bayesian estimation error analysed in Chapter 3.
3. Mean Squared Error (MSE) decay of the Bayesian estimator derived in Chapter 3.
4. Equilibrium curvature of payoff with respect to belief analysed in Chapter 4.


The scripts implement the theoretical expressions derived in the dissertation and generate the figures used in the numerical illustrations.

---

## Repository Contents

### Price Bounds Script

```bash
PriceBounds.py
```

This script computes and plots the state-dependent price-of-information bounds $\Psi_{\min}, \Psi_{\max}$ derived in Chapter 2.

The script:

* simulates a Markov signal process,
* computes state-dependent Bayesian beliefs,
* evaluates equilibrium strategies, and
* computes the resulting price bounds across the belief space.

It classifies outcomes into equilibrium regimes (cooperative, competitive, symbiotic, exploitative, and no-deal) and highlights the regions in which trade is feasible.

Running the script produces:

```bash
PriceBounds.png
```

which corresponds to Figure 2.1.


### Martingale Decomposition Script

```bash
MartingaleDecomposition.py
```

This script simulates the Markov chain and decomposes the Bayesian estimation error into its martingale (stochastic) and bias (deterministic) components.

For a realised path, $w_0, \ldots, w_R$  with true transition parameter $\rho^*$, the decomposition is given by

$$\hat{\rho}_R - \rho^* = \frac{M_R}{S_J + 2} + \frac{1 - 2\rho^*}{S_J + 2}$$

where

$( M_R = \sum_{t=1}^{R}(I_t - J_t \rho^*) ), \quad  S_I = n_{11}$, \quad S_J = n_{11} + n_{10}$

Running the script produces:

```bash
MartingaleDecomposition.png
```

which corresponds to Figure 3.1 in the dissertation.

---
### MSE Decay Script

```bash
MSEdecayScript.py
```

This script generates the figure illustrating the decay of the mean squared error

$$
\delta_R = \mathbb{E}\big[(\hat{\rho}_R - \rho^*)^2\big]
$$

as a function of the sample size R.

The plot illustrates the theoretical result derived in Chapter 3:

$$
\delta_R \approx \frac{C_{\mathrm{Var}}}{R},
$$

showing that estimation error declines at the parametric rate ( 1/R ).

The implementation uses example transition probabilities for the Markov chain:

```bash
p_11 = 0.7
p_01 = 0.3
```

These determine the stationary distribution

$$
\pi_1 = \frac{p_{01}}{p_{01} + (1 - p_{11})}
$$

and the variance constant

$$
C_{\mathrm{Var}} = \frac{p_{11}(1 - p_{11})}{\pi_1}.
$$

Running the script produces:

```bash
MSEdecayFig.png
```

which corrsponds to Figure 3.2.
---

### Equilibrium Curvature Script

```bash
CoefficientPlotScript.py
```

This script evaluates the equilibrium curvature coefficient

$$
P_2(\rho^*)
$$

and the associated payoff coefficient

$$
\frac{3 P_2}{4 \lambda D_Y}.
$$

These quantities determine whether informational precision generates a payoff advantage in equilibrium.

The coefficients are evaluated across the belief space

$$
\rho^* \in [0.1, 0.9]
$$

under the symmetric benchmark where public belief coincides with the true parameter.

Baseline parameters:

```bash
N = 3
lambda = 1
R = 100
```

Running this script produces:

```bash
EquilibriumCurvatureFig.png
```

which corresponds tp Figure 4.1.

---

### Structural Coefficients

The repository also includes documentation of the structural coefficients appearing in the payoff functions of the model. These correspond to:

- the volatility-neutral payoff of player $i$ (Equation 4.8), and
- the volatility-averse payoff of player $i$ (Equation 4.17).

These are provided as Markdown documentation rather than executable scripts.

---

## Requirements

```bash
numpy
matplotlib
scipy
```

Install using:

```bash
pip install numpy matplotlib scipy
```
---

## Running the Scripts

```bash
python PriceBounds.py
python MartingaleDecomposition.py
python MSEdecayScript.py
python CoefficientPlotScript.py
```
---

## Relation to Dissertation

### Chapter 2: Markovian Learning Environment

* State-dependent price-of-information bounds
* Feasibility of information trade across equilibrium regimes

### Chapter 3: Statistical Properties of Bayesian Learning

* $1/R$ decay of mean squared error
* Martingale decomposition of estimation error

### Chapter 4: Strategic Curvature and Informational Dominance

* Curvature of equilibrium payoff
* Numerical evaluation of the curvature coefficient

---

## Author

Alice Warren
MSci Mathematics Dissertation
2026
