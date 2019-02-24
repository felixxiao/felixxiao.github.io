---
layout: post
title: "Variance Estimation in a Simple Kalman Filter"
description: ""
category: 
tags: [time series]
---

$$
  \newcommand{\R}{\mathbb{R}}
  \newcommand{\Prb}{\mathbb{P}}
  \newcommand{\Expect}{\mathrm{E}}
  \newcommand{\Var}{\mathrm{Var}}
  \newcommand{\Cov}{\mathrm{Cov}}
  \newcommand{\indep}{\rotatebox[origin=c]{90}{$\models$}}
  \newcommand{\argmin}[1]{\underset{#1}{\mathrm{arg min}}\;}
  \newcommand{\argmax}[1]{\underset{#1}{\mathrm{arg max}}\;}
  \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
  \newcommand{\diff}{\mathrm{d}}
  \DeclareMathOperator{\Tr}{Tr}
  \let\emptyset\varnothing
$$

Consider the following state space model:

$$
S_{t+1} = S_t + \epsilon_t \\
V_t = S_t + \eta_t \\
\epsilon_t \sim \mathcal{N}(0, \sigma_S^2) \\
\eta_t \sim \mathcal{N}(0, \sigma_V^2)
$$

where $$\{\epsilon_t, \eta_t\}$$ are all iid. Given we only observe $$\{V_t\}$$ from $t=1$ to $n$, how can we infer both variances?

If the $\eta_t$ noise terms were zero (we can observe the latent states), this would be easily solved by taking the lag-1 difference of $V_t$, squaring each difference and computing the sample mean to estimate $\sigma_S^2$. However with the addition of the $\eta_t$ noise, both $\sigma_S$ and $\sigma_V$ are entangled in the resulting expression:

$$
\begin{align}
    \Expect (V_{t+1} - V_t)^2
    &= \Expect (S_{t+1} - S_t + \eta_{t+1} - \eta_t)^2 \\
    &= \Expect (\epsilon_t + \eta_{t+1} - \eta_t)^2 \\
    &= \sigma_S^2 + 2 \sigma_V^2
\end{align}
$$

Now we have one equation and two unknowns. The trick is to take the lag-2 difference of the series to derive a second equation:

$$
\begin{align}
    \Expect (V_{t+2} - V_t)^2
    &= \Expect (S_{t+2} - S_t + \eta_{t+2} - \eta_t)^2 \\
    &= \Expect (S_{t+2} - S_{t+1} + S_{t+1} - S_t + \eta_{t+2} - \eta_t)^2 \\
    &= \Expect (\epsilon_{t+1} + \epsilon_t + \eta_{t+1} - \eta_t)^2 \\
    &= 2 \sigma_S^2 + 2 \sigma_V^2
\end{align}
$$

Define $Y_i$ to be the sample mean of the squared lag-$i$ difference:

$$ Y_i = \frac{1}{n-i} \sum_{t=1}^{n-i} (V_{t+i} - V_t)^2 $$

Our estimator is therefore the solution to the linear system:

$$
\begin{bmatrix}
    1 & 2 \\
    2 & 2
\end{bmatrix}
\begin{bmatrix}
    \hat{\sigma}_S^2 \\
    \hat{\sigma}_V^2
\end{bmatrix}
=
\begin{bmatrix}
    Y_1 \\
    Y_2
\end{bmatrix}
$$

It is unbiased since $\Expect[Y_i] = \Expect[(V_{t+i} - V_t)^2]$.

## The OLS estimator with $k$ lagged differences

It is possible to obtain more equations by taking larger lagged differences of $V_t$. Since this leads to an overdetermined linear system, we use the ordinary least squares estimator to decide the sigmas. The estimator has the traditional analytic form:

$$
\begin{bmatrix}
    \hat{\sigma}_S^2 \\
    \hat{\sigma}_V^2
\end{bmatrix}
=
(X^T X)^{-1} X^T
\begin{bmatrix}
    Y_1 \\
    \vdots \\
    Y_k
\end{bmatrix}
$$

where

$$
X = \begin{bmatrix}
    1 & 2 \\
    \vdots & \vdots \\
    k & 2
\end{bmatrix}
$$

Like the $k=2$ special case, the general OLS estimator is also unbiased for the same reason. The more interesting question is what is the optimal $k$ to use given $n$ and the population parameters, where optimal is defined here to be the estimator that has the minimum variance.

### Variance of the $k$ OLS estimator

The variance is hard to compute for this estimator due to the fact that the $Y_i$ are correlated in a rather complicated way. This stands in contrast with normal OLS regression, where the $Y_i$ are all assumed to be uncorrelated, and the variance of the coefficients estimator would monotonically decrease as the number of rows increases. Here, the optimal $k$ can be thought of as the equilibrium point between two conflicting forces:

* If correlation between successive $Y_i$ is sufficiently low, then having more lags decreases the estimator variance by adding new information.
* The variance of $Y_i$ increases monotically with $i$ because larger lags have more noise terms. This makes the bigger $i$ rows far less useful.[^1]

Here's the analytical form of the $k$ lags OLS estimator covariance matrix:

$$
\Var(\hat{\sigma}^2) = (X^T X)^{-1} X^T \Sigma_Y X (X^T X)^{-1}
$$

where

$$
(X^T X)^{-1} X^T =
\frac{3}{k(k-1)}
\begin{bmatrix}
    \frac{4}{k+1} & -1 \\
    -1 & \frac{2k+1}{6}
\end{bmatrix}
\begin{bmatrix}
    1 & \dots & k \\
    2 & \dots & 2
\end{bmatrix}
$$

and letting $i \geq j$,

$$
\begin{align}
(\Sigma_Y)_{ij}
&= \Cov(Y_i, Y_j) \\
&= \frac{2}{(n-i)(n-j)} [g(n, i, j) \sigma_S^4 + h(n, i, j) \sigma_V^4 ] + \frac{4j}{n-j}\sigma_S^2 \sigma_V^2\\
g(n, i, j) &= (n - i) j
    \left( \frac{(j+1)(2j+1)}{3} + (i-j-1) j \right) - \frac{1}{6} (j+1) j^2 (j - 1) \\
h(n, i, j) &= (4 + 2_{\{i = j\}}) (n-i) + 2 (n - i - j)
\end{align}
$$

Everything except the $\Sigma_Y$ can easily be derived with basic matrix operations. Note the above is valid only for $n > i + j$.

### Covariance matrix of $Y$ derivation

The derivation uses a number of lemmas:

$$
\begin{align}
(V_{t+1} - V_t)^2
&= (\epsilon_t + \dots + \epsilon_{t+i-1} + \eta_{t+i} - \eta_t)^2 \\
&= E_{ti}^2 + 2 E_{ti} H_{ti} + H_{ti}^2
\end{align}
$$

where $E_{ti} = \epsilon_t + \dots + \epsilon_{t+i-1}$ and $H_{ti} = \eta_{t+i} - \eta_t$ are jointly Gaussian. Also,

$$
\begin{align}
\Cov(E_{ti}, E_{sj}) &= |\{t, \dots, t+i-1\} \cap \{s, \dots, s+j-1\}| \sigma_S^2 \\
\Cov(H_{ti}, H_{sj}) &=
(1_{\{t=s\}} + 1_{\{t+i=s+j\}} - 1_{\{t+i=s\}} - 1_{\{t=s+j\}}) \sigma_V^2 \\
\Cov(E_{ti}, H_{sj}) &= 0 \\
\Cov(E_{ti}^2, E_{sj} H_{sj}) &= 0 \\
\Cov(H_{ti}^2, E_{sj} H_{sj}) &= 0 \\
\Cov(E_{ti} H_{ti}, E_{sj} H_{sj}) &= \Cov(E_{ti}, E_{sj}) \Cov(H_{ti}, H_{sj}) \\
&= \min(i,j) (1_{\{s=t\}} + 1_{\{s+j=t+i\}}) \sigma_S^2 \sigma_V^2
\end{align}
$$

If $A$ and $B$ are jointly Gaussian with mean 0 and covariance $c$ then $\Cov(A^2, B^2) = 2 c^2$. Putting all this together:

$$
\begin{align}
\Cov(Y_i, Y_j)
&= \Cov \left(
    \frac{1}{n-i} \sum_{t=1}^{n-i} (V_{t+i} - V_t)^2,
    \frac{1}{n-j} \sum_{s=1}^{n-j} (V_{t+j} - V_t)^2
\right) \\
&= \frac{1}{n-i} \frac{1}{n-j} \sum_{t=1}^{n-i} \sum_{s=1}^{n-j} \Cov(E_{ti}^2, E_{sj}^2) + \Cov(H_{ti}^2, H_{sj}^2) + 2 \Cov(E_{ti} H_{ti}, E_{sj} H_{sj})
\end{align}
$$

By carefully counting the terms in the double summation we arrive at the variance formula above.

## Optimal $k$



___

[^1]: This opens up another interesting question of how the rows might be weighted so as to give precedence for the smaller $i$ rows in computing the estimator.