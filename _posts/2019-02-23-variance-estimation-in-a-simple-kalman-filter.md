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

Now we have one equation and two unknowns. The trick is taking the lag-2 difference of the series to derive a second, independent equation:

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
    \sigma_S^2 \\
    \sigma_V^2
\end{bmatrix}
=
\begin{bmatrix}
    Y_1 \\
    Y_2
\end{bmatrix}
$$
