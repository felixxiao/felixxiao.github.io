---
layout: post
title: "Varying Coefficient Models Via Gaussian Process Priors"
description: ""
category: 
tags: [time series, Bayesian stats]
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

Let us consider the classic problem of univariate regression. The response variable $Y$ depends on the predictor via some unknown function and a noise term $\epsilon$ that is independent of $X$:

$$Y = f(X) + \epsilon$$

To infer the structure of $f$ from observed pairs of $X$ and $Y$, we need to specify two things: a restricted class of functions in which to seek $f$ and a loss function that evaluates different members of that function class. The most common choice for the former is the class of affine functions $f(x) = a + b x$, and for the latter a squared error loss -- together these two constitute simple linear regression.

The time series version of this prediction problem adds an implicit ordering to the data points, a feature that is typically ignored on the assumption that $f$ is invariant across time. However, this assumption is unrealistic in many practical applications. For example in my domain of finance, there is a well-known phenomenon of "alpha decay", in which a signal that has predictive power on asset returns loses its predictability over time. Another example can be found in ecology, where environmental changes may lead to a breakdown or reversal of previously observed patterns.

To model this complexity, statisticians have generalized the static linear model into a varying coefficient model. Older work on such models typically impose some semi-parametric structure on the coefficients, for instance local-linearity ([Fan and Zhang 1999](https://projecteuclid.org/euclid.aos/1017939139)). Here, I will explore a different framework for the varying coefficient model that uses a Gaussian process as a Bayesian prior.

A Gaussian process is an infinite, ordered collection of Gaussian random variables, any finite subset of which is jointly Gaussian. A Gaussian process can be fully characterized by its mean $\mu(t)$ and covariance $k(s, t)$ functions, which take in the time-index and output the first two moments of the corresponding random variables in the set. In most statistics applications, $\mu(t) = 0$ and $k(s, t)$ is a function that depends only on $s - t$. In Gaussian process regression, the index is not time but rather the values of the covariates corresponding to the response. Some parametric function called the kernel is chosen for $k$ and the prediction interval for a new input $x'$ is a Gaussian distribution conditioned on all the observed $(x, y)$.

The approach here is different in that it uses a Gaussian process as a prior model for the coefficient. For simplicity, let's consider the no-intercept case:

$$Y_t = b_t X_t + \epsilon_t$$

$$\{b_t\} \sim \mathcal{GP}(0, K)$$

$$\epsilon_t \sim \mathcal{N}(0, \sigma^2)$$

where $K : \R^n \times \R^m \mapsto \R^{n \times m}$ is a positive-semidefinite function that takes in two (sets of) time-indices and outputs the prior covariance (matrix) of $b$ between them.
Let $\tau$ be the set of time-indices for which we observe data and let $\upsilon$ be the set of time-indices we are inferring $b$ for. Let $x$ and $y$ be the vectors of $X$ and $Y$ realizations and let $D : \R^n \mapsto \R^{n \times n}$ map a vector to a diagonal matrix. The joint distribution of $(b_\upsilon, y_\tau)$ is a multivariate Gaussian with mean vector $\[0 0\]$ and covariance matrix:

$$
\begin{bmatrix}
K(\upsilon, \upsilon) & K(\upsilon, \tau) D(x) \\
D(x) K(\tau, \upsilon) & D(x) K(\tau, \tau) D(x) + I \sigma^2
\end{bmatrix}
$$

And the posterior Gaussian of $b_\upsilon$ has mean and covariance:

$$
\begin{align}
\Expect[b_\upsilon | x, y] &= K(\upsilon, \tau) D(x) \left[ D(x) K(\tau, \tau) D(x) + I_n \sigma^2 \right]^{-1} y \\
\Cov(b_\upsilon | x, y) &= K(\upsilon, \upsilon) - K(\upsilon, \tau) D(x) \left[D(x) K(\tau, \tau) D(x) + I_n \sigma^2\right]^{-1} D(x) K(\tau, \upsilon)
\end{align}
$$

Suppose we allow no time-variance so that $K(s,t) = \theta^2$. Then, this is equivalent to ridge regression with $\lambda = \frac{\sigma^2}{\theta^2}$.
