---
layout: post
title: "The German Tank Problem"
description:
date: 2017-10-16
tags: estimation
comments: true
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

## Problem definition

The tank factory stamps each armored vehicle produced with an integer serial number starting from 1. The serial number of the last tank built equals the number of tanks manufactured in total, $n$.

From the wreckage of tanks captured or destroyed, Allied intelligence observes the serial numbers of a subset of $k$ tanks, among which the largest is denoted $M$. The challenge is to infer $n$, the total number of tanks produced by the Germans.

The maximum likelihood estimator of $n$ is $M$. This is because the likelihood as a function of $n$ is
$$
\mathcal{L}(n) = \begin{cases}
  0,         & \text{if } n < m, \\
  (1 / n)^k, & \text{if }n \geq m
\end{cases}
$$
so $n_{\text{MLE}}$ must be at least $m$ for the likelihood to be nonzero, and smaller $n$ means larger $(1 / n)^k$.

Intuitively, $M$ is a bad estimator because it will underestimate $n$ on average. To understand by how much it underestimates $n$, we need to look at the sampling distribution of $M$, which depends on _how_ the tanks are observed. There are two choices that allow for simple modeling:

1. The tanks do not change locations between observations, so it is impossible for the reconnaissance planes to observe the same tank twice. Moreover, no set of (distinct) tanks has any higher probability of being observed than any other set. In other words, it is not the case that tanks 4 and 5 are always observed together. This corresponds to random sampling without replacement.

2. Each serial number observed does not depend on any previous observations. This would be the choice if there was a protracted time between each tank observation that would allow the tanks to change location. This correponds to random sampling with replacement.

Allied statisticians chose the former as the more realistic model and used it to derive the classic minimum variance unbiased estimator for this problem. I'll reproduce the derivation in this post and then discuss the solution for the latter problem of sampling with replacement.

## Sampling without replacement

An unbiased estimator for $n$ can be derived from the expectation of $M$, the largest serial number.
Let $X_i$ be the $i$th observed serial number. The CDF of $M$ is, for $m \geq k$,

$$
\begin{align}
\Prb(M \leq m) &= \Prb(X_1 \leq m, \ldots, X_k \leq m) \\
  &= \Prb(X_1 \leq m) \Prb(X_2 \leq m | X_1 \leq m) \ldots \Prb(X_k \leq m | X_i \leq m \forall i = 1, \ldots, k - 1) \\
  &= \frac{m}{n} \cdot \frac{m-1}{n-1} \cdot \ldots \cdot \frac{m-k+1}{n-k+1} \\
  &= \dfrac{\frac{m!}{(m - k)!}}{\frac{n!}{(n - k)!}} \\
  &= \dfrac{\binom{m}{k}}{\binom{n}{k}} \\
\end{align}
$$

The probability mass function is

$$
\begin{align}
\Prb(M = m) &= \sum_{i=1}^k \Prb(M = m, X_k = M) \\
  &= k \Prb(M = m, X_1 = m) \\
  &= k \Prb(X_1 = m) \Prb(X_2 < m) \Prb(X_3 < m | X_2 < m) \ldots \Prb(X_k < m | X_i < m \forall i = 2, \ldots, k) \\
  &= k \frac{1}{n} \frac{m-1}{n-1} \ldots \frac{m-k+1}{n-k+1} \\
  &= \frac{k}{m} \dfrac{\frac{m!}{(m - k)!}}{\frac{n!}{(n - k)!}} \\
  &= \frac{k}{m} \dfrac{\binom{m}{k}}{\binom{n}{k}}
\end{align}
$$

Using the binomial identity $\sum_{m=k}^n \binom{m}{k} = \binom{n+1}{k+1}$, the expectation of $M$ is

$$
\begin{align}
\Expect[M] &= \sum_{m=k}^n m \Prb(M = m) \\
  &= \sum_{m=k}^n k \dfrac{\binom{m}{k}}{\binom{n}{k}} \\
  &= \dfrac{k}{\binom{n}{k}} \sum_{m=k}^n \binom{m}{k} \\
  &= \dfrac{k}{\binom{n}{k}} \binom{n+1}{k+1} \\
  &= \frac{k (n - k)! k! (n + 1)!}{n! (k + 1)! (n - k)!} \\
  &= \frac{k (n + 1)}{(k + 1)} \\
\end{align}
$$

If we plug in the sample $M$ and rearrange for $n$ we get the unbiased estimator

$$ \hat{N} = M \left(1 + \frac{1}{k}\right) - 1 $$

## Sampling with replacement

To derive the sampling distribution of $M$ (assuming sampling with replacement), we start with the CDF which is $\Prb(M \leq m) = \Prb(X_1 \leq m, \ldots, X_k \leq m) = (\frac{m}{n})^k$ for $0 \leq m \leq n$.
From this we get the probability mass function which is
$\Prb(M = m) = \Prb(M \leq m) - \Prb(M \leq m - 1) = (\frac{m}{n})^k - (\frac{m-1}{n})^k$ for $1 \leq m \leq n$.

As it turns out, the expectation of $M$ doesn't have a closed form, and the nicest form I could get it in requires a probability trick. Since $M$ is always positive,

$$
\begin{align}
M &= \sum_{m=0}^{\infty} \mathbb{1}_{\{M>n\}} \\
\Expect[M] &= \sum_{m=0}^{\infty} \Expect[\mathbb{1}_{\{M>n\}}] \\
  &= \sum_{m=0}^{\infty} \Prb(M > n) \\
  &= \sum_{m=0}^n 1 - \left(\frac{m}{n}\right)^k \\
  &= n + 1 - \sum_{m=0}^n \left(\frac{m}{n}\right)^k
\end{align}
$$

No [closed form](http://mathworld.wolfram.com/PowerSum.html) exists for the summation.
