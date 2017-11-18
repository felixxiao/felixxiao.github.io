---
layout: post
title: "Natural Selection pt. 3"
description: ""
category: 
tags: []
---

$$
  \newcommand{\R}{\mathbb{R}}
  \newcommand{\Prb}{\mathbb{P}}
  \newcommand{\Expect}{\mathbb{E}}
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

This post continues the discussion of estimating the distribution of unobserved fitness scores from data about the population. To recap our model used in the previous post, we have an initial population of $n_0$
organisms each with a fitness score $F_0$ drawn independently from a $\text{Gamma}(\alpha, \beta)$ distribution. An organism with fitness $f$ will produce $S(f)$ offspring, independent of everything else, and distributed according to $\text{Poisson}(f)$. Each of its offspring will also have fitness score $f$.

We have shown that the fitness score of a randomly sampled individual in the $i$th generation converges in distribution to $\text{Gamma}(\alpha + i, \beta)$ as $n_0 \to \infty$. In the previous post, we have derived maximum likelihood estimation procedures for $\alpha$ and $\beta$ with two different cuts of the data:

1. We have genealogy data for each generation's organisms, which tells us which family each organism belongs to as well as how many children it has.

2. We have offspring counts of randomly sampled individuals in various generations but do not know which individuals share the same fitness score.

Neither of these two cases yield closed-form solutions to the MLE estimators. However, the log-likelihood in both cases is at least twice differentiable in the first quadrant, which allows for a wide variety of numerical optimization methods to compute the estimates. Having solved the estimation problem in cases of relatively detailed data, we will now turn our attention to estimating the 0-generation fitness parameters with the lowest level of data granularity, *i.e.* when we only have the sizes of populations $0$ to $k$.

Let $N_i$ be the random size of the $i$th generation,
and $n_i$ its realization.
Let $\mathcal{F}_i$ be all the information, observed and unobserved, about the $i$th population.
This contains $n_i$, the fitness scores of each individual, and
$p_i(f)$ -- the proportion of individuals with fitness score $f$.
Given all the information about the previous generation, the size of population $i+1$ is the sum of Poisson random variables:

$$
\begin{align*}
N_{i+1} | \mathcal{F}_i
&= \sum_{f \in \mathcal{F}_i} \sum_{j = 1}^{n_i p_i(f)} S(f) \\
&= \text{Pois}\Big( n_i \sum_{f \in \mathcal{F}_i} f p_i(f) \Big) \\
&\approx \text{Pois}( n_i \Expect F_i ) \\
&= \text{Pois}\Big(n_i \frac{\alpha + i}{\beta}\Big)
\end{align*}
$$

By substituting the sample average fitness of generation $i$ with its large sample approximation, we can express $N_{i+1}$
in terms of just the data we have $n_i$ and the parameters we try to estimate. From this we can derive the likelihood

$$
\begin{align*}
\mathcal{L}(\alpha, \beta)
&= \Prb(N_1, ..., N_k | N_0) \\
&= \Prb(N_k | N_{k-1}, ..., N_0) ... \Prb(N_1 | N_0) \\
&= \prod_{i=1}^k \Prb(N_k | N_{k-1}) \\
&= \prod_{i=0}^{k-1} \left( n_i \tfrac{\alpha + i}{\beta} \right)^{n_{i+1}} e^{-n_i \frac{\alpha + i}{\beta}} \frac{1}{n_{i+1}!}
\end{align*}
$$

The first order derivatives of the log-likelihood are:

$$
\begin{align*}
\frac{\partial l}{\partial \alpha} &= \sum_{i=1}^{k-1} \frac{n_{i+1}}{\alpha + i} - \frac{n_i}{\beta} \\
\frac{\partial l}{\partial \beta} &= \sum_{i=1}^{k-1} n_i \frac{\alpha + i}{\beta^2} - \frac{n_{i+1}}{\beta}
\end{align*}
$$
