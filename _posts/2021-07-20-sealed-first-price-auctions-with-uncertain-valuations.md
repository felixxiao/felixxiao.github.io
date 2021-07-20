---
layout: post
title: "Sealed first-price auctions with uncertain valuations"
description: ""
category:
tags: [auctions]
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

Consider a market with one buyer and $k$ sellers of a item whose true value is unknown until after the transaction. The buyer requests quotes from each seller and is constrained to buy from the one with the lowest price. The sellers each have their private estimate of the item's value, generated from estimators which are independent, unbiased, and have finite variances known to all.

The sellers control what price they quote back to the buyer, which is obtained by taking their estimate and adding a constant ("spread"). Their payoff is the difference between their quoted price and the item's intrinsic value if the buyer trades with them, and zero otherwise. Let $X_i$ be this difference for seller $i$, with mean equal to the spread. The buyer buys from $i$ when $X_i$ is the minimum of all the $X_j$. Hence the expected payoff is written:

$$\Pi_i = \Expect[X_i | X_i \leq X_j \forall j] \Prb(X_i \leq X_j \forall j)$$

We can re-write this expression in a nicer way:

$$\Pi_i = \Expect[X_i \prod_{j \neq i} F^C_j(X_i)] \tag{1}$$

where $F^C_j(x) = 1 - F_j(x) = \Prb(X_j \gt x)$ is the complementary CDF of $X_j$. It factors into marginal CDFs because we assumed the $X_j$ to be independent.

The reason these two expressions are equal is that, for any random variable $X$, function $f$, and event $A$:

$$\Expect[f(X) | A] \Prb(A) = \Expect[f(X) P(A|X)] \tag{2}$$

Proof:

$$
\begin{align}
\Expect[f(X) | A] \Prb(A)
&= \int f(x) \Prb(X \in \diff x | A) \Prb(A) \\
&= \int f(x) \Prb(X \in \diff x \cup A) \\
&= \int f(x) \Prb(X \in \diff x) \Prb(A | X = x) \\
&= \Expect[f(X) \Prb(A|X)] \\
\end{align}
$$

With objective functions $(1)$ that define individual payoffs, this $k$-player game can be classified as a first-price sealed-bid auction. First-price means the lowest seller receives the amount she offered. Sealed-bid means the quoted prices are unknown to other sellers. The setup differs from the classic description of this auction in that the player payoffs are with respect to a common intrinsic valuation rather than their private valuation. Private valuations are noisy estimates of the intrinsic value, which creates an adverse selection effect: given you are chosen by the buyer, you were more likely to have undervalued the item and lowballed yourself. The reason each seller charges a spread on top of their private valuation is to mitigate this adverse selection; intuition says the more uncertain our estimate, the larger our spread should ideally be.

The rest of the post will dive into the problem of what each seller's optimal $\mu_i = \Expect[X_i]$ should be, given all the other sellers' spreads and variances. For simplicity we will assume all the $X_i$ belong to the same type of distribution with possibly different parameters. In a small number of cases I will show there is an analytical solution to the optimal spread for each seller that will enable us to find Nash equilibria strategies. In the general case where there is no analytic solution, even for simple distributions, I will describe a fast Newton algorithm for computing the optimal spread.

## Exponential distribution

In the case where each $X_i$ is an exponential random variable with rate parameter $\lambda_i$ the payoff function has an analytic form:

$$ \Pi_i = \frac{\lambda_i}{(\sum_j \lambda_j)^2} $$

## Two-player Gaussian distribution

Suppose there are two sellers with estimation variances $\sigma_1^2$ and $\sigma_2^2$. Let

$$\mu^* = \frac{1}{2} \sqrt{2 \pi (\sigma_1^2 + \sigma_2^2)}$$

Then $(\mu^* , \mu^* )$ is a Nash equilibrium pair of spreads.

Suppose $(\mu_1, \mu_2)$ is a Nash equilibrium pair of spreads. Then $\mu_1 = \mu_2 = \mu^* $.

## Perfect valuations

The following two statements are equivalent:

1. There are at least two players with estimation variance of zero (perfect valuations).

2. In all Nash equilibria, there are at least two optimal spreads of zero.
