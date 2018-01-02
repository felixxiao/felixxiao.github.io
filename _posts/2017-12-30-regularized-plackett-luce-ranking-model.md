---
layout: post
title: "Regularized Plackett-Luce Ranking Model"
description: "How to gamble in horse racing"
category: 
tags: []
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

Suppose you are betting money in a horse race, and you make or lose money depending on whether you can correctly predict which horse will win first, second, or third place. Suppose further you have access to rankings from past races, e.g.:

| Race | Horse | Rank |
| ---- | ----- | ---- |
| A    | 01    | 3    |
| A    | 02    | 1    |
| A    | 03    | 2    |
| B    | 02    | 3    |
| B    | 04    | 1    |

... and so on. How would use this data to predict the results of the upcoming race?

The Plackett-Luce model (or generalized Bradley-Terry model) is the most widely studied probabilistic model for rank data, favored for its mathematical simplicity and intuitive theory. The model is derived from Luce's choice axiom, which requires that the probability of selecting one item over another is independent of all other alternatives:

$$ \frac{P(i | A)}{P(i | A)} = \frac{P(i | B)}{P(j | B)} $$

for any sets $A$ and $B$ of available choices containing $i$ and $j$. The PL model is consistent with this axiom and states that every item in the universe of choices is associated with a nonnegative strength number $w_i$. The probability of choosing item $i$ among alteratives $A \in i$ is

$$ P(i | A) = \frac{w_i}{\sum_{j \in A} w_j} $$

The extension to rank data can be conceptualized thus. The model generates the item in first place by the above distribution. For the second place item, we choose among the alternatives excluding the already-chosen item in first place. The remaining items are drawn in the same manner by excluding from the set of choices items already ranked higher. The joint probability of a ranking $i_1 \succ i_2 \succ ... \succ i_m$ is

$$ \Prb(i_1 \succ i_2 \succ ... \succ i_m | A) =
  \frac{w_{i_1}}{w_{i_1} + w_{i_2} + ... + w_{i_m}}
  \frac{w_{i_2}}{w_{i_2} + ... + w_{i_m}}
  \cdots
  \frac{w_{i_{m-1}}}{w_{i_{m-1}} + w_{i_m}} $$

Hence a single ranking of $m$ items consists of $m-1$ comparisons. Let our data $\mathcal{D} = \\{(i_l, A_l)\\}_{l=1}^n$ be a collection of $n$ comparisons, where the $i_l$ indices the winner of the $l$th comparison and $A_l$ represents the set of alternatives which contains $i_l$. The log-likelihood of this data with nonnegative strengths $w$ is

$$ \log \mathcal{L}(w) = \sum_{l=1}^n \left( \log w_{i_l} - \log \sum_{j \in A_l} w_j \right) $$

There is no closed-form expression for the strengths that maximize the likelihood. The log-likelihood function is concave after a reparametrization with $w_i = \exp(\theta_i)$. There are two types of iterative algorithms for finding its global maximizer:

1. The well-established [Minorization-Maximization](https://en.wikipedia.org/wiki/MM_algorithm) algorithm presented in [Hunter 2004](https://projecteuclid.org/euclid.aos/1079120141). It operates by using a local approximation function to the log-likelihood that has an analytic maximizer and alternates between constructing an approximation function at the current point and maximizing it to get the next point.

2. A more recent spectral algorithm based on the stationary distribution of a Markov chain illustrated in [Maystre and Grossglauser 2015](https://papers.nips.cc/paper/5681-fast-and-accurate-inference-of-plackettluce-models)

Both of these two varieties of algorithms work well when the comparison data meets a particular regularity condition -- Assumption 1 of Hunter 2004:

> In every possible partition of the [set of choices] into two nonempty subsets, some [item] in the second set beats some [item] in the ﬁrst set at least once.

For instance, if there is a strict subset of the choices that are never compared with any choices outside of themselves, then Plackett-Luce cannot infer strength parameters for the whole set. Or, if there are any items that always get chosen in every single comparison or never win any comparison, then this assumption is also violated.

This turns out to be inconvenient for many real-world data sets. The 2002 NASCAR data set that Hunter 2004 uses, for instance, violates the assumption -- the set of drivers had to be cut down to make inference work. [Guiver and Snelson 2009](https://www.microsoft.com/en-us/research/publication/bayesian-inference-for-plackett-luce-ranking-model/) proposed a Bayesian version of Plackett-Luce that handles this issue. Using a Gamma-distributed prior on the choice strengths that pulls their values to be closer to one another, their method no longer requires the regularity condition.

In this blog post, I will present a frequentist version of a regularized Plackett-Luce model that extends the MM algorithm of Hunter 2004. The extended algorithm, which to my knowledge does not yet exist in the literature, maximizes the PL log-likelihood plus a tunable penalty term that encourages the optimal choice strengths to be close to one another. The penalized log-likelihood is

$$ \mathcal{l}_\alpha (w) = \sum_{l=1}^n \left( \log w_{i_l} - \log \sum_{j \in A_l} w_j \right) +\alpha \sum_{k=1}^m \log w_k $$

for $\alpha \geq 0$ and $ \sum_{k=1}^m w_k = 1 $.

Without the equality constraint on the $w_k$, the problem becomes unbounded -- the objective can always be increased by multiplying the current $w$ by a scalar greater than 1. To see why the $\log w_k$ term encourages the strengths to be closer together, see the below plot of the level sets of $\log w_1 + \log w_2$ for $m = 2$. The line segment represents the feasible set of $w$ given by the equality and nonnegative constraints. The level sets that intersect the line segment closer to its mid point have higher function values.

It is interesting to note that this penalized log-likelihood formulation is equivalent to Bayesian inference of the Plackett-Luce model with a symmetric Dirichlet prior on the choice strengths. The penalty term is also concave thanks to $\alpha \geq 0$, so following the technique in Hunter 2004, we obtain the following minorization of the above:

$$ Q(w | w^{(k)}) = \sum_{l=1}^n \left( \log w_{i_l} - \frac{\sum_{j \in A_l} w_j}{\sum_{j \in A_l} w_j^{(k)}} \right) + \alpha \sum_{k=1}^m \log w_k $$

It is easy to see that $Q(w \| w^{(k)})$, plus a constant not depending on $w$, satisfies the required properties.

1. $ Q(w^{(k)} \| w^{(k)}) + C = \mathcal{l}_\alpha (w^{(k)}) $
2. $ Q(w \| w^{(k)}) + C \leq \mathcal{l}_\alpha (w) $

Taking partial derivatives of the Lagrangian of $Q(w \| w^{(k)})$ with respect to $w_k$ gets us the update equation

$$ w_k^{(k+1)} = \frac{c_j + \alpha}{d_j^{(k)} + \lambda} $$

where $\lambda$ is the Lagrange multiplier associated with the sum-to-one constraint, $c_j$ is the number of comparisons item $j$ won in and

$$ d_j^{(k)} = \sum_{l : j \in A_l} \frac{1}{\sum_{i \in A_l} w_i^{(k)}} $$

To find the value of the Lagrange multiplier, apply the sum-to-one constraint and use Newton's method to find a positive root of

$$ g(\lambda) = - 1 + \sum_{j=1}^m \frac{c_j + \alpha}{d_j^{(k)} + \lambda} $$

with $\lambda$ initialized to $\max(0, n + m \alpha - \sum_{j=1}^m d_j^{(k)} / m)$.

[comment]: <> *Claim*: The Lagrange multiplier above is always nonnegative at a critical point.

[comment]: <> Proof:

## Results