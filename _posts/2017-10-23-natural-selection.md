---
layout: post
title: "Natural Selection"
description: ""
date: 2017-10-23
tags: [Bayesian stats]
comments: true
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

A few weeks ago I thought of an interesting and rather elegant application of Bayesian inference for population ecology. In high school biology, our class had spent a good amount of time studying Darwin's theory of evolution by natural selection. The skeleton of the theory is as follows: The survival and reproductive capability of an organism is encoded in a set of innate traits that are passed down from parent to child. For any population of organisms where there is a diversity of such traits, there are differences in the ability of individual organisms to pass on their genes, especially when some traits become critical to survival. As a result, the genes conferring greater reproductivity get passed on to the next generation, where they exist in a higher percentage than in the first. The distribution of traits in the second generation has changed from the first, and over many iterations of this process of natural selection, the genetic composition of the population can drastically change. In this post I will show a way to precisely quantify how fast that change occurs.

## Problem setup

Consider an initial population of $n_0$ organisms, each associated with a number ranging from 0 to 1 that reflects its probability of producing offspring. Suppose we sample an organism from the population and observe its probability of reproducing. Denote this value with the random variable $F_0$ for fitness, whose exact distribution will be decided later.

Let $S_0$ denote the number offspring of a randomly sampled organism. For now, let's set $S$ to be 0-1 Bernoulli random variable parametrized by F. In other words, a randomly sampled organism will have probability $F$ of having one child and probability $1 - F$ of having none. For simplicity, let's also assume that each offspring will have the same fitness score as its parent.
We can now use this framework to say things about, e.g., the probability that a random member of the original population reproduces ($\Prb(S_0 = 1)$), the distribution of fitness in the next population ($F_1$), and the distribution of how large that population is ($N_1$).

[comment]: # Remark about the abstraction from genome

## Bayesian inference

Because each organism can have at most one offspring and because child fitness equals parent fitness, it is the case that

$$ \Prb(F_1 \in \diff f) = \Prb(F_0 \in \diff f | S_0 = 1) $$

where $\Prb(X \in \diff x)$ is shorthand for the pdf of $X$ at $x$ multiplied by an infinitesimal, or equivalently $\diff \Prb(X \leq x)$.
Using Bayes' law the right hand side is proportional to

$$
\begin{align}
  & \Prb(S_0 = 1 | F_0 = f) \Prb(F_0 \in \diff f) \\
= & f \Prb(F_0 \in \diff f)
\end{align}
$$

Because these equations hold for any generation, this forms a recursive equation which can be "flattened out" to

$$ \Prb(F_n \in \diff f) = h(n) f^n \Prb(F_0 \in \diff f) $$

where $h(n)$ is a scaling function that does not depend on $f$.[^1]

This equation is identical to Bayesian inference of the mean of a Bernoulli random variable when all the samples are 1. The prior distribution is the previous generation's fitness distribution, the likelihood is $f$, and the posterior is the next generation's fitness distribution. The [Beta distribution](https://en.wikipedia.org/wiki/Beta_distribution) is typically used as the prior distribution of the Bernoulli parameter because the posterior is [also](https://en.wikipedia.org/wiki/Conjugate_prior) a Beta distribution. If we adopt a $\text{Beta}(\alpha, \beta)$ as the zeroth generation fitness distribution then the fitness distribution of the $n$th generation is

$$ F_n \sim \text{Beta}(\alpha + n, \beta) $$

## Random number of offspring

Now let's consider a more realistic model that allows for each organism to have more than one offspring. The distribution of the number of offspring is parametrized by the organism's fitness, which no longer is required to be in the interval $[0, 1]$. Again, I'll use $S(f)$ or
$S | F = f$ to denote the random number of children of an organism with fitness $f$.

Bayes' rule will not help us get the next generation fitness distribution when there are multiple offspring. We will need a different approach, which requires more precisely defining what a fitness distribution is.

Populations have a finite number of organisms, yet for mathematical simplicity we have pretended that they were of infinite size. Reverting back to the finite population framework, $\Prb(F_1 \in \diff f)$ represents the *number* of organisms in the 1st generation with fitness score in the vicinity of $f$ divided by the *total number* of organisms in the 1st generation. Letting $n_0$ be the population size of the *0th* generation, we can express the former as

$$ \sum_{i=1}^{n_0 \Prb(F_0 \in \diff f)} S(f) $$

which is the sum of numbers of children born to the $n_0 \Prb(F_0 \in \diff f)$ members of the 0th generation with fitness score near $f$. If we take the limit as $n_0$ approaches infinity of this quantity divided by $n_0$ we get

$$
\begin{align}
\lim_{n_0 \to \infty} \frac{1}{n_0} \sum_{i=1}^{n_0 \Prb(F_0 \in \diff f)} S(f)
&= \Prb(F_0 \in \diff f) \lim_{n_0 \to \infty} \frac{1}{n_0 \Prb(F_0 \in \diff f)} \sum_{i=1}^{n_0 \Prb(F_0 \in \diff f)} S(f) \\
&= \Prb(F_0 \in \diff f) \Expect[S | F = f]
\end{align}
$$

[almost surely](https://en.wikipedia.org/wiki/Law_of_large_numbers#Strong_law). Hence we can express the 1st generation fitness pdf as

$$
\begin{align}
\Prb(F_1 \in \diff f)
&= \frac{\Expect[S | F = f] \Prb(F_0 \in \diff f)}{\sum_f \Expect[S | F = f] \Prb(F_0 \in \diff f)} \\
&= \frac{\Expect[S | F = f] \Prb(F_0 \in \diff f)}{\Expect[S(F_0)]}
\end{align}
$$

A very intuitive result! This says that the fitness distribution of the first generation is that of the zeroth generation weighted by the relative average reproductive rates. A quick sanity check affirms that this is a valid pdf that is nonnegative everywhere and integrates to 1. Interestingly, the next generation fitness distribution only depends on the distribution of $S$ via its conditional expectation given $F$.

Let's define some concrete distributions to illustrate this result: let's assume the zeroth generation fitness distribution is uniform over $[0, 1]$
and that $\Expect[S | F = f] = f$,
which would be the case if $S | F \sim \text{Poisson}(F)$. Then

$$ \Prb(F_n \in \diff f) = \begin{cases}
(n + 1) f^n \diff f, & \text{if } f \in [0, 1] \\
0, & \text{otherwise}
\end{cases}
$$

*Proof*: $\Prb(F_1 \in \diff f) = \frac{f}{\Expect[F_0]} \Prb(F_0 \in \diff f) = 2 f \diff f$ checks out for the base case. For the inductive step,

$$
\begin{align}
\Expect[F_n]
&= \int_0^1 (n+1) f^{n+1} \diff f \\
&= \frac{n+1}{n+2}
\end{align}
$$

$$
\begin{align}
\Prb(F_{n+1} \in \diff f)
&= \frac{f}{\Expect[F_n]} \Prb(F_n \in \diff f) \\
&= \frac{f(n+2)}{n+1} (n+1) f^n \diff f \\
&= (n + 2) f^{n+1} \diff f
\end{align}
$$

$\square$

___

[^1]: Because the pdf must integrate to 1, $h(n) = 1 / \int_{-\infty}^\infty f^n \Prb(F_0 \in \diff f)$.
