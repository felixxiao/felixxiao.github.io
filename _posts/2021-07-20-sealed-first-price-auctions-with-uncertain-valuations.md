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

$$\Expect[f(X) | A] \Prb(A) = \Expect[f(X) \Prb(A|X)] \tag{2}$$

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

This is obtained by directly evaluating the integral $(1)$. Of course this distribution isn't ideal for our application, but as far as I am aware it is the only continous distribution with a closed form payoff that doesn't use indicator functions.

## Two-player Gaussian distribution

Now let's consider a more applicable model using Gaussian random variables for $X$. The special case of two players has an interesting result.

---

Suppose there are two sellers with estimation variances $\sigma_1^2$ and $\sigma_2^2$. Let

$$ \mu^* = \frac{1}{2} \sqrt{2 \pi (\sigma_1^2 + \sigma_2^2)} \tag{3} $$

Then $(\mu^* , \mu^* )$ is a Nash equilibrium.

---

Proof:

$$
\begin{align}
\Pi_1
&= \Expect[X_1 F_2^C(X_1)] \\
&= \Expect[X_1] - \Expect[X_1 F_2(X_1)] \\
\Expect[X_1 F_2(X_1)]
&= \int_{-\infty}^{\infty} x F_2(x) p_1(x) \diff x \\
&= \int_{-\infty}^{\infty} (\sigma_1 z + \mu_1) F_2(\sigma_1 z + \mu_1) p_1(\sigma_1 z + \mu_1) \sigma_1 \diff z \\
&= \sigma_1 \int_{-\infty}^{\infty} (\sigma_1 z + \mu_1) \Phi(\frac{\sigma_1}{\sigma_2} z + \frac{\mu_1 - \mu_2}{\sigma_2}) \phi(z) \diff z \\
\end{align}
$$

where $\phi(x)$ is the standard Gaussian PDF and $\Phi(x)$ the CDF. From here we need two lemmas.

$$
\begin{align}
\int_{-\infty}^{\infty} x \phi(x) \Phi(a + bx) \diff x &= \frac{b}{\sqrt{1 + b^2}} \phi(\frac{a}{\sqrt{1 + b^2}}) \\
\int_{-\infty}^{\infty} \phi(x) \Phi(a + bx) \diff x &= \Phi(\frac{a}{\sqrt{1 + b^2}})
\end{align}
$$

To prove the former, it suffices to verify that the integrand has the antiderivative

$$ \frac{b}{t} \phi(\frac{a}{t}) \Phi(x t + \frac{ab}{t}) - \phi(x) \Phi(a + bx) + C$$

where $t = \sqrt{1 + b^2}$.

Proof of the latter:

$$
\begin{align}
\int_{-\infty}^{\infty} \phi(x) \Phi(a + bx) dx
&= E \Phi(a + bX) \\
&= E P(Z \leq a + b X | X) \\
&= P(Z \leq a + b X) \\
&= P(Z - b X \leq a) \\
&= P(\mathcal{N}(0, 1 + b^2) \leq a) \\
&= P(\mathcal{N}(0, 1) \leq \frac{a}{\sqrt{1 + b^2}}) \\
&= \Phi(\frac{a}{\sqrt{1 + b^2}})
\end{align}
$$

where the third equality is an application of $(2)$. After some algebra, we get


$$\Pi_1 = \mu_1 - \frac{\sigma_1^2}{\sqrt{\sigma_1^2 + \sigma_2^2}} \phi(\frac{\mu_1 - \mu_2}{\sqrt{\sigma_1^2 + \sigma_2^2}}) - \mu_1 \Phi(\frac{\mu_1 - \mu_2}{\sqrt{\sigma_1^2 + \sigma_2^2}})$$

and similarly for $\Pi_2$. Letting $d = \frac{\mu_1 - \mu_2}{\sqrt{\sigma_1^2 + \sigma_2^2}}$, the first derivative of payoff for the two players can be written:

$$
\begin{align}
\frac{\partial \Pi_1}{\partial \mu_1} &= 1 - \frac{\sigma_1^2}{\sigma_1^2 + \sigma_2^2} d \phi(d) - \Phi(d) - \frac{\mu_1}{\sqrt{\sigma_1^2 + \sigma_2^2}} \phi(d) \\
\frac{\partial \Pi_2}{\partial \mu_2}
&= 1 - \frac{\sigma_2^2}{\sigma_1^2 + \sigma_2^2} (-d) \phi(-d) - \Phi(-d) - \frac{\mu_2}{\sqrt{\sigma_1^2 + \sigma_2^2}} \phi(-d) \\
&= \frac{\sigma_2^2}{\sigma_1^2 + \sigma_2^2} d \phi(d) + \Phi(d) - \frac{\mu_2}{\sqrt{\sigma_1^2 + \sigma_2^2}} \phi(d) \\
\end{align}
$$

(From the symmetry of the standard Gaussian distribution, we know $\phi(-x) = \phi(x)$ and $\Phi(-x) = 1 - \Phi(x)$).

Plugging $\mu^* $ into $\mu_1$ and $\mu_2$ (so $d = 0$) will make both derivatives zero. $\square$

---

Suppose $(\mu_1, \mu_2)$ is a Nash equilibrium pair of spreads. Then $\mu_1 = \mu_2 = \mu^* $.

---

Proof: $\frac{\partial \Pi_2}{\partial \mu_2} - \frac{\partial \Pi_1}{\partial \mu_1} = 2 d \phi(d) + 2 \Phi(d) - 1 = 0$ only if $d = 0$. Hence $\mu_1 = \mu_2$. To make both derivatives zero, they must equal $\mu^* $. $\square$

## Numerical integration and optimization

In the case of more than two players with arbitrary spreads and valuation variances, there is no closed form expression for the expected profit or its spread derivatives under a Gaussian distribution, so we must use numerical integration. Then we can use the evaluation of the derivatives in an optimization algorithm to compute the best strategy for one player given the others'.

Since we only care about the derivatives of an expected profit function with respect to the player's own spread, we can take advantage of the following fact:

$$ \frac{\partial \Pi_i}{\partial \theta_i} = E[X_i h(X_i) \frac{\partial}{\partial \theta_i}\log(p_i(X_i))] \tag{4} $$

where $h_i(x) = \prod_{j \neq i} F_j^C(x)$ (which doesn't include $\theta_i$). Note this is true for any continuous distribution of $X_i$ parameterized by $\theta_i$. Now suppose $X_i$ is Gaussian, so $\frac{\partial \log p}{\partial u} = \frac{x - u}{\sigma^2}$. Let $g_p(u) = E[X^p h(X)]$ so $\Pi = g_1(u)$. It is also true that for any natural number $p$,

$$\frac{\partial g_p}{\partial u} = \frac{g_{p+1} - g_p u}{\sigma^2}$$

So we can write derivatives of expected profit in terms of the "moments" and evaluate them too using numerical integration.

$$
\begin{align}
\frac{\partial g_1}{\partial u} &= \frac{g_{2} - g_1 u}{\sigma^2} \\
\frac{\partial^2 g_1}{(\partial u)^2} &= \frac{1}{\sigma^2} (\frac{1}{\sigma^2} (g_3 - 2 g_2 u + g_1 u^2) + g_1) \\
\end{align}
$$


## Perfect valuations

The following two statements are equivalent:

1. There are at least two players with estimation variance of zero (perfect valuations).

2. In all Nash equilibria, there are at least two optimal spreads of zero.

