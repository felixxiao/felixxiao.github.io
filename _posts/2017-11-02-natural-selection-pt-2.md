---
layout: post
title: "Natural Selection pt. 2"
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

This post is an extension of the previous one on modeling the dynamics of population fitness over generations. We will focus on the problem of inferring fitness distribution when different aspects of the data are known, starting from the easiest and most complete data to the hardest and most obfuscated.

## Inferrence with genealogies

In our first case, the dataset consists of the genealogical trees of $n$ individuals sampled from the initial population and their descendents up to some number of generations.

Let $S_{ij}$ represent the number of offspring of the $j$th member of the $i$th family.
Let $m_i$ represent the number of individuals in each family tree with known number of children.

The joint likelihood can be written as

$$
\begin{equation}
\mathcal{L} = \prod_{i=1}^n p(F_i = f_i) \prod_{j=1}^{m_i} \Prb(S_{ij} = s_{ij} | F_i = f_i)
\end{equation}
$$

since the $S_{ij}$ are independent of one another given $F_i$. The $F$ in this problem are latent variables. In the case of known genealogies, the latent $F$ is relatively easy to deal with because each $F_i$ is associated with a known partition of the data $S_{ij}$.[^1] A naive way of estimating the distribution of $F$ would be to first get point estimates of the $F_i$ of each family and then use the estimates themselves to back out the distribution. The weakness of this approach is that it doesn't account for the different precision of each point estimate, and would weigh a less precise fitness score estimate from a small family the same as the more precise estimate from a large one. To fit our model, it would be better to maximize the *marginal* likelihood $\Prb(S)$ of the data we observe.

In most cases, the marginal density of a compound distribution does not have an analytic form, as the integral over all the latent variables is too complicated.[^2] However, if one specifies the latent variable to be the conjugate prior of the observed, then the marginal does have a closed form. Therefore we will specify the generative distributions as follows. For number of offspring given fitness we'll use Poisson, and for fitness we'll use its conjugate prior which is Gamma.

$$
S_{ij} | F_i \sim \text{Poisson}(F_i) \\
F_i \sim \text{Gamma}(\alpha, \beta)
$$

The joint density of one family is

$$
p(f_i) \prod_{j=1}^{m_i} p(s_{ij} | f_i)
= \frac{\beta^\alpha}{\Gamma(\alpha)} f_i^{\alpha - 1 + s_i} e^{-(\beta + m_i) f_i} \frac{1}{\prod_{j=1}^{m_i} s_{ij} !}
$$

where $s_i = \sum_{j=1}^{m_i} s_{ij}$.
The marginal likelihood is

$$
\begin{align*}
p(s)
&= \int \cdots \int \prod_{i=1}^n p(f_i) \prod_{j=1}^{m_i} p(s_{ij} | f_i) \diff f_1 \cdots \diff f_n \\
&= \prod_{i=1}^n C_i \frac{\beta^\alpha}{\Gamma(\alpha)} \int_0^\infty f_i^{\alpha + s_i - 1} e^{-(\beta + m_i) f_i} \diff f_i \\
&= \prod_{i=1}^n C_i \frac{\beta^\alpha}{\Gamma(\alpha)} \frac{\Gamma(\alpha + s_i)}{(\beta + m_i)^{\alpha + s_i}} \\
\end{align*}
$$

where $C_i = 1 / \prod_{j=1}^{m_i} s_{ij}!$ can be disregarded as it doesn't depend on $\alpha$ or $\beta$. The first partial derivatives of the log marginal likehood are

$$
\begin{align*}
\frac{\diff l}{\diff \alpha} &=
\sum_{i=1}^n \log \beta - \psi_0(\alpha) + \psi_0(\alpha + s_i) - \log(\beta + m_i) \\

\frac{\diff l}{\diff \beta} &=
\sum_{i=1}^n \frac{\alpha}{\beta} - \frac{\alpha + s_i}{\beta + m_i} \\
\end{align*}
$$

where $\psi_0(\alpha) = \frac{\diff}{\diff \alpha} \log \Gamma(\alpha)$ is the [digamma function](https://en.wikipedia.org/wiki/Digamma_function). There is no analytic solution for the optimal parameters, so we must use an algorithm like gradient ascent or Newton's method to maximize the function. Newton's method, which uses second-order derivatives, converges faster particularly when the number of parameters is small, like in this case. The update equations are

$$
\begin{align}
\alpha_{t+1} &= \alpha_t - \left( l_\alpha^{(t)} l_{\beta \beta}^{(t)} - l_\beta^{(t)} l_{\alpha \beta}^{(t)} \right) / \text{det}(\mathbf{H}^{(t)}) \\
\beta_{t+1} &= \beta_{t+1} - \left( l_\beta^{(t)} l_{\alpha \alpha}^{(t)} - l_\alpha^{(t)} l_{\alpha \beta}^{(t)} \right) / \text{det}(\mathbf{H}^{(t)}) \\
\end{align}
$$

where $l_\alpha^{(t)}$ is the partial derivative of the log likelihood with respect to $\alpha$, evaluated with $\alpha_t$ and $\beta_t$) and $\mathbf{H}^{(t)}$ is the Hessian matrix of $l^{(t)}$.

## Inference with offspring counts

Now suppose we don't have genealogical information, and the only data we have access to is the offspring counts of randomly sampled, unrelated individuals from various generations. Like before, we'll use $s_{ij}$ to denote the number of offspring born to the $j$th individual sampled from generation $i$. Let $\mathcal{I}$ be the set of generations, not necessarily consecutive, from which we have data. We'll use $n$ for the cardinality of $\mathcal{I}$ and $m_i$ the number of individuals from generation $i$th in our dataset.

First we need to model how the distribution of $F_i$ changes over each generation given a starting distribution. As before, we'll let the starting generation fitness $F_0$ be a gamma distribution parametrized by $\alpha$ and $\beta$, and the number of offspring born given fitness $f$ be $\text{Poisson}(f)$. Using a result from the previous post, we get a recursive equation for the distribution of the next generation in terms of the distribution of the previous:

$$ \Prb(F_{n+1} \in \diff f) = \frac{f}{\Expect F_n} \Prb(F_n \in \diff f) $$

*Claim:* $F_n \sim \text{Gamma}(\alpha + n, \beta)$.

*Proof:* The base case of $n = 0$ is trivially true. For the inductive step:

$$
\begin{align}
\frac{f}{\Expect F_n} \Prb(F_n \in \diff f)
&= \frac{f}{\frac{\alpha + n}{\beta}} \frac{\beta^{\alpha + n}}{\Gamma(\alpha + n)} f^{\alpha + n - 1} e^{-\beta f} \diff f \\
\Prb(F_{n+1} \in \diff f)
&= \frac{\beta^{\alpha + n + 1}}{\Gamma(\alpha + n + 1)} f^{\alpha + n} e^{-\beta f} \diff f \\
\end{align}
$$

which is the pdf of $\text{Gamma}(\alpha + n + 1, \beta)$ $\square$.

Armed with this information we can derive the marginal distribution of $S_{ij}$:

$$
\begin{align}
\Prb(S_{ij} = s)
&= \int_0^\infty p(s, f) \diff f \\
&= \int_0^\infty \frac{1}{s!} \frac{\beta^{\alpha + i}}{\Gamma(\alpha + i)} f^{\alpha + i + s - 1} e^{-(\beta + 1) f} \diff f \\
&= \frac{1}{s!} \frac{\beta^{\alpha + i}}{\Gamma(\alpha + i)} \frac{\Gamma(\alpha + i + s)}{(\beta + 1)^{\alpha + i + s}}
\end{align}
$$

The marginal likelihood is maximized when $\alpha$ and $\beta$ satisfy:

$$
\begin{equation}
\sum_{i \in \mathcal{I}} \sum_{j=1}^{m_i} \log \beta - \psi_0(\alpha + i) + \psi_0(\alpha + i + s_{ij}) - \log(\beta + 1) = 0 \\

\beta = \frac{\alpha \sum_{i \in \mathcal{I}} m_i (1 + i)}{\sum_{i \in \mathcal{I}} \sum_{j=1}^{m_i} s_{ij}} 
\end{equation}
$$

___

[^1]: In an opposite case, consider a Gaussian mixture model where each data point is generated by first picking one of $n$ different Gaussian distribution and sampling from that distribution. It is not known whether any two data points are sampled from the same Gaussian.

[^2]: In such cases the approach would be to use an [EM algorithm](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Description)
