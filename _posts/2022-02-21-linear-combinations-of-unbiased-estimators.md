---
layout: post
title: "Linear Combinations of Unbiased Estimators"
description: ""
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

Meta-analysis is a procedure commonly used in the natural sciences that seeks to combine the results of multiple studies on the same research question to produce a better model than what each of the studies could achieve alone. Although meta-analysis deals with a wide array research problems and encompasses a large body of statistical tools, I want to begin today's discussion with a simple example. Suppose there were two studies A and B that have each produced estimates on an unknown quantity $\mu$ and reported values of $A$ and $B$ along with standard errors $\sigma_A$ and $\sigma_B$. Suppose further that these two estimators (random variables) are unbiased (have expectations equalling $\mu$) and, perhaps due to using different data, are independent. How would a meta-analysis combine them to produce the best possible estimator? What does "best" mean?

## Minimum variance unbiased linear combination

We can formulate the problem in terms of optimization. The goal is to find values for the weights of a linear combination $w_A A + w_B B$ that satisfy some conditions. Here are two reasonable ones.

1. The linear combination is unbiased. Since the components are unbiased, this constraint reduces to $w_A + w_B = 1$.

2. Subject to being unbiased, the variance of the linear combination, which is expressed as $w_A^2 \sigma_A^2 + w_B^2 \sigma_B^2$, is the lowest it can be.

A little calculus will show that the optimal weights are proportional to the inverse of their respective estimator variances, and normalized to sum to 1.

## $k$ estimators with known covariance

Generalizing the problem to a vector of $k$ estimators that are potentially correlated, and a known covariance matrix $\Sigma$, the optimization problem can be written as follows (let $\boldsymbol{1}$ denote a $k$-length vector of ones)

$$
\newcommand{\ones}{\boldsymbol{1}}
\begin{align}
    \min && w' \Sigma w \\
    \text{s.t.} && w' \boldsymbol{1} = 1
\end{align}
$$

This is solved using the Lagrange multiplier method to give the optimal weights of

$$ w^* = \frac{\Sigma^{-1} \boldsymbol{1}}{\boldsymbol{1}' \Sigma^{-1} \boldsymbol{1}} $$

for which the objective function, the variance of the linear combination, attains its minimum of

$$ \frac{1}{\ones' \Sigma^{-1} \ones} $$

The inverse of the covariance matrix is sometimes called the *precision* matrix. Take the sum of all the entries of the precision matrix. The reciprocal of this "total precision" is the minimum variance of the unbiased linear combination.

We made an assumption in this solution that the covariance matrix is invertible. Now let's deal with the case where that does not hold.

## Singular covariance matrix

The covariance matrix is a real symmetric positive semidefinite matrix that by the spectral theorem can be diagonalized into an orthogonal basis of eigenvectors and nonnegative eigenvalues. When the matrix is singular then some of the eigenvalues are zero and the corresponding eigenvectors form a basis for the kernel of the matrix, the subspace of vectors that map to the zero vector.

Because the eigenvectors $v$ form a basis for $\R^k$, any weight vector can be expressed in terms of this basis

$$ w = a_1 v_1 + ... + a_k v_k $$

and the variance of the linear combination is equal to

$$ \lambda_1 a_1^2 + ... + \lambda_k a_k^2 $$

Suppose $\lambda_k = 0$ in the singular case. Then if we allocate $w$ purely to $v_k$ and set $a_i = 0$ for all $i \neq k$, the variance of the linear combo becomes zero. However, we still need to deal with the weights-sum-to-one constraint. This is doable if there exists coefficient $a$ such that $a {v_k}' \ones = 1$, which is equivalent to saying that $v_k$ and $\ones$ are not orthogonal.

More generally, if the ones vector is not orthogonal to the kernel of $\Sigma$, then the optimal weight lies in the kernel and gives a variance of zero in the combined estimators. In these special cases it is possible to deduce the unknown quantity exactly.

This was a bit counterintuitive to me so I'll provide an example. Let there be two estimators with unit variance and negative one covariance. The weight vector of $\begin{matrix} [\frac{1}{2} & \frac{1}{2}] \end{matrix}$ gives 1) an unbiased estimator 2) with zero variance. The intuition is that these two estimates perfectly mirror each other around the true quantity, which is just located at their midpoint.

What about the case where the ones vector _is_ orthogonal to the kernel? We can drop the kernel basis completely since it neither affects the objective or the constraint, and focus on the non-zero basis. We can then use the same analytic solution to $w$ but replace $\Sigma^{-1}$ with the pseudoinverse.

Lastly, I leave one exercise for the reader: if $\ones$ is an eigenvector of $\Sigma$, show that an optimal weight vector is $\frac{1}{k} \ones$.

## The bias-variance decomposition

Let's revisit our objective. What if we didn't care about the unbiased-ness constraint and only wanted to minimize expected squared error? For an invertible covariance matrix, does that lead to a different optimal weight vector?

Those familiar with statistical theory know what expected squared error of an estimate can always be decomposed into variance and squared bias.

$$\begin{align}
\Expect[w' X - \mu]^2
&= \Expect[w' X - w' \ones \mu + w' \ones \mu - \mu]^2 \\
&= \Expect[w' (X - \ones \mu)]^2 + (w' \ones - 1)^2 \mu^2 \\
&= w' \Sigma w + \mu^2 (w' \ones - 1)^2 \\
\end{align}
$$

To go further let's say we know the ratio $\Omega = \Sigma / \mu^2$ and scale the objective by $1 / \mu^2$. Then the optimal weights are

$$ w^* = (\Omega + \ones \ones')^{-1} \ones $$

Using my [favorite matrix identity](https://en.wikipedia.org/wiki/Woodbury_matrix_identity), this can be written another way

$$ w^* = \frac{\Omega^{-1} \ones}{1 + \ones' \Omega^{-1} \ones} $$

which is beautifully similar to the unbiased solution and results in the lower expected squared error of

$$ \frac{1}{1 + \ones' \Omega^{-1} \ones} $$


