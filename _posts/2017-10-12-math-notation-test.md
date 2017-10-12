---
layout: post
title: "Math Notation Test"
description: ""
category: 
tags: []
---

$$
  \newcommand{\R}{\mathbb{R}}
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

Let $\varphi_X(s) = \Expect[e^{i s^T X}]$ denote the characteristic
function of a random vector $X$. For some positive weight function
$w : \R^p \times \R^q \mapsto [0, \infty)$ define the norm
$\vert\cdot\vert_w : \{\gamma : \R^p \times \R^q \mapsto \mathbb{C}\} \mapsto [0, \infty)$ as
$\vert\gamma\vert_w^2 = \int_{\R^{p+q}} \vert\gamma(s,t)\vert^2 w(s,t)\diff s \diff t$.

When the subscript is omitted, the $\| \cdot \|$ is meant to denote
the Euclidean norm.

Let $X$ and $Y$ be two $p$ and $q$-dimensional (respectively) random
vectors with finite first moments and characteristic functions
$\varphi_X$ and $\varphi_Y$. Let the weight function
$w(s,t) = \dfrac{c_p c_q}{\|s\|^{1+p} \|t\|^{1+q}}$ where
$c_d = \dfrac{\Gamma(\frac{d+1}{2})}{\pi^{(d+1) / 2}}$
The *squared* distance covariance between $X$ and $Y$ is defined

$$
\begin{align*}
\mathcal{V}^2 (X,Y)
&= \| \varphi_{X,Y}(s,t) - \varphi_X(s)\varphi_Y(t) \|_w^2 \\
&= c_p c_q \int_{\R^{p+q}}
   \dfrac{|\varphi_{X,Y}(s,t) - \varphi_X(s)\varphi_Y(t)|^2}
         {\|s\|^{1+p} \|t\|^{1+q}} \diff s \diff t.
\end{align*}
$$

Distance covariance is defined
$\mathcal{V}(X,Y) = \sqrt{\mathcal{V}^2(X,Y)}$.
