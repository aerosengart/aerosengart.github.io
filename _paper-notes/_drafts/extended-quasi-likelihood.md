---
layout: distill
title: An Extended Quasi-Likelihood Function
description: 
date: 2026-03-18
tabs: true
tags: likelihood theory glmm
toc:
    - name: Background
bibliography: 2026-03-18-extended-quasi-likelihood.bib
---

In this post, I will work through the discussions in Nelder and Pregibon's <i>An Extended Quasi-Likelihood Function</i>.<d-cite key=nelder1987></d-cite> Their work proposes an extension to the quasi-likelihood function introduced by Wedderburn (1974)<d-cite key=wedderburn1974></d-cite> (see <a href="/paper-notes/quasi-likelihood">my post on that paper</a>).

---

## Background
For an observation, $y$, with mean $\mu$ and variance function $v(\mu)$, Wedderburn (1974)<d-cite key=wedderburn1974></d-cite> defined the <i>log quasi-likelihood</i> and <i>quasi-score</i>, respectively, as:

$$
\begin{aligned}
\ell_q(y; \mu) &= \int_{-\infty}^\mu \frac{y - u}{v(u)} du + h(y) \\
\frac{\partial \ell_q(y; \mu)}{\partial \mu_i} &= \frac{y - \mu}{v(\mu)}
\end{aligned}
$$

<aside><p>In Nelder and Pregibon, they use $Q(y; \mu)$ to denote the log quasi-likelihood.</p></aside>

for some function of $y$, $h(\cdot)$. 

We can also define a deviance function via a "log-likelihood ratio"-type statistic as:

$$
D(y; \mu) = - 2 \left(\ell_q(y; \mu) - \ell_q(y; y)\right) = -2 \int_{y}^\mu \frac{y - u}{v(u)} du
$$

Suppose we have $n$ independent observations, $y_i$, with associated $p$-dimensional covariate vectors, $\mathbf{x}_i$. Each observation has a mean, $\mu_i$, and a variance, $v(\mu_i) = \phi V(\mu_i)$. We assume that the mean is some <i>strictly monotone and differentiable</i> function of a linear predictor. That is:

$$
g(\mu_i) = \mathbf{x}_i^\top \beta
$$

for some vector of coefficients, $\beta \in \mathbb{R}^q$. Similar to true likelihood functions, the log quasi-likelihood of the dataset is the sum of the log quasi-likelihoods for each observation:

$$
\ell_q(\mathbf{y}; \mathbf{\mu}) = \sum_{i = 1}^n \ell_q(y_i; \mu_i)
$$

Since the maximum quasi-likelihood estimates are asymptotically normal, a deviance-based test statistic for nested hypotheses will be asymptotically chi-squared. However, this only holds if their variance functions are identical (since $v(\mu)$ appears in the denominator of the deviance expression). The goal of Nelder and Pregiborn's work is to define an <i>extended</i> quasi-likelihood function that permits this comparison. 

---

## An Extended Quasi-Likelihood
The <i>extended quasi-likelihood</i> for an observation, $y$, with mean $\mu$ and variance $v(\mu) = \phi V(\mu)$ is given by:

$$
Q^+(y; \mu) = -\frac{1}{2} \log(2 \pi v(\mu)) - \frac{1}{2} \frac{D(y; \mu)}{\phi}
$$

Since $\boldsymbol{\beta}$ only appears in $D(y; \mu)$, 






