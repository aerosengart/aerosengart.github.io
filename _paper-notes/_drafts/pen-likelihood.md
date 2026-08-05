---
layout: distill
title: Penalized Likelihood for General Semi-Parametric Regression Models
description:
date: 2026-03-16
tabs: true
tags: inference models likelihood
toc:
    - name: Background
    - name: Set-Up
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2

bibliography: 2026-03-16-pen-likelihood.bib
---

In anticipation of my second project write-up, I've decided to go through some of the literature on generalized linear mixed models and the like. In this post, I'll be working through an early paper by Peter J. Green titled <i>Penalized Likelihood for General Semi-parametric Regression Models</i>.<d-cite key=green1987></d-cite> This will provide some foundation for later study of <a href="/stats-ml/glmm/#penalized-quasi-likelihood">penalized quasi-likelihood</a> for generalized linear mixed models.

---

## Background
The main concern of Green's work is a regression model with a composite likelihood function. The model will be parametrized by an $m$-dimensional vector, $\theta$, which describes the systematic sources of variation in the responses, denoted by $\mathbf{y}$, through their dependence on some set(s) of explanatory covariates. $m$ is assumed to be so large that parametric inference is inappropriate.

More specifically, $\theta$ will be assumed to be known except for a $p$-dimensional vector, $\beta$, and some real-valued function $\gamma \in \mathcal{G}$ for some (infinite-dimensional) linear space, $\mathcal{G}$. Such a model will have a <i>composite likelihood function</i>, which we denote with:

$$
\mathcal{L}(\mathbf{y}; \theta(\beta, \gamma))
$$

<aside><p>Here, $\beta$ is the parametric part, and $\gamma$ is the non-parametric part, which gives us the name <i>semi-parametric regression</i>.</p></aside>

An example is:

$$
\theta_i = \mathbf{x}_i^\top \beta + \gamma(\mathbf{t}_i)
$$

which introduces dependence through an additive semi-parametric function.

Both $\beta$ and $\gamma$ are unknown model parameters, but the estimation of $\beta$ is our main goal. Unfortunately, without placing additional constraints on $\gamma$, the model parameters are unidentifiable. Green proposes estimated via penalized likelihood in order to overcome this difficulty.

---

## Set-Up
Because $\gamma$ is a non-parametric component, its values must be discretized at some point in the estimation procedure to make the problem tractable (since $\theta$ depends upon it). Green opts to restrict $\gamma$ to some finite-dimensional, linear subspace $\mathcal{F} \subset \mathcal{G}$ with the form:

$$
\mathcal{F} = \text{span}\left\{ \phi_j \rvert j = 1, 2, \dots, q \right\}
$$

for $q$ known basis functions that are independent of $\mathbf{y}$, but are permitted to depend on any covariates. 

