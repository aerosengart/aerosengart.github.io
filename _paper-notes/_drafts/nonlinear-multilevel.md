---
layout: distill
title: Nonlinear Multilevel Models, with an Application to Discrete Response Data
description:
date: 2026-03-19
tabs: true
tags: glmm inference models
toc:
    - name: Background
bibliography: 2026-03-16-nonlinear-multilevel.bib
---

In anticipation of my second project write-up, I've decided to go through some of the literature on generalized linear mixed models and the like. In this post, I'll be working through the results in <i>Nonlinear Multilevel Models, with an Application to Discrete Response Data</i> by Harvey Goldstein<d-cite key=goldstein1991></d-cite> as a preliminary foundation for <a href="/stats-ml/glmm/#marginal-quasi-likelihood">marginal quasi-likelihood</a> for generalized linear mixed models.

---

## Set-Up
We have $n$ observations of some response, $y_i$, each with an associated $p$-dimensional vector of fixed effects covariates, $\mathbf{x}_i$, and a $q$-dimensional vector of random effects covariates, $\mathbf{z}_i$. We collect these vectors as the rows in the design matrices, $\mathbf{X}$ and $\mathbf{Z}$, which are $n \times p$ and $n \times q$, respectively.

<aside><p>Quantities without a subscript will be arbitrary or the corresponding random variate.</p></aside>

A general multi-level model with $K > 1$ levels is given by:

$$
\begin{equation}
\label{eq:model}
\begin{aligned}
y &= \mathbf{x} \boldsymbol{\alpha} + \mathbf{z} \boldsymbol{\beta} \\
\mathbb{E}\left[ \boldsymbol{\beta} \right] &= \mathbf{0}_{(q \times k)} \\
\mathbb{E}\left[ (\mathbf{z} \boldsymbol{\beta})(\mathbf{z} \boldsymbol{\beta})^\top \right] &= \mathbf{V}
\end{aligned}
\end{equation}
$$

where $\boldsymbol{\alpha}$ is a $q$-dimensional vector of fixed effects, $\boldsymbol{\beta}$ is a $(q \times k)$-dimensional vector of the random effects (i.e. level-specific coefficients), and $\mathbf{V}$ is unknown but usually assumed to have the following structure:

$$
\begin{aligned}
\mathbf{V}_{k + 1} &= \mathbf{V}_k + \mathbf{z}^{k + 1} \left(\mathbb{I}_{n_{k + 1}} \otimes \Omega_{k + 1} \right) (\mathbf{z}^{k + 1})^\top
\end{aligned}
$$

The quantity $\Omega_{k+1}$ is the covariance matrix for the random effects that vary over the subjects in the $(k+1)$-th level, and $\mathbf{z}^{k+1}$ are the covariates associated with those effects. In the base case, $\mathbf{V}_1$ is the matrix of contributions to the total variation made at level $1$ of the model. 

<!-- #region example -->
<div class="example">
<strong>Example.</strong>
A relatively simple example that comes from Goldstein (1986)<d-cite key=goldstein1986></d-cite> is that of children in different classrooms at different schools. We let $y_{k,i,j}$ denote some response from the $j$-th child in the $i$-th classroom at the $k$-th school. We assume:

$$
y_{k,i,j} = \alpha^*_{k,i,j} + \beta^*_{k,i} + \gamma^*_k
$$

Each term in the righthand side of the above is assumed to be a linear combination of fixed and random effects:

$$
\begin{aligned}
\gamma^*_k &= \sum_{l = 0}^q \gamma_l w_{l,k} + v_k \\
\beta^*_{k,i} &= \sum_{l = 0}^p \beta_{l,k} z_{l,k,i} + u_{k,i} \\
\alpha^*_{k,i,j} &= \sum_{l = 0}^r \alpha_{l,k,i} x_{l,k,i,j} + e_{k,i,j}
\end{aligned}
$$

The $\gamma_l$, $\beta_{l,k}$, and $\alpha_{l,k,i}$ are the school-, classroom-, and child-level coefficients (respectively), and $w_{l,k}$, $z_{l,k,i}$, and $x_{l,k,i,j}$ are the associated covariates. The random components are assumed to satisfy:

$$
\begin{aligned}
\mathbb{E}[v_k] &= 0 & \text{Var}(v_k) &= \sigma_v^2 \\
\mathbb{E}[u_{k,i}] &= 0 & \text{Var}(u_{k,i}) &= \sigma_{u,k}^2 \\
\mathbb{E}[e_{k,i,j}] &= 0 & \text{Var}(e_{k,i,j}) \sigma^2_{e_{k,i}} \\
\end{aligned}
$$

Under the additional assumption that the random components are uncorrelated with each other, this set-up implies that:

$$
\text{Var}(y_{k,i,j}) = \sigma_v^2 + \sigma^2_{u_{k}} + \sigma^2_{e_{k,i}}
$$
</div>
<!-- #endregion -->

The multi-level model can be extended to non-linear functions, denoted by $f(\cdot)$, of the linear predictor via:

$$
\begin{equation}
\label{eq:nl-model}
\begin{aligned}
y &= f(\mathbf{x}_\alpha \boldsymbol{\alpha} + \mathbf{z}_\beta \boldsymbol{\beta}) + \mathbf{x}_a \boldsymbol{a} + \mathbf{z}_b \boldsymbol{b}
\end{aligned}
\end{equation}
$$

In this case, $$\boldsymbol{\beta}$$ and $$\boldsymbol{b}$$ are the $q$- and $q'$-dimensional random effects vectors (still with zero means), and $$\mathbf{z}_\beta$$ and $$\mathbf{z}_b$$ are their corresponding design matrices. Similarly, $$\boldsymbol{\alpha}$$ and $$\boldsymbol{a}$$ are the $p$- and $p'$-dimensional fixed effects vectors, and $$\mathbf{x}_\alpha$$ and $$\mathbf{x}_a$$ are their corresponding design matrices.

---

## Estimation

Let $\hat{\theta}_t = (\hat{\boldsymbol{\alpha}}_t, \hat{\boldsymbol{a}}_t, \hat{\boldsymbol{\beta}}_t, \hat{\boldsymbol{b}}_t)^\top$ be the vector of parameter estimates at iteration $t$. The non-linear piece of the model in Eq. \eqref{eq:nl-model} can be rewritten as:

<aside><p>By $\hat{\boldsymbol{\beta}}_t$ and $\hat{\boldsymbol{b}}_t$, we mean the corresponding variances and covariances.</p></aside>

$$
f(\mathbf{x}_\alpha \boldsymbol{\alpha} + \mathbf{z}_\beta \boldsymbol{\beta}) = f\left(\sum_{i = 1}^p x_{\alpha, i} \alpha_i + \sum_{j = 1}^q z_{\beta, j} \beta_j \right)
$$

At iteration $t+1$, we can approximate the non-linear component via a first-order Taylor approximation about the estimates at itertation $t$:

$$
\frac{\partial}{\partial z_{\beta, j'} \beta_{j'}} \left[ \sum_{j = 1}^q z_{\beta, j} \beta_j \right] = z_{\beta, j} \beta_j
$$

$$
f(\mathbf{x}_\alpha \boldsymbol{\alpha} + \mathbf{z}_\beta \boldsymbol{\beta})
= f\left(\sum_{i = 1}^p x_{\alpha, i} \alpha_{t,i} + \sum_{j = 1}^q z_{\beta, j} \beta_{t,j} \right) + \frac{\partial}{\partial } (\theta - \hat{\theta}_t)
$$





---
The model of concern for Goldstein is the two level loglinear model for a discrete response. The mean proportion for the $h$-th level $1$ unit in the $i$-th cell of the $j$-th level $2$ unit is given by:

$$
\begin{equation}
\label{eq:log-linear}
\begin{aligned}
\pi_{h,i,j} &= \exp\left(\sum_{k = 0}^l \alpha_{j,k} x_{h,i,j,k} \right) \\
\sum_{i = 1}^{m_j} \pi_{h,i,j} &= 1
\end{aligned}
\end{equation}
$$
---



