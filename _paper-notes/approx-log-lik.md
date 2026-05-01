---
layout: distill
title: Approximations to the Log-Likelihood Function in the Nonlinear Mixed-Effects Model
description:
date: 2026-03-21
tabs: true
tags: glmm inference models
toc:
    - name: Set-Up
    - name: LME Approximation
    - name: Laplace Approximation
      subsections:
        - name: Example
    - name: Importance Sampling Approximation
    - name: Adaptive Gaussian Quadrature Approximation
      subsections:
        - name: Example
bibliography: paper-notes.bib
---

In <i>Approximations to the Log-Likelihood Function in the Nonlinear Mixed Model</i>, Pinheiro and Bates discuss four different log-likelihood function approximations: the linear mixed effects approximation, the Laplace approximation, importance sampling, and adaptive Gaussian quadrature.<d-cite key=pinheiro1995></d-cite> Throughout the post, I will applying these methods to an example. Let's get the data set-up out of the way first. 

We assume that to have a dataset consisting of $N$ observations of: scalar response, $y$; a $p$-dimensional vector of covariates associated with fixed effects, $\mathbf{x}$; a scalar-valued covariate associated with random effects, $z$. We assume the observations are divided into $K$ independent clusters with $n_k$ observations per cluster such that $\sum_{k = 1}^K n_k = N$.

We assume that the responses follow a negative binomial distribution with conditonal means and variances given by:

$$
\begin{aligned}
\mathbb{E}[y_{k,i} \rvert \boldsymbol{\beta}] &= \mu_{k,i} = \exp(\eta_{k,i}) = \mathbf{x}_{k,i}^\top \boldsymbol{\alpha} + z_{k,i}\beta_k \\
\text{Var}(y_{k,i} \rvert \boldsymbol{\beta}) &= \mu_{k,i} + \frac{1}{\gamma} \mu_{k,i}^2
\end{aligned}
$$

Here, $\gamma > 0$ is an overdispersion parameter, $\boldsymbol{\alpha} \in \mathbb{R}^p$ is a vector of fixed effect coefficients, and $\beta_k \in \mathbb{R}$ is a vector of random effect coefficients satisfying:

$$
\begin{aligned}
\beta_k &\overset{iid}{\sim} \mathcal{N}(0, \tau^2); \hspace{4mm} k \in [K] 
\end{aligned}
$$

where $\tau^2$ is a variance component. We assume that observations within a cluster are conditionally independent.

Omitting the covariates from our notation, the conditional log-likelihood for cluster $k$, $\mathbf{y}_k$, is given by:

$$
\log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k))
= \sum_{i = 1}^{n_k} \left[ \log\left( \frac{\Gamma(y_{k,i} + \gamma + 2)}{y_{k,i}! \Gamma(\gamma + 2)}\right) + y_{k,i} \left[ \log(\mu_{k,i}) - \log(\mu_{k,i} + \gamma) \right] + \gamma\left[ \log(\gamma) - \log(\mu_{k,i} + \gamma)\right] \right]
$$

where $\Gamma(\cdot)$ denotes the Gamma function. We have placed a Gaussian prior on the random effects:

$$
\log(p(\boldsymbol{\beta}; \tau^2)) = -\frac{1}{2} \sum_{k = 1}^K \left[\log(2 \pi \tau^2) + \frac{\beta_k^2}{\tau^2} \right]
$$

Our goal will be to obtain maximum likelihood estimates for $\boldsymbol{\alpha}$, $\gamma$, and $\tau^2$. This requires maximizing:

$$
p(\mathbf{y}; \boldsymbol{\alpha}, \gamma, \tau^2)
= \int p(\mathbf{y}; \boldsymbol{\alpha}, \gamma \rvert \boldsymbol{\beta}) p(\boldsymbol{\beta}; \tau^2) d \boldsymbol{\beta}
$$


---

## Set-Up
The mixed effects model of Pinheiro and Bates has two levels. In the top-most level, we assume the data can be divided into $K$ groups, the $i$-th of which has $n_i$ observations. We also assume that the (conditional) mean of each observation is some non-linear function of a group-specific parameter vector, $\boldsymbol{\phi}$, and covariates, $\boldsymbol{\gamma}$. The $j$-th observation in the $i$-th group, denoted by $y_{i,j}$, is modeled as:

$$
\begin{equation}
\label{eq:model}
\begin{aligned}
y_{i,j} &= f(\boldsymbol{\phi}_{i,j}, \boldsymbol{\gamma}_{i,j}) + \epsilon_{i,j} \\
\boldsymbol{\phi}_{i,j} &= \boldsymbol{\phi}_{i,j}(\boldsymbol{\alpha}, \boldsymbol{\beta}_i) = \mathbf{X}_{i,j} \boldsymbol{\alpha} + \mathbf{Z}_{i,j} \boldsymbol{\beta}_i \\
\epsilon_{i,j} &\overset{iid}{\sim} \mathcal{N}(0, \sigma^2) \\
\boldsymbol{\beta}_i &\sim \mathcal{N}(\mathbf{0}_q, \sigma^2 \mathbf{D})
\end{aligned}
\end{equation}
$$

where $$\boldsymbol{\alpha}$$ is the $p$-dimensional vector of fixed effects, $$\boldsymbol{\beta}_i$$ is the $q$-dimensional vector of random effects for group $i$, an $\mathbf{X}$ and $\mathbf{Z}$ are the fixed effects and random effects design matrices, respectively, for the $j$-th observation in the $i$-th group. The model assumes that the observation errors, $$\epsilon_{i,j}$$ are independent of the random effects, and that the random effects have some general covariance matrix, $$\sigma^2 \mathbf{D}$$. 

Fitting the model in Eq. \eqref{eq:model} can be done in a variety of ways. Here, we are concerned with maximum likelihood estimation, which requires finding the values of the $\boldsymbol{\alpha}$, $\sigma^2$, and $\mathbf{D}$ which maximize the <i>marginal -likelihood</i> given by:

$$
\begin{equation}
\label{eq:marg-lik}
\mathcal{L}(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) = \int p(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D} \rvert \boldsymbol{\beta}) p(\boldsymbol{\beta}; \sigma^2, \mathbf{D}) d \boldsymbol{\beta}
\end{equation}
$$

<!-- #region marg-lik-group-i -->
{% tabs marg-lik-group-i %}
{% tab marg-lik-group-i statement %}
Due to the Gaussianity assumptions, the marginal likelihood for a single group, $i$, can be rewritten as:

$$
\begin{equation}
\label{eq:group-log-lik}
\begin{aligned}
\mathcal{L}(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \int \left[ (2 \pi \sigma^2)^{- \frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-\frac{1}{2\sigma^2}\left[ \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)\rvert \rvert_2^2 + \boldsymbol{\beta}_i \left[ \mathbf{D}\right]^{-1} \boldsymbol{\beta}_i \right] \right) \right] d \boldsymbol{\beta}_i 
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab marg-lik-group-i proof %}
$$
\begin{aligned}
\mathcal{L}(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \int 
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D} \rvert \boldsymbol{\beta}_i) p(\boldsymbol{\beta}_i; \sigma^2, \mathbf{D}) d \boldsymbol{\beta}_i \\
&= \int \left[ (2 \pi)^{- \frac{n_i}{2}}  \rvert \sigma^2 \mathbf{I}_{n_i \times n_i} \rvert^{-\frac{1}{2}} \exp\left(-\frac{1}{2}(\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i))^\top \left[ \sigma^2 \mathbb{I}_{n_i \times n_i} \right]^{-1} (\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)) \right)\right] \left[ (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2} \boldsymbol{\beta}_i \left[ \sigma^2 \mathbf{D}\right]^{-1} \boldsymbol{\beta}_i \right)\right] d \boldsymbol{\beta}_i \\
&= \int \left[ (2 \pi)^{- \frac{n_i + q}{2}} (\sigma^2)^{- \frac{n_i}{2}} (\sigma^2)^{-\frac{q}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-\frac{1}{2\sigma^2}\left[ (\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i))^\top(\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)) + \boldsymbol{\beta}_i \left[ \mathbf{D}\right]^{-1} \boldsymbol{\beta}_i \right] \right) \right] d \boldsymbol{\beta}_i \\
&= \int \left[ (2 \pi \sigma^2)^{- \frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-\frac{1}{2\sigma^2}\left[ \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)\rvert \rvert_2^2 + \boldsymbol{\beta}_i \left[ \mathbf{D}\right]^{-1} \boldsymbol{\beta}_i \right] \right) \right] d \boldsymbol{\beta}_i 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

These integrals are, in general, intractable due to the fact that they may be multi-dimensional and also involve the non-linear function, $f(\cdot)$. Thus, approximations are needed for model fitting and inference. 

---

## LME Approximation
The linear mixed effects (LME) method was developed by Lindstrom and Bates in their work <i>Newton-Raphson and EM Algorithms for Linear Mixed-Effects Models for Repeated Measures Data</i>.<d-cite key=lindstrom1988></d-cite> The log-likelihood approximation that Pinheiro and Bates refer to as the <i>LME approximation</i> is one step of the two step algorithm proposed by Lindstrom and Bates. We detail these steps below, which are alternated until convergence.

### Penalized Non-Linear Least Squares Step
Denote the estimate of $\mathbf{D}$ at iteration $t$ with $\hat{\mathbf{D}}^{(t)}$. We obtain the conditional mode of $\boldsymbol{\beta}$ and the conditional estimate of $\boldsymbol{\alpha}$ as:

$$
\begin{equation}
\label{eq:pnls-optim}
(\hat{\boldsymbol{\alpha}}^{(t)}, \hat{\boldsymbol{\beta}}^{(t)}) = \underset{\boldsymbol{\alpha}, \boldsymbol{\beta}}{\arg \min} \left\{ \sum_{i = 1}^K \left(\rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \rvert \rvert_2^2 + \boldsymbol{\beta}_i^\top \left[\hat{\mathbf{D}}^{(t)}\right]^{-1} \boldsymbol{\beta}_i \right) \right\}
\end{equation}
$$

where $$\mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma})$$ is the vector of $$f(\boldsymbol{\phi}_{i,j}, \boldsymbol{\gamma}_{i,j})$$ for group $i$. 

### Linear Mixed Effects Step
In this step, we find $\hat{\mathbf{D}}^{(t+1)}$ with an approximation of the log-likelihood via a first-order Taylor expansion. The point of expansion is the current estimate/mode, 

Recall that $\boldsymbol{\phi} = \boldsymbol{\phi}(\boldsymbol{\alpha}, \boldsymbol{\beta})$. A first-order Taylor expansion of $\mathbf{f}(\cdot)$ about $\hat{\boldsymbol{\theta}}_i^{(t)} = (\hat{\boldsymbol{\alpha}}^{(t)}, \hat{\boldsymbol{\beta}}^{(t)}_i)$ is given by:

$$
\begin{aligned}
\hat{\mathbf{f}}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)
&= 
\mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \\
&\approx 
\left. \left[ \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}} + (\boldsymbol{\theta}_i - \hat{\boldsymbol{\theta}}_i^{(t)})^\top\left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\theta}_i} \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}} \\
&= \left. \left[ \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}} + (\boldsymbol{\alpha} - \hat{\boldsymbol{\alpha}}^{(t)})^\top\left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\alpha}} \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}} + (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i^{(t)})^\top\left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}}
\end{aligned}
$$

Define the following quantities:

$$
\begin{aligned}
\tilde{\mathbf{f}}^{\boldsymbol{\alpha}}_i(\hat{\boldsymbol{\theta}}_i^{(t)}) &= \left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\alpha}} \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}} \\
\tilde{\mathbf{f}}^{\boldsymbol{\beta}_i}_i(\hat{\boldsymbol{\theta}}_i^{(t)}) &= \left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}} \\
\hat{\mathbf{y}}_i^{(t)} &= \mathbf{y}_i - \left. \left[ \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \right] \right\rvert_{\boldsymbol{\theta}_i = \hat{\boldsymbol{\theta}}_i^{(t)}}  - (\hat{\boldsymbol{\alpha}}^{(t)})^\top \tilde{\mathbf{f}}_i^{\alpha}(\hat{\boldsymbol{\theta}}_i^{(t)}) - (\hat{\boldsymbol{\beta}}_i^{(t)})^\top \tilde{\mathbf{f}}_i^{\alpha}(\hat{\boldsymbol{\theta}}_i^{(t)})
\end{aligned}
$$ 


{% tabs lme-approx %}
{% tab lme-approx statement %}
The approximate log-likelihood is then:

$$
\ell_{LME}(\boldsymbol{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
= - \frac{1}{2} \sum_{i = 1}^K \left[ \log \left( \left\rvert \sigma^2 \left(\mathbb{I}_{n_i \times n_i} + (\tilde{\mathbf{f}}_i^{\boldsymbol{\beta}_i}(\hat{\boldsymbol{\theta}}_i^{(t)}))\mathbf{D} \tilde{\mathbf{f}}_i^{\boldsymbol{\beta}_i}(\hat{\boldsymbol{\theta}}_i^{(t)})^\top  \right) \right\rvert \right) + \frac{1}{\sigma^2} \left(\hat{\mathbf{y}}_i^{(t)} - \tilde{\mathbf{f}}_i^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}}_i^{(t)}) \boldsymbol{\alpha} \right)^\top \left[ \mathbb{I}_{n_i \times n_i} + (\tilde{\mathbf{f}}_i^{\boldsymbol{\beta}_i}(\hat{\boldsymbol{\theta}}_i^{(t)}))^\top \mathbf{D} \tilde{\mathbf{f}}_i^{\boldsymbol{\beta}_i}(\hat{\boldsymbol{\theta}}_i^{(t)}) \right]^{-1}  \left(\hat{\mathbf{y}}_i^{(t)} - \tilde{\mathbf{f}}_i^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}}_i^{(t)}) \boldsymbol{\alpha} \right) \right]
$$
{% endtab %}
{% tab lme-approx proof %}
TO DO
{% endtab %}
{% endtabs %}

Notice that (for a fixed iteration $t$) this approximation is the same as the log-likelihood function for a linear mixed effects model where the response is $\hat{\mathbf{y}}^{(t)} = (\hat{\mathbf{y}}^{(t)}_1, \dots, \hat{\mathbf{y}}_K^{(t)})^\top$, and the fixed and random effects design matrices are given by:

$$
\begin{aligned}
\tilde{\mathbf{f}}^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}}^{(t)}) 
&= 
\begin{bmatrix}
\rvert &  & \rvert \\
\tilde{\mathbf{f}}_1^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}_1}^{(t)}) 
& \dots & \tilde{\mathbf{f}}_K^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}_K}^{(t)}) \\
\rvert & & \rvert
\end{bmatrix}
\\
\tilde{\mathbf{f}}^{\boldsymbol{\beta}}(\hat{\boldsymbol{\theta}}^{(t)})
&= 
\begin{bmatrix}
\rvert &  & \rvert \\
\tilde{\mathbf{f}}_1^{\boldsymbol{\beta}_1}(\hat{\boldsymbol{\theta}_1}^{(t)}) 
& \dots & \tilde{\mathbf{f}}_K^{\boldsymbol{\beta}_K}(\hat{\boldsymbol{\theta}_K}^{(t)}) \\
\rvert & & \rvert
\end{bmatrix}
\end{aligned}
$$

One can also perform the estimation using a <i>restricted</i> log-likelihood approximation as in Lindstrom and Bates (1990)<d-cite key=lindstrom1988></d-cite>:

$$
\ell_{LME}^R(\boldsymbol{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
= - \frac{1}{2} \sum_{i = 1}^K  \log \left( \left\rvert \sigma^2 (\tilde{\mathbf{f}}_i^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}}_i^{(t)}))^\top  \left[\mathbb{I}_{n_i \times n_i} + (\tilde{\mathbf{f}}_i^{\boldsymbol{\beta}_i}(\hat{\boldsymbol{\theta}}_i^{(t)}))\mathbf{D} (\tilde{\mathbf{f}}_i^{\boldsymbol{\beta}_i}(\hat{\boldsymbol{\theta}}_i^{(t)}))^\top  \right]^{-1} \right\rvert \tilde{\mathbf{f}}_i^{\boldsymbol{\alpha}}(\hat{\boldsymbol{\theta}}_i^{(t)}) \right) + \ell_{LME}(\boldsymbol{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
$$

This is equivalent to a Laplacian approximation of the integral in Eq. \eqref{eq:marg-lik} but with a flat prior on the fixed effects.

---

## Laplace Approximation
A Laplace approximation performs a second-order Taylor expansion of the conditional log-likelihood about a maximum a posteriori (MAP) estimate of the random effects with respect to some prior. Since the point of expansion is chosen to be the mode, the first-order term in the expansion disappears, and the integrand resembles a Gaussian kernel.

Define:

$$
\begin{aligned}
g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i)
&= \frac{1}{2 \sigma^2} \left[ \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \rvert \rvert_2^2 + \boldsymbol{\beta}_i^\top \mathbf{D}^{-1} \boldsymbol{\beta}_i \right] \\
\hat{\boldsymbol{\beta}}_i &= \underset{\boldsymbol{\beta}_i}{\arg \min} \left\{ g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i) \right\} \\
g(\hat{\boldsymbol{\beta}}_i) &= \left. \left[  g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i) \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} \\
g'(\hat{\boldsymbol{\beta}_i}) &= \left. \left[ \frac{\partial g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i)}{\partial \boldsymbol{\beta}_i} \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} \\
g''(\hat{\boldsymbol{\beta}_i}) &= \left. \left[ \frac{\partial^2 g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i)}{\partial \boldsymbol{\beta}_i \partial \boldsymbol{\beta}_i^\top} \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i}
\end{aligned}
$$

<aside><p>Note: I have an a factor of $\frac{1}{2\sigma^2}$ in $g(\cdot)$ compared to the publication.</p></aside>

<!-- #region laplace-nlmem -->
{% tabs laplace-nlmem %}
{% tab laplace-nlmem statement %}
$$
\begin{aligned}
\ell_{LAP}(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&\approx - \frac{1}{2} \sum_{i = 1}^K \left[ (n_i + q) \log(\sigma^2)  + n_i \log(2 \pi) + \log\left( \rvert g''(\hat{\boldsymbol{\beta}}_i) \rvert \right) + \log(\rvert \mathbf{D} \rvert) + 2 g(\hat{\boldsymbol{\beta}}_i) \right]
\end{aligned}
$$
{% endtab %}
{% tab laplace-nlmem proof %}
A second-order Taylor expansion of $g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i)$ about $\hat{\boldsymbol{\beta}}_i$ is given by:

$$
\begin{aligned}
g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \boldsymbol{\beta}_i)
&\approx g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \hat{\boldsymbol{\beta}}_i) + (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i)^\top \underbrace{g'(\hat{\boldsymbol{\beta}}_i)}_{= \mathbf{0}_q} + \frac{1}{2} (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i)^\top g''(\hat{\boldsymbol{\beta}}_i) (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i) \\
&=  g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \hat{\boldsymbol{\beta}}_i) + \frac{1}{2} (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i)^\top g''(\hat{\boldsymbol{\beta}}_i) (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i) 
\end{aligned}
$$

Using the second-order Taylor approximation of $g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D} \rvert \boldsymbol{\beta}_i)$, we can write:

$$
\begin{aligned}
p(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \prod_{i = 1}^K \int (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-g(\mathbf{y}_i;, \boldsymbol{\beta}_i, \boldsymbol{\alpha}, \mathbf{D})\right) d \boldsymbol{\beta}_i \\
&\approx \prod_{i = 1}^K \int (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-\left[ g(\hat{\boldsymbol{\beta}}_i) + \frac{1}{2} (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i)^\top g''(\hat{\boldsymbol{\beta}}_i) (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i) \right] \right) d \boldsymbol{\beta}_i \\
&= \prod_{i = 1}^K (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-g(\hat{\boldsymbol{\beta}}_i) \right) \int \frac{(2 \pi)^{\frac{q}{2}} \left\rvert \left[ g''(\hat{\boldsymbol{\beta}}_i) \right]^{-1} \right\rvert^{\frac{1}{2}}}{(2 \pi)^{\frac{q}{2}} \left\rvert \left[ g''(\hat{\boldsymbol{\beta}}_i)\right]^{-1} \right\rvert^{\frac{1}{2}}}\exp\left( -\frac{1}{2} (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i)^\top g''(\hat{\boldsymbol{\beta}}_i) (\boldsymbol{\beta}_i - \hat{\boldsymbol{\beta}}_i) \right) d \boldsymbol{\beta}_i \\
&= \prod_{i = 1}^K  (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} (2 \pi)^{\frac{q}{2}} \left\rvert g''(\hat{\boldsymbol{\beta}}_i) \right\rvert^{-\frac{1}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-g(\hat{\boldsymbol{\beta}}_i) \right)
\end{aligned}
$$

Notice how the third line in the above illustrates how the integrand is approximately proportional to a Gaussian density with mean and covariance (this will be helpful later!):

$$
\begin{equation}
\label{eq:approx-gauss-density}
\begin{aligned}
\boldsymbol{\mu}_i &= \hat{\boldsymbol{\beta}}_i \\
\boldsymbol{\Sigma}_i &= \left[ g''(\hat{\boldsymbol{\beta}}_i) \right]^{-1}
\end{aligned}
\end{equation}
$$

Pinheiro and Bates do one more approximation: the Hessian of $g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D} \rvert \boldsymbol{\beta}_i)$. By the chain rule, we have:

$$
\begin{aligned}
g''(\hat{\boldsymbol{\beta}}_i)
&= \frac{1}{2\sigma^2}  \left. \left[ \left(\frac{\partial}{\partial \boldsymbol{\beta}_i} \left[ -2 \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top} (\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)) + \mathbf{D}^{-1} \boldsymbol{\beta}_i \right] \right) \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} \\
&= \frac{1}{2 \sigma^2} \left. \left[ \left(-2 \left[ \frac{\partial^2 \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i \partial \boldsymbol{\beta}_i^\top} (\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)) - \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top} \right] + \mathbf{D}^{-1} \right) \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} \\
&= \frac{1}{\sigma^2}  \left. \left[ -\frac{\partial^2 \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i \partial \boldsymbol{\beta}_i^\top} (\mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)) + \frac{1}{2} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top} + \frac{1}{2} \mathbf{D}^{-1} \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} \\
&\approx 
\frac{1}{2\sigma^2} \left( \left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top}\right]\right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} + \mathbf{D}^{-1} \right)
\end{aligned}
$$

The last line follows by an argument similar to that used in Gauss-Newton optimization. When evaluated at $\boldsymbol{\hat{\beta}}_i$, the first term should be very, very close to zero (since the fit should be nearly exact).

Plugging this in and taking the natural logarithm yields the Laplace approximation in the original publication:

$$
\begin{aligned}
\ell_{LAP}(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \log\left[ \prod_{i = 1}^K (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} (2 \pi)^{\frac{q}{2}} \left\rvert  g''(\hat{\boldsymbol{\beta}}_i) \right\rvert^{-\frac{1}{2}} \rvert \mathbf{D} \rvert^{-\frac{1}{2}} \exp\left(-g(\hat{\boldsymbol{\beta}}_i) \right)  \right] \\
&= \sum_{i = 1}^K \left[ -\frac{n_i + q}{2} \log(2 \pi \sigma^2) + \frac{q}{2} \log(2 \pi) - \frac{1}{2} \log\left( \rvert g''(\hat{\boldsymbol{\beta}}_i) \rvert \right) - \frac{1}{2} \log(\rvert \mathbf{D} \rvert) - g(\hat{\boldsymbol{\beta}}_i)  \right] \\
&= - \frac{1}{2} \sum_{i = 1}^K \left[ (n_i + q) \log(\sigma^2)  + n_i \log(2 \pi) + \log\left( \rvert g''(\hat{\boldsymbol{\beta}}_i) \rvert \right) + \log(\rvert \mathbf{D} \rvert) + 2 g(\hat{\boldsymbol{\beta}}_i) \right]
\end{aligned}
$$
{% endtab %}
{% endtabs %}

<!-- #endregion -->

<!-- We can pretty easily derive the MLE of $\sigma^2$ based upon this Laplace approximation of the marginal log-likelihood:

$$
\begin{aligned}
\frac{\partial \ell_{LAP}}{\partial \sigma^2}
&= - \frac{1}{2 \sigma^2} \left(\sum_{i = 1}^K n_i \right) - \frac{1}{2} \sum_{i = 1}^K \left[ \frac{1}{\sigma^2} - \frac{1}{(\sigma^2)^2} \left(  \left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} + \mathbf{D}^{-1} \right) \right] \\
\overset{\ell_{LAP} = 0}{\implies}
\frac{1}{2 \sigma^2} \sum_{i = 1}^K \left[ n_i + 1 \right] &= \frac{1}{2 (\sigma^2)^2} \sum_{i = 1}^K \left(  \left.\left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} + \mathbf{D}^{-1} \right) \\
\implies \sigma^2 &= \frac{1}{K + \sum_{i = 1}^K n_i} \sum_{i = 1}^K \left( \left. \left[ \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i^\top} \frac{\partial \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i)}{\partial \boldsymbol{\beta}_i} \right] \right\rvert_{\boldsymbol{\beta}_i = \hat{\boldsymbol{\beta}}_i} + \mathbf{D}^{-1} \right)
\end{aligned}
$$ -->

### Example
For each cluster, we define:

$$
q(\beta_k) = \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2))
$$

Suppose we knew the value of $\beta_k$ that maximized $q(\beta_k)$ (denote this with $\hat{\beta}_k$). We could then construct a second-order Taylor approximation of $q(\beta_k)$ as:

$$
q(\beta_k) \approx q(\hat{\beta}_k) + (\beta_k - \hat{\beta}_k) q'(\hat{\beta}_k) + \frac{1}{2}(\beta_k - \hat{\beta}_k)^2 q''(\hat{\beta}_k)
$$

<!-- #region q-derivs -->
{% tabs q-derivs %}
{% tab q-derivs statement %}
$$
\begin{aligned}
\frac{d q(\beta_k)}{d \beta_k} &= - \frac{\beta_k}{\tau^2} + \sum_{i = 1}^{n_k} z_{k,i} \left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \\
\frac{d^2 q(\beta_k)}{d \beta_k^2} &= - \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \right]
\end{aligned}
$$
{% endtab %}
{% tab q-derivs proof %}
First, we have the following intermediate derivatives:

$$
\begin{aligned}
\frac{d}{d \beta_k} \left[ \log(\mu_{k,j}) \right] 
&= \frac{d}{d \beta_k} \left[ \mathbf{x}_{k,j}^\top \boldsymbol{\alpha} + z_{k,j} \beta_k  \right]
= z_{k,j} \\
\frac{d}{d \beta_k} \left[ \mu_{k,j} \right]
&= \exp(\eta_{k,j}) \frac{d}{d \beta_k} \left[ \mathbf{x}_{k,j}^\top \boldsymbol{\alpha} + z_{k,j} \beta_k \right]
= z_{k,j} \mu_{k,j} \\
\frac{d^2}{d \beta_k^2} \left[ \mu_{k,j} \right] 
&=\frac{d}{d \beta_k} \left[ z_{k,j} \mu_{k,j} \right] 
= z_{k,j}^2 \mu_{k,j}
\end{aligned}
$$

The first derivative of $q(\beta_k)$ with respect to $\beta_k$ is:

$$
\begin{aligned}
\frac{d}{d \beta_k} \left[ \log(p(\beta_k; \tau^2)) \right]
&= \frac{d}{d\beta_k} \left[ - \frac{1}{2} \log(2 \pi \tau^2) - \frac{\beta_k^2}{2 \tau^2}\right] 
= - \frac{\beta_k}{\tau^2} \\
\frac{d}{d\beta_k} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) \right] 
&= \sum_{i = 1}^{n_k} \left[  \frac{d}{d \beta_k} \left[ \log\left( \frac{\Gamma(y_{k,i} + \gamma + 2)}{y_{k,i}! \Gamma(\gamma + 2)}\right) + y_{k,i} \left[ \log(\mu_{k,i}) - \log(\mu_{k,i} + \gamma) \right] + \gamma\left[ \log(\gamma) - \log(\mu_{k,i} + \gamma)\right] \right] \right] \\
&= \sum_{i = 1}^{n_k} \left[ y_{k,i} \left( \frac{d}{d \beta_k} \left[ \log(\mu_{k,i}) \right] - \frac{d}{d \beta_k}\left[ \log(\mu_{k,i} + \gamma) \right] \right) - \gamma \frac{d}{d \beta_k} \left[ \log(\mu_{k,i} + \gamma) \right] \right]  \\
&= \sum_{i = 1}^{n_k} \left[ y_{k,i} \left(z_{k,i} - \frac{1}{\mu_{k,i} + \gamma} (z_{k,i} \mu_{k,i}) \right) - \frac{\gamma}{\mu_{k,i} + \gamma} (z_{k,i} \mu_{k,i}) \right] \\
&= \sum_{i = 1}^{n_k} z_{k,i} \left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \\
\implies \frac{d}{d \beta_k} \left[ q(\beta_k) \right]
&= \frac{d}{d \beta_k} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2)) \right] \\
&= - \frac{\beta_k}{\tau^2} + \sum_{i = 1}^{n_k} z_{k,i} \left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right]
\end{aligned}
$$

And the second derivative is:

$$
\begin{aligned}
\frac{d^2}{d \beta_k^2} \left[ \log(p(\beta_k; \tau^2)) \right] 
&= \frac{d}{\beta_k} \left[ - \frac{\beta_k}{\tau^2} \right]
= - \frac{1}{\tau^2} \\
\frac{d^2}{d\beta_k^2} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) \right] 
&= \sum_{i = 1}^{n_k} z_{k,i} \frac{d}{d \beta_k}\left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \\
&= \sum_{i = 1}^{n_k} z_{k,i} \left[ - \frac{(\mu_{k,i} + \gamma) \frac{d}{d \beta_k} \left[ \mu_{k,i} \right] - \mu_{k,i} \frac{d}{d \beta_k} \left[ \mu_{k,i} + \gamma\right]}{(\mu_{k,i} + \gamma)^2} - \frac{(\mu_{k,i} + \gamma) \frac{d}{d \beta_k} \left[ \gamma \mu_{k,i} \right] - \gamma \mu_{k,i} \frac{d}{d \beta_k} \left[ \mu_{k,i} + \gamma \right]}{(\mu_{k,i} + \gamma)^2} \right] \\
&= - \sum_{i = 1}^{n_k} z_{k,i} \left[ \frac{(\mu_{k,i} + \gamma)(z_{k,i} \mu_{k,i}) - \mu_{k,i}(z_{k,i} \mu_{k,i})}{(\mu_{k,i} + \gamma)^2} + \frac{(\mu_{k,i} + \gamma)(\gamma z_{k,i} \mu_{k,i}) - \gamma \mu_{k,i}(z_{k,i} \mu_{k,i})}{(\mu_{k,i} + \gamma)^2} \right] \\
&= - \sum_{i = 1}^{n_k} z_{k,i} \left[ \frac{z_{k,i} \mu_{k,i}^2 + \gamma z_{k,i} \mu_{k,i} - z_{k,i} \mu_{k,i}^2 + \gamma z_{k,i} \mu_{k,i}^2 + \gamma^2 z_{k,i} \mu_{k,i} - \gamma z_{k,i} \mu_{k,i}^2}{(\mu_{k,i} + \gamma)^2} \right] \\
&= - \sum_{i = 1}^{n_k} z_{k,i} \left[ \frac{\gamma z_{k,i} \mu_{k,i} + \gamma^2 z_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2}\right] \\
&= -  \gamma (\gamma + 1)\sum_{i = 1}^{n_k} \frac{z_{k,i}^2 \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \\
\implies \frac{d^2}{d \beta_k^2} \left[ q(\beta_k) \right] 
&= \frac{d^2}{d \beta_k^2} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2)) \right] \\
&= - \frac{1}{\tau^2} - \gamma (\gamma + 1)\sum_{i = 1}^{n_k} \frac{z_{k,i}^2 \mu_{k,i}}{(\mu_{k,i} + \gamma)^2}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

Because $\hat{\beta}_k$ maximizes $q(\beta_k)$, it should also satisfy $q'(\beta_k) = 0$. Thus, the approximation can be simplified to:

$$
q(\beta_k) \approx  q(\hat{\beta}_k) + \frac{1}{2}(\beta_k - \hat{\beta}_k)^2 q''(\hat{\beta}_k)
$$

If we have the posterior mode for $\boldsymbol{\beta}$, then we can just plug this value into the $q(\beta_k)$ functions and use some arithmetric tricks to get:

<!-- #region ell-laplace-ex -->
{% tabs ell-laplace-ex %}
{% tab ell-laplace-ex statement %}
$$
\ell_{LAP}(\boldsymbol{\alpha}, \gamma, \tau^2; \mathbf{y}_k)  = \sqrt{2 \pi} \exp\left(  \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \hat{\beta}_k)) + \log(p(\hat{\beta}_k; \tau^2)) \right) \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \right]^{-\frac{1}{2}}
$$
{% endtab %}
{% tab ell-laplace-ex proof %}
$$
\begin{aligned}
\ell(\boldsymbol{\alpha}, \gamma, \tau^2; \mathbf{y}_k) 
&=\int p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k) p(\beta_k; \tau^2) d \beta_k \\
&= \int \exp(q(\beta_k)) d\beta_k \\
&\approx \int \exp\left(  q(\hat{\beta}_k) + \frac{1}{2}(\beta_k - \hat{\beta}_k)^2 q''(\hat{\beta}_k) \right) d\beta_k \\
&= \exp(q(\hat{\beta}_k)) \int \frac{(2 \pi (-q''(\hat{\beta}_k))^{-1})^{-\frac{1}{2}}}{(2 \pi (-q''(\hat{\beta}_k))^{-1})^{-\frac{1}{2}}}\exp\left(- \frac{1}{2(-q''(\hat{\beta}_k))^{-1}}(\beta_k - \hat{\beta}_k)^2 \right) d \beta_k \\
&=  \exp(q(\hat{\beta}_k)) (2 \pi (-q''(\hat{\beta}_k))^{-1})^{\frac{1}{2}} \underbrace{\int (2 \pi (-q''(\hat{\beta}_k))^{-1})^{-\frac{1}{2}}\exp\left(- \frac{1}{2(-q''(\hat{\beta}_k))^{-1}}(\beta_k - \hat{\beta}_k)^2 \right) d \beta_k}_{=1} \\
&= \sqrt{2 \pi} \exp(q(\hat{\beta}_k)) (-q''(\hat{\beta}_k))^{-\frac{1}{2}} \\
&= \sqrt{2 \pi} \exp\left(  \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \hat{\beta}_k)) + \log(p(\hat{\beta}_k; \tau^2)) \right) \left(- \left( - \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \right] \right) \right)^{-\frac{1}{2}} \\
&= \sqrt{2 \pi} \exp\left(  \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \hat{\beta}_k)) + \log(p(\hat{\beta}_k; \tau^2)) \right) \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \right]^{-\frac{1}{2}}
\end{aligned} 
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->





---

## Importance Sampling Approximation
An alternative to the previous approximations, which were more analytical in flavor, is to take a sampling-based approach. Once can view the integral in Eq. \eqref{eq:marg-lik} as the expectation of the conditional log-likelihood taken over the distribution of the random effects. In importance sampling, we take repeated samples from a specially chosen distribution in order to approximate this expectation as the average of our samples.

Let $$\mathbf{z}^* \sim q(\cdot)$$ be a latent variable that we introduce. We pick $$\mathbf{z}^*$$ such that $$\boldsymbol{\beta}_{i} = \boldsymbol{\beta}_i(\mathbf{z}^*)$$ (some function, $\boldsymbol{\beta}_i(\cdot)$, of $$\mathbf{z}^*$$) and so that it has a simpler distribution than the random effects prior distribution. Importance sampling makes the following approximation:

$$
\begin{aligned}
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) 
&=  \int p(\mathbf{y}_i, \boldsymbol{\beta}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) d \boldsymbol{\beta}_i \\
&=  \int \frac{p(\mathbf{y}_i, \boldsymbol{\beta}_i(\mathbf{z}^*); \boldsymbol{\alpha}, \sigma^2, \mathbf{D})}{q(\mathbf{z}^*)} q(\mathbf{z}^*) d \mathbf{z}^* \\
&= \mathbb{E}_{\mathbf{z}^*} \left[ \frac{p(\mathbf{y}_i, \boldsymbol{\beta}_i(\mathbf{z}^*); \boldsymbol{\alpha}, \sigma^2, \mathbf{D})}{q(\mathbf{z}^*)} \right] \\
&\approx \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \frac{p(\mathbf{y}_i, \boldsymbol{\beta}_i(\mathbf{z}^*_j); \boldsymbol{\alpha}, \sigma^2, \mathbf{D})}{q(\mathbf{z}^*_j)}
\end{aligned}
$$

where $n^{IS}$ is the number of importance samples taken. When working out the Laplace approximation (see Eq. \eqref{eq:approx-gauss-density}), we found that the integrand for cluster $i$ (i.e. the joint distribution barring constants) is approximately proportional to a Gaussian density with mean and covariance:

$$
\begin{aligned}
\boldsymbol{\mu}_i &= \hat{\boldsymbol{\beta}}_i \\
\boldsymbol{\Sigma}_i &= \left[ g''(\hat{\boldsymbol{\beta}}_i) \right]^{-1}
\end{aligned}
$$

This implies that the posterior distribution should also be proportional to this density. Since this is the optimal proposal distribution (i.e. $q(\cdot)$), we will draw samples as:

$$
\begin{aligned}
\mathbf{z}^* &\sim \mathcal{N}(\mathbf{0}_q, \mathbb{I}_{q \times q}) \\
\boldsymbol{\beta}_i(\mathbf{z}^*) &= \boldsymbol{\mu}_i + \boldsymbol{\Lambda}^{-\frac{1}{2}}\left[\boldsymbol{\Sigma}_i\right] \mathbf{z}^*
\end{aligned}
$$

where $\Lambda^{-\frac{1}{2}}[\cdot]$ denotes the inverse Cholesky factor. Let a quantity with $\mathbf{z}^*$ as its argument denote the quantity evaluated with respect to the importance samples. For example:

$$
\boldsymbol{\phi}_i(\mathbf{z}^*) = \mathbf{X}_{i,j} \boldsymbol{\alpha} + \mathbf{Z}_{i,j} \boldsymbol{\beta}_i(\mathbf{z}^*)
$$

<!-- #region is-nlmem -->
{% tabs is-nlmem %}
{% tab is-nlmem statement %}
$$
\ell_{IS}(\boldsymbol{\alpha}, \sigma^2, \mathbf{D}; \mathbf{y}) = \sum_{i = 1}^K \left[ - \frac{n_i}{2} \log(2 \pi \sigma^2) - \frac{1}{2} \log(\rvert \sigma^2 \mathbf{D} \rvert) + \frac{1}{2} \log(\rvert \boldsymbol{\Sigma}_i \rvert) + \log\left( \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \exp\left( - g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2} \rvert \rvert \mathbf{z}^*_j\rvert \rvert_2^2 \right) \right) \right] 
$$
{% endtab %}
{% tab is-nlmem proof %}
For a single cluster, we obtain:

$$
\begin{aligned}
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&\approx \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \frac{p(\mathbf{y}_i, \boldsymbol{\beta}_i(\mathbf{z}^*_j); \boldsymbol{\alpha}, \sigma^2, \mathbf{D})}{q(\mathbf{z}^*_j)} \\
&= \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \frac{p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D} \rvert \boldsymbol{\beta}_i(\mathbf{z}^*_j)) p(\boldsymbol{\beta}_i(\mathbf{z}^*_j); \sigma^2, \mathbf{D})}{q(\boldsymbol{\beta}_i(\mathbf{z}^*_j))} \\
&= \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \underbrace{\left[ (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \left\rvert \sigma^2 \mathbf{D} \right\rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2 \sigma^2} \left[ \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i(\mathbf{z}^*_j); \boldsymbol{\gamma}_i) \rvert \rvert_2^2 + \boldsymbol{\beta}_i^\top(\mathbf{z}^*_j)\mathbf{D} ^{-1} \boldsymbol{\beta}^\top_i(\mathbf{z}^*_j) \right] \right) \right]}_{= p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D} \rvert \boldsymbol{\beta}_i(\mathbf{z}^*_j)) p(\boldsymbol{\beta}_i(\mathbf{z}^*_j); \sigma^2, \mathbf{D})} \underbrace{\left[ (2 \pi)^{\frac{q}{2}} \left\rvert \boldsymbol{\Sigma}_i \right\rvert^{\frac{1}{2}} \exp\left(\frac{1}{2}(\boldsymbol{\beta}_i(\mathbf{z}^*_j) - \boldsymbol{\mu}_i)^\top \boldsymbol{\Sigma}_i^{-1}(\boldsymbol{\beta}_i(\mathbf{z}^*_j)- \boldsymbol{\mu}_i \right) \right]}_{= \left[ q(\boldsymbol{\beta}_{i}(\mathbf{z}^*_j)) \right]^{-1}} \\
&= \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \left[ (2 \pi \sigma^2)^{-\frac{n_i}{2}} \left\rvert \sigma^2 \mathbf{D} \right\rvert^{-\frac{1}{2}} \left\rvert \boldsymbol{\Sigma}_i \right\rvert^{\frac{1}{2}} \exp\left(- g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2}(\boldsymbol{\mu}_i + \boldsymbol{\Lambda}^{- \frac{1}{2}}\left[ \boldsymbol{\Sigma}_i \right] \mathbf{z}^*_j - \boldsymbol{\mu}_i)^\top \boldsymbol{\Sigma}_i^{-1} (\boldsymbol{\mu}_i + \boldsymbol{\Lambda}^{- \frac{1}{2}}\left[ \boldsymbol{\Sigma}_i \right] \mathbf{z}^*_j - \boldsymbol{\mu}_i) \right) \right] \\
&= \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \left[ (2 \pi \sigma^2)^{- \frac{n_i}{2}} \left\rvert \sigma^2 \mathbf{D} \right\rvert^{-\frac{1}{2}} \left\rvert \boldsymbol{\Sigma}_i \right\rvert^{\frac{1}{2}} \exp\left( - g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2} (\mathbf{z}^*_j)^\top \underbrace{\left( \boldsymbol{\Lambda}^{-\frac{1}{2}} \left[ \boldsymbol{\Sigma}_i \right] \right)^\top  \boldsymbol{\Sigma}_i^{-1} \boldsymbol{\Lambda}^{-\frac{1}{2}} \left[ \boldsymbol{\Sigma}_i \right]}_{= \mathbb{I}_{q \times q}} \mathbf{z}^*_j \right) \right] \\
&= \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \left[ (2 \pi \sigma^2)^{- \frac{n_i}{2}} \left\rvert \sigma^2 \mathbf{D} \right\rvert^{-\frac{1}{2}} \left\rvert \boldsymbol{\Sigma}_i \right\rvert^{\frac{1}{2}} \exp\left( - g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2} \rvert \rvert \mathbf{z}^*_j\rvert \rvert_2^2 \right) \right]
\end{aligned}
$$

We then obtain the importance sampling approximation by extending the above to all of the clusters and taking the natural logarithm:

$$
\begin{aligned}
\log\left( \prod_{i = 1}^K p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})  \right)
&\approx \log\left(  \prod_{i = 1}^K \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \left[ (2 \pi \sigma^2)^{- \frac{n_i}{2}} \left\rvert \sigma^2 \mathbf{D} \right\rvert^{-\frac{1}{2}} \left\rvert \boldsymbol{\Sigma}_i \right\rvert^{\frac{1}{2}} \exp\left( - g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2} \rvert \rvert \mathbf{z}^*_j\rvert \rvert_2^2 \right) \right]  \right) \\
&= \sum_{i = 1}^K \log\left( (2 \pi \sigma^2)^{-\frac{n_i}{2}} \left\rvert \sigma^2 \mathbf{D} \right\rvert^{-\frac{1}{2}} \left\rvert \boldsymbol{\Sigma}_i \right\rvert^{\frac{1}{2}} \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \exp\left( - g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2} \rvert \rvert \mathbf{z}^*_j\rvert \rvert_2^2 \right) \right) \\
&= \sum_{i = 1}^K \left[ - \frac{n_i}{2} \log(2 \pi \sigma^2) - \frac{1}{2} \log(\rvert \sigma^2 \mathbf{D} \rvert) + \frac{1}{2} \log(\rvert \boldsymbol{\Sigma}_i \rvert) + \log\left( \frac{1}{n^{IS}} \sum_{j = 1}^{n^{IS}} \exp\left( - g(\boldsymbol{\beta}_i(\mathbf{z}^*_j)) + \frac{1}{2} \rvert \rvert \mathbf{z}^*_j\rvert \rvert_2^2 \right) \right) \right] 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

The above can be simplified by using the approximation to $g''(\hat{\boldsymbol{\beta}}_i)$ that we found in the preivous section.

Because importance sampling relies on taking random samples, there can be numerical issues if we sampling too far into the tails of the proposal distribution. This can be overcome by replacing the choice of $q(\cdot)$ with something with fatter tails (like a Student's $t$-distribution with the same mean and a large number of degrees of freedom). Alternatively, we can forego sampling altogether and use <i>quadrature</i>. 

---

## Adaptive Gaussian Quadrature Approximation
Because we have a multi-dimensional integral, we have to repeatedly apply one-dimensional quadrature rules, but this can still work well in practice. First, let's cover "standard" Gaussian quadrature where we do not concern ourselves with finding "good" quadrature nodes. 
### Gaussian Quadrature
Recall that <a href="/stats-ml/quadrature/">Gauss-Hermite quadrature</a> makes the following approximation:

$$
\begin{aligned}
\int f(x) \exp(-x^2) dx \approx \sum_{j = 1}^m w_j f(\nu_j)
\end{aligned}
$$

where $\nu_j$ and $w_j$ are the abscissae of the $m$-th order Hermite polynomial and the associated weights, respectively. This approximation is essentially the same as importance sampling where our proposal distribution is a standard Gaussian (and we transform the integrand suitably to match). However, instead of approximating the expectation (with respect to this proposal distribution) as the empirical mean over some random sample, we use a deterministic set of points (i.e. the roots of a Hermite polynomial) and also weight the summands.

To use quadrature, we will need to manipulate the integral into the correct form. As we did in the previous section, we can rewrite the random effects as:

$$
\begin{aligned}
\mathbf{z}^* &\sim \mathcal{N}(\mathbf{0}_q, \mathbb{I}_{q \times q}) \\
\boldsymbol{\beta}_i &= \boldsymbol{\Lambda}^{\frac{1}{2}} \left[ \sigma^2 \mathbf{D} \right] \mathbf{z}^*;
\end{aligned}
$$

To make this one-dimensional, we can equivalently define $$\mathbf{z}^* = (z^*_1, \dots, z^*_q)^\top$$ for $$z^*_j \sim \mathcal{N}(0, 1)$$. 

<!-- #region gq-marg-i -->
{% tabs gq-marg-i %}
{% tab gq-marg-i statement %}
We can then rewrite the marginal for a single cluster as:

$$
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
= \idotsint (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \bar{\mathbf{z}}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{j = 1}^q \exp\left(- (\bar{z}^*_j)^2 \right) d \bar{z}^*_1 \dots d \bar{z}^*_q
$$

where $$\bar{z}^*_j = \sqrt{2} z^*_j$$.
{% endtab %}
{% tab gq-marg-i proof %}
For a single cluster, we then have:

$$
\begin{aligned}
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \int p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D} \rvert \boldsymbol{\beta}_i) p(\boldsymbol{\beta}_i; \sigma^2, \mathbf{D})   d \boldsymbol{\beta}_i \\
&= \int (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left(- \frac{1}{2 \sigma^2} \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i, \boldsymbol{\gamma}_i) \rvert \rvert_2^2 - \frac{1}{2} \boldsymbol{\beta}_i^\top \left[\sigma^2\mathbf{D}\right]^{-1}\boldsymbol{\beta}_i \right) d \boldsymbol{\beta}_i \\
&= \int (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left(- \frac{1}{2 \sigma^2} \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i(\mathbf{z}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2 - \frac{1}{2} \left(\boldsymbol{\Lambda}^{\frac{1}{2}} \left[ \sigma^2 \mathbf{D} \right] \mathbf{z}^*\right)^\top \left[ \sigma^2 \mathbf{D}\right]^{-1} \left(\boldsymbol{\Lambda}^{\frac{1}{2}} \left[ \sigma^2 \mathbf{D} \right] \mathbf{z}^*\right) \right) d \mathbf{z}^* \\
&= \int (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i(\mathbf{z}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2 - \frac{1}{2} (\mathbf{z}^*)^\top \underbrace{\left(\boldsymbol{\Lambda}^{\frac{1}{2}} \left[ \sigma^2 \mathbf{D} \right] \right)^\top \left[ \sigma^2 \mathbf{D}\right]^{-1} \boldsymbol{\Lambda}^{\frac{1}{2}} \left[ \sigma^2 \mathbf{D} \right] }_{= \mathbb{I}_{q \times q}} \mathbf{z}^*\right)  d \mathbf{z}^* \\ 
&= \int (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i(\mathbf{z}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2 - \frac{1}{2} \rvert \rvert \mathbf{z}^* \rvert \rvert_2^2 \right) d \mathbf{z}^*  \\
&= \idotsint (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i(\mathbf{z}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2 - \frac{1}{2} \rvert \rvert \mathbf{z}^* \rvert \rvert_2^2 \right) d z^*_1 \dots d z^*_q \\
&= \idotsint (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i(\mathbf{z}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{j = 1}^q \exp\left(-\frac{(z^*_j)^2}{2} \right) d z^*_1 \dots d z^*_q 
\end{aligned}
$$

Suppose we define $$\bar{z}^*_j = \sqrt{2} z^*_j$$. We can then rewrite the marginal as:

$$
\begin{aligned}
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \idotsint (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \bar{\mathbf{z}}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{j = 1}^q \exp\left(-\frac{(\sqrt{2} \bar{z}^*_j)^2}{2} \right) d \bar{z}^*_1 \dots d \bar{z}^*_q  \\
&= \idotsint (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \bar{\mathbf{z}}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{j = 1}^q \exp\left(- (\bar{z}^*_j)^2 \right) d \bar{z}^*_1 \dots d \bar{z}^*_q 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

We can use quadrature to approximate each dimension of the multiple integral separately.

<!-- #region gh-quad-1 -->
{% tabs gh-quad-1 %}
{% tab gh-quad-1 statement %}
Define $$\boldsymbol{\nu}_{j_1, \dots, j_q} = (\nu_{j_1}, \dots, \nu_{j_q})^\top$$. The Gauss-Hermite quadrature approximation is:

$$
\ell_{GH}(\boldsymbol{\alpha}, \sigma^2, \mathbf{D}; \mathbf{y})
= \sum_{i = 1}^K \left[- \frac{n_i + q}{2} \log(2 \pi \sigma^2) - \frac{1}{2} \log(\rvert \mathbf{D} \rvert) + \log \left(\sum_{j_1 = 1}^{n^{GQ}} \dots \sum_{j_q = 1}^{n^{GQ}} \exp\left( - \frac{1}{2 \sigma^2} \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{l = 1}^q w_{j_l}  \right)  \right]
$$
{% endtab %}
{% tab gh-quad-1 proof %}
Define $$\boldsymbol{\nu}_{j_1, \dots, j_q} = (\nu_{j_1}, \dots, \nu_{j_q})^\top$$. For a single cluster, Gauss-Hermite quadrature yields:

$$
\begin{aligned}
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&= \idotsint (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \bar{\mathbf{z}}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{j = 1}^q \exp\left(- (\bar{z}^*_j)^2 \right) d \bar{z}^*_1 \dots d \bar{z}^*_q \\
&= (2 \pi \sigma^2)^{-\frac{n_i}{2}} (2 \pi)^{- \frac{q}{2}} \rvert \sigma^2 \mathbf{D} \rvert^{- \frac{1}{2}} \underbrace{\idotsint \left[ \left[ \exp\left( - \frac{1}{2 \sigma^2}  \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \bar{\mathbf{z}}^*), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \exp(-(\bar{z}^*_1)^2) d \bar{z}^*_1 \right] \dots \right] \exp(- (\bar{z}^*_q)^2) d \bar{z}^*_q}_{\text{GHQ here}} \\
&\approx (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{- \frac{1}{2}} \sum_{j_1 = 1}^{n^{GQ}} \dots \sum_{j_q = 1}^{n^{GQ}} \exp\left( - \frac{1}{2 \sigma^2} \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{l = 1}^q w_{j_l}
\end{aligned}
$$

Extending to all of the clusters and taking the natural logarithm:

$$
\begin{aligned}
\log(p(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}))
&= \log\left( \prod_{i = 1}^K p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) \right)  \\
&\approx \sum_{i = 1}^K \log\left(  (2 \pi \sigma^2)^{-\frac{n_i + q}{2}} \rvert \mathbf{D} \rvert^{- \frac{1}{2}}  \sum_{j_1 = 1}^{n^{GQ}} \dots \sum_{j_q = 1}^{n^{GQ}} \exp\left( - \frac{1}{2 \sigma^2} \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{l = 1}^q w_{j_l}  \right) \\
&= \sum_{i = 1}^K \left[- \frac{n_i + q}{2} \log(2 \pi \sigma^2) - \frac{1}{2} \log(\rvert \mathbf{D} \rvert) + \log \left(\sum_{j_1 = 1}^{n^{GQ}} \dots \sum_{j_q = 1}^{n^{GQ}} \exp\left( - \frac{1}{2 \sigma^2} \rvert \rvert \mathbf{y}_i - \mathbf{f}_i(\boldsymbol{\phi}_i( 2^{-\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}), \boldsymbol{\gamma}_i) \rvert \rvert_2^2\right) \prod_{l = 1}^q w_{j_l}  \right)  \right]
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->


### Adaptive Gaussian Quadrature
<a href="/paper-notes/gauss-hermite/">Adaptive Gaussian quadrature</a> is very similar. Instead of directly using the weights and abscissae of a Hermite polynomial, we transform them so that the nodes lie in areas of high density (with respect to the posterior). 

For a one-dimensional integral, adaptive Gaussian quadrature makes the following approximation:

$$
\int p(x, y) dy = \int p(x \rvert y) p(y) dy \approx  \sqrt{2 \hat{\sigma}^2}  \sum_{i = 1}^{n^{AGQ}} w_i \exp(\nu_i^2) p(x, \hat{\mu} + \nu_i \sqrt{2 \hat{\sigma}^2})
$$

where $\hat{\mu}$ and $\hat{\sigma}^2$ are the mean and variance of Laplace approximation to the posterior, $p(y \rvert x)$, and $w_i$ and $\nu_i$ are the appropriate weights and abscissae. If, instead, $y$ is $q$-dimensional (we denote this version as $\mathbf{y}$), and we apply the approximation to each coordinate independently, we get:

$$
\int p(x, \mathbf{y}) d\mathbf{y}
\approx  \sqrt{2} \rvert \hat{\boldsymbol{\Sigma}} \rvert^{\frac{1}{2}} \sum_{j_1 = 1}^{n^{AGQ}} \dots \sum_{j_q = 1}^{n^{AGQ}} p(x, \hat{\boldsymbol{\mu}} + \sqrt{2} \rvert \hat{\boldsymbol{\Sigma}} \rvert^{\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}) \prod_{k = 1}^q c_{j_k} \exp(v_{j_q}^2) 
$$

where $$\boldsymbol{\nu}_{j_1, \dots, j_q} = (\nu_{j_1}, \dots, \nu_{j_q})^\top$$. 

<!-- #region gh-quad-2 -->
{% tabs gh-quad-2 %}
{% tab gh-quad-2 statement %}
The adaptive Gaussian quadrature approximation is:

$$
\ell_{AGQ}(\boldsymbol{\alpha}, \sigma^2, \mathbf{D}; \mathbf{y})
= \sum_{i = 1}^K \left[ \frac{1}{2} \log(2) + \frac{1}{2} \log\left( \rvert \boldsymbol{\Sigma}_i \rvert \right) + \log\left( \sum_{j_1 = 1}^{n^{AGQ}} \dots \sum_{j_q = 1}^{n^{AGQ}} \exp(- g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \hat{\boldsymbol{\mu}} + \sqrt{2} \rvert \hat{\boldsymbol{\Sigma}} \rvert^{\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}) + \rvert \rvert \boldsymbol{\nu}_{j_1, \dots, j_q} \rvert \rvert_2^2) \prod_{k = 1}^q c_{j_k}  \right)  \right]
$$
{% endtab %}
{% tab gh-quad-2 proof %}
As we found before, the posterior, $p(\boldsymbol{\beta}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D} \rvert \mathbf{y}_i)$, is (approximately) proportional to a Gaussian density with mean and covariance:

$$
\begin{aligned}
\boldsymbol{\mu}_i &= \hat{\boldsymbol{\beta}}_i \\
\boldsymbol{\Sigma}_i &= \left[ g''(\hat{\boldsymbol{\beta}}_i) \right]^{-1}
\end{aligned}
$$

For a single cluster, adaptive Gaussian quadrature yields:

$$
\begin{aligned}
p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D})
&\approx \sqrt{2} \rvert \boldsymbol{\Sigma}_i \rvert^{\frac{1}{2}} \sum_{j_1 = 1}^{n^{AGQ}} \dots \sum_{j_q = 1}^{n^{AGQ}} p(\mathbf{y}_i, \hat{\boldsymbol{\mu}} + \sqrt{2} \rvert \hat{\boldsymbol{\Sigma}} \rvert^{\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) \prod_{k = 1}^q c_{j_k} \exp(v_{j_q}^2) 
\end{aligned}
$$

Extending to all of the clusters and taking the natural logarithm:

$$
\begin{aligned}
\log(p(\mathbf{y}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}))
&= \log\left( \prod_{i = 1}^K p(\mathbf{y}_i; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) \right)  \\
&\approx \sum_{i = 1}^K \log\left( \sqrt{2} \rvert \boldsymbol{\Sigma}_i \rvert^{\frac{1}{2}} \sum_{j_1 = 1}^{n^{AGQ}} \dots \sum_{j_q = 1}^{n^{AGQ}} p(\mathbf{y}_i, \hat{\boldsymbol{\mu}} + \sqrt{2} \rvert \hat{\boldsymbol{\Sigma}} \rvert^{\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}; \boldsymbol{\alpha}, \sigma^2, \mathbf{D}) \prod_{k = 1}^q c_{j_k} \exp(v_{j_q}^2) \right) \\
&= \sum_{i = 1}^K \left[ \frac{1}{2} \log(2) + \frac{1}{2} \log\left( \rvert \boldsymbol{\Sigma}_i \rvert \right) + \log\left( \sum_{j_1 = 1}^{n^{AGQ}} \dots \sum_{j_q = 1}^{n^{AGQ}} \exp(- g(\mathbf{y}_i; \boldsymbol{\alpha}, \mathbf{D}, \hat{\boldsymbol{\mu}} + \sqrt{2} \rvert \hat{\boldsymbol{\Sigma}} \rvert^{\frac{1}{2}} \boldsymbol{\nu}_{j_1, \dots, j_q}) + \rvert \rvert \boldsymbol{\nu}_{j_1, \dots, j_q} \rvert \rvert_2^2) \prod_{k = 1}^q c_{j_k}  \right)  \right]
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

### Example
We now derive the adaptive Gaussian quadrature approximation to the integral. As in the Laplace approximation, we will define:

$$
\begin{aligned}
q(\beta_k) &= \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2)) \\
\hat{\mu}_k  &= \hat{\beta}_k^{(t+1)} \\
\hat{\sigma}_k^2 &= \left[ - q''(\hat{\beta}_k^{(t+1)}) \right]^{-1}
\end{aligned}
$$

The adaptive Gaussian quadrature approximation for one cluster is:

$$
\begin{aligned}
p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2)
&= \int p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k) p(\beta_k; \tau^2) d \beta_k \\
&= \int \exp(q(\beta_k)) d \beta_k \\
&\approx \sqrt{2 \hat{\sigma}_k^2} \sum_{i = 1}^{n^{AGQ}} w_i \exp(\nu_i^2) \exp(q(\hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2}))
\end{aligned}
$$

where $w_i$ and $\nu_i$ are the standard weights and abscissae of the $n^{AGQ}$-order Hermite polynomial. Thus, for the full data, we have:

$$
\begin{aligned}
\ell_{AGQ}(\boldsymbol{\alpha}, \gamma, \tau^2; \mathbf{y})
&= \log\left( \prod_{k = 1}^K \sqrt{2 \hat{\sigma}_k^2} \sum_{i = 1}^{n^{AGQ}} w_i \exp(\nu_i^2) \exp\left(q\left(\hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2}\right)\right) \right) \\
&= \sum_{k = 1}^K \left[\frac{1}{2}  \log\left( 2 \hat{\sigma}_k^2\right) + \log\left( \sum_{i = 1}^{n^{AGQ}} w_i \exp(\nu_i^2) \exp\left(q\left(\hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2}\right)\right) \right) \right]
\end{aligned}
$$




<!-- 
#### Updating Random Effects
Note that the mode, $\hat{\beta}_k$, is unknown. However, we can use the fact that $q'(\hat{\beta}_k) = 0$ to see that:

$$
\begin{aligned}
0 &= \left. \frac{d q(\beta_k)}{d \beta_k} \right\rvert_{\beta_k = \hat{\beta}_k} \\
\implies 0 &= - \frac{\hat{\beta_k}}{\tau^2} + \sum_{i = 1}^{n_k} z_{k,i} \left[ y_{k,i} - \frac{\hat{\mu}_{k,i}}{\hat{\mu}_{k,i} + \gamma} - \frac{\gamma \hat{\mu}_{k,i}}{\hat{\mu}_{k,i} + \gamma} \right] \\
\implies \hat{\beta}_k &= \tau^2 \sum_{i = 1}^{n_k} z_{k,i} \left[ y_{k,i} - \frac{\hat{\mu}_{k,i}}{\hat{\mu}_{k,i} + \gamma} - \frac{\gamma \hat{\mu}_{k,i}}{\hat{\mu}_{k,i} + \gamma} \right]
\end{aligned}
$$

We can then iteratively compute $\hat{\beta}_k$ by first initializing some $\hat{\beta}_k^{(0)}$ and computing for $t = 0, 1, \dots$ until convergence:

$$
\begin{aligned}
\hat{\mu}_{k,i}^{(t)} &= \exp\left( \mathbf{x}_{k,i}^\top \boldsymbol{\alpha} + z_{k,i} \hat{\beta}_k^{(t)} \right) \\
\hat{\beta}_{k,i}^{(t+1)} &= \tau^2 \sum_{i = 1}^{n_k} z_{k,i} \left[ y_{k,i} - \frac{\hat{\mu}^{(t)}_{k,i}}{\hat{\mu}^{(t)}_{k,i} + \gamma} - \frac{\gamma \hat{\mu}^{(t)}_{k,i}}{\hat{\mu}^{(t)}_{k,i} + \gamma} \right]
\end{aligned}
$$

We can do this simultaneously for all $k = 1, \dots, K$. 

#### Updating Other Parameters
The previous step yields a posterior mode for the random effects assuming we know the values of $\boldsymbol{\alpha}$, $\gamma$, and $\tau^2$. Now we consider the case when we know the value of $\boldsymbol{\beta}$ and update these other parameters.  -->

<!-- We can the update the other parameters, given a previous iteration's posterior mdoe, $\hat{\beta}_k^{(t)}$:

$$
\hat{\boldsymbol{\alpha}}^{(t)}, \hat{\gamma}^{(t)}, {\hat{\tau}^2}^{(t)} = \underset{\boldsymbol{\alpha}, \gamma, \tau^2}{\arg\max} \left\{ \sqrt{2 \pi} \exp\left(  \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \hat{\beta}^{(t)}_k)) + \log(p(\hat{\beta}^{(t)}_k; \tau^2)) \right) \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \hat{\mu}^{(t)}_{k,i}}{(\hat{\mu}^{(t)}_{k,i} + \gamma)^2} \right]^{-\frac{1}{2}} \right\}
$$ -->


<!-- 

### Details
Let $\hat{\boldsymbol{\alpha}}^{(t)}$, $\hat{\gamma}^{(t)}$, $(\hat{\tau}^2)^{(t)}$, and $\hat{\boldsymbol{\beta}}^{(t)}$ be the values of $\boldsymbol{\alpha}$, $\gamma$, $\tau^2$, and $\boldsymbol{\beta}$ at iteration $t$. 

#### Updating Random Effects


Because of our assumptions, we can easily find the value of $\hat{\boldsymbol{\beta}}^{(t+1)}$ as:

$$
\begin{aligned}
\hat{\boldsymbol{\beta}}^{(t+1)} 
&= \underset{\boldsymbol{\beta}}{\arg \max} \left\{ \log( p(\mathbf{y}; \hat{\boldsymbol{\alpha}}^{(t)}, \hat{\gamma}^{(t)} \rvert \boldsymbol{\beta}) p(\boldsymbol{\beta}; (\hat{\tau}^{(t)})^2) ) \right\} \\
&= \underset{\boldsymbol{\beta}}{\arg \max} \left\{ 
\sum_{k = 1}^K \left[ \sum_{i = 1}^{n_k} \left[ \log\left( \frac{\Gamma(y_{k,i} + \hat{\gamma}^{(t)} + 2)}{y_{k,i}! \Gamma(\hat{\gamma}^{(t)} + 2)} \right) + y_{k,i}\left( \log(\hat{\mu}^{(t)}_{k,i}) - \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right) + \hat{\gamma}^{(t)}\left( \log(\hat{\gamma}^{(t)}) - \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right) \right] - \frac{1}{2} \log(2 \pi (\hat{\tau}^2)^{(t)}) - \frac{\beta_k^2}{2 (\hat{\tau}^2)^{(t)}} \right]
\right\}
\end{aligned}
$$

where:

$$
\begin{aligned}
\hat{\mu}_{k,i}^{(t)} &= \exp\left( \mathbf{x}_{k,i}^\top \hat{\boldsymbol{\alpha}}^{(t)} + z_{k,i} \beta_k \right)
\end{aligned}
$$

Any term that does not involve $\boldsymbol{\beta}$ can be omitted from the objective function without changing the solution, so we can instead optimize:

$$
\begin{aligned}
\hat{\boldsymbol{\beta}}^{(t+1)} 
&= \underset{\boldsymbol{\beta}}{\arg \max} \left\{ 
\sum_{k = 1}^K \left[ \sum_{i = 1}^{n_k} \left[ y_{k,i}\left( \log(\hat{\mu}^{(t)}_{k,i}) - \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right) - \hat{\gamma}^{(t)} \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right] - \frac{\beta_k^2}{2 (\hat{\tau}^2)^{(t)}} \right] 
\right\}
\end{aligned}
$$

Furthermore, note that because each $\beta_k$ only appears in one of the summands over $k = 1, \dots, K$, we can perform the optimization over each cluster independently.  -->

<!-- 
#### Summary
Given the parameter values at iteration $t$ (i.e. $\hat{\boldsymbol{\alpha}}^{(t)}$, $\hat{\gamma}^{(t)}$, $(\hat{\tau}^2)^{(t)}$, and $\hat{\boldsymbol{\beta}}^{(t)}$), and the weights and points (i.e. $w_i$ and $\nu_i$ for $i = 1, \dots, n^{AGQ}$ for an $n^{AGQ}$-order quadrature-based approximation), we can then explicitly write out all of the parameter updates at iteration $t+1$.

First, we define:

$$
\begin{aligned}
\hat{\mu}_{k,i}^{(t)} &= \exp\left( \mathbf{x}_{k,i}^\top \hat{\boldsymbol{\alpha}}^{(t)} + z_{k,i} \beta_k \right)
\end{aligned}
$$

The random effects update is given by:

$$
\hat{\boldsymbol{\beta}}^{(t+1)} 
= \underset{\boldsymbol{\beta}}{\arg \max} \left\{ \sum_{k = 1}^K \sum_{i = 1}^{n_k} \left[ y_{k,i}\left( \log(\hat{\mu}^{(t)}_{k,i}) - \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right) - \hat{\gamma}^{(t)} \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right] - \frac{\beta_k^2}{2 (\hat{\tau}^2)^{(t)}} \right\}
$$

Now that we have $\hat{\beta}_k^{(t+1)}$, we define:

$$
\begin{aligned}
\hat{\mu}_k 
&= \hat{\beta}_k^{(t+1)} \\
\hat{\sigma}_k^2
&= \left[ - \left( - \left[\frac{1}{(\hat{\tau}^2)^{(t)}} + \sum_{i = 1}^{n_k} \frac{(1 + \hat{\gamma}^{(t)}) z^2_{k,i} \hat{\mu}^{(t)}_{k,i}}{(\mu_{k,i} + \hat{\gamma}^{(t)})^2} \right]\right) \right]^{-1}
= \frac{(\hat{\tau}^2)^{(t)}}{1 + (\hat{\tau}^2)^{(t)} (1 + \hat{\gamma}^{(t)}) \sum_{i = 1}^{n_k} \frac{z^2_{k,i} \hat{\mu}^{(t)}_{k,i}}{(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)})^2}} \\
q\left(\hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2}\right)
&= \log\left(p\left(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2} \right)\right) + \log\left(p\left(\hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2}; \tau^2 \right) \right)
\end{aligned}
$$

The other parameter updates are then:

$$
\hat{\boldsymbol{\alpha}}^{(t+1)}, \hat{\gamma}^{(t+1)}, (\hat{\tau}^2)^{(t+1)} 
= \underset{\boldsymbol{\alpha}, \gamma, \tau^2}{\arg \max} \left\{ \sum_{k = 1}^K \left[\frac{1}{2}  \log\left( 2 \hat{\sigma}_k^2\right) + \log\left( \sum_{i = 1}^{n^{AGQ}} w_i \exp(\nu_i^2) \exp\left(q\left(\hat{\mu}_k + \nu_i \sqrt{2 \hat{\sigma}_k^2}\right)\right) \right) \right] \right\}
$$ 
 -->
