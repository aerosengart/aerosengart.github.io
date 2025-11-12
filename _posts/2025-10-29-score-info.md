---
layout: post
title:  "Score and Information - GLMMs"
date: 29 October 2025
categories: posts
tags: ["likelihood", "information", "generalized-models", "score-test"]
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

This post is just a catch-all for my derivations for my score test project. Our set-up is as follows. We have $n$ observations coming from $k$ different clusters, each of size $n_t$ for $t \in [k]$. The full data will be denoted by $\mathbf{y}$. Though $\mathbf{y}$ is a vector, we'll denote the $j$-th observation from cluster $i$ with $\mathbf{y}_{i,j}$. For example, $$\mathbf{y}_{i,j}$$ denotes element $$\sum_{l = 1}^{i - 1} n_l + j$$ of $\mathbf{y}$. We'll also denote the $n_i$-dimensional vector of responses for cluster $i$ with $\mathbf{y}_i$. 

For each observation, we will have $p$ fixed effect covariates arranged in a $p$-dimensional vector, $$\mathbf{x}_{i, j}$$, and $q$ random effects covariates in a $q$-dimensional vector, $$\mathbf{z}_{i,j}$$. We'll assume that the observations within the same cluster are independent.

Our model comes in the form of a specification of the conditional mean, $\mu_{i,j} = \mathbb{E}[\mathbf{y}_{i,j} \rvert \beta_i]$ (where we suppress the addition conditioning on the covariates themselves). For a monotonic and differentiable link function (e.g. $\log(\cdot)$ or $\text{logit}(\cdot)$), the conditional mean of the $j$-th observation in group $i$ is assumed to be given by:

$$
\mu_{i,j} = g^{-1}\left(\alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j} \right)
\label{eq:glmm}
$$

We then assume that the observations themselves follow some exponential family distribution with measurement errors, $\epsilon_{i,j}$, which is the deviation of the response from its (unit-specific) conditional mean. These errors are assumed to have mean zero and be independent of each other and of the random effects. We further assume the responses, $\mathbf{y}_{i,j}$, conditional on the random effects (and the covariates), are independent with variances equal to some function of the conditional mean. 

In general, we will assume that:

$$
\beta_i \overset{iid}{\sim} \mathcal{N}\left(\mathbf{0}_q, D(\tau^2) \right)
$$

for some variance component, $\tau^2$. We'll use $[\cdot] \rvert_{H_0}$ to denote evaluation of the function in brackets when setting $\beta$ equal to $\beta_0$. We'll also use a superscript $0$ (e.g. $\mu^0$, $\eta^0$, etc.) to denote the quantity under the null hypothesis (i.e. $\tau^2 = \mathbf{0} \implies \beta = \mathbf{0}$).


---

## Gaussian Case
In this example, we'll have the simple setting of a Gaussian response, which means $g(\cdot)$ is the identity function. We will have a fixed (but cluster-specific) intercept and a single random slope. We will have $k$ clusters and $n$ observations per cluster. We assume:

$$
\mathbf{y}_{i, j} = \alpha_i + \beta_i \mathbf{z}_{i,j} + \epsilon_{i,j}, 
\hspace{8mm} \epsilon_{i,j} \overset{iid}{\sim} \mathcal{N}(0, \sigma^2), 
\hspace{5mm} \beta_i \overset{iid}{\sim} \mathcal{N}(0, \tau^2)
$$

where we also assume the random effects and errors are independent. $$\mathbf{z}_i \in \mathbb{R}^n$$ is the vector of covariate values for the $n$ samples in cluster $i$. We'll denote the vector of responses for cluster $i$ with $\mathbf{y}_i$ so that $$\mathbf{y}_{i,j}$$ denotes the $j$-th component of said vector. Marginally, the response vector $\mathbf{y}_i$ has mean $$\alpha_i \mathbb{1}_n$$ and variance-covariance matrix:

$$
\Sigma_{y_i} =
\sigma^2 \mathbb{I}_{n \times n} + \tau^2 \mathbf{z}_i \mathbf{z}_i^\top
$$

<details>
<summary>Proof.</summary>
For a single cluster:
$$
\begin{aligned}
\mathbb{E}\left[ (\mathbf{y}_{i,j} - \alpha_i)^2\right] 
&= \mathbb{E}\left[ (\beta_i \mathbf{z}_{i,j} + \epsilon_{i,j})^2 \right] \\
&= \mathbb{E}\left[\beta_i^2 \mathbf{z}_{i,j}^2 \right] + 2 \mathbb{E}\left[ \beta_i \mathbf{z}_{i,j} \epsilon_{i,j} \right] + \mathbb{E}\left[ \epsilon_{i,j}^2 \right] \\
&= \tau^2 \mathbf{z}_{i,j}^2 + \sigma^2 \\
\mathbb{E}\left[ (\mathbf{y}_{i,j} - \alpha_i)(\mathbf{y}_{i,j'} - \alpha_i) \right]
&= \mathbb{E}\left[ (\beta_i \mathbf{z}_{i,j} + \epsilon_{i,j})(\beta_i \mathbf{z}_{i,j'} + \epsilon_{i,j'})\right] \\
&= \mathbb{E}\left[ \beta_i^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \right] + \mathbb{E}\left[ \beta_i \mathbf{z}_{i,j} \epsilon_{i,j'}\right] + \mathbb{E}\left[ \beta_i \mathbf{z}_{i,j'} \epsilon_{i,j}\right] + \mathbb{E}\left[ \epsilon_{i,j} \epsilon_{i,j'}\right] \\
&= \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}
\nonumber
\end{aligned}
$$
Thus, the variance-covariance matrix for $\mathbf{y}_i$:
$$
\Sigma_{y_i} = 
\begin{bmatrix}
\sigma^2 + \tau^2 \mathbf{z}_{i,1}^2 & \dots & \tau^2 \mathbf{z}_{i,1} \mathbf{z}_{i,n} \\
\vdots & \ddots & \vdots \\
\tau^2 \mathbf{z}_{i,n} \mathbf{z}_{i, 1} & \dots & \sigma^2 + \tau^2 \mathbf{z}_{i, n}^2
\end{bmatrix} =
\sigma^2 \mathbb{I}_{n \times n} + \tau^2 \mathbf{z}_i \mathbf{z}_i^\top
\nonumber
$$
</details>

Since the $\beta_i$ are independent, observations from different clusters have covariance zero. Let $\mathbf{y} = (\mathbf{y}_1, \dots, \mathbf{y}_k)$ denote the full data, $\alpha = \begin{bmatrix} \alpha_1 & \dots & \alpha_k\end{bmatrix}^\top$, $\beta = \begin{bmatrix} \beta_1 & \dots & \beta_k\end{bmatrix}^\top$, and $\theta = (\alpha, \beta)$. The complete, marginal likelihood and log-likelihood are:

$$
\begin{aligned}
\mathcal{L}(\theta; \mathbf{y}) &= \prod_{i = 1}^k (2 \pi)^{-\frac{n}{2}} \rvert \Sigma_{y_i} \rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right) \\
\ell(\theta; \mathbf{y}) &= \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
$$

### Score and Information
It is easiest to write the score after evaluating it at the MLE of the parameter vector under $H_0$, which we denote with $\hat{\theta}$ (uncollapse the proof to see all of the details). The MLE is given by:

$$
\hat{\theta} = \begin{bmatrix}
    \frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{1,j} \\
    \vdots \\
    \frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{k,j} \\
    \frac{1}{nk} \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    0
\end{bmatrix}
$$

Thus, the score evaluated at $\theta = \hat{\theta}$ is:

$$
\begin{aligned}
U_\theta (\hat{\theta})
&= \begin{bmatrix}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \bigg\rvert_{\theta = \hat{\theta}} \\
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \bigg\rvert_{\theta = \hat{\theta}} \\
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \bigg\rvert_{\theta = \hat{\theta}}
\end{bmatrix} 
= \begin{bmatrix}
\frac{1}{\hat{\sigma}^2} (\mathbf{y}_1 - \hat{\alpha}_1 \mathbf{1}_n)^\top \mathbf{1}_n \\
\vdots \\
\frac{1}{\hat{\sigma}^2}  (\mathbf{y}_k - \hat{\alpha}_k \mathbf{1}_n)^\top \mathbf{1}_n \\
- \frac{1}{2} \sum_{i = 1}^k \left[ \frac{n}{\hat{\sigma}^2} - \frac{1}{(\hat{\sigma}^2)^2} (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n)^\top (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n) \right] \\
-\frac{1}{2} \sum_{i = 1}^k \left[ \frac{\text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right]}{\hat{\sigma}^2} + \frac{1}{(\hat{\sigma}^2)^2}(\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n)^\top \mathbf{z}_i \mathbf{z}_i^\top (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n) \right] 
\end{bmatrix}
\end{aligned}
$$

<details>
<summary>Proof.</summary>
We first find the gradient of the log-likelihood with respect to $\theta$ parameter-wise. Using the Sherman-Morrison formula, we can find $\Sigma_{y_i}^{-1}$ to be:
$$
\Sigma_{y_i}^{-1}  = \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top
\nonumber
$$
<details>
<summary>Proof.</summary>
$$
\begin{aligned}
  \Sigma_{y_i}^{-1} 
    &= \left[ \sigma^2 \mathbb{I}_{n \times n} + \mathbf{z}_i [\tau^2 \mathbb{I}_{n \times n}] \mathbf{z}_i^\top \right]^{-1}  \\
    &= \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \left(1 + \tau^2 \mathbf{z}_i^\top \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right]\mathbf{z}_i  \right)^{-1} \left( \left(\frac{1}{\sigma^2} \mathbb{I}_{n \times n}\right) \left(\tau^2 \mathbf{z}_{i} \mathbf{z}_i^\top \right) \left(\frac{1}{\sigma^2} \mathbb{I}_{n \times n}\right) \right) \\
    &= \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \left(1 + \frac{\tau^2}{\sigma^2} \mathbf{z}_i^\top \mathbf{z}_i \right)^{-1}\left(\frac{\tau^2}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top \right) \\
    &= \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top
\end{aligned}
\nonumber
$$
</details>
First, let's find the derivative with respect to $\sigma^2$:
$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2}
    = -\frac{1}{2}\sum_{i = 1}^k \left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\nonumber
$$
<details>
<summary>Proof.</summary>
$$ 
\begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \log(\rvert \Sigma_y \rvert) \right] 
    &= \text{tr}\left[ \Sigma_y^{-1} \frac{\partial}{\partial \sigma^2} \left[\Sigma_y\right] \right] \\
    &= \text{tr}\left[ \Sigma_y^{-1} \mathbb{I}_{n \times n} \right] \\
    &= \text{tr}\left[ \Sigma^{-1} \right] \\
    &= \text{tr}\left[\frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &=\text{tr}\left[ \frac{1}{\sigma^2}\mathbb{I}_{n \times n} \right] \text{tr}\left[- \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top\right] \\
    &=\frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]
\end{aligned}
\nonumber
$$
$$ 
\begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \Sigma_y^{-1} \right]
    &= - \Sigma_y^{-1} \frac{\partial}{\partial \sigma^2} \left[ \Sigma_y\right] \Sigma_y^{-1} \\
    &= -\left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbb{I}_{n \times n} \left[\frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= - \left[ \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} - \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top + \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top 
\end{aligned}
\nonumber
$$
The above imply:
$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2}
    &= \frac{\partial}{\partial \sigma^2} \left[ \sum_{i = 1}^k - \frac{1}{2} \log(\rvert \Sigma_y \rvert) - \frac{1}{2}(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \Sigma_y^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2}\sum_{i = 1}^k \left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
\nonumber
$$
</details>
We do the same with $\tau^2$:
$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2}
  = - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\nonumber
$$
<details>
<summary>Proof.</summary>
$$
\begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \log(\rvert \Sigma_y \rvert) \right] 
    &= \text{tr}\left[ \Sigma_y^{-1} \frac{\partial}{\partial \tau^2} \left[\Sigma_y\right] \right] \\
    &= \text{tr}\left[ \Sigma_y^{-1} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= \text{tr} \left[ \left( \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top\right) \mathbf{z}_i \mathbf{z}_i^\top\right] \\
    &= \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right]
\end{aligned}
\nonumber
$$
$$ 
\begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \Sigma_y^{-1} \right]
    &= - \Sigma_y^{-1} \frac{\partial}{\partial \tau^2} \left[ \Sigma_y\right] \Sigma_y^{-1} \\
    &= -\left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \left[\frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top
\end{aligned}
\nonumber
$$
$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2}
    &= \frac{\partial}{\partial \tau^2} \left[ \sum_{i = 1}^k - \frac{1}{2} \log(\rvert \Sigma_y \rvert) - \frac{1}{2}(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \Sigma_y^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
\nonumber
$$
</details>
And then take the gradient with respect to $\alpha$:
$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} 
= \begin{bmatrix}
  (\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \Sigma_y^{-1} \mathbf{1}_n \\
  \vdots \\
  (\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \Sigma_y^{-1} \mathbf{1}_n
\end{bmatrix}
\nonumber
$$
<details>
<summary>Proof.</summary>
We do the computations component-wise:
$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} 
    &= \sum_{i = 1}^k - \frac{1}{2} \frac{\partial}{\partial \alpha_j} \left[ (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \Sigma^{-1}_y (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= - \frac{1}{2} \left(2 (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_y^{-1}(- \mathbf{1}_n) \right) \\
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_y^{-1} \mathbf{1}_n
\end{aligned}
$$
So then:
$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} 
    &= \begin{bmatrix}
        (\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \Sigma_y^{-1} \mathbf{1}_n \\
        \vdots \\
        (\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \Sigma_y^{-1} \mathbf{1}_n
    \end{bmatrix}
\end{aligned}
$$
</details>
We can then find the MLE vector, $\hat{\theta}$, by setting the above equations equal to zero, substituting $\tau^2 = 0$, and solving. The MLE vector, $\hat{\theta}$ is:
$$ 
\hat{\theta} = \begin{bmatrix}
    \frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{1,j} \\
    \vdots \\
    \frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{k,j} \\
    \frac{1}{nk} \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    0
\end{bmatrix}
\nonumber
$$
<details>
<summary>Proof.</summary>
We set the derivative with respect to $\sigma^2$ equal to zero, substitute $\tau^2 = 0$ (under $H_0$), and solve for $\sigma^2$:
$$ 
\begin{aligned}
    0 &= \frac{\partial}{\partial \sigma^2} \left[ \ell(\theta; \mathbf{y}) \right] \bigg\rvert_{\theta = \theta_0}\\
    0 &= -\frac{1}{2}\sum_{i = 1}^k \left[ \frac{n}{\sigma^2} - \frac{0}{\sigma^2(\sigma^2 + 0 \cdot\mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\cdot 0}{(\sigma^2)^2(\sigma^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top + \frac{(0)^2}{(\sigma^2)^2(\sigma^2 - 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    0 &= -\frac{1}{2}\sum_{i = 1}^k \left[ \frac{n}{\sigma^2} + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n)\right] \\
    0 &= - \frac{nk}{2 \sigma^2} - \frac{1}{2} \sum_{i = 1}^k - \frac{1}{(\sigma^2)^2} (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    0 &= - \frac{nk}{2\sigma^2} + \frac{1}{2 (\sigma^2)^2} \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    \frac{nk}{2 \sigma^2} &= \frac{1}{2(\sigma^2)^2} \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\ 
    \sigma^2 n k &=  \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    \sigma^2 &= \frac{1}{nk} \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n)
\end{aligned}
\nonumber
$$
We set the gradient w.r.t $\alpha$ equal to zero, substitute $\tau^2 = 0$ (under $H_0$), and solve for $\alpha$.
$$ 
\begin{aligned}
    \mathbf{0} &= \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \bigg\rvert_{\theta = \theta_0} \\
    \mathbf{0} &= \begin{bmatrix}
        (\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \left[\frac{1}{\sigma^2}\mathbb{I}_{n \times n} \right] \mathbf{1}_n \\
        \vdots \\
        (\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \left[\frac{1}{\sigma^2}\mathbb{I}_{n \times n} \right] \mathbf{1}_n
    \end{bmatrix}
\end{aligned}
\nonumber
$$
Since each entry of the gradient only has one component of $\alpha$, we can solve then all separately:
$$ 
\begin{aligned}
    0 &= (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right] \mathbf{1}_n \\
    0 &= \frac{1}{\sigma^2} (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \mathbf{1}_n \\
    0 &= \frac{1}{\sigma^2} \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha_i) \\
    0 &= \frac{1}{\sigma^2} \left( \sum_{j =1 }^n \mathbf{y}_{i,j} - n \alpha_i \right) \\
    n \alpha_i &= \sum_{j =1 }^n \mathbf{y}_{i,j} \\
    \alpha_i &= \frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{i,j}
\end{aligned}
$$
</details>
Plugging $\hat{\theta}$ into the derivatives we found above, we get:
$$
\begin{aligned}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \theta} \bigg\rvert_{\theta = \hat{\theta}}
&= \begin{bmatrix}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \bigg\rvert_{\theta = \hat{\theta}} \\
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \bigg\rvert_{\theta = \hat{\theta}} \\
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \bigg\rvert_{\theta = \hat{\theta}}
\end{bmatrix} 
= \begin{bmatrix}
\frac{1}{\hat{\sigma}^2} (\mathbf{y}_1 - \hat{\alpha}_1 \mathbf{1}_n)^\top \mathbf{1}_n \\
\vdots \\
\frac{1}{\hat{\sigma}^2}  (\mathbf{y}_k - \hat{\alpha}_k \mathbf{1}_n)^\top \mathbf{1}_n \\
- \frac{1}{2} \sum_{i = 1}^k \left[ \frac{n}{\hat{\sigma}^2} - \frac{1}{(\hat{\sigma}^2)^2} (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n)^\top (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n) \right] \\
-\frac{1}{2} \sum_{i = 1}^k \left[ \frac{\text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right]}{\hat{\sigma}^2} + (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n)^\top \left[ \frac{\mathbf{z}_i \mathbf{z}_i^\top}{(\hat{\sigma}^2)^2} \right] (\mathbf{y}_i - \hat{\alpha}_i \mathbf{1}_n) \right] 
\end{bmatrix}
\end{aligned}
\nonumber
$$
where $\hat{\alpha}$ and $\hat{\sigma}$ denote the values of $\alpha$ and $\sigma$ in the MLE vector.
</details>

To find the information, we need to compute the second-order derivatives of the log-likelihood, take the expectation under $H_0$ of minus those quantities, and evaluate them by plugging in $\hat{\theta}$:

$$ 
\begin{aligned}
\mathcal{I}_{\theta, \theta} (\hat{\theta})
    &= -\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \theta \partial \theta^\top}\right]\bigg\rvert_{\theta = \hat{\theta}} 
    = \begin{bmatrix}
        \frac{n}{\hat{\sigma}^2} & \dots & 0 & 0 & 0\\
        \vdots & \ddots & \vdots & \vdots & \vdots \\
        0 & \dots & \frac{n}{\hat{\sigma}^2} & 0 & 0 \\
        0 & \dots & 0 & \frac{nk}{2\hat{\sigma}^2} & \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right]\\
        0 & \dots & 0 & \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] & \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right]
    \end{bmatrix}
\end{aligned}
$$

<details>
<summary>Proof.</summary>
We start by taking the derivative with respect to $\theta$ (component-wise) of the first derivative with respect to $\sigma^2$:
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] &=
- \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] &= 
-\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\sigma^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-2\tau^2(2\sigma^2 + 3\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
&= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ -\frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n 
\end{aligned}
\nonumber
$$
<details>
<summary>Proof.</summary>
$$ 
\begin{aligned}
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k \frac{\partial}{\partial \sigma^2}\left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} - \frac{-\tau^2(2\sigma^2+\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} + \frac{-2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    &= - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right]
    &= - \frac{1}{2} \sum_{i = 1}^k \frac{\partial}{\partial \sigma^2} \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\sigma^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-2\tau^2(2\sigma^2 + 3\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= \frac{\partial}{\partial \sigma^2} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \right] \\ 
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ -\frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n 
\end{aligned}
\nonumber
$$
</details>
Next, we do the same with the derivative with respect to $\tau^2$:
$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[\frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\end{aligned}
\nonumber
$$
<details>
<summary>Proof.</summary>
$$ 
\begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k \frac{\partial}{\partial \tau^2}\left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2}\sum_{i = 1}^k 0 - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[0 + \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[\frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right]
    &= - \frac{1}{2} \sum_{i = 1}^k \frac{\partial}{\partial \tau^2} \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2}\sum_{i = 1}^k 0 - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ 0 + \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= \frac{\partial}{\partial \tau^2} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \right] \\ 
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ 0 - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n 
\end{aligned}
\nonumber
$$
</details>
And finally with $\alpha_j$:
$$
\begin{aligned}
\frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right]
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_j^\top \mathbf{z}_j)} \mathbf{z}_j \mathbf{z}_j^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_j^\top \mathbf{z}_j)^2} \mathbf{z}_j \mathbf{z}_j^\top \mathbf{z}_j \mathbf{z}_j^\top \right]\mathbf{1}_n  \\
\frac{\partial}{\partial \alpha} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right]
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right]  \mathbf{1}_n \\
\frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= - \mathbf{1}_n^\top  \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\frac{\partial}{\partial \alpha_{j'}} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= 0 
\end{aligned}
$$
<details>
<summary>Proof.</summary>
$$ 
\begin{aligned}
\frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k \frac{\partial}{\partial \alpha_j} \left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2}\left[0 - 0 - 2(\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right]\mathbf{1}_n \right] \\
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_j^\top \mathbf{z}_j)} \mathbf{z}_j \mathbf{z}_j^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_j^\top \mathbf{z}_j)^2} \mathbf{z}_j \mathbf{z}_j^\top \mathbf{z}_j \mathbf{z}_j^\top \right]\mathbf{1}_n \\
\frac{\partial}{\partial \alpha} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right]
    &= - \frac{1}{2} \sum_{i = 1}^k \frac{\partial}{\partial \alpha_j} \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2} \left[0 - 0 - 2(\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \right] \\
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= \frac{\partial}{\partial \alpha_j} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \right] \\ 
    &=  \frac{\partial}{\partial \alpha_j} \left[ \sum_{h = 1}^n (\mathbf{y}_{j,h} - \alpha_j) \sum_{l = 1}^n \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right]_{h, l} \right] \\ 
    &= - \sum_{h = 1}^n \sum_{l = 1}^n  \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right]_{h, l} \\
    &= - \mathbf{1}_n^\top  \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\frac{\partial}{\partial \alpha_{j'}} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= \frac{\partial}{\partial \alpha_{j'}} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \right] \\ 
    &= 0 
\end{aligned}
\nonumber
$$
</details>
We take the expectation under the null of all of the terms we found in the previous section. First we do the $\sigma^2$ terms:
$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right] &=  
    - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right]
    &= -\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\sigma^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-2\tau^2(2\sigma^2 + 3\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right]
    &= 0 
\end{aligned}
\nonumber
$$
<details>
<summary>Proof.</summary>
$$ 
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right]
    &= \mathbb{E}\left[ - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \mathbb{E}\left[ (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &=  - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbb{E}\left[ (\mathbf{y}_i - \alpha_i \mathbf{1}_n) (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \right] \right] \\
    &= - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right]
    &= -\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\sigma^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-2\tau^2(2\sigma^2 + 3\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right]
    &= 0 
\end{aligned}
\nonumber
$$
</details>
And the $\tau^2$ terms:
$$ 
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right]
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr} \left[ \left[ \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right]
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] +\text{tr} \left[\left[ \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right] 
    &= 0
\end{aligned}
\nonumber
$$
And finally the $\alpha_j$ ones:
$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right] &= 0 \\
\mathbb{E}\left[ \frac{\partial}{\partial \alpha} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right] &= 0 \\
\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right] &= - \mathbf{1}_n^\top  \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\mathbb{E}\left[\frac{\partial}{\partial \alpha_{j'}} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right] &= 0 
\end{aligned}
$$
<details>
<summary>Proof.</summary>
$$ 
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right]
    &= \mathbb{E}\left[(\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_j^\top \mathbf{z}_j)} \mathbf{z}_j \mathbf{z}_j^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_j^\top \mathbf{z}_j)^2} \mathbf{z}_j \mathbf{z}_j^\top \mathbf{z}_j \mathbf{z}_j^\top \right]\mathbf{1}_n  \right] \\
    &= 0 \\
\mathbb{E}\left[ \frac{\partial}{\partial \alpha} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right]
    &= \mathbb{E}\left[(\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \right] \\
    &= 0 \\
\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right]
    &= - \mathbf{1}_n^\top  \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\mathbb{E}\left[\frac{\partial}{\partial \alpha_{j'}} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right]
    &= 0 
\end{aligned}
\nonumber
$$
</details>
We then evaluate the Fisher information at the MLEs. Note that:
$$ 
\begin{aligned}
    \Sigma_{y_i} \bigg\rvert_{\theta = \hat{\theta}} &= \hat{\sigma}^2 \mathbb{I}_{n \times n}
\end{aligned}
\nonumber
$$
Thus:
$$ 
\begin{aligned}
-\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right] \bigg\rvert_{\theta = \hat{\theta}}
    &= \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\hat{\sigma}^2)^2} + \frac{0 \cdot (2\hat{\sigma}^2 + 0)}{(\hat{\sigma}^2)^2(\hat{\sigma}^2 + 0)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \text{tr}\left[ \left[ \frac{2}{(\hat{\sigma}^2)^3} \mathbb{I}_{n \times n} - 0 \cdot \mathbf{z}_i \mathbf{z}_i^\top + 0 \cdot \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i}\rvert_{\theta = \hat{\theta}}  \right] \\
    &= \frac{1}{2}\sum_{i = 1}^k - \frac{n}{(\hat{\sigma}^2)^2} + \text{tr}\left[  \left[ \frac{2}{(\hat{\sigma}^2)^3} \mathbb{I}_{n \times n} \right] \left[ \hat{\sigma}^2 \mathbb{I}_{n \times n} \right] \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[ -\frac{n}{(\hat{\sigma}^2)^2} +  \frac{2}{\hat{\sigma}^2} \text{tr}\left[ \mathbb{I}_{n \times n} \right]\right] \\
    &= \frac{1}{2}\sum_{i = 1}^k \left[ -\frac{n}{(\hat{\sigma}^2)^2} + \frac{2n}{\hat{\sigma}^2} \right] \\
    &= \sum_{i = 1}^k \frac{n}{2 \hat{\sigma}^2} \\
    &= \frac{nk}{2\hat{\sigma}^2} \\
-\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right] \bigg\rvert_{\theta = \hat{\theta}}
    &=\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\hat{\sigma}^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-0 \cdot (2\hat{\sigma}^2 + 0 \cdot  \mathbf{z}_i^\top \mathbf{z}_i)}{(\hat{\sigma}^2)^2(\hat{\sigma}^2 + 0 \cdot  \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr}\left[ \left[ \frac{2}{(\hat{\sigma}^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-0 \cdot (2\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)}{(\hat{\sigma}^2)^3(\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(0)^2(2\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)}{(\hat{\sigma}^2)^3 (\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \rvert_{\theta = \hat{\theta}}  \right] \\
    &= \frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\hat{\sigma}^2)^2} \text{tr}\left[  \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr}\left[ \left[ \frac{2}{(\hat{\sigma}^2)^3} \mathbf{z}_i \mathbf{z}_i^\top \right] \left[ \hat{\sigma}^2 \mathbb{I}_{n \times n} \right] \right] \\
    &= \frac{1}{2}\sum_{i = 1}^k \left( - \frac{1}{(\hat{\sigma}^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \frac{2}{(\hat{\sigma}^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \right) \\
    &= \sum_{i = 1}^k \frac{1}{2(\hat{\sigma}^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k  \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \\
-\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \right] \right] \bigg\rvert_{\theta = \hat{\theta}} &= 0 \\
- \mathbb{E}\left[ \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right]\bigg\rvert_{\theta = \hat{\theta}}
    &= \frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr} \left[ \left[\frac{2}{\hat{\sigma}^2(\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top - \frac{2 \cdot 0 }{\hat{\sigma}^2(\hat{\sigma}^2 + 0 \cdot  \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i}\rvert_{\theta = \hat{\theta}} \right] \\
    &= \frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\hat{\sigma}^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \text{tr}\left[ \left[\frac{2}{(\hat{\sigma}^2)^3} \mathbf{z}_i \mathbf{z}_i^\top \right] \left[ \hat{\sigma}^2 \mathbb{I}_{n \times n} \right]\right] \\
    &= \frac{1}{2}\sum_{i =1 }^k \left( -\frac{1}{(\hat{\sigma}^2)^2} + \frac{2}{(\hat{\sigma}^2)^2} \right) \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr} \left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \\
-\mathbb{E}\left[ \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right]
    &= \frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\hat{\sigma}^2 + 0 \cdot\mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] +\text{tr} \left[\left[ \frac{2}{\hat{\sigma}^2(\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\cdot 0}{\hat{\sigma}^2(\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\right] \Sigma_{y_i} \rvert_{\theta = \hat{\theta}} \right] \\
    &= \frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\hat{\sigma}^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] +\text{tr} \left[\left[ \frac{2}{(\hat{\sigma}^2)^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \left[\hat{\sigma}^2 \mathbb{I}_{n \times n} \right] \right] \\
    &= \frac{1}{2}\sum_{i = 1}^k \left( -\frac{1}{(\hat{\sigma}^2)^2} + \frac{2}{(\hat{\sigma}^2)^2} \right)\text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  \right] \\
-\mathbb{E}\left[ \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \right] \right] \bigg\rvert_{\theta = \hat{\theta}} &= 0 \\ 
-\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right] \bigg\rvert_{\theta = \hat{\theta}}
    &= 0 \\
-\mathbb{E}\left[ \frac{\partial}{\partial \alpha} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right] \bigg\rvert_{\theta = \hat{\theta}}
    &= 0 \\
- \mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right]
    &= \mathbf{1}_n^\top  \left[ \frac{1}{\hat{\sigma}^2} \mathbb{I}_{n \times n} - \frac{0}{\hat{\sigma}^2(\hat{\sigma}^2 + 0 \cdot \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
    &= \frac{n}{\hat{\sigma}^2} \\
- \mathbb{E}\left[\frac{\partial}{\partial \alpha_{j'}} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right] \bigg\rvert_{\theta = \hat{\theta}}
    &= 0 
\end{aligned}
\nonumber
$$
Putting all of the above together into a matrix yields:
$$ 
\begin{aligned}
    \mathcal{I}_{\theta, \theta} (\hat{\theta})
    &= -\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \theta \partial \theta^\top}\right]\bigg\rvert_{\theta = \hat{\theta}} \\
    &= \begin{bmatrix}
        \frac{n}{\hat{\sigma}^2} & \dots & 0 & 0 & 0\\
        \vdots & \ddots & \vdots & \vdots & \vdots \\
        0 & \dots & \frac{n}{\hat{\sigma}^2} & 0 & 0 \\
        0 & \dots & 0 & \frac{nk}{2\hat{\sigma}^2} & \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right]\\
        0 & \dots & 0 & \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] & \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right]
    \end{bmatrix}
\end{aligned}
\nonumber
$$
</details>

---

## Negative Binomial Case
In this example, we'll let the responses be negative binomial. To keep things simple, we'll say we only have a single fixed intercept and a single random effect. We let $\phi > 0$, denote the <i>known</i> dispersion parameter and assume the conditional mean to be given by:

$$
\mu_{i,j} = \exp\left( \alpha_i + \beta_i \mathbf{z}_{i,j} \right)
\label{eq:neg-bin-mean}
$$

The likelihood based on a single observation, $\mathbf{y}_{i,j}$, is given by:

$$
\mathcal{L}(\mathbf{y}_{i,j}; \alpha_i, \tau^2 \rvert \beta_i) = 
\frac{\Gamma\left(\mathbf{y}_{i,j} + \frac{1}{\phi}\right)}{\Gamma(\mathbf{y}_{i,j} + 1) \Gamma\left(\frac{1}{\phi} \right)}\left(\frac{1}{1 + \phi \mathbf{y}_{i,j}}\right)^{\frac{1}{\phi}} \left( \frac{\phi \mu_{i,j}}{1 + \phi \mu_{i,j}} \right)^{\mathbf{y}_{i,j}}
\label{eq:neg-bin-single-lik}
$$

where $\Gamma(\cdot)$ is the gamma function:

$$
\Gamma(x) = \int_0^\infty t^{x - 1} \exp(-t) dt
$$

The above parametrization of the likelihood implies that the conditional variance of the responses is given by:

$$
V(\mu_{i,j}) = \mu_{i,j} + \frac{1}{\phi} \mu_{i,j}^2
$$

The conditional log-likelihood based on cluster $i$ is:

$$
\ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i) = \sum_{j = 1}^{n_i} \left[ \log \Gamma \left( \mathbf{y}_{i,j} + \frac{1}{\phi} \right)  - \log \Gamma\left(\mathbf{y}_{i,j} + 1\right) - \log\Gamma\left(\frac{1}{\phi} \right) - \frac{1}{\phi} \log\left(1 + \phi \mathbf{y}_{i,j} \right) + \mathbf{y}_{i,j} \left( \log(\phi \mu_{i,j}) - \log(1 + \phi \mu_{i,j}) \right) \right]
\label{eq:neg-bin-full-cond-ll}
$$

### Pseudo-Likelihood Approach
We follow a pseudo-likelihood approach (see <a href="/posts/2025/06/04/glmm.html">here</a>). We assume to have the following generalized linear mixed model:

$$
\mathbf{y}_{i,j} \rvert \beta_i \sim \text{NegBin}(\mu_{i,j}, \phi);
\hspace{10mm} \mu_{i,j} = \exp\left(\eta_{i,j}\right) = \exp\left(\alpha_i + \beta_i \mathbf{z}_{i,j}\right)
\label{eq:glmm-y}
$$

We'll use a superscript $\star$ to denote a quantity evaluated at the parameter estimates made under $H_0$ (i.e. $\tau^2 = \mathbf{0}$). Our <i>working</i> responses and errors are:

$$
\mathbf{y}^\star_{i,j} = \alpha_i + \beta_i \mathbf{z}_{i,j} + \epsilon^\star_{i,j};
\hspace{10mm}
\epsilon^\star_{i,j} \sim \mathcal{N}\left(0, \frac{V(\hat{\mu}_{i,j})}{\delta^2(\hat{\eta}_{i,j})}\right)
$$

where $$\delta(\hat{\eta}_{i,j}) = \frac{\partial g^{-1}(\eta_{i,j})}{\partial \eta_{i,j}}\bigg\rvert_{\eta_{i,j} = \hat{\eta}_{i,j}}$$. We can then just apply all of the results we found in the previous section to this case but make $$\hat{\sigma}^2$$ different for each observation, where $$\hat{\sigma}^2_{i,j} = \text{Var}(\epsilon_{i,j})$$. 

To do so, we need an estimate of $\alpha_i$ under $H_0$. In this case, the model for a cluster $i$ reduces down to:

$$
\mathbf{y}^\star_{i} = \alpha_i \mathbf{1}_n + \epsilon^*_{i}
\hspace{10mm}
\mathbf{y}^\star_{i} = \begin{bmatrix} \mathbf{y}^*_{i, 1} \\ \vdots \\ \mathbf{y}^*_{i, n} \end{bmatrix};
\hspace{2mm}
\epsilon^*_{i} = \begin{bmatrix} \epsilon^*_{i, 1} \\ \vdots \\ \epsilon^*_{i,n} \end{bmatrix}
$$

Thus, a solution can be found in closed form via <a href="/posts/2025/06/03/glm.html#weighted-least-squares">weighted least squares</a> as:

$$
\hat{\alpha}_i = (\mathbf{1}_n^\top \mathbf{W}_i \mathbf{1}_n)^{-1} \mathbf{1}_n^\top \mathbf{W}_i \mathbf{y}_i;
\hspace{8mm}
\mathbf{W}_i =
\begin{bmatrix}
\frac{1}{\hat{\sigma}_{i,1}^2} & \dots & 0 \\
\vdots & \ddots & \vdots \\
0 & \dots & \frac{1}{\hat{\sigma}_{i,n}^2}
\end{bmatrix}
$$

<!-- 
### Derivatives of the Conditional Log-Likelihood
First, note that:

$$
\frac{\partial \mu_{i,j}}{\partial \alpha_i} = \exp(\mathbf{x}_{i,j}^\top \alpha_i  + \mathbf{z}_{i,j}^\top \beta_i ) \mathbf{x}^\top_{i,j} = \mu_{i,j} \mathbf{x}^\top_{i,j};
\hspace{8mm}
\frac{\partial \mu_{i,j}}{\partial \beta_i} = \exp(\mathbf{x}_{i,j}^\top \alpha_i  + \mathbf{z}_{i,j}^\top \beta_i ) \mathbf{z}^\top_{i,j} = \mu_{i,j} \mathbf{z}^\top_{i,j}
\label{eq:deriv-mu}
$$

We have:

$$
\frac{\partial \ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)}{\partial \beta_i} = \sum_{j = 1}^{n_i}  \frac{\mathbf{y}_{i,j}}{1 + \phi \mu_{i,j}}\mathbf{z}^\top_{i,j}; 
\hspace{8mm}
\frac{\partial \ell(\mathbf{y}_i; \alpha_i, \tau^2, \phi \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top}  = -\sum_{j = 1}^{n_i} \frac{\mathbf{y}_{i,j}\phi \mu_{i,j} }{(1  + \phi \mu_{i,j})^2} \mathbf{z}_{i,j}\mathbf{z}_{i,j}^\top
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{\partial \ell(\mathbf{y}_i; \alpha_i, \tau^2\rvert \beta_i)}{\partial \beta_i}
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \left( \frac{\partial}{\partial \beta_i} \left[ \log(\phi \mu_{i,j}) - \log(1 + \phi \mu_{i,j}) \right] \right)  \\
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \left( \frac{\phi\mu_{i,j} }{\phi \mu_{i,j}} - \frac{ \phi \mu_{i,j}}{1 + \phi \mu_{i,j}} \right)\mathbf{z}_{i,j}^\top \\
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \left( 1  - \frac{\phi \mu_{i,j}}{1 + \phi \mu_{i,j}}  \right)  \mathbf{z}_{i,j}^\top \\
&= \sum_{j = 1}^{n_i}  \frac{\mathbf{y}_{i,j}}{1 + \phi \mu_{i,j}}\mathbf{z}^\top_{i,j} \\
\frac{\partial \ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top}
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \frac{\partial}{\partial \beta_i^\top} \left[ \frac{1}{1 + \phi \mu_{i,j}}  \right]  \mathbf{z}^\top_{i,j} \\
&= -\sum_{j = 1}^{n_i} \frac{\mathbf{y}_{i,j}\phi \mu_{i,j} }{(1 + \phi \mu_{i,j})^2} \mathbf{z}_{i,j}\mathbf{z}_{i,j}^\top
\end{aligned}
\nonumber
$$
</details>

$$
\frac{d \ell(\mathbf{y}_i; \alpha_i, \tau^2  \rvert \beta_i)}{d \alpha_i} = \sum_{j = 1}^{n_i}  \frac{\mathbf{y}_{i,j} }{1 + \phi \mu_{i,j}}\mathbf{x}^\top_{i,j} 
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d \ell(\mathbf{y}; \alpha_i, \tau^2 \rvert \beta_i)}{d \alpha_i} 
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \left( \frac{d}{d\alpha_i} \left[ \log(\phi \mu_{i,j}) - \log(1 + \phi \mu_{i,j}) \right] \right) \\
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \left( \frac{\phi \mu_{i,j}}{\phi \mu_{i,j}} - \frac{\phi \mu_{i,j}}{1 + \phi \mu_{i,j}} \right)  \mathbf{x}_{i,j}^\top \\ 
&= \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \left( 1 - \frac{\phi \mu_{i,j}}{1 + \phi \mu_{i,j}} \right)\mathbf{x}^\top_{i,j}   \\
&= \sum_{j = 1}^{n_i} \frac{\mathbf{y}_{i,j}  }{1 + \phi \mu_{i,j}}\mathbf{x}^\top_{i,j}
\end{aligned}
\nonumber
$$
</details>


### Approximating the Marginal Log-Likelihood
We'll assume $\tau^2$ is a scalar, so $D(\tau^2) = \tau^2 \mathbb{I}_{q \times q}$. Let $$\zeta_{i,j} = \phi \mu_{i,j}$$ and $$\xi_{i,j} = \frac{1}{1 + \zeta_{i,j}}$$. We can then define:

$$
\begin{aligned}
\mathbf{K}_i &= D^{-1}(\tau^2) - \left[ \frac{\partial \ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i )}{\partial \beta_i \partial \beta_i^\top} \right]\bigg\rvert_{H_0} = \frac{1}{\tau^2}\mathbb{I}_{q \times q} + \sum_{j = 1}^{n_i}\frac{\mathbf{y}_{i,j}  \phi \mu_{i,j}^0}{(1 + \phi \mu_{i,j}^0)^2}\mathbf{z}_{i,j} \mathbf{z}_{i,j}^\top = \frac{1}{\tau^2}\mathbb{I}_{q \times q} + \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{y}_{i,j} \zeta_{i,j} \mathbf{z}_{i,j} \mathbf{z}_{i,j}^\top \\
\mathbf{h}_i &= \left[ \frac{\partial \ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)}{\partial \beta_i} \right] \bigg\rvert_{H_0} = \sum_{j = 1}^{n_i}  \frac{\mathbf{y}_{i,j} }{1 + \phi \mu^0_{i,j}}\mathbf{z}_{i,j}^\top = \sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{z}^\top_{i,j}
\end{aligned}
$$

Using marginal quasi-likelihood (see <a href="/posts/2025/06/04/glmm.html">here</a>), we approximate the marginal log-likelihood as:

$$
\begin{aligned}
\ell(\mathbf{y}; \alpha, \tau^2)
&\approx \sum_{i = 1}^k \left[ - \frac{1}{2} \log\left( \rvert \tau^2 \mathbf{K}_i  \rvert \right) + [\ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)]\rvert_{H_0} + \frac{1}{2} \mathbf{h}_i^\top \mathbf{K}_i^{-1} \mathbf{h}_i \right] 
\end{aligned}
$$

We can explicitly compute $\mathbf{K}_i^{-1}$ with a little bit of work. Notice that the second term is just the sum of (scaled) rank one matrices, which itself is just a matrix. Let's define:

$$
\mathbf{G}_{i} = \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{y}_{i,j} \zeta_{i,j} \mathbf{z}_{i,j} \mathbf{z}_{i,j}^\top
\implies 
\mathbf{K}_i = \frac{1}{\tau^2} \mathbb{I}_{q \times q} + \mathbf{G}_i 
$$

Using this, we have:

$$
\mathbf{K}_i^{-1} = (\tau^2)^2 \mathbf{K}_i
$$

<details>
<summary>Proof.</summary>
We can then use the <a href="https://en.wikipedia.org/wiki/Woodbury_matrix_identity#Pseudoinverse_with_positive_semidefinite_matrices">Woodbury identity</a> to invert $\mathbf{K}_i$. For two square matrices with the same dimensions, the inverse of their sum is:
$$
(\mathbf{A} + \mathbf{B})^{-1} = \mathbf{A}^{-1} - (\mathbf{A} + \mathbf{A} \mathbf{B}^{-1} \mathbf{A})^{-1}
$$
Then:
$$
\begin{aligned}
\mathbf{K}_i^{-1}
&= \left(\frac{1}{\tau^2} \mathbb{I}_{q \times q} + \mathbf{G}_i \right)^{-1} \\
&= \tau^2 \mathbb{I}_{q \times q} - \left[\frac{1}{(\tau^2)^2} \mathbb{I}_{q \times q} \mathbf{G}_i^{-1} \mathbb{I}_{q \times q}\right]^{-1} \\
&= \tau^2 \mathbb{I}_{q \times q} - \left[ \frac{1}{(\tau^2)^2}\mathbb{I}_{q \times q} \right]^{-1} [\mathbf{G}^{-1}_i]^{-1} \\
&=  \tau^2 \mathbb{I}_{q \times q} - (\tau^2)^2 \mathbb{I}_{q \times q} \mathbf{G} \\
&= \tau^2 \left( \mathbb{I}_{q \times q} - \tau^2 \mathbf{G} \right) \\
&= (\tau^2)^2 \mathbf{K}_i
\end{aligned}
$$
</details>


### Score
We'll derive the score component-by-component. First, let's find the derivatives of $$\zeta_{i,j} = \phi \mu_{i,j}$$ and $$\xi_{i,j} = \frac{1}{1 + \phi \mu_{i,j}}$$ w.r.t the parameters:

$$
\begin{aligned}
\frac{d \zeta_{i,j}}{d \alpha_i} &= \zeta_{i,j} \mathbf{x}_{i,j} 
&
\frac{d \xi_{i,j}}{d \alpha_i} &= -\xi_{i,j}^2\zeta_{i,j} \mathbf{x}_{i,j}
\end{aligned}
$$

<details>
<summary>Proof.</summary>
First, the derivatives w.r.t. $\alpha_i$:
$$
\begin{aligned}
\frac{d \zeta_{i,j}}{d \alpha_i}
&= \frac{d}{\alpha_i} \left[ \phi \mu_{i,j} \right] \\
&= \phi \exp(\alpha_i \mathbf{x}_{i,j} + \beta_i \mathbf{z}_{i,j}) \mathbf{x}_{i,j} \\
&= \phi \mu_{i,j} \mathbf{x}_{i,j} \\
&= \zeta_{i,j} \mathbf{x}_{i,j} \\
\frac{d \xi_{i,j}}{d \alpha_i}
&= \frac{d}{d \alpha_i} \left[ \frac{1}{1 + \phi \mu_{i,j}}  \right] \\
&= -\frac{\phi \mu_{i,j} \mathbf{x}_{i,j}}{(1 + \phi \mu_{i,j})^2} \\
&= -\xi_{i,j}^2\zeta_{i,j} \mathbf{x}_{i,j}
\end{aligned}
\nonumber
$$
</details>

With the above, we can derive:

$$
\begin{aligned}
\frac{d \mathbf{h}_i}{d \alpha_i} &= - \sum_{j = 1}^{n_i} (\xi^0_{i,j})^2  \mathbf{y}_{i,j} \mathbf{z}_{i,j} \mathbf{x}_{i,j} \zeta^0_{i,j}
& 
\frac{d \mathbf{K}_i}{d \alpha_i} &= \sum_{j = 1}^{n_i} (\xi^0_{i,j})^2 \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j} \mathbf{x}_{i,j} \zeta^0_{i,j} \left( 1 - 2 \xi^0_{i,j} \zeta^0_{i,j} \right) \\
\frac{d \mathbf{h}_i}{d \tau^2} &= \mathbf{0}
&
\frac{d \mathbf{K}_i}{d \tau^2} &= -\frac{1}{(\tau^2)^2}\mathbb{I}_{q \times q}
\end{aligned}
$$

<details>
<summary>Proof.</summary>
In the following, I've dropped the superscript $0$ to denote evaluation under $H_0$:
$$
\begin{aligned}
\frac{d \mathbf{h}_i}{d \alpha_i} &= 
\frac{d}{d \alpha_i} \left[ \sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{z}_{i,j} \right] \\
&= - \sum_{j = 1}^{n_i} \xi_{i,j}^2 \zeta_{i,j} \mathbf{y}_{i,j} \mathbf{z}_{i,j} \mathbf{x}_{i,j} \\
\frac{d \mathbf{K}_i}{d \alpha_i} 
&= \frac{d}{d \alpha_i} \left[ \frac{1}{\tau^2} + \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j} \zeta_{i,j} \right] \\
&= \sum_{j = 1}^{n_i} \left[ \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j} \left(\xi_{i,j}^2 \frac{d \zeta_{i,j}}{d \alpha_i} + \zeta_{i,j} \frac{d \xi^2_{i,j}}{d \alpha_i} \right)\right] \\
&= \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j}\left(\xi_{i,j}^2 \zeta_{i,j} \mathbf{x}_{i,j} - 2 \xi^3_{i,j} \zeta^2_{i,j} \mathbf{x}_{i,j} \right) \\
&= \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j} \mathbf{x}_{i,j} \zeta_{i,j} \left( 1 - 2 \xi_{i,j} \zeta_{i,j} \right) 
\end{aligned}
\nonumber
$$
</details>


We then have:


We omit the superscript $0$, but all quantities are evaluated under $H_0$. 
$$
\begin{aligned}
\frac{d \ell(\mathbf{y}; \alpha, \tau^2, \phi)}{d \alpha_i}
&= \frac{d}{d \alpha_i} \left[  \sum_{i = 1}^k \left[ - \frac{1}{2} \log\left( \tau^2 \mathbf{K}_i \right) + [\ell(\mathbf{y}_i; \alpha_i, \tau^2, \phi \rvert \beta_i)]\rvert_{H_0} + \frac{1}{2} \mathbf{h}_i^2 \mathbf{K}_i^{-1}\right]  \right]  \\
&= - \frac{1}{2} \frac{d}{d \alpha_i} \left[ \log\left( \tau^2 \mathbf{K}_i \right) \right] + \left[ \frac{d \ell(\mathbf{y}_i; \alpha_i, \tau^2, \phi \rvert \beta_i)}{d \alpha_i} \right] \bigg\rvert_{H_0} + \frac{1}{2} \frac{d}{d \alpha_i} \left[ \mathbf{h}_i^2 \mathbf{K}_i^{-1} \right] \\
&= -\frac{\tau^2 \frac{d \mathbf{K}_i}{d \alpha_i}}{2 \tau^2 \mathbf{K}_i}  + \sum_{j = 1}^{n_i}  \frac{\mathbf{y}_{i,j} \mathbf{x}_{i,j} }{1 + \phi \mu_{i,j}} + \frac{1}{2} \left(2 \mathbf{h}_i \mathbf{K}_i^{-1} \frac{d \mathbf{h}_i}{d \alpha_i} - \mathbf{h}_i^2 \mathbf{K}_i^{-2} \frac{d \mathbf{K}_i}{d \alpha_i} \right) \\
&= -\frac{\frac{d \mathbf{K}_i}{d \alpha_i}}{2\mathbf{K}_i} \left(1 + \frac{\mathbf{h}_i^2}{\mathbf{K}_i}\right) + \frac{\mathbf{h}_i \frac{d \mathbf{h}_i}{d \alpha_i}}{\mathbf{K}_i} +  \sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{x}_{i,j} \\
&= -\frac{\sum_{j = 1}^{n_i}\xi_{i,j}^2 \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j} \mathbf{x}_{i,j} \zeta_{i,j}(1 - 2 \xi_{i,j} \zeta_{i,j})}{2\left(\frac{1}{\tau^2} + \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j} \mathbf{y}_{i,j} \zeta_{i,j}\right)} \left(1 + \frac{\left(\sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{z}_{i,j} \right)^2}{\frac{1}{\tau^2} + \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j} \mathbf{y}_{i,j} \zeta_{i,j}}\right) - \frac{\left(\sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{z}_{i,j} \right) \left( \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{y}_{i,j} \mathbf{z}_{i,j} \mathbf{x}_{i,j} \zeta_{i,j} \right) }{\frac{1}{\tau^2} + \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j} \mathbf{y}_{i,j} \zeta_{i,j}} + \sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{x}_{i,j} \\
&= -\frac{\tau^2\sum_{j = 1}^{n_i}\xi_{i,j}^2 \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j} \mathbf{x}_{i,j} \zeta_{i,j}(1 - 2 \xi_{i,j} \zeta_{i,j})}{2\left(1 + \tau^2 \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j} \mathbf{y}_{i,j} \zeta_{i,j}\right)} \left(1 + \frac{\tau^2\left(\sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{z}_{i,j} \right)^2}{1 + \tau^2\sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j} \mathbf{y}_{i,j} \zeta_{i,j}}\right) - \frac{\tau^2 \left(\sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{z}_{i,j} \right) \left( \sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{y}_{i,j} \mathbf{z}_{i,j} \mathbf{x}_{i,j} \zeta_{i,j} \right) }{1 + \tau^2\sum_{j = 1}^{n_i} \xi_{i,j}^2 \mathbf{z}_{i,j} \mathbf{y}_{i,j} \zeta_{i,j}} + \sum_{j = 1}^{n_i} \xi_{i,j} \mathbf{y}_{i,j} \mathbf{x}_{i,j}
\end{aligned}
\nonumber
$$ -->
<!-- 

---

## Poisson Case
For this example, we'll assume the responses are Poisson and that we only have a single fixed effect and a single random effect. We can write the log-likelihood for cluster $i$ in vector notation as:

$$
\ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)
= \sum_{j = 1}^{n_i} \left[ \mathbf{y}_{i,j} \log(\mu_{i,j}) - \mu_{i,j} - \log(\mathbf{y}_{i,j}!) \right] 
= \mathbf{y}_i^\top \log(\mu_i) - \mu_i^\top \mathbf{1}_{n_i} - \sum_{j = 1}^{n_i} \log( \mathbf{y}_{i,j}!)
$$

Using marginal quasi-likelihood (see <a href="/posts/2025/06/04/glmm.html">here</a>), we approximate the marginal log-likelihood as:

$$
\begin{aligned}
  \ell(\mathbf{y}; \alpha, \tau^2) 
  &\approx \sum_{i = 1}^k \left[-\frac{1}{2}\log\left(1 + \tau^2 \exp(\alpha_i) \sum_{j= 1}^{n_i} \mathbf{z}_{i,j}^2\right) + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j}\alpha_i - \exp(\alpha_i) - \log(\mathbf{y}_{i,j}!)) + \frac{1}{2} \left( \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right)^2 \left(\frac{1}{\tau^2} +  \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2 \right)^{-1}\right]
\end{aligned}
$$  

<details>
<summary>Proof.</summary>
We have:
$$
\frac{d \ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)}{d \beta_i} = \sum_{j = 1}^{n_i} \left[ \mathbf{y}_{i,j} \mathbf{z}_{i,j} - \exp(\alpha_i \mathbf{x}_{i,j} + \beta_i \mathbf{z}_{i,j})\mathbf{z}_{i,j} \right]
\hspace{5mm}
\frac{d^2 \ell(\mathbf{y}_i; \alpha_i, \tau^2 \rvert \beta_i)}{d \beta_i^2} = - \sum_{j = 1}^{n_i} \exp(\alpha_i \mathbf{x}_{i,j} + \beta_i \mathbf{z}_{i,j}) \mathbf{z}_{i,j}^2
$$
In our case, we only have a random slope with $D(\tau^2) = \tau^2$, so $\mathbf{K}_i$ and $\mathbf{h}_i$ are fairly simple:
$$
\mathbf{K}_i = \frac{1}{\tau^2} - \left(- \sum_{j = 1}^{n_i} \mu_{i,j}^0 \mathbf{z}_{i,j}^2\right);
\hspace{5mm}
\mathbf{h}_i = \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}
$$
We can then just plug this into the formula from my <a href="/posts/2025/06/04/glmm.html">other post</a>.
</details>


---

### Score Functions

Let $$\dot{\mathbf{z}}_i = \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}(\mathbf{y}_{i,j} - \exp(\alpha_i))$$ and define:

$$
\begin{aligned}
\upsilon_i &= \frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)} 
&
\zeta_i &= \frac{\exp(\alpha_i)}{\frac{1}{\tau^2} + \exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}} 
= \frac{\tau^2 \exp(\alpha_i)}{1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}} \\
\xi_i &= \frac{\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)} = 
\frac{\tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}
&
\psi_i &= \frac{1}{1 + 2 \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \left(\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2}
\end{aligned}
$$


We have:

$$
\begin{aligned}
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \bigg\rvert_{H_0} 
&\approx -\upsilon_i + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))
- \zeta_i \dot{\mathbf{z}}_i \left[\xi_i \dot{\mathbf{z}}_i + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right] 
\end{aligned}
$$

<details>
<summary>Proof.</summary>
  We have:
  $$
  \begin{aligned}
  \frac{d}{d \alpha_i} \left[ - \frac{1}{2}\log\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)\right] \bigg\rvert_{H_0} 
  &= -\frac{1}{2}\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{-1} \left(\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \\
  \frac{d}{d \alpha_i} \left[ \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} \alpha_i - \exp(\alpha_i) - \log(\mathbf{y}_{i,j}!)) \right] 
  &= \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \\
  \frac{d}{d \alpha_i} \left[ \frac{1}{2} \left( \sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j} \right)^2 \right]
  &= -\left(\sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right) \left(\exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j} \right)\\
  \frac{d}{d \alpha_i} \left[ \left( \frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{-1} \right] 
  &= - \left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{-2} \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)
  \nonumber
  \end{aligned}
  $$
  Thus, the derivative is:
  $$
  \begin{aligned}
  \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \bigg\rvert_{H_0} 
  &\approx \frac{d}{\alpha_i} \left[ -\frac{1}{2}\log\left(1 + \tau^2 \exp(\alpha_i) \sum_{j= 1}^{n_i} \mathbf{z}_{i,j}^2\right) + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j}\alpha_i - \exp(\alpha_i) - \log(\mathbf{y}_{i,j}!)) + \frac{1}{2} \left( \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right)^2 \left(\frac{1}{\tau^2} +  \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2 \right)^{-1} \right] \\
  &= -\frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}
  + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))
  + \frac{-\left(\sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right) \left(\exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j} \right)}{\frac{1}{\tau^2} + \exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2} 
  + \frac{\frac{1}{2} \left( \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right)^2\left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}{-\left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{2}} \\
  &= -\frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}
  + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))
  - \frac{\left(\sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right) \left(\exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j} \right)}{\frac{1}{\tau^2} + \exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2} 
  - \frac{ \left( \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right)^2\left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}{2\left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{2}}  \\
  &= -\frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}
  + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))
  -\frac{\exp(\alpha_i)\sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}}{\frac{1}{\tau^2} + \exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}} \left[\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}+ \frac{\left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j} \right)\left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)}{2\left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)}\right] \\
  &= -\upsilon_i + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))
  - \zeta_i \left[ \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}(\mathbf{y}_{i,j} - \exp(\alpha_i)) \right]\left[\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} + \xi_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} (\mathbf{y}_{i,j} - \exp(\alpha_i))\right] 
  \end{aligned}
  \nonumber
  $$
</details>

$$
\begin{aligned}
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2} 
&\approx \frac{1}{2}\sum_{i = 1}^k  \left[-\frac{2 \upsilon_i}{\tau^2} +\dot{\mathbf{z}}^2_i \psi_i \right]
\end{aligned}
$$

<details>
<summary>Proof.</summary>
  We have:
  $$
  \begin{aligned}
  \frac{d}{d \tau^2} \left[ - \frac{1}{2}\log\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)\right]
  &= -\frac{1}{2}\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{-1} \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \\
  \frac{d}{d \tau^2} \left[ \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} \alpha_i - \exp(\alpha_i) - \log(\mathbf{y}_{i,j}!)) \right] 
  &= 0\\
  \frac{d}{d \tau^2} \left[ \frac{1}{2} \left( \sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j} \right)^2 \right]
  &= 0 \\
  \frac{d}{d \tau^2} \left[ \left( \frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{-1} \right] 
  &= \left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^{-2} (\tau^2)^{-2}
  \end{aligned}
  \nonumber
  $$
  Thus, the derivative is:
  $$
  \begin{aligned}
  \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2} 
  &\approx \sum_{i = 1}^k \frac{d}{d\tau^2} \left[ -\frac{1}{2}\log\left(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right) + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j}\alpha_i - \exp(\alpha_i) - \log(\mathbf{y}_{i,j}!)) + \frac{1}{2} \left( \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j}\right)^2 \left(\frac{1}{\tau^2} +  \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2 \right)^{-1} \right] \\
  &= \sum_{i = 1}^k \left[-\frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)} + \frac{\left( \sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j} \right)^2}{2 (\tau^2)^2 \left(\frac{1}{\tau^2} + \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2}\right] \\
  &= \sum_{i = 1}^k  \left[-\frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2\left( 1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)} + \frac{\left( \sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j} \right)^2}{2 (\tau^2)^2 \left(\frac{1}{(\tau^2)^2} + 2\frac{1}{\tau^2} \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2 \right) }\right] \\
  &= \frac{1}{2}\sum_{i = 1}^k  \left[-\frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{1 + \tau^2 \exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 } + \frac{\left( \sum_{j = 1}^{n_i}(\mathbf{y}_{i,j} - \exp(\alpha_i))\mathbf{z}_{i,j} \right)^2}{1 + 2\tau^2\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \left(\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2}\right] \\
  &= \frac{1}{2}\sum_{i = 1}^k  \left[-\frac{2 \upsilon_i}{\tau^2} + \psi_i \left( \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}(\mathbf{y}_{i,j} - \exp(\alpha_i)) \right)^2 \right]
  \end{aligned}
  $$
</details>

Notice that under the null hypothesis, the gradient of the log-likelihood with respect to $\alpha$ becomes:

$$
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha} \approx 
\begin{bmatrix}
\sum_{j = 1}^{n_1} (\mathbf{y}_{1, j} - \exp(\alpha_1)) \\
\vdots \\
\sum_{j = 1}^{n_k} (\mathbf{y}_{k, j} - \exp(\alpha_k))
\end{bmatrix}
$$

Setting the above equal to $\mathbf{0}_{k}$ and solving for $\alpha$ gives us the MLE under the null:

$$
\hat{\alpha} = 
\begin{bmatrix}
\log\left(\frac{1}{n_1} \sum_{j = 1}^{n_1} \mathbf{y}_{1, j}\right) \\
\vdots \\
\log\left(\frac{1}{n_k} \sum_{j = 1}^{n_k} \mathbf{y}_{k, j}\right)
\end{bmatrix}
$$


#### A Note On The MLE
Suppose we let $\alpha^* = \alpha + \epsilon$ where $\alpha$ is the <i>true</i> value of $\alpha$ and $\epsilon$ is a little bit of jitter (negative or positive). Consider evaluating the gradient of the marginal log-likelihood w.r.t. $\alpha$ at $(\alpha, \tau^2) = (\alpha^*, 0)$. We know that $\upsilon_i = \zeta_i = \xi_i = 0$ for $\tau^2 = 0$, so the component of the score corresponding to $\alpha_i$ simplifies to:

$$
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \bigg\rvert_{(\alpha, \tau^2) = (\alpha^*, 0)}
= \left(\sum_{j = 1}^{n_i} \mathbf{y}_{i,j}\right) - n_i \exp(\epsilon)\exp(\alpha_i)
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \bigg\rvert_{(\alpha, \tau^2) = (\alpha^*, 0)}
&= \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i + \epsilon)) \\
&= \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i) \exp(\epsilon)) \\
&= \left(\sum_{j = 1}^{n_i} \mathbf{y}_{i,j}\right) - n_i \exp(\epsilon)\exp(\alpha_i)
\end{aligned}
\nonumber
$$
</details>

Since $e^x > 0$ and $n_i > 0$, we see that adding any jitter to $\alpha$ corresponds to scaling the amount we are subtracting off the cluster sum by a factor that is exponential in the amount of jitter. 

Now let's consider the gradient of the marginal log-likelihood w.r.t. $\tau^2$:

$$
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2} \bigg\rvert_{(\alpha, \tau^2) = (\alpha^*, 0)} =
\frac{1}{2}\sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}\right) \left( \exp(\epsilon)\exp(\alpha_i) - \frac{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'} + \frac{1}{2}\mathbf{z}_{i,j}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} \right)^2   + c
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2} \bigg\rvert_{(\alpha, \tau^2) = (\alpha^*, 0)}
&= \frac{1}{2} \sum_{i = 1}^k \left[ -\exp(\alpha_i + \epsilon) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 +  \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} (\mathbf{y}_{i,j} - \exp(\alpha_i + \epsilon))\right)^2 \right] \\
&= \frac{1}{2} \sum_{i = 1}^k \left[ - \exp(\alpha_i + \epsilon) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 (\mathbf{y}_{i,j} - \exp(\alpha_i + \epsilon))^2 + \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j} \mathbf{z}_{i,j'} (\mathbf{y}_{i,j} - \exp(\alpha_i + \epsilon))(\mathbf{y}_{i,j'} - \exp(\alpha_i + \epsilon)) \right] \\
&= \frac{1}{2} \sum_{i = 1}^k \left[ - \exp(\alpha_i + \epsilon) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{y}_{i,j}^2 - 2\exp(\alpha_i + \epsilon)\sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2 \mathbf{y}_{i,j}  + \exp^2(\alpha_i + \epsilon) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j} \mathbf{z}_{i,j'} (\mathbf{y}_{i,j}\mathbf{y}_{i,j'} - \exp(\alpha_i + \epsilon) \mathbf{y}_{i,j} - \exp(\alpha_i + \epsilon)\mathbf{y}_{i,j'} - \exp^2(\alpha_i + \epsilon)) \right] \\
&= \frac{1}{2}\sum_{i = 1}^k \left[ - \exp(\alpha_i + \epsilon) \sum_{j =1 }^{n_i} \mathbf{z}_{i,j}^2 + \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \mathbf{y}_{i,j}\right)^2 - 2\exp(\alpha_i + \epsilon)\left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \mathbf{y}_{i,j}\right)\left(\sum_{j =1 }^{n_i} \mathbf{z}_{i,j} \right) + \exp^2(\alpha_i + \epsilon) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right)^2 \right] \\
&= \frac{1}{2}\sum_{i = 1}^k \sum_{j = 1}^{n_i} \left[ - \exp(\alpha_i + \epsilon) \mathbf{z}^2_{i,j} + \mathbf{z}_{i,j} \mathbf{y}_{i,j}\left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'}\right) - 2\exp(\alpha_i + \epsilon)\mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'}\right)+ \exp^2(\alpha_i + \epsilon) \mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}\right) \right] \\
&= \frac{1}{2}\sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}\right) \left[ - \exp(\alpha_i + \epsilon) \frac{\mathbf{z}_{i,j}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} + \frac{\mathbf{y}_{i,j}\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} - 2\exp(\alpha_i + \epsilon)\frac{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}}+ \exp^2(\alpha_i + \epsilon) \right] \\
&= \frac{1}{2}\sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}\right) \left[ \exp^2(\alpha_i + \epsilon) - 2 \exp(\alpha_i + \epsilon) \left( \frac{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'} + \frac{1}{2}\mathbf{z}_{i,j}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}}  \right) + \frac{\mathbf{y}_{i,j}\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} \right] \\
&= \frac{1}{2}\sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}\right) \left[ \left( \exp(\alpha_i + \epsilon) - \frac{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'} + \frac{1}{2}\mathbf{z}_{i,j}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} \right)^2 - \left(\frac{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'} + \frac{1}{2}\mathbf{z}_{i,j}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}}\right)^2  + \frac{\mathbf{y}_{i,j}\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} \right] \\
&= \frac{1}{2}\sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \left(\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}\right) \left( \exp(\alpha_i + \epsilon) - \frac{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'} \mathbf{y}_{i,j'} + \frac{1}{2}\mathbf{z}_{i,j}}{\sum_{j' = 1}^{n_i} \mathbf{z}_{i,j'}} \right)^2  + c
\end{aligned}
\nonumber
$$
were $c$ is defined appropriately.
</details>

where $c$ is a term that does not depend upon $\epsilon$. For similar reasons as before, we see that all the jitter does is scale the leading term in the quadratic within the sum. Thus, as $\epsilon$ grows in magnitude, the score will also grow in magnitude exponentially in $2\epsilon$. 


---

### Information
We'll derive the information matrix in chunks. First, we have:

$$
\mathbb{E}_{H_0}\left[\dot{\mathbf{z}}_i \right] = 0;
\hspace{5mm}
\mathbb{E}_{H_0}\left[\dot{\mathbf{z}}_i^2 \right] = \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0;
\hspace{5mm}
\mathbb{E}_{H_0}\left[\dot{\mathbf{z}}_i^3 \right] = \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \mu_{i,j}^0;
\hspace{5mm}
\mathbb{E}_{H_0}\left[\dot{\mathbf{z}}_i^4 \right] = 3\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4 ((\mu_{i,j}^0)^2 + \mu_{i,j}^0) + 3\sum_{j = 1}^{n_i} \sum_{j' \neq j}  \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \mu_{i,j}^0 \mu_{i,j'}^0
$$

<details>
<summary>Proof.</summary>
  $$
  \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i\right] = \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right] = 0
  \nonumber
  $$
  $$
  \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right] = \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^2 \right]}_{=\text{Var}(\mathbf{y}_{i,j})} + \sum_{j = 1}^{n_i} \sum_{j' \neq j}  \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right]}_{=0}\underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'} - \exp(\alpha_i) \right]}_{=0} = \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0
  \nonumber
  $$
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^3 \right]
  &= \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^3 \right]}_{= \text{ 3rd central moment } \\ = \mu_{i,j}^0} + 2 \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'} \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^2 \right]}_{= \text{Var}(\mathbf{y}_{i,j})} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'} - \exp(\alpha_i) \right]}_{=0} \\
  &\hspace{5mm} + \sum_{j = 1}^{n_i} \sum_{j' \neq j} \sum_{j'' \neq j, j'} \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \mathbf{z}_{i,j''} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right]}_{=0} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'} - \exp(\alpha_i) \right]}_{=0}\underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j''} - \exp(\alpha_i) \right]}_{=0} \\
  &= \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \mu_{i,j}^0
  \end{aligned}
  \nonumber
  $$
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^4 \right] 
  &= \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4 \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^4 \right]}_{= \text{ 4th central moment } \\ = 3(\mu_{i,j}^0)^2 + \mu_{i,j}^0} + 4 \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j}^3 \mathbf{z}_{i,j'}\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i)^3 \right] \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'} - \exp(\alpha_i) \right]}_{= 0} \\
  &\hspace{5mm} + 3 \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^2 \right]}_{= \text{Var}(\mathbf{y}_{i,j})}\underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^2 \right]}_{= \text{Var}(\mathbf{y}_{i,j'})}
  + 12 \sum_{j = 1}^{n_i} \sum_{j' \neq j} \sum_{j'' \neq j, j'} \mathbf{z}^2_{i,j} \mathbf{z}_{i,j'} \mathbf{z}_{i,j''} \mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^2 \right] \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'} - \exp(\alpha_i)  \right]}_{= 0}\underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j''} - \exp(\alpha_i)  \right]}_{= 0} \\
  &\hspace{5mm} + \sum_{j = 1}^{n_i} \sum_{j' \neq j} \sum_{j'' \neq j, j'} \sum_{j''' \neq j, j', j''} \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \mathbf{z}_{i,j''} \mathbf{z}_{i,j'''} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i)  \right]}_{= 0} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'} - \exp(\alpha_i)  \right]}_{= 0} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j''} - \exp(\alpha_i)  \right]}_{= 0} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j'''} - \exp(\alpha_i)  \right]}_{= 0} \\
  &= \sum_{j = 1}^{n_i} 3 \mathbf{z}_{i,j}^4 ((\mu_{i,j}^0)^2 + \mu_{i,j}^0) + 3\sum_{j = 1}^{n_i} \sum_{j' \neq j}  \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \mu_{i,j}^0 \mu_{i,j'}^0 
  \end{aligned}
  \nonumber
  $$
</details>

We can show that:

$$
\mathbb{E}_{H_0}\left[ \left(\frac{\partial \ell(\mathbf{y}; \alpha, \tau^2)}{\partial \alpha_i}\right)^2 \right]
\approx n_i \exp(\alpha_i)
$$

<details>
<summary>Proof.</summary>
  All expectations are taking under the null hypothesis, so $\mathbb{E}_{H_0}[\mathbf{y}_{i,j}] = \exp(\alpha_i)$. Furthermore, we have that $\upsilon_i = 0$ and $\psi_i = 1$ under the null (also $\zeta_i$ and $\xi_i$ go to $0$ as we take $\tau^2 \rightarrow 0$). The cross terms come out to:
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \upsilon_i \zeta_i \dot{\mathbf{z}}_i \left[\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} + \xi_i \dot{\mathbf{z}}_i \right] \right] 
  &= \upsilon_i \zeta_i \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \underbrace{\mathbb{E}_{H_0}\left[\dot{\mathbf{z}}_i\right]}_{=0} + \xi_i \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right] \right) 
  = \upsilon_i \zeta_i \xi_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0 
  = \upsilon_i \zeta_i \xi_i  \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \\
  \mathbb{E}_{H_0}\left[ \upsilon_i \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right]
  &= \upsilon_i \sum_{j = 1}^{n_i} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right]}_{=0}
  = 0 \\
  \mathbb{E}_{H_0}\left[ \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right)\left( \zeta_i \dot{\mathbf{z}}_i \left[ \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} + \xi_i \dot{\mathbf{z}}_i \right]\right) \right]
  &= \zeta_i \left(\sum_{j = 1}^{n_i}\mathbf{z}_{i,j}\right) \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i))\right] + \zeta_i \xi_i \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}^2_i \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right] \\
  &= \zeta_i \left(\sum_{j = 1}^{n_i}\mathbf{z}_{i,j}\right) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \mu_{i,j}^0 \right) + \zeta_i \xi_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0 \\
  &= \zeta_i \exp(\alpha_i) \left[ \left(\sum_{j = 1}^{n_i}\mathbf{z}_{i,j}\right)^2 + \xi_i \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)\right]
  \end{aligned}
  \nonumber
  $$
  All of the above terms go to $0$ under the null (as we take $\tau^2 \rightarrow 0$), so their sum also goes to $0$. Looking at the squared terms:
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \upsilon_i^2 \right]
  &= \upsilon_i^2 \\
  \mathbb{E}_{H_0}\left[ \left( \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right)^2 \right]
  &= \sum_{j = 1}^{n_i} \mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))^2 \right] + \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i))(\mathbf{y}_{i, j'} - \exp(\alpha_i)) \right] \\
  &= \sum_{j = 1}^{n_i} \text{Var}(\mathbf{y}_{i,j}) + \sum_{j = 1}^{n_i} \sum_{j' \neq j} \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right]}_{=0} \underbrace{\mathbb{E}_{H_0}\left[(\mathbf{y}_{i, j'} - \exp(\alpha_i)) \right]}_{=0} \\
  &= \sum_{j = 1}^{n_i} \mu_{i,j}^0 \\
  &= n_i \exp(\alpha_i) \\
  \mathbb{E}_{H_0}\left[ \left( \zeta_i \dot{\mathbf{z}}_i \left[\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} + \xi_i \dot{\mathbf{z}}_i \right]\right)^2\right] 
  &= \zeta_i^2 \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + \xi_i \dot{\mathbf{z}}_i \right)^2 \right] \\
  &= \zeta_i^2 \left[\left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)^2 \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right]+ 2 \xi^2_i \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}^3_i \right] + \xi_i^2 \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^4 \right] \right] \\
  &= \zeta_i^2 \left[ \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2\left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0\right) + 2 \xi_i^2 \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \mu_{i,j}^0 \right) + 3 \xi_i^2 \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \mu_{i,j}^0 \mu_{i,j'}^0\right] \\
  &= \zeta_i^2 \exp(\alpha_i) \left[ \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^3 + 2 \xi_i^2 \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \right) + 3 \xi_i^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \sum_{j' \neq j} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \right] \\
  \end{aligned}
  \nonumber
  $$
  Only the second term is non-zero under the null hypothesis, so:
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \left(\frac{\partial \ell(\mathbf{y}; \alpha, \tau^2)}{\partial \alpha_i}\right)^2 \right]
  \approx n_i \exp(\alpha_i)
  \end{aligned}
  $$
</details>

$$
\begin{aligned}
\mathbb{E}_{H_0}\left[\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_l} \right] 
&\approx 
  \mathbb{E}_{H_0}\left[\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_l} \right] = 0 
\end{aligned}
$$

<details>
<summary>Proof.</summary>
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \upsilon_i \upsilon_l \right] &= \upsilon_i \upsilon_l \\
  \mathbb{E}_{H_0}\left[ \upsilon_i \sum_{j = 1}^{n_l} (\mathbf{y}_{l, j} - \exp(\alpha_l)) \right]
  &= \upsilon_i \sum_{j = 1}^{n_l} \underbrace{\mathbb{E}_{H_0} \left[ \mathbf{y}_{l,j} - \exp(\alpha_l) \right]}_{=0} = 0 \\
  \mathbb{E}_{H_0}\left[ \upsilon_i \zeta_l \dot{\mathbf{z}}_l \left(\xi_l \dot{\mathbf{z}}_l + \sum_{j = 1}^{n_l} \mathbf{z}_{l, j} \right) \right] 
  &= \upsilon_i \zeta_l \xi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \right] + \upsilon_l \zeta_l \left( \mathbf{z}_{l,j} \right) \underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l \right]}_{=0} \\
  &= \upsilon_i \zeta_l \xi_l \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \mu_{l, j}^0 \\
  &= \upsilon_i \zeta_l \xi_l \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l, j}^2 \\
  \mathbb{E}_{H_0}\left[ \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right) \left(\sum_{j = 1}^{n_l} (\mathbf{y}_{l,j} - \exp(\alpha_l)) \right)\right] 
  &= \sum_{j = 1}^{n_i} \sum_{j = 1}^{n_l} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right]}_{=0} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{l,j} - \exp(\alpha_l) \right]}_{=0} \\
  &= 0 \\
  \mathbb{E}_{H_0}\left[ \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right) \zeta_l \dot{\mathbf{z}}_l\left( \xi_l \dot{\mathbf{z}}_l + \sum_{j = 1}^{n_l} \mathbf{z}_{l,j} \right)\right] 
  &= \zeta_l \xi_l \sum_{j = 1}^{n_i} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i)\right]}_{=0} \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l\right] + \sum_{j = 1}^{n_i} \sum_{j = 1}^{n_l} \mathbf{z}_{l,j} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right]}_{=0} \\
  &= 0 \\
  \mathbb{E}_{H_0}\left[ \zeta_i \zeta_l \dot{\mathbf{z}}_i \dot{\mathbf{z}}_l \left(\xi_i \dot{\mathbf{z}}_i + \sum_{j =1}^{n_i} \mathbf{z}_{i,j}\right)\left(\xi_l \dot{\mathbf{z}}_l + \sum_{j = 1}^{n_l} \mathbf{z}_{l,j} \right)\right] 
  &= \zeta_i \zeta_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i \dot{\mathbf{z}}_l\left(\xi_i \xi_l \dot{\mathbf{z}}_i \dot{\mathbf{z}}_l + \xi_i \dot{\mathbf{z}}_i \sum_{j = 1}^{n_l} \mathbf{z}_{l,j} + \xi_l \dot{\mathbf{z}}_l \sum_{j = 1}^{n_l} \mathbf{z}_{i,j} + \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right)\left(  \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}\right)\right) \right] \\
  &= \zeta_i \zeta_l \left[ \xi_i \xi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right] \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \right] + \xi_i \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j} \right) \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right] \underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l \right]}_{=0} +  \xi_l \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \right] \underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i \right]}_{=0} + \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right)\left(  \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}\right)\underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i \right]}_{=0}\underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l \right]}_{=0}\right] \\
  &= \zeta_i \zeta_l \xi_i \xi_l \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0 \right) \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \mu_{l,j}^0 \right) \\
  &= \zeta_i \zeta_l \xi_i \xi_l \eta(\alpha_i) \eta(\alpha_l)  \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \right)
  \end{aligned}
  \nonumber
  $$
  As $\tau^2 \rightarrow 0$, all of the above terms go to zero, so we get:
  $$
  \mathbb{E}_{H_0}\left[\frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_l} \right] \approx 0 
  $$
  under the null hypothesis.
</details>

$$
\mathbb{E}_{H_0}\left[ \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right] \approx \frac{1}{2} \sum_{l = 1}^k \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2
$$

<details>
<summary>Proof.</summary>
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \upsilon_i \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right]
  &\approx \upsilon_i \mathbb{E}_{H_0}\left[ \frac{1}{2} \sum_{l = 1}^k \left[ - \frac{2\upsilon_l}{\tau^2} +\dot{\mathbf{z}}_l^2  \psi_l \right]\right] \\
  &= \frac{1}{2} \upsilon_i \sum_{l = 1}^{k} - \frac{2\upsilon_l}{\tau^2} + \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \right] \\
  &= \frac{1}{2} \upsilon_i \sum_{l = 1}^{k} \left( - \frac{2\upsilon_l}{\tau^2} + \psi_l \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \mu_{l,j}^0 \right) \\
  &= \frac{1}{2} \upsilon_i \sum_{l = 1}^{k} \left(- \frac{2 \upsilon_l}{\tau^2} + \psi_l \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right) \\
  &= \frac{1}{2} \upsilon_i \sum_{l = 1}^{k} \left(- \frac{\exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2}{1 + \tau^2 \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2}+ \psi_l \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right)\\
  &= \frac{1}{2} \upsilon_i \sum_{l = 1}^{k} \left(\exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right) \left(\psi_l - \frac{1}{1 + \tau^2 \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2} \right)\\
  \mathbb{E}_{H_0}\left[\left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right) \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2} \right]
  &= \frac{1}{2} \sum_{l = 1}^k \mathbb{E}_{H_0}\left[ \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right) \left(- \frac{2\upsilon_l}{\tau^2} + \dot{\mathbf{z}}_l^2  \psi_l \right) \right] \\
  &= \frac{1}{2}\sum_{l = 1}^k \left[ - \frac{2 \upsilon_l}{\tau^2} \sum_{j = 1}^{n_i} \underbrace{\mathbb{E}_{H_0}\left[ \mathbf{y}_{i,j} - \exp(\alpha_i) \right]}_{=0} +  \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \sum_{j = 1}^{n_i}  (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right] \right] \\
  &= \frac{1}{2} \sum_{l = 1}^k \left[ \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \mu_{l,j}^0 + \sum_{j = 1}^{n_i} \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \right] \underbrace{\mathbb{E}_{H_0}\left[ (\mathbf{y}_{i,j} - \exp(\alpha_i)) \right]}_{=0}\right] \\
  &= \frac{1}{2}  \sum_{l = 1}^k \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\\
  \mathbb{E}_{H_0}\left[ \zeta_i \dot{\mathbf{z}}_i \left[ \xi_i \dot{\mathbf{z}}_i + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right] \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right]
  &\approx  \frac{1}{2}\zeta_i \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i \left[ \xi_i \dot{\mathbf{z}}_i + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right] \left(\sum_{l = 1}^k  \left[-\frac{2 \upsilon_l}{\tau^2} +\dot{\mathbf{z}}^2_l \psi_l \right]\right) \right] \\
  &= \frac{1}{2}\zeta_i \sum_{l = 1}^k \mathbb{E}_{H_0}\left[ \left( \xi_i \dot{\mathbf{z}}_i^2 + \dot{\mathbf{z}}_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right)\left(-\frac{2 \upsilon_l}{\tau^2} +\dot{\mathbf{z}}^2_l \psi_l \right)\right] \\
  &= \frac{1}{2}\zeta_i \sum_{l = 1}^k \left[ -\frac{2 \upsilon_l \xi_i}{\tau^2} \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right] - \frac{2 \upsilon_l}{\tau^2} \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \underbrace{\mathbb{E}_{H_0} \left[ \dot{\mathbf{z}}_i \right]}_{=0} + \xi_i \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \dot{\mathbf{z}}_l \right] + \psi_l \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right)\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \dot{\mathbf{z}}_i \right] \right] \\
  &= \frac{1}{2}\zeta_i  \left[ - \frac{2\xi_i}{\tau^2} \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0 \right) \sum_{l = 1}^k \upsilon_l + \xi_i \sum_{l = 1 \\ l \neq i}^k \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right]  \underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l \right]}_{=0} + \xi_i \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^3 \right] + \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \left( \sum_{l = 1\\ l \neq i}^k \psi_l \mathbb{E}_{H_0}\left[  \mathbf{z}_{l}^2\right] \underbrace{\mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i \right]}_{=0} + \psi_i \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^3 \right] \right) \right] \\
  &= \frac{1}{2}\zeta_i  \left[ - \frac{2\xi_i \exp(\alpha_i)}{\tau^2} \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \sum_{l = 1}^k \upsilon_l + \xi_i \psi_i \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^3 \right] + \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \psi_i \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^3 \right] \right] \\
  &= \frac{1}{2}\zeta_i  \left[ - \frac{2\xi_i \exp(\alpha_i)}{\tau^2} \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \sum_{l = 1}^k \upsilon_l + \left( \xi_i \psi_i + \psi_i\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}  \right)\left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \mu_{i,j}^0 \right) \right] \\
  &= \frac{1}{2}\zeta_i \exp(\alpha_i) \left[ - \frac{2\xi_i}{\tau^2} \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \sum_{l = 1}^k \upsilon_l + \psi_i \left( \xi_i + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right)\left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^3 \right) \right] \\
  \end{aligned}
  \nonumber
  $$
  Since $\upsilon_i \rightarrow 0$ and $\zeta_i \rightarrow 0$ as $\tau^2 \rightarrow 0$, the first and last term above go to zero in the limit too. Thus:
  $$
  \mathbb{E}_{H_0}\left[ \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right] \approx \frac{1}{2} \sum_{l = 1}^k \exp(\alpha_l) \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2
  $$
  under the null hypothesis. 
</details>

$$
\mathbb{E}_{H_0}\left[ \left( \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right)^2 \right]
\approx 
\frac{1}{4} \sum_{i =1}^k \sum_{l = 1 \\ l \neq i}^k \exp(\alpha_i) \exp(\alpha_l) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2  \right) \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right)
+ \frac{1}{4}\sum_{i = 1}^{k}\left(3 (\exp^2(\alpha_i) + \exp(\alpha_i)) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4  + \exp^2(\alpha_i)\sum_{j = 1}^{n_i} \sum_{j' = 1 \\ j' \neq j}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \right)
$$

<details>
<summary>Proof.</summary>
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \left( \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right)^2 \right]
  &= \frac{1}{4} \sum_{i = 1}^{k} \sum_{l = 1}^k \mathbb{E}_{H_0}\left[ \left(- \frac{2 \upsilon_i}{\tau^2} + \dot{\mathbf{z}}_i^2 \psi_i\right) \left(- \frac{2 \upsilon_l}{\tau^2} + \dot{\mathbf{z}}_l^2 \psi_l \right)\right] \\
  &= \frac{1}{4} \sum_{i = 1}^k \sum_{l = 1}^k \left[ \frac{4 \upsilon_i \upsilon_l}{(\tau^2)^2} - \frac{2 \upsilon_i \psi_l }{\tau^2} \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2\right] - \frac{2 \upsilon_l \psi_i}{\tau^2} \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \right] + \psi_i \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \dot{\mathbf{z}}_l^2 \right] \right] \\
  &= \frac{1}{4} \sum_{i = 1}^k \sum_{l = 1}^k \left[ \frac{4 \upsilon_i \upsilon_l}{(\tau^2)^2} - \frac{2 \upsilon_i \psi_l }{\tau^2} \left( \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \mu_{l,j}^0 \right) - \frac{2 \upsilon_l \psi_i}{\tau^2} \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0 \right)\right] + \frac{1}{4} \sum_{i =1}^k \sum_{l = 1 \\ l \neq i}^k \psi_i \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2 \dot{\mathbf{z}}_l^2 \right] + \frac{1}{4}\sum_{i = 1}^{k} \psi_i^2 \mathbb{E}_{H_0}\left[\dot{\mathbf{z}}_i^4  \right] \\
  &= \frac{1}{4} \sum_{i = 1}^k \sum_{l = 1}^k \left[ \frac{4 \upsilon_i \upsilon_l}{(\tau^2)^2} - \frac{2 \upsilon_i \psi_l }{\tau^2} \left( \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2 \mu_{l,j}^0 \right) - \frac{2 \upsilon_l \psi_i}{\tau^2} \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \mu_{i,j}^0 \right)\right] + \frac{1}{4} \sum_{i =1}^k \sum_{l = 1 \\ l \neq i}^k \psi_i \psi_l \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_i^2\right] \mathbb{E}_{H_0}\left[ \dot{\mathbf{z}}_l^2 \right] \\
  &\hspace{5mm} + \frac{1}{4}\sum_{i = 1}^{k} \psi_i^2 \left(3 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4 ((\mu_{i,j}^0)^2 + \mu_{i,j}^0) + \sum_{j = 1}^{n_i} \sum_{j' = 1 \\ j' \neq j}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \mu_{i,j}^0 \mu_{i,j'}^0 \right)\\
  &= \sum_{i = 1}^k \sum_{l = 1}^k \left[ \frac{\upsilon_i \upsilon_l}{(\tau^2)^2} - \frac{\upsilon_i \psi_l \exp(\alpha_l)}{2\tau^2} \left( \sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right) - \frac{\upsilon_l \psi_i \exp(\alpha_i)}{2\tau^2} \left( \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)\right] + \frac{1}{4} \sum_{i =1}^k \sum_{l = 1 \\ l \neq i}^k \psi_i \psi_l \exp(\alpha_i) \exp(\alpha_l) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2  \right) \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right) \\
  &\hspace{5mm} + \frac{1}{4}\sum_{i = 1}^{k} \psi_i^2 \left(3 (\exp^2(\alpha_i) + \exp(\alpha_i)) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4  + \exp^2(\alpha_i)\sum_{j = 1}^{n_i} \sum_{j' = 1 \\ j' \neq j}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \right)\\
  \end{aligned}
  \nonumber
  $$
  Evaluating the above under the null hypothesis (using the fact that $\upsilon_i = 0$ and $\psi_i = 1$ as $\tau^2 \rightarrow 0$):
  $$
  \begin{aligned}
  \mathbb{E}_{H_0}\left[ \left( \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2}\right)^2 \right]
  &\approx 
  \frac{1}{4} \sum_{i =1}^k \sum_{l = 1 \\ l \neq i}^k \exp(\alpha_i) \exp(\alpha_l) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2  \right) \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right)
  + \frac{1}{4}\sum_{i = 1}^{k}\left(3 (\exp^2(\alpha_i) + \exp(\alpha_i)) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4  + \exp^2(\alpha_i)\sum_{j = 1}^{n_i} \sum_{j' = 1 \\ j' \neq j}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \right)
  \end{aligned}
  \nonumber
  $$
</details>

Putting the above together into a big matrix, we get:

$$
\mathcal{I} \rvert_{H_0} 
=
\begin{bmatrix}
\mathcal{I}_{\alpha, \alpha} \rvert_{H_0} & \mathcal{I}_{\alpha, \tau^2}\rvert_{H_0}  \\
\mathcal{I}^\top_{\alpha, \tau^2}\rvert_{H_0}  & \mathcal{I}_{\tau^2, \tau^2}\rvert_{H_0} 
\end{bmatrix}
$$

where:

$$
\begin{aligned}
\mathcal{I}_{\alpha, \alpha}\rvert_{H_0}  &= 
\begin{bmatrix}
  n_1 \exp(\alpha_1) & 0 & \dots & 0 \\
  0 & n_2 \exp(\alpha_2) & \dots & 0 \\
  \vdots & \vdots & \ddots & \vdots \\
  0 & 0 & \dots & n_k \exp(\alpha_k)
\end{bmatrix}
\hspace{10mm}
\mathcal{I}_{\alpha, \tau^2}\rvert_{H_0}  = 
\begin{bmatrix}
\frac{1}{2}\sum_{i = 1}^{k} \exp(\alpha_i) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)  \\
\vdots \\
\frac{1}{2}\sum_{i = 1}^{k} \exp(\alpha_i) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) 
\end{bmatrix} \\
\mathcal{I}_{\tau^2, \tau^2} \rvert_{H_0} 
&= \begin{bmatrix}
\frac{1}{4} \sum_{i =1}^k \sum_{l = 1 \\ l \neq i}^k \exp(\alpha_i) \exp(\alpha_l) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2  \right) \left(\sum_{j = 1}^{n_l} \mathbf{z}_{l,j}^2\right)
+ \frac{1}{4}\sum_{i = 1}^{k}\left(3 (\exp^2(\alpha_i) + \exp(\alpha_i)) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^4  + \exp^2(\alpha_i)\sum_{j = 1}^{n_i} \sum_{j' = 1 \\ j' \neq j}^{n_i} \mathbf{z}_{i,j}^2 \mathbf{z}_{i,j'}^2 \right)
\end{bmatrix}
\end{aligned}
$$


$$
\begin{aligned}
\frac{d \upsilon_i}{d \alpha_i} &= \frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} &
\frac{d \xi_i}{d \alpha_i}  &= - \frac{\exp(\alpha_i)(\tau^2 \sum_{j =1}^{n_i} \mathbf{z}_{i,j}^2 )^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} 
\\
\frac{d \zeta_i}{d \alpha_i} &= \frac{\tau^2 \exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} &
\frac{d \dot{\mathbf{z}}_i}{d \alpha_i} &=-\exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}
\end{aligned}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d \upsilon_i}{d \alpha_i} 
&= \frac{d}{d \alpha_i} \left[ \frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)}\right] \\
&= \frac{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)(\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2) - 2\left( \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \right)^2}{(2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2))^2} \\
&= \frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \\
\frac{d \xi_i}{d \alpha_i} 
&= \frac{d}{d \alpha_i} \left[ \frac{\tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)} \right] \\
&= \frac{-2\exp(\alpha_i) (\tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{(2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2))^2} \\
&= - \frac{\exp(\alpha_i)(\tau^2 \sum_{j =1}^{n_i} \mathbf{z}_{i,j}^2 )^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \\
\frac{d \zeta_i}{d \alpha_i} 
&= \frac{d}{d \alpha_i} \left[ \frac{\tau^2 \exp(\alpha_i)}{1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}} \right] \\
&= \frac{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})(\tau^2 \exp(\alpha_i)) - (\tau^2 \exp(\alpha_i))^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} \\
&= \frac{\tau^2 \exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} \\
\frac{d \dot{\mathbf{z}}_i}{d \alpha_i} 
&= \frac{d}{d \alpha_i} \left[\sum_{j =1 }^{n_i} \mathbf{z}_{i,j}(\mathbf{y}_{i,j} - \exp(\alpha_i)) \right] \\
&= - \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}
\end{aligned}
\nonumber
$$
</details>


From these identities, we can derive (component-by-component) the information as the negative expected Hessian:

$$
- \mathbb{E}_{H_0}\left[ \frac{d^2  \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i^2}   \right] 
=  n_i \exp(\alpha_i)
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d}{d \alpha_i} \left[ \zeta_i \dot{\mathbf{z}}_i \right]
&= \dot{\mathbf{z}}_i \frac{d \zeta_i}{\partial \alpha_i} + \zeta_i \frac{d \dot{\mathbf{z}}_i}{d \alpha_i} \\
&= \frac{\tau^2 \exp(\alpha_i) \dot{\mathbf{z}}_i}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j})^2} - \zeta_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \\
\frac{d}{d \alpha_i} \left[ \xi_i \dot{\mathbf{z}}_i \right] 
&= \dot{\mathbf{z}}_i \frac{d \xi_i}{d \alpha_i} + \xi_i \frac{d \dot{\mathbf{z}}_i}{d \alpha_i} \\
&= \frac{\dot{\mathbf{z}}_i  \exp(\alpha_i) (\tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} - \xi_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}
\end{aligned}
$$
$$
\begin{aligned}
\frac{d^2  \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i^2}  
&= \frac{d}{d \alpha_i} \left[ \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \right] \\
&\approx \frac{d}{d \alpha_i} \left[ - \upsilon_i + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) - \zeta_i \dot{\mathbf{z}}_i \left( \xi_i \dot{\mathbf{z}}_i + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \right] \\
&= - \frac{d \upsilon_i}{d \alpha_i} - n_i\exp(\alpha_i) - \xi_i \dot{\mathbf{z}}_i \frac{d}{d \alpha_i} \left[ \zeta_i \dot{\mathbf{z}}_i \right] - \zeta_i \dot{\mathbf{z}}_i \frac{d}{d \alpha_i} \left[ \xi_i \dot{\mathbf{z}}_i \right] - \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \frac{d}{d \alpha_i} \left[ \zeta_i \dot{\mathbf{z}}_i \right] \\
&= - \frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} - n_i \exp(\alpha_i) - \left(\xi_i \dot{\mathbf{z}}_i - \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \left[\frac{\tau^2 \exp(\alpha_i) \dot{\mathbf{z}}_i}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} - \zeta_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right] - \zeta_i \dot{\mathbf{z}}_i \left[ \frac{\dot{\mathbf{z}}_i \exp(\alpha_i) (\tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} - \xi_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right] \\
&= - \frac{\tau^2\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\left(1 + \dot{\mathbf{z}}_i^2 \zeta_i \tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2)^2} - \left(\xi_i \dot{\mathbf{z}}_i - \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \left[\frac{\tau^2 \exp(\alpha_i) \dot{\mathbf{z}}_i}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} - \zeta_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right] + \zeta_i \xi_i  \dot{\mathbf{z}}_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} - n_i \exp(\alpha_i) \\
&= - \frac{\tau^2\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\left(1 + \dot{\mathbf{z}}_i^2 \zeta_i \tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2)^2} -  \frac{\left(\xi_i \dot{\mathbf{z}}_i - \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \tau^2 \exp(\alpha_i) \dot{\mathbf{z}}_i}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} + \left( \zeta_i \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}  \right) \left( 2\xi_i  \dot{\mathbf{z}}_i - \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right)  - n_i \exp(\alpha_i)
\end{aligned}
$$
Recall that $$\mathbb{E}_{H_0}[\dot{\mathbf{z}}_i] = 0$$ and $$\mathbb{E}_{H_0} \left[ \dot{\mathbf{z}}_i^2 \right] = \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2$$, so taking the negative expectation of the above yields:
$$
\begin{aligned}
- \mathbb{E}_{H_0}\left[ \frac{d^2  \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i^2}   \right] 
&= \frac{\tau^2\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\left(1 + \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \zeta_i \tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i}\mathbf{z}_{i,j}^2)^2} + \frac{\tau^2 \exp(\alpha_i)^2 \xi_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} + n_i \exp(\alpha_i) \\
\implies
- \mathbb{E}_{H_0}\left[ \frac{d^2  \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i^2} \right]  \bigg\rvert_{H_0}
&=  n_i \exp(\alpha_i)
\end{aligned} 
$$
</details>

$$
- \mathbb{E}_{H_0} \left[ \frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d (\tau^2)^2} \right] \bigg\rvert_{H_0}
= \frac{3}{2} \sum_{i = 1}^{k} \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2  \\
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d}{d \tau^2} \left[ \frac{2 \upsilon_i}{\tau^2} \right] &= \frac{d}{d \tau^2} \left[ \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}\right] \\
&= -\frac{(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \\
\frac{d \psi_i}{d \tau^2} &= \frac{d}{d \tau^2} \left[ \left(1 + 2 \tau^2 \exp(\alpha_i) \sum_{j =1}^{n_i} \mathbf{z}_{i,j}^2 + (\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2 \right)^{-1} \right] \\
&= -\left(1 + 2 \tau^2 \exp(\alpha_i) \sum_{j =1}^{n_i} \mathbf{z}_{i,j}^2 + (\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2\right)^{-2}\left(2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 + 2 \tau^2 \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)^2 \right) \\
&= - 2\psi_i^2 \exp(\alpha_i) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right) \left(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) 
\end{aligned}
\nonumber
$$
It follows that:
$$
\begin{aligned}
\frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d (\tau^2)^2}
&= \frac{d}{d \tau^2} \left[ \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \tau^2} \right] \\
&= \frac{d}{d \tau^2} \left[ \frac{1}{2} \sum_{i = 1}^k \left( - \frac{2 \upsilon_i}{\tau^2} + \dot{\mathbf{z}}_i^2 \psi_i \right) \right] \\
&= \frac{1}{2} \sum_{i = 1}^k \left( \frac{(\exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{( 1+ \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} - 2 \dot{\mathbf{z}}_i^2 \psi_i^2 \exp(\alpha_i) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right) \left(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \right) 
\end{aligned}
\nonumber
$$
And thus:
$$
\begin{aligned}
- \mathbb{E}_{H_0} \left[ \frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d (\tau^2)^2} \right]
&= \frac{1}{2} \sum_{i = 1}^k  \left[ - \frac{(\exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{( 1+ \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} + 2 \mathbb{E}\left[ \dot{\mathbf{z}}_i^2 \right] \psi_i^2 \exp(\alpha_i) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right) \left(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \right] \\
&=\frac{1}{2} \sum_{i = 1}^k  \left[ - \frac{(\exp(\alpha_i)\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}{( 1+ \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} + 2 \psi_i^2 \left(\exp(\alpha_i)  \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)^2 \left(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right) \right] 
\end{aligned}
\nonumber
$$
Evaluating this under $H_0$, we get:
$$
\begin{aligned}
- \mathbb{E}_{H_0} \left[ \frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d (\tau^2)^2} \right] \bigg\rvert_{H_0}
&= \frac{1}{2} \sum_{i = 1}^{k} \left[ -\left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2 + 2 \left( \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2\right] \\
&= \frac{1}{2} \sum_{i = 1}^{k} \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2 \right)^2  \\
\end{aligned}
\nonumber
$$
</details>

$$
\begin{aligned}
\frac{d \upsilon_i}{d\tau^2} &= \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j =1 }^{n_i} \mathbf{z}_{i,j}^2)^2} 
&
\frac{d \xi_i}{d \tau^2} &= \frac{\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2}  \\
\frac{d \zeta_i}{d \tau^2} &= \frac{\exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} 
&
\frac{d \dot{\mathbf{z}}_i}{d \tau^2} &= 0
\end{aligned}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d \upsilon_i}{d\tau^2}
&= \frac{d}{d \tau^2} \left[ \frac{\tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)} \right] \\
&= \frac{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2) - 2\tau^2 \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)^2}{(2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2))^2} \\
&= \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j =1 }^{n_i} \mathbf{z}_{i,j}^2)^2} \\
\frac{d \xi_i}{d \tau^2}
&= \frac{d}{\tau^2} \left[\frac{\tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)} \right] \\
&= \frac{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2) - 2 \tau^2 \exp(\alpha_i) \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right)^2}{\left(2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2) \right)^2} \\
&= \frac{\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2 (1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \\
\frac{d \zeta_i}{d \tau^2}
&= \frac{d}{d \tau^2} \left[ \frac{\tau^2 \exp(\alpha_i)}{1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}}\right] \\
&= \frac{\exp(\alpha_i)(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})  - \exp^2(\alpha_i) \tau^2 \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}}{(1 + \tau^2 \exp(\alpha_i) \sum_{j =1 }^{n_i} \mathbf{z}_{i,j})^2} \\
&= \frac{\exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2} \\
\frac{d \dot{\mathbf{z}}_i}{d \tau^2} 
&= \frac{d}{\tau^2} \left[ \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}(\mathbf{y}_{i,j} - \exp(\alpha_i)) \right] \\
&= 0
\end{aligned}
\nonumber
$$
</details>

$$
-\mathbb{E}_{H_0}\left[\frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i d \tau^2} \right]  \bigg\rvert_{H_0}
= \frac{1}{2}\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i d \tau^2}
&= \frac{d}{d \tau^2} \left[ \frac{d \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i} \right] \\
&= \frac{d}{d \tau^2} \left[ - \upsilon_i + \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \exp(\alpha_i)) - \zeta_i \dot{\mathbf{z}}_i \left( \xi_i \dot{\mathbf{z}}_i + \sum_{j = 1}^{n_i} \mathbf{z}_{i,j} \right) \right] \\
&= - \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} - \dot{\mathbf{z}}_i^2 \left( \xi_i \frac{d \zeta_i}{d \tau^2} + \zeta_i \frac{d \xi_i}{d \tau^2}\right) - \dot{\mathbf{z}}_i \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right) \frac{d \zeta_i}{d \tau^2} \\
&= - \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} - \dot{\mathbf{z}}_i^2 \left( \frac{\xi_i \exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2}+\frac{\zeta_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \right) - \dot{\mathbf{z}}_i \left(\sum_{j = 1}^{n_i} \mathbf{z}_{i,j}\right) \left(\frac{\exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2}\right)\\
\end{aligned}
\nonumber
$$
Taking the negative expectation of the above under $H_0$:
$$
\begin{aligned}
-\mathbb{E}_{H_0}\left[\frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i d \tau^2} \right] 
&= \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} + \mathbb{E}_{H_0}[\dot{\mathbf{z}}_i^2] \left( \frac{\xi_i \exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2}+\frac{\zeta_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \right) \\
&= \frac{\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} + \left(\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2\right) \left( \frac{\xi_i \exp(\alpha_i)}{(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j})^2}+\frac{\zeta_i \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2}{2(1 + \tau^2 \exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2)^2} \right)
\end{aligned}
\nonumber
$$
And evaluating under $H_0$:
$$
-\mathbb{E}_{H_0}\left[\frac{d^2 \ell(\mathbf{y}; \alpha, \tau^2)}{d \alpha_i d \tau^2} \right]  \bigg\rvert_{H_0}
= \frac{1}{2}\exp(\alpha_i) \sum_{j = 1}^{n_i} \mathbf{z}_{i,j}^2
\nonumber
$$
</details>

---

### Poisson Case - Pseudo-Likelihood Approach -->



---

<!-- 
## References

[^fn-fitzmaurice]: Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). Applied longitudinal analysis (Second edition). Wiley. -->
