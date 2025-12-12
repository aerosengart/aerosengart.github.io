---
layout: post
title: Score and Information
description: Calculations for GLMMs
date: 2025-10-29
tabs: true
tags: glmm information score
categories: glmm
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    beginning: true
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
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

{% tabs covar %}

{% tab covar equation %}

$$
\Sigma_{y_i} =
\sigma^2 \mathbb{I}_{n \times n} + \tau^2 \mathbf{z}_i \mathbf{z}_i^\top
$$

{% endtab %}

{% tab covar proof %}

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

{% endtab %}

{% endtabs %}


Since the $\beta_i$ are independent, observations from different clusters have covariance zero. Let $\mathbf{y} = (\mathbf{y}_1, \dots, \mathbf{y}_k)$ denote the full data, $\alpha = \begin{bmatrix} \alpha_1 & \dots & \alpha_k\end{bmatrix}^\top$, $\beta = \begin{bmatrix} \beta_1 & \dots & \beta_k\end{bmatrix}^\top$, and $\theta = (\alpha, \beta)$. The complete, marginal likelihood and log-likelihood are:

$$
\begin{aligned}
\mathcal{L}(\theta; \mathbf{y}) &= \prod_{i = 1}^k (2 \pi)^{-\frac{n}{2}} \rvert \Sigma_{y_i} \rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right) \\
\ell(\theta; \mathbf{y}) &= \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
$$

### Score
We first find the gradient of the log-likelihood with respect to $\theta$ parameter-wise. Using the Sherman-Morrison formula, we can find $\Sigma_{y_i}^{-1}$ to be:

{% tabs sigma-inv %}

{% tab sigma-inv equation %}

$$
\Sigma_{y_i}^{-1}  = \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top
$$

{% endtab %}

{% tab sigma-inv proof %}

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

{% endtab %}

{% endtabs %}



#### Derivatives
Let's find the derivative with respect to $\sigma^2$:

{% tabs ell-deriv-1 %}

{% tab ell-deriv-1 equation %}

$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2}
    = -\frac{1}{2}\sum_{i = 1}^k \left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
$$

{% endtab %}

{% tab ell-deriv-1 proof %}

$$ 
\begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right] 
    &= \text{tr}\left[ \Sigma_{y_i}^{-1} \frac{\partial}{\partial \sigma^2} \left[\Sigma_{y_i}\right] \right] \\
    &= \text{tr}\left[ \Sigma_{y_i}^{-1} \mathbb{I}_{n \times n} \right] \\
    &= \text{tr}\left[ \Sigma^{-1} \right] \\
    &= \text{tr}\left[\frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &=\text{tr}\left[ \frac{1}{\sigma^2}\mathbb{I}_{n \times n} \right] \text{tr}\left[- \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top\right] \\
    &=\frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]
\end{aligned}
\nonumber
$$

$$ 
\begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \Sigma_{y_i}^{-1} \right]
    &= - \Sigma_{y_i}^{-1} \frac{\partial}{\partial \sigma^2} \left[ \Sigma_{y_i}\right] \Sigma_{y_i}^{-1} \\
    &= -\left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbb{I}_{n \times n} \left[\frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= - \left[ \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} - \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top + \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top 
\end{aligned}
$$

The above imply:

$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2}
    &= \frac{\partial}{\partial \sigma^2} \left[ \sum_{i = 1}^k - \frac{1}{2} \log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2}(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= -\frac{1}{2}\sum_{i = 1}^k \left[ \frac{n}{\sigma^2} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[\mathbf{z}_i \mathbf{z}_i^\top\right]   + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
$$

{% endtab %}

{% endtabs %}

We do the same with $\tau^2$:

{% tabs ell-tau-1 %}

{% tab ell-tau-1 equation %}

$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2}
= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
$$

{% endtab %}

{% tab ell-tau-1 proof %}

$$
\begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right] 
    &= \text{tr}\left[ \Sigma_{y_i}^{-1} \frac{\partial}{\partial \tau^2} \left[\Sigma_{y_i}\right] \right] \\
    &= \text{tr}\left[ \Sigma_{y_i}^{-1} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= \text{tr} \left[ \left( \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top\right) \mathbf{z}_i \mathbf{z}_i^\top\right] \\
    &= \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right]
\end{aligned}
\nonumber
$$

$$ 
\begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i}^{-1} \right]
    &= - \Sigma_{y_i}^{-1} \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i}\right] \Sigma_{y_i}^{-1} \\
    &= -\left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \left[ \mathbf{z}_i \mathbf{z}_i^\top \right] \left[\frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \\
    &= - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top
\end{aligned}
\nonumber
$$

$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2}
    &= \frac{\partial}{\partial \tau^2} \left[ \sum_{i = 1}^k - \frac{1}{2} \log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2}(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{1}{\sigma^2} \text{tr}[\mathbf{z}_i \mathbf{z}_i^\top] - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top  \left[ - \frac{1}{(\sigma^2)^2} \mathbf{z}_i \mathbf{z}_i^\top + \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top  - \frac{(\tau^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top\mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
$$

{% endtab %}

{% endtabs %}

And then take the gradient with respect to $\alpha$:

{% tabs ell-alpha-1 %}

{% tab ell-alpha-1 equation %}

$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} 
= \begin{bmatrix}
(\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} \mathbf{1}_n \\
\vdots \\
(\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} \mathbf{1}_n
\end{bmatrix}
$$

{% endtab %}

{% tab ell-alpha-1 proof %}

We do the computations component-wise:
$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} 
    &= \sum_{i = 1}^k - \frac{1}{2} \frac{\partial}{\partial \alpha_j} \left[ (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \Sigma^{-1}_{y_i} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
    &= - \frac{1}{2} \left(2 (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_i}^{-1}(- \mathbf{1}_n) \right) \\
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} \mathbf{1}_n
\end{aligned}
$$
So then:
$$ 
\begin{aligned}
    \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} 
    &= \begin{bmatrix}
        (\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} \mathbf{1}_n \\
        \vdots \\
        (\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} \mathbf{1}_n
    \end{bmatrix}
\end{aligned}
$$

{% endtab %}

{% endtabs %}

#### MLEs
We can then find the MLE vector, $\hat{\theta}$, by setting the above equations equal to zero, substituting $\tau^2 = 0$, and solving. The MLE vector, $\hat{\theta}$ is:

{% tabs mle-1 %}

{% tab mle-1 equation %}

$$ 
\hat{\theta} = \begin{bmatrix}
\frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{1,j} \\
\vdots \\
\frac{1}{n} \sum_{j = 1}^n \mathbf{y}_{k,j} \\
\frac{1}{nk} \sum_{i = 1}^k (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
0
\end{bmatrix}
$$

{% endtab %}

{% tab mle-1 proof %}

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

{% endtab %}

{% endtabs %}

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

### Information
To find the information, we need to compute the second-order derivatives of the log-likelihood, take the expectation under $H_0$ of minus those quantities, and evaluate them by plugging in $\hat{\theta}$.

#### Derivatives
We start by taking the derivative with respect to $\theta$ (component-wise) of the first derivative with respect to $\sigma^2$:

{% tabs deriv-theta-sigma %}

{% tab deriv-theta-sigma equation %}

$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] &=
- \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] &= 
-\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\sigma^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{(\sigma^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-2\tau^2(2\sigma^2 + 3\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
&= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ -\frac{1}{(\sigma^2)^2} \mathbb{I}_{n \times n} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n 
\end{aligned}
$$

{% endtab %}

{% tab deriv-theta-sigma proof %}

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
$$

{% endtab %}

{% endtabs %}

Next, we do the same with the derivative with respect to $\tau^2$:

{% tabs deriv-theta-tau %}

{% tab deriv-theta-tau equation %}

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[\frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right]
    &= -\frac{1}{2}\sum_{i = 1}^k - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ \frac{2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{2\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\mathbf{z}_i \mathbf{z}_i^\top\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right]
    &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ - \frac{1}{(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\end{aligned}
$$

{% endtab %}

{% tab deriv-theta-tau proof %}

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
$$

{% endtab %}

{% endtabs %}

And finally with $\alpha_j$:

{% tabs deriv-theta-alpha %}

{% tab deriv-theta-alpha equation %}

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

{% endtab %}

{% tab deriv-theta-tau proof %}

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
$$

{% endtab %}

{% endtabs %}


#### Expectations
We take the expectation under the null of all of the terms we found in the previous section. First we do the $\sigma^2$ terms:

{% tabs expectation-sigma-1 %}
{% tab expectation-sigma-1 equation %}

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right] &=  
    - \frac{1}{2}\sum_{i = 1}^k -\frac{n}{(\sigma^2)^2} + \frac{\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top\right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbb{I}_{n \times n} - \frac{2\tau^2(3\sigma^2 + 2\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2}\mathbf{z}_i \mathbf{z}_i^\top + \frac{2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z})^3}\mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right]
    &= -\frac{1}{2}\sum_{i = 1}^k -\frac{1}{(\sigma^2)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \right] - \frac{-\tau^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \text{tr}\left[ \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] + \text{tr}\left[ \left[ \frac{2}{(\sigma^2)^3} \mathbf{z}_i \mathbf{z}_i^\top + \frac{-2\tau^2(2\sigma^2 + 3\tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^2} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top - \frac{-2(\tau^2)^2(2\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)}{(\sigma^2)^3 (\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)^3} \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \Sigma_{y_i} \right] \\
\mathbb{E}\left[ \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right]
    &= 0 
\end{aligned}
$$

{% endtab %}

{% tab expectation-sigma-1 proof %}

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
$$

{% endtab %}

{% endtabs %}

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
$$

And finally the $\alpha_j$ ones:

{% tabs expectation-alpha-1 %}

{% tab expectation-alpha-1 equation %}

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2} \right] \right] &= 0 \\
\mathbb{E}\left[ \frac{\partial}{\partial \alpha} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \right] &= 0 \\
\mathbb{E}\left[ \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right] &= - \mathbf{1}_n^\top  \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 \mathbf{z}_i^\top \mathbf{z}_i)} \mathbf{z}_i \mathbf{z}_i^\top \right] \mathbf{1}_n \\
\mathbb{E}\left[\frac{\partial}{\partial \alpha_{j'}} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \right] &= 0 
\end{aligned}
$$

{% endtab %}

{% tab expectation-alpha-1 proof %}

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
$$

{% endtab %}

{% endtabs %}


We then evaluate the Fisher information at the MLEs:

{% tabs fisher-1 %}

{% tab fisher-1 equation %}

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
$$

{% endtab %}

{% tab fisher-1 proof %}

Note that:

$$ 
\begin{aligned}
    \Sigma_{y_i} \bigg\rvert_{\theta = \hat{\theta}} &= \hat{\sigma}^2 \mathbb{I}_{n \times n}
\end{aligned}
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
$$

{% endtab %}

{% endtabs %}
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

To do so, we need an estimate of $\alpha_i$ under $H_0$, which we can do with  <a href="/posts/2025/06/03/glm.html#weighted-least-squares">iteratively reweighted least squares</a> or some comparable algorithm. With these estimates, we can compute $$\hat{\sigma}_{i,j}^2$$ and $$\mathbf{y}^\star_{i,j}$$ and proceed as in the Gaussian case but without $\sigma^2$ in the parameter vector since we assume them to be fixed at $$\hat{\sigma}^2_{i,j}$$. 

Since $g(\cdot) = \log(\cdot)$, we have that $$\delta(\hat{\eta}_{i,j}) = \exp(\hat{\eta}_{i,j})$$, implying that the working error variances are:

$$
\frac{V(\hat{\mu}_{i,j})}{\delta^2(\hat{\eta}_{i,j})} = \frac{\exp(\hat{\eta}_{i,j}) + \frac{1}{\phi}\exp(\hat{\eta}_{i,j})}{\exp^2(\hat{\eta}_{i,j})} = \frac{1}{\exp(\hat{\eta}_{i,j})}\left(1 + \frac{1}{\phi}\right)
$$

where we recall that we assume $\phi$ is know and $\hat{\eta}_{i,j} = \hat{\alpha}_i$ under $H_0$. Dropping the star superscript, the likelihood and log-likelihood functions we will work with are giving by:

$$
\begin{aligned}
\mathcal{L}(\theta; \mathbf{y}) &= \prod_{i = 1}^k (2 \pi)^{-\frac{n}{2}} \rvert \Sigma_{y_i} \rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right) \\
\ell(\theta; \mathbf{y}) &= \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
$$

but using the working responses and their respective covariances and whatnot. 



### Score
The marginal covariance matrix is very similar to the Gaussian outcome model above. The only thing that has changed is that each error has its own variance:

$$
\Sigma_{y_i} = \text{diag}([\hat{\sigma}^2_{i,1}, \dots, \hat{\sigma}^2_{i, n}]) + \tau^2 \mathbf{z}_i \mathbf{z}_i^\top
$$

Its inverse, $\Sigma^{-1}_{y_i}$, is:

{% tabs sigma-inv-2 %}

{% tab sigma-inv-2 equation %}

$$
\Sigma^{-1}_{y_i} = \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} 
$$

{% endtab %}

{% tab sigma-inv-2 proof %}

First, let's let $$\mathbf{w}_i = \tau^2 \mathbf{z}_i$$, which is the $n$-vector $\mathbf{z}_i$ where each coordinate has been multiplied by $\tau^2$. We'll also let $$\hat{\sigma}^2 = (\hat{\sigma}^2_{i,1}, \dots, \hat{\sigma}^2_{i, n})^\top$$, the vector of the error variances for cluster $i$, and $\frac{1}{\hat{\sigma}^2}$ will be the vector of the reciprocals of the coordinates of $\hat{\sigma}^2$. Using the <a href="https://en.wikipedia.org/wiki/Sherman–Morrison_formula">Sherman-Morrison formula</a>, we have:

$$
\begin{aligned}
\Sigma^{-1}_{y_i} &= \left(\text{diag}(\hat{\sigma}^2) - \mathbf{w}_i \mathbf{v}_i^\top\right)^{-1} \\
&= \text{diag}^{-1}(\hat{\sigma}^2) - \frac{\text{diag}^{-1}(\hat{\sigma}^2) \mathbf{w}_i \mathbf{z}_i^\top \text{diag}^{-1}(\hat{\sigma}^2)}{1 + \mathbf{z}_i^\top \text{diag}^{-1}(\hat{\sigma}^2) \mathbf{w}_i} \\
&= \text{diag}\left(\frac{1}{\hat{\sigma}^2}\right) - \frac{\text{diag}\left(\frac{\tau^2}{(\hat{\sigma}^2)^2}\right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \mathbf{z}_i^\top \text{diag}\left(\frac{\tau^2}{\hat{\sigma}^2}\right) \mathbf{z}_i} \\ 
&= \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma_{i,1}^2}, \dots, \frac{\tau^2}{\sigma_{i,n}^2} \right] \right) \mathbf{z}_i} \\
&= \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} 
\end{aligned}
$$

We also have:

$$
[\Sigma^{-1}_{y_i}]_{j,j'} = - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)}
\hspace{10mm}
\text{and}
\hspace{10mm}
[\Sigma^{-1}_{y_i}]_{j,j} = \frac{1}{\hat{\sigma}^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)}
$$

{% endtab %}

{% endtabs %}


#### Derivatives
We first find the derivative with respect to $\tau^2$:

{% tabs deriv-tau-2 %}

{% tab deriv-tau-2 equation %}

$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} = 
- \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} - 
\sum_{l = 1}^n \sum_{j = 1}^n \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}\right)} - 
\sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2}
\right]
$$

{% endtab %}

{% tab deriv-tau-2 proof %}

First, the log determinant term:

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \log (\rvert \Sigma_{y_i} \rvert) \right]
&= \text{tr}\left[ \Sigma^{-1}_{y_i} \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i} \right] \right] \\
&= \text{tr} \left[ \left( \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right) \mathbf{z}_i \mathbf{z}_i^\top \right] \\
&= \sum_{l = 1}^n \sum_{j = 1}^n \left( \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right)_{l,j} \left(\mathbf{z}_i \mathbf{z}_i^\top \right)_{l,j} \\
&= \sum_{l = 1}^n \left( \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right)_{l,l} \left(\mathbf{z}_i \mathbf{z}_i^\top \right)_{l,l} 
+ \sum_{l = 1}^n \sum_{j \neq l} \left( \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right)_{l,j} \left(\mathbf{z}_i \mathbf{z}_i^\top \right)_{l,j} \\
&= \sum_{l = 1}^n \left( \frac{1}{\hat{\sigma}^2_{i,l}} - \frac{\frac{\tau^2}{(\hat{\sigma}_{i,l}^2)^2}\mathbf{z}_{i,l}^2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right) \mathbf{z}_{i,l}^2
- \sum_{l = 1}^n \sum_{j \neq l} \left( \frac{\frac{\tau^2}{(\hat{\sigma}^2_{i,l})^2} \mathbf{z}_{i,l} \mathbf{z}_{i,j}}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right) \mathbf{z}_{i,l} \mathbf{z}_{i,j} \\
&= \sum_{l = 1}^n \left( \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} - \frac{\frac{\tau^2}{(\hat{\sigma}_{i,l}^2)^2}\mathbf{z}_{i,l}^4}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right) 
- \sum_{l = 1}^n \sum_{j \neq l} \left( \frac{\frac{\tau^2}{(\hat{\sigma}^2_{i,l})^2} \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} \right) \\
&= \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} - 
\sum_{l = 1}^n \sum_{j = 1}^n \frac{\frac{\tau^2}{(\hat{\sigma}_{i,l}^2)^2}\mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}}
\end{aligned}
$$

Next, the quadratic term. We will compute the derivative of $\Sigma^{-1}_{y_i}$ with respect to $\tau^2$ element-wise:

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2}\left[ [\Sigma^{-1}_{y_i}]_{j,j'}\right] 
&= \frac{\partial}{\partial \tau^2} \left[ - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)} \right] \\
&= -\frac{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)(\mathbf{z}_{i,j}\mathbf{z}_{i,j'}) - (\tau^2 \mathbf{z}_{i,j}\mathbf{z}_{i,j'})\left( (\hat{\sigma}_{i,j}^2)^2\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)}{(\hat{\sigma}_{i,j}^2)^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
&= -\frac{\mathbf{z}_{i,j}\mathbf{z}_{i,j'} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}- \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
&= -\frac{\mathbf{z}_{i,j}\mathbf{z}_{i,j'}}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
\frac{\partial}{\partial \tau^2}\left[ [\Sigma^{-1}_{y_i}]_{j,j}\right] 
&= \frac{\partial}{\partial \tau^2} \left[ \frac{1}{\hat{\sigma}^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)} \right] \\
&= - \frac{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)\mathbf{z}_{i,j}^2 - \tau^2 \mathbf{z}_{i,j}^2 (\hat{\sigma}_{i,j}^2)^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}}{(\hat{\sigma}_{i,j}^2)^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
&= - \frac{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)\mathbf{z}_{i,j}^2 - \tau^2 \mathbf{z}_{i,j}^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
&= - \frac{\mathbf{z}_{i,j}^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}- \tau^2\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
&= - \frac{\mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \\
\end{aligned}
$$

In matrix notation, we have:

$$
\frac{\partial}{\partial \tau^2}\left[ \Sigma^{-1}_{y_i} \right]
= \text{diag}\left(\left[ -\frac{1}{(\hat{\sigma}^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2}, \dots, -\frac{1}{(\hat{\sigma}^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top
$$

Then:

$$
\begin{aligned}
(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \frac{\partial}{\partial \tau^2} \left[ \Sigma^{-1}_{y_i} \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) 
&= (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ 
\text{diag}\left(\left[ -\frac{1}{(\hat{\sigma}^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2}, \dots, -\frac{1}{(\hat{\sigma}^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top
\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
&= \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i)  \left[ \text{diag}\left(\left[ -\frac{1}{(\hat{\sigma}^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2}, \dots, -\frac{1}{(\hat{\sigma}^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top
\right]_{a,b} \\
&= - \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2}
\end{aligned}
$$

And thus:

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \ell(\theta; \mathbf{y}) \right]
&= \frac{\partial}{\partial \tau^2} \left[ \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \right] \\
&= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{\partial}{\partial \tau^2} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i}^{-1} \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]\\
&= - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} - 
\sum_{l = 1}^n \sum_{j = 1}^n \frac{\frac{\tau^2}{(\hat{\sigma}_{i,l}^2)^2}\mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}} - 
\sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2}
\right]
\end{aligned}
$$

{% endtab %}

{% endtabs %}

And now we find the gradient with respect to $\alpha$:

{% tabs deriv-alpha-2 %}

{% tab deriv-alpha-2 equation %}

$$
\frac{\partial}{\partial \alpha}[ \ell(\theta; \mathbf{y})] = \begin{bmatrix}
(\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \Sigma_{y_1}^{-1} \mathbf{1}_n \\
\vdots \\
(\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \Sigma_{y_k}^{-1} \mathbf{1}_n 
\end{bmatrix}
$$

{% endtab %}

{% tab deriv-alpha-2 proof %}

We do this component-wise:

$$
\begin{aligned}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} 
&= \frac{\partial}{\partial \alpha_j} \left[ \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \right] \\
&= - \frac{1}{2} \frac{\partial}{\partial \alpha_j} \left[(\mathbf{y}_j -  \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} (\mathbf{y}_j - \alpha_j \mathbf{1}_n)  \right] \\
&= - \frac{1}{2}\left(2 (\mathbf{y}_j -  \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} (-\mathbf{1}_n) \right) \\
&= (\mathbf{y}_i - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\implies
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} 
&= \begin{bmatrix}
(\mathbf{y}_1 - \alpha_1 \mathbf{1}_n)^\top \Sigma_{y_1}^{-1} \mathbf{1}_n \\
\vdots \\
(\mathbf{y}_k - \alpha_k \mathbf{1}_n)^\top \Sigma_{y_k}^{-1} \mathbf{1}_n 
\end{bmatrix}
\end{aligned}
$$

{% endtab %}
{% endtabs %}


#### MLEs
We can find $\hat{\theta}$ by setting the above equal to zero and substitute $\tau^2 = 0$. We get:

{% tabs mle-2 %}

{% tab mle-2 equation %}

$$
\hat{\theta} = \begin{bmatrix} 
\left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{1,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{1,l'}}{\hat{\sigma}^2_{1,l'}}  \\
\vdots \\
\left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{k,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{k,l'}}{\hat{\sigma}^2_{k,l'}}  \\
0
\end{bmatrix}
$$

{% endtab %}

{% tab mle-2 proof %}

We only need to deal with $\alpha$, which we can again do component-wise. First notice that:

$$
\begin{aligned}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n  \\
&= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left(  \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{j, 1}}, \dots, \frac{1}{\hat{\sigma}^2_{j,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{j,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{j,n}^2)^2} \right] \right) \mathbf{z}_j \mathbf{z}_j^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}^2_{j,l}}} \right) \mathbf{1}_n \\
&= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left(\begin{bmatrix} 
\frac{1}{\hat{\sigma}^2_{j, 1}} \\
\vdots \\
\frac{1}{\hat{\sigma}^2_{j,n}}
\end{bmatrix} - \left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right)\text{diag}\left( \left[ \frac{\tau^2}{(\hat{\sigma}_{j,1}^2)^2}, \dots, \frac{\tau^2}{(\hat{\sigma}_{j,n}^2)^2} \right] \right) \begin{bmatrix} 
\sum_{l = 1}^n \mathbf{z}_{j,1} \mathbf{z}_{j,l} \\
\vdots \\
\sum_{l = 1}^n \mathbf{z}_{j,n} \mathbf{z}_{j,l}
\end{bmatrix} \right) \\
&= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\hat{\sigma}^2_{j,l}} \right)- 
\left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \begin{bmatrix} 
\frac{\tau^2}{(\hat{\sigma}_{j,1}^2)^2} \sum_{l = 1}^n \mathbf{z}_{j,1} \mathbf{z}_{j,l} \\
\vdots \\
\frac{\tau^2}{(\hat{\sigma}_{j,n}^2)^2}  \sum_{l = 1}^n \mathbf{z}_{j,n} \mathbf{z}_{j,l}
\end{bmatrix} \\
&= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\hat{\sigma}^2_{j,l}} \right)- 
\left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \sum_{l' = 1}^n (\mathbf{y}_{j,l'} - \alpha_j) 
\frac{\tau^2}{(\hat{\sigma}_{j,l'}^2)^2} \sum_{l = 1}^n \mathbf{z}_{j,l'} \mathbf{z}_{j,l}  \\
&= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\hat{\sigma}^2_{j,l}} \right)- 
\left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \left( \sum_{l' = 1}^n 
\frac{\tau^2 (\mathbf{y}_{j,l'} - \alpha_j) \mathbf{z}_{j,l'} }{(\hat{\sigma}_{j,l'}^2)^2}  \right)
\end{aligned}
$$

So then we solve for $\alpha_j$ in:

$$
\begin{aligned}
0 &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\implies 0 &= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\hat{\sigma}^2_{j,l}} \right)- 
\left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \left( \sum_{l' = 1}^n 
\frac{\tau^2 (\mathbf{y}_{j,l'} - \alpha_j) \mathbf{z}_{j,l'} }{(\hat{\sigma}_{j,l'}^2)^2}  \right)  \\
\implies 0 &= \sum_{l = 1}^n \frac{\mathbf{y}_{j,l}}{\hat{\sigma}^2_{j,l}} - \alpha_j \sum_{l = 1}^n \frac{1}{\hat{\sigma}^2_{j,l}} - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \left[ \sum_{l' = 1}^n \frac{\tau^2 \mathbf{y}_{j,l'} \mathbf{z}_{j,l'}}{(\hat{\sigma}^2_{j,l'})^2} - \alpha_j \sum_{l' =1}^n \frac{\tau^2 \mathbf{z}_{j,l'}}{(\hat{\sigma}^2_{j,l'})} \right] \\
\implies \alpha_j \left[ \sum_{l = 1}^n \frac{1}{\hat{\sigma}^2_{j,l}} - \sum_{l = 1}^n \frac{\tau^2 \mathbf{z}_{j,l'}}{(\hat{\sigma}^2_{j,l'})^2}\right] &= \sum_{l = 1}^n \frac{\mathbf{y}_{j,l}}{\hat{\sigma}^2_{j,l}} - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \sum_{l' = 1}^n\frac{\tau^2 \mathbf{y}_{j,l'} \mathbf{z}_{j,l'}}{(\hat{\sigma}^2_{j,l'})^2} \\
\implies \alpha_j \left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{j,l'}} \left[1 - \frac{\tau^2 \mathbf{z}_{j,l'}}{\hat{\sigma}^2_{j,l'}} \right]\right) &= \sum_{l' = 1}^n  \frac{\mathbf{y}_{j,l'}}{\hat{\sigma}^2_{j,l'}} \left[ 1 - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \frac{\tau^2 \mathbf{z}_{j,l'}}{\hat{\sigma}^2_{j,l'}} \right] \\
\implies \alpha_j &= \left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{j,l'}} \left[1 - \frac{\tau^2 \mathbf{z}_{j,l'}}{\hat{\sigma}^2_{j,l'}} \right]\right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{j,l'}}{\hat{\sigma}^2_{j,l'}} \left[ 1 - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\hat{\sigma}_{j,l}^2}} \right) \frac{\tau^2 \mathbf{z}_{j,l'}}{\hat{\sigma}^2_{j,l'}} \right]
\end{aligned}
$$

Under $H_0$, $\tau^2 = 0$, so we get:

$$
\hat{\alpha}_j = \left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{j,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{j,l'}}{\hat{\sigma}^2_{j,l'}} 
$$

{% endtab %}

{% endtabs %}

Thus, the score evaluated at $\theta = \hat{\theta}$ is then:

$$
U_{\theta}(\hat{\theta}) = \begin{bmatrix}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \bigg\rvert_{\theta = \hat{\theta}} \\
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \bigg\rvert_{\theta = \hat{\theta}} 
\end{bmatrix}
= \begin{bmatrix}
\left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{1,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{1,l'}}{\hat{\sigma}^2_{1,l'}}  \\
\vdots \\
\left(\sum_{l' = 1}^n \frac{1}{\hat{\sigma}^2_{k,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{k,l'}}{\hat{\sigma}^2_{k,l'}}  \\
- \frac{1}{2}\sum_{i = 1}^k \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}_{i,l}^2} + \sum_{l = 1}^n \sum_{l' = 1}^n \frac{\mathbf{z}_{i,l}\mathbf{z}_{i,l'}(\mathbf{y}_{i,l} - \hat{\alpha}_i)(\mathbf{y}_{i,l'} - \hat{\alpha}_i)}{\hat{\sigma}^2_{i,l} \hat{\sigma}^2_{i,l'}}\right]
\end{bmatrix} 
$$

### Information
As before, to find the information, we need to compute the second-order derivatives of the log-likelihood, take the expectation under $H_0$ of minus those quantities, and evaluate them by plugging in $\hat{\theta}$.

#### Derivatives
We'll take all of the derivatives component-wise. We'll start with those with respect to $\tau^2$. 

{% tabs deriv-theta-tau-2 %}

{% tab deriv-theta-tau-2 equation %}

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \theta}
&= \begin{bmatrix}
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{1,l}}{\hat{\sigma}^2_{1,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{1,l'}(\mathbf{y}_{1,l'} - \alpha_1)}{(\hat{\sigma}^2_{1,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{1,l''}\right) \\
    \vdots \\
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{k,l}}{\hat{\sigma}^2_{k,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{k,l'}(\mathbf{y}_{k,l'} - \alpha_k)}{(\hat{\sigma}^2_{k,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{k,l''}\right) \\
    \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\hat{\sigma}^2_{i,a}} \right)^2 \right]
\end{bmatrix}
\end{aligned}
\label{eq:info-tau-tau}
$$

{% endtab %}

{% tab deriv-theta-tau-2 proof %}

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial (\tau^2)^2}
    &= \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \\
    &= \frac{\partial}{\partial \tau^2} \left[ - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} - \sum_{l = 1}^n \sum_{j = 1}^n \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}\right)} - \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right] \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k  \left[\sum_{l = 1}^n \sum_{j = 1}^n \frac{\partial}{\partial \tau^2} \left[  \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}\right)} \right] + \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\partial}{\partial \tau^2} \left[ \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right] \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[ 
        \sum_{l = 1}^n \sum_{j = 1}^n \left(\frac{\mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2 (\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}\right) - \tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2 (\hat{\sigma}^2_{i,l})^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} }{\left((\hat{\sigma}^2_{i,l})^2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}\right)\right)^2}\right) + \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{-2\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}_{i,l}}\right) \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}_{i,l}^2}}{\left((\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)^2\right)^2}\right)
    \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\hat{\sigma}^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)}{(\hat{\sigma}^2_{i,a})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)^3}\right) \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{(\hat{\sigma}^2_{i,l})^2} \right) \left( \sum_{j = 1}^n \mathbf{z}_{i,j}^2 \right)
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{(\hat{\sigma}^2_{i,a})^2} \right) \left( \sum_{b = 1}^n  (\mathbf{y}_{i,b} - \alpha_i) \mathbf{z}_{i,b} \right) \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{(\hat{\sigma}^2_{i,l})^2} \right) \left( \sum_{l = 1}^n \mathbf{z}_{i,l}^2 \right)
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\hat{\sigma}^2_{i,a}} \right) \left( \sum_{a = 1}^n  (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\hat{\sigma}^2_{i,a}} \right) \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\hat{\sigma}^2_{i,a}} \right)^2 \right]
    
\end{aligned}
$$

Next, we take the derivative (with respect to $\tau^2$) of the derivative with respect to $\alpha_j$:

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2  \partial \alpha_j }
    &= \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \\
    &= \frac{\partial}{\partial \tau^2} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \right] \\
    &= \begin{bmatrix}
        (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ \text{diag}\left( \begin{bmatrix} -\frac{1}{(\hat{\sigma}^2_{j,1})^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} & \dots & -\frac{1}{(\hat{\sigma}^2_{j,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} \end{bmatrix} \right) \mathbf{z}_{j} \mathbf{z}_j^\top  \right] \mathbf{1}_n
    \end{bmatrix} \\
    &= \begin{bmatrix}
        -\frac{\mathbf{y}_{j, 1} - \alpha_j}{(\hat{\sigma}^2_{j,1})^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} & \dots & -\frac{\mathbf{y}_{j,n} - \alpha_j}{(\hat{\sigma}^2_{j,n})^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2}
    \end{bmatrix} \mathbf{z}_j \mathbf{z}_j^\top \mathbf{1}_n \\
    &= \left(- \sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'}(\mathbf{y}_{j,l'} - \alpha_j)}{(\hat{\sigma}^2_{j,l'})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 }\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right) \\
    &= - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'}(\mathbf{y}_{j,l'} - \alpha_j)}{(\hat{\sigma}^2_{j,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right)
\end{aligned}
$$

Putting the two together into a matrix yield:

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \theta }
&= \frac{\partial}{\partial \tau^2} \left[ \begin{bmatrix} \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \\ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \end{bmatrix} \right] \\
&= \begin{bmatrix}
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{1,l}}{\hat{\sigma}^2_{1,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{1,l'}(\mathbf{y}_{1,l'} - \alpha_1)}{(\hat{\sigma}^2_{1,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{1,l''}\right) \\
    \vdots \\
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{k,l}}{\hat{\sigma}^2_{k,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{k,l'}(\mathbf{y}_{k,l'} - \alpha_k)}{(\hat{\sigma}^2_{k,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{k,l''}\right) \\
    \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\hat{\sigma}^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)}{(\hat{\sigma}^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)^3}\right) \right]
\end{bmatrix}
\end{aligned}
$$

{% endtab %}

{% endtabs %}


We then do the same but with respect to $\alpha$.

{% tabs deriv-theta-alpha-2 %}

{% tab deriv-theta-alpha-2 equation %}

$$
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \theta}
=
\begin{bmatrix}
0 \\
\vdots \\
- \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\vdots \\
0 \\
-\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right)
\end{bmatrix}
$$

{% endtab %}

{% tab deriv-theta-alpha-2 proof %}

First, we find the gradient (with respect to $\alpha$) of the derivative with respect to $\tau^2$. This should be equivalent to the corresponding components of Eq. \eqref{eq:info-tau-tau}.

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \tau^2 } 
&= \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \\
&= \frac{\partial}{\partial \alpha_j} \left[  - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} - \sum_{l = 1}^n \sum_{j = 1}^n \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\hat{\sigma}_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}}\right)} - \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right]  \right] \\
&= \frac{\partial}{\partial \alpha_j} \left[ \frac{1}{2} \sum_{i = 1}^k \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\hat{\sigma}^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}}\right)^2} \right] \\
&= \frac{1}{2} \sum_{a = 1}^n \sum_{b = 1}^n \frac{\partial}{\partial \alpha_j} \left[ (\mathbf{y}_{j,a} - \alpha_j)(\mathbf{y}_{j,b} - \alpha_j) \frac{\mathbf{z}_{j,a} \mathbf{z}_{j,b}}{(\hat{\sigma}^2_{j,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} \right] \\
&= \frac{1}{2} \sum_{a = 1}^n \sum_{b = 1}^n \frac{\mathbf{z}_{j,a} \mathbf{z}_{j,b}}{(\hat{\sigma}^2_{j,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} \frac{\partial}{\partial \alpha_j} \left[ \mathbf{y}_{j,a} \mathbf{y}_{j,b} - \alpha_j \mathbf{y}_{j,a} - \alpha \mathbf{y}_{j,b} + \alpha_j^2 \right] \\
&= -\frac{1}{2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \sum_{a = 1}^n \sum_{b = 1}^n \left[ \frac{\mathbf{z}_{j,a}\mathbf{z}_{j,b}(\mathbf{y}_{j, a} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} + \frac{\mathbf{z}_{j,a}\mathbf{z}_{j,b}(\mathbf{y}_{j,b} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} \right] \\
&= -\frac{1}{2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \left[ \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} \sum_{b = 1}^n \mathbf{z}_{j,b} + \sum_{b = 1}^n  \frac{\mathbf{z}_{j,b}(\mathbf{y}_{j,b} - \alpha_j)}{(\hat{\sigma}_{j,b}^2)^2}  \sum_{a = 1}^n \mathbf{z}_{j,a}\right] \\
&= -\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right)
\end{aligned}
$$

Now, we find the vector of second derivatives of the log-likelihood with respect to the components of $\alpha$:

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j^2} 
&= \frac{\partial}{\partial \alpha_j} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n\right] \\
&= - \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\frac{\partial^2 \ell(\theta; \mathbf{y})}{ \partial \alpha_{j'} \partial \alpha_j}
&= \frac{\partial}{\partial \alpha_{j'}} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n\right] \\
&= 0 \\
\implies
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \alpha} 
&= \begin{bmatrix}
0 \\
\vdots \\
- \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\vdots \\
0
\end{bmatrix}
\end{aligned}
$$

Putting the two above results together gives us:

$$
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \theta}
=
\begin{bmatrix}
0 \\
\vdots \\
- \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\vdots \\
0 \\
-\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right)
\end{bmatrix}
$$

{% endtab %}

{% endtabs %}

#### Expectations
We now take the expectation of the above vectors. We'll evaluate the second order partial derivatives with respect to $\tau^2$ first.

{% tabs info-tau-tau-2 %}

{% tab info-tau-tau-2 equation %}

$$
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \theta } \right]
= \begin{bmatrix}
0 \\
\vdots \\
0 \\
\frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \right]
\end{bmatrix}
$$

{% endtab %}

{% tab info-tau-tau-2 proof %}

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial (\tau^2)^2 }\right]
&= \mathbb{E}\left[ \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\hat{\sigma}^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)}{(\hat{\sigma}^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)^3}\right) \right] \right] \\
&= \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n \mathbb{E}\left[ (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i)\right] \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\hat{\sigma}^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)}{(\hat{\sigma}^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)^3}\right) \right] \\
&=  \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n \text{Cov}\left(\mathbf{y}_{i,a}, \mathbf{y}_{i,b} \right) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\hat{\sigma}^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)}{(\hat{\sigma}^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\hat{\sigma}^2_{i,l}} \right)^3}\right) \right] \\
&=  \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \right] & \left(\text{obs. ind. under } H_0 \right) \\
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \alpha_j}\right] 
&= \mathbb{E}\left[ - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'}(\mathbf{y}_{j,l'} - \alpha_j)}{(\hat{\sigma}^2_{j,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right) \right] \\
&= - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'} \mathbb{E}\left[ (\mathbf{y}_{j,l'} - \alpha_j)\right] }{(\hat{\sigma}^2_{j,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right)  \\
&= 0
\end{aligned}
$$

{% endtab %}

{% endtabs %}

And then we do the same for $\alpha_j$:

{% tabs info-alpha-alpha-2 %}

{% tab info-alpha-alpha-2 equation %}

$$
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \theta } \right]
= \begin{bmatrix}
0 \\
\vdots \\
- \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n  \\
\vdots \\
0 \\
0
\end{bmatrix}
$$

{% endtab %}

{% tab info-alpha-alpha-2 proof %}

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \tau^2} \right] 
&= \mathbb{E}\left[ -\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\hat{\sigma}_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right) \right] \\
&= -\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\hat{\sigma}^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a} \mathbb{E}\left[\mathbf{y}_{j, a} - \alpha_j\right] }{(\hat{\sigma}_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right) \\
&= 0 \\
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j^2} \right] 
&= \mathbb{E}\left[ - \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \right] \\
&= - \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \alpha_{j'}} \right] 
&= 0
\end{aligned}
$$

{% endtab %}

{% endtabs %}


We then evaluate the Fisher information at the MLE:

{% tabs info-2 %}

{% tab info-2 equation %}

$$
\mathcal{I}_{\theta, \theta}(\hat{\theta})
= \begin{bmatrix}
- \mathbf{1}_n^\top\text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{1,1}}, \dots, \frac{1}{\hat{\sigma}^2_{1,n}}\right]\right) \mathbf{1}_n & 0 & \dots & 0 \\
0 & - \mathbf{1}_n^\top \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{2,1}}, \dots, \frac{1}{\hat{\sigma}^2_{2,n}}\right]\right) \mathbf{1}_n & \dots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \dots & \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2} \right]
\end{bmatrix}
$$

{% endtab %}

{% tab info-2 proof %}

Note that:

$$
\Sigma_{y_i}^{-1} \bigg\rvert_{\theta = \hat{\theta}} = \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{i,1}}, \dots, \frac{1}{\hat{\sigma}^2_{i,n}}\right]\right)
$$

Then we have:

$$
\begin{aligned}
- \mathbb{E}\left. \left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \theta} \right] \right\rvert_{\theta = \hat{\theta}} 
&= - \mathbb{E} \left.\left[ 
    \begin{bmatrix}
    0 \\
    \vdots \\
    0 \\
    \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\hat{\sigma}^2_{i,l}} \right)^2} \right] 
    \end{bmatrix}
\right] \right\rvert_{\theta = \hat{\theta}} \\
&= \begin{bmatrix}
    0 \\
    \vdots \\
    0 \\
    \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2} \right]
\end{bmatrix} \\
- \mathbb{E} \left[  \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \theta} \right] \bigg\rvert_{\theta = \hat{\theta}}  
&= - \mathbb{E} \left. \left[ 
    \begin{bmatrix}
    0 \\
    \vdots \\
    - \mathbf{1}_n^\top \Sigma_{y_j}^{-1} \mathbf{1}_n  \\
    \vdots \\
    0 \\
    0 
    \end{bmatrix} \right] \right\rvert_{\theta = \hat{\theta}} \\
&= \begin{bmatrix}
    0 \\
    \vdots \\
    - \mathbf{1}_n^\top \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{j,1}}, \dots, \frac{1}{\hat{\sigma}^2_{j,n}}\right]\right) \mathbf{1}_n  \\
    \vdots \\
    0 \\
    0
\end{bmatrix}
\end{aligned}
$$

Putting these together into a big matrix:

$$
\begin{aligned}
\mathcal{I}_{\theta, \theta}(\hat{\theta})
&= - \mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \theta \partial \theta^\top} \right] \bigg\rvert_{\theta = \hat{\theta}} \\
&= \begin{bmatrix}
- \mathbf{1}_n^\top\text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{1,1}}, \dots, \frac{1}{\hat{\sigma}^2_{1,n}}\right]\right) \mathbf{1}_n & 0 & \dots & 0 \\
0 & - \mathbf{1}_n^\top \text{diag}\left(\left[ \frac{1}{\hat{\sigma}^2_{2,1}}, \dots, \frac{1}{\hat{\sigma}^2_{2,n}}\right]\right) \mathbf{1}_n & \dots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \dots & \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\hat{\sigma}^2_{i,l})^2} \right]
\end{bmatrix}
\end{aligned}
$$

{% endtab %}

{% endtabs %}


<!-- 

Let $\mathbf{0}_{n \times n}$ denote the $n$ by $n$ matrix of all zeroes, and let $\mathbf{e}_j$ denote the $j$-th standard basis vector in $n$ dimensions (i.e. the $n$-dimensional vector of all zeroes except the $j$-th coordinate that is a $1$). 



First, note that we have:
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \Sigma_{y_i} \right] &= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ 
    \begin{bmatrix}
    \sigma^2_{i,1} & \dots & 0 \\
    \vdots & \ddots & \vdots \\
    0 & \dots & \sigma^2_{i,n}
    \end{bmatrix} 
\right] = 
\begin{bmatrix} 
0 & \dots & 0 & \dots & 0 \\ 
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 1 & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & \dots & 0 \\ 
\end{bmatrix}
= \mathbf{e}_j \mathbf{e}_j^\top \\
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \Sigma_{y_{i'}} \right] &= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ 
    \begin{bmatrix}
    \sigma^2_{i',1} & \dots & 0 \\
    \vdots & \ddots & \vdots \\
    0 & \dots & \sigma^2_{i',n}
    \end{bmatrix} 
\right] = 
\begin{bmatrix} 
0 & \dots & 0 & \dots & 0 \\ 
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & \dots & 0 \\ 
\end{bmatrix}
= \mathbf{0}_{n \times n}
\end{aligned}
\nonumber
$$
Using standard matrix identities, we can compute:
$$
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right] = 
\frac{1}{\sigma^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j}^2}{(\sigma^2_{i,j})^2\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}
\nonumber
$$
<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right]
&= \text{tr}\left[ \Sigma^{-1}_{y_i} \frac{\partial}{\partial \sigma^2_{i,j}} \left[ \Sigma_{y_i} \right] \right] \\
&= \text{tr}\left[ \Sigma^{-1}_{y_i} \mathbf{e}_j \mathbf{e}_j^\top \right] \\
&= \text{tr}\left[ \mathbf{e}_j^\top \Sigma^{-1}_{y_i}  \mathbf{e}_j \right] & \left(\text{linearity of trace}\right) \\
&=  \mathbf{e}_j^\top \Sigma^{-1}_{y_i}  \mathbf{e}_j  & \left(\text{trace of } 1 \times 1 \text{ matrix}\right) \\
&= \left[ \Sigma^{-1}_{y_i} \right]_{j,j} \\
&= \frac{1}{\sigma^2_{i,j}} - \frac{\left[ \text{diag}\left(\left[\frac{\tau^2}{(\sigma^2_{i,1})^2}, \dots, \frac{\tau^2}{(\sigma^2_{i,n})^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top \right]_{j,j}}{1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma^2_{i,1}}, \dots, \frac{\tau^2}{\sigma^2_{i,n}} \right]\right) \mathbf{z}_i} & \left(\text{only numerator is matrix}\right) \\
&= \frac{1}{\sigma^2_{i,j}} - \frac{ \frac{\tau^2}{(\sigma^2_{i,j})^2} \mathbf{z}_{i,j}^2 }{1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma^2_{i,1}}, \dots, \frac{\tau^2}{\sigma^2_{i,n}} \right]\right) \mathbf{z}_i} \\
&= \frac{1}{\sigma^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j}^2}{(\sigma^2_{i,j})^2\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}
\end{aligned}
\nonumber
$$
</details>
$$
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \Sigma^{-1}_{y_i} \right] = \left(\mathbb{I}_{n \times n} - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \right) \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top + \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)^2 \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top 
\nonumber
$$
<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \Sigma^{-1}_{y_i} \right]
&= \Sigma^{-1}_{y_i} \frac{\partial }{\partial \sigma^2_{i,j}} \left[ \Sigma_{y_i} \right] \Sigma^{-1}_{y_i} \\
&= \left( \text{diag}\left(\left[ \frac{1}{\sigma^2_{i, 1}}, \dots, \frac{1}{\sigma^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}} \right) \mathbf{e}_j \mathbf{e}_j^\top \left( \text{diag}\left(\left[ \frac{1}{\sigma^2_{i, 1}}, \dots, \frac{1}{\sigma^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}} \right) \\
&= \left( \text{diag}\left(\left[ \frac{1}{\sigma^2_{i, 1}}, \dots, \frac{1}{\sigma^2_{i,n}} \right]\right)\mathbf{e}_j - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top \mathbf{e}_j}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}} \right) \left( \mathbf{e}_j^\top \text{diag}\left(\left[ \frac{1}{\sigma^2_{i, 1}}, \dots, \frac{1}{\sigma^2_{i,n}} \right]\right) - \frac{\mathbf{e}_j^\top \text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}} \right) \\
&= \left( 
\begin{bmatrix} 0 \\ \vdots \\ \frac{1}{\sigma^2_{i,j}} \\ \vdots \\ 0 \end{bmatrix} 
- \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \begin{bmatrix} 0  \\ \vdots \\ \frac{1}{\sigma^2_{i,j}} \\ \vdots \\ 0 \end{bmatrix} \right) \left(  \begin{bmatrix} 0  \\ \vdots \\ \frac{1}{\sigma^2_{i,j}} \\ \vdots \\ 0 \end{bmatrix}^\top   -  \left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}} \right)\begin{bmatrix} 0 \\ \vdots \\ \frac{\tau^2}{(\sigma^2_{i,j})^2} \\ \vdots \\ 0 \end{bmatrix}^\top \mathbf{z}_i \mathbf{z}_i^\top \right) \\
&= \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top - \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top  - \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)\text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top  \mathbf{z}_i \mathbf{z}_i^\top + \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)^2 \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \\
&= \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top  + \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)^2 \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \\
&= \left(\mathbb{I}_{n \times n} - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \right) \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top + \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)^2 \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top 
\end{aligned}
\nonumber
$$
</details>
We can then find the derivative of the marginal log-likelihood with respect to each $\sigma^2_{i,j}$ as:
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \ell(\theta; \mathbf{y})\right] 
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[\sum_{l = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_l} \rvert) - \frac{1}{2} (\mathbf{y}_l -  \alpha_l \mathbf{1}_n)^\top \Sigma_{y_l}^{-1} (\mathbf{y}_l - \alpha_l \mathbf{1}_n) \right] \right] \\
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) \right] - \frac{\partial}{\partial \sigma^2_{i,j}} \left[ \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \\
&= -\frac{1}{2} \left[ \frac{1}{\sigma^2_{i,j}} - \frac{\tau^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,j})^2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)} + (\mathbf{y}_i - \alpha_i\mathbf{1}_n)^\top \left[ \left(\mathbb{I}_{n \times n} - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \right) \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top + \left(\frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)^2 \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
\nonumber
$$

Note that:

$$
\begin{aligned}
\left(\mathbb{I}_{n \times n} - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_i \mathbf{z}_i^\top \right) \text{diag}\left( \left[ \frac{1}{(\sigma^2_{i,1})^2}, \dots, \frac{1}{(\sigma^2_{i,n})^2} \right]\right)\mathbf{e}_j \mathbf{e}_j^\top
&= \begin{bmatrix}
1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,1}^2 & \dots & - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)  \mathbf{z}_{i,1}\mathbf{z}_{i,n} \\
\vdots & \ddots & \vdots \\
- \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)\mathbf{z}_{i,n} \mathbf{z}_{i,1} & \dots & 1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,n}^2
\end{bmatrix} \mathbf{e}_j \mathbf{e}_j^\top \\
&= \begin{bmatrix}
- \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)  \mathbf{z}_{i,1}\mathbf{z}_{i,j} \\
\vdots \\
1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,j}^2
\\
\vdots \\
- \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)  \mathbf{z}_{i,n}\mathbf{z}_{i,j} \\
\end{bmatrix} \mathbf{e}_j^\top \\
&= \begin{bmatrix}
0 & \dots & 1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,j}^2 & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,j}^2 & \dots & 0
\end{bmatrix}
\end{aligned}
\nonumber
$$

Thus:

$$
\begin{aligned}
(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \begin{bmatrix}
0 & \dots & 1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,j}^2 & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 1 - \left(\frac{2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right) \mathbf{z}_{i,j}^2 & \dots & 0
\end{bmatrix} (\mathbf{y}_i - \alpha_i \mathbf{1}_n)
&= \left(1 - \frac{2\mathbf{z}_{i,j}^2}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}}\right)  \sum_{l = 1}^n (\mathbf{y}_{i,l} - \alpha_i)^2
\end{aligned}
$$


 -->

<!-- 

First note that:

$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \text{diag} \left( \left[ \frac{1}{\sigma^2_{i,1}}, \dots, \frac{1}{\sigma^2_{i,n}} \right] \right)\right] 
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ 
    \begin{bmatrix}
    \frac{1}{\sigma^2_{i,1}} & \dots & 0 \\
    \vdots & \ddots & \vdots \\
    0 & \dots & \frac{1}{\sigma^2_{i,n}}
    \end{bmatrix}
\right] 
= \begin{bmatrix} 
0 & \dots & 0 & \dots & 0 \\ 
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & - \frac{1}{(\sigma^2_{i,j})^2} & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & \dots & 0 \\ 
\end{bmatrix} 
=  \text{diag}\left(\left[- \frac{1}{(\sigma^2_{i,1})^2}, \dots, -\frac{1}{(\sigma^2_{i,n})^2} \right] \right) \mathbf{e}_j \mathbf{e}_j^\top \\
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \text{diag} \left( \left[ \frac{\tau^2}{\sigma^2_{i,1}}, \dots, \frac{\tau^2}{\sigma^2_{i,n}} \right] \right)\right] 
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ 
    \begin{bmatrix}
    \frac{\tau^2}{\sigma^2_{i,1}} & \dots & 0 \\
    \vdots & \ddots & \vdots \\
    0 & \dots & \frac{\tau^2}{\sigma^2_{i,n}}
    \end{bmatrix}
\right] 
= \begin{bmatrix} 
0 & \dots & 0 & \dots & 0 \\ 
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & - \frac{\tau^2}{(\sigma^2_{i,j})^2} & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & \dots & 0 \\ 
\end{bmatrix} 
= \text{diag}\left(\left[- \frac{\tau^2}{(\sigma^2_{i,1})^2}, \dots, -\frac{\tau^2}{(\sigma^2_{i,n})^2} \right] \right) \mathbf{e}_j  \mathbf{e}_j^\top \\
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \text{diag} \left( \left[ \frac{\tau^2}{(\sigma^2_{i,1})^2}, \dots, \frac{\tau^2}{(\sigma^2_{i,n})^2} \right] \right)\right] 
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ 
    \begin{bmatrix}
    \frac{\tau^2}{(\sigma^2_{i,1})^2} & \dots & 0 \\
    \vdots & \ddots & \vdots \\
    0 & \dots & \frac{\tau^2}{(\sigma^2_{i,n})^2}
    \end{bmatrix}
\right] 
= \begin{bmatrix} 
0 & \dots & 0 & \dots & 0 \\ 
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & - \frac{2\tau^2}{(\sigma^2_{i,j})^3} & \dots & 0 \\
\vdots & \ddots & \vdots & \ddots & \vdots \\
0 & \dots & 0 & \dots & 0 \\ 
\end{bmatrix} 
= \text{diag}\left(\left[- \frac{2\tau^2}{(\sigma^2_{i,1})^3}, \dots, -\frac{2\tau^2}{(\sigma^2_{i,n})^3} \right] \right) \mathbf{e}_j \mathbf{e}_j^\top
\end{aligned}
\nonumber
$$

So then:

$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right]
&= -\frac{\tau^2 \mathbf{z}_{i,j}^2}{(\sigma^2_{i,j})^2} \\
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top \right]
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ \text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \right] \mathbf{z}_i \mathbf{z}_i^\top 
= \text{diag}\left( \left[-\frac{2\tau^2}{(\sigma^2_{i,1})^3}, \dots, -\frac{2\tau^2}{(\sigma^2_{i,n})^3} \right] \right) \mathbf{e}_j \mathbf{e}_j^\top \mathbf{z}_i \mathbf{z}_i^\top
\end{aligned}
\nonumber
$$ 




$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2_{i,j}} \left[ \Sigma^{-1}_{y_i} \right]
&= \frac{\partial}{\partial \sigma^2_{i,j}} \left[ \text{diag}\left(\left[ \frac{1}{\sigma^2_{i, 1}}, \dots, \frac{1}{\sigma^2_{i,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma_{i,1}^2}, \dots, \frac{\tau^2}{\sigma_{i,n}^2} \right] \right) \mathbf{z}_i} \right] \\
&= \text{diag}\left(\left[- \frac{1}{(\sigma^2_{i,1})^2}, \dots, -\frac{1}{(\sigma^2_{i,n})^2} \right] \right) \mathbf{e}_j\mathbf{e}_j^\top - \frac{\partial}{\partial \sigma^2_{i,j}} \left[ \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma_{i,1}^2}, \dots, \frac{\tau^2}{\sigma_{i,n}^2} \right] \right) \mathbf{z}_i} \right] \\
&= \text{diag}\left(\left[- \frac{1}{(\sigma^2_{i,1})^2}, \dots, -\frac{1}{(\sigma^2_{i,n})^2} \right] \right) \mathbf{e}_j\mathbf{e}_j^\top - \frac{\left(1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma_{i,1}^2}, \dots, \frac{\tau^2}{\sigma_{i,n}^2} \right] \right) \mathbf{z}_i\right) \text{diag}\left( \left[-\frac{2\tau^2}{(\sigma^2_{i,1})^3}, \dots, -\frac{2\tau^2}{(\sigma^2_{i,n})^3} \right] \right) \mathbf{e}_j \mathbf{e}_j^\top \mathbf{z}_i \mathbf{z}_i^\top - \text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{i,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{i,n}^2)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top }{}  \\
\end{aligned}
\nonumber
$$



 -->








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
