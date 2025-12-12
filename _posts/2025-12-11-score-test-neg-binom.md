---
layout: distill
title: Score Calcuations
description: The Negative Binomial Case
date: 2025-12-11
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
    - name: Introduction
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
---

## Introduction

This post is just a catch-all for my derivations for my score test project with negative binomial outcomes. Our set-up is as follows. We have $n$ observations coming from $k$ different clusters, each of size $n_t$ for $t \in [k]$. The full data will be denoted by $\mathbf{y}$. Though $\mathbf{y}$ is a vector, we'll denote the $j$-th observation from cluster $i$ with $\mathbf{y}_{i,j}$. For example, $$\mathbf{y}_{i,j}$$ denotes element $$\sum_{l = 1}^{i - 1} n_l + j$$ of $\mathbf{y}$. We'll also denote the $n_i$-dimensional vector of responses for cluster $i$ with $\mathbf{y}_i$. 

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

## Negative Binomial Case

### Set-Up
In this example, we'll let the responses be negative binomial. To keep things simple, we'll say we only have a single fixed intercept and a single random slope. We let $\phi > 0$ denote the <i>unknown</i> dispersion parameter and assume the conditional mean to be given by:

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

We assume to have the following generalized linear mixed model:

$$
\begin{equation}
\label{eq:glmm-y}
\begin{aligned}
\mathbf{y}_{i,j} \rvert \beta_i &\sim \text{NegBin}(\mu_{i,j}, \phi) \\
\mu_{i,j} &= \exp\left(\eta_{i,j}\right) = \exp\left(\alpha_i + \beta_i \mathbf{z}_{i,j}\right)
\end{aligned}
\end{equation}
$$

### Approximation 
We follow a pseudo-likelihood approach (see <a href="/blog/2025/glmm.html">here</a>) which permits a Gaussian approximation of the outcome distribution via a linearization. Suppose we find the MLEs of the $\alpha_i$ and $\phi$ terms under $H_0$, which we can do with  <a href="/blog/2025/glm.html#weighted-least-squares">iteratively reweighted least squares</a> or some comparable algorithm (basically just fitting the null model). Denote these estimates with $\tilde{\alpha_i}$ and $\tilde{\phi}$, respectively. We can compute our <i>working</i> responses and errors as:

$$
\begin{equation}
\label{eq:working-responses}
\mathbf{y}^\star_{i,j} = \alpha_i + \beta_i \mathbf{z}_{i,j} + \epsilon^\star_{i,j}; \hspace{10mm} 
\epsilon^\star_{i,j} \sim \mathcal{N}\left(0, \frac{V(\tilde{\mu}_{i,j})}{\delta^2(\tilde{\eta}_{i,j})}\right)
\end{equation}
$$

where 

$$
\tilde{\eta}_{i,j} = \tilde{\alpha}_i; \hspace{10mm}
\delta(\tilde{\eta}_{i,j}) = \frac{\partial g^{-1}(\eta_{i,j})}{\partial \eta_{i,j}}\bigg\rvert_{\eta_{i,j} = \tilde{\eta}_{i,j}}
$$

Since $g(\cdot) = \log(\cdot)$, we have that $$\delta(\tilde{\eta}_{i,j}) = \exp(\tilde{\eta}_{i,j})$$, implying that the working error variances are:

$$
\frac{V(\tilde{\mu}_{i,j})}{\delta^2(\tilde{\eta}_{i,j})} = \frac{\exp(\tilde{\eta}_{i,j}) + \frac{1}{\phi}\exp(\tilde{\eta}_{i,j})}{\exp^2(\tilde{\eta}_{i,j})} = \frac{1}{\exp(\tilde{\eta}_{i,j})}\left(1 + \frac{1}{\tilde{\phi}}\right)
$$

where we recall that $\tilde{\eta}_{i,j} = \tilde{\alpha}_i$ under $H_0$. Eq. \eqref{eq:working-responses} describes a linear mixed model (i.e. Gaussian outcomes), defined in terms of the working responses and errors. The only big difference is that we force the error variances to all be different and fix them at:

$$
\sigma_{i,j}^2 = \frac{V(\tilde{\mu}_{i,j})}{\delta^2(\tilde{\eta}_{i,j})} = \frac{1}{\exp(\tilde{\alpha}_i)} \left(1 + \frac{1}{\tilde{\phi}} \right)
$$

We have simplified some of the calculations since we can now use the Gaussian likelihood instead:

$$
\begin{aligned}
\mathcal{L}(\theta; \mathbf{y}^\star) &= \prod_{i = 1}^k (2 \pi)^{-\frac{n}{2}} \rvert \Sigma_{y_i} \rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right) \\
\ell(\theta; \mathbf{y}^\star) &= \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y^\star_i} \rvert) - \frac{1}{2} (\mathbf{y}^\star_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y^\star_i}^{-1} (\mathbf{y}^\star_i - \alpha_i \mathbf{1}_n) \right]
\end{aligned}
$$


In the rest of the post, I'll drop the $\star$ superscript for readability. Just be sure to remember that we are not dealing with the original observations but with their transformations.


### Score
The marginal covariance matrix is very similar to the Gaussian outcome model above. The only thing that has changed is that each error has its own variance:

$$
\Sigma_{y_i} = \text{diag}\left(\begin{bmatrix} \sigma^2_{i,1} & \dots & \sigma^2_{i, n}\end{bmatrix}\right) + \tau^2 \mathbf{z}_i \mathbf{z}_i^\top
$$

Its inverse, $\Sigma^{-1}_{y_i}$, is:

{% tabs sigma-inv-2 %}

{% tab sigma-inv-2 equation %}

$$
\Sigma^{-1}_{y_i} = \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) - \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) \mathbf{z}_i \mathbf{z}_i^\top
$$

{% endtab %}

{% tab sigma-inv-2 proof %}

First, let's let $$\mathbf{w}_i = \tau^2 \mathbf{z}_i$$, which is the $n$-vector $\mathbf{z}_i$ where each coordinate has been multiplied by $\tau^2$. We'll also let $$\sigma^2 = (\sigma^2_{i,1}, \dots, \sigma^2_{i, n})^\top$$, the vector of the error variances for cluster $i$, and $\frac{1}{\sigma^2}$ will be the vector of the reciprocals of the coordinates of $\sigma^2$. Using the <a href="https://en.wikipedia.org/wiki/Sherman–Morrison_formula">Sherman-Morrison formula</a>, we have:

$$
\begin{aligned}
\Sigma^{-1}_{y_i} &= \left(\text{diag}(\sigma^2) - \mathbf{w}_i \mathbf{v}_i^\top\right)^{-1} \\
&= \text{diag}^{-1}(\sigma^2) - \frac{\text{diag}^{-1}(\sigma^2) \mathbf{w}_i \mathbf{z}_i^\top \text{diag}^{-1}(\sigma^2)}{1 + \mathbf{z}_i^\top \text{diag}^{-1}(\sigma^2) \mathbf{w}_i} \\
&= \text{diag}\left(\frac{1}{\sigma^2}\right) - \frac{\text{diag}\left(\frac{\tau^2}{(\sigma^2)^2}\right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \mathbf{z}_i^\top \text{diag}\left(\frac{\tau^2}{\sigma^2}\right) \mathbf{z}_i} \\ 
&= \text{diag}\left(\begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix}  \right) - \frac{\text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2} \end{bmatrix} \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \mathbf{z}_i^\top \text{diag}\left( \left[ \frac{\tau^2}{\sigma_{i,1}^2}, \dots, \frac{\tau^2}{\sigma_{i,n}^2} \right] \right) \mathbf{z}_i} \\
&= \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) - \frac{\text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2} \end{bmatrix} \right) \mathbf{z}_i \mathbf{z}_i^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}} \\
&= \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) - \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} \end{bmatrix} \right) \mathbf{z}_i \mathbf{z}_i^\top
\end{aligned}
$$

We also have:

$$
[\Sigma^{-1}_{y_i}]_{j,j'} = - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)}
\hspace{10mm}
\text{and}
\hspace{10mm}
[\Sigma^{-1}_{y_i}]_{j,j} = \frac{1}{\sigma^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j}^2}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)}
$$

{% endtab %}

{% endtabs %}


#### Derivatives
We first find the derivative with respect to $\tau^2$:

{% tabs deriv-tau-2 %}

{% tab deriv-tau-2 equation %}

$$
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} = 
 - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{j = 1}^n \left[ \frac{\mathbf{z}_{i,j}^2}{\sigma^2_{i,j}} - \sum_{h = 1}^n \frac{\tau^2}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)} \mathbf{z}_{i,j}^2 \mathbf{z}^2_{i,h}\right]  - 
\sum_{j = 1}^n \sum_{h = 1}^n (\mathbf{y}_{i,j} - \alpha_i)(\mathbf{y}_{i,h} - \alpha_i) \frac{\mathbf{z}_{i,j} \mathbf{z}_{i,h}}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2}
\right]
$$

{% endtab %}

{% tab deriv-tau-2 proof %}

First, the log determinant term:

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \log (\rvert \Sigma_{y_i} \rvert) \right]
&= \text{tr}\left[ \Sigma^{-1}_{y_i} \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i} \right] \right] \\
&= \text{tr} \left[ \left( \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) - \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) \mathbf{z}_i \mathbf{z}_i^\top \right) \mathbf{z}_i \mathbf{z}_i^\top \right] \\
&= \text{tr} \left[ \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) \mathbf{z}_{i} \mathbf{z}_i^\top \right] - \text{tr} \left[ \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) \mathbf{z}_i \mathbf{z}_i^\top \mathbf{z}_i \mathbf{z}_i^\top \right] \\
&= \sum_{j = 1}^n \sum_{k = 1}^n \left( \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) \right)_{j,k} \left(\mathbf{z}_i \mathbf{z}_{i}^\top \right)_{k,j} - \sum_{j = 1}^n \sum_{k = 1}^n \left( \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) \mathbf{z}_i \mathbf{z}_i^\top \right)_{j,k} \left( \mathbf{z}_i \mathbf{z}_i^\top \right)_{k,j} \\
&= \sum_{j = 1}^n  \frac{\mathbf{z}_{i,j}^2}{\sigma^2_{i,j}} - \sum_{j = 1}^n \sum_{k = 1}^n \frac{\tau^2}{(\sigma^2_{i,j})^2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)} \mathbf{z}_{i,j} \mathbf{z}_{i,k} \mathbf{z}_{i,k} \mathbf{z}_{i,j} \\
&= \sum_{j = 1}^n \left[ \frac{\mathbf{z}_{i,j}^2}{\sigma^2_{i,j}} - \sum_{k = 1}^n \frac{\tau^2}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)} \mathbf{z}_{i,j}^2 \mathbf{z}^2_{i,k}\right]
\end{aligned}
$$

Next, the quadratic term. We will compute the derivative of $\Sigma^{-1}_{y_i}$ with respect to $\tau^2$ element-wise:

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2}\left[ [\Sigma^{-1}_{y_i}]_{j,j'}\right] 
&= \frac{\partial}{\partial \tau^2} \left[ - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)} \right] \\
&= -\frac{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)(\mathbf{z}_{i,j}\mathbf{z}_{i,j'}) - (\tau^2 \mathbf{z}_{i,j}\mathbf{z}_{i,j'})\left( (\sigma_{i,j}^2)^2\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)}{(\sigma_{i,j}^2)^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
&= -\frac{\mathbf{z}_{i,j}\mathbf{z}_{i,j'} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}- \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
&= -\frac{\mathbf{z}_{i,j}\mathbf{z}_{i,j'}}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
\frac{\partial}{\partial \tau^2}\left[ [\Sigma^{-1}_{y_i}]_{j,j}\right] 
&= \frac{\partial}{\partial \tau^2} \left[ \frac{1}{\sigma^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j}^2}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)} \right] \\
&= - \frac{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)\mathbf{z}_{i,j}^2 - \tau^2 \mathbf{z}_{i,j}^2 (\sigma_{i,j}^2)^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}}{(\sigma_{i,j}^2)^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
&= - \frac{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)\mathbf{z}_{i,j}^2 - \tau^2 \mathbf{z}_{i,j}^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
&= - \frac{\mathbf{z}_{i,j}^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}- \tau^2\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
&= - \frac{\mathbf{z}_{i,j}^2}{(\sigma_{i,j}^2)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \\
\end{aligned}
$$

In matrix notation, we have:

$$
\frac{\partial}{\partial \tau^2}\left[ \Sigma^{-1}_{y_i} \right]
= \text{diag}\left(\left[ -\frac{1}{(\sigma^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2}, \dots, -\frac{1}{(\sigma^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top
$$

Then:

$$
\begin{aligned}
(\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \frac{\partial}{\partial \tau^2} \left[ \Sigma^{-1}_{y_i} \right](\mathbf{y}_i - \alpha_i \mathbf{1}_n) 
&= (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \left[ 
\text{diag}\left(\left[ -\frac{1}{(\sigma^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2}, \dots, -\frac{1}{(\sigma^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top
\right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \\
&= \sum_{j = 1}^n \sum_{h = 1}^n (\mathbf{y}_{i,j} - \alpha_i)(\mathbf{y}_{i,h} - \alpha_i)  \left[ \text{diag}\left(\left[ -\frac{1}{(\sigma^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2}, \dots, -\frac{1}{(\sigma^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right] \right) \mathbf{z}_i \mathbf{z}_i^\top
\right]_{j,k} \\
&= - \sum_{j = 1}^n \sum_{h = 1}^n (\mathbf{y}_{i,j} - \alpha_i)(\mathbf{y}_{i,h} - \alpha_i) \frac{\mathbf{z}_{i,j} \mathbf{z}_{i,h}}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2}
\end{aligned}
$$

And thus:

$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \ell(\theta; \mathbf{y}) \right]
&= \frac{\partial}{\partial \tau^2} \left[ \sum_{i = 1}^k \left[ -\frac{n}{2} \log(2 \pi) - \frac{1}{2}\log(\rvert \Sigma_{y_i} \rvert) - \frac{1}{2} (\mathbf{y}_i -  \alpha_i \mathbf{1}_n)^\top \Sigma_{y_i}^{-1} (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right] \right] \\
&= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{\partial}{\partial \tau^2} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right] + (\mathbf{y}_i - \alpha_i \mathbf{1}_n)^\top \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i}^{-1} \right] (\mathbf{y}_i - \alpha_i \mathbf{1}_n) \right]\\
&= - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{j = 1}^n \left[ \frac{\mathbf{z}_{i,j}^2}{\sigma^2_{i,j}} - \sum_{k = 1}^n \frac{\tau^2}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)} \mathbf{z}_{i,j}^2 \mathbf{z}^2_{i,k}\right]  - 
\sum_{j = 1}^n \sum_{k = 1}^n (\mathbf{y}_{i,j} - \alpha_i)(\mathbf{y}_{i,k} - \alpha_i) \frac{\mathbf{z}_{i,j} \mathbf{z}_{i,k}}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2}
\right]
\end{aligned}
$$

{% endtab %}

{% endtabs %}

Next, we find the derivative with respect to $\phi$:


{% tabs deriv-phi-1 %}
{% tab deriv-phi-1 equation %}
Coming soon...
{% endtab %}
{% tab deriv-phi-1 proof %}
First, let's find the derivative of the diagonal components of $\Sigma_{y_i}$:

$$
\begin{aligned}
\frac{\partial}{\partial \phi} \left[ (\Sigma_{y_i})_{j,j} \right]
&= \frac{\partial}{\partial \phi} \left[ \sigma^2_{i,j} + \tau^2 \mathbf{z}_{i,j}^2 \right] \\
&= \frac{\partial}{\partial \phi} \left[ \frac{1}{\tilde{\mu}_{i,j}} \left(1 + \frac{1}{\phi} \right) + \tau^2 \mathbf{z}_{i,j}^2  \right] \\
&= \frac{1}{\tilde{\mu}_{i,j}}  \frac{\partial}{\partial \phi} \left[ \frac{1}{\phi} \right] \\
&=  -\frac{1}{\tilde{\mu}_{i,j} \phi^2 } 
\end{aligned}
$$

This implies that $$\frac{\partial}{\partial \phi}\left[ \Sigma_{y_i} \right] = \text{diag}\left( \begin{bmatrix} -\frac{1}{\tilde{\mu}_{i,1} \phi^2 } & \dots & -\frac{1}{\tilde{\mu}_{i,n} \phi^2 } \end{bmatrix} \right)$$.

$$
\begin{aligned}
\frac{\partial}{\partial \phi} \left[ \log(\rvert \Sigma_{y_i} \rvert) \right]
&= \text{tr} \left[ \Sigma_{y_i}^{-1} \frac{\partial}{\partial \tau^2} \left[ \Sigma_{y_i} \right] \right] \\
&= \text{tr} \left[ \left( \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i, 1}} & \dots & \frac{1}{\sigma^2_{i,n}} \end{bmatrix} \right) - \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) \mathbf{z}_i \mathbf{z}_i^\top \right) \text{diag}\left( \begin{bmatrix} -\frac{1}{\tilde{\mu}_{i,1} \phi^2 } & \dots & -\frac{1}{\tilde{\mu}_{i,n} \phi^2 } \end{bmatrix} \right) \right] \\
&= \text{tr} \left[ \text{diag}\left(\begin{bmatrix} - \frac{1}{\tilde{\mu}_{i,1} \sigma^2_{i,1} \phi^2} & \dots & - \frac{1}{\tilde{\mu}_{i,n} \sigma^2_{i,n} \phi^2} \end{bmatrix}\right) \right] - 
\text{tr} \left[ \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) (\mathbf{z}_i \mathbf{z}_i^\top)  \text{diag}\left( \begin{bmatrix} -\frac{1}{\tilde{\mu}_{i,1} \phi^2 } & \dots & -\frac{1}{\tilde{\mu}_{i,n} \phi^2 } \end{bmatrix} \right) \right]  \\
&= \sum_{j = 1}^n \sum_{h = 1}^n \left[ \text{diag}\left(\begin{bmatrix} - \frac{1}{\tilde{\mu}_{i,1} \sigma^2_{i,1} \phi^2} & \dots & - \frac{1}{\tilde{\mu}_{i,n} \sigma^2_{i,n} \phi^2} \end{bmatrix}\right)\right]_{j,h} - \text{tr}\left[ \mathbf{z}_i^\top  \text{diag}\left( \begin{bmatrix} -\frac{1}{\tilde{\mu}_{i,1} \phi^2 } & \dots & -\frac{1}{\tilde{\mu}_{i,n} \phi^2 } \end{bmatrix} \right) \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix} \right) \mathbf{z}_i   \right] \\
&= \sum_{j = 1}^n -\frac{1}{\tilde{\mu}_{i,j} \sigma^2_{i,j} \phi^2} - \mathbf{z}_i^\top  \text{diag}\left( \begin{bmatrix} -\frac{1}{\tilde{\mu}_{i,1} \phi^2 } & \dots & -\frac{1}{\tilde{\mu}_{i,n} \phi^2 } \end{bmatrix} \right) \text{diag}\left( \begin{bmatrix} \frac{\tau^2}{(\sigma_{i,1}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots & \frac{\tau^2}{(\sigma_{i,n}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix} \right) \mathbf{z}_i \\
&= \sum_{j = 1}^n -\frac{1}{\tilde{\mu}_{i,j} \sigma^2_{i,j} \phi^2} - \mathbf{z}_{i}^\top \text{diag} \left( \begin{bmatrix} - \frac{\tau^2}{\tilde{\mu}_{i,1} \phi^2 (\sigma^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)} & \dots & - \frac{\tau^2}{\tilde{\mu}_{i,n} \phi^2 (\sigma^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}  \end{bmatrix} \right) \mathbf{z}_i \\
&= - \sum_{j = 1}^n \frac{1}{\tilde{\mu}_{i,j} \left(\frac{1}{\tilde{\mu}_{i,j}} \left(1 + \frac{1}{\phi}\right) \right) \phi^2} - \sum_{j = 1}^n \sum_{h = 1}^n \mathbf{z}_{i,j} \mathbf{z}_{i,h} \left[ \text{diag} \left( \begin{bmatrix} - \frac{\tau^2}{\tilde{\mu}_{i,1} \phi^2 (\sigma^2_{i,1})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)} & \dots & - \frac{\tau^2}{\tilde{\mu}_{i,n} \phi^2 (\sigma^2_{i,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}  \end{bmatrix} \right) \right]_{j,h} \\
&= - \sum_{j = 1}^n \frac{1}{\phi^2 + \phi} - \sum_{j = 1}^n \mathbf{z}_{i,j}^2 \left(- \frac{\tau^2}{\tilde{\mu}_{i,j} \phi^2 (\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}  \right) \\
&= - \sum_{j = 1}^n \left[ \frac{1}{\phi^2 + \phi} - \frac{\tau^2 \mathbf{z}_{i,j}^2 }{\tilde{\mu}_{i,n} \phi^2 \left( \frac{1}{\tilde{\mu}_{i,j}} \left(1 + \frac{1}{\phi}\right) \right)^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right) }\right] \\
&= - \sum_{j = 1}^n \left[ \frac{1}{\phi^2 + \phi} - \frac{\tau^2 \mathbf{z}_{i,j}^2 }{ \sigma^2_{i,j} (\phi^2 + \phi) \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right) }\right]  \\ 
&= -\sum_{j = 1}^n \frac{\sigma^2_{i,j}  \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right) - \tau^2 \mathbf{z}_{i,j}^2 }{\sigma^2_{i,j} (\phi^2 + \phi) \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right) } 
&= \sum_{j = 1}^n \frac{\tau^2 \mathbf{z}_{i,j}^2  -  \left( \sigma^2_{i,j} + \sigma^2_{i,j} \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}{\sigma^2_{i,j} (\phi^2 + \phi) \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right) } 
\end{aligned}
$$

Next, the quadratic term. This requires the following:

$$
\begin{aligned}
\frac{\partial}{\partial \phi} \left[ [ \Sigma^{-1}_{y_i}]_{j,j'} \right]
&= \frac{\partial}{\partial \phi} \left[ - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}\right] \\
&= - \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \frac{\partial}{\partial \phi} \left[ (\sigma^2_{i,j})^{-2} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1} \right] \\
&= - \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \left[ -2(\sigma^2_{i,j})^{-3} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1} \frac{\partial}{\partial \phi} \left[ \sigma^2_{i,j} \right] - \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1}  (\sigma^2_{i,j})^{-2} \frac{\partial}{\partial \phi} \left[ 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right] \right] \\
&= - \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \left[ -2(\sigma^2_{i,j})^{-3} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1} \frac{\partial}{\partial \phi} \left[ \frac{1}{\tilde{\mu}_{i,j}}\left(1 + \frac{1}{\phi} \right) \right] - \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1}  (\sigma^2_{i,j})^{-2}  \tau^2 \sum_{l = 1}^n \mathbf{z}_{i,l}^2 \frac{\partial}{\partial \phi} \left[ \left( \sigma^2_{i,l} \right)^{-1}\right] \right] \\
&= - \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \left[ 2(\sigma^2_{i,j})^{-3} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1} \tilde{\mu}^{-1}_{i,j} \phi^{-2} - \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1}  (\sigma^2_{i,j})^{-2}  \tau^2 \sum_{l = 1}^n - \mathbf{z}^2_{i,l} (\sigma^2_{i,l})^{-2} \frac{\partial}{\partial \phi} \left[  \sigma^2_{i,l} \right] \right] \\
&= - \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \left[ 2(\sigma^2_{i,j})^{-3} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1} \tilde{\mu}^{-1}_{i,j} \phi^{-2} + \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1}  (\sigma^2_{i,j})^{-2}  \tau^2 \sum_{l = 1}^n \mathbf{z}^2_{i,l} (\sigma^2_{i,l})^{-2} \tilde{\mu}_{i,l}^{-1} \phi^{-2} \right] \\
&= - \tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^{-1} (\sigma^2_{i,j})^{-2}  \phi^{-2} \left[ 2(\sigma^2_{i,j})^{-1} \tilde{\mu}^{-1}_{i,j} +  \tau^2 \sum_{l = 1}^n \mathbf{z}^2_{i,l} (\sigma^2_{i,l})^{-2} \tilde{\mu}_{i,l}^{-1} \right] \\
&= - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^{2}  \phi^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right) }  \left[ \frac{2}{1 + \frac{1}{\phi}} + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}^2_{i,l}}{ \sigma^2_{i,l} \left(1 + \frac{1}{\phi} \right)} \right] \\
&= - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^{2}  \phi^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right) \left(1 + \frac{1}{\phi} \right)}  \left[ 2 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}^2_{i,l}}{ \sigma^2_{i,l} } \right] \\
\end{aligned}
$$

And:

$$
\begin{aligned}
\frac{\partial}{\partial \phi} \left[ [ \Sigma^{-1}_{y_i}]_{j,j} \right]
&= \frac{\partial}{\partial \phi} \left[  \frac{1}{\sigma^2_{i,j}} - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}\right] \\
&= -\frac{1}{(\sigma^2_{i,j})^2 } \frac{\partial}{\partial\phi} \left[ \sigma^2_{i,j} \right] - \frac{\partial}{\partial \phi} \left[ \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}\right] \\
&= \frac{1}{(\sigma^2_{i,j})^2 \tilde{\mu}_{i,j} \phi^2} - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^{2}  \phi^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right) \left(1 + \frac{1}{\phi} \right)}  \left[ 2 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}^2_{i,l}}{ \sigma^2_{i,l} } \right] \\
&=  \frac{1}{\sigma^2_{i,j} \phi^2 \left(1 + \frac{1}{\phi} \right)} - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^{2}  \phi^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right) \left(1 + \frac{1}{\phi} \right)}  \left[ 2 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}^2_{i,l}}{ \sigma^2_{i,l} } \right] \\
&= \frac{1}{\sigma^2_{i,j}\left( \phi^2  + \phi\right)} - \frac{\tau^2 \mathbf{z}_{i,j} \mathbf{z}_{i,j'}}{(\sigma^2_{i,j})^{2} \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right) \left( \phi^2 + \phi\right)}  \left[ 2 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}^2_{i,l}}{ \sigma^2_{i,l} } \right] 
\end{aligned}
$$

In matrix notation:

$$
\begin{aligned}
\frac{\partial}{\partial \phi} \left[ \Sigma^{-1}_{y_i} \right]
&=  \text{diag}\left( \begin{bmatrix} \frac{1}{\sigma^2_{i,1} (\phi^2 + \phi)} & \dots & \frac{1}{\sigma^2_{i,n} (\phi^2 + \phi)} \end{bmatrix} \right) - \text{diag} \left( \begin{bmatrix} \frac{ \tau^2 \left(2 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,1})^2 (\phi^2 + \phi) \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} & \dots &  \frac{ \tau^2 \left(2 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,n})^2 (\phi^2 + \phi) \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)}  \end{bmatrix}  \right) \mathbf{z}_i \mathbf{z}_i^\top
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
\left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{1,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{1,l'}}{\sigma^2_{1,l'}}  \\
\vdots \\
\left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{k,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{k,l'}}{\sigma^2_{k,l'}}  \\
0
\end{bmatrix}
$$

{% endtab %}

{% tab mle-2 proof %}

We only need to deal with $\alpha$, which we can again do component-wise. First notice that:

$$
\begin{aligned}
\frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n  \\
&= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left(  \text{diag}\left(\left[ \frac{1}{\sigma^2_{j, 1}}, \dots, \frac{1}{\sigma^2_{j,n}} \right]\right) - \frac{\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{j,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{j,n}^2)^2} \right] \right) \mathbf{z}_j \mathbf{z}_j^\top}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma^2_{j,l}}} \right) \mathbf{1}_n \\
&= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left(\begin{bmatrix} 
\frac{1}{\sigma^2_{j, 1}} \\
\vdots \\
\frac{1}{\sigma^2_{j,n}}
\end{bmatrix} - \left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right)\text{diag}\left( \left[ \frac{\tau^2}{(\sigma_{j,1}^2)^2}, \dots, \frac{\tau^2}{(\sigma_{j,n}^2)^2} \right] \right) \begin{bmatrix} 
\sum_{l = 1}^n \mathbf{z}_{j,1} \mathbf{z}_{j,l} \\
\vdots \\
\sum_{l = 1}^n \mathbf{z}_{j,n} \mathbf{z}_{j,l}
\end{bmatrix} \right) \\
&= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\sigma^2_{j,l}} \right)- 
\left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \begin{bmatrix} 
\frac{\tau^2}{(\sigma_{j,1}^2)^2} \sum_{l = 1}^n \mathbf{z}_{j,1} \mathbf{z}_{j,l} \\
\vdots \\
\frac{\tau^2}{(\sigma_{j,n}^2)^2}  \sum_{l = 1}^n \mathbf{z}_{j,n} \mathbf{z}_{j,l}
\end{bmatrix} \\
&= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\sigma^2_{j,l}} \right)- 
\left( \frac{1}{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \sum_{l' = 1}^n (\mathbf{y}_{j,l'} - \alpha_j) 
\frac{\tau^2}{(\sigma_{j,l'}^2)^2} \sum_{l = 1}^n \mathbf{z}_{j,l'} \mathbf{z}_{j,l}  \\
&= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\sigma^2_{j,l}} \right)- 
\left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \left( \sum_{l' = 1}^n 
\frac{\tau^2 (\mathbf{y}_{j,l'} - \alpha_j) \mathbf{z}_{j,l'} }{(\sigma_{j,l'}^2)^2}  \right)
\end{aligned}
$$

So then we solve for $\alpha_j$ in:

$$
\begin{aligned}
0 &= (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \\
\implies 0 &= \left( \sum_{l = 1}^n \frac{\mathbf{y}_{j,l} - \alpha_j}{\sigma^2_{j,l}} \right)- 
\left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \left( \sum_{l' = 1}^n 
\frac{\tau^2 (\mathbf{y}_{j,l'} - \alpha_j) \mathbf{z}_{j,l'} }{(\sigma_{j,l'}^2)^2}  \right)  \\
\implies 0 &= \sum_{l = 1}^n \frac{\mathbf{y}_{j,l}}{\sigma^2_{j,l}} - \alpha_j \sum_{l = 1}^n \frac{1}{\sigma^2_{j,l}} - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \left[ \sum_{l' = 1}^n \frac{\tau^2 \mathbf{y}_{j,l'} \mathbf{z}_{j,l'}}{(\sigma^2_{j,l'})^2} - \alpha_j \sum_{l' =1}^n \frac{\tau^2 \mathbf{z}_{j,l'}}{(\sigma^2_{j,l'})} \right] \\
\implies \alpha_j \left[ \sum_{l = 1}^n \frac{1}{\sigma^2_{j,l}} - \sum_{l = 1}^n \frac{\tau^2 \mathbf{z}_{j,l'}}{(\sigma^2_{j,l'})^2}\right] &= \sum_{l = 1}^n \frac{\mathbf{y}_{j,l}}{\sigma^2_{j,l}} - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \sum_{l' = 1}^n\frac{\tau^2 \mathbf{y}_{j,l'} \mathbf{z}_{j,l'}}{(\sigma^2_{j,l'})^2} \\
\implies \alpha_j \left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{j,l'}} \left[1 - \frac{\tau^2 \mathbf{z}_{j,l'}}{\sigma^2_{j,l'}} \right]\right) &= \sum_{l' = 1}^n  \frac{\mathbf{y}_{j,l'}}{\sigma^2_{j,l'}} \left[ 1 - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \frac{\tau^2 \mathbf{z}_{j,l'}}{\sigma^2_{j,l'}} \right] \\
\implies \alpha_j &= \left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{j,l'}} \left[1 - \frac{\tau^2 \mathbf{z}_{j,l'}}{\sigma^2_{j,l'}} \right]\right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{j,l'}}{\sigma^2_{j,l'}} \left[ 1 - \left( \frac{\sum_{l = 1}^n \mathbf{z}_{j,l} }{1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}^2}{\sigma_{j,l}^2}} \right) \frac{\tau^2 \mathbf{z}_{j,l'}}{\sigma^2_{j,l'}} \right]
\end{aligned}
$$

Under $H_0$, $\tau^2 = 0$, so we get:

$$
\hat{\alpha}_j = \left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{j,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{j,l'}}{\sigma^2_{j,l'}} 
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
\left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{1,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{1,l'}}{\sigma^2_{1,l'}}  \\
\vdots \\
\left(\sum_{l' = 1}^n \frac{1}{\sigma^2_{k,l'}} \right)^{-1} \sum_{l' = 1}^n  \frac{\mathbf{y}_{k,l'}}{\sigma^2_{k,l'}}  \\
- \frac{1}{2}\sum_{i = 1}^k \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma_{i,l}^2} + \sum_{l = 1}^n \sum_{l' = 1}^n \frac{\mathbf{z}_{i,l}\mathbf{z}_{i,l'}(\mathbf{y}_{i,l} - \hat{\alpha}_i)(\mathbf{y}_{i,l'} - \hat{\alpha}_i)}{\sigma^2_{i,l} \sigma^2_{i,l'}}\right]
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
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{1,l}}{\sigma^2_{1,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{1,l'}(\mathbf{y}_{1,l'} - \alpha_1)}{(\sigma^2_{1,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{1,l''}\right) \\
    \vdots \\
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{k,l}}{\sigma^2_{k,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{k,l'}(\mathbf{y}_{k,l'} - \alpha_k)}{(\sigma^2_{k,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{k,l''}\right) \\
    \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\sigma^2_{i,a}} \right)^2 \right]
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
    &= \frac{\partial}{\partial \tau^2} \left[ - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} - \sum_{l = 1}^n \sum_{j = 1}^n \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\sigma_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} - \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\sigma^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right] \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k  \left[\sum_{l = 1}^n \sum_{j = 1}^n \frac{\partial}{\partial \tau^2} \left[  \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\sigma_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} \right] + \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\partial}{\partial \tau^2} \left[ \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\sigma^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right] \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[ 
        \sum_{l = 1}^n \sum_{j = 1}^n \left(\frac{\mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2 (\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right) - \tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2 (\sigma^2_{i,l})^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} }{\left((\sigma^2_{i,l})^2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)\right)^2}\right) + \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{-2\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\sigma^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma_{i,l}}\right) \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma_{i,l}^2}}{\left((\sigma^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)^2\right)^2}\right)
    \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\sigma^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,a})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)^3}\right) \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{(\sigma^2_{i,l})^2} \right) \left( \sum_{j = 1}^n \mathbf{z}_{i,j}^2 \right)
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{(\sigma^2_{i,a})^2} \right) \left( \sum_{b = 1}^n  (\mathbf{y}_{i,b} - \alpha_i) \mathbf{z}_{i,b} \right) \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{(\sigma^2_{i,l})^2} \right) \left( \sum_{l = 1}^n \mathbf{z}_{i,l}^2 \right)
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\sigma^2_{i,a}} \right) \left( \sum_{a = 1}^n  (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\sigma^2_{i,a}} \right) \right] \\
    &= \frac{1}{2} \sum_{i = 1}^k \left[  \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \left( \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2
    - \frac{2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^3} \left( \sum_{a = 1}^n (\mathbf{y}_{i,a} - \alpha_i) \frac{\mathbf{z}_{i,a}}{\sigma^2_{i,a}} \right)^2 \right]
    
\end{aligned}
$$

Next, we take the derivative (with respect to $\tau^2$) of the derivative with respect to $\alpha_j$:

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2  \partial \alpha_j }
    &= \frac{\partial}{\partial \tau^2} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha_j} \right] \\
    &= \frac{\partial}{\partial \tau^2} \left[ (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \Sigma_{y_j}^{-1} \mathbf{1}_n \right] \\
    &= \begin{bmatrix}
        (\mathbf{y}_j - \alpha_j \mathbf{1}_n)^\top \left[ \text{diag}\left( \begin{bmatrix} -\frac{1}{(\sigma^2_{j,1})^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} & \dots & -\frac{1}{(\sigma^2_{j,n})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} \end{bmatrix} \right) \mathbf{z}_{j} \mathbf{z}_j^\top  \right] \mathbf{1}_n
    \end{bmatrix} \\
    &= \begin{bmatrix}
        -\frac{\mathbf{y}_{j, 1} - \alpha_j}{(\sigma^2_{j,1})^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} & \dots & -\frac{\mathbf{y}_{j,n} - \alpha_j}{(\sigma^2_{j,n})^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2}
    \end{bmatrix} \mathbf{z}_j \mathbf{z}_j^\top \mathbf{1}_n \\
    &= \left(- \sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'}(\mathbf{y}_{j,l'} - \alpha_j)}{(\sigma^2_{j,l'})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 }\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right) \\
    &= - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'}(\mathbf{y}_{j,l'} - \alpha_j)}{(\sigma^2_{j,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right)
\end{aligned}
$$

Putting the two together into a matrix yield:

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \theta }
&= \frac{\partial}{\partial \tau^2} \left[ \begin{bmatrix} \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha} \\ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \end{bmatrix} \right] \\
&= \begin{bmatrix}
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{1,l}}{\sigma^2_{1,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{1,l'}(\mathbf{y}_{1,l'} - \alpha_1)}{(\sigma^2_{1,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{1,l''}\right) \\
    \vdots \\
    - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{k,l}}{\sigma^2_{k,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{k,l'}(\mathbf{y}_{k,l'} - \alpha_k)}{(\sigma^2_{k,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{k,l''}\right) \\
    \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\sigma^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)^3}\right) \right]
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
-\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\sigma_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right)
\end{bmatrix}
$$

{% endtab %}

{% tab deriv-theta-alpha-2 proof %}

First, we find the gradient (with respect to $\alpha$) of the derivative with respect to $\tau^2$. This should be equivalent to the corresponding components of Eq. \eqref{eq:info-tau-tau}.

$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \alpha_j \partial \tau^2 } 
&= \frac{\partial}{\partial \alpha_j} \left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2} \right] \\
&= \frac{\partial}{\partial \alpha_j} \left[  - \frac{1}{2} \sum_{i = 1}^k  \left[ \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} - \sum_{l = 1}^n \sum_{j = 1}^n \frac{\tau^2 \mathbf{z}_{i,l}^2 \mathbf{z}_{i,j}^2}{(\sigma_{i,l}^2)^2\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}}\right)} - \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\sigma^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right]  \right] \\
&= \frac{\partial}{\partial \alpha_j} \left[ \frac{1}{2} \sum_{i = 1}^k \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \frac{\mathbf{z}_{i,a} \mathbf{z}_{i,b}}{(\sigma^2_{i,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}}\right)^2} \right] \\
&= \frac{1}{2} \sum_{a = 1}^n \sum_{b = 1}^n \frac{\partial}{\partial \alpha_j} \left[ (\mathbf{y}_{j,a} - \alpha_j)(\mathbf{y}_{j,b} - \alpha_j) \frac{\mathbf{z}_{j,a} \mathbf{z}_{j,b}}{(\sigma^2_{j,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} \right] \\
&= \frac{1}{2} \sum_{a = 1}^n \sum_{b = 1}^n \frac{\mathbf{z}_{j,a} \mathbf{z}_{j,b}}{(\sigma^2_{j,a})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} \frac{\partial}{\partial \alpha_j} \left[ \mathbf{y}_{j,a} \mathbf{y}_{j,b} - \alpha_j \mathbf{y}_{j,a} - \alpha \mathbf{y}_{j,b} + \alpha_j^2 \right] \\
&= -\frac{1}{2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \sum_{a = 1}^n \sum_{b = 1}^n \left[ \frac{\mathbf{z}_{j,a}\mathbf{z}_{j,b}(\mathbf{y}_{j, a} - \alpha_j)}{(\sigma_{j,a}^2)^2} + \frac{\mathbf{z}_{j,a}\mathbf{z}_{j,b}(\mathbf{y}_{j,b} - \alpha_j)}{(\sigma_{j,a}^2)^2} \right] \\
&= -\frac{1}{2 \left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \left[ \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\sigma_{j,a}^2)^2} \sum_{b = 1}^n \mathbf{z}_{j,b} + \sum_{b = 1}^n  \frac{\mathbf{z}_{j,b}(\mathbf{y}_{j,b} - \alpha_j)}{(\sigma_{j,b}^2)^2}  \sum_{a = 1}^n \mathbf{z}_{j,a}\right] \\
&= -\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\sigma_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right)
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
-\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\sigma_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right)
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
\frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \right]
\end{bmatrix}
$$

{% endtab %}

{% tab info-tau-tau-2 proof %}

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial (\tau^2)^2 }\right]
&= \mathbb{E}\left[ \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\sigma^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)^3}\right) \right] \right] \\
&= \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n \mathbb{E}\left[ (\mathbf{y}_{i,a} - \alpha_i)(\mathbf{y}_{i,b} - \alpha_i)\right] \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\sigma^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)^3}\right) \right] \\
&=  \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} - 2 \sum_{a = 1}^n \sum_{b = 1}^n \text{Cov}\left(\mathbf{y}_{i,a}, \mathbf{y}_{i,b} \right) \left(\frac{\mathbf{z}_{i,a}\mathbf{z}_{i,b}(\sigma^2_{i,a})^2 \left(\sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)}{(\sigma^2_{i,l})^4 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}}{\sigma^2_{i,l}} \right)^3}\right) \right] \\
&=  \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \right] & \left(\text{obs. ind. under } H_0 \right) \\
\mathbb{E}\left[ \frac{\partial^2 \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \alpha_j}\right] 
&= \mathbb{E}\left[ - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'}(\mathbf{y}_{j,l'} - \alpha_j)}{(\sigma^2_{j,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right) \right] \\
&= - \frac{1}{\left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2} \left(\sum_{l' = 1}^n \frac{\mathbf{z}_{j,l'} \mathbb{E}\left[ (\mathbf{y}_{j,l'} - \alpha_j)\right] }{(\sigma^2_{j,l'})^2}\right)\left(\sum_{l'' = 1}^n \mathbf{z}_{j,l''}\right)  \\
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
&= \mathbb{E}\left[ -\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a}(\mathbf{y}_{j, a} - \alpha_j)}{(\sigma_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right) \right] \\
&= -\frac{1}{\left( 1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{j,l}}{\sigma^2_{j,l}}\right)^2 } \sum_{a = 1}^n \frac{\mathbf{z}_{j,a} \mathbb{E}\left[\mathbf{y}_{j, a} - \alpha_j\right] }{(\sigma_{j,a}^2)^2} \left(\sum_{b = 1}^n \mathbf{z}_{j,b} \right) \\
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
- \mathbf{1}_n^\top\text{diag}\left(\left[ \frac{1}{\sigma^2_{1,1}}, \dots, \frac{1}{\sigma^2_{1,n}}\right]\right) \mathbf{1}_n & 0 & \dots & 0 \\
0 & - \mathbf{1}_n^\top \text{diag}\left(\left[ \frac{1}{\sigma^2_{2,1}}, \dots, \frac{1}{\sigma^2_{2,n}}\right]\right) \mathbf{1}_n & \dots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \dots & \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2} \right]
\end{bmatrix}
$$

{% endtab %}

{% tab info-2 proof %}

Note that:

$$
\Sigma_{y_i}^{-1} \bigg\rvert_{\theta = \hat{\theta}} = \text{diag}\left(\left[ \frac{1}{\sigma^2_{i,1}}, \dots, \frac{1}{\sigma^2_{i,n}}\right]\right)
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
    \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2 \left(1 + \tau^2 \sum_{l = 1}^n \frac{\mathbf{z}_{i,l}^2}{\sigma^2_{i,l}} \right)^2} \right] 
    \end{bmatrix}
\right] \right\rvert_{\theta = \hat{\theta}} \\
&= \begin{bmatrix}
    0 \\
    \vdots \\
    0 \\
    \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2} \right]
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
    - \mathbf{1}_n^\top \text{diag}\left(\left[ \frac{1}{\sigma^2_{j,1}}, \dots, \frac{1}{\sigma^2_{j,n}}\right]\right) \mathbf{1}_n  \\
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
- \mathbf{1}_n^\top\text{diag}\left(\left[ \frac{1}{\sigma^2_{1,1}}, \dots, \frac{1}{\sigma^2_{1,n}}\right]\right) \mathbf{1}_n & 0 & \dots & 0 \\
0 & - \mathbf{1}_n^\top \text{diag}\left(\left[ \frac{1}{\sigma^2_{2,1}}, \dots, \frac{1}{\sigma^2_{2,n}}\right]\right) \mathbf{1}_n & \dots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \dots & \frac{1}{2} \sum_{i = 1}^k \left[ \sum_{l = 1}^n \sum_{j = 1}^n \frac{\mathbf{z}_{i,l}^2\mathbf{z}_{i,j}^2}{(\sigma^2_{i,l})^2} \right]
\end{bmatrix}
\end{aligned}
$$

{% endtab %}

{% endtabs %}


---

<!-- 
## References

[^fn-fitzmaurice]: Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). Applied longitudinal analysis (Second edition). Wiley. -->
