---
layout: distill
title: Computational Methods For Mixed Models
description: lme4
date: 2026-03-19
tabs: true
tags: methods computation glmm lmm
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Mixed Effects Models
      subsections:
        - name: Spherical Random Effects
        - name: Density and Likelihood
    - name: Model Fitting
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: stats-ml.bib
---

I thought it would be wise to return to some of the basics and study how (generalized) linear mixed models are fit in standard packages. This way, I will have a sturdy foundation upon which to build my understanding of the more "non-standard" cases. In this post, I'll mostly be focusing on the theory behind the `lme4` package in `R`, which can be used to fit (generalized) linear mixed models of all different types.<d-cite key=bates2026></d-cite>

---

## Mixed Effects Models
We have some $n$-dimensional response vector, $\boldsymbol{\mathcal{Y}}$, and a random effects vector, $\boldsymbol{\mathcal{B}} \in \mathbb{R}^q$. An observation of $\boldsymbol{\mathcal{Y}}$ is denoted by $\mathbf{y}$. For simplicity, we assume that the unconditional distribution of $\boldsymbol{\mathbf{B}}$ is (multivariate) Gaussian:

$$
\boldsymbol{\mathcal{B}} \sim \mathcal{N}\left(\mathbf{0}_q, \sigma^2 \boldsymbol{\Lambda}(\boldsymbol{\theta}) \boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \right)
$$

though this normality assumption has been relaxed in other literature. $\boldsymbol{\Lambda}(\boldsymbol{\theta})$ is the $q \times q$ left factor of the covariance matrix of the random effects that is parametrized by some $m$-dimensional vector, $\boldsymbol{\theta}$, where $m < < q$. $\boldsymbol{\Lambda}(\boldsymbol{\theta})$ is invertible when $\boldsymbol{\theta}$ is not on the boundary (e.g. the coordinates of $\boldsymbol{\theta}$ are constrained to some particular intervals). The $\sigma^2$ is the <strong>common scale parameter</strong> that is frequently set to $1$. 

We do not observe $\boldsymbol{\mathcal{B}}$, but our model involves it through the conditional distribution of $\boldsymbol{\mathcal{Y}}$. We assume that:

$$
\begin{equation}
\label{eq:conditional-mean}
\boldsymbol{\mu}_{\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{B}}}(\boldsymbol{\beta}) = \mathbb{E}\left[\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{B}} = \boldsymbol{\beta} \right] = g\left(\mathbf{X} \boldsymbol{\alpha} + \mathbf{Z} \boldsymbol{\beta}\right)
\end{equation}
$$

where $\mathbf{X}$ is an $n \times p$ design matrix of fixed effects covariates, $\mathbf{Z}$ is an $n \times q$ design matrix of random effects covariates, and $\boldsymbol{\alpha}$ is a $p$-dimensional vector fo fixed effects. 

Eq. \eqref{eq:conditional-mean} gives the <i>conditional mean</i> of our response given the random effects. The <strong>link function</strong>, $g(\cdot)$, (applied element-wise) related the conditional mean to the <strong>linear predictor</strong>:

$$
\boldsymbol{\eta}(\boldsymbol{\beta}, \boldsymbol{\theta}, \boldsymbol{\alpha}) = \mathbf{X} \boldsymbol{\alpha} + \mathbf{Z} \boldsymbol{\beta}
$$

We make the additional assumptions that the components of $\boldsymbol{\mathcal{Y}}$ are independent <i>conditional on $\boldsymbol{\mathcal{B}}$</i> and that they all have the same distributional forms that are fully determined by Eq. \eqref{eq:conditional-mean} and $\sigma$.

<div class="example">
<strong>Example (Linear Mixed Model).</strong>
<br>
A simple example that is very frequently used is the <i>linear mixed model</i>. The conditional distribution of the responses given the random effects is assumed to be multivariate Gaussian:

$$
(\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{B}} = \boldsymbol{\beta}) \sim \mathcal{N}\left( \mathbf{X} \boldsymbol{\alpha} + \mathbf{Z} \boldsymbol{\beta}, \sigma^1 \mathbb{I}_{n \times n} \right)
$$

In this case, $g(\cdot)$ is the <i>identity link function</i>. 
</div>

### Spherical Random Effects
Let us define a new random variable:

$$
\boldsymbol{\mathcal{U}} \sim \mathcal{N}(\mathbf{0}_q, \sigma^2 \mathbb{I}_{q \times q})
$$

Bates terms $\boldsymbol{\mathcal{U}}$ a <strong>"spherical" Gaussian</strong> random variable because its distribution is symmetrical about its mean (i.e. constant diagonal covariance matrix). Thus, the probability contours will be <i>spheres</i> centered at $\mathbf{0}_q$. 

We can rewrite $\boldsymbol{\mathcal{B}}$ as:

$$
\begin{aligned}
\mathbb{E}\left[\boldsymbol{\Lambda}(\boldsymbol{\theta}) \boldsymbol{\mathcal{U}}\right] &= \mathbf{0}_q \\
\text{Cov}(\boldsymbol{\Lambda}(\boldsymbol{\theta}) \boldsymbol{\mathcal{U}}) &= \boldsymbol{\Lambda}(\boldsymbol{\theta}) \text{Cov}(\boldsymbol{\mathcal{U}}) \boldsymbol{\Lambda}^\top(\boldsymbol{\theta})
= \sigma^2 \boldsymbol{\Lambda}(\boldsymbol{\theta})  \boldsymbol{\Lambda}^\top(\boldsymbol{\theta})  \\
\implies \boldsymbol{\mathcal{B}} &= \boldsymbol{\Lambda}(\boldsymbol{\theta}) \boldsymbol{\mathcal{U}} \\
\end{aligned}
$$

If we let $\mathbf{u}$ be a particular realization of $\boldsymbol{\mathcal{U}}$, then we can rewrite the linear predictor as a function of $\mathbf{u}$:

$$
\boldsymbol{\eta}(\mathbf{u}) := \boldsymbol{\eta}(\mathbf{u}, \boldsymbol{\theta}, \boldsymbol{\alpha}) = \mathbf{X} \boldsymbol{\alpha} + \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta})\mathbf{u}
$$

### Density and Likelihood
Our goal is statistical inference on the model parameters; however, we never obtain observations of $\boldsymbol{\mathcal{B}}$ or $\boldsymbol{\mathcal{U}}$. Let $f_{X \rvert Y}(X \rvert Y)$ denote the conditional probability density function of some random variable $X$ given $Y$.

If we wish to find the <i>maximum likeihood estimates</i> (MLEs) of $\boldsymbol{\alpha}$, $\boldsymbol{\theta}$, and $\sigma$, then we need to maximize:

$$
\begin{equation}
\label{eq:marg-lik}
\mathcal{L}(\boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{y}) = \int_{\mathbb{R}^q} f_{\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{U}}}(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{u}) f_{\boldsymbol{\mathcal{U}}}(\mathbf{u}; \sigma) d \mathbf{u}
\end{equation}
$$

which is the <i>marginal</i> likelihood of $\boldsymbol{\alpha}$, $\boldsymbol{\theta}$, and $\sigma$ given the data, $\mathbf{y}$. This is equal to the marginal density of $\boldsymbol{\mathcal{Y}}$, $f_{\boldsymbol{\mathcal{Y}}}(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma)$.

The integrand in Eq. \eqref{eq:marg-lik} is equal to the <i>joint</i> density of $\boldsymbol{\alpha}$, $\boldsymbol{\theta}$, $\sigma$, and $\boldsymbol{u}$, and the joint density is proportional to the posterior density:

$$
\begin{aligned}
f_{\boldsymbol{\mathcal{Y}}, \boldsymbol{\mathcal{U}}}(\mathbf{y}, \mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma) &=  f_{\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{U}}}(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{u}) f_{\boldsymbol{\mathcal{U}}}(\mathbf{u}; \sigma) \\
\implies f_{\boldsymbol{\mathcal{U}} \rvert \boldsymbol{\mathcal{Y}}}(\mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{y}) &\propto  f_{\boldsymbol{\mathcal{Y}}, \boldsymbol{\mathcal{U}}}(\mathbf{y}, \mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma) \\
=: h(\mathbf{y}, \mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma)
\end{aligned}
$$

<aside><p>We just need to divide by the marginal density of $\boldsymbol{\mathcal{Y}}$ to obtain equality via Bayes' Theorem.</p></aside>

One thing to note is that Bates emphasizes that the conditional distribution of $(\boldsymbol{\mathcal{U}} \rvert \boldsymbol{\mathcal{Y}} = \mathbf{y})$ is <i>always</i> continuous with density $f_{\boldsymbol{\mathcal{U}} \rvert \boldsymbol{\mathcal{Y}}}(\mathbf{u};\boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma\rvert \mathbf{y})$. The same goes for $(\boldsymbol{\mathcal{B}} \rvert \boldsymbol{\mathcal{Y}} = \mathbf{y})1$. 

---

## Model Fitting
Finding the MLEs can be difficult, especially in the case of non-identity link functions, because it involves maximizing evaluating the possible intractable integral in Eq. \eqref{eq:marg-lik}. As a work-around, we can use an iterative process (one of which is covered in <a href="/stats-ml/glmm/">this post</a>). 

Suppose we are given some value of $\boldsymbol{\theta}$. We can use this to find the <strong>conditional mode</strong> of the random effects, denoted by $\tilde{\mathbf{u}}(\boldsymbol{\theta})$, and the <strong>conditional estimate</strong> of the fixed effects, denoted by $\tilde{\boldsymbol{\alpha}}$, which are defined as:

<aside><p>$\sigma$ is omitted from the notation because the conditional mode and estimate are independent of its value.</p></aside>

$$
\begin{equation}
\label{eq:obj-func}
\begin{aligned}
\begin{bmatrix}
\tilde{\mathbf{u}}(\boldsymbol{\theta}) \\
\tilde{\boldsymbol{\alpha}}(\boldsymbol{\theta})
\end{bmatrix}
&= \underset{\mathbf{u}, \boldsymbol{\alpha}}{\arg \max} \left\{ f_{\boldsymbol{\mathcal{U}} \rvert \boldsymbol{\mathcal{Y}}}(\mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{y}) \right\} \\
&= \underset{\mathbf{u}, \boldsymbol{\alpha}}{\arg \min} \left\{ -2 \log\left( f_{\boldsymbol{\mathcal{U}} \rvert \boldsymbol{\mathcal{Y}}}(\mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{y}) \right) \right\}
\end{aligned}
\end{equation}
$$

<aside><p>The objective in the second line of Eq. \eqref{eq:obj-func} is on the <strong>deviance</strong> scale.</p></aside>






### Derivations
Suppose we are dealing with a linear mixed model. All of the distributions are continuous, so we can deal with their densities and write these on the deviance scale:

$$
\begin{aligned}
-2 \log(f_{\boldsymbol{\mathcal{U}}})(\mathbf{u}) 
    &= -2 \left[-\frac{q}{2}\log(2 \pi) - \frac{1}{2} \log \left( \rvert \sigma^2 \mathbb{I}_{q \times q} \rvert \right) - \frac{1}{2}\mathbf{u}^\top \left[ \sigma^2 \mathbb{I}_{q \times q} \right]^{-1} \mathbf{u} \right]\\
    &= q \log(2 \pi \sigma^2) + \frac{\rvert \rvert \mathbf{u} \rvert \rvert_2^2}{\sigma^2} \\
-2 \log(f_{\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{U}}}(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{u}))
    &= -2 \left[ - \frac{n}{2}\log(2 \pi) - \frac{1}{2}\log\left(\rvert \sigma^2 \mathbb{I}_{n \times n} \rvert\right) - \frac{1}{2}(\mathbf{y} - \boldsymbol{\eta}(\mathbf{u}, \boldsymbol{\theta}, \boldsymbol{\alpha}))^\top \left[ \sigma^2 \mathbb{I}_{n \times n}\right]^{-1}(\mathbf{y} - \boldsymbol{\eta}(\mathbf{u}, \boldsymbol{\theta}, \boldsymbol{\alpha})) \right] \\
    &= n \log(2 \pi \sigma^2) + \frac{\rvert \rvert \mathbf{y} - \boldsymbol{\eta}(\mathbf{u}, \boldsymbol{\theta}, \boldsymbol{\alpha}) \rvert \rvert_2^2}{\sigma^2} \\
-2 \log(h(\mathbf{y}, \mathbf{u}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma)) 
    &= -2 \log(f_\boldsymbol{\mathcal{U}}(\mathbf{u})) -2 \log(f_{\boldsymbol{\mathcal{Y}} \rvert \boldsymbol{\mathcal{U}}}(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta}, \sigma \rvert \mathbf{u})) \\
    &= (n + q) \log(2 \pi \sigma^2) + \frac{\rvert \rvert \mathbf{y} - \boldsymbol{\eta}(\mathbf{u}, \boldsymbol{\alpha}, \boldsymbol{\theta}) \rvert \rvert_2^2 + \rvert \rvert \mathbf{u} \rvert \rvert_2^2}{\sigma^2}
\end{aligned}
$$

Here, we can define the <strong>discrepancy function</strong>:

$$
\begin{equation}
\label{eq:discrepancy}
d(\mathbf{u}; \mathbf{y} \rvert \boldsymbol{\alpha}, \boldsymbol{\theta}) = \underbrace{\rvert \rvert \mathbf{y} - \boldsymbol{\eta}(\mathbf{u}, \boldsymbol{\alpha}, \boldsymbol{\theta}) \rvert \rvert_2^2}_{(\star)} + \underbrace{\rvert \rvert \mathbf{u} \rvert \rvert_2^2}_{(\dagger)}
\end{equation}
$$

The first term, $(\star)$, is the sum of squared residuals based upon the linear predictor, and the second term, $(\dagger)$, is just the squared $\ell^2$-norm of $\mathbf{u}$. Thus, Eq. \eqref{eq:discrepancy} is essentially a penalized residual sum of squares where the penalty is on the size (in terms of $\ell^2$-norm) of $\mathbf{u}$. 

The discrepancy function can be rewritten in the following way:

$$
\begin{aligned}
d(\mathbf{u}; \mathbf{y} \rvert \boldsymbol{\alpha}, \boldsymbol{\theta}) 
&= 
\rvert \rvert  
\mathbf{y} - \mathbf{X} \boldsymbol{\alpha} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta})\mathbf{u} 
\rvert \rvert_2^2
+ \rvert \rvert \mathbf{u} - \mathbf{0}_q \rvert \rvert_2^2 \\
&= \left\rvert \left\rvert  
\mathbf{y} - 
\begin{bmatrix}
\mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) & \mathbf{X}
\end{bmatrix}
\begin{bmatrix} 
\mathbf{u} \\ \boldsymbol{\alpha}
\end{bmatrix}
\right\rvert \right\rvert_2^2
-
\left\rvert \left\rvert \mathbf{0}_{q} -
\begin{bmatrix}
\mathbb{I}_{q \times q} & \mathbf{0}_{q \times p}
\end{bmatrix}
\begin{bmatrix} 
\mathbf{u} \\
\boldsymbol{\alpha}
\end{bmatrix}
\right\rvert \right\rvert_2^2 \\
&= 
\left\rvert \left\rvert  
\begin{bmatrix}
\mathbf{y} \\ \mathbf{0}_{q}
\end{bmatrix} 
- 
\begin{bmatrix}
\mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) & \mathbf{X} \\
\mathbb{I}_{q \times q} & \mathbf{0}_{q \times p}
\end{bmatrix}
\begin{bmatrix} 
\mathbf{u} \\ \boldsymbol{\alpha}
\end{bmatrix}
\right\rvert \right\rvert_2^2
\end{aligned}
$$

This is simply the residual sum of squares for a standard linear model for an augmented dataset with $q$ additional <i>pseudo-observations</i> that are constructed such that the above holds. Notice that the discrepancy is quadratic in $\mathbf{u}$ and $\boldsymbol{\beta}$, so we minimize it by taking the derivative and equating that with zero. The minimizers of the discrepancy will then be given by $\begin{bmatrix} \tilde{\mathbf{u}}(\boldsymbol{\theta}) & \tilde{\boldsymbol{\alpha}}(\boldsymbol{\theta}) \end{bmatrix}^\top$ that satisfy:

<aside><p>Why?</p></aside>

$$
\begin{aligned}
\begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) + \mathbb{I}_{q \times q} &
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{X} \\
\mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) &
\mathbf{X}^\top \mathbf{X}
\end{bmatrix} 
\begin{bmatrix}
\tilde{\mathbf{u}}(\boldsymbol{\theta}) \\ 
\tilde{\boldsymbol{\alpha}}(\boldsymbol{\theta})
\end{bmatrix}
&= \begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{y} \\
\mathbf{X}^\top \mathbf{y}
\end{bmatrix}
\end{aligned}
$$

<details>
<summary>Proof.</summary>
The components of the vector of partial derivatives are:
$$
\begin{aligned}
\frac{\partial d(\mathbf{u}; \mathbf{y}; \boldsymbol{\alpha, \boldsymbol{\theta}})}{\partial \mathbf{u}}
&= \frac{\partial}{\partial \mathbf{u}} \left[ 
\left\rvert \left\rvert  
\begin{bmatrix}
\mathbf{y} \\ \mathbf{0}_{q}
\end{bmatrix} 
- 
\begin{bmatrix}
\mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) & \mathbf{X} \\
\mathbb{I}_{q \times q} & \mathbf{0}_{q \times p}
\end{bmatrix}
\begin{bmatrix} 
\mathbf{u} \\ \boldsymbol{\alpha}
\end{bmatrix}
\right\rvert \right\rvert_2^2 \right] \\
&= 2 \frac{\partial}{\partial \mathbf{u}^\top} \left[ \begin{bmatrix}
\mathbf{y} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} - \mathbf{X} \boldsymbol{\alpha} \\ \mathbf{0}_q - \mathbb{I}_{q \times q} \mathbf{u} - \mathbf{0}_{q \times p} \boldsymbol{\alpha}
\end{bmatrix} \right] \begin{bmatrix}
\mathbf{y} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} - \mathbf{X} \boldsymbol{\alpha} \\ \mathbf{0}_q - \mathbb{I}_{q \times q} \mathbf{u} - \mathbf{0}_{q \times p} \boldsymbol{\alpha}
\end{bmatrix}
\\
&= 2 \begin{bmatrix}
- \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \\
- \mathbb{I}_{q \times q}
\end{bmatrix}^\top  
\begin{bmatrix}
\mathbf{y} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} - \mathbf{X} \boldsymbol{\alpha} \\ \mathbf{0}_q - \mathbb{I}_{q \times q} \mathbf{u} - \mathbf{0}_{q \times p} \boldsymbol{\alpha}
\end{bmatrix}\\
&= 2 \left(   - \boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{y} + \boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} + \boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top\mathbf{X} \boldsymbol{\alpha} + \mathbf{u}\right) \\
\frac{\partial d(\mathbf{u}; \mathbf{y}; \boldsymbol{\alpha, \boldsymbol{\theta}})}{\partial \boldsymbol{\alpha}}
&= \frac{\partial}{\partial \boldsymbol{\alpha}} \left[ 
\left\rvert \left\rvert  
\begin{bmatrix}
\mathbf{y} \\ \mathbf{0}_{q}
\end{bmatrix} 
- 
\begin{bmatrix}
\mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) & \mathbf{X} \\
\mathbb{I}_{q \times q} & \mathbf{0}_{q \times p}
\end{bmatrix}
\begin{bmatrix} 
\mathbf{u} \\ \boldsymbol{\alpha}
\end{bmatrix}
\right\rvert \right\rvert_2^2 \right] \\
&= 2\frac{\partial}{\partial \boldsymbol{\alpha}^\top} \left[ \begin{bmatrix}
\mathbf{y} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} - \mathbf{X} \boldsymbol{\alpha} \\ \mathbf{0}_q - \mathbb{I}_{q \times q} \mathbf{u} - \mathbf{0}_{q \times p} \boldsymbol{\alpha}
\end{bmatrix} \right]  \begin{bmatrix}
\mathbf{y} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} - \mathbf{X} \boldsymbol{\alpha} \\ \mathbf{0}_q - \mathbb{I}_{q \times q} \mathbf{u} - \mathbf{0}_{q \times p} \boldsymbol{\alpha}
\end{bmatrix} \\
&= 2
\begin{bmatrix}
- \mathbf{X} \\
- \mathbf{0}_{q \times p}
\end{bmatrix}^\top \begin{bmatrix}
\mathbf{y} - \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} - \mathbf{X} \boldsymbol{\alpha} \\ \mathbf{0}_q - \mathbb{I}_{q \times q} \mathbf{u} - \mathbf{0}_{q \times p} \boldsymbol{\alpha}
\end{bmatrix} \\
&= 2 \left(  -\mathbf{X}^\top \mathbf{y} + \mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} + \mathbf{X}^\top \mathbf{X} \boldsymbol{\alpha} \right) 
\end{aligned}
$$
Setting this equal to zero and rearranging yields the system of equations:
$$
\begin{aligned}
\begin{bmatrix}
\frac{\partial d(\mathbf{u}; \mathbf{y}; \boldsymbol{\alpha, \boldsymbol{\theta}})}{\partial \mathbf{u}} \\
\frac{\partial d(\mathbf{u}; \mathbf{y}; \boldsymbol{\alpha, \boldsymbol{\theta}})}{\partial \boldsymbol{\alpha}}
\end{bmatrix}
&= \begin{bmatrix}
\mathbf{0}_q \\
\mathbf{0}_p
\end{bmatrix} \\
\implies
\begin{bmatrix}
2 \left( - \boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{y} + \boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} + \boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top\mathbf{X} \boldsymbol{\alpha} + \mathbf{u}\right)
\\
2 \left(  -\mathbf{X}^\top \mathbf{y} + \mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) \mathbf{u} + \mathbf{X}^\top \mathbf{X} \boldsymbol{\alpha} \right) 
\end{bmatrix}
&= \begin{bmatrix}
\mathbf{0}_q \\
\mathbf{0}_p
\end{bmatrix} \\
\implies
\begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) + \mathbb{I}_{q \times q} \\
\mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta})
\end{bmatrix} \mathbf{u} +
\begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{X} \\
\mathbf{X}^\top \mathbf{X} 
\end{bmatrix} \boldsymbol{\alpha}
&= \begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{y} \\
\mathbf{X}^\top \mathbf{y}
\end{bmatrix} \\
\implies 
\begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) + \mathbb{I}_{q \times q} &
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{X} \\
\mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) &
\mathbf{X}^\top \mathbf{X}
\end{bmatrix} 
\begin{bmatrix}
\mathbf{u} \\
\boldsymbol{\alpha}
\end{bmatrix}
&= \begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{y} \\
\mathbf{X}^\top \mathbf{y}
\end{bmatrix}
\end{aligned}
$$
</details>


Notice that the leading matrix is symmetric, and recall that we (usually) assume $\mathbf{X}$ to be full rank. This implies that $\mathbf{X}^\top \mathbf{X}$ should be invertible. We can then partition the matrix as:

$$
\begin{bmatrix}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) + \mathbb{I}_{q \times q} &
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{X} \\
\mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) &
\mathbf{X}^\top \mathbf{X} 
\end{bmatrix} 
= 
\begin{bmatrix}
\mathbf{A} & \mathbf{B} \\
\mathbf{B}^\top & \mathbf{C}
\end{bmatrix}
$$

with invertible sub-matrix $\mathbf{C}$. Using <a href="https://en.wikipedia.org/wiki/Schur_complement#Conditions_for_positive_definiteness_and_semi-definiteness">some matrix facts</a>, we have that the leading matrix will be positive definite if the Schur complement of $\mathbf{C}$ is also positive definite. That is, if:

$$
\mathbf{A} - \mathbf{B} \mathbf{C}^{-1} \mathbf{B}^\top \succ 0
$$

$$
\begin{aligned}
\mathbf{A} - \mathbf{B} \mathbf{C}^{-1} \mathbf{B}^{-1}  
&= \boldsymbol{\Lambda}^\top(\boldsymbol{\theta}) \mathbf{Z}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta}) + \mathbb{I}_{q \times q} - 
\mathbf{X}^\top \mathbf{Z} \boldsymbol{\Lambda}(\boldsymbol{\theta})
\left(\mathbf{X}^\top \mathbf{X} \right)^{-1}
\boldsymbol{\Lambda}^\top(\boldsymbol{\theta})\mathbf{Z}^\top \mathbf{X}
\end{aligned}
$$
