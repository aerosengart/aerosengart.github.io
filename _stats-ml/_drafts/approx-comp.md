---
layout: distill
title: Approximation Methods Comparison
description:
date: 2026-03-09
tabs: true
tags: theory research
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Background
      subsections:
        - name: Rewriting The Likelihood
    - name: Quadrature
      subsections:
        - name: Gaussian Case
        - name: Poisson Case
    - name: Marginal Quasi-Likelihood
      subsections:
        - name: Gaussian Case
        - name: Poisson Case
    - name: Laplace Approximation
      subsections:
        - name: Gaussian Case
        - name: Poisson Case
bibliography: stats-ml.bib
---

The question of how certain marginal likelihood approximations compare has come up in my research. While I await a response from the university library on acquiring <i>Comparison of Computationally Efficient Approximate Methods for Nonlinear and Generalized Linear Mixed Effects Models</i> by Sihaoyu Gao and Lang Wu, I've decided to try to do some preliminary investigation. 

---

## Background
Our setting will be a simple <a href="/stats-ml/glmm">(generalized) linear mixed model</a> with a fixed intercept and a random intercept. We'll assume our $N$ datapoints come from $K$ clusters, and the $i$-th response in the $k$-th cluster will be denoted by $y_i^k$. The clusters will be independent, and observations within the same cluster will be independent conditional on the random effects and independent of the random effects. For simplicity, we'll assume the clusters are all of size $n$. We'll further assume the random effects (Gaussian), so our set-up is:

$$
\begin{aligned}
\mathbb{E}\left[ y_i^k \rvert \beta_k \right] &= \mu_i^k \\
g(\mu_i^k) &= \alpha + \beta_k \\
\beta_k &\overset{iid}{\sim} \mathcal{N}(0, \tau^2)
\end{aligned}
$$

where $g(\cdot)$ is a differentiable and invertible link function. In this post, it will sometimes be $\log(\cdot)$ and sometimes be the identity, which correspond to the Poisson and Gaussian response distributions, respectively. These corresponds to the following conditional likelihoods:

$$
\begin{equation}
\label{eq:poi-gauss}
\begin{aligned}
\mathcal{L}_p(y_i^k; \theta \rvert \beta_k) &= \frac{(\mu_i^k)^{y_i^k} \exp\left(-\mu_i^k\right)}{y_i^k!} & \mu_i^k = \exp\left( \alpha + \beta_k \right) \\
\mathcal{L}_g(y_i^k; \theta \rvert \beta_k) &= \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(- \frac{(y_i^k - \mu_i^k)^2}{2 \sigma^2}\right) & \mu_i^k = \alpha + \beta_k
\end{aligned}
\end{equation}
$$

where $\sigma^2$ is the same across all observations. Let $\theta = (\alpha, \tau^2)^\top$ be the model parameter vector. The marginal likelihood can be obtained by integrating the random effects out of the joint likelihood:

$$
\begin{equation}
\label{eq:eq-setup}
\begin{aligned}
\mathcal{L}(\mathbf{y}; \theta)
    &= \int \mathcal{L}(\mathbf{y}; \theta \rvert \boldsymbol\beta) \mathcal{L}(\boldsymbol\beta) d \boldsymbol\beta
\end{aligned}
\end{equation}
$$

However, Eq. \eqref{eq:eq-setup} can be difficult to integrate (perhaps due to it being multidimensional or the conditional likelihood being cumbersome to work with), which leads us to the need for approximations. 

We will be concerned with three different methods: <a href="/stats-ml/quadrature">quadrature</a>, <a href="/stats-ml/glmm#marginal-quasi-likelihood">marginal quasi-likelihood</a>, and a <a href="https://en.wikipedia.org/wiki/Laplace%27s_approximation">Laplace approximation</a>.

In all that follows, we will assume that the fixed effects and any other nuisance parameters have been estimated under the null hypothesis that $\tau^2 = 0$. These estimates will be denoted with a hat.

### Rewriting The Likelihood
Let's first do some rewriting of the likelihood so that the later approximations are a bit easier. By the independence assumptions, we have:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}; \theta \rvert \boldsymbol \beta) 
    &= \prod_{k = 1}^K \prod_{i = 1}^n \mathcal{L}(y_i^k; \theta \rvert \beta_k) \\
\mathcal{L}(\mathbf{y}; \theta, \boldsymbol \beta)
    &= \prod_{k = 1}^K \prod_{i = 1}^n \mathcal{L}(y_i^k; \theta, \beta_k)
\end{aligned}
$$

We can then write:

$$
\begin{equation}
\label{eq:marg-lik}
\begin{aligned}
\mathcal{L}(\mathbf{y}; \theta) 
    &= \prod_{k = 1}^K \int \mathcal{L}(\mathbf{y}^k; \theta, \beta_k) d \beta_k \\
    &= \prod_{k = 1}^K \int \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \mathcal{L}(\beta_k) d \beta_k \\
    &= \prod_{k = 1}^K \int \frac{1}{\sqrt{2 \pi \tau^2}} \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \exp\left(- \frac{\beta_k^2}{2 \tau^2}\right) d \beta_k \\
\end{aligned}
\end{equation}
$$

Note that, in the Gaussian case, we can actually evaluate the integral.

{% tabs gauss-case %}
{% tab gauss-case statement %}
$$
\begin{aligned}
\mathcal{L}(\mathbf{y}; \theta)
    &=  \prod_{k = 1}^K \prod_{i = 1}^n \frac{1}{\sqrt{2\pi(\tau^2 + \sigma^2)}} \exp\left(- \frac{(y_i^k - \alpha)^2}{2(\tau^2 + \sigma^2)}\right)
\end{aligned}
$$
{% endtab %}
{% tab gauss-case proof %}
$$
\begin{aligned}
\mathcal{L}(y_i^k; \theta)
    &= \int \frac{1}{\sqrt{2 \pi \tau^2}} \mathcal{L}(y_i^k; \theta \rvert \beta_k) \exp\left( - \frac{\beta_k^2}{2 \tau^2} \right) d\beta_k \\
    &= \frac{1}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \int \exp\left(- \frac{(y_i^k - \mu_i^k)^2}{2 \sigma^2} \right) \exp\left(- \frac{\beta_k^2}{2 \tau^2} \right) d\beta_k \\
    &= \frac{1}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \int \exp\left(- \frac{y_i^k - 2(\alpha + \beta_k)y_i^k + (\alpha + \beta_k)^2}{2 \sigma^2} - \frac{\beta_k^2}{2 \tau^2} \right) d\beta_k \\    
    &= \frac{1}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \int \exp\left(- \frac{(y_i^k - \alpha)^2}{2 \sigma^2}\right) \exp\left( - \frac{(2\alpha - 2 y_i^k)\beta_k + \beta_k^2}{2 \sigma^2} - \frac{\beta_k^2}{2 \tau^2} \right) d\beta_k \\
    &= \frac{1}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \exp\left(- \frac{(y_i^k - \alpha)^2}{2 \sigma^2}\right) \int  \exp\left( - \frac{(\tau^2 + \sigma^2) \beta_k^2 - 2 \tau^2(y_i^k - \alpha) \beta_k}{2 \sigma^2 \tau^2}  \right) d\beta_k \\
    &= \frac{1}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \exp\left(- \frac{(y_i^k - \alpha)^2}{2 \sigma^2}\right) \int  \exp\left( - \frac{\beta_k^2 - 2 \frac{\tau^2(y_i^k - \alpha)}{\tau^2 + \sigma^2} \beta_k}{2\frac{\sigma^2 \tau^2}{\tau^2 + \sigma^2}}  \right) d\beta_k \\
    &= \frac{1}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \exp\left(- \frac{(y_i^k - \alpha)^2}{2 \sigma^2}\right) \int  \exp\left( - \frac{\left(\beta_k - \frac{\tau^2(y_i^k - \alpha)}{\tau^2 + \sigma^2}\right)^2 - \left(\frac{\tau^2(y_i^k - \alpha)}{\tau^2 + \sigma^2}\right)^2}{2 \frac{\sigma^2 \tau^2}{\tau^2 + \sigma^2}}\right) d\beta_k \\
    &= \frac{\sqrt{2 \pi \frac{\sigma^2 \tau^2}{\tau^2 + \sigma^2}}}{(\sqrt{2 \pi \tau^2})(\sqrt{2 \pi \sigma^2})} \exp\left(- \frac{(y_i^k - \alpha)^2}{2 \sigma^2} + \frac{\left(\frac{\tau^2(y_i^k - \alpha)}{\tau^2 + \sigma^2}\right)^2}{\frac{2 \sigma^2 \tau^2}{\tau^2 + \sigma^2}}\right) \underbrace{\int \frac{1}{\sqrt{2 \pi \frac{\sigma^2 \tau^2}{\tau^2 + \sigma^2}}}\exp\left( - \frac{\left(\beta_k - \frac{\tau^2(y_i^k - \alpha)}{\tau^2 + \sigma^2}\right)^2}{2\frac{ \sigma^2 \tau^2}{(\tau^2 + \sigma^2)}}\right) d\beta_k}_{=1} \\
    &= \frac{1}{\sqrt{2\pi(\tau^2 + \sigma^2)}} \exp\left(- \frac{(y_i^k - \alpha)^2}{2 \sigma^2} + \frac{\tau^2(y_i - \alpha)^2}{2 \sigma^2(\tau^2 + \sigma^2)}\right) \\
    &= \frac{1}{\sqrt{2\pi(\tau^2 + \sigma^2)}} \exp\left(- \frac{(\tau^2 + \sigma^2)(y_i^k - \alpha)^2 - \tau^2(y_i - \alpha)^2}{2 \sigma^2(\tau^2 + \sigma^2)}\right) \\
    &= \frac{1}{\sqrt{2\pi(\tau^2 + \sigma^2)}} \exp\left(- \frac{(y_i^k - \alpha)^2}{2(\tau^2 + \sigma^2)}\right) \\
\end{aligned}
$$
{% endtab %}
{% endtabs %}



---

## Quadrature
Quadrature allows us to approximate one-dimensional integrals of functions that follow particular forms. We will use <a href="/stats-ml/quadrature/#gauss-hermite-quadrature">Gauss-Hermite quadrature</a>, which approximates integrals of functions of the form $\exp\left(-x^2\right) f(x)$. We thus need to manipulate the joint likelihood into this form.

Let's define $x_k = \frac{\beta_k}{\sqrt{2 \tau^2}}$, which implies that $\beta_k = \sqrt{2 \tau^2 x_k^2}$. We then have:

$$
\begin{aligned}
x_k &= \frac{\beta_k}{\sqrt{2 \tau^2}} \\
\implies \frac{d x_k}{d \beta_k} &=  \frac{1}{\sqrt{2 \tau^2}} \\
\implies d \beta_k &= \sqrt{2 \tau^2} d x_k
\end{aligned}
$$

We substitute and write the marginal likelihood for cluster $k$ and $\beta_k$ as:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}; \theta)
    &= \prod_{k = 1}^K \int \frac{\sqrt{2 \tau^2}}{\sqrt{2 \pi \tau^2}} \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \exp\left(- x_k^2 \right) d x_k \\
    &=  \frac{1}{\sqrt{\pi}} \prod_{k = 1}^K \int \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \exp\left(- x_k^2 \right) d x_k
\end{aligned}
$$

Let's define:

$$
f(x_k) = \left. \left[ \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = \sqrt{2 \tau^2 x_k^2}}
$$

which is the conditional likelihood of cluster $k$ given $\beta_k = \sqrt{2 \tau^2 x_k^2}$. This gives us the form we need to use Gauss-Hermite quadrature!

The order $m$ approximation of the marginal likelihood is then given by:

$$
\begin{equation}
\label{eq:ghq-approx}
\begin{aligned}
\mathcal{L}_{ghq}(\mathbf{y}; \theta)
&= \pi^{-\frac{K}{2}} \prod_{k = 1}^K  \left( \sum_{j = 1}^m \frac{2^{m-1} m! \sqrt{\pi}}{m^2 [H_{m-1}(z_j)]^2} f(z_j) \right)
\end{aligned}
\end{equation}
$$

where $z_j$ are the roots of the $m$-th (physicists') Hermite polynomial. 

### Gaussian Case
In the Gaussian case, we have:

$$
\begin{aligned}
f(z_j) 
&= \left. \left[ \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = \sqrt{2 \tau^2 z_j^2}} \\
&= \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(- \frac{\left(y_i^k - \left(\alpha + \sqrt{2 \tau^2 z_j^2}\right)\right)^2}{2 \sigma^2}\right)
\end{aligned} 
$$

So then:

$$
\mathcal{L}_{ghq}(\mathbf{y}; \theta)
= \pi^{-\frac{K}{2}} \prod_{k = 1}^K  (2 \pi \sigma^2)^{-\frac{n}{2}} \sum_{j = 1}^m  \left( \frac{2^{m-1} m! \sqrt{\pi}}{m^2 [H_{m-1}(z_j)]^2} \exp\left(- \sum_{i = 1}^n \frac{\left(y_i^k - \left(\alpha + \sqrt{2 \tau^2 z_j^2}\right)\right)^2}{2 \sigma^2}\right) \right)
$$

### Poisson Case
In the Poisson case, we have:

$$
\begin{aligned}
f(z_j) 
&= \left. \left[ \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = \sqrt{2 \tau^2 z_j^2}} \\
&= \prod_{i = 1}^n \frac{(\exp(\alpha + \sqrt{2 \tau^2 z_j^2}))^{y_i^k} \exp\left(-\left(\exp\left(\alpha + \sqrt{2 \tau^2 z_j^2}\right)\right)\right)}{y_i^k!}
\end{aligned} 
$$

So then:

$$
\mathcal{L}_{ghq}(\mathbf{y}; \theta)
= (\pi)^{-\frac{K}{2}} \prod_{k = 1}^K  \sum_{j = 1}^m \left(  \frac{2^{m-1} m! \sqrt{\pi}}{m^2 [H_{m-1}(z_j)]^2} \prod_{i = 1}^n \left(\frac{(\exp(\alpha + \sqrt{2 \tau^2 z_j^2}))^{y_i^k} \exp\left(-\left(\exp\left(\alpha + \sqrt{2 \tau^2 z_j^2}\right)\right)\right)}{y_i^k!}\right) \right)
$$


---

## Marginal Quasi-Likelihood
The method of marginal quasi-likelihood begins with approximating the <i>conditional</i> log quasi-likelihood with a Taylor expansion about the random effects mean, $\boldsymbol \beta = \mathbf{0}_k$. For cluster $k$, this is given by:

<aside><p>We'll use $\ell$ to denote the log quasi-likelihood.</p></aside>

$$
\ell(\mathbf{y}^k; \theta \rvert \beta_k)
    \approx \left. \left[ \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} 
    + \beta_k \left. \left[ \frac{\partial \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right\rvert_{\beta_k = 0} 
    + \frac{1}{2} \beta_k^2 \left. \left[ \frac{\partial^2 \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = 0} 
$$

Exponentiating this and plugging it into Eq. \eqref{eq:marg-lik}:

$$
\begin{equation}
\label{eq:mql-1}
\begin{aligned}
\mathcal{L}_{mql}(\mathbf{y}^k; \theta)
    &= \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \tau^2}} \int \exp\left( \ell_{mql}(y_i^k; \theta \rvert \beta_k) \right) \exp\left(- \frac{\beta_k^2}{2 \tau^2} \right) d \beta_k \\
    &= \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \tau^2}} \int \exp\left( \left. \left[\ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} + \beta_k \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} + \frac{1}{2} \beta_k^2\left.\left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right \rvert_{\beta_k = 0} - \frac{\beta_k^2}{2 \tau^2} \right) d\beta_k \\
    &= \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \tau^2}} \exp\left(\left. \left[\ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} \right) \int \exp\left( \underbrace{\beta_k \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} - \frac{1}{2} \beta_k^2 \left( \frac{1}{\tau^2} - \left.\left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right \rvert_{\beta_k = 0} \right)}_{(\star)} \right) d\beta_k \\
\end{aligned}
\end{equation}
$$

Define:

$$
\sigma^2_{i,k} = \left( \frac{1}{\tau^2} - \left.\left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right \rvert_{\beta_k = 0} \right) ^{-1} = \frac{\tau^2}{1 - \tau^2 \left.\left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right \rvert_{\beta_k = 0}}
$$

We complete the square in $(\star)$:

$$
\begin{aligned}
(\star)
    &=  - \frac{1}{2 \sigma^2_{i,k}}\beta_k^2 + \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \beta_k  \\
    &= -\frac{1}{2\sigma^2_{i,k}} \left[ \beta_k^2 - 2 \sigma^2_{i,k} \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \beta_k  \right] \\
    &= \underbrace{- \frac{1}{2\sigma^2_{i,k}} \left(\beta_k - \sigma^2_{i,k} \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2}_{(*)} + \underbrace{\frac{1}{2\sigma^2_{i,k}} \left( \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2}_{(**)} \\
\end{aligned}
$$

Notice that $(*)$ is the unnormalized probability density function for a random variable with distribution:

$$
\mathcal{N}\left(\sigma^2_{i,k} \left. \left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k}\right] \right\rvert_{\beta_k = 0}, \sigma^2_{i,k} \right)
$$

and $(**)$ does not involve $\beta_k$. We can therefore simplify Eq. \eqref{eq:mql-1} to:

$$
\begin{aligned}
\mathcal{L}_{mql}(\mathbf{y}^k; \theta)
    &= \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \tau^2}} \exp\left(\left. \left[\ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} \right) \int \exp\left(\beta_k \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} - \frac{1}{2} \beta_k^2 \left( \frac{1}{\tau^2} - \left.\left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right \rvert_{\beta_k = 0} \right)\right) d\beta_k \\
    &= \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \tau^2}} \exp\left(\left. \left[\ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} \right) \frac{\sqrt{2 \pi \sigma^2_{i,k}}}{\sqrt{2 \pi \sigma^2_{i,k}}} \int \exp\left( - \frac{1}{2\sigma^2_{i,k}} \left(\beta_k - \sigma^2_{i,k} \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2 + \frac{1}{2\sigma^2_{i,k}} \left( \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2\right) d\beta_k \\
    &= \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left(\left. \left[\ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} + \frac{1}{2\sigma^2_{i,k}} \left( \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2\right) \underbrace{\int \frac{1}{\sqrt{2 \pi \sigma^2_{i,k}}}\exp\left( - \frac{1}{2\sigma^2_{i,k}} \left(\beta_k - \sigma^2_{i,k} \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2\right) d\beta_k}_{=1} \\
    &=  \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left(\left. \left[\ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} + \frac{1}{2\sigma^2_{i,k}} \left( \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2\right)
\end{aligned}
$$

Thus, the marginal quasi-likelihood approximation of the marginal likelihood is given by:

$$
\begin{equation}
\label{eq:mql-approx}
\mathcal{L}_{mql}(\mathbf{y}; \theta)
    =\prod_{k = 1}^K \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left(\left. \left[ \ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = 0} + \frac{1}{2\sigma^2_{i,k}} \left( \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = 0} \right)^2 \right)
\end{equation}
$$

### Gaussian Case
In the Gaussian case, we have:

$$
\begin{aligned}
\mathcal{L} (y_i^k; \theta \rvert \beta_k) 
&= \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(- \frac{(y_i^k - (\alpha + \beta_k))^2}{2 \sigma^2}\right) \\
\ell(y_i^k; \theta \rvert \beta_k)
&= -\frac{1}{2} \log(2 \pi \sigma^2) - \frac{(y_i^k - (\alpha + \beta_k))^2}{2 \sigma^2}
\end{aligned}
$$

We then have:

$$
\begin{equation}
\label{eq:beta-derivs-gauss}
\begin{aligned}
\frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k}
&= - \frac{1}{2\sigma^2}\left(-2(y_i^k - (\alpha + \beta_k))\right) \\
&= \frac{y_i^k - (\alpha + \beta_k)}{\sigma^2} \\
\frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2}
&= - \frac{1}{\sigma^2}
\end{aligned}
\end{equation}
$$

Plugging in $\boldsymbol \beta = \mathbf{0}_k$ yields:

$$
\begin{aligned}
\sigma^2_{i,k}
&= \left( \frac{1}{\tau^2} + \frac{1}{\sigma^2} \right)^{-1} 
= \frac{\tau^2 \sigma^2}{\tau^2 + \sigma^2} \\
\mathcal{L}_{mql}(\mathbf{y}; \theta)
&= \prod_{k = 1}^K \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left( -\frac{1}{2} \log(2 \pi \sigma^2) - \frac{(y_i^k - \alpha)^2}{2 \sigma^2} + \frac{1}{2 \sigma^2_{i,k}} \left( \frac{y_i^k - \alpha}{\sigma^2} \right)^2  \right)
\end{aligned}
$$

### Poisson Case
In the Poisson case, we have:

$$
\begin{aligned}
\mathcal{L}(y_i^k; \theta \rvert \beta_k) 
&= \frac{(\exp\left(\alpha + \beta_k \right))^{y_i^k} \exp\left(-\exp\left(\alpha + \beta_k \right) \right)}{y_i^k!} \\
\ell(y_i^k; \theta \rvert \beta_k) 
&= y_i^k \left(\alpha + \beta_k \right) - \exp\left(\alpha + \beta_k \right) - \log(y_i^k!) \\
\end{aligned}
$$

And thus:

$$
\begin{equation}
\label{eq:beta-derivs-poi}
\begin{aligned}
\frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k}
&= y_i^k - \exp\left(\alpha + \beta_k \right)\\
\frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2}
&= - \exp\left(\alpha + \beta_k \right)
\end{aligned}
\end{equation}
$$

Plugging in $\boldsymbol \beta = \mathbf{0}_k$ yields:

$$
\begin{aligned}
\sigma^2_{i,k}
&= \left( \frac{1}{\tau^2} + \exp(\alpha) \right)^{-1} 
= \frac{\tau^2}{1 + \tau^2 \exp(\alpha)} \\
\mathcal{L}_{mql}(\mathbf{y}; \theta)
&= \prod_{k = 1}^K \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left( y_i^k \alpha - \exp(\alpha) - \log(y_i^k!) + \frac{1}{2 \sigma^2_{i,k}} \left( y_i^k - \exp(\alpha) \right)^2  \right)
\end{aligned}
$$

---

## Laplace Approximation
In a <a href="https://en.wikipedia.org/wiki/Laplace%27s_approximation">Laplace approximation</a>, we do a Taylor approximation of the log-likelihood about the maximum a posterior estimate of $\boldsymbol \beta$, which we denote with $\hat{\boldsymbol \beta}$. All of the derivations are the same as in <a href="#marginal-quasi-likelihood">the previous section</a> except we use $\boldsymbol \beta = \hat{\boldsymbol \beta}$ instead of $\boldsymbol \beta = \mathbf{0}_k$. 

<aside><p>If $\ell$ is the log quasi-likelihood, then this is equivalent to <a href="/stats-ml/glmm/#penalized-quasi-likelihood">penalized quasi-likelihood</a> where $\hat{\boldsymbol \beta}$ is a current estimate.</p></aside>

$$
\begin{equation}
\label{eq:lap-approx}
\mathcal{L}_{lap}(\mathbf{y}; \theta)
    = \prod_{k = 1}^K \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left(\left. \left[ \ell(y_i^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = \hat{\beta}_k} + \frac{1}{2\sigma^2_{i,k}} \left( \left.\left[ \frac{\partial \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k} \right] \right \rvert_{\beta_k = \hat{\beta}_k} \right)^2 \right)
\end{equation}
$$

where:

$$
\sigma^2_{i,k} = \left(\frac{1}{\tau^2} - \left. \left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2}\right] \right\rvert_{\beta_k = \hat{\beta}_k} \right)^{-1} = \frac{\tau^2}{1 - \tau^2 \left. \left[ \frac{\partial^2 \ell(y_i^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \hat{\beta}_k}}
$$

### Gaussian Case
We use the derivatives we found in Eq. \eqref{eq:beta-derivs-gauss} and plug in $\boldsymbol \beta = \hat{\boldsymbol \beta}$ to get:

$$
\begin{aligned}
\sigma^2_{i,k}
&= \left( \frac{1}{\tau^2} + \frac{1}{\sigma^2} \right)^{-1} 
= \frac{\tau^2 \sigma^2}{\tau^2 + \sigma^2} \\
\mathcal{L}_{mql}(\mathbf{y}; \theta)
&= \prod_{k = 1}^K \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left( -\frac{1}{2} \log(2 \pi \sigma^2) - \frac{(y_i^k - (\alpha + \hat{\beta}_k))^2}{2 \sigma^2} + \frac{1}{2 \sigma^2_{i,k}} \left( \frac{y_i^k - (\alpha + \hat{\beta}_k)}{\sigma^2} \right)^2  \right)
\end{aligned}
$$

### Poisson Case
In the Poisson case, we have:

$$
\begin{aligned}
\sigma^2_{i,k}
&= \left( \frac{1}{\tau^2} + \exp(\alpha + \hat{\beta}_k) \right)^{-1} 
= \frac{\tau^2}{1 + \tau^2 \exp(\alpha + \hat{\beta}_k)} \\
\mathcal{L}_{mql}(\mathbf{y}; \theta)
&= \prod_{k = 1}^K \prod_{i = 1}^n \sqrt{\frac{\sigma^2_{i,k}}{\tau^2}} \exp\left( y_i^k (\alpha + \hat{\beta}_k) - \exp(\alpha + \hat{\beta}_k) - \log(y_i^k!) + \frac{1}{2 \sigma^2_{i,k}} \left( y_i^k - \exp(\alpha + \hat{\beta}_k) \right)^2  \right)
\end{aligned}
$$

---

## Likelihood Expansion
One last (bonus) method is to perform a Taylor approximation of the conditional <i>likelihood</i> directly rather than the exponentiating an approximation of the log-likelihood. We can then substitute in some value $\dot{\boldsymbol \beta}$ for $\boldsymbol \beta$ (likely either $\hat{\boldsymbol \beta}$ or $\mathbf{0}_k$). For cluster $k$, we have:

$$
\mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k)
\approx 
\left. \left[ \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \right] \right\rvert_{\beta_k = \dot{\beta}_k} + \beta_k \left. \left[ \frac{\partial \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k}\right] \right\rvert_{\beta_k = \dot{\beta}_k} + \frac{1}{2} \beta_k^2 \left. \left[ \frac{\partial^2 \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k} 
$$

Since $\mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) = \exp(\ell(\mathbf{y}^k; \theta \rvert \beta_k))$, we can write:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k)
&\approx 
\left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)  \right] \right\rvert_{\beta_k = \dot{\beta}_k} + \beta_k \left. \left[ \frac{\partial \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right) }{\partial \beta_k}\right] \right\rvert_{\beta_k = \dot{\beta}_k} + \frac{1}{2} \beta_k^2 \left. \left[ \frac{\partial^2 \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k}  \\
&= \left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)  \right] \right\rvert_{\beta_k = \dot{\beta}_k} + \beta_k \left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)  \right] \right\rvert_{\beta_k = \dot{\beta}_k}\left. \left[ \frac{\partial \ell(\mathbf{y}^k; \theta \rvert \beta_k) }{\partial \beta_k}\right] \right\rvert_{\beta_k = \dot{\beta}_k} + \frac{1}{2} \beta_k^2 \left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right) \right] \right\rvert_{\beta_k = \dot{\beta}_k} \left. \left[ \frac{\partial^2 \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k}  \\
&= \left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)  \right] \right\rvert_{\beta_k = \dot{\beta}_k} \left(1 + \beta_k \left. \left[ \frac{\partial \ell(\mathbf{y}^k; \theta \rvert \beta_k) }{\partial \beta_k}\right] \right\rvert_{\beta_k = \dot{\beta}_k} + \frac{1}{2} \beta_k^2 \left. \left[ \frac{\partial^2 \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k} \right) 
\end{aligned}
$$

The marginal likelihood can then be computed as the expectation of the above expression with respect to the distribution of $\boldsymbol \beta$:

$$
\begin{aligned}
\mathbb{E}\left[ \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k) \right]
&\approx  \left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)  \right] \right\rvert_{\beta_k = \dot{\beta}_k} \mathbb{E}\left[ 1 + \beta_k \left. \left[ \frac{\partial \ell(\mathbf{y}^k; \theta \rvert \beta_k) }{\partial \beta_k}\right] \right\rvert_{\beta_k = \dot{\beta}_k} + \frac{1}{2} \beta_k^2 \left. \left[ \frac{\partial^2 \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k} \right] \\
&= \left. \left[ \exp\left( \ell(\mathbf{y}^k; \theta \rvert \beta_k) \right)  \right] \right\rvert_{\beta_k = \dot{\beta}_k} \left(1 + \frac{1}{2}\text{Var}(\beta_k) \left. \left[ \frac{\partial^2 \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k} \right) & \left(\mathbb{E}[\beta_k] = 0 \right) \\
&= \left. \left[ \mathcal{L}(\mathbf{y}^k; \theta \rvert \beta_k)  \right] \right\rvert_{\beta_k = \dot{\beta}_k} \left(1 + \frac{\tau^2}{2} \left. \left[ \frac{\partial^2 \ell(\mathbf{y}^k; \theta \rvert \beta_k)}{\partial \beta_k^2} \right] \right\rvert_{\beta_k = \dot{\beta}_k} \right)
\end{aligned}
$$

### Gaussian Case
We use the derivatives we found in Eq. \eqref{eq:beta-derivs-gauss} to get:

$$
\begin{aligned}
\mathcal{L}_{le}(\mathbf{y}; \theta)
&= \prod_{k = 1}^K \left[ \left(\prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(- \frac{(y_i^k -\mu_i^k)^2}{2 \sigma^2} \right)\right) \left(1 - \frac{\tau^2}{2 \sigma^2} \right) \right]
\end{aligned}
$$

### Poisson Case
In the Poisson case, we have:

$$
\begin{aligned}
\mathcal{L}_{le}(\mathbf{y}; \theta)
&= \prod_{k = 1}^K \left[ \left( \prod_{i = 1}^n \frac{(\exp(\alpha + \dot{\beta}_k))^{y_i^k} \exp(-\exp(\alpha + \dot{\beta}_k))}{y_i^k!}\right) \left(1 - \frac{\tau^2}{2} \exp(\alpha + \dot{\beta}_k) \right)\right]
\end{aligned}
$$

