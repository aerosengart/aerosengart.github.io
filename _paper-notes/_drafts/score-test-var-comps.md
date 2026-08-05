---
layout: distill
title: A Score Test for Variance Components in Generalized Linear Models
description:
date: 2026-01-05
tabs: true
tags: glmm likelihood theory paper-review
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Set-Up
      subsections:
          - name: Some Non-Standard Notation
          - name: The Problem
    - name: Quasi-Likelihood
      subsections:
          - name: Intermediate Quantities
          - name: Rewriting The Likelihood
    - name: Score
    - name: Information
      subsections:
          - name: Intermediate Quantities
          - name: Information Components
    - name: Individual Variance Component Tests
      subsections:
          - name: Intermediate Quantities
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: 2026-01-05-score-test-lin.bib
---

In this page, I will be trying to derive the results in Lin (1997)<d-cite key=lin1997></d-cite> myself as a way to make sure I understand what's going on.


## Set-Up
We will assume to have $N$ total observations. For each, we have some response, $y_i$, a $p \times 1$ fixed effects covariate vector, $\mathbf{x}_i$, and a $q \times 1$ random effects covariate vector, $\mathbf{z}_i$. We assume that, conditional on the random effects, $\beta$, the responses are independent and have means and variances:

$$
\begin{aligned}
\mathbb{E}[y_j \rvert \beta] &= \mu_j \\
\text{Var}(y_j \rvert \beta) &= V(\mu_j) = \phi v(\mu_j)
\end{aligned}
$$

for some scale parameter, $\phi$. 

We will assume to have a clustered design, so each observation will belong to some cluster $k \in [K]$. We further assume a generalized linear model with monotonic and differentiable $g(\cdot)$, and $p \times 1$ vector of fixed effects is appropriate. We assume that $\mathbf{D}(\tau^2)$ is a $q \times q$ diagonal matrix parametrized by the $m$-dimensional vector-valued variance component, $\tau^2 = (\tau^2_1, \dots, \tau^2_m)^\top$. Our model is:

<aside><p>The observations in cluster $k$ will be denoted with a superscript $k$ (i.e. $y_j^k$), and we'll let $n_k$ denote the number of observations in cluster $k$ (such that $\sum_{k = 1}^K n_k = N$).</p></aside>

$$
\begin{aligned}
g(\mu^k_i) &= \eta^k_i = \mathbf{x}_i^\top \alpha + \mathbf{z}_i^\top \beta_k \\
\beta &\sim \mathcal{N}(\mathbf{0}_q, \mathbf{D}(\tau^2))
\end{aligned}
$$

where $\beta_k$ is the $h$-dimensional sub-vector of $\beta$ corresponding to cluster $k$. We've kept things general so that we match the notation in the original article (for the most part), but it's important to note that when we are dealing with clustered data with independent clusters, we have that:

$$
\begin{aligned}
\beta_i &\in \mathbb{R}^h \\
\beta &= (\beta_1, \dots, \beta_K) \\
\mathbf{D}(\tau^2) &= \mathbf{M}(\tau^2) \otimes \mathbb{I}_{k \times k}
\end{aligned}
$$

where $\mathbf{M}(\tau^2)$ is the $h \times h$ covariance matrix for a single cluster. In this case, $q = h \times k$, and the random effects covariates vector for a single observation, $\mathbf{z}_i$, is a vector of zeros with the appropriate coordinates replaced with the $h$ covariate values. 

We'll construct vectors of our data for each cluster:

$$
\begin{aligned}
\mathbf{y}^k &= (y^k_1, \dots, y^k_{n_k})^\top \\
\mu^k &= (\mu^k_1, \dots, \mu^k_{n_k})^\top \\
\eta^k &= (\eta^k_1, \dots, \eta^k_{n_k})^\top
\end{aligned}
$$

And we can concatenate these vectors to construct large vectors for the complete data:

$$
\begin{aligned}
\mathbf{y} &= (\mathbf{y}^1, \dots, \mathbf{y}^K)^\top = (y_1, \dots, y_N)^\top \\
\mu &= (\mu^1, \dots \mu^K)^\top = (\mu_1, \dots, \mu_N)^\top \\
\eta &= (\eta^1, \dots, \eta^K)^\top = (\eta_1, \dots, \eta_N)^\top
\end{aligned}
$$

Furthermore, we will let $\mathbf{V}(\mu)$ denote the diagonal matrix with the $V(\mu_i) = \phi v(\mu_i)$ terms along its main diagonal, and $\mathbf{X}$ will be the $N \times p$ matrix with rows equal to the $\mathbf{x}_i^\top$ vectors.

<div class="example">
  <strong>Example.</strong>
  <br>
    In this example, we'll assume we have negative binomial data with a scalar variance component (i.e. $m = 1$), a fixed (cluster-specific) intercept, a fixed slope, and a corresponding random slope (so $z_i^k = x_{i,0}^k$): 
    $$
    g(\mu_i^k) = \alpha^0_k + \alpha_1 x_{i,0}^k+ \beta_k z_i^k; 
    \hspace{5mm}
    \beta_k \overset{iid}{\sim} \mathcal{N}(0, \tau^2)
    $$
    The link and variance functions are:
    $$
    g(\mu^k_i) = \log(\mu^k_i);
    \hspace{5mm}
    v(\mu_i^k) = \mu_i^k + \frac{1}{\gamma} (\mu_i^k)^2
    $$
    where we assume $\gamma$ is a known dispersion parameter. We are interested in testing:
    $$
    H_0: \alpha_1 = 0, \beta_k = 0 \hspace{1mm} \text{for } k \in [K]
    \hspace{5mm} \text{vs.} \hspace{5mm} 
    H_1: \alpha_1 \neq 0, \beta_k \neq 0 \hspace{1mm} \text{for } k \in [K]
    $$
</div>

### Some Non-Standard Notation
We'll use a bar (e.g. $$\bar{\mu}_i$$) to denote a term evaluated under the null hypothesis that $$\tau^2 = \mathbf{0}_{m}$$, and a hat (e.g. $$\hat{\mu}_i$$) to denote a term evaluated at the parameter maximum likelihood estimates under the null hypothesis as well.

### The Problem
We are interesting in assessing the need for heterogeneity in the effects of our model. We test this by considering the null hypothesis that $$\tau^2 = \mathbf{0}_m$$. Notice that, under $H_0$, the random effects have zero variance. Since they are assumed to have mean zero, this implies that $$\beta = \mathbf{0}_{q}$$, and the null hypothesis can be changed to this. 

---

## Quasi-Likelihood
The model has parameters $\alpha$, $\phi$, $\tau^2$, and random effects, $\beta$. Let $\theta = (\alpha, \phi, \tau^2)$. The conditional log quasi-likelihood is:

$$

\ell_q(y_i; \theta \rvert \beta) = \int_{y_i}^{\mu_i} \frac{y_i - u}{\phi v(u)} du
$$ 

Let $f(\beta; \tau^2)$ and $\mathcal{F}(\beta; \tau^2)$ denote the log density and distribution functions of $\beta$. Under the assumption that $\mathbf{D}(\tau^2)$ is diagonal (i.e. the random effects are independent), the integrated (marginal) quasi-likelihood is:

$$
\begin{eqn}
\mathcal{L}_q(\mathbf{y}; \theta) 
&= \int \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) d \mathcal{F}(\beta; \tau^2) 
\end{eqn}
$$

And the marginal log quasi-likelihood is the logarithm of this function:

$$
\begin{eqn}
\mathcal{L}_q(\mathbf{y}; \theta) 
= \log \left[ \int \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta)  \right) d \mathcal{F}(\beta; \tau^2) \right]
\end{eqn}
$$

This integral is unwieldy, so we use a second-order Taylor approximation of the conditional quasi-likelihood about the random effects mean to simplify our calculations. We therefore construct an approximation via Taylor expansion. 

<!-- #region log-marg-quasi-lik -->
{% tabs log-marg-quasi-lik %}
{% tab log-marg-quasi-lik statement %}
$$
\begin{aligned}
\ell_q(\mathbf{y}; \theta) 
&\approx \sum_{i = 1}^N \left. \ell_q(y_i; \theta \rvert \beta)\right\rvert_{\tau^2 = \mathbf{0}_{m}}+ \frac{1}{2} \text{tr} \left[ \mathbf{Z}^\top \left( \left. \left[ \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta} \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta^\top} + \frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta )}{\partial \eta \partial \eta^\top} \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} \right) \mathbf{Z} \mathbf{D}(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert)
\end{aligned}
$$

where $$\frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta}$$ is the $N$-dimensional vector whose $i$-th element is equal to $$\frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i}$$, and $$\frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta )}{\partial \eta \partial \eta^\top}$$ is the $N \times N$ dimensional diagonal matrix whose $i$-th diagonal element is equal to $$\frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2}$$. 
{% endtab %}
{% tab log-marg-quasi-lik proof %}
The integrated (conditional) quasi-likelihood can be written via a Taylor expansion:

$$
\begin{aligned}
\mathcal{L}_q\left( \mathbf{y}; \theta \rvert \beta \right)
&= \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \\
&= \left. \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \right\rvert_{\tau^2 = \mathbf{0}_{m}} + \left. \frac{\partial}{\partial \beta} \left[ \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \right]\right\rvert_{\tau^2 = \mathbf{0}_{m}} \beta + \frac{1}{2}\beta^\top  \left[ \left. \frac{\partial}{\partial \beta \partial \beta^\top} \left[ \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right)  \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} \right]  \beta + \dots \\
&= \left. \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \right\rvert_{\tau^2 = \mathbf{0}_{m}} \left( 1 + \left. \frac{\partial}{\partial \beta} \left[ \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right]\right\rvert_{\tau^2 = \mathbf{0}_{m}} \beta + \frac{1}{2}\beta^\top \left[  \left.\frac{\partial}{\partial \beta \partial \beta^\top} \left[ \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} \right]  \beta + e \right)
\end{aligned}
$$

where $e$ is the residual that holds the higher order terms. Taking the expectation and applying the chain rule yields an approximate marginal quasi log-likelihood:

$$
\mathcal{L}_q(\mathbf{y}; \theta) 
\approx \left. \exp\left(\sum_{i = 1}^N \ell_q(y_i; \theta \rvert \beta) \right)\right\rvert_{\tau^2 = \mathbf{0}_{m}} \left(  1 + \frac{1}{2} \text{tr} \left[ \left( \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \mathbf{z}_i \right)  \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta )}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \mathbf{z}_i^\top \right) + \sum_{i = 1}^N \left. \frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2}\right\rvert_{\tau^2 = \mathbf{0}_{m}} \mathbf{z}_i \mathbf{z}_i^\top \right) \mathbf{D}(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert) \right)
$$

Taking the logarithm and using the fact that $\log(x + 1) \approx x$ for small $x$, we can write the marginal log quasi-likelihood approximation as:

$$
\begin{aligned}
\ell_q(\mathbf{y}; \theta) 
&\approx 
\sum_{i = 1}^N \left. \ell_q(y_i; \theta \rvert \beta)\right\rvert_{\tau^2 = \mathbf{0}_{m}} + \frac{1}{2} \text{tr} \left[ \left( \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \mathbf{z}_i \right)  \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta )}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \mathbf{z}_i^\top \right) + \sum_{i = 1}^N \left. \frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \mathbf{z}_i \mathbf{z}_i^\top \right) \mathbf{D}(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert) \\
&= \sum_{i = 1}^N \left. \ell_q(y_i; \theta \rvert \beta)\right\rvert_{\tau^2 = \mathbf{0}_{m}}+ \frac{1}{2} \text{tr} \left[ \mathbf{Z}^\top \left( \left. \left[ \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta} \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta^\top} + \frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta )}{\partial \eta \partial \eta^\top} \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} \right) \mathbf{Z} \mathbf{D}(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert)
\end{aligned}
$$

where $$\frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta}$$ is the $N$-dimensional vector whose $i$-th element is equal to $$\frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i}$$, and $$\frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta )}{\partial \eta \partial \eta^\top}$$ is the $N \times N$ dimensional diagonal matrix whose $i$-th diagonal element is equal to $$\frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2}$$. 
{% endtab %}
{% endtabs %}
<!-- #endregion -->

We will drop the $o(\rvert \rvert \tau^2 \rvert \rvert)$ in what follows as we have assumed moment conditions hold that make this term negligible.

### Intermediate Quantities
We will define some helpful quantities that will make the notation easier as we continue. Recall that a single subscript (with no superscript) directly indexes the <i>entire</i> vector (over all clusters). We should also note that we further assume that $g(\cdot)$ is <i>strictly monotone</i> and <i>continuous</i> so that:

$$
\begin{equation}
\label{eq:deriv-assumption}
\frac{\partial g^{-1}(\eta_i)}{\partial \eta_i} = \frac{\partial \mu_i}{\partial \eta_i} = \left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{-1} = \left(\frac{\partial g(\mu_i)}{\partial \mu_i}\right)^{-1} 
\end{equation}
$$

This implies that $g(\cdot)$ is invertible, and its derivative is never $0$. This holds for many standard choices of link function (e.g. $\log$, $\text{logit}$, the identity). 

We now define:

$$
\begin{equation}
\label{eq:intermediate-terms}
\begin{aligned}
\delta_i &= \left[ \frac{\partial g(\mu_i)}{\partial\mu_i}\right]^{-1} \\
\omega_i &= \left[ V(\mu_i) \left(\frac{\partial g(\mu_i)}{\partial\mu_i}\right)^2 \right]^{-1} = \left[ \phi v(\mu_i) \left(\frac{\partial g(\mu_i)}{\partial\mu_i}\right)^2 \right]^{-1} = \frac{\delta_i^2}{\phi v(\mu_i)} \\
e_i &= \frac{\frac{\partial V(\mu_i)}{\partial \mu_i} \frac{\partial g(\mu_i)}{\partial \mu_i} + V(\mu_i) \frac{\partial^2 g(\mu_i)}{\partial \mu_i^2}}{V^2(\mu_i)\left(\frac{\partial g(\mu_i)}{\partial \mu_i}\right)^3} = \frac{\phi\left(\frac{\partial v(\mu_i)}{\partial\mu_i}\right)\left(\frac{\partial g(\mu_i)}{\partial\mu_i}\right) + \phi v(\mu_i)\left(\frac{\partial^2 g(\mu_i)}{\partial \mu_i^2}\right)}{\phi^2 v^2(\mu_i) \left(\frac{\partial g(\mu_i)}{\partial\mu_i}\right)^3} \\
\xi_i &= \omega_i + e_i(y_i - \mu_i)
\end{aligned}
\end{equation}
$$

<aside><p>In the original publication, $\omega_i$ is denoted by $w_i$, and $\xi_i$ is denoted by $w_{oi}$.</p></aside>

We also define the $N \times N$ diagonal matrices $\Delta$, $\Omega$, and $\Xi$, which have the $\delta_i$, $\omega_i$, and $\xi_i$ terms along their main diagonals (respectively).

<!-- 
<div class="example">
  <strong>Example.</strong>
  <br>
  As another reminder, the link function is $\log(\cdot)$, and the variance function is $v(\mu_i) = \mu_i + \frac{1}{\gamma} (\mu_i)^2$. Thus:
  $$
  \begin{aligned}
  \frac{\partial g(\mu_i)}{\partial \mu_i} &= \frac{\partial}{\partial \mu_i} \left[ \log(\mu_i) \right] = \frac{1}{\mu_i} \\
  \frac{\partial^2 g(\mu_i)}{\partial \mu_i^2} &= \frac{\partial}{\partial \mu_i} \left[ \frac{1}{\mu_i} \right] = -\frac{1}{\mu_i^2} \\
  \frac{\partial v(\mu_i)}{\partial \mu_i} &= \frac{\partial}{\partial \mu_i}\left[ \mu_i + \frac{1}{\gamma}\mu_i^2 \right] = 1 + \frac{2}{\gamma} \mu_i
  \end{aligned}
  $$
  The intermediate terms are then:
  $$
  \begin{aligned}
  \delta_i &= \left[ \frac{1}{\mu_i} \right]^{-1} = \mu_i \\
  \omega_i &= \frac{\mu_i^2}{\phi \left(\mu_i + \frac{1}{\gamma} (\mu_i)^2\right)} = \frac{\mu_i}{\phi \left(1 + \frac{1}{\gamma} \mu_i\right)} \\
  e_i &= \frac{\phi \left(1 + \frac{2}{\gamma} \mu_i\right) \left(\frac{1}{\mu_i}\right) + \phi\left(\mu_i + \frac{1}{\gamma} \mu_i^2 \right)\left(- \frac{1}{\mu_i^2}\right)}{\phi^2 \left(\mu_i + \frac{1}{\gamma}\mu_i^2 \right)^2\left(\frac{1}{\mu_i^3}\right)} = \frac{\frac{1}{\mu_i} + \frac{2}{\gamma} - \frac{1}{\mu_i} - \frac{1}{\gamma}}{\phi \left(1+ \frac{1}{\gamma}\mu_i\right)^2 \frac{1}{\mu_i}} = \frac{\mu_i}{\phi \gamma\left(1+ \frac{1}{\gamma} \mu_i \right)^2}\\
  \xi_i &= \frac{\mu_i}{\phi\left(1 + \frac{1}{\gamma} \mu_i \right)} + \frac{\mu_i(y_i - \mu_i)}{\phi \gamma\left(1+ \frac{1}{\gamma} \mu_i \right)^2} = \frac{\gamma \mu_i \left(1 + \frac{1}{\gamma} \mu_i \right) + \mu_i}{\phi \gamma \left(1 + \frac{1}{\gamma} \mu_i \right)^2} = \frac{\mu \left( 1 + \gamma + \mu_i \right)}{\phi \gamma \left(1 + \frac{1}{\gamma} \mu_i \right)^2}
  \end{aligned}
  $$
</div> 
-->


### Rewriting The Likelihood
Notice that our expression for the approximate marginal log quasi-likelihood involves the first and second order partial derivatives of the conditional log quasi-likelihood. We will use the terms we introduced above to rewrite the equation. We have:

<!-- #region deriv-1-eta -->
{% tabs deriv-1-eta %}
{% tab deriv-1-eta equation %}
$$
\begin{equation}
\label{eq:deriv-1-eta}
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} 
&= \bar{\omega}_i \bar{\delta}_i^{-1} \left(y_i - \bar{\mu}_i\right)
&\implies
\left. \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta} \right\rvert_{\tau^2 = \mathbf{0}_{m}}
&= \bar{\Omega}\bar{\Delta}^{-1} (\mathbf{y} - \bar{\mu})
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-1-eta proof %}
We first take the partial derivative with respect to $\mu_i$:

$$
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \mu_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} 
&= \left. \frac{\partial}{\partial \mu_i} \left[ \int_{y_i}^{\mu_i} \frac{y_i - u}{\phi v(u)} du \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} \\
&= \left. \left[ \frac{y_i - \mu_i}{\phi v(\mu_i)} \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} & \left(\text{FTC}\right) \\
&= \frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)}
\end{aligned}
$$

We then have, by the chain rule:

$$
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} 
&= \left. \left[ \frac{\partial \mu_i}{\partial \eta_i} \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \mu_i} \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} & \left(\text{chain rule} \right)\\
&= \left. \left[ \frac{\partial \mu_i}{\partial \eta_i} \frac{\partial}{\partial \mu_i} \left[ \int_{y_i}^{\mu_i} \frac{y_i - u}{\phi v(u)} du \right] \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}} \\
&= \left. \left[ \left(\frac{\partial \mu_i}{\partial \eta_i}\right) \left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right) \right] \right\rvert_{\tau^2 = \mathbf{0}_{m}}  & \left( \text{FTC} \right) \\
&= \left. \frac{\partial \mu_i}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \left(\frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)}\right)
\end{aligned}
$$

Note that, by Eq. $\eqref{eq:deriv-assumption}$:

$$
\begin{aligned}
\delta_i &= \left[\frac{\partial g(\mu_i)}{\partial \mu_i}\right]^{-1} 
= \left[\frac{\partial \eta_i}{\partial \mu_i}\right]^{-1} 
= \frac{\partial \mu_i}{\partial \eta_i}
\end{aligned}
$$

So we have:

$$
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} 
&= \left. \frac{\partial \mu_i}{\partial \eta_i} \right\rvert_{\tau^2 = \mathbf{0}_{m}} \left(\frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)}\right) \\
&= \left. \frac{\left(\frac{\partial \mu_i}{\partial \eta_i}\right)^2}{\phi v(\mu_i) \frac{\partial \mu_i}{\partial \eta_i}}\right\rvert_{\tau^2 = \mathbf{0}_{m}} \left(y_i - \bar{\mu}_i\right) \\
&= \bar{\omega}_i \bar{\delta}_i^{-1} \left(y_i - \bar{\mu}_i\right)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

Similarly, we can find the second-order partial derivatives with respect to the linear predictor:

<!-- #region deriv-2-eta -->
{% tabs deriv-2-eta %}
{% tab deriv-2-eta equation %}
$$
\begin{equation}
\label{eq:deriv-2-eta}
\begin{aligned}
\left. \frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2} \right\rvert_{\tau^2 = \mathbf{0}_{m}}
&= -\bar{\xi}_i 
&\implies \left. \frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta \partial \eta^\top} \right\rvert_{\tau^2 = \mathbf{0}_{m}} 
&= -\bar{\Xi}
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-2-eta proof %}
First, we note that:

$$
\begin{equation}
\label{eq:temp}
\begin{aligned}
\frac{\partial^2 \mu_i}{\partial \eta_i^2} 
&= \left(\frac{\partial \mu_i}{\partial \eta_i} \right) \frac{\partial}{\partial \mu_i}\left[ \frac{\partial \mu_i}{\partial \eta_i}\right] &\left(\text{chain rule}\right) \\
&= \left(\frac{\partial \mu_i}{\partial \eta_i} \right) \frac{\partial}{\partial \mu_i}\left[ \left( \frac{\partial \eta_i}{\partial \mu_i} \right)^{-1} \right] &\left(\text{Eq. } \eqref{eq:deriv-assumption} \right) \\
&= - \left(\frac{\partial \mu_i}{\partial \eta_i} \right) \left(\frac{\partial \eta_i}{\partial \mu_i} \right)^{-2} \frac{\partial}{\partial \mu_i}\left[ \frac{\partial \eta_i}{\partial \mu_i} \right] \\
&= -  \left(\frac{\partial \eta_i}{\partial \mu_i} \right)^{-3} \left(\frac{\partial^2 \eta_i}{\partial \mu_i^2} \right) & \left(\text{Eq. } \eqref{eq:deriv-assumption}\right)\\
\implies \frac{\partial^2 \eta_i}{\partial \mu_i^2} &= -\left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{3} \left(\frac{\partial^2 \mu_i}{\partial \eta_i^2}\right)
\end{aligned}
\end{equation}
$$

It follows that:

$$
\begin{aligned}
\frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta_i^2}
&= \frac{\partial}{\partial \eta_i}\left[ \omega_i \delta_i^{-1} \left(\frac{y_i - \mu_i}{\phi v(\mu_i)} \right)\right] \\
&= \frac{\partial}{\partial \eta_i}\left[\left(\frac{\partial \mu_i}{\partial \eta_i}\right) \left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right)\right] & \left(\text{see proof of Eq. \eqref{eq:deriv-1-eta}}\right) \\
&= \frac{\partial}{\partial \eta_i} \left[ \frac{\partial \mu_i}{\partial \eta_i}\right]\left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right) + \left(\frac{\partial \mu_i}{\partial \eta_i}\right) \frac{\partial}{\partial \eta_i}\left[\frac{y_i - \mu_i}{\phi v(\mu_i)}\right] & \left(\text{product rule}\right)\\
&= \left(\frac{\partial^2 \mu_i}{\partial \eta_i^2}\right)\left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right) + \left(\frac{\partial \mu_i}{\partial \eta_i}\right)^2 \frac{\partial}{\partial \mu_i}\left[\frac{y_i - \mu_i}{\phi v(\mu_i)}\right] & \left(\text{chain rule}\right)\\
&= \left(\frac{\partial^2 \mu_i}{\partial \eta_i^2}\right)\left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right) + \left(\frac{\partial \mu_i}{\partial \eta_i}\right)^2 \left(\frac{- \phi v(\mu_i) - \phi (y_i - \mu_i)\frac{\partial v(\mu_i)}{\partial \mu_i}}{(\phi v(\mu_i))^2}\right) & \left(\text{quotient rule}\right)\\
&= - \left(\frac{\partial \eta_i}{\partial \mu_i}\right)^3 \left(\frac{\partial^2 \eta_i}{\partial \mu_i^2}\right)\left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right) + \left(\frac{\partial \mu_i}{\partial \eta_i}\right)^2 \left(\frac{- \phi v(\mu_i) - \phi (y_i - \mu_i)\frac{\partial v(\mu_i)}{\partial \mu_i}}{(\phi v(\mu_i))^2}\right) & \left(\text{Eq. \eqref{eq:temp}}\right) \\
&= - \frac{\frac{\partial^2 \eta_i}{\partial \mu_i^2} \left( \frac{y_i - \mu_i}{\phi v(\mu_i)}\right)}{\left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{3}}  + \frac{-\phi v(\mu_i) - \phi \frac{\partial v(\mu_i)}{\partial \mu_i}(y_i - \mu_i)}{(\phi v(\mu_i))^2 \left( \frac{\partial \eta_i}{\partial \mu_i}\right)^2 } & \left( \frac{\partial \mu_i}{\partial \eta_i} = \left(\frac{\partial \eta_i}{\partial \mu_i} \right)^{-1} \text{ by Eq. \eqref{eq:deriv-assumption}}\right) \\
&= -(y_i - \mu_i) \left(\frac{\frac{\partial^2 \eta_i}{\partial \mu_i^2}}{\phi v(\mu_i) \left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{3}}\right) - \frac{1}{\phi v(\mu_i) \left(\frac{\partial \eta_i}{\partial \mu_i} \right)^2} - (y_i - \mu_i) \left(\frac{\phi \frac{\partial v(\mu_i)}{\partial \mu_i}}{(\phi v(\mu_i))^2 \left( \frac{\partial \eta_i}{\partial \mu_i}\right)^2}\right) \\
&= -(y_i - \mu_i) \left[\frac{\frac{\partial^2 \eta_i}{\partial \mu_i^2}}{\phi v(\mu_i) \left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{3}}+ \frac{\phi \frac{\partial v(\mu_i)}{\partial \mu_i}}{(\phi v(\mu_i))^2 \left( \frac{\partial \eta_i}{\partial \mu_i}\right)^2} \right] - \frac{1}{\phi v(\mu_i) \left(\frac{\partial \eta_i}{\partial \mu_i} \right)^2}\\
&= -(y_i - \mu_i) \left[\frac{\phi v(\mu_i) \frac{\partial^2 \eta_i}{\partial \mu_i^2}}{(\phi v(\mu_i))^2 \left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{3}}+ \frac{\phi \left(\frac{\partial v(\mu_i)}{\partial \mu_i}\right)\left(\frac{\partial \eta_i}{\partial \mu_i}\right)}{(\phi v(\mu_i))^2 \left( \frac{\partial \eta_i}{\partial \mu_i}\right)^3} \right] - \omega_i & \left(\text{Eq. \eqref{eq:deriv-assumption}}\right)\\
&= -(y_i - \mu_i) \left[\frac{\phi v(\mu_i) \frac{\partial^2 \eta_i}{\partial \mu_i^2} + \phi \left(\frac{\partial v(\mu_i)}{\partial \mu_i}\right)\left(\frac{\partial \eta_i}{\partial \mu_i}\right)}{(\phi v(\mu_i))^2 \left( \frac{\partial \eta_i}{\partial \mu_i}\right)^3} \right] - \omega_i & \left(\text{Eq. \eqref{eq:deriv-assumption}}\right)\\
&= -\omega_i - e_i(y_i - \mu_i) \\
&= -\xi_i
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

Let's also define the following quantities:

$$
\begin{aligned}
    \tilde{\mathbf{y}} &= \bar{\Omega} \bar{\Delta}^{-1} (\mathbf{y} - \bar{\mu});
    &
    \mathbf{M} &= \mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top;
    &
    (\mathbf{A}_l) &= \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top
\end{aligned}
$$

We can rewrite the marginal log quasi-likelihood as:

$$
\begin{equation}
\label{eq:rewrite-mlql}
\begin{aligned}
    \ell_q(\mathbf{y}; \theta) &\approx \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) + \frac{1}{2} \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} - \frac{1}{2} \text{tr}\left[ \bar{\Xi} \mathbf{M} \right]
\end{aligned}
\end{equation}
$$

<div class="example">
    <strong>Example.</strong>
    <br>
    For a negative binomial response,  we have the following identities:
    $$
    \begin{aligned}
    g(\mu^k_i) &= \eta^k_i    &    V(\mu^k_i) &= v(\mu_i^k) = \mu^k_i + \frac{1}{\gamma} (\mu^k_i)^2 \\
    \frac{\partial g(\mu^k_i)}{\partial \mu_i^k} &= \frac{1}{\mu_i^k}    &    \frac{\partial V(\mu^k_i)}{\partial \mu_i^k} &= 1 + \frac{2}{\gamma} \mu_i^k \\
    \frac{\partial^2 g(\mu^k_i)}{\partial (\mu_i^k)^2} &= -\frac{1}{(\mu_i^k)^2}    &    \frac{\partial^2 V(\mu^k_i)}{\partial (\mu_i^k)^2} &= \frac{2}{\gamma}
    \end{aligned}
    $$
    Thus, the intermediate terms are:
    $$
    \begin{aligned}
    \delta_i^k &= \mu_i^k   &    \omega_i^k &= \frac{\mu_i^k \gamma}{\mu_i^k + \gamma} \\
    e_i^k &= \frac{\mu_i^k \gamma}{(\mu_i^k + \gamma)^2}   &  
    \xi_i^k &= \frac{\mu_i^k \gamma}{\mu_i^k + \gamma} \left( \frac{y_i^k - \mu_i^k}{\mu_i^k + \gamma}\right)
    \end{aligned}
    $$
</div>


---

## Score

The score vector, $U(\theta) = (U^\top_{\alpha}(\theta), U_{\phi}(\theta), U^\top_{\tau^2}(\theta))^\top$, has components:

$$
\begin{align}
U_{\alpha_j}(\theta) 
    &= \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \alpha_j} + \frac{1}{2} \frac{\partial}{\partial \alpha_j} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right] - \frac{1}{2}  \frac{\partial}{\partial \alpha_j} \left[ \text{tr}\left[ \bar{\Xi} \mathbf{M} \right] \right] \label{eq:score-alpha-j} \\
U_{\phi}(\theta) 
    &= \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \phi} + \frac{1}{2} \frac{\partial}{\partial \phi} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right] - \frac{1}{2}  \frac{\partial}{\partial \phi} \left[ \text{tr}\left[ \bar{\Xi} \mathbf{M} \right] \right] \label{eq:score-phi} \\ 
U_{\tau^2}(\theta) 
    &= \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \tau^2} + \frac{1}{2} \frac{\partial}{\partial \tau^2} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right] - \frac{1}{2} \frac{\partial}{\partial \tau^2} \left[ \text{tr}\left[  \bar{\Xi} \mathbf{M} \right] \right] \label{eq:score-tau-sqr} 
\end{align}
$$

We'll start with the partial derivatives of the conditional log quasi-likelihood (the first term in each equation above).

<!-- #region deriv-1-clql -->
{% tabs deriv-1-clql %}
{% tab deriv-1-clql statement %}
$$
\begin{equation}
\begin{aligned}
\frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \alpha_j}
  &= \tilde{\mathbf{y}}^\top \mathbf{X}_{\cdot, j} \\
\frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \phi} 
  &= - \frac{1}{\phi} \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) \\
\frac{\partial \ell_q(\mathbf{y} \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \tau^2_l}
  &= 0 
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-1-clql proof %}
$$
\begin{aligned}
\frac{\partial \ell_q(y_i; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \alpha_j}
    &= \frac{\partial \bar{\eta}_i}{\partial \alpha_j}\frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial}{\partial \bar{\mu}_i} \left[ \int_{y_i}^{\bar{\mu}_i} \frac{y_i - u}{V(u)} du \right] \\
    &= \mathbf{x}_{i,j} \left(\frac{\partial \bar{\eta}_i}{\partial \bar{\mu}_i}\right)^{-1} \left(\frac{y_i - \bar{\mu}_i}{V(\bar{\mu}_i)}\right) & \left(\text{Eq. } \eqref{eq:deriv-assumption} \right) \\
    &= \mathbf{x}_{i,j} \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^{-1} \left(\frac{y_i - \bar{\mu}_i}{V(\bar{\mu}_i)}\right) & \left(\text{FTC}\right)\\
    &= \mathbf{x}_{i,j} (y_i - \bar{\mu}_i) \frac{\bar{\delta}_i}{V(\bar{\mu}_i)} \\
    &= \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i) \\
    &= \mathbf{x}_{i,j} \tilde{y}_i \\
\implies \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \alpha_j}
    &= \tilde{\mathbf{y}}^\top \mathbf{X}_{\cdot, j} \\
\frac{\partial \ell_q(y_i; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \phi}
    &= \frac{\partial}{\partial \phi} \left[ \int_{y_i}^{\bar{\mu}_i} \frac{y_i - u}{V(u)} du \right] \\
    &= \frac{\partial}{\partial \phi} \left[ \int_{y_i}^{\bar{\mu}_i} \frac{y_i - u}{\phi v(u)} du \right] \\
    &= \frac{\partial}{\partial \phi} \left[ \frac{1}{\phi} \right]\int_{y_i}^{\bar{\mu}_i} \frac{y_i - u}{v(u)} du  \\
    &= - \frac{1}{\phi^2} \int_{y_i}^{\bar{\mu}_i} \frac{y_i - u}{v(u)} du \\
    &= - \frac{1}{\phi} \ell_q(y_i; \theta \rvert \beta = \mathbf{0}_{m \times k}) \\
\implies \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \phi}
  &= - \frac{1}{\phi} \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) \\
\frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k})}{\partial \tau^2_l} 
    &= 0
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

Next, we compute the partial derivatives of the quadratic form (middle term).

<!-- #region deriv-1-quad -->
{% tabs deriv-1-quad %}
{% tab deriv-1-quad statement %}
$$
\begin{equation}
\label{eq:deriv-1-quad}
\begin{aligned}
\frac{\partial}{\partial \alpha_j} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right]
    &=   -2 \left( \bar{\Xi} \mathbf{X}_{\cdot, j} \right)^\top \left(\mathbf{M} \tilde{y} \right) \\
\frac{\partial}{\partial \phi} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right]
    &= - \frac{2}{\phi} \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \\
\frac{\partial}{\partial \tau^2_l} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right] 
    &=  \tilde{\mathbf{y}}^\top (\mathbf{A}_l) \tilde{\mathbf{y}} 
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-1-quad proof %}
First, we will find the partial derivatives of $\mathbf{M}$. 

$$
\begin{aligned}
\frac{\partial \mathbf{M}_{i,i'}}{\partial \alpha_j}
    &= \frac{\partial}{\partial \alpha_j} \left[ \left(\mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top \right)_{i, i'} \right] \\
    &= 0 \\
\frac{\partial \mathbf{M}_{i,i'}}{\partial \phi}
    &= \frac{\partial}{\partial \phi} \left[ \left(\mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top \right)_{i, i'} \right] \\ 
    &= 0 \\
\frac{\partial \mathbf{M}_{i,i'}}{\partial \tau^2_l}
    &= \frac{\partial}{\partial \tau^2_l} \left[ \left(\mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top \right)_{i, i'} \right] \\
    &= \left( \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top\right)_{i,i'} \\
    &= \left( \mathbf{Z} (\mathbf{A}_l) \mathbf{Z}^\top\right)_{i,i'}
\end{aligned}
$$

Next, we'll find the partial derivatives of $\tilde{\mathbf{y}}$:

$$
\begin{aligned}
\frac{\partial \tilde{y}_i}{\partial \alpha_j}
    &= \frac{\partial}{\partial \alpha_j} \left[ \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i) \right] \\
    &= (y_i - \bar{\mu}_i) \frac{\partial}{\partial \alpha_j} \left[ \frac{\bar{\delta}_i}{V(\bar{\mu}_i)}\right] + \bar{\omega}_i \bar{\delta}_i^{-1} \frac{\partial}{\partial \alpha_j} \left[ y_i - \bar{\mu}_i \right] \\
    &= (y_i - \bar{\mu}_i) \left( \frac{1}{V(\bar{\mu}_i)} \frac{\partial}{\partial \alpha_j} \left[ \left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^{-1} \right] + \bar{\delta}_i \frac{\partial}{\partial \alpha_j} \left[ (V(\bar{\mu}_i))^{-1} \right] \right) + \bar{\omega}_i \bar{\delta}_i^{-1} \left(- \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \right) \\
    &= (y_i - \bar{\mu}_i) \left(\frac{1}{V(\bar{\mu}_i)} \left(- \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^{-2} \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} \right) + \bar{\delta}_i \left(- \frac{1}{V^2(\bar{\mu}_i)} \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial V(\bar{\mu})}{\partial \bar{\mu}_i} \right) \right) - \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} \left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^{-1} \\
    &= (y_i - \bar{\mu}_i) \left(- \frac{\mathbf{x}_{i,j} \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^{-1}\frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{V(\bar{\mu}_i)\left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^2} - \frac{\mathbf{x}_{i,j} \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^{-1} \frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i}}{V^2(\bar{\mu}_i) \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}} \right) - \mathbf{x}_{i,j} \bar{\omega}_i & \left(\text{Eq. } \eqref{eq:deriv-assumption} \right) \\
    &= - \mathbf{x}_{i,j} \left[ (y_i - \bar{\mu}_i) \left( \frac{V(\bar{\mu}_i) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{V^2(\bar{\mu}_i) \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^3} + \frac{\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i}}{V^2(\bar{\mu}_i) \left(\frac{ \partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^3} \right) + \bar{\omega}_i  \right]\\
    &=-\mathbf{x}_{i,j} \left(\bar{\omega}_i  + \bar{e}_i (y_i - \bar{\mu}_i) \right) \\
    &= -\mathbf{x}_{i,j} \bar{\xi}_i \\
\frac{\partial \tilde{y}_i}{\partial \phi}
    &= \frac{\partial}{\partial \phi} \left[ \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i) \right] \\ 
    &= (y_i - \bar{\mu}_i) \frac{\partial}{\partial \phi} \left[ \frac{\bar{\delta}_i}{\phi v(\bar{\mu}_i)}\right] \\
    &= - \frac{\bar{\delta}_i (y_i - \bar{\mu}_i)}{\phi^2 v(\bar{\mu}_i) } \\
    &= - \frac{1}{\phi} \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i)  \\
    &= - \frac{\tilde{y}_i}{\phi} \\
\frac{\partial \tilde{y}_i}{\partial \tau^2}
    &= \frac{\partial}{\partial \tau^2} \left[ \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i) \right] \\ 
    &= 0 
\end{aligned}
$$

And thus:

$$
\begin{aligned}
\frac{\partial}{\partial \alpha_j} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right]
    &= \sum_{i = 1}^N \sum_{i' = 1}^N  \mathbf{M}_{i,i'} \frac{\partial}{\partial \alpha_j} \left[ \tilde{y}_i \tilde{y}_{i'}  \right] \\
    &= \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \left(\tilde{y}_{i'} \frac{\partial \tilde{y}_i}{\partial \alpha_j} + \tilde{y}_i \frac{\partial \tilde{y}_{i'}}{\partial \alpha_j} \right) \\
    &= \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \left(\tilde{y}_{i'} \left(- x_{i,j} \bar{\xi}_i \right) + \tilde{y}_i \left( - x_{i', j} \bar{\xi}_i' \right) \right) \\
    &= -\sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \tilde{y}_{i'} x_{i,j} \bar{\xi}_i - \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \tilde{y}_i x_{i', j} \bar{\xi}_i' \\
    &= - 2 \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \tilde{y}_{i'} x_{i,j} \bar{\xi}_i & \left(\mathbf{M} \text{ sym.} \right)\\
    &= - 2 \sum_{i = 1}^N  x_{i,j} \bar{\xi}_i  \sum_{i' = 1}^N \mathbf{M}_{i,i'} \tilde{y}_{i'} \\
    &= -2 \sum_{i = 1}^N x_{i,j} \bar{\xi}_i \left(\mathbf{M}^\top \tilde{y} \right)_i \\
    &= -2 \left( \bar{\Xi} \mathbf{X}_{\cdot, j} \right)^\top \left(\mathbf{M} \tilde{y} \right) & \left(\mathbf{M} \text{ sym.}\right) \\
\frac{\partial}{\partial \phi} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right]
    &= 2 \tilde{\mathbf{y}}^\top \mathbf{M} \frac{\partial \tilde{\mathbf{y}}}{\partial \phi} \\
    &= - \frac{2}{\phi} \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \\
\frac{\partial}{\partial \tau^2_l} \left[ \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} \right] 
    &= \tilde{\mathbf{y}}^\top \frac{\partial \mathbf{M}}{\partial \tau^2_l} \tilde{\mathbf{y}} \\
    &=  \tilde{\mathbf{y}}^\top \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top \tilde{\mathbf{y}}  \\
    &= \tilde{\mathbf{y}}^\top (\mathbf{A}_l) \tilde{\mathbf{y}}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

Finally, we can take the derivatives of the trace (last term):

<!-- #region deriv-1-trace -->
{% tabs deriv-1-trace %}
{% tab deriv-1-trace statement %}
$$
\begin{equation}
\begin{aligned}
\frac{\partial}{\partial \alpha_j} \left[ \text{tr}\left[ \bar{\Xi} \mathbf{M} \right] \right] 
    &= \sum_{i = 1}^N \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - 3 \bar{\omega}_i \bar{\delta}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} \right)  \\
\frac{\partial}{\partial \phi}\left[ \text{tr} \left[ \bar{\Xi} \mathbf{M} \right]  \right]
    &= - \frac{1}{\phi} \sum_{i = 1}^N \mathbf{M}_{i,i} \bar{\xi}_i \\
\frac{\partial}{\partial \tau^2} \left[ \text{tr} \left[ \bar{\Xi} \mathbf{M} \right]  \right] 
    &= \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \bar{\xi}_i 
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-1-trace proof %}
$$
\begin{aligned}
\frac{\partial \bar{\xi}_i}{\partial \alpha_j} 
    &= \frac{\partial}{\partial \alpha_j} \left[ \bar{\omega}_i + \bar{e}_i(y_i - \bar{\mu}_i) \right] \\
    &= \frac{\partial}{\partial \alpha_j} \left[ \frac{\bar{\delta}_i^2}{V(\bar{\mu}_i)} \right] + y_i \frac{\partial}{\partial \alpha_j} \left[ \bar{e}_i\right] - \frac{\partial}{\partial \alpha_j} \left[ \bar{e}_i \bar{\mu}_i \right] \\
    &= \frac{2 \bar{\delta}_i}{V(\bar{\mu}_i)} \frac{\partial \bar{\delta}_i}{\partial \alpha_j} + \frac{\bar{\delta}_i^2}{V^2(\bar{\mu}_i)} \frac{\partial V(\bar{\mu}_i)}{\partial \alpha_j} + y_i \frac{\partial}{\partial \alpha_j} \left[ \bar{e}_i\right] - \frac{\partial}{\partial \alpha_j} \left[ \bar{e}_i \bar{\mu}_i \right] \\ 
    &= \frac{2 \bar{\delta}_i}{V(\bar{\mu}_i)} \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial \bar{\delta}_i}{\partial \bar{\mu}_i} + \frac{\bar{\delta}_i^2}{V^2(\bar{\mu}_i)} \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i} + y_i \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} -  \frac{\partial \bar{\eta}_i}{\partial \alpha_j} \frac{\partial \bar{\mu}_i}{\partial \bar{\eta}_i} \frac{\partial}{\partial \bar{\mu}_i}\left[ \bar{e}_i \bar{\mu}_i \right] \\ 
    &= \mathbf{x}_{i,j} \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^{-1} \left( \frac{2 \bar{\delta}_i}{V(\bar{\mu}_i)}  \frac{\partial \bar{\delta}_i}{\partial \bar{\mu}_i} + \frac{\bar{\delta}_i^2}{V^2(\bar{\mu}_i)}\frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i} + y_i \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \frac{\partial}{\partial \bar{\mu}_i}\left[ \bar{e}_i \bar{\mu}_i \right] \right) & \left(\text{Eq. } \eqref{eq:deriv-assumption}\right)\\  
    &= \mathbf{x}_{i,j} \bar{\delta}_i \left( \frac{2 \bar{\delta}_i}{V(\bar{\mu}_i)}  \frac{\partial}{\partial \bar{\mu}_i} \left[ \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^{-1} \right]+ \frac{\bar{\delta}_i^2}{V^2(\bar{\mu}_i)}\frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i} + y_i \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \frac{\partial}{\partial \bar{\mu}_i}\left[ \bar{e}_i \bar{\mu}_i \right] \right) & \left(\text{Eq. } \eqref{eq:deriv-assumption}\right) \\  
     &= \mathbf{x}_{i,j} \bar{\delta}_i \left(- \frac{2 \bar{\delta}_i}{V(\bar{\mu}_i)\left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^2} \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} + \frac{\bar{\delta}_i^2}{V^2(\bar{\mu}_i)}\frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i} + y_i \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \frac{\partial}{\partial \bar{\mu}_i}\left[ \bar{e}_i \bar{\mu}_i \right] \right) & \left(\text{Eq. } \eqref{eq:deriv-assumption}\right) \\  
     &= \mathbf{x}_{i,j} \bar{\delta}_i \left(- \frac{2 V(\bar{\mu}_i) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{V^2(\bar{\mu}_i)\left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^3} + \frac{ \frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i}\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}}{V^2(\bar{\mu}_i)\left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^3}  + (y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \bar{e}_i \right) \\
     &= \mathbf{x}_{i,j} \bar{\delta}_i \left(- \frac{3 V(\bar{\mu}_i) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{V^2(\bar{\mu}_i)\left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^3} + \frac{ V(\bar{\mu}_i) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} + \frac{\partial V(\bar{\mu}_i)}{\partial \bar{\mu}_i}\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}}{V^2(\bar{\mu}_i)\left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^3}  + (y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \bar{e}_i \right) \\
     &= \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \frac{3 \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{V(\bar{\mu}_i)\left( \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} \right)^3} \right) \\
     &= \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - \frac{3 \bar{\delta}_i^3 \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{V(\bar{\mu}_i)} \right) \\
     &= \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - 3 \bar{\omega}_i \bar{\delta}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} \right) \\ 
\implies \frac{\partial}{\partial \alpha_j} \text{tr}\left[ \bar{\Xi} \mathbf{M} \right]
    &= \sum_{i = 1}^N \sum_{j = 1}^N \frac{\partial}{\partial \alpha_j} \left[ \bar{\Xi}_{i,j} \right] \mathbf{M}_{i,j} \\
    &= \sum_{i = 1}^N\frac{\partial}{\partial \alpha_j} \left[ \bar{\xi}_i \right] \mathbf{M}_{i,i} & \left(\bar{\Xi} \text{ diag.}\right) \\
    &= \sum_{i = 1}^N \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial\bar{e}_i}{\partial \bar{\mu}_i} - 3 \bar{\omega}_i \bar{\delta}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} \right) \\
\frac{\partial \bar{\xi}_i}{\partial \phi}
    &= \frac{\partial}{\partial \phi} \left[ \bar{\omega}_i + \bar{e}_i(y_i - \bar{\mu}_i) \right] \\
    &= \frac{\partial \bar{\omega}_i}{\partial \phi} + (y_i - \bar{\mu}_i) \frac{\partial \bar{e}_i}{\partial \phi} \\
    &= \frac{\partial}{\partial \phi} \left[\frac{\bar{\delta}_i^2}{\phi v(\bar{\mu}_i) } \right] + (y_i - \bar{\mu}_i) \frac{\partial}{\partial \phi} \left[\frac{\phi \frac{\partial v(\bar{\mu}_i)}{\partial \bar{\mu}_i} \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} + \phi v(\bar{\mu_i}) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{\phi^2 v^2(\bar{\mu}_i) \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^3 } \right] \\
    &= - \frac{1}{\phi^2} \frac{\bar{\delta}^2}{v(\bar{\mu}_i)} - \frac{y_i - \bar{\mu}_i}{\phi^2} \frac{\frac{\partial v(\bar{\mu}_i)}{\partial \bar{\mu}_i} \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} + v(\bar{\mu_i}) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{v^2(\bar{\mu}_i) \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^3 } \\
    &= - \frac{\bar{\omega}_i}{\phi} - \frac{y_i - \bar{\mu}_i}{\phi^2} \frac{\phi \frac{\partial v(\bar{\mu}_i)}{\partial \bar{\mu}_i} \frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i} + \phi v(\bar{\mu_i}) \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}}{\phi v^2(\bar{\mu}_i) \left(\frac{\partial g(\bar{\mu}_i)}{\partial \bar{\mu}_i}\right)^3 } \\
    &= - \frac{1}{\phi}\left(\bar{\omega}_i + \bar{e}_i (y_i - \bar{\mu}_i)\right) \\
    &= - \frac{\bar{\xi}_i}{\phi} \\
\implies \frac{\partial}{\partial \phi} \left[ \text{tr}\left[ \bar{\Xi} \mathbf{M} \right] \right]
    &= - \frac{1}{\phi} \text{tr}\left[ \bar{\Xi} \mathbf{M} \right] \\
    &= - \frac{1}{\phi} \sum_{i = 1}^N \sum_{j = 1}^N \bar{\Xi}_{i,j} \mathbf{M}_{i,j} \\
    &= - \frac{1}{\phi} \sum_{i = 1}^N \mathbf{M}_{i,i} \bar{\xi}_i  & \left(\bar{\Xi} \text{ diag.}\right)\\
\frac{\partial}{\partial \tau^2_l} \left[ \bar{\Xi} \mathbf{M} \right] 
    &= \bar{\Xi} \frac{\partial \mathbf{M}}{\partial \tau^2} \\
    &= \bar{\Xi} \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2} \mathbf{Z}^\top \\
    &= \bar{\Xi} (\mathbf{A}_l) \\
\implies \frac{\partial}{\partial \tau^2_l} \left[ \text{tr}\left[ \bar{\Xi} \mathbf{M} \right] \right] 
    &= \text{tr}\left[ \bar{\Xi} (\mathbf{A}_l) \right] \\
    &= \sum_{i = 1}^N \sum_{j = 1}^N \bar{\Xi}_{i,j} (\mathbf{A}_l)_{i,j} \\
    &= \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \bar{\xi}_i  & \left(\bar{\Xi} \text{ diag.}\right)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->


Let's define the follow terms, which do not depend on the responses:

$$
\begin{aligned}
\rho_{i,j} &= \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} - \frac{1}{2} \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \frac{\partial \bar{e}_i}{\partial \bar{\mu}_i} + \mathbf{x}_{i,j} \bar{\xi}_i \bar{\omega}_i \bar{\delta}_i^{-1} \left(\sum_{i' = 1}^N \mathbf{M}_{i,i'} \right) 
\\
\sigma_{i,j} &=  \frac{3}{2} \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i^2 \bar{\omega}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}
\end{aligned}
$$

<!-- #region score-rewrite-2 -->
{% tabs score-rewrite-2 %}
{% tab score-rewrite-2 statement %}
The score components can be simplified to:

$$
\begin{align}
U_{\alpha_j}(\theta) 
&= \sum_{i = 1}^N \left[ (y_i - \bar{\mu}_i) \rho_{i,j} + \sigma_{i,j} \right] \label{eq:score-alpha-j-2}\\
U_\phi(\theta)
&= \frac{1}{\phi} \left( - \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) - \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \bar{\omega}_i \bar{\omega}_{i'} \bar{\delta}_i^{-1} \bar{\delta}_{i'}^{-1} (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) + \frac{1}{2} \sum_{i = 1}^N \mathbf{M}_{i,i} \bar{\xi}_i \right) \label{eq:score-phi-2} \\
U_{\tau^2_l}(\theta)
&= \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i,i'} \bar{\omega}_i \bar{\omega}_{i'} \bar{\delta}_i^{-1} \bar{\delta}_{i'}^{-1} (y_i - \bar{\mu}_i) (y_{i'} - \bar{\mu}_{i'}) - \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \bar{\xi}_i \label{eq:score-tau-sqr-2}
\end{align}
$$
{% endtab %}
{% tab score-rewrite-2 proof %}
For the fixed effects:

$$
\begin{aligned}
U_{\alpha_j}(\theta) 
&= \tilde{\mathbf{y}}^\top \mathbf{X}_{\cdot, j} + \frac{2}{2} \left(\bar{\Xi} \mathbf{X}_{\cdot, j} \right)^\top (\mathbf{M} \tilde{\mathbf{y}}) - \frac{1}{2} \sum_{i = 1}^N \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial \bar{e}_i}{\partial \bar{\mu}_i} - 3 \bar{\omega}_i \bar{\delta}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2} \right) \\
&= \sum_{i = 1}^N \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i)
+ \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \mathbf{x}_{i,j} \bar{\xi}_i \bar{\omega}_i \bar{\delta}_i^{-1} (y_i - \bar{\mu}_i)
- \frac{1}{2} \sum_{i = 1}^N \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \left((y_i - \bar{\mu}_i) \frac{\partial \bar{e}_i}{\partial \bar{\mu}_i} - 3 \bar{\omega}_i \bar{\delta}_i \frac{\partial^2 g(\bar{mu}_i)}{\partial \bar{\mu}_i^2} \right) \\
&= \sum_{i = 1}^N \left[ (y_i - \bar{\mu}_i) \left(\mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} - \frac{1}{2} \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \frac{\partial \bar{e}_i}{\partial \bar{\mu}_i} + \mathbf{x}_{i,j} \bar{\xi}_i \bar{\omega}_i \bar{\delta}_i^{-1} \left(\sum_{i' = 1}^N \mathbf{M}_{i,i'} \right)\right) + \frac{3}{2} \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i^2 \bar{\omega}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}\right]  \\
&= \sum_{i = 1}^N \left[ (y_i - \bar{\mu}_i) \rho_{i,j} + \sigma_{i,j} \right]
\end{aligned}
$$

For $\phi$:

$$
\begin{aligned}
U_{\phi}(\theta)
&= - \frac{1}{\phi} \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) - \frac{2}{2 \phi} \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} + \frac{1}{2\phi} \sum_{i = 1}^N \mathbf{M}_{i,i} \bar{\xi}_i \\
&= \frac{1}{\phi} \left( - \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) - \sum_{i = 1}^N \sum_{i' = 1}^N \mathbf{M}_{i,i'} \bar{\omega}_i \bar{\omega}_{i'} \bar{\delta}_i^{-1} \bar{\delta}_{i'}^{-1} (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) + \frac{1}{2} \sum_{i = 1}^N \mathbf{M}_{i,i} \bar{\xi}_i \right)
\end{aligned} 
$$

For $\tau^2$:

$$
\begin{aligned}
U_{\tau^2_l}(\theta)
&= 0 + \frac{1}{2} \tilde{\mathbf{y}}^\top (\mathbf{A}_l) \tilde{\mathbf{y}} - \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \bar{\xi}_i \\
&= \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i,i'} \bar{\omega}_i \bar{\omega}_{i'} \bar{\delta}_i^{-1} \bar{\delta}_{i'}^{-1} (y_i - \bar{\mu}_i) (y_{i'} - \bar{\mu}_{i'}) - \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \bar{\xi}_i
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->


<div class="example">
    <strong>Example.</strong>
    <br>
    Our example is a bit different because the hypothesis test concerns a variance component <i>and</i> a fixed effect. The variance component and $\alpha_1$ pose no problems; $U_{\tau^2}(\boldsymbol{\theta})$ is exactly as in Eq. \eqref{eq:score-tau-sqr-2}, and $U_{\alpha_1}(\boldsymbol{\theta})$ is given in Eq. \eqref{eq:score-alpha-j-2}. 
    <br>
    However, the intercept terms are cluster-specific, so the derivative of the marginal likelihood w.r.t to a particular $\alpha_k^0$ will only involve terms from cluster $k$. Since these terms are intercepts, the corresponding covariate values are $0$, so we can drop $\mathbf{x}$ from our expressions. Thus, $U_{\alpha_k^0}(\boldsymbol{\theta})$ is as in Eq. \eqref{eq:score-alpha-j-2} with:

    $$
    \begin{aligned}
    &\rho_{i, 0}^{k'} = 0 &
    &\rho_{i, 0}^k = \bar{\omega}_i^k \left(\bar{\delta}_i^k\right)^{-1} - \frac{1}{2} \left[ \mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top \right]_{i,i} \bar{\delta}_i^k \frac{\partial \bar{e}_i^k}{\partial \bar{\mu}_i^k} + \bar{\xi}_i^k \bar{\omega}_i^k \left( \bar{\delta}_i^j\right)^{-1} \left( \sum_{i' = 1}^{n_k} \left[ \mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top\right]_{i,i'} \right) \\
    &\sigma_{i,0}^{k'} = 0 &
    &\sigma_{i,0}^k = \frac{3}{2} \left[ \mathbf{Z} \mathbf{D}(\tau^2) \mathbf{Z}^\top \right]_{i,i} \left(\bar{\delta}_i^k\right)^2 \bar{\omega}_i^k \frac{\partial^2 g(\bar{\mu}_i^k)}{\partial \left(\bar{\mu}_i^k\right)^2}
    \end{aligned}
    $$

    Since we assume that $\mathbf{D}(\tau^2) = \mathbf{0}_{q \times q}$ when $\tau^2 = \mathbf{0}_q$, the non-zero components will simplify to $\rho_{i,0}^k = \bar{\omega}_i^k \left(\bar{\delta}_i^k\right)^{-1}$ and $\sigma_{i,0}^k = 0$.
</div>


---

## Information

With expectations taken with respect to $\mathbf{y}$ under $H_0$, the information matrix is a function of $\theta$ and has blocks:

$$
\begin{aligned}
\mathcal{I}(\theta) 
    &= \mathbb{E}\left[ U(\theta) U^\top (\theta) \right] \\
    &= \begin{bmatrix}
        \mathbb{E}\left[ U_{\alpha}(\theta) U^\top_{\alpha}(\theta) \right] & \mathbb{E}\left[ U_{\alpha}(\theta) U_\phi(\theta) \right] & \mathbb{E}\left[ U_\alpha(\theta) U_{\tau^2}(\theta) \right] \\
        \mathbb{E}\left[  U_{\phi}(\theta) U^\top_{\alpha}(\theta) \right] & \mathbb{E}\left[ U_\phi(\theta)U_\phi(\theta) \right] & \mathbb{E}\left[ U_\phi(\theta)U_{\tau^2}(\theta) \right] \\
        \mathbb{E}\left[ U_{\tau^2}(\theta)U^\top_\alpha(\theta) \right] & \mathbb{E}\left[ U_{\phi}(\theta) U_{\tau^2}(\theta) \right] & \mathbb{E}\left[ U_{\tau^2}(\theta) U_{\tau^2}(\theta)\right]
    \end{bmatrix}
\end{aligned}
$$

### Intermediate Quantities
We will let $m_{r,i}$ denote the $r$-th central moment of $y_i$ and $\kappa_{r, i}$ the $r$-th cumulant of $y_i$. Recall that the second and third cumulants equal the second and third central moments, respectively. Furthermore, we have (as a fact) that $\kappa_{4,i} = m_{4,i} - 3m_{2,i}^2$ and that:

$$
\begin{equation}
\label{eq:omega-facts}
\begin{aligned}
\bar{\delta}_i^{-2} m_{2,i}
  &= \frac{V(\bar{\mu}_i)}{\bar{\delta}_i^2}
  = \bar{\omega}_i^{-1} \\
\bar{\omega}_i \bar{\delta}_i^{-1} m_{2,i}
  &= \frac{\bar{\delta}^2}{V(\bar{\mu}_i)} \bar{\delta}_i^{-1} V(\bar{\mu}_i) 
  = \bar{\delta}_i
\end{aligned}
\end{equation}
$$

If we assume that the following relationship holds for $r = 2$ and $r = 3$:

$$
\begin{equation}
\label{eq:kappa-assumption}
\begin{aligned}
\kappa_{2, i} &= \phi v(\mu_i) &
\kappa_{(r+1), i} &= \kappa_{2, i} \frac{\partial \kappa_{r, i}}{\partial \mu_i}
\end{aligned}
\end{equation}
$$

then we have: 

$$
\begin{aligned}
\kappa_{3, i} &= \phi^2 v(\mu_i) \frac{\partial v(\mu_i)}{\partial\mu_i} &
\kappa_{4, i} &= \phi^3 v(\mu_i) \left( v(\mu_i) \frac{\partial^2 v(\mu_i)}{\partial\mu_i^2} + \left(\frac{\partial v(\mu_i)}{\partial \mu_i} \right)^2 \right)
\end{aligned}
$$


<!-- 
<div class="example">
  <strong>Example.</strong>
  <br>
  $$
  \begin{aligned}
  \kappa_{2, i} &= \phi \left(\mu_i + \frac{1}{\gamma} \mu_i^2 \right) \\
  \kappa_{3, i} &= \phi^2 \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right)\left(1 + \frac{2}{\gamma} \mu_i \right) \\
  &= \phi^2 \mu_i \left(1 + \frac{1}{\gamma} \mu_i \right)\left(1 + \frac{2}{\gamma} \mu_i \right)  \\
  \kappa_{4, i} &= \phi^3 \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right) \left( \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right) \frac{\partial}{\partial \mu_i} \left[1 + \frac{2}{\gamma} \mu_i \right] + \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right)\\
  &= \phi^3\left[ \frac{2}{\gamma} \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right)^2 + \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right) \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right]
  \end{aligned}
  $$
</div> -->


We will also let $\mathbf{A}$ and $\mathbf{B}$ be $N \times N$ symmetrical matrices that do not depend on $\mathbf{y}$. We have the following results.

<!-- #region y-tilde-exp -->
{% tabs y-tilde-exp %}
{% tab y-tilde-exp statement %}

$$
\begin{align}
\sum_{\substack{i = 1 \\ j = 1 \\ i' = 1 \\ j' = 1}}^N \mathbf{A}_{i,j} \mathbf{B}_{i',j'} \mathbb{E}\left[ \tilde{y}_i \tilde{y}_{j} \tilde{y}_{i'} \tilde{y}_{j'} \right] 
&= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
    + \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_j 
    + 2  \sum_{i = 1}^N\sum_{j = 1}^N \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i \bar{\omega}_j \label{eq:y-four} \\
\sum_{\substack{i = 1 \\ j = 1 \\ i' = 1 \\ j' = 1}}^N \mathbf{A}_{i,j} \mathbf{B}_{i',j'} \mathbb{E}\left[ \bar{\Xi}_{i,j} \tilde{y}_{i'} \tilde{y}_{j'} \right] 
&= \sum_{i = 1}^n \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i}
    + \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j, j} \bar{\omega}_i \bar{\omega}_{j}  \label{eq:xi-y-two} \\
\sum_{\substack{i = 1 \\ j = 1 \\ i' = 1 \\ j' = 1}}^N  \mathbf{A}_{i,j} \mathbf{B}_{i',j'}\mathbb{E}\left[ \bar{\Xi}_{i,j} \bar{\Xi}_{i',j'} \right]
&= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{e}_i^2 \kappa_{2,i} 
        + \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_j  \label{eq:xi-two}
\end{align}
$$

{% endtab %}
{% tab y-tilde-exp proof %}
The first equation is as follows:

$$
\begin{aligned}
    \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \mathbf{A}_{i,j} \mathbf{B}_{i',j'} \mathbb{E}\left[ \tilde{y}_i \tilde{y}_{j} \tilde{y}_{i'} \tilde{y}_{j'} \right] 
    &= \sum_{i = 1}^N\mathbf{A}_{i,i} \mathbf{B}_{i,i} \mathbb{E}\left[ \tilde{y}_i^4 \right] 
        + \underbrace{\sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i',i'} \mathbb{E}\left[ \tilde{y}_i^2 \tilde{y}_{i'}^2 \right]}_{i = j; i' = j'; i \neq i'} 
        + \underbrace{\sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{i,j} \mathbb{E}\left[ \tilde{y}_i^2 \tilde{y}_{j}^2 \right]}_{i = i'; j = j'; i \neq j} 
        + \underbrace{\sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{j,i} \mathbb{E}\left[ \tilde{y}_i^2 \tilde{y}_{j}^2 \right]}_{i = j'; j = i'; i \neq j} \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} m_{4,i} 
        + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i',i'} \bar{\omega}_i^2 \bar{\omega}_j^2 \bar{\delta}_i^{-2} \bar{\delta}_j^{-2} m_{2,i} m_{2,j}  
        + \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i^2 \bar{\omega}_j^2 \bar{\delta}_i^{-2} \bar{\delta}_j^{-2} m_{2,i} m_{2,j}  
        + \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{j,i} \bar{\omega}_i^2 \bar{\omega}_j^2 \bar{\delta}_i^{-2} \bar{\delta}_j^{-2} m_{2,i} m_{2,j} \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} (\kappa_{4,i} + 3\kappa_{2,i}^2) 
        + \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_{j} 
        + 2 \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i \bar{\omega}_j 
        & \left(\mathbf{A}, \mathbf{B} \text{ sym.}\right)\\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
        + 3 \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} m_{2,i}^2 
        + \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_{j} 
        + 2 \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i \bar{\omega}_j \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
        + 3 \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^2 
        + \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_{j} 
        + 2 \sum_{i = 1}^N \sum_{j \neq i} \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i \bar{\omega}_j 
        & \left(\text{Eq. } \eqref{eq:omega-facts}\right) \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
        + 3 \underbrace{\sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^2}_{(a)} 
        + \left(\sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_j - \underbrace{\sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^2}_{(a)}\right) 
        + 2 \left( \sum_{i = 1}^N\sum_{j = 1}^N \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i \bar{\omega}_j - \underbrace{\sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i}\bar{\omega}_i^2}_{(a)}\right) \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
        + \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_j 
        + 2  \sum_{i = 1}^N\sum_{j = 1}^N \mathbf{A}_{i,j} \mathbf{B}_{i,j}\bar{\omega}_i \bar{\omega}_j 
        &\left( \text{cancel } (a)\right)\\
\end{aligned}
$$

Now the second equation:

$$
\begin{aligned}
  \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \mathbf{A}_{i,j} \mathbf{B}_{i',j'} \mathbb{E}\left[ \bar{\Xi}_{i,j} \tilde{y}_{i'} \tilde{y}_{j'} \right] 
    &= \sum_{i = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N  \mathbf{A}_{i,i} \mathbf{B}_{i',j'} \mathbb{E}\left[ \bar{\xi}_i \tilde{y}_{i'} \tilde{y}_{j'} \right]  
        & \left(\bar{\Xi} \text{ sym.}\right) \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \mathbb{E}\left[ \bar{\xi}_i \tilde{y}_i^2 \right] 
        + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i',i'} \mathbb{E}\left[ \bar{\xi}_i \tilde{y}_{i'}^2 \right] & \left(\mathbb{E}\left[ \bar{\xi}_i \tilde{y}_{i'} \tilde{y}_{j} \right] = 0 \text{ for } i' \neq j \right)\\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \mathbb{E}\left[ (\bar{\omega}_i + \bar{e}_i(y_i - \bar{\mu}_i)) \bar{\omega}_i^2 \bar{\delta}_i^{-2} (y_i - \bar{\mu}_i)^2 \right] 
        + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i',i'} \mathbb{E}\left[ (\bar{\omega}_i + \bar{e}_i(y_i - \bar{\mu}_i)) \bar{\omega}_{i'}^2 \bar{\delta}_{i'}^{-2} (y_{i'} - \bar{\mu}_{i'})^2 \right] \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \left( \bar{\omega}_i^3 \bar{\delta}_i^{-2} m_{2,i} + \bar{e}_i \bar{\omega}_i^2 \bar{\delta}_i^2 m_{3,i}\right) 
        + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i', i'} \bar{\omega}_i \bar{\omega}_{i'}^2 \bar{\delta}_{i'}^{-2} m_{2,i'}
        & \left( \text{Eq. } \eqref{eq:omega-facts} \right) \\
    &= \sum_{i = 1}^n \mathbf{A}_{i,i} \mathbf{B}_{i,i} \left(\bar{\omega}_i^2 + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right)
        + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i', i'} \bar{\omega}_i \bar{\omega}_{i'}
        & \left( \text{Eq. } \eqref{eq:omega-facts} \right) \\
    &= \sum_{i = 1}^n \mathbf{A}_{i,i} \mathbf{B}_{i,i} \left(\bar{\omega}_i^2 + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
        + \left( \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j, j} \bar{\omega}_i \bar{\omega}_{j} - \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{\omega}_i^2 \right)\\
    &= \sum_{i = 1}^n \mathbf{A}_{i,i} \mathbf{B}_{i,i}\bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} 
        + \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j, j} \bar{\omega}_i \bar{\omega}_{j}
\end{aligned}
$$

And then the third:

$$
\begin{aligned}
    \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \mathbf{A}_{i,j} \mathbf{B}_{i',j'}\mathbb{E}\left[ \bar{\Xi}_{i,j} \bar{\Xi}_{i',j'} \right]
    &= \sum_{i = 1}^N  \sum_{i' = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i',i'}\mathbb{E}\left[ \bar{\xi}_i \bar{\xi}_{i'} \right] 
        & \left(\bar{\Xi} \text{ is diag.}\right) \\
    &= \sum_{i = 1}^N  \mathbf{A}_{i,i} \mathbf{B}_{i,i}\mathbb{E}\left[ \bar{\xi}_i^2 \right] 
        + \sum_{i = 1}^N  \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i',i'}\mathbb{E}\left[ \bar{\xi}_i \bar{\xi}_{i'} \right] \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \mathbb{E}\left[ \bar{\omega}_i^2 
        + 2 \bar{\omega}_i \bar{e}_i(y_i -\bar{\mu}_i) + \bar{e}_i^2(y_i - \bar{\mu}_i)^2\right] + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i', i'}  \mathbb{E}\left[ \bar{\omega}_i \bar{\omega}_{i'} + \bar{\omega}_{i'} \bar{e}_i(y_i -\bar{\mu}_i)+ \bar{\omega}_{i} \bar{e}_{i'}(y_{i'} -\bar{\mu}_{i'}) + \bar{e}_i \bar{e}_{i'}(y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'})\right] \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \left(\bar{\omega}_i^2 + \bar{e}_i^2 m_{2,i} \right) 
        + \sum_{i = 1}^N \sum_{i' \neq i} \mathbf{A}_{i,i} \mathbf{B}_{i',i'} \bar{\omega}_i \bar{\omega}_{i'} \\
    &= \sum_{i = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{i,i} \bar{e}_i^2 \kappa_{2,i} 
        + \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{A}_{i,i} \mathbf{B}_{j,j} \bar{\omega}_i \bar{\omega}_j \\
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->


### Information Components
We define:

$$
\begin{aligned}
r_{i,i} &= \bar{\omega}_i^4 \bar{\delta}_i^{-4} \bar{\kappa}_{4, i} + 2 \bar{\omega}_i^2 + \bar{e}^2_i \bar{\kappa}_{2, i} - 2 \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \bar{\kappa}_{3, i} \\
r_{i, i'} &= 2 \bar{\omega}_i \bar{\omega}_{i'} \\
c_i &= \bar{\omega}_i^3 \bar{\delta}_i^{-3} \bar{\kappa}_{3,i} - \bar{\omega}_i \bar{\delta}_i^{-1} \bar{e}_i \bar{\kappa}_{2,i} \\
\mathbf{A}_j &= \mathbf{Z} \left. \left[ \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_j} \right] \right\rvert_{\tau^2 = \mathbf{0}_m} \mathbf{Z}^\top
\end{aligned}
$$

We also define $$\mathbf{a}_j = \text{diag}(\mathbf{A}_j)$$ (which is $n$-dimensional) and has elements denoted by $$a^j_{i,i}$$. 

The covariance of the score can be derived with a lot of algebra.

<!-- #region info-alpha-alpha -->
{% tabs info-alpha-alpha %}
{% tab info-alpha-alpha statement %}
$$
\begin{aligned}
\left. \mathcal{I}_{\alpha_j, \alpha_{j'}}(\theta) \right\rvert_{H_0}
&= \sum_{i = 1}^N \mathbf{x}_{i,j} \mathbf{x}_{i,j'} \bar{\omega}_i
\end{aligned}
$$
{% endtab %}
{% tab info-alpha-alpha proof %}
Recall the following definitions:

$$
\begin{aligned}
\rho_{i,j} &= \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} - \frac{1}{2} \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i \frac{\partial \bar{e}_i}{\partial \bar{\mu}_i} + \mathbf{x}_{i,j} \bar{\xi}_i \bar{\omega}_i \bar{\delta}_i^{-1} \left(\sum_{i' = 1}^N \mathbf{M}_{i,i'} \right) 
\\
\sigma_{i,j} &=  \frac{3}{2} \mathbf{M}_{i,i} \mathbf{x}_{i,j} \bar{\delta}_i^2 \bar{\omega}_i \frac{\partial^2 g(\bar{\mu}_i)}{\partial \bar{\mu}_i^2}
\end{aligned}
$$

Then: 

$$
\begin{aligned}
\mathcal{I}_{\alpha_j, \alpha_{j'}}(\theta)
&= \mathbb{E}\left[ \left(\sum_{i = 1}^N \left[ (y_i - \bar{\mu}_i) \rho_{i,j} +  \sigma_{i,j} \right]  \right) \left(\sum_{i = 1}^N \left[ (y_i - \bar{\mu}_i) \rho_{i,j'} +  \sigma_{i,j'} \right]  \right)\right] \\
&= \sum_{i = 1}^N \mathbb{E}\left[ ((y_i - \bar{\mu}_i) \rho_{i,j} +  \sigma_{i,j} )((y_i - \bar{\mu}_i) \rho_{i,j'} +  \sigma_{i,j'} )\right] + \sum_{i = 1}^N \sum_{i' \neq i} \mathbb{E}\left[ (y_i - \bar{\mu}_i) \rho_{i,j} +  \sigma_{i,j}  \right] \mathbb{E}\left[ (y_{i'} - \bar{\mu}_{i'}) \rho_{i',j'} +  \sigma_{i',j'}  \right] \\
&= \sum_{i = 1}^N \left[ \rho_{i,j} \rho_{i,j'} \mathbb{E}\left[ (y_i - \bar{\mu}_i)^2 \right] + \sigma_{i,j} \sigma_{i,j'} \right] + \sum_{i = 1}^N \sum_{i' \neq i} \sigma_{i,j} \sigma_{i',j'} \\
&= \sum_{i = 1}^N \rho_{i,j} \rho_{i,j'} m_{2,i}  + \sum_{i = 1}^N \sum_{i' = 1}^N \sigma_{i,j} \sigma_{i',j'} 
\end{aligned}
$$

Under $H_0$, $\tau^2 = \mathbf{0}_q$, so all of the terms involving $\mathbf{M}$ will be zero. Thus:

$$
\begin{aligned}
\left. \mathcal{I}_{\alpha_j, \alpha_{j'}}(\theta) \right\rvert_{H_0}
&= \sum_{i = 1}^N \mathbf{x}_{i,j} \mathbf{x}_{i,j'} \left( \bar{\omega}_i \bar{\delta}_i^{-1} \right)^2  m_{2,i}  \\
&= \sum_{i = 1}^N \mathbf{x}_{i,j} \mathbf{x}_{i,j'} \bar{\omega}_i \bar{\delta}_i^{-1} \bar{\delta}_i & \left(\text{Eq. } \eqref{eq:omega-facts} \right) \\
&= \sum_{i = 1}^N \mathbf{x}_{i,j} \mathbf{x}_{i,j'} \bar{\omega}_i
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

<!-- #region i-tau-tau -->
{% tabs i-tau-tau %}
{% tab i-tau-tau statement %}
$$
\begin{equation}
\label{eq:i-tau-tau}
\left. \mathcal{I}_{\tau^2_l, \tau^2_{l'}}(\theta) \right\rvert_{H_0}
= \frac{1}{4} \left[ \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} r_{i,i} + \sum_{i = 1}^N \sum_{j \neq i} (\mathbf{A}_{l})_{i,j} (\mathbf{A}_{l'})_{i,j} r_{i,j} \right] 
\end{equation}
$$
{% endtab %}
{% tab i-tau-tau proof %}
$$
\begin{aligned}
  \mathcal{I}_{\tau^2_l, \tau^2_{l'}}(\theta)
    &= \mathbb{E}\left[ U_{\tau^2_l}(\theta) U_{\tau^2_{l'}}(\theta) \right] \\
    &= \mathbb{E}\left[ \left( \frac{1}{2} \tilde{\mathbf{y}}^\top \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top \tilde{\mathbf{y}} - \frac{1}{2} \text{tr}\left[ \bar{\Xi} \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top \right] \right) \left( \frac{1}{2} \tilde{\mathbf{y}}^\top \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_{l'}} \mathbf{Z}^\top - \frac{1}{2} \text{tr}\left[ \bar{\Xi} \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_{l'}} \mathbf{Z}^\top \right] \right)\right] \\
    &= \frac{1}{4} \mathbb{E}\left[ \left(\sum_{i = 1}^N \sum_{j = 1}^N \tilde{y}_i \tilde{y}_j \dot{\mathbf{M}}_{i,j}^l - \sum_{i = 1}^N \sum_{j = 1}^N \bar{\Xi}_{i,j} \dot{\mathbf{M}}_{i,j}^l \right) \left(\sum_{i = 1}^N \sum_{j = 1}^N \tilde{y}_i \tilde{y}_j \dot{\mathbf{M}}_{i,j}^{l'} - \sum_{i = 1}^N \sum_{j = 1}^N \bar{\Xi}_{i,j} \dot{\mathbf{M}}_{i,j}^{l'} \right) \right] \\
    &= \frac{1}{4}\left[ \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[ \tilde{y}_i \tilde{y}_{j} \tilde{y}_{i'} \tilde{y}_{j'} \right]  
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[ \bar{\Xi}_{i,j} \tilde{y}_{i'} \tilde{y}_{j'} \right] 
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N  \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[ \bar{\Xi}_{i',j'} \tilde{y}_i \tilde{y}_j \right] 
        + \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[\bar{\Xi}_{i,j} \bar{\Xi}_{i',j'} \right]\right] \\
    &= \frac{1}{4} \left[ \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
        + \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{j,j} \bar{\omega}_i \bar{\omega}_j 
        + 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} (\mathbf{A}_{l'})_{i,j}\bar{\omega}_i \bar{\omega}_j 
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[ \bar{\Xi}_{i,j} \tilde{y}_{i'} \tilde{y}_{j'} \right] 
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N  \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[ \bar{\Xi}_{i',j'} \tilde{y}_i \tilde{y}_j \right] 
        + \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \bar{e}_i^2 \kappa_{2,i} 
        + \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{j,j} \bar{\omega}_i \bar{\omega}_j \right] 
        & \left(\text{Eq. } \eqref{eq:y-four}, \eqref{eq:xi-two}\right)\\
    &= \frac{1}{4} \left[ 
      \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \left(\bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + \bar{e}_i^2 \kappa_{2,i} \right)
        + 2\sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{j,j} \bar{\omega}_i \bar{\omega}_j 
        + 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} (\mathbf{A}_{l'})_{i,j}\bar{\omega}_i \bar{\omega}_j 
        - 2 \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \dot{\mathbf{M}}_{i',j'}^{l'} \mathbb{E}\left[ \bar{\Xi}_{i,j} \tilde{y}_{i'} \tilde{y}_{j'} \right]
    \right] \\
    &= \frac{1}{4} \left[ 
      \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \left(\bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + \bar{e}_i^2 \kappa_{2,i} \right)
        + 2\sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{j,j} \bar{\omega}_i \bar{\omega}_j 
        + 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} (\mathbf{A}_{l'})_{i,j}\bar{\omega}_i \bar{\omega}_j 
        - 2\sum_{i = 1}^n (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \left(\bar{\omega}_i^2 + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
        - 2\sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{j, j} \bar{\omega}_i \bar{\omega}_{j} 
        + 2\sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \bar{\omega}_i^2 
    \right] 
    & \left(\text{Eq. } \eqref{eq:xi-y-two}\right) \\
    &= \frac{1}{4} \left[ 
      \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \left(\bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + \bar{e}_i^2 \kappa_{2,i} -2\bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i}\right)
        + 2  \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,j} (\mathbf{A}_{l'})_{i,j}\bar{\omega}_i \bar{\omega}_j 
    \right] \\
    &= \frac{1}{4} \left[ 
      \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} \left(\bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + 2\bar{\omega}_i^2 + \bar{e}_i^2 \kappa_{2,i} -2\bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i}\right)
        + 2  \sum_{i = 1}^N \sum_{j \neq i} (\mathbf{A}_l)_{i,j} (\mathbf{A}_{l'})_{i,j}\bar{\omega}_i \bar{\omega}_j 
    \right] \\
    &= \frac{1}{4} \left[ \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} (\mathbf{A}_{l'})_{i,i} r_{i,i} + \sum_{i = 1}^N \sum_{j \neq i} (\mathbf{A}_{l})_{i,j} (\mathbf{A}_{l'})_{i,j} r_{i,j} \right] 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

<!-- #region i-tau-alpha -->
{% tabs i-tau-alpha %}
{% tab i-tau-alpha statement %}
$$
\begin{aligned}
\left. \mathcal{I}_{\alpha_j, \tau^2_l}(\theta)\right\rvert_{H_0}
&= \frac{1}{2} \mathbf{X}_{\cdot, j}^\top \mathbf{C} \mathbf{a}_l 
\end{aligned}
$$
{% endtab %}
{% tab i-tau-alpha proof %}
$$
\begin{aligned}
\mathcal{I}_{\alpha_j, \tau^2_l}(\theta)
&= \frac{1}{2} \mathbb{E}\left[ \left(\sum_{i = 1}^N \left[ (y_i - \bar{\mu}_i) \rho_{i,j} + \sigma_{i,j} \right] \right) \left( \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i,i'} \tilde{y}_i \tilde{y}_{i'} - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \bar{\xi}_i \right) \right] \\
&= \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N \sum_{i'' = 1}^N \mathbb{E}\left[ (\mathbf{A}_l)_{i',i''} \tilde{y}_{i'} \tilde{y}_{i''} \left[ (y_i - \bar{\mu}_i) \rho_{i,j} + \sigma_{i,j} \right] \right] + \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N \mathbb{E}\left[ (\mathbf{A}_l)_{i',i'} \bar{\xi}_{i'}  \left[ (y_i - \bar{\mu}_i) \rho_{i,j} + \sigma_{i,j} \right] \right] \\
&= \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N \sum_{i'' = 1}^N (\mathbf{A}_l)_{i',i''} \rho_{i,j}  \mathbb{E}\left[\tilde{y}_{i'} \tilde{y}_{i''}  (y_i - \bar{\mu}_i) \right] 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N \sum_{i'' = 1}^N (\mathbf{A}_l)_{i',i''} \sigma_{i,j}  \mathbb{E}\left[\tilde{y}_{i'} \tilde{y}_{i''}\right]
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'}  \rho_{i,j}  \mathbb{E}\left[ \bar{\xi}_{i'} (y_i - \bar{\mu}_i)\right]
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \mathbb{E}\left[\bar{\xi}_{i'} \right] \\
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \mathbb{E}\left[ \tilde{y}_i^2 (y_i - \bar{\mu}_i) \right] +  \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N \sum_{i'' = 1}^N (\mathbf{A}_l)_{i',i''} \sigma_{i,j}  \mathbb{E}\left[\tilde{y}_{i'} \tilde{y}_{i''}\right]
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'}  \rho_{i,j}  \mathbb{E}\left[ (\bar{\omega}_{i'} + \bar{e}_{i'}(y_{i'} - \bar{\mu}_{i'})) (y_i - \bar{\mu}_i)\right]
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \bar{\omega}_{i'} & \left(\mathbb{E}[\tilde{y}_{i'} \tilde{y}_{i''} (y_i - \bar{\mu}_i)] = 0 \text{ unless } i = i' = i''\right)  \\
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{\omega}_i^2 \left(\bar{\delta}_i\right)^{-2} m_{3,i} 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \mathbb{E}\left[ \tilde{y}_{i'}^2 \right] 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \rho_{i,j} \left(\bar{\omega}_{i'} \mathbb{E}\left[ y_i - \bar{\mu}_i \right] + \bar{e}_{i'} \mathbb{E}\left[ (y_{i'} - \bar{\mu}_{i'})(y_i - \bar{\mu}_i) \right]\right) 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i', i'} \sigma_{i,j} \bar{\omega}_{i'} & \left( \mathbb{E}[\tilde{y}_{i'} \tilde{y}_{i''}] 0 \text{ unless } i' = i'' \right) \\
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{\omega}_i^2 \left(\bar{\delta}_i\right)^{-2} m_{3,i} 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \bar{\omega}_{i'}^2 \left(\bar{\delta}_{i'}\right)^{-2} m_{2,i'} 
+ \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{e}_{i} \mathbb{E}\left[ (y_i - \bar{\mu}_i)^2 \right] 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i', i'} \sigma_{i,j} \bar{\omega}_{i'}& \left(\mathbb{E}\left[ (y_{i'} - \bar{\mu}_{i'})(y_i - \bar{\mu}_i) \right] = 0 \text{ unless } i = i'\right) \\
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{\omega}_i^2 \left(\bar{\delta}_i\right)^{-2} m_{3,i} 
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \bar{\omega}_{i'}
+ \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{e}_{i} m_{2,i}
+ \frac{1}{2} \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i', i'} \sigma_{i,j} \bar{\omega}_{i'}  & \left(\text{Eq. } \eqref{eq:omega-facts}\right) \\ 
&=  \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{\omega}_i^2 \left(\bar{\delta}_i\right)^{-2} m_{3,i} 
+ \sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \bar{\omega}_{i'}
+ \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{e}_{i} m_{2,i}
\end{aligned}
$$

Under $H_0$, we have $\rho_{i,j} = \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1}$ and $\sigma_{i,j} = 0$. Thus:

$$
\begin{aligned}
\left. \mathcal{I}_{\alpha_j, \tau^2_l}(\theta)\right\rvert_{H_0}
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{\omega}_i^2 \left(\bar{\delta}_i\right)^{-2} m_{3,i} 
+ \underbrace{\sum_{i = 1}^N \sum_{i' = 1}^N (\mathbf{A}_l)_{i',i'} \sigma_{i,j} \bar{\omega}_{i'}}_{=0}
+ \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \rho_{i,j} \bar{e}_{i} m_{2,i} \\
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{x}_{i,j} \bar{\omega}_i \bar{\delta}_i^{-1} \left[ \bar{\omega}_i^2 \left(\bar{\delta}_i\right)^{-2} m_{3,i} + \bar{e}_{i} m_{2,i} \right]\\ 
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{x}_{i,j} \left[ \bar{\omega}_i^3 \left(\bar{\delta}_i\right)^{-3} m_{3,i} + \bar{e}_{i} \bar{\omega}_i \bar{\delta}_i^{-1} m_{2,i} \right]
\end{aligned}
$$

And then:

$$
\begin{aligned}
\left. \mathcal{I}_{\alpha_j, \tau^2_l}(\theta)\right\rvert_{H_0}
&= \frac{1}{2} \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{x}_{i,j} c_i \\
&= \frac{1}{2} \mathbf{X}_{\cdot, j}^\top \mathbf{C} \mathbf{a}_l 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion i-tau-alpha -->

<div class="example">
    <strong>Example.</strong>
    <br>
    We can directly use the information calculations we did in the most general case for the components involving any shared fixed effects and the variance component. However, the cluster-specific intercept must be treated differently.
    <br>
    For arbitrary $j, j' \in [K]$, we have:
    $$
    \begin{aligned}
    \mathcal{I}_{\alpha_j^0, \alpha_{j'}^0}(\boldsymbol{\theta})
    &= \mathbb{E}\left[ \left(\sum_{i = 1}^{n_j} \left[ (y_i^j - \bar{\mu}_i^j) \rho_{i,0}^j + \sigma_{i,0}^j \right] \right) \left(\sum_{i = 1}^{n_{j'}}  \left[ (y_i^{j'} - \bar{\mu}_i^{j'}) \rho_{i,0}^{j'} + \sigma_{i,0}^{j'} \right] \right) \right] \\
    &= \sum_{i = 1}^{n_j} \sum_{i' = 1}^{n_{j'}} \mathbb{E} \left[(y_i^j - \bar{\mu}_i^j)(y_{i'}^{j'} - \bar{\mu}_{i'}^{j'}) \right] \rho_{i,0}^j \rho_{i', 0}^{j'} & \left(\sigma_{i,0}^k = 0 \text{ under } H_0\right)  \\
    &= \begin{cases}
    0 & \text{ if } j \neq j' \\
    \sum_{i = 1}^{n_j} (\rho_{i,0}^j)^2 m_{2,i}^j & \text{ if } j = j'
    \end{cases} 
    \end{aligned}
    $$
    For the cross terms involving the intercept and slope:
    $$
    \begin{aligned}
    \mathcal{I}_{\alpha_j^0, \alpha_1}(\boldsymbol{\theta})
    &= \mathbb{E}\left[ \left(\sum_{i = 1}^{n_j} \left[ (y_i^j - \bar{\mu}_i^j) \rho_{i,0}^j + \sigma_{i,0}^j \right] \right) \left(\sum_{i = 1}^{n_k} \sum_{k = 1}^K  \left[ (y_i^{k} - \bar{\mu}_i^{k}) \rho_{i,1}^{k} + \sigma_{i,1}^{k} \right] \right) \right] \\
    &= \sum_{i = 1}^{n_j} \sum_{i' = 1}^{n_k} \sum_{k = 1}^K \mathbb{E} \left[(y_i^j - \bar{\mu}_i^j)(y_{i'}^{k} - \bar{\mu}_{i'}^{k}) \right] \rho_{i,1}^k \rho_{i', 1}^{k} & \left(\sigma_{i,0}^k = 0 \text{ under } H_0\right)  \\
    &= \sum_{i = 1}^{n_j} \sum_{i' = 1}^{n_j} \mathbb{E}\left[ (y_i^j - \bar{\mu}_i^j)(y_{i'}^j - \bar{\mu}_{i'}^j) \right] \rho_{i,1}^j \rho_{i',1}^j & \left( 0 \text{ if } j \neq k\right) \\
    &= \sum_{i = 1}^{n_j} \left(\rho_{i,1}^j \right)^2 m_{2,i} & \left( 0 \text{ if } i \neq i' \right)
    \end{aligned}
    $$
    Under $H_0$ and using Eq. \eqref{eq:omega-facts}, we see:
    $$
    \begin{aligned}
    \mathcal{I}_{\alpha_j^0, \alpha_{j}^0}(\boldsymbol{\theta})
    &= \sum_{i = 1}^{n_j} \left[ \bar{\omega}_i^j \left( \bar{\delta}_i^j\right)^{-1} \right]^2 m_{2,i} = \sum_{i = 1}^{n_j} \bar{\omega}_i^j  \\
    \mathcal{I}_{\alpha_j^0, \alpha_1}(\boldsymbol{\theta}) 
    &= \sum_{i = 1}^{n_j} \mathbf{x}_{i,j}  \left[ \bar{\omega}_i^j \left(\bar{\delta}_i^j\right)^{-1} \right]^2 m_{2,i} 
    = \sum_{i = 1}^{n_j} \mathbf{x}_{i,j} \bar{\omega}_i^j 
    \end{aligned}
    $$
</div>



































































<!-- 
---

$$
\begin{aligned}
\mathcal{I}_{\tau^2_l, \phi}(\theta) 
    &= \mathbb{E}\left[ U_{\tau^2_l}(\theta) U_\phi (\theta) \right] \\
    &= \mathbb{E}\left[ \left( \frac{1}{2} \tilde{\mathbf{y}}^\top \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top - \frac{1}{2} \text{tr}\left[ \bar{\Xi} \mathbf{Z} \frac{\partial \mathbf{D}(\tau^2)}{\partial \tau^2_l} \mathbf{Z}^\top \right] \right) \left( \frac{1}{\phi} \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) - \frac{1}{\phi} \tilde{\mathbf{y}}^\top \mathbf{M} \tilde{\mathbf{y}} + \frac{1}{2} \text{tr}\left[ \frac{\partial \bar{\Xi}}{\partial \phi} \mathbf{M} \right] \right)\right] \\
    &= \frac{1}{2} \mathbb{E}\left[ \left(\sum_{i = 1}^N \sum_{j = 1}^N \tilde{y}_i \tilde{y}_j \dot{\mathbf{M}}_{i,j}^l - \sum_{i = 1}^N \sum_{j = 1}^N \bar{\Xi}_{i,j} \dot{\mathbf{M}}_{i,j}^l \right) \left( \frac{1}{\phi} \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) - \frac{1}{\phi} \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{M}_{i,j} \ \tilde{y}_i \tilde{y}_j + \frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \frac{\partial \bar{\Xi}_{i,j}}{\partial \phi} \mathbf{M}_{i,j} \right)\right] \\
    &=\frac{1}{2} \mathbb{E}\left[ \left(\sum_{i = 1}^N \sum_{j = 1}^N \dot{\mathbf{M}}_{i,j}^l (\tilde{y}_i \tilde{y}_j -  \bar{\Xi}_{i,j})  \right) \left( \frac{1}{\phi} \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) - \frac{1}{\phi} \sum_{i = 1}^N \sum_{j = 1}^N \mathbf{M}_{i,j}  \tilde{y}_i \tilde{y}_j + \frac{1}{2} \sum_{i = 1}^N \mathbf{M}_{i,i} \frac{\partial \bar{\xi}_i}{\partial \phi}  \right) \right] & \left(\bar{\Xi} \text{ is diagonal}\right) \\
    &= \frac{1}{2} \left[ \sum_{i = 1}^N \sum_{j = 1}^N \frac{\dot{\mathbf{M}}_{i,j}^l}{\phi} \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\Xi}_{i,j}) \right] 
        - \frac{1}{\phi} \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',j'} \mathbb{E}\left[ (\tilde{y}_i \tilde{y}_j - \bar{\Xi}_{i,j}) \tilde{y}_{i'} \tilde{y}_{j'} \right] 
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',i'} \mathbb{E}\left[ \bar{\xi}_{i'} (\tilde{y}_i \tilde{y}_j - \bar{\Xi}_{i,j})\right] \right] \\
    &= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i', j'} \mathbb{E}\left[ \tilde{y}_i \tilde{y}_j \tilde{y}_{i'} \tilde{y}_{j'} \right] 
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i', j'} \mathbb{E}\left[ \bar{\Xi}_{i,j} \tilde{y}_{i'} \tilde{y}_{j'} \right]
        - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',i'} \mathbb{E}\left[ \bar{\xi}_{i'} (\tilde{y}_i \tilde{y}_j - \bar{\Xi}_{i,j})\right] \right]  \\
    &=  \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} 
      - \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j 
      - \sum_{i = 1}^n (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left(\bar{\omega}_i^2 + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j, j} \bar{\omega}_i \bar{\omega}_{j} 
      + \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \bar{\omega}_i^2 
      - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',i'} \mathbb{E}\left[ \bar{\xi}_{i'} (\tilde{y}_i \tilde{y}_j - \bar{\Xi}_{i,j})\right] \right] & \left(\text{Eq. } \eqref{eq:y-four}, \eqref{eq:xi-y-two}\right) \\
    &= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 2 \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j 
      - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',i'} \mathbb{E}\left[ \bar{\xi}_{i'} \tilde{y}_i \tilde{y}_j \right] +  \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',i'} \mathbb{E}\left[\bar{\xi}_{i'} \bar{\Xi}_{i,j}\right] \right] \\
  &= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 2 \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j 
      - \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',j'} \mathbb{E}\left[ \bar{\Xi}_{i',j'} \tilde{y}_i \tilde{y}_j \right] +  \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',j'} \mathbb{E}\left[\bar{\Xi}_{i', j'} \bar{\Xi}_{i,j}\right] \right]
       & \left(\bar{\Xi} \text{ diag.}\right) \\
\end{aligned}
$$


$$
\begin{aligned}
&= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 2 \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j 
      - \sum_{i = 1}^n (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left(\bar{\omega}_i^2 + \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j, j} \bar{\omega}_i \bar{\omega}_{j} 
      + \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \bar{\omega}_i^2  +  \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',j'} \mathbb{E}\left[\bar{\Xi}_{i', j'} \bar{\Xi}_{i,j}\right] \right]
       & \left(\text{Eq. } \eqref{eq:xi-y-two}\right) \\
&= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} - \bar{\omega}_i^2 + 2 \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 4 \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j  
      +  \sum_{i = 1}^N \sum_{j = 1}^N \sum_{i' = 1}^N \sum_{j' = 1}^N \dot{\mathbf{M}}_{i,j}^l \mathbf{M}_{i',j'} \mathbb{E}\left[\bar{\Xi}_{i', j'} \bar{\Xi}_{i,j}\right] \right] \\
  &= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} - \bar{\omega}_i^2 + 2 \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 4 \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j  
      + \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \bar{e}_i^2 \kappa_{2,i} 
      + \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j  \right] & \left(\text{Eq. } \eqref{eq:xi-two}\right)\\
&= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} - \bar{\omega}_i^2 -\bar{e}_i^2 \kappa_{2,i} + 2 \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 3 \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \bar{\omega}_i^2 
      - 3 \sum_{i = 1}^N \sum_{j \neq i} (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N\sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j  \right] \\
&= \frac{1}{2\phi} \left[ \sum_{i = 1}^N  \dot{\mathbf{M}}_{i,i}^l \mathbb{E}\left[ \ell_q(\mathbf{y}; \theta \rvert \beta = \mathbf{0}_{m \times k}) (\tilde{y}_i \tilde{y}_j - \bar{\xi}_i) \right] 
      - \sum_{i = 1}^N (\mathbf{A}_l)_{i,i} \mathbf{M}_{i,i} \left( \bar{\omega}_i^4 \bar{\delta}_i^{-4} \kappa_{4,i} + 2 \bar{\omega}_i^2 - \bar{e}_i^2 \kappa_{2,i} + 2 \bar{\omega}_i^2 \bar{\delta}_i^{-2} \bar{e}_i \kappa_{3,i} \right) 
      - 3 \sum_{i = 1}^N \sum_{j \neq i} (\mathbf{A}_l)_{i,i} \mathbf{M}_{j,j} \bar{\omega}_i \bar{\omega}_j 
      - 2  \sum_{i = 1}^N \sum_{j = 1}^N (\mathbf{A}_l)_{i,j} \mathbf{M}_{i,j}\bar{\omega}_i \bar{\omega}_j  \right] \\
\end{aligned}
$$

 -->


---





<!-- ---

## Individual Variance Component Tests
We may also be interested in testing whether a single coordinate of the variance component is zero while not impposing that the other are. This is testing the hypotheses:

$$
H_0: \tau^2_j = 0 \hspace{15mm} \text{vs.} \hspace{15mm} H_1: \tau^2_j > 0
$$

where we notice that the alternative is restricted to the positive reals.






We'll use a subscript $-j$ to denote a vector with the $j$-th component removed; i.e. $$\mathbf{v}_{-j} = (\mathbf{v}_1, \dots, \mathbf{v}_{j - 1}, \mathbf{v}_{j + 1}, \dots, \mathbf{v}_m)^\top$$. We'll let $f(\beta_j)$and $f(\beta_{-j})$ denote the log density functions of $\beta_j$ and $\beta_{-j}$, respectively. Define:

$$
\ell_q(\mathbf{y}; \beta_j) = \log \int \exp \left( \sum_{i = 1}^N \ell_q(y_i; \theta \rvert \beta) + f(\beta_{-j}) \right) d\beta_{-j}
$$

and rewrite the marginal log quasi-likelihood as:

$$
\begin{aligned}
\ell_q(\mathbf{y}; \theta) 
&= \log \int \exp\left( \ell_q(\mathbf{y}; \theta \rvert \beta_j) + f(\beta_j) \right) d \beta_j \\
&= \log \int \exp\left( \log \int \exp \left( \sum_{i = 1}^N \ell_q(y_i; \theta \rvert \beta) + f(\beta_{-j}) \right) d\beta_{-j} + f(\beta_j) \right) d \beta_j \\
&= \log \int 
\end{aligned}
$$
 -->
