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
    - name: Quasi-Likelihood
    - name: Score
      subsections:
          - name: Intermediate Quantities
          - name: Derivative Terms
          - name: Computing The Score
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
We will assume to have $N$ total observations. For each, we have some response, $y_i$, a $p \times 1$ fixed effects covariate vector, $\mathbf{x}_i$, and a $q \times 1$ random effects covariate vector, $\mathbf{z}_i$. 

We'll assume to have a clustered design, so each observation will belong to some cluster $k \in [K]$. The observations in cluster $k$ will be denoted with a superscript $k$ (i.e. $y_j^k$), and we'll let $n_k$ denote the number of observations in cluster $k$ (such that $\sum_{k = 1}^K n_k = N$). 

We assume that, conditional on the $q \times 1$ vector of random effects, $\beta^k$, the responses are independent and have means and variances:

$$
\begin{aligned}
\mathbb{E}[y^k_j \rvert \beta_k] &= \mu^k_j \\
\text{Var}(y^k_j \rvert \beta_k) &= V(\mu^k_j) = \phi v(\mu^k_j)
\end{aligned}
$$

for some scale parameter, $\phi$. We further assume a generalized linear model is appropriate:

$$
\begin{aligned}
g(\mu^k_j) &= \eta^k_j = \mathbf{x}_i^\top \alpha + \mathbf{z}_i^\top \beta_k \\
\beta_k &\overset{iid}{\sim} \mathcal{N}(\mathbf{0}_q, D(\tau^2))
\end{aligned}
$$

with monotonic and differentiable $g(\cdot)$, and $p \times 1$ vector of fixed effects. We assume that $D(\tau^2)$ is a $q \times q$ diagonal matrix parametrized by the $m$-dimensional vector-valued variance component, $\tau^2 = (\tau^2_1, \dots, \tau^2_m)^\top$. 

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
\mathbf{y} &= (\mathbf{y}^1, \dots, \mathbf{y}^K)^\top \\
\mu &= (\mu^1, \dots \mu^K)^\top \\
\eta &= (\eta^1, \dots, \eta^K)^\top
\end{aligned}
$$

Furthermore, we will let $\mathbf{V}(\mu)$ denote the diagonal matrix with the $V(\mu_i) = \phi v(\mu_i)$ terms along its main diagonal, and $\mathbf{X}$ will be the $N \times p$ matrix with rows equal to the $\mathbf{x}_i^\top$ vectors.

<div class="example">
  <strong>Example.</strong>
  <br>
  In this example, we'll assume we have negative binomial data with a scalar variance component (i.e. $m = 1$), a single random effect, and a fixed, cluster-specific intercept:
  $$
  \begin{aligned}
  g(\mu_j^k) &= \alpha_k + \beta_k z_i \\
  \beta &\sim \mathcal{N}(\mathbf{0}_K, \tau^2 \mathbb{I}_{K \times K})
  \end{aligned}
  $$
  The link and variance functions are:
   $$
   \begin{aligned}
   g(\mu^k_i) &= \log(\mu^k_i) \\
   v(\mu_i^k) &= \mu_i^k + \frac{1}{\gamma} (\mu_i^k)^2
   \end{aligned}
   $$
   where we assume $\gamma$ is a known dispersion parameter.
</div>

### Some Non-Standard Notation
Let $$\tilde{\mathbf{z}}_i$$ be the $(q \times k)$-dimensional vector of zeros where the coordinates corresponding to observation $i$ are replaced with $\mathbf{z}_i$. We'll use $$\tilde{\mathbf{Z}}$$ to denote the $$N \times (q \times K)$$ matrix whose $i$-th row is $$\tilde{\mathbf{z}}_i^\top$$. Similarly, we'll use $$\tilde{D(\tau^2)} = D(\tau^2) \otimes \mathbb{I}_{k \times k}$$ to denote the $(q \times K) \times (q \times K)$ block covariance matrix for <i>all</i> of the random effects. 

We'll also use a bar (e.g. $$\bar{\mu}_i^k$$) to denote a term evaluated under the null hypothesis that $$\beta = \mathbf{0}_{q \times k}$$, and a hat (e.g. $$\hat{\mu}_i^k$$) to denote a term evaluated at the parameter maximum likelihood estimates under the null hypothesis as well. 


---

## Quasi-Likelihood
The model has parameters $\alpha$, $\phi$, $\tau^2$, and random effects $\beta = (\beta_1^\top, \dots, \beta_K^\top)^\top$. Let $\theta = (\alpha, \phi, \tau^2)$. The conditional log quasi-likelihood is:

$$
\ell_q(y^k_i; \theta \rvert \beta_k) = \int_{y^k_i}^{\mu_i^k} \frac{y^k_i - u}{\phi v(u)} du
$$ 

Let $f(\beta; \tau^2)$ and $\mathcal{F}(\beta; \tau^2)$ denote the log density and distribution functions of $\beta$. Under the assumption that $D(\tau^2)$ is diagonal (i.e. the random effects are independent), the integrated (marginal) quasi-likelihood is:

$$
\begin{aligned}
\mathcal{L}_q(\mathbf{y}; \theta) 
&= \int \exp\left( \sum_{k = 1}^K \sum_{i = 1}^{n_k} \ell_q(y^k_i; \theta \rvert \beta_k) \right) d \mathcal{F}(\beta; \tau^2)  \\
&= \int \exp\left(  \sum_{k = 1}^K \sum_{i = 1}^{n_k} \ell_q(y^k_i; \theta \rvert \beta_k)  \right) \exp(f(\beta; \tau^2)) d \beta 
\end{aligned}
$$

And the marginal log quasi-likelihood is the logarithm of this function:

$$
\mathcal{L}_q(\mathbf{y}; \theta) 
= \log \left[ \int \exp\left(  \sum_{k = 1}^K \sum_{i = 1}^{n_k} \ell_q(y^k_i; \theta \rvert \beta_k)  \right) \exp(f(\beta; \tau^2)) d \beta  \right]
$$

This integral is unwieldy, so we use a second-order Taylor approximation of the conditional quasi-likelihood about the random effects mean to simplify our calculations. If observation $i$ is a member of cluster $k$, we will let $\tilde{\mathbf{z}}_i$ denote the $(q \times K)$-dimensional whose only non-zero entries are the $qk$ to $q(k + 1)$ (those corresponding to the cluster of observation $i$). The integrated (conditional) quasi-likelihood can be written via a Taylor expansion:

$$
\begin{aligned}
\mathcal{L}_q\left( \mathbf{y}; \theta \rvert \beta \right)
&= \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \\
&= \left. \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \right\rvert_{\beta = \mathbf{0}_{q \times k}} + \left. \frac{\partial}{\partial \beta} \left[ \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \right]\right\rvert_{\beta = \mathbf{0}_{q \times k}} \beta + \frac{1}{2}\beta^\top  \left[ \left. \frac{\partial}{\partial \beta \partial \beta^\top} \left[ \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right)  \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} \right]  \beta + \dots \\
&=  \left. \exp\left( \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right) \right\rvert_{\beta = \mathbf{0}_{q \times k}} \left( 1 + \left. \frac{\partial}{\partial \beta} \left[ \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right]\right\rvert_{\beta = \mathbf{0}_{q \times k}} \beta + \frac{1}{2}\beta^\top \left[  \left.\frac{\partial}{\partial \beta \partial \beta^\top} \left[ \sum_{i = 1}^{N} \ell_q(y_i; \theta \rvert \beta) \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} \right]  \beta + \epsilon \right) 
\end{aligned}
$$

where $\epsilon$ is the residual that holds the higher order terms. Taking the expectation and applying the chain rule yields an approximate marginal quasi log-likelihood:

$$
\mathcal{L}_q(\mathbf{y}; \theta) 
\approx \left. \exp\left(\sum_{i = 1}^N \ell_q(y_i; \theta \rvert \beta) \right)\right\rvert_{\beta = \mathbf{0}_{q \times k}} \left(  1 + \frac{1}{2} \text{tr} \left[ \left( \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \tilde{\mathbf{z}}_i \right)  \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta )}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \tilde{\mathbf{z}}_i^\top \right) + \sum_{i = 1}^N \left. \frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2}\right\rvert_{\beta = \mathbf{0}_{q \times k}} \tilde{\mathbf{z}}_i \tilde{\mathbf{z}}_i^\top \right) D(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert) \right)
$$

Taking the logarithm and using the fact that $\log(x + 1) \approx x$ for small $x$, we can write the marginal log quasi-likelihood approximation as:

$$
\begin{aligned}
\ell_q(\mathbf{y}; \theta) 
&\approx 
\sum_{i = 1}^N \left. \ell_q(y_i; \theta \rvert \beta)\right\rvert_{\beta = \mathbf{0}_{q \times k}} + \frac{1}{2} \text{tr} \left[ \left( \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \tilde{\mathbf{z}}_i \right)  \left(\sum_{i = 1}^N \left. \frac{\partial \ell_q(y_i; \theta \rvert \beta )}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \tilde{\mathbf{z}}_i^\top \right) + \sum_{i = 1}^N \left. \frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \tilde{\mathbf{z}}_i \tilde{\mathbf{z}}_i^\top \right) \tilde{D}(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert) \\
&= \sum_{i = 1}^N \left. \ell_q(y_i; \theta \rvert \beta)\right\rvert_{\beta = \mathbf{0}_{q \times k}}+ \frac{1}{2} \text{tr} \left[ \tilde{\mathbf{Z}}^\top \left( \left. \left[ \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta} \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta^\top} + \frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta )}{\partial \eta \partial \eta^\top} \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} \right) \tilde{\mathbf{Z}} \tilde{D}(\tau^2) \right] + o(\rvert \rvert \tau^2 \rvert \rvert)
\end{aligned}
$$

where $$\frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta}$$ is the $N$-dimensional vector whose $i$-th element is equal to $$\frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i}$$, and $$\frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta )}{\partial \eta \partial \eta^\top}$$ is the $N \times N$ dimensional diagonal matrix whose $i$-th diagonal element is equal to $$\frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2}$$. 

We will drop the $o(\rvert \rvert \tau^2 \rvert \rvert)$ in what follows as we have assumed moment conditions hold that make this term negligible. 

---

## Score
### Intermediate Quantities
We will define some helpful quantities that will make the notation easier as we continue. Recall that a single subscript (with no superscript) directly indexes the <i>entire</i> vector (over all clusters). We should also note that we further assume that $g(\cdot)$ is <i>strictly monotone</i> and <i>continuous</i> so that:

$$
\begin{equation}
\label{eq:deriv-assumption}
\frac{\partial g^{-1}(\eta_i)}{\partial \eta_i} = \frac{\partial \mu_i}{\partial \eta_i} = \left(\frac{\partial \eta_i}{\partial \mu_i}\right)^{-1}
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
e_i &= \frac{\frac{\partial V(\mu_i)}{\partial \mu_i} \frac{\partial g(\mu_i)}{\partial \mu_i} + V(\mu_i) \frac{\partial^2 g(\mu_i)}{\partial \mu_i^2}}{(V(\mu_i)^2)\left(\frac{\partial g(\mu_i)}{\partial \mu_i}\right)^3} = \frac{\phi\left(\frac{\partial v(\mu_i)}{\partial\mu_i}\right)\left(\frac{\partial g(\mu_i)}{\partial\mu_i}\right) + \phi v(\mu_i)\left(\frac{\partial^2 g(\mu_i)}{\partial \mu_i^2}\right)}{\phi^2 v^2(\mu_i) \left(\frac{\partial g(\mu_i)}{\partial\mu_i}\right)^3} \\
\xi_i &= \omega_i + e_i(y_i - \mu_i)
\end{aligned}
\end{equation}
$$

<aside><p>In the original publication, $\omega_i$ is denoted by $w_i$, and $\xi_i$ is denoted by $w_{oi}$.</p></aside>

We also define the $N \times N$ diagonal matrices $\Delta$, $\Omega$, and $\Xi$, which have the $\delta_i$, $\omega_i$, and $\xi_i$ terms along their main diagonals (respectively).

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

### Derivative Terms
Notice that the expression above all involve the partial derivatives of $$\ell_q(\mathbf{y}; \theta \rvert \beta)$$ with respect to the components of the linear predictor, $\eta$. With an application of the chain rule, we can derive the first order partial derivatives. Recall that a bar above a variable indicates its evaluation at the true parameter values under $H_0$ (i.e. $$\beta = \mathbf{0}_{q \times k}$$). To streamline our analysis later, let's do some initial derivations. We have:

<!-- #region deriv-1-mu -->
{% tabs deriv-1-mu %}
{% tab deriv-1-mu equation %}
$$
\begin{equation}
\label{eq:deriv-1-mu}
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \mu_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
&=  \frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)} \\
\implies
\left. \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \mu} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
&= \mathbf{V}^{-1}(\bar{\mu})(\mathbf{y}_i - \mu_i)
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-1-mu proof %}
$$
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \mu_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
&= \left. \frac{\partial}{\partial \mu_i} \left[ \int_{y_i}^{\mu_i} \frac{y_i - u}{\phi v(u)} du \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} \\
&= \left. \left[ \frac{y_i - \mu_i}{\phi v(\mu_i)} \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} & \left(\text{FTC}\right) \\
&= \frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

<!-- #region deriv-1-eta -->
{% tabs deriv-1-eta %}
{% tab deriv-1-eta equation %}
$$
\begin{equation}
\label{eq:deriv-1-eta}
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
&= \bar{\omega}_i \bar{\delta}_i^{-1} \left(y_i - \bar{\mu}_i\right) \\
\implies
\left. \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta} \right\rvert_{\beta = \mathbf{0}_{q \times k}}
&= \bar{\Omega}\bar{\Delta}^{-1} (\mathbf{y} - \bar{\mu})
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab deriv-1-eta proof %}
$$
\begin{aligned}
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
&= \left. \left[ \frac{\partial \mu_i}{\partial \eta_i} \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \mu_i} \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} & \left(\text{chain rule} \right)\\
&= \left. \left[ \frac{\partial \mu_i}{\partial \eta_i} \frac{\partial}{\partial \mu_i} \left[ \int_{y_i}^{\mu_i} \frac{y_i - u}{\phi v(u)} du \right] \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}} \\
&= \left. \left[ \left(\frac{\partial \mu_i}{\partial \eta_i}\right) \left(\frac{y_i - \mu_i}{\phi v(\mu_i)}\right) \right] \right\rvert_{\beta = \mathbf{0}_{q \times k}}  & \left( \text{FTC} \right) \\
&= \left. \frac{\partial \mu_i}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \left(\frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)}\right)
\end{aligned}
$$

Note that, by Eq. \eqref{eq:deriv-assumption}:

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
\left. \frac{\partial \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
&= \left. \frac{\partial \mu_i}{\partial \eta_i} \right\rvert_{\beta = \mathbf{0}_{q \times k}} \left(\frac{y_i - \bar{\mu}_i}{\phi v(\bar{\mu}_i)}\right) \\
&= \left. \frac{\left(\frac{\partial \mu_i}{\partial \eta_i}\right)^2}{\phi v(\mu_i) \frac{\partial \mu_i}{\partial \eta_i}}\right\rvert_{\beta = \mathbf{0}_{q \times k}} \left(y_i - \bar{\mu}_i\right) \\
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
\left. \frac{\partial^2 \ell_q(y_i; \theta \rvert \beta)}{\partial \eta_i^2} \right\rvert_{\beta = \mathbf{0}_{q \times k}}
&= -\bar{\xi}_i \\
\implies \left. \frac{\partial^2 \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \eta \partial \eta^\top} \right\rvert_{\beta = \mathbf{0}_{q \times k}} 
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

We can then rewrite the approximate marginal log quasi-likelihood as:

$$
\begin{equation}
\label{eq:rewrite-mlql}
\begin{aligned}
\ell_q(\mathbf{y}; \theta)
&\approx \left. \ell_q(\mathbf{y}; \theta \rvert \beta)\right\rvert_{\beta = \mathbf{0}_{q \times k}}+ \frac{1}{2} \text{tr} \left[ \tilde{\mathbf{Z}}^\top \left( \bar{\Omega}\bar{\Delta}^{-1}(\mathbf{y} - \bar{\mu})(\mathbf{y} - \bar{\mu})^\top \bar{\Delta}^{-1} \bar{\Omega} - \bar{\Xi} \right) \tilde{\mathbf{Z}} \tilde{D}(\tau^2) \right] 
\end{aligned}
\end{equation}
$$

<div class="example">
  <strong>Example.</strong>
  <br>
  Recall that in our negative binomial setting, the link function is $\log(\cdot)$, and the variance function is $v(\mu_i^k) = \mu_i^k + \frac{1}{\gamma} (\mu_i^k)^2$. Thus:
  $$
  \begin{aligned}
  \frac{\partial \mu_i^k}{\partial \eta_i^k} &= \frac{\partial}{\partial \eta_i^k} \left[ \exp(\eta_i^k) \right] = \exp(\eta_i^k) = \mu_i^k \\
  \frac{\partial v(\mu_i^k)}{\partial \mu_i^k} &= \frac{\partial}{\partial \mu_i^k}\left[ \mu_i^k + \frac{1}{\gamma}(\mu_i^k)^2 \right] = 1 + \frac{2}{\gamma} \mu_i^k
  \end{aligned}
  $$
  Plugging these values into Eqs. \eqref{eq:deriv-1-eta} and \eqref{eq:deriv-2-eta} gives us:
  $$
  \begin{aligned}
  \frac{\partial \ell_q(y_i^k; \theta \rvert \beta_i = 0)}{\partial \eta_i^k} 
  &= \left(\frac{y_i^k - \bar{\mu}_i^k}{\phi\left(\bar{\mu}_i^k + \frac{1}{\gamma}(\bar{\mu}_i^k)^2\right)}\right)\left(\bar{\mu}_i^k \right) \\
  &= \frac{\bar{\mu}_i^k (y_i^k - \bar{\mu}_i^k)}{\phi\left(\bar{\mu}_i^k + \frac{1}{\gamma}(\bar{\mu}_i^k)^2\right)} \\
  \frac{\partial^2 \ell_q(y_i^k; \theta \rvert \beta_i)}{\partial (\eta^i_k)^2} 
  &= \left(- \frac{\bar{\mu}_i + \frac{1}{\gamma}(\bar{\mu}_i^k)^2 + \left(1 + \frac{2}{\gamma} \bar{\mu}_i^k \right)(y_i^k - \bar{\mu}_i^k)}{\phi \left(\bar{\mu}_i^k + \frac{1}{\gamma}(\bar{\mu}_i^k)^2\right)^2}\right)\left(\bar{\mu}_i^k \right)^2 + \left(\frac{y_i^k - \bar{\mu}_i^k}{\phi \left(\bar{\mu}_i^k + \frac{1}{\gamma}(\bar{\mu}_i^k)^2 \right)}\right) \left(\bar{\mu}_i^k \right) \\
  &= - \frac{\bar{\mu}^k_i + \frac{1}{\gamma}(\bar{\mu}_i^k)^2 + \left(1 + \frac{2}{\gamma} \bar{\mu}_i^k \right)(y_i^k - \bar{\mu}_i^k)}{\phi \left(1 + \frac{1}{\gamma}\bar{\mu}_i^k\right)^2} + \frac{y_i^k - \bar{\mu}_i^k}{\phi \left(1 + \frac{1}{\gamma}\bar{\mu}_i^k \right)}
  \end{aligned}
  $$
</div>

### Computing The Score
This brings us to the actual score, which is the gradient of the log-likelihood (or of the approximate marginal log quasi-likelihood in this case). Unfortunately, the derivative with respect to the $\alpha_j$ terms is very difficult to compute, and we don't actually need it explicitly for the test statistic. We'll just focus on the derivative with respect to the variance component, $\tau^2$, here. 

Notice that $\tau^2$ only appears in the covariance matrix for the random effects within the approximate marginal log quasi-likelihood. Thus:

$$
\begin{equation}
\label{eq:score-tau}
\begin{aligned}
U_{\tau^2_j}(\hat{\theta}) &= \left. \frac{\partial \ell_q(\mathbf{y}; \theta \rvert \beta)}{\partial \tau^2_j} \right\rvert_{\theta = \hat{\theta}} \\
&= \frac{1}{2} \text{tr} \left[\tilde{\mathbf{Z}}^\top \left(\bar{\Omega} \bar{\Delta}^{-1}(\mathbf{y} - \bar{\mu})(\mathbf{y} - \bar{\mu})^\top \bar{\Delta}^{-1} \bar{\Omega} - \bar{\Xi}\right) \tilde{\mathbf{Z}} \left. \left[ \frac{\partial \tilde{D}(\tau^2)}{\partial \tau^2_j} \right] \right\rvert_{\theta = \hat{\theta}} \right] 
\end{aligned}
\end{equation}
$$



---

## Information

### Intermediate Quantities

Let $\kappa_{r, i}$ denote the $r$-th cumulant of $y_i$. If we assume that the following relationship holds for $r = 2$ and $r = 3$:

$$
\begin{aligned}
\kappa_{2, i} &= \phi v(\mu_i) \\
\kappa_{(r+1), i} &= \kappa_{2, i} \frac{\partial \kappa_{r, i}}{\partial \mu_i}
\end{aligned}
$$

then we have: 

$$
\begin{aligned}
\kappa_{3, i} &= \phi^2 v(\mu_i) \frac{\partial v(\mu_i)}{\partial\mu_i} \\
\kappa_{4, i} &= \phi^3 v(\mu_i) \left( v(\mu_i) \frac{\partial^2 v(\mu_i)}{\partial\mu_i^2} + \left(\frac{\partial v(\mu_i)}{\partial \mu_i} \right)^2 \right)
\end{aligned}
$$

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
</div>

We define:

$$
\begin{aligned}
r_{i,i} &= \omega_i^4 \delta_i^{-4} \kappa_{4, i} + 2 \omega_i^2 + e^2_i \kappa_{2, i} - 2 \omega_i^2 \delta_i^{-2} e_i \kappa_{3, i} \\
r_{i, i'} &= 2 \omega_i \omega_{i'} \\
c_i &= \omega_i^3 \delta_i^{-3} \kappa_{3,i} - \omega_i \delta_i^{-1} e_i \kappa_{2,i} \\
\mathbf{A}_j &= \tilde{\mathbf{Z}} \left. \left[ \frac{\partial D(\tau^2)}{\partial \tau^2_j} \right] \right\rvert_{\tau^2 = \mathbf{0}_m} \tilde{\mathbf{Z}}^\top
\end{aligned}
$$


<!-- #region example-r -->
<div class="example">
  <strong>Example.</strong>
  <br>
  First note that:
  $$
  \omega_i \delta_i^{-1} = \frac{\mu_i}{\mu_i \phi \left(1 + \frac{1}{\gamma} \mu_i\right)} = \frac{1}{\phi\left(1 + \frac{1}{\gamma} \mu_i \right)}
  $$
  Then:
  $$
  \begin{aligned}
  r_{i,i} &= \frac{\mu_i \left[ \frac{2}{\gamma} \mu_i \left(1 + \frac{1}{\gamma} \mu_i\right) + \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right]}{\phi\left(1 + \frac{1}{\gamma} \mu_i \right)^3} +  \frac{2 \mu_i^2}{\phi^2 \left( 1 + \frac{1}{\gamma} \mu_i\right)^2} + \frac{\mu_i^3}{\gamma^2 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^3} - \frac{2 \mu_i^2 \left(1 + \frac{2}{\gamma} \mu_i \right)}{\gamma^3 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^5}  \\
  r_{i,i'} &= \frac{2 \mu_i \mu_{i'}}{\phi^2\left(1 + \frac{1}{\gamma} \mu_i \right)\left(1 + \frac{1}{\gamma} \mu_{i'} \right)} \\
  c_i &= \frac{ \mu_i \left(1 + \frac{2}{\gamma} \mu_i \right)}{\phi \left(1 + \frac{1}{\gamma} \mu_i\right)^2}  \\
  \mathbf{A}_j &= \tilde{\mathbf{Z}} \tilde{\mathbf{Z}}^\top
  \end{aligned}
  $$

  <details>
  <summary>Proof.</summary>
  $$
  \begin{aligned}
  r_{i,i} &= \frac{\kappa_{4, i}}{\phi^4\left(1 + \frac{1}{\gamma} \mu_i \right)^4} +  \frac{2 \mu_i^2}{\phi^2 \left( 1 + \frac{1}{\gamma} \mu_i\right)^2} + \frac{\kappa_{2,i} \mu_i^2}{\phi^2 \gamma^2 \left(1 + \frac{1}{\gamma} \mu_i\right)^4} - \frac{2 \kappa_{3,i} \mu_i}{\phi^3 \gamma^3\left(1 + \frac{1}{\gamma} \mu_i\right)^6} \\
  &= \frac{\phi^3\left[ \frac{2}{\gamma} \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right)^2 + \left(\mu_i + \frac{1}{\gamma} \mu_i^2\right) \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right]}{\phi^4\left(1 + \frac{1}{\gamma} \mu_i \right)^4} +  \frac{2 \mu_i^2}{\phi^2 \left( 1 + \frac{1}{\gamma} \mu_i\right)^2} + \frac{\phi \left(\mu_i + \frac{1}{\gamma} \mu_i^2 \right) \mu_i^2}{\phi^2 \gamma^2 \left(1 + \frac{1}{\gamma} \mu_i\right)^4} - \frac{2 \phi^2 \mu_i \left(1 + \frac{1}{\gamma} \mu_i \right)\left(1 + \frac{2}{\gamma} \mu_i \right) \mu_i}{\phi^3 \gamma^3\left(1 + \frac{1}{\gamma} \mu_i\right)^6}  \\
  &= \frac{\left[ \frac{2}{\gamma} \mu_i^2 \left(1 + \frac{1}{\gamma} \mu_i\right)^2 + \mu_i \left(1 + \frac{1}{\gamma} \mu_i\right) \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right]}{\phi\left(1 + \frac{1}{\gamma} \mu_i \right)^4} +  \frac{2 \mu_i^2}{\phi^2 \left( 1 + \frac{1}{\gamma} \mu_i\right)^2} + \frac{\left(1 + \frac{1}{\gamma} \mu_i\right) \mu_i^3}{\gamma^2 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^4} - \frac{2 \mu_i^2 \left(1 + \frac{2}{\gamma} \mu_i \right)}{\gamma^3 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^5}  \\
  &= \frac{\mu_i \left(1 + \frac{1}{\gamma} \mu_i\right)\left[ \frac{2}{\gamma} \mu_i \left(1 + \frac{1}{\gamma} \mu_i\right) + \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right]}{\phi\left(1 + \frac{1}{\gamma} \mu_i \right)^4} +  \frac{2 \mu_i^2}{\phi^2 \left( 1 + \frac{1}{\gamma} \mu_i\right)^2} + \frac{\mu_i^3}{\gamma^2 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^3} - \frac{2 \mu_i^2 \left(1 + \frac{2}{\gamma} \mu_i \right)}{\gamma^3 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^5} \\
  &= \frac{\mu_i \left[ \frac{2}{\gamma} \mu_i \left(1 + \frac{1}{\gamma} \mu_i\right) + \left(1 + \frac{2}{\gamma} \mu_i \right)^2 \right]}{\phi\left(1 + \frac{1}{\gamma} \mu_i \right)^3} +  \frac{2 \mu_i^2}{\phi^2 \left( 1 + \frac{1}{\gamma} \mu_i\right)^2} + \frac{\mu_i^3}{\gamma^2 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^3} - \frac{2 \mu_i^2 \left(1 + \frac{2}{\gamma} \mu_i \right)}{\gamma^3 \phi \left(1 + \frac{1}{\gamma} \mu_i\right)^5} \\
  c_i &= \frac{\kappa_{3,i}}{\phi^3 \left(1 + \frac{1}{\gamma} \mu_i\right)^3} \\
  &= \frac{\phi^2 \mu_i \left(1 + \frac{1}{\gamma} \mu_i \right) \left(1 + \frac{2}{\gamma} \mu_i \right)}{\phi^3 \left(1 + \frac{1}{\gamma} \mu_i\right)^3} \\
  &= \frac{ \mu_i \left(1 + \frac{2}{\gamma} \mu_i \right)}{\phi \left(1 + \frac{1}{\gamma} \mu_i\right)^2} 
  \end{aligned}
  $$
  </details>
</div>
<!-- #endregion -->


We also define $$\mathbf{a}_j = \text{diag}(\mathbf{A}_j)$$ (which is $n$-dimensional) and has elements denoted by $$a^j_{i,i}$$. 


### Information
The Fisher information is the variance of the score.

$$
\begin{aligned}
\mathcal{I}_{\tau^2_j, \tau^2_k} &= \frac{1}{4} \left( \sum_{i = 1}^n \left[ a_{i,i}^j a_{i,i}^k r_{i,i} + \sum_{i \neq i'}^n a_{i,i'}^j a_{i,i'}^k r_{i,i'} \right] \right) \\
\mathcal{I}_{\alpha, \tau^2_j} &= \frac{1}{2} \sum_{i = 1}^n a_{i,i}^j c_i \mathbf{x}_i \\
\mathcal{I}_{\alpha, \alpha} &= \sum_{i = 1}^n \omega_i \mathbf{x}_i \mathbf{x}_i^\top
\end{aligned}
$$



Recall that the trace of the product of two $n \times n$ matrices can be written as a double summation:

$$
\text{tr}[AB] = \sum_{i = 1}^n \sum_{i' = 1}^n A_{i,i'} B_{i,i'}
$$

$$
\begin{aligned}
U_{\tau^2_j}(\hat{\theta}) &= \frac{1}{2} \text{tr}\left[ \tilde{\mathbf{Z}}^\top \left( \bar{\Omega}\bar{\Delta}^{-1}(\mathbf{y} - \bar{\mu})(\mathbf{y} - \bar{\mu})^\top \bar{\Delta}^{-1} \bar{\Omega} - \bar{\Xi} \right) \tilde{\mathbf{Z}} \left. \left[ \frac{\partial \tilde{D}(\tau^2)}{\partial \tau^2_j} \right] \right\rvert_{\theta = \hat{\theta}} \right] \\
&= \frac{1}{2} \text{tr}\left[\left( \bar{\Omega}\bar{\Delta}^{-1}(\mathbf{y} - \bar{\mu})(\mathbf{y} - \bar{\mu})^\top \bar{\Delta}^{-1} \bar{\Omega} - \bar{\Xi} \right) \tilde{\mathbf{Z}} \left. \left[ \frac{\partial \tilde{D}(\tau^2)}{\partial \tau^2_j} \right] \right\rvert_{\theta = \hat{\theta}}  \tilde{\mathbf{Z}}^\top  \right] \\
&= \frac{1}{2} \sum_{i = 1}^n \sum_{i' = 1}^n \left( \bar{\Omega}\bar{\Delta}^{-1}(\mathbf{y} - \bar{\mu})(\mathbf{y} - \bar{\mu})^\top \bar{\Delta}^{-1} \bar{\Omega} - \bar{\Xi} \right)_{i,i'} \left( \tilde{\mathbf{Z}} \left. \left[ \frac{\partial \tilde{D}(\tau^2)}{\partial \tau^2_j} \right] \right\rvert_{\theta = \hat{\theta}}  \tilde{\mathbf{Z}}^\top \right)_{i,i'}
\end{aligned}
$$

Let $$\zeta_{i,i'} = \left( \tilde{\mathbf{Z}} \left[ \frac{\partial \tilde{D}(\tau^2)}{\partial \tau^2_j} \right]\tilde{\mathbf{Z}}^\top \right)_{i,i'}$$, and let $$\bar{\psi}_i = \bar{\omega}_i \bar{\delta}_i^{-1}$$. We can then write this as:

$$
\begin{aligned}
U_{\tau^2_j}(\hat{\theta}) &= \frac{1}{2} \sum_{i = 1}^n \sum_{i' = 1}^n \left( \bar{\Omega}\bar{\Delta}^{-1}(\mathbf{y} - \bar{\mu})(\mathbf{y} - \bar{\mu})^\top \bar{\Delta}^{-1} \bar{\Omega} - \bar{\Xi} \right)_{i,i'} \left( \tilde{\mathbf{Z}} \left[ \frac{\partial \tilde{D}(\tau^2)}{\partial \tau^2_j} \right] \tilde{\mathbf{Z}}^\top \right)_{i,i'} \\
&= \frac{1}{2}\sum_{i = 1}^n \sum_{i' = 1}^n \left( \bar{\omega}_i \bar{\omega}_{i'} \bar{\delta}_i^{-1} \bar{\delta}_{i'}^{-1} (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) - \bar{\xi}_{i,i'} \right) \zeta_{i,i'} \\
&= \frac{1}{2} \left[ \sum_{i = 1}^n \sum_{i' = 1}^n \bar{\psi}_i \bar{\psi}_{i'} (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) \zeta_{i,i'} - \sum_{i = 1}^n \bar{\xi}_i \zeta_{i,i} \right] & \left(\Xi \text{ is diagonal}\right) \\
&= \frac{1}{2} \left[ \sum_{i = 1}^n \sum_{i' = 1}^n \bar{\psi}_i \bar{\psi}_{i'} (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) \zeta_{i,i'} - \sum_{i = 1}^n (\bar{\omega}_i + \bar{e}_i (y_i - \bar{\mu}_i)) \zeta_{i,i} \right] 
\end{aligned}
$$

This leads us to:

$$
\begin{aligned}
\mathcal{I}_{\tau^2_r, \tau^2_s} 
&= \mathbb{E}_{H_0} \left[ U_{\tau^2_r} U_{\tau^2_s} \right] \\
&= \frac{1}{4} \mathbb{E}_{H_0} \left[ \left( \sum_{i = 1}^n \sum_{i' = 1}^n \bar{\psi}_i \bar{\psi}_{i'} (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) \zeta_{i,i'} - \sum_{i = 1}^n (\bar{\omega}_i + \bar{e}_i (y_i - \bar{\mu}_i)) \zeta_{i,i} \right) \left( \sum_{j = 1}^n \sum_{j' = 1}^n \bar{\psi}_j \bar{\psi}_{j'} (y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) \zeta_{j,j'} - \sum_{j = 1}^n (\bar{\omega}_j + \bar{e}_j (y_j - \bar{\mu}_j)) \zeta_{j,j} \right) \right] \\
&= \frac{1}{4} \left[
  \underbrace{\sum_{i = 1}^n \sum_{i' = 1}^n \sum_{j = 1}^n  \sum_{j' = 1}^n \bar{\psi}_i \bar{\psi}_{i'} \bar{\psi}_j \bar{\psi}_{j'} \zeta_{i,i'} \zeta_{j, j'} \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'})(y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) \right]}_{(i)} \right. \\
&\hspace{5mm} -\left.  
  \underbrace{\sum_{i = 1}^n \sum_{i' = 1}^n \sum_{j = 1}^n  \bar{\psi}_i \bar{\psi}_{i'} \zeta_{i,i'} \zeta_{j,j} \left ( \bar{\omega}_j \mathbb{E}_{H_0} \left[ (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) \right] + \bar{e}_i \mathbb{E}_{H_0} \left[ (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'})(y_j - \bar{\mu}_j) \right] \right)}_{(ii)}  \right. \\
&\hspace{5mm} - \left.
  \underbrace{\sum_{j = 1}^n \sum_{j' = 1}^n \sum_{i = 1}^n \bar{\psi}_j \bar{\psi}_{j'} \zeta_{j,j'} \zeta_{i,i} \left( \bar{\omega}_i \mathbb{E}_{H_0} \left[ (y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) \right] + \bar{e}_i \mathbb{E}_{H_0} \left[ (y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) (y_i - \bar{\mu}_i) \right] \right)}_{(iii)} \right. \\
&\hspace{5mm} + \left. 
  \underbrace{\sum_{i = 1}^n \sum_{j = 1}^n \zeta_{i,i} \zeta_{j,j}\mathbb{E}_{H_0}\left[ (\bar{\omega}_i + \bar{e}_i(y_i - \bar{\mu}_i))(\bar{\omega}_j + \bar{e}_j(y_j - \bar{\mu}_j)) \right]}_{(iv)}
\right] \\
&= \frac{1}{4} \left[ 
  \sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 \nu_i^4 + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,i} \zeta_{j,j} \nu_i^2 \nu_j^2 + 2 \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,j}^2 \nu_i^2 \nu_j^2
\right. \\
&\hspace{5mm} - \left. 
  \sum_{i = 1}^n \bar{\psi}^2 \zeta_{i,i}^2 \left(\bar{\omega}_i \nu_i^2 + \bar{e}_i \nu_i^3 \right) + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \nu_i^2
\right. \\
&\hspace{5mm} - \left.
  \sum_{i = 1}^n \bar{\psi}_i^2 \zeta_{i,i}^2 \left(\bar{\omega}_i \nu_i^2 + \bar{e}_i \nu_i^3 \right) + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \nu^2_i
\right. \\
&\hspace{5mm} + \left.
  \sum_{i = 1}^n \zeta_{i,i}^2 \left(\omega_i^2 + \bar{e}_i^2 \nu_i^2 \right)
\right] \\
&= \frac{1}{4} \left[ 
  \sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 \nu_i^4 + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,i} \zeta_{j,j} \nu_i^2 \nu_j^2 + 2 \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,j}^2 \nu_i^2 \nu_j^2
\right. \\
&\hspace{5mm} - \left. 
  2\sum_{i = 1}^n \bar{\psi}^2 \zeta_{i,i}^2 \left(\bar{\omega}_i \nu_i^2 + \bar{e}_i \nu_i^3 \right)
- 2\sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \nu_i^2 
+ \sum_{i = 1}^n \zeta_{i,i}^2 \left(\omega_i^2 + \bar{e}_i^2 \nu_i^2 \right)
\right] \\
\end{aligned}
$$


Now, notice that $\nu_i^2 = V(\mu_i) = \kappa_{2,i}$, $\nu_i^3 = \kappa_{3,i}$, and $\nu_i^4 = \kappa_{4,i} + 3 \kappa_{2,i}^2$. Furthermore, we have that:

$$
\bar{\psi}_i^2 \nu_i^2 = \bar{\omega}_i^2 \bar{\delta}^{-2}_i V(\mu_i) = \bar{\omega}_i^2 \frac{V(\mu_i)}{\bar{\delta}_i^2} = \frac{\bar{\omega}_i^2}{\bar{\omega}_i} = \bar{\omega}_i
$$

Thus:

$$
\begin{aligned}
\mathcal{I}_{\tau^2_r, \tau^2_s} 
&= \frac{1}{4} \left[ 
  \sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 \nu_i^4 + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,i} \zeta_{j,j} \nu_i^2 \nu_j^2 + 2 \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,j}^2 \nu_i^2 \nu_j^2
\right. \\
&\hspace{5mm} - \left. 
  2\sum_{i = 1}^n \bar{\psi}^2 \zeta_{i,i}^2 \left(\bar{\omega}_i \nu_i^2 + \bar{e}_i \nu_i^3 \right)
- 2\sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \nu_i^2 
+ \sum_{i = 1}^n \zeta_{i,i}^2 \left(\omega_i^2 + \bar{e}_i^2 \nu_i^2 \right)
\right] \\
&= \frac{1}{4} \left[ \sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 (\kappa_{4,i} + 3 (V(\mu_i))^2) + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\omega}_i \bar{\omega}_j \zeta_{i,i} \zeta_{j,j} + 2 \sum_{i = 1}^n \sum_{j = 1}^n \bar{\omega}_i \bar{\omega}_j \zeta_{i,j}^2 \right. \\
&\hspace{5mm} - \left. 
  2\sum_{i = 1}^n \bar{\psi}^2 \zeta_{i,i}^2 \left(\bar{\omega}_i V(\mu_i) + \bar{e}_i \kappa_{3,i} \right)
- 2\sum_{i = 1}^n \sum_{j = 1}^n \bar{\omega}_i \bar{\omega}_j \zeta_{i,i} \zeta_{j,j} 
+ \sum_{i = 1}^n \zeta_{i,i}^2 \left(\omega_i^2 + \bar{e}_i^2 V(\mu_i) \right)
\right] \\
\end{aligned}
$$





---


Notice that we have four different indices: $i$, $i'$, $j$, and $j'$. There are four possible cases: all indices are equal, three are equal and one is different, two are equal and the others are different, two pairs are equal, and all are different. If we have a case where we are taking the expectation of the product of residuals and a single index is different from all of the others (the second, third, and last case), then the entire expectation will equal zero by the independence of observations. 
<br>
Though we only need to consider the remaining cases when evaluating $(i)$, $(ii)$, $(iii)$, and $(iv)$, we must be sure to remember that there are three possible ways to choose two pairs. 
<br>
Let $\nu_i^j$ denote the $j$-th central moment for observation $i$. We have:

$$
\begin{aligned}
(i) 
&= \sum_{i = 1}^n \sum_{i' = 1}^n \sum_{j = 1}^n \sum_{j' = 1}^n \bar{\psi}_i \bar{\psi}_{i'} \bar{\psi}_j \bar{\psi}_{j'} \zeta_{i,i'} \zeta_{j, j'} \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i) (y_{i'} - \bar{\mu}_{i'}) (y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) \right] \\
&= \underbrace{\sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i)^4 \right]}_{\text{i = i' = j = j'}} \\
&\hspace{5mm} + \underbrace{\sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,i} \zeta_{j,j} \mathbb{E}_{H_0} \left[(y_i - \bar{\mu}_i)^2 \right] \mathbb{E}_{H_0} \left[ (y_j - \bar{\mu}_j)^2 \right]}_{i = i'; j = j'} \\
&\hspace{5mm} + \underbrace{\sum_{i = 1}^n \sum_{i' = 1}^n \bar{\psi}_i^2 \bar{\psi}_{i'}^2 \zeta_{i,i'}^2 \mathbb{E}_{H_0} \left[(y_i - \bar{\mu}_i)^2 \right] \mathbb{E}_{H_0} \left[ (y_{i'} - \bar{\mu}_{i'})^2 \right]}_{i = j; i' = j'} \\
&\hspace{5mm} + \underbrace{\sum_{i = 1}^n \sum_{i' = 1}^n \bar{\psi}_i^2 \bar{\psi}_{i'}^2 \zeta_{i,i'} \zeta_{i',i} \mathbb{E}_{H_0} \left[(y_i - \bar{\mu}_i)^2 \right] \mathbb{E}_{H_0} \left[ (y_{i'} - \bar{\mu}_{i'})^2 \right]}_{i = j'; i' = j} \\
&= \sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i)^4 \right] \\
&\hspace{5mm} + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,i} \zeta_{j,j} \mathbb{E}_{H_0} \left[(y_i - \bar{\mu}_i)^2 \right] \mathbb{E}_{H_0} \left[ (y_j - \bar{\mu}_j)^2 \right]\\
&\hspace{5mm} + 2\sum_{i = 1}^n \sum_{i' = 1}^n \bar{\psi}_i^2 \bar{\psi}_{i'}^2 \zeta_{i,i'}^2 \mathbb{E}_{H_0} \left[(y_i - \bar{\mu}_i)^2 \right] \mathbb{E}_{H_0} \left[ (y_{i'} - \bar{\mu}_{i'})^2 \right] & \left(\zeta \text{ is symmetric}\right)\\
&= \sum_{i = 1}^n \bar{\psi}_i^4 \zeta_{i,i}^2 \nu_i^4 + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,i} \zeta_{j,j} \nu_i^2 \nu_j^2 + 2 \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \bar{\psi}_j^2 \zeta_{i,j}^2 \nu_i^2 \nu_j^2 & \left(\text{switch index}\right)
\end{aligned}
$$


$$
\begin{aligned}
(ii) 
&= \sum_{i = 1}^n \sum_{i' = 1}^n \sum_{j = 1}^n  \bar{\psi}_i \bar{\psi}_{i'} \zeta_{i,i'} \zeta_{j,j} \left ( \bar{\omega}_j \mathbb{E}_{H_0} \left[ (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'}) \right] + \bar{e}_i \mathbb{E}_{H_0} \left[ (y_i - \bar{\mu}_i)(y_{i'} - \bar{\mu}_{i'})(y_j - \bar{\mu}_j) \right] \right) \\
&= \underbrace{\sum_{i = 1}^n \bar{\psi}_i^2 \zeta_{i,i}^2 \left(\bar{\omega}_i \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i)^2 \right] + \bar{e}_i \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i)^3 \right]\right)}_{i=i'=j} \\
&\hspace{5mm} + \underbrace{\sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \mathbb{E}_{H_0}\left[(y_i - \bar{\mu}_i)^2 \right] }_{i = i' \neq j} \\
&= \sum_{i = 1}^n \bar{\psi}^2 \zeta_{i,i}^2 \left(\bar{\omega}_j \nu_i^2 + \bar{e}_i \nu_i^3 \right) + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \nu_i^2
\end{aligned}
$$


$$
\begin{aligned}
(iii) 
&= \sum_{j = 1}^n \sum_{j' = 1}^n \sum_{i = 1}^n \bar{\psi}_j \bar{\psi}_{j'} \zeta_{j,j'} \zeta_{i,i} \left( \bar{\omega}_i \mathbb{E}_{H_0} \left[ (y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) \right] + \bar{e}_i \mathbb{E}_{H_0} \left[ (y_j - \bar{\mu}_j)(y_{j'} - \bar{\mu}_{j'}) (y_i - \bar{\mu}_i) \right] \right) \\
&= \underbrace{\sum_{j = 1}^n \bar{\psi}_j^2 \zeta_{j,j}^2 \left(\bar{\omega}_j \mathbb{E}_{H_0}\left[ (\mathbf{y}_j - \bar{\mu}_j)^2 \right] + \bar{e}_j \mathbb{E}_{H_0}\left[ (y_j - \bar{\mu}_j)^3 \right] \right)}_{j = j' = i} \\
&\hspace{5mm} + \underbrace{\sum_{j = 1}^n \sum_{i = 1}^n \bar{\psi}_j^2 \zeta_{j,j}\zeta_{i,i} \bar{\omega}_i \mathbb{E}_{H_0}\left[ (y_j - \bar{\mu}_j)^2 \right]}_{j = j' \neq i} \\
&= \sum_{i = 1}^n \bar{\psi}_i^2 \zeta_{i,i}^2 \left(\bar{\omega}_i \nu_i^2 + \bar{e}_i \nu_i^3 \right) + \sum_{i = 1}^n \sum_{j = 1}^n \bar{\psi}_i^2 \zeta_{i,i} \zeta_{j,j} \bar{\omega}_j \nu^2_i
\end{aligned}
$$

$$
\begin{aligned}
(iv) 
&= \sum_{i = 1}^n \sum_{j = 1}^n \zeta_{i,i} \zeta_{j,j}\mathbb{E}_{H_0}\left[ (\bar{\omega}_i + \bar{e}_i(y_i - \bar{\mu}_i))(\bar{\omega}_j + \bar{e}_j(y_j - \bar{\mu}_j)) \right] \\
&= \underbrace{\sum_{i = 1}^n \zeta_{i,i}^2 \left( \bar{\omega}_i^2 + \bar{e}_i^2 \mathbb{E}_{H_0}\left[ (y_i - \bar{\mu}_i)^2 \right]\right)}_{i = j} \\
&= \sum_{i = 1}^n \zeta_{i,i}^2 \left(\omega_i^2 + \bar{e}_i^2 \nu_i^2 \right)
\end{aligned}
$$







---

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

