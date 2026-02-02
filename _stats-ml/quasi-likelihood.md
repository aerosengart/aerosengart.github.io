---
layout: distill
title: Quasi-Likelihood
description: A Primer
date: 2025-05-30
tabs: true
tags: glmm likelihood theory primer
toc:
    - name: Set-Up
    - name: Connections To Likelihood
    - name: Estimation
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: 2025-05-30-quasi-likelihood.bib
---

This post is a review of quasi-likelihood theory and mostly relies upon Wedderburn<d-cite key=wedderburn1974></d-cite> and Breslow & Clayton<d-cite key=breslow1993></d-cite>. Quasi-likelihood functions provide an alternative and less demanding way of characterizing the distribution of observations compared to specifying a true likelihood function. In essence, we simply assume a particular relationship between the mean and the variance rather than a particular distributional family. Then a so-called <i>quasi-likelihood</i> function can be defined and used for parameter estimation.

---
## Set-Up
Let $x_1, \dots, x_n$ denote our (independent) observations. Denote the expectation and variance of an arbitrary $x_i$ with $\mathbb{E}[x_i] = \mu_i$ and $\text{var}(x_i) = \phi V(\mu_i)$ where $\phi > 0$ is some scale parameter. $V(\mu_i)$ is some (known) function of the mean. 

It's important to note that the expectations and variances need not be identical, but we enforce that the variances are proportional to some (shared) function of the expectations. We assume that $\mu_i = g(\beta_1, \dots, \beta_m)$. That is, the expectations are some (known) function of parameters $\beta_1, \dots, \beta_m$. 

The quasi-likelihood is easier to explain after defining the <i>quasi-score</i> function.

<div class="definition">
<strong>Definition (Quasi-Score).</strong>
<br>
  The <i>quasi-score</i> of $x_i$ is given by:
  $$
  \begin{equation}
  \label{eq:quasi-score}
  U(x_i; \mu_i) = \frac{x_i - \mu_i}{\phi V(\mu_i)}
  \end{equation}
  $$
</div>

In likelihood theory, the score function is the gradient of the log-likelihood function with respect to the parameters. In a similar fashion, the quasi-score is the gradient of the quasi-likelihood function with respect to the mean. 

<div class="definition">
<strong>Definition (Quasi-Likelihood).</strong>
<br>
  The <i>quasi-likelihood</i> (or, more precisely, the <i>quasi log-likelihood</i>) of $x_i$ is given by:
  $$
  \begin{equation}
  \label{eq:quasi-likelihood}
  \begin{aligned}
  \ell_q(x_i; \mu_i) &= \int_{x_i}^{\mu_i} \frac{x_i - z}{\phi V(z)} dz \\
  &\iff \\
  \frac{\partial}{\partial \mu_i} [\ell_q(x_i; \mu_i)] &= U(x_i; \mu_i)
  \end{aligned}
  \end{equation}
  $$
</div>

The quasi-score satisfies several of the properties that the score in likelihood theory satisfies, which justifies its name as a <i>quasi</i> score. 

<div id="theorem-1"></div>
<div class="theorem">
<strong>Theorem 1.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs theorem-1-wedderburn %}
{% tab theorem-1-wedderburn theorem %}
Let $x$ be some observation with expectation $\mu$ and variance $\phi V(\mu)$ for some $\phi > 0$. Suppose $\mu = g(\beta_1, \dots, \beta_m)$ for some (continuous and differentiable) funciton $g(\cdot)$. The quasi-score and quasi-likelihood, as defined in Eqs. \eqref{eq:quasi-score} and \eqref{eq:quasi-likelihood}, satisfy the following properties:

<ol>
<li>$\mathbb{E}\left[ U(x; \mu) \right] = 0$</li>
<li>$\mathbb{E}\left[ \frac{\partial \ell_q(x; \mu) }{\partial \beta_i} \right] = 0$ for all $i \in [m]$</li>
<li>$\text{Var}(U(x; \mu)) = - \mathbb{E}\left[ \frac{\partial^2 \ell_q(x;\mu)}{\partial \mu^2}\right] = \frac{1}{\phi V(\mu)}$</li>
<li>$\mathbb{E}\left[ \frac{\partial \ell_q(x; \mu)}{\partial \beta_i} \frac{\partial \ell_q(x; \mu)}{\partial \beta_j}\right] = - \mathbb{E}\left[ \frac{\partial^2 \ell_q(x; \mu)}{\partial \beta_i \partial \beta_j} \right] = \frac{1}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}$</li>
</ol>
{% endtab %}
{% tab theorem-1-wedderburn proof %}
Proving (1) is simple:

$$
\mathbb{E}\left[ U(x; \mu) \right] 
= \mathbb{E}\left[ \frac{x - \mu}{\phi V(\mu)} \right] 
= \frac{1}{\phi V(\mu)} \mathbb{E}\left[ x  - \mu\right] 
= \frac{1}{\phi V(\mu)} (\mu - \mu)
= 0
$$

To show (2), we note that:

$$
\frac{\partial \ell_q(x; \mu)}{\partial \beta_i} = \frac{\partial \ell_q(x; \mu)}{\mu} \frac{\partial \mu}{\beta_j}
\hspace{5mm} \implies \hspace{5mm}
\mathbb{E}\left[ \frac{\partial \ell_q(x; \mu) }{\partial \beta_i}  \right]
= \mathbb{E}\left[ \frac{\partial \ell_q(x; \mu)}{\partial \mu} \frac{\partial \mu}{\beta_i} \right]
= \frac{\partial \mu}{\beta_i} \mathbb{E}\left[ U(x; \mu) \right] 
= 0
\nonumber
$$

Showing (3) is also relatively easy. We first show that $\text{Var}(U(x;\mu)) = \frac{1}{\phi V(\mu)}$:

$$
\text{Var}(U(x; \mu)) 
= \text{Var}\left(\frac{x - \mu}{\phi V(\mu)}\right)
= \frac{1}{\phi^2 V^2(\mu)} \text{Var}(x - \mu) 
= \frac{1}{\phi^2 V^2(\mu)} \text{Var}(x)
= \frac{\phi V(\mu)}{\phi^2 V^2(\mu)}
= \frac{1}{\phi V(\mu)}
$$

Next, we show that $- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x;\mu)}{\partial \mu^2}\right] = \frac{1}{\phi V(\mu)}$:

$$
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x;\mu)}{\partial \mu^2}\right]
= - \mathbb{E}\left[ \frac{\partial U(x; \mu)}{\partial \mu} \right] 
= - \mathbb{E}\left[ \frac{\partial}{\partial \mu} \left[ \frac{x - \mu}{\phi V(\mu)}\right] \right]
= - \mathbb{E}\left[ \frac{\phi V(\mu)(-1) - (x- \mu)(\phi V'(\mu))}{\phi^2 V^2(\mu)} \right]
= \frac{1}{\phi V(\mu)} + \frac{(\mathbb{E}[x]- \mu)(\phi V'(\mu))}{\phi^2 V^2(\mu)}
= \frac{1}{\phi V(\mu)}
$$

For (4), we first show that $\mathbb{E}\left[ \frac{\partial \ell_q(x; \mu)}{\partial \beta_i} \frac{\partial \ell_q(x; \mu)}{\partial \beta_j}\right] = \frac{1}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}$:

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial \ell_q(x; \mu)}{\partial \beta_i} \frac{\partial \ell_q(x; \mu)}{\partial \beta_j}\right]
&= \mathbb{E}\left[ \frac{\partial \ell_q(x;\mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \ell_q(x;\mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_j}\right] \\
&= \mathbb{E}\left[ \left( \frac{\partial \ell_q(x;\mu)}{\partial \mu} \right)^2 \right]  \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \mathbb{E}\left[ U^2(x; \mu) \right] \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \mathbb{E}\left[ \frac{(x - \mu)^2}{\phi^2 V^2(\mu)} \right] \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \frac{\phi V(\mu)}{\phi^2 V(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \frac{1}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}
\end{aligned}
$$

Then we show that $- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x; \mu)}{\partial \beta_i \partial \beta_j} \right] = \frac{1}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}$:

$$
\begin{aligned}
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x; \mu)}{\partial \beta_i \partial \beta_j} \right] 
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[ \frac{\partial \ell_q(x; \mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_i}\right] \right] \\
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[ U(x; \mu) \frac{\partial \mu}{\partial \beta_i}\right] \right] \\
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[ \frac{x - \mu}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_i}\right] \right] \\
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[\frac{x-\mu}{\phi V(\mu)} \right] \frac{\partial \mu}{\partial \beta_i} \right] + \underbrace{\mathbb{E}\left[ \frac{x - \mu}{\phi V(\mu)}\frac{\partial}{\partial \beta_j} \left[\frac{\partial \mu}{\partial \beta_i} \right] \right]}_{=0}  \\
&= - \mathbb{E}\left[ \frac{\phi V(\mu) \frac{\partial}{\partial \beta_j}[x - \mu] - (x-\mu) \frac{\partial}{\partial \beta_j}[\phi V(\mu)]}{\phi^2 V^(\mu)}\right] \frac{\partial \mu}{\partial \beta_i} \\
&= \left(- \mathbb{E}\left[ - \frac{1}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_j} \right] + \underbrace{\mathbb{E}\left[ x- \mu\right]}_{=0} \frac{\frac{\partial}{\partial \beta_j}[\phi V(\mu)]}{\phi^2 V^2(\mu)} \right) \frac{\partial \mu}{\partial \beta_i} \\
&= \frac{1}{\phi V(\mu)} \frac{\partial \mu}{\partial \beta_j} \frac{\partial \mu}{\partial \beta_i}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

---

## Connections To Likelihood

Suppose that the distribution of $x$ is a function of $\mu$ such that a log-likelihood can be explicitly written. Let $\ell(z; \mu)$ denote this log-likelihood. The following property is due to the above theorem:

<div class="theorem">
<strong>Claim.</strong>
{% tabs claim-lik1 %}
{% tab claim-lik1 statement %}

$$
\begin{equation}
\label{eq:corollary-1}
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x; \mu)}{\partial \mu^2} \right] \leq - \mathbb{E}\left[ \frac{\partial^2 \ell(x; \mu)}{\partial \mu^2} \right]
\end{equation}
$$
{% endtab %}
{% tab claim-lik1 proof %}
By (4) of <a href="#theorem-1">Theorem 1</a>:

$$
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x; \mu)}{\partial \mu^2}\right] = \frac{1}{\phi V(\mu)}
$$

Our problem then becomes showing that:

$$
\frac{1}{\phi V(\mu)} \leq - \mathbb{E}\left[\frac{\partial^2 \ell(x; \mu)}{\partial \mu^2} \right]
\hspace{5mm} \iff \hspace{5mm}
\phi V(\mu) \geq - \frac{1}{\mathbb{E}\left[\frac{\partial^2 \ell(x; \mu)}{\partial \mu^2} \right]}
$$

Under certain regularity conditions (see <a href="/posts/2025/02/03/likelihood-theory.html">my likelihood post</a>), we have that $-\mathbb{E}\left[\frac{\partial^2 \ell(x; \mu)}{\partial \mu^2} \right]$ is the Fisher information. The result follows directly from the <a href="https://en.wikipedia.org/wiki/Cramér–Rao_bound">Cramér-Rao bound</a>.
{% endtab %}
{% endtabs %}
</div>

Wedderburn provides an additional connection between quasi-likelihood and likelihood functions for one-parameter distributions specified by the mean. 

<div id="theorem-2"></div>
<div class="theorem">
<strong>Theorem 2.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs theorem-2-wedderburn %}
{% tab theorem-2-wedderburn statement %}
Let $x$ be some observation with expectation $\mu$ and variance $\phi V(\mu)$ for some $\phi > 0$. Suppose $\mu = g(\beta_1, \dots, \beta_m)$ for some (continuous and differentiable) function $g(\cdot)$. The log-likelihood function, $\ell(x; \mu)$, for $x$ satisfies:

$$
\begin{equation}
\label{eq:ll-condition}
\frac{\partial}{\partial \mu} \left[ \ell(x; \mu) \right] = \frac{x - \mu}{\phi V(\mu)}
\end{equation}
$$

if and only if the density function of $x$ can be written, with respect to some measure, as:

$$
f_x = \exp\left( x \theta - h(\theta) \right)
$$
{% endtab %}
{% tab theorem-2-wedderburn proof %}
We first prove the forwards direction. Assume the log-likelihood satisfies Eq. \eqref{eq:ll-condition}. We integrate with respect to $\mu$:

$$
\begin{aligned}
&\int \frac{\partial}{\partial \mu} \left[ \ell(x; \mu) \right] d\mu = \int \frac{x - \mu}{\phi V(\mu)} d\mu  \\
\implies
&\ell(x; \mu) = \frac{x}{\phi} \int \frac{1}{V(\mu)} d \mu - \frac{1}{\phi}\int \frac{\mu}{V(\mu)} d\mu
\end{aligned}
$$

Substituting in $\theta = \frac{1}{\phi} \int \frac{1}{V(\mu)} d\mu$:

$$
\ell(x; \mu) = x \theta - ?
$$
{% endtab %}
{% endtabs %}
</div>

The theorem can be summarized quite nicely: the quasi-likelihood function will equal the log-likelihood function <i>if and only if</i> the distribution comes from an exponential family. 

Extending the previous corollary, we see that for a one-parameter exponential family, Eq. \eqref{eq:corollary-1} obtains equality. Under certain regularity conditions (see <a href="/stats-ml/likelihood-theory">my likelihood post</a>), $-\mathbb{E}\left[\frac{\partial^2 \ell(x; \mu)}{\partial \mu^2} \right]$ is the Fisher information, which describes the amount of information about $\mu$ that is held in $x$.

Since equality is obtained, we can also think of $-\mathbb{E}\left[\frac{\partial^2 \ell_q(x; \mu)}{\partial \mu^2} \right]$ as describing the amount of information about $\mu$ that is held in $x$. In addition, the difference between the former and the latter can be thought of as the amount of information gained by knowing, specifically, the distribution of $z$. 

---

## Estimation
Let $x_{1:n}= (x_1, \dots, x_n)$ for independent observations $x_1, \dots, x_n$, and let $\mu_{1:n} = (\mu_1, \dots, \mu_n)$. We'll denote the gradient of the (full) quasi-likelihood with respect to the parameters $\beta_1, \dots, \beta_m$ with:

$$
\mathbf{u} = \frac{\partial \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta} = 
\begin{bmatrix}
\frac{\partial \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta_1} \\
\vdots \\
\frac{\partial \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta_m}
\end{bmatrix}
\label{eq:gradient-ql}
$$

By <a href="#theorem-1">Theorem 1</a>, $\mathbf{u}$ has mean vector:

$$
\begin{equation}
\label{eq:u-mean}
\begin{aligned}
\mathbb{E}[\mathbf{u}]
&= \mathbb{E}\left[ \frac{\partial \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta} \right]  \\
&= \mathbb{E}\left[ 
    \begin{bmatrix}
    \sum_{i = 1}^n  \frac{\partial \ell_q(x_i; \mu_i)}{\partial \beta_1} \\
    \vdots \\
    \sum_{i = 1}^n  \frac{\partial \ell_q(x_i; \mu_i)}{\partial \beta_m} \\
    \end{bmatrix}
\right] \\
&= \sum_{i = 1}^n \mathbb{E}\left[ 
    \begin{bmatrix}
    \frac{\partial \ell_q(x_i; \mu_i)}{\partial \beta_1} \\
    \vdots \\
    \frac{\partial \ell_q(x_i; \mu_i)}{\partial \beta_m} \\
    \end{bmatrix}
\right]\\
&= \mathbf{0}
\end{aligned}
\end{equation}
$$

and covariance matrix: 

$$
\begin{equation}
\label{eq:u-cov}
\begin{aligned}
\text{Cov}(\mathbf{u})
&= \mathbb{E}\left[ \frac{\partial \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta} \frac{\partial \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta^\top}\right]  \\
&= - \mathbb{E}\left[\frac{\partial^2 \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta \partial \beta^\top}\right]
\end{aligned}
\end{equation}
$$

The <i>maximum quasi-likelihood estimates</i> of $\beta$, denoted by $\hat{\beta}$, are found by setting $\mathbf{u}$ equal to $\mathbf{0}$ and solving for $\beta$, just like we would do for maximum likelihood estimation. 

<div id="theorem-3"></div>
<div class="theorem">
<strong>Theorem 3.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs theorem-3-wedderburn %}
{% tab theorem-3-wedderburn theorem %}
Let $x$ be some observation with expectation $\mu$ and variance $\phi V(\mu)$ for some $\phi > 0$. Suppose $\mu = g(\beta_1, \dots, \beta_m)$ for some (continuous and differentiable) function $g(\cdot)$. 

Denote the gradient of the (full) quasi-likelihood with respect to the parameters $\beta_1, \dots, \beta_m$ with $\mathbf{u}$, and let $\hat{\beta}$ be the maximum quasi-likelihood estimates of $\beta$. The mean of $\hat{\beta}$ is approximately $\mathbf{0}$, and the covariance of $\hat{\beta}$ is approximately:

$$
\text{Cov}(\hat{\beta}) \approx \text{Cov}^{-1}(\mathbf{u}) = \left[ -\mathbb{E}\left[ \frac{\partial^2 \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta \partial \beta^\top} \right] \right]^{-1}
$$

if $\phi = 1$.
{% endtab %}
{% tab theorem-3-wedderburn proof %}
Let $\hat{\mathbf{u}}$ denote the gradient vector evaluated at the maximum quasi-likelihood estimate, $\hat{\beta}$. Since $\hat{\beta}$ is the value of $\beta$ such that $\mathbf{u}$ equals $\mathbf{0}$, a first-order Taylor approximation of $\mathbf{u}$ gives us:

$$
\begin{aligned}
\mathbf{u} &\approx \hat{\mathbf{u}} + \frac{\partial \mathbf{u}}{\partial \beta} (\beta - \hat{\beta}) \\
&= \frac{\partial^2 \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta \partial \beta^\top} (\beta - \hat{\beta}) \\
\implies
\beta - \hat{\beta} &\approx \left[ \frac{\partial^2 \ell_q(x_{1:n}; \mu_{1:n})}{\partial \beta \partial \beta^\top} \right]^{-1} \mathbf{u}
\end{aligned}
$$

If we approximate the inverted matrix by its expectation, whose elements are given in <a href="#theorem-1">Theorem 1</a>, we get:

$$
\begin{aligned}
&\beta - \hat{\beta} \approx  -\text{Cov}^{-1}(\mathbf{u})\mathbf{u} \\
\implies
&\hat{\beta} \approx \beta + \text{Cov}^{-1}(\mathbf{u})\mathbf{u} 
\end{aligned}
$$

Since $\mathbf{u}$ has expectation zero, it is clear that $\mathbb{E}[\hat{\beta}] \approx \mathbf{0}$ as well.

The first term on the right-hand side of the above expression is fixed (the true parameter value), so the <i>approximate</i> covariance matrix of $\hat{\beta}$ is:

$$
\begin{aligned}
\text{Cov}(\hat{\beta}) &\approx \text{Cov}\left( \beta + \text{Cov}^{-1}(\mathbf{u})\mathbf{u} \right) \\
&= \mathbb{E}\left[ \text{Cov}^{-1}(\mathbf{u})\mathbf{u} \left(\text{Cov}^{-1}(\mathbf{u})\mathbf{u}\right)^\top \right] \\
&= \text{Cov}^{-1}(\mathbf{u})\mathbb{E}\left[ \mathbf{u} \mathbf{u}^\top\right] \text{Cov}^{-1}(\mathbf{u}) \\
&= \text{Cov}^{-1}(\mathbf{u})
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

The above theorem holds only for a scale parameter $\phi = 1$. If we relax this assumption to $\phi > 0$, the expectation does not change, but we need to estimate $\phi$ before we can approximate the covariance of the maximum quasi-likelihood estimates. 
