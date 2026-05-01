---
layout: distill
title: Quasi-Likelihood Functions, Generalized Linear Models, and the Gauss-Newton Method
description: 
date: 2026-03-16
tabs: true
tags: likelihood theory
toc:
    - name: Quasi-Likelihood
      subsections:
        - name: Properties
    - name: Estimation
      subsection:
        - name: Weighted Least Squares

bibliography: paper-notes.bib
---

(Approximate) inference for generalized linear mixed model often relies upon quasi-likelihood functions for the given model. In this post, I go through <i>Quasi-likelihood functions, generalized linear models, and the Gauss-Newton method</i> by R. W. M. Wedderburn to get a better sense of the theory behind these functions.<d-cite key=wedderburn1974></d-cite>

The motivation behind Wedderburn's work is estimation of regression model parameters under relatively minimal assumptions on the distribution of the responses. Wedderburn introduces the idea of the <strong>quasi-likelihood function</strong>, which is an analogue of the likelihood function for fitting models in which only the relationship between the mean and variance is specified. 

Note that I have already covered much of this paper in my <a href="/stats-ml/quasi-likelihood/">quasi-likelihood post</a>, but I've reiterated the results here. 

---

## Quasi-Likelihood
We assume to have $n$ independent observations, $\mathbf{z} = (z_1, \dots, z_n)^\top$, each having expectation, $\mu_i$, and variance, $v(\mu_i) = \phi V(\mu_i)$, with $\phi > 0$. We assume $V(\cdot)$ is a known function of the mean and that $\mu_i = g(\beta_1, \dots, \beta_m)$. That is, the expectations are some (known) function of parameters $\beta_1, \dots, \beta_m$. 

It's important to note that the expectations and variances need not be identical, but we enforce that the variances are proportional to some (shared) function of the expectations. 

The quasi-likelihood funtion, or, perhaps, more accurately the <i>log</i> quasi-likelihood function, can be defined as followed.

<div class="definition">
<strong>Definition (Quasi-Likelihood).</strong>
<br>
  The <i>quasi-likelihood</i> of $z_i$, denoted by $\ell_q(z_i; \mu_i)$, satisfies:
  $$
  \begin{equation}
  \label{eq:quasi-score}
  \frac{\partial \ell_q(z_i; \mu_i)}{\partial \mu_i} = \frac{z_i - \mu_i}{v(\mu_i)}
  \end{equation}
  $$
</div>

<aside><p>In the original publication, $\ell_q(z_i; \mu_i)$ is denoted by $K(z_i, \mu_i)$.</p></aside>

Eq. \eqref{eq:quasi-score} is often referred to as the <strong>quasi-score</strong>, which we denote with $U_q(z_i; \mu_i)$. It implies that the quasi-likelihood also satisfies:

$$
\ell_q(z_i; \mu_i) = \int_{-\infty}^{\mu_i} \frac{z_i - u}{v(u)} du + h(z_i)
$$

for some function of $z_i$, $h(\cdot)$. 

### Properties 
The quasi-likelihood satisfies many properties that are also satisfied by log-likelihood functions.

<!-- #region ql-props -->
<div id="theorem-1"></div>
<div class="theorem">
<strong>Theorem 1.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs ql-props %}
{% tab ql-props statement %}
Let $z$ be an arbitrary observation, and let $\mu = g(\beta_1, \dots, \beta_m)$ be its expectation. The quasi-likelihood satisfies the following properties:

<ol>
<li>$\mathbb{E}\left[\frac{\partial \ell_q(z; \mu)}{\partial \mu} \right] = 0$</li>
<li>$\mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \beta_i} \right] = 0$</li>
<li>$\left( \mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \mu}\right] \right)^2 = -\mathbb{E}\left[ \frac{\partial^2 \ell_q(z; \mu)}{\partial \mu^2}\right] = \frac{1}{v(\mu)}$</li>
<li>$\mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \beta_i} \frac{\partial \ell_q(z; \mu)}{\partial \beta_j} \right] = - \mathbb{E}\left[ \frac{\partial^2 \ell_q(z; \mu)}{\partial \beta_i \partial \beta_j} \right] = \frac{1}{v(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}$</li>
</ol>

where all quantities have been evaluated at the true parameter values.
{% endtab %}
{% tab ql-props proof %}
Proving (1) is simple:

$$
\mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \mu} \right] 
= \mathbb{E}\left[ \frac{z - \mu}{v(\mu)} \right] 
= \frac{1}{v(\mu)} \mathbb{E}\left[ z  - \mu\right] 
= 0
$$

To show (2), we note that:

$$
\begin{aligned}
&\frac{\partial \ell_q(z; \mu)}{\partial \beta_i} = \frac{\partial \ell_q(z; \mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_i} \\
\implies &\mathbb{E}\left[ \frac{\partial \ell_q(z; \mu) }{\partial \beta_i}  \right]
= \mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \mu} \right] \frac{\partial \mu}{\beta_i} 
= 0
\end{aligned}
$$

Showing (3) is also relatively easy. We first show that $\text{Var}\left( \frac{\partial \ell_q(z; \mu)}{\partial \mu}\right) = \frac{1}{v(\mu)}$:

$$
\text{Var}\left( \frac{\partial \ell_q(z; \mu)}{\partial \mu}\right)
= \text{Var}\left(\frac{z - \mu}{v(\mu)}\right)
= \frac{\text{Var}(z)}{v^2(\mu)} 
= \frac{1}{v(\mu)}
$$

Next, we show that $- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x;\mu)}{\partial \mu^2}\right] = \frac{1}{v(\mu)}$:

$$
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(x;\mu)}{\partial \mu^2}\right]
= - \mathbb{E}\left[ \frac{\partial}{\partial \mu} \left[ \frac{z - \mu}{v(\mu)}\right] \right]
= - \mathbb{E}\left[ \frac{-v(\mu)-(z - \mu)v'(\mu)}{v^2(\mu)}\right] \\
= \frac{1}{v(\mu)}
$$

For (4), we first show that $\mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \beta_i} \frac{\partial \ell_q(z; \mu)}{\partial \beta_j}\right] = \frac{1}{v(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}$:

$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial \ell_q(z; \mu)}{\partial \beta_i} \frac{\partial \ell_q(z; \mu)}{\partial \beta_j}\right]
&= \mathbb{E}\left[ \frac{\partial \ell_q(z;\mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \ell_q(z;\mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_j}\right] \\
&= \mathbb{E}\left[ \left( \frac{\partial \ell_q(z;\mu)}{\partial \mu} \right)^2 \right]  \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \mathbb{E}\left[ \frac{(z - \mu)^2}{v^2(\mu)} \right] \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \frac{Var(z)}{v^2(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j} \\
&= \frac{1}{v(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}
\end{aligned}
$$

Then we show that $- \mathbb{E}\left[ \frac{\partial^2 \ell_q(z; \mu)}{\partial \beta_i \partial \beta_j} \right] = \frac{1}{v(\mu)} \frac{\partial \mu}{\partial \beta_i} \frac{\partial \mu}{\partial \beta_j}$:

$$
\begin{aligned}
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(z; \mu)}{\partial \beta_i \partial \beta_j} \right] 
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[ \frac{\partial \ell_q(z; \mu)}{\partial \mu} \frac{\partial \mu}{\partial \beta_i}\right] \right] \\
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[ \frac{z - \mu}{v(\mu)} \frac{\partial \mu}{\partial \beta_i}\right] \right] \\
&= - \mathbb{E}\left[ \frac{\partial}{\partial \beta_j} \left[\frac{z-\mu}{v(\mu)} \right] \frac{\partial \mu}{\partial \beta_i} \right] + \underbrace{\mathbb{E}\left[ \frac{z - \mu}{v(\mu)}\frac{\partial}{\partial \beta_j} \left[\frac{\partial \mu}{\partial \beta_i} \right] \right]}_{=0}  \\
&= - \mathbb{E}\left[ \frac{v(\mu) \frac{\partial}{\partial \beta_j}[z - \mu] - (z-\mu) \frac{\partial}{\partial \beta_j}[v(\mu)]}{v^2(\mu)}\right] \frac{\partial \mu}{\partial \beta_i} \\
&= \left(- \mathbb{E}\left[ - \frac{1}{v(\mu)} \frac{\partial \mu}{\partial \beta_j} \right] + \underbrace{\mathbb{E}\left[ z- \mu\right]}_{=0} \frac{\frac{\partial}{\partial \beta_j}[v(\mu)]}{v^2(\mu)} \right) \frac{\partial \mu}{\partial \beta_i} \\
&= \frac{1}{v(\mu)} \frac{\partial \mu}{\partial \beta_j} \frac{\partial \mu}{\partial \beta_i}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Suppose that the distribution of $z$ is a function of $\mu$ such that a log-likelihood can be explicitly written. Let $\ell(z; \mu)$ denote this log-likelihood. The following property is due to the above theorem:

<!-- #region claim-lik1 -->
<div class="theorem">
<strong>Claim.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs claim-lik1 %}
{% tab claim-lik1 statement %}

$$
\begin{equation}
\label{eq:corollary-1}
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(z; \mu)}{\partial \mu^2} \right] \leq - \mathbb{E}\left[ \frac{\partial^2 \ell(z; \mu)}{\partial \mu^2} \right]
\end{equation}
$$
{% endtab %}
{% tab claim-lik1 proof %}
By (4) of <a href="#theorem-1">Theorem 1</a>:

$$
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(z; \mu)}{\partial \mu^2}\right] = \frac{1}{v(\mu)}
$$

Our problem then becomes showing that:

$$
\frac{1}{v(\mu)} \leq - \mathbb{E}\left[\frac{\partial^2 \ell(z; \mu)}{\partial \mu^2} \right]
\hspace{5mm} \iff \hspace{5mm}
v(\mu) \geq - \frac{1}{\mathbb{E}\left[\frac{\partial^2 \ell(z; \mu)}{\partial \mu^2} \right]}
$$

Under certain regularity conditions (see <a href="/stats-ml/likelihood-theory">my likelihood post</a>), we have that $-\mathbb{E}\left[\frac{\partial^2 \ell(z; \mu)}{\partial \mu^2} \right]$ is the Fisher information. The result follows directly from the <a href="https://en.wikipedia.org/wiki/Cramér–Rao_bound">Cramér-Rao bound</a>.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

One can also show that, for certain one-parameter distributions, the quasi-likelihod is equivalent to the log-likelihood. 

<!-- #region ql-lik -->
<div id="theorem-2"></div>
<div class="theorem">
<strong>Theorem 2.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs theorem-2-wedderburn %}
{% tab theorem-2-wedderburn statement %}
Let $z$ be some observation with expectation $\mu$ and variance $v(\mu) = \phi V(\mu)$ for some $\phi > 0$. Suppose $\mu = g(\beta_1, \dots, \beta_m)$ for some (continuous and differentiable) function $g(\cdot)$. The log-likelihood function, $\ell(z; \mu)$, for $z$ satisfies:

$$
\begin{equation}
\label{eq:ll-condition}
\frac{\partial \ell(z; \mu)}{\partial \mu} = \frac{z - \mu}{v(\mu)}
\end{equation}
$$

if and only if the density function of $z$ can be written, with respect to some measure, as:

$$
f_z = \exp\left( z \theta - h(\theta) \right)
$$
{% endtab %}
{% tab theorem-2-wedderburn proof %}
We first prove the forwards direction. Assume the log-likelihood satisfies Eq. \eqref{eq:ll-condition}. We integrate with respect to $\mu$:

$$
\begin{aligned}
&\int \frac{\partial \ell(z; \mu)}{\partial \mu} d\mu = \int \frac{z - \mu}{v(\mu)} d\mu  \\
\implies
&\ell(z; \mu) = z \int \frac{1}{v(\mu)} d \mu - \int \frac{\mu}{v(\mu)} d\mu
\end{aligned}
$$

Substituting in $\theta = \int \frac{1}{v(\mu)} d\mu$:

$$
\ell(z; \mu) = z \theta - \int \frac{\mu}{v(\mu)} d\mu
$$

Now, the backwards direction. Suppose that, for measure $m$ on the real line, the density function of $z$ can be written as $\exp\left( z \theta - g(\theta) \right)$. Then:

$$
\begin{aligned}
1 &= \int \exp\left( z \theta - g(\theta) \right) dm(z) \\
\implies \exp\left(g(\theta)\right) &= \int \exp\left(z \theta \right) dm(z)
\end{aligned}
$$

The moment generating function of $z$ is then given by:

$$
\begin{aligned}
M_z(t) &= \mathbb{E}\left[ \exp\left( tz \right) \right] \\
&= \int \exp\left(tz\right) \exp\left(z\theta - g(\theta) \right) dm(z) \\
&= \exp\left(-g(\theta)\right) \int \exp\left(z (\theta + t) \right) dm(z) \\
&= \exp\left(g(\theta + t) -g(\theta)\right)
\end{aligned}
$$

The cumulant generating function of $z$ is then:

$$
K_z(t) = \log(M_z(t)) = g(\theta + t) - g(\theta)
$$

The first- and second-order derivatives of $K_z(t)$ give us the first and second central moments (i.e. the mean and variance).

$$
\begin{aligned}
\mathbb{E}\left[ z - \mu \right] &= \left. \frac{\partial K_z(t)}{\partial t} \right\rvert_{t = 0} = \left. \frac{\partial g(\theta + t)}{\partial t} \right\rvert_{t=0} \\
\mathbb{E}\left[ (z - \mu)^2 \right]  &= \left. \frac{\partial^2 K_z(t)}{\partial t^2} \right\rvert_{t = 0 } = \left. \frac{\partial^2 g(\theta + t)}{\partial t^2} \right\rvert_{t = 0} 
\end{aligned}
$$

Letting $s = \theta + t$, we see that:

$$
\begin{aligned}
\left. \frac{\partial g(\theta + t)}{\partial t} \right\rvert_{t=0} 
&= \left. \frac{\partial s}{\partial t} \frac{\partial g(\theta + t)}{\partial s} \right\rvert_{t = 0} \\
&= \left. \frac{\partial g(s)}{\partial s} \right \rvert_{s = \theta} \\
&= g'(\theta) \\
\left. \frac{\partial^2 g(\theta + t)}{\partial t^2} \right\rvert_{t = 0} 
&= \left. \frac{\partial}{\partial t} \left[ \frac{\partial g(\theta + t)}{\partial t} \right] \right\rvert_{t = 0}  \\
&= \left. \frac{\partial s}{\partial t} \frac{\partial}{\partial s}\left[ \frac{\partial s}{\partial t} \frac{\partial g(s)}{\partial s} \right] \right\rvert_{s = \theta}  \\
&= \left. \frac{\partial^2 g(s)}{\partial s^2} \right \rvert_{s = \theta} \\
&= g''(\theta)
\end{aligned}
$$

It also follows that:

$$
\begin{aligned}
\frac{\partial \mu}{\partial \theta} &= \frac{\partial}{\partial \theta} \left[ \left. \frac{\partial g(t)}{\partial t} \right\rvert_{t=\theta} \right]\\
&= g''(\theta) \\
&= \mathbb{E}\left[ (z - \mu)^2 \right] \\
&= v(\mu)
\end{aligned}
$$

Using the above results, we see that:

$$
\begin{aligned}
\frac{\partial \ell(z; \mu)}{\partial \mu}
&= \frac{\partial \theta}{\partial \mu} \frac{\partial \ell(z; \mu)}{\partial \theta}  \\
&= \frac{\partial \theta}{\partial \mu} \frac{\partial}{\partial \theta} \left[ z\theta - g(\theta) \right] \\
&= \left( z - g'(\theta) \right) \frac{\partial \theta}{\partial \mu} \\
&\overset{(i)}{=} \frac{z - \mu}{v(\mu)}
\end{aligned}
$$

where $(i)$ holds if $\frac{\partial \mu}{\partial \theta} = \frac{1}{\frac{\partial \theta}{\partial \mu}}$, which will be true if $g'(\theta)$ is invertible (from $g'(\theta) = \mu$ above). 
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

The theorem can be summarized quite nicely: the quasi-likelihood function will equal the log-likelihood function <i>if and only if</i> the distribution comes from an exponential family. 

Extending the previous corollary, we see that for a one-parameter exponential family, Eq. \eqref{eq:corollary-1} obtains equality. Under certain regularity conditions (see <a href="/stats-ml/likelihood-theory">my likelihood post</a>), $-\mathbb{E}\left[\frac{\partial^2 \ell(z; \mu)}{\partial \mu^2} \right]$ is the Fisher information, which describes the amount of information about $\mu$ that is held in $z$.

Since equality is obtained, we can also think of $-\mathbb{E}\left[\frac{\partial^2 \ell_q(z; \mu)}{\partial \mu^2} \right]$ as describing the amount of information about $\mu$ that is held in $z$. In addition, the difference between the former and the latter can be thought of as the amount of information gained by knowing, specifically, the distribution of $z$. 

---

## Estimation
Let $\mathbf{z} = (z_1, \dots, z_n)^\top$ and $\boldsymbol \mu = (\mu_1, \dots, \mu_n)^\top$. We'll denote the gradient of the (full) quasi-likelihood with respect to the parameters $\beta_1, \dots, \beta_m$ with:

$$
\begin{equation}
\label{eq:gradient-ql}
\mathbf{u} = \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta} = 
\begin{bmatrix}
\frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_1} \\
\vdots \\
\frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_m}
\end{bmatrix}
\end{equation}
$$

By <a href="#theorem-1">Theorem 1</a>, $\mathbf{u}$ satisfied:

$$
\begin{equation}
\label{eq:u-mean-var}
\begin{aligned}
\mathbb{E}[\mathbf{u}] &= \mathbf{0}_m \\
\text{Cov}(\mathbf{u}) &= - \mathbb{E}\left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta \partial \beta^\top} \right]
\end{aligned}
\end{equation}
$$

<details>
<summary>Proof.</summary>
By (1) in <a href="#theorem-1">Theorem 1</a>:
$$
\begin{aligned}
\mathbb{E}[\mathbf{u}]
&= \mathbb{E}\left[ \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta} \right] \\
&= \mathbb{E}\left[ 
    \begin{bmatrix}
    \sum_{i = 1}^n  \frac{\partial \ell_q(z_i; \mu_i)}{\partial \beta_1} \\
    \vdots \\
    \sum_{i = 1}^n  \frac{\partial \ell_q(z_i; \mu_i)}{\partial \beta_m} \\
    \end{bmatrix}
\right] \\
&= \sum_{i = 1}^n \mathbb{E}\left[ 
    \begin{bmatrix}
    \frac{\partial \ell_q(z_i; \mu_i)}{\partial \beta_1} \\
    \vdots \\
    \frac{\partial \ell_q(z_i; \mu_i)}{\partial \beta_m} \\
    \end{bmatrix}
\right] \\
&= \mathbf{0}_m
\end{aligned}
$$
And:
$$
\begin{aligned}
\text{Cov}\left(\frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i}, \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_j}\right)
&= \mathbb{E}\left[ \left( \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i} - \mathbb{E}\left[ \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i} \right] \right) \left(\frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_j}  - \mathbb{E}\left[ \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_j} \right] \right)\right] \\
&= \mathbb{E}\left[ \left( \frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i} \right) \left(\frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_j} \right)\right] \\
&= - \mathbb{E}\left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i \partial \beta_j} \right]
\end{aligned}
$$
</details>

The <i>maximum quasi-likelihood estimates</i> of $\beta$, denoted by $\hat{\beta}$, are found by setting $\mathbf{u}$ equal to $\mathbf{0}_m$ and solving for $\beta$, just like we would do for maximum likelihood estimation. 

<!-- #region theorem-3-wedderburn -->
<div id="theorem-3"></div>
<div class="theorem">
<strong>Theorem 3.<d-cite key=wedderburn1974></d-cite></strong>
{% tabs theorem-3-wedderburn %}
{% tab theorem-3-wedderburn theorem %}
The mean of $\hat{\beta}$ is approximately $\beta$, and the covariance of $\hat{\beta}$ is approximately:

$$
\text{Cov}(\hat{\beta}) \approx \text{Cov}^{-1}(\mathbf{u}) = \left[ -\mathbb{E}\left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta \partial \beta^\top} \right] \right]^{-1}
$$

if $\phi = 1$ (i.e. the mean-variance relationship is known entirely).
{% endtab %}
{% tab theorem-3-wedderburn proof %}
Let $\hat{\mathbf{u}}$ denote the gradient vector evaluated at the maximum quasi-likelihood estimate, $\hat{\beta}$. Since $\hat{\beta}$ is the value of $\beta$ such that $\mathbf{u}$ equals $\mathbf{0}$, a first-order Taylor approximation of $\mathbf{u}$ gives us:

$$
\begin{aligned}
\mathbf{u} &\approx \hat{\mathbf{u}} + \frac{\partial \mathbf{u}}{\partial \beta} (\beta - \hat{\beta}) \\
&= \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta \partial \beta^\top} (\beta - \hat{\beta}) \\
\implies
\beta - \hat{\beta} &\approx \left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta \partial \beta^\top} \right]^{-1} \mathbf{u}
\end{aligned}
$$

If we approximate the inverted matrix by its expectation, whose elements are given in <a href="#theorem-1">Theorem 1</a>, we get:

$$
\begin{aligned}
\beta - \hat{\beta} &\approx  -\text{Cov}^{-1}(\mathbf{u})\mathbf{u} \\
\implies
\hat{\beta} &\approx \beta + \text{Cov}^{-1}(\mathbf{u})\mathbf{u} 
\end{aligned}
$$

Since $\mathbf{u}$ has expectation zero, it is clear that $\mathbb{E}[\hat{\beta}] \approx \beta$.

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
<!-- #endregion -->

It can be shown that maximum quasi-likelihood estimates are asymptotically normal (I believe by an appeal to the weak law of large numbers similar to the argument used to show that maximum likelihood estimates are asymptotically normal). 

The above theorem holds only for a scale parameter $\phi = 1$. If we relax this assumption to $\phi > 0$, the expectation does not change, but we need to estimate $\phi$ before we can approximate the covariance of the maximum quasi-likelihood estimates. 

Wedderburn suggests the following estimate:

$$
\hat{\phi} = \frac{1}{n - m} \sum_{i = 1}^n \frac{(z_i - \hat{\mu}_i)^2}{V(\hat{\mu}_i)}
$$

where $\hat{\mu}_i$ is the mean evaluated at the maximum quasi-likelihood estimate of $\boldsymbol{\beta}$. 

### Weighted Least Squares
Suppose we wish to use Newton-Raphson to estimate $\boldsymbol \beta$. Newton-Raphson is an algorithm that can be used for optimization via a second-order approximation of the objective function. In its basic form, Newton-Raphson is used as a root-finding algorithm; however, if we let our objective function be a derivative (and everything behaves nicely enough), then finding the root of the derivative is equivalent to finding a local optimum. 

Our goal is to find:

$$
\hat{\boldsymbol{\beta}} = \underset{\boldsymbol \beta}{\arg \max} \left\{ \ell_q(\mathbf{z}; \boldsymbol{\mu}) \right\}
$$

which will be the solution to:

$$
\frac{\partial \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta}} = \mathbf{0}_m
$$

The Newton-Raphson updates are given by:

$$
\begin{equation}
\label{eq:n-r}
\hat{\boldsymbol{\beta}}_{t + 1}
= \hat{\boldsymbol{\beta}}_t - \left[ \left. \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta} \partial \boldsymbol{\beta}^\top} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t} \right]^{-1} \left(\left. \frac{\partial \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right)
\end{equation}
$$

Let $\mathbf{v}(\boldsymbol{\mu}) = \text{diag}(v(\boldsymbol{\mu}))$. Notice that:

$$
\begin{aligned}
\frac{\partial \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta}} &= \left(\frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}\right)^\top \mathbf{v}^{-1}(\boldsymbol{\mu})  \left(\mathbf{z} - \boldsymbol{\mu}\right) \\
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta} \partial \boldsymbol{\beta}^\top} \right] 
&= \left(\frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}\right)^\top \mathbf{v}^{-1}(\boldsymbol{\mu}) \left(\frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}\right)
\end{aligned}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\frac{\partial \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i} 
&= \sum_{j = 1}^n \frac{\partial}{\partial \beta_i} \left[  \ell_q(\mathbf{z}; \boldsymbol \mu) \right] \\
&= \sum_{j = 1}^n \frac{\partial \mu_j}{\partial \beta_i} \frac{\partial}{\partial \mu_j} \left[  \ell_q(\mathbf{z}; \boldsymbol \mu) \right] \\
&= \sum_{j  = 1}^n \frac{\partial \mu_j}{\partial \beta_i} \frac{z_j - \mu_j}{v(\mu_j)} \\
&= \left(\frac{\partial \boldsymbol \mu}{\partial \beta_i} \right)^\top \mathbf{v}^{-1}(\boldsymbol{\mu}) \left( \mathbf{z} - \boldsymbol \mu\right) \\
\implies
\frac{\partial \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta}} &= \left(\frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}\right)^\top \mathbf{v}^{-1}(\boldsymbol{\mu})  \left(\mathbf{z} - \boldsymbol{\mu}\right)
\end{aligned}
$$
By <a href="#theorem-1">Theorem 1</a>, we have:
$$
\begin{aligned}
- \mathbb{E}\left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol \mu)}{\partial \beta_i \beta_j} \right]
&= - \sum_{k = 1}^n \mathbb{E}\left[ \frac{\partial^2 \ell_q(z_k; \mu_k)}{\partial \beta_i \beta_j} \right]  \\
&= \sum_{k = 1}^n \frac{1}{v(\mu_k)} \frac{\partial \mu_k}{\partial \beta_i} \frac{\partial \mu_k}{\partial \beta_j} \\
&= \left(\frac{\partial \boldsymbol{\mu}}{\partial \beta_i}\right)^\top \mathbf{v}^{-1}(\boldsymbol{\mu}) \left(\frac{\partial \boldsymbol{\mu}}{\partial \beta_j} \right) \\
\implies - \mathbb{E}\left[ \frac{\partial^2 \ell_q(\mathbf{z}; \boldsymbol{\mu})}{\partial \boldsymbol{\beta} \partial \boldsymbol{\beta}^\top} \right] 
&= \left(\frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}\right)^\top \mathbf{v}^{-1}(\boldsymbol{\mu}) \left(\frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}\right)
\end{aligned}
$$
</details>

Define the <i>residuals</i> for our setting as:

$$
\mathbf{r}(\boldsymbol \beta) = \mathbf{z} - \boldsymbol{\mu} = \mathbf{z} - \mathbf{g}(\boldsymbol{\beta})
$$

Since the data are considered fixed in this setting, we have:

$$
\begin{aligned}
\frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} &= \frac{\partial \mathbf{z}}{\partial \boldsymbol{\beta}} - \frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}
= - \frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}}
\end{aligned}
$$

Let $$\hat{\boldsymbol{\mu}}_t = \boldsymbol{\mu} \rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}$$. We can then rewrite the Newton-Raphson updates as:

$$
\begin{aligned}
\hat{\boldsymbol \beta}_{t + 1}
&= \hat{\boldsymbol \beta}_t - \left[ \left(\left. \frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t} \right)^\top \mathbf{v}^{-1}(\hat{\boldsymbol{\mu}}_t) \left(\left. \frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right) \right]^{-1} \left( \left. \frac{\partial \boldsymbol{\mu}}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right)^\top \mathbf{v}^{-1}(\hat{\boldsymbol{\mu}}_t) (\mathbf{z} - \hat{\boldsymbol{\mu}}_t)  \\
&= \hat{\boldsymbol \beta}_t - \left[ \left(-\left. \frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t} \right)^\top \mathbf{v}^{-1}(\hat{\boldsymbol{\mu}}_t) \left(\left. -\frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right) \right]^{-1} \left( - \left. \frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right)^\top \mathbf{v}^{-1}(\hat{\boldsymbol{\mu}}_t) \mathbf{r}(\hat{\boldsymbol{\beta}}_t) \\
&= \hat{\boldsymbol \beta}_t - \left[ \left(\left. \frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t} \right)^\top \mathbf{v}^{-1}(\hat{\boldsymbol{\mu}}_t) \left(\left. \frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right) \right]^{-1} \left(\left. \frac{\partial \mathbf{r}(\boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right\rvert_{\boldsymbol{\beta} = \hat{\boldsymbol{\beta}}_t}\right)^\top \mathbf{v}^{-1}(\hat{\boldsymbol{\mu}}_t) (-\mathbf{r}(\hat{\boldsymbol{\beta}}_t))
\end{aligned}
$$

Notice that the last line in the above are exactly the <a href="https://en.wikipedia.org/wiki/Weighted_least_squares">weighted least squares</a> solution for regressing the (negative) residuals on the gradient (with respect to $\boldsymbol{\beta}$) of the residuals with weights given by $\frac{1}{v((\hat{\mu}_i)_t)}$ (the inverse of the variance function of observation $i$ computed at the current estimate $\hat{\boldsymbol{\beta}}_t$). 

The above equivalence shows that one can use the <a href="https://en.wikipedia.org/wiki/Gauss–Newton_algorithm">Gauss-Newton method</a> with an added weight matrix ($\mathbf{v}(\boldsymbol{\mu})$) to obtain maximum quasi-likelihood estimates of $\boldsymbol{\beta}$.



