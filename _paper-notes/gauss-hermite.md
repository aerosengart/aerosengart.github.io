---
layout: distill
title: A Note on Gauss-Hermite Quadrature
description:
date: 2026-03-14
tabs: true
tags: quadrature integration numerical
toc:
    - name: Background
    - name: Quadrature Procedure
      subsections:
        - name: Set-Up
        - name: Quadrature Details
        - name: Intuition
    - name: Connection To The Laplace Approximation
    - name: Example
      subsections:
        - name: Data Set-Up
        - name: An Approximation
        - name: Quadrature
bibliography: paper-notes.bib
---

In a <a href="/stats-ml/approx-comp/">recent post</a>, I compared a handful of different integral approximation methods. One of these, Gauss-Hermite quadrature, ended up not being very accurate in simulations, and I've decided to do a little bit more digging into the choice of knots in the approximation. In this post, I'll go through <i>A Note on Gauss-Hermite Quadrature</i> by Liu and Pierce.<d-cite key=liu1994></d-cite>

---

## Background
I've already done a basic overview of <a href="/stats-ml/quadrature/">quadrature</a> in general, but I'll reiterate some of the main points here for convenience. 

Gauss-Hermite quadrature is a method of numerical integration using the <a href="https://en.wikipedia.org/wiki/Hermite_polynomials">physicists' Hermite polynomials</a>. It is used to approximate integrals of the following form as a weighted sum:

$$
\begin{equation}
\label{eq:g-h}
\int_{-\infty}^{\infty} f(x) \exp(-x^2) dx \approx \sum_{i = 1}^m \omega(z_i) f(z_i)
\end{equation}
$$

where the $z_i$ are the roots of the $m$-th order Hermite polynomial, and $\omega(z_i)$ are associated weights. Both the roots and weights are deterministic, and the roots are symmetric about $0$, so the function of interest may end up being evaluated at "uninteresting" points in its domain. In Lin and Pierce's work, they describe a way to sidestep this issue and effectively use Gauss-Hermite quadrature to approximate indefinite integrals of functions, $g(t) > 0$, that satisfy some light conditions. 

---

## Quadrature Procedure
More specifically, for a function, $g(t)$, that is positive, unimodal, and that the ratio of $g(t)$ to a Gaussian curve is a relatively smooth function, our goal is to approximate an integral of the form:

$$
\int_{-\infty}^\infty g(t) dt
$$ 

such that the nodes fall in the region where $g(t)$ is "interesting". 


### Set-Up
Let $\phi(t; \mu, \sigma)$ be some Gaussian density. The integral form in Eq. \eqref{eq:g-h} can be rewritten, via suitable scaling and centering, as:

$$
\begin{equation}
\label{eq:g-h-rewrite}
\int_{-\infty}^\infty f(t) \phi(t; \mu, \sigma) dt = \sum_{i = 1}^m \omega(z'_i) f(z'_i)
\end{equation}
$$

The quadrature nodes then can be chosen as $z'_i = \mu + \sqrt{2 \sigma^2} z_i,$, and the weights as $\omega(z'_i) = \frac{\omega(z_i)}{\sqrt{\pi}}$. 

Suppose we select $\hat{\mu}$ be the mode of $g(t)$ and let:

$$
\begin{equation}
\label{eq:hat-sigma}
\begin{aligned}
\hat{\sigma} = \frac{1}{\sqrt{- \left.\frac{\partial^2 \log(g(t))}{\partial t^2} \right\rvert_{t = \hat{\mu}} }}
\implies
\hat{\sigma}^2 = \left[ - \left.\frac{\partial^2 \log(g(t))}{\partial t^2} \right\rvert_{t = \hat{\mu}} \right]^{-1}
\end{aligned}
\end{equation}
$$

These choices of parameter values make it so that, at the mode of $g(t)$, the first and second derivatives of $\phi(t; \hat{\mu}, \hat{\sigma})$ and $\log(g(t))$ are identical. 

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\left. \frac{\partial \log(g(t))}{\partial t} \right \rvert_{t = \hat{\mu}} &= 0 \\
\left. \frac{\partial}{\partial t} \left[ \log (\phi(t; \hat{\mu}, \hat{\sigma})) \right] \right \rvert_{t = \hat{\mu}}
&= \left. \frac{\partial}{\partial t} \left[ -\frac{1}{2} \log(2 \pi \hat{\sigma}^2) - \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right] \right\rvert_{t = \hat{\mu}}
= - \frac{\hat{\mu} - \hat{\mu}}{\hat{\sigma}^2} = 0 \\
\left. \frac{\partial^2}{\partial t^2} \left[ \log (\phi(t; \hat{\mu}, \hat{\sigma})) \right] \right \rvert_{t = \hat{\mu}}
&= - \frac{1}{\hat{\sigma}^2} 
= \left.  \frac{\partial^2 \log( g(t) ) }{\partial t^2} \right\rvert_{t = \hat{\mu}}
\end{aligned}
$$
</details>

We can then define a function, $h(t)$, and write:

$$
\begin{equation}
\label{eq:h-t}
\begin{aligned}
h(t) &= \frac{g(t)}{\phi(t; \hat{\mu}, \hat{\sigma})} \\
\implies
\int_{-\infty}^\infty g(t) dt &= \int_{- \infty}^\infty h(t) \phi(t; \hat{\mu}, \hat{\sigma}) dt
\end{aligned}
\end{equation}
$$

We then apply Gauss-Hermite quadrature to Eq. \eqref{eq:h-t}:

$$
\begin{equation}
\label{eq:g-h-3}
\begin{aligned}
\int_{-\infty}^\infty g(t) dt 
&= \int_{-\infty}^\infty h(t) \phi(t; \hat{\mu}, \hat{\sigma}) dt \\
&\approx \sum_{i = 1}^m \omega(z'_i) g(z'_i) \\
&= \sum_{i = 1}^m \frac{\omega(z_i)}{\sqrt{\pi}} h(\hat{\mu} + z_i \sqrt{2 \hat{\sigma}^2}) \\
&= \sum_{i = 1}^m \frac{\omega(z_i)}{\sqrt{\pi}} \frac{g(\hat{\mu} + z_i \sqrt{2 \hat{\sigma}^2})}{\frac{1}{\sqrt{2 \pi \hat{\sigma}^2}}\exp\left(- \frac{(\hat{\mu} + z_i \sqrt{2 \hat{\sigma}^2} - \hat{\mu})^2}{2 \hat{\sigma}^2}\right)} \\
&= \sqrt{2 \hat{\sigma}^2} \sum_{i = 1}^m \omega(z_i)\exp\left( z_i^2 \right) g(\hat{\mu} + z_i \sqrt{2 \hat{\sigma}^2})  \\
\end{aligned}
\end{equation}
$$

### Quadrature Details
First, note that:

$$
h(\hat{\mu}) = \frac{g(\hat{\mu})}{\phi(\hat{\mu}; \hat{\mu}, \hat{\sigma})} = \frac{g(\hat{\mu})}{\frac{1}{\sqrt{2 \pi \hat{\sigma}^2}} \exp\left(- \frac{(\hat{\mu} - \hat{\mu})^2}{2 \hat{\sigma}^2}\right)} = g(\hat{\mu})\sqrt{2 \pi \hat{\sigma}^2}
$$


Define:

$$
\begin{equation}
\label{eq:c-k}
c_k = \frac{1}{k!} \left. \left[ \frac{\partial^{(k)}}{\partial t^{(k)}}\left[ \exp\left(\log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \right] \right] \right\rvert_{t = \hat{\mu}}
\end{equation}
$$

A Taylor expansion of $h(t)$ about $t = \hat{\mu}$ is given by:

$$
\begin{equation}
\label{eq:h-t-exp}
\begin{aligned}
h(t)
&= h(\hat{\mu}) \left(1 + \sum_{k = 1}^\infty c_k (t-\hat{\mu})^k \right) 
\end{aligned}
\end{equation}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
h(t) 
&= \sum_{k = 0}^\infty \left[ \frac{\partial^{(k)} h(t)}{\partial t^{(k)}}\right] \frac{(t-\hat{\mu})^k}{k!} \\
&= h(\hat{\mu}) + \sum_{k = 1}^\infty \left. \left[ \frac{\partial^{(k)} h(t)}{\partial t^{(k)}}\right] \right\rvert_{t = \hat{\mu}} \frac{(t-\hat{\mu})^k}{k!}  \\
&= h(\hat{\mu}) \left(1 + \frac{1}{h(\hat{\mu})} \sum_{k = 1}^\infty \left. \left[ \frac{\partial^{(k)} h(t)}{\partial t^{(k)}}\right] \right\rvert_{t = \hat{\mu}}\frac{(t-\hat{\mu})^k}{k!}  \right)  \\
&= h(\hat{\mu}) \left(1 + \sum_{k = 0}^\infty \left. \left[ \frac{\partial^{(k)}}{\partial t^{(k)}}\left[ \frac{h(t)}{h(\hat{\mu})}\right] \right] \right\rvert_{t = \hat{\mu}} \frac{(t-\hat{\mu})^k}{k!} \right) \\
&= h(\hat{\mu}) \left(1 + \sum_{k = 0}^\infty \left. \left[ \frac{\partial^{(k)}}{\partial t^{(k)}}\left[ \frac{g(t)}{\phi(t; \hat{\mu}, \hat{\sigma}) g(\hat{\mu}) \sqrt{2 \pi \hat{\sigma}^2}}\right] \right] \right\rvert_{t = \hat{\mu}} \frac{(t-\hat{\mu})^k}{k!} \right) \\
&= h(\hat{\mu}) \left(1 + \sum_{k = 0}^\infty \left. \left[ \frac{\partial^{(k)}}{\partial t^{(k)}}\left[ \frac{g(t)}{g(\hat{\mu})}  \exp\left(\frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \right] \right] \right\rvert_{t = \hat{\mu}} \frac{(t-\hat{\mu})^k}{k!} \right) \\
&= h(\hat{\mu}) \left(1 + \sum_{k = 0}^\infty \left. \left[ \frac{\partial^{(k)}}{\partial t^{(k)}}\left[ \exp\left(\log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \right] \right] \right\rvert_{t = \hat{\mu}} \frac{(t-\hat{\mu})^k}{k!} \right)
\end{aligned}
$$
</details>

It can be shown that, by the choice of $\hat{\mu}$ and $\hat{\sigma}^2$, $c_1 = c_2 = 0$.

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
c_1 &= \frac{1}{1} \left. \left[ \frac{\partial}{\partial t}\left[ \exp\left(\log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \right] \right] \right\rvert_{t = \hat{\mu}} \\
&= \exp\left(\log(g(\hat{\mu})) - \log(g(\hat{\mu})) + \frac{(\hat{\mu} - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \left( \left. \left[ \frac{\partial \log(g(t))}{\partial t} \right]\right\rvert_{t = \hat{\mu}} + \exp\left(\frac{(\hat{\mu} - \hat{\mu})^2}{2 \hat{\sigma}^2}\right)\frac{\partial}{\partial t} \left[ \frac{(t-\hat{\mu})^2}{2 \hat{\sigma}^2} \right]\right) \\
&= \exp(0) \left(0 + \exp(0) \frac{\hat{\mu} - \hat{\mu}}{\hat{\sigma}^2}\right)  & \left(\hat{\mu} \text{ is mode}\right) \\
&= 0  \\
c_2 &= \frac{1}{2} \left. \left[ \frac{\partial^2}{\partial t^2}\left[ \exp\left(\log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \right] \right] \right\rvert_{t = \hat{\mu}} \\
&= \frac{1}{2} \left. \frac{\partial}{\partial t} \left[ \exp\left(\log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \frac{\partial}{\partial t}\left[ \log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2\hat{\sigma}^2}\right] \right] \right\rvert_{t = \hat{\mu}} \\
&= \frac{1}{2} \left. \frac{\partial}{\partial t} \left[ \exp\left(\log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \right] \right\rvert_{t = \hat{\mu}} \left. \left[ \frac{\partial}{\partial t} \left[ \log(g(t)) - \log(g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2\hat{\sigma}^2}\right] \right]  \right\rvert_{t = \hat{\mu}} +
\frac{1}{2} \exp\left(\log(g(\hat{\mu})) - \log(g(\hat{\mu})) + \frac{(\hat{\mu} - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \left. \frac{\partial^2}{\partial t^2} \left[ \log(g(t)) - \log(g(\hat{\mu})) + \frac{(\hat{\mu} - \hat{\mu})^2}{2 \hat{\sigma}^2} \right] \right\rvert_{t = \hat{\mu}} \\
&= \frac{c_1}{2} \left( \left. \left[ \frac{ \partial \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}}  + \frac{\hat{\mu} - \hat{\mu}}{\hat{\sigma}^2}\right) + \frac{\exp(0)}{2} \left( \left. \left[ \frac{\partial^2 \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} + \frac{1}{\hat{\sigma}^2} \right) \\
&= 0 + \frac{1}{2} \left. \left[ \frac{\partial^2 \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} + \frac{1}{2\hat{\sigma}^2}  & \left(c_1 = 0 \right) \\
&= \frac{1}{2}\left. \left[ \frac{\partial^2 \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} + \frac{1}{\left[- 2\left. \frac{\partial^2 \log(g(t))}{\partial t^2} \right\rvert_{t = \hat{\mu}} \right]^{-1}} & \left(\text{Eq. } \eqref{eq:hat-sigma} \right)\\
&= \frac{1}{2} \left. \left[ \frac{\partial^2 \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} - \frac{1}{2}\left. \left[ \frac{\partial^2 \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} \\
&= 0 
\end{aligned}
$$
</details>

Now, we place $h(t)$ with the Taylor expansion: 

$$
\begin{equation}
\label{eq:g-h-4}
\begin{aligned}
\int_{-\infty}^\infty g(t) dt 
&= \int_{-\infty}^\infty h(t) \phi(t; \hat{\mu}, \hat{\sigma}) dt  \\
&= \int_{-\infty}^\infty  \phi(t; \hat{\mu}, \hat{\sigma})  h(\hat{\mu})\left(1 + \sum_{k = 1}^\infty c_k(t - \hat{\mu})^k \right) dt  \\
&= \mathbb{E}_{t \sim \phi(\hat{\mu}, \hat{\sigma})}\left[ h(\hat{\mu}) \left(1 + \sum_{k = 1}^\infty c_k(t - \hat{\mu})^k\right) \right] \\
&= h(\hat{\mu}) \left(1 + \sum_{k = 1}^\infty c_k  \mathbb{E}_{t \sim \phi(\hat{\mu}, \hat{\sigma})}\left[  (t - \hat{\mu})^k \right]  \right)\\
\end{aligned}
\end{equation}
$$

In the last line above, the terms involving odd values of $k$ disappear because $\phi(\hat{\mu}, \hat{\sigma})$ is symmetric (so the odd central moments are zero). 

Thus, if it were the case that all of the coefficients beyond $c_{2(m+1)}$ were $0$, the approximation would be exact. In other words, $m$-order Gaussian quadrature will perform well if $h(t)$ can be approximated well by a polynomial of order $2m + 1$ in the domain where $g(t)$ is large. 

### Intuition
If we take a Bayesian perspective, we can think of $g(\cdot)$ as a joint distribution and the integral we are interested in approximating will yield a marginal distribution. In this case, $\phi(\cdot; \cdot, \cdot)$ serves as the prior over the latent variable we are integrating over, and $h(\cdot)$ is the conditional distribution of the latent variable and the variable of interest. Since the posterior is proportional to the joint (by Bayes' Theorem), the choices of $\hat{\mu}$ and $\hat{\sigma}$ are the posterior mode and standard deviation of a second-order Taylor approximation to the posterior. 


---

## Connection To The Laplace Approximation
Consider the Laplace approximation to:

$$
\int_{-\infty}^\infty g(t) dt = \int_{-\infty}^\infty \exp(\log(g(t))) dt
$$

Again, we let $\hat{\mu}$ be the mode of $g(t)$ and $\hat{\sigma}$ be as defined in Eq. \eqref{eq:hat-sigma}. A second-order Taylor expansion of $\log(g(t))$ about $t = \hat{\mu}$ yields the approximation:

$$
\int_{-\infty}^\infty g(t) dt
\approx
g(\hat{\mu}) \sqrt{2 \pi \hat{\sigma}^2}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\int_{-\infty}^\infty g(t) dt
&\approx \int_{-\infty}^\infty \exp\left(\log(g(\hat{\mu})) + \left. \left[ \frac{\partial \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} (t - \hat{\mu})  + \frac{1}{2}   \left. \left[ \frac{\partial \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} (t - \hat{\mu})^2 \right) dt \\
&= \exp(\log(g(\hat{\mu}))) \int_{-\infty}^\infty \exp\left(\frac{1}{2}   \left. \left[ \frac{\partial \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}} (t - \hat{\mu})^2 \right) dt & \left(\hat{\mu} \text{ is mode}\right)\\
&= g(\hat{\mu}) \int_{-\infty}^\infty \exp\left( - \frac{(t - \hat{\mu})^2}{2 \left[ - \left. \left[ \frac{\partial \log(g(t))}{\partial t} \right] \right\rvert_{t = \hat{\mu}}\right]^{-1}} \right) dt \\
&= g(\hat{\mu}) \int_{-\infty}^\infty \exp\left( - \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2} \right) dt & \left(\text{Eq. } \eqref{eq:hat-sigma} \right)\\
&= g(\hat{\mu}) \sqrt{2 \pi \hat{\sigma}^2} \int_{-\infty}^\infty \frac{1}{\sqrt{2 \pi \hat{\sigma}^2}} \exp\left( - \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2} \right) dt \\
&= g(\hat{\mu}) \sqrt{2 \pi \hat{\sigma}^2}
\end{aligned}
$$
</details>

Now, suppose we perform the approximation specified in Eq. \eqref{eq:g-h-3} but with only one node. The Hermine polynomial, $H_n(z)$, of order $n$ has weights given by:

$$
\omega(z_i) = \frac{2^{n-1} n! \sqrt{\pi}}{n^2 \left[ H_{n-1}(z_i)\right]^2}
$$

A single node corresponds to $H_1(z) = 2z$, implying that $z_1 = 0$. Recall that $H_0(z)  = 1$. Thus, the corresponding weight is given by:

$$
\omega(z_1) = \frac{2^{1-1} 1! \sqrt{\pi}}{1^2 \left[ H_0(z_1) \right]^2} = \sqrt{\pi}
$$

The approximation can then be written as:

$$
\begin{aligned}
\int_{\infty}^\infty g(t) dt 
&\approx \sqrt{2 \hat{\sigma}^2} \omega(z_1) \exp(z_1^2) g(\hat{\mu} + z_1 \sqrt{2 \hat{\sigma}^2}) \\
&= \sqrt{2 \pi \hat{\sigma}^2}  g(\hat{\mu})
\end{aligned}
$$

which is exactly the same as the Laplace approximation!


---

## Example
I'll do a brief example related to my score test project. 

### Data Set-Up
We assume that to have a dataset consisting of $N$ observations of: scalar response, $y$; a $p$-dimensional vector of covariates associated with fixed effects, $\mathbf{x}$; a scalar-valued covariate associated with random effects, $z$. We assume the observations are divided into $K$ independent clusters with $n_k$ observations per cluster such that $\sum_{k = 1}^K n_k = N$.

We assume that the responses follow a negative binomial distribution with conditonal means and variances given by:

$$
\begin{aligned}
\mathbb{E}[y_{k,i} \rvert \boldsymbol{\beta}] &= \mu_{k,i} = \exp(\eta_{k,i}) = \mathbf{x}_{k,i}^\top \boldsymbol{\alpha} + \mathbf{z}_{k,i}^\top \boldsymbol{\beta}_k \\
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
\log(p(\mathbf{y} \rvert \boldsymbol{\alpha}, \gamma, \beta_k))
= \sum_{i = 1}^{n_k} \left[ \log\left( \frac{\Gamma(y_{k,i} + \gamma + 2)}{y_{k,i}! \Gamma(\gamma + 2)}\right) + y_{k,i} \left[ \log(\mu_{k,i}) - \log(\mu_{k,i} + \gamma) \right] + \gamma\left[ \log(\gamma) - \log(\mu_{k,i} + \gamma)\right] \right]
$$

where $\Gamma(\cdot)$ denotes the Gamma function. We have placed a Gaussian prior on the random effects:

$$
\log(p(\boldsymbol{\beta};\tau^2)) = -\frac{1}{2} \sum_{k = 1}^K \left[\log(2 \pi \tau^2) + \frac{\beta_k^2}{\tau^2} \right]
$$

Our goal will be to approximate the marginal likelihood, which can be obtained by integrating over the random effects with respect to their prior distribution:

$$
p(\mathbf{y}; \boldsymbol{\alpha}, \gamma, \tau^2)
= \int p(\mathbf{y}; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \boldsymbol{\beta}) p(\boldsymbol{\beta}; \tau^2) d \boldsymbol{\beta}
$$


### An Approximation
We define:

$$
q(\beta_k) := \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \beta_k)) - \frac{\beta_k^2}{2 \tau^2} \\
$$

We will do a second-order Laplace approximation of each $q(\beta_k)$ about the posterior mode of $\beta_k$, which we denote with $\hat{\beta}_k$. Derivations of the first- and second-order derivatives of $q(\beta_k)$ are deferred to the end of this post (see <a href="#derivations">here</a>).

$$
q(\beta_k) 
\approx q(\hat{\beta}_k) + (\beta_k - \hat{\beta}_k) q'(\hat{\beta}_k) + \frac{1}{2}(\beta_k - \hat{\beta}_k)^2 q''(\hat{\beta}_k)
$$

Recall that the posterior is proportional to the joint distribution (see <a href="https://en.wikipedia.org/wiki/Bayes%27_theorem">Bayes Theorem</a>), and the constant of proportionality (for fixed data) is just the marginal probability which does not depend on the random effects. Thus, the value of $\boldsymbol{\beta}$ that maximize the posterior (i.e. the posterior mode) will also be the maximizing value of the joint probability. That is, $q'(\hat{\beta}_k) = 0$, so we can write:

$$
\begin{aligned}
p(\mathbf{y}_k, \beta_k; \boldsymbol{\alpha}, \gamma, \tau^2)
&= \exp(\log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \beta_k))) (2 \pi \tau^2)^{-\frac{1}{2}} \exp\left(- \frac{\beta_k^2}{2 \tau^2} \right) \\
&= \underbrace{(2 \pi \tau^2)^{-\frac{1}{2}}}_{=: C} \exp\left( \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \beta_k)) - \frac{\beta_k^2 }{2 \tau^2} \right) \\
&= C \exp\left( \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \beta_k)) - \frac{\beta_k^2}{2\tau^2} \right) \\
&= C \exp\left( \sum_{k = 1}^K q(\beta_k) \right) \\
&\approx C \exp\left( \sum_{k = 1}^K \left[ q(\hat{\beta}_k) + \frac{1}{2}(\beta_k - \hat{\beta}_k)^2 q''(\hat{\beta}_k) \right] \right) \\
&= C \exp\left( q(\hat{\beta}_k)\right) \exp\left(\frac{-(\beta_k - \hat{\beta}_k)^2}{2 (-q''(\hat{\beta}_k))^{-1}} \right)
\end{aligned}
$$

Notice that the last line is proportional to a Gaussian density with mean $\hat{\beta}_k$ and variance $-\frac{1}{q''(\hat{\beta}_k)}$. 


### Quadrature
Since quadrature is for one-dimensional integrals, we need to reformulate the goal. I believe that, since the random effects are independent, and the observations within a cluster are also independent, the clustered observations should also be marginally independent. This allows us to write:

$$
\begin{aligned}
p(\mathbf{y}; \boldsymbol{\alpha}, \gamma, \tau^2)
&= \prod_{k = 1}^K p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2) \\
&= \prod_{k = 1}^K \int p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \beta_k) p(\beta_k; \tau^2) d \beta_k 
\end{aligned}
$$

We define:

$$
\begin{aligned}
&\hat{\mu}_k = \hat{\boldsymbol{\beta}}_k;& 
&\hat{\sigma}_k^2 = \left[ -q''(\hat{\beta}_k) \right]^{-1};& 
&h(\beta_k) = \exp\left( q(\hat{\beta}_k) \right);&
&g(\beta_k) = h(\beta_k) \phi(\beta_k; \hat{\mu}_k, \hat{\sigma}_k)
\end{aligned}
$$

To get:

$$
\begin{aligned}
p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\alpha}, \gamma, \tau^2)
&= \prod_{k = 1}^K \int p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma, \tau^2 \rvert \beta_k) p(\beta_k; \tau^2) d \beta_k \\
&\approx \prod_{k = 1}^K \int C \exp \left(q(\hat{\beta}_k)\right) \exp\left(\frac{-(\beta_k - \hat{\beta}_k)^2}{2 (-q''(\hat{\beta}_k))^{-1}}\right) d \beta_k  \\
&= \prod_{k = 1}^K \int (2 \pi \tau^2)^{-\frac{1}{2}} \frac{h(\hat{\beta}_k)}{ (2 \pi \hat{\sigma}_k^2)^{\frac{1}{2}} \phi(\hat{\beta}_k; \hat{\mu}_k, \hat{\sigma}_k)} d \beta_k \\
&= \prod_{k = 1}^K (2 \pi \tau^2)^{-\frac{1}{2}} (2 \pi \hat{\sigma}_k^2)^{-\frac{1}{2}} \underbrace{\int h(\hat{\beta}_k) \phi(\hat{\beta}_k; \hat{\mu}_k, \hat{\sigma}_k) d \beta_k}_{(\star)} 
\end{aligned}
$$

We notice now that $(\star)$ is the form used in Liu and Pierce's adaptive GH quadrature procedure. Their approximation is:

$$
\begin{aligned}
\int h(\hat{\beta}_k) \phi(\hat{\beta}_k; \hat{\mu}, \hat{\sigma}) d \beta_k
&\approx \sqrt{2 \hat{\sigma}^2_k} \sum_{i = 1}^m w_i \exp(x_i^2) g(\hat{\mu}_k + x_i \sqrt{2 \hat{\sigma}^2_k})
\end{aligned}
$$

where $x_i$ is the $i$-th root of the $m$-order Hermite polynomial, and $w_i$ is the corresponding weight.


### Derivations
We have the following intermediate terms:

$$
\begin{aligned}
\frac{d \eta_{k,i}}{d \beta_k} 
&= \frac{d}{d \beta_k} \left[ \mathbf{x}_{k,i}^\top \boldsymbol{\alpha} + \beta_k z_{k,i} \right] =
 z_{k,i} \\
\frac{d \mu_{k,i}}{d \beta_k} 
&= \frac{d}{d \beta_k} \left[ \exp(\eta_{k,i}) \right] 
= \exp(\eta_{k,i}) z_{k,i}
= \mu_{k,i} z_{k,i} \\
\frac{d}{d \beta_k} \left[ \frac{\mu_{k,i}}{\mu_{k,i} + \gamma}\right] 
&= \frac{(\mu_{k,i} + \gamma) \frac{d \mu_{k,i}}{d \beta_k} - \mu_{k,i}\frac{d}{d \beta_k} \left[ \mu_{k,i} + \gamma \right]}{(\mu_{k,i} + \gamma)^2} 
= \frac{(\mu_{k,i} + \gamma) \mu_{k,i} z_{k,i} - \mu_{k,i}^2 z_{k,i}}{(\mu_{k,i} + \gamma)^2}
= \frac{\gamma \mu_{k,i} z_{k,i}}{(\mu_{k,i} + \gamma)^2}
\end{aligned}
$$

The first derivative is then:

$$
\begin{aligned}
\frac{d q(\beta_k)}{d \beta_k}
&= \left( \sum_{i = 1}^{n_k} \frac{d}{d\beta_k} \left[  \log\left( \frac{\Gamma(y_{k,i} + \gamma + 2)}{y_{k,i}! \Gamma(\gamma + 2)}\right) + y_{k,i} \left[ \log(\mu_{k,i}) - \log(\mu_{k,i} + \gamma) \right] + \gamma\left[ \log(\gamma) - \log(\mu_{k,i} + \gamma)\right]\right] \right) - \frac{\beta_k}{\tau^2}\\
&= \left( \sum_{i = 1}^{n_k} \left[  y_{k,i} \left(z_{k,i} - \frac{\frac{d \mu_{k,i}}{d \beta_k}}{\mu_{k,i} + \gamma} \right) + \gamma \left( -\frac{\frac{d \mu_{k,i}}{d \beta_k}}{\mu_{k,i} + \gamma} \right) \right]\right)  - \frac{\beta_k}{\tau^2} \\ 
&= \left( \sum_{i = 1}^{n_k} \left[  y_{k,i} \left(z_{k,i} - \frac{\mu_{k,i} z_{k,i}}{\mu_{k,i} + \gamma} \right) + \gamma \left( -\frac{\mu_{k,i} z_{k,i}}{\mu_{k,i} + \gamma} \right) - \frac{\beta_k}{\tau^2}\right] \right) - \frac{\beta_k}{\tau^2}\\
&=\left( \sum_{i = 1}^{n_k} z_{k,i} \left[ y_{k,i} - \frac{(1 + \gamma) \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \right)- \frac{\beta_k}{\tau^2}
\end{aligned}
$$

And the second derivative is:

$$
\begin{aligned}
\frac{d^2 q(\beta_k)}{d \beta_k^2}
&= \left( \sum_{i = 1}^{n_k} z_{k,i} \frac{d}{d \beta_k} \left[ y_{k,i} - \frac{(1 + \gamma) \mu_{k,i}}{\mu_{k,i} + \gamma} \right]  \right) - \frac{1}{\tau^2}\\ 
&= \left( \sum_{i = 1}^{n_k} z_{k,i} \left[ - \frac{(1 + \gamma)\mu_{k,i} z_{k,i}}{(\mu_{k,i} + \gamma)^2} \right]  \right) - \frac{1}{\tau^2} \\
&= - \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \right]
\end{aligned}
$$



<!-- ## Asymptotics
Liu and Pierce also include an asymptotic analysis of their proposed variation on Gauss-Hermite quadrature. The consider approximating integrals of the form:

$$
\int_{-\infty}^\infty \exp\left(n l(t) \right) dt
$$

for some unimodal function, $l(t)$, with $m$-order Gauss-Hermite quadrature. In this case, we can treat $g(t) = \exp(n l(t))$ and plug-in to the expressions we found before. 

As shown <a href="#quadrature-details">above</a>, we need only consider the terms after the first $2m + 1$ in Eq. \eqref{eq:g-h-4} because of our Taylor expansion of $h(t)$. It thus suffices to consider the integral of just the $2(m + 1)$-th term, since the coefficients in the Taylor expansion decay. That is, the error will be of the same order as this integral.

Define $\hat{\mu}$ as the mode of $l(t)$, and let $\hat{\sigma} = \left[ -\left. \left[ \frac{d^2 l(t)}{d t^2} \right] \right\rvert_{t = \hat{\mu}}\right]^{-\frac{1}{2}}$. Define, too, $\pi(t)$ to satisfy $h(t) = h(\hat{\mu}) \pi(t)$. It follows from Eq. \eqref{eq:h-t-exp} and Eq. \eqref{eq:c-k} that:

$$
\begin{aligned}
\pi(t) &= 1 + \sum_{k = 1}^{\infty} c_k(t - \hat{\mu})^k \\
&= 1 + \sum_{k = 1}^{\infty} \left. \left[ \frac{\partial^{(k)}}{\partial t^{(k)}} \left[ \exp\left(\log(n q(t)) - \log(n g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2} \right) \right] \right] \right\rvert_{t = \hat{\mu}} \frac{(t - \hat{\mu})^k}{k!} \\
&= \exp\left(\log(n g(t)) - \log (n g(\hat{\mu})) + \frac{n(t - \hat{\mu})^2}{2 \hat{\sigma}^2}\right) \\
\implies c_{2(m + 1)} &= \frac{1}{(2(m + 1))!}  \left. \left[ \frac{\partial^{(2(m+1))}}{\partial t^{(2(m+1))}} \left[ \exp\left(\log(n g(t)) - \log(n g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2} \right) \right] \right] \right\rvert_{t = \hat{\mu}}
\end{aligned}
$$

It is not shown in the paper, but it is possible to prove:

$$
\left. \left[ \frac{\partial^{(2(m+1))}}{\partial t^{(2(m+1))}} \left[ \exp\left(\log(n g(t)) - \log(n g(\hat{\mu})) + \frac{(t - \hat{\mu})^2}{2 \hat{\sigma}^2} \right) \right] \right] \right\rvert_{t = \hat{\mu}} = O\left( n^{\lfloor \frac{2(m + 1)}{3} \rfloor } \right) 
$$

Furthermore, it is <a href="https://en.wikipedia.org/wiki/Normal_distribution#Moments">known</a> that the $k$-th central moment (for even $k$) for a Gaussian distribution is proportional to the $k$-th power of the standard deviation. Thus:

$$
\begin{aligned}
\int_{-\infty}^\infty \exp\left( \phi(t; \hat{\mu}, \hat{\sigma}) h(\hat{\mu}) c_{2(m + 1)} (t - \hat{\mu})^{2(m+1)}  \right) dt 
&= h(\hat{\mu}) \mathbb{E}_{t \sim \phi(\hat{\mu}, \hat{\sigma})} \left[ c_{2(m + 1)} (t - \hat{\mu})^{2(m + 1)} \right] \\
&\propto c_{2(m+1)} \mathbb{E}_{t \sim \phi(\hat{\mu}, \hat{\sigma})} \left[ (t - \hat{\mu})^{2(m + 1)} \right] \\
&\propto c_{2(m+1)} \hat{\sigma}^{m+1} \\
&= c_{2(m+1)} 
\end{aligned}
$$
 -->
