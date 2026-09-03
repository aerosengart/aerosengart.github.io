---
layout: distill
title: Expectation-Maximization
description:
date: 2026-04-17
tabs: true
tags: numerical optimization
toc:
    - name: General Algorithm
      subsections:
        - name: Intuition
        - name: Algorithm
    - name: Extensions
    - name: Example
      subsections:
        - name: Data Set-Up
bibliography: stats-ml.bib
---

In this post, I will briefly cover what goes into the expectation-maximization algorithm and work through an example relevant to my score test project. 

---

## General Algorithm
The expectation-maximization (EM) algorithm is used to obtain maximum likelihood or maximum a posterior estimates of parameters. We assume to have observed data, $\mathbf{y}$, and some latent variables/missing data, $\mathbf{z}$, along with a model parametrized by an unknown vector, $\boldsymbol{\theta}$. Perhaps there are some nuisance parameters and additional data, but perhaps not. It will not affect the algorithm (just the notation).


### Intuition 
Our goal is to find the value of $\boldsymbol{\theta}$ that maximizes the <i>marginal</i> log-likelihood:

$$
\ell(\boldsymbol{\theta}; \mathbf{y}) 
= \log \left[ \int p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta})d \mathbf{z}  \right]
$$

This can be very hard (or impossible) to evaluate, let alone maximize. However, let $q(\mathbf{z})$ be some arbitrary distribution for the latent variables. We can instead maximize the equivalent form:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{y}) 
&= \log\left[ \int \frac{p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta}) q(\mathbf{z})}{q(\mathbf{z})} d \mathbf{z}  \right] \\
&= \log\left[ \mathbb{E}_{\mathbf{z} \sim q( \cdot; \boldsymbol{\theta})} \left[ \frac{p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta})}{q(\mathbf{z})}  \right] \right] 
\end{aligned}
$$

An application of Jensen's inequality yields:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{y})
&= \log\left[ \mathbb{E}_{\mathbf{z} \sim q( \cdot; \boldsymbol{\theta})} \left[ \frac{p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta})}{q(\mathbf{z})}  \right] \right] \\
&\geq \mathbb{E}_{\mathbf{z} \sim q( \cdot; \boldsymbol{\theta})} \left[ \log\left[  \frac{p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta})}{q(\mathbf{z})} \right] \right] \\
&= \mathbb{E}_{\mathbf{z} \sim q( \cdot; \boldsymbol{\theta})} \left[ \log\left(p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta}) \right) \right] - \underbrace{\mathbb{E}_{\mathbf{z} \sim q( \cdot; \boldsymbol{\theta})} \left[ \log( q(\mathbf{z})) \right]}_{=: -\mathbb{H}(q)}  \\
&= \mathbf{L}(\boldsymbol{\theta}, q(\mathbf{z}); \mathbf{y})
\end{aligned}
$$

<aside><p>$\mathbb{H}(q)$ is the <strong>entropy</strong> of $q$.</p></aside>

We call $\mathbf{L}_q (\boldsymbol{\theta}; \mathbf{y})$ the <strong>ELBO</strong> (evidence lower bound) because it lower bounds the <i>evidence</i> (the marginal log-likelihood). We can manipulate the ELBO as follows:

$$
\begin{aligned}
\mathbf{L}(\boldsymbol{\theta}, q(\mathbf{z}); \mathbf{y})
&= \int \log\left[ \frac{p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta})}{q(\mathbf{z})}\right] q(\mathbf{z}) d \mathbf{z} \\
&= \int \log\left[ \frac{p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y}) p(\mathbf{y}; \boldsymbol{\theta})}{q(\mathbf{z})}\right] q(\mathbf{z}) d \mathbf{z} \\
&= \int \log\left[ \frac{p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y})}{q(\mathbf{z})}\right] q(\mathbf{z})  d \mathbf{z} + \int \log\left[ p(\mathbf{y}; \boldsymbol{\theta}) \right] q(\mathbf{z}) d \mathbf{z} \\
&= - \int \log\left[ \frac{q(\mathbf{z})}{p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y})}\right] q(\mathbf{z})  d \mathbf{z} + \log\left( p(\mathbf{y}; \boldsymbol{\theta}) \right)  \underbrace{\int q(\mathbf{z}) d \mathbf{z}}_{= 1} \\
&= -D_{KL}(q(\mathbf{z}) \rvert \rvert p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y})) + \log(p(\mathbf{y}; \boldsymbol{\theta}))
\end{aligned}
$$

where $D_{KL}(q(\mathbf{z}) \rvert \rvert p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y}))$ is the Kullback-Leibler divergence between $q(\mathbf{z})$ and $p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y})$. Since this divergence is always non-negative and will be zero if, and only if, $q(\mathbf{z}) = p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y})$, the ELBO can be maximized over choices of $q$ by choosing $q(\mathbf{z}) = p(\mathbf{z}; \boldsymbol{\theta} \rvert \mathbf{y})$. With this choice, the ELBO becomes equal to the marginal log-likelihood. 

Expectation-maximization algorithms are iterative procedures. Suppose we have some intermediate guess as to the value of $\boldsymbol{\theta}$ (denote this with $\hat{\boldsymbol{\theta}}^{(t)}$). If, at the next iteration, we update our guess to $\hat{\boldsymbol{\theta}}^{(t+1)}$ such that:

$$
\mathbf{L}(\hat{\boldsymbol{\theta}}^{(t)}, p(\mathbf{z}; \boldsymbol{\theta}^{(t)} \rvert \mathbf{y}); \mathbf{y}) \leq \mathbf{L}(\hat{\boldsymbol{\theta}}^{(t+1)}, p(\mathbf{z}; \boldsymbol{\theta}^{(t+1)} \rvert \mathbf{y}); \mathbf{y}) 
$$

Then we see that:

$$
\begin{aligned}
0 &\leq \mathbf{L}(\hat{\boldsymbol{\theta}}^{(t+1)}, p(\mathbf{z}; \hat{\boldsymbol{\theta}}^{(t+1)} \rvert \mathbf{y}); \mathbf{y}) - \mathbf{L}(\hat{\boldsymbol{\theta}}^{(t)}, p(\mathbf{z}; \boldsymbol{\theta}^{(t)} \rvert \mathbf{y}); \mathbf{y}) \\
\implies 0 &\leq 
\left[ \underbrace{-D_{KL}(p(\mathbf{z}; \hat{\boldsymbol{\theta}}^{(t+1)} \rvert \mathbf{y}) \rvert \rvert p(\mathbf{z}; \hat{\boldsymbol{\theta}}^{(t+1)} \rvert \mathbf{y}))}_{=0} + \log(p(\mathbf{y}; \hat{\boldsymbol{\theta}}^{(t+1)})) \right] - \left[ \underbrace{-D_{KL}(p(\mathbf{z}; \hat{\boldsymbol{\theta}}^{(t)} \rvert \mathbf{y}) \rvert \rvert p(\mathbf{z}; \hat{\boldsymbol{\theta}}^{(t)} \rvert \mathbf{y}))}_{=0} + \log(p(\mathbf{y}; \hat{\boldsymbol{\theta}}^{(t)}))  \right] \\
\implies 0 &\leq 
\log(p(\mathbf{y}; \hat{\boldsymbol{\theta}}^{(t+1)})) - \log(p(\mathbf{y}; \hat{\boldsymbol{\theta}}^{(t)})) 
\end{aligned}
$$

This shows that updating our parameter vector values such that the ELBO increases will also increase the marginal log-likelihood.

### Algorithm
The EM algorithm interates the following steps until convergence. Let $\boldsymbol{\theta}^{(t)}$ denote the value of $\boldsymbol{\theta}$ at iteration $t$. 

#### Expectation Step
The algorithm begins with computing the ELBO. 

$$
\begin{aligned}
\mathbf{L}(\boldsymbol{\theta}, q(\mathbf{z}); \mathbf{y}) 
&= \int \log\left[ \frac{p(\mathbf{y}, \mathbf{z}; \boldsymbol{\theta})}{q(\mathbf{z})}\right] q(\mathbf{z}) d \mathbf{z}  \\
&= \mathbb{E}_{\mathbf{z} \sim q(\cdot)} \left[ \log\left[ \frac{p(\mathbf{y}, \mathbf{z};  \boldsymbol{\theta})}{q(\mathbf{z})} \right] \right]
\end{aligned}
$$

The ELBO is just an integral over the latent variables, which gives us the <i>expectation</i> part of the algorithm's name.

#### Maximization Step
We then update our parameter values by maximizing the ELBO:

$$
\hat{\boldsymbol{\theta}}^{(t + 1)} = \underset{\boldsymbol{\theta}}{\arg \max} \left\{ \mathbb{E}_{\mathbf{z} \sim q(\cdot)} \left[ \log\left[ \frac{p(\mathbf{y}, \mathbf{z};  \boldsymbol{\theta})}{q(\mathbf{z})} \right] \right] \right\}
$$

---

## Extensions
Sometimes the expectation step is still very hard. There are many extensions that aim to overcome this.

### Hard EM
If we don't really want to think about any uncertainty with respect to the latent variables, we can replace the average over their values by using a point mass at a MAP estimate for the posterior distribution. That is, we pick:

$$
\begin{aligned}
q(\mathbf{z}) &= \delta_{\hat{\mathbf{z}}}(\mathbf{z}) \\
\hat{\mathbf{z}} &= \underset{\mathbf{z}}{\arg \max} \left\{ p(\mathbf{z}; \hat{\boldsymbol{\theta}}^{(t)} \rvert \mathbf{y}) \right\}
\end{aligned}
$$


### Variational EM
Sometimes we do not have access to the posterior. We can instead replace it with some restricted class of proposal distributions. If this class is rich enough to include the true posterior, then the above intuition holds. However, if it does not, then EM is not guaranteed to increase the marginal log-likelihood at each step; the ELBO will still (monotonically) increase. 

In variational EM, we introduce an intermediate step where we optimize over choices of $q$ and then proceed with the EM steps as above. At iteration $t$, we pick:

$$
q^* = \underset{q \in \mathcal{Q}}{\arg \min} \left\{ D_{KL}(q(\mathbf{z}) \rvert \rvert p(\mathbf{z}; \boldsymbol{\theta}^{(t)} \rvert \mathbf{y})) \right\}
$$

for some class $\mathcal{Q}$ (e.g. Gaussian distributions or something). 

---

## Example

### Data Set-Up
We assume that to have a dataset consisting of $N$ observations of: scalar response, $y$; a $p$-dimensional vector of covariates associated with fixed effects, $\mathbf{x}$; a scalar-valued covariate associated with random effects, $z$. We assume the observations are divided into $K$ independent clusters with $n_k$ observations per cluster such that $\sum_{k = 1}^K n_k = N$. Where convenient, we will collect the parameter into a single vector $\boldsymbol{\theta} = (\boldsymbol{\alpha}^\top, \gamma, \tau^2)^\top$ to simplify notation.

We assume that the responses follow a negative binomial distribution with conditonal means and variances given by:

$$
\begin{aligned}
\mathbb{E}[y_{k,i} \rvert \boldsymbol{\beta}] &= \mu_{k,i} = \exp(\eta_{k,i}) = \mathbf{x}_{k,i}^\top \boldsymbol{\alpha} + z_{k,i}\beta_k \\
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
\log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k))
= \sum_{i = 1}^{n_k} \left[ \log\left( \frac{\Gamma(y_{k,i} + \gamma + 2)}{y_{k,i}! \Gamma(\gamma + 2)}\right) + y_{k,i} \left[ \log(\mu_{k,i}) - \log(\mu_{k,i} + \gamma) \right] + \gamma\left[ \log(\gamma) - \log(\mu_{k,i} + \gamma)\right] \right]
$$

where $\Gamma(\cdot)$ denotes the Gamma function. We have placed a Gaussian prior on the random effects:

$$
\log(p(\boldsymbol{\beta}; \tau^2)) = -\frac{1}{2} \sum_{k = 1}^K \left[\log(2 \pi \tau^2) + \frac{\beta_k^2}{\tau^2} \right]
$$

In this case, the posterior is not known in closed form, so we will use a Laplace approximation for our posterior distribution (similar to what we do with adaptive Gaussian quadrature). 

### Approximating The Posterior
Let a hat and a superscript $(t)$ denote a parameter estimate at iteration $t$. We will use a Laplace approximation to the posterior. 

#### Finding The MAP
We will first just consider a single cluster, $k$, and define:

$$
\begin{aligned}
\zeta(\beta_k; \boldsymbol{\theta}) &:= \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2))
\end{aligned}
$$

By Bayes' Theorem, we have:

$$
p(\beta_k; \hat{\boldsymbol{\theta}}^{(t)} \rvert \mathbf{y}_k) 
= \frac{p(\mathbf{y}_k, \beta_k; \hat{\boldsymbol{\theta}}^{(t)})}{p(\mathbf{y}_k; \hat{\boldsymbol{\theta}}^{(t)})} 
= \frac{p(\mathbf{y}_k; \hat{\boldsymbol{\alpha}}^{(t)}, \hat{\gamma}^{(t)} \rvert \beta_k) p(\beta_k; {\hat{\tau}^2}^{(t)})}{p(\mathbf{y}_k; \hat{\boldsymbol{\theta}}^{(t)})}
= \frac{\exp(\zeta(\beta_k; \hat{\boldsymbol{\theta}}^{(t)}))}{p(\mathbf{y}_k; \hat{\boldsymbol{\theta}}^{(t)})}
$$

Since the denominator does not involve $\boldsymbol{\beta}$, we can find the MAP by maximizing the numerator alone, which is the same as maximizing $\zeta(\beta_k)$ for all $k = 1, \dots, K$. Given $\hat{\boldsymbol{\theta}}^{(t)}$, this is:

$$
\begin{aligned}
\beta_k^{(t+1)} 
&= \underset{\beta_k}{\arg \max} \left\{ \exp(\zeta(\beta_k; \hat{\boldsymbol{\theta}}^{(t)})) \right\}
= \underset{\beta_k}{\arg \max} \left\{  \zeta(\beta_k; \hat{\boldsymbol{\theta}}^{(t)}) \right\}
\end{aligned}
$$

<!-- &= \underset{\beta_K}{\arg \max} \left\{ \sum_{i = 1}^{n_k} \left[ y_{k,i}\left( \log(\hat{\mu}^{(t)}_{k,i}) - \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right) - \hat{\gamma}^{(t)} \log(\hat{\mu}^{(t)}_{k,i} + \hat{\gamma}^{(t)}) \right] - \frac{\beta_k^2}{2 {\hat{\tau}^2}^{(t)}} \right\} -->

where we've dropped any term that does not involve $\boldsymbol{\beta}$. This steps gives us the MAP estimate of $\boldsymbol{\beta}$ at iteration $(t+1)$.

#### Laplace Approximation
We now construct an approximation to the posterior distribution of $\beta_k$ given $\boldsymbol{\theta}^{(t)}$. We use a Gaussian distribution with the same mean and variance as a Laplace approximation to the posterior. This requires the first and second derivatives of the posterior with respect to $\beta_k$, which we denote with apostrophes. Since the MAP estimate is the posterior mode, we have:

<aside><p>This is technically not a true Laplace approximation to the posterior because we have to omit the denominator (the marginal distribution, which is unknown).</p></aside>

$$
\begin{aligned}
p(\beta_k; {\hat{\tau}^2}^{(t)} \rvert \mathbf{y}_k)
&\propto \exp(\zeta(\beta_k; \hat{\boldsymbol{\theta}}^{(t)})) \\
&\approx \exp\left(
  \zeta(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) + (\beta_k - \hat{\beta}_k) \zeta'(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) + \frac{1}{2} (\beta_k - \hat{\beta}_k)^2 \zeta''(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)})
\right) \\
&= \exp\left( \zeta(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) \right) \exp\left(\frac{1}{2} (\beta_k - \hat{\beta}_k)^2 \zeta''(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) \right) \\
&= \sqrt{2 \pi {\hat{\tau}^2}^{(t)}} \exp\left( \zeta(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) \right) \phi(\beta_k; \hat{\beta}_k^{(t+1)}, {\hat{\sigma}^2_k}^{(t)})
\end{aligned}
$$

where $\phi(\cdot; \hat{\mu}, \hat{\sigma}^2)$ denotes a Gaussian density with mean $\hat{\mu}$ and variance $\hat{\sigma}^2$ and where ${\hat{\sigma}^2}_k^{(t)}$ is given by:

$$
\begin{aligned}
{\hat{\sigma}^2}_k^{(t)} &= \left[ - \zeta''(\hat{\beta}_k^{(t+1)}; \hat{\boldsymbol{\theta}}^{(t)}) \right]^{-1}
\end{aligned}
$$

<!-- #region q-derivs -->
{% tabs q-derivs %}
{% tab q-derivs statement %}
The derivatives are given by:

$$
\begin{aligned}
\frac{d \zeta(\beta_k; \boldsymbol{\theta})}{d \beta_k} &= - \frac{\beta_k}{2 \tau^2} + \sum_{i = 1}^{n_k} z_{k,i} \left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \\
\frac{d^2 \zeta(\beta_k; \boldsymbol{\theta})}{d \beta_k^2} &= - \left[\frac{1}{\tau^2} + \sum_{i = 1}^{n_k} \frac{(1 + \gamma) z^2_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \right]
\end{aligned}
$$
{% endtab %}
{% tab q-derivs proof %}
First, we have the following intermediate derivatives:

$$
\begin{aligned}
\frac{d}{d \beta_k} \left[ \log(\mu_{k,j}) \right] 
&= \frac{d}{d \beta_k} \left[ \mathbf{x}_{k,j}^\top \boldsymbol{\alpha} + z_{k,j} \beta_k  \right]
= z_{k,j} \\
\frac{d}{d \beta_k} \left[ \mu_{k,j} \right]
&= \exp(\eta_{k,j}) \frac{d}{d \beta_k} \left[ \mathbf{x}_{k,j}^\top \boldsymbol{\alpha} + z_{k,j} \beta_k \right]
= z_{k,j} \mu_{k,j} \\
\frac{d^2}{d \beta_k^2} \left[ \mu_{k,j} \right] 
&=\frac{d}{d \beta_k} \left[ z_{k,j} \mu_{k,j} \right] 
= z_{k,j}^2 \mu_{k,j}
\end{aligned}
$$

The first derivative of $\zeta(\beta_k; \boldsymbol{\theta})$ with respect to $\beta_k$ is:

$$
\begin{aligned}
\frac{d}{d \beta_k} \left[ \log(p(\beta_k; \tau^2)) \right]
&= \frac{d}{d\beta_k} \left[ - \frac{1}{2} \log(2 \pi \tau^2) - \frac{\beta_k^2}{2 \tau^2}\right] 
= - \frac{\beta_k}{2 \tau^2} \\
\frac{d}{d\beta_k} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) \right] 
&= \sum_{i = 1}^{n_k} \left[  \frac{d}{d \beta_k} \left[ \log\left( \frac{\Gamma(y_{k,i} + \gamma + 2)}{y_{k,i}! \Gamma(\gamma + 2)}\right) + y_{k,i} \left[ \log(\mu_{k,i}) - \log(\mu_{k,i} + \gamma) \right] + \gamma\left[ \log(\gamma) - \log(\mu_{k,i} + \gamma)\right] \right] \right] \\
&= \sum_{i = 1}^{n_k} \left[ y_{k,i} \left( \frac{d}{d \beta_k} \left[ \log(\mu_{k,i}) \right] - \frac{d}{d \beta_k}\left[ \log(\mu_{k,i} + \gamma) \right] \right) - \gamma \frac{d}{d \beta_k} \left[ \log(\mu_{k,i} + \gamma) \right] \right]  \\
&= \sum_{i = 1}^{n_k} \left[ y_{k,i} \left(z_{k,i} - \frac{1}{\mu_{k,i} + \gamma} (z_{k,i} \mu_{k,i}) \right) - \frac{\gamma}{\mu_{k,i} + \gamma} (z_{k,i} \mu_{k,i}) \right] \\
&= \sum_{i = 1}^{n_k} z_{k,i} \left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \\
\implies \frac{d}{d \beta_k} \left[ \zeta(\beta_k; \boldsymbol{\theta}) \right]
&= \frac{d}{d \beta_k} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2)) \right] \\
&= - \frac{\beta_k}{2 \tau^2} + \sum_{i = 1}^{n_k} z_{k,i} \left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right]
\end{aligned}
$$

And the second derivative is:

$$
\begin{aligned}
\frac{d^2}{d \beta_k^2} \left[ \log(p(\beta_k; \tau^2)) \right] 
&= \frac{d}{\beta_k} \left[ - \frac{\beta_k}{2 \tau^2} \right]
= - \frac{1}{2 \tau^2} \\
\frac{d^2}{d\beta_k^2} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) \right] 
&= \sum_{i = 1}^{n_k} z_{k,i} \frac{d}{d \beta_k}\left[y_{k,i} - \frac{\mu_{k,i}}{\mu_{k,i} + \gamma} - \frac{\gamma \mu_{k,i}}{\mu_{k,i} + \gamma} \right] \\
&= \sum_{i = 1}^{n_k} z_{k,i} \left[ - \frac{(\mu_{k,i} + \gamma) \frac{d}{d \beta_k} \left[ \mu_{k,i} \right] - \mu_{k,i} \frac{d}{d \beta_k} \left[ \mu_{k,i} + \gamma\right]}{(\mu_{k,i} + \gamma)^2} - \frac{(\mu_{k,i} + \gamma) \frac{d}{d \beta_k} \left[ \gamma \mu_{k,i} \right] - \gamma \mu_{k,i} \frac{d}{d \beta_k} \left[ \mu_{k,i} + \gamma \right]}{(\mu_{k,i} + \gamma)^2} \right] \\
&= - \sum_{i = 1}^{n_k} z_{k,i} \left[ \frac{(\mu_{k,i} + \gamma)(z_{k,i} \mu_{k,i}) - \mu_{k,i}(z_{k,i} \mu_{k,i})}{(\mu_{k,i} + \gamma)^2} + \frac{(\mu_{k,i} + \gamma)(\gamma z_{k,i} \mu_{k,i}) - \gamma \mu_{k,i}(z_{k,i} \mu_{k,i})}{(\mu_{k,i} + \gamma)^2} \right] \\
&= - \sum_{i = 1}^{n_k} z_{k,i} \left[ \frac{z_{k,i} \mu_{k,i}^2 + \gamma z_{k,i} \mu_{k,i} - z_{k,i} \mu_{k,i}^2 + \gamma z_{k,i} \mu_{k,i}^2 + \gamma^2 z_{k,i} \mu_{k,i} - \gamma z_{k,i} \mu_{k,i}^2}{(\mu_{k,i} + \gamma)^2} \right] \\
&= - \sum_{i = 1}^{n_k} z_{k,i} \left[ \frac{\gamma z_{k,i} \mu_{k,i} + \gamma^2 z_{k,i} \mu_{k,i}}{(\mu_{k,i} + \gamma)^2}\right] \\
&= -  \gamma (\gamma + 1)\sum_{i = 1}^{n_k} \frac{z_{k,i}^2 \mu_{k,i}}{(\mu_{k,i} + \gamma)^2} \\
\implies \frac{d^2}{d \beta_k^2} \left[ \zeta(\beta_k; \boldsymbol{\theta}) \right] 
&= \frac{d^2}{d \beta_k^2} \left[ \log(p(\mathbf{y}_k; \boldsymbol{\alpha}, \gamma \rvert \beta_k)) + \log(p(\beta_k; \tau^2)) \right] \\
&= - \frac{1}{2 \tau^2} - \gamma (\gamma + 1)\sum_{i = 1}^{n_k} \frac{z_{k,i}^2 \mu_{k,i}}{(\mu_{k,i} + \gamma)^2}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->


### Algorithm
At iteration $t + 1$, we begin by taking the previous iteration's parameter values and computing an approximation posterior distribution for each cluster:

$$
p(\beta_k; {\hat{\tau}^2}^{(t)} \rvert \mathbf{y}_k)
\approx
\sqrt{2 \pi {\hat{\tau}^2}^{(t)}} \exp(\zeta(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)})) \phi(\beta_k; \hat{\beta}_k^{(t+1)}, {\hat{\sigma}^2}_k^{(t)})
$$

Let $\phi\left(\boldsymbol{\beta}; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)$ denote the <i>multivariate</i> Gaussian density for the joint posterior distribution of the random effects. 

We then find $\hat{\boldsymbol{\theta}}^{(t+1)}$ as:

$$
\begin{aligned}
\hat{\boldsymbol{\theta}}^{(t+1)}
&= \underset{\boldsymbol{\theta}}{\arg \max} \left\{ \int \log\left( \frac{p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta})}{p(\boldsymbol{\beta}; {\hat{\tau}^2}^{(t)} \rvert \mathbf{y})} \right) p(\boldsymbol{\beta}; {\hat{\tau}^2}^{(t)} \rvert \mathbf{y}) d \boldsymbol{\beta}  \right\} \\
&\approx \underset{\boldsymbol{\theta}}{\arg \max} \left\{ \int \log\left( \frac{p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta})}{p(\boldsymbol{\beta}; {\hat{\tau}^2}^{(t)} \rvert \mathbf{y})} \right) (2 \pi {\hat{\tau}^2}^{(t)})^{-\frac{K}{2}} \exp\left(\sum_{k = 1}^K \zeta(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) \right) \phi(\boldsymbol{\beta}; \hat{\boldsymbol{\beta}}^{(t+1)}, {\hat{\boldsymbol{\sigma}}^2}^{(t)})  d \boldsymbol{\beta}  \right\} \\
&= \underset{\boldsymbol{\theta}}{\arg \max} \left\{ (2 \pi {\hat{\tau}^2}^{(t)})^{-\frac{K}{2}} \exp\left(\sum_{k = 1}^K \zeta(\hat{\beta}_k; \hat{\boldsymbol{\theta}}^{(t)}) \right)\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[ \log\left( \frac{p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta})}{p(\boldsymbol{\beta}; {\hat{\tau}^2}^{(t)} \rvert \mathbf{y})} \right) \right]  \right\} \\
&= \underset{\boldsymbol{\theta}}{\arg \max} \left\{\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[ \log\left( p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta}) \right) \right] \right\} \\
\end{aligned}
$$

where in the last line we drop the leading factors and the second term as they are constant with respect to $\boldsymbol{\theta}$ and do not affect the optimization problem. We can compute the expectation as:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[ \log\left( p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta}) \right) \right]
&= \mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[ \log\left( p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right)  \right] + \mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[  \log\left( p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta}) \right) \right] 
\end{aligned}
$$

The second term has a closed form since the prior is Gaussian:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[  \log\left( p(\mathbf{y}, \boldsymbol{\beta}; \boldsymbol{\theta}) \right) \right] 
&= \mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)} \right)} \left[ \log\left( \prod_{k = 1}^K (2 \pi \tau^2)^{-\frac{1}{2}} \exp(-\frac{\beta_k^2}{2 \tau^2}) \right) \right] \\
&= \sum_{k = 1}^K \left[ -\frac{1}{2} \log(2 \pi \tau^2) - \frac{1}{2 \tau^2} \mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)} \right)} \left[ \beta_k^2 \right] \right] \\
&= \sum_{k = 1}^K \left[ -\frac{1}{2} \log(2 \pi \tau^2) - \frac{1}{2 \tau^2} \left( \text{Var}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)} \right)} \left[ \beta_k \right] + \left(\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)} \right)} \left[ \beta_k \right] \right)^2 \right)  \right] \\
&=  \sum_{k = 1}^K \left[ -\frac{1}{2} \log(2 \pi \tau^2) - \frac{1}{2 \tau^2} \left( {\hat{\sigma}^2}_k^{(t)} + (\hat{\beta}_k^{(t+1)})^2 \right)  \right] \\
&=  -\frac{K}{2} \log(2 \pi \tau^2) - \frac{1}{2 \tau^2}\sum_{k = 1}^K \left[ {\hat{\sigma}^2}_k^{(t)} + (\hat{\beta}_k^{(t+1)})^2 \right] \\
\end{aligned}
$$

The first term requires a little more care. We will define:

$$
x_k^2 := \frac{(\beta_k - \hat{\beta}_k^{(t+1)})^2}{2 {\hat{\sigma}^2}_k^{(t)}}; \hspace{5mm} 
\mathbf{x} := (x_1, \dots, x_K)^\top
$$

This implies that the random effects can be rewritten as:

$$
\begin{aligned}
\beta_k &= \sqrt{2 {\hat{\sigma}^2}_k^{(t)}} x_k + \hat{\beta}_k^{(t+1)} \\
\boldsymbol{\beta}^* &:= (\sqrt{2 {\hat{\sigma}^2}_1^{(t)}} x_1 + \hat{\beta}_1^{(t+1)}, \dots, \sqrt{2 {\hat{\sigma}^2}_K^{(t)}} x_K + \hat{\beta}_K^{(t+1)})^\top
\end{aligned}
$$

We also have:

$$
\begin{aligned}
d \beta_k &= \frac{d}{d x_k} \left[  \sqrt{2 {\hat{\sigma}^2}_k^{(t)}} x_k + \hat{\beta}_k^{(t+1)} \right] \\
\implies d \beta_k &= \sqrt{2 {\hat{\sigma}^2}_k^{(t)}} dx_k
\end{aligned}
$$

We can rewrite the expectation as:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[ \log\left( p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right)  \right] 
&= \int \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta})) \phi(\boldsymbol{\beta}; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}) d \boldsymbol{\beta} \\
&= (2 \pi)^{-\frac{K}{2}} \rvert \hat{\boldsymbol{\Sigma}}^{(t)} \rvert^{-\frac{1}{2}} \int \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta})) \exp\left( -\frac{1}{2} (\boldsymbol{\beta} - \hat{\boldsymbol{\beta}}^{(t + 1)})^\top \left[ \hat{\boldsymbol{\Sigma}}^{(t)} \right]^{-1} (\boldsymbol{\beta} - \hat{\boldsymbol{\beta}}^{(t + 1)}) \right) d \boldsymbol{\beta} \\
&= (2 \pi)^{-\frac{K}{2}} \rvert \hat{\boldsymbol{\Sigma}}^{(t)} \rvert^{-\frac{1}{2}} \int \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta})) \exp\left( - \frac{1}{2} \sum_{k = 1}^K \frac{(\beta_k - \hat{\beta}_k^{(t+1)})^2}{ {\hat{\sigma}^2}_k^{(t)} } \right) d \boldsymbol{\beta} 
\\
&= (2 \pi)^{-\frac{K}{2}} \rvert \hat{\boldsymbol{\Sigma}}^{(t)} \rvert^{-\frac{1}{2}} \idotsint \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}^*)) \left( \prod_{k = 1}^K \exp\left( - x_k^2 \right) \sqrt{2 {\hat{\sigma}^2}_k^{(t)}} \right) d x_1 \dots d x_K \\
&=  (2 \pi)^{-\frac{K}{2}} \left( \prod_{k = 1}^K {\hat{\sigma}^2}_k^{(t)} \right)^{-\frac{1}{2}} \left(\prod_{k = 1}^K \sqrt{2 {\hat{\sigma}^2}_k^{(t)}} \right)  \idotsint \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}^*)) \left( \prod_{k = 1}^K \exp\left( - x_k^2 \right)  \right) d x_1 \dots d x_K \\
&=  (2 \pi)^{-\frac{K}{2}} 2^{\frac{K}{2}} \left( \prod_{k = 1}^K {\hat{\sigma}^2}_k^{(t)} \right)^{-\frac{1}{2}} \left(\prod_{k = 1}^K {\hat{\sigma}^2}_k^{(t)} \right)^{\frac{1}{2}}  \idotsint \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}^*)) \left( \prod_{k = 1}^K \exp\left( - x_k^2 \right)  \right) d x_1 \dots d x_K \\
&=  \pi^{-\frac{K}{2}} \idotsint \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}^*)) \left( \prod_{k = 1}^K \exp\left( - x_k^2 \right)  \right) d x_1 \dots d x_K \\
\end{aligned}
$$

The multiple integral in the last line is intractable. Thus, we must use $K$ successive applications of Gauss-Hermite quadrature in order to approximate it. Let $\omega_{i_k}$ and $\nu_{i_k}$ (for $i = 1, \dots, m$) denote the weights and nodes from the Hermite polynomial of order $m$ using in the $k$-th approximation. We then have:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\beta} \sim \phi\left(\cdot; \hat{\boldsymbol{\beta}}^{(t+1)}, \hat{\boldsymbol{\Sigma}}^{(t)}\right)} \left[ \log\left( p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right)  \right] 
&= \pi^{-\frac{K}{2}}\idotsint \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}^*)) \prod_{k = 1}^K \exp\left( - x_k^2 \right) d x_1 \dots d x_K \\
&= \pi^{-\frac{K}{2}} \idotsint \left[ \int \left[ \int \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert \boldsymbol{\beta}^*))  \exp\left( - x_1^2 \right) d x_1 \right] \exp\left(-x_2^2\right) d x_2 \right] \dots \exp\left(-x_K^2\right) d x_K \\
&\approx \pi^{-\frac{K}{2}} \idotsint \left[ \int  \sum_{i_1 = 1}^m \omega_{i_1} \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert (\nu_{i, 1}, \beta_2, \dots, \beta_K))) \exp(-x_2^2) dx_2 \right] \dots \exp(-x_K^2) dx_K \\
&\approx \pi^{-\frac{K}{2}}  \idotsint \sum_{i_2 = 1}^m  \omega_{i_2} \left( \sum_{i_1 = 1}^m \omega_{i_1} \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert (\nu_{i_1}, \nu_{i_2}, \dots, \beta_K))) \right) \dots \exp(-x_K^2) dx_K \\ 
&\approx \pi^{-\frac{K}{2}}\sum_{i_K = 1}^m \omega_{i_K} \left(\dots  \left( \sum_{i_2 = 1}^m  \omega_{i_2} \left( \sum_{i_1 = 1}^m \omega_{i_1} \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert (\nu_{i_1}, \nu_{i_2}, \dots, \nu_{i_K}))) \right) \right) \right) \\
&= \pi^{-\frac{K}{2}}  \sum_{i_K = 1}^m \dots \sum_{i_1 = 1}^m \left( \prod_{k = 1}^K \omega_{i_k}\right) \log(p(\mathbf{y}; \boldsymbol{\theta} \rvert (\nu_{i_1}, \nu_{i_2}, \dots, \nu_{i_K}))) 
\end{aligned}
$$

Unfortunately, quadrature results in having to compute $m^K$ sets of points because we have to apply it $K$ times. This is leading to memory issues on my computer <3 so I am going to just try sampling instead.
