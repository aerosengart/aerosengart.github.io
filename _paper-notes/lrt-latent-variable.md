---
layout: distill
title: A Note on Likelihood Ratio Tests for Models with Latent Variables
description: 
date: 2026-05-04
tabs: true
tags: likelihood testing glmm
toc:
  - name: Background
  - name: General Problem
  - name: Regularity Conditions
  - name: Example
bibliography: paper-notes.bib
---

In this post I cover an overview paper<d-cite key=chen2020></d-cite> of using likelihood ratio tests in non-standard settings, specifically for latent variable models. I'm mostly interested in the mixed model case, so I will focus on those sections. However, the paper has broader applications.


## Background
We will assume to have a sample of $J$-dimensional observations denoted by $\mathbf{y}_1$, $\dots$, $\mathbf{y}_n$, which are independent and identically distributed observations such that $$\mathbf{y}_i = (y_{i,1}, \dots, y_{i,J})^\top$$. The observations will have a distribution from some family parametrized by $k$-dimensional vector, $\boldsymbol{\theta}$:

$$
\mathcal{P}_{\boldsymbol{\Theta}} = \{P_\boldsymbol{\theta} \rvert \boldsymbol{\theta} \in \boldsymbol{\Theta} \subset \mathbb{R}^k \}
$$

<aside><p>I believe $\boldsymbol{\Theta}$ must be open as well.</p></aside>

We assume that all $$P_{\boldsymbol{\theta}} \in \mathcal{P}_{\boldsymbol{\Theta}}$$ are dominated by some $\sigma$-finite measure, $\nu$. We will also assume that each $$P_{\boldsymbol{\theta}}$$ has a probability density function, $$p_{\boldsymbol{\theta}}: \mathbb{R}^J \rightarrow [0, \infty)$$ defined with respect to $\nu$.

<details>
<summary>What does that mean?</summary>
We can consider our observations from a measure theoretic perspective. 
<br>
Our probability space will be denoted by $(\Omega, \mathcal{F}, P_{\boldsymbol{\theta}})$, where $\Omega$ is our sample space, $\mathcal{F}$ is a $\sigma$-algebra over $\Omega$, and $P_{\boldsymbol{\theta}}$ is a probability measure over $(\Omega, \mathcal{F})$. 
Each observation can be thought of as a realization of some random variable, $\mathbf{y}: \Omega \rightarrow \mathbb{R}^J$, which is an $\mathcal{F}$-measureable function.
<br>
Since each $P_{\boldsymbol{\theta}}$ is a probability measure, it is finite and therefore a finite measure. This implies it is also <strong>$\sigma$-finite</strong>, which means that there exist countably man $F_1, F_2, \dots \in \mathcal{F}$ such that $P_{\boldsymbol{\theta}}(F_n) < \infty$ for all $n \in \mathbb{N}$ such that $\Omega = \bigcup_{n \in \mathbb{N}} F_n$. 
<br>
A $\sigma$-finite measure, $\nu$, is said to <strong>dominate</strong> all $P_{\boldsymbol{\theta}}$ if $\nu << P_{\boldsymbol{\theta}}$ for all $\boldsymbol{\theta}$. This implies that, for any $F \in \mathcal{F}$, $P_{\boldsymbol{\theta}}(F) = 0$ implies $\nu(F) = 0$ as well.
<br>
The <strong>density</strong> is defined as the Radon-Nikodym derivative, which we will not cover in detail here. Briefly, a density function, $p_{\boldsymbol{\theta}}$, is one that satisfies:
$$
\mathbb{P}(\mathbf{y} \in F) - \int_{\mathbf{y}^{-1}(F)} dP_{\boldsymbol{\theta}} = \int_F p_\boldsymbol{\theta} d \nu
$$
for any $F \in \mathcal{F}$.
</details>

We will use $\boldsymbol{\Sigma}_{pd}^{J \times J}$ and $\boldsymbol{\Sigma}_d^{J \times J}$ to denote the spaces of $J \times J$ strictly positive definite and diagonal matrices, respectively. Furthermore, define the mappings $$\rho: \boldsymbol{\Sigma}_{pd}^{J \times J} \mapsto \mathbb{R}^{\frac{J(J+1)}{2}}$$ and $$\mu: \boldsymbol{\Sigma}_d^{J \times J} \mapsto \mathbb{R}^J$$ such that:

$$
\begin{aligned}
\rho(\boldsymbol{\Sigma}) &= (\Sigma_{1,1}, \Sigma_{1,2}, \dots, \Sigma_{1,J}, \Sigma_{2,2}, \dots, \Sigma_{2,J}, \dots, \Sigma_{J,J})^\top &\text{ for } \boldsymbol{\Sigma} \in \boldsymbol{\Sigma}_{pd}^{J \times J} \\
\mu(\boldsymbol{\Sigma}) &= (\Sigma_{1,1}, \Sigma_{2,2}, \dots, \Sigma_{J,J})^\top &\text{ for } \boldsymbol{\Sigma} \in \boldsymbol{\Sigma}_{d}^{J \times J}
\end{aligned}
$$

That is, $\rho(\boldsymbol{\Sigma})$ returns the vector of upper triangular entries (row-major order), and $\mu(\boldsymbol{\Sigma})$ returns the vector of diagonal entries.

---

## General Problem 
In general, we are interested in testing hypotheses of the form:

$$
\begin{equation}
\label{eq:general-hypothesis}
H_0: \boldsymbol{\theta} \in \boldsymbol{\Theta}_0 \hspace{5mm} \text{ vs. } \hspace{5mm} H_1: \boldsymbol{\theta} \in \boldsymbol{\Theta} \setminus \boldsymbol{\Theta}_0
\end{equation}
$$

We consider testing via the likelihood ratio test (LRT) statistic. Let $\boldsymbol{\theta}^* \in \boldsymbol{\Theta}_0$ be the true value of $\boldsymbol{\theta}$. The log-likelihood function, due to the i.i.d. assumption, is given by:

$$
\begin{equation}
\label{eq:loglik}
\ell(\boldsymbol{\theta}) = \sum_{i = 1}^n \log (p_{\boldsymbol{\theta}}(\mathbf{y}_i))
\end{equation}
$$

The LRT statistic is then defined as:

$$
\begin{equation}
\label{eq:lrt}
\lambda = 2 \left[ \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \ell(\boldsymbol{\theta}) \right\} - \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\sup} \left\{ \ell(\boldsymbol{\theta}) \right\} \right]
\end{equation}
$$

---

## Regularity Conditions
The asymptotic distribution of $\lambda$ traditionally relies upon the idea of <i>local asymptotic normality</i>.<d-cite key=vaart2012></d-cite> In essence, we replace our hypothesis by constructing a new one that involves a Gaussian approximation via Taylor approximation. This allows us to use the nice Gaussian distribution to derive asymptotic results. To use this idea, however, we need to establish that the maximum likelihood estimator of our parameter is asymptotically normal. That is, that:

$$
\sqrt{n} (\hat{\boldsymbol{\theta}}_{\boldsymbol{\Theta}} - \boldsymbol{\theta}^*) \rightsquigarrow \mathcal{N}\left(\mathbf{0}_k, \mathcal{I}^{-1}(\boldsymbol{\theta}^*)\right)
$$

The above holds if Conditions 1-4 hold. 

#### Condition 1
The true parameter, $\boldsymbol{\theta}^* \in \boldsymbol{\Theta}_0$, lies in the interior of $\boldsymbol{\Theta}$. 

#### Condition 2a
The family $\mathcal{P}_{\boldsymbol{\theta}}$ is quadratic mean differentiable at $\boldsymbol{\theta}^*$. 

<div id="qmd"></div>
<div class=definition>
<strong>Definition (Quadratic Mean Differentiability).<d-cite key=vaart2012></d-cite></strong>
<br>
Suppose we are in the setting explained in the <a href="#background">Background</a> section. A family of distributions, $\mathcal{P}_{\boldsymbol{\theta}} = \left\{ P_{\boldsymbol{\theta}} \rvert \boldsymbol{\theta} \in \boldsymbol{\Theta} \subset \mathbb{R}^k \right\}$, is called <i>quadratic mean differentiable (QMD)</i> at $\boldsymbol{\theta}^*$ if there exists a vector of real-valued functions, $\boldsymbol{\eta}(\cdot, \boldsymbol{\theta}^*) = (\eta_1(\cdot, \boldsymbol{\theta}^*), \dots, \eta_k(\cdot, \boldsymbol{\theta}^*))^\top$, satisfying:
$$
\underset{\mathbf{h} \rightarrow \mathbf{0}}{\lim} \left\{ 
   \int_{\mathbb{R}^J} \left( \sqrt{p_{\boldsymbol{\theta}^* + \mathbf{h}}(\mathbf{y})} - \sqrt{p_{\boldsymbol{\theta}^*}(\mathbf{y})} - \frac{1}{2} \mathbf{h}^\top \boldsymbol{\eta}(\mathbf{y}, \boldsymbol{\theta}^*) \sqrt{p_{\boldsymbol{\theta}^*}(\mathbf{y})}\right)^2 d\nu(\mathbf{y})
\right\} = o(\rvert \rvert \mathbf{h} \rvert \rvert^2_2)
$$
The function, $\boldsymbol{\eta}(\cdot, \boldsymbol{\theta}^*)$ is called the <i>quadratic mean derivative of $\mathcal{P}_{\boldsymbol{\theta}}$ at $\boldsymbol{\theta}^*$</i>. 
<br>
If $\mathcal{P}_{\boldsymbol{\theta}}$ is QMD at all $\boldsymbol{\theta}^*$, then it is just called <i>quadratic mean differentiable</i>.
</div>

The above condition can be more readily seen as relevant to estimation settings if we define $\mathbf{h} = \boldsymbol{\theta} - \boldsymbol{\theta}^*$. 

#### Condition 2b
The Fisher information matrix for $\mathcal{P}_\boldsymbol{\theta}$ is positive definite. 

<div id="fi-mat"></div>
<div class=definition>
<strong>Definition (Fisher Information Matrix).<d-cite key=lehmann2008></d-cite></strong>
<br>
Let $\mathcal{P}_\boldsymbol{\theta}$ be a QMD family of distributions with derivative $\boldsymbol{\eta}(\cdot, \boldsymbol{\theta})$. The <i>Fisher information matrix</i> for $\mathcal{P}_\boldsymbol{\theta}$, denoted by $\mathcal{I}(\boldsymbol{\theta})$, is the $k \times k$ matrix with entries:
$$
\mathcal{I}_{i,j}(\boldsymbol{\theta}) = 4 \int_{\mathbb{R}^J} \eta_i(\mathbf{y}, \boldsymbol{\theta}) \eta_j(\mathbf{y}, \boldsymbol{\theta}) d \nu(\mathbf{y})
$$
It is independent of the choice of $nu$.
</div>

This condition concerns the convexity of the likelihood function near the MLE. A function is strongly convex if, and only if, its Hessian is positive definite. If the Fisher information matrix is positive definite, then the likelihood function will be strongly convex near the MLE. This means that problems that rely on numerical optimization will converge quickly to the MLE, and we will be able to compute the likelihood ratio accurately (since we will be able to actually <i>get</i> to the maxima). 



#### Condition 3
There exists a neighborhood of $$\boldsymbol{\theta}^*$$, $$U_{\boldsymbol{\theta}^*} \subset \boldsymbol{\Theta}$$, and a measurable function, $M(\mathbf{y})$, satisfying:

$$
\int_{\mathbb{R}^J} M^2(\mathbf{y})  dP_{\boldsymbol{\theta}^*}(\mathbf{y}) < \infty
$$

<aside><p>$M(\mathbf{y})$ is <strong>square-integrable</strong>.</p></aside>

such that, for any $$\boldsymbol{\theta}_1, \boldsymbol{\theta}_2 \in U_{\boldsymbol{\theta}^*}$$, the following holds:

$$
\rvert \log(p_{\boldsymbol{\theta}_1}(\mathbf{y})) - \log(p_{\boldsymbol{\theta}_2}(\mathbf{y}))
\rvert \leq M(\mathbf{y}) \rvert \rvert \boldsymbol{\theta}_1 - \boldsymbol{\theta}_2 \rvert \rvert_2
$$

#### Condition 4
The following maximum likelihood estimators (MLEs) are consistent under $P_{\boldsymbol{\theta}^*}$:

$$
\begin{aligned}
\hat{\boldsymbol{\theta}}_{\boldsymbol{\Theta}} &= \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\arg \max} \left\{ \ell(\boldsymbol{\theta}) \right\} \\
\hat{\boldsymbol{\theta}}_{\boldsymbol{\Theta}_0} &= \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\arg \max} \left\{  \ell(\boldsymbol{\theta}) \right\}
\end{aligned}
$$

### Additional Conditions
Suppose Condition 1 or 2 is violated by our setting. In order to adapt the LRT to these kinds of cases, we must reformulate our test. We instead suppose that we are testing nested sub-models of a full/saturated model. We assume, now, that $\boldsymbol{\theta}^*$ is in the <i>interior</i> of $\boldsymbol{\Theta}$ and consider testing:

$$
H_0: \boldsymbol{\theta} \in \boldsymbol{\Theta} \hspace{5mm} \text{vs.} \hspace{5mm} H_1: \boldsymbol{\theta} \in \boldsymbol{\Theta}_1 \setminus \boldsymbol{\Theta}_0
$$

where $\boldsymbol{\Theta}_0 \subset \boldsymbol{\Theta}_1 \subset \boldsymbol{\Theta} \subset \mathbb{R}^k$. 

The asymptotic distribution can be derived under a few additional conditions.

#### Condition 5
For any $$\mathbf{t} \in \mathcal{T}_{\boldsymbol{\Theta}_0}(\boldsymbol{\theta}^*)$$, there exist $\epsilon > 0$ and $$\boldsymbol{\alpha}: [0, \epsilon) \rightarrow \boldsymbol{\Theta}_0$$ where $$\boldsymbol{\alpha}(0) = \boldsymbol{\theta}^*$$, such that:

$$
\mathbf{t} = \underset{t \rightarrow 0^+}{\lim} \frac{\boldsymbol{\alpha}(t) - \boldsymbol{\alpha}(0)}{t}
$$

where $$\mathcal{T}_{\boldsymbol{\Theta}_0}(\boldsymbol{\theta}^*)$$ is the <strong>tangent cone</strong> of $\boldsymbol{\Theta}_0$ at $$\boldsymbol{\theta}^*$$ defined as:

$$
\mathcal{T}_{\boldsymbol{\Theta}_0}(\boldsymbol{\theta}^*) = \left\{ \mathbf{v} \in \mathbb{R}^k \rvert \mathbf{v} = \underset{n \rightarrow \infty}{\lim} \left\{ r_n (\boldsymbol{\theta}_n - \boldsymbol{\theta}^*) \right\}; r_n \in \mathbb{R}_{>0}; \boldsymbol{\theta}_n \in \boldsymbol{\Theta}_0 \text{ with } \boldsymbol{\theta}_n \rightarrow \boldsymbol{\theta}^* \right\}
$$

This condition ensures that the shape of the null parameter space around the true parameter value behaves nicely even when on the boundary.

#### Condition 6
The following maximum likelihood estimator is consistent under $P_{\boldsymbol{\theta}^*}$:

$$
\begin{aligned}
\hat{\boldsymbol{\theta}}_{\boldsymbol{\Theta}_1} &= \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_1}{\arg \max} \left\{ \ell(\boldsymbol{\theta}) \right\} \\
\end{aligned}
$$

#### Condition 7
For any $$\mathbf{t} \in \mathcal{T}_{\boldsymbol{\Theta}_1}(\boldsymbol{\theta}^*)$$, there exist $\epsilon > 0$ and $$\boldsymbol{\alpha}: [0, \epsilon) \rightarrow \boldsymbol{\Theta}_1$$ where $$\boldsymbol{\alpha}(0) = \boldsymbol{\theta}^*$$, such that:

$$
\mathbf{t} = \underset{t \rightarrow 0^+}{\lim} \frac{\boldsymbol{\alpha}(t) - \boldsymbol{\alpha}(0)}{t}
$$

where $$\mathcal{T}_{\boldsymbol{\Theta}_1}(\boldsymbol{\theta}^*)$$ is the <strong>tangent cone</strong> of $\boldsymbol{\Theta}_1$ at $$\boldsymbol{\theta}^*$$. 

<div id="chen-theorem-2"></div>
<div class=theorem>
<strong>Theorem 2.<d-cite key=chen2020></d-cite></strong>
<br>
Suppose Conditions 1-7 hold. Let $\boldsymbol{\theta}^* \in \boldsymbol{\Theta}_0$ be the true parameter vector. The likelihood ratio statistic based upon a sample size of $n$, $\lambda_n$, satisfies:
$$
\lambda \rightsquigarrow \underset{\mathbf{t} \in \mathcal{T}_{\boldsymbol{\Theta}_0}(\boldsymbol{\theta}^*)}{\min} \left\{ \rvert \rvert \mathbf{z} - \mathcal{I}^{\frac{1}{2}}(\boldsymbol{\theta}^*) \mathbf{t} \rvert \rvert_2^2 \right\} - \underset{\mathbf{t} \in \mathcal{T}_{\boldsymbol{\Theta}_1}(\boldsymbol{\theta}^*)}{\min} \left\{ \rvert \rvert \mathbf{z} - \mathcal{I}^{\frac{1}{2}}(\boldsymbol{\theta}^*) \mathbf{t} \rvert \rvert_2^2 \right\}
$$
where $\mathbf{z} = (z_1, \dots, z_k)^\top$ with $z_i \sim \mathcal{N}(0, 1)$, and $\mathcal{I}^{\frac{1}{2}} (\boldsymbol{\theta}^*)$ satisfies:
$$
\mathcal{I}(\boldsymbol{\theta}^*) = \mathcal{I}^{\frac{1}{2}} (\boldsymbol{\theta}^*) \left[\mathcal{I}^{\frac{1}{2}} (\boldsymbol{\theta}^*) \right]^\top
$$
</div>

---

## Example
As a concrete example, we will assume a simple linear mixed effects model where the responses come from $K$ equally sized groups. We let $y_{i,j}$ denote the response for the $j$-th individual in group $i$. Our model is:

$$
\begin{equation}
\label{eq:ri-model}
\begin{aligned}
y_{i,j} &= \alpha^0_i + \alpha_1 x_{i,j} + \beta_i x_{i,j} + \epsilon_{i,j} \\
\beta_i &\sim \mathcal{N}(0, \tau^2) \\
\epsilon_{i,j} &\sim \mathcal{N}(0, \sigma^2)
\end{aligned}
\end{equation}
$$

We are interested in testing the following hypothesis:

$$
\begin{equation}
\label{eq:h0}
H_0: \alpha_1 = 0 \text{ and } \tau^2 = 0 \hspace{5mm} \text{ vs. } \hspace{5mm} H_1: \alpha_1 \neq 0 \text{ and } \tau^2 > 0
\end{equation}
$$

We reframe our hypothesis test in terms of the parameter spaces. Let $\boldsymbol{\alpha}^0 = (\alpha_1^0, \dots, \alpha_K^0)^\top$. Note that the marginal covariance matrix is $\boldsymbol{\Sigma} = \sigma^2 \mathbb{I}_{J \times J} + \tau^2 \mathbf{1}_J \mathbf{1}_j^\top$. Thus:

$$
\begin{aligned}
\boldsymbol{\Theta} &= \left\{ (\boldsymbol{\alpha}^0, \alpha_1, \rho^\top(\boldsymbol{\Sigma}))^\top \rvert \boldsymbol{\alpha}^0 \in \mathbb{R}^K, \alpha_1 \in \mathbb{R}, \boldsymbol{\Sigma} \in \mathbb{R}_{pd}^{J \times J} \right\} \\
\boldsymbol{\Theta}_0 &= \left\{ (\boldsymbol{\alpha}^0, \alpha_1, \rho^\top(\boldsymbol{\Sigma}))^\top \rvert \boldsymbol{\alpha}^0 \in \mathbb{R}^K, \alpha_1 = 0, \boldsymbol{\Sigma} = \sigma^2 \mathbb{I}_{J \times J}, \sigma^2 > 0 \right\} \\
\boldsymbol{\Theta}_1 &= \left\{ (\boldsymbol{\alpha}^0, \alpha_1, \rho^\top(\boldsymbol{\Sigma}))^\top \rvert \boldsymbol{\alpha}^0 \in \mathbb{R}^K, \alpha_1 = 0, \boldsymbol{\Sigma} = \sigma^2 \mathbb{I}_{J \times J} + \tau^2 \mathbf{1}_J \mathbf{1}_J^\top, \tau^2 \geq 0 , \sigma^2 > 0 \right\}
\end{aligned}
$$

Let $\boldsymbol{\theta}^* \in \boldsymbol{\Theta}_0$ denote the true parameter vector. Clearly, Condition 1 holds because, in this case, the marginal covariance matrix is ${\sigma^2}^* \mathbb{I}_J$ with ${\sigma^2}^* > 0$, which lies in the interior of the set of all positive definite matrices. We can then define the tangent cones for the null and alternative parameter spaces.

$$
\begin{aligned}
\mathcal{T}_{\boldsymbol{\Theta}_0}(\boldsymbol{\theta}^*) 
    &= \left\{ (b^0_1, \dots, b^0_K, b_1, \rho(\boldsymbol{\Sigma}))^\top \rvert b^0_1, \dots, b^0_K, b_1, b_2 \in \mathbb{R}, \boldsymbol{\Sigma} = b_2 \mathbb{I}_{J \times J} \right\} \\
\mathcal{T}_{\boldsymbol{\Theta}_1}(\boldsymbol{\theta}^*) 
    &= \left\{ (b^0_1, \dots, b^0_K, b_1, \rho(\boldsymbol{\Sigma}))^\top \rvert b^0_1, \dots, b^0_K, b_1, b_2 \in \mathbb{R}, \boldsymbol{\Sigma} = b_2 \mathbb{I}_{J \times J} + b_3 \mathbf{1}_J \mathbf{1}_J^\top, b_3 \geq 0 \right\}
\end{aligned}
$$

We can get the asymptotic distribution of $\lambda$, the LRT statistic, via Theorem 2. In addition, it can be shown that this is a mixture of $\chi^2$ variables. 



