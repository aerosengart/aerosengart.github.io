---
layout: distill
title: Constrained Statistical Inference (Chapter 4)
description: Tests in General Parametric Models
date: 2026-07-13
tabs: true
tags: likelihood theory testing constraints
toc:
  - name: Background
    subsections:
        - name: Notation
        - name: Assumptions
  - name: No Constraints
    subsections:
        - name: Simple Case
        - name: Nuisance Parameters
  - name: Linear Constraints
    subsections:
        - name: Likelihood Ratio Statistic
        - name: Wald-Type Statistic
        - name: Distance Statistic
        - name: Score Statistic
        - name: Extension to Non-IID Case
  - name: Non-Linear Constraints
    subsections:
        - name: Cones
bibliography: stats-ml.bib
---

In this post, I will be going through some of the content in Chapter 4 of Silvapulle and Sen's <i>Constrained Statistical Inference: Inequality, Order, and Shape Restrictions</i>.<d-cite key=silvapulle2005></d-cite>. The chapter concerns testing parameter vectors more generally (i.e. not mean vectors for multivariate Gaussian data).


## Background
### Notation
We will assume to have independent and identically distributed $\mathbf{x}_1, \dots, \mathbf{x}_n$ with common distribution $F(\mathbf{x}; \boldsymbol{\theta})$ and probability density function $f(\mathbf{x}; \boldsymbol{\theta})$ parametrized by $\boldsymbol{\theta} \in \boldsymbol{\Theta} \subset \mathbb{R}^p$. 

We write the log-likelihood, score function, and information matrix for one observation as:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{x}) &= \log f(\mathbf{x}; \boldsymbol{\theta}) \\
\mathbf{S}(\boldsymbol{\theta}) &= \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x})}{\partial \boldsymbol{\theta}} \\
\mathcal{I}(\boldsymbol{\theta}) &= \mathbb{E}_{\boldsymbol{\theta}} \left[ \frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}^\top}\right]
\end{aligned}
$$

We will let $$\boldsymbol{\theta}_0$$ denote the true value of the parameter vector, and we will let $$\hat{\boldsymbol{\theta}}$$ denote the global maximum lieklihood estimate of the parameter vector over all of $\boldsymbol{\Theta}$. For some $$\mathbf{A} \in \boldsymbol{\Theta}$$ such that $$\boldsymbol{\theta}_0 \in \mathbf{A}$$, we define the MLE over $\mathbf{A}$ as the value of the parameter vector, $$\hat{\boldsymbol{\theta}}_\mathbf{A}$$, satisfying:

$$
\ell(\hat{\boldsymbol{\theta}}_\mathbf{A}; \mathbf{x})
  = \underset{\boldsymbol{\theta} \in \mathbf{A}}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{x}) \right\} + o_p(1)
$$

For some sequence $\{ \mathbf{x}_n \}_n$, we will denote its probability limit when $\boldsymbol{\theta}$ is true value with $$\underset{\boldsymbol{\theta}}{\text{plim}} [\mathbf{x}_n]$$, so:

$$
\underset{n \rightarrow \infty}{\lim} \left\{ \rvert \mathbf{x}_n - \underset{\boldsymbol{\theta}}{\text{plim}} [\mathbf{x}_n] \rvert > \epsilon \right\} = 0; \hspace{2mm} \forall \epsilon > 0 
$$

### Assumptions
We assume that $\hat{\boldsymbol{\theta}}_\mathbf{A}$ is consistent. This means that it converges in probability to the true value $\boldsymbol{\theta}_0$.

<aside><p><strong>Assumption A1.</strong></p></aside>

In what follows, we assume that $\hat{\boldsymbol{\theta}}_\mathbf{A}$ is consistent when $\boldsymbol{\theta}_0 \in \mathbf{A}$. We also assume that the following regularity conditions hold:

<div id="condition-q"></div>
<ol>
<li>If $F(\mathbf{x}; \boldsymbol{\theta}_1) = F(\mathbf{x}; \boldsymbol{\theta}_2)$ for $\boldsymbol{\theta}_1, \boldsymbol{\theta}_2 \in \boldsymbol{\Theta}$, then $\boldsymbol{\theta}_1 = \boldsymbol{\theta}_2$.</li>
<li>$\frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}$, $\frac{\partial^2 \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}$, and $\frac{\partial^3 \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_k}$ exist almost everywhere for all $i, j, k = 1, \dots, p$.</li>
<li>These exists $G(\mathbf{x})$ such that $\int G(\mathbf{x}) d\mathbf{x} < \infty$ and the absolute values of the partial derivatives of $\log f(\mathbf{x}; \boldsymbol{\theta})$ with respect to $\boldsymbol{\theta}$ up to order three are bounded by $G(\mathbf{x})$ in a neighborhood of the true parameter value, $\boldsymbol{\theta}_0$.</li>
<li>$\mathcal{I}(\boldsymbol{\theta})$ has finite entries and is positive definite.</li>
</ol>

<aside><p><strong>Condition Q.</strong></p></aside>

The conditions above ensure that the null distribution of the likelihood ratio test statistic for $\boldsymbol{\theta} = \boldsymbol{\theta}_0$ against $\boldsymbol{\theta} \neq \boldsymbol{\theta}_0$ is asymptotically $\chi^2$ and that the global MLE is (suitably centered and scaled) asymptotically normal.

These conditions imply the following results.

<div id="prop-4-2-1"></div>
<div class="theorem">
  <strong>Proposition 4.2.1.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  {% tabs p-4-2-1 %}
  {% tab p-4-2-1 statement %}
  Assume the above regularity conditions hold. Then:
  <ol>
  <li>$\frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \mathcal{I}(\boldsymbol{\theta}_0))$</li>
  <li>$\frac{1}{\sqrt{n}} \mathcal{I}^{-1}(\boldsymbol{\theta}_0)\mathbf{S}(\boldsymbol{\theta}_0) \rightsquigarrow \sqrt{n}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0) + o_p(1)$</li>
  <li>$\sqrt{n}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}) \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \mathcal{I}^{-1}(\boldsymbol{\theta}_0))$</li>
  <li>Under $H_0$, the likelihood ratio test statistic for testing $\mathbf{R} \boldsymbol{\theta} = \mathbf{0}$ against $\mathbf{R}\boldsymbol{\theta} \neq \mathbf{0}$ is asymptotically $\chi_r^2$ distributed where $r$ is the rank of $\mathbf{R}$.</li>
  </ol>
  {% endtab %}
  {% tab p-4-2-1 proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

As noted by Silvapulle and Sen, the proofs of the above statements are done with approximations to the log-likelihood function and gradient. In what follows, we often assume that the MLE is $\sqrt{n}$-consistent, which means that $$\sqrt{n}(\hat{\boldsymbol{\theta}}_\mathbf{A} - \boldsymbol{\theta}_0) = O_p(1)$$ for $$\boldsymbol{\theta}_0 \in \mathbf{A} \subset \boldsymbol{\Theta}$$. However, if the regularity conditions hold, and $$\hat{\boldsymbol{\theta}}_\mathbf{A}$$ is consistent, then it will be $\sqrt{n}$-consistent as well.

<div id="prop-4-2-2"></div>
<div class="theorem">
  <strong>Proposition 4.2.2.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  {% tabs p-4-2-2 %}
  {% tab p-4-2-2 statement %}
  Assume the above regularity conditions in Condition Q hold and that we have $\mathbf{u} = \sqrt{n}(\boldsymbol{\theta} - \boldsymbol{\theta}_0)$ and $K > 0$. Then:

  <strong>(A)</strong>
  For $\underset{\rvert \rvert \mathbf{u} \rvert \rvert_2 < K}{\sup} \left\\{ \rvert r_n(\mathbf{u}) \rvert \right\\} = o_p(1)$:

  $$
  \begin{equation}
  \label{eq:quad-approx-1}
  \ell(\boldsymbol{\theta}; \mathbf{x}) = \ell(\boldsymbol{\theta}_0; \mathbf{x}) + \frac{1}{\sqrt{n}} \mathbf{u}^\top \mathbf{S}(\boldsymbol{\theta}_0) - \frac{1}{2}\mathbf{u}^\top \mathcal{I}(\boldsymbol{\theta}_0) \mathbf{u} + r_n(\mathbf{u}) 
  \end{equation}
  $$

  <strong>(B)</strong>
  For $\mathbf{z}_n = \frac{1}{\sqrt{n}} \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0)$, $\underset{\rvert \rvert \mathbf{u} \rvert \rvert_2 < K}{\sup} \left\\{ \rvert \delta_n(\mathbf{u}) \rvert \right\\} = o_p(1)$, and $\delta_n(\mathbf{u}) = \frac{1}{\sqrt{n}} \rvert \rvert \mathbf{u} \rvert \rvert_2^3 O_p(1)$:

  $$
  \begin{equation}
  \label{eq:quad-approx-2} 
  \ell(\boldsymbol{\theta}; \mathbf{x}) = \ell(\boldsymbol{\theta}_0; \mathbf{x}) + \frac{1}{2n} \mathbf{S}^\top(\boldsymbol{\theta}_0) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0) - \frac{1}{2}(\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0) (\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u}) 
  \end{equation}
  $$

  <strong>(C)</strong>
  For $\mathbf{z}_n = \sqrt{n}(\hat{\boldsymbol{\theta} - \boldsymbol{\theta}}_0)$ where $\hat{\boldsymbol{\theta}}$ is $\sqrt{n}$-consistent and $\underset{\rvert \rvert \mathbf{u} \rvert \rvert < K}{\sup} \left\\{ \rvert \epsilon_n(\mathbf{u}) \rvert \right\\} = o_p(1)$:

  $$
  \begin{equation}
  \label{eq:quad-approx-3} 
  \ell(\boldsymbol{\theta}; \mathbf{x}) = \ell(\hat{\boldsymbol{\theta}}; \mathbf{x}) - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \epsilon_n(\mathbf{u}) 
  \end{equation}
  $$
  {% endtab %}
  {% tab p-4-2-2 proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

Silvapulle and Sen note that the $\sqrt{n}$-consistency comes in handy when we consider the above quadratic approximations. By definition, $$\hat{\boldsymbol{\theta}}_\mathbf{A}$$ satisfies $$\sqrt{n}(\hat{\boldsymbol{\theta}}_\mathbf{A} - \boldsymbol{\theta}_0) = O_p(1)$$ if it is $\sqrt{n}$-consistent. This means that there exists some $$M \in \mathbb{R}_{> 0}$$ (and therefore for some $$K \in \mathbb{R}_{>0}$$) such that:

$$
\sqrt{n} (\hat{\boldsymbol{\theta}}_{\mathbf{A}} - \boldsymbol{\theta}_0) < M
\implies
\sqrt{n} \rvert \rvert\hat{\boldsymbol{\theta}}_{\mathbf{A}} - \boldsymbol{\theta}_0\rvert\rvert_2 < K
$$

It then suffices to consider only those values of $\boldsymbol{\theta}$ that are in this neighborhood of $\boldsymbol{\theta}_0$ when considering the global maximum of the log-likelihood, while allows us to ignore the remainder terms as they will vanish uniformly on this neighborhood. 

<aside><p>Why?</p></aside>

---

## No Constraints
We will assume the data are i.i.d. and that <a href="#assumptions">Condition Q</a> holds. Let $\mathbf{x}$ denote the collection of $n$ i.i.d. datapoints $\mathbf{x}_1, \dots, \mathbf{x}_n$, so $\ell(\boldsymbol{\theta}; \mathbf{x})$ is the <i>full</i> data log-likelihood (the sum of the log-likelihood values for each observation). We will also assume $\boldsymbol{\Theta} \subset \mathbb{R}^p$ is open. 

### Simple Case
Consider testing:

$$
\begin{equation}
\label{eq:hypotheses-0}
H_0: \boldsymbol{\theta} = \boldsymbol{\theta}_0
\hspace{5mm} \text{vs} \hspace{5mm}
H_1: \boldsymbol{\theta} \neq \boldsymbol{\theta}_0
\end{equation}
$$

Since $\log(x) \leq x - 1$, we have:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}}\left[ \log\left( \frac{f(\mathbf{x}; \boldsymbol{\theta}^*)}{f(\mathbf{x}; \boldsymbol{\theta})} \right) \right] 
&\leq \mathbb{E}_{\boldsymbol{\theta}}\left[ \frac{f(\mathbf{x}; \boldsymbol{\theta}^*)}{f(\mathbf{x}; \boldsymbol{\theta})} \right] 
= 0 \\
\implies \mathbb{E}_{\boldsymbol{\theta}}\left[ \log(f(\mathbf{x}; \boldsymbol{\theta}^*))\right]  - \mathbb{E}_{\boldsymbol{\theta}}\left[\log(f(\mathbf{x}; \boldsymbol{\theta})) \right] &\leq 0 \\
\implies \mathbb{E}_{\boldsymbol{\theta}}\left[ \log(f(\mathbf{x}; \boldsymbol{\theta}^*))\right] &\leq \mathbb{E}_{\boldsymbol{\theta}}\left[\log(f(\mathbf{x}; \boldsymbol{\theta})) \right]
\end{aligned}
$$

<aside><p>I think $\boldsymbol{\theta}^*$ is an arbitrary value of $\boldsymbol{\theta}$ and $\boldsymbol{\theta}$ is the true value (that is not necessarily equal to $\boldsymbol{\theta}_0$)?</p></aside>

An application of the law of large numbers yields:

$$
\begin{aligned}
\underset{\boldsymbol{\theta}}{\text{plim}} \left[\frac{1}{n} \ell(\boldsymbol{\theta}^*; \mathbf{x}) \right]
  &= \mathbb{E}_{\boldsymbol{\theta}}\left[ \log(f(\mathbf{x}; \boldsymbol{\theta}^*)) \right] 
  \leq \mathbb{E}_{\boldsymbol{\theta}}\left[\log(f(\mathbf{x}; \boldsymbol{\theta})) \right] 
  = \underset{\boldsymbol{\theta}}{\text{plim}} \left[\frac{1}{n} \ell(\boldsymbol{\theta}; \mathbf{x}) \right]
\end{aligned}
$$

The above inequality implies that the probability limit of $$\frac{1}{n} \ell(\boldsymbol{\theta}^*; \mathbf{x})$$ will be maximized by $\boldsymbol{\theta}$ (i.e. at $\boldsymbol{\theta}^* = \boldsymbol{\theta}$). If we assume the likelihood function is nice, then this point is also a stationary point. Thus, if the slope of the log-likelihood function is far from zero at $\boldsymbol{\theta}_0$, then this provides evidence against $H_0$. 

The following hold:

$$
\begin{align}
&\underset{\boldsymbol{\theta}_0}{\text{plim}} \left[ \frac{1}{n} \mathbf{S}(\boldsymbol{\theta}_0) \right] = \mathbf{0}_p & (\text{under } H_0 )\\
&\underset{\boldsymbol{\theta} \neq \boldsymbol{\theta}_0}{\text{plim}} \left[ \frac{1}{n} \mathbf{S}(\boldsymbol{\theta}_0) \right] = \mathbb{E}_{\boldsymbol{\theta}}\left[ \frac{\partial \ell(\boldsymbol{\theta}_0; \mathbf{x})}{\partial \boldsymbol{\theta}} \right] & ( \text{under } H_1 )\\
&\frac{1}{\sqrt{n}}\mathbf{S}(\boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \mathcal{I}(\boldsymbol{\theta}_0)) & \text{ under } H_0
\end{align}
$$

<aside><p>Include proof?</p></aside>

The <strong>likelihood ratio statistic</strong>, <strong>score statistic</strong>, and <strong>Wald statistic</strong> for testing Eq. \eqref{eq:hypotheses-0} are, respectively:

$$
\begin{align}
\text{LRT} &= -2 \left[ \ell(\hat{\boldsymbol{\theta}}_0; \mathbf{x}) - \ell(\hat{\boldsymbol{\theta}}_1; \mathbf{x}) \right] \\
\text{S}_L &= \frac{1}{n} \mathbf{S}^\top(\boldsymbol{\theta}_0) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0) \\
\text{T}_W &= n (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0)^\top \hat{\mathcal{I}} (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0)
\end{align}
$$

Let $$\boldsymbol{\theta}_n = \left(\boldsymbol{\theta}_0 + \frac{1}{\sqrt{n}} \boldsymbol{\delta}\right)$$ be a sequence of parameter values for some fixed $\boldsymbol{\delta}$ and fixed test level. We will let $H_n: \boldsymbol{\theta} = \boldsymbol{\theta}_n$ be a corresponding sequence of local alternatives. It can be shown that:

$$
\begin{aligned}
\text{LRT} &= \text{S}_L + o_p(1) = \text{T}_W + o_p(1) \hspace{4mm} \text{under } H_0 \text{ and } H_n\\
\end{aligned}
$$

This implies that all three test statistics have an asymptotic $\chi^2(\boldsymbol{\Delta})$ distribution under $H_n$, where $\boldsymbol{\Delta} = \boldsymbol{\delta}^\top \mathcal{I}(\boldsymbol{\theta}_0) \boldsymbol{\delta}$. 

### Nuisance Parameters
We will now assume that $\boldsymbol{\theta}$ can be partitioned into the $k$-dimensional nuisance and $(p - k)$-dimensional structural parameters:

$$
\boldsymbol{\theta} = (\boldsymbol{\lambda}^\top, \boldsymbol{\psi}^\top)^\top =: (\boldsymbol{\lambda} : \boldsymbol{\psi})
$$

We will also partition the score function, information matrix, and inverse information matrix conformably with $\boldsymbol{\theta}$. Let:

$$
(\mathcal{I}^{-1}(\boldsymbol{\theta}))_{\boldsymbol{\psi}, \boldsymbol{\psi}} = \mathcal{I}^{-1}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}) = \left[ \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\psi}}(\boldsymbol{\theta}) - \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathcal{I}^{-1}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathcal{I}_{\boldsymbol{\lambda}, \boldsymbol{\psi}}(\boldsymbol{\theta}) \right]^{-1}
$$

Consider testing:

$$
\begin{equation}
\label{eq:hypotheses-0b}
H_0: \boldsymbol{\psi} = \boldsymbol{\psi}_0
\hspace{5mm} \text{vs} \hspace{5mm}
H_1: \boldsymbol{\psi} \neq \boldsymbol{\psi}_0
\end{equation}
$$

Denote the true value of the nuisance parameters with $\boldsymbol{\lambda}_0$, and let $\boldsymbol{\theta}_0 = (\boldsymbol{\lambda}_0 : \boldsymbol{\psi}_0)$ be the true value of $\boldsymbol{\theta}$ under $H_0$. Let $\hat{\boldsymbol{\lambda}}_0$ and $\hat{\boldsymbol{\theta}}_0 = (\hat{\boldsymbol{\lambda}}_0: \boldsymbol{\psi}_0)$ denote the MLE of $\boldsymbol{\lambda}$ and $\boldsymbol{\theta}$ under $H_0$, respectively. Let $\hat{\boldsymbol{\theta}}_1 = (\hat{\boldsymbol{\lambda}}_1 : \hat{\boldsymbol{\psi}}_1)$ denote the same under $H_1$. 

The <strong>likeihood ratio statistic</strong>, <strong>score statistic</strong>, and <strong>Wald statistic</strong> for testing Eq. \eqref{eq:hypotheses-0b} are, respectively:

$$
\begin{align}
\text{LRT} &= -2 \left[ \ell(\hat{\boldsymbol{\theta}}_0; \mathbf{x}) - \ell(\hat{\boldsymbol{\theta}}_1; \mathbf{x}) \right]\\
\text{S}_L &= \frac{1}{n} \mathbf{S}^\top_{\boldsymbol{\psi}}(\hat{\boldsymbol{\theta}}_0) \mathcal{I}^{-1}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\hat{\boldsymbol{\theta}}_0) \mathbf{S}_{\boldsymbol{\psi}}(\hat{\boldsymbol{\theta}}_0) \\
\text{T}_W &= n (\hat{\boldsymbol{\psi}}_1 - \boldsymbol{\psi}_0)^\top \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\hat{\boldsymbol{\theta}}_1) (\hat{\boldsymbol{\psi}}_1 - \boldsymbol{\psi}_0)
\end{align}
$$

For some fixed $\boldsymbol{\delta} \in \mathbb{R}^k$, we define the sequence of local alternatives:

$$
H_n: \boldsymbol{\psi}_n = \boldsymbol{\psi}_0 + \frac{1}{\sqrt{n}}\boldsymbol{\delta}
$$

It can be shown that:

<aside><p>How?</p></aside>

$$
\begin{align}
&\sqrt{n}(\hat{\boldsymbol{\psi}} - \boldsymbol{\psi}_0) \rightsquigarrow \mathcal{N}(\boldsymbol{\delta}, \mathcal{I}^{-1}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0)) & \text{under } H_n \\
&\frac{1}{\sqrt{n}} \mathbf{S}_{\boldsymbol{\psi}}(\hat{\boldsymbol{\theta}}_0) \rightsquigarrow \mathcal{N}\left(\mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0) \boldsymbol{\delta}, \mathcal{I}_{\boldsymbol{\psi},\boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0) \right) & \text{under } H_n
\end{align}
$$

It can also be shown that:

$$
\begin{aligned}
\text{LRT} &= \text{S}_L + o_p(1) = \text{T}_W + o_p(1) \hspace{4mm} \text{under } H_0 \text{ and } H_n\\
\end{aligned}
$$

and that all three statistics are asymptotically $\chi_k^2(\boldsymbol{\Delta})$ under $H_n$ where $$\boldsymbol{\Delta} = \boldsymbol{\delta}^\top \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0) \boldsymbol{\delta}$$. 

### Extension to Estimating Equations
The above results can be extended to cases where we do not use the score function as $\mathbf{S}(\boldsymbol{\theta})$ in our test statistic. Instead, we let $\mathbf{S}(\boldsymbol{\theta}) = \mathbf{0}_p$ be an <strong>estimating equation</strong> where $\mathbf{S}(\boldsymbol{\theta})$ is our <strong>estimating function</strong>.

To make sure our estimating function is chosen suitably, we must make sure it satisfies a couple of conditions. We assume that there exist non-singular $\mathbf{G}(\boldsymbol{\theta})$ and $\mathbf{V}(\boldsymbol{\theta})$ satisfying for $a > 0$:

<ol>
<li>$\frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \mathbf{V}(\boldsymbol{\theta}_0))$</li>
<li>$\underset{\rvert \rvert \mathbf{u} \rvert \rvert \leq a}{\sup} \left\{ \left\rvert \left\rvert \frac{1}{\sqrt{n}} \left[ \mathbf{S}\left(\boldsymbol{\theta}_0 + \frac{1}{\sqrt{n}} \mathbf{u}\right) - \mathbf{S}(\boldsymbol{\theta}_0)\right] + \mathbf{G}(\boldsymbol{\theta}_0)\mathbf{u} \right\rvert \right\rvert \right\} = o_p(1)$</li>
</ol>

<aside><p><strong>Condition A.</strong></p></aside>

uniformly in $\boldsymbol{\theta}_0$ as $n \rightarrow \infty$.

Suppose we partition $\mathbf{G}(\boldsymbol{\theta})$ and $\mathbf{V}(\boldsymbol{\theta})$ to conform with the partitioning of $\boldsymbol{\theta} = (\boldsymbol{\lambda} : \boldsymbol{\psi})$. If the above conditions hold, then:

$$
\sqrt{n}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}\left(\mathbf{0}_p, \mathbf{G}^{-1}(\boldsymbol{\theta}_0) \mathbf{V}(\boldsymbol{\theta}_0) (\mathbf{G}^{-1}(\boldsymbol{\theta}_0))^\top\right)
$$

where $\hat{\boldsymbol{\theta}}$ is the solution to $\mathbf{S}(\hat{\boldsymbol{\theta}}) = \mathbf{0}_p$, and $\boldsymbol{\theta}_0$ is the true value of the parameter vector under $H_0$. 


<div id="prop-4-5-1"></div>
<div class="theorem">
  <strong>Proposition 4.5.1.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  {% tabs p-4-5-1 %}
  {% tab p-4-5-1 statement %}
  Define:

  $$
  \begin{aligned}
  \mathbf{C}(\boldsymbol{\theta})
    &= [\mathbf{V}_{\boldsymbol{\psi}, \boldsymbol{\psi}}(\boldsymbol{\theta}) - \mathbf{G}_{\boldsymbol{\psi}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathbf{G}^{-1}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathbf{V}_{\boldsymbol{\lambda}, \boldsymbol{\psi}}(\boldsymbol{\theta})] - [ \mathbf{V}^\top_{\boldsymbol{\lambda}, \boldsymbol{\psi}}(\boldsymbol{\theta}) - \mathbf{G}_{\boldsymbol{\psi}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathbf{G}^{-1}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathbf{V}_{\boldsymbol{\lambda},\boldsymbol{\lambda}}(\boldsymbol{\theta})](\mathbf{G}^{-1}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}}(\boldsymbol{\theta}))^\top \mathbf{G}^\top_{\boldsymbol{\psi}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \\
  \mathbf{Z}(\boldsymbol{\theta})
    &= \frac{1}{\sqrt{n}} \left[ \mathbf{S}_{\boldsymbol{\psi}}(\boldsymbol{\theta}) - \mathbf{G}_{\boldsymbol{\psi}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathbf{G}^{-1}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}}(\boldsymbol{\theta}) \mathbf{S}_{\boldsymbol{\lambda}}(\boldsymbol{\theta}) \right]
  \end{aligned}
  $$

  Let $$\hat{\boldsymbol{\lambda}}_0$$ be the estimator of $\boldsymbol{\lambda}$ under $H_0$ that satisfies $$\mathbf{S}_\boldsymbol{\lambda}( (\hat{\boldsymbol{\lambda}}_0 : \boldsymbol{\psi}_0) ) = \mathbf{0}_p$$. Let $$\hat{\boldsymbol{\theta}}_0 = (\hat{\boldsymbol{\theta}}_0 ; \boldsymbol{\psi})$$. Also let $$H_n : \boldsymbol{\psi} = \boldsymbol{\psi}_0 + \frac{1}{\sqrt{n}} \boldsymbol{\delta}$$ define a squence of local alternatives for some fixed $\boldsymbol{\delta} \in \mathbb{R}^k$. Then:

  <ol>
  <li>$\mathbf{z}(\hat{\boldsymbol{\theta}}_0) = \mathbf{z}(\boldsymbol{\theta}_0) + o_p(1)$ under $H_0$</li>
  <li>$\mathbf{z}(\boldsymbol{\theta}_0)$ and $\mathbf{z}(\hat{\boldsymbol{\theta}}_0) \rightsquigarrow \mathcal{N}\left(\mathbf{G}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0) \boldsymbol{\delta}, \mathbf{C}(\boldsymbol{\theta}_0)\right)$ under $H_n$</li>
  </ol>

  And a <strong>local score-type statistic</strong> for testing Eq. \eqref{eq:hypotheses-0b} is:

  $$
  \text{S}_L = \mathbf{z}^\top(\hat{\boldsymbol{\theta}}_0) \mathbf{C}^{-1}(\hat{\boldsymbol{\theta}}_0) \mathbf{z}(\hat{\boldsymbol{\theta}}_0)
  $$ 

  which is asymptotically $$\chi_k^2(\boldsymbol{\Delta})$$ under $H_n$ where $$\boldsymbol{\Delta} = \boldsymbol{\delta}^\top \mathbf{G}^\top_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0) \mathbf{C}^{-1}(\boldsymbol{\theta}_0) \mathbf{G}_{\boldsymbol{\psi}, \boldsymbol{\psi} \cdot \boldsymbol{\lambda}}(\boldsymbol{\theta}_0) \boldsymbol{\delta}$$.
  {% endtab %}
  {% tab p-4-5-1 proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

We can also define a Wald-type statistic:

$$
\text{T}_W = n (\hat{\boldsymbol{\psi}} - \boldsymbol{\psi}_0)^\top (\hat{\mathbf{A}}_{\boldsymbol{\psi}, \boldsymbol{\psi}}(\boldsymbol{\theta}))^{-1}(\hat{\boldsymbol{\psi}} - \boldsymbol{\psi}_0);
\hspace{5mm} \mathbf{A}(\boldsymbol{\theta}) = \mathbf{G}^{-1}(\boldsymbol{\theta})\mathbf{V}(\boldsymbol{\theta}) \mathbf{G}^\top(\boldsymbol{\theta})
$$

where $\hat{\mathbf{A}}(\boldsymbol{\theta})$ is a consistent estimator of the Fisher information. Under the sequence of local alternatives given in <a href="prop-4-5-1">Proposition 4.5.1</a>, $\text{T}_W$ also is asymptotically $\chi_k^2(\boldsymbol{\Delta})$.


---

## Linear Constraints
In this section, we assume that $\boldsymbol{\Theta}$ is open and consider testing:

$$
\begin{equation}
\label{eq:hypotheses-1}
H_0: \mathbf{R} \boldsymbol{\theta} = \mathbf{0}_p
\hspace{5mm} \text{vs} \hspace{5mm}
H_1: \mathbf{R} \boldsymbol{\theta} \geq \mathbf{0}_p, \mathbf{R} \boldsymbol{\theta} \neq \mathbf{0}_p
\end{equation}
$$

where $\mathbf{R}$ is a full rank, fixed, $r \times p$ matrix. Note that the null hypothesis states that $\boldsymbol{\theta}$ lies in the null space of $\mathbf{R}$ (equivalently, the orthogonal complement of the row space). 

In what follows, we will also assume the set-up given previously and that Condition Q is satisfied. We will also let $\boldsymbol{\Theta}_0$ denote the part of the parameter space containing values of $\boldsymbol{\theta}$ that are valid under $H_0$ and $\boldsymbol{\Theta}_1$ denotes the same for $H_1$.

Under Condition Q, all of the following test statistics are asymptotically equivalent for testing Eq. \eqref{eq:hypotheses-1}:

$$
\text{LRT} = \text{W} + o_p(1) = \text{S}_G + o_p(1) = \text{D} + o_p(1)
$$

Their asymptotic null distributions are also the same (given in Eq. \eqref{eq:lrt-dist-asymp}).

### Likelihood Ratio Statistic
Suppose $H_0$ is true, and let $\mathbf{z} \in \mathbb{R}^p$. Define the quantity:

$$
D(\mathbf{z}) = \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\min} \left\{ (\mathbf{z} - \boldsymbol{\theta})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z} - \boldsymbol{\theta}) \right\} - \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_1}{\min} \left\{ (\mathbf{z} - \boldsymbol{\theta})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z} - \boldsymbol{\theta}) \right\} 
$$

Consider the quadratic approximation given in Eq. \eqref{eq:quad-approx-2}. For $\mathbf{u} = \sqrt{n}(\boldsymbol{\theta} - \boldsymbol{\theta}_0)$, $K>0$, and:

$$
\mathbf{z}_n = \frac{1}{\sqrt{n}} \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0);
\hspace{5mm}
\underset{\rvert \rvert \mathbf{u} \rvert \rvert < K}{\sup} \left\{ \delta_n(\mathbf{u}) \right\} = o_p(1)
$$

we can write:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{z})
  &= \underbrace{\ell(\boldsymbol{\theta}_0; \mathbf{z}) + \frac{1}{2 n} \mathbf{S}^\top(\boldsymbol{\theta}_0) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0)}_{=: A_n} - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u}) \\
  &= A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u})
\end{aligned}
$$

where $A_n$ is independent of $\mathbf{u}$. Now, we have:

$$
\begin{aligned}
\mathbf{R}\boldsymbol{\theta}_0 &= \mathbf{R} \mathbf{0}_p = \mathbf{0}_p \\
\implies \mathbf{R} \mathbf{u} &= \mathbf{R} \left(\sqrt{n} (\boldsymbol{\theta} - \boldsymbol{\theta}_0) \right) = \sqrt{n} \mathbf{R} \boldsymbol{\theta} 
\end{aligned}
$$

This implies that $$\mathbf{R} \boldsymbol{\theta} \geq \mathbf{0}_p$$ if, and only if, $$\mathbf{R}\mathbf{u} \geq \mathbf{0}_p$$ and $$\mathbf{R} \boldsymbol{\theta} = \mathbf{0}_p$$ if, and only if, $$\mathbf{R} \mathbf{u} = \mathbf{0}_p$$. Then:

$$
\text{LRT} = D(\mathbf{z}_n) + o_p(1)
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\text{LRT} &= -2 \left[ \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{z}) \right\} - \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_1}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{z}) \right\}  \right] \\
&= 2 \left[ \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_1}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{z}) \right\} - \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{z}) \right\}  \right] \\
&= 2 \left[ \underset{\mathbf{u} : \mathbf{R} \mathbf{u} \geq \mathbf{0}_p}{\sup} \left\{ A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u}) \right\} - \underset{\mathbf{u} : \mathbf{R} \mathbf{u} =\mathbf{0}_p}{\sup} \left\{ A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u}) \right\}  \right] \\
&= 2 A_n - 2 A_n - \underset{\mathbf{u} : \mathbf{R} \mathbf{u} \geq \mathbf{0}_p}{\inf} \left\{ A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u}) \right\} + \underset{\mathbf{u} : \mathbf{R} \mathbf{u} =\mathbf{0}_p}{\inf} \left\{ A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u}) \right\} \\
&= \underset{\mathbf{u} : \mathbf{R} \mathbf{u} =\mathbf{0}_p}{\sup} \left\{ A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) \right\} - \underset{\mathbf{u} : \mathbf{R} \mathbf{u} \geq \mathbf{0}_p}{\sup} \left\{ A_n - \frac{1}{2} (\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0)(\mathbf{z}_n - \mathbf{u}) \right\} + o_p(1) \\
&= D(\mathbf{z}_n) + o_p(1)
\end{aligned}
$$
</details>

Since $\mathbf{z}_n \rightsquigarrow \mathbf{x}$ and $D(\mathbf{z})$ is continuous in $\mathbf{z}$, $D(\mathbf{z}_n)$ (and therefore the LRT statistic) have the same asymptotic distribution as $D(\mathbf{z})$. For $c > 0$, this distribution is:

$$
\begin{equation}
\label{eq:lrt-dist-asymp}
\underset{n \rightarrow \infty}{\lim} \mathbb{P}_{\boldsymbol{\theta}_0} (\text{LRT} \geq c \rvert H_0) = \sum_{i = 1}^r \omega_i(r, \mathbf{R} \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{R}^\top, \mathcal{O}^q) \mathbb{P}(\chi_i^2 \geq c)
\end{equation}
$$

where $\mathcal{O}^q$ is the non-negative orthant in $q$ dimensions. This statistic requires estimating the MLE under $H_1$. 

### Wald-Type Statistic
We can define a Wald statistic for Eq. \eqref{eq:hypotheses-1}:

$$
\text{W} = n (\mathbf{R} \hat{\boldsymbol{\theta}}_1)^\top (\mathbf{R} \hat{\mathcal{I}}_1^{-1} \mathbf{R}^\top)^{-1}(\mathbf{R} \hat{\boldsymbol{\theta}}_1)
$$

where $\hat{\boldsymbol{\theta}}_1$ is the MLE of $\boldsymbol{\theta}$ under $H_1$, and $\hat{\mathcal{I}}_1$ is an estimator of the Fisher information under $H_1$ (e.g. $\mathcal{I}(\hat{\boldsymbol{\theta}}_1)$ or $\mathcal{I}(\hat{\boldsymbol{\theta}})$). 

One way to interpret this Wald statistic is that it kind of measures the length of $\mathbf{R} \hat{\boldsymbol{\theta}}_1$ (w.r.t the inner product induced by $\mathbf{R} \hat{\mathcal{I}}^{-1}_1 \mathbf{R}^\top$). Since we estimate $\hat{\boldsymbol{\theta}}_1$ under $H_1$, we would expect $$\rvert \rvert \mathbf{R} \hat{\boldsymbol{\theta}}_1 \rvert \rvert_{\mathbf{R} \hat{\mathcal{I}}^{-1}_1 \mathbf{R}^\top}^2$$ to be small when $H_0$ is true ($\mathbf{R} \boldsymbol{\theta} = \mathbf{0}_p$) and to be large when it is not true. This statistic requires estimating the MLE under $H_1$. 

### Distance Statistic
We can define a statistic that is based upon comparing the distance between the global MLE and the MLEs under $H_0$ and $H_1$. 

$$
\text{D} = \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\inf} \left\{ n (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta})^\top \hat{\mathcal{I}} (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}) \right\} - \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_1}{\inf}\left\{ n (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta})^\top \hat{\mathcal{I}} (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}) \right\}
$$

where $\hat{\mathcal{I}}$ is consistent estimator of the Fisher information. This statistic only requires computing the global MLE. 

### Score Statistic
Consider the quantity:

$$
\mathbf{S}(\hat{\boldsymbol{\theta}}_0) - \mathbf{S}(\hat{\boldsymbol{\theta}}_1)
$$

If $H_0$ is true, then we would expect $\hat{\boldsymbol{\theta}}_0$ and $\hat{\boldsymbol{\theta}}_1$ to be close, thus making the above quantity close to zero. However, if $H_0$ is not true, we would expect the $\mathbf{S}(\hat{\boldsymbol{\theta}}_1)$ to be close to zero (since we will aim to find a maximizer of the log-likelihood over $\boldsymbol{\Theta}_1$, i.e. a stationary point), but $\mathbf{S}(\hat{\boldsymbol{\theta}}_0)$ should be far from zero. This would make the quantity far from zero. Rather than directly use the score function, we use the <i>effective score</i>. A <strong>global</strong> score test statistic is given by:

<aside><p>This statistic is called "global" because we consider the shape over all of $\boldsymbol{\Theta}$.</p></aside>

$$
\text{S}_G = (\mathbf{R} \hat{\mathcal{I}}^{-1} [\mathbf{S}(\hat{\boldsymbol{\theta}}_0 - \mathbf{S}(\hat{\boldsymbol{\theta}}_1))])^\top (\mathbf{R} \hat{\mathcal{I}}^{-1} \mathbf{R}^\top)^{-1} (\mathbf{R} \hat{\mathcal{I}}^{-1}(\mathbf{S}(\hat{\boldsymbol{\theta}}_0) - \mathbf{S}(\hat{\boldsymbol{\theta}}_1)))
$$

where, as before, $\hat{\mathcal{I}}$ is a consistent estimator of the Fisher information. This statistic requires estimating the MLE under $H_1$. 

### Extension to Non-IID Case
Suppose that we have independent observations, but they are not identically distributed. In this case, let $f_i(\mathbf{x}_i; \boldsymbol{\theta})$ denote the density function for observation $i$. Let $\mathbf{X}$ denote the collection of all $n$ observations. We assume that <a href="#assumptions">Assumption A1</a> holds.

<aside><p>This setting includes testing coefficients in linear regression.</p></aside>

If there exists a positive definite matrix $\mathcal{V}(\boldsymbol{\theta})$ such that:

$$
\frac{1}{\sqrt{n}}\mathbf{S}(\boldsymbol{\theta}) \rightsquigarrow \mathcal{N}\left(\mathbf{0}_p, \mathcal{V}(\boldsymbol{\theta})\right);
\hspace{10mm}
\frac{1}{n} \frac{\partial^2 \ell(\boldsymbol{\theta}; \mathbf{X})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \overset{as}{\rightarrow} - \mathcal{V}(\boldsymbol{\theta})
$$

such that, for $$\mathbf{z}_n = \frac{1}{\sqrt{n}} \mathcal{V}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0)$$ and $\mathbf{u} = \sqrt{n}(\boldsymbol{\theta} - \boldsymbol{\theta}_0)$ with $$\underset{\rvert \rvert \mathbf{u} \rvert \rvert < K}{\sup} \left\{ \rvert \delta_n(\mathbf{u}) \rvert \right\} = o_p(1)$$ for some $K > 0$:

$$
\ell(\boldsymbol{\theta}; \mathbf{X})
  = \ell(\boldsymbol{\theta}_0; \mathbf{X}) + \frac{1}{2n} \mathbf{S}^\top(\boldsymbol{\theta}_0)\mathcal{V}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0) - \frac{1}{2}(\mathbf{z}_n - \mathbf{u})^\top \mathcal{I}(\boldsymbol{\theta}_0) (\mathbf{z}_n - \mathbf{u}) + \delta_n(\mathbf{u})
$$

and $\mathcal{V}(\boldsymbol{\theta})$ is continuous in $\boldsymbol{\theta}$, and the MLEs under $H_0$ and $H_1$ are $\sqrt{n}$-consistent, then the following test statistics for testing Eq. \eqref{eq:hypotheses-1} are asymptotically equivalent:

$$
\begin{align}
\text{LRT} &= -2 \left[ \ell(\hat{\boldsymbol{\theta}}_0; \mathbf{X}) - \ell(\hat{\boldsymbol{\theta}}_1; \mathbf{X}) \right] \\
\text{W} &=  n(\mathbf{R} \hat{\boldsymbol{\theta}}_1)^\top (\mathbf{R} \hat{\mathcal{V}}^{-1} \mathbf{R}^\top)^{-1} (\mathbf{R} \hat{\boldsymbol{\theta}}_1) \\
\text{S}_G &= \frac{1}{n} \left[ \mathbf{R} \hat{\mathcal{V}}^{-1} (\mathbf{S}(\hat{\boldsymbol{\theta}}_1) - \mathbf{S}(\hat{\boldsymbol{\theta}}_0)) \right]^\top (\mathbf{R} \hat{\mathcal{V}}^{-1} \mathbf{R}^\top)^{-1} \left[ \mathbf{R} \hat{\mathcal{V}}^{-1} (\mathbf{S}(\hat{\boldsymbol{\theta}}_1) - \mathbf{S}(\hat{\boldsymbol{\theta}}_0)) \right] \\
\text{D} &= \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\inf} \left\{ n (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta})^\top \hat{\mathcal{V}}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}) \right\} - \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\inf} \left\{ n (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta})^\top \hat{\mathcal{V}}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}) \right\} 
\end{align}
$$

<aside><p>I think $\hat{\mathcal{V}}$ is a consistent estimator fo $\mathcal{V}(\boldsymbol{\theta})$.</p></aside>

Furthermore:

$$
\text{LRT} = \text{W} + o_p(1) = \text{S}_G + o_p(1) = \text{D} + o_p(1)
$$

and as $n \rightarrow \infty$:

$$
\mathbb{P}_{\boldsymbol{\theta}_0}\left(\text{T} \geq c \rvert H_0 \right) \rightarrow \sum_{i = 0}^r \omega_i(r, \mathbf{R} \mathcal{V}^{-1}(\boldsymbol{\theta}_0) \mathbf{R}^\top, \mathcal{O}^p) \mathbb{P}(\chi_i^2 \geq c)
$$

<aside><p>See Proposition 4.3.4.<d-cite key=silvapulle2005></d-cite></p></aside>

where $\text{T}$ is any one of the test statistics. This result follows from the fact that the asymptotic distribution of $\text{LRT}$ under $H_0$ can be shown to be the same as the distribution of:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\inf} \left\{ (\mathbf{z} - \boldsymbol{\theta})^\top \mathcal{V}(\boldsymbol{\theta}_0)(\mathbf{z} - \boldsymbol{\theta}) \right\} -
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_1}{\inf} \left\{ (\mathbf{z} - \boldsymbol{\theta})^\top \mathcal{V}(\boldsymbol{\theta}_0)(\mathbf{z} - \boldsymbol{\theta}) \right\}
$$

with $$\mathbf{z} \sim \mathcal{N}(\mathbf{0}_p, \mathcal{V}^{-1}(\boldsymbol{\theta}_0))$$.

---

## Non-Linear Constraints

### Closed and Convex Cone
First, let's consider the i.i.d. setting introduced in previous sections, and let's test hypotheses of the form:
$$
\begin{equation}
\label{eq:hypotheses-2}
H_0: \boldsymbol{\theta} = \mathbf{0}_p
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\theta} \in \mathcal{C}
\end{equation}
$$

where $\mathcal{C}$ is a closed, convex cone in $\mathbb{R}^p$. We will assume that $\mathbf{0}_p \in \text{int}(\mathcal{C})$. Let us define the quantity:

$$
\mathbf{U} = \frac{1}{\sqrt{n}} \mathcal{I}^{-1}(\mathbf{0}_p)\mathbf{S}(\boldsymbol{0}_p)
$$

<aside><p>$\mathbf{S}(\boldsymbol{\theta}_0) = \mathbf{S}(\mathbf{0}_p)$.</p></aside>

where $\mathbf{S}(\boldsymbol{\theta})$ is the score function. 

<!-- 
Fix $\boldsymbol{\delta} \in \mathbb{R}^p$, and define $\boldsymbol{\theta}_n = \frac{1}{\sqrt{n}} \boldsymbol{\delta}$. Via a first order Taylor expansion of $\mathbf{S}(\boldsymbol{\theta})$ at $\boldsymbol{\theta}_0$ about $\boldsymbol{\theta}_n$, we can show:

$$
\begin{aligned}
\frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_0) 
  &= \frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_n) + \left[ \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x})}{\partial \boldsymbol{\theta}} \right]_{\boldsymbol{\theta} = \boldsymbol{\theta}_n} (\boldsymbol{\theta}_0 - \boldsymbol{\theta}_n) + o_p(1) \\
\implies \frac{1}{\sqrt{n}} \mathbf{S}(\mathbf{0}_p) 
  &= \frac{1}{\sqrt{n}} \mathbf{S}\left(\boldsymbol{\theta}_n \right) - \frac{1}{\sqrt{n}} \mathcal{I}\left( \boldsymbol{\theta}_n \right) \boldsymbol{\delta} + o_p(1)
\end{aligned}
$$

Since $\boldsymbol{\theta}_n \rightarrow \mathbf{0}_p$ as $n \rightarrow \infty$, we should have that $\mathcal{I}(\boldsymbol{\theta}_n) \boldsymbol{\delta} \overset{p}{\rightarrow} \mathcal{I}(\mathbf{0}_p) \boldsymbol{\delta}$.  -->

We can define a <strong>local score test statistic</strong> and a <strong>distance statistic</strong> for testing Eq. \eqref{eq:hypotheses-2}:

$$
\begin{align}
\text{LRT} &= -2 \left[ \ell(\hat{\boldsymbol{\theta}}_0; \mathbf{x}) - \ell(\hat{\boldsymbol{\theta}}_1; \mathbf{x}) \right] \\
\text{S}_L &= \mathbf{U}^\top \mathcal{I}(\mathbf{0}_p) \mathbf{U} - \underset{\mathbf{c} \in \mathcal{C}}{\min} \left\{ (\mathbf{U} - \mathbf{c})^\top \mathcal{I}(\mathbf{0}_p)(\mathbf{U} - \mathbf{c}) \right\} \\
\text{D} &= n \left[ \hat{\boldsymbol{\theta}}_1^\top \hat{\mathcal{I}}_0 \hat{\boldsymbol{\theta}}_1 - \underset{\mathbf{c} \in \mathcal{C}}{\min} \left\{ (\hat{\boldsymbol{\theta}}_1 - \mathbf{c})^\top \hat{\mathcal{I}}_0(\hat{\boldsymbol{\theta}}_1 - \mathbf{c}) \right\} \right]
\end{align}
$$

where $$\hat{\mathcal{I}}_0$$ is a consistent estimator of $$\mathcal{I}(\mathbf{0}_p)$$. It can be shown that:

$$
\text{LRT} = \text{S}_L + o_p(1) = \text{D} + o_p(1)
$$

<aside><p>See Proposition 4.6.1.<d-cite key=silvapulle2005></d-cite></p></aside>

under $H_0$, and their asymptotic null distributions are all $$\bar{\chi}^2(\mathcal{I}^{-1}(\mathbf{0}_p), \mathcal{C})$$. 



