---
layout: distill
title: Constrained Statistical Inference (Chapter 3)
description: Tests on Multivariate Normal Mean
date: 2026-07-13
tabs: true
tags: likelihood theory testing constraints
toc:
  - name: Background
    subsections:
        - name: Linear Algebra and Geometry
        - name: Other Stuff
  - name: Set-Up
  - name: Likelihood Ratio Test in Two Dimensions
  - name: Chi-Bar-Square Distribution
    subsections:
        - name: Geometric Interpretation
        - name: Results
        - name: Extensions
  - name: Unknown Covariance
bibliography: stats-ml.bib
---

In this post, I will be going through some of the content in Chapter 3 of Silvapulle and Sen's <i>Constrained Statistical Inference: Inequality, Order, and Shape Restrictions</i>.<d-cite key=silvapulle2005></d-cite>. The chapter concerns testing the mean of a multivariate Normal variate.

---

## Background

### Linear Algebra and Geometry
For now, we will concern ourselves with Euclidean space, $\mathbb{R}^p$.

We will let $\mathcal{C}_a$, $\mathcal{C}_b$, $\mathcal{C}$, and $\mathcal{M}$ to be subsets of $\mathbb{R}^p$ where $\mathcal{C}$ is a closed, <a href="/stats-ml/definitions#convex-set">convex</a> <a href = "/stats-ml/definitions#cone">cone</a>, $\mathcal{M}$ is a <a href="/stats-ml/definitions#vector-space">vector space</a>, $\mathcal{M} \subset \mathcal{C}$, and $\mathcal{C}_a \subset \mathcal{C}_b$. We will also define a <strong>polyhedral</strong> in $\mathbb{R}^p$ as a set:

$$
\mathcal{P} = \left\{ \boldsymbol{\theta} \in \mathbb{R}^p : \mathbf{a}_1^\top \boldsymbol{\theta} \geq 0, \dots, \mathbf{a}_k^\top \boldsymbol{\theta} \geq 0 \right\}
$$

for some $\mathbf{a}_1, \dots, \mathbf{a}_k \in \mathbb{R}^p$. 

Positive definite matrices induce inner products on vector spaces. Let $\mathbf{M}$ denote a positive definite matrix. The induced inner product, norm (called the <i>$\mathbf{M}$-norm</i>), and distance function (called <i>$\mathbf{M}$-distance</i>) are given by, for $\mathbf{x}, \mathbf{y} \in \mathbb{R}^p$:

$$
\begin{aligned}
\langle \mathbf{x}, \mathbf{y} \rangle_\mathbf{M} &:= \mathbf{x}^\top \mathbf{M}^{-1} \mathbf{y} \\
\rvert \rvert \mathbf{x} \rvert \rvert_\mathbf{M} &:= \sqrt{\mathbf{x}^\top \mathbf{M}^{-1} \mathbf{x}} \\
d_\mathbf{M}(\mathbf{x}, \mathbf{y}) &:= \rvert \rvert \mathbf{x} - \mathbf{y} \rvert \rvert_{\mathbf{M}} 
\end{aligned}
$$

We will denote vectors that are orthogonal with respect to the induced inner product as $\mathbf{x} \perp_{\mathbf{M}} \mathbf{y}$ and will call them <i>$\mathbf{M}$-orthogonal</i>. Similarly, the $\mathbf{M}$-projection of a vector $\mathbf{x}$ onto some suitable subset will be denoted by $\Pi_\mathbf{M}(\mathbf{x})$ and will be the point on the subset that is closest to $\mathbf{x}$ in terms of the $\mathbf{M}$-distance. 

We can define the orthogonal projection of $\mathbf{x} \in \mathbb{R}^p$ onto a vector subspace $\mathcal{W} \subset \mathbb{R}^p$ as:

$$
\begin{equation}
\label{eq:proj-subsp}
\Pi_{\mathbf{V}}(\mathbf{x} \rvert \mathcal{W}) = \underset{\mathbf{a} \in \mathcal{W}}{\arg \min} \left\{ (\mathbf{x} - \mathbf{a})^\top \mathbf{V}^{-1} (\mathbf{x} - \mathbf{a}) \right\}
\end{equation}
$$

Such a projection will exist and be unique. This implies that any $\mathbf{x} \in \mathbf{R}^p$ has a unique <i>orthogonal decomposition</i> given by:

$$
\mathbf{x} = \Pi_{\mathbf{V}}(\mathbf{x} \rvert \mathcal{W})  + \Pi_{\mathbf{V}}(\mathbf{x} \rvert \mathcal{W}^\perp) 
$$

where $\mathcal{W}^\perp$ is the <a href="/stats-ml/definitions/orthogonal-complement">orthogonal complement</a> of $\mathcal{W}$. 

Projections onto vector spaces satisfy several additional nice properties. If we have another vector space $\mathcal{M} \subset \mathcal{W}$, then projecting $\mathbf{x}$ onto $\mathcal{W}$ and then projecting that projection onto $\mathcal{M}$ will result in the same vector as if we had directly projected onto $\mathcal{M}$. That is:

$$
\begin{aligned}
\Pi_{\mathbf{V}}( \Pi_{\mathf{V}}(\mathbf{x} \rvert \mathcal{C}) \rvert \mathcal{M}) = \Pi_{\mathbf{V}}(\mathbf{x} \rvert \mathcal{C})
\end{aligned}
$$

A similar result holds for a closed convex cone $\mathcal{C}$ satisfying $\mathcal{M} \subset \mathcal{C} \subset \mathcal{W} \subset \mathbb{R}^p$. 

Similarly, for some closed, convex set $\mathbf{A} \in \mathbb{R}^p$ and some $\mathbf{x} \in \mathbb{R}^p$, we can define the <strong>projection</strong> of $\mathbf{x}$ onto $\mathbf{A}$ as the closest point on $\mathbf{A}$ to $\mathbf{x}$ in terms of $\mathbf{V}$-distance:

$$
\begin{equation}
\label{eq:proj-set}
\Pi_{\mathbf{V}}(\mathbf{x} \rvert \mathbf{A}) = \underset{\mathbf{a} \in \mathbf{A}}{\arg \min} \left\{ (\mathbf{x} - \mathbf{a})^\top \mathbf{V}^{-1} (\mathbf{x} - \mathbf{a}) \right\}
\end{equation}
$$

The following proposition provides some useful facts about the geometry of projections.

<div id="prop-3-12-3"></div>
<div class="theorem">
  <strong>Proposition 3.12.3<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  {% tabs prop-3-12-3 %}
  {% tab prop-3-12-3 statement %}
  Let $\mathbf{A}$ be a non-empty, closed, convex set in $\mathbb{R}^p$, and let $\mathbf{a} \in \mathbb{R}^p \setminus \mathbf{A}$. Then:

  $$
  \begin{equation}
  \label{eq:p-3-12-3-c}
  \exists \mathbf{p} \in \mathbb{R}^p \text{ and } \alpha \in \mathbb{R}  \text{ s.t. } \mathbf{p}^\top \mathbf{a} > \alpha \text{ and } \mathbf{p}^\top \mathbf{x} \leq \alpha \text{ for each } \mathbf{x} \in \mathbf{A}
  \end{equation}
  $$

  And:

  $$
  \begin{equation}
  \label{eq:p-3-12-3-a} 
  \exists \text{ unique } \mathbf{a}^* \in \mathbf{A} \text{ s.t. } \mathbf{a}^* = \underset{\mathbf{x} \in \mathbf{A}}{\arg \min} \left\{ \rvert \rvert \mathbf{a} - \mathbf{x} \rvert \rvert \right\} 
  \end{equation}
  $$

  Furthermore, $\mathbf{a}^* \in \mathbf{A}$ is the unique closest point in the previous equation if, and only if:
  
  $$
  \begin{equation}
  \label{eq:p-3-12-3-b}
  (\mathbf{a} - \mathbf{a}^*)^\top(\mathbf{x} - \mathbf{a}^*) \leq 0 
  \hspace{2mm} 
  \forall \mathbf{x} \in \mathbf{A}
  \end{equation}
  $$
  {% endtab %}
  {% tab prop-3-12-3 proof %}
  Proof to be included.
  {% endtab %}
  {% endtabs %}
</div>

There is also a unique orthogonal decomposition of vectors through their projection onto a cone. Let $\mathcal{C} \subset \mathbb{R}^p$, and let $\mathcal{C}^0$ denote its <a href="/stats-ml/definitions/polar-cone">polar cone</a>. We have:

$$
\begin{equation}
\label{eq:o-decomp-cone}
\mathbf{x} = \Pi_\mathbf{V}(\mathbf{x} \rvert \mathcal{C}) + \Pi_\mathbf{V}(\mathbf{x} \rvert \mathcal{C}^0)
\end{equation}
$$

<aside><p>Include proof?</p></aside>

where $$\Pi_\mathbf{V}^\top(\mathbf{x} \rvert \mathcal{C})\Pi_\mathbf{V}(\mathbf{x} \rvert \mathcal{C}^) = 0$$. 

### Other Stuff
Include proofs here.

---

## Set-Up
We will consider the setting of a $k$-dimensional multivariate Gaussian population with mean vector $\boldsymbol{\theta} = (\theta_1, \theta_2, \dots, \theta_k)^\top$ and <i>known</i> (and non-singular) covariance matrix, $\boldsymbol{\Sigma}$. We have $n$ independent samples, denoted by $\mathbf{x}_1, \dots, \mathbf{x}_n$, from this population.

We are concerned with what Silvapulle and Sen call <strong>Type A Tests</strong>, which are hypothesis testing problems on $\boldsymbol{\theta}$ of the form:

$$
\begin{equation}
\label{eq:hypotheses}
H_0: \boldsymbol{\theta} \in \mathcal{M}
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\theta} \in \mathcal{C} \setminus \mathcal{M}
\end{equation}
$$

In most other settings, we assume to have $n$ independent and identically distributed observations, $\mathbf{x}_1, \dots, \mathbf{x}_n$, with which we use to conduct our hypothesis test. However, because the sample mean is <a href="">sufficient</a> for $\boldsymbol{\theta}$, we can instead frame our discussions in terms of a single "observation", $\mathbf{x}$, that is basically a sample mean.

### Geometric Interpretation
We mostly will be concerned with the likelihood ratio test. Recall that the probability density function (PDF) for a multivariate Normal random vector is:

$$
\begin{aligned}
f(\mathbf{x}; \boldsymbol{\theta}, \boldsymbol{\Sigma}) &= (2 \pi)^{-\frac{k}{2}} \rvert \boldsymbol{\Sigma} \rvert^{-\frac{1}{2}} \exp\left(- \frac{1}{2} (\mathbf{x} - \boldsymbol{\theta})^\top \boldsymbol{\Sigma}^{-1}(\mathbf{x} - \boldsymbol{\theta}) \right) \\
&\propto \exp\left( - \frac{1}{2} (\mathbf{x} - \boldsymbol{\theta})^\top \boldsymbol{\Sigma}^{-1}(\mathbf{x} - \boldsymbol{\theta}) \right)
\end{aligned}
$$

<aside><p>The second line is called the <strong>kernel</strong>.</p></aside>

The log-likelihood is just the natural logarithm of the PDF reframed as a function of $\boldsymbol{\theta}$:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{x}) = - \frac{k}{2} \log(2 \pi) - \frac{1}{2} \log\left(\rvert \boldsymbol{\Sigma} \rvert\right) - \frac{1}{2} (\mathbf{x} - \boldsymbol{\theta})^\top \boldsymbol{\Sigma}^{-1}(\mathbf{x} - \boldsymbol{\theta}) 
\end{aligned}
$$

We can find the kernel of the log-likelihood by removing any constants and scaling factors. This yields:

$$
\begin{aligned}
K(\boldsymbol{\theta}) = (\mathbf{x} - \boldsymbol{\theta})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \boldsymbol{\theta}) 
\end{aligned}
$$

Now, since we assumed that $\boldsymbol{\Sigma}$ is positive definite, there exists some <i>invertible</i> matrix $\mathbf{A}$ such that $\mathbf{A}^\top \mathbf{A} = \boldsymbol{\Sigma}^{-1}$ (see <a href="https://en.wikipedia.org/wiki/Definite_matrix#Decomposition">here</a>). Define:

$$
\begin{equation}
\label{eq:transform-unit-var}
\mathbf{y} = \mathbf{A} \mathbf{x}
\iff
\mathbf{x} = \mathbf{A}^{-1} \mathbf{y}
\end{equation}
$$

As a linear transformation of a Gaussian random vector, $\mathbf{y}$ is also Gaussian. We also see that:

$$
\begin{aligned}
\boldsymbol{\mu} &:=  \mathbb{E}[\mathbf{y}] = \mathbf{A} \boldsymbol{\theta} \\
\mathbf{V} &:= \text{Var}(\mathbf{y}) = \mathbf{A}\boldsymbol{\Sigma} \mathbf{A}^\top = \mathbb{I}_{k \times k}
\end{aligned}
$$

<aside><p>The coordinates of $\mathbf{y}$ are uncorrelated and have unit variance!</p></aside>

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\mathbb{E}\left[ \mathbf{y} \right] 
    &= \mathbb{E}\left[ \mathbf{A} \mathbf{x} \right] \\
    &= \mathbf{A} \mathbb{E}\left[\mathbf{x} \right] \\
    &= \mathbf{A} \boldsymbol{\theta}  \\
\text{Var}(\mathbf{y}) 
    &= \text{Var}(\mathbf{A} \mathbf{x} ) \\
    &= \mathbf{A} \text{Var}(\mathbf{x}) \mathbf{A}^\top \\
    &= \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top
\end{aligned}
$$
Multiplying by $\mathbb{I}_{k \times k} = (\mathbf{A} \mathbf{A}^{-1})^\top$, we see:
$$
\begin{aligned}
\text{Var}(\mathbf{y})
    &= \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top \\
    &= (\mathbf{A} \mathbf{A}^{-1})^\top \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top \\
    &= (\mathbf{A}^{-1})^\top \mathbf{A}^\top \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top \\
    &= (\mathbf{A}^\top)^{-1} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Sigma} \mathbf{A}^\top \\
    &= (\mathbf{A}^\top)^{-1} \mathbf{A}^\top \\
    &= \mathbb{I}_{k \times k}
\end{aligned}
$$
</details>

We can then rewrite the kernel of the log-likelihood as:

$$
\begin{equation}
\label{eq:rewrite-kernel}
K(\boldsymbol{\theta}) 
    = K(\mathbf{A}^{-1} \boldsymbol{\mu})
    = (\mathbf{y} - \boldsymbol{\mu})^\top (\mathbf{y} - \boldsymbol{\mu}) \\
    = \rvert \rvert \mathbf{y} - \boldsymbol{\mu} \rvert \rvert_2^2
\end{equation}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
K(\boldsymbol{\theta}) 
    &= K(\mathbf{A}^{-1} \boldsymbol{\mu}) \\
    &= (\mathbf{A}^{-1} \mathbf{y} - \mathbf{A}^{-1} \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{A}^{-1}\mathbf{y} - \mathbf{A}^{-1} \boldsymbol{\mu}) \\
    &= (\mathbf{A}^{-1} (\mathbf{y} - \boldsymbol{\mu}))^\top \boldsymbol{\Sigma}^{-1} (\mathbf{A}^{-1} (\mathbf{y} - \boldsymbol{\mu})) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top(\mathbf{A}^{-1})^\top \boldsymbol{\Sigma}^{-1} \mathbf{A}^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top(\mathbf{A}^\top)^{-1} \boldsymbol{\Sigma}^{-1} \mathbf{A}^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top \left[\mathbf{A} \boldsymbol{\Sigma}\mathbf{A}^\top \right]^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top \left[(\mathbf{A}^\top)^{-1} \mathbf{A}^\top \mathbf{A} \boldsymbol{\Sigma}\mathbf{A}^\top \right]^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top \left[(\mathbf{A}^\top)^{-1} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Sigma}\mathbf{A}^\top \right]^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top \left[(\mathbf{A}^\top)^{-1} \mathbf{A}^\top \right]^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top \mathbf{V}^{-1} (\mathbf{y} - \boldsymbol{\mu}) \\
    &= (\mathbf{y} - \boldsymbol{\mu})^\top(\mathbf{y} - \boldsymbol{\mu}) \\
    &= \rvert \rvert \mathbf{y} - \boldsymbol{\mu} \rvert \rvert_2^2
\end{aligned}
$$
</details>

Thus, the kernel of the log-likelihood can be thought of as the (Euclidean) distance between the transformed data, $\mathbf{y}$, and the transformed mean vector, $\boldsymbol{\mu}$. 

<aside><p>Add geometric interpretation?</p></aside>

---

## Likelihood Ratio Test in Two Dimensions

### Set-Up
Suppose we have $\mathbf{x} = (x_1, x_2)^\top$ with $\mathbf{x} \sim \mathcal{N}(\boldsymbol{\theta}, \mathbb{I}_{2 \times 2})$. We wish to perform a likelihood ratio test of:

$$
\begin{equation}
\label{eq:lrt-ex-1}
H_0: \boldsymbol{\theta} = \mathbf{0}_2 
\hspace{5mm}
\text{vs.}
\hspace{5mm}
H_1: \boldsymbol{\theta} \in \mathcal{C} 
\end{equation}
$$

where $$\mathcal{C} = \left\{ \mathbf{c} \in \mathbb{R}^2 : \mathbf{R} \mathbf{c} \geq \mathbf{0}_2, \hspace{2mm} \mathbf{R} \mathbf{c} \neq \mathbf{0}_2 \right\}$$ for some non-singular, $2 \times 2$ matrix $\mathbf{R}$. $\mathcal{C}$ is clearly a closed, convex cone, and thus we can consider its <a href="/stats-ml/definitions#polar-cone">polar cone</a>, $\mathcal{C}^0$. Since we are in $\mathbb{R}^2$, we will let $\mathbf{l}$ and $\mathbf{u}$ denote unit vectors along the lower and upper boundaries of $\mathcal{C}$, and we will let $\mathbf{l}^0$ and $\mathbf{u}^0$ denote the same for $\mathcal{C}^0$. 

### Test Statistic
By the definition of the polar cone, we know that $\mathbf{l} \perp \mathbf{l}^0$ and $\mathbf{u} \perp \mathbf{u}^0$. Thus, these four vectors partition the space into four cones, $\mathcal{S}_1, \mathcal{S}_2, \mathcal{S}_3, \mathcal{S}_4$, which are analogous to the four quadrants of the plane. WLOG, we will assume $\mathcal{C} = \mathcal{S}_1$, $\mathcal{C}^0 = \mathcal{S}_3$, and the cones between $\mathbf{u}$ and $\mathbf{u}^0$ and between $\mathbf{l}$ and $\mathbf{l}^0$ are $\mathcal{S}_2$ and $\mathcal{S}_4$, respectively.

As we saw in Eq. \eqref{eq:rewrite-kernel}, the kernel of the log-likelihood function in this case is:

$$
K(\boldsymbol{\theta}) = \rvert \rvert \mathbf{x} - \boldsymbol{\theta} \rvert \rvert_2^2
$$

because the covariance matrix is the identity matrix. The maximum likelihood estimate (MLE) of $\boldsymbol{\theta}$, which we will denote with $\hat{\boldsymbol{\theta}}$, will then be the value of $\boldsymbol{\theta}$ that <i>minimizes</i> $K(\boldsymbol{\theta})$ (because we multiplied by $-2$ to obtain the kernel). When we consider the constraint imposed by $H_1$, we see that the value that minimizes $K(\boldsymbol{\theta})$ will be the closest point to $\mathbf{x}$ that is on $\mathcal{C}$ (i.e. the <i>orthogonal projection</i> of $\mathbf{x}$ onto $\mathcal{C}$). Let $\Pi_{\mathcal{C}}(\mathbf{x})$ denote this projection.

We now consider the value of $\hat{\boldsymbol{\theta}}_1$, the MLE under $H_1$, in cases.

#### Case 1: $$\mathbf{x} \in \mathcal{S}_1$$
If $\mathbf{x} \in \mathcal{S}_1$, then we need not do any projection. Thus:

$$
\hat{\boldsymbol{\theta}}_1 = \mathbf{x}
$$

#### Case 2: $$\mathbf{x} \in \mathcal{S}_2$$
If $\mathbf{x} \in \mathcal{S}_2$, then the projection of $\mathcal{C}$ is given by the projection of $\mathbf{x}$ onto $\mathbf{u}$: 

$$
\hat{\boldsymbol{\theta}}_1 = \frac{\langle \mathbf{x}, \mathbf{u} \rangle}{\rvert \rvert \mathbf{u} \rvert \rvert_2^2} \mathbf{u} = \langle \mathbf{x}, \mathbf{u} \rangle \mathbf{u}
$$


#### Case 3: $$\mathbf{x} \in \mathcal{S}_2$$
If $\mathbf{x} \in \mathcal{S}_3$, then the projection will be the origin:

$$
\hat{\boldsymbol{\theta}}_3 = \mathbf{0}_2
$$


#### Case 4: $$\mathbf{x} \in \mathcal{S}_4$$
Similar to Case 2, if $\mathbf{x} \in \mathcal{S}_4$, then the projection of $\mathcal{C}$ is given by the projection of $\mathbf{x}$ onto $\mathbf{l}$: 

$$
\hat{\boldsymbol{\theta}}_1 = \frac{\langle \mathbf{x}, \mathbf{l} \rangle}{\rvert \rvert \mathbf{l} \rvert \rvert_2^2} \mathbf{l} = \langle \mathbf{x}, \mathbf{l} \rangle \mathbf{l}
$$

In all cases, we have that $$(\mathbf{x} - \Pi_{\mathcal{C}}(\mathbf{x})) \perp \Pi_{\mathcal{C}}(\mathbf{x})$$:

$$
\begin{aligned}
\langle \mathbf{x} - \Pi_{\mathcal{C}}(\mathbf{x}),\Pi_{\mathcal{C}}(\mathbf{x}) \rangle 
    &= \begin{cases}
    (\mathbf{x} - \mathbf{x})^\top \mathbf{x} & \mathbf{x} \in \mathcal{S}_1 \\
    (\mathbf{x} - \langle \mathbf{x}, \mathbf{u} \rangle \mathbf{u})^\top (\langle \mathbf{x}, \mathbf{u} \rangle \mathbf{u}) & \mathbf{x} \in \mathcal{S}_2 \\
    (\mathbf{x} - \mathbf{0}_2)^\top \mathbf{0}_2 & \mathbf{x} \in \mathcal{S}_3 \\
    (\mathbf{x} - \langle \mathbf{x}, \mathbf{l} \rangle \mathbf{l})^\top (\langle \mathbf{x}, \mathbf{l} \rangle \mathbf{l}) & \mathbf{x} \in \mathcal{S}_4
    \end{cases} \\
    &= \begin{cases}
    0 & \mathbf{x} \in \mathcal{S}_1 \\
    (\mathbf{x}^\top \mathbf{u})^2 - (\mathbf{x}^\top \mathbf{u})^2 \mathbf{u}^\top \mathbf{u} & \mathbf{x} \in \mathcal{S}_2 \\
    0 & \mathbf{x} \in \mathcal{S}_3 \\
    (\mathbf{x}^\top \mathbf{l})^2 - (\mathbf{x}^\top \mathbf{l})^2 \mathbf{l}^\top \mathbf{l} & \mathbf{x} \in \mathcal{S}_4
    \end{cases} \\
    &= 0
\end{aligned} 
$$

The likelihood ratio test statistic is then:

$$
\begin{aligned}
\text{LRT} 
    &= -2 (\ell(\hat{\boldsymbol{\theta}}_0; \mathbf{x}) - \ell(\hat{\boldsymbol{\theta}}_1; \mathbf{x})) \\
    &= -2 \left[ \left(-\frac{2}{2} \log(2 \pi) - \frac{1}{2}(\mathbf{x} - \hat{\boldsymbol{\theta}}_0)^\top(\mathbf{x} - \hat{\boldsymbol{\theta}}_0)\right) - \left(-\frac{2}{2} \log(2 \pi) - \frac{1}{2}(\mathbf{x} - \hat{\boldsymbol{\theta}}_1)^\top(\mathbf{x} - \hat{\boldsymbol{\theta}}_1)\right) \right] \\
    &= \rvert \rvert \mathbf{x}\rvert \rvert_2^2 - \rvert \rvert \mathbf{x} - \hat{\boldsymbol{\theta}}_1 \rvert \rvert_2^2 \\
    &= \rvert \rvert \mathbf{x} \rvert \rvert_2^2 - \left( \rvert \rvert \hat{\boldsymbol{\theta}}_1 \rvert \rvert_2^2 - 2 \langle \mathbf{x}, \hat{\boldsymbol{\theta}}_1 \rangle + \rvert \rvert \mathbf{x} \rvert \rvert_2^2 \right) \\
    &= -\rvert \rvert \hat{\boldsymbol{\theta}}_1 \rvert \rvert_2^2 + 2 \langle \mathbf{x}, \hat{\boldsymbol{\theta}}_1 \rangle - 2 \underbrace{\langle \mathbf{x} - \hat{\boldsymbol{\theta}}_1, \hat{\boldsymbol{\theta}}_1 \rangle}_{=0} \\
    &= -\rvert \rvert \hat{\boldsymbol{\theta}}_1 \rvert \rvert_2^2 + 2 \langle \mathbf{x}, \hat{\boldsymbol{\theta}}_1 \rangle - 2\langle \mathbf{x}, \hat{\boldsymbol{\theta}}_1\rangle + 2 \langle \hat{\boldsymbol{\theta}}_1, \hat{\boldsymbol{\theta}}_1 \rangle\\
    &= \rvert \rvert \hat{\boldsymbol{\theta}}_1 \rvert \rvert_2^2
\end{aligned}
$$

### Null Distribution
Since the $\mathcal{S}_1$, $\mathcal{S}_2$, $\mathcal{S}_3$, and $\mathcal{S}_4$ partition the plane, we can use the law of total probability to derive a tail bound. For $c > 0$, we see that:

$$
\begin{aligned}
\mathbb{P}(\text{LRT} \leq c ) &= \sum_{i = 1}^4 \mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_i) \mathbb{P}(\mathbf{x} \in \mathcal{S}_i)
\end{aligned}
$$

Note that, under $H_0$, the probability of $\mathbf{x}$ falling in each quadrant is $\frac{1}{4}$ because of its symmetry about the origin and the fact that the quadrants partition the plane into four equal sized pieces. Similarly, the probability of $\mathbf{x}$ falling in any $\mathcal{S}_i$ is equal to the proportion of the circumference of the unit circle that is covered by the arc contained in $\mathcal{S}_i$. We now look at each $\mathcal{S}_i$ separately.

Since the upper and lower boundaries of $\mathcal{S}_1$ and $\mathcal{S}_3$ are orthogonal, the angles at the vertices of $\mathcal{S}_2$ and $\mathcal{S}_4$ are both $\frac{\pi}{2}$ radians. We will let $\gamma$ denote the angle (in radians) at the vertex of $\mathcal{S}_3 = \mathcal{C}^0$. Then the angle at the vertex of $\mathcal{S}_1$ is $\pi - \gamma$ radians.

Thus:

$$
\mathbb{P}\left(\mathbf{x} \in \mathcal{S}_i \right)
= \begin{cases}
\frac{\pi - \gamma}{2 \pi} & i = 1\\
\frac{1}{4} & i = 2 \\
\frac{\gamma}{2 \pi} & i = 3\\
\frac{1}{4} & i = 4
\end{cases}
$$

We can then proceed in cases to find $\mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_i)$ for each $i$.

#### Case 1: $$\mathbf{x} \in \mathcal{S}_1$$
Under $H_0$, $$\mathbf{x} \sim \mathcal{N}(\mathbf{0}_2, \mathbb{I}_{2 \times 2})$$. In polar coordinates, we have $r^2 = x_1^2 + x_2^2$ and $\alpha = \arctan\left(frac{x_2}{x_1}\right)$, which can be shown to be independent. This implies that the event that $x_1^2 + x_2^2 = \rvert \rvert \mathbf{x} \rvert \rvert_2^2 \leq c$ is independent of the event that $\mathbf{x} \in \mathcal{S}_1$. Thus:

<aside><p>Include proof of independence?</p></aside>

$$
\begin{aligned}
\mathbb{P}\left(\text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_1 \right) 
    &= \mathbb{P}\left(\rvert \rvert \hat{\boldsymbol{\theta}}_1 \rvert \rvert_2^2 \leq c \rvert \mathbf{x} \in \mathcal{S}_1 \right) \\
    &= \mathbb{P}\left(\rvert \rvert \mathbf{x} \rvert \rvert_2^2 \leq c \rvert \mathbf{x} \in \mathcal{S}_1 \right) \\
    &= \mathbb{P}\left( x_1^2 + x_2^2 \leq c \right) & \left(r \perp \alpha \right)\\
    &= \mathbb{P}(\chi_2^2 \leq c ) & \left( x_1, x_2 \overset{iid}{\sim} \mathcal{N}(0, 1) \right)
\end{aligned}
$$

#### Case 2: $$\mathbf{x} \in \mathcal{S}_2$$
Consider $\mathbf{u}$ and $\mathbf{u}^0$, which we know satisfy $\mathbf{u} \perp \mathbf{u}^0$. We can consider a new coordinate system with axes given by the directions of $\mathbf{u}$ and $\mathbf{u}^0$. This is basically just rotating the standard axes to align with $\mathbf{u}$ and $\mathbf{u}^0$, which essentially aligns $\mathcal{S}_2$ with the first quadrant of the new coordinate plane, $\mathcal{Q}_1'$. 

Since $\hat{\boldsymbol{\theta}}_1$ is the projection of $\mathbf{x}$ onto $\mathbf{u}$, which is our new $x$-axis, the LRT statistic in the new coordinate system only concerns the first coordinate of $\mathbf{x}$, while the other has been projected to zero. Let $\text{LRT}'$, $\mathbf{x}'$, $x_1'$, and $x_2'$ denote the test statistic and the data in the new coordinate system. We see that:

$$
\begin{aligned}
\mathbb{P}\left( \text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_2 \right)
    &= \mathbb{P}\left(\text{LRT}' \leq c \rvert \mathbf{x}' \in \mathcal{Q}'_1 \right) \\
    &= \mathbb{P}\left((x'_1)^2 \leq c \rvert x'_1 \geq 0, x'_2 \geq 0 \right) \\
    &= \mathbb{P}\left((x'_1)^2 \leq c \rvert x'_1 \geq 0\right) & \left((x'_1)^2 \perp x'_2 \right) \\
    &= \mathbb{P}\left( (x'_1)^2 \leq c \right) & \left(x'_1 \text{ sym. about } 0 \right) \\
    &= \mathbb{P}\left( \chi_1^2 \leq c \right) & \left(x'_1 \sim \mathcal{N}(0, 1)\right)
\end{aligned}
$$

#### Case 3: $$\mathbf{x} \in \mathcal{S}_3$$
In this case, $\hat{\boldsymbol{\theta}}_1 = \mathbf{0}_2$ no matter where $\mathbf{x}$ falls in $\mathcal{S}_3$. Thus:

$$
\begin{aligned}
\mathbb{P}\left(\text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_3\right)
    &= \mathbb{P}\left( \rvert \rvert \mathbf{0}_2 \rvert \rvert_2^2 \leq c \rvert \mathbf{x} \in \mathcal{S}_3 \right) \\
    &= \mathbb{P}\left(0 \leq c \right) \\
    &= \mathbb{P}\left(\chi^2_0 \leq c \right)
\end{aligned}
$$

Note that the last line equals one since $c > 0$.

#### Case 4: $$\mathbf{x} \in \mathcal{S}_4$$
This case is the same as Case 2 but with the lower boundaries of $\mathcal{C}$ and $\mathcal{C}^0$. We get:

$$
\mathbb{P}\left( \text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_4 \right) = \mathbb{P}(\chi_1^2 \leq c)
$$

We then combine the above results to get:

$$
\begin{aligned}
\mathbb{P}(\text{LRT} \leq c ) 
    &= \sum_{i = 1}^4 \mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x} \in \mathcal{S}_i) \mathbb{P}(\mathbf{x} \in \mathcal{S}_i) \\
    &= \frac{\pi - \gamma}{2 \pi} \mathbb{P}\left(\chi_2^2 \leq c \right) 
    + \frac{1}{4} \mathbb{P} \left(\chi_1^2 \leq c \right)
    + \frac{\gamma}{2 \pi} \mathbb{P} \left(\chi_0^2 \leq c \right)
    + \frac{1}{4} \mathbb{P} \left(\chi_1^2 \leq c \right)  \\
    &=  \frac{\gamma}{2 \pi} \mathbb{P}\left(\chi_0^2 \leq c \right) + \pi \mathbb{P} \left(\chi_1^2 \leq c \right)  + \frac{\pi - \gamma}{2 \pi} \mathbb{P}\left(\chi_2^2 \leq c \right) 
\end{aligned}
$$

### Correlated Setting
Now, suppose $\mathbf{x} \sim \mathcal{N}\left(\boldsymbol{\theta}, \boldsymbol{\Sigma}\right)$ for some positive definite covariance matrix $\boldsymbol{\Sigma} \neq \mathbb{I}_{2 \times}$. Consider testing:

$$
\begin{equation}
\label{eq:lrt-ex-2}
H_0: \boldsymbol{\theta} = \mathbf{0}_2 
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\theta} \in \mathcal{C}
\end{equation}
$$

Let $\mathbf{A}$ be the matrix such that $\mathbf{A}^\top \mathbf{A} = \boldsymbol{\Sigma}^{-1}$. As discussed <a href="#a-bit-of-geometry">earlier</a>, the transformed data vector $\mathbf{y} = \mathbf{A} \mathbf{x}$ satisfies $\mathbf{y} \sim \mathcal{N}\left( \boldsymbol{\mu}, \mathbf{V}\right)$ with $\boldsymbol{\mu} = \mathbf{A} \boldsymbol{\theta}$ and $\mathbf{V} = \mathbb{I}_{2 \times 2}$. If we also define:

$$
\mathcal{P} = \left\{ \mathbf{c} : \mathbf{R} \mathbf{A}^{-1} \mathbf{c} \geq \mathbf{0}_2 \right\} 
$$

then testing problem described in Eq. \eqref{eq:lrt-ex-2} is equivalent to using $\mathbf{y}$ to test:

$$
H_0: \boldsymbol{\mu} = \mathbf{0} 
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\mu} \in \mathcal{P}
$$

which is the same as the problem in Eq. \eqref{eq:lrt-ex-1}.

---

## Chi-Bar-Square Distribution
The mixture of $\chi^2$ distributions derived in the previous section is an example of a $\bar{\chi}^2$ (chi-bar-square) distribution. We formally define this concept now.

<div id="chi-bar-square"></div>
<div class="definition">
  <strong>Definition (Chi-Bar-Square Distribution).</strong>
  <br>
  Let $\mathcal{C} \subset \mathbb{R}^p$ be a closed, convex cone, and let $\mathbf{z} \sim \mathcal{N}(\mathbf{0}_p, \mathbf{V})$ for positive definite $\mathbf{V}$. 
  The <i>chi-bar-square distribution</i> is the distribution of the quantity:
  $$
  \mathbf{z}^\top \mathbf{V}^{-1} \mathbf{z} - \underset{\mathbf{c} \in \mathcal{C}}{\min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\}
  $$
  We let $\bar{\chi}^2(\mathbf{V}, \mathcal{C})$ denote such a random variable.
</div>

### Geometric Interpretation
Let $$\Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C})$$ denote the $\mathbf{V}$-projection of $\mathbf{z}$ onto $\mathcal{C}$, and let $$\Pi_{\mathbf{V}}(\mathbf{z} \rvert \mathcal{C}^0)$$ denote its $\mathbf{V}$-projection onto $\mathcal{C}^0$, the polar cone of $\mathcal{C}$. By definition, these projections satisfy:

$$
\begin{aligned}
\Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) &= \underset{\mathbf{c} \in \mathcal{C}}{\arg \min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\} \\
\Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) &= \underset{\mathbf{c} \in \mathcal{C}^0}{\arg \min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\}
\end{aligned}
$$

It can be shown formally that $$\Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) \perp \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0)$$. Intuitively, these vectors lie on the boundaries of $\mathcal{C}$ and $\mathcal{C}^0$, respectively. Since every vector in $\mathcal{C}^0$ must form a non-acute angle with every vector in $\mathcal{C}$, the result must be true. It can also be shown that $$\mathbf{z} = \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) + \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0)$$.

<aside><p>See Proposition 3.12.4 for a proof.</p></aside>

These vectors describe a rectangle with one vertex at the point where $\mathcal{C}$ and $\mathcal{C}^0$ meet. This rectangle has a diagonal equal to $\mathbf{z}$ and sides (denoted by $s_i$ for $i = 1, 2, 3, 4$) given by:

<aside><p>Fig. 3.9 on pg. 76 is helpful here.</p></aside>

$$
\begin{aligned}
s_1 &= \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) \\
s_2 &= \mathbf{z} - \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) = \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0)\\
s_3 &= \mathbf{z} - \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) = \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C})\\
s_4 &= \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0)
\end{aligned}
$$

Applying th Pythagorean theorem to $\mathbf{z}$, $s_1$ and $s_2$, and also to $\mathbf{z}$, $s_3$, and $s_4$, we see that:

$$
\begin{aligned}
&\rvert \rvert \mathbf{z} \rvert \rvert_{\mathbf{V}}^2 
    = \rvert \rvert s_1 \rvert \rvert_{\mathbf{V}}^2 + \rvert \rvert s_2 \rvert \rvert^2_\mathbf{V}  \\
\implies & \mathbf{z} \mathbf{V}^{-1} \mathbf{z} = \Pi^\top_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) \mathbf{V}^{-1} \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) + \rvert \rvert \mathbf{z} - \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) \rvert \rvert_\mathbf{V}^2 \\
\implies & \mathbf{z} \mathbf{V}^{-1} \mathbf{z} = \underset{\mathbf{c} \in \mathcal{C}}{\min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\} + \rvert \rvert \mathbf{z} - \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) \rvert \rvert_\mathbf{V}^2 \\
\implies & \bar{\chi}^2(\mathbf{V}, \mathcal{C}) = \rvert \rvert \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) \rvert \rvert_\mathbf{V}^2 \\
&\rvert \rvert \mathbf{z} \rvert \rvert_{\mathbf{V}}^2 
    = \rvert \rvert s_3 \rvert \rvert_{\mathbf{V}}^2 + \rvert \rvert s_2 \rvert \rvert^4_\mathbf{V}  \\
\implies & \mathbf{z} \mathbf{V}^{-1} \mathbf{z} = \Pi^\top_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) \mathbf{V}^{-1} \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) + \rvert \rvert \mathbf{z} - \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) \rvert \rvert_\mathbf{V}^2 \\
\implies & \mathbf{z} \mathbf{V}^{-1} \mathbf{z} = \underset{\mathbf{c} \in \mathcal{C}^0}{\min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\} + \rvert \rvert \mathbf{z} - \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}^0) \rvert \rvert_\mathbf{V}^2 \\
\implies & \bar{\chi}^2(\mathbf{V}, \mathcal{C}^0) = \rvert \rvert \Pi_\mathbf{V}(\mathbf{z} \rvert \mathcal{C}) \rvert \rvert_\mathbf{V}^2 
\end{aligned}
$$

Thus, the length of the orthogonal projection of $\mathbf{z}$ onto a cone $\mathbf{C}$ with respect to a positive definite $\mathbf{V}$ is a chi-bar-square random variable. 

### Results
We finally come to the main theorem of the section.

<div id="theorem-3-4-2"></div>
<div class="theorem">
  <strong>Theorem 3.4.2.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  {% tabs th-3-4-2 %}
  {% tab th-3-4-2 statement %}
  Let $\mathcal{C} \subset \mathbb{R}^p$ be a closed, convex cone, and let $\mathbf{V}$ be a $p \times p$ positive definite matrix. The distribution of the quantity $\bar{\chi}^2(\mathbf{V}, \mathcal{C})$ is:

  $$
  \mathbb{P}\left(\bar{\chi}^2(\mathbf{V}, \mathcal{C}) \leq c \right)
    = \sum_{i = 0}^p \omega_i(p, \mathbf{V}, \mathcal{C}) \mathbb{P}\left(\chi_i^2 \leq c \right)
  $$
  {% endtab %}
  {% tab th-3-4-2 proof %}
  Proof to be included.
  {% endtab %}
  {% endtabs %}
</div>

Recall that in the <a href="#likelihood-ratio-test-in-two-dimensions">previous section</a> we saw that the likelihood ratio test statistic for testing hypotheses of the form given in Eq. \eqref{eq:lrt-ex-1} or Eq. \eqref{eq:lrt-ex-2} has a null distribution given by a mixture of $\chi^2$ distributions. We also saw that the LRT statistic is equal to the squared length of the maximum likelihood estimate under $H_1$ and that this MLE will be the projection of the data point onto the alternative parameter space. Since the alternative parameter space is a closed, convex cone, the LRT statistic is a chi-bar-square random variable (as we saw <a href="#geometric-interpretation">above</a>). 

An alternative way to show this is to look at the LRT statistic more directly. Let $\mathbf{x} \sim \mathcal{N}\left(\boldsymbol{\theta}, \boldsymbol{\Sigma}\right)$ and consider testing $\boldsymbol{\theta} = \mathbf{0}_k$ against $\boldsymbol{\theta} \in \mathcal{C} \setminus \mathbf{0}_k$. We have:

$$
\begin{aligned}
\text{LRT}
    &= -2 \left[ - \frac{1}{2} \left(\mathbf{x} - \hat{\boldsymbol{\theta}}_0 \right)^\top \boldsymbol{\Sigma}^{-1} \left(\mathbf{x} - \hat{\boldsymbol{\theta}}_0 \right) + \frac{1}{2} \left(\mathbf{x} - \hat{\boldsymbol{\theta}}_1 \right)^\top \boldsymbol{\Sigma}^{-1} \left(\mathbf{x} - \hat{\boldsymbol{\theta}}_1 \right)\right] \
    &=  \left(\mathbf{x} -  \mathbf{0}_k \right)^\top \boldsymbol{\Sigma}^{-1} \left(\mathbf{x} - \mathbf{0}_k \right) + \underset{\boldsymbol{\theta} \in \mathcal{C}}{\max} \left\{ - \left(\mathbf{x} - \boldsymbol{\theta} \right)^\top \boldsymbol{\Sigma}^{-1} \left(\mathbf{x} - \boldsymbol{\theta}  \right) \right\} \\
    &= \mathbf{x}^\top \boldsymbol{\Sigma}^{-1} \mathbf{x}  + \underset{\boldsymbol{\theta} \in \mathcal{C}}{\min} \left\{ \left(\mathbf{x} - \boldsymbol{\theta} \right)^\top \boldsymbol{\Sigma}^{-1} \left(\mathbf{x} - \boldsymbol{\theta}  \right) \right\} 
\end{aligned}
$$

which is exactly the form of $\bar{\chi}^2$.

The weights $\omega_i(p, \mathbf{V}, \mathcal{C})$ can be difficult to compute exactly and sometimes even approximately. However, we do have some nice properties in certain cases. Let $\mathbf{V}$ be a non-singular convariance matrix, and let $\mathcal{C}$ be a closed, convex cone. 

The weights will satisfy $$\sum_{i = 0}^p (-1)^i \omega_i(p, \mathbf{V}, \mathcal{C}) = 0$$ and $$0 \leq \omega_i(p, \mathbf{V}, \mathcal{C}) \leq \frac{1}{2}$$. In addition, $$\omega_i(p, \mathbf{V}, \mathcal{C}^0) = \omega_{p - i}(p, \mathbf{V}, \mathcal{C})$$. This implies that if $\mathcal{C}$ is the non-negative orthant, then $$\omega_i(p, \mathbf{V}, \mathcal{C}) = \omega_{p-i}(p, \mathbf{V}^{-1}, \mathcal{C})$$. If, in addition, $\mathbf{z} \sim \mathcal{N}(\mathbf{0}_k, \mathbf{V})$, then $\omega_i(p, \mathbf{V}, \mathcal{C})$ is the probability that $\Pi_{\mathbf{V}}(\mathbf{z} \rvert \mathcal{C})$ has exactly $i$ positive coordinates. 

<aside><p>See Proposition 3.6.1.</p></aside>

If the constraints that define the closed, convex cone $\mathcal{C}$ are linear and independent (e.g. $\boldsymbol{\theta}_j \geq 0$ for all $j$), then we can simplify our setting slightly. We can transform our setting so that we instead consider $$\omega_i(k', \mathbf{W}, \mathcal{O}^{k'})$$ for some $k'$, $\mathbf{W}$, and non-negative orthant in $k'$ dimensions $\mathcal{O}^{k'}$. 

### Extensions
The previous sub-section covered the case of a simple null hypothesis. However, the results still hold when $$H_0: \boldsymbol{\theta} \in \mathcal{M}$$ for some vector space $\mathcal{M} \subset \mathcal{C}$. 

<div id="theorem-3-7-1"></div>
<div class="theorem">
  <strong>Theorem 3.7.1.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  {% tabs th-3-7-1 %}
  {% tab th-3-7-1 statement %}
  Let $\mathcal{C} \subset \mathbb{R}^p$ be a closed, convex cone, and let $\mathbf{V}$ be a $p \times p$ positive definite matrix. Let $\mathbf{x} \in \mathcal{N}(\boldsymbol{\theta}, \mathbf{V})$ and consider testing:
  
  $$
    H_0: \boldsymbol{\theta} \in \mathcal{M} 
    \hspace{5mm} \text{vs.} \hspace{5mm}
    H_1: \boldsymbol{\theta} \in \mathcal{C} \setminus \mathcal{M}
  $$

  for some vector space $\mathcal{M} \subset \mathcal{C}$. The likelihood ratio test statistic has a null distribution:

  $$
  \mathbb{P}\left(\text{LRT} \leq c \rvert H_0 \right) = \sum_{i = 0}^p \omega_i\left(p, \mathbf{V}, \mathcal{C} \cap \mathcal{M}^\perp \right) \mathbb{P}\left(\chi_i^2 \leq c\right)
  $$

  where $\mathcal{M}^\perp$ is the orthogonal complement of $\mathcal{M}$ with respect to $\langle \cdot, \cdot, \rangle_\mathbf{V}$. Also, this null distribution is the same for every value $\boldsymbol{\theta}$ can take on under $H_0$ (i.e. it is <strong>similar</strong>).
  {% endtab %}
  {% tab th-3-7-1 proof %}
  Proof to be included.
  {% endtab %}
  {% endtabs %}
</div>

---

## Unknown Covariance
The previous results required the variance-covariance matrix of the data to be known. This allows us to focus our attention on <i>only</i> the quadratic form in the log-likelihood function (we ignored the determinant term). There are several ways to go about accounting for an unknown covariance matrix, and readers should look at sub-section 3.10 for details.




