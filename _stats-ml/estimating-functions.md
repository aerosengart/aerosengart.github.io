---
layout: distill
title: Inference Functions
description: A Primer
date: 2026-08-14
tabs: true
tags: theory estimation primer
toc:
  - name: Set-Up
  - name: Characteristics
  - name: Properties
bibliography: stats-ml.bib
---

In this post, I'll cover some of the main ideas about inference functions, which are kind of like a parent concept to estimating equations and estimates. This is by no means exhaustive or even complete; see the end of the post for a couple of other sources on the topic.

---

## Set-Up
Let's start with a quick overview of what <i>inference functions</i> are. Chapter 3 of Song (2007)<d-cite key=song2007></d-cite> has a very good one, which we follow here. 

Inference functions provide a method of parameter estimation without the hassle of dealing with the specifics of a probability distribution. We let $\mathcal{X}$ denote our sample space, and we will consider the class, $\mathcal{P}$, of distributions parametrized by $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ with $\boldsymbol{\Theta} \subseteq \mathbb{R}^p$.

An inference function is just a special function of our data and the parameter upon which the sample depends. 

<div class="definition">
<strong>Definition (Inference Function).</strong>
<br>
An <i>inference function</i> (or <i>estimating function</i>) is any function $\boldsymbol{\Psi}: \mathcal{X} \times \boldsymbol{\Theta} \rightarrow \mathbb{R}^p$ such that:

<ul>
<li> $\boldsymbol{\Psi}(\cdot; \boldsymbol{\theta})$ is measurable for any $\boldsymbol{\theta} \in \boldsymbol{\Theta}$</li>
<li>$\boldsymbol{\Psi}(\mathbf{x}; \cdot)$ is continuous in a compact subspace of $\boldsymbol{\Theta}$ containing the true parameter $\boldsymbol{\theta}_0$ for any $\mathbf{x} \in \mathcal{X}$</li>
</ul>
</div>

Clearly, $\boldsymbol{\Psi}$ must have at least $p$ linearly independent components if $\boldsymbol{\theta}$ is $p$-dimensional. However, it has $q > p$ such components, then we say that $\boldsymbol{\theta}$ is <strong>over-identified</strong>. For inference function $\boldsymbol{\Psi}$, the corresponding <storng>estimating equation</strong> is given by:

$$
\boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta}) = \mathbf{0}_p
$$

The solution to the above equation, $\hat{\boldsymbol{\theta}}$, is the <strong>estimate</strong>. Two inference functions are called <strong>equivalent</strong> (denoted as $\boldsymbol{\Psi} \sim \boldsymbol{\Phi}$) if they yield the same estimate for any sample $\mathbf{x} \in \mathcal{X}$. 


---

## Characteristics
We can describe inference functions with several different nice properties.

<div class="definition">
<strong>Definition (Unbiased).</strong>
<br>
An inference function $\boldsymbol{\Psi}$ is called <i>unbiased</i> if it satisfies, for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$:

$$
\mathbb{E}_{\boldsymbol{\theta}}\left[ \boldsymbol{\Psi}(\mathbf{X}; \boldsymbol{\theta}) \right] = \mathbf{0}_p
$$
</div>

An unbiased inference function will have expectation zero for any value of the parameter $\boldsymbol{\theta}$ (where the expectation is with respect to $\mathbf{X} \sim P_{\boldsymbol{\theta}}$). We can also define an <strong>additive inference function</strong>.

<div class="definition">
<strong>Definition (Additive).</strong>
<br>
Let $\mathbf{x} = (x_1, \dots, x_K)^\top$ be a sample where $x_i \overset{iid}{\sim} p_{\boldsymbol{\theta}}$ which are random variables or vectors (in the case of vectors, we permit the coordinates to be correlated). 
<br>
An <i>additive inference function</i> is given by:

$$
\boldsymbol{Psi}_K(\mathbf{x}; \boldsymbol{\theta}) = \sum_{i = 1}^K \boldsymbol{\Psi}(x_i; \boldsymbol{\theta})
$$

$\boldsymbol{\Psi}_K(\cdot; \cdot)$ is the <strong>kernel inference function</strong>.
</div>

An inference function can also be described as being <strong>regular</strong>.

<div class="definition">
<strong>Definition (Regular).</strong>
<br>
Let $\boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta})$ be a $p$-dimensional inference function, and let $\phi_i(\cdot, \cdot)$ denote its $i$-th component. We call $\boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta})$ <i>regular</i> if it satisfies:

<ul>
<li>$\mathbf{E}_{\boldsymbol{\theta}}\left[ \boldsymbol{\Psi}(\mathbf{X}; \boldsymbol{\theta}) \right] = \mathbf{0}_p$, for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$</li>
<li>$\frac{\partial \boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_j}$ exists, for all $\mathbf{x} \in \mathcal{X}$ and $j = 1, \dots, p$</li>
<li>For any $j = 1, \dots, p$ and bounded function, $f(\mathbf{x})$, that is independent of $\boldsymbol{\theta}$, $\frac{\partial}{\partial \boldsymbol{\theta}_j}\left[ \int_{\mathcal{X}} f(\mathbf{x}) \boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta}) p(\mathbf{x}; \boldsymbol{\theta}) d \mathbf{x} \right] = \int_{\mathcal{X}} f(\mathbf{x}) \frac{\partial}{\partial \boldsymbol{\theta}_j} \left[ \boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta}) p(\mathbf{x}; \boldsymbol{\theta}) \right] d\mathbf{x}$</li>
<li>$\mathbf{E}_{\boldsymbol{\theta}}\left[ \phi_j(\mathbf{X}; \boldsymbol{\theta}) \phi_k(\mathbf{X}; \boldsymbol{\theta}) \right]$ exists for all $j, k = 1, \dots, p$, and the $p \times p$ matrix (the <strong>variability matrix</strong>) $\mathbf{V}_{\boldsymbol{\Psi}}(\boldsymbol{\theta}) = \mathbb{E}_{\boldsymbol{\theta}}\left[ \boldsymbol{\Psi}(\mathbf{X}; \boldsymbol{\theta}) \boldsymbol{\Psi}^\top(\mathbf{X}; \boldsymbol{\theta}) \right]$ is positive-definite</li>
<li>The $p \times p$ matrix (the <strong>sensitivity matrix</strong>) $\mathbf{S}_{\boldsymbol{\Psi}}(\boldsymbol{\theta}) = \mathbb{E}_{\boldsymbol{\theta}}\left[ \frac{\partial \boldsymbol{\Psi}(\mathbf{X}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \right]$ is non-singular</li>
</ul>
</div>

Similarly, a statistical model is called regular if its score function is regular. We will use $\mathcal{G}$ to denote the class of $p$-dimensional, regular inference functions.

We will also define the <strong>Crowder class</strong> of regular inference functions, $\mathcal{G}_C$, as the subclass of $\mathcal{G}$ consisting of members satisfying:

$$
\boldsymbol{\Psi}_C(\boldsymbol{\theta}) = \sum_{i = 1}^K \mathbf{C}_i(\boldsymbol{\theta}) \boldsymbol{\Psi}_i(\mathbf{x}_i; \boldsymbol{\theta}); \hspace{5mm} \boldsymbol{\theta} \in \boldsymbol{\Theta} \subseteq \mathbb{R}^p
$$

where $\boldsymbol{\Psi}_i \in \mathcal{G}$ and $\mathbf{C}_i(\boldsymbol{\theta})$ is a non-random matrix function of $\boldsymbol{\theta}$ chosen such that the sequence of roots $$\{ \hat{\boldsymbol{\theta}}_K \}_{K \geq 1}$$ to $\boldsymbol{\Psi}_C(\boldsymbol{\theta}) = \mathbf{0}_p$ is consistent.

---

## Properties
In what follows, it will be helpful to define the quantity:

$$
\lambda(\boldsymbol{\theta}) = \mathbb{E}_{\boldsymbol{\theta}_0}\left[ \boldsymbol{\Psi}(\mathbf{X}; \boldsymbol{\theta})\right] = \int \boldsymbol{\Psi}(x; \boldsymbol{\theta}) p(x; \boldsymbol{\theta}_0) dx
$$

where $\boldsymbol{\theta}_0 \in \boldsymbol{\Theta}$ is the true parameter value. We now state a consistency theorem and an asymptotic normality theorem for additive inference functions.

<div class="theorem">
<strong>Theorem 3.4 (Consistency).</strong>
<br>
Let $\boldsymbol{\Psi}_K$ be an additive inference function that is unbiased at $\boldsymbol{\theta}_0$. If $\lambda(\boldsymbol{\theta})$ has a unique zero at $\boldsymbol{\theta}_0$, then there exists a sequence of roots, $\{ \hat{\boldsymbol{\theta}}_K \}$ to the estimating equation, $\boldsymbol{\Psi}_K(\mathbf{x}; \boldsymbol{\theta}) = \mathbf{0}_p$, such that:

$$
\hat{\boldsymbol{\theta}}_K \overset{p}{\rightarrow} \boldsymbol{\theta}_0; \hspace{5mm} \text{under } P_{\boldsymbol{\theta}_0}
$$
</div>

For asymptotic normality, we restrict our attention to regular inference functions. For a regular inference function, we can define the <strong>Godambe information</strong>:

$$
\mathbf{J}_{\boldsymbol{\Psi}}(\boldsymbol{\theta}) = \mathbf{S}^\top_{\boldsymbol{\Psi}}(\boldsymbol{\theta}) \mathbf{V}^{-1}_{\boldsymbol{\Psi}}(\boldsymbol{\theta}) \mathbf{S}_{\boldsymbol{\Psi}}(\boldsymbol{\theta})
$$

Now we come to the normality result. We omit the proof, but it is similar to the proof of asymptotic normality for maximum likelihood estimators (i.e. relies upon a Taylor expansion and controlling higher order terms with the boundedness condition).

<div class="theorem">
<strong>Theorem 3.11 (Asymptotic Normality).</strong>
<br>
Let $$\{ \hat{\boldsymbol{\theta}}_K\}_{K \geq 1}$$ be a sequence of roots of the estimating equations:

$$
\boldsymbol{\Psi}_K(\boldsymbol{\theta}) = \sum_{i = 1}^K \boldsymbol{\Psi}(\mathbf{x}_i; \boldsymbol{\theta}); \hspace{5mm} K \geq 1
$$

where $\boldsymbol{\Psi} \in \mathcal{G}$. Suppose $\hat{\boldsymbol{\theta}}_K$ is consistent and satisfies:

$$
\left\rvert \left\rvert \frac{\partial^2 \boldsymbol{\Psi}(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right\rvert \right\rvert < M(\mathbf{x}); \hspace{5mm} \boldsymbol{\theta} \in \mathbf{N}(\boldsymbol{\theta}_0)
$$

for a neighborhood, $\mathbf{N}(\boldsymbol{\theta}_0)$, centered at the true parameter value $\boldsymbol{\theta}_0$ and for a $P_{\boldsymbol{\theta}}$-measurable function, $M(\mathbf{x})$, satisfying $\mathbb{E}_{\boldsymbol{\theta}_0}[M(\mathbf{X})] < \infty$. Then:

$$
\sqrt{K}(\hat{\boldsymbol{\theta}}_K - \boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}\left(\mathbf{0}_p, \mathbf{J}^{-1}_{\boldsymbol{\Psi}}(\boldsymbol{\theta}_0)\right); \hspace{5mm} \text{under } P_{\boldsymbol{\theta}_0}
$$
</div>

Let $\mathbf{U}(\mathbf{x}; \boldsymbol{\theta}) = \frac{\partial \log(p(\mathbf{x}; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}}$ denote the score function. If $\mathbf{U} \in \mathcal{G}$, then the Fisher information matrix (for a single observation) is given by:

$$
\mathcal{I}(\boldsymbol{\theta}) = \mathbb{E}_{\boldsymbol{\theta}}\left[ \mathbf{U}(\mathbf{X}; \boldsymbol{\theta})\mathbf{U}^\top(\mathbf{X}; \boldsymbol{\theta}) \right] = - \mathbb{E}_{\boldsymbol{\theta}} \left[ \frac{\partial}{\partial \boldsymbol{\theta}^\top}\left[ \mathbf{U}(\mathbf{X}; \boldsymbol{\theta}) \right] \right]
$$

Regular inference functions also satisfy the <strong>Godambe inequality</strong>.

<div class="theorem">
<strong>Theorem 3.12 (Godambe Inequality).</strong>
<br>
Let $\boldsymbol{\Psi} \in \mathcal{G}$. Then:

$$
\mathcal{I}(\boldsymbol{\theta}) - \mathbf{J}_{\boldsymbol{\Psi}}(\boldsymbol{\theta})
$$

is positive semi-definite $\forall \boldsymbol{\theta} \in \boldsymbol{\Theta}$ and positive definite if, and only if, $\boldsymbol{\Psi} \sim \mathbf{U}$. 
</div>

We also have the <strong>Crowder optimality theorem</strong>, which tells us which inference function will be optimal (in terms of efficiency) over $\mathcal{G}_C$.

<div class="theorem">
<strong>Theorem 3.13 (Crowder Optimality).</strong>
<br>
The optimal inference function in the Crowder class, $\mathcal{G}_C$, is $\boldsymbol{\Psi}_K(\boldsymbol{\theta}) \in \mathcal{G}_C$ with matrix $\mathbf{C}_i(\cdot)$ defined as:

$$
\mathbf{C}_i(\boldsymbol{\theta}) = \mathbb{E}_{\boldsymbol{\theta}}\left[\left(\frac{\partial \boldsymbol{\Psi}_i(\mathbf{X}_i; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \right)^\top \text{Var}_{\boldsymbol{\theta}}^{-1}\left(\boldsymbol{\Psi}_i(\mathbf{X}_i; \boldsymbol{\theta}) \right)\right]
$$
</div>

---

There is a lot more to say about inference functions. See Heyde's <i>Quasi-Likelihood and Its Applications: A General Approach to Optimal Parameter Estimation</i><d-cite key=heyde1997></d-cite> and McLeish's <i>The Theory and Application of Statistical Inference Functions</i><d-cite key=mcleish1988></d-cite>.