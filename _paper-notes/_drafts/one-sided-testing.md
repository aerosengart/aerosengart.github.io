---
layout: distill
title: One-Sided Testing Problems in Multivariate Analysis
description: 
date: 2026-02-18
tabs: true
tags: likelihood theory testing paper-review
toc:
  - name: Notation
  - name: Definitions
bibliography: 2026-02-18-one-sided-testing.bib
---

Perlman's 1969<d-cite key=perlman1969></d-cite> paper on one-sided testing problems is cited in Silvapulle and Silvapulle's work<d-cite key=silvapulle1997></d-cite> on score-type tests for one-sided alternatives that I worked through (exactly!) one year ago (see <a href="/paper-notes/score-test-one-sided/">here</a>). I've decided to go through some of this earlier work in order to understand where the $p$-value bound comes from for Silvapulle and Silvapulle's test statistics.

<aside><p>The key contribution of Perlman is the study of likelihood ratio tests against one-sided alternatives with $\Sigma$ <i>unknown</i>.</p></aside>

---

## Notation
We will assume to have $N$ independent $p$-dimensional multivariate Gaussian observations with mean $\mu$ and covariance matrix $\Sigma$. First, I'll reiterate Perlman's notation conventions and set-up for easy reference.


The variables $c$ and $d$ will be used to denote arbitrary positive constants. $\mathcal{E}_p$ will be used to denote $p$-dimensional Euclidean space, and inequalities for vectors will be applied element-wise (e.g. $\mathbf{x} \geq 0$ denotes a vector with non-negative coordinates). Arbitrary $(p-1)$-dimensional subspaces in $\mathcal{E}_p$ are denoted by $\mathcal{K}$, and $\mathcal{O}$ is the closed positive orthant. 

It is a fact that there exist two closed halfspaces such that their intersection equals $\mathcal{K}$, and we denote either one of these with $\mathcal{K}^+$. $\mathcal{L}$ will denote a ray from the origin, and $A(\mathcal{A})$ will denote the image of any set $\mathcal{A}$ under the linear transformation $A: \mathcal{E}_p \rightarrow \mathcal{E}_p$. 

Similar to the case with vectors, if matrix $\mathbf{M}$ is square, then $\mathbf{M} > 0$ will specify that the matrix is symmetric positive definite. Its unique square root will be denoted by $$\mathbf{M}^{\frac{1}{2}}$$. Such matrices can be partitioned into $k^2$ blocks, which we denote by $$\mathbf{M}_{i,j}$$ for $$i, j = 1, \dots, k$$. We let $$\mathbf{M}_{i,i \cdot j} = \mathbf{M}_{i,i} - \mathbf{M}_{i,j} \mathbf{M}_{j,j}^{-1} \mathbf{M}_{j,i}$$, the <a href="https://en.wikipedia.org/wiki/Schur_complement">Schur complement</a> of the block $\mathbf{M}_{j,j}$. 

We denote the mutivariate distribution with mean $\mu$ and covariance matrix $\Sigma$ with $\mathcal{N}(\mu, \Sigma)$, and we denote the <a href="https://en.wikipedia.org/wiki/Wishart_distribution">Wishart distribution</a> with $n$ degrees of freedom and expectation $n \Sigma$ with $\mathcal{W}(n, \Sigma)$. The chi-squared distribution with $n$ degrees of freedom will be denoted by $\chi_n^2$.

<aside><p>The Wishart distribution is a multivariate generalization of the Gamma distribution.</p></aside>

---

## Definitions
Next, we'll cover some of the important definitions for sets needed for later results.

<div id="positively-homogeneous"></div>
<div class="definition">
<strong>Definition (Positively Homogeneous).<d-cite key=perlman1969></d-cite></strong>
<br>
A set $\mathcal{P} \in \mathcal{E}_p$ is called <i>positively homogeneous</i> if:

$$
\mathbf{x} \in \mathcal{P} \implies c \mathbf{x} \in \mathcal{P} \hspace{5mm} \forall c \in \mathbb{R}^+
$$
</div>

In what follows, we assume that positively homogeneous sets are closed.

<div id="one-sided"></div>
<div class="definition">
<strong>Definition (One-sided).<d-cite key=perlman1969></d-cite></strong>
<br>
A set $\mathcal{A}$ with at least one non-zero point is called <i>one-sided</i> (with respect to the origin) if there exists at least one non-zero point $\mathbf{z}$ such that $\mathbf{a}^\top \mathbf{z} > 0$ for all non-zero $\mathbf{a} \in \mathcal{A}$.
</div>

We call a closed, positively homogeneous, one-sided set a <strong>cone</strong> and denote it by $\mathcal{C}$. 

---

## Preliminary Results
Below we list a few properties of cones that will be helpful later on.

<div id="cone-properties"></div>
<div class="theorem">
<strong>Lemmata (Cone Properties).<d-cite key=perlman1969></d-cite></strong>
<br>
Let $\mathcal{C}$ be any arbitrary cone.
<ol>
<li>For $\mathcal{C}$, there exists a ray, $\mathcal{L}$, and a half-space, $\mathcal{K}^+$ such that $\mathcal{L} \subset \mathcal{C} \subset \mathcal{K}^+$.</li>
<li>For $\mathcal{C}$, there exists a right circular cone $\mathcal{C}_\lambda$ with $0 < \lambda < 1$ such that $\mathcal{C} \subset \mathcal{C}_\lambda$.</li>
<li>If $\mathcal{C}$ contains an open set, then there exists a right circular cone $\mathcal{C}_\nu$ with $0 < \nu < 1$ such that $\mathcal{C}_\nu \subset \mathcal{C}$.</li>
<li>For any right circular cone $\mathcal{C}_\lambda$ with $0 < \lambda < 1$, there exists a ray, $\mathcal{L}$, a half-space, $\mathcal{K}^+$, and sequences $\{ A_n \}$ and $\{ B_n \}$ of invertible linear transformations that satisfy:
$$
\begin{aligned}
&A_n(\mathcal{C}_\lambda) \supset A_{n + 1}(\mathcal{C}_\lambda) 
&\bigcap_{n = 1}^\infty A_n(\mathcal{C}_n) = \mathcal{L} \\
&B_n(\mathcal{C}_\lambda) \subset B_{n+1}(\mathcal{C}_\lambda) 
&\bigcup_{n = 1}^\infty B_n(\mathcal{C}_\lambda) = \text{int}(\mathcal{K}^+)
\end{aligned}
$$
</li>
</ol>
</div>

<aside><p>A <strong>right circular cone</strong> is a cone with a circular base that is perpendicular to its axis (the line passing from its base to its apex).
<br><br>
We use $\text{int}(\mathcal{A})$ to denote the <strong>interior</strong> of $\mathcal{A}$, which is the set of points in $\mathcal{A}$ such that the open balls centered at those points are completely contained in $\mathcal{A}$.</p></aside>

We also note the following decomposition of a quadratic form constructed from a partitioned vector and a partitioned matrix.

<div id="partitioned-quad-form"></div>
<div class="theorem">
<strong>Lemma (Partitioned Vector and Matrix).<d-cite key=perlman1969></d-cite></strong>
<br>
Let $\mathbf{z}$ be a $p$-dimensional vector, and let $\mathbf{V}$ be a $p \times p$ matrix satisfying $\mathbf{V} > 0$. Suppose $\mathbf{z}^\top = (\mathbf{z}_1^\top, \mathbf{z}_2^\top)$ with $q$-dimensional $\mathbf{z}_1$, and suppose $\mathbf{V}$ is partitioned into four blocks with $\mathbf{V}_{1,1}$ being $q \times q$. We have:

$$
\begin{aligned}
\mathbf{z}^\top \mathbf{V}^{-1} \mathbf{z} 
&= (\mathbf{z}_1 - \mathbf{V}_{1,2} \mathbf{V}_{2,2}^{-1} \mathbf{z}_2)^\top \mathbf{V}_{1,1 \cdot 2}^{-1}(\mathbf{z}_1 - \mathbf{V}_{1,2} \mathbf{V}_{2,2}^{-1} \mathbf{z}_2) + \mathbf{z}_2^\top \mathbf{V}_{2,2}^{-1} \mathbf{z}_2 \\
&= (\mathbf{z}_2 - \mathbf{V}_{2,1}\mathbf{V}_{1,1}^{-1}\mathbf{z}_1)^\top\mathbf{V}_{2,2 \cdot 1}^{-1} (\mathbf{z}_2 - \mathbf{V}_{2,1} \mathbf{V}_{1,1}^{-1} \mathbf{z}_1) + \mathbf{z}_1^\top \mathbf{V}_{1,1}^{-1} \mathbf{z}_1
\end{aligned}
$$
</div>

We also have the follows results, which we will not prove. Perlman states that they can be found in a mimeograph from Stein (1966), but I could not find them online.


<div id="multivariate-distribution"></div>
<div class="theorem">
<strong>Lemmata (Multivariate Distribution Theory).<d-cite key=perlman1969></d-cite></strong>
<br>
Let $\mathbf{x} \sim \mathcal{N}(\mu, \Sigma)$ be a $p$-dimensional multivariate Gaussian random variable, and let $\mathbf{S} \sim \mathcal{W}(n, \Sigma)$ be an $(p \times p)$-dimensional random variable independent of $\mathbf{x}$. We construct the partition $\mathbf{x}^\top = (\mathbf{x}_1^\top, \mathbf{x}_2^\top)$ with $k$-dimensional $\mathbf{x}_1$ where $1 \leq k < p$. We partition $\mathbf{S}$, $\mu$, and $\Sigma$ similarly. Define:

$$
\begin{aligned}
\gamma &= \frac{\mathbf{x}_1 - \mathbf{S}_{1,2} \mathbf{S}_{2,2}^{-1}\mathbf{x}_2}{\left(1 + \mathbf{x}_2^\top \mathbf{S}_{2,2}^{-1}\mathbf{x}_2\right)^{\frac{1}{2}}} \\
\delta &= \frac{\gamma^\top \mathbf{S}^{-1}_{1,1} \gamma}{\gamma^\top \Sigma^{-1}_{1,1 \cdot 2}\gamma}
\end{aligned}
$$

The following results hold:

<ol>
<li>$\mathbf{x}_1 \mid \mathbf{x}_2 \sim \mathcal{N}(\mu_1 + \Sigma_{1,2}\Sigma_{2,2}^{-1}(\mathbf{x}_2 - \mu_2), \Sigma_{1,1 \cdot 2})$</li>
<li>$\mathbf{S}_{1,1 \cdot 2} \sim \mathcal{W}(n - p + k, \Sigma_{1,1 \cdot 2})$ and is independent of $\mathbf{S}_{1,2}$ and $\mathbf{S}_{2,2}$</li>
<li>$\mathbf{S}_{1,2}\mathbf{S}_{2,2}^{-\frac{1}{2}} \mid \mathbf{S}_{2,2} \sim \mathcal{N}(\Sigma_{1,2}\Sigma_{2,2}^{-1}\mathbf{S}_{2,2}^{\frac{1}{2}}, \Sigma_{1,1 \cdot 2} \otimes \mathbb{I}_{k \times k})$</li>
<li>$\frac{\mathbf{x}^\top \mathbf{S}^{-1}\mathbf{x}}{\mathbf{x}^\top \Sigma^{-1}\mathbf{x}} \sim \frac{1}{\chi^2_{n - p + 1}}$ and is independent of $\mathbf{x}$</li>
<li>$\mathbf{x}^\top \mathbf{S}^{-1}\mathbf{x} \sim \frac{\chi_p^2(\mu^\top \Sigma^{-1}\mu)}{\chi^2_{n - p + 1}}$ which is the ratio of two independent chi-squared variables with the numerator being non-central</li>
<li>$\delta \rvert \mathbf{S}_{1,2}, \mathbf{S}_{2,2}, \mathbf{x}_1, \mathbf{x}_2 \sim \frac{1}{\chi^2_{n - p + 1}}$, implying that $\delta$ is independent of $\gamma$ as well as $\mathbf{S}_{1,2}$, $\mathbf{S}_{2,2}$, $\mathbf{x}_1$ and $\mathbf{x}_2$</li>
<li>$\mathbf{x}_1 \rvert \mathbf{S}_{2,2}, \mathbf{x}_2 \sim \mathcal{N}(\mu_1 + \Sigma_{1,2} \Sigma_{2,2}^{-1}(\mathbf{x}_2 - \mu_2), \Sigma_{1,1 \cdot \mu_2})$ and $\mathbf{S}_{1,2}\mathbf{S}_{2,2}^{-1} \mathbf{x}_2 \rvert \mathbf{S}_{2,2}, \mathbf{x}_2 \sim \mathcal{N}(\Sigma_{1,2} \Sigma_{2,2}^{-1} \mathbf{x}_2, (\mathbf{x}_2^\top \mathbf{S}_{2,2}^{-1} \mathbf{x}_2)\Sigma_{1,1 \cdot 2})$ and the two are independent (conditionally)</li>
<li>$\gamma \rvert \mathbf{S}_{2,2} \mathbf{x}_2 \sim \mathcal{N}\left(\frac{\mu_1 - \Sigma_{1,2} \Sigma_{2,2}^{-1}\mu_2}{(1 + \mathbf{x}_2^\top \mathbf{S}_{2,2}^{-1}\mathbf{x}_2)^{\frac{1}{2}}}, \Sigma_{1,1 \cdot 2}\right)$</li>
<li>If $\mu = \mathbf{0}$, then $\gamma \rvert \mathbf{S}_{2,2}, \mathbf{x}_2 \sim \mathcal{N}(\mathbf{0}, \Sigma_{1,1 \cdot 2})$ and is independent of $\mathbf{S}_{2,2}$ and $\mathbf{x}_2$, and $\gamma^\top \Sigma_{1,1 \cdot 2}^{-1} \gamma \rvert \mathbf{S}_{2,2}, \mathbf{x}_2 \sim \chi_k^2$</li>
</ol>
</div>

And finally, we have results concerning <strong>orthogonally invariant</strong> random variables.

<div id="orthogonally-invariant"></div>
<div class="theorem">
<strong>Lemma 3.1.<d-cite key=perlman1969></d-cite></strong>
<br>
Let $\mathbf{z} \in \mathcal{E}_p$ with $p \geq 2$ be a random vector with a density with respect to the Lebesgue measure. Let $\mathcal{S} = \{ \mathbf{x} \mid \mathbf{x}^\top \mathbf{x} = 1 \}$, the unit sphere in $\mathcal{E}_p$, and let $m$ denote the Lebesgue measure. Assume $\mathbf{z}$ is <i>orthogonally invariant</i>; i.e. for any orthogonal matrix $\Gamma$, $\mathbf{z}$ and $\Gamma \mathbf{z}$ are identically distributed. We have:
<ol>
<li>$\mathbf{z}^\top \mathbf{z}$ is independent of $\frac{\mathbf{z}}{\rvert \rvert \mathbf{z} \rvert \rvert_2}$</li>
<li>$\frac{\mathbf{z}}{\rvert \rvert \mathbf{z} \rvert \rvert_2}$ is uniformly distributed over $\mathcal{S}$</li>
<li>For any measurable, positively homogeneous set $\mathcal{P}$, $\mathbb{P}(\mathbf{z} \in \mathcal{P}) = \frac{m(\mathcal{P} \cap \mathcal{S})}{m(\mathcal{S})}$</li>
<li>For any invertible linear transformation $\mathbf{A}: \mathcal{E}_p \rightarrow \mathcal{E}_p$, $\mathbb{P}(\mathbf{A}\mathbf{z} \geq 0) = \frac{m(\mathbf{A}^{-1}(\mathcal{O}) \cap \mathcal{S})}{m(\mathcal{S})}$</li>
</ol>
</div>

---






