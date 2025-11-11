---
layout: post
title:  "One-Sided Score Test"
date: 18 February 2025
categories: posts
tags: ["theory", "likelihood", "score-test"]
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

In many cases, we may want to test the null hypothesis that a parameter is zero against a one-sided alternative (e.g. the parameter is non-negative). In this setting, we are constraining the alternative parameter space, and for some parameters (such as variance components), the value of the parameter under the null may be on the boundary. 

In this post, I'll cover some of the literature on tests against one-sided alternatives. We'll mostly be in the mindset of tests of homogeneity in mixed models (i.e. testing whether the random effects variance is zero), but a lot of these results are generally applicable. 

## Some Intuition
First, let's discuss some of the intuition behind the score test (also called Rao's test (after C. R. Rao) and the Lagrange multiplier test). Suppose we have some model or data-generating process parametrized by $\theta$. Our goal will be to test:

$$
H_0: \theta = \theta_0
$$

Let $\ell(\theta; \mathbf{X})$ denote the log-likelihood function for parameter $\theta$ given data $\mathbf{X}$. We can imagine that if our maximum likelihood estimate, $\hat{\theta}$, is far from $\theta_0$, then our data provide evidence against $H_0$.

The score test uses the slope of the log-likelihood (i.e. the derivative), called the <i>score</i> to determine what it means for $\hat{\theta}$ to be far from $\theta_0$. If the derivative is quite large (in absolute value) at $\theta_0$, then that implies that we have moved quite far away from the root of the log-likelihood, $\hat{\theta}$. 

Under the assumption that the log-likelihood is partially differentiable w.r.t. each component of $\theta$ and the Fisher information exists and is invertible at $\theta_0$, the standard score test statistic for i.i.d. sample $\mathbf{X} = (\mathbf{X}_1, \dots, \mathbf{X}_n)$ is computed as:

$$
t = \frac{1}{n} U^\top(\theta_0) \mathcal{I}^{-1} (\theta_0) U(\theta_0)
$$

where 

$$
U(\theta^*) = \sum_{i = 1}^n \frac{\partial}{\partial \theta} \ell(\theta; \mathbf{X}_i) \rvert_{\theta = \theta^*}; 
\hspace{8mm}
\mathcal{I}(\theta^*) = \mathbb{E}\left[ U(\theta) U^\top(\theta) \right] \bigg\rvert_{\theta = \theta^*}
$$

Assuming $U(\theta)$ has finite variance and mean zero (which it will under certain conditions — see <a href="/posts/2025/02/02/likelihood-theory.html">this post on likelihood theory</a>), it is asymptotically multivariate Gaussian when suitably centered and scaled (by the central limit theorem).[^fn-dasgupta] 

<details>
<summary>Proof.</summary>
Notice that $U(\theta_0)$ is the sum of (functions of) i.i.d. random variables. As we noted above, $U(\theta_0)$ has mean zero under certain regularity conditions. The CLT states that:
$$
\frac{\sqrt{n}\left(\frac{1}{n} \sum_{i = 1}^n X_i - \mu \right) }{\sigma} \rightsquigarrow \mathcal{N}(0, 1)
\hspace{4mm} \iff \hspace{4mm}
(\frac{1}{n} \sum_{i = 1}^n X_i - \mu_X) \rightsquigarrow \mathcal{N}\left(0, \frac{\sigma^2}{n}\right)
\nonumber
$$
If $\mathcal{I}(\theta)$ is the covariance of $U(\theta)$, then $n^2 \mathcal{I}(\theta)$ is the covariance of $\frac{1}{n}U(\theta)$. Thus:
$$
\sqrt{n} \left(\frac{1}{n} U^\top(\theta)\right) \rightsquigarrow \mathcal{N}\left(\mathbf{0}, \mathcal{I}(\theta) \right)
\hspace{4mm} \iff \hspace{4mm}
\sqrt{n} \left(\frac{1}{n} U^\top(\theta)\right)\mathcal{I}^{-1/2}(\theta) \rightsquigarrow \mathcal{N}\left(\mathbf{0}, \mathbb{I}\right)
\nonumber
$$
This can then be used to derive the asymptotic null distribution of the score test statistic, since the distribution of a squared Gaussian random variable is $\chi^2$:
$$
\left( \sqrt{n} \left(\frac{1}{n} U^\top(\theta)\right)\mathcal{I}^{-1/2}(\theta)\right)^2 \rightsquigarrow \chi^2_k
\hspace{4mm} \iff \hspace{4mm}
\frac{1}{n} U^\top(\theta) \mathcal{I}^{-1}(\theta) U(\theta) \rightsquigarrow \chi^2_k
$$
where $k$ is the dimension of $\theta$.
</details>

The alternative for the above test is implicitly two-sided: $H_1: \theta \neq \theta_0$.


<!-- ---
## All About Cones (and some linear algebra)

A _cone_, $C$, in $\mathbb{R}^m$ is defined as the set $$C := \{ t\mathbf{x} \in C \rvert \mathbf{x} \in C \}$$ for any $t > 0$. $C$ is called a _pointed cone_ if $\mathbf{x} \in C$ and $-\mathbf{x} \in C$ implies $\mathbf{x} = \mathbf{0}$, the zero vector. 

The _cone generated_ by $S$ for non-empty subset, $S$, of some real vector space, $X$, is defined as $$\text{cone}(S) := \{ a \mathbf{x} \rvert a \geq 0, \mathbf{x} \in S \}$$. We'll denote the _closure_ of $S$ with $\text{cl}(S)$, which is the set of all points in $S$ together with all of the limit points of $S$.


Let $$\bar{\mathbf{x}} \in \text{cl}(S)$$. A _tangent vector_ to $S$ at $\bar{\mathbf{x}}$ is any vector $\mathbf{h} \in X$ such that there exists a sequence $$(\mathbf{x}_n)_{n \in \mathbb{N}}$$ with $$\mathbf{x}_n \in S$$ and a sequence $$(a_n)_{n \in \mathbf{N}}$$ with $$a_n > 0$$ satisfying:

$$
\bar{\mathbf{x}} = \underset{n \rightarrow \infty}{\lim} x_n \hspace{8mm} \text{ and } \hspace{8mm}  \mathbf{h} = \underset{n \rightarrow \infty}{\lim} a_n(x_n - \bar{\mathbf{x}})
\nonumber
$$

We can then define the _sequential Bouligand tangent cone_ to $S$ at $\bar{\mathbf{x}}$, $T(S, \bar{\mathbf{x}})$, as the set of all such tangent vectors. This is also called the _contingent cone_ to $S$ at $\bar{\mathbf{x}}$.

<details>
<summary>Another Cone.</summary>
The <i>sequential Clarke tangent cone</i> is defined as the set of vectors, $\mathbf{h} \in X$, such that, for every sequence $(\mathbf{x}_n)_{n in \mathbb{N}}$ with $\mathbf{x}_n \in S$ satisfying $\bar{\mathbf{x}} = \underset{n \rightarrow \infty}{\lim} \mathbf{x}_n$ and every sequences of reals $(a_n)_{n \in \mathbb{N}}$ such that $a_n > 0$ and $0 = \underset{n \rightarrow \infty}{\lim} a_n$, there exists some sequence of vectors $(\mathbf{h}_n)_{n \in \mathbb{N}}$ such that $\mathbf{h} = \underset{n \rightarrow \infty}{\lim} \mathbf{h}_n$ and $\mathbf{x}_n + a_n \mathbf{h}_n \in S$ for all $n \in \mathbb{N}$.

The Clarke tangent cone is closed and convex, always. Furthermore, the Clarke tangent cone is contained within the Bouligand tangent cone when $\bar{\mathbf{x}} \in S$, and the two are the same if $X$ is convex. 
</details>

In geometry, if we assume that $S$ is closed and convex, then the _solid tangent cone_ to $S$ at $\bar{\mathbf{x}} \in b(S)$ (where $b(S)$ is the boundary of $S$) is the closure, $\text{cl}(C)$, of the cone, $C$, made of all rays that start at $\bar{\mathbf{x}}$ and intersect with $S$ at at least one $\mathbf{y} \in S$ where $\mathbf{y} \neq \bar{\mathbf{x}}$. That is, the set $$\{ \bar{\mathbf{x}} + t \mathbf{y} \rvert \mathbf{y} \in K, t \in [0, +\infty) \}$$. -->

---

## Introduction

In our setting, it is important to consider a one-sided alternative because variances are non-negative. As pointed out by Hall and Praestgaard (2001)[^fn-hall], omnibus tests like that of Lin (1997)[^fn-lin], which implicitly test against a two-sided alternative, have power against these impossible cases. They also explain that the claim that Lin's global test for homogeneity is locally asymptotically most stringest does not hold in particular boundary cases. Thus, it seems worthwhile to pursue explicitly one-sided tests, most of which are related to cones. 

### All About Cones
As an introduction into the theoretical/geometric setting we are interested in, I will cover some of the background and results in Shapiro[^fn-shapiro] This will give us a better foundation for what's to come.

Shapiro restricts his attention to 
<span class="popup" onclick="PopupFunc('pop1')">
    closed
    <span class="popuptext" id="pop1">
     A set is <i>closed</i> if it contains all of its boundary points.
    </span>
</span>
and 
<span class="popup" onclick="PopupFunc('pop2')">
    convex
    <span class="popuptext" id="pop2">
    A set, $S$, is <i>convex</i> if $tx + (1-t)y \in S$ for $x,y \in S$ and $t \in [0, 1]$.
    </span>
</span>
cones, so I'll do the same for the rest of this section.

A _cone_, $C$, in $\mathbb{R}^m$ is defined as the set $$C := \{ t\mathbf{x} \in C \rvert \mathbf{x} \in C \}$$ for any $t > 0$. $C$ is called a _pointed cone_ if $\mathbf{x} \in C$ and $-\mathbf{x} \in C$ implies $\mathbf{x} = \mathbf{0}$, the zero vector. 

Shapiro denotes the orthogonal projection of a point onto $C$ with 

$$
P(\mathbf{x}, C) = \underset{\eta \in C}{\arg\min} \left\{ (\mathbf{x} - \eta)^\top \mathbf{U} (\mathbf{x} - \eta) \right\}
$$ 

where $\mathbf{U}$ is any positive-definite matrix. This orthogonal projection maps the input point, $\mathbf{x}$, to the closest point on $C$, $\eta$, where closeness is defined by the norm associated with the matrix $\mathbf{U}$. We'll denote this norm with $\rvert \rvert \mathbf{x} \rvert \rvert = \sqrt{\mathbf{x}^\top \mathbf{U} \mathbf{x}}$, and we'll use $\langle\mathbf{x}, \mathbf{y}\rangle = \mathbf{x}^\top \mathbf{U}\mathbf{y}$ to denote the inner product of $\mathbf{x}$ and $\mathbf{y}$ associated with $\mathbf{U}$.

This brings us to the <i>dual</i> cone, which I have very little intuition...

<div id="dual-cone"></div>
<div class="definition">
  <body>
  <strong>Definition (Dual Cone).</strong>
  <br>
 The <i>dual cone</i> is the set $C^0= \{ \mathbf{y} \rvert \langle \mathbf{x}, \mathbf{y} \rangle \leq 0 \hspace{2mm} \forall \mathbf{x} \in C \}$. If $C$ is a vector subspace, then $C^0$ is the orthogonal complement of $C$, and if $C$ is closed and convex, then $(C^0)^0 = C$. 
  </body>
</div>

If we have another convex cone $K$, then:

$$
\rvert \rvert \mathbf{x} - P(\mathbf{x}, C) \rvert \rvert^2 = \rvert \rvert \mathbf{x} - P(\mathbf{x}, K) \rvert \rvert^2 + \rvert \rvert P(\mathbf{x}, K) - P(\mathbf{x}, C) \rvert \rvert^2
\nonumber
$$

if $C$ or $K$ is a lienar space and $C \subset K$. 

Furthermore, if $C$ is a vector subspace, then $\mathbf{x} - P(\mathbf{x}, C) = P(\mathbf{x}, C^0)$. That is, the difference between $\mathbf{x}$ and the orthogonal projection of $\mathbf{x}$ onto $C$ is equivalent to its orthogonal projection onto the dual cone $C^0$.

<div id="chi-bar-squared"></div>
<div class="definition">
  <body>
  <strong>Definition ($\bar{\chi}^2$-Statistic).</strong>
  <br>
  Let $\mathbf{y} \sim \mathcal{N}(\mathbf{0}, \mathbf{V})$ be a Gaussian random vector of $m$ dimensions with some covariance matrix, $\mathbf{V}$, and let $C$ be a convex cone. A <i>$\bar{\chi}^2$-statistic</i> is given by the following:
  $$
  \begin{aligned}
  \bar{\chi}^2 &= \mathbf{y}^\top \mathbf{V}^{-1} \mathbf{y} - \underset{\eta \in C}{\min} \left\{ (\mathbf{y} - \eta)^\top \mathbf{V}^{-1}(\mathbf{y} - \eta) \right\} \\
  &= \mathbf{y}^\top \mathbf{V}^{-1} \mathbf{y} - (\mathbf{y} - P(\mathbf{y}, C))^\top \mathbf{V}^{-1}(\mathbf{y} - P(\mathbf{y}, C)) \\
  &= (\mathbf{y} - \mathbf{y} + P(\mathbf{y}, C))^\top \mathbf{V}^{-1}(\mathbf{y} - \mathbf{y} + P(\mathbf{y}, C)) \\
  &= \rvert \rvert P(\mathbf{y}, C) \rvert \rvert^2
  \end{aligned}
  \nonumber
  $$
  where in the above, the inner products/norms are taken using the matrix $\mathbf{V}^{-1}$.
  </body>
</div>

The chi-bar-squared statistic follows a mixture of $\chi^2$ distributions:

$$
\mathbb{P}(\bar{\chi}^2 \geq c) = \sum_{i = 1}^m w_i \mathbb{P}(\chi_i^2 \geq c)
\nonumber
$$

where $\chi_i^2$ is a $\chi^2$ random variable with $i$ degrees of freedom, and $w_i$ are individuals weights that sum to $1$. As is standard, we let $\chi^2_0$ be a point mass at $0$. We'll denote this mixture distribution as $\mathcal{\bar{X}}^2(\mathbf{V}, C)$, since it depends on both $\mathbf{V}$ and the cone, $C$.

Note that, if we have the dual cone to $C$, then we can ignore the first term in the test statistic equation to get:

$$
\bar{\chi}^2 = \underset{\eta \in D}{\min} \left\{ (\mathbf{y} - \eta)^\top \mathbf{V}^{-1} (\mathbf{y} - \eta) \right\}
\nonumber
$$

which follows a $\mathcal{\bar{X}}^2(\mathbf{V}, C^0)$ distribution. We can also note that the mixture's weights satisfy:

$$
w_i(m, \mathbf{V}, C^0) = w_{m - i}(m, \mathbf{V}, C) \hspace{15mm} i = 0, \dots, m
\nonumber
$$


<!-- #### Other Stuff
Suppose $C$ is in a vector space, $L$, that is generated by the vectors in an $m \times k$ matrix, $\Delta$, that has rank $k$, and let $K$ be some convex cone in $\mathbb{R}^k$ such that we can write $C$ as:

$$
C = \Delta K = \{ \eta \rvert \eta = \Delta \beta, \beta \in K \}
\nonumber
$$

and therefore we can rewrite the chi-bar-squared statistic as:

$$
\mathcal{\bar{X}}^2 = \underset{\beta \in K}{\min} \left\{ (\mathbf{y} - \Delta \beta)^\top \mathbf{V}^{-1} (\mathbf{y} - \Delta \beta) \right\}
\nonumber
$$ -->

---

## A Multivariate One-Sided Test
Some of the later literature is based upon mid-century work by Kudô (1963)[^fn-kudo] for testing the mean vector of a multivariate Gaussian distribution. Kudô's work has some very nice geometric interpretations that permit the derivation of the limiting distribution of his test statistic, and all of it begins with the likelihood ratio. 

### Set-Up
We have a multivariate Gaussian population with mean vector $\theta^\top = (\theta_1, \dots, \theta_k)$ and known, non-singular, positive definite variance-covariance matrix $\Sigma$. Given a sample of $n$ i.i.d. observations $\mathbf{X} = \{ \mathbf{X}^{(1)}, \dots, \mathbf{X}^{(n)} \}$, we wish to test:

$$
H_0: \theta_i = 0 \hspace{2mm} i = 1, 2, \dots, k
\hspace{8mm} \text{vs.} \hspace{8mm}
H_1: \theta_i \geq 0 \hspace{2mm} i = 1, 2, \dots, k
$$

where at least one inequality in the alternative hypothesis setting is strict. 

Letting $\bar{\mathbf{X}} = \frac{1}{n} \sum_{i = 1}^n \mathbf{X}^{(i)}$, the sample mean vector, we can rewrite the likelihood of the sample under $H_1$ as:

$$
\begin{aligned}
\mathcal{L}(\theta, \Sigma; \mathbf{X}) 
&= \frac{1}{(\sqrt{2 \pi})^{kn} \rvert \Sigma \rvert^n} \exp\left(- \frac{1}{2} \sum_{i = 1}^n (\mathbf{X}^{(i)} - \theta)^\top \Sigma^{-1} (\mathbf{X}^{(i)} - \theta) \right)  \\
&= \frac{1}{(\sqrt{2 \pi})^{kn} \rvert \Sigma \rvert^n} \exp\left(- \frac{1}{2} \sum_{i = 1}^n \left[ (\mathbf{X}^{(i)} - \bar{\mathbf{X}})^\top \Sigma^{-1} (\mathbf{X}^{(i)} - \bar{\mathbf{X}}) + n(\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right] \right)  \\
\end{aligned}
$$

The standard likelihood ratio test statistic is given by:

$$
\begin{aligned}
t_{LRT} &= \frac{\underset{\theta_i = 0 \\ i = 1, \dots, k}{\max} \{ \mathcal{L}(\theta, \Sigma; \mathbf{X}) \}}{\underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\max} \{ \mathcal{L}(\theta, \Sigma; \mathbf{X}) \}} 
= \frac{\exp\left(-\frac{1}{2} n \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}}\right)}{\underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\max} \{ \exp\left(-\frac{1}{2}n(\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta)\right) \}}
\end{aligned}
$$

Notice that the argument maximum of the above is equivalent to:

$$
\bar{\chi}^2 = n \left[ \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}} - \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\min} \left\{ (\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right\} \right] 
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\theta^* &= \frac{\exp\left(-\frac{1}{2} n \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}}\right)}{\underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\max} \{ \exp\left(-\frac{1}{2}n(\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta)\right) \}} \\
&= \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\min} \left\{  \frac{\exp\left(-\frac{1}{2} n \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}}\right)}{\exp\left(-\frac{1}{2}n(\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta)\right)} \right\} & \left(\text{maximize denom. = minimize quotient} \right) \\
&= \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\min} \left\{  \exp\left( -\frac{1}{2} n \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}} + \frac{1}{2}n(\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right) \right\} \\
&= \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\min} \left\{   -\frac{1}{2} n \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}} + \frac{1}{2}n(\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right\} & \left(\exp \text{ monotonic} \right)\\
&= n \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\min} \left\{ -\bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}} + (\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right\} \\
&= n \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\max} \left\{ \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}} - (\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right\} & \left( \text{minimize = maximize negative} \right)\\
&= n \left[ \bar{\mathbf{X}}^\top \Sigma^{-1}\bar{\mathbf{X}} - \underset{\theta_i \geq 0 \\ i = 1, \dots, k}{\arg\min} \left\{ (\bar{\mathbf{X}} - \theta)^\top \Sigma^{-1}(\bar{\mathbf{X}} - \theta) \right\} \right] & \left(\text{ maximize quantity = minimize positive subtraction} \right)
\end{aligned}
\nonumber
$$
</details>

so we can just look at this quantity. 


---

## A One-Sided Score Test
I'll next dip into the work of Silvapulle and Silvapulle[^fn-silvapulle], who present a score-type test statistic for one-side alternative hypotheses based upon estimating functions instead of the true score function.

### Set-Up
The authors present a fairly general setting where we do not assume we know the exact form of the distribution of the observations, only that they depend on some $k \times 1$-dimensional vector-valued parameter, $\theta$, that is partitioned into the nuisance parameters, $\lambda$, and the components of interest, $\psi$. We'll write the partitioned parameter vector as $(\lambda : \psi) = (\lambda^\top, \psi^\top)^\top$ where $\lambda$ is $(k - q) \times 1$ and $\psi$ is $q \times 1$.

We're interested in testing hypotheses of the form:

$$
H_0: \psi = \mathbf{0}_q \hspace{10mm} H_A: \psi \in \mathcal{C}
\nonumber
$$

where $\mathbf{0}_q$ is a $q$-dimensional vector of zeros and $\mathcal{C}$ is the $q$-dimensional Euclidean space, or a closed, convex cone in $q$-dimensional Euclidean space with its vertex at the origin. The latter case encompasses alternatives of the form $\psi \geq \mathbf{0}_q$, which is the alternative hypothesis I am interested in for variance component testing.

Define $\mathbf{D}(\lambda)$ as some (fixed) matrix-valued function of the nuisance parameter, $\lambda$, that is independent of some (potentially vector-valued) random variable, $\delta$. Also define the $q \times 1$ vector, $\mathbf{U}_0$, as some function of the data. 

Suppose under the sequence of alternatives $K_n: \psi = n^{-1/2}\delta$, $\mathbf{U}_0$ satisfies:

$$
\mathbf{U}_0 \rightsquigarrow \mathcal{N}(\delta, \mathbf{D}(\lambda))
\label{eq:U-condition}
$$ 

as our sample size $n \rightarrow \infty$. Notice that testing the alternative hypothesis that $\psi \geq 0$ is equivalent to testing $\delta \geq \mathbf{0}_q$, which implies that the null hypothesis is equivalent to $\delta = \mathbf{0}_q$. 

---

### A General Test Statistic
The authors define a very general test statistic as the following.

<div class="definition">
  <body>
  <strong>Test Statistic.</strong>
  <br>
  Let $\tilde{\mathbf{D}}(\lambda)$ be any consistent estimator under the null hypothesis of $\mathbf{D}(\lambda)$, the asymptotic covariance matrix of $\mathbf{U}_0$. 
  <br>
  The test statistic for $H_0: \psi = \mathbf{0}$ against $H_A: \psi \in \mathcal{C}$ has the form:
  $$
  T = \mathbf{U}_0^\top \tilde{\mathbf{D}}(\lambda)^{-1}\mathbf{U}_0 - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{U}_0 - \mathbf{b})^\top \tilde{\mathbf{D}}(\lambda)^{-1}(\mathbf{U}_0 - \mathbf{b})\right\}
  \label{eq:test-stat-1}
  $$
  </body>
</div>

A $p$-value for large sample sizes can be found by defining $\mathbf{Z} \sim \mathcal{N}(\mathbf{0}, \mathbf{D}(\lambda))$ and:

$$
\xi(t, \mathbf{D}(\lambda), \mathcal{C}) = \mathbb{P}\left( \left[ \mathbf{Z}^\top \mathbf{D}(\lambda)^{-1} \mathbf{Z} - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{Z} - \mathbf{b})^\top \mathbf{D}(\lambda)^{-1}(\mathbf{Z} - \mathbf{b}) \right\} \right] \geq t \right)
\label{eq:xi-defn}
$$

The quantity $1 - \xi(t, \mathbf{D}(\lambda), \mathcal{C})$ follows a chi-bar-squared distribution; that is, a mixture of chi-squared distributions as we introduced in the previous section.[^fn-shapiro] The weights for the mixture can be hard to find, but we can get around this using the fact that, for large enough $n$ and under $H_0$ (i.e. $\delta = \mathbf{0}$), $\mathbf{U}_0$ is approximately $\mathcal{N}(\mathbf{0}, \mathbf{D}(\lambda))$. Thus, $\mathbb{P}(T \geq t; \lambda) \approx \xi(t, \mathbf{D}(\lambda), \mathcal{C})$. 

Suppose we observe a value of $T$, $t^*$. Define $$\mathbf{D}^*(\lambda)$$ as a consistent estimator of $\mathbf{D}(\lambda)$ for any $\lambda$. Then it follows that:

$$
p \approx \underset{\lambda}{\sup} \left\{ \xi(t^*, \mathbf{D}^*(\lambda), \mathcal{C}) \right\}
\nonumber
$$

for large enough $n$ because $\lambda$ is a nuisance parameter, so we can take the "best" probability over all of its values. 

---

### Use
How do we use this test statistic in practice? This is pretty much just a question of what function, $\mathbf{U}_0$, of our data we want to pick. Silvapulle and Silvapulle explain how to construct a score-type test statistic using their general statistic.  

Let's define $\mathbf{S}_n(\theta)$ as any $k \times 1$ vector estimating equation (so it should have expectation zero) for $\theta$ (e.g. the score function or something else). We need this vector to satisfy a couple of conditions.

<div class="definition">
  <body>
  <strong>Conditions.</strong>
  <br>
  Suppose $\mathbf{S}_n(\theta)$ is such that there exist non-singular $\mathbf{G}(\theta)$ and $\mathbf{V}(\theta)$ satisfying for any $a > 0$:
  $$
  \frac{1}{\sqrt{n}}\mathbf{S}_n(\theta) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathbf{V}(\theta))
  \label{eq:condition-a1}
  $$
  and
  $$
  \underset{\rvert \rvert \mathbf{h} \rvert \rvert \leq a}{\sup} \left\{ \frac{1}{\sqrt{n}} \left( \mathbf{S}_n\left(\theta + \frac{1}{\sqrt{n}} \mathbf{h}\right) - \mathbf{S}_n(\theta) \right)  + \mathbf{G}(\theta) \mathbf{h} \right\} = o_p(1)
  \label{eq:condition-a2}
  $$
  where $o_p(1)$ is stochastic order notation for convergence in probability to $0$.
  </body>
</div>

The first condition basically states that $\mathbf{S}_n(\theta)$ is asymptotically Gaussian when suitably scaled. Let's take a closer look at the second condition (this will be kind of hand-wavy).

Suppose we fix $\mathbf{h}$. The directional derivative of $\mathbf{S}_n(\theta)$ at $\theta$ along $\mathbf{h}$ is given by the limit:

$$
\nabla_{\mathbf{h}} \mathbf{S}_n(\theta) = \underset{s \rightarrow 0}{\lim} \left[ \frac{\mathbf{S}_n(\theta + s \mathbf{h}) - \mathbf{S}_n(\theta)}{s \rvert \rvert \mathbf{h} \rvert \rvert} \right]
$$

Technically, this is the definition for a scalar function, but we can just use the above notation to mean the limits are taken element-wise to get the result for a vector-valued function.

If we let $s = \frac{1}{\sqrt{n}}$, then we can rewrite the above as the limit as $n \rightarrow \infty$:

$$
\nabla_{\mathbf{h}} \mathbf{S}_n(\theta) = \underset{n \rightarrow \infty}{\lim} \left[  \frac{\sqrt{n}}{\rvert \rvert \mathbf{h} \rvert \rvert} \left( \mathbf{S}_n \left(\theta + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\theta) \right) \right]
$$

Scaling by $- \frac{1}{n}$, we get:

$$
\nabla_{\mathbf{h}} \left[ -\frac{1}{n} \mathbf{S}_n(\theta)\right] = -\frac{1}{\rvert \rvert \mathbf{h} \rvert \rvert} \underset{n \rightarrow \infty}{\lim} \left[  \frac{1}{\sqrt{n} } \left( \mathbf{S}_n \left(\theta + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\theta) \right) \right]
$$

Recall that for differentiable functions, the directional derivative is equal to the dot product between the gradient and the normalized direction vector. That is:

$$
\begin{aligned}
&\nabla_{\mathbf{h}} \left[ -\frac{1}{n} \mathbf{S}_n(\theta)\right] = \frac{\partial}{\partial \theta} \left[ - \frac{1}{n} \mathbf{S}_n(\theta) \right] \cdot \frac{\mathbf{h}}{\rvert \rvert \mathbf{h}\rvert \rvert} \\
\implies &-\underset{n \rightarrow \infty}{\lim} \left[  \frac{1}{\sqrt{n} } \left( \mathbf{S}_n \left(\theta + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\theta) \right) \right] = \frac{\partial}{\partial \theta} \left[ - \frac{1}{n} \mathbf{S}_n(\theta) \right] \cdot \mathbf{h}
\end{aligned}
$$

Let's subtract $\mathbf{G}(\theta) \mathbf{h}$ from both sides.

$$
\begin{aligned}
&- \mathbf{G}(\theta) \mathbf{h} - \underset{n \rightarrow \infty}{\lim} \left[  \frac{1}{\sqrt{n} } \left( \mathbf{S}_n \left(\theta + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\theta) \right) \right] = \left( \frac{\partial}{\partial \theta} \left[ - \frac{1}{n} \mathbf{S}_n(\theta) \right] - \mathbf{G}(\theta) \right) \mathbf{h} \\
\implies & - \left(\underset{n \rightarrow \infty}{\lim} \left[  \frac{1}{\sqrt{n} } \left( \mathbf{S}_n \left(\theta + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\theta) \right) \right] + \mathbf{G}(\theta) \mathbf{h}  \right)  = \left( \frac{\partial}{\partial \theta} \left[ - \frac{1}{n} \mathbf{S}_n(\theta) \right] - \mathbf{G}(\theta) \right) \mathbf{h}
\end{aligned}
$$

Notice that both $\mathbf{G}(\theta)$ and $\mathbf{h}$ are independent of $n$, so we can take them inside the limit:

$$
- \left(\underset{n \rightarrow \infty}{\lim} \left[  \frac{1}{\sqrt{n} } \left( \mathbf{S}_n \left(\theta + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\theta) \right) + \mathbf{G}(\theta) \mathbf{h}  \right] \right)  = \left( \frac{\partial}{\partial \theta} \left[ - \frac{1}{n} \mathbf{S}_n(\theta) \right] - \mathbf{G}(\theta) \right) \mathbf{h}
$$

Condition \eqref{eq:condition-a2} basically implies that the lefthand side will approach $0$, which itself implies that $\mathbf{G}(\theta) = \frac{\partial}{\partial \theta} \left[ - \frac{1}{n} \mathbf{S}_n(\theta) \right]$ in the limit. 

Let's partition our vectors and matrices in the following ways:

$$
\mathbf{S}_n(\theta) =
\begin{bmatrix}
\mathbf{S}^\top_{n, \lambda}(\theta) & 
\mathbf{S}^\top_{n, \psi}(\theta)
\end{bmatrix}^\top
\hspace{2mm} \text{ and } \hspace{2mm}
\mathbf{G}(\theta) =
\begin{bmatrix}
\mathbf{G}_{\lambda, \lambda}(\theta) & 
\mathbf{G}_{\lambda, \psi}(\theta) \\ 
\mathbf{G}_{\psi, \lambda}(\theta) &
\mathbf{G}_{\psi, \psi}(\theta)
\end{bmatrix}
\hspace{2mm} \text{ and } \hspace{2mm}
\mathbf{V}(\theta) =
\begin{bmatrix}
\mathbf{V}_{\lambda, \lambda}(\theta) & 
\mathbf{V}_{\lambda, \psi}(\theta) \\ 
\mathbf{V}_{\psi, \lambda}(\theta) &
\mathbf{V}_{\psi, \psi}(\theta)
\end{bmatrix}
$$

Let $$\theta_0 = (\lambda : \mathbf{0})$$ denote the value of $\theta$ under the null hypothesis, and suppose the null is true. Define the quantities:

$$
\begin{aligned}
\mathbf{Z}_n(\theta_0) &= n^{-1/2} \left( \mathbf{S}_{n, \psi}(\theta_0) - \mathbf{G}_{\psi, \lambda}(\theta_0) \mathbf{G}_{\lambda, \lambda}^{-1}(\theta_0) \mathbf{S}_{n, \lambda}(\theta_0) \right) \\
\mathbf{C}(\theta_0) &= \mathbf{V}_{\psi, \psi}(\theta_0) - \mathbf{G}_{\psi, \lambda}(\theta_0) \mathbf{G}^{-1}_{\lambda, \lambda}(\theta_0) \mathbf{V}_{\lambda, \psi}(\theta_0) - \left( \mathbf{V}_{\psi, \lambda}(\theta_0) - \mathbf{G}_{\psi, \lambda}(\theta_0) \mathbf{G}_{\lambda, \lambda}^{-1}(\theta_0) \mathbf{V}_{\lambda, \lambda}(\theta_0) \right)\left(\mathbf{G}^{-1}_{\lambda, \lambda}(\theta_0)\right)^\top \mathbf{G}^\top_{\psi, \lambda}(\theta_0)
\end{aligned}
$$

Since the above condition (Eq. \eqref{eq:condition-a1}) is assumed to be satisfied, and $\mathbf{Z}_n$ is just a function of $\mathbf{S}_n$:

$$
n^{-1/2} \mathbf{S}_n(\theta_0) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathbf{V}(\theta_0)) 
\hspace{2mm} \implies \hspace{2mm}
\mathbf{Z}_n(\theta_0) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathbf{C}(\theta_0))
$$

Denote consistent estimators for $\mathbf{G}(\theta_0)$ and $\mathbf{V}(\theta_0)$ with $\tilde{\mathbf{G}}(\theta_0)$ and $\tilde{\mathbf{V}}(\theta_0)$, respectively. Furthermore, let $\tilde{\lambda}$ denote a "suitable" estimator for $\lambda$ (where suitable is not really specific, but examples are given in the text), and let $\tilde{\theta}_0 = (\tilde{\lambda} : \mathbf{0})$. Define:

$$
\begin{aligned}
\mathbf{G}^{\psi, \psi}(\theta) &= \left( \mathbf{G}_{\psi, \psi}(\theta) - \mathbf{G}_{\psi, \lambda}(\theta) \mathbf{G}^{-1}_{\lambda, \lambda}(\theta) \mathbf{G}_{\lambda, \psi}(\theta) \right)^{-1} \\
\mathbf{A}(\theta) &= \left( \mathbf{G}^\top(\theta) \mathbf{V}^{-1}(\theta)\mathbf{G}(\theta) \right)^{-1} \\
\tilde{\mathbf{Z}}_n(\tilde{\theta}_0) &= n^{-1/2} \left( \mathbf{S}_{n, \psi}(\tilde{\theta}_0) - \tilde{\mathbf{G}}_{\psi, \lambda}(\tilde{\theta}_0) \tilde{\mathbf{G}}_{\lambda, \lambda}^{-1}(\tilde{\theta}_0) \mathbf{S}_{n, \lambda}(\tilde{\theta}_0) \right) 
\end{aligned}
$$

Using these, define:

$$
\mathbf{U}(\tilde{\theta}_0) = \tilde{\mathbf{G}}^{\psi, \psi}(\theta_0) \tilde{\mathbf{Z}}_n(\tilde{\theta}_0)
$$

where $\tilde{\mathbf{G}}^{\psi, \psi}$ is a consistent estimator for $\mathbf{G}^{\psi, \psi}(\theta_0)$.

Let's partition $\mathbf{A}(\theta)$ in the same way that we did with $\mathbf{V}(\theta)$ and $\mathbf{G}(\theta)$. With some work, we can see that for a fixed $$\delta \in \mathcal{C}$$, $$\mathbf{U}(\tilde{\theta}_0) \rightsquigarrow \mathcal{N}(\delta, \mathbf{A}_{\psi, \psi}(\theta_0))$$ under the sequence of alternatives $H_n: \psi = n^{-1/2} \delta$ as we take $n \rightarrow \infty$. Thus, $\mathbf{U}$ is a function of the data satisfying the condition in Eq. \eqref{eq:U-condition} and can be used in our test statistic construction.

<div id="test-statistic"></div>
<div class="definition">
  <body>
  <strong>Definition (One-Sided Test Statistic).</strong>
  <br>
  The test statistic for $H_0: \psi = \mathbf{0}$ against $H_A: \psi \in \mathcal{C}$ is given by:
  $$
  T_s = \mathbf{U}^\top(\tilde{\theta}_0) \tilde{\mathbf{A}}_{\psi, \psi}^{-1}(\tilde{\theta}_0) \mathbf{U}(\tilde{\theta}_0) - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{U}(\tilde{\theta}_0) - \mathbf{b})^\top \tilde{\mathbf{A}}_{\psi, \psi}^{-1}(\tilde{\theta}_0) (\mathbf{U}(\tilde{\theta}_0) - \mathbf{b}) \right\}
  \label{eq:test-stat-2}
  $$
  where $\tilde{\mathbf{A}}_{\psi, \psi}^{-1}(\tilde{\theta}_0)$ is the partition of $\mathbf{A}_{\psi, \psi}^{-1}(\tilde{\theta}_0)$ corresponding to $(\psi, \psi)$ and constructed using $\tilde{\mathbf{G}}(\tilde{\theta}_0)$ and $\tilde{\mathbf{V}}(\tilde{\theta}_0)$ (I think...the authors never define $\tilde{\mathbf{A}}$).
  </body>
</div>

A large sample $p$-value can be obtained as:

$$
p \approx \underset{\lambda}{\sup} \left\{ \xi(t^*, \mathbf{A}_{\psi, \psi}(\theta_0), \mathcal{C}) \right\}
$$

where $\xi(\cdot, \cdot, \cdot)$ is defined as in Eq. \eqref{eq:xi-defn} and $t^*$ is the observed value of $T_s$.

<details>
<summary>A Long Example.</summary>
Let's try going through a simple example situation where we can apply this one-sided test statistic. We'll do a random intercepts model.
<br>
Our goal is to test:
$$
H_0: \tau^2 = 0 
\hspace{5mm} \text{vs.} \hspace{5mm} 
H_A: \tau^2 > 0
$$
Since $\tau^2$ is a scalar, the alternative hypothesis is equivalent to $\tau^2 \in \mathcal{C}$ where $\mathcal{C}$ is the non-negative ray in $\mathbb{R}$. We will choose $\mathbf{S}_n(\theta)$ to be the score, so we need to do some derivations.
<br>
We will use the maximum likelihood estimator (assuming $H_0$ is true) to get our estimates for the nuisance parameters. Let $\theta_0$ denote our parameter vector under the null hypothesis, and let $\tilde{\theta}_0$ denote this same vector with the nuisance parameters replaced with their MLEs. The influence function for a maximum likelihood estimator is the score function, and I am pretty sure that maximum likelihood estimators are regular, so we should have that:
$$
\mathbf{S}(\theta_0) \rightsquigarrow \mathcal{N}\left(\mathbf{0}, \mathcal{I}^{-1}(\theta) \right)
$$
where $$\mathbf{S}(\theta_0)$$ is the gradient of the log-likelihood of $\theta$ given $\mathbf{y}$ (the score) evaluated under the null, and $\mathcal{I}(\theta) = \mathbb{E}_{\theta_0}\left[ \mathbf{S}(\theta_0) \mathbf{S}^\top(\theta_0) \right]$, the Fisher information. For some more explanation, see my <a href="/posts/2025/06/17/estimation-theory.html">estimation theory post</a>. 
To match the notation in the previous sections, we have $\mathbf{D}(\lambda) = \mathcal{I}^{-1}(\theta_0)$, and $\mathbf{U}_0 = \mathbf{S}(\tilde{\theta}_0)$. Since we don't know the true values of the other components  We can then use a consistent estimate of the Fisher information, denote it by $\tilde{\mathcal{I}}(\tilde{\theta}_0)$, for $\tilde{\mathbf{D}}(\tilde{\theta}_0)$.
<!-- 
<div class="example">
  <body>
  <strong>Example (Set-Up).</strong>
  <br>
  We will have $N$ total observations divided into $k$ different clusters with each clustering having $n_i$ observations in it (i.e. $\sum_{i=1}^k n_i = N$).
  <br>
  Let $\mathbf{y}_i \in \mathbb{R}^{n_i}$ denote the vector of responses for cluster $i$, let $\alpha$ be our fixed intercept, and let $\beta \in \mathbb{R}^{k \times 1}$ denote the random intercepts vector. We will use $\mathbf{y}_{i,j}$ to denote the $j$-th response in cluster $i$. 
  <br>
  We assume that the random effects are i.i.d. Gaussian with mean zero and variance $\tau^2$. 
  $$
  \beta_{i} \overset{iid}{\sim} \mathcal{N}\left(0, \tau^2 \right), \hspace{2mm} \forall i \in [k]
  $$
  Define the linear predictors and means as follows:
  $$
  \eta_{i,j} = \alpha + \beta_i, \hspace{5mm} \mu_{i,j} = \exp(\eta_{i,j})
  $$
  We assume that the data are drawn:
  $$
  \mathbf{y}_{i,j} \rvert \beta_i \overset{iid}{\sim} \text{Poi}(\mu_{i,j})
  $$
  The (total) conditional log-likelihood can be written as:
  $$
  \ell(\alpha, \tau^2; \mathbf{y} \rvert \beta) = \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} \log(\mu_{i,j}) - \mu_{i,j} - \log(\mathbf{y}_{i,j}!) \right)
  $$
  However, this can cause some issues down the road, so we'll use the log quasi-likelihood instead:
  $$
  \ell_q(\mathbf{y}; \alpha, \tau^2 \rvert \beta) = \sum_{i = 1}^k \sum_{j =1}^{n_i} \left( \mathbf{y}_{i,j} \log(\mu_{i,j}) - \mathbf{y}_{i,j} \log(\mathbf{y}_{i,j}) - \mu_{i,j} + \mathbf{y}_{i,j} \right)
  $$
  </body>
</div> -->
We assume to have $N$ total observations coming from $k$ different clusters. The $i$-th cluster will have $n_i$ observations, and we'll later assume that $n_i = n$ for all $i \in [k]$. We will assume independent clusters and that observations from the same cluster are i.i.d. Gaussian with mean $(\alpha + \beta_i)\mathbf{1}_{n_i}$ and variance $\sigma^2$:
$$
\mathbf{y}_i \rvert \beta_i \sim \mathcal{N}\left( (\alpha + \beta_i)\mathbf{1}_{n_i}, \sigma^2 \mathbb{I}_{n_i \times n_i}\right)
$$
Integrating out the random effects, which we assume to be i.i.d. Gaussian with mean zero and variance $\tau^2$, we get the marginal distribution of the observations as:
$$
\mathbf{y}_i \sim \mathcal{N}\left( (\alpha + \beta_i)\mathbf{1}_{n_i}, \sigma^2 \mathbb{I}_{n_i \times n_i} + \tau^2 \mathbf{1}_{n_i}\mathbf{1}_{n_i}^\top \right)
$$
<!-- #### Marginal Likelihood
We need the <i>marginal</i> log-likelihood, which means we need to integrate out the random effects. This can be kind of difficult, so instead we will take the expectation (with respect to the random effects distribution) of a Taylor approximation of the conditional log-likelihood. For a deeper discussion, see <a href="/posts/2025/06/04/glmm">this post</a>. -->
<!-- 
<div class="example">
  <body>
  <strong>Example (Marginal Likelihood).</strong>
  <br>
  The Taylor expansion of the conditional log-likelihood centered at $\beta_0$ is:
  $$
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \approx \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\rvert_{\beta = \beta_0} + (\beta - \beta_0)^\top \left[ \frac{\partial \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta } \right] \bigg\rvert_{\beta = \beta_0} + \frac{1}{2} (\beta - \beta_0)^\top \left[ \frac{\partial^2 \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta \partial \beta^\top} \right] \bigg\rvert_{\beta = \beta_0} (\beta - \beta_0)
  $$
  <details>
  <summary>First Derivatives.</summary>
  $$
  \begin{aligned}
  \frac{\partial \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta_h}
  &= \frac{\partial}{\partial \beta_h} \left[ \exp\left( \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} \log(\mu_{i,j}) - \mu_{i,j} - \log(\mathbf{y}_{i,j}!) \right) \right) \right] \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \frac{\partial}{\partial \beta_h} \left[ \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} \log(\mu_{i,j}) - \mu_{i,j} - \log(\mathbf{y}_{i,j}!) \right)\right] \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \sum_{j = 1}^{n_h} \left( \mathbf{y}_{h,j} \frac{\partial}{\partial \beta_h} \left[ \alpha + \beta_h \right]  - \frac{\partial}{\partial \beta_h} \left[ \exp(\alpha + \beta_h) \right]\right) \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \sum_{j = 1}^{n_h} \left(\mathbf{y}_{h,j} - \exp(\alpha + \beta_h) \right) \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \sum_{j = 1}^{n_h} \left( \mathbf{y}_{h,j} - \mu_{h,j} \right)
  \end{aligned}
  \nonumber
  $$
  Thus, the gradient is:
  $$
  \frac{\partial \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta} 
  = \begin{bmatrix}
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)\sum_{j = 1}^{n_1} \left( \mathbf{y}_{1,j} - \mu_{1,j} \right)\\
  \vdots \\
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)\sum_{j = 1}^{n_k} \left( \mathbf{y}_{k,j} - \mu_{k,j} \right)
  \end{bmatrix}
  \nonumber
  $$
  </details>
  <details>
  <summary>Second Derivatives.</summary>
  First, the diagonal terms:
  $$
  \begin{aligned}
  \frac{\partial^2 \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta_h^2} 
  &= \frac{\partial}{\partial \beta_h} \left[  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \sum_{j = 1}^{n_h} \left( \mathbf{y}_{h,j} - \mu_{h,j} \right) \right] \\
  &=  \frac{\partial \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta_h} \sum_{j = 1}^{n_h} \left( \mathbf{y}_{h,j} - \mu_{h,j} \right) + \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \frac{\partial}{\partial \beta_h} \left[ \sum_{j = 1}^{n_h} \left( \mathbf{y}_{h,j} - \mu_{h,j} \right) \right] \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left(\sum_{j = 1}^{n_h} \left(\mathbf{y}_{h, j} - \mu_{h, j} \right) \right)^2 + \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \sum_{j = 1}^{n_h} \frac{\partial}{\partial \beta_h} \left[ \mathbf{y}_{h, j} - \exp(\alpha + \beta_h) \right] \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left(\sum_{j = 1}^{n_h} \left(\mathbf{y}_{h, j} - \mu_{h, j} \right) \right)^2
  - \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \sum_{j = 1}^{n_h} \exp(\alpha + \beta_h)  \\
  &= \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left[ \left(\sum_{j = 1}^{n_h} \left(\mathbf{y}_{h, j} - \mu_{h, j} \right) \right)^2
  -\sum_{j = 1}^{n_h} \mu_{h,j} \right] 
  \end{aligned} 
  \nonumber
  $$
  Next we do the off-diagonal terms:
  $$
  \begin{aligned}
  \frac{\partial^2 \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta_h \partial \beta_{h'}} 
  &= \frac{\partial}{\partial \beta_{h'}} \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)\left(\sum_{j = 1}^{n_h} \left(\mathbf{y}_{h, j} - \mu_{h, j} \right) \right) \right] \\
  &=  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left(\sum_{j = 1}^{n_{h'}} \left(\mathbf{y}_{h', j} - \mu_{h', j} \right) \right)\left(\sum_{j = 1}^{n_h} \left(\mathbf{y}_{h, j} - \mu_{h, j} \right) \right)
  \end{aligned}
  \nonumber
  $$
  Putting the two pieces together:
  $$
  \frac{\partial^2 \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \beta \partial \beta^\top}
  = \begin{bmatrix}
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left[ \left(\sum_{j = 1}^{n_1} \left(\mathbf{y}_{1, j} - \mu_{1, j} \right) \right)^2
  -\sum_{j = 1}^{n_1} \mu_{1,j} \right] 
  & \dots &  
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left(\sum_{j = 1}^{n_{1}} \left(\mathbf{y}_{1, j} - \mu_{1, j} \right) \right)\left(\sum_{j = 1}^{n_k} \left(\mathbf{y}_{k, j} - \mu_{k, j} \right) \right) \\
  \vdots & \ddots & \vdots \\
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left(\sum_{j = 1}^{n_{k}} \left(\mathbf{y}_{k, j} - \mu_{k, j} \right) \right)\left(\sum_{j = 1}^{n_1} \left(\mathbf{y}_{1, j} - \mu_{h, j} \right) \right) 
  & \dots & 
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \left[ \left(\sum_{j = 1}^{n_k} \left(\mathbf{y}_{k, j} - \mu_{k, j} \right) \right)^2
  -\sum_{j = 1}^{n_k} \mu_{k,j} \right] 
  \end{bmatrix}
  \nonumber
  $$
  </details>
  We will do our expansion about $\beta_0 = \mathbf{0}$, the value under the null hypothesis. 
  $$
  \begin{aligned}
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)
  &\approx 
  \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}}
  + \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \sum_{i = 1}^k \beta_i \left(\sum_{j = 1}^{n_i} \left(\mathbf{y}_{i, j} - \mu_{i,j}^0 \right) \right) \\
  &\hspace{5mm} + \frac{1}{2} \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \sum_{i = 1}^k \left[ \beta_i^2  \left[\left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i, j} - \mu^0_{i,j} \right) \right)^2 - \sum_{j = 1}^{n_i} \mu^0_{i, j} \right] + \sum_{i' \neq i} \beta_i \beta_{i'} \left(\sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu^0_{i,j} \right)\right) \left(\sum_{j = 1}^{n_{i'}} \left(\mathbf{y}_{i',j} - \mu^0_{i',j} \right)\right) \right]
  \end{aligned}
  $$
  Since $\beta \sim \mathcal{N}(\mathbf{0}, \tau^2 \mathbb{I}_{k \times k})$, taking the expectation w.r.t. $\beta$ gives us:
  $$
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y})
  \approx 
  \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \left( 1 + \frac{\tau^2}{2} \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right)
  $$
  <details>
  <summary>Derivations.</summary>
  $$
  \begin{aligned}
  \mathcal{L}_q(\alpha, \tau^2; \mathbf{y})
  &= \mathbb{E}_\beta \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right] \\
  &\approx 
  \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}}
  + \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \sum_{i = 1}^k \mathbb{E}_\beta \left[\beta_i\right] \left(\sum_{j = 1}^{n_i} \left(\mathbf{y}_{i, j} - \mu_{i,j}^0 \right) \right) \\
  &\hspace{5mm} + \frac{1}{2} \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \sum_{i = 1}^k \left[ \mathbb{E}\left[\beta_i^2\right]  \left[\left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i, j} - \mu^0_{i,j} \right) \right)^2 - \sum_{j = 1}^{n_i} \mu^0_{i, j} \right] + \sum_{i' \neq i} \mathbb{E}\left[\beta_i\right] \mathbb{E}\left[\beta_{i'}\right] \left(\sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu^0_{i,j} \right)\right) \left(\sum_{j = 1}^{n_{i'}} \left(\mathbf{y}_{i',j} - \mu^0_{i',j} \right)\right) \right] \\
  &=  \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \left( 1 + \frac{1}{2} \sum_{i = 1}^k \tau^2 \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right) \\
  &= \left[ \mathcal{L}_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} \left( 1 + \frac{\tau^2}{2} \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right)
  \end{aligned}
  \nonumber
  $$
  </details>
  And taking the natural logarithm gives us the marginal log quasi-likelihood:
  $$
  \ell_q(\alpha, \tau^2; \mathbf{y})
  \approx
  \left[ \ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} + \log\left(1 + \frac{\tau^2}{2} \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right)
  $$
  Since $\log(1 + x) \approx x$ for small enough $x$:
  $$
  \ell_q(\alpha, \tau^2; \mathbf{y}) 
  \approx
  \left[ \ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} +  \frac{\tau^2}{2}\sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] 
  $$
  </body>
</div> -->
<!-- 
#### Score
We now need to derive the score vector, $\mathbf{S}(\theta)$. This is straightforward; we just take the gradient with respect to $\theta = (\alpha, \tau^2)$.  -->
<!-- 
<div class="example">
  <body>
  <strong>Example (Score).</strong>
  <br>
  We'll do the differentiation component-by-component. 
  $$
  \begin{aligned}
  \frac{\partial \ell_q(\alpha, \tau^2; \mathbf{y})}{\partial \alpha}
  &\approx 
  \underbrace{\frac{\partial}{\partial \alpha} \left[ [\ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)]\bigg\rvert_{\beta = \mathbf{0}} \right]}_{(a)}
  + \underbrace{\frac{\partial}{\partial \alpha} \left[ \frac{\tau^2}{2} \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right]}_{(b)} \\
  &= \left( \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} - \mu^0_{i,j} \right) \right)
  -\frac{\tau^2}{2} \sum_{i = 1}^k \left(\sum_{j = 1}^{n_i} \mu_{i,j}^0 \right) \left[2 \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0 \right) \right) + 1\right]
  \end{aligned}
  \nonumber
  $$
  <details>
  <summary>Derivations Of $(a)$ And $(b)$.</summary>
  $$
  \begin{aligned}
  (a)
  &= \sum_{i = 1}^k \sum_{j = 1}^{n_i} \frac{\partial}{\partial \alpha} \left[ \mathbf{y}_{i,j} \log(\mu^0_{i,j}) - \mathbf{y}_{i,j} \log(\mathbf{y}_{i,j}) - \mu^0_{i,j} + \mathbf{y}_{i,j} \right] \\
  &= \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} \frac{\partial}{\partial \alpha} \left[ \alpha \right] - \frac{\partial}{\partial \alpha} \left[ \exp(\alpha) \right] \right) \\
  &= \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} - \exp(\alpha) \right) \\
  &= \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left( \mathbf{y}_{i,j} - \mu^0_{i,j} \right)
  \end{aligned}
  \nonumber
  $$
  $$
  \begin{aligned}
  (b)
  &= \frac{\tau^2}{2} \sum_{i = 1}^k  \frac{\partial}{\partial \alpha} \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0  \right] \\
  &= \frac{\tau^2}{2} \sum_{i = 1}^k \left[ 2 \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0 \right) \right)\left(\sum_{j = 1}^{n_i} \frac{\partial}{\partial \alpha} \left[\mathbf{y}_{i,j} - \mu^0_{i,j}\right] \right)- \sum_{j = 1}^{n_i} \frac{\partial}{\partial \alpha} [\mu^0_{i,j}] \right] \\
  &= \frac{\tau^2}{2} \sum_{i = 1}^k \left[ 2 \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0 \right) \right)\left(\sum_{j = 1}^{n_i}  - \mu^0_{i,j} \right) - \sum_{j = 1}^{n_i} \mu^0_{i,j} \right]  \\
  &= -\frac{\tau^2}{2} \sum_{i = 1}^k \left(\sum_{j = 1}^{n_i} \mu_{i,j}^0 \right) \left[2 \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0 \right) \right) + 1\right]
  \end{aligned}
  \nonumber
  $$
  </details>
  $$
  \begin{aligned}
  \frac{\partial \ell_q(\alpha, \tau^2; \mathbf{y})}{\partial \tau^2}
  &\approx \frac{1}{2} \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right]
  \end{aligned}
  \nonumber
  $$
  </body>
</div> -->
<!-- #### Information
This become a bit tricky here as we need to compute the variance of the score under the null hypothesis. Unfortunately, we can't use the negative expected Hessian because there is no $\tau^2$ term in the second component of the score, so that would be zero. We'll have to take the variance directly. Luckily, the score has expectation zero, so we don't have to do any centering.
<div class="example">
  <body>
  <strong>Example (Information).</strong>
  <br>
  $$
  \begin{aligned}
  \mathcal{I}_{\alpha, \alpha} 
  &= \mathbb{E}\left[ \frac{\partial \ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \alpha}\frac{\partial \ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta)}{\partial \alpha} \right] \\
  &\approx \mathbb{E}\left[ \left( \left[ \ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) \right]\bigg\rvert_{\beta = \mathbf{0}} +  \frac{\tau^2}{2}\sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right)^2 \right] \\
  &= \underbrace{\mathbb{E}\left[ \left( [ \ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) ] \bigg\rvert_{\beta = \mathbf{0}} \right)^2 \right]}_{(a)} \\
  &\hspace{5mm} + \frac{2 \tau^2}{2} \underbrace{\mathbb{E}\left[ [\ell_q(\alpha, \tau^2; \mathbf{y} \rvert \beta) ]\bigg\rvert_{\beta = \mathbf{0}} \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right]}_{(b)} \\
  &\hspace{5mm} + \frac{(\tau^2)^2}{4} \underbrace{\mathbb{E}\left[ \left( \sum_{i = 1}^k \left[ \left( \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} - \mu_{i,j}^0\right) \right)^2 - \sum_{j = 1}^{n_i} \mu_{i,j}^0 \right] \right)^2 \right]}_{(c)}
  \end{aligned}
  \nonumber
  $$
  $$
  \begin{aligned}
  (a) &= \mathbb{E}\left[ \left( \sum_{i = 1}^k \sum_{j = 1}^{n_i} \left(\mathbf{y}_{i,j} \log(\mu_{i,j}^0) - \mathbf{y}_{i,j} \log(\mathbf{y}_{i,j}) - \mu^0_{i,j} + \mathbf{y}_{i,j} \right) \right)^2\right]
  \end{aligned}
  \nonumber
  $$
  </body>
</div>
 -->
The complete marginal likelihood is given by:
$$
\ell(\theta; \mathbf{y}) = -\frac{1}{2} \sum_{i = 1}^k \left[ n_i \log(2 \pi) + \log(\rvert \Sigma \rvert) + (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right]
$$
where $$\Sigma = \sigma^2 \mathbb{I}_{n_i \times n_i} + \tau^2 \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top$$.
Note that the Sherman-Morrison formula gives us:
$$
\Sigma^{-1} = \frac{1}{\sigma^2} \mathbb{I}_{n \times n} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n)} \mathbf{1}_{n}  \mathbf{1}_{n}^\top
$$
<details>
<summary>Proof.</summary>
$$
\begin{aligned}
  \Sigma^{-1}
  &= \left[ \sigma^2 \mathbb{I}_{n_i \times n_i} + \tau^2 \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right]^{-1} \\
  &= \left[ \sigma^2 \mathbb{I}_{n_i \times n_i} \right]^{-1} - \frac{\left[ \sigma^2 \mathbb{I}_{n_i \times n_i} \right]^{-1} \left(\tau^2 \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right) \left[ \sigma^2 \mathbb{I}_{n_i \times n_i} \right]^{-1}}{1 + \tau^2 \mathbf{1}_{n_i}^\top\left[ \sigma^2 \mathbb{I}_{n_i \times n_i} \right]^{-1} \mathbf{1}_{n_i} } \\
  &= \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} - \frac{\frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} \left(\tau^2 \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top\right) \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} }{1 + \tau^2 \mathbf{1}_{n_i}^\top \left(\frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} \right) \mathbf{1}_{n_i}} \\
  &= \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} - \frac{\left(\frac{\tau^2}{\sigma^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top\right) \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} }{1 + \frac{\tau^2}{\sigma^2} \mathbf{1}_{n_i}^\top  \mathbf{1}_{n_i}} \\
  &= \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} - \frac{\frac{\tau^2}{(\sigma^2)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top }{1 + \frac{\tau^2 n_i}{\sigma^2}} \\
  &= \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top
\end{aligned}
\nonumber
$$
</details>
We have:
$$
\mathbf{S}(\theta) = \frac{\partial}{\partial \theta} \left[ \ell(\theta; \mathbf{y})\right] 
= 
\begin{bmatrix}
\frac{\partial}{\partial \alpha} \left[ \ell(\theta; \mathbf{y})\right] \\
\frac{\partial}{\partial \sigma^2} \left[ \ell(\theta; \mathbf{y})\right] \\
\frac{\partial}{\partial \tau^2} \left[ \ell(\theta; \mathbf{y})\right] 
\end{bmatrix} 
=
\begin{bmatrix}
- \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j}) \\
- \frac{1}{2} \left[\frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)}  + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right] \\
- \frac{1}{2} \left[\frac{kn}{\sigma^2 + \tau^2 n} \right] +  \frac{1}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2 
\end{bmatrix}
$$
<details>
<summary>Proof.</summary>
$$
\frac{\partial}{\partial \tau^2} \left[ \ell(\theta; \mathbf{y})\right] 
= - \frac{1}{2} \left[\frac{kn}{\sigma^2 + \tau^2 n} \right] +  \frac{1}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2 
$$
<details>
<summary>Proof.</summary> 
  Using identities from <a href="https://en.wikipedia.org/wiki/Matrix_calculus#:~:text=Identities%3A%20scalar%2Dby%2Dscalar%2C%20with%20matrices%20involved">this table</a> on Wikipedia, we have:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \log(\rvert \Sigma \rvert) \right] 
    &= \text{tr} \left[ \Sigma^{-1} \frac{\partial \Sigma}{\partial \tau^2} \right] \\
    &= \text{tr} \left[ \Sigma^{-1} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right] \\
    &= \text{tr} \left[  \left( \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top  \right) \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right] \\
    &= \text{tr} \left[ \frac{1}{\sigma^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top - \frac{\tau^2 n_i}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right] \\
    &= \text{tr} \left[ \left( \frac{1}{\sigma^2} - \frac{\tau^2 n_i}{\sigma^2(\sigma^2 + \tau^2 n_i)} \right) \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right] \\
    &= \left( \frac{1}{\sigma^2} - \frac{\tau^2 n_i}{\sigma^2(\sigma^2 + \tau^2 n_i)} \right) \text{tr}\left[ \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right] \\
    &= \frac{n_i}{\sigma^2} - \frac{\tau^2 n_i^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \\
    &= \frac{n_i\sigma^2 + \tau^2 n_i^2 - \tau^2 n_i^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \\
    &= \frac{n_i}{\sigma^2 + \tau^2 n_i}
  \end{aligned}
  \nonumber
  $$
  $$
  \begin{aligned}
    \frac{\partial \Sigma^{-1}}{\partial \tau^2} 
    &= - \Sigma^{-1} \frac{\partial \Sigma}{\partial \tau^2} \Sigma^{-1} \\
    &= - \left( \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top \right) \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \left(\frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top\right) \\
    &= - \left( \frac{1}{\sigma^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top - \frac{\tau^2 n_i}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top\right)\left(\frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top\right) \\
    &= - \left( \left( \frac{1}{(\sigma^2)^2} - \frac{\tau^2 n_i}{(\sigma^2)^2 (\sigma^2 + \tau^2 n_i)}\right)  \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top - \left(\frac{\tau^2 n_i}{(\sigma^2)^2 (\sigma^2 + \tau^2 n_i)} - \frac{(\tau^2)^2 n_i^2}{(\sigma^2)^2 (\sigma^2 + \tau^2 n_i)^2} \right) \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right) \\
    &= -\left( \frac{1}{(\sigma^2)^2} - \frac{2\tau^2 n_i}{(\sigma^2)^2 (\sigma^2 + \tau^2 n_i)} + \frac{(\tau^2)^2 n_i^2}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2}\right) \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \\
    &= -\left(\frac{\tau^2 n_i}{\sigma^2(\sigma^2 + \tau^2 n_i)} - \frac{1}{\sigma^2} \right)^2 \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \\
    &= - \left(- \frac{\sigma^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \right)^2 \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \\
    &= -\frac{1}{(\sigma^2 + \tau^2 n_i)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top
  \end{aligned}
  $$
  The above equations can be used to derive:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right]
    &= \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \frac{\partial \Sigma^{-1}}{\tau^2} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \\
    &= \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \left(- \frac{1}{(\sigma^2 + \tau^2 n_i)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top \right) (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \\
    &= \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)(\mathbf{y}_{i, j'} - \alpha) \left(- \frac{1}{(\sigma^2 + \tau^2 n_i)^2}\right) \\
    &= -  \sum_{i = 1}^k \frac{1}{(\sigma^2 + \tau^2 n_i)^2}\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha) \sum_{j' = 1}^{n_i} (\mathbf{y}_{i, j'} - \alpha) \\
    &= -  \sum_{i = 1}^k  \frac{1}{(\sigma^2 + \tau^2 n_i)^2}\left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha) \right) \left(\sum_{j' = 1}^{n_i} (\mathbf{y}_{i, j'} - \alpha) \right) \\
    &= -  \sum_{i = 1}^k  \frac{1}{(\sigma^2 + \tau^2 n_i)^2} \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha) \right)^2
  \end{aligned}
  \nonumber
  $$
  Putting all of the above together, we get:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \ell(\alpha, \sigma^2, \tau^2; \mathbf{y})\right] 
    &= -\frac{1}{2} \sum_{i = 1}^k \left[  \frac{\partial}{\partial \tau^2} \left[ \log\left( \rvert \Sigma \rvert \right) \right] + \frac{\partial}{\partial \tau^2} \left[ (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right] \right] \\
    &= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{n_i}{\sigma^2 + \tau^2 n_i} \right] - \frac{1}{2} \frac{\partial}{\partial \tau^2} \left[ \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right] \\
    &= - \frac{1}{2} \sum_{i = 1}^k \left[\frac{n_i}{\sigma^2 + \tau^2 n_i} \right] +  \sum_{i = 1}^k \left[ \frac{1}{2(\sigma^2 + \tau^2 n_i)^2} \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha) \right)^2 \right]
  \end{aligned}
  \nonumber
  $$
  Since $n_i = n$ for all $i \in [k]$, this simplifies to:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \tau^2} \left[ \ell(\alpha, \sigma^2, \tau^2; \mathbf{y})\right] 
    &= - \frac{1}{2} \sum_{i = 1}^k \left[\frac{n}{\sigma^2 + \tau^2 n} \right] +  \sum_{i = 1}^k \left[ \frac{1}{2(\sigma^2 + \tau^2 n)^2} \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2 \right] \\
    &= - \frac{1}{2} \left[\frac{kn}{\sigma^2 + \tau^2 n} \right] +  \frac{1}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2  \\
  \end{aligned}
  \nonumber
  $$
</details>
$$
\frac{\partial}{\partial \sigma^2} \left[ \ell(\theta; \mathbf{y})\right] 
= - \frac{1}{2} \left[\frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)}  + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right]
$$
<details>
<summary>Proof.</summary>
  We repeat the process with $\sigma^2$.
  $$
  \begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \log(\rvert \Sigma \rvert) \right]
    &= \text{tr} \left[ \Sigma^{-1} \frac{\partial \Sigma}{\partial \sigma^2} \right] \\
    &=  \text{tr} \left[ \Sigma^{-1} \mathbb{I}_{n_i \times n_i} \right] \\
    &= \text{tr} \left[ \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top \right] \\
    &= \frac{n_i}{\sigma^2} - \frac{\tau^2 n_i}{\sigma^2(\sigma^2 + \tau^2 n_i)} \\
    &= \frac{n_i(\sigma^2 - \tau^2 + \tau^2 n_i)}{\sigma^2(\sigma^2 + \tau^2 n_i)}
  \end{aligned}
  $$
  $$
  \begin{aligned}
    \frac{\partial \Sigma^{-1}}{\partial \sigma^2} 
    &= - \Sigma^{-1} \frac{\partial \Sigma}{\partial \sigma^2} \Sigma^{-1} \\
    &= - \left( \frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top \right) \mathbb{I}_{n_i \times n_i} \left(\frac{1}{\sigma^2} \mathbb{I}_{n_i \times n_i} -\frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i}  \mathbf{1}_{n_i}^\top\right) \\
    &= - \left(\frac{1}{(\sigma^2)^2} \mathbb{I}_{n_i \times n_i} - \frac{2\tau^2}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top + \frac{(\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top\right) \\
    &= - \left(\frac{1}{(\sigma^2)^2} \mathbb{I}_{n_i \times n_i} - \left(\frac{2\tau^2(\sigma^2 + \tau^2 n_i)}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} - \frac{(\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \right) \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top\right) \\
    &= - \left(\frac{1}{(\sigma^2)^2} \mathbb{I}_{n_i \times n_i} - \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top\right) \\
    &= \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n_i \times n_i}
  \end{aligned}
  $$
  The above equations can be used to derive:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right]
    &= \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \frac{\partial \Sigma^{-1}}{\sigma^2} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \\
    &= \sum_{i = 1}^k (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \left( \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \mathbf{1}_{n_i} \mathbf{1}_{n_i}^\top - \frac{1}{(\sigma^2)^2} \mathbb{I}_{n_i \times n_i} \right) (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \\
    &= \sum_{i = 1}^k \left[ \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} -\alpha ) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2 \right]
  \end{aligned}
  $$
  Putting all of the above together, we get:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \ell(\theta; \mathbf{y})\right] 
    &= -\frac{1}{2} \sum_{i = 1}^k \left[  \frac{\partial}{\partial \sigma^2} \left[ \log\left( \rvert \Sigma \rvert \right) \right] + \frac{\partial}{\partial \sigma^2} \left[ (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right] \right] \\
    &= - \frac{1}{2}  \sum_{i = 1}^k \left[ \frac{n_i(\sigma^2 - \tau^2 + \tau^2 n_i)}{\sigma^2(\sigma^2 + \tau^2 n_i)} + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n_i}{(\sigma^2)^2(\sigma^2 + \tau^2 n_i)^2} \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} -\alpha ) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2 \right] 
  \end{aligned}
  $$
  Since $n_i = n$ for all $i \in [k]$, this simplifies to:
  $$
  \begin{aligned}
    \frac{\partial}{\partial \sigma^2} \left[ \ell(\theta; \mathbf{y})\right] 
    &= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{n(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2(\sigma^2 + \tau^2 n)} + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right] \\
    &= - \frac{1}{2} \left[\frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)}  + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right]
  \end{aligned}
  \nonumber
  $$
</details>
$$
\frac{\partial}{\partial \alpha} \left[ \ell(\theta; \mathbf{y})\right] 
= - \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j})
$$
<details>
<summary>Proof.</summary>
$$
\begin{aligned}
  \frac{\partial}{\partial \alpha} \left[\ell(\alpha, \sigma^2, \tau^2; \mathbf{y}) \right] 
  &= - \frac{1}{2} \sum_{i = 1}^k \frac{\partial}{\partial \alpha} \left[(\mathbf{y}_i - \alpha \mathbf{1}_{n_i})^\top \Sigma^{-1} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i}) \right] \\
  &= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{\partial}{\partial \alpha} \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})_{j} (\mathbf{y}_i - \alpha \mathbf{1}_{n_i})_{j'} \left[ \Sigma^{-1} \right]_{j, j'} \right] \\
  &= - \frac{1}{2} \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i}  \left[ \Sigma^{-1} \right]_{j, j'} \frac{\partial}{\partial \alpha}  \left[  (\mathbf{y}_{i,j} - \alpha) (\mathbf{y}_{i,j'} - \alpha) \right]\\
  &= - \frac{1}{2} \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i}  \left[ \Sigma^{-1} \right]_{j, j'} \left( -\mathbf{y}_{i,j} - \mathbf{y}_{i,j'} + 2 \alpha\right)\\
  &= - \frac{1}{2} \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} \left( (\alpha - \mathbf{y}_{i,j}) + (\alpha - \mathbf{y}_{i, j'}) \right) \\
  &= - \frac{1}{2} \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} (\alpha - \mathbf{y}_{i,j}) - \frac{1}{2} \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'}(\alpha - \mathbf{y}_{i, j'}) \\
  &=  - \sum_{i = 1}^k \sum_{j = 1}^{n_i} (\alpha - \mathbf{y}_{i,j}) \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} 
\end{aligned} 
$$
Substituting for $\Sigma^{-1}$ and simplifying due to the assumption that $n_i = n$:
$$
\begin{aligned}
\frac{\partial}{\partial \alpha} \left[\ell(\alpha, \sigma^2, \tau^2; \mathbf{y}) \right] 
&= - \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j}) \sum_{j' = 1}^{n} \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{\tau^2}{\sigma^2(\sigma^2 + \tau^2 n)} \mathbf{1}_n \mathbf{1}_n^\top \right]_{j, j'} \\
&= - \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j})  \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \\
&= - \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j})
\end{aligned}
$$
</details>
</details>
We will use the maximum likelihood estimates under $H_0$ as our "suitable" estimates for the nuisance parameters. Thus:
$$
\tilde{\theta}_0 
=
\begin{bmatrix}
\hat{\alpha} \\
\hat{\sigma}^2 \\
0
\end{bmatrix}
=
\begin{bmatrix}
\frac{1}{kn} \sum_{j = 1}^{n}  \sum_{i = 1}^k  \sum_{j = 1}^n \mathbf{y}_{i,j}  \\
\frac{1}{k n} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2\\
0
\end{bmatrix}
$$
<details>
<summary>Proof.</summary>
$$
\hat{\sigma}^2 = \frac{\sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2}{k n}
$$
<details>
<summary>Proof.</summary>
Under the null hypothesis, $\tau^2 = 0$. Plugging this into the score for $\sigma^2$:
$$
\begin{aligned}
  \frac{\partial}{\partial \sigma^2} \left[ \ell(\alpha, \sigma^2, \tau^2; \mathbf{y})\right]\bigg\rvert_{\tau^2 = 0} 
  &= - \frac{1}{2} \sum_{i = 1}^k \left[ \frac{n_i(\sigma^2 - 0 + 0)}{\sigma^2(\sigma^2 + 0)} + \frac{0 + 0}{(\sigma^2)^2(\sigma^2 + 0)^2} \left(\sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} -\alpha ) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2 \right] \\
  &= -\frac{1}{2} \sum_{i = 1}^k \left[ \frac{n_i}{\sigma^2} - \frac{1}{(\sigma^2)^2} \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2 \right]
\end{aligned}
$$
Equating the above with $0$ and solving for $\sigma^2$ yields:
$$
\begin{aligned}
  0 &= -\frac{1}{2} \sum_{i = 1}^k \left[ \frac{n_i}{\sigma^2} - \frac{1}{(\sigma^2)^2} \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2 \right] \\
  \frac{1}{\sigma^2}\sum_{i = 1}^k n_i &= \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2 \\
  \sigma^2 &= \frac{\sum_{i = 1}^k \sum_{j = 1}^{n_i} (\mathbf{y}_{i,j} - \alpha)^2}{\sum_{i = 1}^k n_i}
\end{aligned}
$$
Simplifying with $n_i = n$:
$$
\begin{aligned}
\hat{\sigma}^2 &= \frac{\sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2}{\sum_{i = 1}^k n} \\
&= \frac{\sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2}{k n}
\end{aligned}
$$
</details>
$$
\hat{\alpha} = \frac{1}{kn} \sum_{j = 1}^{n}  \sum_{i = 1}^k  \mathbf{y}_{i,j} 
$$
<details>
<summary>Proof.</summary>
Setting the derivative equal to $0$ and solving for $\alpha$:
$$
\begin{aligned}
0 &= \frac{\partial}{\partial \alpha} \left[ \ell(\theta; \mathbf{y})\right] \\
0 &= - \sum_{i = 1}^k \sum_{j = 1}^{n_i} (\alpha - \mathbf{y}_{i,j}) \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} \\
0 &= - \alpha \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} + \sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} \\
\alpha \sum_{i = 1}^k \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} &= \sum_{i = 1}^k \sum_{j = 1}^{n_i} \mathbf{y}_{i,j} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'}  \\
\alpha \left(k \mathbf{1}_{n_i}^\top \Sigma^{-1} \mathbf{1}_{n_i} \right) &= \sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} \sum_{i = 1}^k  \mathbf{y}_{i,j} \\
\alpha &= \frac{\sum_{j = 1}^{n_i} \sum_{j' = 1}^{n_i} \left[ \Sigma^{-1} \right]_{j, j'} \sum_{i = 1}^k  \mathbf{y}_{i,j}}{k \mathbf{1}_{n_i}^\top \Sigma^{-1} \mathbf{1}_{n_i} }
\end{aligned}
\nonumber
$$
Substituting in for $\Sigma^{-1}$ under the null and simplifying, since $n_i = n$:
$$
\begin{aligned}
\hat{\alpha} &= \frac{\sum_{j = 1}^{n} \sum_{j' = 1}^{n} \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{0}{\sigma^2 (\sigma^2 + 0)} \mathbf{1}_n \mathbf{1}_n^\top \right]_{j, j'} \sum_{i = 1}^k  \mathbf{y}_{i,j}}{k \mathbf{1}_{n}^\top \left( \frac{1}{\sigma^2} \mathbb{I}_{n \times n} - \frac{0}{\sigma^2 (\sigma^2 + 0)} \mathbf{1}_n \mathbf{1}_n^\top \right) \mathbf{1}_{n} }\\
&= \frac{\sum_{j = 1}^{n} \sum_{j' = 1}^{n} \left[ \frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right]_{j, j'} \sum_{i = 1}^k  \mathbf{y}_{i,j}}{\frac{kn}{\sigma^2}}\\
&= \frac{\sigma^2 \sum_{j = 1}^{n} \frac{1}{\sigma^2} \sum_{i = 1}^k  \mathbf{y}_{i,j}}{kn} \\
&= \frac{1}{kn} \sum_{j = 1}^{n}  \sum_{i = 1}^k  \mathbf{y}_{i,j} 
\end{aligned}
$$
</details>
</details>
Thus:
$$
\mathbf{S}(\tilde{\theta}_0) = 
\begin{bmatrix}
  \frac{1}{\hat{\sigma}^2} \sum_{i = 1}^k (\mathbf{y}_{i,j} - \hat{\alpha}) \\
  - \frac{kn}{2\hat{\sigma}^2} + \frac{1}{2(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \hat{\alpha})^2 \\
  - \frac{kn}{2 \hat{\sigma}^2} + \frac{1}{2(\hat{\sigma}^2)^2} \sum_{i = 1}^k \left(\sum_{j = 1}^n (\mathbf{y}_{i,j} - \hat{\alpha}) \right)^2
\end{bmatrix}
$$
We'll derive the components of the Hessian by each variable:
$$
\frac{\partial}{\partial \alpha} \left[ \frac{\partial}{\partial \theta}\left[\ell (\theta; \mathbf{y}) \right] \right] = 
\begin{bmatrix}
-\frac{nk}{\sigma^2} + \frac{k n^2 \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)} \\
- \frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \\
-\frac{n}{(\sigma^2 + \tau^2 n)^2}\sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)
\end{bmatrix}
$$
<details>
<summary>Second Derivatives w.r.t. $\alpha$.</summary>
$$
\begin{aligned}
  \frac{\partial^2}{\partial \alpha^2} \left[ \ell(\theta; \mathbf{y}) \right] 
&= \frac{\partial}{\partial \alpha} \left[ - \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j}) \right] \\
&= - \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} 1 \\
&= - nk \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \\
&= -\frac{nk}{\sigma^2} + \frac{k n^2 \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
  \frac{\partial^2}{\partial \alpha \partial \sigma^2} \left[ \ell(\theta; \mathbf{y}) \right] 
&= - \frac{1}{2} \left[\frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)}  + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right] \\
&= - \frac{1}{2} \left[\frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \frac{\partial}{\partial \alpha} \left[ \left(\sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right)^2 \right]  - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^n \frac{\partial}{\partial \alpha} \left[(\mathbf{y}_{i,j} - \alpha)^2 \right] \right] \\
&=  - \frac{1}{2} \left[\frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k 2\left(\sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right) \frac{\partial}{\partial \alpha} \left[ \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)  \right]   - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^n -2(\mathbf{y}_{i,j} - \alpha)  \right] \\ 
&=  - \frac{1}{2} \left[\frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k 2\left(\sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right)  \sum_{j = 1}^n (-1) + \frac{2}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)  \right] \\ 
&=  - \left[-\frac{2\tau^2 \sigma^2n + (\tau^2)^2 n^2}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)  + \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)  \right] \\ 
&=  \frac{2\tau^2 \sigma^2n + (\tau^2)^2 n^2}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) - \frac{(\sigma^2 + \tau^2 n )^2}{(\sigma^2)^2 (\sigma^2 + \tau^2 n )^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)  \\ 
&=  \frac{2\tau^2 \sigma^2n + (\tau^2)^2 n^2}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) - \frac{(\sigma^2)^2 + 2 \sigma^2 \tau^2 n + (\tau^2)^2 n^2}{(\sigma^2)^2 (\sigma^2 + \tau^2 n )^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)  \\ 
&=  \frac{-(\sigma^2)^2}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \\
&= - \frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \\
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
  \frac{\partial^2}{\partial \alpha \partial \tau^2} \left[ \ell(\theta; \mathbf{y}) \right] 
&= \frac{\partial}{\partial \alpha} \left[ - \frac{1}{2} \left[\frac{kn}{\sigma^2 + \tau^2 n} \right] +  \frac{1}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2  \right] \\
&= \frac{2}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j}- \alpha)\frac{\partial}{\partial \alpha} \left[ \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right] \\
&= \frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \left[ \sum_{j = 1}^n (-1) \right]\\
&= -\frac{n}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha)
\end{aligned}
\nonumber
$$
</details>
$$
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial}{\partial \theta}\left[\ell (\theta; \mathbf{y}) \right] \right] = 
\begin{bmatrix}
-\frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j}- \alpha) \\
- \frac{1}{2} \left( \left[\frac{k - kn}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2}\right] + \frac{2}{n} \left[\frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 + \frac{2}{(\sigma^2)^3} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2  \right) \\
\frac{kn}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2n)^3}\sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2
\end{bmatrix}
$$
<details>
<summary>Second Derivatives w.r.t. $\sigma^2$.</summary>
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial}{\partial \alpha} \ell(\theta; \mathbf{y}) \right]
&= \frac{\partial}{\partial \sigma^2} \left[ - \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j})  \right] \\
&= \frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j}) \\
&= -\frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
  \frac{\partial}{\partial \sigma^2} \left[ \frac{\partial}{\partial \sigma^2} \ell(\theta; \mathbf{y}) \right]
&=   \frac{\partial}{\partial \sigma^2} \left[ - \frac{1}{2} \left[\frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)}  + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right] \right] \\
&= - \frac{1}{2}\left( \frac{\partial}{\partial \sigma^2} \left[ \frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)} \right] + \frac{\partial}{\partial \sigma^2} \left[ \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2\right] - \frac{\partial}{\partial \sigma^2}\left[ \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right] \right) \\
&= - \frac{1}{2} \left( \left[\frac{k - kn}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2}\right] + \frac{2}{n} \left[\frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 + \frac{2}{(\sigma^2)^3} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2  \right)
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
\frac{\partial}{\partial \sigma^2} \left[ \frac{\partial}{\partial \tau^2} \ell(\theta; \mathbf{y}) \right]
&= \frac{\partial}{\partial \sigma^2} \left[ - \frac{1}{2} \left[\frac{kn}{\sigma^2 + \tau^2 n} \right] +  \frac{1}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2  \right] \\
&= - \frac{1}{2}\left(\frac{-kn}{(\sigma^2 + \tau^2 n)^2}\right) + \frac{1}{(\sigma^2 + \tau^2n)^3}\sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2 \\
&= \frac{kn}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2n)^3}\sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2
\end{aligned}
\nonumber
$$
</details>
$$
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial}{\partial \theta}\left[\ell (\theta; \mathbf{y}) \right] \right] = 
\begin{bmatrix}
-\frac{n}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \\
\frac{kn}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3}\sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 \\
\frac{kn^2}{2(\sigma^2 + \tau^2n)^2} - \frac{n}{(\sigma^2 + \tau^2 n)^3}\sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2
\end{bmatrix}
$$
<details>
<summary>Second Derivatives w.r.t. $\tau^2$.</summary>
$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial}{\partial \alpha} \ell(\theta; \mathbf{y}) \right] 
&= \frac{\partial}{\partial \tau^2} \left[ - \left( \frac{1}{\sigma^2} - \frac{n \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)}\right) \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j}) \right]  \\
&= \frac{n}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\alpha - \mathbf{y}_{i,j})
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial}{\partial \sigma^2} \ell(\theta; \mathbf{y}) \right] 
&= \frac{\partial}{\partial \tau^2} \left[ - \frac{1}{2} \left[\frac{kn(\sigma^2 - \tau^2 + \tau^2 n)}{\sigma^2( \sigma^2 + \tau^2 n)}  + \frac{2\tau^2 \sigma^2 + (\tau^2)^2 n}{(\sigma^2)^2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 - \frac{1}{(\sigma^2)^2} \sum_{i = 1}^k \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha)^2 \right] \right] \\
&= - \frac{1}{2}\left(
  - \frac{kn}{(\sigma^2 + \tau^2 n)^2} + \frac{2}{(\sigma^2 + \tau^2 n)^3}\sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 
\right) \\
&= \frac{kn}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3}\sum_{i = 1}^k  \left(\sum_{j = 1}^{n} (\mathbf{y}_{i,j} -\alpha) \right)^2 
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
\frac{\partial}{\partial \tau^2} \left[ \frac{\partial}{\partial \tau^2} \ell(\theta; \mathbf{y}) \right] 
&= \frac{\partial}{\partial \tau^2} \left[ - \frac{1}{2} \left[\frac{kn}{\sigma^2 + \tau^2 n} \right] +  \frac{1}{2(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2  \right] \\
&= - \frac{1}{2}\frac{-kn^2}{(\sigma^2 + \tau^2n)^2} - \frac{n}{(\sigma^2 + \tau^2 n)^3}\sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2 \\
&= \frac{kn^2}{2(\sigma^2 + \tau^2n)^2} - \frac{n}{(\sigma^2 + \tau^2 n)^3}\sum_{i = 1}^k  \left( \sum_{j = 1}^{n} (\mathbf{y}_{i,j} - \alpha) \right)^2
\end{aligned}
\nonumber
$$
</details>
Taking the expectation, we get:
$$
\mathbb{E}\left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \theta \partial \theta^\top} \right] =
\begin{bmatrix}
- \frac{nk}{\sigma^2} + \frac{k n^2 \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)} & 0 & \\ 
0 & - \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2nk(\sigma^2 + \tau^2 n)}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right]+ \frac{2 n k (\sigma^2 + \tau^2)}{(\sigma^2)^3} \right) & -\frac{nk}{2(\sigma^2 + \tau^2 n)^2} \\
0 & -\frac{nk}{2(\sigma^2 + \tau^2 n)^2} & - \frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2}
\end{bmatrix}
$$
<details>
<summary>Derivations.</summary>
$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \alpha \partial \theta} \right]
&= \begin{bmatrix}
- \frac{nk}{\sigma^2} + \frac{k n^2 \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)} \\
- \frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^n \sum_{j = 1}^n \mathbb{E}\left[(\mathbf{y}_{i,j} - \alpha) \right] \\
- \frac{n}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^n \sum_{j = 1}^n \mathbb{E}\left[(\mathbf{y}_{i,j} - \alpha) \right]
\end{bmatrix} \\
&= \begin{bmatrix}
- \frac{nk}{\sigma^2} + \frac{k n^2 \tau^2}{\sigma^2(\sigma^2 + \tau^2 n)} \\
0 \\
0
\end{bmatrix}
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \sigma^2 \partial \theta} \right]
&= \begin{bmatrix}
-\frac{1}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n \mathbb{E}\left[ (\mathbf{y}_{i,j} - \alpha)\right] \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \sum_{i = 1}^k \mathbb{E}\left[ \left( \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right)^2 \right] + \frac{2}{(\sigma^2)^3} \sum_{i = 1}^n \sum_{j = 1}^n \mathbb{E}\left[ (\mathbf{y}_{i,j} - \alpha)^2 \right] \right) \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \mathbb{E}\left[  \left( \sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right)^2 \right]
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \sum_{i = 1}^k \sum_{j = 1}^n \sum_{j'  = 1}^n\text{Cov}\left(\mathbf{y}_{i,j}, \mathbf{y}_{i,j'} \right)+ \frac{2}{(\sigma^2)^3} \sum_{i = 1}^n \sum_{j = 1}^n \text{Var}(\mathbf{y}_{i,j}) \right) \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n \sum_{j' = 1}^n \text{Cov}\left(\mathbf{y}_{i,j}, \mathbf{y}_{i,j'} \right)
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \sum_{i = 1}^k  \sum_{j = 1}^n \left[\text{Var}(\mathbf{y}_{i,j}) + \sum_{j' \neq j} \text{Cov}\left(\mathbf{y}_{i,j}, \mathbf{y}_{i,j'} \right) \right] + \frac{2}{(\sigma^2)^3} \sum_{i = 1}^n \sum_{j = 1}^n (\sigma^2 + \tau^2)\right) \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n\left[ \text{Var}(\mathbf{y}_{i,j}) + \sum_{j' \neq j}\text{Cov}\left(\mathbf{y}_{i,j}, \mathbf{y}_{i,j'} \right) \right]
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \sum_{i = 1}^k  \sum_{j = 1}^n \left[\sigma^2 + \tau^2 + \sum_{j' \neq j} \tau^2 \right] + \frac{2 n k (\sigma^2 + \tau^2)}{(\sigma^2)^3} \right) \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n\left[\sigma^2 + \tau^2 + \sum_{j' \neq j} \tau^2 \right] 
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right] \left[nk(\sigma^2 + \tau^2) + kn(n-1)\tau^2\right] + \frac{2 n k (\sigma^2 + \tau^2)}{(\sigma^2)^3} \right) \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \left(nk(\sigma^2 + \tau^2) + nk(n-1) \tau^2\right)
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2nk(\sigma^2 + \tau^2 n)}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right]+ \frac{2 n k (\sigma^2 + \tau^2)}{(\sigma^2)^3} \right) \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{nk(\sigma^2 + \tau^2 n)}{(\sigma^2 + \tau^2 n)^3} 
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
- \frac{1}{2}\left(\left[\frac{k - nk}{(\sigma^2)^2} - \frac{k}{(\sigma^2 + \tau^2 n)^2} \right] + \frac{2nk(\sigma^2 + \tau^2 n)}{n} \left[ \frac{1}{(\sigma^2 + \tau^2 n)^3} - \frac{1}{(\sigma^2)^3} \right]+ \frac{2 n k (\sigma^2 + \tau^2)}{(\sigma^2)^3} \right) \\
-\frac{nk}{2(\sigma^2 + \tau^2 n)^2}
\end{bmatrix}
\end{aligned}
\nonumber
$$
$$
\begin{aligned}
\mathbb{E}\left[ \frac{\partial \ell(\theta; \mathbf{y})}{\partial \tau^2 \partial \theta} \right] 
&= \begin{bmatrix}
- \frac{n}{(\sigma^2 + \tau^2 n)^2} \sum_{i = 1}^k \sum_{j = 1}^n \mathbb{E}\left[ (\mathbf{y}_{i,j} - \alpha) \right] \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \mathbb{E}\left[ \left(\sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right)^2\right] \\
\frac{n^2k}{2(\sigma^2 + \tau^2 n)^2} - \frac{n}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \mathbb{E}\left[ \left(\sum_{j = 1}^n (\mathbf{y}_{i,j} - \alpha) \right)^2\right] \\
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n \left[ \text{Var}(\mathbf{y}_{i,j}) + \sum_{j' \neq j} \text{Cov}(\mathbf{y}_{i,j}, \mathbf{y}_{i,j'}) \right] \\
\frac{n^2k}{2(\sigma^2 + \tau^2 n)^2} - \frac{n}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n \left[ \text{Var}(\mathbf{y}_{i,j}) + \sum_{j' \neq j} \text{Cov}(\mathbf{y}_{i,j}, \mathbf{y}_{i,j'}) \right]
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{1}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n \left[\sigma^2 + \tau^2 + \sum_{j' \neq j} \tau^2 \right] \\
\frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2} - \frac{n}{(\sigma^2 + \tau^2 n)^3} \sum_{i = 1}^k \sum_{j = 1}^n \left[\sigma^2 + \tau^2 + \sum_{j' \neq j} \tau^2 \right]
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{nk(\sigma^2 + \tau^2 + (n-1)\tau^2)}{(\sigma^2 + \tau^2 n)^3}  \\
\frac{n^2k}{2(\sigma^2 + \tau^2 n)^2} - \frac{n(nk(\sigma^2 + \tau^2 + (n-1)\tau^2))}{(\sigma^2 + \tau^2 n)^3}
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
\frac{nk}{2(\sigma^2 + \tau^2 n)^2} - \frac{nk(\sigma^2 + \tau^2n)}{(\sigma^2 + \tau^2 n)^3} \\
\frac{n^2 k}{2(\sigma^2 + \tau^2 n)^3} - \frac{n^2 k(\sigma^2 + \tau^2 n)}{(\sigma^2 + \tau^2 n)^3}
\end{bmatrix} \\
&= \begin{bmatrix}
0 \\
-\frac{nk}{2(\sigma^2 + \tau^2 n)^2} \\
- \frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2}
\end{bmatrix}
\end{aligned}
\nonumber
$$
</details>
Multiplying by $-1$ and evaluating at $\tilde{\theta}_0 = (\hat{\alpha}, \hat{\sigma}^2, 0)$ gives us the information matrix:
$$
\tilde{\mathcal{I}}(\tilde{\theta}_0) =
\begin{bmatrix}
\frac{nk}{\hat{\sigma}^2} & 0 & 0 \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2} & \frac{nk}{2(\hat{\sigma}^2)^2} \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2} & \frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2}-
\end{bmatrix}
$$
<details>
<summary>Derivations.</summary>
$$
\begin{aligned}
\tilde{\mathcal{I}}(\tilde{\theta}_0) &= 
\begin{bmatrix}
\frac{nk}{\hat{\sigma}^2} & 0 & 0 \\
0 & \frac{1}{2}\left(\frac{k - nk - k}{(\hat{\sigma}^2)^2} + \frac{2nk}{(\hat{\sigma}^2)^2} \right) & \frac{nk}{2(\hat{\sigma}^2)^2} \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2} & \frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2}
\end{bmatrix} \\
&= \begin{bmatrix}
\frac{nk}{\hat{\sigma}^2} & 0 & 0 \\
0 & \frac{- nk}{2(\hat{\sigma}^2)^2} + \frac{nk}{(\hat{\sigma}^2)^2} & \frac{nk}{2(\hat{\sigma}^2)^2} \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2} & \frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2}
\end{bmatrix} \\
&= \begin{bmatrix}
\frac{nk}{\hat{\sigma}^2} & 0 & 0 \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2} & \frac{nk}{2(\hat{\sigma}^2)^4} \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2} & \frac{n^2 k}{2(\sigma^2 + \tau^2 n)^2}
\end{bmatrix}
\end{aligned}
$$
</details>
We'll partition $\tilde{\mathcal{I}}(\tilde{\theta}_0)$ into the nuisance parameter ($\lambda = (\alpha, \sigma^2)$) and the parameter of interest, $\tau^2$. 
$$
\tilde{\mathcal{I}}(\tilde{\theta}_0) =
\begin{bmatrix}
  \tilde{\mathcal{I}}_{\lambda, \lambda}(\tilde{\theta}_0) & \tilde{\mathcal{I}}_{\lambda, \tau^2}(\tilde{\theta}_0) \\
  \tilde{\mathcal{I}}_{\tau^2, \lambda}(\tilde{\theta}_0) & \tilde{\mathcal{I}}_{\tau^2, \tau^2}(\tilde{\theta}_0)
\end{bmatrix}
$$
where we let:
$$
\begin{aligned}
\tilde{\mathcal{I}}_{\lambda, \lambda}(\tilde{\theta}_0) &= 
\begin{bmatrix}
  \tilde{\mathcal{I}}_{\alpha, \alpha}(\tilde{\theta}_0) & \tilde{\mathcal{I}}_{\alpha, \sigma^2}(\tilde{\theta}_0) \\
  \tilde{\mathcal{I}}_{\sigma^2, \alpha}(\tilde{\theta}_0) & \tilde{\mathcal{I}}_{\sigma^2, \sigma^2}(\tilde{\theta}_0)
\end{bmatrix} =
\begin{bmatrix}
\frac{nk}{\hat{\sigma}^2} & 0 \\
0 & \frac{nk}{2(\hat{\sigma}^2)^2}
\end{bmatrix} \\
\tilde{\mathcal{I}}_{\tau^2, \lambda}(\tilde{\theta}_0) &= 
\begin{bmatrix}
  \tilde{\mathcal{I}}_{\alpha, \tau^2}(\tilde{\theta}_0) & \tilde{\mathcal{I}}_{\sigma^2, \tau^2}(\tilde{\theta}_0)
\end{bmatrix} =
\begin{bmatrix}
  0 & \frac{nk}{2(\hat{\sigma}^2)^2}
\end{bmatrix} \\
\tilde{\mathcal{I}}_{\tau^2, \tau^2}(\tilde{\theta}_0) &= 
\begin{bmatrix}
\frac{n^2 k}{2(\hat{\sigma}^2)^2}
\end{bmatrix} 
\end{aligned}
$$
We only really need the block of the inverse Fisher information that corresponds to $(\tau^2, \tau^2)$. Using the <a href="https://chrisyeh96.github.io/2021/05/19/schur-complement.html#schur-complements">matrix inverse/Schur complement theorem</a>, we can find this pretty easily:
$$
\begin{aligned}
\tilde{\mathcal{I}}^{-1}_{\tau^2, \tau^2} (\tilde{\theta}_0)
&= \left[\tilde{\mathcal{I}}_{\tau^2, \tau^2}(\tilde{\theta}_0) - \tilde{\mathcal{I}}_{\tau^2, \lambda}(\tilde{\theta}_0) \tilde{\mathcal{I}}^{-1}_{\lambda, \lambda}(\tilde{\theta}_0) \tilde{\mathcal{I}}_{\lambda, \tau^2}(\tilde{\theta}_0)\right]^{-1} \\
&= \left[\frac{n^2 k}{2(\hat{\sigma}^2)^2} - 
\begin{bmatrix}
0 & \frac{n^2 k^2}{4(\hat{\sigma}^2)^2}
\end{bmatrix}
\begin{bmatrix}
0 \\
\frac{nk}{2(\hat{\sigma}^2)^2}
\end{bmatrix}\right]^{-1} \\
&= \left[\frac{n^2 k}{2(\hat{\sigma}^2)^2} - \frac{n^3 k^3}{8 (\hat{\sigma}^2)^6}\right]^{-1}
\end{aligned}
$$
Following Eq. \eqref{eq:test-stat-1}, we get our one-sided test statistic:
$$
T = \mathbf{S}^\top(\tilde{\theta}_0) \tilde{\mathcal{I}}_{\tau^2, \tau^2}^{-1}(\tilde{\theta}_0) \mathbf{S}(\tilde{\theta}_0) - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ \left(\mathbf{S}(\tilde{\theta}_0) - \mathbf{b}\right)^\top \tilde{\mathcal{I}}_{\tau^2, \tau^2}^{-1}(\tilde{\theta}_0) \left( \mathbf{S}(\tilde{\theta}_0)- \mathbf{b}\right) \right\}
$$
Since we are comparing a model with one random effect to one with zero, the large sample $p$-value is a fifty-fifty mixture of a $\chi^2_0$ and a $\chi^2_1$ distribution.<span markdown="1">[^fn-shapiro]</span>
</details>

---

### Results

<div class="theorem">
  <body>
  <strong>Lemma 1 (Silvapulle and Silvapulle<span markdown="1">[^fn-silvapulle]</span>).</strong>
  <br>
  Let $\hat{\theta}$ be an estimator of $\theta$ using the entire parameter space (no restrictions imposed). Let $\mathcal{P}$ denote a closed and convex cone with its vertex at the origin. Let $\mathbf{B}$ be a positive definite matrix independent of $\theta$, and let $\mathbf{B}_{\psi, \psi \cdot \lambda} = \mathbf{B}_{\psi, \psi} - \mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi}$. 
  <br>
  Note that:

  $$
  \begin{aligned}
  (\hat{\theta} - \theta)^\top \mathbf{B} (\hat{\theta} - \theta) &= 
  \left( \begin{bmatrix} \hat{\psi} \\ \hat{\lambda} \end{bmatrix} - \begin{bmatrix} \psi \\ \lambda \end{bmatrix} \right)^\top  
  \begin{bmatrix} \mathbf{B}_{\psi, \psi} & \mathbf{B}_{\psi, \lambda} \\ \mathbf{B}_{\lambda, \psi} & \mathbf{B}_{\lambda, \lambda} \end{bmatrix}
  \left(  \begin{bmatrix} \hat{\psi} \\ \hat{\lambda} \end{bmatrix} - \begin{bmatrix} \psi \\ \lambda \end{bmatrix} \right)
  \end{aligned}
  \nonumber
  $$

  The minimum of the above expression over just $\psi \in \mathcal{P}$ is equivalent to $\underset{\psi \in \mathcal{P}}{\min} \left\{ (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \psi \cdot \lambda}(\hat{\psi} - \psi) \right\}$ where we get $\mathbf{B}_{\psi, \psi \cdot \lambda}$ by adjusting $\mathbf{B}_{\psi, \psi}$ for the uncertainty in $\hat{\lambda}$. 
  
  Let $(\bar{\psi} : \bar{\lambda}) := \underset{\psi \in \mathcal{P}}{\arg \min} \left\{ (\hat{\theta} - \theta)^\top \mathbf{B}(\hat{\theta} - \theta)\right\}$. Then $\bar{\lambda} = \hat{\lambda} + \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi}(\hat{\psi} - \bar{\psi})$.

  <details>
  <summary>Proof.</summary>
  There are two steps that go into the proof. We need to show that:

  $$
  \underset{\psi \in \mathcal{P}}{\min} \left\{ (\hat{\theta} - \theta)^\top \mathbf{B} (\hat{\theta} - \theta) \right\} = \underset{\psi \in \mathcal{P}}{\min} \left\{ (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \psi \cdot \lambda}(\hat{\psi} - \psi) \right\}
  \nonumber
  $$

  We also need to show that $\bar{\lambda} = \hat{\lambda} + \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi}(\hat{\psi} - \bar{\psi})$. 
  <br>
  The second claim can be shown by supposing we know $\bar{\psi}$ and finding $\bar{\lambda}$ by optimizing over $\lambda$. Plugging in our expression for $\bar{\lambda}$ into the objective function, we can decompose it into two pieces, one of which is the objective function in the RHS of the first claim. Then all we need to do is show that the minimizer of that piece is the same as $\bar{\psi}$.
  <br>
  <hr>
  <br>
  Since the minimization is over $\psi \in \mathcal{P}$, we'll assume $\psi \in \mathcal{P}$. To save space, define $J(\theta) = (\hat{\theta} - \theta)^\top \mathbf{B}(\hat{\theta} - \theta)$. 
  <br>
  Notice that $\underset{\theta}{\min} \left\{ J(\theta) \right\}$ is the same as $\underset{\psi}{\min} \left\{ J(\psi : \bar{\lambda}) \right\}$ and also $\underset{\lambda}{\min} \left\{ J(\bar{\psi} : \lambda) \right\}$. That is, minimizing $J(\theta)$ over all values of $\theta$ is the same as if we plugged in the minimizing value of $\psi$ ($\bar{\psi}$) and then just minimized over values of $\lambda$ (and vice versa).
  <br>
  Since $J(\bar{\psi} : \lambda)$ has a quadratic form, we can set the derivative equal to $0$ and solve to get a minimizing value for $\lambda$:

  $$
  \begin{aligned}
  \frac{\partial}{\partial \lambda} \left[ J(\bar{\psi} : \lambda) \right]  &=  
  \frac{\partial}{\partial \lambda} \left[ 
    \left( \begin{bmatrix} \hat{\psi} \\ \hat{\lambda} \end{bmatrix} - \begin{bmatrix} \psi \\ \lambda \end{bmatrix} \right)^\top  
  \begin{bmatrix} \mathbf{B}_{\psi, \psi} & \mathbf{B}_{\psi, \lambda} \\ \mathbf{B}_{\lambda, \psi} & \mathbf{B}_{\lambda, \lambda} \end{bmatrix}
  \left(  \begin{bmatrix} \hat{\psi} \\ \hat{\lambda} \end{bmatrix} - \begin{bmatrix} \psi \\ \lambda \end{bmatrix} \right) \right]  \\
  &= \frac{\partial}{\partial \lambda} \left[ 
    \begin{bmatrix}
    (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \psi} + (\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \psi} \\
    (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} + (\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \lambda} 
    \end{bmatrix}^\top
    \left(  \begin{bmatrix} \hat{\psi} \\ \hat{\lambda} \end{bmatrix} - \begin{bmatrix} \psi \\ \lambda \end{bmatrix} \right) 
  \right] \\
  &= \frac{\partial}{\partial \lambda} \left[
    (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \psi}(\hat{\psi} - \psi) + (\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \psi) + 
    (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda}(\hat{\lambda} - \lambda) + (\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \lambda} (\hat{\lambda} - \lambda)
  \right] \\
  &= 2(\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \lambda} - (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} - \left( \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \psi) \right)^\top\\
  &= -2(\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \lambda} - 2(\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda}
  \end{aligned}
  \nonumber
  $$

  where in the last line we assume the authors mean <i>symmetric</i> positive definite when they say positive definite. 
  <br>

  Setting this equal to $0$ yields:

  $$
  \begin{aligned}
  &0 = \frac{\partial}{\partial \lambda} \left[ J(\bar{\psi} : \lambda) \right] \\
  \implies &0 = -2(\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \lambda} - 2(\hat{\psi} - \bar{\psi})^\top \mathbf{B}_{\psi, \lambda} \\
  \implies &-(\hat{\psi} - \bar{\psi})^\top\mathbf{B}_{\psi, \lambda} = (\hat{\lambda} - \lambda)^\top \mathbf{B}_{\lambda, \lambda} \\
  \implies &\lambda^\top \mathbf{B}_{\lambda, \lambda} = (\hat{\psi} - \bar{\psi})^\top \mathbf{B}_{\psi, \lambda} + \hat{\lambda}^\top \mathbf{B}_{\lambda, \lambda} \\
  \implies &\lambda^\top = (\hat{\psi} - \bar{\psi})^\top\mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} + \hat{\lambda}^\top \\
  \implies &\bar{\lambda} = \hat{\lambda} + \mathbf{B}_{\lambda, \lambda}^{-1}\mathbf{B}_{\lambda, \psi} (\hat{\psi} - \bar{\psi})
  \end{aligned}
  \nonumber
  $$

  This shows that that there exists a point $\bar{\lambda}$ satisfying the form specified in the Lemma. Now we want to show that this corresponds to the minimum over $\psi \in \mathcal{P}$.
  <br>
  Since this value of $\lambda$ is optimal for any fixed $\psi$, we can just plug this into $J(\psi : \bar{\lambda})$ and minimize to find a value for $\bar{\psi}$:

  $$
  \begin{aligned}   
  J(\psi : \bar{\lambda})  &=  
   \begin{bmatrix} \hat{\psi} - \psi \\ -\mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \bar{\psi}) \end{bmatrix}^\top  
  \begin{bmatrix} \mathbf{B}_{\psi, \psi} & \mathbf{B}_{\psi, \lambda} \\ \mathbf{B}_{\lambda, \psi} & \mathbf{B}_{\lambda, \lambda} \end{bmatrix}
  \begin{bmatrix}  \hat{\psi} - \psi \\ -\mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \bar{\psi}) \end{bmatrix} \\
  &= \begin{bmatrix}
  (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \psi} - (\hat{\psi} - \bar{\psi})^\top  \mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi} \\
  (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} - (\hat{\psi} - \bar{\psi})^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \lambda}
  \end{bmatrix}^\top
  \begin{bmatrix}  \hat{\psi} - \psi \\ -\mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \bar{\psi}) \end{bmatrix} \\
  &= 
  (\hat{\psi} - \psi)^\top\mathbf{B}_{\psi, \psi}(\hat{\psi} - \psi) \underbrace{
  - (\hat{\psi} - \bar{\psi})^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi}(\hat{\psi} - \psi) 
  - (\hat{\psi} - \psi)^\top\mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \bar{\psi}) 
  + (\hat{\psi} - \bar{\psi})^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}_{\lambda, \lambda}^{-1} \mathbf{B}_{\lambda, \psi} (\hat{\psi} - \bar{\psi})}_{(a)}\\
  &= (\hat{\psi} - \psi)^\top\mathbf{B}_{\psi, \psi}(\hat{\psi} - \psi) - (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}^{-1}_{\lambda, \lambda} \mathbf{B}_{\lambda, \psi}(\hat{\psi} - \psi) +  (\bar{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}^{-1}_{\lambda, \lambda} \mathbf{B}_{\lambda, \psi}(\bar{\psi} - \psi) \\
  &= \underbrace{(\hat{\psi} - \psi)^\top\left( \mathbf{B}_{\psi, \psi} - \mathbf{B}_{\psi, \lambda} \mathbf{B}^{-1}_{\lambda, \lambda} \mathbf{B}_{\lambda, \psi} \right) (\hat{\psi} - \psi)}_{=: f(\psi)} + \underbrace{(\bar{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}^{-1}_{\lambda, \lambda} \mathbf{B}_{\lambda, \psi}(\bar{\psi} - \psi)}_{=: g(\psi)}
  \end{aligned}
  \nonumber
  $$    

  <details>
  <summary>Details Of $(a)$.</summary>
  In $(a)$, we do the following simplification. Let $\bar{\mathbf{B}} = \mathbf{B}_{\psi, \lambda} \mathbf{B}^{-1}_{\lambda, \lambda} \mathbf{B}_{\lambda, \psi}$. Then:

  $$
  \begin{aligned}
  (a) &= -(\hat{\psi} - \bar{\psi})^\top \bar{\mathbf{B}} (\hat{\psi} - \psi) - (\hat{\psi} - \psi)^\top \bar{\mathbf{B}}  (\hat{\psi} - \bar{\psi}) + (\hat{\psi} - \bar{\psi})^\top \bar{\mathbf{B}} (\hat{\psi} - \bar{\psi}) \\
  &= \left(-\hat{\psi}^\top\bar{\mathbf{B}}\hat{\psi}+ \bar{\psi}^\top \bar{\mathbf{B}} \hat{\psi} + \hat{\psi}^\top \bar{\mathbf{B}} \psi - \bar{\psi}^\top \bar{\mathbf{B}} \psi \right)
  + \left( - \hat{\psi}^\top \bar{\mathbf{B}} \hat{\psi} + \psi^\top \bar{\mathbf{B}} \hat{\psi} + \hat{\psi}^\top \bar{\mathbf{B}} \bar{\psi} - \psi^\top \bar{\mathbf{B}} \bar{\psi} \right) 
  + \left( \hat{\psi}^\top \bar{\mathbf{B}} \hat{\psi} - \bar{\psi}^\top \bar{\mathbf{B}} \hat{\psi} - \hat{\psi}^\top \bar{\mathbf{B}} \bar{\psi} + \bar{\psi}^\top \bar{\mathbf{B}} \bar{\psi} \right) \\
  &= \underbrace{-\hat{\psi}^\top \bar{\mathbf{B}}\hat{\psi} + \hat{\psi}^\top \bar{\mathbf{B}} \psi + \psi^\top \bar{\mathbf{B}} \hat{\psi}}_{(i)} + \underbrace{\bar{\psi}^\top \bar{\mathbf{B}} \bar{\psi} -\bar{\psi}^\top \bar{\mathbf{B}} \psi - \psi^\top \bar{\mathbf{B}} \bar{\psi}}_{(ii)} \\
  &= \underbrace{(\bar{\psi} - \psi)^\top \bar{\mathbf{B}} (\bar{\psi} - \psi) - \psi^\top \bar{\mathbf{B}} \psi}_{=(ii)} - \underbrace{\left((\hat{\psi} - \psi)^\top \bar{\mathbf{B}}(\hat{\psi} - \psi) + \psi^\top \bar{\mathbf{B}}\psi\right)}_{=(i)} \\
  &= (\bar{\psi} - \psi)^\top \bar{\mathbf{B}} (\bar{\psi} - \psi) - (\hat{\psi} - \psi)^\top \bar{\mathbf{B}}(\hat{\psi} - \psi)
  \end{aligned}
  \nonumber
  $$
  </details>

  We can see that minimizing $J(\psi: \bar{\lambda})$ over values of $\psi$ is equivalent to minimizing $f(\psi) + g(\psi)$, where:

  $$
  f(\psi) = (\hat{\psi} - \psi)^\top \mathbf{B}_{\psi, \psi \cdot \lambda} (\hat{\psi} - \psi)\hspace{25mm} g(\psi) = (\bar{\psi} - \psi)^\top \mathbf{B}_{\psi, \lambda} \mathbf{B}^{-1}_{\lambda, \lambda} \mathbf{B}_{\lambda, \psi}(\bar{\psi} - \psi)
  \nonumber
  $$

  We also know that $f(\psi)$ and $g(\psi)$ are strictly convex and convex, respectively.

  <details>
  <summary>Proofs Of Convexity.</summary>
  First, let's look at $f(\psi)$, in which the middle portion is the Schur complement of $\mathbf{B}_{\lambda, \lambda}$. The positive definiteness of $\mathbf{B}$ implies that the Schur complement is also positive definite (since we also assume $\mathbf{B}_{\lambda, \lambda}$ is invertible). Since $f(\psi)$ has a quadratic form associated with a positive definite matrix, it is strictly convex.
  <br>
  <br>
  Looking at $g(\psi)$, we see that the middle portion, $\mathbf{B}_{\lambda, \lambda}^{-1}$ is positive definite due to the fact that the principal sub-matrices of a positive definite matrix are also positive definite and the inverse of a positive definite matrix is also positive definite. Secondly, we know that pre- and post-multiplying a positive definite matrix by another matrix will yield a positive semi-definite matrix. Thus, $g(\psi)$ is positive semi-definite, which implies that it is convex since it is a quadratic form associated with a positive semi-definite matrix.
  </details>
  
  Denote the minimizer of $f(\psi)$ with $\psi^*$ and consider the line segment $[ \bar{\psi}, \psi^* ]$, which connects the minimizing value of $\psi$ for $J(\theta)$ with the minimizing value of $\psi$ for $f(\psi)$. 
  <br>
  Supposedly we can show that $\psi^* = \bar{\psi}$, which would complete the proof. However, I am unsure of how to show this. Here is my best attempt so far, which is an adaptation of AlephZero's proof on Math StackExchange:
  <br>
  <i> Assume $(\psi^*, \bar{\lambda})$ is a solution to the minimization problem (I don't know how to show this). </i> Suppose $\psi^* \neq \bar{\psi}$. Since both $\bar{\psi}$ and $\psi^*$ are within $\mathcal{P}$, which is closed and convex, the entire line segment is also contained in $\mathcal{P}$ and any convex combination of $\bar{\psi}$ and $\psi^*$ is also a solution. Let $a := J(\bar{\psi} : \bar{\lambda}) = J(\psi^* : \bar{\lambda})$, and let $t \in (0, 1)$. Let $t(\psi : \lambda)$ denote the scaling of the respecting entries in the concatenated vector of parameter values. We then have:

  $$
  \begin{aligned}
  J(t(\bar{\psi} : \bar{\lambda}) + (1 - t)(\psi^* : \bar{\lambda})) &= J(t\bar{\psi} + (1-t)\psi^* : \bar{\lambda}) \\
  &= f(t\bar{\psi} + (1-t)\psi^*) + g(t\bar{\psi} + (1-t) \psi^*) \\
  &< t f(\bar{\psi}) + (1-t)f(\psi^*) + tg(\bar{\psi}) + (1-t)g(\psi^*) \hspace{15mm} f \text{ is strictly convex and } g \text{ is convex} \\
  &= t J(\bar{\psi} : \bar{\lambda}) + (1 - t) J(\psi^* : \bar{\lambda}) \\
  &= a
  \end{aligned}
  \nonumber
  $$

  This implies any convex combination of $\bar{\psi}$ and $\psi^*$ achieves a value smaller than the minimum, which is a contradiction. Thus, $\bar{\psi} = \psi^*$. 

  <br>
  <br>
  <p style="color:red;">TODO: FINISH PROOF</p>
  </details>
  </body>
</div>

<br>
With the above lemma, we can prove the following theorem.

<div class="theorem">
  <body>
  <strong>Theorem 1 (Silvapulle and Silvapulle<span markdown="1">[^fn-silvapulle]</span>).</strong>
  <br>
  Define $\mathbf{S}_n(\theta) = \frac{\partial}{\partial \theta} \left[ \ell(\theta) \right]$ as the score function (the derivative of the log-likelihood), and assume that it satisfies Condition \eqref{eq:condition-a}. Suppose we are testing $H_0: \psi = \mathbf{0}$  against $H_A: \psi \in \mathcal{C}$ for $\mathcal{C}$ as defined above. As $n \rightarrow \infty$, the likelihood ratio test statistic, $LR = -2 \left(\ell(\theta_0) - \ell(\hat{\theta}) \right)$ where $\hat{\theta}$ is the MLE of $\theta$ over the entire parameter space, satisfies:

  $$
  LR = T_s + o_p(1)
  \nonumber
  $$

  under the null hypothesis.

  <details>
  <summary>Proof.</summary>
  We'll assume we're in the setting where the null hypothesis is true ($\theta_0 = (\lambda : \mathbf{0})$). Since we've taken our estimating equation, $\mathbf{S}_n(\theta)$, as the score, we can take $\mathbf{G}(\theta_0)$ to be the Fisher information under the null (WHY?).
  <br>
  First, let's rewrite the log-likelihood, $\ell(\theta)$, with a Taylor expansion about the unrestricted estimator, $\hat{\theta}$:

  $$
  \ell(\theta) = \ell(\hat{\theta}) + \frac{n}{2}(\hat{\theta} - \theta)^\top \mathbf{G}(\theta_0) (\hat{\theta} -\theta) + \Delta(\theta)
  \nonumber
  $$

  where, for any $a > 0$, $\underset{\rvert \rvert \theta - \theta_0 \rvert \rvert \leq \frac{a}{\sqrt{n}}}{\sup} \left\{ \rvert \Delta(\theta) \rvert \right\} = o_p(1)$. This last term is just a bound on all of the terms in the Taylor expansion that are of a higher order than two, written with stochastic order notation. A related proof can be found in my post "A Score Test Primer". The multiplication by $n$ is just because we assume we have $n$ i.i.d. observations, so we can just multiply the value for a single observation by the number we have.
  <br>
  From this expansion, we can rewrite the likelihood ratio test statistic as:

  $$
  LR = n\left[ \underset{\psi = \mathbf{0}}{\min} \left\{ (\hat{\theta} - \theta)^\top \mathbf{G}(\theta_0) (\hat{\theta} - \theta)\right\} - \underset{\psi \in \mathcal{C}}{\min}\left\{ (\hat{\theta} - \theta)^\top \mathbf{G}(\theta_0) (\hat{\theta} - \theta) \right\} \right] + o_p(1)
  \nonumber
  $$

  Using Lemma 1, we can see that:

  $$
  \begin{aligned}
  LR &= n\left[ \underset{\psi = \mathbf{0}}{\min} \left\{ (\hat{\psi} - \psi)^\top \mathbf{G}_{\psi, \psi \cdot \lambda}(\theta_0) (\hat{\psi} - \psi) \right\} - \underset{\psi \in \mathcal{C}}{\min}\left\{(\hat{\psi} - \psi)^\top \mathbf{G}_{\psi, \psi \cdot \lambda}(\theta_0) (\hat{\psi} - \psi) \right\} \right] + o_p(1) \\
  &= n \left[ \underset{\psi = \mathbf{0}}{\min} \left\{ \hat{\psi}^\top \mathbf{G}_{\psi, \psi \cdot \lambda}(\theta_0) \hat{\psi}\right\} + \underset{\psi \in \mathcal{C}}{\min}\left\{(\hat{\psi} - \psi)^\top \mathbf{G}_{\psi, \psi \cdot \lambda}(\theta_0) (\hat{\psi} - \psi) \right\}\right] + o_p(1) \\
  &= n \left[ \hat{\psi}^\top \mathbf{G}_{\psi, \psi \cdot \lambda}(\theta_0) \hat{\psi} + \underset{\psi \in \mathcal{C}}{\inf}\left\{(\hat{\psi} - \psi)^\top \mathbf{G}_{\psi, \psi \cdot \lambda}(\theta_0) (\hat{\psi} - \psi) \right\}\right] + o_p(1)  \hspace{15mm} \mathcal{C} \text{ is closed, so } \min \text{ is same as } \inf \\
  \end{aligned}
  \nonumber
  $$

  Notice that $\sqrt{n} \mathbf{G}(\theta_0) (\hat{\theta} - \theta_0) = \frac{1}{\sqrt{n}} \mathbf{S}_n(\theta_0) + o_p(1)$ because 
$$
T_s = \mathbf{U}^\top \tilde{\mathbf{A}}_{\psi, \psi}^{-1}\mathbf{U} - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{U} - \mathbf{b})^\top \tilde{\mathbf{A}}_{\psi, \psi}^{-1}(\mathbf{U} - \mathbf{b}) \right\}
$$
    <p style="color:red;">TODO: FINISH PROOF</p>
  </details>
  </body>
</div>


---
## References

[^fn-aleph]: AlephZero. “Does This Special Case of Convex Quadratic Programming Have a Partially-Unique Solution?” Mathematics Stack Exchange, February 12, 2019. https://math.stackexchange.com/questions/3108304/does-this-special-case-of-convex-quadratic-programming-have-a-partially-unique-s. 

[^fn-shapiro]: Shapiro, A. (1988). Towards a Unified Theory of Inequality Constrained Testing in Multivariate Analysis. International Statistical Review / Revue Internationale de Statistique, 56(1), 49. https://doi.org/10.2307/1403361

[^fn-silvapulle]: Silvapulle, M. J., & Silvapulle, P. (1995). A Score Test Against One-Sided Alternatives.

[^fn-coxhinkley]: Cox, D. R., & Hinkley, D. V. (2017). Theoretical statistics. CRC Press, Taylor & Francis Group.

[^fn-hall]: Hall, D. B., & Praestgaard, J. T. (2001). Order-restricted score tests for homogeneity in generalised linear and nonlinear mixed models. Biometrika, 88(3), 739–751. https://doi.org/10.1093/biomet/88.3.739

[^fn-kudo]: Kudo, A. (1963). A Multivariate Analogue of the One-Sided Test. Biometrika, 50, 403–418.

[^fn-lin]: Lin, X. (1997). Variance component testing in generalised linear models with random effects. Biometrika, 84(2), 309–326. https://doi.org/10.1093/biomet/84.2.309

[^fn-dasgupta]: DasGupta, A. (2008). Asymptotic Theory of Statistics and Probability. Springer New York. https://doi.org/10.1007/978-0-387-75971-5.
