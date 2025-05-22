---
layout: post
title:  "One-Sided Score Test"
date: 18 February 2025
categories: posts
tags: ["theory", "likelihood", "score test"]
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

For some reason, I have not been able to exactly identify what has been going on in my simulations to lead to the strange distributions. My gut tells me it is either a mistake in my derivations or something that I'm not realizing in theory (perhaps related to the fact that the null hypothesis parameter value is on the boundary of the parameter space?).

Instead of trying to derive the score test statistic again (for the umpteenth time...), I'm going to try a different angle. This post will cover comments and results related to one-sided score-type tests.

Note: Not all of the proofs are finished/included. I am hoping to find the time to return to thist post and complete them.


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
## Shapiro (1988)

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

Shapiro denotes the orthogonal projection of a point onto $C$ with $$P(\mathbf{x}, C) = \underset{\eta \in C}{\min} \left\{ (\mathbf{x} - \eta)^\top \mathbf{U} (\mathbf{x} - \eta) \right\}$$ where $\mathbf{U}$ is any positive-definite matrix. We'll use $\rvert \rvert \mathbf{x} \rvert \rvert = \sqrt{\mathbf{x}^\top \mathbf{U} \mathbf{x}}$ to denote the norm of $\mathbf{x}$ associated with $\mathbf{U}$, and we'll use $\langle\mathbf{x}, \mathbf{y}\rangle = \mathbf{x}^\top \mathbf{U}\mathbf{y}$ to denote the inner product of $\mathbf{x}$ and $\mathbf{y}$ associated with $\mathbf{U}$.

<div id="dual-cone"></div>
<div class="definition">
  <body>
  <strong>Definition (Dual Cone).</strong>
  <br>
 The <i>dual</i> or <i>polar cone</i> is the set $C^0= \{ \mathbf{y} \rvert \langle \mathbf{x}, \mathbf{y} \rangle \leq 0 \hspace{2mm} \forall \mathbf{x} \in C \}$. If $C$ is a vector subspace, then $C^0$ is the orthogonal complement of $C$, and if $C$ is closed and convex, then $(C^0)^0 = C$. 
  </body>
</div>

If we have another convex cone $K$, then:

$$
\rvert \rvert \mathbf{x} - P(\mathbf{x}, C) \rvert \rvert^2 = \rvert \rvert \mathbf{x} - P(\mathbf{x}, K) \rvert \rvert^2 + \rvert \rvert P(\mathbf{x}, K) - P(\mathbf{x}, C) \rvert \rvert^2
\nonumber
$$

if $C$ or $K$ is a vector space. Furthermore, if $C$ is a vector subspace, then $\mathbf{x} - P(\mathbf{x}, C) = P(\mathbf{x}, C^0)$. That is, the difference between $\mathbf{x}$ and the orthogonal projection of $\mathbf{x}$ onto $C$ is equivalent to its orthogonal projection onto the dual cone $C^0$.


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

where $\chi_i^2$ is a $\chi^2$ random variable with $i$ degrees of freedom, and $w_i$ are individuals wegiths that sum to $0$. As is standard, we let $\chi^2_0$ be a point mass at $0$. We'll denote this mixture distribution as $\mathcal{\bar{X}}^2(\mathbf{V}, C)$, since it depends on both $\mathbf{V}$ and the cone, $C$.

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
## Silvapulle And Silvapulle

I'll next dip into the work of Silvapulle and Silvapulle[^fn-silvapulle], who present a score-type test statistic for one-side alternative hypotheses based upon estimating functions instead of the true score function.

#### Set-Up
The authors present a fairly general setting where we do not assume we know the exact form of the distribution of the observations, only that they depend on some $k \times 1$-dimensional vector-valued parameter, $\theta$, that is partitioned into the nuisance parameters, $\lambda$, and the components of interest, $\psi$. We'll write the partitioned parameter vector as $(\lambda : \psi) = (\lambda^\top, \psi^\top)^\top$ where $\lambda$ is $(k - q) \times 1$ and $\psi$ is $q \times 1$.

We're interested in testing hypotheses of the form:

$$
H_0: \psi = 0 \hspace{10mm} H_A: \psi \in \mathcal{C}
\nonumber
$$

where $\mathcal{C}$ is some $q$-dimensional Euclidean space, or a closed, convex cone with its vertex at the origin. Define the $q \times 1$ vector, $\mathbf{U}_0$, as some function of the data satisfying, under the sequence of alternatives $K_n: \psi = n^{-1/2}\delta$, $U_0 \rightsquigarrow \mathcal{N}(\delta, \mathbf{D})$ as our sample size $n \rightarrow \infty$. In the previous, $\mathbf{D} = \mathbf{D}(\lambda)$ is some (fixed) matrix-valued function of the nuisance parameters.

Notice that testing the alternative hypothesis that $\psi \geq 0$ is equivalent to testing $\delta \geq 0$, which implies that the null hypothesis is equivalent to $\delta = 0$. We will keep this in mind in the followign steps.

The authors define an intermediate test statistic with:

$$
T = \mathbf{U}_0^\top \tilde{\mathbf{D}}^{-1}\mathbf{U}_0 - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{U}_0 - \mathbf{b})^\top \tilde{\mathbf{D}}^{-1}(\mathbf{U}_0 - \mathbf{b})\right\}
\label{eq:test-stat-1}
$$

where $\tilde{\mathbf{D}}$ is any consistent estimator of $\mathbf{D}$ under the null hypothesis.

A $p$-value for large sample sizes can be found by defining $\mathbf{Z} \sim \mathcal{N}(\mathbf{0}, \mathcal{D})$ and:

$$
\xi(t, \mathbf{D}, \mathcal{C}) = \mathbb{P}\left( \left[ \mathbf{Z}^\top \mathbf{D}^{-1} \mathbf{Z} - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{Z} - \mathbf{b})^\top \mathbf{D}^{-1}(\mathbf{Z} - \mathbf{b}) \right\} \right] \geq t \right)
\nonumber
$$

The quantity $1 - \xi(t, \mathbf{D}, \mathcal{C})$ follows a chi-bar-squared distribution; that is, a mixture of chi-squared distributions as we introduced in the previous section (Shapiro, 1988). The weights for the mixture can be hard to find, but we can get around this using the fact that, for large enough $n$ and under $H_0$ (i.e. $\delta = 0$), $\mathbf{U}_0$ is approximately $\mathcal{N}(\mathbf{0}, \mathbf{D}(\lambda))$. Thus, $\mathbb{P}(T \geq t; \lambda) \approx \xi(t, \mathbf{D}(\lambda), \mathcal{C})$. 

Suppose we observe a value of $T$, $t^*$. Define $$\mathbf{D}^*(\lambda)$$ as a consistent estimator of $\mathbf{D}(\lambda)$ for any $\lambda$. Then it follows that:

$$
p \approx \underset{\lambda}{\sup} \left\{ \xi(t^*, \mathbf{D}^*(\lambda), \mathcal{C}) \right\}
\nonumber
$$

for large enough $n$ because $\lambda$ is a nuisance parameter, so we can take the "best" probability over all of its values. 

#### Test Statistic
Define $\mathbf{S}_n(\theta) = \mathbf{0}$ as any $k \times 1$ vector estimating equation for $\theta$ (e.g. the score function or something else!) such that there exist non-singular $\mathbf{G}(\theta)$ and $\mathbf{V}(\theta)$ satisfying for any $a > 0$:

$$
\frac{1}{\sqrt{n}}\mathbf{S}_n(\theta) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathbf{V}(\theta)) \hspace{10mm} \text{and} \hspace{10mm} \underset{\rvert \rvert \mathbf{h} \rvert \rvert \leq a}{\sup} \left\{ \frac{1}{\sqrt{n}} \left( \mathbf{S}_n\left(\theta + \frac{1}{\sqrt{n}} \mathbf{h}\right) - \mathbf{S}_n(\theta) \right)  + \mathbf{G}(\theta) \mathbf{h} \right\} = o_p(1)
\label{eq:condition-a}
$$

(Recall $o_p$ is stochastic order notation for convergence in probability to $0$). Partition $\mathbf{S}_n$ into $$(S_{n, \lambda} : S_{n, \psi})$$ and $$\mathbf{G}$$ into $$\mathbf{G}_{\lambda, \lambda}$$, $$\mathbf{G}_{\lambda, \psi}$$, $$\mathbf{G}_{\psi, \lambda}$$, and $$\mathbf{G}_{\psi, \psi}$$. Let $$\theta_0 = (\lambda : \mathbf{0})$$ denote the value of $\theta$ under the null hypothesis.

Suppose the null hypothesis is true. Define $$\mathbf{Z}_n(\theta_0) = n^{-1/2} \left( \mathbf{S}_{n, \psi}(\theta_0) - \mathbf{G}_{\psi, \lambda}(\theta_0) \mathbf{G}_{\lambda, \lambda}^{-1}(\theta_0) \mathbf{S}_{n, \lambda}(\theta_0) \right)$$. Since the above condition (Eq. \eqref{eq:condition-a}) is assumed to be satisfied, we know that $$n^{-1/2} \mathbf{S}_n(\theta_0) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathbf{V}(\theta_0))$$. Since $\mathbf{Z}_n$ is just a function of $\mathbf{S}_n$, this implies that $$\mathbf{Z}_n(\theta_0) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathbf{C}(\theta_0))$$, where $$\mathbf{C} = \mathbf{V}_{\psi, \psi} - \mathbf{G}_{\psi, \lambda} \mathbf{G}^{-1}_{\lambda, \lambda} \mathbf{V}_{\lambda, \psi} - \left( \mathbf{V}_{\psi, \lambda} - \mathbf{G}_{\psi, \lambda} \mathbf{G}_{\lambda, \lambda}^{-1} \mathbf{V}_{\lambda, \lambda} \right)\left(\mathbf{G}^{-1}_{\lambda, \lambda}\right)^\top \mathbf{G}^\top_{\psi, \lambda}$$ (in which we have partitioned $\mathbf{V}(\theta)$ in the same way as we did with $\mathbf{G}(\theta)$).

Denote consistent estimators for $\mathbf{G}(\theta_0)$ and $\mathbf{V}(\theta_0)$ with $\tilde{\mathbf{G}}$ and $\tilde{\mathbf{V}}$, respectively. Furthermore, let $\tilde{\lambda}$ denote a "suitable" estimator for $\lambda$ (where suitable is not really specific, but examples are given in the text). Let $\tilde{\theta}_0 = (\tilde{\lambda} : \mathbf{0})$.

Define the following:

$$
\mathbf{G}^{\psi, \psi} = \left( \mathbf{G}_{\psi, \psi} - \mathbf{G}_{\psi, \lambda} \mathbf{G}^{-1}_{\lambda, \lambda} \mathbf{G}_{\lambda, \psi} \right)^{-1} \hspace{15mm} \mathbf{A}(\theta) = \left( \mathbf{G}^\top(\theta) \mathbf{V}^{-1}(\theta)\mathbf{G}(\theta) \right)^{-1} \hspace{15mm} \mathbf{U} = \tilde{\mathbf{G}}^{\psi, \psi} \tilde{\mathbf{Z}}_n
$$

where $\tilde{\mathbf{G}}^{\psi, \psi}$ is a consistent estimator for $\mathbf{G}^{\psi, \psi}$, and $\tilde{\mathbf{Z}}_n$ is $\mathbf{Z}_n$ found using $\tilde{\theta}_0$, $\tilde{\mathbf{G}}$, and $\tilde{\mathbf{V}}$.

Let's partition $\mathbf{A}(\theta)$ in the same way that we did with $\mathbf{V}(\theta)$ and $\mathbf{G}(\theta)$. With some work, we can see that for a fixed $\delta \in \mathcal{C}$, $\mathbf{U} \rightsquigarrow \mathcal{N}(\delta, \mathbf{A}_{\psi, \psi}(\theta_0))$ under the sequence of alternatives $H_n: \psi = n^{-1/2} \delta$ as we take $n \rightarrow \infty$. 

The test statistic for one-sided alternatives is then:

$$
T_s = \mathbf{U}^\top \tilde{\mathbf{A}}_{\psi, \psi}^{-1} \mathbf{U} - \underset{\mathbf{b} \in \mathcal{C}}{\inf} \left\{ (\mathbf{U} - \mathbf{b})^\top \tilde{\mathbf{A}}_{\psi, \psi}^{-1} (\mathbf{U} - \mathbf{b}) \right\}
\label{eq:test-stat-2}
$$

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