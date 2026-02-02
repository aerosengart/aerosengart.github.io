---
layout: distill
title: Estimation Theory
description: A Primer
date: 2025-06-17
tabs: true
tags: theory estimation primer
toc:
  - name: Set-Up
    subsections:
        - name: Notation
  - name: Hilbert Spaces
  - name: Extremum Estimators
    subsections:
        - name: Influence Functions
  - name: Consistency
    subsections:
        - name: Showing Consistency
  - name: Asymptotic Normality
    subsections:
        - name: Estimating Asymptotic Variance
  - name: Efficiency
    subsections:
        - name: Efficient Influence Functions
        - name: Efficient Score
  - name: Appendix
    subsections:
        - name: Auxiliary Definitions
        - name: Auxiliary Results
bibliography: 2025-06-17-estimation-theory.bib
---


Recently, I've run into some confusion in parameter estimation. In this post, I'll cover the relevant sections in Tsiatis' <i>Semiparametric Theory and Missing Data</i><d-cite key=tsiatis2006></d-cite> and Newey and McFadden's <i>Large Sample Estimation and Hypothesis Testing</i><d-cite key=newey1994></d-cite> in my journey to answering some of my questions:

<ol>
<li>What is the <i>efficient information matrix</i>?</li>
<li>When exactly are the quasi-likelihood estimators asymptotically normal and consistent?</li>
<li>Are any regularity conditions violated in the variance component estimation scheme for generalized linear mixed models?</li>
</ol>

---

## Set-Up
We assume we have independent and identically distributed random variables, $X_1, \dots, X_n$, arranged into a vector $\mathbf{X}$ defined on the probability space $(\mathcal{X}, \mathcal{A}, P)$. We assume $X_i \sim P$ where $P \in \mathcal{P}$ for some family of distributions, $\mathcal{P}$. We also assume that $P$ is defined by the values of parameters, $\theta$, where some coordinates are of interest and some are a nuisance (i.e. we can partition the parameter space). We use $\beta$ and $\eta$ to denote the parameters of interest and the nuisance parameters, respectively. We have:

$$
\begin{equation}
\label{eq:parameter-vector}
\theta = (\beta^\top, \eta^\top)^\top \in \mathbb{R}^{p \times 1} \hspace{2mm} (p = q + r),
\hspace{8mm}
\beta \in \mathbb{R}^{q \times 1},
\hspace{8mm}
\eta \in \mathbb{R}^{r \times 1}
\end{equation}
$$

As in Tsiatis' book, we'll restrict the parameter space, $\Omega$, to subspaces of linear vector spaces, which themselves could be finite- or infinite-dimensional (depending on whether there exist a finite number of elements in it that span it). In most applications/examples, this assumption will not be too restrictive since we will usually be working in Euclidean spaces. 

### Notation
Uppercase letters denote random variables, and lowercase letters denote realizations of these variables. Boldface indicates a vector or matrix, exactly which will be clear from context. 

We will indicate parameters that are functions with parenthesis (e.g. $\gamma(\cdot)$), and we'll use $\hat{\gamma}_n$ to denote an estimator of parameter $\gamma$. We'll use a subscript of $0$ to denote the true parameter values (i.e. $\theta_0 = (\beta_0^\top, \eta_0^\top)^\top$). 

$\nu_X(x)$ denotes the <i>dominating measure</i> on which densities for $\mathbf{X}$ are defined. That is, the dominating measure is any measure $\mu(\cdot)$ such that the density $f_\mathbf{X}(\mathbf{x})$ exists and such that $\mu'(A) = 0$ implies $\mu(A) = 0$. For continuous $X$, we usually use the Lebesgue measure, and for discrete $X$, we use the counting measure.

I am pretty sure all of the expectations are taking with respect to the true data-generating distribution, $P_X(x, \theta_0)$. I've denoted this with a subscript of $\theta_0$ on the expectation, but I may have missed a few. 

Convergence in probability and convergence in distribution when the density of the random variable $X$ is equal to $p_X(x; \beta, \eta)$ are denoted by, respectively:

$$
\xrightarrow{P\{ \beta, \eta \}},
\hspace{10mm}
\xrightarrow{D\{ \beta, \eta \}}
$$

General convergence in probability and distribution are denoted by:

$$
\overset{p}{\rightarrow},
\hspace{10mm}
\rightsquigarrow
$$

---

## Hilbert Spaces
Before we dive into our discussion of semiparametric inference, we need to first look at some geometric definitions. 

Consider the space of functions of the form $h(X)$ where $h: \mathcal{X} \rightarrow \mathbb{R}^q$, are measurable (w.r.t. the probability space we described earlier), and are such that:

$$
\mathbb{E}[h(X)] = 0, 
\hspace{10mm}
\mathbb{E}\left[ h^\top(X) h(X)\right] < \infty
$$

That is, we consider the space of all measurable, $q$-dimensional random functions with mean zero and finite variance. Notice that this space is a linear space. Let $\mathbf{0}$ (the constant function outputting a $q$-dimensional vector of zeros) denote the origin. Though these are random functions, we will sometimes drop the $(X)$ from our notation to just write $h$. 

Later we will concern ourselves with Hilbert spaces of these types of functions. We define these formally for completeness. Recall that a Hilbert space is a complete normed linear vector space with an inner product. 

<div id="def-hilbert-space"></div>
<div class="definition">
  <strong>Definition (Hilbert Space of Mean-Zero Random Functions).</strong>
  <br>
    Consider the space of $q$-dimensional random functions $h: \mathcal{X} \rightarrow \mathbb{R}^q$ that are measurable (w.r.t. probability space $(\mathcal{X}, \mathcal{A}, P)$) and satisfy:
    $$
    \mathbb{E}[h(X)] = 0, 
    \hspace{10mm}
    \mathbb{E}\left[ h^\top(X) h(X)\right] < \infty
    $$
    Define the <i>covariance inner product</i> as:
    $$
    \langle h_1, h_2 \rangle = \mathbb{E}\left[ h_1^\top(X) h_2(X) \right],
    \hspace{10mm}
    h_1, h_2 \in \mathcal{H}
    $$
    and denote the norm induced by this inner product with:
    $$
    \rvert \rvert h \rvert \rvert = \langle h, h \rangle^{1/2}
    $$
    This space, denoted by $\mathcal{H}$, is a normed linear vector space with an inner product. By the $L_2$-completeness theorem, it is also complete, so it is a Hilbert space.
</div>

Our Hilbert spaces of interest (defined above) are dependent upon the true value $\theta_0$, and this is the value with respect to which we take the expectation in the inner product. Thus, the space will change if $\theta_0$ changes. 

We come to a theorem that will be helpful later in some definitions and results.

<div id="projection-theorem"></div>
<div class="theorem">
  <strong>Theorem (Projection Theorem).</strong>
  {% tabs project-theorem %}
  {% tab project-theorem statement %}
  Let $\mathcal{H}$ denote the Hilbert space <a href="#def-hilbert-space">defined above</a>, and let $\mathcal{U}$ denote a closed linear subspace. For any $h \in \mathcal{H}$, there exists a unique $u_0 \in \mathcal{U}$ such that:

  $$
  \rvert \rvert h - u_0 \rvert \rvert \leq \rvert \rvert h - u \rvert \rvert,
  \hspace{10mm}
  \forall u \in \mathcal{U}
  $$

  We call $u_0$ the <i>projection of $h$ onto $\mathcal{U}$</i>, and we usually use $\Pi (h \rvert \mathcal{U})$ to denote it. $\Pi(h \rvert \mathcal{U})$ satisfies:

  $$
  \langle \Pi(h \rvert \mathcal{U}), u \rangle = 0,
  \hspace{10mm} 
  \forall u \in \mathcal{U}
  $$
  {% endtab %}
  {% tab project-theorem proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

A version of the Pythagorean theorem and the Cauchy-Schwarz inequality can be derived as a result of the Projection theorem. 

---

## Extremum Estimators
Before we dive into more specifics, let's define a few different classes of estimators.

<div class="definition">
  <strong>Definition (Extremum Estimator).</strong>
  <br>
   Let $\hat{Q}_n(X, \theta)$ be some objective function the depends on some data, $X$, and a sample size, $n$. An <i>extremum estimator</i> is $\hat{\theta}_n$ such that:
   $$
   \hat{\theta}_n = \underset{\theta \in \Omega}{\arg \max} \left\{ \hat{Q}_n(X, \theta) \right\}
   $$
</div>

The name comes from the fact that we are maximizing the objective (i.e. finding an extreme value). An example is the maximum likelihood estimator (MLE), which has the (normalized) log-likelihood as its objective function:

$$
\hat{Q}_n(X, \theta) = \frac{1}{n} \sum_{i = 1}^n \log (p_X(x_i; \theta))
$$

In certain cases, the maximization can be done by taking the derivative and setting this equal to zero. However, we should be cautious as there are cases when there are many such solutions. In this case, a solution to this first-order condition may not be the global maximum of $\hat{Q}_n(X, \theta)$ (i.e. a solution may be a local maximum). Luckily, if the extremum estimator is consistent and $\theta_0$ is on the interior of $\Omega$, then the extremum estimator will be included in the set of solutions to setting the derivative equal to zero.

The extremum estimator class is quite broad. A subclass of extremum estimator is the $m$-estimator. 

<div id="m-estimator"></div>
<div class="definition">
  <strong>Definition ($m$-Estimator).</strong>
  <br>
    Let $m(X, \theta)$ be a $p \times 1$-dimensional function of $X$ satisfying:
    $$
    \mathbb{E}_\theta[m(X, \theta)] = \mathbf{0}_{p \times 1},
    \hspace{10mm}
    \mathbb{E}_\theta\left[ m^\top(X, \theta) m(X, \theta)\right] < \infty,
    \hspace{10mm}
    \mathbb{E}_\theta \left[ m(X, \theta) m^\top(X, \theta) \right] \text{ positive definite } \forall \theta \in \Omega
    $$
    The $m$-estimator of $\theta$ is the solution (if it exists) of:
    $$
    \sum_{i = 1}^n m(X_i, \hat{\theta}_n) = 0,
    \hspace{10mm}
    X_1, \dots, X_n \overset{iid}{\sim} p_X(x; \theta), 
    \hspace{2mm}
    \theta \in \Omega \subset \mathbb{R}^p
    $$
    Note: if $m(X, \theta)$ is differentiable, then take the derivative, set equal to $\mathbf{0}$, and solve. In this case, we call $\hat{\theta}_n$ an $m$-estimator of <i>$\psi$-type</i>. Otherwise, it is of <i>$\rho$-type</i>(see <a href="https://en.wikipedia.org/wiki/M-estimator">here</a>).
</div>

Put simply, an extremum estimator that maximizes a sample average is an $m$-estimator[^fn-newey]. Before we look at an example, let's define the <i>score</i>.

<div class="definition">
  <body>
  <strong>Definition (Score).</strong>
  <br>
  Suppose $X \sim p_X(x; \theta)$ where $\theta = (\beta^\top, \eta^\top)^\top$. The <i>score function</i>, $U_\theta(X, \theta)$, is defined as the vector of first order partial derivatives of the log-likelihood with respect to the parameter vector:
  $$
  U_\theta(x, \theta^*) = \frac{\partial \log(p_X(x, \theta))}{\partial \theta} \bigg\rvert_{\theta = \theta^*}
  $$
  where $\theta^*$ is some value of $\theta$. When evaluated at the true value of $\theta$ (the value that generated the data), we call $U_\theta(x, \theta_0)$ the <i>score vector</i>.
  </body>
</div>

The MLE is also an example of an $m$-estimator. Consider the function $m(X, \theta) = \log (p_X(x; \theta))$, the log-likelihood. We can maximize the log-likelihood by taking the derivative with respect to $\theta$ and setting this equal to zero. This is called the <i>score equation in $\theta$</i>:

$$
\sum_{i = 1}^n U_\theta(X_i, \theta) = 0
$$

### Influence Functions
Some estimators (asymptotically linear ones...to be defined later) can be analyzed with respect to a particular defining function, called the <i>influence function</i>.

<div class="definition">
  <body>
  <strong>Definition (Influence Function).</strong>
  <br>
    Let $\hat{\theta}_n$ be an estimator for parameter $\theta$. Define the $q$-dimensional measurable random function (i.e. a random vector) $\psi(X, \theta_0)$ such that $\mathbb{E}_{\theta_0}[\psi(X, \theta_0)] = \mathbf{0}$, $\mathbb{E}_{\theta_0}[\psi(X, \theta_0)\psi^\top(X, \theta_0)]$ exists, and:
    $$
    \sqrt{n}(\hat{\theta}_n - \theta) = \frac{1}{\sqrt{n}}\sum_{i = 1}^n \psi(X_i, \theta_0) + o_p(1)
    $$
    where $o_p(1)$ indicates a term that converges in probability to $0$. We call $\psi(X_i, \theta_0)$ the <i>$i$-th influence function of $\hat{\theta}_n$</i> or the <i>influence function of the $i$-th observations of $\hat{\theta}_n$</i>. 
  </body>
</div>

Note that the influence function is always defined with respect to the true data-generating distribution; it is evaluated at $\theta_0$, and the expectations are always taken with respect to the parameter being $\theta_0$. Thus, we don't really need to write $\psi(X, \theta_0)$, since the influence function is not a function of $\theta$ (it is just dependent on $\theta_0$). However, I keep the full notation to avoid confusion.

We can also define special properties of estimators to help us decribe the ones that are desirable. We say that $\hat{\theta}_n$ is <i>asymptotically linear</i> if it has an influence function(s). Thus, the maximum likelihood estimator is asymptotically linear. 

The influence function(s) provides us a clear way to analyze the behavior of our (asymptotically linear) estimator as $n \rightarrow \infty$:

$$
\begin{equation}
\label{eq:al-dist}
\begin{aligned}
\frac{1}{\sqrt{n}}\sum_{i = 1}^n \psi(X_i, \theta_0) 
&\rightsquigarrow \mathcal{N}(\mathbf{0}_{q \times q}, \mathbb{E}_{\theta_0}\left[ \psi(X, \theta_0) \psi(X, \theta_0)^\top \right]) & \left(\text{CLT}\right) \\
\sqrt{n}(\hat{\theta}_n - \theta_0)
&\rightsquigarrow \mathcal{N}(\mathbf{0}_{q \times q}, \mathbb{E}_{\theta_0}\left[ \psi(X, \theta_0) \psi(X, \theta_0)^\top \right]) & \left(\text{Slutsky's theorem}\right)
\end{aligned}
\end{equation}
$$

Eq. \eqref{eq:al-dist} implies that the asymptotic variance of the estimator is the variance of its influence function. Furthermore, an asymptotically linear estimator is (effectively) uniquely identified by its influence function.

<div class="theorem">
<strong>Theorem 3.1.<d-cite key=tsiatis2006></d-cite></strong>
{% tabs tsiatis-3-1 %}
{% tab tsiatis-3-1 statement %}
An asymptotically linear estimator has an almost surely unique influence function.
{% endtab %}
{% tab tsiatis-3-1 proof %}
Tsiatis proceeds by contradiction. Let our estimator be denoted by $\hat{\theta}_n$ with true value $\theta_0$, and let $n$ be our sample size. Let $\psi(X, \theta_0)$ be an influence function, and assume that there exists another influence function $\psi^*(X, \theta_0)$. Thus:

$$
\mathbb{E}_{\theta_0}[\psi(X, \theta_0)] = \mathbb{E}_{\theta_0}[\psi^*(X, \theta_0)] = 0, 
\hspace{10mm}
\sqrt{n}(\hat{\theta}_n - \theta_0) = \frac{1}{\sqrt{n}}\sum_{i = 1}^n \psi(X_i, \theta_0) + o_p(1) = \frac{1}{\sqrt{n}}\sum_{i = 1}^n \psi^*(X_i, \theta_0) + o_p(1)
$$

Recall that $X_1, \dots, X_n$ are i.i.d., so, by the central limit theorem:

$$
\frac{1}{\sqrt{n}}\sum_{i = 1}^n (\psi(X_i, \theta_0) - \psi^*(X_i, \theta_0)) \rightsquigarrow \mathcal{N}\left(\mathbf{0}, \mathbb{E}_{\theta_0}\left[ (\psi(X, \theta_0) - \psi^*(X, \theta_0))(\psi(X, \theta_0) - \psi^*(X, \theta_0))^\top\right] \right)
\nonumber
$$

However, by the <a href="https://en.wikipedia.org/wiki/Continuous_mapping_theorem">continuous mapping theorem</a>:

$$
\begin{aligned}
&\frac{1}{\sqrt{n}}\sum_{i = 1}^n (\psi(X_i, \theta_0) - \psi^*(X_i, \theta_0)) = o_p(1) \\
\implies 
&\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \bigg\rvert \frac{1}{\sqrt{n}}\sum_{i = 1}^n (\psi(X_i, \theta_0) - \psi^*(X_i, \theta_0)) \bigg\rvert \geq \epsilon \right) = 0, \hspace{5mm} \forall \epsilon > 0 
\end{aligned}
$$

For both of the above to be true, we need:

$$
\mathbb{E}_{\theta_0}\left[ (\psi(X, \theta_0) - \psi^*(X, \theta_0))(\psi(X, \theta_0) - \psi^*(X, \theta_0))^\top\right] = \mathbf{0}_{q \times q}
$$

which implies $\psi(X, \theta_0) = \psi^*(X, \theta_0)$ almost surely.
{% endtab %}
{% endtabs %}
</div>

The score vector also satisfies nice properties, which we summarize in the following theorem. First, let's define what it means for an estimator to be "regular".

<div class="definition">
  <strong>Definition (Regular Estimator).</strong>
  <br>
    Suppose we are considering a family of distributions, $\mathcal{P}$, indexed by parameter $\theta = (\beta^\top, \eta^\top)^\top$. Let $\theta^*$ be some fixed value for $\theta$. We refer to this set-up as a <i>local data generating process (LDGP)</i>.
    <br>
    For each value of $n$, we assume we have $n$ data points:

    $$
    X_{1, n}, X_{2, n}, \dots, X_{n, n} \overset{iid}{\sim} P(X, \theta_n)
    $$

    where $P(X, \theta_n) \in \mathcal{P}$ is the distribution from the family where $\theta = \theta_n$ and such that: 

    $$
    \underset{n \rightarrow \infty}{\lim} \sqrt{n}(\theta_n - \theta^*)  \rightarrow \mathbf{c}
    $$ 

    where $\mathbf{c}$ is some constant vector. This condition implies the two values of $\theta$ are "close".
    <br>
    An estimator $\hat{\theta}_n$ of $\theta$ is called <i>regular</i> if, for each $\theta^*$, the asymptotic distribution of $\sqrt{n}(\hat{\theta}_n - \theta)$ does not depend on the LDGP. That is, if the limiting distribution does not depend on $\theta_n$.
</div>

In many cases, the maximum likelihood estimator is regular. When an estimator is both asymptotically linear and regular, we say it is <i>RAL</i>. 

<div id="theorem-3-2"></div>
<div class="theorem">
<strong>Theorem 3.2, Corollary 1.<d-cite key=tsiatis2006></d-cite></strong>
{% tabs t-3-2-cor-1 %}
{% tab t-3-2-cor-1 statement %}
Let $\beta(\theta)$ be a $q$-dimensional function of $p$-dimensional parameter $\theta$ such that $q < p$. Assume the following exists:

$$
\begin{equation}
\label{eq:gamma}
\Gamma(\theta) = \frac{\partial \beta(\theta)}{\partial \theta^\top}
\end{equation}
$$

which is the $q \times p$ matrix of first order partial derivatives of vector $\beta(\theta)$ with respect to $\theta$, and assume it is continuous in $\theta$ in a neighborhood of $\theta_0$, the true parameter value. Let $$\hat{\beta}_n$$ denote an asymptotically linear estimator with influence function $\psi(X, \theta_0)$ such that $$\mathbb{E}_\theta[\psi^\top(X, \theta_0) \psi(X, \theta_0)]$$ exists and is also continuous in $\theta$ in a neighborhood of $\theta_0$. If $\hat{\beta}_n$ is regular, then:

$$
\begin{equation}
\label{eq:condition-3-2}
\mathbb{E}_{\theta_0}\left[ \psi(X,\theta_0) U_\theta^\top(X, \theta_0) \right] = \Gamma(\theta_0)
\end{equation}
$$

If the parameter space can be partitioned as $\theta = (\beta^\top, \eta^\top)^\top$ where $\beta$ is $q$-dimensional and $\eta$ is $r$-dimensional, then:

$$
\begin{equation}
\label{eq:corollary-1}
\mathbb{E}_{\theta_0}\left[ \psi(X, \theta_0) U_\beta^\top(X, \theta_0)\right] = \mathbb{I}_{q \times q},
\hspace{10mm}
\mathbb{E}_{\theta_0}\left[ \psi(X, \theta_0) U_\eta^\top(X, \theta_0) \right] = \mathbf{0}_{q \times r}
\end{equation}
$$

where:

$$
U_\beta(x, \theta_0) = \frac{\partial \log(p_X(x, \theta))}{\partial \beta} \bigg\rvert_{\theta = \theta_0}, 
\hspace{10mm}
U_\eta(x, \theta_0) = \frac{\partial \log(p_X(x, \theta))}{\partial \eta} \bigg\rvert_{\theta = \theta_0}, 
$$
{% endtab %}
{% tab t-3-2-cor-1 proof %}
See pg. 34 in <d-cite key=tsiatis2006></d-cite>.
{% endtab %}
{% endtabs %}
</div>

Note that the influence functions for RAL estimators are in the subspace of elements of $\mathcal{H}$ (the Hilbert space of mean-zero measurable random functions) satisfying Eq. \eqref{eq:corollary-1}. Similarly, each element in that subspace is the influence function for some RAL estimator. 

From this geometric perspective, we can see that the asymptotic variance of an RAL estimator is basically the squared distance between the origin and its influence function in our special Hilbert space.

---

## Consistency
We can show under what conditions extremum estimators are consistent ($\hat{\theta}_n \overset{p}{\rightarrow} \theta_0$). This is called the <i>Basic Consistency Theorem</i> by Newey and McFadden. In what follows, we introduce the probability limit, $Q_0(X, \theta)$, of $\hat{Q}_n(X, \theta)$, which is the quantity described in <a href="#uniform-convergence">the appendix</a>.

<div id="theorem-2-1"></div>
<div class="theorem">
<strong>Theorem 2.1.<d-cite key=newey1994></d-cite></strong>
{% tabs newey-2-1 %}
{% tab newey-2-1 statement %}
Let $\hat{Q}_n(X,\theta)$ be the objective function for an extremum estimator, $\hat{\theta}_n$. If there exists a function $Q_0(\theta)$ satisfying:

<ol>
<li>Identification: $Q_0(X,\theta)$ has a unique maximum at $\theta_0$, the true value of $\theta$</li>
<li>Boundedness: $\Omega$, the parameter space, is compact</li>
<li>Continuity: $Q_0(X,\theta)$ is continuous</li>
<li>Uniform Convergence: $\hat{Q}_n(X, \theta)$ converges uniformly in probability to $Q_0(X, \theta)$</li>
</ol>

then $\hat{\theta}_n$ is consistent; i.e.:

$$
  \hat{\theta}_n \overset{p}{\rightarrow} \theta_0
$$
{% endtab %}
{% tab newey-2-1 proof %}
See pg. 2121 in <d-cite key=newey1994></d-cite>
{% endtab %}
{% endtabs %}
</div>

The authors note that the some of the conditions can be relaxed. Instead of assuming that $\hat{\theta}_n$ maximizes $\hat{Q}_n(X, \theta)$, we can instead assume that it "nearly" maximizes it:

$$
\hat{Q}_n(\hat{\theta}_n) \geq \underset{\theta \in \Omega}{\sup} \hat{Q}_n(\theta) + o_p(1)
$$

The second condition can be relaxed if the objective function, $\hat{Q}_n(X, \theta)$, is concave. Then the assumption of compactness of the parameter space, $\Omega$, can be exchanged for just convexity. 

In addition, the third condition can be relaxed to just upper semi-continuity rather than continuity proper. That is, we assume that, for any $\theta \in \Omega$ and any $\epsilon > 0$, there exists an open subset $\mathcal{B} \subset \Omega$ such that $\theta \in \mathcal{B}$ and such that:

$$
Q_0(X, \theta') < Q_0(X, \theta) + \epsilon \hspace{5mm} \forall \theta' \in \mathcal{B}
$$

The fourth condition can be changed to just require that:

$$
\hat{Q}_n(X, \theta_0) \overset{p}{\rightarrow} Q_0(\theta_0), 
\hspace{5mm} \text{and} \hspace{5mm}
\hat{Q}_n(X, \theta) < Q_0(X, \theta) + \epsilon \hspace{5mm} \forall \epsilon > 0, \hspace{1mm} \forall \theta \in \Omega
$$

with probability approaching $1$. If we make the stronger assumption that:

$$
\underset{\theta \in \Omega}{\sup} \big\rvert \hat{Q}_n(X, \theta) - Q_0(X, \theta) \big\rvert \overset{as}{\rightarrow} 0
$$

instead of the fourth condition, then we have that $\hat{\theta}_n \overset{as}{\rightarrow} \theta_0$ (i.e. $\hat{\theta}_n$ is <i>strongly consistent</i>).


### Showing Consistency
In order to use <a href="#theorem-2-1">Theorem 2.1</a>, one must be able to show that the conditions (or their relaxations) hold. This can be difficult in practice, we often try to show some other property that are sufficient for the conditions. Newey and McFadden call these <i>primitive conditions</i>.

<div id="lemma-2-4"></div>
<div class="theorem">
<strong>Lemma 2.4.<d-cite key=newey1994></d-cite></strong>
{% tabs newey-2-4 %}
{% tab newey-2-4 statement %}
Let $A(X, \theta)$ be a matrix of functions of observation $X$ and parameter $\theta$. Let $\rvert \rvert A \rvert \rvert$ denote the Euclidean norm. Let $\Omega$ be a compact parameter space, and suppose our data are independent and identically distributed.

Suppose that each element of $A(X, \theta)$ is continuous at each $\theta \in \Omega$ with probability one, and suppose that there exists a function $d(X)$ such that $\rvert \rvert A(X, \theta) \rvert \rvert \leq d(X)$ for all $\theta \in \Omega$. Assume $\mathbb{E}[d(X)] < \infty$. Then:

$$
\mathbb{E}[A(X, \theta)] \text{ is continuous}
$$

and:

$$
\underset{\theta \in \Omega}{\sup} \left\vert \left\vert\frac{1}{n} \sum_{i = 1}^n \left(A(X_i, \theta) - \mathbb{E}[A(X, \theta)] \right) \right\vert \right\vert \overset{p}{\rightarrow} 0
$$
{% endtab %}
{% tab newey-2-4 proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>

The above lemma can be used for sample averages, which is exactly what we deal with in maximum likelihood estimation (and several other estimation settings). This leads us to the next theorem that states that, under certain conditions, the MLE is consistent:

<div id="lemma-2-5"></div>
<div class="theorem">
<strong>Lemma 2.5.<d-cite key=newey1994></d-cite></strong>
{% tabs newey-2-5 %}
{% tab newey-2-5 statement %}
Let $X_1, X_2, \dots$ are i.i.d. data and let $p_X(x_i, \theta_0)$ denote their probability density function. If the following conditions are met, then $\hat{\theta}_n \overset{p}{\rightarrow} \theta_0$ (the MLE is consistent):

<ol>
<li>$\theta \neq \theta_0$ implies $p_X(x_i, \theta) \neq p_X(x_i, \theta_0)$ for any $\theta \in \Omega$</li>
<li>$\theta_0 \in \Omega$ where $\Omega$ is compact</li>
<li>$\log(p_X(x_i, \theta))$ is continuous for all $\theta \in \Omega$ with probability one</li>
<li> $\mathbb{E}_{\theta_0}\left[ \underset{\theta \in \Omega}{\sup} \left\{ \left\vert \log(p_X(x, \theta)) \right\vert \right\} \right] < \infty$</li>
</ol>
{% endtab %}
{% tab newey-2-5 proof %}
The result follows from ensuring the conditions of <a href="#theorem-2-1">Theorem 2.1</a> are satisfied and then applying <a href="#lemma-2-4">Lemma 2.4</a>.
{% endtab %}
{% endtabs %}
</div>

---

## Asymptotic Normality
To construct confidence intervals, we often rely upon asymptotical normality of an estimator. 

<div id="theorem-3-1"></div>
<div class="theorem">
<strong>Theorem 3.1.<d-cite key=newey1994></d-cite></strong>
{% tabs newey-3-1 %}
{% tab newey-3-1 statement %}
Let $\hat{\theta}_n$ be an extremum estimator; that is, it maximizes some objective function $\hat{Q}_n(X, \theta)$ subject to $\theta \in \Omega$ give a sample size of $n$. If the following conditions are satisfied:

<ol>
<li>$\hat{\theta}_n \overset{p}{\rightarrow} \theta_0$</li>
<li>$\theta_0 \in \text{ interior of } \Omega$</li>
<li>$\hat{Q}_n(X, \theta)$ is twice continuously differentiable in a neighborhood $\mathcal{B}$ of $\theta_0$</li>
<li>$\sqrt{n} \frac{\partial \hat{Q}_n(X,\theta)}{\partial \theta} \bigg\rvert_{\theta = \theta_0} \rightsquigarrow \mathcal{N}(\mathbf{0}, \Sigma)$</li>
<li>There exists matrix function $H(\theta)$ that is continuous at $\theta_0$ and satisfies $\underset{\theta \in \mathcal{B}}{\sup} \left\vert \left\vert \frac{\partial^2 \hat{Q}_n(X, \theta)}{\partial \theta \partial \theta^\top} - H(\theta) \right\vert\right\vert \overset{p}{\rightarrow} 0$</li>
<li>$H = H(\theta_0)$ is non-singular</li>
</ol>

Then:

$$
\begin{equation}
\label{eq:asymptotic-dist}
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightsquigarrow \mathcal{N}\left(\mathbf{0}, H^{-1} \Sigma H^{-1}\right)
\end{equation}
$$

where $\Sigma$ is the asymptotic variance of $$\sqrt{n} \left[ \frac{\partial \hat{Q}_n(X, \theta)}{\partial \theta} \right] \bigg\rvert_{\theta = \theta_0}$$, and $$H = \underset{n \rightarrow \infty}{\text{plim}} \left[ \left. \frac{\partial^2 \hat{Q}_n(\theta)}{\partial \theta \partial \theta^\top} \right\vert_{\theta = \theta_0} \right]$$.
{% endtab %}
{% tab newey-3-1 proof %}
See Section 3.5 in <d-cite key=newey1994></d-cite>.
{% endtab %}
{% endtabs %}
</div>

Applied to maximum likelihood estimators, <a href="#theorem-3-1">Theorem 3.1</a> gives:

$$
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightsquigarrow \mathcal{N}(\mathbf{0}, \mathcal{I}^{-1}(\theta_0))
$$

<details>
<summary>Intuition.</summary>
The more expansive explanation of the intuition is provided at the beginning of Section 3 in Newey and McFadden, but we summarize it here for convenience.
<br>
For MLE, $\hat{Q}_n(\theta) = \frac{1}{n} \sum_{i = 1}^n \log(p_X(x_i, \theta))$. Suppose that the log-likelihood is differentiable and the maximum likelihood estimate, $\hat{\theta}$, is on the interior of the parameter space. This implies that the gradient of the log-likelihood will be equal to $0$ when evaluated at $\hat{\theta}_n$:
$$
\left. \left[ \frac{1}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta}\right] \right\vert_{\theta = \hat{\theta}_n} = 0
\nonumber
$$
Suppose also that the log-likelihood is twice differentiable and the second derivative is continuous. We can apply the <a href="https://en.wikipedia.org/wiki/Mean_value_theorem">mean value theorem</a> (MVT) to each element of the gradient to rewrite the first-order condition. Let $\bar{\theta}$ be the vector of intermediate points between each coordinate of $\hat{\theta}_n$ and $\theta_0$ that are used in the MVT. The $(i,j)$-th element is:
$$
\left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta_i \partial \theta_j} \right] \right\vert_{\theta = \bar{\theta}} 
= \frac{1}{(\hat{\theta}_n)_j - (\theta_0)_j} \left( \left. \left[ \frac{1}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \hat{\theta}_n}
- \left. \left[ \frac{1}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \theta_0}\right)
\nonumber
$$
Using the fact that the gradient of the log-likelihood evaluated at the maximum likelihood estimate is $0$:
$$
((\hat{\theta}_n)_j - (\theta_0)_j) 
\left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta_i \partial \theta_j} \right] \right\vert_{\theta = \bar{\theta}}  
= - \left. \left[ \frac{1}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \theta_0}
\nonumber
$$
And rearranging terms:
$$
(\hat{\theta_j} - (\theta_0)_j) 
\left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta_i \partial \theta_j} \right] \right\vert_{\theta = \bar{\theta}}  + \left. \left[ \frac{1}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \theta_0} = 0
\nonumber
$$
Rewriting in matrix notation:
$$
\left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta \partial \theta^\top} \right] \right\vert_{\theta = \bar{\theta}}(\hat{\theta}_n - (\theta_0))  + \left. \left[ \frac{1}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta}\right] \right\vert_{\theta = \theta_0} = \mathbf{0}
\nonumber
$$
Multiplying by $\sqrt{n}$ and rearranging again:
$$
\sqrt{n}(\hat{\theta_j} - (\theta_0)_j) 
= - \left[ \left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta \partial \theta^\top} \right] \right\vert_{\theta = \bar{\theta}} \right]^{-1}
\left. \left[ \frac{\sqrt{n}}{n}\sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \theta_0} 
\nonumber
$$
As the score, the second term has zero mean. In addition, it is the sum of (functions of) i.i.d. random variables, so we can apply the central limit theorem to see that:
$$
\left. \left[ \frac{\sqrt{n}}{n} \sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \theta_0} 
\rightsquigarrow
\mathcal{N}\left(\mathbf{0}, \mathcal{I}(\theta_0)\right)
\nonumber
$$
where $\mathcal{I}(\theta_0) = \text{Var}\left( \left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial \log(p_X(x_i, \theta))}{\partial \theta_i}\right] \right\vert_{\theta = \theta_0} \right)$.
<br>
Since $\bar{\theta}$ is just a function of $\hat{\theta}_n$ (since $\theta_0$ is fixed), if $\hat{\theta}_n$ is consistent, then $\bar{\theta}$ should also be. Supposing we can use a uniform law of large numbers (see <a href="#lemma-2-4">Lemma 2.4</a>), then:
$$
\left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta \partial \theta^\top} \right] \right\vert_{\theta = \bar{\theta}}
\overset{p}{\rightarrow} \mathbb{E}_{\theta_0}\left[ \left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta \partial \theta^\top} \right] \right\vert_{\theta = \theta_0} \right] = H
\nonumber
$$
where the righthand side is the expected Hessian, denoted by $H$. 
<br>
Since we assumed the log-likelihood is twice continuously differentiable, the <a href="https://en.wikipedia.org/wiki/Inverse_function_theorem">inverse function theorem</a> gives us:
$$
-\left[\left. \left[ \frac{1}{n} \sum_{i = 1}^n \frac{\partial^2 \log(p_X(x_i, \theta))}{\partial \theta \partial \theta^\top} \right] \right\vert_{\theta = \bar{\theta}} \right]^{-1}
\overset{p}{\rightarrow} -H^{-1}
\nonumber
$$
Since $-H^{-1}$ is a constant, by <a href="https://en.wikipedia.org/wiki/Slutsky%27s_theorem">Slutsky's theorem</a>:
$$
\sqrt{n}(\hat{\theta} - \theta_0) \rightsquigarrow \mathcal{N}\left(\mathbf{0}, H^{-1}\mathcal{I}(\theta_0) H^{-1}\right)
$$
Under our regularity conditions, the negative expected Hessian is the variance of the score, so:
$$
-H^{-1} = \mathcal{I}^{-1}(\theta_0) \implies \mathcal{N}\left(\mathbf{0}, H^{-1}\mathcal{I}(\theta_0) H^{-1}\right) = \mathcal{N}\left(\mathbf{0}, \mathcal{I}^{-1}(\theta_0)\right)
\nonumber
$$
</details>

Newey and McFadden make some important points about the conditions stated in <a href="#theorem-3-1">Theorem 3.1</a>. Perhaps most notable is the need for the true parameter value to be located in the interior of the parameter space. When this condition is not met ($\hat{\theta}_n$ is on the boundary even as $n \rightarrow \infty$), the asymptotic normality result is no longer guaranteed (though it could still be true). This comes up in variance estimation, since we usually constrain our space to positive reals. We also need the average score over the sample to satisfy a central limit theorem (this gives us the "base" Normal distribution) and the inverse Hessian needs to converge to a constant (so we can apply Slutsky's theorem). 

### Estimating Asymptotic Variance
For confidence intervals, we need a consistent estimator of the asymptotic variance, $H^{-1} \Sigma H^{-1}$, of the estimator. Usually we do this by estimating the components and then plugging these in; i.e. we find $\hat{H}^{-1}$ and $\hat{\Sigma}$ and use $\hat{H}^{-1} \hat{\Sigma} \hat{H}^{-1}$

<div id="theorem-4-1"></div>
<div class="theorem">
<strong>Theorem 4.1.<d-cite key=newey1994></d-cite></strong>
{% tabs newey-4-1 %}
{% tab newey-4-1 statement %}
Under the conditions of <a href="#theorem-3-1">Theorem 3.1</a>, if $$\hat{H} = \left. \left[ \frac{\partial \hat{Q}_n(\theta)}{\partial \theta}\right]\right\vert_{\theta = \hat{\theta}}$ and $\hat{\Sigma} \overset{p}{\rightarrow} \Sigma$$, then:

$$
\hat{H}^{-1} \hat{\Sigma} \hat{H}^{-1} \overset{p}{\rightarrow} H^{-1} \Sigma H^{-1}
$$
{% endtab %}
{% tab newey-4-1 proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>

This estimator is sometimes called a <i>sandwich estimator</i> since we have the asymptotic variance "sandwiched" between the same term twice.

For maximum likelihood estimators, one can estimate the asymptotic variance with any consistent estimator of the inverse Fisher information matrix. Some examples provided by Newey and McFadden are to use the negative Hessian (if regularity conditions are met) or the sample average of the outer product of the score. However, one should use caution because, under  model misspecification, estimators of the inverse information matrix may not be consistent. 

---

## Efficiency
We can compare different asymptotically normal estimators by their <i>efficiency</i>. But before we begin, we need to define a few new quantities.

<div id="tangent-space"></div>
<div class="definition">
  <strong>Definition (Tangent Space).<d-cite key=tsiatis2006></d-cite></strong>
  <br>
  Let $\mathcal{H}$ be the Hilbert space of $q$-dimensional measurable functions of $X$ with mean zero and finite variance and equipped with inner product $\langle h_1, h_2 \rangle = \mathbb{E}_{\theta_0}[ h_1^\top h_2]$ as we defined <a href="#def-hilbert-space">above</a>.
  <br>
  The <i>tangent space</i> is the linear subspace of $\mathcal{H}$ spanned by the score vector $U_\theta(X, \theta_0)$:

  $$
  \mathcal{T} = \left\{ B U_\theta(X, \theta_0) : B  \in \mathbb{R}^{q \times p} \right\}
  $$
</div>

We can define a similar space when the parameter vector can be partitioned as $\theta = (\beta^\top, \eta^\top)^\top$. The <i>nuisance tangent space</i> is the linear subspace of $\mathcal{H}$ spanned by the nuisance score vector $U_\eta(X, \theta_0)$:

$$
\begin{equation}
\label{eq:nuisance-space}
\Lambda = \left\{ B U_\eta(X, \theta_0): B \in \mathbb{R}^{q \times r}\right\}
\end{equation}
$$

The tangent space generated by $U_\beta(X, \theta_0)$ is:

$$
\begin{equation}
\label{eq:interest-space}
\mathcal{T}_\beta = \left\{ B U_\beta(X, \theta_0): B \in \mathbb{R}^{q \times p} \right\}
\end{equation}
$$

Notably, the direct sum of these two spaces equals the tangent space generated by the entire score vector:

$$
\mathcal{T} = \mathcal{T}_\beta \oplus \Lambda
$$

We can also construct the set of elements that are orthogonal to $\Lambda$ are then given by $h - \Pi(h \rvert \Lambda)$ for all $h \in \mathcal{H}$ where $\Pi(h \rvert \Lambda)$ is the <i>residual of $h$ after projecting onto $\Lambda$.

<div id="residual"></div>
<div class="definition">
  <strong>Definition (Residual).<d-cite key=newey1994></d-cite></strong>
  <br>
  Let $h \in \mathcal{H}$ be an element in the <a href="#def-hilbert-space">Hilbert space defined above</a>. By the <a href="#projection-theorem">projection theorem</a>, there exists a unique element $a_0 \in \Lambda$ in the <a href="#nuisance-tangent-space">nuisance tangent space</a> such that:

  $$
  \rvert \rvert h - a_0 \rvert \rvert \leq \rvert \rvert h - a \rvert \rvert \hspace{3mm} \text{ and } \hspace{3mm} \langle h - a_0, a \rangle \hspace{10mm} \forall a \in \Lambda
  $$

  The element $h - a_0$ is called the <i>residual of $h$ after projecting onto $\Lambda$</i> and is equal to:

  $$
  h - a_0 = \Pi(h \rvert \Lambda^\perp)
  $$
</div>

### Efficient Influence Functions
We call the influence function with the smallest variance matrix the <i>efficient influence function</i>. In the case of RAL estimators, we say that $\psi^{(1)}(X, \theta_0)$ has <i>smaller (asymptotic) variance</i> than $\psi^{(2)}(X, \theta_0)$ if, for all $q \times 1$ constant vectors $a$:

$$
\text{Var}(\psi^{(1)}(X, \theta_0)) \leq \text{Var}(\psi^{(2)}(X, \theta_0)) \iff \text{Var}(a^\top \psi^{(1)}(X, \theta_0)) \leq \text{Var}(a^\top \psi^{(2)}(X, \theta_0))
$$

<details>
<summary>Implications.</summary>
Since influence functions have zero mean, the righthand side of the above condition is equivalent to:

$$
a^\top \mathbb{E}_{\theta_0}[\psi^{(1)}(X, \theta_0) (\psi^{(1)}(X, \theta_0))^\top]a 
\leq a^\top \mathbb{E}_{\theta_0}[\psi^{(2)}(X, \theta_0) (\psi^{(2)}(X, \theta_0))^\top]a
$$

for all $q \times 1$ constant vectors $a$, which itself is equivalent to:

$$
a^\top \left[ \mathbb{E}_{\theta_0}[\psi^{(2)}(X, \theta_0) (\psi^{(2)}(X, \theta_0))^\top] - \mathbb{E}_{\theta_0}[\psi^{(1)}(X, \theta_0) (\psi^{(1)}(X, \theta_0))^\top] \right] a \geq 0
$$
</details>

This brings us to a theorem that defines the subspace of influence functions. 

<div id="theorem-3-4"></div>
<div class="theorem">
<strong>Theorem 3.4.<d-cite key=tsiatis2006></d-cite></strong>
{% tabs tsiatis-3-4 %}
{% tab tsiatis-3-4 statement %}
The set of all influence functions is the linear variety $$\psi^*(X, \theta_0) + \mathcal{T}^\perp$$ where $$\psi^*(X, \theta_0)$$ is any influence function, and $$\mathcal{T}^\perp$$ is the space perpendicular to the tangent space.
{% endtab %}
{% tab tsiatis-3-4 proof %}
See pg. 45-46 of <d-cite key=tsiatis2006></d-cite>.
{% endtab %}
{% endtabs %}
</div>

This result states that we can construct all influence functions of RAL estimators by taking an arbitrary influence function and adding any element from the orthogonal complement of the tangent space to it. <a href="#theorem-3-4">Theorem 3.4</a> can be used to define the efficient influence function.

<div id="theorem-3-5"></div>
<div class="theorem">
<strong>Theorem 3.5.<d-cite key=tsiatis2006></d-cite></strong>
{% tabs tsiatis-3-5 %}
{% tab tsiatis-3-5 statement %}
Let $$\psi^*(X, \theta_0)$$ be any influence function, and let $\mathcal{T}$ be the tangent space generated by the score vector. The <i>efficient influence function</i> is given by the projection of $\psi^*(X, \theta_0)$ onto the tangent space:

$$
\psi_{\text{eff}}(X, \theta_0) = \psi^*(X, \theta_0) - \Pi(\psi^*(X, \theta_0) \rvert \mathcal{T}^\perp) = \Pi(\psi^*(X, \theta_0) \rvert \mathcal{T})
$$

The efficient influence function can be written as:

$$
\psi_{\text{eff}}(X, \theta_0) = \Gamma(\theta_0) \mathcal{I}^{-1}(\theta_0) U_\theta(X, \theta_0)
$$

where $\Gamma(\theta_0)$ is the matrix defined in Eq. \eqref{eq:gamma} and $\mathcal{I}(\theta_0)$ is the information matrix.
{% endtab %}
{% tab tsiatis-3-5 proof %}
See pg. 46-47 in <d-cite key=tsiatis2006></d-cite>.
{% endtab %}
{% endtabs %}
</div>

### Efficient Score
We can define an efficient version of the score if we are in the setting where we can partition $\theta = (\beta^\top, \eta^\top)^\top$.

<div id="efficient-score"></div>
<div class="definition">
  <strong>Definition (Efficient Score).<d-cite key=tsiatis2006></d-cite></strong>
  <br>
  The <i>efficient score</i> is the residual of the score with respect to $\beta$ after projecting onto the nuisance tangent space:

  $$
  \begin{aligned}
  U_{\text{eff}}(X, \theta_0) &= \Pi(U_\beta(X, \theta_0) \rvert \Lambda^{\perp}) \\
  &= U_\beta(X, \theta_0) - \Pi(U_\beta(X, \theta_0) \rvert \Lambda) \\
  &= U_\beta(X, \theta_0) - \mathbb{E}_{\theta_0}\left[ U_\beta(X, \theta_0) U^\top_\eta(X, \theta_0)\right] \left[ \mathbb{E}\left[ U_\eta(X, \theta_0) U^\top_\eta(X, \theta_0) \right] \right]^{-1} U_\eta(X, \theta_0)
  \end{aligned}
  $$ 
</div>

In this setting, we can construct the efficient influence function as:

$$
\begin{equation}
\label{eq:partitioned-efficient-influence}
\psi_\text{eff}(X, \theta_0) = \left[ \mathbb{E}_{\theta_0}\left[ U_\text{eff}(X, \theta_0) U^\top_\text{eff}(X, \theta_0) \right] \right]^{-1} U_\text{eff}(X, \theta_0)
\end{equation}
$$

and has variance equal to:

$$
\text{Var}(\psi_{\text{eff}}(X, \theta_0)) = \left[ \mathbb{E}_{\theta_0} \left[ U_\text{eff}(X, \theta_0) U^\top(X, \theta_0) \right] \right]^{-1}
$$

which is the inverse variance matrix of the efficient score. If we partition the variance of the score vector as:

$$
\text{Var}(U_\theta(X, \theta_0)) = \mathcal{I} = 
\begin{bmatrix}
\mathcal{I}_{\beta, \beta} = \mathbb{E}_{\theta_0}\left[ U_\beta(X, \theta_0) U^\top_\beta(X, \theta_0) \right]
& \mathcal{I}_{\beta, \eta} = \mathbb{E}_{\theta_0}\left[ U_\beta(X, \theta_0) U^\top_\eta(X, \theta_0) \right] \\
\mathcal{I}_{\eta, \beta} = \mathbb{E}_{\theta_0}\left[ U_\eta(X, \theta_0) U^\top_\beta(X, \theta_0) \right] 
& \mathcal{I}_{\eta, \eta} = \mathbb{E}_{\theta_0}\left[ U_\eta(X, \theta_0) U^\top_\eta(X,  \theta_0) \right]
\end{bmatrix}
$$

Then we can use the <a href="https://en.wikipedia.org/wiki/Schur_complement">Schur complement formula</a> to get the variance of the efficient influence function is:

$$
\begin{equation}
\label{eq:var-eff-influence}
\text{Var}(\psi_{\text{eff}}(X, \theta_0)) = \left[\mathcal{I}_{\beta, \beta} - \mathcal{I}_{\beta, \eta} \mathcal{I}^{-1}_{\eta, \eta}  \mathcal{I}^\top_{\beta, \eta} \right]^{-1}
\end{equation}
$$

---

## Appendix 

### Auxiliary Definitions

<div id="uniform-convergence"></div>
<div class="definition">
<strong>Definition (Uniform Convergence In Probability).<d-cite key=newey1994></d-cite></strong>
<br>
Let $\hat{Q}_n(X,\theta)$ be the objective function for an extremum estimator, $\hat{\theta}_n$. We say that $\hat{Q}_n(X, \theta)$ converges <i>uniformly in probability</i> to $Q_0(X,\theta)$ if:

$$
\underset{\theta \in \Omega}{\sup} \big\rvert \hat{Q}_n(X, \theta) - Q_0(X,\theta) \rvert \overset{p}{\rightarrow} 0 
$$
</div>


<div id="direct-sum"></div>
<div class="definition">
<strong>Definition (Direct Sum).<d-cite key=tsiatis2006></d-cite></strong>
<br>
Let $M, N \subset \mathcal{H}$ be linear subspaces of Hilbert space $\mathcal{H}$. Their <i>direct sum</i>, denoted by $M \oplus N$, is the linear subspace in $\mathcal{H}$ such that every $x \in M \oplus N$ can be uniquely represented as $x = m + n$ for some $m \in M$ and $n \in N$. 
</div>

<div id="orthogonal-complement"></div>
<div class="definition">
<strong>Definition (Orthogonal Complement).<d-cite key=tsiatis2006></d-cite></strong>
<br>
  Let $M \subset \mathcal{H}$ be a linear subspace of Hilbert space $\mathcal{H}$. The <i>orthogonal complement</i> of $M$, denoted by $M^\perp$, is the linear subspace of elements in $\mathcal{H}$ that are orthogonal to $M$. It also holds that:
  
  $$
  \mathcal{H} = M \oplus M^\perp
  $$
</div>

<div id="linear-variety"></div>
<div class="definition">
<strong>Definition (Linear Variety).<d-cite key=tsiatis2006></d-cite></strong>
<br>
Let $M$ be a linear subspace of Hilbert space $\mathcal{H}$. A <i>linear variety</i> or <i>affine space</i> is the set:
$$
V = x_0 + M
$$
for all $x_0 \in \mathcal{H}$ such that $x_0 \not\in M$ and $\rvert \rvert x_0 \rvert \rvert \neq 0$. 
<br>
A linear variety is simply a shifting of a linear subspace.
</div>

### Auxiliary Results

The following result is called the <i>Information inequality</i> and states that if $\theta_0$ can be identified using a MLE, then the limiting objective function will have a unique maximum at the true value.

<div id="information-inequality"></div>
<div class="theorem">
<strong>Lemma 2.2. <d-cite key=newey1994></d-cite></strong>
{% tabs newey-2-2 %}
{% tab newey-2-2 statement %}
  If, for any $\theta \in \Omega$ such that $\theta \neq \theta_0$ implies that $p_X(x, \theta) \neq p_X(x, \theta_0)$ (i.e. $\theta_0$ is identified) and $\mathbb{E}[\rvert \log p_X(x, \theta) \rvert ] < \infty$ for all $\theta$, then the limiting objective function $Q_0(\theta) = \mathbb{E}\left[ \log(p_X(x, \theta)) \right]$ has a unique maximum at $\theta_0$.
{% endtab %}
{% tab newey-2-2 proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>
