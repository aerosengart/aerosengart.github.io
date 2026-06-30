---
layout: distill
title: Maximum Likelihood Estimation
description: A Primer
date: 2026-06-30
tabs: true
tags: theory likelihood primer
toc:  
  - name: Background
    subsections:
        - name: Notation
        - name: Definitions
  - name: Regularity Conditions
  - name: Results
bibliography: stats-ml.bib
---

This post mostly covers Chapter 7 of <i>Elements of Large-Sample Theory</i> by Erich Lehmann<d-cite key=lehmann2004></d-cite> with some supplemental material from <i>Theory of Statistics</i> by Mark Schervish,<d-cite key=schervish1995></d-cite> and <i>Theory of Point Estimation</i> by Lehmann and Casella.<d-cite key=lehmann2005></d-cite>


---

## Background 

### Notations
Suppose we have some probability space $(S, \mathcal{A}, \mu)$. A random variable is some function $X: S \rightarrow \mathcal{X}$ where $\mathcal{X}$ is the sample space (which we only require to be a <a href="/measure-theory/chapter2/#borel-set">Borel space</a> but is usually a subset of Euclidean space) with <a href="/maesure-theory/chapter2/#sigma-field">$\sigma$-field</a> $\mathcal{B}$. Individual elements of $\mathcal{X}$ are denoted with $x$.

Let $$\mathcal{P}_\Theta$$ be a family of distributions for $X$ parametrized by $\Theta: S \rightarrow \Omega$ where $\Omega$ is the parameter space with $\sigma$-field $\tau$. We allow $\Theta$ to be scalar- or vector-valued. Here, we assume it is $k$-dimensional. Individual elements of $\Omega$ are denoted with $\theta$. Denote the conditional distribution of $X$ given $\Theta = \theta$ with $P_\theta$ (which is a distribution on $(\mathcal{X},\mathcal{B})$). 

To make our notation match less measure theoretical texts, we'll use $X$ to denote a random variable with realizations denoted with the lowercase $x$. A parameter will be denoted with $\Theta$ with particular values denoted by $\theta$ and its true value by $\theta^*$. The probability density or mass function of $X$ given parameter $\Theta$ evaluated at a particular $x$ and $\theta$ will be denoted by $f_{X \rvert \Theta}(x; \theta)$ or, more compactly, $f(x; \theta)$. We will also use the notation $\mathbf{X} = (X_1, \dots, X_n)$ and $\mathbf{x} = (x_1, \dots, x_n)$ to denote a collection of $n$ random variables or vectors $X$ and observations of each, respectively. Unless otherwise stated, we assume $X_1, \dots, X_n$ are i.i.d.

### Definitions
First, we need to lay out some definitions and background information for our discussion. The <strong>likelihood</strong> of $\theta$ given $\mathbf{X} = \mathbf{x}$ is $f(\mathbf{x}; \theta)$ viewed as a function of $\theta$:

$$
\mathcal{L}(\theta; \mathbf{x}) = f_{\mathbf{X} \rvert \Theta}(\mathbf{x}; \theta) = \prod_{i = 1}^n f_{X_i \rvert \Theta}(x_i; \theta)
$$

where the last equality follows from the independence assumption. We also often work with the natural logarithm of $\mathcal{L}(\theta; \mathbf{x})$, denoted by $\ell(\theta; \mathbf{x})$, and called the <strong>log-likelihood</strong>. 

<aside><p>Note: the likelihood is also a function of the data!</p></aside>

We can estimate $\theta$ by finding the value that maximizes $\mathcal{L}(\theta; \mathbf{x})$, which is called the <strong>maximum likelihood estimate (MLE)</strong>. Intuitively, the MLE is the value that is most probable when it comes to having generated the data at hand.

<div id="ml-estimator"></div>
<div class="definition">
<strong>Definition (Maximum Likelihood Estimator; Schervish, pg. 3).</strong>
<br>
Let $\phi$ be an estimator of $\Theta$. We call $\phi$ a <i>maximum likelihood estimator (MLE)</i> if it satisfies, for all $x \in \mathcal{X}$:

$$
\underset{\theta \in \Omega}{\sup} \left\{ \mathcal{L}(\theta; x) \right\} = \mathcal{L}(\phi(x); x)
$$
</div>

Though we do not know the true value of $\Theta$, we can try to find the MLE. This choice of estimator makes sense because of the following theorem.

<div class="theorem">
  <strong>Theorem 3.2.<d-cite key=lehmann2005></d-cite></strong>
  {% tabs thrm-3-2 %}
  {% tab thrm-3-2 statement %}
  Assume Conditions 1-4 (given <a href="#regularity-conditions">below</a>) hold. Then:

  $$
  \mathbb{P}_{\theta^*}\left( \mathcal{L}(\theta^*; \mathbf{X}) > \mathcal{L}(\theta; \mathbf{X}) \right) \rightarrow 1 \hspace{5mm} \text{ as } \hspace{5mm} n \rightarrow \infty
  $$

  for any fixed $\theta \neq \theta^*$.
  {% endtab %}
  {% tab thrm-3-2 proof %}
  Fix $\theta \in \Omega$. Note that:

  $$
  \begin{aligned}
  &\mathcal{L}(\theta^*; \mathbf{X}) > \mathcal{L}(\theta; \mathbf{X})  \\
  \implies &\log(f(\mathbf{X}; \theta^*)) > \log(f(\mathbf{X}; \theta)) \\ 
  \implies &\log(f(\mathbf{X}; \theta)) - \log(f(\mathbf{X}; \theta^*)) < 0 \\
  \implies &\log\left(\prod_{i = 1}^n f(X_i; \theta) \right) - \log\left( \prod_{i = 1}^n f(X_i; \theta^*) \right) < 0 \\
  \implies &\sum_{i = 1}^n \left[ \log\left( \frac{f(X_i; \theta)}{f(X_i; \theta^*)} \right) \right] < 0 \\
  \implies &\frac{1}{n} \sum_{i = 1}^n \log\left( \frac{f(X_i; \theta)}{f(X_i; \theta^*)} \right) < 0
  \end{aligned}
  $$

  By the <a href="https://en.wikipedia.org/wiki/Law_of_large_numbers">law of large numbers</a>:

  $$
  \frac{1}{n} \sum_{i = 1}^n \log\left( \frac{f(X_i; \theta)}{f(X_i; \theta^*)} \right) \overset{p}{\rightarrow} \mathbb{E}_{\theta^*}\left[ \log\left( \frac{f(X; \theta)}{f(X; \theta^*)} \right) \right]
  $$

  Since the logarithm is strictly concave, an application of <a href="https://en.wikipedia.org/wiki/Jensen%27s_inequality">Jensen's inequality</a> yields:

  $$
  \mathbb{E}_{\theta^*}\left[ \log\left( \frac{f(X; \theta)}{f(X; \theta^*)} \right) \right] < \log\left(  \mathbb{E}_{\theta^*}\left[\frac{f(X; \theta)}{f(X; \theta^*)} \right] \right) \overset{(i)}{=} 0
  $$

  I do not immediately see how to get $(i)$, but if the expectation is $1$, then it follows.
  {% endtab %}
  {% endtabs %}
</div>

The MLE is a function $\hat{\theta}: \mathcal{X} \rightarrow \Omega$ mapping from the sample space to the parameter space. 
MLEs exhibit the <i>invariance property</i>, which is, in words, that a function of an MLE is the MLE of that function.

<div class="theorem">
<strong>Invariance Property of Maximum Likelihood Estimators.</strong>
{% tabs invariance-prop %}
{% tab invariance-prop statement %}
Let $\hat{\theta}$ be an MLE of $\Theta$, and let $g$ be some function of $\theta$. Then $g(\hat{\theta})$ is an MLE of $g(\Theta)$. 
{% endtab %}
{% tab invariance-prop proof %}
Define the <i>induced likelihood function</i>:
  
$$
\mathcal{L}^*(\eta; x) = \underset{\theta: g(\theta) = \eta}{\sup} \left\{ \mathcal{L}(\theta; x) \right\}
$$

which is a function of $\eta$ equal to the maximum value of the likelihood function over all values of $\theta$ such that $g(\theta) = \eta$. Let:

$$
\hat{\eta} = \underset{\eta}{\arg\sup}\left\{ \mathcal{L}^*(\eta; x) \right\};
\hspace{5mm}
\hat{\theta} = \underset{\theta}{\arg\sup}\left\{ \mathcal{L}(\theta; x) \right\}
$$

We have:

$$
\begin{aligned}
  \mathcal{L}^*(\hat{\eta}; x)
  &= \underset{\eta}{\sup}\left\{ \mathcal{L}^*(\eta; x) \right\} \\
  &= \underset{\eta}{\sup}\left\{ \underset{\theta: g(\theta) = \eta}{\sup} \left\{ \mathcal{L}(\theta; x) \right\} \right\} \\
  &= \underset{\theta}{\sup}\left\{ \mathcal{L}(\theta; x) \right\}  \\
  &= \mathcal{L}(\hat{\theta}; x) \\
  &= \underset{\theta: g(\theta) = g(\hat{\theta})}{\sup} \left\{ \mathcal{L}(\theta; x) \right\} \\
  &= \mathcal{L}^*(\hat{\theta}; x)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Unfortunately, the MLE is not guaranteed to be unique or even to exist (see Ex. 7.1.2<d-cite key=lehmann2004></d-cite>). To make the problem a bit easier, we relax our request to just a <i>local</i> maximum of the likelihood, but we desire our estimator to have some nice properties. One of which is <i>consistency</i>. 

<div class="definition">
<strong>Definition (Consistent Estimator).</strong>
<br>
Let $X_1, X_2, \dots$ be an infinite sample where $X_i \sim P_{\theta}$ with $\theta \in \Omega$, and let $\{ T_n(\theta) \}_{n = 1}^\infty$ be a sequence of estimators of some function of $\theta$, $g(\theta)$ ($g(\cdot)$ could be the identity function). The sequence of estimators is called <i>consistent</i> if, for all $\theta \in \Omega$:

$$
\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \rvert \rvert T_n(\theta) - g(\theta) \rvert \rvert_2 > \epsilon \right) = \mathbf{0}_k
$$

for all $\epsilon > 0$. In other words, $T_n(\theta) = g(\theta)+ o_p(1)$. (Technically, this is a <i>weakly</i> consistent estimator).
</div>

In order to establish the existence of a sequence of estimates that are local maxima of the likelihood, we require several conditions to be met. We begin with a couple of definitions.

<div class="definition">
<strong>Definition (Score; Schervish, pg. 111).</strong>
<br>
The <i>score function</i> (or <i>score statistic</i> or <i>efficient score</i>) is the $k$-dimensional vector with $i$-th coordinate given by:

$$
\begin{equation}
\label{eq:score}
    U_{\Theta_i}(\theta; x) = \left. \frac{\partial \ell(\Theta; x)}{\partial \Theta_i} \right\rvert_{\Theta = \theta}
\end{equation}
$$
</div>

We also have the <i>information matrix</i>, which quantifies how much our data tells us about the parameter.

<div class="definition">
<strong>Definition (Fisher Information Matrix; Schervish, pg. 111).</strong>
<br>
The <i>Fisher information matrix</i> is defined as the $k \times k$ matrix with $(i,j)$-th element given by:

$$
\begin{equation}
\label{eq:info-mat}
\mathcal{I}_{i,j}(\theta) 
    = \text{Cov}_{\theta}\left( \frac{\partial \ell(\Theta; X)}{\partial \Theta_i}, \frac{\partial \ell(\Theta; X)}{\partial \Theta_j} \right)
    = \text{Cov}_{\theta}\left( U_{\Theta_i}(\theta; X), U_{\Theta_j}(\theta; X) \right)
\end{equation}
$$

where the expectation is with respect to $X \sim P_{\theta}$ (i.e. assuming $\Theta = \theta$). 
</div>

Both definitions can be extended to the case where we condition on a statistic. These are called the <i>conditional score function</i> and <i>conditional Fisher information</i>, and are denoted by $U_{\Theta \rvert T}(\theta; x)$ and $\mathcal{I}_{\Theta \rvert T}(\theta)$ for statistic $T$, respectively (see pg. 111 and Section 2.3.3 in Schervish).

---

## Regularity Conditions
Below are the conditions needed for stating some of the results discussed later in this post. Note that not all of them are needed for each result. 

#### Condition 1
If $$P_{\theta_1} = P_{\theta_2}$$ for $\theta_1, \theta_2 \in \Omega$, then $\theta_1 = \theta_2$.

In words, this implies that the parameter is identifiable and that the distributions in the class $\mathcal{P}_{\Theta}$ are distinct.

#### Condition 2
For some $\delta > 0$, there exists a set, $N$, defined as:

$$
N = \left\{\theta \in \mathbb{R}^k \hspace{1mm} : \hspace{2mm} \rvert \rvert \theta - \theta^* \rvert \rvert_2 < \delta \right\}
$$

such that $N \subseteq \Omega$. 

<aside><p>$N$ is called a <strong>neighborhood</strong>.</p></aside>

This condition can be replaced by the stronger one (but easier to check) that $\Omega$ is open. In other words, we check that it holds for all $\theta \in \Omega$ since we do not know what $\theta^*$ is. 

#### Condition 3
We observe $n$ i.i.d. random variables (or vectors), $X_1, \dots, X_n$, which have <i>continuous</i> probability density or mass function $f_{X \rvert \Theta}(x; \theta) = f(x; \theta)$. 

#### Condition 4
The support of $f(x; \theta)$, defined as:

$$
\mathbf{A} = \left\{ x \hspace{1mm} : \hspace{1mm} f(x; \theta) > 0 \right\}
$$

is independent of $\theta$ (and the same for all $X_i$?).

#### Condition 5
The first-order partial derivatives of $f(x; \theta)$ with respect to $\Theta$, denoted by $\frac{\partial f(x; \theta)}{\partial \Theta_i}$, exist for all $x \in \mathbf{A}$ and $i = 1, \dots, k$. 

Note: Lehmann and Casella only require that the partial derivatives (in 5, 5a, and 5b) exist <i>almost everywhere</i> (pg. 462).<d-cite key=lehmann2005></d-cite>

#### Condition 5a
The second-order partial derivatives of $f(x; \theta)$ with respect to $\theta$, denoted by $\frac{\partial^2 f(x; \theta)}{\partial \Theta_i \partial \Theta_j}$, exist for all $x \in \mathbf{A}$, $\theta \in \Omega$, and $i,j = 1, \dots, k$.

#### Condition 5b
The third-order partial derivatives of $f(x; \theta)$ with respect to $\theta$, denoted by $\frac{\partial^3 f(x; \theta)}{\partial \Theta_i \partial \Theta_j \partial \Theta_l}$, exist for all $x \in \mathbf{A}$, $\theta \in \Omega$, and $i,j,l = 1, \dots, k$.

#### Condition 6
The order of differentiation with respect to $\Theta$ and integration can be exchanged when differentiating $\int f(x; \theta) dx$. That is:

$$
\frac{\partial}{\partial \Theta_i} \int f(x; \theta) dx = \int \frac{\partial f(x; \theta)}{\partial \Theta_i}  dx
$$

for $i = 1, \dots, k$. 

<details>
<summary>When will this hold?</summary>
Since differentiation is basically just a particular limit, we can use results about the interchanging of the integral and limit to get results about interchanging the integral with differentiation.

<div class="theorem">
  <strong>Dominated Convergence Theorem.</strong>
  {% tabs dom-conv %}
  {% tab dom-conv statement %}
  For a sequence of measurable functions $\left\{f_n\right\}_{n = 1}^\infty$ and measurable functions $f$ and $g$ satisfying $f_n(x) \rightarrow f(x)$ almost everywhere, $\rvert f_n(x) \rvert \leq g(x)$ almost everywhere, and $\int g(x) d\mu(x) < \infty$:

  $$
  \underset{n \rightarrow \infty}{\lim} \int f_n(x) d\mu(x) = \int f(x) d\mu(x)
  $$

  The dominated convergence theorem states that the integral of the limit of a sequence of measurable functions equals the limit of the integral of each element in the sequence.
  {% endtab %}
  {% tab dom-conv proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

The Dominated Convergence Theorem states that we can interchange the order of limits and integrals for (certain) functions that are always smaller than (in absolute value) some other function with finite integral. If we define a function that mimics the form of the derivative as a limit (something along the lines of $h(x) = \frac{f(x + \delta) - f(x)}{\delta}$), then we can use this theorem to get results for derivatives and integrals. This is the basic idea behind the <a href="https://en.wikipedia.org/wiki/Leibniz_integral_rule#">Leibniz integral rule</a>.

<div class="theorem">
  <strong>Leibniz Integral Rule.</strong>
  {% tabs leibniz %}
  {% tab leibniz statement %}
  Let $\Omega$ be an open subset of $\mathbb{R}$ and $\mathcal{M} = (S, \mathcal{A}, \mu)$ be a measure space. Let $f: \Omega \times \mathcal{M} \rightarrow \mathbb{R}$ be a function that satisfies:
  <ul>
  <li>$f(x; \Theta)$ is Lebesgue-integrable in $x$ for all $\theta \in \Omega$</li>
  <li>$\frac{\partial f(x; \theta)}{\partial \Theta}$ exists for all $\theta \in \Omega$ and for <i>almost all</i> $x \in \mathcal{M}$</li>
  <li>There exists integrable function $g: \mathcal{M} \rightarrow \mathbb{R}$ that is integrable and satisfies $\big\rvert \frac{\partial f(x; \theta)}{\partial \Theta} \big\rvert \leq g(x)$ for all $\theta \in \Omega$ and almost every $x \in \mathcal{M}$</li>
  </ul>
  Then, for all $\theta \in \Omega$:
  $$
  \frac{\partial}{\partial \Theta} \int_\mathcal{M} f(x; \theta) dx = \int_{\mathcal{M}} \frac{\partial}{\partial \Theta} f(x; \theta) dx
  $$
  {% endtab %}
  {% tab leibniz proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>
In summary, if our log-likelihood satisfies the (Lebesgue version of the) Leibniz Rule conditions, then it will satisfy Condition 2. 
</details>

#### Condition 6a
The order of differentiation with respect to $\Theta$ and integration can be exchanged when taking the second-order partial derivatives of $\int f(x; \theta) dx$.  That is:

$$
\frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \int f(x; \theta) dx = \int \frac{\partial^2 f(x; \theta)}{\partial \Theta_i \partial \Theta_j}  dx
$$

for $i,j = 1, \dots, k$. 

#### Condition 6b
Torder of differentiation with respect to $\theta$ and integration can be exchanged when taking the second-order partial derivatives of $\int f(x; \theta) dx$. That is:

$$
\frac{\partial^3}{\partial \Theta_i \partial \Theta_j \partial \Theta_l} \int f(x; \theta) dx = \int \frac{\partial^3 f(x; \theta)}{\partial \Theta_i \partial \Theta_j \partial \Theta_l}  dx
$$

for $i,j,l = 1, \dots, k$. 

#### Condition 7
There exist $M_{i,j,l}(x)$ and $c(\theta^*) > 0$ such that:

$$
\begin{aligned}
&\left\rvert \frac{\partial^3 \ell(\theta; x)}{\partial \Theta_i \partial \Theta_j \partial \Theta_l} \right\rvert \leq M_{i,j,l}(x) \\
\text{for all } &\theta \in \left\{ \theta \in \Omega \hspace{1mm} : \hspace{1mm} \rvert \rvert \theta_i - \theta_0 \rvert \rvert_2^2 < c(\theta^*) \right\}
\end{aligned}
$$

with $$\mathbb{E}_{\theta^*}\left[ M_{i,j,l}(X) \right] < \infty$$ for all $i,j,l = 1, \dots, k$.

#### Condition 8
For all $i,j = 1, \dots, k$, $\mathcal{I}_{i,j}(\theta) < \infty$, and $\mathcal{I}(\theta)$ is positive definite.

Note: Lehmann and Casella only require $\mathcal{I}(\theta)$ to be positive definite for all $\theta \in \mathbf{N}$ (pg. 463).<d-cite key=lehmann2005></d-cite> 


---

## Results
We now can state several results that arise when these conditions hold. The following tells us when the score has expectation zero and when the information matrix can be found as the negative expected Hessian.

<div class="theorem">
<strong>Theorem 7.5.1. (Lehmann, pg. 500)</strong>
{% tabs thrm-7-5-1 %}
{% tab thrm-7-5-1 statement %}
Assume Conditions 1-5, and 6 hold. Then:

$$
\begin{equation}
\label{eq:7-5-1a}
\begin{aligned}
\mathbb{E}_{\theta^*}\left[ U_{\Theta}(\theta^*; X) \right] &= \mathbf{0}_k \\
\mathcal{I}_{i,j}(\theta^*) &= \mathbb{E}_{\theta^*} \left[ U_{\Theta_i}(\theta^*; X) U_{\Theta_j}(\theta^*; X) \right]
\end{aligned}
\end{equation}
$$

If 5a and 6a also hold, then:

$$
\begin{equation}
\label{eq:7-5-1b}
\begin{aligned}
\mathbb{E}_{\theta^*}\left[ U_{\Theta}(\theta^*; X) \right] &= \mathbf{0}_k \\
\mathcal{I}_{i,j}(\theta^*) &= - \mathbb{E}_{\theta^*}\left[ \frac{\partial^2 \log(f(X; \theta))}{\partial \Theta_i \partial \Theta_j} \right]\\
\end{aligned}
\end{equation}
$$

and $\mathcal{I}(\theta)$ is positive semi-definite.
{% endtab %}
{% tab thrm-7-5-1 proof %}
To show the first claim in Eqs. \eqref{eq:7-5-1a} and \eqref{eq:7-5-1b}; 

$$
\begin{aligned}
\mathbb{E}_{\theta^*}\left[ U_{\Theta_i}(x; \theta) \right]
&= \int f(x; \theta^*) \frac{\partial \log f(x; \theta^*)}{\partial \Theta_i} d x \\
&= \int \frac{\partial f(x; \theta^*)}{\partial \Theta_i} d x \\
&= \frac{\partial}{\partial \Theta_i} \left[ \int f(x; \theta^*) d x \right]\\
&= \frac{\partial}{\partial \Theta_i} [1]\\
&= 0 
\end{aligned}
$$

The second claim in Eq. \eqref{eq:7-5-1a} follows from the definition of the Fisher information matrix and the fact that the score function has mean zero (shown above). To show the second claim in Eq. \eqref{eq:7-5-1b}, note that:

$$
\begin{aligned}
&\int f(x; \theta) dx = 1 \\
\implies &\frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[ \int f(x; \theta) dx \right] = 0 
\end{aligned}
$$

By Conditions <a href="#condition-5a">5a</a> and <a href="#condition-6a">6a</a>, we have:

$$
\begin{aligned}
&\frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[ \int f(x; \theta^*) dx \right] = \int \frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[  f(x; \theta^*) \right] dx  \\
\implies &\int \frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[  f(x; \theta^*) \right] dx = 0 \\
\implies &\int \frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[  f(x; \theta^*) \right] \frac{f(x; \theta^*)}{f(x; \theta^*)} dx = 0 \\
\implies &\mathbb{E}_{\theta^*} \left[ \frac{\frac{\partial^2 f(X; \theta^*) }{\partial \Theta_i \partial \Theta_j}}{f(X; \theta^*)}\right] = 0
\end{aligned}
$$

By the chain and quotient rules, we have that:

$$
\begin{aligned}
\frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[  \log(f(x; \theta)) \right] 
&= \frac{\partial}{\partial \Theta_i} \left[ \frac{\frac{\partial f(x; \theta)}{\partial \Theta_i}}{f(x; \theta)}\right] \\
&= \frac{\frac{\partial^2 f(x; \theta)}{\partial \Theta_i \partial \Theta_j}f(x; \theta) - \frac{\partial f(x; \theta)}{\partial \Theta_j}\frac{\partial f(x; \theta)}{\partial \Theta_i}}{(f(x; \theta))^2} \\
&= \frac{\frac{\partial^2 f(x; \theta)}{\partial \Theta_i \partial \Theta_j}}{f(x; \theta)} - \frac{\frac{\partial f(x; \theta)}{\partial \Theta_j}\frac{\partial f(x; \theta)}{\partial \Theta_i}}{(f(x; \theta))^2} 
\end{aligned}
$$

Thus:

$$
\begin{aligned}
&\int \frac{\partial^2}{\partial \Theta_i \partial \Theta_j} \left[  \log(f(x; \theta^*)) \right] f(x; \theta^*) dx = \int\left[\frac{\frac{\partial^2 f(x; \theta^*)}{\partial \Theta_i \partial \Theta_j}}{f(x; \theta^*)} - \frac{\frac{\partial f(x; \theta^*)}{\partial \Theta_j}\frac{\partial f(x; \theta^*)}{\partial \Theta_i}}{(f(x; \theta^*))^2}  \right] f(x; \theta^*) dx\\
\implies &\mathbb{E}_{\theta^*} \left[ \frac{\partial^2 \log(f(x; \theta^*))}{\partial \Theta_i \partial \Theta_j}  \right] = \mathbb{E}_{\theta^*} \left[ \frac{\frac{\partial^2 f(x; \theta^*)}{\partial \Theta_i \partial \Theta_j}}{f(x; \theta^*)} - \frac{\frac{\partial f(x; \theta^*)}{\partial \Theta_j}\frac{\partial f(x; \theta^*)}{\partial \Theta_i}}{(f(x; \theta^*))^2}  \right] \\
\implies &\mathbb{E}_{\theta^*} \left[ \frac{\partial^2 \log(f(x; \theta^*))}{\partial \Theta_i \partial \Theta_j} \right] = 0 - \mathbb{E}\left[ \frac{\frac{\partial f(x; \theta^*)}{\partial \Theta_j}}{f(x; \theta^*)} \frac{\frac{\partial f(x; \theta^*)}{\partial \Theta_i}}{f(x; \theta^*)} \right] \\
\implies &\mathbb{E}_{\theta^*} \left[ \frac{\partial^2 \log(f(x; \theta^*))}{\partial \Theta_i \partial \Theta_j} \right] = - \mathbb{E}\left[ \frac{\partial \log( f(x; \theta^*))}{\partial \Theta_j} \frac{\partial \log(f(x; \theta^*))}{\partial \Theta_i} \right] \\
\implies &-\mathbb{E}_{\theta^*} \left[ \frac{\partial^2 \log(f(x; \theta^*))}{\partial \Theta_i \partial \Theta_j} \right] = \text{Cov}\left(U_{\Theta_i}(\theta^*; X), U_{\Theta_j}(\theta^*; X) \right)
\end{aligned}
$$

where the last line follows from the fact that the score has mean zero (shown above). Thus:

$$
\mathcal{I}_{i,j}(\theta^*) = -\mathbb{E}_{\theta^*} \left[ \frac{\partial^2 \log(f(x; \theta^*))}{\partial \Theta_i \partial \Theta_j} \right]
$$

By the above proof, we have shown that $\mathcal{I}(\theta^*)$ is the variance-covariance matrix of $U_{\Theta}(\theta^*; X)$. Let $\mathbf{v} = (v_1, \dots, v_k)^\top$ be a constant (non-zero) vector. We then see that:

$$
\begin{aligned}
\mathbf{v}^\top \mathcal{I}(\theta^*) \mathbf{v}
    &= \sum_{i = 1}^k \sum_{j = 1}^k \mathcal{I}_{i,j}(\theta^*) v_i v_j \\
    &= \sum_{i = 1}^k \mathcal{I}_{i,i}(\theta^*) v_i^2 + \sum_{i = 1}^k \sum_{i' \neq i} \mathcal{I}_{i,i'}(\theta^*) v_i v_{i'} \\
    &= \sum_{i = 1}^k \text{Var}(v_i U_{\Theta_i}(\theta^*; X)) + \sum_{i = 1}^k \sum_{i' \neq i} v_i v_{i'} \text{Cov}(U_{\Theta_i}(\theta^*; X), U_{\Theta_j}(\theta^*; X)) \\
    &= \text{Var}\left( \sum_{i = 1}^k v_i U_{\Theta_i}(\theta^*; X) \right) \\
    &\geq 0
\end{aligned}
$$

where the last equality follows from the equation for the variance of a linear combination.
{% endtab %}
{% endtabs %}
</div>

If we make a few more assumptions, then we can show that there exists a solution to $U_{\theta}(\theta; \mathbf{x}) = \mathbf{0}_k$ that is asymptotically normal and consistent.

<div class="theorem">
  <strong>Theorem 7.5.2.</strong>
  {% tabs thrm-7-5-2 %}
  {% tab thrm-7-5-2 statement %}
  Assume that Conditions 1-8 hold. Then there exists a solution, $\hat{\theta}$, to $U_{\Theta}(\hat{\theta}; \mathbf{x}) = \mathbf{0}_k$ that is consistent. Furthermore, any solution will satisfy:

  $$
  \sqrt{n} (\hat{\theta} - \theta^*) \rightsquigarrow \mathcal{N}\left(\mathbf{0}_k, \mathcal{I}^{-1}(\theta^*)\right)
  $$

  (Technically, the asymptotic variance holds everywhere except for a set of values of $\theta$ with measure $0$.)
  {% endtab %}
  {% tab thrm-7-5-2 proof %}
    See pg. 463-465 of Lehmann and Casella for a proof.<d-cite key=lehmann2005></d-cite>
  {% endtab %}
  {% endtabs %}
</div>

<aside><p>Note: $\mathcal{I}^{-1}(\theta^*)$ is the <i>variance of the asymptotic distribution</i>, which is not necessarily equal to the limit of the variance! It is also the Fisher information of a single arbitrary $X$ since they are i.i.d. (Schervish uses the notation $\mathcal{I}^{-1}_{X_1}(\theta^*)$ to emphasize this).</p></aside>


We can also define the <i>efficiency</i> of an estimator. 

<div class="definition">
<strong>Definition (Efficiency).</strong>
<br>
Let $\boldsymbol{\delta}_1$ and $\boldsymbol{\delta}_2$ be two estimators of $\theta$ such that:

$$
\begin{aligned}
\sqrt{n}(\boldsymbol{\delta}_1 - \theta) &\rightsquigarrow \mathcal{N}\left(\mathbf{0}_k, \boldsymbol{\Sigma}_1(\theta) \right) \\
\sqrt{n}(\boldsymbol{\delta}_2 - \theta) &\rightsquigarrow \mathcal{N}\left(\mathbf{0}_k, \boldsymbol{\Sigma}_2(\theta) \right) 
\end{aligned}
$$

We call $\boldsymbol{\delta}_1$ <i>more efficient</i> than $\boldsymbol{\delta}_2$ if: 

$$\boldsymbol{\Sigma}_2(\theta) - \boldsymbol{\Sigma}_1(\theta)$$

is positive semi-definite for all $\theta \in \Omega$. 
</div>

Similarly, we call $\hat{\theta}$ an <strong>efficient estimator</strong> of $\theta$ if it is asymptotically normal with mean zero and non-singular covariance matrix $\boldsymbol{\Sigma}(\theta) = \mathcal{I}^{-1}(\theta)$. An example is the estimator discussed in Theorem 7.5.2. (<a href="https://en.wikipedia.org/wiki/Cramér–Rao_bound">Cramér-Rao Bound</a>!). 

<aside><p>Only if we also assume that $\boldsymbol{\Sigma}_{i,j}(\theta)$ and $\mathcal{I}_{i,j}^{-1}(\theta)$ are continuous functions of $\theta$ for all $i,j$.</p></aside>

The results above can easily be extended to cases where we have multiple independent samples or if we have nuisance parameters. In the latter case, the efficiency will generally decrease with the number of nuisance parameters unless the parameter of interest is independent of them.
