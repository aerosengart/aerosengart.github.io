---
layout: distill
title: Maximum Likelihood Estimation (Standard)
description: A Primer
date: 2026-07-27
tabs: true
tags: theory likelihood primer
toc:  
  - name: Background
    subsections:
        - name: Notation and Set-Up
        - name: Definitions
        - name: Properties
  - name: Regularity Conditions
  - name: Properties
  - name: Estimation
bibliography: stats-ml.bib
---

This post mostly covers Chapter 7 of <i>Elements of Large-Sample Theory</i> by Erich Lehmann<d-cite key=lehmann2004></d-cite> with some supplemental material from <i>Theory of Statistics</i> by Mark Schervish,<d-cite key=schervish1995></d-cite> and <i>Theory of Point Estimation</i> by Lehmann and Casella.<d-cite key=lehmann2005></d-cite>

---

## Background 

### Notation and Set-Up
Suppose we have some probability space $(S, \mathcal{A}, \mu)$. A random variable is some function $X: S \rightarrow \mathcal{X}$ where $\mathcal{X}$ is the sample space (which we only require to be a <a href="/measure-theory/chapter2/#borel-set">Borel space</a> but is usually a subset of Euclidean space) with <a href="/measure-theory/chapter2/#sigma-field">$\sigma$-field</a> $\mathcal{B}$. Individual elements of $\mathcal{X}$ are denoted with $x$.

Let $$\mathcal{P}_\boldsymbol{\theta}$$ be a family of distributions for $X$ parametrized by $\boldsymbol{\theta}: S \rightarrow \boldsymbol{\Theta}$ where $\boldsymbol{\Theta}$ is the parameter space with $\sigma$-field $\tau$. We allow $\boldsymbol{\Theta}$ to be scalar- or vector-valued. Here, we assume it is $p$-dimensional. Individual elements of $\boldsymbol{\Theta}$ are denoted with $\boldsymbol{\theta}$. Denote the conditional distribution of $X$ given $$\boldsymbol{\theta} = \boldsymbol{\theta}^*$$ with $P_\boldsymbol{\theta}$ (which is a distribution on $(\mathcal{X},\mathcal{B})$). 

To make our notation match less measure theoretical texts, we'll use $X$ to denote a random variable with realizations denoted with the lowercase $x$. A parameter will be denoted with $\boldsymbol{\theta}$ with particular values denoted by superscripts and subscripts and its true value by $$\boldsymbol{\theta}^*$$. The probability density or mass function of $X$ given parameter $\boldsymbol{\theta}$ evaluated at a particular $x$ and $$\boldsymbol{\theta}_0$$ will be denoted by $$f_{X \rvert \boldsymbol{\theta}}(x; \boldsymbol{\theta}_0)$$ or, more compactly, $$f(x; \boldsymbol{\theta}_0)$$. We will also use the notation $\mathbf{X} = (X_1, \dots, X_n)$ and $\mathbf{x} = (x_1, \dots, x_n)$ to denote a collection of $n$ random variables or vectors $X$ and observations of each, respectively. Unless otherwise stated, we assume $X_1, \dots, X_n$ are i.i.d.

### Definitions
First, we need to lay out some definitions and background information for our discussion. The <i>score function</i> will be a star player in later sections.

<div class="definition">
<strong>Definition (Score; Schervish, pg. 111).</strong>
<br>
The <i>score function</i> (or <i>score statistic</i> or <i>efficient score</i>) is the $k$-dimensional vector with $i$-th coordinate given by:

$$
\begin{equation}
\label{eq:score}
    U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}_0; x) = \left. \frac{\partial \log f(x; \boldsymbol{\theta}_0)}{\partial \boldsymbol{\theta}_i} \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0}
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
\mathcal{I}_{i,j}(\boldsymbol{\theta}; X) 
    = \text{Cov}_{\boldsymbol{\theta}}\left( \frac{\partial \log f(X; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}, \frac{\partial \log f(X; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_j} \right)
    = \text{Cov}_{\boldsymbol{\theta}}\left( U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}; X), U_{\boldsymbol{\theta}_j}(\boldsymbol{\theta}; X) \right)
\end{equation}
$$

where the expectation is with respect to $X \sim P_{\boldsymbol{\theta}}$.
</div>

Both definitions can be extended to the case where we condition on a statistic. These are called the <i>conditional score function</i> and <i>conditional Fisher information</i>, and are denoted by $U_{\boldsymbol{\theta} \rvert T}(\boldsymbol{\theta}; x)$ and $\mathcal{I}_{\boldsymbol{\theta}; X \rvert T}(\boldsymbol{\theta})$ for statistic $T$, respectively (see pg. 111 and Section 2.3.3 in Schervish).

<div class="definition">
<strong>Definition (Likelihood).</strong>
The <i>likelihood</i> of $\boldsymbol{\theta}$ given $\mathbf{X} = \mathbf{x}$ is the probabilityt mass/density function $f(\mathbf{x}; \boldsymbol{\theta})$ viewed as a function of $\boldsymbol{\theta}$:

$$
\mathcal{L}(\boldsymbol{\theta}; \mathbf{x}) = f_{\mathbf{X} \rvert \boldsymbol{\theta}}(\mathbf{x}; \boldsymbol{\theta}) = \prod_{i = 1}^n f_{X_i \rvert \boldsymbol{\theta}}(x_i; \boldsymbol{\theta})
$$

where the last equality follows from the independence assumption. 
</div>

<aside><p>Note: the likelihood is also a function of the data!</p></aside>

We also often work with the natural logarithm of $\mathcal{L}(\boldsymbol{\theta}; \mathbf{x})$, denoted by $\ell(\boldsymbol{\theta}; \mathbf{x})$, and called the <strong>log-likelihood</strong>.

In later sections, we will rely on assumptions on the uniform convergence and consistency of estimators.

<div class="definition">
<strong>Definition (Uniform Convergence in Probability).<d-cite key=newey1994></d-cite></strong>
Let $Y^\gamma$ be a random variable, let $\{ Y_n^\gamma \}_{n \in \mathbb{N}}$ be a sequence of random variables, and let $\gamma \in \boldsymbol{\Gamma}$. We say that the sequence <i>converges uniformly in probability</i> to $Y^\gamma$ if, for any $\epsilon > 0$:

$$
\underset{n \rightarrow \infty}{\lim} \left\{ \mathbb{P}\left( \underset{\gamma \in \boldsymbol{\Gamma}}{\sup} \left\{ Y_n^\gamma - Y^\gamma \right\} > \epsilon \right) \right\} = 0
$$
</div>

If $\{Y^\gamma_n\}_{n \in \mathbb{N}}$ is a sequence of estimators, then uniform convergence in probability to $Y^\gamma$ means that, in the limit, our estimators will be arbitrarily "good".

By <i>uniformly consistent</i>, I believe we mean that a sequence of estimators converges uniformly in probability to the true parameter value.


---

## Regularity Conditions
Below are the conditions needed for stating some of the results discussed later in this post. The statements come from Lehmann (2005) and Lehmann and Casella (2005).

<aside><p>These are often call the Cramér conditions.</p></aside>

Note that not all of the conditions are needed for each result,and some can be relaxed slightly. For example, Conditions 6 and 6A could be exchanged for requiring Eq. \eqref{eq:5-3-b} to hold.

#### Condition 1
If $$P_{\boldsymbol{\theta}_1} = P_{\boldsymbol{\theta}_2}$$ for $\boldsymbol{\theta}_1, \boldsymbol{\theta}_2 \in \boldsymbol{\Theta}$, then $\boldsymbol{\theta}_1 = \boldsymbol{\theta}_2$.

In words, this implies that the parameter is identifiable and that the distributions in the class $\mathcal{P}_{\boldsymbol{\theta}}$ are distinct.

#### Condition 2
For some $\delta > 0$, there exists a set, $\boldsymbol{\omega}$, defined as:

$$
\boldsymbol{\omega} = \left\{\boldsymbol{\theta} \in \mathbb{R}^k \hspace{1mm} : \hspace{2mm} \rvert \rvert \boldsymbol{\theta} - \boldsymbol{\theta}^* \rvert \rvert_2 < \delta \right\}
$$

such that $\boldsymbol{\omega} \subseteq \boldsymbol{\Theta}$. 

<aside><p>$\boldsymbol{\omega}$ is called a <strong>neighborhood</strong>.</p></aside>

This condition can be replaced by the stronger one (but easier to check) that $\boldsymbol{\Theta}$ is open. In other words, we check that it holds for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ since we do not know what $\boldsymbol{\theta}^*$ is. 

#### Condition 3
We observe $X_1, \dots, X_n$, which are i.i.d. and have probability density or mass function $f_{X \rvert \boldsymbol{\theta}}(x; \boldsymbol{\theta}) = f(x; \boldsymbol{\theta})$ with respect to measure $\mu$.

#### Condition 4
The support of $f(x; \boldsymbol{\theta})$, defined as:

$$
\mathbf{A} = \left\{ x \hspace{1mm} : \hspace{1mm} f(x; \boldsymbol{\theta}) > 0 \right\}
$$

is the same for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ (i.e. independent of $\boldsymbol{\theta}$).

#### Condition 5
The first-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to $\Theta$, denoted by $\frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}$, exist for all $x \in \mathbf{A}$ and $i = 1, \dots, k$. 

Note: Lehmann and Casella only require that the partial derivatives (in 5, 5a, and 5b) exist for <i>almost all</i> $x \in \mathcal{X}$ and only all $\boldsymbol{\theta} \in \boldsymbol{\omega}$ (pg. 462).<d-cite key=lehmann2005></d-cite>

#### Condition 5a
The second-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to $\boldsymbol{\theta}$, denoted by $\frac{\partial^2 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}$, exist for all $x \in \mathbf{A}$, $\boldsymbol{\theta} \in \boldsymbol{\Theta}$, and $i,j = 1, \dots, k$.

#### Condition 5b
The third-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to $\boldsymbol{\theta}$, denoted by $\frac{\partial^3 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_l}$, exist for all $x \in \mathbf{A}$, $\boldsymbol{\theta} \in \boldsymbol{\Theta}$, and $i,j,l = 1, \dots, k$.

#### Condition 6
The order of differentiation with respect to $\Theta$ and integration can be exchanged when differentiating $\int f(x; \boldsymbol{\theta}) dx$. That is:

$$
\frac{\partial}{\partial \boldsymbol{\theta}_i} \int f(x; \boldsymbol{\theta}) dx = \int \frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}  dx
$$

for $i = 1, \dots, k$. 

<details>
<summary>When will this hold?</summary>
Since differentiation is basically just a particular limit, we can use results about the interchanging of the integral and limit to get results about interchanging the integral with differentiation.

<div class="theorem">
  <strong>Dominated Convergence Theorem.</strong>
  {% tabs dom-conv %}
  {% tab dom-conv statement %}
  For a sequence of measurable functions $$\left\{f_n\right\}_{n = 1}^\infty$$ and measurable functions $f$ and $g$ satisfying $f_n(x) \rightarrow f(x)$ almost everywhere, $\rvert f_n(x) \rvert \leq g(x)$ almost everywhere, and $\int g(x) d\mu(x) < \infty$:

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
  Let $\boldsymbol{\Theta}$ be an open subset of $\mathbb{R}$ and $\mathcal{M} = (S, \mathcal{A}, \mu)$ be a measure space. Let $f: \boldsymbol{\Theta} \times \mathcal{M} \rightarrow \mathbb{R}$ be a function that satisfies:
  <ul>
  <li>$f(x; \Theta)$ is Lebesgue-integrable in $x$ for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$</li>
  <li>$\frac{\partial f(x; \boldsymbol{\theta})}{\partial \Theta}$ exists for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ and for <i>almost all</i> $x \in \mathcal{M}$</li>
  <li>There exists integrable function $g: \mathcal{M} \rightarrow \mathbb{R}$ that is integrable and satisfies $\big\rvert \frac{\partial f(x; \boldsymbol{\theta})}{\partial \Theta} \big\rvert \leq g(x)$ for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ and almost every $x \in \mathcal{M}$</li>
  </ul>
  Then, for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$:
  $$
  \frac{\partial}{\partial \Theta} \int_\mathcal{M} f(x; \boldsymbol{\theta}) dx = \int_{\mathcal{M}} \frac{\partial}{\partial \Theta} f(x; \boldsymbol{\theta}) dx
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
The order of differentiation with respect to $\boldsymbol{\theta}$ and integration can be exchanged when taking the second-order partial derivatives of $\int f(x; \boldsymbol{\theta}) dx$.  That is:

$$
\frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \int f(x; \boldsymbol{\theta}) dx = \int \frac{\partial^2 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}  dx
$$

for $i,j = 1, \dots, k$. 

#### Condition 6b
Torder of differentiation with respect to $\boldsymbol{\theta}$ and integration can be exchanged when taking the second-order partial derivatives of $\int f(x; \boldsymbol{\theta}) dx$. That is:

$$
\frac{\partial^3}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_l} \int f(x; \boldsymbol{\theta}) dx = \int \frac{\partial^3 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_l}  dx
$$

for $i,j,l = 1, \dots, k$. 

#### Condition 7
There exist $M_{i,j,l}(x)$ and $c(\boldsymbol{\theta}^*) > 0$ such that:

$$
\begin{aligned}
&\left\rvert \frac{\partial^3 \ell(\boldsymbol{\theta}; x)}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_l} \right\rvert \leq M_{i,j,l}(x) \\
\text{for all } &\boldsymbol{\theta} \in \boldsymbol{\omega}
\end{aligned}
$$

with $$\mathbb{E}_{\boldsymbol{\theta}^*}\left[ M_{i,j,l}(X) \right] < \infty$$ for all $i,j,l = 1, \dots, k$.

#### Condition 8
For all $i,j = 1, \dots, k$, $\mathcal{I}_{i,j}(\boldsymbol{\theta}; X) < \infty$, and $\mathcal{I}(\boldsymbol{\theta}; X)$ is positive definite.

Note: Lehmann and Casella only require $\mathcal{I}(\boldsymbol{\theta}; X)$ to be positive definite for all $\boldsymbol{\theta} \in \boldsymbol{\omega}$ (pg. 463).<d-cite key=lehmann2005></d-cite> 

---

## Properties
Before we dive into maximum likelihood estimation, we can show a few properties of the score and information that will be helpful for later proofs.

<div id="lemma-5-3-l"></div>
<div class="theorem">
<strong>Lemma 5.3 (Lehmann and Casella, pg. 116)</strong>
{% tabs thrm-7-5-1 %}
{% tab thrm-7-5-1 statement %}
Suppose $\boldsymbol{\Theta}$ is open, and Conditions 3, 4, 5, 6 hold. Then:

$$
\begin{equation}
\label{eq:5-3-a}
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}^*}\left[ U_{\boldsymbol{\theta}}(\boldsymbol{\theta}^*; X) \right] &= \mathbf{0}_p \\
\mathcal{I}_{i,j}(\boldsymbol{\theta}^*; X) &= \mathbb{E}_{\boldsymbol{\theta}^*} \left[ U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; X) U_{\theta_j}(\boldsymbol{\theta}^*; X) \right]
\end{aligned}
\end{equation}
$$

Suppose Conditions 5A and 6A also hold. Then:

$$
\begin{equation}
\label{eq:5-3-b}
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}^*}\left[ U_{\theta}(\boldsymbol{\theta}^*; X) \right] &= \mathbf{0}_k \\
\mathcal{I}_{i,j}(\boldsymbol{\theta}^*; X) &= - \mathbb{E}_{\boldsymbol{\theta}^*}\left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}))}{\partial \theta_i \partial \theta_j} \right]\\
\end{aligned}
\end{equation}
$$

and $\mathcal{I}(\boldsymbol{\theta}; X)$ is positive semi-definite.
{% endtab %}
{% tab thrm-7-5-1 proof %}
To show the first claim in Eqs. \eqref{eq:5-3-a} and \eqref{eq:5-3-b}; 

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}^*}\left[ U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}; X) \right]
&= \int f(X; \boldsymbol{\theta}^*) \frac{\partial \log f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i} d X \\
&= \int \frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i} d X \\
&= \frac{\partial}{\partial \boldsymbol{\theta}_i} \left[ \int f(X; \boldsymbol{\theta}^*) d X \right]\\
&= \frac{\partial}{\partial \boldsymbol{\theta}_i} [1]\\
&= 0 
\end{aligned}
$$

The second claim in Eq. \eqref{eq:5-3-a} follows from the definition of the Fisher information matrix and the fact that the score function has mean zero (shown above). To show the second claim in Eq. \eqref{eq:5-3-b}, note that:

$$
\begin{aligned}
&\int f(X; \boldsymbol{\theta}) dX = 1 
\implies \frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[ \int f(X; \boldsymbol{\theta}) d \right] = 0 
\end{aligned}
$$

Then by Conditions 5A and 6A:

$$
\begin{aligned}
&\frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[ \int f(X; \boldsymbol{\theta}^*) dx \right] = \int \frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[  f(X; \boldsymbol{\theta}^*) \right] dX \\
\implies &\int \frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[  f(X; \boldsymbol{\theta}^*) \right] dX = 0 \\
\implies &\int \frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[  f(X; \boldsymbol{\theta}^*) \right] \frac{f(x; \boldsymbol{\theta}^*)}{f(X; \boldsymbol{\theta}^*)} dX = 0 \\
\implies &\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\frac{\partial^2 f(X; \boldsymbol{\theta}^*) }{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}}{f(X; \boldsymbol{\theta}^*)}\right] = 0
\end{aligned}
$$

By the chain and quotient rules, we have that:

$$
\begin{aligned}
\frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[  \log(f(x; \boldsymbol{\theta})) \right] 
&= \frac{\partial}{\partial \boldsymbol{\theta}_i} \left[ \frac{\frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}}{f(x; \boldsymbol{\theta})}\right] \\
&= \frac{\frac{\partial^2 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}f(x; \boldsymbol{\theta}) - \frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_j}\frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}}{(f(x; \boldsymbol{\theta}))^2} \\
&= \frac{\frac{\partial^2 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}}{f(x; \boldsymbol{\theta})} - \frac{\frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_j}\frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}}{(f(x; \boldsymbol{\theta}))^2} 
\end{aligned}
$$

Thus:

$$
\begin{aligned}
&\int \frac{\partial^2}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \left[  \log(f(X; \boldsymbol{\theta}^*)) \right] f(X; \boldsymbol{\theta}^*) dX = \int\left[\frac{\frac{\partial^2 f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}}{f(X; \boldsymbol{\theta}^*)} - \frac{\frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_j}\frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i}}{(f(X; \boldsymbol{\theta}^*))^2}  \right] f(X; \boldsymbol{\theta}^*) dx\\
\implies &\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}  \right] = \mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\frac{\partial^2 f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}}{f(X; \boldsymbol{\theta}^*)} - \frac{\frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_j}\frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i}}{(f(X; \boldsymbol{\theta}^*))^2}  \right] \\
\implies &\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right] = 0 - \mathbb{E}\left[ \frac{\frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_j}}{f(X; \boldsymbol{\theta}^*)} \frac{\frac{\partial f(X; \boldsymbol{\theta}^*)}{\partial \boldsymbol{\theta}_i}}{f(X; \boldsymbol{\theta}^*)} \right] \\
\implies &\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right] = - \mathbb{E}\left[ \frac{\partial \log( f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_j} \frac{\partial \log(f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_i} \right] \\
\implies &-\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right] = \text{Cov}\left(U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; X), U_{\boldsymbol{\theta}_j}(\boldsymbol{\theta}^*; X) \right)
\end{aligned}
$$

where the last line follows from the fact that the score has mean zero (shown above). Thus:

$$
\mathcal{I}_{i,j}(\boldsymbol{\theta}^*; X) = -\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}^*))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right]
$$

By the above proof, we have shown that $$\mathcal{I}(\boldsymbol{\theta}^*; X)$$ is the variance-covariance matrix of $$U_{\boldsymbol{\theta}}(\boldsymbol{\theta}^*; X)$$. Let $$\mathbf{v} = (v_1, \dots, v_k)^\top$$ be a constant (non-zero) vector. We then see that:

$$
\begin{aligned}
\mathbf{v}^\top \mathcal{I}(\boldsymbol{\theta}^*; X) \mathbf{v}
    &= \sum_{i = 1}^k \sum_{j = 1}^k \mathcal{I}_{i,j}(\boldsymbol{\theta}^*) v_i v_j \\
    &= \sum_{i = 1}^k \mathcal{I}_{i,i}(\boldsymbol{\theta}^*; X) v_i^2 + \sum_{i = 1}^k \sum_{i' \neq i} \mathcal{I}_{i,i'}(\boldsymbol{\theta}^*) v_i v_{i'} \\
    &= \sum_{i = 1}^k \text{Var}(v_i U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; X)) + \sum_{i = 1}^k \sum_{i' \neq i} v_i v_{i'} \text{Cov}(U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; X), U_{\boldsymbol{\theta}_j}(\boldsymbol{\theta}^*; X)) \\
    &= \text{Var}\left( \sum_{i = 1}^k v_i U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; X) \right) \\
    &\geq 0
\end{aligned}
$$

where the last equality follows from the equation for the variance of a linear combination.
{% endtab %}
{% endtabs %}
</div>

We can also show that the information is additive for independent variables.

<div class="theorem">
<strong>Theorem 5.8. (Lehmann and Casella, pg. 119)</strong>
{% tabs thrm-5-8-l %}
{% tab thrm-5-8-l statement %}
Let $X$ and $Y$ be independent with respective densities $p_{\boldsymbol{\theta}}$ and $q_{\boldsymbol{\theta}}$ with respect to measures $\mu$ and $\nu$. Suppose Conditions 1-4 in Lemma 5.3 hold so that the score has mean zero. 
Let $\mathcal{I}_X(\boldsymbol{\theta})$, $\mathcal{I}_Y(\boldsymbol{\theta})$, and $\mathcal{I}(\boldsymbol{\theta})$ denote the information for $X$, $Y$, and $(X,Y)$, respectively. Then:

$$
\mathcal{I}(\boldsymbol{\theta}) = \mathcal{I}_X(\boldsymbol{\theta}) + \mathcal{I}_Y(\boldsymbol{\theta})
$$
{% endtab %}
{% tab thrm-5-8-l proof %}
$$
\begin{aligned}
\mathcal{I}(\boldsymbol{\theta})
    &= \text{Var}_{\boldsymbol{\theta}}\left( \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}} + \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}} \right) \\
    &= \mathbb{E}\left[ \left(  \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}} + \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}}  \right) \left(  \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}} + \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}}  \right)^\top \right] \\
    &= \mathbb{E}\left[ \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}} \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}^\top} + 2 \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}} \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}^\top} + \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}}\frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}^\top}\right] \\
    &= \text{Var}\left( \frac{\partial \log p_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}} \right) + \text{Var}\left( \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}} \right) + 2 \underbrace{\mathbb{E}\left[ \frac{\partial \log q_{\boldsymbol{\theta}}(X)}{\partial \boldsymbol{\theta}}  \right]}_{=\mathbf{0}_p} \underbrace{\mathbb{E}\left[ \frac{\partial \log q_{\boldsymbol{\theta}}(Y)}{\partial \boldsymbol{\theta}^\top} \right]}_{=\mathbf{0}_p} \\
    &= \mathcal{I}_X(\boldsymbol{\theta}) + \mathcal{I}_Y(\boldsymbol{\theta})
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

<aside><p>Applying Theorem 5.8 to i.i.d. $X_1, \dots, X_n$ gives us the result that the information in $\mathbf{X} = (X_1, \dots, X_n)$ about $\boldsymbol{\theta}$ is $n \mathcal{I}_X(\boldsymbol{\theta}; X)$.</p></aside>

We can also show that the score function evaluated at the true parameter value is asymptotically normal.

<div class="theorem">
<strong>Lemma (shown in Lehmann, pg. 470)</strong>
{% tabs lemma-470 %}
{% tab lemma-470 statement %}
Suppose the conditions of Lemma 5.3 are satisfied. Then:

$$
\frac{1}{\sqrt{n}} U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; \mathbf{X}) \rightsquigarrow \mathcal{N}(0, \mathcal{I}_{i,i}(\boldsymbol{\theta}^*; \mathbf{X}))
$$
{% endtab %}
{% tab lemma-470 proof %}
$$
\begin{aligned}
\frac{1}{\sqrt{n}} U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; \mathbf{X})
  &= \sqrt{n} \frac{1}{n} \left. \frac{\partial}{\partial \boldsymbol{\theta}_i} \left[ \log f(\mathbf{X}; \boldsymbol{\theta}) \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}^*} \\
  &= \sqrt{n} \frac{1}{n} \left. \frac{\partial}{\partial \boldsymbol{\theta}_i} \left[ \log \left( \prod_{j = 1}^n f(X_j; \boldsymbol{\theta})  \right) \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}^*} \\
  &= \sqrt{n} \frac{1}{n} \sum_{j = 1}^n \left. \frac{\partial}{\partial \boldsymbol{\theta}_i} \left[ \log \left(  f(X_j; \boldsymbol{\theta})  \right) \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}^*} \\
  &= \sqrt{n} \underbrace{\frac{1}{n} \sum_{j = 1}^n U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; X_j)}_{(a)}\\
\end{aligned}
$$

$(a)$ is the average of i.i.d. random variables with mean zero and variance $\mathcal{I}_{i,i}(\boldsymbol{\theta}; X)$ by Eq. \eqref{eq:5-3-a}. Since $\text{Var}(c Y) = c^2 \text{Var}(Y)$ for random variable $Y$ and constant $c$, we get:

$$
\frac{1}{\sqrt{n}} U_{\boldsymbol{\theta}_i}(\boldsymbol{\theta}^*; \mathbf{X}) \rightsquigarrow \mathcal{N}(0, n \mathcal{I}_{i,i}(\boldsymbol{\theta}^*; X))
$$

via the central limit theorem.
{% endtab %}
{% endtabs %}
</div>


<div class="theorem">
<strong>Theorem 6.6. (Lehmann and Casella, pg. 127)</strong>
{% tabs thrm-6-6-l %}
{% tab thrm-6-6-l statement %}
TODO
{% endtab %}
{% tab thrm-6-6-l proof %}
TODO
{% endtab %}
{% endtabs %}
</div>

---

## Estimation
The following result provides motivation for maximum likelihood estimation since the value of $f(X; \boldsymbol{\theta})$ will be at its greatest (with high probability) at the true value. 

<div class="theorem">
<strong>Theorem 3.2. (Lehmann and Casella, pg. 444)</strong>
{% tabs thrm-3-2 %}
{% tab thrm-3-2 statement %}
Suppose Conditions 1, 2, 3, and 4 hold with $\boldsymbol{\omega}$ open and $\boldsymbol{\theta}^* \in \text{int}(\boldsymbol{\omega})$. Let $\boldsymbol{\theta}^*$ denote the true value of $\boldsymbol{\theta}$. Then, for any other (fixed) value $\boldsymbol{\theta}_0$ of $\boldsymbol{\theta}$:

$$
\mathbb{P}_{\boldsymbol{\theta}^*} \left(\mathcal{L}(\boldsymbol{\theta}^*; \mathbf{x}) > \mathcal{L}(\boldsymbol{\theta}_0; \mathbf{x}) \right) \rightarrow 1 \hspace{3mm} \text{as} \hspace{3mm} n \rightarrow \infty
$$
{% endtab %}
{% tab thrm-3-2 proof %}
Notice that the event we are considering can be rewritten to:

$$
\begin{aligned}
&\mathcal{L}(\boldsymbol{\theta}^*; \mathbf{X}) > \mathcal{L}(\boldsymbol{\theta}_0; \mathbf{X}) \\
\implies &\prod_{i = 1}^n f(\boldsymbol{\theta}^*; X_i) > \prod_{i = 1}^n f(\boldsymbol{\theta}_0; X_i) \\
\implies & \log\left( \prod_{i = 1}^n f(\boldsymbol{\theta}^*; X_i) \right) > \log\left( \prod_{i = 1}^n f(\boldsymbol{\theta}; X_i) \right) \\
\implies & \sum_{i = 1}^n \log \left(f(X_i; \boldsymbol{\theta}^*)\right) > \sum_{i = 1}^n \log\left( f(X_i; \boldsymbol{\theta}_0) \right) \\
\implies & 0 > \sum_{i = 1}^n \left[ \log\left( f(X_i; \boldsymbol{\theta}_0)  \right) - \log \left(f(X_i; \boldsymbol{\theta}^*) \right) \right] \\
\implies &\sum_{i =1}^n \log\left(\frac{f(X_i; \boldsymbol{\theta}_0)}{f(X_i; \boldsymbol{\theta}^*)} \right) < 0 \\
\implies &\frac{1}{n}\sum_{i =1}^n \log\left(\frac{f(X_i; \boldsymbol{\theta}_0)}{f(X_i; \boldsymbol{\theta}^*)} \right) < 0 
\end{aligned}
$$

By the <a href="https://en.wikipedia.org/wiki/Law_of_large_numbers">law of large numbers</a>:

$$
\frac{1}{n}\sum_{i =1}^n \log\left(\frac{f(X_i; \boldsymbol{\theta}_0)}{f(X_i; \boldsymbol{\theta}^*)} \right) \rightarrow \mathbb{E}_{\boldsymbol{\theta}^*} \left[ \log\left(\frac{f(X_i; \boldsymbol{\theta}_0)}{f(X_i; \boldsymbol{\theta}^*)} \right) \right]
\hspace{3mm} \text{as} \hspace{3mm}
n \rightarrow \infty
$$

Since $\log$ is strictly concave, $-\log$ is strictly convex. Thus, an application of <a href="https://en.wikipedia.org/wiki/Jensen%27s_inequality">Jensen's inequality</a> yields:

$$
\begin{aligned}
&\mathbb{E}_{\boldsymbol{\theta}^*} \left[ -\log\left(\frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)} \right) \right] > -\log\left( \mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)}  \right] \right)  \\
\implies &-\mathbb{E}_{\boldsymbol{\theta}^*} \left[\log\left(\frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)} \right) \right] > -\log\left( \mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)}  \right] \right) \\
\implies &\mathbb{E}_{\boldsymbol{\theta}^*} \left[\log\left(\frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)} \right) \right] < \log\left( \mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)}  \right] \right)
\end{aligned}
$$

We have:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)} \right]
&= \int \frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)} f(X; \boldsymbol{\theta}^*) d X 
= \int f(X; \boldsymbol{\theta}_0) d X 
= 1 \\
\implies \mathbb{E}_{\boldsymbol{\theta}^*} \left[\log\left(\frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)} \right) \right] &< \log\left( \mathbb{E}_{\boldsymbol{\theta}^*} \left[ \frac{f(X; \boldsymbol{\theta}_0)}{f(X; \boldsymbol{\theta}^*)}  \right] \right) = 0
\end{aligned}
$$

Thus, in the limit, the event is always true, which yields the desired result.
{% endtab %}
{% endtabs %}
</div>

This brings us to the definition of the <i>maximum likelihood estimator</i>.

<div id="ml-estimator"></div>
<div class="definition">
<strong>Definition (Maximum Likelihood Estimator; Schervish, pg. 3).</strong>
<br>
Let $\phi$ be an estimator of $\boldsymbol{\theta}$. We call $\phi$ a <i>maximum likelihood estimator (MLE)</i> if it satisfies, for all $x \in \mathcal{X}$:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\} = \mathcal{L}(\phi(x); x)
$$
</div>

<aside><p>Technically, we need $\phi$ to be unique.</p></aside>

Some very nice properties are had by the MLE under the conditions set forth above. 

<div class="theorem">
<strong>Theorem 5.1. (Lehmann and Casella, pg. 463)</strong>
{% tabs thrm-5-1-l %}
{% tab thrm-5-1-l statement %}
Suppose Conditions 1-8 hold with $\boldsymbol{\omega}$ open and $\boldsymbol{\theta}^* \in \boldsymbol{\omega}$ in Condition 2. With probability going to 1 as $n \rightarrow \infty$, there exists $\hat{\boldsymbol{\theta}}_n$ that is the solution to the likelihood equations:

$$
U_{\boldsymbol{\theta}}(\hat{\boldsymbol{\theta}}_n; \mathbf{x}) = \mathbf{0}_p
$$

which satisfies:

<ol>
<li>$\hat{\boldsymbol{\theta}}_i$ is consistent for $\boldsymbol{\theta}_i$</li>
<li>$\sqrt{n}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta})$ is asymptotically normal with mean $\mathbf{0}_p$ and covariance matrix $\mathcal{I}^{-1}(\boldsymbol{\theta})$</li>
<li>$\hat{\boldsymbol{\theta}}_i$ is asymptotically efficient</li>
</ol>
{% endtab %}
{% tab thrm-5-1-l proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>

<!-- ---

## Non-Standard Conditions
The above discussion of maximum likelihood estimation is mostly built off of Abraham Wald's paper <i>Tests of Statistical Hypotheses Concerning Several Parameters When the Number of Observation is Large</i>.<d-cite key=wald1943></d-cite> Notably, there is the major assumption that the true value of $\boldsymbol{\theta}$ lies in an open subset of $\boldsymbol{\Theta}$, which does not have to be the case for certain testing problems. For example, perhaps we are testing variance components in a mixed effects model:

$$
H_0: \sigma^2 = 0 
\hspace{5mm} \text{vs.} \hspace{5mm} 
H_1: \sigma^2 > 0
$$

Because negative variances are not permitted, the value of the parameter of interest lies on the boundary of the natural parameter space (i.e. the space parameter values that yield a valid distribution for our observations). P. A. P. Moran approaches this problem and derives the asymptotic distribution of maximum likelihood estimators.<d-cite key=moran1971></d-cite>

### Set-Up
We will basically assume the same set-up as before with the additional condition that $P_{\boldsymbol{\theta}}(X)$ is either absolutely continuous in $\mathcal{X}$ or purely discrete with a countable number of points of increase. For simplicity, we will concern ourselves with a closed and bounded subset $\boldsymbol{\Theta}_0 \subseteq \boldsymbol{\Theta}$ defined as:

$$
\boldsymbol{\Theta}_0 = \left\{ \boldsymbol{\theta} \in \boldsymbol{\Theta} : 0 \leq \boldsymbol{\theta}_i \leq b_i; \hspace{1mm} b_i > 0 \hspace{1mm}; i = 1, \dots, p \right\}
$$

We will consider the testing problem:

$$
\begin{equation}
\label{eq:hypotheses-2}
H_0: \boldsymbol{\theta}_1 = 0, \hspace{1mm} 0 < \boldsymbol{\theta}_i < b_i \hspace{2mm} i = 2, \dots, p
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: 0 < \boldsymbol{\theta}_i \leq b_i \hspace{2mm} \forall i
\end{equation}
$$

We will need to define a couple of functions before we begin.

$$
\begin{align}
\psi_{i,j}(x, \boldsymbol{\theta}', \delta) &= \underset{\boldsymbol{\theta} : \rvert \rvert \boldsymbol{\theta} - \boldsymbol{\theta}' \rvert \rvert_{\infty} \leq \delta}{\inf} \left\{ \frac{\partial^2 \log (f(\mathbf{x}; \boldsymbol{\theta}))}{\boldsymbol{\theta}_i \boldsymbol{\theta}_j} \right\} \label{eq:inf-func} \\
\phi_{i,j}(x, \boldsymbol{\theta}', \delta) &= \underset{\boldsymbol{\theta} : \rvert \rvert \boldsymbol{\theta} - \boldsymbol{\theta}' \rvert \rvert_{\infty} \leq \delta}{\sup} \left\{ \frac{\partial^2 \log (f(\mathbf{x}; \boldsymbol{\theta}))}{\boldsymbol{\theta}_i \boldsymbol{\theta}_j} \right\} \label{eq:sup-func}
\end{align}
$$

### Regularity Assumptions
Below are Moran's assumptions needed for his results.

#### Assumption 1
$f(x; \boldsymbol{\theta})$ is continuous in $\boldsymbol{\theta}$ on $\boldsymbol{\Theta}_0$.

This assumption ensures that at least one maximum likelihood estimator exists. 

#### Assumption 2
Let $\mathbf{D}_n \subseteq \mathcal{X}$ be the set of points such that $\frac{\partial^2 f(x_l; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}$ is continuous in $\boldsymbol{\theta}$ for $l = 1, \dots, n$ and $i, j = 1, \dots, p$. The probability of $\mathbf{D}_n$ is one.

This assumptions states that the second-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to all coordinates of $\boldsymbol{\theta}$ are continuous in $\boldsymbol{\theta}$ almost everywhere. 

<aside><p>I guess almost surely since we are dealing with a probability space...</p></aside>

#### Assumption 3
The maximum likelihood estimator $\hat{\boldsymbol{\theta}}$ is consistent uniformly for $\boldsymbol{\Theta}_0$.

See Moran (1971) for a discussion of what this means and when it occurs.<d-cite key=moran1971b></d-cite>

<aside><p><i>Add definition?</i></p></aside>

#### Assumption 4
For any two sequences of values in $\boldsymbol{\Theta}_0$, $\{ \boldsymbol{\theta}^1(n) \}$ and $\{ \boldsymbol{\theta}^2 \}$, and and any sequence, $\{ \delta(n)\}$, such that:

$$
\begin{aligned}
\underset{n \rightarrow \infty}{\lim} \boldsymbol{\theta}^1(n) &= \underset{n \rightarrow \infty}{\lim} \boldsymbol{\theta}^2(n) = \boldsymbol{\theta} \in \boldsymbol{\Theta}_0 \\
\underset{n \rightarrow \infty}{\lim} \delta(n) &= 0 
\end{aligned}
$$

we have, uniformly in $\boldsymbol{\theta}$:

$$
\underset{n \rightarrow \infty}{\lim} \mathbb{E}_{\boldsymbol{\theta}^1(n)} \left[ \psi_{i,j}(X, \boldsymbol{\theta}^2(n), \delta(n)) \right] = \underset{n \rightarrow \infty}{\lim} \mathbb{E}_{\boldsymbol{\theta}^1(n)}\left[ \phi_{i,j}(X, \boldsymbol{\theta}^2(n), \delta(n)) \right] = \mathbb{E}\left[ \frac{\partial^2 \log (f(X; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right]
$$

This assumption has something to do with completeness/Cauchy-ness, but I don't feel like figuring it out right now.

#### Assumption 5
There exists $\epsilon > 0$ such that $$\mathbb{E}_{\boldsymbol{\theta}^1} \left[ \psi^2_{i,j}(X, \boldsymbol{\theta}^2, \delta) \right]$$ and $$\mathbb{E}_{\boldsymbol{\theta}^1} \left[ \phi^2_{i,j}(X, \boldsymbol{\theta}^2, \delta) \right]$$ are bounded functions of $\boldsymbol{\theta}^1$, $\boldsymbol{\theta}^2$, and $\delta$ over $\boldsymbol{\Theta}'_0$, the subset defined as:

$$
\boldsymbol{\Theta}'_0 = \left\{ \boldsymbol{\theta}^1, \boldsymbol{\theta}^2 \in \boldsymbol{\Theta}_0 : \rvert \rvert \boldsymbol{\theta}^1 - \boldsymbol{\theta}^2 \rvert \rvert_\infty \leq \epsilon; \hspace{1mm} 0 \leq \delta \leq \epsilon \right\}
$$

#### Assumption 6
The following holds:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}_0}{\inf} \left\{ \left\rvert - \mathbb{E}_{\boldsymbol{\theta}}\left[ \frac{\partial^2 \log f(X; \boldsymbol{\theta})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right] \right\rvert \right\} > 0
$$

#### Assumption 7
For all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ and $i, j= 1, \dots, p$:

$$
\int \frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i} dx = \int \frac{\partial^2 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} dx = 0
$$

<aside><p>These are summations if $P_{\boldsymbol{\theta}}(X)$ is discrete.</p></aside>

This assumption is basically that of 6, 6a, and 6b above (exchanging the order of differentiation and integration). It ensures that the score will have mean zero and that the Fisher information can be computed as the negative expected Hessian.

#### Assumption 8
There exists $\eta > 0$ such that:

$$
\mathbb{E}_\boldsymbol{\theta} \left[ \left\rvert \frac{\partial \log(f(X; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i} \right\rvert^{2 + \eta} \right];
\hspace{5mm}
i = 1, \dots, p
$$

<aside><p>Lyapunov CLT?</p></aside>

are bounded functions of $\boldsymbol{\theta} \in \boldsymbol{\Theta}_0$.

### Main Results
There are a few different results in Moran's paper, but we only cover Theorem 1 here. Let $\hat{\boldsymbol{\theta}}$ and $$\boldsymbol{\theta}^*$$ denote the likelihood estimator and the true value of $\boldsymbol{\theta}$, respectively.

<div class="theorem">
<strong>Theorem 1 (Moran, 1971).</strong>
{% tabs thrm-1-m %}
{% tab thrm-1-m statement %}
Define:

$$
\begin{aligned}
c_{i,j} &= - \mathbb{E}_{\boldsymbol{\theta}}\left[ \frac{\partial^2 \log(f(X; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right]; \hspace{2mm} i,j = 1, \dots, p \\
z^n_i &= \sqrt{n}(\hat{\boldsymbol{\theta}}_i - \boldsymbol{\theta}^*_i); \hspace{2mm} i = 1, \dots, p
\end{aligned}
$$

Under $H_0$ described in Eq. \eqref{eq:hypotheses-2}:

$$
\mathbb{P}_{\boldsymbol{\theta}^*}\left( z_n < t \right) \rightarrow \mathbb{P}_{\boldsymbol{\theta}^*}\left( z < t \right) = \frac{1}{2}\mathcal{N}() + \frac{1}{2}
$$

uniformly in $t$ and $\boldsymbol{\theta}$.
{% endtab %}
{% tab thrm-1-m proof %}
In what follows, we assume $H_0$ holds. We will denote the complete data with $\mathbf{x}^n = (x_1, \dots, x_n)$.
<br><br>
By Assumption 7, the score has expectation zero, and the Fisher information (which has elements equal to $c_{i,j}$) can be computed as the negative expected Hessian. This matrix, by Assumption 6, is positive definite, and any sub-matrix obtained by removing any number of rows and the columns with the same indices <a href="https://ximera.osu.edu/oerlinalg/LinearAlgebra/RTH-0045/main#thm:024830">will also be positive definite</a>.
<br><br>
Assumption 3 states the the MLE $\hat{\boldsymbol{\theta}}$ is consistent uniformly on $\boldsymbol{\Theta}_0$. This means that, as we take $n \rightarrow \infty$, the probability that $$\rvert \rvert \hat{\boldsymbol{\theta}} - \boldsymbol{\theta}^* \rvert \rvert_\infty$$ is arbitrarily small will go to one. This implies that, with probability going to one, $\hat{\boldsymbol{\theta}}$ will lie in any neighborhood of $$\boldsymbol{\theta}^*$$ that is open relative to $\boldsymbol{\Theta}_0$. The reason we have open neighborhood is because any point arbitrarily close to $$\boldsymbol{\theta}^*$$ that is not in a given neighborhood will also not be in $\boldsymbol{\Theta}_0$.
<br><br>
By Assumption 2, the first- and second-order partial derivatives of $f(x; \boldsymbol{\theta})$ exist almost everywhere. Since $\hat{\boldsymbol{\theta}}$ maximizes the likelihood, it will be a stationary point of the log-likelihood. Thus:

$$
\begin{aligned}
\left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i} \right\rvert_{\boldsymbol{\theta} = \hat{\boldsymbol{\theta}}} &= 0; \hspace{2mm} i = 2, \dots, p
\end{aligned}
$$

If $\hat{\boldsymbol{\theta}}_1 = 0$, then the derivative $$\left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_1} \right\rvert_{\boldsymbol{\theta} = \hat{\boldsymbol{\theta}}}$$ will be negative since we must take the limit from the right heading towards zero. The log-likelihood should be increasing as we approach zero, which means that it is decreasing generally. If $\hat{\boldsymbol{\theta}}_1 > 0$, then the derivative will be zero like in the cases of $i = 2, \dots, p$.
<br><br>
We now proceed in cases. First suppose $\hat{\boldsymbol{\theta}}_1 > 0$. We now define:

$$
y_i^n = \frac{1}{\sqrt{n}} \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i} \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}^*}
$$

By Assumption 2 and the <a href="https://en.wikipedia.org/wiki/Mean_value_theorem">mean value theorem</a>, there exists $\bar{\boldsymbol{\theta}}$ lying somewhere between $\boldsymbol{\theta}$ and $\hat{\boldsymbol{\theta}}$ satisfying:

$$
\begin{aligned}
\left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i} \right\rvert_{\boldsymbol{\theta} = \hat{\boldsymbol{\theta}}} 
  &= \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i} \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}^*} + \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}^\top} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}^*) \\
\implies
\left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}} \right\rvert_{\boldsymbol{\theta} = \hat{\boldsymbol{\theta}}} 
  &= \sqrt{n} y_i^n + \frac{1}{\sqrt{n}} \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} \mathbf{z}^n \\
\implies 
\frac{1}{\sqrt{n}} \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}} \right\rvert_{\boldsymbol{\theta} = \hat{\boldsymbol{\theta}}} 
  &= y_i^n + \frac{1}{n} \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} \mathbf{z}^n 
\end{aligned}
$$

for all $x_l \in \mathbf{D}_n$. However, we know that the LHS of the previous equation equals $\mathbf{0}_p$ by the argument made earlier in the proof. Thus:

$$
\begin{aligned}
0 &= y_i^n + \frac{1}{n} \left. \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} \mathbf{z}^n  \\
\implies
\mathbf{z}^n &= - y_i^n \left[\left. \frac{1}{n} \frac{\partial \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} \right]^{-1} \\
\implies 
\mathbf{z}^n_j &= -
\end{aligned}
$$

We will pick some very small $\nu > 0$ and define:

$$
\mathbf{Q}_n(\boldsymbol{\theta}) = \left\{ \mathbf{x}^n \in \mathbf{D}_n : \left\rvert \frac{1}{n} \left. \frac{\partial^2 \log(f(\mathbf{x}^n; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} + c_{i,j} \right\rvert < \nu \right\}
$$

It is shown in Wald (1943), pg. 431<d-cite key=wald1943></d-cite> that, uniformly in $\boldsymbol{\theta}$:

$$
\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left(\mathbf{Q}_n(\boldsymbol{\theta}) \rvert \boldsymbol{\theta} \right) = 1
$$

Notice that we can rewrite the condition defining $\mathbf{Q}_n(\boldsymbol{\theta})$ as:

$$
\begin{aligned}
&\left\rvert \frac{1}{n} \left. \sum_{l = 1}^n \frac{\partial^2 \log(f(x_l; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} + c_{i,j} \right\rvert < \nu  \\
&= 
\end{aligned}
$$

It then follows that, for samples in $\mathbf{Q}_n(\boldsymbol{\theta})$ and conditional on $$z^n_1 = \sqrt{n}(\hat{\boldsymbol{\theta}}_1 - \boldsymbol{\theta}^*_1) = \sqrt{n}\hat{\boldsymbol{\theta}}_1 > 0$$:

$$
\begin{aligned}
0 &= y_i^n + \sum_{j = 1}^p z_j^n \frac{1}{n}  \left( \left. \sum_{l = 1}^n \frac{\partial \log(f(x_l; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} \right) \\
\implies 0 &= y_i^n 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

$$
\begin{aligned}
0 &= y_i^n + \sum_{j = 1}^p \frac{1}{n} z_j^n \left. \sum_{l = 1}^n \frac{\partial \log(f(x_l; \boldsymbol{\theta}))}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right\rvert_{\boldsymbol{\theta} = \bar{\boldsymbol{\theta}}} \\

\end{aligned}
$$ -->




































































---

## Results


We can also define the <i>efficiency</i> of an estimator. 

<div class="definition">
<strong>Definition (Efficiency).</strong>
<br>
Let $\boldsymbol{\delta}_1$ and $\boldsymbol{\delta}_2$ be two estimators of $\boldsymbol{\theta}$ such that:

$$
\begin{aligned}
\sqrt{n}(\boldsymbol{\delta}_1 - \boldsymbol{\theta}) &\rightsquigarrow \mathcal{N}\left(\mathbf{0}_k, \boldsymbol{\Sigma}_1(\boldsymbol{\theta}) \right) \\
\sqrt{n}(\boldsymbol{\delta}_2 - \boldsymbol{\theta}) &\rightsquigarrow \mathcal{N}\left(\mathbf{0}_k, \boldsymbol{\Sigma}_2(\boldsymbol{\theta}) \right) 
\end{aligned}
$$

We call $\boldsymbol{\delta}_1$ <i>more efficient</i> than $\boldsymbol{\delta}_2$ if: 

$$\boldsymbol{\Sigma}_2(\boldsymbol{\theta}) - \boldsymbol{\Sigma}_1(\boldsymbol{\theta})$$

is positive semi-definite for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$. 
</div>

Similarly, we call $\hat{\boldsymbol{\theta}}$ an <strong>efficient estimator</strong> of $\boldsymbol{\theta}$ if it is asymptotically normal with mean zero and non-singular covariance matrix $\boldsymbol{\Sigma}(\boldsymbol{\theta}) = \mathcal{I}^{-1}(\boldsymbol{\theta})$. An example is the estimator discussed in Theorem 7.5.2. (<a href="https://en.wikipedia.org/wiki/Cramér–Rao_bound">Cramér-Rao Bound</a>!). 

<aside><p>Only if we also assume that $\boldsymbol{\Sigma}_{i,j}(\boldsymbol{\theta})$ and $\mathcal{I}_{i,j}^{-1}(\boldsymbol{\theta})$ are continuous functions of $\boldsymbol{\theta}$ for all $i,j$.</p></aside>

The results above can easily be extended to cases where we have multiple independent samples or if we have nuisance parameters. In the latter case, the efficiency will generally decrease with the number of nuisance parameters unless the parameter of interest is independent of them.






---

## Newey and McFadden
In this section, we work up to the asymptotic normality of the maximum likelihood estimator as shown in Newey and McFadden.<d-cite key=newey1994></d-cite>

<div class="definition">
<strong>Definition (Extremum Estimator).<d-cite key=newey1994></d-cite></strong>
Let $\boldsymbol{\Theta}$ denote a parameter space, let $\boldsymbol{\theta}$, and let $\hat{Q}_n(\boldsymbol{\theta})$ be an objective function. An estimator, $\hat{\boldsymbol{\theta}}$ is called an <i>extremum estimator</i> if there exists a $\hat{Q}_n(\boldsymbol{\theta})$ such that: 

$$
\begin{equation}
\label{eq:ext-est}
\hat{\boldsymbol{\theta}} = \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\arg \max} \left\{ \hat{Q}_n(\boldsymbol{\theta}) \right\}
\end{equation}
$$
</div>

One example of an extremum estimator is the maximum likelihood estimator. The corresponding objective function is the likelihood (or log-likelihood). 


### Consistency
Before we begin to show when the maximum likelihood estimator is consistent, we must define what it means for an objective function for an extremum estimator to <strong>uniformly converge in probability</strong>.

<div class="definition">
<strong>Definition (Uniform Convergence in Probability).<d-cite key=newey1994></d-cite></strong>
An objective function $\hat{Q}_n(\boldsymbol{\theta})$ <i>converges uniformly in probability</i> to $Q_0(\boldsymbol{\theta})$ if:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \hat{Q}_n(\boldsymbol{\theta}) - Q_0(\boldsymbol{\theta}) \right\} \overset{p}{\rightarrow} 0
$$
</div>

With this property of $\hat{Q}_n(\boldsymbol{\theta})$ (among others), we can show consistency of an extremum estimator on any topological space. 

<div id="theorem-2-1"></div>
<div class="theorem">
<strong>Theorem 2.1. (Newey and McFadden, pg. 2121)</strong>
{% tabs thrm-2-1 %}
{% tab thrm-2-1 statement %}
Let $\hat{\boldsymbol{\theta}}$ be an extremum estimator that maximizes objective function $\hat{Q}_n(\boldsymbol{\theta})$. If there exists a function $Q_0(\boldsymbol{\theta})$ such that:

<ol>
<li>$Q_0(\boldsymbol{\theta})$ has a unique maximizer, $\boldsymbol{\theta}_0$</li>
<li>The parameter space $\boldsymbol{\Theta}$ is compact</li>
<li>$Q_0(\boldsymbol{\theta})$ is continuous</li>
<li>$\hat{Q}_n(\boldsymbol{\theta})$ converges uniformly in probability to $Q_0(\boldsymbol{\theta})$</li>
</ol>

Then $\hat{\boldsymbol{\theta}}$ is consistent for $\boldsymbol{\theta}$; i.e. $\hat{\boldsymbol{\theta}} \overset{p}{\rightarrow} \boldsymbol{\theta}$.

<i>Technically, we also need $\hat{\boldsymbol{\theta}}$, $\hat{Q}_n(\hat{\boldsymbol{\theta}})$, and $\hat{Q}_n(\hat{\boldsymbol{\theta}}_0)$ to be measurable.</i>
{% endtab %}
{% tab thrm-2-1 proof %}
By definition, $\hat{\boldsymbol{\theta}} = \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\arg \max} \left\\{ \hat{Q}_n(\boldsymbol{\theta}) \right\\}$. Then, for any $\epsilon > 0$ (WHY?):

$$
\begin{aligned}
\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \hat{Q}_n(\hat{\boldsymbol{\theta}}) > \hat{Q}_n(\boldsymbol{\theta}_0) - \frac{\epsilon}{3} \right) &= 1 \\
\implies \underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \hat{Q}_n(\hat{\boldsymbol{\theta}}) - \frac{\epsilon}{3} > \hat{Q}_n(\boldsymbol{\theta}_0) - \frac{2 \epsilon}{3} \right) &= 1
\end{aligned}
$$

By uniform convergence in probability of $\hat{Q}_n(\boldsymbol{\theta})$, we also have:

$$
\begin{aligned}
\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( Q_0(\hat{\boldsymbol{\theta}}) > \hat{Q}_n(\hat{\boldsymbol{\theta}}) - \frac{\epsilon}{3} \right) &= 1 
\end{aligned}
$$

And: 

$$
\begin{aligned}
\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \hat{Q}_n(\boldsymbol{\theta}_0) > Q_0(\boldsymbol{\theta}_0) - \frac{\epsilon}{3} \right) &= 1 \\
\implies \underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \hat{Q}_n(\boldsymbol{\theta}_0) - \frac{2\epsilon}{3} > Q_0(\boldsymbol{\theta}_0) - \epsilon \right) &= 1 \\
\end{aligned}
$$

Then, with probability approaching 1, we have:

$$
\begin{aligned}
Q_0(\hat{\boldsymbol{\theta}}) &> \hat{Q}_n(\hat{\boldsymbol{\theta}}) - \frac{\epsilon}{3} > \hat{Q}_n(\boldsymbol{\theta}_0) - \frac{2 \epsilon}{3} > Q_0(\boldsymbol{\theta}_0) - \epsilon
\end{aligned}
$$

Since $\boldsymbol{\Theta}$ is compact, let us choose some open subset $\mathbf{C} \subset \boldsymbol{\Theta}$ such that $\boldsymbol{\theta}_0 \in \mathbf{C}$. The complement of $\mathbf{C}$, denoted by $\mathbf{C}^c$, is closed (by definition of open/closed). Thus, $\boldsymbol{\Theta} \cap \mathbf{C}^c$ will also be compact (see <a href="https://math.stackexchange.com/questions/35038/is-the-intersection-of-a-closed-set-and-a-compact-set-always-compact">here</a>). 

Since $Q_0(\boldsymbol{\theta})$ is continuous and $\boldsymbol{\Theta} \cap \mathbf{C}^c$ is compact, $Q_0(\boldsymbol{\theta})$ will achieve its supremum at some $\boldsymbol{\theta}^* \in \boldsymbol{\Theta} \cap \mathbf{C}^c$ and:

$$
Q_0(\boldsymbol{\theta}^*) < Q_0(\boldsymbol{\theta}_0)
$$

Let:

$$
\epsilon = Q_0(\boldsymbol{\theta}_0) - Q_0(\boldsymbol{\theta}^*)
$$

Then, with probability approaching 1:

$$
Q_0(\hat{\boldsymbol{\theta}}) &> Q_0(\boldsymbol{\theta}_0) - \epsilon  \\
\implies
Q_0(\hat{\boldsymbol{\theta}}) &> Q_0(\boldsymbol{\theta}_0) - \left[ Q_0(\boldsymbol{\theta}_0) - Q_0(\boldsymbol{\theta}^*) \right] \\ 
\implies
Q_0(\hat{\boldsymbol{\theta}}) &> Q_0(\boldsymbol{\theta}^*)
$$

This implies $\hat{\boldsymbol{\theta}} \notin \mathbf{C}^c$, which implies it is in $\mathbf{C}$. 
{% endtab %}
{% endtabs %}
</div>

<aside><p>Why do we need continuity of the objective function limit ($Q_0(\boldsymbol{\theta})$)?</p></aside>

We call Conditions 3 and 4 (continuity and uniform convergence in probability) in the previous theorem <i>regularity conditions</i> while 1 and 2 (identification and boundedness) are called <i>substantive</i>. In order to apply the above result, the conditions must be verified, which can be difficult to do. Instead, we often verify simpler properties. 

Regarding Condition 1, the maximum likelihood estimator has the additional nice property that satisfying the identification condition implies a <i>unique maximum</i>.

<div class="theorem">
<strong>Lemma 2.2 (Newey and McFadden, pg. 2124)</strong>
{% tabs lemma-2-2 %}
{% tab lemma-2-2 statement %}
If $\boldsymbol{\theta} \neq \boldsymbol{\theta}_0$ for $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ implies that $f(x; \boldsymbol{\theta}) \neq f(x; \boldsymbol{\theta}_0)$ (i.e. $\boldsymbol{\theta}_0$ is <i>identified</i>) and $\mathbb{E}\left[ \rvert \log(f(x; \boldsymbol{\theta})) \rvert \right] < \infty$ for all $\boldsymbol{\theta}$, then $Q_0(\boldsymbol{\theta}) = \mathbb{E}[\log(f(x; \boldsymbol{\theta}))]$ has a unique maximum at $\boldsymbol{\theta}$. 
{% endtab %}
{% tab lemma-2-2 proof %}
For a non-constant, positive random variable, $-\log(\cdot)$ is a strictly convex function. Clearly $\frac{f(x; \boldsymbol{\theta})}{f(x; \boldsymbol{\theta}_0)} \in [0, 1]$ and will not equal 1 if $\boldsymbol{\theta} \neq \boldsymbol{\theta}_0$. In thise case, we can apply Jensen's inequality to see:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}_0}\left[ - \log\left( \frac{f(x; \boldsymbol{\theta})}{f(x; \boldsymbol{\theta}_0)}\right) \right]
  &> - \log\left( \mathbb{E}_{\boldsymbol{\theta}_0}\left[ \frac{f(x; \boldsymbol{\theta})}{f(x; \boldsymbol{\theta}_0)}\right] \right) \\
  &= - \log \left( \int \frac{f(x; \boldsymbol{\theta})}{f(x; \boldsymbol{\theta}_0)} f(x; \boldsymbol{\theta}_0) dx \right) \\
  &= - \log \left( \int f(x; \boldsymbol{\theta}) dx \right) \\
  &= -\log(1) \\
  &= 0
\end{aligned}
$$

Thus, $Q_0(\boldsymbol{\theta}_0) > Q_0(\boldsymbol{\theta})$, so $\boldsymbol{\theta}_0$ is a unique maximizer.
{% endtab %}
{% endtabs %}
</div>

<aside><p>The above result is often called the <strong>information inequality</strong>.</p></aside>

The above lemma shows the needed conditions for Condition 1 in <a href="#theorem-2-1">Theorem 2.1</a> to be satisfied when $\hat{\boldsymbol{\theta}}$ is the maximum likelihood estimator.

Condition 2 is often left unchecked because there are many cases where it is not actually needed for proving consistency in all cases. 

For Conditions 3 and 4, we use the following lemma which allows us to simply inspect $\mathbf{A}(x, \boldsymbol{\theta})$ and construct an upper bound on its Euclidean norm that may depend on the data.

<div class="theorem">
<strong>Lemma 2.4. (Newey and McFadden, pg. 2129)</strong>
{% tabs thrm-2-4 %}
{% tab thrm-2-4 statement %}
Let $\mathbf{A}(x, \boldsymbol{\theta})$ be a matrix whose elements, denoted by $\mathbf{A}_{i,j}(x, \boldsymbol{\theta})$, are functions of an observation, $x$, and the parameter, $\boldsymbol{\theta}$. Suppose:

<ol>
<li>$X_1, \dots, X_n$ are i.i.d.</li>
<li>$\boldsymbol{\Theta}$ is compact</li>
<li>For all $i,j$, $\mathbf{A}_{i,j}(x, \boldsymbol{\theta})$ is continuous at each $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ with probability one</li>
<li>There exists $d(x)$ with $\mathbb{E}[d(x)] < \infty$ such that $\rvert \rvert \mathbf{A}_{i,j}(x, \boldsymbol{\theta}) \rvert \rvert_2 \leq d(x)$ for all $\boldsymbol{\theta}$</li>
</ol>

Then $\mathbb{E}[\mathbf{A}(x, \boldsymbol{\theta})]$ is continuous and:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \rvert \rvert \frac{1}{n} \sum_{i = 1}^n \left[ \mathbf{A}(x_i, \boldsymbol{\theta}) - \mathbb{E}\left[ \mathbf{A}(x, \boldsymbol{\theta})\right] \right] \rvert \rvert_2 \right\} \overset{p}{\rightarrow} 0
$$
{% endtab %}
{% tab thrm-2-4 proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>

<aside><p>This is a uniform law of large numbers.</p></aside>

However, it may be illustrative to explore when uniform convergence in probability holds. First, we must introduce the concept of <i>stochastic equicontinuity</i>.

<div class="definition">
<strong> Definition (Stochastic Equicontinuity).<d-cite key=newey1994></d-cite></strong>
For every $\epsilon, \eta > 0$, there exists a sequence of random variable $\{ \hat{\Delta} \}_{n = 0}^\infty$ and sample size $n_0$ such that, for $n \geq n_0$, $\mathbb{P}(\rvert \hat{\boldsymbol{\Delta}}_n \rvert > \epsilon) < \eta$ and for each $\boldsymbol{\theta}$, there exists an open set $\mathbf{C}$ such that $\boldsymbol{\theta} \in \mathbf{C}$ an:

$$
\underset{\tilde{\boldsymbol{\theta}} \in \mathbf{C}}{\sup} \left\{ \rvert \hat{Q}_n(\tilde{\boldsymbol{\theta}}) - \hat{Q}_n(\boldsymbol{\theta}) \right\} \leq \hat{\boldsymbol{\Delta}}_n
\hspace{5mm} \text{for } n \geq n_0
$$
</div>

Stochastic equicontinuity is an extension of the concept of equicontinuity for non-random quantities to random ones. Intuitively, equicontinuity means that the values of a function at two different points can be made arbitrarily close by picking the two points to be close enough together. Similarly, stochastic equicontinuity can be interpreted as putting a bound on how much the output of $\hat{Q}_n(\cdot)$ can change when we consider $\boldsymbol{\theta}$ and points close to $\boldsymbol{\theta}$. However, the "stochastic" part of the name means that with arbitrarily high probability, this bound will arbitrarity small. 

This brings us to our next result.

<div class="theorem">
<strong>Lemma 2.8 (Newey and McFadden, pg. 2137)</strong>
{% tabs thrm-2-8 %}
{% tab thrm-2-8 statement %}
Let $\hat{\boldsymbol{\theta}}$ be an extremum estimator with objective function $\hat{Q}_n(\boldsymbol{\theta})$. Suppose $\boldsymbol{\Theta}$ is compact, and $Q_0(\boldsymbol{\theta})$ is continuous. Then:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \rvert \hat{Q}_n(\boldsymbol{\theta}) - Q_0(\boldsymbol{\theta}) \right\} \overset{p}{\rightarrow} 0
$$

if, and only if, $\hat{Q}_n(\boldsymbol{\theta}) \overset{p}{\rightarrow} Q_0(\boldsymbol{\theta})$ for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ and $\hat{Q}_n(\boldsymbol{\theta})$ is stochastically equicontinuous.
{% endtab %}
{% tab thrm-2-8 proof %}
See Newey (1991).<d-cite key=newey1991></d-cite>
{% endtab %}
{% endtabs %}
</div>



We can instead show a summary result that states "easy to check" conditions that are needed for consistency of the maximum likelihood estimator.

<aside><p>Newey and McFadden refer to such conditions as <i>primitive</i> (easy to interpret).</p></aside>

<div class="theorem">
<strong>Theorem 2.5. (Newey and McFadden, pg. 2131)</strong>
{% tabs thrm-2-5 %}
{% tab thrm-2-5 statement %}
Let $\hat{\boldsymbol\boldsymbol{\theta}}$ is the MLE of $\boldsymbol{\theta}$ based upon i.i.d. $X_1, X_2, \dots$ with probability density function $f(x_i; \boldsymbol{\theta}_0)$. Suppose the following hold:

<ol>
<li>$\boldsymbol{\theta} \neq \boldsymbol{\theta}_0 \implies f(x_i; \boldsymbol{\theta}) \neq f(x_i; \boldsymbol{\theta}_0)$</li>
<li>$\boldsymbol{\theta}_0 \in \boldsymbol{\Theta}$ and $\boldsymbol{\Theta}_0$ is compact</li>
<li>$\log(f(x_i; \boldsymbol{\theta}))$ is continuous at each $\boldsymbol{\theta} \in \boldsymbol{\Theta}$ with probability 1</li>
<li>$\mathbb{E}\left[ \underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \rvert \log(f(x; \boldsymbol{\theta})) \rvert \right\} \right] < \infty$</li>
</ol>

Then $\hat{\boldsymbol{\theta}}$ is consistent for $\boldsymbol{\theta}$; that is, $\hat{\boldsymbol{\theta}} \overset{p}{\rightarrow} \boldsymbol{\theta}_0$. 
{% endtab %}
{% tab thrm-2-5 proof %}
The conditions listed ensure that the conditions of Theorem 2.1, Lemma 2.2, and Lemma 2.4 hold. Thus, we can claim the result in Theorem 2.1.
{% endtab %}
{% endtabs %}
</div>

A nice fact is that the compactness condition can be relaxed for the MLE because the objective function is concave.

<div class="theorem">
<strong>Theorem 2.7 (Newey and McFadden, pg. 2133)</strong>
{% tabs thrm-2-7 %}
{% tab thrm-2-7 statement %}
Let $\hat{\boldsymbol{\theta}}$ be an extremum estimator that maximizes objective function $\hat{Q}_n(\boldsymbol{\theta})$. If there exists a function $Q_0(\boldsymbol{\theta})$ such that:

<ol>
<li>$Q_0(\boldsymbol{\theta})$ has a unique maximizer, $\boldsymbol{\theta}_0$</li>
<li>$\boldsymbol{\theta}_0 \in \text{int}(\boldsymbol{\Theta})$ with convex $\boldsymbol{\Theta}$</li>
<li>$\hat{Q}_n(\boldsymbol{\theta})$ is concave</li>
<li>$\hat{Q}_n(\boldsymbol{\theta}) \overset{p}{\rightarrow} Q_0(\boldsymbol{\theta})$ for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$</li>
</ol>

Then $\hat{\boldsymbol{\theta}}$ will exist with probability going to 1, and $\hat{\boldsymbol{\theta}} \overset{p}{\rightarrow} \boldsymbol{\theta}$ with probability going to 1.
{% endtab %}
{% tab thrm-2-7 proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>


### Asymptotic Normality


<div class="theorem">
<strong>Theorem 3.1. (Newey and McFadden, pg. 2143)</strong>
{% tabs thrm-3-1 %}
{% tab thrm-3-1 statement %}
Let $\hat{\boldsymbol{\theta}}$ satisfy Eq. \eqref{eq:ext-est} and the following conditions:

<ol>
<li>$\hat{\boldsymbol{\theta}} \overset{p}{\rightarrow} \boldsymbol{\theta}$</li>
<li>$\boldsymbol{\theta}_0 \in \text{int}(\boldsymbol{\Theta})$</li>
<li>$\hat{Q}_n(\boldsymbol{\theta})$ is twice continuously differntiable in a neighborhood $\mathcal{K}$ of $\boldsymbol{\theta}_0$</li>
<li>$\sqrt{n} \left. \frac{\partial\hat{Q}_n(\boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \boldsymbol{\Sigma})$</li>
<li>There exists a matrix $\mathbf{H}(\boldsymbol{\theta})$ that is continuous at $\boldsymbol{\theta}_0$ and satisfies $\underset{\boldsymbol{\theta} \in \mathcal{K}}{\sup} \left\{ \rvert \rvert\frac{\partial^2 \hat{Q}_n(\boldsymbol{\theta})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} - \mathbf{H}(\boldsymbol{\theta}) \rvert \rvert \right\} \overset{p}{\rightarrow} 0$</li>
<li>$\mathbf{H}(\boldsymbol{\theta}_0)$ is non-singular</li>
</ol>

Then $\hat{\boldsymbol{\theta}}$ is asymptotically normal. That is:

$$
\sqrt{n}(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \mathbf{H}^{-1}(\boldsymbol{\theta}_0) \boldsymbol{\Sigma} \mathbf{H}^{-1}(\boldsymbol{\theta}_0))
$$
{% endtab %}
{% tab thrm-3-1 proof %}
TO DO: See 3.5 in Newey.
{% endtab %}
{% endtabs %}
</div>

Condition 4 is the asymptotic normality of the score function. This can be ensured under the conditions discussed in <a href="#theorem-7-5-1">Theorem 7.5.1</a> above. Condition 5 is the continuity of the matrix of second order partial derivatives at the true value and its uniform convergence over a neighborhood around the true parameter value. Condition 6 ensures identifiability of $\boldsymbol{\theta}_0$ as a unique local maximum of the objective function. If <a href="#condition-6">Condition 6</a> and <a href="#condition-6a">6a</a> above hold, then we can exchange the order of integration and differentiation so that $\mathbf{H}(\boldsymbol{\theta}_0) = \frac{}{}.




---

<div class="theorem">
<strong>Theorem 3.2. (Lehmann and Casella, pg. 444)</strong>
{% tabs thrm-3-2 %}
{% tab thrm-3-2 statement %}
{% endtab %}
{% tab thrm-3-2 proof %}
{% endtab %}
{% endtabs %}
</div>




---
We can estimate $\boldsymbol{\theta}$ by finding the value that maximizes $\mathcal{L}(\boldsymbol{\theta}; \mathbf{x})$, which is called the <strong>maximum likelihood estimate (MLE)</strong>. Intuitively, the MLE is the value that is most probable when it comes to having generated the data at hand.

<div id="ml-estimator"></div>
<div class="definition">
<strong>Definition (Maximum Likelihood Estimator; Schervish, pg. 3).</strong>
<br>
Let $\phi$ be an estimator of $\Theta$. We call $\phi$ a <i>maximum likelihood estimator (MLE)</i> if it satisfies, for all $x \in \mathcal{X}$:

$$
\underset{\boldsymbol{\theta} \in \boldsymbol{\Theta}}{\sup} \left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\} = \mathcal{L}(\phi(x); x)
$$
</div>

Though we do not know the true value of $\Theta$, we can try to find the MLE. This choice of estimator makes sense because of the following theorem.

<div class="theorem">
  <strong>Theorem 3.2.<d-cite key=lehmann2005></d-cite></strong>
  {% tabs thrm-3-2 %}
  {% tab thrm-3-2 statement %}
  Assume Conditions 1-4 (given <a href="#regularity-conditions">below</a>) hold. Then:

  $$
  \mathbb{P}_{\boldsymbol{\theta}^*}\left( \mathcal{L}(\boldsymbol{\theta}^*; \mathbf{x}) > \mathcal{L}(\boldsymbol{\theta}; \mathbf{x}) \right) \rightarrow 1 \hspace{5mm} \text{ as } \hspace{5mm} n \rightarrow \infty
  $$

  for any fixed $\boldsymbol{\theta} \neq \boldsymbol{\theta}^*$.
  {% endtab %}
  {% tab thrm-3-2 proof %}
  Fix $\boldsymbol{\theta} \in \boldsymbol{\Theta}$. Note that:

  $$
  \begin{aligned}
  &\mathcal{L}(\boldsymbol{\theta}^*; \mathbf{x}) > \mathcal{L}(\boldsymbol{\theta}; \mathbf{x})  \\
  \implies &\log(f(\mathbf{x}; \boldsymbol{\theta}^*)) > \log(f(\mathbf{x}; \boldsymbol{\theta})) \\ 
  \implies &\log(f(\mathbf{x}; \boldsymbol{\theta})) - \log(f(\mathbf{x}; \boldsymbol{\theta}^*)) < 0 \\
  \implies &\log\left(\prod_{i = 1}^n f(x_i; \boldsymbol{\theta}) \right) - \log\left( \prod_{i = 1}^n f(x_i; \boldsymbol{\theta}^*) \right) < 0 \\
  \implies &\sum_{i = 1}^n \left[ \log\left( \frac{f(x_i; \boldsymbol{\theta})}{f(x_i; \boldsymbol{\theta}^*)} \right) \right] < 0 \\
  \implies &\frac{1}{n} \sum_{i = 1}^n \log\left( \frac{f(x_i; \boldsymbol{\theta})}{f(x_i; \boldsymbol{\theta}^*)} \right) < 0
  \end{aligned}
  $$

  By the <a href="https://en.wikipedia.org/wiki/Law_of_large_numbers">law of large numbers</a>:

  $$
  \frac{1}{n} \sum_{i = 1}^n \log\left( \frac{f(X_i; \boldsymbol{\theta})}{f(X_i; \boldsymbol{\theta}^*)} \right) \overset{p}{\rightarrow} \mathbb{E}_{\boldsymbol{\theta}^*}\left[ \log\left( \frac{f(X; \boldsymbol{\theta})}{f(X; \boldsymbol{\theta}^*)} \right) \right]
  $$

  Since the logarithm is strictly concave for positive values, an application of <a href="https://en.wikipedia.org/wiki/Jensen%27s_inequality">Jensen's inequality</a> yields:

  $$
  \begin{aligned}
  \mathbb{E}_{\boldsymbol{\theta}^*}\left[ \log\left( \frac{f(X; \boldsymbol{\theta})}{f(X; \boldsymbol{\theta}^*)} \right) \right] 
  &< \log\left(  \mathbb{E}_{\boldsymbol{\theta}^*}\left[\frac{f(X; \boldsymbol{\theta})}{f(X; \boldsymbol{\theta}^*)} \right] \right) \\
  &= \log\left( \int  \frac{f(X; \boldsymbol{\theta})}{f(X; \boldsymbol{\theta}^*)} f(X; \boldsymbol{\theta}^*) dX \right)
  \end{aligned}
  $$

  {% endtab %}
  {% endtabs %}
</div>

The MLE is a function $\hat{\boldsymbol{\theta}}: \mathcal{X} \rightarrow \boldsymbol{\Theta}$ mapping from the sample space to the parameter space. 
MLEs exhibit the <i>invariance property</i>, which is, in words, that a function of an MLE is the MLE of that function.

<div class="theorem">
<strong>Invariance Property of Maximum Likelihood Estimators.</strong>
{% tabs invariance-prop %}
{% tab invariance-prop statement %}
Let $\hat{\boldsymbol{\theta}}$ be an MLE of $\Theta$, and let $g$ be some function of $\boldsymbol{\theta}$. Then $g(\hat{\boldsymbol{\theta}})$ is an MLE of $g(\Theta)$. 
{% endtab %}
{% tab invariance-prop proof %}
Define the <i>induced likelihood function</i>:
  
$$
\mathcal{L}^*(\eta; x) = \underset{\boldsymbol{\theta}: g(\boldsymbol{\theta}) = \eta}{\sup} \left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\}
$$

which is a function of $\eta$ equal to the maximum value of the likelihood function over all values of $\boldsymbol{\theta}$ such that $g(\boldsymbol{\theta}) = \eta$. Let:

$$
\hat{\eta} = \underset{\eta}{\arg\sup}\left\{ \mathcal{L}^*(\eta; x) \right\};
\hspace{5mm}
\hat{\boldsymbol{\theta}} = \underset{\boldsymbol{\theta}}{\arg\sup}\left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\}
$$

We have:

$$
\begin{aligned}
  \mathcal{L}^*(\hat{\eta}; x)
  &= \underset{\eta}{\sup}\left\{ \mathcal{L}^*(\eta; x) \right\} \\
  &= \underset{\eta}{\sup}\left\{ \underset{\boldsymbol{\theta}: g(\boldsymbol{\theta}) = \eta}{\sup} \left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\} \right\} \\
  &= \underset{\boldsymbol{\theta}}{\sup}\left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\}  \\
  &= \mathcal{L}(\hat{\boldsymbol{\theta}}; x) \\
  &= \underset{\boldsymbol{\theta}: g(\boldsymbol{\theta}) = g(\hat{\boldsymbol{\theta}})}{\sup} \left\{ \mathcal{L}(\boldsymbol{\theta}; x) \right\} \\
  &= \mathcal{L}^*(\hat{\boldsymbol{\theta}}; x)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Unfortunately, the MLE is not guaranteed to be unique or even to exist (see Ex. 7.1.2<d-cite key=lehmann2004></d-cite>). To make the problem a bit easier, we relax our request to just a <i>local</i> maximum of the likelihood, but we desire our estimator to have some nice properties. One of which is <i>consistency</i>. 

<div class="definition">
<strong>Definition (Consistent Estimator).</strong>
<br>
Let $X_1, X_2, \dots$ be an infinite sample where $X_i \sim P_{\boldsymbol{\theta}}$ with $\boldsymbol{\theta} \in \boldsymbol{\Theta}$, and let $\{ T_n(\boldsymbol{\theta}) \}_{n = 1}^\infty$ be a sequence of estimators of some function of $\boldsymbol{\theta}$, $g(\boldsymbol{\theta})$ ($g(\cdot)$ could be the identity function). The sequence of estimators is called <i>consistent</i> if, for all $\boldsymbol{\theta} \in \boldsymbol{\Theta}$:

$$
\underset{n \rightarrow \infty}{\lim} \mathbb{P}\left( \rvert \rvert T_n(\boldsymbol{\theta}) - g(\boldsymbol{\theta}) \rvert \rvert_2 > \epsilon \right) = \mathbf{0}_k
$$

for all $\epsilon > 0$. In other words, $T_n(\boldsymbol{\theta}) = g(\boldsymbol{\theta})+ o_p(1)$. (Technically, this is a <i>weakly</i> consistent estimator).
</div>

In order to establish the existence of a sequence of estimates that are local maxima of the likelihood, we require several conditions to be met. We begin with a couple of definitions.

