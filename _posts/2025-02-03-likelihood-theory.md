---
layout: distill
title: Likelihood and Large-Sample Theory
description: A Primer
date: 2025-02-02
tabs: true
tags: theory likelihood primer
toc:  
  - name: Background
  - name: Score and Information
    subsections:
        - name: Fisher Information Conditions
        - name: Implications
  - name: Maximum Likelihood
    subsections:
        - name: Finding MLEs
        - name: Asymptotic Normality
bibliography: 2025-02-02-likelihood-theory.bib
---

The score test in non-standard conditions has been the motivation for much of my reading these past few months. However, it has led me to wonder about the small details of the test in standard conditions. What <i>exactly</i> are the regularity conditions and when are they usually satisfied? When can we appeal to large-sample theory for the score test? It is slightly more challenging than anticipated to get a straight answer to these questions, and this is the purpose of this post. 

This post relies on some measure theory, which I've covered in <a href="/posts/2025/measure-theory">another post</a>. Most of the content comes from Moran (1971)<d-cite key=moran1971></d-cite> and Schervish (1995)<d-cite key=schervish1995></d-cite>.

---

## Background
Suppose we have some probability space $(S, \mathcal{A}, \mu)$. A random variable is some function $X: S \rightarrow \mathcal{X}$ where $\mathcal{X}$ is the sample space (which we only require to be a Borel space but is usually a subset of Euclidean space) with $\sigma$-field $\mathcal{B}$. Individual elements of $\mathcal{X}$ are denoted with $x$.

Let $$\mathcal{P}_\Theta$$ be a family of distributions for $X$ parametrized by $\Theta: S \rightarrow \Omega$ where $\Omega$ is the parameter space with $\sigma$-field $\tau$. Individual elements of $\Omega$ are denoted with $\theta$. Denote the conditional distribution of $X$ given $\Theta = \theta$ with $P_\theta$ (which is a distribution on $(\mathcal{X},\mathcal{B})$). 

To make our notation match less measure theoretical texts, we'll use $X$ to denote a random variable with realizations denoted with the lowercase $x$. A parameter will be denoted with $\Theta$ with particular values denoted by $\theta$ and its true value by $\theta^*$. The density of $X$ given parameter $\Theta$ evaluated at a particular $x$ and $\theta$ will be denoted by $f_{X \rvert \Theta}(x; \theta)$ or, more compactly, $f(x; \theta)$.


---

## Score and Information

We'll first define a <strong>very</strong> important quantity in likelihood theory: the <i>score function</i>.

<div class="definition">
<strong>Definition (Score).</strong>
<br>
Suppose $\Theta$ is $k$-dimensional, let $f_{X \rvert \Theta}(x; \theta)$ denote the density of $X$ given $\Theta = \theta$. The <i>score function</i> (or <i>score statistic</i>) is given by:

$$
\begin{equation}
\label{eq:score}
U_\Theta(x; \theta) = \frac{\partial \log f_{X \rvert \Theta}(x; \Theta)}{\partial \Theta} \bigg\rvert_{\Theta = \theta} = \frac{\partial \log f_{X \rvert \Theta}(x; \theta)}{\partial \theta}
\end{equation}
$$
</div>

The score function is the gradient of the log density of the data with respect to the parameter. It describes the curvature of the log density at a particular value of the parameter $\Theta$. 

<div class="definition">
<strong>Definition (Fisher Information).</strong>
<br>
Suppose $\Theta$ is $k$-dimensional, let $f_{X \rvert \Theta}(x; \theta)$ denote the density of $X$ given $\Theta = \theta$. The <i>Fisher information</i> is the expectation of the squared gradient of the log density of the data:

$$
\begin{equation}
\label{eq:information}
\mathcal{I}_X(\theta) = \mathbb{E}_\Theta \left[ U_\Theta(x; \theta) U_\theta(x; \theta)^\top \right]
\end{equation}
$$
</div>

The Fisher information describes the amount of information about $\Theta$ held by $X$. It is also the variance of the score function under conditions when the score has expectation zero.


### Fisher Information Conditions

Schervish outlines several regularity conditions, which he terms the <i>Fisher Information (FI) conditions</i>, that are needed for the definition of the Fisher information and some nice results about the properties of the score.

#### Condition 1

> There exists a subset of the sample space, $B$, with measure $0$ (i.e. $\mu(B) = 0$) such that $\frac{\partial f_{X \rvert \Theta}(x; \theta)}{\partial \theta_i}$ exists for any $x \notin B$ and all values (and coordinates) of $\theta$.

Condition 1 requires that the partial derivatives with respect to all coordinates of $\theta$ (for all values of $\theta$) exists almost surely.

Intuitively (and a bit hand-wavily), this means that the derivatives must exist for all possible values of $\theta$ for pretty much any sample. This implies that log-likelihood functions that have cusps or points will not be differentiable at the particular value of $\theta$ where the feature occurs, implying a violation of this condition. 

<div class="example">
  <strong>Example.</strong>
  Suppose we have $n$ i.i.d. samples $x_1, \dots, x_n \sim Unif(0, \theta)$ for some value $\theta$. The log-likelihood is given by:

  $$
  \begin{aligned}
  \ell(\theta; x_1, \dots, x_n) &= \log \left( \prod_{i = 1}^n  \frac{1}{\theta} \mathbb{I}(0 \leq x_i \leq \theta) \right) \\
  &= \log \left( \frac{1}{\theta^n} \mathbb{I}(0 \leq x_i \leq \theta; \forall i) \right) \\
  &= \log \left(\frac{1}{\theta^n}\mathbb{I}(\min_i x_i \geq 0) \mathbb{I}(\max_i x_i \leq \theta) \right)
  \end{aligned}
  $$

  The MLE is $\max_i x_i$. Clearly, the above is not differentiable at this point.
</div>

#### Condition 2
> The order of integration and differentiation can be exchanged for all coordinates of $\theta$. That is:
> $$
> \frac{\partial}{\partial \theta_i} \int f_{X \rvert \Theta}(x; \theta) d\mu(x) = \int \frac{\partial f_X(x; \theta)}{\partial \theta_i} d \mu(x)
> $$

Condition 2 states the order of integration of differentiation can be exchanged. Since differentiation is basically just a particular limit, we can use results about the interchanging of the integral and limit to get results about interchaing the integral with differentiation.

<div class="theorem">
  <strong>Theorem 1 (Dominated Convergence Theorem).</strong>
  {% tabs dom-conv %}
  {% tab dom-conv statement %}
  For a sequence of measurable functions $\{f_n\}_{n = 1}^\infty$ and measurable functions $f$ and $g$ satisfying $f_n(x) \rightarrow f(x)$ almost everywhere, $\rvert f_n(x) \rvert \leq g(x)$ almost everywhere, and $\int g(x) d\mu(x) < \infty$:

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

The Dominated Convergence Theorem states that we can interchange the order of limits and integrals for (certain) functions that are always smaller than (in absolute value) some other function with finite integral. If we define a function that mimics the form of the derivative as a limit (something along the lines of $h(x) = \frac{f(x + \delta) - f(x)}{\delta}$), then we can use this theorem to get results for derivatives and integrals.

This is the basic idea of the Leibniz integral rule:

<div class="theorem">
  <strong>Theorem 2 (Leibniz Integral Rule).</strong>
  {% tabs leibniz %}
  {% tab leibniz statement %}
  Let $\Omega$ be an open subset of $\mathbb{R}$ and $(S, \mathcal{A}, \mu) =: \mathcal{M}$ be a measure space. Let $f: \Omega \times \mathcal{M} \rightarrow \mathbb{R}$ be a function that satisfies:
  <ul>
  <li>$f(x, \theta)$ is Lebesgue-integrable in $x$ for all $\theta \in \Omega$</li>
  <li>$\frac{\partial f(x, \theta)}{\partial \theta}$ exists for all $\theta \in \Omega$ and for <i>almost all</i> $x \in \mathcal{M}$</li>
  <li>There exists integrable function $g: \mathcal{M} \rightarrow \mathbb{R}$ that is integrable and satisfies $\big\rvert \frac{\partial f(x, \theta)}{\partial \theta} \big\rvert \leq g(x)$ for all $\theta \in \Omega$ and almost every $x \in \mathcal{M}$</li>
  </ul>
  Then for all $\theta \in \Omega$:
  $$
  \frac{d}{d\theta} \int_\mathcal{M} f(x, \theta) dx = \int_{\mathcal{M}} \frac{\partial}{\partial \theta} f(x, \theta) dx
  $$
  {% endtab %}
  {% tab leibniz proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

In summary, if our log-likelihood/density satisfies the (Lebesgue version of the) Leibniz Rule conditions, then it will satisfy Condition 2. 

#### Condition 3
> The set $$C = \{ x: f_X(x \rvert \theta) > 0 \}$$ is the same $\forall \theta$.

Condition 3 states that the support of $f_X(x; \theta)$ should not depend on $\theta$. This is fairly easy to verify because we usually assume we know the family of distributions that our data are drawn from. I won't go into any more details than this.

### Implications
The FI conditions allow us to obtain the following results that are pretty fundamental for later likelihood theory. 

<div class="theorem">
<strong>Claim (Score Expectation).</strong>
{% tabs score-exp-0 %}
{% tab score-exp-0 statement %}
When these conditions hold, the score has expectation $0$. 
{% endtab %}
{% tab score-exp-0 proof %}
$$
\begin{aligned}
  \mathbb{E}_{\Theta}\left[ U_\Theta(x; \theta) \right]
  &= \int f_{X \rvert \Theta}(x; \theta) \frac{\partial \log f_{X \rvert \Theta}(x; \theta)}{\partial \Theta} d x \\
  &= \int \frac{\partial f_{X \rvert \Theta}(x; \theta)}{\partial \Theta} d x \\
  &= \frac{\partial}{\partial \Theta} \int f_{X \rvert \Theta}(x; \theta) d x \\
  &= \frac{\partial}{\partial \Theta} [1]\\
  &= 0 
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Furthermore, we have the following result under additional constraints.

<div class="theorem">
<strong>Claim (Fisher Information).</strong>
{% tabs fisher-info-1 %}
{% tab fisher-info-1 statement %}
If the log-likelihood is twice differentiable with respect to $\theta$, then the Fisher information is equal to:

$$
\mathcal{I}_X(\theta) = - \mathbb{E}_\Theta \left[ \frac{\partial^2 \log f_{X \rvert \Theta}(x; \Theta)}{\partial \Theta \partial \Theta^\top} \bigg\rvert \theta \right]
$$
{% endtab %}
{% tab fisher-info-1 proof %}
$$
\begin{aligned}
\frac{\partial^2}{\partial \Theta \partial \Theta^\top} \left[ \log f_{X \rvert \Theta}(x; \Theta)\right] 
&= \frac{\partial}{\partial \Theta} \left[ \frac{1}{f_{X \rvert \Theta}(x; \Theta)} \frac{\partial}{\partial \Theta^\top} \left[ f_{X \rvert \Theta}(x; \Theta) \right] \right] \\
&= \frac{\partial}{\partial \Theta} \left[ \frac{1}{f_{X\rvert \Theta}(x; \Theta)} \right] \frac{\partial}{\partial \Theta^\top} \left[ f_{X \rvert \Theta}(x; \Theta) \right] + \frac{1}{f_{X\rvert \Theta}(x; \Theta)} \frac{\partial^2}{\partial \Theta \partial \Theta^\top} \left[ f_{X \rvert \Theta}(x; \Theta) \right] \\
&= - \frac{\frac{\partial}{\partial \Theta} \left[ f_{X \rvert \Theta}(x; \Theta) \right] \frac{\partial}{\partial \Theta^\top} \left[ f_{X \rvert \Theta}(x; \Theta) \right]}{f^2_{X \rvert \Theta}(x; \Theta)} + \frac{\frac{\partial^2}{\partial \Theta \partial \Theta^\top}[ f_{X \rvert \Theta}(x; \Theta)]}{f_{X \rvert \Theta}(x; \Theta)} \\
&= \frac{\frac{\partial^2}{\partial \Theta \partial \Theta^\top}[ f_{X \rvert \Theta}(x; \Theta)]}{f_{X \rvert \Theta}(x; \Theta)} - \left(\frac{\partial}{\partial \Theta} [ \log f_{X \rvert \Theta}(x; \Theta)]\right)^2 
\end{aligned}
$$

Taking the expected value:

$$
\begin{aligned}
\mathbb{E}_\Theta \left[ \frac{\partial^2}{\partial \Theta \partial \Theta^\top} \left[ \log f_{X \rvert \Theta}(x; \Theta)\right] \bigg\rvert \theta \right]
&= \mathbb{E}_\Theta \left[ \frac{\frac{\partial^2}{\partial \theta \partial \theta^\top}[ f_{X \rvert \Theta}(x; \theta)]}{f_{X \rvert \Theta}(x; \theta)} - \left(\frac{\partial}{\partial \theta} [ \log f_{X \rvert \Theta}(x; \theta)]\right)^2  \right] \\
&= \mathbb{E}_\Theta \left[ \frac{\frac{\partial^2}{\partial \theta \partial \theta^\top}[ f_{X \rvert \Theta}(x; \theta)]}{f_{X \rvert \Theta}(x; \theta)} \right] - \mathbb{E}_\Theta \left[ \frac{\partial}{\partial \theta} [ \log f_{X \rvert \Theta}(x; \theta)] \right] \\
&= \mathbb{E}_\Theta \left[ \frac{\frac{\partial^2}{\partial \theta \partial \theta^\top}[ f_{X \rvert \Theta}(x; \theta)]}{f_{X \rvert \Theta}(x; \theta)} \right] - \mathcal{I}_X(\theta) \\
&= \int_\mathbb{R} \left( \frac{\frac{\partial^2}{\partial \theta \partial \theta^\top}[ f_{X \rvert \Theta}(x; \theta)]}{f_{X \rvert \Theta}(x; \theta)} \right) f(x; \theta) dx - \mathcal{I}_X(\theta) \\
&= \int_\mathbb{R} \left( \frac{\partial^2}{\partial \theta \partial \theta^\top}[ f_{X \rvert \Theta}(x; \theta)]  \right) dx - \mathcal{I}_X(\theta) \\
&\overset{(i)}{=} \frac{\partial^2}{\partial \theta \partial \theta^\top} \left[ \underbrace{\int_\mathbb{R} f_{X \rvert \Theta}(x; \theta) dx}_{= 1} \right] - \mathcal{I}_X \\
&= - \mathcal{I}_X 
\end{aligned}
$$

In $(i)$, we rely on the regularity conditions (specifically number 2 above) so we can interchange the order of differentiation and integration.
{% endtab %}
{% endtabs %}
</div>

---

## Maximum Likelihood
If we consider the observations $x$ as fixed, then we can define the <i>likelihood function</i> as a function of $\Theta$:

$$
\begin{equation}
\label{eq:lik-func}
\mathcal{L}(\theta; x) = f_{X \rvert \Theta}(x; \theta)
\end{equation}
$$

One of the most common settings in which the likelihood function will be useful is in statistical inference. A good starting point is in point estimation. Intuitively, it seems reasonable to judge the quality of a parameter estimate by how probable it is one would observe the sample at hand under the assumption that the estimate is the true parameter value. Or, in another way, we might think that the best estimate we could come up with is the one that is most likely to result in the observations we have. Thus, maximum likelihood estimation is born. 

<div class="definition">
  <strong>Definition (Maximum Likelihood Estimator).</strong>
  <br>
  Let $X$ be a random variable with density $f_{X \rvert \Theta}(x; \theta)$, and let $x$ be some realization of $X$. A <i>maximum likelihood estimator (MLE)</i> is any random quantity:

  $$
  \hat{\theta} = \underset{\theta \in \Omega}{\arg\max}\left\{ f_{X \rvert \Theta}(x; \theta) \right\} 
               = \underset{\theta \in \Omega}{\arg \max}\left\{ \mathcal{L}(\theta; x) \right\}
  $$

  The MLE is a function $\hat{\theta}: \mathcal{X} \rightarrow \Omega$ mapping from the sample space to the parameter space. 
</div>

If the parameter space $\Omega$ is compact and the likelihood function is continuous over $\Omega$, then maximum likelihood estimate will exist for a given sample (i.e. the supremum of the maximum likelihood estimator will be achieved in $\Omega$). If the parameter space is open, then the likelihood may increase and never reach a supremum. 

<div class="example">
  <strong>Example.</strong>
  Consider $x$ distributed uniformly on the open interval $(0, \theta)$. The likelihood is:

  $$
  \mathcal{L}(\theta; x) = \frac{1}{\theta} \mathbb{I}(x > 0) \mathbb{I}(x < \theta)
  $$ 

  This function is decreasing on the interval $(0, \theta)$, and the maximum is never achieved. 
</div>

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

### Finding MLEs
The easiest way to find an MLE is to use set the log-likelihood equal to zero (since monotonic transformations will not affect the $\arg \max$ or $\arg \min$). 

<div class="theorem">
  <strong>Theorem 3 (Fermat's Interior Extremem Theorem).</strong>
  {% tabs fermat %}
  {% tab fermat statement %}
  For function $f: A \rightarrow \mathbb{R}$, let $x_0 \in A$ be a local extremum of $f$. If $f$ is differentiable at $x_0$, then $f'(x_0) = 0$.

  A simple corollary states that a global extremum $x_0$ of $f$ must fall into one of the following cases: (1) $x_0$ is on the boundary of $A$; (2) $f$ is not differentiable at $x_0$; (3) $x_0$ is a stationary point of $f$.
  {% endtab %}
  {% tab fermat proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

For Fermat's Theorem to apply, the log-likelihood must be differentiable in $\Omega$. The previous example is one such where this is not the case. This theorem also implies that this $x_0$ occurs on the _interior_ of the domain of $f$. We would have to also check the boundary points if $A$ were closed.

### Asymptotic Normality 
The MLE is asymptotically normal under suitable conditions.

<div class="theorem">
  <strong>Theorem 4 (Asymptotic Normality of MLEs).</strong>
  {% tabs asymp-norm-mle %}
  {% tab asymp-norm-mle statement %}
  Let $\Omega \subseteq \mathbb{R}^p$ and $\{ X_n \}_{i = 1}^\infty$ be conditionally i.i.d. random variables given $\Theta = \theta^*$ with distribution $P_{\theta^*}$. Let $\hat{\theta}_n$ be an MLE for $\Theta$ and assume $\hat{\theta}_n \overset{p}{\rightarrow} \theta^*$. 

  Further assume that the second partial derivatives with respect to $\Theta$ of the densities are continuous and that the order of differentiation and integration can be exchanged. Suppose there exists function $H_r(x, \theta)$ such that, for $\theta^*$ in the interior of $\Omega$ and each $k,j$, the following is satisfied with $\underset{r \rightarrow 0}{\lim} \mathbb{E}_{\theta^*} \left[ H_r(X, \theta^*) \right] = 0$:

  $$
  \underset{\rvert\rvert \theta - \theta^* \rvert \rvert \leq r}{\sup} \left\{  \frac{\partial^2 \ell(\theta^*; x)}{\partial \theta_k \partial \theta_j} - \frac{\partial^2 \ell(\theta; x)}{\partial \theta_k \partial \theta_j} \right\} \leq H_r(x, \theta^*)
  $$

  And finally, assume the Fisher information for a single data point, $\mathcal{I}_X(\theta)$, is finite and non-singular. Then, under $P_{\theta^*}$:

  $$
  \sqrt{n}\left( \hat{\theta}_n - \theta^* \right) \rightsquigarrow \mathcal{N}\left(0, \mathcal{I}^{-1}_{X}(\theta^*) \right)
  $$

  In words, this theorem states that the MLE (suitably centered and scaled) is asymptotically normal.
  {% endtab %}
  {% tab asymp-norm-mle %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

The requirements are that the the true parameter value is on the interior of the parameter space (if it is restricted); the MLE is consistent; the density is nice enough (continuous second derivatives); the order of integration and differentiation can be exchanged; and that there is a function with finite mean that bounds the difference between the second derivatives of the log-likelihoods for two values of $\Theta$. This last condition is a uniform law of large numbers.

<!-- FINISH PROOF?
<details>
<summary>Proof.</summary>
Since we have an i.i.d. sample, we can write the likelihood (scaled by $\frac{1}{n}$) as:

$$
\ell(\theta) := \ell_\Theta(\theta; x_1, \dots, x_n) = \frac{1}{n} \sum_{i = 1}^n \log f_{X_i \rvert \Theta}(x_i; \theta)
\nonumber
$$

Since $\theta^*$ is in the interior of the parameter space, there is an open neighborhood centered at $\theta^*$ that is contained in the parameter space. By the consistency assumption on $\hat{\theta}_n$, we have that:

$$
Z_n I_{int(\Omega)^c}(\hat{\theta}_n) = o_P(\frac{1}{\sqrt{n}})
\nonumber
$$
</details> -->

We often supply $$\hat{\theta}_n$$ for (usually unknown) $$\theta^*$$ in the Fisher information, which yields the _expected Fisher information_. The expected Fisher information converges in probability to the Fisher Information given $$\Theta = \theta^*$$.

We could also instead use the scaled matrix of second partial derivaties of the log-likelihood function near to $\hat{\theta}_n$:

$$
- \frac{1}{n} \left\{ \frac{\partial^2 \ell(\Theta; x)}{\partial \Theta_i \Theta_j} \bigg\rvert_{\Theta = \hat{\theta}_n} \right\}
$$

which is called the <i>observed Fisher information</i>. It's basically the sample analog of the expected information.

