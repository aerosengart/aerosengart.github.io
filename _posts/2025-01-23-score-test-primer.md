---
layout: post
title:  "A Score Test Primer"
date: 23 January 2025
categories: posts
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

Recently, I've been doing a lot of exploring in the literature related to score tests and variance component testing. A few questions have come up regarding how a score test can work when the true value of the parameter under scrutiny is on the boundary of the parameter space. 

Let's first take a quick detour into the score test in general before diving into its use for variance component testing in mixed models. For sake of space, I'll use the notation $$\mathbb{E}_{\theta^*}\left[ \cdot \right]$$ to denote the expectation taken when the true value of $\theta$ is $\theta^*$.

## Score Test

#### Score

Suppose we have random vector $Y$ with probability density function $f_Y(y; \theta)$. Let $\ell(\theta; y) = \log(f_Y(y; \theta))$ denote its log-likelihood. Assuming that $f_Y(y; \theta)$ is partially differentiable for $y$ and all components of $\theta$, the _(efficient) score statistic_ (since it is a function of data $y$) for $\theta$, which we denote with $$U_\theta(y)$$, is given by the gradient of the log-likelihood with respect to $\theta$:

$$
U_\theta(y) = \frac{\partial \ell(\theta; y)}{\partial \theta} 
\nonumber
$$

The score describes the curvature of the log-likelihood at a particular value of a parameter. 

Following Cox & Hinkley (1974), suppose we are interested in testing the simple null hypothesis $H_0: \theta = \theta_0$ against the two-sided alternative $H_1: \theta \neq \theta_0$. If we assume that the order of integration and differentiation can be exchanged (Leibniz rule), then we have the following result (where $\theta_0$ is the true value of $\theta$):

$$
\mathbb{E}_{\theta_0} \left[ U_\theta(y) \right] = 0
\label{eq:score-zero-mean}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\mathbb{E}\left[ U_\theta(y; \theta_0) \bigg \rvert \theta_0 \right] &= \mathbb{E}\left[ \frac{\partial \ell(\theta; y)}{\partial \theta}\bigg\rvert_{\theta = \theta_0} \right] \\
&= \left[ \int  \frac{\partial}{\partial \theta} \left[ \log\left( f_Y(y; \theta) \right)  \right]  f_Y(y; \theta) dy \right] \bigg\rvert_{\theta = \theta_0} \\
&= \left[ \int  \frac{\partial f_Y(y; \theta) }{\partial \theta}\frac{1}{f_Y(y; \theta)} f_Y(y; \theta) dy \right] \bigg\rvert_{\theta = \theta_0} \\
&= \left[ \frac{\partial }{\partial \theta} \int   f_Y(y; \theta)  dy \right] \bigg\rvert_{\theta = \theta_0} \\
&\overset{(i)}{=} \frac{\partial }{\partial \theta} 1 \\
&= 0
\end{aligned}
\nonumber
$$

In $(i)$, we use the fact that densities integrate to $1$.
</details>

If the data are independent, then the total log-likelihood is the sum of the log-likelihoods for each observation. Thus, the score will also be the sum of the individual score statistics for each observation.


#### Fisher Information
The Fisher information is the amount of information about the parameter $\theta$ that is held by the data, $y$. It is found as the expectation of the squared derivative of the log-likelihood:

$$
\mathcal{I}_\theta = \mathbb{E}_{\theta^*} \left[ \left( \frac{\partial}{\partial \theta} \ell(\theta; y) \right)^2 \right]  \overset{(i)}{=} \text{Var}_{\theta^*} \left[ U_\theta(y) \right]  \overset{(ii)}{=} -\mathbb{E}_{\theta^*} \left[ \frac{\partial^2}{\partial \theta^2} \ell(y; \theta) \right]
\label{eq:fisher-information}
$$

where $(i)$ and $(ii)$ hold under certain regularity conditions (see below).

<details>
<summary>Proof.</summary>
In the following, I omit the subscript $\theta^*$ and just use $\theta$ to denote the true value.
For $(i)$:
$$
\begin{aligned}
\text{Var} \left[ U_\theta(y) \right] &= \mathbb{E}\left[ U_\theta^2(y) \right] - \underbrace{\mathbb{E}\left[ U_\theta(y) \right]^2}_{= 0} = \mathbb{E}\left[ U_\theta^2(y) \right]
\end{aligned}
\nonumber
$$
<br>
For $(ii)$:
$$
\begin{aligned}
\frac{\partial^2 \ell(\theta; y)}{\partial \theta^2} &= \frac{\partial}{\partial \theta} \left[ \frac{\frac{\partial f(y; \theta)}{\partial \theta} }{f(y; \theta)} \right] = \frac{\frac{\partial^2 f(y; \theta)}{\partial \theta^2}}{f(y; \theta)} - \left(\frac{\frac{\partial f(y; \theta)}{\partial \theta}}{f(y; \theta)}\right)^2 = \frac{\frac{\partial^2 f(y; \theta)}{\partial \theta^2}}{f(y; \theta)} - \left(\frac{\partial \ell(\theta; y)}{\partial \theta}\right)^2 \\
\implies - \mathbb{E}\left[ \frac{\partial^2 \ell(\theta; y)}{\partial \theta^2} \right] &= \mathbb{E}\left[ -\frac{\frac{\partial^2 f(y; \theta)}{\partial \theta^2}}{f(y; \theta)} + \left(\frac{\partial \ell(\theta; y)}{\partial \theta}\right)^2\right] = - \int \frac{\partial^2 f(y; \theta)}{\partial \theta^2} dy + \int \left( \frac{\partial \ell(\theta; y)}{\partial \theta}\right)^2 dy \overset{(a)}{=} \int \left( \frac{\partial \ell(\theta; y)}{\partial \theta}\right)^2 dy = \mathbb{E}\left[ U_\theta^2(y) \right]
\end{aligned}
\nonumber
$$
where $(a)$ follows from:
$$
\int \left( \frac{\partial f(\theta; y)}{\partial \theta}\right)^2 dy  = \frac{\partial^2}{\partial \theta^2} \int f(\theta; y) dy = 0
\nonumber
$$
</details>

The Fisher Information will be non-negative in the scalar case and positive semidefinite/definite in the vector case (unless the scores for each component of the parameter vector are linearly dependent).

#### Conditions
We assume some regularity conditions hold so that the score and Fisher Information have several desirable properties. In order for Eq. \eqref{eq:score-zero-mean} and and for $(i)$ in Eq. \eqref{eq:fisher-information} to hold, we need the following to be met, as outlined in Schervish (1995):

1) The parameter space is an open interval

2) The support of the density of $y$ is independent of the parameter $\theta$

3) The first partial derivatives exist and are finite for all $y$ and $\theta$ 

4) The order of integration and differentiation can be exchanged

For $(ii)$ to hold, we additionally need the log-likelihood to be twice differentiable for all $y$ and $\theta$. In addition, in the case of vector-valued parameters, (3) must hold with respect to each component of $\theta$. 

If these conditions are not met, then several problems can occur including the score vanishing identically (leading to the Fisher Information also vanishing) or the score being undefined at the null value. 

#### Test Statistic

In the case of scalar-valued parameters, the score test statistic is given by:

$$
S_\theta(y; \theta_0) = \frac{U^2_\theta(y)}{\mathcal{I}_\theta} \bigg\rvert_{\theta = \theta_0}
\nonumber
$$

For vectors, it can be written in matrix form:

$$
S_\theta(y; \theta_0) = U^\top_\theta(y) \mathcal{I}_\theta^{-1} U_\theta(y)\bigg\rvert_{\theta = \theta_0}
\nonumber
$$

<details>
<summary>Derivation Of Score Test As A Locally Most Powerful Test.</summary>
For simplicity, assume we have a one-dimensional parameter of interest, though the results are extendable to the multi-dimensional case. Let $H_0: \theta = \theta_0$ and $H_A: \theta > \theta_0$, which is a composite alternative. We want to find a locally most powerful test, so we will consider the specific alternative $H_A: \theta = \theta_0 + \delta$ for some small $\delta$. 
<br>
Recall that by the Neyman-Pearson Lemma, the likelihood ratio is a uniformly most powerful test for a simple null against a simple alternative. Thus, to find a locally most powerful test, we can use the likelihood ratio corresponding to our specific $H_A$. 
<br>
Let $f_Y(y; \theta)$ denote the density of $Y$ and let $\ell(\theta; y)$ denote the log-likelihood function of $\theta$ given $y$. We can write the log-likelihood of $\theta_0 + \delta$ using a Taylor approximation (for simplicity, we'll just do first order): 

$$
\ell(\theta_0 + \delta; y) = \ell(\theta_0; y) + \delta \frac{\partial \ell(\theta; y)}{\partial \theta} \bigg\rvert_{\theta = \theta_0} + R
\nonumber
$$

where the remainder, $R = o(\delta)$.

<details>
<summary>Proof.</summary>
The remainder can be written as:

$$
\bigg \rvert \ell(\theta_0 + \delta; y) -  \left(\ell(\theta_0; y) + \delta \frac{\partial \ell(\theta; y)}{\partial \theta}\bigg\rvert_{\theta = \theta_0} \right)  \bigg \rvert 
= \bigg \rvert \ell(\theta_0 + \delta; y) -  \ell(\theta_0; y)  -  \delta\frac{\partial \ell(\theta; y)}{\partial \theta}\bigg\rvert_{\theta = \theta_0} \bigg\rvert
= \bigg \rvert \log\left( \frac{f_Y(y; \theta_0 + \delta)}{f_Y(y; \theta_0)}\right) -  \delta\frac{\partial \ell(\theta; y)}{\partial \theta}\bigg\rvert_{\theta = \theta_0} \bigg\rvert
\nonumber
$$

As $\delta \rightarrow 0$, the RHS of the above will go to 0, which implies it is $o(\delta)$.
</details>

The above equation implies:

$$
\ell(\theta_0 + \delta; y) - \ell(\theta_0; y) = \delta \frac{\partial \ell(\theta; y)}{\partial \theta} \bigg\rvert_{\theta = \theta_0} + o(\delta) \hspace{3mm} \implies \hspace{3mm} \log\left( \frac{f_Y(y; \theta_0 + \delta)}{f_Y(y; \theta_0)} \right) = \delta \frac{\partial \ell(\theta; y)}{\partial \theta} \bigg\rvert_{\theta = \theta_0} + o(\delta)
\nonumber
$$

Thus, for very small $\delta$, the log-likelihood ratio $\log LR_{A; 0}(y) = \frac{\ell(\theta_0 + \delta; y)}{\ell(\theta_0; y)}$ will be achieved for large values of:

$$
U_\theta(y; \theta_0) = \frac{\partial \ell(\theta; y)}{\partial \theta}\bigg\rvert_{\theta = \theta_0}
\nonumber
$$

 Intuitively, we can obtain the critical region associated with the likelihood ratio ($c_{LR} := \{ y \big\rvert LR_{A; 0}(y) > t \}$ for some $t$) through the log-likelihood ratio, since the logarithm is a monotonically increasing function. As $\delta$ shrinks, large values of $LR_{A; 0}(y)$ will occur only when the efficient score is also suitably large (since our approximation error shrinks as well); thus, we can recover the critical region associated with the LR using $U_\theta(y; \theta_0)$. 
<br>
<br>
If the data are i.i.d., then we can decompose the total log-likelihood into the sum of log-likelihoods for each data point, thereby decomposing the score into a sum of independent components. This means we can appeal to central limit theory and use a Gaussian approximation for the distribution of $U_\theta$ under the null (the mean is $0$, and the variance is the total Fisher information). The efficient score is therefore a locally most powerful test for alternatives of the form $\theta = \theta_0 + \delta$, where significance is determined by larger values for $\delta > 0$ and smaller values for $\delta < 0$. 
</details>

The score test has the benefit of only requiring a model fit under the null hypothesis, which can greatly simplify the computational requirements of calculating the test statistic. 

<details>
<summary>A Note On The Score Test.</summary>
In econometrics, the score test is often called the Lagrange multiplier test due to its relationshiop with the use of Lagrange multipliers for optimization.
<br>
Our null hypothesis is that $H_0: \theta = \theta_0$. If this is true, then we would expect the derivative of the log-likelihood to be close to $0$ because we could expect $\ell(y; \theta_0)$ to be close to that of the (unrestricted) MLE of $\theta$. 
<br>
Define $h(\theta) = \theta - \theta_0$. Wee can reformulate the null hypothesis as $H_0: h(\theta) = 0$. If the null hypothesis is true, then we would expect $\theta_0$ to be close to a solution to the constrained maximization of $\ell(y; \theta)$ subject to $h(\theta) = 0$. The Lagrangian function for this case is given by $L(y; \theta) = \ell(y; \theta) + \lambda h(\theta)$, and a solution will occur at a stationary point. In other words, a solution to $\frac{d L(y; \theta)}{d \theta} = 0$ will coincide with a solution, $(\theta^*, \lambda^*)$, to the following system:

$$
\begin{array}{l}
\frac{d \ell(y; \theta^*)}{d \theta^*} + \lambda^* \frac{d h(\theta^*)}{d \theta^*} = 0 \\
h(\theta^*) = 0
\end{array}
\nonumber
$$

The Lagrange multiplier, $\lambda$, can be thought of as the amount of change we would see in the solution when the constraint value changes. In our context, a large $\lambda$ implies that $\theta_0$ is not the true value of $\theta$. 
<br>
A test statistic can be formed using the above as:
$$
S_\theta(y; \theta_0) = \frac{\left(\frac{d \ell(y; \theta)}{d \theta}\right)^2}{\mathcal{I}_\theta}\bigg\rvert_{\theta = \theta_0} = \frac{\left(\lambda \frac{d h(\theta)}{d \theta}\right)^2}{\mathcal{I}_\theta} \bigg\rvert_{\theta = \theta_0}
\nonumber
$$

Somewhat loosely, we've squared the score because we do not care if the deviations are negative or positive, and the scaling essentially adjusts for the average curvature of the log-likelihood for this dataset. The latter is needed because we want the test statistic to be larger for a $\theta_0$ that is farther from the maximizer of the log-likelihood. In the case where two log-likelihoods have the same slope at $\theta_0$, $\theta_0$ will closer to the maximizer for the log-likelihood that has a larger second derivative (i.e. the "curvier" log-likelihood) at $\theta_0$.
</details>

For situations satisfying our regularity conditions (and a few others — see "Some Asymptotic Theory" below), we can derive the limiting distribution of the score test statistic for simple null hypotheses as chi-squared.


<details>
<summary>Some Asymptotic Theory.</summary>
Conducting a hypothesis test requires knowing the distribution of the test statistic under the null hypothesis. This can sometimes be derived using standard theoretical results so long as the sample size grows large and we satisfy some regularity conditions. Cox & Hinkley (1997) list the following conditions on the density $f_Y(y; \theta)$:
<br>
1) The parameter space is closed, compact, and has finite dimension. The true parameter value must also lie on the interior of the space. 
<br>
2) The probability distributions defined by distinct values of $\theta$ are distinct.
<br>
3) The first, second, and third order derivatives of the log-likelihood $\ell(\theta; y)$ with respect to $\theta$ exist in a neighborhood of the true value almost surely.
<br>
4) In the neighborhood defined in 3), $\frac{1}{n}\big\rvert \frac{\partial^3 \ell(\theta; y)}{\partial \theta^3} \big\rvert \leq g(Y)$ for some function of $Y$, $g(Y)$, and $\mathbb{E}[g(Y)]$ exists.
<br>
5) The total (for all of the data) Fisher Information is finite and positive definite in a neighborhood of the true parameter value, and $(i)$ and $(ii)$ of Eq. \eqref{eq:fisher-information} hold.
<br>
<i>I am unsure if these are necessary only for their discussion in order to simplify things for the student reader or if they are requirements standard for asymptotic theory.</i>

<details>
<summary>A Quick Aside</summary>
In general, for a test statistic $T_n$, we want to show that $T_n = T + o_p(1)$ for random variable $T$ for which we know the probability distribution. Here, $o_p(1)$ is some variable $Z_n$ satisfying $\underset{n \rightarrow infty}{\lim} \mathbb{P}\left( \rvert Z_n \rvert \geq \epsilon \right) = 0$ for any $\epsilon > 0$ (i.e. $T_n$ converges in probability to $T$). 

Cox and Hinkley also use the notation $\mathbb{E}[T_n]$ to denote the mean of some limiting distribution, which is kind of hand-wavy. In this case, they mean $\mathbb{E}[T_n]$ is an approximation that converges to the mean of the limiting distribution as the sample size grows large or a truncated $T_n$ also has the same limiting distribution.
</details>

First, we must detour into maximum likelihood theory. For a general log-likelihood function, $\ell(\theta; y)$, the maximum likelihood estimate of $\theta$, which we denote with $\hat{\theta}$, will satisfy $\ell(\hat{\theta}; y) \geq \ell(\theta^*; y)$ for all other $\theta^* \in \Theta$.
<br>
Finding this value can be done by identifying the values of $\theta$ that satisfy $\frac{\partial \ell(\theta; y)}{\partial \theta} = 0$ (for all componenets of $\theta$ in the vector-valued parameter case) and where $\frac{\partial^2 \ell(\theta; y)}{\partial \theta^2} < 0$ (or negative definite in the vector-valued parameter case). Solutions will be local maxima, and one can then evaluate the likelihood at each solution, and the maximizing value can be chosen as the MLE. The likelihood may not be maximized by a unique value, and it may not occur at a stationary point (see uniform example below).
<br>
Note that, if the parameter space is compact and the likelihood is continuous, then the MLE will exists. Further, if the parameter space is convex and the likelihood is strictly concave in the parameter, then the MLE will be unique (if it exists).

<br>
<br>
When doing estimation, we want our estimator to be consistent. Cox & Hinkley do a good job of providing an intuitive explanation for this property: if we had access to as many samples as we wanted, then our estimate should be exactly correct.

Weak consistency is the most useful for our purposes. with the last being most useful for our purposes: If $T_n = \theta + o_p(1)$ (i.e. it converges in probability to $\theta$), then $T_n$ is called weakly consistent for $\theta$. In this case, $T_n$ will satisfy a weak law of large numbers.

<details>
<summary>Others Types Of Consistent.</summary>
Let $T_n$ be our estimates made with i.i.d. random variables $Y_1, \dots, Y_n$. Define $\tilde{F}_n(y) = \frac{1}{n} \sum H(y - Y_{(j)})$ as the empirical c.d.f. where $H(\cdot)$ denotes the Heaviside function. If $T_n$ can be written as some function of the empirical c.d.f. ($T_n = t(\tilde{F}_n(\cdot))$) and $t(F(\cdot, \theta)) = \theta$, then it is called Fisher consistent. Fisher consistency implies weak consistency.
<br>
If $\mathbb{P}\left( \underset{n \rightarrow \infty}{\lim} T_n = \theta\right) = 1$ ($T_n$ converges almost surely to $\theta$), then it is called strongly consistent and will satisfy the strong law of large numbers. 
</details>

Under certain conditions, the MLE is a consistent estimator of the true parameter values:
<br>
1) If $\theta \neq \theta'$, then $f(y; \theta) \neq f(y; \theta')$. 
<br>
2) The parameter space is compact. Alternatively, we could require that there is just a compact neighborhood around the true parameter value and that the likelihood is less than the maximum achieved there.
<br>
3) The log-likelihood is continuous in the parameter for almost all $y$.
<br>
4) There is some $g(y)$ that is integrable with respect to $f(y; \theta^*)$ for true value $\theta^*$ such that $\rvert \ell(\theta; y) \rvert < g(y)$ for all $\theta$ in the parameter space.
<br>
(3) and (4) together implies that a uniform law of large numbers holds, namely: $\underset{\theta \in \Theta}{\sup} \rvert \ell(\theta; y) - \mathbb{E}_{\theta^*} \left[ \ell(\theta; y) \right] \rvert \overset{p}{\rightarrow} 0$. 

<br>
Suppose the following hold:
<br>
1) The dimension of the parameter space does not change with $n$
<br>
2) The density $f(y; \theta)$ is third order differentiable
<br>
3) Integration and differentiable can be exchanged
<br>
4) If $\theta \neq \theta'$, then $f(y; \theta) \neq f(y; \theta')$
<br>
5) $\theta$ is in the interior of $\Theta$
<br>
Then $\sqrt{n}(\hat{\theta} - \theta^*) \rightsquigarrow \mathcal{N}(0, \mathcal{I}_\theta^{-1})$. 
<br>
The proof is a bit involved, so I will omit the details, but there are plenty of references that cover it (see Wasserman). Of note is the fact that the proof relies upon a Taylor series expansion at $\hat{\theta}$ about the true value $\theta^*$. This is where the assumption that $\theta$ is on the interior of the parameter space comes in. Without it, we cannot write the Taylor approximation! Moran (1971) covers some other results in these non-standard situations.
<br>
Something else of importance is that the Fisher information may depend on the true value of $\theta$, which is often unknown. If we use a consistent estimate of it, then the distributional results should still hold. Examples include plugging in $\hat{\theta}$ for $\theta$ or taking the negative Hessian and evaluating that at $\hat{\theta}$.
</details>


<div class="example">
  <body>
  <strong>Example.</strong>
  Suppose we have $n$ observations that are i.i.d. Gaussian with known variance $\sigma^2$. Denote the collection of data points are $y$, where a single point is $y_i$ for $i = 1:n$. 
  <br>
  We want to test $H_0: \mu = 0$ against $H_A: \mu \neq 0$. Their log-likelihood is:

  $$
  \ell(\mu, \sigma; y) = -\frac{n\log(2 \pi \sigma^2)}{2}- \frac{1}{2\sigma^2}\sum_{i = 1}^n (y_i - \mu)^2
  \nonumber
  $$

  It follows that:

  $$
  \begin{aligned}
  U_\mu(y) &= \frac{1}{\sigma^2} \sum_{i = 1}^n (y_i - \mu) \overset{\mu = 0}{=} \frac{1}{\sigma^2} \sum_{i = 1}^n y_i \\
  \mathcal{I}_\mu &= \mathbb{E}\left[ \left( \frac{1}{\sigma^2} \sum_{i = 1}^n (y_i - \mu) \right)^2\right] = \frac{1}{\sigma^4} \sum_{i = 1}^n\left( \mathbb{E}\left[ (y_i - \mu)^2 \right] + 2\sum_{i' < i} \mathbb{E}\left[ y_i - \mu \right]\mathbb{E}\left[ y_{i'} - \mu \right] \right) = \frac{n}{\sigma^2}
  \end{aligned}
  \nonumber
  $$

  The score test statistic is:

  $$
  S_\mu(y; 0) = \frac{\left( \frac{1}{\sigma^2} \sum_{i = 1}^n y_i \right)^2}{\frac{n}{\sigma^2} } = \frac{\left( \sum_{i = 1}^n y_i \right)^2}{\sigma^2 n}
  \nonumber
  $$

  This should be asymptotically $\chi^2_1$ since we satisfy the regularity conditions.
  </body>
</div>


## Variance Component Testing

Variance component testing for linear mixed models, as explored in some prior posts of mine and in many publications from the last century, can be approached from two perspectives that are outlined by Verbeke and Molenberghs in their 2003 paper "The Use of Score Tests for Inference on Variance Components". 

Suppose we have $i = 1:n_j$ observations grouped into $j = 1:k$ clusters, each of which has $n_j$ observations. We assume that we have fixed effects, $\alpha \in \mathbb{R}^{p}$, with the fixed effects covariates for observation $i$ in cluster $j$ denoted by the $p$-dimensional vector $x_{i,j}$. We have random intercepts, $\beta_i \in \mathbb{R}$, that are assumed to be $\beta_i \overset{iid}{\sim} \mathcal{N}(0, \tau^2)$ for variance component $\tau^2$. The error terms are assumed to be $\epsilon_{i,j} \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$ and independent of the random effects. This yields the mixed model:

$$
y_{i,j} \overset{ind}{\sim} \mathcal{N}(x_{i,j}^\top \alpha + \beta_i, \sigma^2)
\label{eq:condition-model}
$$

#### Approach 1

In Approach 1, we are not so interested in the random effects themselves but want to include them to have some additional flexibility in our model. Instead, we want a marginal model (rather than the conditonal model above), which means we integrate out the random effects:

$$
y_{\cdot, j} \overset{ind}{\sim} \mathcal{N}(x_{\cdot, j}\alpha, \tau^2 \mathbf{1}_{n_j} \mathbf{1}_{n_j}^\top + \sigma^2 \mathbb{I}_{n_j})
\label{eq:marg-model}
$$

where $y_{\cdot, j}$ denotes the $n_j$-dimensional vector of reponses for cluster $j$; $x_{\cdot, j}$ denotes the $n_j \times p$ matrix of fixed effects covariates for cluster $j$; and the $n_j$ subscript denote the dimensionality of a vector of ones ($\mathbf{1}$) and the identity matrix ($\mathbb{I}$). 


In this scenario, it is reasonable to allow $\tau^2 < 0$ (as long as $\tau^2 + \sigma^2 > 0$) as it corresponds to a negative ICC (the intracluster correlation coefficient is given by $\frac{\tau^2}{\tau^2 + \sigma^2}$).

#### Approach 2

In Approach 2, we are interested in the conditional model (Eq. \eqref{eq:condition-model}), which means $\tau^2$ must be positive. In this case, classical theory for hypothesis testing will not always hold due to the assumption that the true parameter value lies in an open ball contained in the parameter space, which is not true when $\tau^2 = 0$. 

<br>

In both cases, there are some values that the parameter of interest is permitted to take but cannot be meaningfully interpreted. In Approach 1, these values are dependent upon $\sigma^2$, whereas in Approach 2, it is simply that $\tau^2 > 0$. 



---

## An Illustrative Example

To try to wrap my head around this scenario, I've gone through the calculations for the score test statistic in a (somewhat) simple Gaussian case.

Suppose we have $N = \sum_{j = 1}^k n$ Gaussian observations that can be clustered into $k$ different groups. We'll denote the $i$-th observation in cluster $j$ with $y_{i,j}$, and we'll assume that observations in different clusters are independent.

We'll assume that our observations have conditional means $\mu_{i,j} = \alpha + \beta_j$ where $\alpha$ is some overall intercept and each $\beta_j$ is an intercept associated with cluster $j$. Thus, we have:

$$
\begin{aligned}
y_{\cdot,j} &= \alpha\cdot \mathbf{1}_{n} + \beta_j \cdot \mathbf{1}_{n} + e_{\cdot,j} \\ 
e_{\cdot,j} &\sim \mathcal{N}(0 \cdot \mathbf{1}_n, \sigma^2 \cdot \mathbb{I}_{n \times n}) \\
\beta &\sim \mathcal{N}(0 \cdot \mathbf{1}_k, \tau^2 \cdot \mathbb{I}_{k \times k})
\end{aligned}
\nonumber
$$

where we assume the errors and random effects are independent. Let $$\Sigma_{e_j} = \sigma^2 \mathbb{I}_{n \times n}$$ and $$\Sigma_\beta = \tau^2 \mathbb{I}_{k \times k}$$. To save space, we'll let $$\tilde{y}:= y - \alpha \cdot \mathbf{1}_n \mathbf{1}_k^\top$$, the centered data matrix under the null.

By our independence assumptions, observations within the same group are also independent.

<details>
<summary>Proof Of Independence.</summary>
$$
\begin{aligned}
\text{Cov}(y_{i,j}, y_{i', j}) &= \mathbb{E}\left[ (y_{i,j} - \mu_{i,j})(y_{i', j} - \mu_{i',j})\right] - \\
&= \mathbb{E}\left[(\alpha + \beta_j + e_{i,j} - \alpha - \beta_j)(\alpha + \beta_j + e_{i', j}) \right] \\
&= \mathbb{E}\left[ e_{i,j} e_{i', j}\right] \\
&= 0
\end{aligned}
\nonumber
$$
</details>

#### Likelihood Derivations

The likelihood of $(\alpha, \tau^2)$ given $y$ and conditional on $\beta$ and the density of $\beta$ are:

$$
\mathcal{L}(\alpha, \tau^2; y \rvert \beta) = \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}} \exp \left( - \frac{1}{2}  (y_{\cdot, j} - \mu_{\cdot, j})^\top \Sigma_{e_j}^{-1} (y_{\cdot, j} - \mu_{\cdot, j}) \right) 
\hspace{20mm}
f(\beta) = \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \exp \left( -\frac{1}{2} \beta^\top \Sigma^{-1}_\beta \beta \right) 
\nonumber
$$

Their product is the joint likelihood of $(\alpha, \tau^2, \beta)$ given $y$ is 
<span class="popup" onclick="PopupFunc('pop1')">
    proportional to
    <span class="popuptext" id="pop1">
        The constant of proportionality is $$\frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}$$
    </span>
</span>:


$$
\mathcal{L}(\alpha, \tau^2, \beta; y) \propto \exp \left( - \frac{1}{2} \left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)^\top \left( \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k} \right)\left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right) \right) 
\nonumber
$$

<details>
<summary>Joint Likelihood Derivation.</summary>
$$
\begin{aligned}
\mathcal{L}(\alpha, \tau^2, \beta; y) &= \left(   \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}} \exp \left( - \frac{1}{2}  (y_{\cdot, j} - \mu_{\cdot, j})^\top \Sigma_{e_j}^{-1} (y_{\cdot, j} - \mu_{\cdot, j}) \right)  \right) \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \exp \left( -\frac{1}{2} \beta^\top \Sigma^{-1}_\beta \beta \right)  \right) \\
&= \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2} \left[ \beta^\top \Sigma_\beta^{-1} \beta + \sum_{j = 1}^k (y_{\cdot, j} - \mu_{\cdot, j})^\top \Sigma_{e_j} (y_{\cdot, j} - \mu_{\cdot, j}) \right] \right) \\
&\overset{(i)}{=} \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2} \left[ \beta^\top \Sigma_\beta^{-1} \beta + \sum_{j = 1}^k (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n} - \beta_j \cdot \mathbf{1}_{n})^\top \Sigma_{e_j}^{-1} (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n} - \beta_j \cdot \mathbf{1}_{n}) \right] \right) \\
&= \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2}  \left[ \beta^\top \Sigma_\beta^{-1} \beta + \sum_{j = 1}^k \left[ (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n})^\top \Sigma_{e_j}^{-1} (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n}) - 2 (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n})^\top \Sigma_{e_j}^{-1} \left( \beta_j \cdot \mathbf{1}_{n}\right)  +  \left( \beta_j \cdot \mathbf{1}_{n}\right)^\top \Sigma_{e_j}^{-1} \left( \beta_j \cdot \mathbf{1}_{n}\right) \right] \right] \right) \\
&\overset{(ii)}{=} \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2}  \left[ \beta^\top \left( \Sigma_\beta^{-1} + \frac{n}{\sigma^2} \mathbb{I}_{k \times k}\right) \beta + \sum_{j = 1}^k \left[ (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n})^\top \Sigma_{e_j}^{-1} (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n}) - 2 (y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n})^\top \Sigma_{e_j}^{-1} \left( \beta_j \cdot \mathbf{1}_{n}\right) \right] \right] \right) \\ 
&\overset{(iii)}{=} \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2}  \left[ \beta^\top \left( \left[ \frac{1}{\tau^2} + \frac{n}{\sigma^2} \right] \mathbb{I}_{k \times k}\right) \beta + \sum_{j = 1}^k \left[ \tilde{y}_{\cdot, j}^\top\left(\frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right) \tilde{y}_{\cdot, j} - 2 \sum_{i = 1}^n \frac{\beta_j}{\sigma^2} \left(\tilde{y}_{i, j} \right) \right] \right] \right) \\ 
&=\left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2}  \left[ \beta^\top \left( \left[ \frac{1}{\tau^2} + \frac{n}{\sigma^2} \right] \mathbb{I}_{k \times k}\right) \beta - 2 \sum_{j = 1}^k \left( \frac{\beta_j}{\sigma^2} \sum_{i = 1}^n \tilde{y}_{i, j} \right) + \sum_{j = 1}^k \left[ \tilde{y}_{\cdot, j}^\top\left(\frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right) \tilde{y}_{\cdot, j} \right]\right]  \right) \\ 
&= \left( \frac{1}{\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2}  \left[ \beta^\top \left( \left[ \frac{1}{\tau^2} + \frac{n}{\sigma^2} \right] \mathbb{I}_{k \times k}\right) \beta - \frac{2}{\sigma^2} \beta^\top\left( \tilde{y} \mathbf{1}_k \right) + \sum_{j = 1}^k \left[ \tilde{y}_{\cdot, j}^\top\left(\frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right) \tilde{y}_{\cdot, j} \right]\right]  \right) \\ 
&\overset{(iv)}{=} \left( \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \exp \left( - \frac{1}{2} \left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)^\top \left( \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k} \right)\left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)\right)
\end{aligned}
\nonumber
$$

In $(i)$, we substitute in $$\mu_{\cdot, j} = \alpha \cdot \mathbf{1}_{n} + \beta_j \cdot \mathbf{1}_{n}$$.
<br>
In $(ii)$, we use the following:

$$
\begin{aligned}
\sum_{j = 1}^k (\beta_j \cdot \mathbf{1}_{n})^\top \Sigma_{e_j}^{-1} (\beta_j \cdot \mathbf{1}_{n}) 
&= \sum_{j = 1}^k \sum_{i = 1}^n \sum_{i' = 1}^n \beta^j \left(\Sigma_{e_j}^{-1}\right)_{i,i'}  \\
&= \sum_{j = 1}^k  \beta_j \sum_{i = 1}^{n} \sum_{i' = 1}^{n} \left(\frac{1}{\sigma^2} \mathbb{I}_{n, n}\right)_{i,i'} \\
&= \left( \sum_{i = 1}^{n} \sum_{i' = 1}^{n} \left(\Sigma_{e_j}^{-1}\right)_{i,i'} \right) \sum_{j = 1^k} \beta_j^2 \\
&= \mathbf{1}_{n}^\top \left( \frac{1}{\sigma^2} \mathbb{I}_{n \times n} \right) \mathbf{1}_{n} \left( \beta^\top \beta\right) \\
&= \left( \frac{n}{\sigma^2} \right) \beta^\top \beta \\
&= \beta^\top \left( \frac{n}{\sigma^2} \mathbb{I}_{k \times k} \right)  \beta^\top
\end{aligned}
\nonumber
$$

<br>
In $(iii)$, we define $\tilde{y}_{\cdot, j} := y_{\cdot, j} - \alpha \cdot \mathbf{1}_{n}$ and substitute in for $\Sigma_\beta^{-1}$ and $\Sigma_{e_j}^{-1}$.
<br>
In $(iv)$, we complete the square. Let $A := \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k}$ and $b := -\frac{2}{\sigma^2}(\tilde{y}\mathbf{1}_k)$. Define the following:

$$
\begin{aligned}
h &:= -\frac{1}{2}A^{-1}b = -\frac{1}{2} \left(\frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2 } \mathbb{I}_{k \times k} \right) \left( -\frac{2}{\sigma^2} (\tilde{y}\mathbf{1}_k) \right) = \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \\
g &:= -\frac{1}{4} b^\top A^{-1} b = -\frac{1}{4}\left( -\frac{2}{\sigma^2} (\tilde{y} \mathbf{1}_n) \right)^\top  \left( \frac{\tau^2\sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k} \right) \left( -\frac{2}{\sigma^2} (\tilde{y} \mathbf{1}_n) \right) = - \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left( \frac{\tau^2}{\sigma^2(\sigma^2 + n \tau^2)} \mathbb{I}_{k \times k}\right) \left( \tilde{y} \mathbf{1}_k \right)
\end{aligned}
\nonumber
$$

Then:

$$
\begin{aligned}
\beta^\top \left( \left[ \frac{1}{\tau^2} + \frac{n}{\sigma^2} \right] \mathbb{I}_{k \times k}\right) \beta - \frac{2}{\sigma^2} \beta^\top\left( \tilde{y} \mathbf{1}_n \right) &= (\beta - h)^\top A (\beta - h) + k \\
&= \left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)^\top \left( \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k} \right)\left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right) - \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left( \frac{\tau^2}{\sigma^2(\sigma^2 + n \tau^2)} \mathbb{I}_{k \times k}\right) \left( \tilde{y} \mathbf{1}_k \right)
\end{aligned}
\nonumber
$$
</details>

The likelihood looks very close to a Gaussian density (just without the normalization constant). If we multiply and divide it by $\sqrt{\rvert 2\pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k} \rvert}$ and rearrange some terms, we can then integrate out $\beta$ to get the likelihood of just $(\alpha, \tau^2)$:

$$
\mathcal{L}(\alpha, \tau^2; y) = \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\left(\sigma^2\right)^{k/2} \left( 2 \sigma^2 \pi \right) ^{n/2}} {\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \sigma^2 + n \tau^2\right)^{k/2}}
\nonumber
$$

<details>
<summary>Integration Details.</summary>
Multiplying and dividing by the specified value yields:

$$
\left( \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right)  \underbrace{\frac{1}{\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}  \exp \left( - \frac{1}{2} \left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)^\top \left( \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k} \right)\left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)\right)}_{(*)}
\nonumber
$$

where $(*)$ is a Gaussian density and will integrate to $(1)$:

$$
\begin{aligned}
\mathcal{L}(\alpha, \tau^2; y) &= \int \mathcal{L}(\alpha, \tau^2, \beta; y) d\beta \\
&= \int \left( \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \frac{1}{\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}  \exp \left( - \frac{1}{2} \left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)^\top \left( \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k} \right)\left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)\right)d\beta \\
&= \left( \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}}  \right) \int  \frac{1}{\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}  \exp \left( - \frac{1}{2} \left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)^\top \left( \frac{\sigma^2 + n \tau^2 }{\tau^2 \sigma^2}  \mathbb{I}_{k \times k} \right)\left( \beta - \frac{\tau^2}{\sigma^2 + n \tau^2} (\tilde{y}\mathbf{1}_k) \right)\right) d\beta \\
&= \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \Sigma_\beta \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \Sigma_{e_j} \rvert}} \\
&= \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\sqrt{\big \rvert 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n \tau^2} \mathbb{I}_{k \times k}\big \rvert}}{\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\sqrt{\rvert 2 \pi \left( \tau^2 \mathbb{I}_{k \times k} \right) \rvert}} \prod_{j = 1}^k \frac{1}{\sqrt{\rvert 2 \pi \left( \sigma^2 \mathbb{I}_{n \times n} \right) \rvert}} \\
&= \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\left( 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n\tau^2} \big\rvert \mathbb{I}_{k \times k} \big \rvert \right)^{k/2}} {\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\left( 2\tau^2 \pi \rvert \mathbb{I}_{k \times k} \rvert \right)^{k/2}} \left(  2 \sigma^2 \pi \rvert \mathbb{I}_{n \times n} \rvert\right) ^{n/2}\\
&= \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\left( 2 \pi \frac{\tau^2 \sigma^2}{\sigma^2 + n\tau^2} \right)^{k/2}} {\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}  \right)\left( 2\tau^2 \pi  \right)^{k/2}} \left(  2 \sigma^2 \pi \right) ^{n/2}\\
&= \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\left(\tau^2 \sigma^2\right)^{k/2} \left( 2 \sigma^2 \pi \right) ^{n/2}} {\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left(\tau^2 \right)^{k/2} \left( \sigma^2 + n \tau^2\right)^{k/2}} \\
&= \frac{\exp\left( \frac{\tau^2}{2\sigma^2(\sigma^2 + n \tau^2)} (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)\left(\sigma^2\right)^{k/2} \left( 2 \sigma^2 \pi \right) ^{n/2}} {\exp\left(\frac{1}{2\sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \sigma^2 + n \tau^2\right)^{k/2}} \\
\end{aligned}
\nonumber
$$
</details>

Taking the natural logarithm gives the log-likelihood:

$$
\ell(\alpha, \tau^2; y) = \frac{k}{2} \log(\sigma^2) + \frac{n}{2} \log(2 \sigma^2 \pi) - \frac{k}{2} \log(\sigma^2 + n \tau^2) + \frac{\tau^2}{2 \sigma^2(\sigma^2 + n \tau^2)}\left(\mathbf{1}_k^\top \tilde{y}^\top \right) \left(\tilde{y} \mathbf{1}_k \right) - \frac{1}{2 \sigma^2} \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} 
\nonumber
$$

#### Score Function

We take the gradient of the log-likelihood with respect to $\tau^2$ to get the score function for $\tau^2$:

$$
\begin{aligned}
U_{\tau^2} &=  - \frac{nk}{2(\sigma^2 + n \tau^2)} + \frac{1}{2 \sigma^2(\sigma^2 + n \tau^2)^2}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k) \\
\overset{\tau^2 = 0}{\implies} U_{\tau^2}(y; 0) &= - \frac{nk}{2\sigma^2} + \frac{(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k)}{2 \sigma^4} 
\end{aligned}
\nonumber
$$

<details>
<summary>Score Derivation.</summary>
$$
\begin{aligned}
U_{\tau^2} &= \frac{\ell(\alpha, \tau^2; y)}{\partial \tau^2} \\
&= \frac{\partial}{\partial \tau^2} \left[ -\frac{k}{2}\log(\sigma^2 + n \tau^2) + \frac{\tau^2}{2 \sigma^2(\sigma^2 + n \tau^2)}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k)\right] \\
&= -\frac{n k}{2(\sigma^2 + n \tau^2)} + \frac{2\sigma^2(\sigma^2 + n\tau^2)- 2\sigma^2 \tau^2 n}{4 \sigma^4(\sigma^2 + n \tau^2)^2}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k) \\
&= - \frac{nk}{2(\sigma^2 + n \tau^2)} + \frac{2\sigma^4 + 2 \sigma^2 \tau^2 n - 2 \sigma^2 \tau^2 n}{4 \sigma^4 (\sigma^2 + n \tau^2)^2}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k) \\
&= - \frac{nk}{2(\sigma^2 + n \tau^2)} + \frac{1}{2 \sigma^2(\sigma^2 + n \tau^2)^2}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k)
\end{aligned}
\nonumber
$$
</details>

If we evaluate at $\tau^2 = 0$, then the score has expectation $0$.

<details>
<summary>Proof.</summary>

$$
\begin{aligned}
\mathbb{E}\left[ U_{\tau^2} \right]  \bigg\rvert_{\tau^2 = 0}  &= \mathbb{E}\left[ - \frac{nk}{2(\sigma^2 + n \tau^2)} + \frac{1}{2 \sigma^2(\sigma^2 + n \tau^2)^2}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k) \right]\bigg\rvert_{\tau^2 = 0} \\
&= -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \mathbb{E}\left[ \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j}\right)^\top \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j}\right) \right] \\
&=  -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \mathbb{E}\left[ \sum_{i = 1}^n \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j}\right)_i^2 \right] \\
&=  -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \mathbb{E}\left[ \sum_{i = 1}^n \left( \sum_{j = 1}^k \tilde{y}_{i, j}\right)^2 \right] \\
&=  -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \mathbb{E}\left[ \sum_{i = 1}^n \left( \sum_{j = 1}^k (y_{i, j} - \alpha) \right)^2 \right] \\
&=  -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \sum_{i = 1}^n \mathbb{E}\left[ \sum_{j = 1}^k \sum_{j' = 1}^k (y_{i, j} - \alpha)(y_{i, j'} - \alpha) \right] \\
&=  -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \sum_{i = 1}^n\sum_{j = 1}^k \mathbb{E}\left[  (y_{i, j} - \alpha)^2 + 2\sum_{j' < j} (y_{i, j} - \alpha)(y_{i, j'} - \alpha)\right] \\
&\overset{(i)}{=}  -\frac{nk}{2\sigma^2} + \frac{1}{2 \sigma^4} \sum_{i = 1}^n\sum_{j = 1}^k \sigma^2 \\
&= -\frac{nk}{2\sigma^2} + \frac{nk \sigma^2}{2 \sigma^4}  \\
&= -\frac{nk}{2\sigma^2} + \frac{nk}{2\sigma^2} \\
&= 0
\end{aligned}
\nonumber
$$

In $(i)$, we use the fact that observations in different groups are independent, so their covariance is $0$.
</details>


#### Information

If we assume that we know $\alpha$ and $\sigma^2$, then the Fisher information is the expectation of the squared gradient of the log-likelihood, and we set $\tau^2 = 0$:

$$
\mathcal{I}_{\tau^2}(0) =  \frac{4 \sigma^4 n^2 k^2 + 2 n^2 k + 2 k^2 n - 3nk}{4 \sigma^8}
\nonumber
$$

In most cases, we must estimate $\alpha$ and $\sigma^2$ and account for the uncertainty in the estimation, which we have not done here.

<details>
<summary>Fisher Information Derivations.</summary>

$$
\begin{aligned}
\mathcal{I}_{\tau^2} & = \mathbb{E}\left[ \left( \frac{\partial \ell(\alpha, \tau^2; y)}{\partial \tau} \right)^2 \right] \\
&= \mathbb{E}\left[ \left( \frac{nk}{2(\sigma^2 + n \tau^2)} + \frac{\left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left(\tilde{y} \mathbf{1}_k\right)}{2\sigma^2(\sigma^2 + n \tau^2)^2} \right)^2 \right]  \\
&= \frac{n^2k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{nk}{(\sigma^2 + n \tau^2)^3} \mathbb{E}\left[ \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left(\tilde{y} \mathbf{1}_k\right) \right] + \frac{1}{4 \sigma^4(\sigma^2 + n \tau^2)^4}\mathbb{E}\left[ \left( \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left(\tilde{y} \mathbf{1}_k\right) \right)^2 \right] \\
&\overset{(i)}{=} \frac{n^2k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n^2k^2 \sigma^2}{(\sigma^2 + n \tau^2)^3} + \frac{1}{4 \sigma^4(\sigma^2 + n \tau^2)^4}\mathbb{E}\left[ \left( \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left(\tilde{y} \mathbf{1}_k\right) \right)^2 \right] \\
&\overset{(ii)}{=} \frac{n^2k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n^2k^2 \sigma^2}{(\sigma^2 + n \tau^2)^3} + \frac{\sigma^4nk(2n + 2k - 3)}{4 \sigma^4(\sigma^2 + n \tau^2)^4} \\
&= \frac{n^2k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n^2k^2 \sigma^2}{(\sigma^2 + n \tau^2)^3} + \frac{nk(2n + 2k - 3)}{4(\sigma^2 + n \tau^2)^4} \\
&= \frac{n^2k^2(\sigma^2 + n \tau^2)^2}{4(\sigma^2 + n \tau^2)^4} + \frac{4n^2k^2 \sigma^2(\sigma^2 + n \tau^2)}{4(\sigma^2 + n \tau^2)^4} + \frac{nk(2n + 2k - 3)}{4(\sigma^2 + n \tau^2)^4} \\
&\overset{\tau^2 = 0}{=} \frac{\sigma^4n^2k^2 + 4 \sigma^4n^2 k^2 + 2n^2k + 2k^2 n - 3nk}{4\sigma^8} \\
&= \frac{4 \sigma^4 n^2 k^2 + 2 n^2 k + 2 k^2 n - 3nk}{4 \sigma^8}
\end{aligned}
\nonumber
$$

Under the null hypothesis that $\tau^2 = 0$, $(i)$ follows from:

$$
\begin{aligned}
\mathbb{E}\left[ \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left(\tilde{y} \mathbf{1}_k\right) \right]  &= \mathbb{E}\left[ \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \right) \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j} \right)\right] \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \sum_{j' = 1}^k \mathbb{E}\left[ \tilde{y}_{i, j}^2 \right] \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left( \mathbb{E}\left[ \tilde{y}_{i,j}^2 \right] + 2\sum_{j' < j} \underbrace{\mathbb{E}\left[ \tilde{y}_{i,j} \tilde{y}_{i, j'} \right]}_{=0}  \right) \\
&= nk\sigma^2
\end{aligned}
\nonumber
$$


And $(ii)$ follows from:

$$
\begin{aligned}
\mathbb{E}\left[ \left( \left( \mathbf{1}_k^\top \tilde{y}^\top \right)\left(\tilde{y} \mathbf{1}_k\right) \right)^2 \right] &= \mathbb{E}\left[ \left[ \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j}^\top \right) \left( \sum_{j = 1}^k \tilde{y}_{\cdot, j} \right) \right]^2 \right] \\
&= \mathbb{E}\left[ \left( \sum_{j = 1}^k \sum_{j' = 1}^k \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)^2 \right] \\
&=  \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' = 1}^k \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right] \\
&=  \underbrace{\sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right]}_{(a)} + \underbrace{2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' < l}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right]}_{(b)} \\
&= \sigma^4nk(2n - 1) + 2\sigma^4nk(k - 1) \\
&= 2 \sigma^4 n^2 k - \sigma^4nk + 2\sigma^4 k^2 n - 2 \sigma^4 nk \\
&= 2 \sigma^4 n^2 k + 2\sigma^4 k^2 n - 3\sigma^4nk \\
&= \sigma^4nk(2n + 2k - 3)
\end{aligned}
\nonumber
$$

<details>
<summary>Details Of $(a)$ And $(b)$.</summary>
$$
\begin{aligned}
(a) &= \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right] \\
&= \underbrace{\sum_{j = 1}^k \sum_{j' = 1}^k  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, j'}^\top \tilde{y}_{\cdot, j'}\right) \right]}_{(a1)}  + \underbrace{2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right]}_{(a2)} \\
&= \sigma^4nk(2n - 1) + 2\sigma^4nk^2k(k - 1)
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
(a1) &= \sum_{j = 1}^k \sum_{j' = 1}^k  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, j'}^\top \tilde{y}_{\cdot, j'}\right) \right] \\
&= \sum_{j = 1}^k \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j}\right) \right] + 2\sum_{j' < j} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, j'}^\top \tilde{y}_{\cdot, j'}\right) \right] \\
&= \sum_{j = 1}^k \mathbb{E}\left[ \left( \sum_{i = 1}^n \tilde{y}_{i, j}^2\right)^2  \right] \\
&= \sum_{j = 1}^k \mathbb{E}\left[ \sum_{i = 1}^n \sum_{i' = 1}^n \tilde{y}_{i, j}^2 y_{i', j}^2 \right] \\
&= \sum_{j = 1}^k \sum_{i = 1}^n \left( \mathbb{E}\left[ \tilde{y}_{i, j}^4 \right] + 2\sum_{i' < i} \mathbb{E}\left[\tilde{y}_{i, j}^2 \tilde{y}_{i', j}^2 \right]  \right) \\
&= nk \gamma^4 + 2 \sum_{j = 1}^k \sum_{i = 1}^n \sum_{i' < i}\mathbb{E}\left[\tilde{y}_{i, j}^2 \tilde{y}_{i', j}^2 \right]  \hspace{15mm} \text{ where } \gamma^4  \text{ is the 4th central moment of an arbitrary } y_{i,j} \\
&= nk \gamma^4 + 2 \sum_{j = 1}^k \sum_{i = 1}^n \sum_{i' < i}\mathbb{E}\left[\tilde{y}_{i, j}^2 \right] \mathbb{E} \left[ \tilde{y}_{i', j}^2 \right]  \hspace{15mm} \text{ variables equal to functions of independent variables are independent } \\
&= nk \gamma^4 + 2 \sigma^4 n (n - 1) k  \\
&= 3 \sigma^4 nk + 2\sigma^4 n^2 k  - 2 \sigma^4 n k  \hspace{15mm} \text{ known form of 4th central moment for Gaussian variables} \\
&= 2 \sigma^4 n^2 k - \sigma^4 nk \\
&= \sigma^4 nk (2n - 1) \\
(a2) &= 2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right] \\
&= 2 \sum_{j = 1}^k \sum_{l < j}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right] + 4 \sum_{j = 1}^k \sum_{j' < j} \sum_{l < j'}   \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right]  \\
&= 2 \sum_{j = 1}^k \sum_{l < j} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l}\right) \right] \\
&= 2 \sum_{j = 1}^k \sum_{l < j} \mathbb{E}\left[ \left( \sum_{i = 1}^n \tilde{y}_{i,j}^2 \right) \left( \sum_{i = 1}^n \tilde{y}_{i,l}^2 \right)\right] \\
&= 2 \sum_{j = 1}^k \sum_{l < j} \sum_{i = 1}^n \sum_{i' = 1}^n \mathbb{E}\left[ \tilde{y}_{i,j}^2 \tilde{y}_{i', l}\right] \\ 
&= 2 \sum_{j = 1}^k \sum_{l < j} \sum_{i = 1}^n \sum_{i' = 1}^n \mathbb{E}\left[ \tilde{y}_{i,j}^2 \right] \mathbb{E}\left[ \tilde{y}_{i', l}\right] \\ 
&= 2\sigma^4 n^2 k (k - 1)  
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
(b) &= 2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' < l}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right] \\
&= \underbrace{2 \sum_{j = 1}^k \sum_{j' = 1}^k  \sum_{l' < j'}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, j'}^\top \tilde{y}_{\cdot, l'}\right) \right]}_{(b1)} + \underbrace{4 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'} \sum_{l' < l} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right]}_{(b2)}
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
(b1) &= 2 \sum_{j = 1}^k \sum_{j' = 1}^k  \sum_{l' < j'}  \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, j'}^\top \tilde{y}_{\cdot, l'}\right) \right] = 2 \sum_{j = 1}^k \sum_{l' < j} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \tilde{y}_{\cdot, j'}^\top \tilde{y}_{\cdot, l'}\right) \right] + 2 \sum_{j = 1}^l \sum_{j' < j} \sum_{l' < j'} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, l'}\right) \right] = 0 \\
(b2) &= 4 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'} \sum_{l' < l} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right] = 4 \sum_{j = 1}^k \sum_{l < j} \sum_{l' < l} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right]  + 8 \sum_{j = 1}^k \sum_{j' < j} \sum_{l < j'} \sum_{l' < l} \mathbb{E}\left[ \left( \tilde{y}_{\cdot, j}^\top \tilde{y}_{\cdot, j'} \right)\left( \tilde{y}_{\cdot, l}^\top \tilde{y}_{\cdot, l'}\right) \right] = 0
\end{aligned}
\nonumber
$$
</details>
</details>


#### Test Statistic

The test statistic is the score squared and divided by the Fisher information, all under the null: 

$$
\begin{aligned}
S_{\tau^2}(y; 0) &= \frac{U_{\tau^2}^2(y; 0)}{\mathcal{I}_{\tau^2}(0)} = \frac{\left( -\frac{nk}{2\sigma^2} + \frac{\left(\mathbf{1}_k^\top \tilde{y}^\top \right)\left( \tilde{y} \mathbf{1}_k \right)}{2 \sigma^4} \right)^2}{\frac{4 \sigma^4 n^2 k^2 + 2 n^2 k + 2 k^2 n - 3nk}{4 \sigma^8}} = \frac{\left( \left(\mathbf{1}_k^\top \tilde{y}^\top \right)\left( \tilde{y} \mathbf{1}_k \right) - \sigma^2nk\right)^2}{4 \sigma^4 n^2 k^2 + 2 n^2 k + 2 k^2 n - 3nk}
\end{aligned}
\nonumber
$$



---


## References

Breusch, T. S., & Pagan, A. R. (1980). The Lagrange Multiplier Test and its Applications to Model Specification in Econometrics. The Review of Economic Studies, 47(1), 239. https://doi.org/10.2307/2297111

Buse, A. (1982). The Likelihood Ratio, Wald, and Lagrange Multiplier Tests: An Expository Note. The American Statistician, 36(3a), 153–157. https://doi.org/10.1080/00031305.1982.10482817

Cox, D.R., & Hinkley, D.V. (1974). Theoretical Statistics (1st ed.). Chapman and Hall/CRC. https://doi.org/10.1201/b14832.

Moran, P. A. P. (1971). Maximum-likelihood estimation in non-standard conditions. Mathematical Proceedings of the Cambridge Philosophical Society, 70(3), 441–450. https://doi.org/10.1017/S0305004100050088.

Schervish, M. J. (1995). Theory of Statistics. Springer New York. https://doi.org/10.1007/978-1-4612-4250-5

Verbeke, G., & Molenberghs, G. (2003). The Use of Score Tests for Inference on Variance Components. Biometrics, 59(2), 254–262. https://doi.org/10.1111/1541-0420.00032.

Wasserman, L. https://www.stat.cmu.edu/~larry/=stat705/Lecture15.pdf.



<!-- $$
\begin{aligned}
\mathcal{I}(0) &= \mathbb{E}\left[ U_{\tau^2}^2 \right] \bigg \rvert_{\tau^2 = 0} \\
&= \sigma^4 \left( 4n^2k^2 - 2n^2k - 2k^2n + 3nk\right)
\end{aligned}
\nonumber
$$

<details>
<summary>Details.</summary>
$$
\begin{aligned}
\mathcal{I}(0) &= \mathbb{E}\left[ U^2_{\tau^2} \right] \bigg \rvert_{\tau^2 = 0} \\
&= \mathbb{E}\left[ \left(\frac{nk}{2(\sigma^2 + n \tau^2)} + \frac{1}{2 \sigma^2(\sigma^2 + n \tau^2)^2}(\mathbf{1}_k^\top \tilde{y}^\top)(\tilde{y}\mathbf{1}_k)\right)^2\right] \\
&= \frac{n^2 k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n k}{\sigma^2 (\sigma^2 + n \tau^2)^3} \mathbb{E}\left[ (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k )\right] + \frac{1}{4 \sigma^4 (\sigma^2 + n \tau^2)^4} \mathbb{E}\left[ \left((\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)^2\right]\\
&\overset{(i)}{=} \frac{n^2 k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n^2 k^2}{(\sigma^2 + n \tau^2)^3} + \frac{1}{4 \sigma^4 (\sigma^2 + n \tau^2)^4} \mathbb{E}\left[ \left((\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right)^2\right]\\
&= \frac{n^2 k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n^2 k^2}{(\sigma^2 + n \tau^2)^3} + \frac{1}{4 \sigma^4 (\sigma^2 + n \tau^2)^4} \mathbb{E}\left[ \left( \sum_{i = 1}^n \sum_{j = 1}^k \sum_{j' = 1}^k (y_{i,j} - \alpha)(y_{i, j'} - \alpha) \right)^2\right]\\
&= \frac{n^2 k^2}{4(\sigma^2 + n \tau^2)^2} + \frac{n^2 k^2}{(\sigma^2 + n \tau^2)^3} + \frac{1}{4 \sigma^4 (\sigma^2 + n \tau^2)^4} \underbrace{\mathbb{E}\left[ \left( \sum_{i = 1}^n \sum_{j = 1}^k \sum_{j' = 1}^k (y_{i,j} - \alpha)(y_{i, j'} - \alpha) \right) \left( \sum_{i' = 1}^n \sum_{l = 1}^k \sum_{l' = 1}^k (y_{i', l} - \alpha)(y_{i', l'} - \alpha) \right) \right]}_{(*)} \\
\end{aligned}
\nonumber
$$

In $(i)$, we find $\mathbb{E}\left[ (\mathbf{1}_k^\top \tilde{y}^\top) (\tilde{y} \mathbf{1}_k)\right]$ as we did when finding the expectation of the score above.
<br>
In $(ii)$, we use the following. Let $$\tilde{y}_{i,j} := y_{i,j} - \alpha$$ and $$\gamma^t_{i,j} := \mathbb{E}\left[ \tilde{y}_{i,j}^t \right]$$, the $$t$$-th central moment of $$y_{i,j}$$. Note that $$\gamma^2_{i,j} = \sigma^2$$ and $$\gamma^4_{i,j} = 3\sigma^4$$ for all $$i, j$$. 
Then:

$$
\begin{aligned}
(*) &= \mathbb{E}\left[ \left( \sum_{i = 1}^n \sum_{j = 1}^k \sum_{j' = 1}^k \tilde{y}_{i, j} \tilde{y}_{i, j'} \right) \left( \sum_{i' = 1}^n \sum_{l = 1}^k \sum_{l' = 1}^k \tilde{y}_{i', l} \tilde{y}_{i', l'} \right) \right] \\
&= \sum_{i = 1}^n \sum_{i' = 1}^n \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' = 1}^k \mathbb{E}\left[  \tilde{y}_{i, j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'}  \right] \\
&\overset{(1)}{=}  \sum_{i = 1}^n \sum_{i' = 1}^n\left( \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{l < j} \mathbb{E}\left[ \tilde{y}_{i,j}^2 \tilde{y}_{i', l}^2\right] \right) \\
&= \sum_{i = 1}^n \left( \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^4 \right] + 2 \sum_{j = 1}^k \sum_{l < j} \mathbb{E}\left[ \tilde{y}_{i,j}^2 \tilde{y}_{i, l}^2\right]  + 2 \sum_{i' < i} \left( \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{l < j} \mathbb{E}\left[ \tilde{y}_{i,j}^2 \tilde{y}_{i', l}^2\right] \right)\right) \\
&= \sum_{i = 1}^n \left( \sum_{j = 1}^k \gamma^4_{i,j} + 2 \sum_{j = 1}^k \sum_{l < j} \gamma^2_{i,j} \gamma_{i,l}^2 +2 \sum_{i' < i} \left( \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{l < j} \gamma^2_{i,j} \gamma^2_{i', l} \right)\right) \\
&= \sum_{i = 1}^n \left( \sum_{j = 1}^k 3\sigma^4 + 2 \sum_{j = 1}^k \sum_{l < j} \sigma^4 +2 \sum_{i' < i} \left( \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{l < j} \sigma^4 \right)\right) \\
&= 3nk\sigma^4 + 2nk(k-1) \sigma^4 + 4n(n-1)k(k-1)\sigma^4 + 2 \sum_{i = 1}^n \sum_{i' < i} \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] \\
&\overset{(2)}{=} 3nk\sigma^4 + 2nk(k-1) \sigma^4 + 4n(n-1)k(k-1)\sigma^4 + 2 \sum_{i = 1}^n \sum_{i' < i} \sum_{j = 1}^k \gamma^2_{i, j} \gamma^2_{i', j} \\
&= 3nk\sigma^4 + 2nk(k-1) \sigma^4 + 4n(n-1)k(k-1)\sigma^4 + 2 n(n-1)k \sigma^4 \\
&= \sigma^4\left( 3nk + 2k^2n - 2nk + 4(n^2 -n)(k^2 - k) + 2n^2k - 2nk \right) \\
&= \sigma^4 \left( nk + 2k^2n + 4(n^2k^2 - n^2k - k^2 n + nk) + 2n^2k - 2nk \right) \\
&= \sigma^4 \left( 4n^2k^2 - 2n^2k - 2k^2n + 3nk\right)
\end{aligned}
\nonumber
$$


<details>
<summary>Details Of $(1)$.</summary>
$$
\begin{aligned}
\sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' = 1}^k \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'} \right]
&= \underbrace{\sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l}^2 \right]}_{(a)} + \underbrace{2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'} \right]}_{(b)}  \\
&= \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{l < j} \mathbb{E}\left[ \tilde{y}_{i,j}^2 \tilde{y}_{i', l}^2\right]
\end{aligned}
\nonumber
$$

<details>
<summary>Details Of $(a)$ and $(b)$.</summary>
$$
\begin{aligned}
(a) &= \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l}^2 \right] \\
&= \underbrace{\sum_{j = 1}^k \sum_{j' = 1}^k \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', j'}^2 \right]}_{(a_1)} + \underbrace{2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l}^2 \right]}_{(a_2)} \\
(b) &= 2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l = 1}^k \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'} \right] \\
&= \underbrace{2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l' < j'} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', j'} \tilde{y}_{i', l'} \right]}_{(b_1)}  + \underbrace{4 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'} \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'} \right]}_{(b_2)}
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
(a_1) &= \sum_{j = 1}^k \sum_{j' = 1}^k \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', j'}^2 \right] \\
&= \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{j' < j} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', j'}^2 \right] \\
&= \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] + 2 \sum_{j = 1}^k \sum_{j' < j} \underbrace{\mathbb{E}\left[\tilde{y}_{i,j}\right]}_{=0} \mathbb{E}\left[ \tilde{y}_{i, j'} \tilde{y}_{i', j'}^2 \right] \\
&= \sum_{j = 1}^k \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j}^2 \right] \\
(a_2) &= 2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l}^2 \right] \\
&= 2 \sum_{j = 1}^k \sum_{l < j}\mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', l}^2 \right]  + 4 \sum_{j = 1}^k \sum_{j' < j} \sum_{l < j'} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l}^2 \right] \\
&= 2 \sum_{j = 1}^k \sum_{l < j}\mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', l}^2 \right]  + 4 \sum_{j = 1}^k \sum_{j' < j} \sum_{l < j'} \underbrace{\mathbb{E}\left[\tilde{y}_{i,j} \right] \mathbb{E}\left[\tilde{y}_{i, j'}\right]}_{=0} \mathbb{E}\left[ \tilde{y}_{i', l}^2 \right] \\
&= 2 \sum_{j = 1}^k \sum_{l < j}\mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', l}^2 \right] 
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
(b_1) &= 2 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l' < j'} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', j'} \tilde{y}_{i', l'} \right] \\
&= 2 \sum_{j = 1}^k \sum_{l' < j}  \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j} \tilde{y}_{i', l'} \right] + 4\sum_{j = 1}^k \sum_{j' < j} \sum_{l' < j}  \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', j'} \tilde{y}_{i', l'} \right] \\
&= 2 \sum_{j = 1}^k \sum_{l' < j}  \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', j} \right] \underbrace{\mathbb{E}\left[\tilde{y}_{i', l'} \right]}_{=0} + 4\sum_{j = 1}^k \sum_{j' < j} \sum_{l' < j}  \mathbb{E}\left[\tilde{y}_{i,j} \right] \mathbb{E}\left[ \tilde{y}_{i, j'} \tilde{y}_{i', j'} \right] \underbrace{\mathbb{E}\left[ \tilde{y}_{i', l'} \right]}_{=0} \\
&= 0 \\
(b_2) &= 4 \sum_{j = 1}^k \sum_{j' = 1}^k \sum_{l < j'} \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'} \right] \\
&= 4 \sum_{j = 1}^k \sum_{l < j} \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j}^2 \tilde{y}_{i', l} \tilde{y}_{i', l'} \right] + 8 \sum_{j = 1}^k \sum_{j' < j} \sum_{l < j'} \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j} \tilde{y}_{i, j'} \tilde{y}_{i', l} \tilde{y}_{i', l'} \right]  \\
&= 4 \sum_{j = 1}^k \sum_{l < j} \sum_{l' < l} \mathbb{E}\left[\tilde{y}_{i,j}^2 \right] \mathbb{E}\left[ \tilde{y}_{i', l} \right] \underbrace{\mathbb{E}\left[\tilde{y}_{i', l'} \right]}_{=0} + 8 \sum_{j = 1}^k \sum_{j' < j} \sum_{l < j'} \sum_{l' < l} \underbrace{\mathbb{E}\left[\tilde{y}_{i,j}\right] \mathbb{E}\left[ \tilde{y}_{i, j'}\right] \mathbb{E}\left[ \tilde{y}_{i', l} \right] \mathbb{E}\left[\tilde{y}_{i', l'} \right]}_{=0}  \\
&= 0
\end{aligned}
\nonumber
$$
</details>
</details>

In $(2)$ we use the fact that, under the null, $\tau^2 = 0$, so observations in the same group are independent. (I hope this is true!)
</details> -->



<!-- For simplicity, assume our data have p.d.f. $f_Y(y; \theta)$ for parameter $\theta \in \Theta$ and that $\theta$ lies in the interior of $\Theta$. We want to test $H_0: \theta = \theta_0$ against $H_A: \theta \in \Theta_A \subset \Theta$. The maximum likelihood ratio critical region will be the part of the sample space achieving large $W' := 2\left( \underset{\theta' \in \Theta_A}{\sup} \left[ \ell(\theta'; y) \right]- \ell(\theta_0; y) \right)$. If $\theta_0 \notin \Theta_A$, then we use $W := 2 \left( \ell(\hat{\theta}; y) - \ell(\theta_0; y) \right) $.

We can use a Taylor expansion at $\hat{\theta}$ about $\theta_0$ of the log-likelihood to rewrite $W$. For some $\theta'$ between $\hat{\theta}$ and $\theta_0$:

$$
\begin{aligned}
\ell(\theta_0; y) &= \ell(\hat{\theta}; y) + (\theta_0 - \hat{\theta})\ell'(\hat{\theta}; y) + \frac{1}{2}(\theta_0 - \hat{\theta})^2 \ell''(\theta'; y) \\
\implies  \ell(\theta_0; y) &= \ell(\hat{\theta}; y) - (\hat{\theta} - \theta_0)U_\theta(y) \bigg\rvert_{\theta = \hat{\theta}} - \frac{1}{2}(\hat{\theta} - \theta_0)^2 U'_\theta(y) \bigg\rvert_{\theta = \theta'} \\
\implies W &= 2 \left( \ell(\hat{\theta}; y) - \ell(\hat{\theta}; y) + (\hat{\theta} - \theta_0) U_\theta(y) \bigg\rvert_{\theta = \hat{\theta}}  + \frac{1}{2}(\hat{\theta} - \theta_0)^2 U'_\theta(y) \bigg\rvert_{\theta = \theta'} \right) \\
\end{aligned}
$$

Under suitable conditions, the score evaluated at $$\hat{\theta}$$ will be zero. Furthermore, $$\theta'$$ must also be consistent (assuming the MLE is consistent). Let $$U(\theta^*) := U_\theta(y)\big\rvert_{\theta = \theta^*}$$ and let $$\mathcal{I}(\theta^*)$$ be the same for the Fisher information. We can do a similar Taylor expansion of $$U(\hat{\theta})$$ about $$\theta_0$$:

$$
\begin{aligned}
U(\hat{\theta}) &= U(\theta_0) + (\hat{\theta} - \theta_0)U'(\theta_0) + \frac{1}{2}(\hat{\theta} - \theta_0)^2 U''(\theta') \\
\implies 0 &= U(\theta_0) + (\hat{\theta} - \theta_0)\mathcal{I}(\theta_0) + \frac{1}{2}(\hat{\theta} - \theta_0)^2 U''(\theta') \\
\implies -U(\theta_0) &=  + (\hat{\theta} - \theta_0)\mathcal{I}(\theta_0) + \frac{1}{2}(\hat{\theta} - \theta_0)^2 U''(\theta') \\
\end{aligned}
\nonumber
$$ -->
