---
layout: distill
title: Rao Score Test
description:
date: 2026-04-07
tabs: true
tags: likelihood theory paper-review
toc:
    - name: Set-Up
      subsections:
         - name: Some More Assumptions
         - name: Distribution Result
    - name: Hypothesis Testing
      subsections:
         - name: Simple
         - name: Composite
    - name: Example
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: 2026-04-07-score-test.bib

---

I've been reading about the score test, but I've found the [Wikipedia page](https://en.wikipedia.org/wiki/Score_test) to be a bit lacking. 

What better way to clarify my understanding than to go back to the beginning and read C. R. Rao's original work from 1948?

## Set-Up
Let $\mathbf{x} = (x_1, \dots, x_p)^\top$, $\mathbf{y} = (y_1, \dots, y_q)^\top$, $\dots$ be observations with probability densities $f_1(x; \boldsymbol{\theta})$, $f_2(y; \boldsymbol{\theta})$, $\dots$ where $\boldsymbol{\theta} \in \mathbb{R}^k$ is some vector-valued parameter. Further assume that each density function $f_j(\cdot; \boldsymbol{\theta})$ depends on at least one component of $\boldsymbol{\theta}$. 

Assuming the sets of observations are independent, the likelihood is given by:

$$
\mathcal{L}(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots) = f_1(x; \boldsymbol{\theta};) \times f_2(y; \boldsymbol{\theta};) \times \dots$$

Taking the natural logarithm yields the log-likelihood, which is the logarithm of the joint density: 

$$
\ell(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots) = \log \left( \mathcal{L}(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots) \right) = \log \left( f(\mathbf{x}, \mathbf{y}, \dots; \boldsymbol{\theta}) \right)
$$

We define the <strong>score</strong> as the vector of first-order partial derivatives of the log-likelihood whose $i$-th component is given by:

$$
\phi_i(\boldsymbol{\theta}) = \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots)}{\partial \theta_i}; \hspace{5mm} i = 1, \dots, k
$$

Let $$\mathcal{I}(\boldsymbol{\theta}_0)$$ denote the variance-covariance matrix of the scores evalued at $\boldsymbol{\theta}_0$, and let $$\mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta}_0)$$ denote the $(i,j)$-th entry of its inverse. It can be shown that $$\left. \mathbb{E}[\phi_i(\boldsymbol{\theta})] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}^*} = 0$$ for all $i$ where $$\boldsymbol{\theta}^*$$ $is the <i>true</i> parameter value.

<aside><p>In the original text, $\mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta}_0)$ is denoted by $\alpha^{i,j}$.</p></aside>

<details>
  <summary>Proof.</summary>
  There's a little bit of abuse of notation...sorry.
  $$
  \begin{align}
  \mathbb{E}[\phi_i] &= \mathbb{E}\left[ \frac{\partial \log \mathcal{L}}{\partial \theta_i} \right] \\
  				     &= \int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm} 
                  \frac{1}{f(x, y, \dots; \theta)} \frac{\partial}{\partial \theta_i} \log \mathcal{L}\theta \hspace{2mm} 
                  dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &= \int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm}
                  \frac{1}{f(x, y, \dots; \theta)} \frac{\partial f(x, y, \dots; \theta)}{\partial \theta_i} 
                  f(x, y, \dots; \theta) \hspace{2mm} 
                  dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &= int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm}
                  \frac{\partial f(x, y, \dots; \theta)}{\partial \theta_i} \hspace{2mm} dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &\overset{(i)}{=} \frac{\partial}{\partial \theta_i} \int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm} f(x, y, \dots; \theta) \hspace{2mm} dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &\overset{(ii)}{=} \frac{\partial}{\partial \theta_i} 1\\
               &= 0
  \end{align}
  $$
  In $(i)$ we use the fact that, under [certain regularity condition](https://en.wikipedia.org/wiki/Leibniz_integral_rule), we can interchange the order of integration and differentiation. In $(ii)$ we use the fact that the integral of a probability density function evaluates to $1$.
</details>

Under <a href="/stats-ml/likelihood-theory/">certain conditions</a> we can also show that $$\boldsymbol{\phi}(\boldsymbol{\theta}^*)$$ goes to a multivariate Gaussian distribution with mean zero and covariance matrix, $$\mathcal{I}(\boldsymbol{\theta}^*)$$. 

<aside><p>Condition 1.</p></aside>

<details>
  <summary>Proof.</summary>
  As a function of independent random variables, the score is also a random variable. In fact, it is the sum of independent random variables. Let $f_j(\boldsymbol{\theta})$ denote the density of random variable $j$ (i.e. for $j = 1$, we have $f_j(\boldsymbol{\theta}) = f_1(\mathbf{x}; \boldsymbol{\theta})$). Let $\ell_j(\boldsymbol{\theta})$ denote the log-likelihood of $\boldsymbol{\theta}$ using only observations of random variable $j$. We then have:
  $$
  \begin{aligned}
  \phi_i &= \frac{\partial \log \mathcal{L}(\boldsymbol{\theta})}{\partial \theta_i} \\
         &= \frac{\partial \log \prod_j f_j(\boldsymbol{\theta})}{\partial \theta_i} \\
         &= \frac{\partial \sum_j \log f_j(\boldsymbol{\theta})}{\partial \theta_i} \\
         &= \sum_j \frac{\partial \log f_j(\boldsymbol{\theta})}{\partial \theta_i} \\
         &= \sum_j \frac{\partial \ell_j(\boldsymbol{\theta})}{\partial \theta_i}
  \end{aligned}
  $$
  As shown above, $\mathbb{E}[\phi_i] = 0 < \infty$. Assume that $A$ has finite diagonal elements.
  If we further assume that, for some $\delta > 0$, 
  $$
  \underset{n \rightarrow \infty}{\lim} \frac{1}{\alpha_{i,i}} \sum_j \mathbb{E}[\rvert \phi_i - \mathbb{E}[\phi_i] \rvert^{2+\delta}] = 0
  $$
  then the Lyapunov CLT states that, as $n \rightarrow \infty$, $\sum_{j} (\phi_i - \mathbb{E}[\phi_i])$ converges in distribution to a multivariate Normal distribution with mean zero and covariance $\mathbf{A}$.
</details>

We also assume that there exists $\eta > 0$ such that, for $j = 1, \dots, k$:

$$
\begin{equation}
\label{eq:eta-exist}
\begin{aligned}
\mathbb{E}\left[ \left( \frac{\partial \log f_1(\mathbf{x} \rvert \boldsymbol{\theta})}{\partial \theta_j} \right)^{2 + \eta} \right] &< \infty; \\
\mathbb{E}\left[ \left(\frac{\partial \log f_2(\mathbf{y} \rvert \boldsymbol{\theta})}{\partial \theta_j}  \right)^{2 + \eta}\right] &< \infty; \\
\dots &< \infty
\end{aligned}
\end{equation}
$$

<aside><p>Condition 2.</p></aside>

If the set of non-vanishing components of the sequence:

$$
\frac{\partial \log f_1(\mathbf{x} \rvert \boldsymbol{\theta})}{\partial \theta_j}, \frac{\partial \log f_2(\mathbf{y} \rvert \boldsymbol{\theta})}{\partial \theta_j}, \dots; \hspace{5mm} i = 1, 2, \dots
$$

<aside><p>Condition 3.</p></aside>

for <i>some</i> $j = 1, \dots, k$ is sufficiently large, then it can be shown that the statistic:

$$
\begin{equation}
\label{eq:ts-1}
\begin{aligned}
\chi^2(\boldsymbol{\theta}^*) 
&= \sum_{i = 1}^k \sum_{j = 1}^k \mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta}^*) \phi_i(\boldsymbol{\theta}^*) \phi_j(\boldsymbol{\theta}^*) \\
&= \boldsymbol{\phi}^\top(\boldsymbol{\theta}^*) \mathcal{I}^{-1}(\boldsymbol{\theta}^*) \boldsymbol{\phi}(\boldsymbol{\theta}^*)
\end{aligned}
\end{equation}
$$ 

is asymptotically $\chi^2$ with $k$ degrees of freedom (since $\boldsymbol{\theta}$ is $k$-dimensional). 

### Some More Assumptions
The above is not that helpful unless we know that the true value of $\boldsymbol{\theta}$. Suppose instead we impose $s < k$ restrictions on our parameters of the form:

$$
\begin{equation}
\label{eq:lin-const}
\psi_l(\theta_1, \dots, \theta_k) = 0; \hspace{5mm} l = 1, 2, \dots, s
\end{equation}
$$

Assume that Eq. \eqref{eq:lin-const} are such that $s$ of the parameters are determined as functions of the $k-s$ parameters. Without loss of generality, we will say that $$\theta_{k - s + 1}, \dots, \theta_k$$ are determined by $$\theta_1, \dots, \theta_{k-s}$$. This means that the likelihood is only really a function of $$\theta_1, \dots, \theta_{k-s}$$. 

We also assume that the joint distribution of $$\hat{\theta}_1, \dots, \hat{\theta}_{k - s}$$ tends to a multivariate Gaussian with variances/covariances of order $$O\left(\frac{1}{n}\right)$$ for large enough sample sizes, $n$.

Rao mentions that this assumption (joint distribution of MLEs) holds if Eq. \eqref{eq:eta-exist} is satisfied and the MLEs are <strong>uniformly consistent</strong>.

<aside><p>Condition 4.</p></aside>

<div class="definition">
<strong>Definition (Uniformly Consistent).<d-cite key=wald1943></d-cite></strong>
<br>
Let $\mathcal{X}$ be our sample space, and let $\hat{\boldsymbol{\theta}}_m(\mathbf{x})$ denote the maximum likelihood estimate of some parameter vector $\boldsymbol{\theta}$ made with $m$ samples $\mathbf{x} \in \mathcal{X}$. 
<br>
We call $\hat{\boldsymbol{\theta}}$ <i>uniformly consistent</i> if, for any $\epsilon > 0$, there exists $M \in \mathbb{N}$ such that for all $m \geq M$ and for all $\mathbf{x} \in \mathcal{X}$:
$$
\underset{n \rightarrow \infty}{\lim} \left\{ \mathbb{P}\left(\rvert \rvert \hat{\boldsymbol{\theta}}_m(\mathbf{x}) - \boldsymbol{\theta} \rvert \rvert < \epsilon \right) \right\} = 1
$$
</div>

This is shown in the proof of Proposition I of Wald (1943).<d-cite key=wald1943></d-cite>

### Distribution Result
Using the method of <a href="/stats-ml/constrained-optim/">Lagrange multipliers</a>, the maximum likelihood estimates of the parameter, denoted by $\hat{\boldsymbol{\theta}} = (\hat{\theta}_1, \dots, \hat{\theta}_k)^\top$ will be the solution to the system of equations:

$$
\begin{equation}
\label{eq:lagrange}
\begin{aligned}
\phi_i + \sum_{j = 1}^s \lambda_j \frac{\partial \psi_j}{\partial \theta_i} &= 0; \hspace{5mm} &i = 1, 2, \dots, k \\
\psi_l &= 0; \hspace{5mm} &l = 1, 2, \dots, s
\end{aligned}
\end{equation}
$$

It then follows that the test statistic in Eq. \eqref{eq:ts-1} evaluated at $\hat{\boldsymbol{\theta}}$ (instead of the true values):

$$
\begin{equation}
\label{eq:ts-2}
\chi^2(\hat{\boldsymbol{\theta}}) = \boldsymbol{\phi}^\top(\hat{\boldsymbol{\theta}}) \mathcal{I}^{-1}(\hat{\boldsymbol{\theta}}) \boldsymbol{\phi}(\hat{\boldsymbol{\theta}})
\end{equation}
$$

is asymptotically $\chi^2$ with $s$ degrees of freedom if the assumptions hold.

---

## Hypothesis Testing

### Simple
Suppose we wish to test a hypothesis concerning $\boldsymbol{\theta}$:

$$
H_0: \boldsymbol{\theta} = \boldsymbol{\theta}_0 \hspace{5mm} \text{vs.} \hspace{5mm} H_1: \boldsymbol{\theta} = \boldsymbol{\theta}_1
$$

An "ideal" test, according to Rao, should be based upon a statistic that is "maximally discriminatory" when $\boldsymbol{\theta}^1$ differs only by a little bit from $\boldsymbol{\theta}^0$. One way to restate this is that we want our test to have the best power (compared to other tests) when we deviate only slightly from $H_0$.

From Calculus I, we know that we can approximate the change in the log-likelihood as we move from $\boldsymbol{\theta}^0$ to $\boldsymbol{\theta}^1$ as:

$$
\begin{equation}
\label{eq:deriv-approx}
\begin{aligned}
\ell(\boldsymbol{\theta}^1; \mathbf{x}, \mathbf{y}, \dots) = \ell(\boldsymbol{\theta}^0 + \mathbf{h}; \mathbf{x}, \mathbf{y}, \dots) 
&\approx \ell(\boldsymbol{\theta}^0; \mathbf{x}, \mathbf{y}, \dots) + (\boldsymbol{\theta}^0 + \mathbf{h} - \boldsymbol{\theta}^0)^\top \left. \left[\frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots)}{\partial \boldsymbol{\theta}} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \\
&= \ell(\boldsymbol{\theta}^0; \mathbf{x}, \mathbf{y}, \dots) + \mathbf{h}^\top \boldsymbol{\phi}(\boldsymbol{\theta}_0)
\end{aligned}
\end{equation}
$$

As $h_1, \dots, h_k \rightarrow 0$, the approximation becomes better and better. We also know that, by the <a href="https://en.wikipedia.org/wiki/Neyman–Pearson_lemma"></a>, a most powerful test for small deviations will be based on the likelihood ratio. That is, for $c > 0$, a level $\alpha$ test where we reject $H_0$ if the following holds:

$$
\begin{aligned}
\frac{\mathcal{L}(\boldsymbol{\theta}^0; \mathbf{x}, \mathbf{y}, \dots)}{\mathcal{L}(\boldsymbol{\theta}^1; \mathbf{x}, \mathbf{y}, \dots)} \geq c \\
\iff \ell(\boldsymbol{\theta}^0; \mathbf{x}, \mathbf{y}, \dots) - \ell(\boldsymbol{\theta}^1; \mathbf{x}, \mathbf{y}, \dots) &\geq \log(c)
\end{aligned}
$$

will have power no smaller than any other level $\alpha$ test. By Eq. \eqref{eq:deriv-approx}, we have that:

$$
\ell(\boldsymbol{\theta}^1; \mathbf{x}, \mathbf{y}, \dots) - \ell(\boldsymbol{\theta}^0; \mathbf{x}, \mathbf{y}, \dots) \approx \mathbf{h}^\top\boldsymbol{\phi}(\boldsymbol{\theta}_0)
$$

This motivates a test based upon the statistic:

$$
U(\boldsymbol{\theta}_0) = \mathbf{h}^\top \boldsymbol{\phi}(\boldsymbol{\theta}_0) = \sum_{i = 1}^k h_i\phi_i(\boldsymbol{\theta}_0)
$$

In the previous section, we showed that the asymptotic distribution of the statistic:

$$
\chi^2(\boldsymbol{\theta}^*) = \boldsymbol{\phi}^\top(\boldsymbol{\theta}^*) \mathcal{I}^{-1}(\boldsymbol{\theta}^*) \boldsymbol{\phi}(\boldsymbol{\theta}^*)
$$

is $\chi^2$ with $k - s$ degrees of freedom. Thus, under $H_0$ the test statistic:

$$
\begin{equation}
\label{eq:test-1}
\begin{aligned}
\chi^2(\boldsymbol{\theta}^0) &= \frac{(U(\boldsymbol{\theta}^0))^2}{\mathcal{I}_\mathbf{h}(\boldsymbol{\theta}^0)} \\
\mathcal{I}_\mathbf{h}(\boldsymbol{\theta}) &= \text{Cov}\left(U(\boldsymbol{\theta})\right) = \text{Cov}\left(\mathbf{h}^\top \boldsymbol{\phi}(\boldsymbol{\theta}) \right) = \mathbf{h}^\top \text{Cov}(\boldsymbol{\phi}(\boldsymbol{\theta})) \mathbf{h} = \mathbf{h}^\top \mathcal{I}(\boldsymbol{\theta}) \mathbf{h}
\end{aligned}
\end{equation}
$$

should be $\chi^2$ with $1$ degree of freedom (since we do not impose restrictions on the components of $\boldsymbol{\theta}$).

<aside><p>Why $1$ degree of freedom...?</p></aside>

Eq. \eqref{eq:test-1} is nice, but it requires us to specify the direction of deviation from the null (i.e. we have to know $\mathbf{h}$). Instead, Rao proposes maximizing $\chi^2$ over all possible $\mathbf{h}$. This maximum is given by:

$$
\chi^2(\boldsymbol{\theta}^0) = \boldsymbol{\phi}^\top(\boldsymbol{\theta}^0) \mathcal{I}^{-1}(\boldsymbol{\theta}^0) \boldsymbol{\phi}(\boldsymbol{\theta}^0)
$$

<aside><p>How?</p></aside>

which we showed is asymptotically $\chi^2$ with $k$ degrees of freedom (see Eq. \eqref{eq:ts-1}). 

### Composite
Suppose we now wish to test a hypothesis concerning $\boldsymbol{\theta}$:

$$
H_0: \boldsymbol{\theta} \in \boldsymbol{\Theta}^0 \hspace{5mm} \text{vs.} \hspace{5mm} H_1: \boldsymbol{\theta} \notin \boldsymbol{\Theta}^0
$$

for some constrained parameter space, $\boldsymbol{\Theta}^0$. We can rewrite the constrained parameter space with some set of restrictions on the values that $\theta_1, \dots, \theta_k$ can take on. These restrictions can then be restated in the form:

$$
\psi_l(\theta_1, \dots, \theta_k) = 0; \hspace{5mm} l = 1, \dots, s
$$

Notice that there is a single set of values, $$\{\theta^*_1, \dots, \theta^*_k\}$$, that are responsible for generating our data, $\mathbf{x}$, $\mathbf{y}$, $\dots$. Rather than test $H_0$ directly, we can instead convert it into a test of a simple hypothesis: whether the estimate of $\boldsymbol{\theta}$ under the restrictions are plausible given the data. This leads to the statistic:

$$
\chi^2(\hat{\boldsymbol{\theta}}) = \boldsymbol{\phi}^\top(\hat{\boldsymbol{\theta}}) \mathcal{I}^{-1}(\hat{\boldsymbol{\theta}}) \boldsymbol{\phi}(\hat{\boldsymbol{\theta}})
$$

which will be $\chi^2$ with $s$ degrees of freedom (see Eq. \eqref{eq:ts-2}) under suitable conditions.

---

## Example
Rao provides the following two-dimensional example that hopefully will provide a bit more intuition on the results above.

We have $\boldsymbol{\theta} = (\theta_1, \theta_2)^\top$, and let $\boldsymbol{\theta}^*$ denote the true value of $\boldsymbol{\theta}$. We assume we have the restriction that $\theta_2 = \gamma(\theta_1)$. We denote its derivative with:

$$
\frac{d \theta_2}{d \theta_1} = \lambda(\theta_1)
$$

The MLE, $\hat{\boldsymbol{\theta}}$, will then satisfy:

$$
\begin{aligned}
\phi_1(\hat{\boldsymbol{\theta}}) + \lambda(\hat{\theta}_1) \phi_2(\hat{\boldsymbol{\theta}}) &= 0 \\
\hat{\theta}_2 - \gamma(\hat{\theta}_1) &= 0  
\end{aligned}
$$

Define:

$$
\chi_1^2 := \boldsymbol{\phi}^\top(\hat{\boldsymbol{\theta}}) \mathcal{I}^{-1}(\hat{\boldsymbol{\theta}}) \boldsymbol{\phi}(\hat{\boldsymbol{\theta}})
$$

Recall that:

$$
\begin{aligned}
\phi_i(\boldsymbol{\theta}) &= \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots)}{\partial \theta_i} \\
\frac{\partial \phi_i(\boldsymbol{\theta})}{\partial \theta_j} &= \frac{\partial^2 \ell(\boldsymbol{\theta}; \mathbf{x}, \mathbf{y}, \dots)}{\partial \theta_i \partial \theta_j} = \mathcal{I}_{i,j}
\end{aligned}
$$

A second-order Taylor approximation of Eq. \eqref{eq:ts-1} about $\hat{\theta}_1$ yields:

$$
\begin{aligned}
\chi_0^2 
&= \boldsymbol{\phi}^\top(\boldsymbol{\theta}_0) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \boldsymbol{\phi}(\boldsymbol{\theta}_0) +  (\boldsymbol{\theta} - \boldsymbol{\theta}_0)^\top \left. \frac{\partial}{\partial \boldsymbol{\theta}} \left[ \boldsymbol{\phi}^\top(\boldsymbol{\theta}) \mathcal{I}^{-1}(\boldsymbol{\theta}) \boldsymbol{\phi}(\boldsymbol{\theta})  \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0}  \\
\end{aligned}
$$

Let $\boldsymbol{\theta}_0 = (\hat{\theta}_1, \theta_2)^\top$.  

$$
\begin{aligned}
\chi_0^2 
&\approx \boldsymbol{\phi}^\top(\boldsymbol{\theta}_0) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \boldsymbol{\phi}(\boldsymbol{\theta}_0) +  \underbrace{(\theta_1 - \hat{\theta}_1) \left. \frac{\partial}{\partial \theta_1} \left[ \boldsymbol{\phi}^\top(\boldsymbol{\theta}) \mathcal{I}^{-1}(\boldsymbol{\theta}) \boldsymbol{\phi}(\boldsymbol{\theta})  \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0}}_{(i)}  + \frac{1}{2} \underbrace{(\theta_1 - \hat{\theta}_1)^2 \left. \frac{\partial^2}{\partial \theta_1^2} \left[ \boldsymbol{\phi}^\top(\boldsymbol{\theta}) \mathcal{I}^{-1}(\boldsymbol{\theta}) \boldsymbol{\phi}(\boldsymbol{\theta})  \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0}}_{=(ii)}
\end{aligned}
$$

If we assume $\frac{\partial \phi_i(\boldsymbol{\theta})}{\partial \theta_j} = \mathcal{I}_{i,j}(\boldsymbol{\theta})$ (which will hold true )

$$
\begin{aligned}
(i)
&= (\theta_1 - \hat{\theta}_1) \left. \frac{\partial}{\partial \theta_1} \left[ \boldsymbol{\phi}^\top(\boldsymbol{\theta}) \mathcal{I}^{-1}(\boldsymbol{\theta}) \boldsymbol{\phi}(\boldsymbol{\theta})  \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \\
&= (\theta_1 - \hat{\theta}_1) \sum_{i = 1}^2 \sum_{j = 1}^2 \left( \mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta}_0) \left. \left[ \frac{\partial}{\partial \theta_1}\left[ \phi_i(\boldsymbol{\theta}) \phi_j(\boldsymbol{\theta}) \right] \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} + \phi_i(\boldsymbol{\theta}) \phi_j(\boldsymbol{\theta}) \left. \left[ \frac{\partial \mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta})}{\partial \theta_1} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \right) \\
&= (\theta_1 - \hat{\theta}_1) \sum_{i = 1}^2 \sum_{j = 1}^2 \left(\mathcal{I}_{i,j}^{-1} \left[ \mathcal{I}_{i,1}(\boldsymbol{\theta}_0) \phi_j(\boldsymbol{\theta}_0) + \mathcal{I}_{j,1}(\boldsymbol{\theta}_0) \phi_i(\boldsymbol{\theta}_0)\right] + \phi_i(\boldsymbol{\theta}_0) \phi_j(\boldsymbol{\theta}_0) \left. \left[ \frac{\partial \mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta})}{\partial \theta_1} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \right) &\left( \frac{\partial \phi_i(\boldsymbol{\theta})}{\partial \theta_j} \overset{?}{=} \mathcal{I}_{i,j}(\boldsymbol{\theta}) \right) \\
&=  \sum_{i = 1}^2 \sum_{j = 1}^2 \left(\mathcal{I}_{i,j}^{-1} \left[ \mathcal{I}_{i,1}(\boldsymbol{\theta}_0) \phi_j(\boldsymbol{\theta}_0) + \mathcal{I}_{j,1}(\boldsymbol{\theta}_0) \phi_i(\boldsymbol{\theta}_0)\right] + \phi_i(\boldsymbol{\theta}_0) \phi_j(\boldsymbol{\theta}_0) \left. \left[ \frac{\partial \mathcal{I}^{-1}_{i,j}(\boldsymbol{\theta})}{\partial \theta_1} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \right) \\
\end{aligned}
$$








