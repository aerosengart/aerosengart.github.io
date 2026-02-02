---
layout: distill
title: Generalized Linear Models
description: A Primer
date: 2025-06-03
tabs: true
tags: regression likelihood models primer
toc:
  - name: Set-Up
  - name: Linear Models
    subsections:
        - name: Least Squares
        - name: Likelihood
        - name: Uncertainty Quantification
  - name: Generalized Linear Models
    subsections:
        - name: Model
        - name: Estimation
        - name: Uncertainty Quantification

bibliography: 2025-06-03-glm.bib
---

This post is a primer on generalized linear models and their associated estimation procedures.

## Set-Up

Let's assume we have covariate matrix, $\mathbf{X}$, response vector $\mathbf{y}$, and parameter (coefficient) vector, $\beta$, and error vector, $\epsilon$, given by:

$$
\mathbf{X} = 
\begin{bmatrix}
1 & x_{1,1} & \dots & x_{1, m} \\
1 & x_{2, 1} & \dots & x_{2, m} \\
\vdots & \vdots & \ddots & \vdots \\
1 & x_{n, 1} & \dots & x_{n, m}
\end{bmatrix},
\hspace{8mm}
\mathbf{y} = 
\begin{bmatrix}
y_1 \\
\vdots \\
y_n
\end{bmatrix},
\hspace{8mm} 
\beta = 
\begin{bmatrix}
\beta_0 \\
\beta_1 \\
\vdots \\
\beta_m
\end{bmatrix},
\hspace{8mm} 
\epsilon = 
\begin{bmatrix}
\epsilon_1 \\
\vdots \\
\epsilon_n
\end{bmatrix}
$$

Assume that the following holds:

$$
\begin{equation}
\label{eq:assumption-1}
g(\mathbf{y}) = \mathbf{X} \beta + \epsilon; \hspace{5mm} \epsilon_i \overset{iid}{\sim} F
\end{equation}
$$

where $g(\cdot)$ is some function (called a <i>link</i> function) and $F$ is some noise distribution. Equivalently, we can assume:

$$
\begin{equation}
\label{eq:assumption-2}
\mathbb{E}\left[ g(\mathbf{y}) \rvert \mathbf{X} \right] = \mathbf{X}\beta + \mathbb{E}[\epsilon]
\end{equation}
$$

For the rest of the post, we'll use $\mathbf{x}_i$ to denote:

$$
(1, x_{i, 1}, \dots, x_{i, m})^\top
$$

---

## Linear Models
Let's first cover the basics of linear regression, which will be a good base upon which to build the theory of generalized linear models. For basic linear regression, we assume that $g(\cdot)$ is thhe identity function (i.e. $g(\mathbf{y}) = \mathbf{y}$) and that, conditional on $\mathbf{X}$, $F = \mathcal{N}(0, \sigma^2)$. 


### Least Squares
We usually estimate $\beta$ via <i>ordinary least squares (OLS)</i>, which provides a closed form solution. The least squares objective is to minimize the <i>residual sum of squares</i>:

$$
RSS(\beta) = \sum_{i = 1}^n (y_i - \mathbf{x}_i^\top \beta )^2 = (\mathbf{y} - \mathbf{X} \beta)^\top (\mathbf{y} - \mathbf{X} \beta)
$$

The above function is quadratic in $\beta$, so we can differentiate with respect to $\beta$, set that equal to $0$, and solve for $\beta$ to find a minimizer:

$$
\begin{aligned}
&\mathbf{X}^\top(\mathbf{y} - \mathbf{X} \beta) = 0 \\
\implies
&\mathbf{X}^\top \mathbf{y} - \mathbf{X}^\top \mathbf{X} \beta = 0 \\
\implies
&\hat{\beta}_{OLS} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}
\end{aligned}
$$

This solution is unique if $\mathbf{X}^\top \mathbf{X}$ is non-singular. 

#### Generalized Least Squares
<a href="https://en.wikipedia.org/w/index.php?title=Generalized_least_squares&oldid=1292245629">Generalized least squares (GLS)</a> is a method to estimate the parameters in a linear regression model when there is correlation between the errors. Suppose we are in the linear regression setting as before except we assume that we have a <i>known</i> and non-singular conditional covariance matrix, $\Omega$, for the errors. That is:

$$
\text{Cov}(\epsilon \rvert \mathbf{X}) = \Omega
$$

GLS estimates $\beta$ through the following optimization problem:

$$
\begin{aligned}
\hat{\beta}_{GLS} &= \underset{\beta}{\arg \min} \left\{ (\mathbf{y} - \mathbf{X} \beta)^\top  \Omega^{-1} (\mathbf{y} - \mathbf{X}\beta) \right\} \\
&= \underset{\beta}{\arg \min} \left\{ -2 \beta^\top \mathbf{X}^\top \Omega^{-1} \mathbf{y} + \beta^\top \mathbf{X}^\top \Omega^{-1} \mathbf{X} \beta \right\}
\end{aligned}
$$

The above expression is quadratic in $\beta$, so we take the gradient, set it equal to zero, and solve for $\beta$:

$$
\begin{aligned}
0 &= -2\mathbf{X}^\top \Omega^{-1} \mathbf{y} + 2\mathbf{X}^\top \Omega^{-1} \mathbf{X} \beta   \\
\implies
\hat{\beta}_{GLS} &= (\mathbf{X}^\top \Omega{-1} \mathbf{X})^{-1} \mathbf{X}^\top \Omega{-1} \mathbf{y}
\end{aligned}
$$

Just like in OLS, the GLS estimator of $\beta$ is unbiased, and its conditional covariance matrix also has a similiar form:

$$
\mathbb{E}[\hat{\beta}_{GLS} \rvert \mathbf{X}] = \beta,
\hspace{8mm}
\text{Cov}(\hat{\beta}_{GLS} \rvert \mathbf{X}) = (\mathbf{X}^\top \Omega^{-1} \mathbf{X})^{-1}
$$

#### Weighted Least Squares
<a href="https://en.wikipedia.org/w/index.php?title=Weighted_least_squares&oldid=1279139018">Weighted least squares (WLS)</a> is equivalent to generalized least squares when $\Omega$ is a diagonal matrix. However, we'll discuss it in more detail here for completeness. 

The idea behind WLS is that we may have more or less confidence in the reliability of our observations, so we want to weight their contributions to the objective function during estimation. Let $\mathbf{W}$ be an $n \times n$ diagonal matrix of weights. The WLS goal is to find:

$$
\hat{\beta}_{WLS} = \underset{\beta}{\arg \min} \left\{ (\mathbf{y} - \mathbf{X} \beta)^\top \mathbf{W} (\mathbf{y} - \mathbf{X} \beta) \right\}
$$

Since this is quadratic in $\beta$, we take the derivative with respect to $\beta$. Solving the following for $\beta$ will give us our estimate:

$$
\begin{equation}
\label{eq:wls-obj}
\frac{\partial}{\partial \beta} \left[ (\mathbf{y} - \mathbf{X} \beta)^\top \mathbf{W} (\mathbf{y} - \mathbf{X} \beta) \right] = -2 \mathbf{X}^\top \mathbf{W} (\mathbf{y} - \mathbf{X}\beta) = 0
\end{equation}
$$

As we showed in the <a href="#generalized-least-squares">generalized least squares section</a>, the solution to this problem is to set:

$$
\hat{\beta}_{WLS} = (\mathbf{X}^\top \mathbf{W} \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{W} \mathbf{y}
$$

However, this solution isn't finished because we don't know what $\mathbf{W}$ is! By the same argument in the OLS case, $\hat{\beta}_{WLS}$ is unbiased for $\beta$ for any $\mathbf{W}$, so how can we determine what the best weight matrix is? 

One way is to minimize the (weighted) mean squared error. For any vector of weights $\lambda = (\lambda_0, \lambda_1, \dots, \lambda_m)^\top$, we want to find the $\mathbf{W}$ such that $\hat{\beta}_{WLS}$ minimizes the expression:

$$
\mathbb{E}\left[ (\lambda^\top(\hat{\beta}_{WLS} - \beta))^2 \right] = \text{Var}(\lambda^\top \hat{\beta}_{WLS})
$$

where the equality follows from the unbiasedness of the WLS estimator.

<!-- Let $\hat{\beta} = \Lambda\mathbf{y}$ be some other estimator where $\Lambda = (\mathbf{X}^\top (\mathbf{W} + )$

ADD PROOF OF WLS ESTIMATOR SOLUTION?? -->


### Likelihood
Since we have assumed a particular distribution for the errors, we can also estimate the coefficient vector with maximum likelihood (see <a href="/posts/2025/02/03/likelihood-theory.html">my likelihood post</a> for details). Using Eq. \eqref{eq:assumption-2}, our assumption is:

$$
\mathbb{E}\left[ \mathbf{y} \rvert \mathbf{X} \right] = \mathbf{X}\beta
$$

Since our observations are i.i.d., the (full) log-likelihood (technically conditional on $\mathbf{X}$...) is:

$$
\begin{aligned}
\ell(\mathbf{y}; \mathbf{X}, \beta) &= \log\left( \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(- \frac{(y_i - \mathbf{x}_i \beta)^2}{2 \sigma^2} \right)\right) \\
&= \frac{n}{\sqrt{2\pi \sigma^2}} - \frac{1}{2 \sigma^2}\sum_{i = 1}^n (y_i - \mathbf{x}_i \beta)^2
\end{aligned}
$$

Taking the derivative with respect to $\beta$, setting this equal to zero, and solving for $\beta$ gives us our maximum likelihood estimate:

$$
\begin{aligned}
\frac{\partial \ell(\mathbf{y}; \mathbf{X}, \beta)}{\partial \beta} &= - \frac{2}{2 \sigma^2}\sum_{i = 1}^n \mathbf{x}_i^\top (y_i - \mathbf{x}_i \beta) \\
&= 0  \\
\implies
\mathbf{X}^\top (\mathbf{y} - \mathbf{X} \beta) &= 0 \\
\implies 
\hat{\beta}_{MLE} &= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}
\end{aligned}
$$

We've just shown that the OLS and ML estimates of the parameters are equivalent!


### Uncertainty Quantification
Often a simple estimate of $\beta$ is not enough; we want to be able to say how confident we are in our estimation. This is usually done through confidence intervals, which require estimates of the mean and variance of our estimates. 

Let $\hat{\beta}$ denote the OLS (and ML since they are equivalent) estimator for $\beta$. The mean vector of $\hat{\beta}$ is given by:

$$
\begin{aligned}
\mathbb{E}\left[ \hat{\beta} \rvert \mathbf{X} \right] 
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}  \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \mathbf{y} \rvert \mathbf{X} \right]  \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \beta  \\
&= \beta
\end{aligned}
$$

<div class="theorem">
<strong>Claim (Estimator Covariance).</strong>
{% tabs covar-claim %}
{% tab covar-claim statement %}
We can similarly derive the covariance matrix of $\hat{\beta}$ as:
$$
\text{Cov}(\hat{\beta} \rvert \mathbf{X}) = \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}
$$
{% endtab %}
{% tab covar-claim proof %}
Let $\mathbf{X}$ be fixed, so we can drop the conditioning. Using the fact that $\mathbb{E}[\hat{\beta}] = \beta$, we see that:

$$
\begin{aligned}
\text{Cov}(\hat{\beta}) 
&= \mathbb{E}\left[ (\hat{\beta} - \beta) (\hat{\beta} - \beta)^\top \right]  \\
&= \mathbb{E}\left[(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \right] - 2 \mathbb{E}\left[ \hat{\beta} \right]\beta^\top + \beta \beta^\top\\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \mathbf{y} \mathbf{y}^\top \right] \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} - \beta \beta^\top \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ (\mathbf{X} \beta + \epsilon)(\mathbf{X} \beta + \epsilon)^\top \right] \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} - \beta \beta^\top\\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \beta \beta^\top \mathbf{X}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} - 2 (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \beta \mathbb{E}\left[ \epsilon \right] \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} + (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \epsilon \epsilon^\top \right]\mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}  - \beta \beta^\top \\
&= \beta \beta^\top  + (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top (\sigma^2 \mathbb{I}_{n \times n}) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}  - \beta \beta^\top \\
&= \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Putting the above results together, we see that the estimator has a Gaussian distribution:

$$
\hat{\beta} \sim \mathcal{N}(\beta, \sigma^2 (\mathbf{X}^\top \mathbf{X}){-1})
$$

We can use this information for additional inference about the parameters (e.g. hypothesis testing).

---

## Generalized Linear Models
Generalized linear models are not that different from classical linear regression models. In this case, we assume that $\mathbf{y}$ follow a distribution from an overdispersed exponential family (e.g. Poisson, Binomial, exponential, Gaussian, etc.).<d-cite key=mccullagh1989></d-cite>

<div id="exponential-family"></div>
<div class="definition">
  <strong>Definition (Overdispersed Exponential Family).</strong>
  <br>
    Let $\theta \in \mathbb{R}^{m \times 1}$ be a parameter vector. We say that the distribution of random vector $\mathbf{y} \in \mathbb{R}^{n \times 1}$ comes from an <i>overdispersed exponential family</i> if the probability density/mass function can be written as:
    $$
    \begin{equation}
    \label{eq:exp-fam}
    \begin{aligned}
    f_Y(\mathbf{y} \rvert \theta, \tau) = \prod_{i = 1}^n h(\mathbf{y}_i, \tau_i) \exp\left(\frac{\mathbf{b}(\theta_i) \mathbf{T}(\mathbf{y}_i) - A(\theta_i)}{d(\tau_i)} \right) 
    \end{aligned}
    \end{equation}
    $$
    <ul>
        <li>$\tau$ is the <i>dispersion parameter</i>, usually related to the variance somehow.</li>
        <li>$d(\tau_i)$ is some function of $\tau_i$ used to introduce the overdispersion, often taken to be $\frac{\tau}{\omega_i}$.</li>
        <li>$h(\mathbf{y}, \tau)$ is a non-negative function. </li>
        <li>The vector $\mathbf{b}(\theta)$ is a function of the parameter vector. If $\mathbf{b}(\theta) = \theta$, then we say the distribution is in <i>canonical form</i>. Any distribution can be written in canonical form by letting $\theta' = \theta$ and then using $\mathbf{b}(\theta') = \theta$.</li>
        <li>The vector $\mathbf{T}(\mathbf{y})$ is a sufficient statistic, meaning $f_Y(\mathbf{f} \rvert \mathbf{T}(\mathbf{y}))$ is independent of $\theta$.</li>
        <li>The scalar $A(\theta)$ is the <i>log-partition function</i> and is the natural logarithm of the normalization constant needed so that the above function integrates to one.</li>
    </ul>
    All $h(\mathbf{y}, \tau)$, $\mathbf{b}(\theta)$, $\mathbf{T}(\mathbf{y})$, $A(\theta)$, and $d(\tau)$ are assumed to be known. However, these functions are non-unique! For example, you can easily scale $\mathbf{T}(\theta)$ one by a constant and scale $\mathbf{b}(\theta)$ by the reciprocal.
</div>

### Model
We'll assume the response, $\mathbf{y}$, comes from an exponential family distribution so that the likelihood has the form in Eq. \eqref{eq:exp-fam}. We'll also take $\mathbf{T}(\mathbf{y}) = \mathbf{y}$, $\mathbf{b}(\theta) = \theta$, and we'll rewrite $h(\mathbf{y}, \tau)$ as $\exp(\log(h(\mathbf{y}, \tau)))$ so we can bring it inside the exponential function. Taking the natural logarithm gives us the log-likelihood:

$$
\begin{equation}
\label{eq:log-lik-glm}
\begin{aligned}
\ell(\mathbf{y}; \theta, \tau) &= \log(f_Y(\mathbf{y} \rvert \theta, \tau)) \\
&= \sum_{i = 1}^n \left[ \frac{\theta_i \mathbf{y}_i - A(\theta_i)}{d(\tau_i)} + \log(h(\mathbf{y}_i, \tau_i))\right]
\end{aligned}
\end{equation}
$$ 

Taking the first- and second-order derivatives of the log-likelihood with respect to the parameter $\theta$ gives us the gradient and Hessian:

$$
\begin{equation}
\label{eq:grad-hessian}
\begin{aligned}
\frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \theta} &= 
\begin{bmatrix}
\frac{\mathbf{y}_1 - \frac{d A(\theta_1)}{d \theta_1}}{d(\tau_1)} \\
\vdots \\
\frac{\mathbf{y}_n - \frac{d A(\theta_n)}{d \theta_1}}{d(\tau_n)} \\
\end{bmatrix} \\
&= \text{diag}\left( \frac{1}{d(\tau)} \right) \left[ \mathbf{y} - \frac{\partial A(\theta)}{\partial \theta}\right] 
\\
\frac{\partial^2 \ell(\mathbf{y}; \theta, \tau)}{\partial \theta \partial \theta^\top} &=
\text{diag}\left( \frac{1}{d(\tau)} \right)
\begin{bmatrix}
- \frac{d^2(A(\theta_1))}{d \theta_1 d \theta_1} & \dots & - \frac{d^2(A(\theta_1))}{d \theta_1 d \theta_n} \\
\vdots & \ddots & \vdots \\
- \frac{d^2(A(\theta_n))}{d \theta_n d \theta_1} & \dots & - \frac{d^2(A(\theta_n))}{d \theta_n d \theta_n}  
\end{bmatrix} \\
&= -\text{diag}\left( \frac{1}{d(\tau)} \right)\frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top}
\end{aligned}
\end{equation}
$$

where we let $\tau = (\tau_1, \dots, \tau_n)^\top$. 

<div class="theorem">
<strong>Claim (Mean and Variance).</strong>
{% tabs mean-var-glm %}
{% tab mean-var-glm statement %}
Using properties of the likelihood function, we can derive the mean (denoted by $\mu$) and variance of $\mathbf{y}$ as:

$$
\begin{equation}
\label{eq:mean-var-glm}
\begin{aligned}
\mu &= \mathbb{E}[\mathbf{y}] = \frac{\partial A(\theta)}{\partial \theta} \\
\text{Cov}(\mathbf{y}) &= \text{diag}\left( d(\tau) \right) \underbrace{\frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top}}_{= V(\mu)}
\end{aligned}
\end{equation}
$$
{% endtab %}
{% tab mean-var-glm proof %}
Under certain regularity conditions, which we'll assume hold here, the expectation of the score function (the gradient of the log-likelihood) should equal zero (see <a href="/blog/2025/likelihood-theory">my likelihood post</a> for a discussion on this). Thus:

$$
\begin{aligned}
&\mathbb{E}\left[ \frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \theta}  \right] = \mathbf{0} \\
\implies
&\text{diag}\left( \frac{1}{d(\tau)} \right) \left(\mathbb{E}[\mathbf{y}] - \frac{\partial A(\theta)}{\partial \theta} \right) = \mathbf{0} \\
\implies
&\mathbb{E}[\mathbf{y}] =  \frac{\partial A(\theta)}{\partial \theta}
\end{aligned}
$$

Under these conditions, we also have that the Fisher information (the variance of the score) is the negative expectation of the Hessian. Thus:

$$
\begin{aligned}
&\mathbb{E}\left[  \frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \theta}  \frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \theta^\top} \right] = -\mathbb{E}\left[\frac{\partial^2 \ell(\mathbf{y}; \theta, \tau)}{\partial \theta \partial \theta^\top}  \right] \\
\implies &\mathbb{E}\left[ \text{diag}\left( \frac{1}{d^2(\tau)} \right) \left(\mathbf{y}  - \frac{\partial A(\theta)}{\partial \theta}\right)\left(\mathbf{y} - \frac{\partial A(\theta)}{\partial \theta}\right)^\top\right] = - \mathbb{E}\left[ -\text{diag}\left( \frac{1}{d(\tau)} \right)\frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top} \right] \\
\implies & \text{diag}\left( \frac{1}{d^2(\tau)} \right)\mathbb{E}\left[ \left(\mathbf{y} -\mathbb{E}[\mathbf{y}]\right)\left(\mathbf{y} - \mathbb{E}[\mathbf{y}]\right)^\top\right] = \text{diag}\left( \frac{1}{d(\tau)} \right)\frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top} \\
\implies &\text{Cov}(\mathbf{y}) = \text{diag}\left(d(\tau)\right) \frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Eq. \eqref{eq:mean-var-glm} shows that the covariance matrix of $\mathbf{y}$ is the product of a function of the dispersion parameter and a function of the parameter vector $\theta$. In the literature, $\frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top}$ is often referred to as the <i>variance function</i> and denoted by $V(\mu)$ (where $\mu$ is the mean $\mathbb{E}[\mathbf{y}]$) since a function of the parameters is just a function of the mean as $\mathbf{X}$ is fixed. 

Arguably more importantly, the righthand side of $\mu = \frac{\partial A(\theta)}{\partial \theta}$ is some function of $\theta$, which we'll denote by $s^{-1}(\theta)$. This lets us link the mean with the canonical parameter through the expression:

$$
\begin{equation}
\label{eq:link-mean-param}
\theta = s(\mu)
\end{equation}
$$

Recall that we assumed that $\mathbb{E}\left[ g(\mathbf{y}) \rvert \mathbf{X} \right] = \mathbf{X} \beta$, which is equivalent to the assumption that $\mu = g^{-1}(\mathbf{X}\beta)$. Thus, we can relate the likelihood parameter, $\theta$, to the regression parameters, $\beta$, by:

$$
\begin{equation}
\label{eq:link-param-beta}
\theta = s(g^{-1}(\mathbf{X}\beta))
\end{equation}
$$

Now we have everything we need for our generalized linear model! 

### Estimation
The parameter estimates for generalized linear models are generally found through maximum likelihood methods. Since we assume a distribution from an exponential family, we can (usually) write a likelihood function as in Eq. \eqref{eq:log-lik-glm}:

$$
\ell(\mathbf{y}; \theta, \tau) = \sum_{i = 1}^n \left[ \frac{\theta_i^\top \mathbf{y}_i - A(\theta_i)}{d(\tau_i)} + \log(h(\mathbf{y}_i, \tau_i)) \right]
$$

As we usually do, we want to take the gradient with respect to the parameters of interest, $\beta$. However, this expression is not in terms of $\beta$, so we use the chain rule. The maximum likelihood estimate of $\beta$ is a solution to:

$$
\begin{equation}
\label{eq:score-beta}
\begin{aligned}
\frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \beta} &= \frac{\partial \theta}{\partial \beta} \frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \theta} = \mathbf{0} \\
\implies \frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \beta} &= \text{diag}\left( d(\tau) \right)\frac{\partial \theta}{\partial \beta} \left(\mathbf{y} - \frac{\partial A(\theta)}{\partial \theta}\right) = \mathbf{0} & \left(\text{Eq. }\eqref{eq:grad-hessian} \right) \\
\implies \frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \beta} &= 
\text{diag}\left( d(\tau) \right)\frac{\partial \theta}{\partial \beta} \left(\mathbf{y} - g^{-1}(\mathbf{X} \beta) \right) = \mathbf{0} & \left(\text{Eq. }\eqref{eq:mean-var-glm} \right)
\end{aligned}
\end{equation}
$$

<div class="theorem">
<strong>Claim (Derivatives).</strong>
{% tabs d-theta-d-beta %}
{% tab d-theta-d-beta statement %}
To find $\frac{\partial \theta}{\partial \beta}$, we use the chain rule again:

$$
\begin{equation}
\label{eq:d-theta-d-beta}
\frac{\partial \theta}{\partial \beta} = \frac{\partial \eta}{\partial \beta} \frac{\partial \theta}{\partial \eta} = \text{diag}\left( \frac{1}{d(\tau)} \right)\mathbf{X}^\top \text{Cov}^{-1}(\mathbf{y})\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-1}
\end{equation}
$$
{% endtab %}
{% tab d-theta-d-beta proof %}
The first derivative is easy:

$$
\frac{\partial \eta}{\partial \beta} = \frac{\partial}{\partial \beta} [\mathbf{X} \beta] = \mathbf{X}^\top
$$

The second is a bit trickier. By the chain rule:

$$
\frac{\partial \theta}{\partial \eta} = \frac{\partial \theta}{\partial \mu} \frac{\partial \mu}{\partial \eta}
$$

Using the fact that $\theta = s(\mu)$ and $\mu = \frac{\partial A(\theta)}{\partial \theta}$, we get:

$$
\begin{aligned}
&\theta = s(\mu) \\
\implies 
&\frac{\partial \theta}{\partial \theta} = \frac{\partial s(\mu)}{\partial \theta} = \frac{\partial s(\mu)}{\partial \mu} \frac{\partial \mu}{\partial \theta} \\
\implies
&\mathbb{I}_{n \times n} = \frac{\partial s(\mu)}{\partial \mu} \frac{\partial}{\partial \theta} \left[ \frac{\partial A(\theta)}{\partial \theta} \right] = \frac{\partial s(\mu)}{\partial \mu} \frac{\partial^2 A(\theta)}{\partial \theta \partial \theta^\top} 
\end{aligned}
$$

Using Eq. \eqref{eq:mean-var-glm}, we get:

$$
\begin{aligned}
&\mathbb{I}_{n \times n} = \frac{\partial s(\mu)}{\partial \mu}  \\
\implies
&\text{diag}\left( \frac{1}{d(\tau)} \right) \text{Cov}(\mathbf{y})  \\
\implies
&\text{diag}\left( d(\tau) \right) \text{Cov}^{-1}(\mathbf{y}) = \frac{\partial s(\mu)}{\partial \mu} = \frac{\partial \theta}{\partial \mu}
\end{aligned}
$$

We can do the same sort of trick for $\frac{\partial \mu}{\partial \eta}$:

$$
\begin{aligned}
&\mu = g^{-1}(\eta)  \\
\implies
&\eta = g(\mu) \\
\implies
&\frac{\partial \eta}{\partial \eta} = \frac{\partial g(\mu)}{\partial \eta} = \frac{\partial g(\mu)}{\partial \mu} \frac{\partial \mu}{\partial \eta} \\
\implies
&\mathbb{I}_{n \times n} = \frac{\partial g(\mu)}{\partial \mu} \frac{\partial \mu}{\partial \eta} \\
\implies 
&\frac{\partial \mu}{\partial \eta} = \left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-1}
\end{aligned}
$$

Putting the two results together, we get:
$$
\frac{\partial \theta}{\partial \eta} = 
\text{diag}\left( d(\tau) \right) \text{Cov}^{-1}(\mathbf{y}) \left[ \frac{\partial g(\mu)}{\mu} \right]^{-1}
$$
{% endtab %}
{% endtabs %}
</div>

Using the fact that each component of $g(\mu)$ uses only the corresponding component in $\mu$, we see that $\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-1}$ is diagonal, so the inverse is the matrix where we take the reciprocal of the diagonal elements. The same goes for $\text{Cov}^{-1}(\mathbf{y})$ since we assumed independent data (i.i.d. errors). 

Plugging Eq. \eqref{eq:d-theta-d-beta} into Eq. \eqref{eq:score-beta}, our objective becomes solving:

$$
\begin{equation}
\label{eq:new-obj}
\frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \beta} =  \mathbf{X}^\top \text{Cov}^{-1}(\mathbf{y})\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-1} \left(\mathbf{y} - g^{-1}(\mathbf{X} \beta) \right) = \mathbf{0}
\end{equation}
$$

Unfortunately, there may not be a closed form solution to the above equation since $g^{-1}(\cdot)$ may be non-linear. To get around this problem, we do a linear approximation of $g^{-1}(\mathbf{X}\beta)$ and optimize this instead.

Let $\hat{\beta}^{(t)}$ be a guess for the optimal value of $\beta$. Also let $$\eta^{(t)} = \mathbf{X}\beta^{(t)}$$ be the linear predictor evaluated at this value. We can then write the first order Taylor approximation of $$g^{-1}(\eta)$$ as:

$$
\begin{equation}
\label{eq:taylor-mu}
g^{-1}(\eta) 
\approx g^{-1}(\eta^{(t)}) + \frac{\partial g^{-1}(\eta)}{\partial \eta} \bigg\rvert_{\eta = \eta^{(t)}} (\eta - \eta^{(t)}) = g^{-1}(\eta^{(t)}) + \mathbf{W}^{(t)}(\eta - \eta^{(t)})
\end{equation}
$$

where $\mathbf{W}^{(t)} = \frac{\partial g^{-1}(\eta)}{\partial \eta} \bigg\rvert_{\eta = \eta^{(t)}}$ is the matrix of first order partial derivatives of $\mu$ with respect to $\eta$ evaluated at our guess for $\eta$, $\eta^{(t)}$. This matrix is diagonal since any component, $\mu_i$, of the mean vector is only a function of one corresponding component of the linear predictor $\eta_i$. 

We now use this approximation into Eq. \eqref{eq:new-obj}:

$$
\begin{aligned}
\frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \beta} &\approx \mathbf{X}^\top \text{Cov}^{-1}(\mathbf{y})\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-1}\bigg\rvert_{\eta = \eta^{(t)}} \left(\mathbf{y} - g^{-1}(\eta^{(t)}) + \mathbf{W}^{(t)} \left( \eta^{(t)} - \eta \right) \right) \\
&= \mathbf{X}^\top \text{Cov}^{-1}(\mathbf{y})\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-1}\bigg\rvert_{\eta = \eta^{(t)}} \mathbf{W}^{(t)} \left([\mathbf{W}^{(t)}]^{-1} (\mathbf{y} - g^{-1}(\eta^{(t)})) + \eta^{(t)} - \eta \right) \\
&= \mathbf{X}^\top \text{Cov}^{-1}(\mathbf{y})\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-2} \bigg\rvert_{\eta = \eta^{(t)}} \left(\frac{\partial g(\mu)}{\partial \mu}\bigg\rvert_{\eta = \eta^{(t)}} (\mathbf{y} - g^{-1}(\eta^{(t)})) + \eta^{(t)} - \eta \right) & \left(\frac{\partial \mu}{\partial \eta} \bigg\rvert_{\eta = \eta^{(t)}} = \frac{\partial g^{-1}(\eta)}{\partial \eta} \bigg\rvert_{\eta = \eta^{(t)}}  = \mathbf{W}^{(t)} = \left[\frac{\partial g(\mu)}{\partial \mu} \right]^{-1} \bigg\rvert_{\eta = \eta^{(t)}} \right)
\end{aligned}
$$

Letting $\mathbf{z}^{(t)} = \frac{\partial g(\mu)}{\partial \mu} \bigg\rvert_{\eta = \eta^{(t)}} (\mathbf{y} - g^{-1}(\eta^{(t)})) + \eta^{(t)}$, and defining $\tilde{\mathbf{W}}^{(t)} = \text{Cov}^{-1}(\mathbf{y})\left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-2}\bigg\rvert_{\eta = \eta^{(t)}}$, we can rewrite this as:

$$
\frac{\partial \ell(\mathbf{y}; \theta, \tau)}{\partial \beta} \approx \mathbf{X}^\top \tilde{\mathbf{W}}^{(t)}  \left(\mathbf{z}^{(t)} - \eta \right)
$$

Notice that setting the above equal to zero is exactly Eq. \eqref{eq:wls-obj} with $\mathbf{z}^{(t)}$ instead of $\mathbf{y}$, so we can use weighted least squares. Though we have a closed form solution, since we are using an approximation, we need to iterate the above process until convergence. So at iteration $t + 1$, we will set:

$$
\begin{aligned}
\eta^{(t+1)} &= \mathbf{X} \beta^{(t)} \\
\mu^{(t+1)} &= g^{-1} \eta^{(t+1)} \\
\mathbf{W}^{(t+1)} &= \text{Cov}^{-1}(\mathbf{y}) \left[ \frac{\partial g(\mu)}{\partial \mu} \right]^{-2} \bigg\rvert_{\eta = \eta^{(t+1)}} \\
\mathbf{z}^{(t+1)} &= \eta^{(t+1)} + \frac{ \partial g(\mu)}{\partial \mu}\bigg\rvert_{\eta = \eta^{(t+1)}}(\mathbf{y} - g^{-1}(\eta^{(t+1)})) \\
\beta^{(t+1)} &= (\mathbf{X}^\top \tilde{\mathbf{W}}^{(t+1)} \mathbf{X})^{-1} \mathbf{X}^\top \tilde{\mathbf{W}}^{(t+1)} \mathbf{z}^{(t+1)} 
\end{aligned}
$$

This process is called <a href="https://en.wikipedia.org/w/index.php?title=Iteratively_reweighted_least_squares&oldid=1279139010">iteratively reweighted least squares</a>.

### Uncertainty Quantification 
Due to the iterative nature of the estimation procedure, it is not so clear how to quantify our uncertainty in our estimates. Since there is no closed form solution to the MLEs, <a href=" https://statisticaloddsandends.wordpress.com/2020/11/20/variance-of-coefficients-for-linear-models-and-generalized-linear-models/">there is also no closed form for the variance of the MLEs</a>.

One can, instead, derive the asymptotic variance of the estimator using <a href="/posts/2025/02/03/likelihood-theory.html">likelihood theory</a>. As a maximum likelihood estimator, it should be asymptotically normal with mean $\beta$ and covariance matrix $-\left[ \mathbb{E}\left[ \frac{\partial^2 \ell(\mathbf{y}; \theta, \tau)}{\partial \beta \beta^\top} \right] \right]^{-1}$ <i>if the model is correct</i>.