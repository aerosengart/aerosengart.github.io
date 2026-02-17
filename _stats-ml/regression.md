---
layout: distill
title: Linear Regression
description: A Primer
date: 2026-01-19
tabs: true
tags: regression models
toc:
  - name: Background
  - name: Estimation
  - name: Prediction
  - name: Inference
  - name: Best Linear Unbiased Estimator
bibliography: stats-ml.bib
---

Linear regression is a method for predicting some outcome using information from other quantitative variables. It's one of the most important models in statistics/machine learning due to its ubiquity. Most of this post follows Chapter 3 in <i>The Elements of Statistical Learning</i>.<d-cite key=hastie2017></d-cite>

---

## Background
Let $X = (X_1, X_2, \dots, X_p)^\top$ be an input vector of $p$ <i>predictors</i>, also called <i>regressors</i> or <i>features</i>, and let $Y$ be a real-valued outcome variable. The goal of linear regression is to predict $Y$ using the information in $X$ by assuming that the mean of $Y$ is a function of the predictors with the form:

$$
\begin{equation}
\label{eq:lin-reg}
f(X) = \beta_0 + \sum_{j = 1}^p X_j \beta_j
\end{equation}
$$

The $\beta = (\beta_0, \beta_1, \dots, \beta_p)^\top$ are the unknown model parameters (called <i>coefficients</i>) which we must estimate by minimizing a chosen loss function over a sample (i.e. a <i>training set</i>). We'll denote a sample of $n$ observations of the predictors and outcome with ordered pairs $(\mathbf{x}_1, y_1), \dots, (\mathbf{x}_n, y_n)$ where $\mathbf{x}_i = (x_1, \dots, x_p)^\top$. 

We'll also use $\mathbf{X} = (\mathbf{1}_n, \mathbf{x}_1, \dots, \mathbf{x}_n)^\top$ to denote the $n \times (p + 1)$ of regressors (plus a prepended $n$-vector of ones) and $\mathbf{y} = (y_1, \dots, y_n)^\top$ to denote the $n$-vector of the sample observations.

---

## Estimation 
Usually, when someone refers to linear regression, they have estimated their parameters using <i>least squares</i>, which uses the <i>residual sum of squares (RSS)</i> as its loss function:

$$
\begin{equation}
\label{eq:rss}
\begin{aligned}
RSS(\beta) &= \sum_{i = 1}^n \left(y_i - f(\mathbf{x}_i)\right)^2 \\
           &= \sum_{i = 1}^n \left(y_i - \beta_0 - \sum_{j = 1}^p \mathbf{x}_{i,j} \beta_j \right)^2 \\
           &= (\mathbf{y} - \mathbf{X} \beta)^\top (\mathbf{y} - \mathbf{X} \beta)
\end{aligned}
\end{equation}
$$

Least squares makes no assumptions about the data distributions; it simply tries to minimize the average squared difference between the predicted and true values where the predicted values <i>must be</i> linear in the regressors. Since Eq. \eqref{eq:rss} is quadratic in $\beta$, we can (under certain conditions) minimize it by taking the gradient, setting that equal to zero, and solving for $\beta$:

$$
\begin{equation}
\begin{aligned}
\frac{\partial}{\partial \beta} \left[ RSS(\beta) \right] 
&= - 2 (\mathbf{y} - \mathbf{X} \beta)^\top \mathbf{X} = -2 \mathbf{X}^\top (\mathbf{y} - \mathbf{X} \beta) \\
\frac{\partial^2}{\partial \beta \partial \beta^\top} \left[ RSS(\beta) \right]
&= 2 \mathbf{X}^\top \mathbf{X}
\end{aligned}
\end{equation}
$$

Assuming that $\mathbf{X}^\top \mathbf{X}$ is invertible (i.e. the Hessian is positive-definite):

$$
\begin{aligned}
& &\frac{\partial}{\partial \beta} \left[ RSS(\beta) \right] &= \mathbf{0} \\
&\implies &-2 \mathbf{X}^\top (\mathbf{y} - \mathbf{X} \beta) &= \mathbf{0} \\
&\implies &\mathbf{X}^\top (\mathbf{y} - \mathbf{X} \beta) &= \mathbf{0} \\
&\implies &\mathbf{X}^\top \mathbf{y} &= \mathbf{X}^\top \mathbf{X} \beta \\
&\implies &(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} &= \beta
\end{aligned}
$$

---

## Prediction
The predicted values are then given by:

$$
\hat{\mathbf{y}} = \mathbf{X} \hat{\beta} = \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} 
$$

The matrix $\mathbf{H} = \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top$ is often called the <i>hat matrix</i> because it puts a "hat" on $\mathbf{y}$ to create the predicted value.

If we also notice that the least squares objective yields a predicted value that is the <i>orthogonal projection</i> of the outcome vector onto the <i>column space</i> of the predictors (i.e. the span of the columns of $\mathbf{X}$), we can see why $\mathbf{H}$ is also called the <i>projection matrix</i>.

---

## Inference
It's important to keep in mind that the least squares parameter estimates, $\hat{\beta}$, are functions of the sample, and we can therefore try to characterize its sampling distribution. 

We will make the following assumptions:

<ul>
    <li>$\mathbf{x}_1, \dots, \mathbf{x}_n$ are fixed</li>
    <li><strong>Independence</strong>: $y_1, \dots, y_n$ are uncorrelated</li>
    <li><strong>Homoscedasticity</strong>: $y_1, \dots, y_n$ have constant variance, $\sigma^2$</li>
    <li><strong>Linearity</strong>: $\mathbb{E}[Y \rvert X] = X\beta$</li>
    <li><strong>Gaussianity</strong>: $Y = \mathbb{E}[Y \rvert X] + \epsilon$ with $\epsilon \sim \mathcal{N}(0, \sigma^2)$</li>
</ul>

We can then derive the mean of $\hat{\beta}$:

<!-- #region mean-ls -->
<div class="theorem">
  <strong>Claim (Mean of Least Squares Estimate).</strong>
  <br>
{% tabs mean-ls %}
{% tab mean-ls statement %}
$$
\mathbb{E}[\hat{\beta}] = \beta
$$
{% endtab %}
{% tab mean-ls proof %}
Here, let $\epsilon = (\epsilon_1, \dots, \epsilon_n)^\top$ with $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$. The expectations below are taken conditional on $\mathbf{X}$ (i.e. with $\mathbf{X}$ fixed):

$$
\begin{aligned}
\mathbb{E}\left[ \hat{\beta} \right] 
&= \mathbb{E} \left[ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} \right] \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \mathbf{y} \right] & \left(\mathbf{X} \text{ fixed}\right) \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \mathbf{X} \beta + \epsilon \right] \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left( \mathbf{X} \beta + \mathbb{E}\left[ \epsilon \right] \right)
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X}\beta & \left(\mathbb{E}[\epsilon] = 0 \right) \\
&= \beta
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

We can also derive its variance-covariance matrix:

<!-- #region var-ls -->
<div class="theorem">
  <strong>Claim (Covariance of Least Squares Estimate).</strong>
  <br>
{% tabs var-ls %}
{% tab var-ls statement %}
$$
\text{Var}(\hat{\beta}) = \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}
$$
{% endtab %}
{% tab var-ls proof %}
Here, let $\epsilon = (\epsilon_1, \dots, \epsilon_n)^\top$ with $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$. The expectations below are taken conditional on $\mathbf{X}$ (i.e. with $\mathbf{X}$ fixed):

$$
\begin{aligned}
\text{Var} \left( \hat{\beta} \right)
&= \mathbb{E} \left[ \left(\hat{\beta} - \mathbb{E}[\hat{\beta}]\right) \left( \hat{\beta} - \mathbb{E}[\hat{\beta}]\right)^\top \right] \\
&= \mathbb{E} \left[ \left(\hat{\beta} - \beta \right)\left(\hat{\beta} - \beta\right)^\top  \right] \\
&= \mathbb{E}\left[ \hat{\beta} \hat{\beta}^\top - 2 \beta\hat{\beta}^\top  - \beta \beta^\top \right] \\
&= \mathbb{E}\left[ \left( (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right) \left( (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right)^\top  \right] -2 \beta \mathbb{E}\left[\hat{\beta}^\top\right] - \beta\beta^\top & \left(\text{linearity of expectation}\right) \\
&= \mathbb{E}\left[ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \right] - 2 \beta \beta^\top + \beta \beta^\top & \left(\text{previous proof}\right) \\
&= \mathbb{E}\left[  (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top (\mathbf{X} \beta + \epsilon) (\mathbf{X} \beta + \epsilon)^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \right] - \beta\beta^\top   \\
&= \mathbb{E}\left[  (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left(\mathbf{X} \beta \beta^\top \mathbf{X}^\top -2 \epsilon \beta^\top \mathbf{X}^\top + \epsilon \epsilon^\top \right) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \right] - \beta \beta^\top  \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left(\mathbf{X} \beta \beta^\top \mathbf{X}^\top - 2 \mathbb{E}\left[  \epsilon \right] \beta^\top \mathbf{X}^\top + \mathbb{E}\left[  \epsilon \epsilon^\top  \right] \right) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top - \beta \beta^\top  \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left(\mathbf{X} \beta \beta^\top \mathbf{X}^\top + \sigma^2 \mathbb{I}_{n \times n} \right) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top - \beta^\top \beta \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \beta \beta^\top \mathbf{X}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \sigma^2  (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{I}_{n \times n} \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top - \beta\beta^\top  \\
&= \beta \beta^\top + \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1} - \beta \beta^\top \\
&= \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Under the above assumptions and given the above derivations, we conclude that:

$$
\hat{\beta} \sim \mathcal{N}\left(\beta, \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}\right)
$$

We can also form an estimate of $\sigma^2$ as:

$$
\hat{\sigma}^2 = \frac{1}{n - p - 1} \sum_{i = 1}^n (y_i - \hat{y}_i)^2
$$

which is unbiased (shown below).

<!-- #region exp-var -->
<div class="theorem">
  <strong>Claim (Mean of Least Squares Estimate).</strong>
  <br>
{% tabs var-est %}
{% tab var-est statement %}
$$
\mathbb{E}[\hat{\sigma}^2] = \sigma^2
$$
{% endtab %}
{% tab var-est proof %}
Here, let $\epsilon = (\epsilon_1, \dots, \epsilon_n)^\top$ with $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$. Let $\tilde{\mathbf{x}}_i = (1, \mathbf{x}_i^\top)^\top$. The expectations below are taken conditional on $\mathbf{X}$ (i.e. with $\mathbf{X}$ fixed):

$$
\begin{aligned}
\mathbb{E}\left[ \hat{\sigma}^2 \right]
&= \mathbb{E}\left[ \frac{1}{n - p - 1} \sum_{i = 1}^n (y_i - \hat{y}_i)^2 \right] \\
&=  \frac{1}{n - p - 1}  \mathbb{E}\left[(\mathbf{y} - \hat{\mathbf{y}})^\top(\mathbf{y} - \hat{\mathbf{y}}) \right] \\
&=  \frac{1}{n - p - 1}  \mathbb{E}\left[\mathbf{y}^\top \mathbf{y} - 2\mathbf{y}^\top \hat{\mathbf{y}} + \hat{\mathbf{y}}^\top \hat{\mathbf{y}} \right] \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[\mathbf{y}^\top \mathbf{y} \right] - 2 \mathbb{E}\left[ \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right] + \mathbb{E}\left[ (\mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y})^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y}  \right] \right) \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[\mathbf{y}^\top \mathbf{y} \right] - 2 \mathbb{E}\left[ \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right]  + \mathbb{E}\left[ \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right] \right) \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[\mathbf{y}^\top \mathbf{y} \right] - 2 \mathbb{E}\left[ \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right]  + \mathbb{E}\left[ \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right] \right) \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[\mathbf{y}^\top \mathbf{y} \right] - \mathbb{E}\left[ \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right] \right) \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[(\mathbf{X} \beta + \epsilon)^\top (\mathbf{X}\beta + \epsilon) \right] -\mathbb{E}\left[ (\mathbf{X} \beta + \epsilon)^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top (\mathbf{X} \beta + \epsilon) \right] \right) \\
&=  \frac{1}{n - p - 1} \mathbb{E}\left[ (\mathbf{X} \beta + \epsilon)^\top\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] (\mathbf{X} \beta + \epsilon) \right] \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[ (\mathbf{X} \beta)^\top \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \beta \right] + 2 \mathbb{E}\left[ \epsilon^\top \left[ \mathbb{I}_{n \times n} - \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \beta \right] + \mathbb{E}\left[ \epsilon^\top \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \epsilon \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \beta^\top \mathbf{X}^\top \left[ \mathbb{I}_{n \times n} - \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \beta  + \mathbb{E}\left[ \text{tr}\left[ \epsilon^\top  \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \epsilon \right] \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \beta^\top \mathbf{X}^\top \left[ \mathbb{I}_{n \times n} - \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \beta  + \mathbb{E}\left[ \text{tr}\left[ \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \epsilon \epsilon^\top  \right] \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \beta^\top \mathbf{X}^\top \mathbf{X} \beta - \beta^\top \mathbf{X}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{X} \beta + \text{tr}\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top  \mathbb{E}\left[ \epsilon \epsilon^\top \right] \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \beta^\top \mathbf{X}^\top \mathbf{X} \beta - \beta^\top \mathbf{X}^\top \mathbf{X} \beta + \text{tr}\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top (\sigma^2 \mathbb{I}_{n \times n}) \right]  \right) 
&=  \frac{1}{n - p - 1} \text{tr}\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top (\sigma^2 \mathbb{I}_{n \times n}) \right] 
&=  \frac{1}{n - p - 1} \text{tr}\left[ \sigma^2 \mathbb{I}_{n \times n} -  \sigma^2 \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \\
&= \frac{1}{n - p - 1} \left( \sigma^2 \text{tr}\left[ \mathbb{I}_{n \times n} \right] -  \sigma^2  \text{tr}\left[ \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \right)\\
&= \frac{1}{n - p - 1} \left( n \sigma^2 - \sigma^2  \text{tr}\left[  (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\mathbf{X} \right] \right)\\
&= \frac{1}{n - p - 1} \left( n \sigma^2 - \sigma^2  \text{tr}\left[\mathbb{I}_{p + 1}\right] \right)\\
&= \frac{1}{n - p - 1} \left( n \sigma^2 - (p + 1)\sigma^2 \right) \\
&= \sigma^2
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

---

## Best Linear Unbiased Estimator
Recall that an unbiased estimator will have mean equal to the true parameter value, and a linear estimator will be a linear combination of the sample responses, $\mathbf{y}$. 

The least squares estimator is a linear, unbiased estimator. As we showed in the previous section, $\mathbb{E}[\hat{\beta}] = \beta$ (unbiased), and it is a linear estimator because, if we let $\mathbf{A} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top$, we have:

$$
\begin{aligned}
\hat{\beta} &= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} = \mathbf{A}\mathbf{y}
\end{aligned}
$$

The least squares estimate is considered the <i>best linear unbiased estimator (BLUE)</i> because it has the smallest variance of all linear unbiased estimates. This result is summarized in the <i>Gauss-Markov Theorem</i>.

<!-- #region gm-theorem -->
<div class="theorem">
  <strong>Gauss-Markov Theorem.</strong><d-cite key=gm2025></d-cite>
  <br>
  {% tabs gm-theorem %}
  {% tab gm-theorem statement %}
  For $\mathbf{y}, \epsilon \in \mathbb{R}^n$, $\mathbf{X} \in \mathbb{R}^{n \times (p + 1)}$, and $\beta \in \mathbb{R}^{p + 1}$, assume $\mathbf{y} = \mathbf{X} \beta + \epsilon$ where all $\epsilon_i$ are independent with mean $0$ and variance $\sigma^2$ (but are not necessarily Gaussian). 
  <br>
  Let $$\hat{\beta}_{OLS}$$ be the ordinary least squares estimator, and let $$\tilde{\beta} = \mathbf{A} \mathbf{y}$$ be some other unbiased linear estimator. The Gauss-Markov Theorem states that the OLS estimator minimizes the mean squared error criterion. That is, for any set of coefficients $$\lambda_1, \dots, \lambda_{p+1}$$:

  $$
    \underset{\hat{\beta}}{\arg \min} \left[ \mathbb{E}\left[ \left( \sum_{j = 1}^{p + 1} \lambda_j (\hat{\beta}_j - \beta_j) \right)^2 \right] \right] = \hat{\beta}_{OLS}
  $$

  which is equivalent to:

  $$  
    \text{Var}(\hat{\beta}) - \text{Var}(\hat{\beta}_{OLS})
  $$

  being positive semi-definite for all other linear unbiased estimators, $\hat{\beta}$. 

  <details>
  <summary>Proof.</summary>
  $$
  \begin{aligned}
    \mathbb{E}\left[ \left( \sum_{j = 1}^{p + 1} \lambda_j (\hat{\beta}_j - \beta_j) \right)^2\right] 
    &= \mathbb{E}\left[ \left((\hat{\beta} - \beta)^\top \lambda \right)^2 \right] \\
    &= \mathbb{E}\left[ \lambda^\top\left( \hat{\beta} - \mathbb{E}\left[ \hat{\beta} \right]\right)\left( \hat{\beta} - \mathbb{E}\left[ \hat{\beta} \right]\right)^\top \lambda \right] \\
    &= \lambda^\top \mathbb{E}\left[ \left( \hat{\beta} - \mathbb{E}\left[ \hat{\beta} \right]\right)\left( \hat{\beta} - \mathbb{E}\left[ \hat{\beta} \right]\right)^\top \right] \lambda \\
    &= \lambda^\top \text{Var}(\hat{\beta}) \lambda \\
    &= \text{Var}(\lambda^\top \hat{\beta})
  \end{aligned}
  $$
  </details>
  {% endtab %}
  {% tab gm-theorem proof %}
  Note that $\tilde{\beta}$ can be rewritten as:

  $$
  \mathbf{A} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}
  $$

  for some $(p + 1) \times n$ matrix, $\mathbf{D}$. Deriving the mean of $\tilde{\beta}$, which we know to be $\mathbf{0}_{p + 1}$:

  $$
  \begin{aligned}
  \mathbb{E}\left[ \tilde{\beta} \right]
  &= \mathbb{E}\left[ \mathbf{A} \mathbf{y} \right] \\
  &= \mathbb{E}\left[ \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D} \right) (\mathbf{X} \beta = \epsilon) \right] \\
  &= \mathbb{E}\left[ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \beta + (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top\epsilon + \mathbf{D} \mathbf{X} \beta + \mathbf{D} \epsilon \right] \\
  &= \mathbb{I}_{(p + 1) \times (p + 1)} \beta + (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[  \epsilon \right] + \mathbf{D} \mathbf{X} \beta + \mathbf{D} \mathbb{E}\left[ \epsilon \right] \\
  &= \beta + \mathbf{D} \mathbf{X} \beta \\
  &= (\mathbb{I}_{(p + 1) \times (p + 1)} + \mathbf{D}\mathbf{X}) \beta
  \end{aligned}
  $$

  Since we assumed that $\tilde{\beta}$ is unbiased:

  $$
  \begin{aligned}
  &\mathbb{E}\left[ \tilde{\beta} \right] = \mathbf{0}_{p + 1} \\
  \implies 
  &(\mathbb{I}_{(p + 1) \times (p + 1)} + \mathbf{D}\mathbf{X}) = \beta  \\
  \implies
  &\mathbf{D}\mathbf{X} = \mathbb{0}_{(p + 1) \times (p + 1)}
  \end{aligned}
  $$

  Now, consider the variance of $\tilde{\beta}$:

  $$
  \begin{aligned}
  \text{Var}(\tilde{\beta})
  &= \text{Var}(\mathbf{A}\mathbf{y}) \\
  &= \mathbf{A} \text{Var}(\mathbf{y}) \mathbf{A}^\top \\
  &= \mathbf{A} \left[\sigma^2 \mathbb{I}_{n \times n} \right] \mathbf{A}^\top \\
  &= \sigma^2 \mathbf{A}\mathbf{A}^\top \\
  &= \sigma^2 \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right)\left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right)^\top \\
  &= \sigma^2 \left[ \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right) \mathbf{X}(\mathbf{X}^\top \mathbf{X})^{-1} + \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right)\mathbf{D}^{-1} \right] \\
  &= \sigma^2 \left[(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X}(\mathbf{X}^\top \mathbf{X})^{-1} + \mathbf{D}\mathbf{X}(\mathbf{X}^\top \mathbf{X})^{-1} + (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{D}^\top + \mathbf{D}\mathbf{D}^\top \right] \\
  &= \sigma^2 \left[(\mathbf{X}^\top \mathbf{X})^{-1} +  \mathbf{0}_{(p + 1) \times (p + 1)}(\mathbf{X}^\top \mathbf{X})^{-1} + (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{0}_{(p + 1) \times (p + 1)} + \mathbf{D}\mathbf{D}^\top \right] \\
  &= \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1} + \sigma^2 \mathbf{D}\mathbf{D}^\top  \\
  &= \text{Var}(\hat{\beta}_{OLS}) + \sigma^2 \mathbf{D} \mathbf{D}^\top
  \end{aligned}
  $$
  
  Because $\sigma^2 > 0$ and $\mathbf{D} \mathbf{D}^\top$ is positive semi-definite, we have the desired result.
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endregion -->