---
layout: distill
title: Linear Regression
description: A Primer
date: 2026-01-19
tabs: true
tags: regression models
toc:
  - name: Background
  - name: Ordinary Least Squares
  - name: Inference
  - name: Best Linear Unbiased Estimator
bibliography: stats-ml.bib
---

Linear regression is a method for predicting some outcome using information from other quantitative variables. It's one of the most important models in statistics/machine learning due to its ubiquity. Most of this post follows Chapter 3 in <i>The Elements of Statistical Learning</i><d-cite key=hastie2017></d-cite> and the ordinary least square Wikipedia article.<d-cite key=ols2026></d-cite>

---

## Background
Let $X = (X_1, X_2, \dots, X_p)^\top$ be an input vector of $p$ <i>predictors</i>, also called <i>regressors</i> or <i>features</i>, and let $Y$ be a real-valued outcome variable. The goal of linear regression is to predict $Y$ using the information in $X$ by assuming that the mean of $Y$ is a function of the predictors with the form:

$$
\begin{equation}
\label{eq:lin-reg}
f(X) = \beta_0 + \sum_{j = 1}^p X_j \beta_j
\end{equation}
$$

The $\boldsymbol{\beta} = (\beta_0, \beta_1, \dots, \beta_p)^\top$ are the unknown model parameters (called <i>coefficients</i>) which we must estimate by minimizing a chosen loss function over a sample (i.e. a <i>training set</i>). We'll denote a sample of $n$ observations of the predictors and outcome with ordered pairs $(\mathbf{x}_1, y_1), \dots, (\mathbf{x}_n, y_n)$ where $\mathbf{x}_i = (x_1, \dots, x_p)^\top$. We'll also use $\mathbf{X} = (\mathbf{1}_n, \mathbf{x}_1, \dots, \mathbf{x}_n)^\top$ to denote the $n \times (p + 1)$ of regressors (plus a prepended $n$-vector of ones) and $\mathbf{y} = (y_1, \dots, y_n)^\top$ to denote the $n$-vector of the sample observations.

We do not directly observe the mean relationship; there is some error in each response that is unobservable, $\epsilon_i$, and our model becomes:

$$
y_i = \mathbf{x}_i^\top \boldsymbol{\beta} + \epsilon_i
$$

We will make the following assumptions:

<ul>
<li><strong>Exogeneity</strong>: $\mathbb{E}[\epsilon_i x_i] = 0$ for all $i$</li>
<li><strong>Homoscedasticity</strong>: $\mathbb{E}[\epsilon_i^2 \rvert x_i] = \sigma^2$ for all $i$</li>
<li><strong>Linearity</strong>: $\mathbb{E}[y_i \rvert \mathbf{x}_i] = \mathbf{x}_i^\top \boldsymbol{\beta}$ for some $\boldsymbol{\beta}$ for all $i$</li>
<li><strong>Uncorrelatedness</strong>: $\mathbb{E}[\epsilon_i \epsilon_j \rvert \mathbf{x}_i, \mathbf{x}_j] = 0$ for all $i, j$</li>
<li><strong>No Perfect Multicollinearity</strong>: $\mathbb{P}\left(\text{rank}(\mathbf{X}) = p \right) = 1$</li>
</ul>

---

## Ordinary Least Squares
Usually, when someone refers to linear regression, they have estimated their parameters using <i>least squares</i>, which uses the <i>residual sum of squares (RSS)</i> as its loss function:

$$
\begin{equation}
\label{eq:rss}
\begin{aligned}
RSS(\beta) &= \sum_{i = 1}^n \left(y_i - f(\mathbf{x}_i)\right)^2 \\
           &= \sum_{i = 1}^n \left(y_i - \beta_0 - \sum_{j = 1}^p \mathbf{x}_{i,j} \beta_j \right)^2 \\
           &= (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})^\top (\mathbf{y} - \mathbf{X} \boldsymbol{\beta}) \\
           &= \rvert \rvert \mathbf{y} - \mathbf{X} \boldsymbol{\beta} \rvert \rvert_2^2
\end{aligned}
\end{equation}
$$

### Optimization Perspective 
Least squares makes no assumptions about the data distributions; it simply tries to minimize the average squared difference between the predicted and true values where the predicted values <i>must be</i> linear in the regressors. Since Eq. \eqref{eq:rss} is quadratic in $\beta$, we can (under certain conditions) minimize it by taking the gradient, setting that equal to zero, and solving for $\beta$:

$$
\begin{equation}
\begin{aligned}
\frac{\partial}{\partial\boldsymbol{\beta}} \left[ RSS(\beta) \right] 
&= - 2 (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})^\top \mathbf{X} = -2 \mathbf{X}^\top (\mathbf{y} - \mathbf{X} \boldsymbol{\beta}) \\
\frac{\partial^2}{\partial \boldsymbol{\beta}\partial \boldsymbol{\beta}^\top} \left[ RSS(\beta) \right]
&= 2 \mathbf{X}^\top \mathbf{X}
\end{aligned}
\end{equation}
$$

Assuming that $\mathbf{X}^\top \mathbf{X}$ is invertible (i.e. the Hessian is positive-definite, which it will be under the "no perfect multicollinearity" assumption):

<aside><p>$\mathbf{X}^\top \mathbf{X}$ is called the <strong>Gram matrix</strong>.</p></aside>

$$
\begin{aligned}
& &\frac{\partial}{\partial\boldsymbol{\beta}} \left[ RSS(\beta) \right] &= \mathbf{0} \\
&\implies &-2 \mathbf{X}^\top (\mathbf{y} - \mathbf{X} \boldsymbol{\beta}) &= \mathbf{0} \\
&\implies &\mathbf{X}^\top (\mathbf{y} - \mathbf{X} \boldsymbol{\beta}) &= \mathbf{0} \\
&\implies &\mathbf{X}^\top \mathbf{y} &= \mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}\\
&\implies &(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} &=\boldsymbol{\beta}
\end{aligned}
$$

<aside><p>This solution is a unique minimizer!</p></aside>

### Prediction
The predicted values are then given by:

$$
\hat{\mathbf{y}} = \mathbf{X} \hat{\boldsymbol{\beta}} = \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} 
$$

The matrix $\mathbf{H} = \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top$ is often called the <i>hat matrix</i> because it puts a "hat" on $\mathbf{y}$ to create the predicted value.


### Geometric Perspective
As stated before, the OLS estimate is given by:

$$
\hat{\boldsymbol{\beta}}_{\text{OLS}} = \underset{\boldsymbol{\beta}}{\arg \min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X} \boldsymbol{\beta} \rvert \rvert^2_2 \right\}
$$

Notice that, for any value of $\hat{\boldsymbol{\beta}}$, we obtain the predicted values by computing:

$$
\hat{\mathbf{y}} = \mathbf{X} \hat{\boldsymbol{\beta}} = \begin{bmatrix}
    \sum_{j = 1}^p x_{1,j} \hat{\beta}_j \\
    \vdots \\
    \sum_{j = 1}^p x_{n, j} \hat{\beta}_j
\end{bmatrix}
= \hat{\beta}_1 \begin{bmatrix}
x_{1,1} \\
\vdots \\
x_{n, 1} 
\end{bmatrix}
+ \dots +
\hat{\beta}_p \begin{bmatrix}
x_{1, p} \\
\vdots \\
x_{n, p}
\end{bmatrix}
$$

The predicted values are linear combinations of the columns of $\mathbf{X}$ (i.e. they lie in the column space of $\mathbf{X}$). Thus, OLS finds the closest points in the column space of $\mathbf{X}$ to $\mathbf{y}$. 

To put it another way, since $\hat{\mathbf{y}} = \mathbf{X} \hat{\boldsymbol{\beta}}$ defines a hyperplane, and $\rvert \rvert \mathbf{y} - \hat{\mathbf{y}} \rvert \rvert_2^2$ is the Euclidean distance between $\mathbf{y}$ and $\hat{\mathbf{y}}$, OLS finds the $\hat{\boldsymbol{\beta}}$ which defines the hyperplane (based on the columns of $\mathbf{X}$) that is closest (in terms of Euclidean distance) to $\mathbf{y}$. Putting together the predicted values and the OLS solution, we see:

<aside><p>We call the difference between $y_i$ and $\hat{y}_i$ a <strong>residual</strong>.</p></aside>

$$
\begin{aligned}
\hat{\mathbf{y}}_{\text{OLS}} &= \mathbf{X} \hat{\boldsymbol{\beta}}_{\text{OLS}} 
= \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}
\end{aligned}
$$

The hat matrix $\mathbf{H}$ is also called the <i>projection matrix</i> because it projects $\mathbf{y}$ onto the column space of $\mathbf{X}$. 

### Statistical Estimation Perspective
Under the additional (and commonly made) assumption that $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$, the OLS estimator can be shown to be identical to the maximum likelihood estimator. 

Under the Gaussianity assumption, $\mathbf{y} \rvert \mathbf{X} \sim \mathcal{N}\left(\mathbf{X} \boldsymbol{\beta}, \sigma^2 \mathbf{I}_{n \times n}\right)$. The log-likelihood function is then:

$$
\begin{aligned}
\ell(\boldsymbol{\beta}, \sigma^2; \mathbf{y} \rvert \mathbf{X})
  &= - \frac{n}{2} \log(2 \pi) - \frac{1}{2} \rvert \sigma^2 \mathbb{I}_{n \times n} \rvert - \frac{1}{2} (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})^\top \left[ \sigma^2 \mathbb{I}_{n \times n} \right]^{-1}(\mathbf{y} - \mathbf{X} \boldsymbol{\beta})
\end{aligned}
$$

Differentiating with respect to $\boldsymbol{\beta}$:

$$
\begin{aligned}
\frac{\partial \ell(\boldsymbol{\beta}, \sigma^2; \mathbf{y} \rvert \mathbf{X})}{\partial \boldsymbol{\beta}}
  &= (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})^\top \left[ \sigma^2 \mathbb{I}_{n \times n} \right]^{-1} \mathbf{X} 
\end{aligned}
$$

Equating with $\mathbf{0}_p$ and solving for $\boldsymbol{\beta}$:

$$
\begin{aligned}
\mathbf{0}_p &= (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})^\top \left[ \sigma^2 \mathbb{I}_{n \times n} \right]^{-1} \mathbf{X} \\
\mathbf{0}_p &= \mathbf{y}^\top \left[ \sigma^2 \mathbb{I}_{n \times n} \right]^{-1} \mathbf{X} - \boldsymbol{\beta}^\top \mathbf{X}^\top \left[ \sigma^2 \mathbb{I}_{n \times n} \right]^{-1} \mathbf{X} \\
\mathbf{0}_p &= \frac{1}{\sigma^2} \mathbf{y}^\top \mathbf{X} - \frac{1}{\sigma^2} \boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} \\
\boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} &= \mathbf{y}^\top \mathbf{X} \\
\boldsymbol{\beta}^\top &= \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \\
\boldsymbol{\beta} &= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}
\end{aligned}
$$

---

## Inference
It's important to keep in mind that the least squares parameter estimates, $\hat{\boldsymbol{\beta}}$, are functions of the sample, and we can therefore try to characterize its sampling distribution. With the assumptions listed above, we can derive the mean of $\hat{\boldsymbol{\beta}}$:

<!-- #region mean-ls -->
<div class="theorem">
  <strong>Claim (Mean of Least Squares Estimate).</strong>
  <br>
{% tabs mean-ls %}
{% tab mean-ls statement %}
$$
\mathbb{E}[\hat{\boldsymbol{\beta}}] =\boldsymbol{\beta}
$$
{% endtab %}
{% tab mean-ls proof %}
Here, let $\epsilon = (\epsilon_1, \dots, \epsilon_n)^\top$ with $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$. The expectations below are taken conditional on $\mathbf{X}$ (i.e. with $\mathbf{X}$ fixed):

$$
\begin{aligned}
\mathbb{E}\left[ \hat{\boldsymbol{\beta}} \right] 
&= \mathbb{E} \left[ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} \right] \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \mathbf{y} \right] & \left(\mathbf{X} \text{ fixed}\right) \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[ \mathbf{X} \boldsymbol{\beta}+ \epsilon \right] \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left( \mathbf{X} \boldsymbol{\beta}+ \mathbb{E}\left[ \epsilon \right] \right) \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X}\boldsymbol{\beta}& \left(\mathbb{E}[\epsilon] = 0 \right) \\
&=\boldsymbol{\beta}
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
\text{Var}(\hat{\boldsymbol{\beta}}) = \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}
$$
{% endtab %}
{% tab var-ls proof %}
Here, let $\epsilon = (\epsilon_1, \dots, \epsilon_n)^\top$ with $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$. The expectations below are taken conditional on $\mathbf{X}$ (i.e. with $\mathbf{X}$ fixed):

$$
\begin{aligned}
\text{Var} \left( \hat{\boldsymbol{\beta}} \right)
&= \mathbb{E} \left[ \left(\hat{\boldsymbol{\beta}} - \mathbb{E}[\hat{\boldsymbol{\beta}}]\right) \left( \hat{\boldsymbol{\beta}} - \mathbb{E}[\hat{\boldsymbol{\beta}}]\right)^\top \right] \\
&= \mathbb{E} \left[ \left(\hat{\boldsymbol{\beta}} - \boldsymbol{\beta}\right)\left(\hat{\boldsymbol{\beta}} -\boldsymbol{\beta}\right)^\top  \right] \\
&= \mathbb{E}\left[ \hat{\boldsymbol{\beta}} \hat{\boldsymbol{\beta}}^\top - 2\boldsymbol{\beta}\hat{\boldsymbol{\beta}}^\top  - \boldsymbol{\beta}\boldsymbol{\beta}^\top \right] \\
&= \mathbb{E}\left[ \left( (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right) \left( (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{y} \right)^\top  \right] -2 \boldsymbol{\beta}\mathbb{E}\left[\hat{\boldsymbol{\beta}}^\top\right] - \boldsymbol{\beta}\boldsymbol{\beta}^\top & \left(\text{linearity of expectation}\right) \\
&= \mathbb{E}\left[ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} \mathbf{y}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \right] - 2 \boldsymbol{\beta}\boldsymbol{\beta}^\top + \boldsymbol{\beta}\boldsymbol{\beta}^\top & \left(\text{previous proof}\right) \\
&= \mathbb{E}\left[  (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top (\mathbf{X} \boldsymbol{\beta}+ \epsilon) (\mathbf{X} \boldsymbol{\beta}+ \epsilon)^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \right] - \boldsymbol{\beta}\boldsymbol{\beta}^\top   \\
&= \mathbb{E}\left[  (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left(\mathbf{X} \boldsymbol{\beta}\boldsymbol{\beta}^\top \mathbf{X}^\top -2 \epsilon \boldsymbol{\beta}^\top \mathbf{X}^\top + \epsilon \epsilon^\top \right) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \right] - \boldsymbol{\beta}\boldsymbol{\beta}^\top  \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left(\mathbf{X} \boldsymbol{\beta}\boldsymbol{\beta}^\top \mathbf{X}^\top - 2 \mathbb{E}\left[  \epsilon \right] \boldsymbol{\beta}^\top \mathbf{X}^\top + \mathbb{E}\left[  \epsilon \epsilon^\top  \right] \right) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top - \boldsymbol{\beta}\boldsymbol{\beta}^\top  \\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \left(\mathbf{X} \boldsymbol{\beta}\boldsymbol{\beta}^\top \mathbf{X}^\top + \sigma^2 \mathbb{I}_{n \times n} \right) \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top - \boldsymbol{\beta}^\top \boldsymbol{\beta}\\
&= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}\boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \sigma^2  (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{I}_{n \times n} \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top - \boldsymbol{\beta}\boldsymbol{\beta}^\top  \\
&= \boldsymbol{\beta}\boldsymbol{\beta}^\top + \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1} - \boldsymbol{\beta}\boldsymbol{\beta}^\top \\
&= \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Under the above assumptions and given the above derivations, we conclude that:

$$
\hat{\boldsymbol{\beta}} \sim \mathcal{N}\left(\beta, \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1}\right)
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
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[(\mathbf{X} \boldsymbol{\beta}+ \epsilon)^\top (\mathbf{X}\boldsymbol{\beta}+ \epsilon) \right] -\mathbb{E}\left[ (\mathbf{X} \boldsymbol{\beta}+ \epsilon)^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top (\mathbf{X} \boldsymbol{\beta}+ \epsilon) \right] \right) \\
&=  \frac{1}{n - p - 1} \mathbb{E}\left[ (\mathbf{X} \boldsymbol{\beta}+ \epsilon)^\top\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] (\mathbf{X} \boldsymbol{\beta}+ \epsilon) \right] \\
&=  \frac{1}{n - p - 1} \left( \mathbb{E}\left[ (\mathbf{X} \boldsymbol{\beta})^\top \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \boldsymbol{\beta}\right] + 2 \mathbb{E}\left[ \epsilon^\top \left[ \mathbb{I}_{n \times n} - \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \boldsymbol{\beta}\right] + \mathbb{E}\left[ \epsilon^\top \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \epsilon \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \boldsymbol{\beta}^\top \mathbf{X}^\top \left[ \mathbb{I}_{n \times n} - \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \boldsymbol{\beta} + \mathbb{E}\left[ \text{tr}\left[ \epsilon^\top  \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \epsilon \right] \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \boldsymbol{\beta}^\top \mathbf{X}^\top \left[ \mathbb{I}_{n \times n} - \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \mathbf{X} \boldsymbol{\beta} + \mathbb{E}\left[ \text{tr}\left[ \left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top\right] \epsilon \epsilon^\top  \right] \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}- \boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}+ \text{tr}\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top  \mathbb{E}\left[ \epsilon \epsilon^\top \right] \right]  \right) \\
&=  \frac{1}{n - p - 1} \left( \boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}- \boldsymbol{\beta}^\top \mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}+ \text{tr}\left[ \mathbb{I}_{n \times n} -  \mathbf{X} (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top (\sigma^2 \mathbb{I}_{n \times n}) \right]  \right) 
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

The least squares estimator is a linear, unbiased estimator. As we showed in the previous section, $\mathbb{E}[\hat{\boldsymbol{\beta}}] =\boldsymbol{\beta}$ (unbiased), and it is a linear estimator because, if we let $\mathbf{A} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top$, we have:

$$
\begin{aligned}
\hat{\boldsymbol{\beta}} &= (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y} = \mathbf{A}\mathbf{y}
\end{aligned}
$$

The least squares estimate is considered the <i>best linear unbiased estimator (BLUE)</i> because it has the smallest variance of all linear unbiased estimates. This result is summarized in the <i>Gauss-Markov Theorem</i>.

<!-- #region gm-theorem -->
<div class="theorem">
  <strong>Gauss-Markov Theorem.</strong><d-cite key=gm2025></d-cite>
  <br>
  {% tabs gm-theorem %}
  {% tab gm-theorem statement %}
  For $\mathbf{y}, \epsilon \in \mathbb{R}^n$, $\mathbf{X} \in \mathbb{R}^{n \times (p + 1)}$, and $\boldsymbol{\beta}\in \mathbb{R}^{p + 1}$, assume $\mathbf{y} = \mathbf{X} \boldsymbol{\beta}+ \epsilon$ where all $\epsilon_i$ are independent with mean $0$ and variance $\sigma^2$ (but are not necessarily Gaussian). 
  <br>
  Let $$\hat{\boldsymbol{\beta}}_{OLS}$$ be the ordinary least squares estimator, and let $$\tilde{\boldsymbol{\beta}} = \mathbf{A} \mathbf{y}$$ be some other unbiased linear estimator. The Gauss-Markov Theorem states that the OLS estimator minimizes the mean squared error criterion. That is, for any set of coefficients $$\lambda_1, \dots, \lambda_{p+1}$$:

  $$
    \underset{\hat{\boldsymbol{\beta}}}{\arg \min} \left[ \mathbb{E}\left[ \left( \sum_{j = 1}^{p + 1} \lambda_j (\hat{\beta}_j - \beta_j) \right)^2 \right] \right] = \hat{\boldsymbol{\beta}}_{OLS}
  $$

  which is equivalent to:

  $$  
    \text{Var}(\hat{\boldsymbol{\beta}}) - \text{Var}(\hat{\boldsymbol{\beta}}_{OLS})
  $$

  being positive semi-definite for all other linear unbiased estimators, $\hat{\boldsymbol{\beta}}$. 

  <details>
  <summary>Proof.</summary>
  $$
  \begin{aligned}
    \mathbb{E}\left[ \left( \sum_{j = 1}^{p + 1} \lambda_j (\hat{\beta}_j - \beta_j) \right)^2\right] 
    &= \mathbb{E}\left[ \left((\hat{\boldsymbol{\beta}} -\boldsymbol{\beta})^\top \lambda \right)^2 \right] \\
    &= \mathbb{E}\left[ \lambda^\top\left( \hat{\boldsymbol{\beta}} - \mathbb{E}\left[ \hat{\boldsymbol{\beta}} \right]\right)\left( \hat{\boldsymbol{\beta}} - \mathbb{E}\left[ \hat{\boldsymbol{\beta}} \right]\right)^\top \lambda \right] \\
    &= \lambda^\top \mathbb{E}\left[ \left( \hat{\boldsymbol{\beta}} - \mathbb{E}\left[ \hat{\boldsymbol{\beta}} \right]\right)\left( \hat{\boldsymbol{\beta}} - \mathbb{E}\left[ \hat{\boldsymbol{\beta}} \right]\right)^\top \right] \lambda \\
    &= \lambda^\top \text{Var}(\hat{\boldsymbol{\beta}}) \lambda \\
    &= \text{Var}(\lambda^\top \hat{\boldsymbol{\beta}})
  \end{aligned}
  $$
  </details>
  {% endtab %}
  {% tab gm-theorem proof %}
  Note that $\tilde{\boldsymbol{\beta}}$ can be rewritten as:

  $$
  \mathbf{A} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}
  $$

  for some $(p + 1) \times n$ matrix, $\mathbf{D}$. Deriving the mean of $\tilde{\boldsymbol{\beta}}$, which we know to be $\mathbf{0}_{p + 1}$:

  $$
  \begin{aligned}
  \mathbb{E}\left[ \tilde{\boldsymbol{\beta}} \right]
  &= \mathbb{E}\left[ \mathbf{A} \mathbf{y} \right] \\
  &= \mathbb{E}\left[ \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D} \right) (\mathbf{X} \boldsymbol{\beta}= \epsilon) \right] \\
  &= \mathbb{E}\left[ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X} \boldsymbol{\beta}+ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top\epsilon + \mathbf{D} \mathbf{X} \boldsymbol{\beta}+ \mathbf{D} \epsilon \right] \\
  &= \mathbb{I}_{(p + 1) \times (p + 1)} \boldsymbol{\beta}+ (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbb{E}\left[  \epsilon \right] + \mathbf{D} \mathbf{X} \boldsymbol{\beta}+ \mathbf{D} \mathbb{E}\left[ \epsilon \right] \\
  &= \boldsymbol{\beta}+ \mathbf{D} \mathbf{X} \boldsymbol{\beta}\\
  &= (\mathbb{I}_{(p + 1) \times (p + 1)} + \mathbf{D}\mathbf{X})\boldsymbol{\beta}
  \end{aligned}
  $$

  Since we assumed that $\tilde{\boldsymbol{\beta}}$ is unbiased:

  $$
  \begin{aligned}
  &\mathbb{E}\left[ \tilde{\boldsymbol{\beta}} \right] = \mathbf{0}_{p + 1} \\
  \implies 
  &(\mathbb{I}_{(p + 1) \times (p + 1)} + \mathbf{D}\mathbf{X}) = \boldsymbol{\beta} \\
  \implies
  &\mathbf{D}\mathbf{X} = \mathbb{0}_{(p + 1) \times (p + 1)}
  \end{aligned}
  $$

  Now, consider the variance of $\tilde{\boldsymbol{\beta}}$:

  $$
  \begin{aligned}
  \text{Var}(\tilde{\boldsymbol{\beta}})
  &= \text{Var}(\mathbf{A}\mathbf{y}) \\
  &= \mathbf{A} \text{Var}(\mathbf{y}) \mathbf{A}^\top \\
  &= \mathbf{A} \left[\sigma^2 \mathbb{I}_{n \times n} \right] \mathbf{A}^\top \\
  &= \sigma^2 \mathbf{A}\mathbf{A}^\top \\
  &= \sigma^2 \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right)\left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right)^\top \\
  &= \sigma^2 \left[ \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right) \mathbf{X}(\mathbf{X}^\top \mathbf{X})^{-1} + \left((\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top + \mathbf{D}\right)\mathbf{D}^{-1} \right] \\
  &= \sigma^2 \left[(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{X}(\mathbf{X}^\top \mathbf{X})^{-1} + \mathbf{D}\mathbf{X}(\mathbf{X}^\top \mathbf{X})^{-1} + (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{X}^\top \mathbf{D}^\top + \mathbf{D}\mathbf{D}^\top \right] \\
  &= \sigma^2 \left[(\mathbf{X}^\top \mathbf{X})^{-1} +  \mathbf{0}_{(p + 1) \times (p + 1)}(\mathbf{X}^\top \mathbf{X})^{-1} + (\mathbf{X}^\top \mathbf{X})^{-1}\mathbf{0}_{(p + 1) \times (p + 1)} + \mathbf{D}\mathbf{D}^\top \right] \\
  &= \sigma^2 (\mathbf{X}^\top \mathbf{X})^{-1} + \sigma^2 \mathbf{D}\mathbf{D}^\top  \\
  &= \text{Var}(\hat{\boldsymbol{\beta}}_{OLS}) + \sigma^2 \mathbf{D} \mathbf{D}^\top
  \end{aligned}
  $$
  
  Because $\sigma^2 > 0$ and $\mathbf{D} \mathbf{D}^\top$ is positive semi-definite, we have the desired result.
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endregion -->