---
layout: distill
title: Shrinkage In Regression
description: 
date: 2026-02-16
tabs: true
tags: shrinkage regression
toc:
  - name: Background
  - name: Ridge Regression
    subsections:
        - name: Linear Algebra Perspective
        - name: Bayesian Perspective
  - name: LASSO Regression
    subsections:
        - name: Bayesian Perspective
        - name: Grouped LASSO
  - name: Elastic-Net Regression
  - name: Others
bibliography: stats-ml.bib
---

Sometimes when doing linear regression, we want to balance model complexity with performance to achieve a model that is useful but more interpretable.

<aside><p>That is, we favor <i>parsimony</i>.</p></aside>

One way to achieve a smaller model is with <strong>shrinkage</strong>, which is the practice of attenuating the effect of certain regressors through their coefficients. This is usually performed by adding a penalty to the objective function that <i>shrinks</i> the model coefficients towards zero.<d-cite key=hastie2017></d-cite>

--- 

## Background
Let's first consider the standard linear regression objective function (see <a href="/stats-ml/regression/">this post</a> for a more in-depth discussion of linear regression). Suppose our inputs are centered and scaled (so we do not have $\beta_0$). We choose the coefficient vector, $\hat{\beta}$, as the one that minimizes the residual sum of squares (i.e. the $\ell_2$-norm of the residuals):

$$
\begin{equation}
\label{eq:ols}
\hat{\beta} = \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 \right\}
\end{equation}
$$

If we want to penalize large coordinates in $\beta$, we can add a term on:

$$
\begin{equation}
\label{eq:shrink}
\hat{\beta}_{S} = \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 + \lambda \rvert \rvert \beta \rvert \rvert \right\}
\end{equation}
$$

In Eq. \eqref{eq:shrink}, $\lambda \geq 0$ is a tuning parameter, usually called the <strong>complexity</strong>, that controls how much shrinkage we desire. It is usually chosen by cross-validation (or something similar) because the optimal value is unknown. The term $\rvert \rvert \beta \rvert \rvert$ is some norm of $\beta$ and measures the length/magnitude of the vector. The choice of norm dictates what type of shrinkage we are performing.

<aside><p>It's important to note that the predictions will be dependent upon the scaling and centering of $\mathbf{y}$ and $\mathbf{X}$, which is why we have standardized the inputs.</p></aside>

---

## Ridge Regression
Suppose we choose the (squared) $\ell_2$-norm for our penalty:

$$
\begin{equation}
\label{eq:ridge}
\hat{\beta}_{R} = \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 + \lambda \rvert \rvert \beta \rvert \rvert_2^2 \right\}
\end{equation}
$$

This is called <strong>ridge regression</strong> and shrinks coefficients towards zero and each other.<d-cite key=hastie2017></d-cite>

<aside><p>This is also called <strong>Tikhonov regularization</strong>.</p></aside>

### Linear Algebra Perspective
The objective can be rewritten in a way that makes the restriction on the size of $\beta$ explicit:

$$
\begin{equation}
\label{eq:ridge-2}
\hat{\beta}_{R} = \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 + \lambda \left( \rvert \rvert \beta \rvert \rvert_2^2 - c\right) \right\}
\end{equation}
$$

where $c \geq 0$ is the constraint value. Taking the derivative of Eq. \eqref{eq:ridge-2} with respect to $\beta$, setting equal to zero, and solving yields the closed form solution:

$$
\begin{aligned}
\frac{\partial}{\partial \beta} \left[ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 + \lambda \left( \rvert \rvert \beta \rvert \rvert_2^2 - c\right) \right] 
&= \frac{\partial}{\partial \beta} \left[ (\mathbf{y} - \mathbf{X} \beta)^\top (\mathbf{y} - \mathbf{X} \beta) + \lambda \left( \beta^\top \beta - c \right) \right] \\
&= -2\mathbf{X}^\top(\mathbf{y} - \mathbf{X} \beta) + 2 \lambda \beta \\
\implies \mathbf{0}_p &= -2\mathbf{X}^\top(\mathbf{y} - \mathbf{X} \beta) + 2 \lambda \beta \\
\implies \lambda \beta &= \mathbf{X}^\top\mathbf{y} - \mathbf{X}^\top \mathbf{X} \beta \\
\implies \mathbf{X}^\top \mathbf{y} &= \left( \lambda \mathbb{I}_{p \times p} + \mathbf{X}^\top \mathbf{X}\right) \beta \\
\implies \hat{\beta}_R &= \left(\mathbf{X}^\top \mathbf{X} + \lambda \mathbb{I}_{p \times p}\right)^{-1} \mathbf{X}^\top \mathbf{y}
\end{aligned}
$$

Now, consider the ridge regression estimate of the model coefficients:

$$
\hat{\beta}_R = \left(\mathbf{X}^\top \mathbf{X} + \lambda \mathbb{I}_{p \times p}\right)^{-1} \mathbf{X}^\top \mathbf{y}
$$

This can be interpreted as a perturbation of the matrix $\mathbf{X}^\top \mathbf{X}$. In the event that this matrix is ill-conditioned, computing the OLS solution may be problematic because it involves $(\mathbf{X}^\top \mathbf{X})^{-1}$. By adding a positive constant to its main diagonal, we make it non-singular and thus avoid issues with inversion. 

### Bayesian Perspective
The ridge regression solution can alternatively be derived as the maximum a posteriori estimate when a Gaussian prior is assumed for $\beta$. Let's assume the following:

$$
\begin{aligned}
\mathbf{y}_i &\overset{iid}{\sim} \mathcal{N}\left(\mathbf{x}_i^\top \beta, \sigma^2 \right) \\
\beta_j &\overset{iid}{\sim} \mathcal{N}\left(0, \tau^2\right)
\end{aligned}
$$

Assuming $\tau^2$ and $\sigma^2$ are known, the posterior density of $\beta$ can be derived as:

$$
\begin{aligned}
f(\beta \rvert \mathbf{y}) 
&= \frac{f(\mathbf{y} \rvert \beta) f(\beta)}{f(\mathbf{y})} \\
&\propto \left[ \prod_{i = 1}^n \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(- \frac{\left( \mathbf{y}_i - \mathbf{x}_i^\top \beta \right)^2}{2 \sigma^2} \right) \right] \left[ \prod_{j = 1}^p \frac{1}{\sqrt{2 \pi \tau^2}} \exp\left(- \frac{\beta_j^2}{2 \tau^2} \right) \right] \\
&= \frac{1}{(2 \pi \sigma^2)^{\frac{n}{2}} (2 \pi \tau^2)^{\frac{p}{2}}} \exp\left( - \frac{1}{2 \sigma^2} \sum_{i = 1}^n \left(\mathbf{y}_i - \mathbf{x}_i^\top \beta \right)^2 - \frac{1}{2 \tau^2} \sum_{j = 1}^p \beta_j^2 \right)
\end{aligned}
$$

We can then take the derivative with respect to $\beta$ of the log-posterior density, we get:

$$
\begin{aligned}
\frac{\partial}{ \partial \beta} \left[ \log f(\beta \rvert \mathbf{y}) \right] &\propto  \frac{\partial}{\partial \beta} \left[ - \log\left( (2 \pi \sigma^2)^{\frac{n}{2}} (2 \pi \tau^2)^{\frac{p}{2}} \right) - \left[ \frac{1}{2 \sigma^2} \sum_{i = 1}^n \left(\mathbf{y}_i - \mathbf{x}_i^\top \beta \right)^2 + \frac{1}{2 \tau^2} \sum_{j = 1}^p \beta_j^2  \right]\right] \\
&= - \left(- \frac{2}{2\sigma^2} \sum_{i = 1}^n \mathbf{x}_i^\top \left(\mathbf{y}_i - \mathbf{x}_i^\top \beta\right) +  \frac{2}{2\tau^2} \beta\right) \\
&= \frac{1}{\sigma^2} \mathbf{X}^\top \left(\mathbf{y} - \mathbf{X} \beta \right) - \frac{1}{\tau^2} \beta
\end{aligned}
$$

Setting this equal to $\mathbf{0}_p$ and solvng for $\beta$ yields the posterior mode:

$$
\begin{aligned}
\mathbf{0}_p &= \frac{1}{\sigma^2} \mathbf{X}^\top \left(\mathbf{y} - \mathbf{X} \beta \right) - \frac{1}{\tau^2} \beta \\
\implies - \frac{1}{\tau^2} \beta &= - \frac{1}{\sigma^2} \mathbf{X}^\top \left(\mathbf{y} - \mathbf{X} \beta \right) \\
\implies -\frac{\sigma^2}{\tau^2} \beta &= -\mathbf{X}^\top \mathbf{y} + \mathbf{X}^\top \mathbf{X} \beta  \\
\implies \mathbf{X}^\top \mathbf{y} &= \left(\mathbf{X}^\top \mathbf{X} + \frac{\sigma^2}{\tau^2} \right) \beta \\
\implies \hat{\beta}_R &=  \left(\mathbf{X}^\top \mathbf{X} + \frac{\sigma^2}{\tau^2} \right)^{-1}\mathbf{X}^\top \mathbf{y} 
\end{aligned}
$$

We see that this is the same as the solution in Eq. \eqref{eq:ridge} but with $\lambda = \frac{\sigma^2}{\tau^2}$. Since the posterior distribution is also Gaussian, the mode equals the mean, so this is also the maximum a posteriori estimate. 

---

## LASSO Regression
We can alternatively choose the $\ell_1$-norm for our penalty:

$$
\begin{equation}
\label{eq:lasso}
\hat{\beta}_{L} = \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 + \lambda \rvert \rvert \beta \rvert \rvert_1 \right\}
\end{equation}
$$

This is called <strong>LASSO regression</strong>, which stands for <i>least absolute shrinkage and selection operator</i>.<d-cite key=hastie2017></d-cite> Eq. \eqref{eq:lasso} is equivalent to:

$$
\begin{equation}
\label{eq:lasso-2}
\begin{aligned}
\hat{\beta}_{L} &= \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 \right\} \\
&\text{subject to} \rvert \rvert \beta \rvert \rvert_1 \leq c
\end{aligned}
\end{equation}
$$

If $c$ is chosen to be small enough, then some coordinates of $\beta$ with be <i>shrunken to exactly zero</i>. This is in contrast to ridge regression which does not fully eliminate any coordinates. Because of this feature, LASSO can be used as a feature selection method. However, there is no closed form solution.

### Bayesian Perspective
The LASSO can also be though of from a Bayesian perspective. We have the same setting as with ridge regression but with a Laplace prior with location parameter $\mu = 0$ and scale parameter $\tau > 0$:

$$
f(\beta_j) = \frac{1}{2 \tau} \exp\left(- \frac{\rvert \beta_j \rvert}{\tau}\right)
$$

The LASSO estimate can then be derived as the posterior mode (but not necessarily the posterior mean because that requires a symmetric, unimodal distribution).


### Grouped LASSO
LASSO regression has also been extended to penalties that affect groups of predictors in the same way. Suppose we have $Q$ groups of predictors with $m_q$ predictors in group $q$. The <strong>grouped LASSO</strong> is the solution to:

$$
\begin{equation}
\label{eq:group-lasso}
\hat{\beta}_{GL} = \underset{\beta}{\arg \min} \left\{  \rvert \rvert \mathbf{y} - \sum_{q = 1}^Q \mathbf{X}_q \beta_q \rvert \rvert_2^2 = \lambda \sum_{q = 1}^Q \sqrt{m_q} \rvert \rvert \beta_q \rvert \rvert_2 \right\}
\end{equation}
$$

where $\mathbf{X}_q$ is the matrix containing the observed values for all samples of the predictors in group $q$. Since $\rvert \rvert \beta_q \rvert \rvert_2 = 0$ if, and only if, $\beta_q = \mathbf{0}$, this shrinkage method aims for group-wise sparsity in addition to shrinking individual coordinates towards $0$.<d-cite key=hastie2017></d-cite>


---

## Elastic-Net Regression
The <strong>elastic-net</strong> estimate is a compromise between ridge and LASSO regression:

$$
\begin{equation}
\hat{\beta}_{EN} =  \underset{\beta}{\arg\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{X}\beta \rvert \rvert_2^2 + \lambda \left(\alpha \rvert \rvert \beta \rvert \rvert_2^2 + (1 - \alpha) \rvert \rvert \beta \rvert \rvert_1 \right)\right\} 
\end{equation}
$$

Here, $\alpha$ is a tuning parameter that controls how much we lean towards ridge regression. That is, $\alpha = 1$ will yield $\hat{\beta}_R$ and $\alpha = 0$ will yield $\hat{\beta}_L$.<d-cite key=hastie2017></d-cite>


---

## Others
There are many other alternative shrinkage programs. Ridge and LASSO regression are both types of <strong>$\ell_p$ regularization</strong>, where the penalty term is some $\ell_p$ norm of $\beta$ (for $p > 0$). Another alternative is <a href="https://en.wikipedia.org/wiki/Least-angle_regression"><strong>least angle regression</strong></a>. 


