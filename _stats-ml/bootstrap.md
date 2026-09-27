---
layout: distill
title: Bootstrap
description: A Primer
date: 2026-09-20
tabs: true
tags: regression models resampling
toc:
  - name: Set-Up
  - name: Bootstrap
bibliography: stats-ml.bib
---

Welcome to another post on basic (but important!) concepts in statistics and machine learning. This post is on the bootstrap and follows Chapter 8 in <i>All of Statistics</i><d-cite key=wasserman2004></d-cite>.

---

## Set-Up
We will, in general, assume to have some training dataset, $\mathcal{T} = \{ (x_i, y_i) \}_{i = 1}^N$, which consists of points $(x,y) \overset{iid}{\sim} \mathcal{F}$ where $\mathcal{F}$ is some joint distribution. Each $x$ is a vector or scalar containing input information, and $y$ is the response variable that we are interested in learning about using the information in $x$. We will fit some prediction model, $\hat{f}$, to $\mathcal{T}$ and compare its predictions to the true $Y$ values using a <i>loss function</i>, which we denote with $\mathcal{L}(y, \hat{f}(x))$. 

In what follows, we will let $z_i = (x_i, y_i)$, the training points. There is a slight abuse of notation where we use lowercase $x$ and $y$ to denote both arbitrary random variables from the joint distribution in question as well as realizations of them.

---

## Boostrap
The bootstrap is a resampling method that aims to estimate aspects of a sampling distribution or model. Th basic idea behind it is that we use our training dataset (our sample) as a population and then take many samples from it. In essence, we sample <i>with replacement</i> $N$ points $B$ times to create $B$ bootstrap datasets. On each bootstrap dataset, we refit the model to investigate how the fit behaves over our $B$ datasets. 

This process can be adjusted and altered in order to improve the bootstrap results or extend the method to more complicated cases, etc.

### Conditional Test Error
For the case of estimating the conditional test error, we can combine the bootstrap idea with the cross-validation idea. The estimate is given by:

$$
\hat{\text{Err}}_{\text{LOOB}} = \frac{1}{N} \sum_{i = 1}^N \frac{1}{\rvert C^{-i} \rvert} \sum_{b \in C^{-i}} \mathcal{L}(y_i, \hat{f}^b(x_i))
$$

Here, $C^{-i}$ is the set of indices of the bootstrap datasets ($b = 1, \dots, N$) that <i>do not</i> contain observation $i$ from $\mathcal{T}$. $\hat{f}^b(x_i)$ is the predicted value of $y_i$ for observation $x_i$ using the model fit using bootstrap sample $b$. This avoids a data leakage problem that occurs when we simply take the average loss over all samples. However, this involves reducing the training set size, which introduces a different source of bias. There are adjustments that can be made to account for this ($0.632$ estimator, $0.632+$ estimator, etc.).

### Variance of an Estimator
For any statistic, $S(\mathcal{T})$, we can estimate aspects of its distribution by computing the empirical version over our bootstrap replicates. For example, if we are in a maximum likelihood setting, we could estimate the variance of the parameter of interest, $\theta$, with the empirical variance over the bootstrap samples:

$$
\hat{\text{Var}}(\theta) = \frac{1}{B - 1} \sum_{i = 1}^B (\hat{\theta}_i - \bar{\theta})^2
$$

where $\hat{\theta}_i$ is the estimate of $\theta$ made from boostrap sample $i$, and $$\bar{\theta} = \frac{1}{B} \sum_{i = 1}^B \hat{\theta}_i$$. 

One could also estimate a confidence interval for the estimator using bootstrapping. A fairly easy way is the <strong>bootstrap percentile interval</strong> where we use the lower and upper quantiles of the values computed over the bootstrap replicates:

$$
C = \left(\theta^*_{\frac{\alpha}{2}}, \theta^*_{1 - \frac{\alpha}{2}} \right)
$$

where $\theta^*_{s}$ is the $s$-th sample quantile of the bootstrap replicates. This can be justified with a monotone transformation-based argument. If one exists, call it $m(\cdot)$, and it is Gaussian with mean equal to $m(\theta)$, then the procedure above will yield an exact $1 - \alpha$ confidence interval (see pg. 116 of Wasserman<d-cite key=wasserman2004></d-cite>).

An alternative way is via <strong>pivotal quantities</strong>. Let $\hat{\theta}$ denote the estimate for $\theta$ made from the original sample and define $R = \hat{\theta} - \theta$. Let $H(\cdot)$ denote the CDF of $R$. We can then define:

$$
C^* = \left(\hat{\theta} - H^{-1}(1 - \alpha/2), \hat{\theta} - H^{-1}(\alpha/2) \right)
$$

This forms an <i>exact</i> $1 - \alpha$ confidence interval for $\theta$. 

<details>
<summary>Proof.</summary>
We have:

$$
\begin{aligned}
\mathbb{P}(a \leq \theta \leq b)
  &= \mathbb{P}(a - \hat{\theta} \leq \theta - \hat{\theta} \leq b - \hat{\theta}) \\
  &= \mathbb{P}(\hat{\theta} - b \leq \hat{\theta} - \theta \leq \hat{\theta} - a) \\
  &= \mathbb{P}(\hat{\theta} - b \leq R \leq \hat{\theta} - a) \\
  &= -H(\hat{\theta} - b) + H(\hat{\theta}-a) \\
  &= -H(\hat{\theta} - \hat{\theta} + H^{-1}(\alpha/2)) + H(\hat{\theta} - \hat{\theta} + H^{-1}(1 - \alpha/2)) \\
  &= -\frac{\alpha}{2} + 1 - \frac{\alpha}_2 \\
  &= 1 - \alpha
\end{aligned}
$$
</details>

Since we don't know $H(\cdot)$, we estimate it with the bootstrap as:

$$
\hat{H}(r) = \frac{1}{B} \sum_{b = 1}^B \mathbf{1}\left\{ R^*_b \leq r \right\}; \hspace{5mm} R^*_b = \hat{\theta}_b - \hat{\theta}
$$

Then we approximate the lower ($a$) and upper ($b$) limits of $C^*$ as:

$$
\begin{aligned}
\hat{a} &= \hat{\theta} - \hat{H}^{-1}(1 - \alpha/2) \\
        &= \hat{\theta} - r^*_{1 - \frac{\alpha}{2}} \\
        &= 2 \hat{\theta} - \theta^*_{1 - \frac{\alpha}{2}}
\hat{b} &= \hat{\theta} - \hat{H}^{_1}(\alpha / 2) \\
        &= \hat{\theta} - r^*_{\frac{\alpha}{2}} \\
        &= 2 \hat{\theta} - \theta^*_{\frac{\alpha}{2}}
\end{aligned}
$$

where $$r^*_s = \theta^*_{s} - \hat{\theta}$$ is the $s$-th sample quantile of the bootstrap pivotal quantities, and $$\theta^*_s$$ is the $s$-th sample quantile of the bootstrap estimates. 
