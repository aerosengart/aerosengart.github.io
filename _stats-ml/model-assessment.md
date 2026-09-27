---
layout: distill
title: Model Assessment and Selection
description: A Primer
date: 2026-09-14
tabs: true
tags: regression models
toc:
  - name: Background
  - name: Model Assessment
  - name: Bias-Variance Trade-Off
bibliography: stats-ml.bib
---

Welcome to another post on basic (but important!) concepts in statistics and machine learning. This post is on model assessment and follows Chapter 7 in <i>The Elements of Statistical Learning</i><d-cite key=hastie2017></d-cite> very closely. I'll also add some of my own thoughts and maybe a bit from other references as well.

---

## Background
We will, in general, assume to have some training dataset, $\mathcal{T} = \{ (x_i, y_i) \}_{i = 1}^N $, which consists of points $(x,y) \overset{iid}{\sim} \mathcal{F}$ where $\mathcal{F}$ is some joint distribution. Each $x$ is a vector or scalar containing input information, and $y$ is the response variable that we are interested in learning about using the information in $x$. We will fit some prediction model, $\hat{f}$, to $\mathcal{T}$ and compare its predictions to the true $Y$ values using a <i>loss function</i>, which we denote with $\mathcal{L}(y, \hat{f}(x))$. 

In what follows, we will assume that $\mathcal{T}$ consists of $N$ training points. There is a slight abuse of notation where we use lowercase $x$ and $y$ to denote both arbitrary random variables from the joint distribution in question as well as realizations of them.

### Definitions
We'll start with several types of errors that quantify how well a model performs. First is the <i>training error</i>, which is easy to compute but doesn't tell us more about how the model will work in general.

<div class="definition">
<strong>Definition (Training Error).</strong> <br>
The <i>training error</i> of our model is defined as the average loss over a training sample:

$$
\text{Err}_{train} = \frac{1}{N} \sum_{i = 1}^N \mathcal{L}(y_i, \hat{f}(x_i))
$$
</div>

To get an idea of how the model performs for data we haven't seen before, we define the <i>test error</i>.

<div class="definition">
<strong>Definition (Test Error).</strong> <br>
The <i>test error</i> of our model is defined as the expected loss on a new point $(x^0, y^0) \sim \mathcal{F}$ for a fixed training set:

$$
\text{Err}_{\text{test}} = \mathbb{E}_{x^0, y^0}\left[ \mathcal{L}(y^0, \hat{f}(x^0)) \rvert \mathcal{T} \right]
$$
</div>

<aside><p>This is also called the <i>generalization error</i>.</p></aside>

This quantifies how <i>this particular fit of our model</i> is expected to perform on new data. Letting the training set vary yields the <i>expected test error</i>. 

<div class="definition">
<strong>Definition (Expected Test Error).</strong> <br>
The <i>expected test error</i> of our model is defined as the expected loss on a new point $(x^0, y^0) \sim \mathcal{F}$ for a fixed training set:

$$
\text{Err} = \mathbb{E}_{\mathcal{T}} \left[ \mathbb{E}_{x^0, y^0} \mathcal{L}(y^0, \hat{f}(x^0)) \rvert \mathcal{T} \right]
$$
</div>

<aside><p>This is also called the <i>expected prediction error</i>.</p></aside>

By the law of total expectation, $\text{Err} = \mathbb{E}[\text{Err}_{\text{test}}]$. In words, this means that the expected test error is the average test error across all possible training sets and testing points. 

This leads us to the <i>in-sample error</i>.

<div class="definition">
<strong>Definition (In-Sample Error).</strong> <br>
Let $y^0_i$ denote a new response observation for training point $x_i$. The <i>in-sample error</i> of our model is defined as:

$$
\text{Err}_{\text{in-sample}} = \frac{1}{N} \sum_{i = 1}^N \mathbb{E}_{y^0}\left[ \mathcal{L}(y^0_i, \hat{f}(x_i)) \rvert \mathcal{T} \right]
$$
</div>

In words, the in-sample error gives us a measure of how well this kind of model can be expected to perform given it was fit to $\mathcal{T}$, where the expectation is with respect to the possible response values associated with each training point.

<div class="definition">
<strong>Definition (Optimism).</strong> <br>
The <i>optimism</i> of our model is defined as:

$$
\text{optimism} = \text{Err}_{\text{in-sample}} - \text{Err}_{\text{train}}
$$

The <i>average optimism</i> of our model is:

$$
\omega = \mathbb{E}_{\mathbf{y}}[\text{Err}_{\text{in-sample}} - \text{Err}_{\text{train}}]
$$
</div>

The optimism of a model is the difference between the in-sample error and the training error. It basically quantifies how much worse (or better) our particular model fit is performing compared to what we expected for this type of model and this training set. The average optimism is simply the expectation of the optimism taken with respect to the responses, $\mathbf{y}$, in $\mathcal{T}$.

<aside><p>Note: the higher the optimism, the greater the difference between the expected in-sample and training errors.</p></aside> 

---

## Model Assessment
How do we judge whether a model is "good"? A reasonable choice is to consider how well we would expect it to perform on data like the ones we used to train it. Due to the relationship between the quantities we defined in the previous section, we can estimate the in-sample error by computing the training error and estimating the expected optimism:

$$
\hat{\text{Err}}_{\text{in-sample}}
 = \text{Err}_{\text{train}} + \hat{\omega} 
$$

We can do this estimation by saving some of the training data to use for this estimation, effectively making it "unseen". 

### Statistics
For certain situations, we can select the "best" model based upon a choice of statistic that is related to the expected optimism is some way. For now, we will consider the standard choices of loss  function(0-1 loss, squared error, entropy loss) For these, the following identity holds:

$$
\begin{aligned}
\label{eq:exp-opt}
\omega = \frac{2}{N} \sum_{i = 1}^N \text{Cov}(\hat{y}_i, y_i)
\end{aligned}
$$

The interpretation here is that the higher the covariance between the fitted and true values are, the higher the expected optimism. In other words, the more $y_i$ affects its fitted value, $\hat{y}_i$, the greater the difference between the training and in-sample errors.

For models linear in the predictors (with, say, $d$ parameters) and fit with squared loss, the expected optimism can be expressed as:

<aside><p>That is, models of the form $y = f(x) + \epsilon$.</p></aside>

$$
\begin{aligned}
\label{eq:lin-sqr-loss}
\omega = \frac{2}{N} \sum_{i = 1}^N \text{Cov}(\hat{y}_i, y_i) = \frac{2}{N} d \sigma_\epsilon^2
\end{aligned}
$$

This shows that the expected optimism increases with the number of parameters (with $d$) and decreases with the training set size (with $N$). Intuitively, this makes sense. As we increase the number of parameters, our model becomes more and more flexible, which sends the training error to zero. 

Our first estimate of the in-sample error is the <i>$C_p$ statistic</i>, which holds when we have a model that is linear in the predictors and fit using squared loss.

<div class="definition">
<strong>Definition ($C_p$ Statistic).</strong><br>
Suppose we have $y = f(x) + \epsilon$ with $\epsilon$ having mean zero and variance $\sigma^2_\epsilon$. We fit a model, $\hat{f}$, with $d$ parameters using squared loss. The $C_p$ statistic is given by:

$$
C_p = \hat{\text{Err}}_{\text{train}} + \frac{2d}{N} \hat{\sigma}^2_{\epsilon}
$$

where $\hat{\sigma}_{\epsilon}^2$ is a suitable estimate of $\sigma^2_{\epsilon}$.
</div>

If we instead fit the model by maximizing the log-likelihood, then we obtain the <i>Akaike information criterion</i>.

<div class="definition">
<strong>Definition (Akaike Information Criterion).</strong><br>
Suppose $y$ has a distribution parameterized by $\theta$. We fit a model, $\hat{f}$, by obtaining a maximum likelihood estimate of $\theta$, $\hat{\theta}$. Suppose the following holds as $N \rightarrow \infty$:

$$
-2 \mathbb{E}\left[ \log\left( P_{\hat{\theta}}(y) \right) \right] \approx - \frac{2}{N} \mathbb{E}\left[ \ell_N(\theta; y) \right] + \frac{2d}{N}
$$

where $P_{\hat{\theta}} \in \mathcal{P}_{\theta}$, the family of densities for $y$ and $\ell_N(\hat{\theta}; y) = \sum_{i = 1}^N \log\left( P_{\hat{\theta}}(y)\right)$ is the maximized log-likelihood function. The <i>Akaike information criterion</i> (AIC) is given by:

$$
\text{AIC} = - \frac{2}{N} \ell_N(\hat{\theta}; y)  + 2 \frac{d}{N}
$$
</div>

<aside><p>Note that this definition only holds for models with $d$ parameters. If we had a non-linear model, we would need to substitute a different measure of model complexity.</p></aside>

Similar to the AIC is the <i>Bayesian information criterion</i> (BIC). It differs slightly in its form but largely in its motivation/derivation.

<div class="definition">
<strong>Definition (Bayesian Information Criterion).</strong><br>
Suppose the same setting as in the definition of AIC. The <i>Bayesian information criterion</i> (BIC) is defined as:

$$
\text{BIC} = - 2 \ell_N(\hat{\theta}; y) + d \log(N) 
$$
</div>

<aside><p>The <strong>Schwarz criterion</strong> is given by $2 \text{BIC}$.</p></aside>

A nice property of the BIC is that it is <i>asymptotically consistent</i>, meaning as $N \rightarrow \infty$, it will select the correct model from a family of models with probability going to $1$. AIC does not exhibit this property and tends to select models with higher complexity. However, this property doesn't mean either one is better than the other. For finite samples, BIC can choose a model that is too simple. 

### Effective Number of Parameters
How do we measure the complexity of a model? For models that are linear in the predictors, one way is with the <i>effective number of parameters</i>.

<div class="definition">
<strong>Definition (Effective Number of Parameters).</strong> <br>
Let $\mathbf{y} = (y_1, y_2, \dots, y_N)^\top$ be the vector of responses in our training set, and let $\hat{\mathbf{y}}$ be the same for a model's predictions. Suppose our model is a linear fitting model; that is, it satisfies:

$$
\hat{\mathbf{y}} = \mathbf{S} \mathbf{y}
$$

where $\mathbf{S}$ is an $N \times N$ matrix that depends on the features $x_1, \dots, x_n$ but not on the responses. The <i>effective number of parameters</i> (also called the <i>effective degrees of freedom</i>) is defined as:

$$
df(\mathbf{S}) = \text{tr}\left[ \mathbf{S} \right]
$$
</div>

For the $C_p$ statistic, this is what we use for $d$!

---

## Bias-Variance Trade-Off
Let's say that we assumed that we are using the <strong>squared-error</strong> loss, $\mathcal{L}(Y, \hat{f}(X)) = (Y - \hat{f}(X))^2$. We'll also assume the following data-generating mechanism:

$$
Y = f(X) + \epsilon; \hspace{5mm} \mathbb{E}[\epsilon] = 0, \hspace{2mm} \text{Var}(\epsilon) = \sigma^2_{\epsilon}
$$

The expected testing error of some regression model $\hat{f}(X)$ at a given point $X = x_0$ can be decomposed as:

<!-- #region bd-trade-off -->
{% tabs bd-trade-off %}
{% tab bd-trade-off statement %}
$$
\begin{aligned}
\text{Err}(x_0)
    &= \sigma^2_\epsilon +  \left(f(x_0) - \mathbb{E}[\hat{f}(x_0)] \right)^2 + \text{Var}(\hat{f}(x_0))
\end{aligned}
$$
{% endtab %}
{% tab bd-trade-off proof %}
$$
\begin{aligned}
\text{Err}(x_0) 
    &= \mathbb{E}\left[ (Y - \hat{f}(x_0))^2 \right] \\
    &= \mathbb{E}\left[ Y^2 - 2 Y \hat{f}(x_0) + \hat{f}^2(x_0) \right] \\
    &= \mathbb{E}\left[ f^2(x_0) - 2f(x_0) \epsilon + \epsilon^2 \right] - 2 \mathbb{E}\left[(f(x_0) + \epsilon) \hat{f}(x_0) \right] + \mathbb{E}\left[ \hat{f}^2(x_0) \right] \\
    &= \mathbb{E}\left[f^2(x_0) \right] + \sigma^2_\epsilon - 2 \mathbb{E}\left[f(x_0) \hat{f}(x_0) \right] + \mathbb{E}\left[ \hat{f}^2(x_0)\right] \\
    &= \sigma^2_\epsilon + \mathbb{E}\left[ (f(x_0) - \hat{f}(x_0))^2\right]
\end{aligned}
$$

Adding and subtracting $\mathbb{E}[\hat{f}(x_0)]$ within the expectation in the last line yields:

$$
\begin{aligned}
\text{Err}(x_0)
    &= \sigma^2_\epsilon + \mathbb{E}\left[ (f(x_0) - \mathbb{E}[\hat{f}(x_0)] + \mathbb{E}[\hat{f}(x_0)] -  \hat{f}(x_0))^2 \right] \\
    &= \sigma^2_\epsilon + \underbrace{\mathbb{E}\left[ (f(x_0) - \mathbb{E}[\hat{f}(x_0)])^2 \right]}_{(a)}  + 2 \underbrace{\mathbb{E}\left[ (f(x_0) - \mathbb{E}[\hat{f}(x_0)])(\mathbb{E}[\hat{f}(x_0)] -  \hat{f}(x_0)) \right]}_{(b)} + \underbrace{\mathbb{E}\left[ (\mathbb{E}[\hat{f}(x_0)] -  \hat{f}(x_0))^2  \right]}_{(c)}
\end{aligned}
$$

Because $f(x_0)$ is a fixed quantity, we have:

$$
\begin{aligned}
(a) &= \mathbb{E}\left[ (f(x_0) - \mathbb{E}[\hat{f}(x_0)])^2 \right] \\
    &= \mathbb{E}\left[ f^2(x_0) - f(x_0) \mathbb{E}[\hat{f}(x_0)] + \left( \mathbb{E}[\hat{f}(x_0)] \right)^2 \right] \\
    &= f^2(x_0) - f(x_0)\mathbb{E}[\hat{f}(x_0)] + \left( \mathbb{E}[\hat{f}(x_0)] \right)^2 \\
    &= \left(f(x_0) - \mathbb{E}[\hat{f}(x_0)] \right)^2
$$

Furthermore, $\mathbb{E}[f(x_0)]$ is also fixed, so:

$$
\begin{align}
(b) &= \mathbb{E}\left[ (f(x_0) - \mathbb{E}[\hat{f}(x_0)])(\mathbb{E}[\hat{f}(x_0)] -  \hat{f}(x_0)) \right] \\
    &= \mathbb{E}\left[ f(x_0)\mathbb{E}[\hat{f}(x_0)]  - f(x_0) \hat{f}(x_0) - \left(\mathbb{E}[\hat{f}(x_0)]\right)^2 + \hat{f}(x_0) \mathbb{E}[\hat{f}(x_0)] \right] \\
    &= f(x_0) \mathbb{E}[\hat{f}(x_0)] - f(x_0) \mathbb{E}[\hat{f}(x_0)] - \left(\mathbb{E}[\hat{f}(x_0)]\right)^2 + \left( \mathbb{E}[\hat{f}(x_0)] \right)^2 \\
    &= 0
\end{align}
$$

Finally:

$$
\begin{align}
(c) &= \mathbb{E}\left[ (\mathbb{E}[\hat{f}(x_0)] -  \hat{f}(x_0))^2  \right] \\
    &= \mathbb{E}\left[ (\hat{f}(x_0) - \mathbb{E}[\hat{f}(x_0)] )^2 \right] \\
    &= \text{Var}(\hat{f}(x_0))
\end{align}
$$

Putting everything together:

$$
\begin{aligned}
\text{Err}(x_0)
    &= \sigma^2_\epsilon +  \left(f(x_0) - \mathbb{E}[\hat{f}(x_0)] \right)^2 + \text{Var}(\hat{f}(x_0))
\end{aligned}
$$
{% endtab %}
{% endtabs %}
<!-- #endregion -->

In the above decomposition, $\sigma^2_\epsilon$ is the <strong>irreducible error</strong> and represents the variation of $Y$ about its mean, $f(X)$. It is irreducible because we cannot get rid of this term no matter how good our model is. The second term, $\left(f(x_0) - \mathbb{E}[\hat{f}(x_0)]\right)^2$, is called the squared <strong>bias</strong>. It represents the difference between the average of our estimate, $\hat{f}(x_0)$, and the true mean, $f(x_0)$. The last term, $\text{Var}(\hat{f}(x_0))$, is the variance of our model. 

