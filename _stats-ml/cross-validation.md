---
layout: distill
title: Cross-Validation
description: A Primer
date: 2026-09-19
tabs: true
tags: regression models resampling
toc:
  - name: Set-Up
  - name: $K$-Fold Cross-Validation
  - name: Notes
bibliography: stats-ml.bib
---

Welcome to another post on basic (but important!) concepts in statistics and machine learning. This post is on cross-validation nd follows Chapter 7 in <i>The Elements of Statistical Learning</i><d-cite key=hastie2017></d-cite>.

In my post on <a href="/stats-ml/model-assessment">model assessment</a>, I mentioned that we can split the training data into subsets so that we have "new" data with which we can evaluate performance. This post is on <i>cross-validation</i>, which is arguably the most popular method of doing this sample splitting. 

---

## Set-Up
We will, in general, assume to have some training dataset, $\mathcal{T} = \{ (x_i, y_i) \}_{i = 1}^N $, which consists of points $(x,y) \overset{iid}{\sim} \mathcal{F}$ where $\mathcal{F}$ is some joint distribution. Each $x$ is a vector or scalar containing input information, and $y$ is the response variable that we are interested in learning about using the information in $x$. We will fit some prediction model, $\hat{f}$, to $\mathcal{T}$ and compare its predictions to the true $Y$ values using a <i>loss function</i>, which we denote with $\mathcal{L}(y, \hat{f}(x))$. 

In what follows, we will assume that $\mathcal{T}$ consists of $N$ training points. There is a slight abuse of notation where we use lowercase $x$ and $y$ to denote both arbitrary random variables from the joint distribution in question as well as realizations of them.

---

## $K$-Fold Cross-Validation
Cross-validation, in general, aims to estimate the conditional test error:

$$
\text{Err}_{\mathcal{T}} = \mathbb{E}_{x^0, y^0}\left[ \mathcal{L}(y^0, \hat{f}(x^0)) \rvert \mathcal{T} \right] 
$$

which measures how well our model is expected to perform on previously unseen data. 

In the simplest case, we could split $\mathcal{T}$ into $$\mathcal{T}_{\text{train}}$$ and $$\mathcal{T}_{\text{test}}$$ containing, say, $70\%$ and $30\%$ of the $N$ points in $\mathcal{T}$, respectively. We would then fit $\hat{f}$ on $$\mathcal{T}_\text{train}$$ and then estimate $$\text{Err}_{\text{test}}$$ by computing the loss on $$\mathcal{T}_{\text{test}}$$. 

In $K$-fold cross validation, we increase the number of splits from $2$ to $K$. 

### Algorithm
Let $\kappa: \{1, 2, \dots, N \} \mapsto \{ 1, 2, \dots, K\}$ be an indexing function that maps from the training set sample indices to the fold indices. 

This will be used to partition $\mathcal{T}$ into $$\mathcal{T}_1, \dots, \mathcal{T}_K$$. Let $\hat{f}^{-k}(x)$ denote the model fitted to $\mathcal{T} \setminus \mathcal{T}_k$. 

The cross-validation estimate of the test error is given by:

$$
\text{CV}(\hat{f}) = \frac{1}{N} \sum_{i = 1}^N \mathcal{L}(y_i, \hat{f}^{-\kappa(i)}(x_i))
$$

### Choice of $K$
A quick and easy way to pick $K$ is simply to look at the size of $\mathcal{T}$ and select a value that would yield reasonably sized folds for training. Common choices are $5$ and $10$. In these cases, we can run into biased estimates of the expected prediction error due to the fact that we are fitting to much smaller training sets. However, we gain lower variance in the estimate because the training sets differ a good deal. 

Another option is to pick $K = N$, which is called <strong>leave-one-out cross-validation (LOOCV)</strong>. In this case, the training sets are so similar to each other (only one point differs!) that the estimate of the expected prediction error can be highly variable, though nearly unbiased. This can also be computationally expensive since it requires fitting the model $N$ times.

---

## Notes
Cross-validation can also be used to construct confidence intervals for parameter estimation. The same procedure is followed, but instead of computing $\hat{y}_i$, we would estimate the quantity of interest from each training subset, $$\mathcal{T} \setminus \mathcal{T}_k$$.

Another thing to note is that cross-validation can be a bit finicky. As we noted above, it comes with its own bias-variance trade-off that depends upon the size of $K$. There is also a distinction to be made about the quantity it is actually estimating. Cross-validation (see Section 7.12<d-cite key=hastie2017></d-cite>), generally, is used for approximating the <i>expected test error</i>. That is, it is not very reliably good for estimating $$\text{Err}_{\text{test}}$$, which is conditional on the training data $\mathcal{T}$.