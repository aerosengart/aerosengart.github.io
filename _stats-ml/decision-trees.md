---
layout: distill
title:  Decision Trees
description: A Primer
date: 2026-10-02
tabs: true
tags: trees ml prediction methods
bibliography: stats-ml.bib

---

In this post, I'm going to review <i>decision trees</i> (also called <i>classification and regression trees (CART)</i>). This will mostly follow Chapter 18 in <i>Probabilistic Machine Learning</i> by Kevin Murphy.<d-cite key=murphy2025></d-cite>

---

## Decision Trees

### Regression Trees
Regression trees are a particular type of decision tree where the features are exclusively real-valued. It works by partitioning the space of the feature space according to a series of simple decision rules. 

Suppose we have $N$ training samples consisting of feature vectors $\mathbf{x} = (x_1, \dots, x_d)$ with $x_i \in \mathbb{R}$ for $i = 1, \dots, d$ and some real-valued response, $y$. One such decision rule could be defined by whether $x_2 > t$ for some threshold $t$. This tells us whether we should move to the left or right branch of the tree. At the next branch, we find another decision rule that tells us which branch to follow. Thus, the tree consists of some number of nodes, each of which specifies a decision rule. Each leaf node specifies our decision (prediction, classification, etc.).

Because the rules have the form $x_i > t_i$, they are parallel to the axes in $\mathbb{R}^d$. If we have $M-1$ total nodes, then the tree defines $M$ regions, $R_1, \dots, R_M$, which partition the input space. Thus, in mathematical notation, a decision tree has the form:

$$
f(\mathbf{x}; \boldsymbol{\theta}) = \sum_{m = 1}^M w_m \mathbb{I}\left\{ \mathbf{x} \in R_m \right\}
$$

In the above, $R_m$ is the region defined by the $m$-th leaf node, $w_m$ is the decision specified by the $m$-th leaf node, and $$\boldsymbol{\theta} = \left\{ (R_m, w_m) : m = 1, \dots, M \right\}$$. Because the responses are also real-value, we take $w_m$ to be the average response value of the training samples that fall in $R_m$:

$$
w_m = \frac{\sum_{i = 1}^N y_n \mathbb{I}\left\{ \mathbf{x}_i \in R_m \right\}}{\sum_{i = 1}^N \mathbb{I} \left\{ \mathbf{x}_i \in R_m \right\}}
$$

### Classification Trees
Classification trees are essentially the same as regression trees but the nodes must be tailored to categorical features. For example, suppose $y$ is binary, $x_1$ specifies the shape of an object, and $x_2$ specifies its color. We could use "$x_1$ is a circle" and "$x_2$ is red" as two of our decision rules. 

<aside><p>Note: these are no longer axis aligned boundaries.</p></aside>

We must also change the leaf nodes to reflect the fact that we have binary output. A simple way is to take the majority vote of the training samples that fall in the given region. Another way is to assign labels by minimizing some training loss. 

### Tree Fitting
In general, fitting $f(\mathbf{x}; \boldsymbol{\theta})$ by minimizing some loss function is NP-complete. In practice, popular implementations use a greedy approach. The basic idea is that, at a given node, we want to minimize the error achieved by both the left and right sub-trees simultaneously. Let $\mathcal{D}_i = \{ (\mathbf{x}_l, y_l) \in N_i \}$ denote the collection of training samples that reach node $i$, $N_i$, and that we are splitting based on feature $j$.

If feature $j$ is real-valued, then we can define the possible left and right splits as $$\mathcal{D}^L_i(j, t) = \{ (\mathbf{x}_l, y_l) \in N_i : x_{l, j} \leq t \}$$ and $$\mathcal{D}^R_i(j, t) = \{ (\mathbf{x}_l, y_l) \in N_i: x_{l,j} > t \}$$ where $t$ is any possible threshold that divides the ordered list of unique values of the $j$-th feature across all $\mathbf{x}_l$ at $N_i$ into two parts. If feature $j$ is categorical, then we would choose one of the unique values taken on by the training data. This gives left and right splits defined as $$\mathcal{D}^L_i(j, t) = \{ (\mathbf{x}_l, y_l) \in N_i : x_{l, j} = t \}$$ and $$\mathcal{D}^R_i(j, t) = \{ (\mathbf{x}_l, y_l) \in N_i: x_{l,j} \neq t \}$$. 

To choose the best feature, $j_i$, and best threshold value, $t_i$, to split on, we optimize:

$$
(j_i, t_i) = \underset{j \in \{ 1, \dots, D \}}{\arg \min} \underset{t}{\min} \left\{  \frac{\rvert \mathcal{D}_i^L(j, t)}{\rvert \mathcal{D}_i \rvert} c(\mathcal{D}^L_i(j, t)) + \frac{\rvert \mathcal{D}_i^R(j, t)}{\rvert \mathcal{D}_i \rvert} c(\mathcal{D}_i^R(j, t))\right\} 
$$

for some choice of cost function, $c(\cdot)$. The optimal choice is then the minimizer over the sum of the error (cost) of the left sub-tree and right sub-tree over all possible features and thresholds at node $N_i$. The sum is weighted by the number of samples that fall into each sub-tree as well. 

For continuous $y$, we might use the mean squared error for $c(\cdot)$:

$$
c(\mathcal{D}_i) = \frac{1}{\rvert \mathcal{D}_i \rvert} \sum_{n \in \mathcal{D}_i}(y_n - \bar{y})^2; \hspace{5mm} \bar{y} = \frac{1}{\rvert \mathcal{D}_i \rvert} \sum_{n \in \mathcal{D}_i} y_n
$$

For categorical $y$, we might use the <strong>Gigi index</strong>. Let $\mathcal{Z}$ denote the set of class labels. The Gini index is defined as:

$$
G_i = \sum_{z \in \mathcal{Z}} \hat{\pi}_{i,z} (1 - \hat{\pi}_{i,z}); \hspace{5mm} \hat{\pi}_{i,z} = \frac{1}{\rvert \mathcal{D}_i \rvert} \sum_{n \in \mathcal{D}_i} \mathbb{I}\left\{ y_n = z \right\}
$$

Another choice is the <strong>entropy</strong> or <strong>deviance</strong>:

$$
H_i = \mathbb{H}(\hat{\boldsymbol{\pi}}_i) = - \sum_{z \in \mathcal{Z}} \hat{\pi}_{i,z} \log(\hat{\pi}_{i,z})
$$

<aside><p>We call a node with $H_i = 0$ <strong>pure</strong>.</p></aside>

We must be careful to not overfit! If we allow the tree to have as many nodes as possible, then we can achieve $0$ training loss by putting each training point in its own leaf (given no nois). One way to avoid this is to limit the depth of the tree beforehand or to grow to the maximum depth and then merge leaf nodes in a process called <strong>pruning</strong>.

---

## Ensemble Methods

Decision trees are really nice because they are interpretable, robust to outliers, can deal with both discrete and continuous data, pretty fast, and automatically do variable selection! They are also invariant under monotone transformations of the features (if you shift all the $x$'s up but maintain their order, then the splits don't change). However, they can be very sensitive to small variation in the input features and often have worse prediction performance to other machine learning methods. This is part of the motivation behind ensemble methods like random forests, model stacking, and boosting. 

In general, ensemble learning is where we fit many models and average over their outputs. If we have $f_m(y \vert \mathbf{x})$ for $m \in \mathcal{M}$, a collection of base models, then the ensemble model will be:

$$
f(y \rvert \mathbf{x}) = \frac{1}{\rvert \mathcal{M} \rvert} \sum_{m \in \mathcal{M}} f_m(y \rvert \mathbf{x})
$$

Alternatively, we can use a <i>committee</i> and take the majority vote of the base models (like with classification trees). In general, ensemble learners will have similar bias to the base models but lower variance. 

### Stacking
<i>Stacking (stacked generalization)</i> creates an ensemble model of the form:

$$
f(y \rvert \mathbf{x}) = \frac{1}{\rvert \mathcal{M} \rvert} \sum_{m \in \mathcal{M}} w_m f_m(y \rvert \mathbf{x})
$$

where the $w_m$ are weights that must be learned from some holdout set of data.

### Bagging
<i>Bagging (bootstrap aggregating)</i> creates an ensemble model whose base models come from fitting on different bootstrap samples of the training data (for more information on bootstrap, see <a href="/stats-ml/bootstrap">my post</a>). Bagging improves model robustness because the ensemble model does not rely so heavily on any one particular sample (recall decision trees can be quite volatile when it comes to small perturbations in the training data). However, for very stable models (ones that don't change so much when the training data is perturbed slightly), bagging does not really help much. Another limitation of bagging is that each base model only uses $63\%$ of the observations, on average, because the bootstrap samples are made <i>with replacement</i>. 

### Random Forests
<i>Random forests</i> basically creates an ensemble model in the same way as stacking, but each base model is only given a random subset of the training data and a random subset of the features to use for creating its nodes. This decorrelates the individual base models, which increasees their diversity and thereby the performance of the ensemble model.
