---
layout: distill
title: Principal Components Analysis
description: A Primer
date: 2026-03-02
tabs: true
tags: linear-algebra primer dimension-reduction
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Set-Up
    - name: PCA Solution
      subsections:
        - name: Induction
    - name: Maximizing Variance
bibliography: stats-ml.bib
---

Principal components analysis is one of the first dimensionality reduction methods taught in introductory machine learning/statistics courses. It's fairly simple and has numerous extensions in classification and prediction.

---

## Set-Up
We will assume to have some (high-dimensional, i.e. $d$ is large) features, and our dataset consists of $n$ observations $$\mathbf{x}_i \in \mathbb{R}^d$$. We'll collect these $$\mathbf{x}_1, \dots, \mathbf{x}_d$$ into an $n \times d$ matrix denoted by $\mathbf{X}$. We'll assume that the $\mathbf{x}$ have been centered so that $$\bar{\mathbf{x}} = \frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i = \mathbf{0}_d$$. 

Our goal is to project these $\mathbf{x}_i$ into some lower dimensional space, $p < < d$, that preserves as much of the "information" contained in the dataset as possible. We do this by assuming that a linear combination of basis functions is a suitable approximation to each $\mathbf{x}_i$. We will denote these basis functions are $$\mathbf{w}_1, \dots, \mathbf{w}_p \in \mathbb{R}^d$$, and the weights for the linear combination will be denoted by $$\mathbf{z}_1, \dots, \mathbf{z}_n \in \mathbb{R}^p$$. Our assumption can be stated as:

<aside><p>The weights are also called the <strong>latent vectors</strong> or <strong>latent factors</strong>.</p></aside>

$$
\mathbf{x}_i \approx \sum_{j = 1}^p \mathbf{z}_{i,j} \mathbf{w}_j
$$

We'll collect the $\mathbf{z}_i$ vectors and $\mathbf{w}_j$ vectors into matrices $\mathbf{Z} \in \mathbb{R}^{n \times p}$ and $\mathbf{W} \in \mathbb{R}^{p \times d}$, respectively. Since we could pick anything we wanted for the weights and basis vectors, we will narrow down our search. Principal components analysis aims to minimize the <strong>reconstruction error</strong>; that is, we want to minimize:

<aside><p>The name comes from the fact that we are <i>reconstructing</i> each $\mathbf{x}_i$ with $\mathbf{W} \mathbf{z}_i$.</p></aside>.

$$
\mathcal{L}(\mathbf{W}, \mathbf{Z}) = \frac{1}{n} \rvert \rvert \mathbf{X} - \mathbf{Z} \mathbf{W}^\top \rvert \rvert_F^2 = \frac{1}{n} \rvert \rvert \mathbf{X}^\top - \mathbf{W} \mathbf{Z}^\top \rvert \rvert_F^2 = \frac{1}{n} \sum_{i = 1}^n \rvert \rvert \mathbf{x}_i - \mathbf{W} \mathbf{z}_i \rvert \rvert^2_2
$$

where $\mathbf{W}$ is restricted to be orthogonal. 

---

## PCA Solution
Suppose, first, that $d = 1$, and let $$\tilde{\mathbf{z}}_1 = (\mathbf{z}_{1,1}, \mathbf{z}_{2,1}, \dots, \mathbf{z}_{n,1})^\top$$ be the vector of the weights associated with $\mathbf{w}_1$ for all of the observations in the dataset. Using the fact that $\mathbf{W}$ is assumed to be orthogonal, the optimization problem is:

$$
\begin{equation}
\label{eq:loss-1}
\begin{aligned}
\underset{\mathbf{w}_1, \tilde{\mathbf{z}}_1}{\arg \min} \left\{ \mathcal{L}(\mathbf{w}_1, \tilde{\mathbf{z}}_1) \right\} 
    &= \underset{\mathbf{w}_1, \tilde{\mathbf{z}}_1}{\arg \min} \left\{ 
        \frac{1}{n} \sum_{i = 1}^n \rvert \rvert \mathbf{x}_i - \tilde{\mathbf{z}}_{i,1}\mathbf{w}_1 \rvert \rvert^2_2 
    \right\} \\
    &= \underset{\mathbf{w}_1, \tilde{\mathbf{z}}_1}{\arg \min} \left\{ 
        \frac{1}{n} \sum_{i = 1}^n (\mathbf{x}_i - \tilde{\mathbf{z}}_{i,1}\mathbf{w}_1)^\top (\mathbf{x}_i - \tilde{\mathbf{z}}_{i,1} \mathbf{w}_1) 
    \right\} \\
    &= \underset{\mathbf{w}_1, \tilde{\mathbf{z}}_1}{\arg \min} \left\{ 
        \frac{1}{n} \sum_{i = 1}^n \left[ \mathbf{x}_i^\top \mathbf{x}_i -2 \tilde{\mathbf{z}}_{i,1} \mathbf{w}_1^\top \mathbf{x}_i + \tilde{\mathbf{z}}_{i,1}^2 \mathbf{w}_1^\top \mathbf{w}_1 \right] 
    \right\} \\
    &= \underset{\mathbf{w}_1, \tilde{\mathbf{z}}_1}{\arg \min} \left\{ 
        \frac{1}{n} \sum_{i = 1}^n \left[ \mathbf{x}_i^\top \mathbf{x}_i -2 \tilde{\mathbf{z}}_{i,1} \mathbf{w}_1^\top \mathbf{x}_i + \tilde{\mathbf{z}}_{i,1}^2 \right] 
    \right\} & \left(\mathbf{w}_1^\top \mathbf{1} = 1 \right)
\end{aligned}
\end{equation}
$$

Let's consider finding the optimal weights. Fixing $\mathbf{w}_1$, we take the derivative of the loss, set that equal to zero, and solve. 

$$
\begin{aligned}
\frac{\partial \mathcal{L}(\mathbf{w}_1, \tilde{\mathbf{z}}_1)}{\partial \tilde{\mathbf{z}}_{i,1}} 
    &= \frac{\partial}{\partial \mathbf{z}_{i,1}} \frac{1}{n} \sum_{j = 1}^n \mathbf{x}_j^\top \mathbf{x}_j - 2 \tilde{\mathbf{z}}_{j,1} \mathbf{w}_1^\top \mathbf{x}_j + \tilde{\mathbf{z}}_{j,1}^2
    = \frac{1}{n} \left(-2 \mathbf{w}_1^\top \mathbf{x}_i + 2 \tilde{\mathbf{z}}_{i,1} \right) \\
\implies \underset{\tilde{\mathbf{z}}_1}{\arg \min} \left\{ \mathcal{L}(\mathbf{w}_1, \tilde{\mathbf{z}}_1) \right\}
    &= \mathbf{w}_1^\top \mathbf{x}_i
\end{aligned}
$$

Let $\tilde{\mathbf{z}}^*_1$ denote this optimal weight vector. Plugging it into Eq. \eqref{eq:loss-1} gives us its loss:

$$
\begin{equation}
\label{eq:loss-2}
\begin{aligned}
\mathcal{L}(\mathbf{w}_1, \tilde{\mathbf{z}}_1^*)
    &= \frac{1}{n} \sum_{i = 1}^n \left[ \mathbf{x}_i^\top \mathbf{x}_i - 2 (\tilde{\mathbf{z}}^*_{i,1})^2 + (\tilde{\mathbf{z}}^*_{i,1})^2 \right] \\
    &= \frac{1}{n} \sum_{i = 1}^n \left[ \mathbf{x}_i^\top \mathbf{x}_i - (\tilde{\mathbf{z}}^*_{i,1})^2 \right] \\
    &\propto - \frac{1}{n} \sum_{i = 1}^n (\mathbf{w}_1^\top \mathbf{x}_i)(\mathbf{w}_1^\top \mathbf{x}_i) \\
    &\propto -\mathbf{w}_1^\top \hat{\Sigma}\mathbf{w}_1 & \left( \frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i = \mathbf{0}_d \right)
\end{aligned}
\end{equation}
$$

where $\hat{\Sigma}$ is the empirical covariance matrix. 

Since the reconstruction error can be made arbitrarily small by giving $\mathbf{w}_1$ arbitrarily large magnitude, we restrict $\rvert \rvert \mathbf{w}_1 \rvert \rvert_2^2 = 1$. Thus, we want to solve the following auxiliary problem:

<aside><p>The restriction is equivalent to $\rvert \rvert \mathbf{w}_1 \rvert \rvert_2^2 -1 = 0$$.</p></aside>

$$
\begin{aligned}
\underset{\mathbf{w}_1: \rvert \rvert \mathbf{w}_1 \rvert \rvert_2^2 = 1}{\arg \min} \left\{ \tilde{\mathbf{L}}(\mathbf{w}_1, \tilde{\mathbf{z}}_1^*) \right\} 
    &= \underset{\mathbf{w}_1: \rvert \rvert \mathbf{w}_1 \rvert \rvert_2^2 = 1}{\arg \min} \left\{ -\mathbf{w}_1^\top \hat{\Sigma}\mathbf{w}_1 \right\} 
\end{aligned}
$$

Via <a href="/stats-ml/constrained-optimization">Lagrange multipliers</a>, we can make this an unconstrained optimization problem. We can then just take the gradient, set that equal to zero, and solve.

$$
\begin{aligned}
\underset{\mathbf{w}_1: \rvert \rvert \mathbf{w}_1 \rvert \rvert_2^2}{\arg \min} \left\{ \tilde{\mathbf{L}}(\mathbf{w}_1, \tilde{\mathbf{z}}_1^*) \right\} 
    &= \underset{\mathbf{w}_1}{\arg \min} \left\{ -
        \mathbf{w}_1^\top \hat{\Sigma} \mathbf{w}_1 + \lambda_1 (\mathbf{w}_1^\top \mathbf{w}_1 - 1) 
    \right\} \\
\frac{\partial \tilde{\mathbf{L}}(\mathbf{w}_1, \tilde{\mathbf{z}}_1^*)}{\partial \mathbf{w}_1}
    &= -2 \hat{\Sigma} \mathbf{w}_1 + 2 \lambda_1 \mathbf{w}_1 = \mathbf{0}_d \\
    \implies \hat{\Sigma} \mathbf{w}_1 &= \lambda_1 \mathbf{w}_1
\end{aligned}
$$

Notice that the last line above implies that $\mathbf{w}_1$ is an <a href="https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors">eigenvector</a> of the sample covariance matrix. Since we wish to minimize $\tilde{\mathbf{L}}(\mathbf{w}_1, \tilde{\mathbf{z}}_1^*) = -\mathbf{w}_1^\top \hat{\Sigma} \mathbf{w}_1$, we want to maximize $\mathbf{w}_1^\top \hat{\Sigma} \mathbf{w}_1$. Left multiplying the last line by $\mathbf{w}_1^\top$, we see that:

$$
\mathbf{w}_1^\top \hat{\Sigma} \mathbf{w}_1 = \lambda_1 \mathbf{w}_1^\top \mathbf{w}_1 = \lambda_1
$$

so we want to make $\lambda_1$ as large as possible. This corresponds to selecting $\mathbf{w}_1$ as the eigenvector of $\hat{\Sigma}$ with the largest eigenvalue. 

### Induction
Now, suppose we want to do the same thing but with $d = 2$. That is, we want to minimize the reconstruction error subject to the constraints that $\mathbf{w}_1^\top \mathbf{w}_2 = 0$ and $\rvert \rvert \mathbf{w}_1 \rvert \rvert_2^2 = \rvert \rvert \mathbf{w}_2 \rvert \rvert_2^2 = 1$. 

Let $$\tilde{\mathbf{z}}_2 = (\mathbf{z}_{1,2}, \mathbf{z}_{2,2}, \dots, \mathbf{z}_{n, 2})^\top$$, and let $\tilde{\mathbf{z}}_1$ be as before. Since $\mathbf{w}_1$ and $\mathbf{w}_2$ are orthogonal, we can first find the optimal $\mathbf{w}_1$ (which we did in the previous section) and then find the optimal $\mathbf{w}_2$. 

We write the loss function, take the derivative, and set that equal to zero:

$$
\begin{aligned}
\mathcal{L}(\mathbf{w}_1, \mathbf{w}_2, \tilde{\mathbf{z}}_1, \tilde{\mathbf{z}}_2)
    &= \frac{1}{n} \sum_{i = 1}^n \rvert \rvert \mathbf{x}_i - \tilde{\mathbf{z}}_{i,1} \mathbf{w}_1 - \tilde{\mathbf{z}}_{i,2} \mathbf{w}_2 \rvert \rvert_2^2 \\
\frac{\partial \mathcal{L}(\mathbf{w}_1, \mathbf{w}_2, \tilde{\mathbf{z}}_1, \tilde{\mathbf{z}}_2)}{\partial \tilde{\mathbf{z}}_{i,2}}
    &= \frac{\partial}{\partial \tilde{\mathbf{z}}_{i,2}} \left[ \frac{1}{n} \sum_{j = 1}^n (\mathbf{x}_j - \tilde{\mathbf{z}}_{j,1} \mathbf{w}_1 - \tilde{\mathbf{z}}_{j,2} \mathbf{w}_2)^\top (\mathbf{x}_j - \tilde{\mathbf{z}}_{j,1} \mathbf{w}_1 - \tilde{\mathbf{z}}_{j,2} \mathbf{w}_2) \right] \\
    &= \frac{1}{n} \left(-2 \mathbf{x}_i^\top \mathbf{w}_2 +2 \tilde{\mathbf{z}}_{i,2} \mathbf{w}_2^\top \right) \\
\implies \tilde{\mathbf{z}}^*_{i,2} &= \mathbf{w}_2^\top \mathbf{x}_i
\end{aligned}
$$

We find the loss for the weights and basis we've found so far (an asterisk denotes an optimal value):

$$
\begin{aligned}
\mathcal{L}(\mathbf{w}_1^*, \mathbf{w}_2, \tilde{\mathbf{z}}^*_1, \tilde{\mathbf{z}}^*_{2})
    &= \frac{1}{n} \sum_{j = 1}^n\rvert \rvert\mathbf{x}_j - \tilde{\mathbf{z}}_{j,1}^* \mathbf{w}_1^* - \tilde{\mathbf{z}}_{j,2}^* \mathbf{w}_2 \rvert\rvert_2^2 \\
    &= \frac{1}{n} \sum_{j = 1}^n \left[ \mathbf{x}_j^\top \mathbf{x}_j - 2 \tilde{\mathbf{z}}_{j,1}^* \mathbf{x}_j^\top \mathbf{w}_1^*  - 2 \tilde{\mathbf{z}}_{j,2}^* \mathbf{x}_j^\top \mathbf{w}_2 + (\tilde{\mathbf{z}}^*_{j,1})^2 (\mathbf{w}_1^*)^\top \mathbf{w}_1  + (\tilde{\mathbf{z}}^*_{j,1})^2 (\mathbf{w}_1^*)^\top \mathbf{w}_1^* \right] & \left((\mathbf{w}_1^*)^\top \mathbf{w}_2 = 0 \right) \\
    &= \frac{1}{n} \sum_{j = 1}^n \left[ \mathbf{x}_j^\top \mathbf{x}_j - 2 \tilde{\mathbf{z}}_{j,1}^* \mathbf{x}_j^\top \mathbf{w}_1^*  - 2 \tilde{\mathbf{z}}_{j,2}^* \mathbf{x}_j^\top \mathbf{w}_2 + (\tilde{\mathbf{z}}^*_{j,1})^2 + (\tilde{\mathbf{z}}^*_{j,1})^2\right] & \left( \rvert \rvert \mathbf{w}_k \rvert \rvert_2^2 = 1 \right) \\
    &= \frac{1}{n} \sum_{j = 1}^n \left[ \mathbf{x}_j^\top \mathbf{x}_j^\top - (\mathbf{w}_1^*)^\top \mathbf{x}_j \mathbf{x}_j^\top \mathbf{w}_1^* - \mathbf{w}_2^\top \mathbf{x}_j \mathbf{x}_j^\top \mathbf{w}_2 \right] \\
    &\propto - \mathbf{w}_2^\top \hat{\Sigma} \mathbf{w}_2 & \left(\mathbf{w}_1^*, \tilde{\mathbf{z}}_1^*, \tilde{\mathbf{z}}_2^* \text{ fixed}\right)
\end{aligned}
$$

We again introduce the Lagrangian:

$$
\begin{aligned}
\underset{\substack{\mathbf{w}_2: &\rvert \rvert \mathbf{w}_2 \rvert \rvert_2^2 = 1, \\ &(\mathbf{w}_1^*)^\top \mathbf{w}_2 = 0}}{\arg \min} \left\{ \tilde{\mathbf{L}}(\mathbf{w}_1^*, \mathbf{w}_2, \tilde{\mathbf{z}}_1^*, \tilde{\mathbf{z}}_2^*) \right\} 
    &= \underset{\mathbf{w}_2}{\arg \min} \left\{ 
        -\mathbf{w}_2^\top \hat{\Sigma} \mathbf{w}_2 - \lambda_{2,1}(\mathbf{w}_2^\top \mathbf{w}_2 - 1) -  \lambda_{2,2} (\mathbf{w}_1^\top \mathbf{w}_1 - 0)
    \right\} \\
\frac{\partial \tilde{\mathbf{L}}(\mathbf{w}_1^*, \mathbf{w}_2, \tilde{\mathbf{z}}_1^*, \tilde{\mathbf{z}}_2^*)}{\partial \mathbf{w}_2}
    &= -2 \hat{\Sigma} \mathbf{w}_2 + 2 \lambda_{2,1} \mathbf{w}_2 + \lambda_{2,2} \mathbf{w}_1^* = \mathbf{0}_d \\
\implies \hat{\Sigma} \mathbf{w}_2 &= \lambda_{2,1} \mathbf{w}_2 + \frac{1}{2} \lambda_{2,2} \mathbf{w}_1^*
\end{aligned}
$$

Left multiplying both sides of the above by $(\mathbf{w}_1^*)^\top$ yields:

$$
\begin{aligned}
(\mathbf{w}_1^*)^\top \hat{\Sigma} \mathbf{w}_2 &= \lambda_{2,1} (\mathbf{w}_1^*)^\top \mathbf{w}_2 + \frac{1}{2} \lambda_{2,2} (\mathbf{w}_1^*)^\top\mathbf{w}_1^* \\
\mathbf{w}_2^\top \hat{\Sigma} \mathbf{w}_1^* &= \lambda_{2,1} + \frac{1}{2} \lambda_{2,2}  & \left(\mathbf{w}_2^\top\mathbf{w}_1^* = 0; \hspace{1mm} \rvert \rvert \mathbf{w}_1^* \rvert \rvert_2^2 = 1\right) \\
\lambda_1 \mathbf{w}_2^\top \mathbf{w}_1^* &= \frac{1}{2} \lambda_{2,2} & \left(\hat{\Sigma} \mathbf{w}_1^* = \lambda_1 \mathbf{w}_1^* \right) \\
\lambda_{2,2} &= 0
\end{aligned}
$$

Plugging this back into our equation, we see that $\mathbf{w}_2$ is an eigenvector of $\hat{\Sigma}$:

$$
\begin{aligned}
\hat{\Sigma} \mathbf{w}_2 &= \lambda_{2,1} \mathbf{w}_2 + \frac{1}{2} \underbrace{\lambda_{2,2}}_{=0} \mathbf{w}_1^* \\
\hat{\Sigma} \mathbf{w}_2 &= \lambda_{2,1} \mathbf{w}_2
\end{aligned}
$$

By the same argument as in the previous section, we see that $\mathbf{w}_2$ should be chosen as the eigenvector corresponding to the second largest eigenvalue of $\hat{\Sigma}$. 

<aside><p>The above can be repeated to show that $\mathbf{W}$ corresponds to the matrix of eigenvectors of $\hat{\Sigma}$.</p></aside>

---

## Maximizing Variance
We often talk about principal components analysis as projecting the data onto a space where the axes are aligned in order of the directions of maximal variation in the data. 

To see why, consider the one-dimensional case, and let $\mathbf{w}_1^*$ and $\tilde{\mathbf{z}}_1^*$ be the optimal basis vector and weights. Since we centered the data, we have:

$$
\begin{aligned}
\mathbb{E}\left[ \tilde{\mathbf{z}}_{i,1}^* \right]
&= \mathbb{E}\left[ (\mathbf{w}_1^*)^\top \mathbf{x}_i \right] 
&= (\mathbf{w}_1^*)^\top \mathbb{E}\left[ \mathbf{x}_i \right]
&= 0
\end{aligned}
$$

Thus, the variance of the projected data is:

$$
\begin{aligned}
\hat{\text{Var}}(\tilde{\mathbf{z}}_{i,1}^*)
&= \mathbb{E}\left[ (\tilde{\mathbf{z}}_{i,1}^*)^2 \right] - \left(\mathbb{E}\left[ \tilde{\mathbf{z}}_{i,1}^* \right]\right)^2 
&= 
\end{aligned}
$$