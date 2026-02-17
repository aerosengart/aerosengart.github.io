---
layout: distill
title: Support Vector Machines
description: A Primer
date: 2026-02-17
tabs: true
tags: classification models
toc:
  - name: Hyperplanes
  - name: Support Vector Classifier
    subsections:
        - name: Optimization
        - name: KKT Conditions
  - name: Support Vector Machine
    subsections:
        - name: Extensions
bibliography: stats-ml.bib
---

There are a multitude of ways to perform classification. One of the foundational techniques for this is the <strong>support vector machine</strong> (<strong>SVM</strong>). SVMs are extensions of separating hyperplanes to classes that cannot be perfectly linearly separated.<d-cite key=hastie2017></d-cite>

---

## Hyperplanes
For any $\beta \in \mathbb{R}^p$, we can represent all hyperplanes in $\mathbb{R}^p$ as the set of vectors satisfying:

$$
\{\mathbf{x} \mid \mathbf{x}^\top \beta = 0 \}
$$

<aside><p>Hyperplanes are also called <strong>affine sets</strong>.</p></aside>

We can then consider the hyperplane $L$ defined by $\beta_L = (\beta_0, \beta^\top)^\top$:

$$
L = \{ \mathbf{x} \mid \mathbf{x}^\top \beta + \beta_0 = 0 \}
$$

Note that, by definition, $\mathbf{x}^\top \beta = -\beta_0$ for any $\mathbf{x} \in L$. Furthermore, for any $\mathbf{x}_1, \mathbf{x}_2 \in L$, we have that:

$$
(\mathbf{x}_1 - \mathbf{x}_2)^\top \beta = \mathbf{x}_1^\top \beta - \mathbf{x}_2^\top \beta = - \beta_0 + \beta_0 = 0
$$

This implies that $\beta$ is orthogonal to the surface of $L$, and $\beta^* = \frac{\beta}{\rvert \rvert \beta \rvert \rvert_2}$ is the vector normal to $L$. 

<aside><p>This also follows from the definition of the <a href="https://en.wikipedia.org/wiki/Normal_(geometry)#Normal_to_planes_and_polygons">normal</a>.</p></aside>

Let $\mathbf{x}$ be an arbitrary point in $L$ and consider, now, the (signed) distance between $\mathbf{x}_0 \in \mathbb{R}^{p}$ and $L$:

$$
\begin{aligned}
\frac{\beta^* \cdot (-(\mathbf{x} - \mathbf{x}_0))}{\rvert \rvert \beta^* \rvert \rvert_2} 
&= \frac{1}{\rvert \rvert \frac{\beta}{\rvert \rvert \beta \rvert \rvert_2} \rvert \rvert_2} (\beta^*)^\top (\mathbf{x}_0 - \mathbf{x}) \\
&=  (\beta^*)^\top \mathbf{x}_0 - (\beta^*)^\top \mathbf{x} \\
&= \frac{1}{\rvert \rvert \beta \rvert \rvert_2} \left[ \beta^\top (\mathbf{x}_0 - \mathbf{x}) \right]\\
&= \frac{1}{\rvert \rvert \beta \rvert \rvert_2} \left[ \beta^\top \mathbf{x}_0 - (- \beta_0) \right]\\
&= \frac{1}{\rvert \rvert \beta \rvert \rvert_2} \left[ \beta^\top \mathbf{x}_0 + \beta_0 \right]
\end{aligned}
$$

This shows that the function $f(\mathbf{x}) = \beta^\top \mathbf{x} + \beta_0 = 0$ that defines $L$ above is proportional to the signed distance between $\mathbf{x}$ and $L$.

---

## Support Vector Classifier
Now, we'll discuss the <strong>support vector classifier</strong>, which is a classifier that constructs the optimal separating hyperplane between two perfectly linearly separable classes. By <i>optimal</i>, we mean that the distance between the closest point to the hyperplane and the hyperplane is maximized. 

<aside><p>This distance is called the <strong>margin</strong>.</p></aside>

Let's say we have a training set of $N$ observations with class labels $y_i \in \\{ -1, 1 \\}$ and $p$-dimensional predictors $\mathbf{x}_i \in \mathbb{R}^{p}$. Let $\mathbf{X}$ be the $n \times p$ matrix with rows equal to the $\mathbf{x}_i$ vectors.

Suppose we have the hyperplane $L$ again, and we classify points with positive distance to $L$ to one class and those with negative distance to $L$ to another. We can represent this binary classifier as:

$$
g(\mathbf{x}) = \text{sign}\left(\mathbf{x}^\top \beta + \beta_0 \right) 
= \begin{cases}
1 & \text{ if } \mathbf{x}^\top \beta + \beta_0 \geq 0 \\
-1 & \text{ else }
\end{cases}
$$

<aside><p>For $0-1$ classes, this classifier is called a <strong>perceptron</strong> in the machine learning literature.</p></aside>

### Optimization 
To find the optimal hyperplane, we want to maximize the distance between the points and $L$. Notice that if $y_i$ is correctly classified, then:

$$
y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq 0
$$

Thus, optimization consists of solving:

$$
\begin{equation}
\label{eq:first-opt}
\begin{aligned}
&\underset{\beta, \beta_0, \rvert \rvert \beta \rvert \rvert_2 = 1}{\max} \left\{ M \right \} \\
&\text{ subject to } 
y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq M 
\text{  } \forall i = 1, \dots, N
\end{aligned}
\end{equation}
$$

Via some rescaling and manipulation, this is equivalent to:

$$
\begin{aligned}
&\underset{\beta, \beta_0}{\min} \left\{ \frac{1}{2}\rvert \rvert \beta \rvert \rvert_2^2 \right \} \\
&\text{ subject to } 
0 \geq 1 - y_i(\mathbf{x}_i^\top \beta + \beta_0)
\text{  } \forall i = 1, \dots, N \\
\end{aligned}
$$

<details>
<summary>Proof.</summary>
Moving the norm requirement into the constraint (and thus rescaling $\beta_0$), this is equivalent to:
$$
\begin{equation}
\label{eq:second-opt}
\begin{aligned}
&\underset{\beta, \beta_0}{\max} \left\{ M \right \} \\
&\text{ subject to } 
\frac{1}{\rvert \rvert \beta \rvert\rvert_2} \left[ y_i(\mathbf{x}_i^\top \beta + \beta_0)\right] \geq M 
\text{  } \forall i = 1, \dots, N \\
\iff
&\underset{\beta, \beta_0}{\max} \left\{ M \right \} \\
&\text{ subject to } 
 y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq \rvert \rvert \beta \rvert\rvert_2 M
\text{  } \forall i = 1, \dots, N \\
\end{aligned}
\end{equation}
$$
Again rescaling (since any scalar multiple of the optimal $\beta$ and $\beta_0$ will also satisfy the constraints) by setting $\rvert \rvert \beta \rvert \rvert_2 = \frac{1}{M}$, we can write:
$$
\begin{equation}
\label{eq:third-opt}
\begin{aligned}
&\underset{\beta, \beta_0}{\max} \left\{ \frac{1}{\rvert \rvert \beta \rvert \rvert_2} \right \} \\
&\text{ subject to } 
y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq 1
\text{  } \forall i = 1, \dots, N \\
\iff
&\underset{\beta, \beta_0}{\min} \left\{ \rvert \rvert \beta \rvert \rvert_2 \right \} \\
&\text{ subject to } 
y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq 1
\text{  } \forall i = 1, \dots, N \\
\iff
&\underset{\beta, \beta_0}{\min} \left\{ \frac{1}{2}\rvert \rvert \beta \rvert \rvert_2^2 \right \} \\
&\text{ subject to } 
y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq 1
\text{  } \forall i = 1, \dots, N \\
\end{aligned}
\end{equation}
$$
The margin is now encoded in the norm of $\beta$ (inversely).We'll rewrite the last line of Eq. \eqref{eq:third-opt} as:
$$
\begin{aligned}
&\underset{\beta, \beta_0}{\min} \left\{ \frac{1}{2}\rvert \rvert \beta \rvert \rvert_2^2 \right \} \\
&\text{ subject to } 
0 \geq 1 - y_i(\mathbf{x}_i^\top \beta + \beta_0)
\text{  } \forall i = 1, \dots, N \\
\end{aligned}
$$
</details>

Note that is a convex optimization problem (quadratic objective with linear inequality constraints), so we can use <a href="/stats-ml/constrained-optim/">duality</a> to rewrite this as:

$$
\begin{equation}
\label{eq:lagrange-dual}
\begin{aligned}
&\underset{\lambda}{\max} \left\{ \underset{\beta, \beta_0}{\inf} \left\{ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + \sum_{i = 1}^N \lambda_i \left(1 - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) \right\} \right\} \\
&\text{ subject to }
\lambda_i \geq 0 \text{  } \forall i = 1, \dots, N
\end{aligned}
\end{equation}
$$

<aside><p>This is called the <strong>primal form</strong>.</p></aside>

Assuming the objective and all of the constraints are convex and continuously differentiable functions, we can recast the problem in terms of the gradients. We have:

$$
\begin{aligned}
\frac{\partial }{\partial \beta} \left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + \sum_{i = 1}^N \lambda_i \left(1 - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) \right] &= \beta - \sum_{i= 1}^N \lambda_i y_i \mathbf{x}_i \\
\implies
\beta &= \sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i \\
\frac{\partial}{\partial \beta_0} \left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + \sum_{i = 1}^N \lambda_i \left(1 - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) \right] &= - \sum_{i = 1}^N \lambda_i y_i \\
\implies 
0 &= \sum_{i = 1}^N \lambda_i y_i
\end{aligned}
$$

We substitute these into the objective function in Eq. \eqref{eq:lagrange-dual} to get:

$$
\begin{equation}
\label{eq:wolfe-dual}
\begin{aligned}
&\underset{\lambda}{\max} \left\{ \sum_{i = 1}^N \lambda_i - \frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i^\top \mathbf{x}_j \right\} \\
&\text{ subject to } \sum_{i = 1}^N \lambda_i y_i = 0 \text{ and } \lambda_i \geq 0 \text{  } \forall i = 1, \dots, N
\end{aligned}
\end{equation}
$$

<aside><p>This is called the <strong>dual form</strong>.</p></aside>

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\mathcal{L}_D &= \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + \sum_{i = 1}^N \lambda_i \left(1 - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right)  \\
&= \frac{1}{2} \left\rvert\left\rvert \sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i \right\rvert\right\rvert_2^2 + \sum_{i = 1}^N \lambda_i \left(1 - y_i\left(\mathbf{x}_i^\top \sum_{j = 1}^N \lambda_j y_j \mathbf{x}_j + \beta_0 \right) \right) \\
&= \frac{1}{2}\left(\sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i \right)^\top \left(\sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i \right) + \sum_{i = 1}^N \lambda_i \left(1 - \sum_{j = 1}^N \lambda_j y_i y_j \mathbf{x}_i^\top \mathbf{x}_j - y_i \beta_0 \right) \\
&= \frac{1}{2}\sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i^\top \mathbf{x}_j - \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i^\top \mathbf{x}_j + \sum_{i = 1}^N \lambda_i y_i \beta_0 + \sum_{i = 1}^N \lambda_i \\
&= \sum_{i = 1}^N \lambda_i - \frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i^\top \mathbf{x}_j - \beta_0 \sum_{i = 1}^N \lambda_i y_i 
\end{aligned}
$$
Under the constraint that $0 = \sum_{i = 1}^N \lambda_i y_i$, the objective then becomes:
$$
\mathcal{L}_D = \sum_{i = 1}^N \lambda_i - \frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i^\top \mathbf{x}_j 
$$
</details>

This is a quadratic optimization problem with a linear constraint and a cone (the non-negative orthant) constraint, and it can be solved with standard packages (e.g. `cvxopt`).


### KKT Conditions
In order to be optimal, the solution must satisfy the <a href="https://en.wikipedia.org/wiki/Karush–Kuhn–Tucker_conditions">Karush-Kuhn-Tucker conditions</a>. That is, they must satisfy:

<ol>
<li><strong>Stationarity</strong>: $$\frac{\partial}{\partial \beta_L} \left[ \frac{1}{2} \rvert\rvert \beta \rvert \rvert_2^2 + \sum_{i = 1}^N \lambda_i \left(1 - y_i(\mathbf{x}_i^\top \beta + \beta_0)\right) \right] = \mathbf{0}_{p + 1}$$</li>
<li><strong>Primal Feasibility</strong>: $$y_i(\mathbf{x}_i^\top \beta + \beta_0) - 1 \geq 0 \text{  } \forall i = 1, \dots, N$$</li>
<li><strong>Dual Feasibility</strong>: $$\lambda_i \geq 0 \text{  } \forall i = 1, \dots, N$$</li>
<li><strong>Complementary Slackness</strong>: $$\lambda_i \left( y_i(\mathbf{x}_i^\top \beta + \beta_0) - 1\right) = 0 \text{  } \forall i = 1, \dots, N$$</li>
</ol>

The stationarity condition states that:

$$
\beta = \sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i
$$

which shows that the optimal $\beta$ must be a linear combination of the $\mathbf{x}_i$ vectors. The complementary slackness condition implies that if $\lambda_i > 0$, then $y_i(\mathbf{x}_i^\top \beta + \beta_0) = 1$. On the other hand, if $\lambda_i = 0$, then the $\mathbf{x}_i$ does not contribute to $\beta$. Thus, the optimal $\beta$ is a linear combination of some of the $\mathbf{x}_i$ vectors, which are called <strong>support vectors</strong>. 

<aside><p>Looking at the earlier iterations of our problem, we see that points that lie exactly on the boundary of the margin will be these points.</p></aside>

Thus, the optimal separating hyperplane (i.e. the support vector classifier) is the one represented by:

$$
\begin{aligned}
G(\mathbf{x}) &= \text{sign}\left(\mathbf{x}^\top \hat{\beta} + \hat{\beta}_0 \right) \\
&= \text{sign}\left( \mathbf{x}^\top  \sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i + \hat{\beta}_0 \right) \\
&= \text{sign}\left( \sum_{i = 1}^N \lambda_i y_i \mathbf{x}^\top \mathbf{x}_i + \hat{\beta}_0 \right) \\
\end{aligned}
$$

where $\hat{\beta}$ and $\hat{\beta}_0$ are the solutions to the dual problem in Eq. \eqref{eq:wolfe-dual}.

---

## Support Vector Machine
The <strong>support vector machine</strong> extends the support vector classifier to classes that are not linearly separable (in the feature space). In this case, we must adjust the optimization problem because there is no longer a separating hyperplane (let alone an optimal one). 

To get around this fact, we introduce <strong>slack variables</strong> $\xi = (\xi_1, xi_2, \dots, \xi_N)^\top$. We adjust the constraint to include these slack variables:

$$
\begin{equation}
\label{eq:svm-1}
y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq M(1 - \xi_i)
\end{equation}
$$

where $\xi_i \geq 0$ for all $i = 1, \dots, N$ and $\sum_{i = 1}^N \xi_i \leq c$ for some constant $c$. By incorporating $\xi$ in this way, we permit each observation to cross the appropriate margin by an amount proportional to $\xi_i$ (i.e. only by $M \xi_i$). The constraint that the sum of slack variables does not exceed $c$ forces the total amount of "crossing" to only be so large. Since misclassifications occur if $\xi_i > 1$, the number of misclassifications cannot exceed $c$. 

We can write the optimization problem as:

$$
\begin{equation}
\label{eq:svm-2}
\begin{aligned}
&\underbrace{\lambda, \mu}{\max} \left\{ \sum_{i = 1}^N \lambda_i -\frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i \mathbf{x}_j  \right\} \\
&\text{ subject to } 
\sum_{i = 1}^N \lambda_i y_i = 0 \text{ and } 0 \leq \lambda_i \leq C \text{  } \forall i = 1, \dots, N
\end{aligned}
\end{equation}
$$

<details>
<summary>Proof.</summary>
As we showed in the SVC case, we can rescale (since any scalar multiple of the optimal $\beta$ and $\beta_0$ will also satisfy the constraints) by setting $\rvert \rvert \beta \rvert \rvert_2 = \frac{1}{M}$ to get:
$$
\begin{aligned}
&\underset{\beta, \beta_0}{\min} \left\{ \frac{1}{2}\rvert \rvert \beta \rvert \rvert_2^2 \right \} &\\
&\text{ subject to } 
&y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq 1 - \xi_i \\
& &\sum_{i = 1}^N \xi_i \leq c \\
& &\xi_i \geq 0
\text{  } \forall i = 1, \dots, N \\
\end{aligned}
$$
and we've just added in the constraints on the slack variables. We can reparametrize the constraint on the sum of the slack variables with $C$ and put it in the objective function to get the equivalent form:
$$
\begin{aligned}
&\underset{\beta, \beta_0}{\min} \left\{ \frac{1}{2}\rvert \rvert \beta \rvert \rvert_2^2  + C \sum_{i = 1}^N \xi_i \right \} &\\
&\text{ subject to } 
&y_i(\mathbf{x}_i^\top \beta + \beta_0) \geq 1 - \xi_i \\
& &\xi_i \geq 0
\text{  } \forall i = 1, \dots, N \\
\end{aligned}
$$
Similar to the SVC case, we rewrite this as the Lagrangian dual problem:
$$
\begin{equation}
\label{eq:svm-lagrange}
\begin{aligned}
&\underset{\lambda, \mu}{\max} \left\{ \underset{\beta, \beta_0}{\inf} \left\{ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) - \sum_{i = 1}^N \mu_i \xi_i \right\} \right\}\\
&\text{ subject to }
\lambda_i, \mu_i \geq 0 \text{  } \forall i = 1, \dots, N
\end{aligned}
\end{equation}
$$
Taking the derivative with respect to $\beta$, $\beta_0$, and $\xi$, and setting equal to zero, we get:
$$
\begin{aligned}
\frac{\partial}{\partial \beta} \left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right)  - \sum_{i = 1}^N \mu_i \xi_i \right] &= \beta - \sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i \\
\implies
\beta &= \sum_{i = 1}^N \lambda_i y_i \mathbf{x} \\
\frac{\partial}{\partial \beta_0}\left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) - \sum_{i = 1}^N \mu_i \xi_i\right] 
&= - \sum_{i = 1}^N \lambda_i y_i \\
\implies
0 &= \sum_{i = 1}^N \lambda_i y_i \\
\frac{\partial}{\partial \xi_i} \left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) - \sum_{i = 1}^N \mu_i \xi_i\right]
&= C - \lambda_i + \mu_i \\
\implies
\lambda_i &= C - \mu_i 
\end{aligned}
$$
Plugging these into the objective funciton of Eq. \eqref{eq:svm-lagrange}, we get:
$$
\begin{aligned}
\mathcal{L}_D 
&= \frac{1}{2} \rvert \rvert \beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) - \sum_{i = 1}^N \mu_i \xi_i \\
&= \frac{1}{2} \left(\sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i \right)^2 + C \sum_{i = 1}^N \xi_i - \sum_{i = 1}^N \lambda_i \left(\xi_i + y_i \left( \mathbf{x}_i^\top \sum_{j = 1}^N \lambda_j y_j \mathbf{x}_j + \beta_0 \right) - 1\right) - \sum_{i = 1}^N \mu_i \xi_i \\
&= \frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i \mathbf{x}_j + C \sum_{i = 1}^N \xi_i - \sum_{i = 1}^N (C - \mu_i) \xi_i + \sum_{i = 1}^N \lambda_i - \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j  \mathbf{x}_i^\top \mathbf{x}_j - \sum_{i = 1}^N \mu_i \xi_i \\
&= \sum_{i = 1}^N \lambda_i -\frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i \mathbf{x}_j 
\end{aligned}
$$
The optimization problem is then:
$$
\begin{aligned}
&\underbrace{\lambda, \mu}{\max} \left\{ \sum_{i = 1}^N \lambda_i -\frac{1}{2} \sum_{i = 1}^N \sum_{j = 1}^N \lambda_i \lambda_j y_i y_j \mathbf{x}_i \mathbf{x}_j  \right\} \\
&\text{ subject to } 
\sum_{i = 1}^N \lambda_i y_i = 0 \text{ and } 0 \leq \lambda_i \leq C \text{  } \forall i = 1, \dots, N
\end{aligned}
$$
</details>

The KKT conditions are:

<ol>
<li><strong>Stationarity</strong>: 
$$
\begin{aligned}
&\frac{\partial}{\partial \beta_L} \left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right) - \sum_{i = 1}^N \mu_i \xi_i \right] = \mathbf{0}_{p + 1} \\
&\frac{\partial}{\partial \xi} \left[ \frac{1}{2} \rvert \rvert\beta \rvert \rvert_2^2 + C \sum_{i = 1}^N \xi_i + \sum_{i = 1}^N \lambda_i \left(1 - \xi_i - y_i (\mathbf{x}_i^\top \beta + \beta_0) \right)  - \sum_{i = 1}^N \mu_i \xi_i  \right] = C\mathbf{1}_N - \mu
\end{aligned}
$$
</li>
<li><strong>Primal Feasibility</strong>: $$y_i(\mathbf{x}_i^\top \beta + \beta_0) - ( 1- \xi_i ) \geq 0 \text{  } \forall i = 1, \dots, N$$</li>
<li><strong>Dual Feasibility</strong>: $$\mu_i \xi_i \geq 0 \text{  } \forall i = 1, \dots, N$$</li>
<li><strong>Complementary Slackness</strong>: $$\lambda_i \left( y_i(\mathbf{x}_i^\top \beta + \beta_0) - (1 - \xi_i)\right) = 0 \text{  } \forall i = 1, \dots, N$$</li>
</ol>

By a similar argument to the SVC case, the optimal $\beta$ is given by:

$$
\hat{\beta} = \sum_{i = 1}^N \lambda_i y_i \mathbf{x}_i
$$

By the complementary slackness requirement, $\lambda_i > 0$ only for those points $\mathbf{x}_i$ that satisfy:

$$
y_i(\mathbf{x}_i^\top \beta + \beta_0) = 1 - \xi_i
$$

These are called <strong>support vectors</strong>, and those with $\xi_i = 0$ will lie exactly on the margin. If $\lambda_i = C$, then if $\xi_i \leq 1$, the point will lie within the margin but be classified correctly. Otherwise, it is incorrectly classified.

The support vector machine is the classifier:

$$
\begin{aligned}
G(\mathbf{x}) &= \text{sign} \left( \mathbf{x}^\top \hat{\beta} + \hat{\beta}_0 \right) \\
&= \text{sign} \left( \mathbf{x}^\top \sum_{i = 1}^N  \lambda_i y_i\mathbf{x}_i + \hat{\beta}_0 \right) \\
&= \text{sign} \left( \sum_{i = 1}^N  \lambda_i y_i \mathbf{x}^\top \mathbf{x}_i + \hat{\beta}_0 \right) \\
\end{aligned}
$$

where $C$ is left as a tuning parameter, and $\hat{\beta}$ and $\hat{\beta}_0$ are the solutions to Eq. \eqref{eq:svm-2}.


### Extensions
The SVM can be extended beyond the inner product with respect to the Euclidean norm via the <a href="/stats-ml/kernel.md">kernel trick</a>. The classifier then becomes:

$$
G(\mathbf{x}) = \text{sign} \left( \sum_{i = 1}^N \lambda_i y_i \mathcal{K}(h(\mathbf{x}), h(\mathbf{x}_i)) + \beta_0 \right)
$$

for some choice of symmetric positive semi-definite kernel function $\mathcal{K}(\cdot, \cdot)$.