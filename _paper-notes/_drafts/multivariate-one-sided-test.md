---
layout: distill
title: A Multivariate Analogue of the One-Sided Test
description: 
date: 2026-02-18
tabs: true
tags: likelihood theory testing paper-review
toc:
  - name: Background
    subsections:
        - name: An Illustrative Example
        - name: Another Example
  - name: Maximum Likelihod Estimation
bibliography: 2026-02-18-multivariate-one-sided-testing.bib
---

The $\bar{\chi}^2$ distribution comes up in the theory behind hypothesis tests against one-sided alternatives. It seems like it would be helpful to have some more intuition, particularly from a more geometric perspective, about how this distribution comes about. In this post, I work through most of Kudô's 1963 paper <i>A Multivariate Analogue of the One-Sided Test</i>.<d-cite key=kudo1963></d-cite>

---

## Background
We will consider the setting of a $k$-dimensional multivariate Gaussian population with mean vector $\boldsymbol{\mu} = (\mu_1, \mu_2, \dots, \mu_k)^\top$ and <i>known</i> (and non-singular) covariance matrix, $\boldsymbol{\Sigma}$. We have $n$ independent samples, denoted by $\mathbf{x}_1, \dots, \mathbf{x}_n$ from this population

For now, we will concern ourselves with Euclidean space. We will let $\mathcal{C}_a$, $\mathcal{C}_b$, $\mathcal{C}$, and $\mathcal{M}$ to be subsets of $\mathbb{R}^p$ where $\mathcal{C}$ is a closed, convex cone, $\mathcal{M}$ is a linear space, $\mathcal{M} \subset \mathcal{C}$, and $\mathcal{C}_a \subset \mathcal{C}_b$. We will also define a <strong>polyhedral</strong> in $\mathbb{R}^p$ as a set:

$$
\mathcal{P} = \left\{ \boldsymbol{\theta} \in \mathbb{R}^p : \mathbf{a}_1^\top \boldsymbol{\theta} \geq 0, \dots, \mathbf{a}_k^\top \boldsymbol{\theta} \geq 0 \right\}
$$

for some $\mathbf{a}_1, \dots, \mathbf{a}_k \in \mathbb{R}^p$. 

We are concerned with what Silvapulle and Sen<d-cite key=silvapulle2005></d-cite> call <strong>Type A Tests</strong>, which are hypothesis testing problems of the form:

$$
\begin{equation}
\label{eq:hypotheses}
H_0: \boldsymbol{\mu} \in \mathcal{M}
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\mu} \in \mathcal{C}, \boldsymbol{\mu} \notin \mathcal{M}
\end{equation}
$$

for some fixed matrix $\mathbf{R}$.


### An Illustrative Example
We first consider the case when $\mathbf{x} = (x_1, x_2)^\top \sim \mathcal{N}\left(\boldsymbol{\mu}, \mathbb{I}_{2 \times 2}\right)$ with $\boldsymbol{\mu} = (\mu_1, \mu_2)^\top$. Suppose $\mathcal{C}$ is the non-negative orthant, and consider the Type A testing problem:

<aside><p>This is Example 3.3.1 in Silvapulle and Sen.<d-cite key=silvapulle2005></d-cite></p></aside>

$$
H_0: \mu_1 = \mu_2 = 0 
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \mu_1 \geq 0, \mu_2 \geq 0; (\mu_1, \mu_2) \neq (0, 0)
$$

Let $\hat{\boldsymbol{\mu}}_0$ and $\hat{\boldsymbol{\mu}}_1$ denote the MLE of $\boldsymbol{\mu}$ under $H_0$ and $H_1$, respectively. We have that the log-likelihood satisfies:

$$
-2 \ell(\boldsymbol{\mu}) = -2 \left[ 
  - \frac{2}{2} \log(2 \pi) - \frac{1}{2} (\mathbf{x} - \boldsymbol{\mu})^\top (\mathbf{x} - \boldsymbol{\mu})
\right] \propto \rvert \rvert \mathbf{x} - \boldsymbol{\mu} \rvert \rvert_2^2
$$

Since the MLE of $\boldsymbol{\mu}$ is the minimizer of the negative log-likelihood, we have that $\hat{\boldsymbol{\mu}}_1$ is the projection of $\mathbf{x}$ onto $\mathcal{C}$ (the closest point in $\mathcal{C}$ to $\mathbf{x}$). We notate the projection as:

$$
\hat{\boldsymbol{\mu}}_1 = \Pi(\mathbf{x} \rvert \mathcal{C})
$$

Thus, if $\mathbf{x}$ falls in the first quadrant ($Q_1$), then $\hat{\boldsymbol{\mu}}_1 = \mathbf{x}$, since the coordinates are all non-negative. If $\mathbf{x}$ falls in the second quadrant ($Q_2$), then $\hat{\boldsymbol{\mu}}_1 = (0, x_2)$. Similarly, if $\mathbf{x}$ falls in $Q_3$ or $Q_4$, $\hat{\boldsymbol{\mu}}_1 = (x_1, 0)$ and $(0, 0)$, respectively.

The likelihood ratio test statistic is given by:

$$
\begin{equation}
\label{eq:lrt-ex}
\begin{aligned}
\text{LRT} 
  &= 2( \ell(\hat{\boldsymbol{\mu}}_1) - \ell(\hat{\boldsymbol{\mu}}_0) ) 
  = -2 (\ell(\hat{\boldsymbol{\mu}}_0)- \ell(\hat{\boldsymbol{\mu}}_1) )
  \propto \rvert \rvert \mathbf{x} \rvert \rvert_2^2 - \rvert \rvert \mathbf{x} - \hat{\boldsymbol{\mu}}_1 \rvert \rvert_2^2
  = \rvert \rvert \hat{\boldsymbol{\mu}}_1 \rvert \rvert_2^2
\end{aligned}
\end{equation}
$$

where the last equality follows from the fact that $\hat{\boldsymbol{\mu}}_1$ is the orthogonal projection of $\mathbf{x}$ onto $\mathcal{C}$, so $\mathbf{x} - \hat{\boldsymbol{\mu}}_1$ is orthogonal to $\hat{\boldsymbol{\mu}}_1$. 

Under $H_0$, we can decompose $\mathbb{P}(\text{LRT} \leq c)$ for $c > 0$ as:

$$
\begin{aligned}
\mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x})
  &= \sum_{i = 1}^4 \mathbb{P}\left(\text{LRT}  \leq c \cap \mathbf{x} \in Q_i\right) \\
  &= \sum_{i = 1}^4 \mathbb{P}\left(\text{LRT}  \leq c \rvert \mathbf{X} \in Q_i\right) \mathbb{P}(\mathbf{x} \in Q_i ) \\
  &= \frac{1}{4} \mathbb{P}(\chi_2^2 \leq c) + \frac{1}{2} \mathbb{P}(\chi_1^2 \leq c) + \frac{1}{4} \\
  &= \sum_{i = 0}^2 \omega_i \mathbb{P}\left(\chi_i^2 \leq c\right) & \left(\boldsymbol{\omega} = (0.25, 0.5, 0.25) \right)
\end{aligned}
$$

<details>
<summary>Proof.</summary>
Under $H_0$, we can show that $\frac{\mathbf{x}}{\rvert \rvert \mathbf{x} \rvert \rvert_2}$ is independent of $\rvert \rvert \mathbf{x} \rvert \rvert_2$, which implies that the distribution of $x_1^2 + x_2^2$ conditional on $\frac{\mathbf{x}}{\rvert \rvert \mathbf{x} \rvert \rvert_2}$ falling in $Q_1$ (the first quadrant), is the same as its unconditional distribution. Also note that $\hat{\boldsymbol{\mu}}_1 = \mathbf{x}$ if $\mathbf{x} \in Q_1$ since it satisfies the constraints. Thus, for some value $c$, we have:
<aside><p>HOW?</p></aside>
$$
\begin{aligned}
\mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x})
  &= \sum_{i = 1}^4 \mathbb{P}\left(\text{LRT}  \leq c \cap \mathbf{x} \in Q_i\right) \\
  &= \sum_{i = 1}^4 \mathbb{P}\left(\text{LRT}  \leq c \rvert \mathbf{X} \in Q_i\right) \mathbb{P}(\mathbf{x} \in Q_i )
\end{aligned}
$$
We then see that:
$$
\begin{aligned}
\mathbb{P}(\text{LRT}  \leq c \rvert \mathbf{x} \in Q_1)
  &= \mathbb{P}\left( \rvert \rvert \hat{\boldsymbol{\mu}}_1 \rvert \rvert^2_2 \leq c \rvert \mathbf{x} \in Q_1 \right) \\
  &= \mathbb{P}\left( \rvert \rvert \mathbf{x} \rvert \rvert_2^2 \leq c \right) & \left(\hat{\boldsymbol{\mu}}_1 = \mathbf{x} \text{ if } \mathbf{x} \in Q_2 \right)\\
  &= \mathbb{P}\left(\chi^2_2 \leq c \right) & \left(x_1 \parallel x_2 \right)
\end{aligned}
$$
We also have:
$$
\begin{aligned}
\mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x} \in Q_2) 
  &= \mathbb{P}\left(\rvert \rvert \hat{\boldsymbol{\mu}}_1 \rvert \rvert^2_2 \leq c \rvert  x_2 \geq 0, x_1 \leq 0  \right) \\
  &= \mathbb{P}\left( x_2^2 \leq c \rvert x_2 \geq 0, x_1 \leq 0 \right) & \left( \hat{\boldsymbol{\mu}}_1 = (0, x_2) \text{ if } \mathbf{x} \in Q_2 \right) \\
  &= \mathbb{P}\left( x_2^2 \leq c \right) & \left(x_2 \text{ sym. about } 0\right) \\
  &= \mathbb{P}(\chi_1^2 \leq c) & \left(x_2 \sim \mathcal{N}(0, 1)\right)
\end{aligned}
$$
By the same argument, $\mathbb{P}(LRT \leq c \rvert \mathbf{x} \in Q_4) = \mathbb{P}(\chi_1^2 \leq c)$. Finally, we see that:
$$
\begin{aligned}
\mathbb{P}\left(LRT \leq c \rvert \mathbf{x} \in Q_3 \right) 
  &= \mathbb{P}\left(\rvert \rvert \hat{\boldsymbol{\mu}}_1 \rvert \rvert_2^2 \leq c \rvert x_1, x_2 \leq 0 \right) \\
  &= \mathbb{P}\left(0 \leq c \right) \\
  &= 1
\end{aligned}
$$
We then plug these results into the expression to obtain the desired result.
</details>

Something of note is that the weights $\boldsymbol{\omega}$ correspond to the probability of $\mathbf{x}$ falling in $Q_3$, $Q_2 \cup Q_4$, and $Q_1$, respectively. 

### Another Example
Suppose we expand our setting to be $\mathbf{x} = (x_1, x_2)^\top \sim \mathcal{N}(\boldsymbol{\mu}, \mathbf{V})$ for some positive definite matrix, $\mathbf{V}$. We now consider the case of testing:

$$
H_0: \boldsymbol{\mu} = \mathbf{0}_2
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\mu} \in \mathcal{C} 
$$

<aside><p>This is Example 3.3.4 in Silvapulle and Sen.<d-cite key=silvapulle2005></d-cite></p></aside>

where we have the closed, convex cone $$\mathcal{C} = \left\{ \mathbf{c} : \mathbf{R} \mathbf{c} \geq \mathbf{0}_2 \right\}$$ for some non-singular $\mathbf{R} \in \mathbb{R}^{2 \times 2}$.

As a positive definite matrix, $\mathbf{V}$ is invertible. We can also factorize the inverse into $\mathbf{V}^{-1} = \mathbf{A}^\top \mathbf{A}$ where $\mathbf{A}$ is invertible (see <a href="https://en.wikipedia.org/wiki/Definite_matrix#Decomposition">this Wikipedia page</a>). We can then define:

$$
\begin{equation}
\label{eq:example-2}
\mathbf{x}^\star = \mathbf{A} \mathbf{x};
\hspace{4mm}
\boldsymbol{\mu}^\star = \mathbf{A} \boldsymbol{\mu};
\hspace{4mm}
\mathcal{C}^\star = \left\{ \mathbf{A} \mathbf{c}: \mathbf{c} \in \mathcal{C} \right\} 
\end{equation}
$$

Notice that:

$$
\mathcal{C}^\star = \left\{ \mathbf{v}: \mathbf{R} \mathbf{A}^{-1} \mathbf{v} \geq \mathbf{0}_2 \right\};
\hspace{5mm}
\mathbf{x}^\star \sim \mathcal{N}\left(\boldsymbol{\mu}^\star, \mathbb{I}_{2 \times 2}\right)
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\mathbb{E}\left[ \mathbf{x}^\star  \right] 
  &= \mathbb{E}\left[ \mathbf{A} \mathbf{x} \right] \\
  &= \mathbf{A} \boldsymbol{\mu} \\
  &= \boldsymbol{\mu}^\star \\
\text{Var}\left(\mathbf{x}^\star\right) 
  &= \text{Var}(\mathbf{A} \mathbf{x}) \\
  &= \mathbf{A} \text{Var}(\mathbf{x}) \mathbf{A}^\top \\
  &= \mathbf{A} \mathbf{V} \mathbf{A}^\top \\
  &= (\mathbf{A}^\top)^{-1} \mathbf{A}^\top \mathbf{A} \mathbf{V} \mathbf{A}^\top \\
  &= (\mathbf{A}^\top)^{-1} \mathbf{V}^{-1} \mathbf{V} \mathbf{A}^\top \\
  &= (\mathbf{A}^\top)^{-1} \mathbf{A}^\top \\
  &= \mathbb{I}_{2 \times 2}
\end{aligned}
$$
</details>

Thus, the problem can be restated equivalently as:

$$
H_0: \boldsymbol{\mu}^\star = \mathbf{0}_2
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\mu}^\star \in \mathcal{C}^\star
$$

where we conduct our tests with the $\mathbf{x}^\star$, which has uncorrelated coordinates with unit variance, instead of $\mathbf{x}$. This problem can also be restated as:

$$
H_0: \boldsymbol{\mu}^\star = \mathbf{0}_2
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \mathbf{R} \boldsymbol{\mu}^\star \geq \mathbf{0}_2
$$

Recall the definition of the <i>polar cone</i> to $\mathcal{C}$ (see <a href="#paper-notes/score-test-one-sided#dual-cone">here</a>). If we define the inner product to be $\langle \mathbf{a}, \mathbf{b} \rangle = \mathbf{a}^\top \mathbf{b}$, then the polar cone, which we denote by $\mathcal{C}^\star_0$, is set of vectors which are at a non-acute angle to the vectors in $\mathcal{C}^\star$. That is:

$$
\mathcal{C}_0^\star = \left\{ \mathbf{d} : \mathbf{d}^\top \mathbf{c} \leq 0 \text{ for } \mathbf{c} \in \mathcal{C}^\star \right\}
$$

This is the space of vectors that will form non-acute angles with any $\mathbf{c} \in \mathcal{C}^\star$. Let $\mathbf{u}$ and $\mathbf{v}$ be the unit vectors that are parallel to the upper and lower boundaries of $\mathcal{C}^\star$ (i.e. the two edges of the cone). Since $\mathbf{R} \mathbf{c} \geq \mathbf{0}_2$, we have that:

$$
\begin{aligned}
\mathbf{R}^\top_{1, \cdot} \mathbf{u} &= 0 \implies
-\mathbf{R}^\top_{1, \cdot} \mathbf{u} = 0 \\
\mathbf{R}^\top_{2, \cdot} \mathbf{v} &= 0 \implies
-\mathbf{R}^\top_{2, \cdot} \mathbf{v} = 0
\end{aligned}
$$

In other words, $-\mathbf{R}^\top_{1, \cdot}$ and $-\mathbf{R}^\top_{2, \cdot}$ form the boundaries of $\mathcal{C}^0$. 

We can partition then plane into four pieces $\{S_1, S_2, S_3, S_4\}$ that mimic the quadrants $\{Q_1, Q_2, Q_3, Q_4\}$. We let $\mathcal{C}^\star = S_1$ and $$\mathcal{C}_0^\star = S_3$$, and then we let $S_2$ and $S_4$ be the cones that fall between $\mathbf{u}$ and $-\mathbf{R}^\top_{1, \cdot}$ and between $\mathbf{v}$ and $-\mathbf{R}^\top_{2, \cdot}$. 

As we saw in the earlier example, the log-likelihood is proportional to the squared Euclidean distance between the data, $\mathbf{x}^\star$, and the mean vector, $\boldsymbol{\mu}^\star$. Thus, the MLE of $\boldsymbol{\mu}^\star$ subject to $\boldsymbol{\mu}^\star \in \mathcal{C}^\star$ will be the point in $\mathcal{C}^\star$ that is closest to $\mathbf{x}^\star$. 

If $$\mathbf{x}^\star \in S_1$$, then $$\hat{\boldsymbol{\mu}}^\star_1 = \mathbf{x}^\star$$ since $$\mathbf{x}^\star \in \mathcal{C}^\star$$. If $$\mathbf{x}^\star \in S_2$$ or $S_4$, then $$\hat{\boldsymbol{\mu}}^\star_1$$ is the projection of $$\mathbf{x}^\star$$ onto $\mathbf{u}$ or $\mathbf{v}$, respectively. Lastly, if $$\mathbf{x}^\star \in S_3$$, then the closest point is the origin. Since:

$$
\mathbf{x}^\star - \hat{\boldsymbol{\mu}}^\star_1 \perp \hat{\boldsymbol{\mu}}^\star_1
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
(\mathbf{x}^\star - \hat{\boldsymbol{\mu}}^\star_1)^\top \hat{\boldsymbol{\mu}}^\star_1
  &= \begin{cases}
  (\mathbf{x}^\star - \mathbf{x}^\star)^\top \mathbf{x}^\star & \mathbf{x}^\star \in S_1 \\
  (\mathbf{x}^\star - (\mathbf{u}^\top \mathbf{x}^\star) \mathbf{u})^\top (\mathbf{u}^\top \mathbf{x}^\star)\mathbf{u} & \mathbf{x}^\star \in S_2 \\
  0 & \mathbf{x}^\star \in S_3 \\
  (\mathbf{x}^\star - (\mathbf{v}^\top \mathbf{x}^\star) \mathbf{v})^\top (\mathbf{v}^\top \mathbf{x}^\star)\mathbf{v} & \mathbf{x}^\star \in S_2 \\
  \end{cases} \\
  &= \begin{cases}
  0 & \mathbf{x}^\star \in S_1 \\
  (\mathbf{u}^\top \mathbf{x}^\star)^2 - (\mathbf{u}^\top \mathbf{x}^\star)^2 \mathbf{u}^\top \mathbf{u} & \mathbf{x}^\star \in S_2 \\
  0 & \mathbf{x}^\star \in S_3 \\
  (\mathbf{v}^\top \mathbf{x}^\star)^2 - (\mathbf{v}^\top \mathbf{x}^\star)^2 \mathbf{v}^\top \mathbf{v} & \mathbf{x}^\star \in S_2 \\& \mathbf{x}^\star \in S_4 \\
  \end{cases} \\
  &= \begin{cases}
  0 & \mathbf{x}^\star \in S_1 \\
  0 & \mathbf{x}^\star \in S_2 \\
  0 & \mathbf{x}^\star \in S_3 \\
  0 & \mathbf{x}^\star \in S_2 \\
  \end{cases} 
\end{aligned}
$$
</details>

the LRT statistic is the same as in Eq. \eqref{eq:lrt-ex}. We can then repeat the same decomposition of $\mathbb{P}(\text{LRT} \leq c)$ for $c > 0$ as before to get:

$$
\begin{aligned}
\mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x} )
  &= \sum_{i = 1}^4 \mathbb{P}(\text{LRT} \leq c \rvert \mathbf{x} \in S_i) \mathbb{P}(\mathbf{x} \in S_i) \\
  &= q \mathbb{P}(\chi_0^2 \leq c) + \frac{1}{2} \mathbb{P}(\chi_1^2 \leq c) + (\frac{1}{2} - q) \mathbb{P}(\chi_2^2 \leq c)
\end{aligned}
$$

<details>
<summary>Proof.</summary>
The case for $\mathbf{x}^\star \in S_1$ is basically the same as in the previous example since the length and angle of $\mathbf{x}^\star$ are independent under $H_0$. 
<br>
For the second case, we consider a rotation so that the axes are defined by $-\mathbf{R}_{1, \cdot}^\top$ and $\mathbf{u}$ rather than the standard basis vectors. With this transformation, if $\mathbf{x}^\star \in S_2$, then we can write the LRT statistic as the Euclidean distance between $\mathbf{x}^\star$ and our new $x$-axis. This is just the square of the first coordinate of $\mathbf{x}^\star$, which is $\chi_1^2$ by the proof in the previous example. A similar argument yields the same result for if $\mathbf{x}^\star \in S_4$. 
<br>
The case of $\mathbf{x}^\star \in S_3$ is identical to before. 
</details>

In the above, $q = \mathbb{P}(\mathbf{x}^\star \in S_3 \rvert H_0) = \mathbb{P}(\mathbf{x}^\star \in \mathcal{C}_0^\star \rvert H_0)$, the probability of the data falling in the dual cone under $H_0$. Since $\mathbf{x}^\star$ is standard normal under $H_0$, it is uniformly distributed (radially) about the origin. Thus, the probability that $\mathbf{x}^\star$ falls in any given $S_i$ is just the proportion of the circumference of the unit circle that is "chopped out" by the cone, $S_i$. If we let $\gamma$ denote the angle of the vertex of $S_3$, then $$q = \frac{\gamma}{2 \pi}$$. Since $$-\mathbf{R}_{1, \cdot}^\top$$ and $$-\mathbf{R}_{2, \cdot}^\top$$ are vectors that define the boundaries of $S_3$, we can find $\gamma$ as:

$$
\gamma = \cos^{-1}\left(\frac{\mathbf{R}_{1, \cdot}^\top \mathbf{R}_{2 \cdot}}{\sqrt{(\mathbf{R}_{1, \cdot}^\top \mathbf{R}_{2, \cdot}) (\mathbf{R}_{2, \cdot}^\top \mathbf{R}_{2, \cdot})}} \right)
$$

The probabilities remain the same, but we can back-transform $\hat{\boldsymbol{\mu}}_1$ using the relations given in Eq. \eqref{eq:example-2} if desired (invariance property of MLE). 



---

The joint distribution of our sample under $H_1$ is:

$$
\begin{aligned}
f(\mathbf{x}_1, \dots, \mathbf{x}_n)
&= \frac{1}{(\sqrt{2 \pi})^{kn} \rvert \boldsymbol{\Sigma} \rvert^n} \exp\left(- \frac{1}{2} \sum_{i = 1}^n (\mathbf{x}_i - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{x}_i - \boldsymbol{\mu}) \right) \\
&= \frac{1}{(\sqrt{2 \pi})^{kn} \rvert \boldsymbol{\Sigma} \rvert^n} \exp\left(- \frac{1}{2} \left[ \sum_{i = 1}^n (\mathbf{x}_i - \bar{\mathbf{x}})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{x}_i - \bar{\mathbf{x}}) + n (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right]\right) \\
\end{aligned}
$$

<aside><p>The second line follows by adding and subtracting $\bar{\mathbf{x}}$, the sample mean.</p></aside>

The likelihood ratio criterion is given by the ratio of the likelihood function evaluated under $H_0$ and $H_1$: 

$$
\begin{aligned}
\lambda_{LR} 
&= \frac{\underset{\substack{\mu_i = 0 \\ i = 1, \dots, k}}{\max} \left\{ f(\mathbf{x}_1, \dots, \mathbf{x}_n) \right\}}{\underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\max} \left\{ f(\mathbf{x}_1, \dots, \mathbf{x}_n) \right\}} \\
&= \frac{\exp\left(- \frac{1}{2}\sum_{i = 1}^n (\mathbf{x}_i - \bar{\mathbf{x}})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{x}_i - \bar{\mathbf{x}})\right)\exp\left( - \frac{n}{2} \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1}\bar{\mathbf{x}}\right)}{\underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\max} \left\{ \exp\left(- \frac{1}{2}  \sum_{i = 1}^n (\mathbf{x}_i - \bar{\mathbf{x}})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{x}_i - \bar{\mathbf{x}})\right)\exp\left(-\frac{n}{2}(\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu})\right)\right\}} \\
&= \frac{\exp\left( - \frac{n}{2} \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1}\bar{\mathbf{x}}\right)}{\underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\max} \left\{\exp\left(-\frac{n}{2}(\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu})\right)\right\}}
\end{aligned}
$$

Notice that maximizing the denominator is equivalent to maximizing the piece inside $\exp(\cdot)$ because the function is monotonic increasing. The likelihood ratio then can be rewritten as:

$$
\lambda_{LR} = \exp\left(- \frac{n}{2} \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1}\bar{\mathbf{x}} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\max} \left\{ -\frac{n}{2}(\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\} \right)
$$

By the same argument, analyzing this test statistic is equivalent to examining just the part within the exponential function. Since maximizing a quantity is equivalent to minimizing its negative, we can instead analyze the <strong>chi-bar squared</strong> statistic:

$$
\begin{equation}
\label{eq:chi-bar-squared}
\bar{\chi}^2 = n\left[ \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1}\bar{\mathbf{x}} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\}\right]
\end{equation}
$$

Note that the minimization portion can be completed with standard convex optimization packages (e.g. `CVXOPT` quadratic cone program) because it is a quadratic form with inequality constraints. 


---

## Distribution
Define the <i>maximum likelihood estimate</i> of $\boldsymbol{\mu}$ given $\mathbf{x}$ is given by:

$$
\hat{\boldsymbol{\mu}}^* = \underset{\begin{split}
\mathbf{m}_i \geq 0 \\
i = 1, \dots, k
\end{split}
}{\arg \min} \left\{ 
(\mathbf{x} - \mathbf{m})^\top \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \mathbf{m})
\right\}
$$

<aside><p>We've omitted the parts that do not involve $\boldsymbol{m}$, taken the natural logarithm, and multiplied by $-1$ to find the minimum.</p></aside>

For now, assume $$\hat{\boldsymbol{\mu}}^*$$ exists.

### Rewriting The Problem
For a symmetric and positive (semi-)definite matrix $\mathbf{D}$, there exists a symmetric and positive (semi-)definite matrix, which we denote by $\mathbf{B}$,<d-cite key=squareroot2026></d-cite> satisfying:

<aside><p>This is called the <strong>principal square root</strong>.</p></aside>

$$
\mathbf{D} = \mathbf{B} \mathbf{B} = \mathbf{B} \mathbf{B}^\top
$$

<aside><p>Recall $\mathbf{B}$ is symmetric.</p></aside>

This decomposition is unique up to orthogonal transformations. That is, for any orthogonal matrix $\mathbf{T}$, $\mathbf{TB}$ will also satisfy the above equation.<d-cite key=definitematrix2026></d-cite> Furthermore, if $\mathbf{D}$ is positive definite, then it is invertible, and so is $\mathbf{B}$. 

Let $\mathbf{A}$ denote such a matrix for $\boldsymbol{\Sigma}$ so that it satisfies:

$$
\begin{equation}
\label{eq:sigma-identity}
\frac{1}{n} \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top = \mathbb{I}_{k \times k}
\end{equation}
$$

Now, suppose we construct $\mathbf{y} = \mathbf{A} \bar{\mathbf{x}}$, which is clearly still multivariate Gaussian as a linear transformation of multivariate Gaussian random variables. It also has the following properties under the null hypothesis.

<div class="theorem">
<strong>Properties. </strong>
{% tabs kudo-mean-covar-y %}
{% tab kudo-mean-covar-y statement %}
For $\mathbf{y} = \mathbf{A} \bar{\mathbf{x}}$ with $\mathbf{A}$ defined as in the above section, we have, under $H_0$:

$$
\begin{aligned}
\mathbb{E}\left[ \mathbf{y} \right] &= \mathbf{0}_k \\
\text{Cov}\left(\mathbf{y}, \mathbf{y} \right) &= \mathbb{I}_{k \times k}
\end{aligned}
$$
{% endtab %}
{% tab kudo-mean-covar-y proof %}
We prove the zero expectation statement first.

$$
\begin{aligned}
\mathbb{E}_{H_0}\left[ \mathbf{y} \right] 
  &= \mathbb{E}_{H_0}\left[ \mathbf{A} \bar{\mathbf{x}} \right] \\
  &= \mathbf{A} \mathbb{E}_{H_0} \left[ \frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i \right] \\
  &= \mathbf{A} \mathbb{E}_{H_0} \left[ \boldsymbol{\mu} \right] \\
  &= \mathbf{A} \mathbf{0}_k \\
  &= \mathbf{0}_k 
\end{aligned}
$$

To show the second statement: 

$$
\begin{aligned}
\text{Cov}(\mathbf{y}, \mathbf{y})
&= \mathbb{E}_{H_0}\left[\left(\mathbf{y} - \mathbb{E}\left[ \mathbf{y} \right]\right) \left(\mathbf{y} - \mathbb{E}\left[ \mathbf{y}\right]\right)^\top \right] \\
&= \mathbb{E}_{H_0}\left[\left(\mathbf{y} - \mathbf{A} \left( \frac{1}{n} \sum_{i = 1}^n \mathbb{E}\left[ \mathbf{x}_i \right] \right) \right) \left(\mathbf{y} -  \mathbf{A} \left( \frac{1}{n} \sum_{i = 1}^n \mathbb{E}\left[ \mathbf{x}_i \right] \right) \right)^\top \right] \\
&= \mathbb{E}_{H_0}\left[\left(\mathbf{y} -  \mathbf{A} \boldsymbol{\mu} \right) \left(\mathbf{y} -  \mathbf{A}\boldsymbol{\mu} \right)^\top \right] \\
&= \mathbb{E}_{H_0}\left[  \left(\mathbf{A} \left(\frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i \right) - \mathbf{A} \boldsymbol{\mu} \right)\left(\mathbf{A} \left(\frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i \right) - \mathbf{A} \boldsymbol{\mu} \right)^\top \right] \\
&= \mathbb{E}_{H_0}\left[ \mathbf{A} \left( \frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i - \boldsymbol{\mu} \right)\left(\frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i - \boldsymbol{\mu} \right)^\top \mathbf{A}^\top \right]\\
&=  \mathbf{A} \mathbb{E}_{H_0}\left[ \frac{1}{n^2} \sum_{i = 1}^n \sum_{j = 1}^n  \left( \mathbf{x}_i - \boldsymbol{\mu} \right) \left(\mathbf{x}_j - \boldsymbol{\mu} \right)^\top\right]  \mathbf{A}^\top \\
&=  \frac{1}{n^2} \sum_{i = 1}^n \sum_{j = 1}^n  \mathbf{A} \mathbb{E}\left[ \left( \mathbf{x}_i - \boldsymbol{\mu} \right) \left(\mathbf{x}_j - \boldsymbol{\mu} \right)^\top\right]  \mathbf{A}^\top \\
&=  \frac{1}{n^2} \sum_{i = 1}^n \sum_{j = 1}^n  \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top \\
&= \mathbf{A} \boldsymbol{\Sigma} \mathbf{A}^\top\\
&= \mathbb{I}_{k \times k}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>


Let $\mathbf{m} = \mathbf{A} \boldsymbol{\mu} = \mathbf{A} \mathbb{E}\left[ \mathbf{x} \right]$, and let $\mathbf{l} = \mathbf{A}^{-1}\boldsymbol{\mu}$. We can then rewrite Eq. \eqref{eq:chi-bar-squared} in terms of $\mathbf{y}$:

$$
\begin{equation}
\label{eq:rewrite-chi-bar}
\begin{aligned}
\bar{\chi}^2 &= n \left[ \mathbf{y}^\top \mathbf{y} - \underset{\substack{\mathbf{l}_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\mathbf{y}- \mathbf{m})^\top (\mathbf{y} - \mathbf{m}) \right\}\right] \\
&= n \left[ \rvert \rvert \mathbf{y} \rvert \rvert_2^2 -\underset{\substack{\mathbf{l}_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ \rvert \rvert \mathbf{y} - \mathbf{m} \rvert \rvert_2^2 \right\} \right]
\end{aligned}
\end{equation}
$$

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\bar{\chi}^2
&= n\left[ \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1}\bar{\mathbf{x}} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\}\right] \\
&= n\left[ \bar{\mathbf{x}}^\top \mathbf{A} \mathbf{A}^{-1} \boldsymbol{\Sigma}^{-1} \mathbf{A}^{-1} \mathbf{A} \bar{\mathbf{x}} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \mathbf{A} \mathbf{A}^{-1} \boldsymbol{\Sigma}^{-1} \mathbf{A}^{-1} \mathbf{A} (\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\}\right] \\
&= n\left[ \mathbf{y}^\top \mathbf{A}^{-1} \boldsymbol{\Sigma}^{-1} \mathbf{A}^{-1} \mathbf{y} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\mathbf{y}- \mathbf{A} \boldsymbol{\mu})^\top \mathbf{A}^{-1} \boldsymbol{\Sigma}^{-1} \mathbf{A}^{-1} (\mathbf{y} - \mathbf{A} \boldsymbol{\mu}) \right\}\right] \\
&= n\left[ \mathbf{y}^\top \left( \mathbf{A} \boldsymbol{\Sigma}  \mathbf{A}\right)^{-1} \mathbf{y} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\mathbf{y}- \mathbf{A} \boldsymbol{\mu})^\top \left( \mathbf{A} \boldsymbol{\Sigma}  \mathbf{A} \right)^{-1} (\mathbf{y} - \mathbf{A} \boldsymbol{\mu}) \right\}\right] \\
&= n \left[ \mathbf{y}^\top \mathbf{y} - \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\mathbf{y}- \mathbf{A} \boldsymbol{\mu})^\top (\mathbf{y} - \mathbf{A} \boldsymbol{\mu}) \right\}\right] & \left( \text{Eq. } \eqref{eq:sigma-identity}\right)\\
&= n \left[ \mathbf{y}^\top \mathbf{y} - \underset{\substack{(\mathbf{A}^{-1} \mathbf{m})_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\mathbf{y}- \mathbf{m})^\top (\mathbf{y} - \mathbf{m}) \right\}\right] & \left( \mathbf{m} = \mathbf{A} \boldsymbol{\mu} = \mathbf{A} \mathbb{E}\left[ \mathbf{x} \right]\right) \\
&= n \left[ \mathbf{y}^\top \mathbf{y} - \underset{\substack{\mathbf{l}_i \geq 0 \\ i = 1, \dots, k}}{\min} \left\{ (\mathbf{y}- \mathbf{m})^\top (\mathbf{y} - \mathbf{m}) \right\}\right] & \left( \mathbf{l} = \mathbf{A}^{-1} \boldsymbol{\mu} \right)
\end{aligned}
$$
</details>

Eq. \eqref{eq:rewrite-chi-bar} provides us insights into the geometric interpretation of $\bar{\chi}^2$. By the definition of $\mathbf{l}$, the feasible set can be rewritten as:

$$
\begin{equation}
\label{eq:cone-c}
\begin{aligned}
\mathcal{C} &= \left\{ \mathbf{m} \in \mathbb{R}^k : \mathbf{A}^{-1} \mathbf{m} \geq \mathbf{0}_k \right\} \\
&= \left\{ \mathbf{l} \in \mathbb{R}^k : \mathbf{l} \geq \mathbf{0}_k \right\}
\end{aligned}
\end{equation}
$$

The set, $\mathcal{C}$, is a <strong>polyhedral cone</strong><d-cite key=cone2025></d-cite>, and, clearly, this set is convex. It is also <a href="https://math.stackexchange.com/questions/1831401/how-do-you-prove-that-ax-mid-x-geq-0-is-closed">closed</a>. Thus, $\bar{\chi}^2$ can be interpreted as the difference (scaled by $n$) between the squared length of $\mathbf{y}$ and the distance between $\mathbf{y}$ and the cone, $\mathcal{C}$. 

<aside><p>Of course, with distance measured in the Euclidean sense.</p></aside>

In what follows, we will require a partitioning of this $\mathcal{C}$. Let $\mathbf{K} = \{1, 2, \dots, k \}$. We will divide $\mathcal{C}$ as:

$$
\mathcal{C} = \underset{\mathbf{M} \in \mathcal{P}(\mathbf{K})}{\bigcup} \mathcal{C}_\mathbf{M};
\hspace{8mm}
\mathcal{C}_{\mathbf{M}} = \left\{ 
  \mathbf{c} \in \mathcal{C} \hspace{2mm} : \hspace{2mm} 
  \begin{split}
    \mathbf{l}_i = 0 \text{ if }  i \notin \mathbf{M} \\
    \mathbf{l}_i > 0 \text{ if } i  \in    \mathbf{M}
  \end{split}
\right\}
$$

Note that $$\mathcal{C}_\mathbf{K}$$ consists of all points on the interior of $\mathcal{C}$, and $$\mathcal{C}_\emptyset$$ is just the singleton of $\mathbf{0}_k$. 

<aside><p>This is called the <a href="https://en.wikipedia.org/wiki/Linear_span"><strong>linear span</strong></a> or <strong>linear hull</strong> of $\mathcal{C}_{\mathbf{M}}$.</p></aside>

We also let $$\mathcal{B}_\mathbf{M}$$ denote the smallest <a href="https://en.wikipedia.org/wiki/Vector_space">linear space</a> such that $$\mathcal{C}_{\mathbf{M}} \subseteq \mathcal{B}_{\mathbf{M}}$$. By definition, $$\mathcal{B}_{\mathbf{M}}$$ will contain all finite linear combinations of elements in $$\mathcal{C}_{\mathbf{M}}$$. That is:

$$
\begin{aligned}
\mathcal{B}_{\mathbf{M}} &= \left\{ \left. \sum_{i = 1}^n \nu_i \mathbf{c}_i \right\rvert
  n \in \mathbb{N}; \mathbf{c}_1, \dots, \mathbf{c}_n \in \mathcal{C}_{\mathbf{M}}; \nu_1, \dots, \nu_n \in \mathbb{R}
\right\} \\
  &= \left\{ \mathbf{v} \in \mathbb{R}^k : \mathbf{v}_i = 0 \text{ for } i \notin \mathbf{M} \right\}
\end{aligned}
$$

<aside><p>Assuming the field that the space is defined over is the reals.</p></aside>

<details>
<summary>Linear Algebra Review.</summary>
Let $\mathcal{V}$ be a subspace of $\mathbb{R}^k$, and let $\{ \mathbf{v}_1, \dots, \mathbf{v}_k \}$ be a set of basis vectors for $\mathcal{V}$. Any other vector, $\mathbf{z} \in \mathbb{R}^k$, can be written <i>uniquely</i> as the sum of a vector in $\mathcal{V}$ and a vector that is orthogonal to $\mathcal{V}$:
$$
\mathbf{z} = \mathbf{z}_{\parallel \mathcal{V}} + \mathbf{z}_{\perp \mathcal{V}}
$$
This is called the <strong>orthogonal decomposition</strong> of $\mathbf{z}$ with respect to $\mathcal{V}$, and $\mathbf{z}_{\parallel \mathcal{V}}$ is the closest (in the Euclidean sense) vector to $\mathbf{z}$ that is in $\mathcal{V}$. This vector is computed as:
$$
\mathbf{z}_{\parallel \mathcal{V}} = \sum_{i = 1}^k \frac{\langle \mathbf{z}, \mathbf{v}_i \rangle}{\langle \mathbf{v}_i, \mathbf{v}_i \rangle} \mathbf{v}_i
$$
For a matrix $\mathbf{V} \in \mathbb{R}^k$, we define its <strong>column space</strong>, denoted by $\text{col}(\mathbf{V})$, as the set of all finite linear combinations of the column vectors of $\mathbf{V}$. Similarly, its <strong>row space</strong>, denoted by $\text{row}(\mathbf{V})$, is the set of all finite linear combinations of the row vectors of $\mathbf{V}$. The <strong>null space</strong> $\mathbf{V}$ is the set of vectors in $\mathbb{R}^k$ that are mapped to zero:
$$
\text{null}(\mathbf{V}) = \left\{ \mathbf{z} \in \mathbb{R}^k : \mathbf{V} \mathbf{z} = 0 \right\}
$$
An orthogonal complement of a subspace $\mathcal{V}$ of some other space $\mathcal{W}$ is the set of all vectors in the $\mathcal{W}$ that are orthogonal to every vector in $\mathcal{V}$. We have the following relationship between the row space, column space, and null space of $\mathbf{V}$:
$$
(\text{row}(\mathbf{V}))^\perp = \text{null}(\mathbf{V}) 
\hspace{5mm} \text{and} \hspace{5mm}
(\text{col}(\mathbf{V}))^\perp = \text{null}(\mathbf{V}^\top)
$$
where a superscript $\perp$ denotes the orthogonal complement.
<br>
Now, consider the column space of $\mathbf{V}$ and suppose $\text{col}(\mathbf{V}) = \mathcal{V}$. A vector $\mathbf{p}$ is in $\text{col}(\mathbf{V})$ if $\mathbf{A}\mathbf{c} = \mathbf{p}$ for some $\mathbf{c}$. It follows that $\mathbf{z}_{\parallel \mathcal{V}} = \mathbf{V} \mathbf{c}$ for some $\mathbf{c}$.
<br>
As dictated by the orthogonal decomposition of $\mathbf{z}$, we know that $\mathbf{z}_{\perp \mathcal{V}}$ must be orthogonal to $\mathcal{V}$, which is equivalent to being orthogonal to the column space of $\mathbf{V}$. As we just saw, the orthogonal complement of the column space of $\mathbf{V}$ is equal to the null space of $\mathbf{V}^\top$. Thus:
$$
\mathbf{V}^\top \mathbf{z}_{\perp \mathcal{V}} = \mathbf{0}_k
$$
The orthogonal decomposition (with some rearranging) states that $\mathbf{z}_{\perp \mathcal{V}} = \mathbf{z} - \mathbf{z}_{\parallel \mathcal{V}}$. Using this fact together with the two previous equations, we get:
$$
\begin{aligned}
\mathbf{0}_k &= \mathbf{V}^\top \left(\mathbf{z} - \mathbf{z}_{\parallel \mathcal{V}}\right) \\
\implies \mathbf{0}_k &= \mathbf{V}^\top \left(\mathbf{z} - \mathbf{V}\mathbf{c}\right) \\
\implies \mathbf{0}_k &= \mathbf{V}^\top \mathbf{z} - \mathbf{V}^\top \mathbf{V} \mathbf{c}\\
\implies \mathbf{V}^\top \mathbf{z} &= \mathbf{V}^\top \mathbf{V} \mathbf{c} \\
\implies \mathbf{c} &= \left(\mathbf{V}^\top \mathbf{V}\right)^{-1} \mathbf{V}^\top \mathbf{z}
\end{aligned}
$$
This brings us to the conclusion that:
$$
\mathbf{z}_{\parallel \mathbf{z}} = \mathbf{V} \left(\mathbf{V}^\top \mathbf{V} \right)^{-1} \mathbf{V}^\top \mathbf{z}
$$
The matrix $\mathbf{P} = \mathbf{V} \left(\mathbf{V}^\top \mathbf{V} \right)^{-1} \mathbf{V}^\top$ is sometimes called the <strong>projection matrix</strong>. If a vector $\mathbf{c} \in \text{col}(\mathbf{V})$, then $\mathbf{P}\mathbf{c} = \mathbf{c}$.
</details>

We will similarly partition the alternative parameter space, $\mathcal{R}$ as:

$$
\mathcal{R} = \underset{\mathbf{M} \subseteq \mathbf{K}}{\bigcup} \mathcal{R}_\mathbf{M};
\hspace{8mm}
\mathcal{R}_{\mathbf{M}} = \left\{ 
  \mathbf{r} \in \mathcal{R} \hspace{2mm} : \hspace{2mm} 
  \begin{split}
    \hat{\boldsymbol{\mu}}^*_i = 0 \text{ if }  i \notin \mathbf{M} \\
    \hat{\boldsymbol{\mu}}^*_i > 0 \text{ if } i  \in    \mathbf{M}
  \end{split}
\right\}
$$

$$\mathcal{R}_{\mathbf{K}}$$ is the set of all points where all components of the MLE, $$\hat{\boldsymbol{\mu}}^*$$, are positive, and $$\mathcal{R}_{\emptyset}$$ is the singleton $$\hat{\boldsymbol{\mu}}^* = \mathbf{0}_k$$ (implying that $\bar{\chi}^2 = 0$).

Let $$\mathbf{A}^{-1}_{i, \cdot}$$ denote the $i$-th row of $\mathbf{A}^{-1}$, and let $$\mathbf{A}^{-1}_{\mathbf{M}, \cdot}$$ denote the sub-matrix of rows indexed $i \in \mathbf{M}$ of $\mathbf{A}^{-1}$ for some set of indices $\mathbf{M}$. Without loss of generality, since we can reorder the indices of $\mathbf{K}$, suppose $$\mathbf{M} = \{m+1, \dots, k\}$$ for some $m > 0$. We will define $$\mathbf{a}_{i, \cdot}$$ for $i \in \mathbf{M}$ as the $k$-dimensional vector such that:

<ul>
<li>$\mathbf{a}_{i, \cdot} \in \text{span}(\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot})$</li>
<li>$\mathbf{a}_{i, \cdot} \perp \mathbf{A}^{-1}_{j, \cdot}$ for $j = 1, \dots, m$ and $i \neq j$</li>
<li>$\langle \mathbf{a}_{i, \cdot}, \mathbf{A}^{-1}_{i, \cdot} \rangle = 1$</li>
</ul>

We will also let $\mathbf{a}$ denote the matrix with rows equal to $$\mathbf{a}_{1, \cdot}, \dots, \mathbf{a}_{m, \cdot}$$. Since $\mathbf{a}_{i, \cdot}$ is in the row space of $\mathbf{A}^{-1}$, we can write $$\mathbf{a}^\top_{i, \cdot} = (\mathbf{A}^{-1})^\top \boldsymbol{\alpha}_i$$ for some coefficient vector $$\boldsymbol{\alpha}_i$$. 

The above conditions imply that $$\mathbf{A}^{-1}_{\mathbf{M}, \cdot} \mathbf{a}^\top_{i, \cdot} = \mathbf{e}_i$$ where $\mathbf{e}_i$ is the $i$-th standard basis vector in $m$ dimensions (all zeroes except the $i$-th coordinate which is equal to one). 

We also note that:

$$
\text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{a}_i, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right) = \text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)
$$

<details>
<summary>Proof.</summary>
By construction, $\mathbf{a}_{i, \cdot}$ is in the span of the first $m$ rows of $\mathbf{A}^{-1}$ but is orthogonal to all of its rows except the $i$-th one with which its inner product is one. Since it is orthogonal:
$$
\text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right) 
\subseteq 
\text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{a}_{i, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)
$$
Consider arbitrary row vector $\mathbf{v} \in \text{span}(\{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{a}_{i, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \})$. We can write $\mathbf{v}$ as the linear combination:
$$
\begin{aligned}
\mathbf{v}
  &= \nu_1 \mathbf{A}^{-1}_{1, \cdot} + \dots + \nu_i \mathbf{a}_{i, \cdot} + \dots + \nu_m \mathbf{A}^{-1}_{m, \cdot}
\end{aligned}
$$
However, because $\mathbf{a}_{i, \cdot} \in \text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)$, we can write:
$$
\begin{aligned}
\mathbf{v} 
  &= \nu_1 \mathbf{A}^{-1}_{1, \cdot} + \dots + \nu_i \mathbf{a}_{i, \cdot} + \dots + \nu_m \mathbf{A}^{-1}_{m, \cdot} \\
  &= \nu_1 \mathbf{A}^{-1}_{1, \cdot} + \dots + \nu_i \left(\lambda_1 \mathbf{A}^{-1}_{1, \cdot} + \dots + \lambda_i \mathbf{A}^{-1}_{i, \cdot} + \dots + \lambda_m \mathbf{A}^{-1}_{m, \cdot} \right) + \dots + \nu_m \mathbf{A}^{-1}_{m, \cdot} \\
  &= (\nu_1 + \nu_i \lambda_1) \mathbf{A}^{-1}_{1, \cdot} + \dots + \nu_i \lambda_i \mathbf{A}^{-1}_{i, \cdot} + \dots + (\nu_m + \nu_i \lambda_m) \mathbf{A}^{-1}_{m, \cdot}\\
\implies \mathbf{v} 
  &\in \text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)
\end{aligned}
$$
This implies that:
$$
\text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{a}_{i, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right) 
\subseteq 
\text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)
$$
And therefore that:
$$
\text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{a}_{i, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right) = \text{span}\left( \{\mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)
$$
</details>

In the next result, we use the following definitions. We will define an orthonormal system $\{ \mathbf{i}_1, \dots, \mathbf{i}_k \}$ as follows. Since any vector can be written as the sum of two orthogonal vectors, we we successively define the $\mathbf{i}_j$ vectors in the system to be those unit vectors that are orthogonal and satisfy:

$$
\begin{aligned}
(\mathbf{A}^{-1}_{1, \cdot})^\top &= c_{1,1} \mathbf{i}_1 \\
(\mathbf{A}^{-1}_{2, \cdot})^\top &= c_{2,1} \mathbf{i}_1 + c_{2,2} \mathbf{i}_2 \\
&\vdots \\
(\mathbf{A}^{-1}_{k, \cdot})^\top &= c_{k,1} \mathbf{i}_1 + c_{k,2} \mathbf{i}_2 + \dots + c_{k,k} \mathbf{i}_k
\end{aligned}
$$

We will also define:

$$
\begin{equation}
\label{eq:lemma-3-1-defns}
\begin{aligned}
\mathbf{I} &= \underbrace{\begin{bmatrix} 
\rvert & & \rvert \\
\mathbf{i}_1 & \dots & \mathbf{i}_k \\
\rvert & & \rvert
\end{bmatrix}}_{k \times k} 
&\mathbf{C} &= \underbrace{
  \begin{bmatrix}
  c_{1,1} & \dots & c_{k, 1} \\
  \vdots & \ddots & \vdots \\
  c_{1,k} & \dots & c_{k,k}
  \end{bmatrix}
}_{k \times k} \\
\mathbf{C}_{\cdot, 1:m} &= \underbrace{\begin{bmatrix}
  \rvert & & \rvert \\
  \mathbf{C}_{\cdot, 1} & \dots & \mathbf{C}_{\cdot, m} \\
  \rvert & & \rvert
\end{bmatrix}}_{k \times m} 
& \mathbf{C}_{1:m, 1:m} &= \underbrace{\begin{bmatrix}
  c_{1,1} & \dots & c_{m,1} \\
  \vdots & \ddots & \vdots \\
  c_{1,m} & \dots & c_{m,m}
  \end{bmatrix}
}_{m \times k}
\end{aligned}
\end{equation}
$$

<div id="kudo-lemma-3-1"></div>
<div class="theorem">
<strong>Lemma 3.1. </strong>
{% tabs kudo-3-1-lemma %}
{% tab kudo-3-1-lemma statement %}
For $$\mathbf{A}_{i, \cdot}^{-1}$$ and $$\mathbf{a}_{i, \cdot}$$ as defined above, the following holds:

$$
\left[\mathbf{A}_{\mathbf{M}, \cdot}^{-1}(\mathbf{A}_{\mathbf{M}, \cdot}^{-1})^\top \right]^{-1} = \mathbf{a}_{\mathbf{M}, \cdot} \mathbf{a}^\top_{\mathbf{M}, \cdot}
$$
{% endtab %}
{% tab kudo-3-1-lemma proof %}
We can write:

$$
\begin{aligned}
\begin{bmatrix}
\rvert & & \rvert \\
(\mathbf{A}^{-1}_{1, \cdot})^\top & \dots & (\mathbf{A}^{-1}_{m, \cdot})^\top \\
\rvert & & \rvert
\end{bmatrix} 
&= \begin{bmatrix}
\rvert & & \rvert \\
\mathbf{i}_1 & \dots & \mathbf{i}_k \\
\rvert & & \rvert
\end{bmatrix}
\begin{bmatrix}
c_{1,1} & \dots & c_{m,1} \\
\vdots & \ddots & \vdots \\
c_{1,k} & \dots & c_{m,k}
\end{bmatrix}
\\
&= \underbrace{\mathbf{I} \mathbf{C}_{\cdot, 1:m}}_{k \times m}
\end{aligned}
$$

Since all of the $$\mathbf{a}_{j, \cdot} \in \text{span}\left( \{ \mathbf{A}^{-1}_{1, \cdot}, \dots, \mathbf{A}^{-1}_{m, \cdot} \} \right)$$, we can also write each $$\mathbf{a}_{j, \cdot}$$ as a linear combination of the rows of $\mathbf{I}$. Let $\mathbf{D}$ be the matrix of the corresponding coefficients. Thus:

$$
\begin{bmatrix}
\rvert & & \rvert \\
\mathbf{a}_{1, \cdot}^\top & \dots & \mathbf{a}^\top_{m, \cdot}\\
\rvert & & \rvert
\end{bmatrix} 
= 
\underbrace{\mathbf{I} \mathbf{D}_{\cdot, 1:m}}_{k \times m}
$$

By the construction of our orthonormal system $$\{ \mathbf{i}_1, \dots, \mathbf{i}_k \}$$, the upper right entries of $\mathbf{C}$ are all zeros (i.e. $c_{i,j} = 0$ for $j > i$). This means that:

$$
\mathbf{C}^\top_{\cdot, 1:m} \mathbf{C}_{\cdot, 1:m}
=
\mathbf{C}^\top_{1:m, 1:m} \mathbf{C}_{1:m, 1:m}
$$

The same fact holds for $\mathbf{D}$. We then have:

$$
\left[ \mathbf{C}^\top_{\cdot, 1:m} \mathbf{C}_{\cdot, 1:m} \right]^{-1}
= 
\left[ \mathbf{C}^\top_{1:m, 1:m} \mathbf{C}_{1:m, 1:m} \right]^{-1} 
= \mathbf{C}_{1:m, 1:m}^{-1} [\mathbf{C}^{-1}_{1:m, 1:m}]^\top
$$

We then have:

$$
\begin{aligned}
\left[ \begin{bmatrix}
  - & \mathbf{A}^{-1}_1 & - \\
  & \vdots & \\
  - & \mathbf{A}^{-1}_m  & - 
\end{bmatrix}
\begin{bmatrix}
\rvert & & \rvert \\
(\mathbf{A}^{-1}_{1, \cdot})^\top & \dots & (\mathbf{A}^{-1}_{m, \cdot})^\top \\
\rvert & & \rvert
\end{bmatrix}
\right]^{-1}
= \left[ \left( \mathbf{I} \mathbf{C}_{\cdot, 1:m}  \right)^\top \mathbf{I} \mathbf{C}_{\cdot, 1:m} \right]^{-1} 
= \left[ \mathbf{C}_{\cdot, 1:m}^\top \underbrace{\mathbf{I}^\top \mathbf{I}}_{=\mathbb{I}_{k \times k}} \mathbf{C}_{\cdot, 1:m} \right]^{-1} 
= \left[ \mathbf{C}^\top_{1:m, 1:m} \mathbf{C}_{1:m,1:m}\right]^{-1}
\end{aligned}
$$

Recall that, by construction, $$\langle \mathbf{A}^{-1}_{j, \cdot}, \mathbf{a}_{l, \cdot} \rangle = 0$$ if $j \neq l$ and $1$ if $j = l$. Thus:

$$
\begin{bmatrix}
  - & \mathbf{A}^{-1}_{1, \cdot} & - \\
  & \vdots & \\
  - & \mathbf{A}^{-1}_{m, \cdot} & - 
\end{bmatrix}
\begin{bmatrix}
\rvert & & \rvert \\
\mathbf{a}_{1, \cdot}^\top & \dots & \mathbf{a}_{m, \cdot}^\top \\
\rvert & & \rvert
\end{bmatrix} 
= \mathbb{I}_{m \times m} 
= \left[  \mathbf{I} \mathbf{C}_{\cdot, 1:m} \right]^\top \mathbf{I} \mathbf{D}_{1:m} 
$$

Since $\mathbf{I}$ has orthonormal rows, $\mathbf{I}^\top \mathbf{I} = \mathbf{I} \mathbf{I}^\top = \mathbb{I}_{k \times k}$. Thus:

$$
\left[  \mathbf{I} \mathbf{C}_{\cdot, 1:m} \right]^\top \mathbf{I} \mathbf{D}_{\cdot, 1:m}
= \mathbf{C}_{\cdot, 1:m}^\top \mathbf{D}_{\cdot, 1:m}
= \mathbf{C}_{1:m, 1:m}^\top \mathbf{D}_{1:m, 1:m}
$$

However, we know that the last line must equal $$\mathbb{I}_{m \times m}$$, which implies that $$\mathbf{C}_{1:m, 1:m}^\top = \mathbf{D}^{-1}_{1:m, 1:m}$$. Thus:

$$
\left[\mathbf{A}^{-1}_{\mathbf{M}, \cdot} (\mathbf{A}^{-1}_{\mathbf{M}, \cdot})^\top\right]^{-1}
=
\left[\mathbf{C}^\top_{1:m, 1:m} \mathbf{C}_{1:m, 1:m} \right]^{-1}
=
\left[ \mathbf{D}^{-1}_{1:m, 1:m} (\mathbf{D}^{-1}_{1:m, 1:m})^\top \right]^{-1}
= \mathbf{D}^\top_{1:m, 1:m} \mathbf{D}_{1:m, 1:m}
$$

And then:

$$
\mathbf{D}^\top_{1:m, 1:m} \mathbf{D}_{1:m, 1:m}
= \mathbf{D}^\top_{\cdot, 1:m} \mathbf{D}_{\cdot, 1:m}
= \mathbf{D}^\top_{\cdot, 1:m} \underbrace{\mathbf{I}^\top \mathbf{I}}_{=\mathbb{I}_{k \times k}} \mathbf{D}_{\cdot, 1:m} 
= \left[ \mathbf{I} \mathbf{D}_{\cdot, 1:m} \right]^\top \mathbf{I} \mathbf{D}_{\cdot, 1:m} 
= \begin{bmatrix}
- & \mathbf{a}_{1, \cdot} & - \\
& \vdots & \\
- & \mathbf{a}_{m, \cdot} & - 
\end{bmatrix}
\begin{bmatrix}
\rvert & & \rvert \\
\mathbf{a}^\top_{1, \cdot} & \dots & \mathbf{a}^\top_{m, \cdot} \\
\rvert & & \rvert
\end{bmatrix}
= \mathbf{a}_{\mathbf{M}, \cdot} \mathbf{a}_{\mathbf{M}, \cdot}^\top
$$
{% endtab %}
{% endtabs %}
</div>

Let $$\mathbf{M}^c = \{ 1, \dots, k \} \setminus \mathbf{M}$$. We note that the nullspace of $$\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}$$ is equal to $$\mathcal{B}_{\mathbf{M}}$$. That is, every vector $$\mathbf{b} \in \mathcal{B}_{\mathbf{M}}$$ is orthogonal to the rows of $$\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}$$.

<details>
<summary>Proof.</summary>
Consider $\mathbf{a} \in \text{row}(\mathbf{A}_{\mathbf{M}^c, \cdot}^{-1})$ and $\mathbf{b} \in \mathcal{B}_{\mathbf{M}}$. We have:
$$
\begin{aligned}
\mathbf{a} &= 
  \sum_{i \in \mathbf{M}^c} \alpha_i (\mathbf{A}_{i, \cdot}^{-1})^\top 
  = \begin{bmatrix} 
  | & & | \\
  (\mathbf{A}_{1, \cdot}^{-1})^\top & \dots & (\mathbf{A}^{-1}_{m, \cdot})^\top \\
  | & & |
\end{bmatrix} \begin{bmatrix}
  \alpha_{1} \\
  \vdots \\
  \alpha_{m}
\end{bmatrix}\\
\mathbf{b} &= \sum_{i \in \mathbf{M}} \nu_i \mathbf{c}_i
 = \begin{bmatrix} 
  | & & | \\
  \mathbf{c}_{m+1} & \dots & \mathbf{c}_k \\
  | & & |
\end{bmatrix} \begin{bmatrix}
  \nu_{m+1} \\
  \vdots \\
  \nu_k
\end{bmatrix}
\end{aligned}
$$
Since $\mathbf{c}_{m+1}, \dots, \mathbf{c}_k \in \mathcal{C}_{\mathbf{M}}$, $\mathbf{b}$ satisfies:
$$
\begin{aligned}
\mathbf{A}_{\mathbf{M}, \cdot}^{-1}\mathbf{b}
  &= \mathbf{A}_{\mathbf{M}^c, \cdot}^{-1} \begin{bmatrix} 
  | & & | \\
  \mathbf{c}_{m+1} & \dots & \mathbf{c}_k \\
  | & & |
\end{bmatrix} \begin{bmatrix}
  \nu_{m+1} \\
  \vdots \\
  \nu_{k}
\end{bmatrix} \\
  &= \mathbf{0}_m
\end{aligned}
$$
Thus, $\mathbf{b}$ is in the null space of $\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}$. 
</details>

Because of the previous fact, we can project a vector $\mathbf{v} \in \mathbb{R}^k$ onto $$\mathcal{B}_{\mathbf{M}}$$ by projecting it onto the null space of $$\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}$$ (i.e. the orthogonal complement of the row space). For a subspace $\mathbf{W}$ and projection matrix $\mathbf{P}$, the projection matrix onto the orthogonal complement is given by $$\mathbf{P}^\perp = \mathbb{I} - \mathbf{P}$$. Since the row space of $$\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}$$ is equivalent to the column space of $$\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}$$, we can compute the projection matrix for $$\mathcal{B}_{\mathbf{M}}$$ as:

$$
\begin{aligned}
\mathbf{P}_{\mathbf{M}} 
  &= \mathbb{I}_{k \times k} - \mathbf{P}_{\text{col}([\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}]^\top)} \\
  &= \mathbb{I}_{k \times k} - (\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot})^\top [\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot} (\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot})^\top]^{-1} \mathbf{A}^{-1}_{\mathbf{M}^c, \cdot}
\end{aligned}
$$

For an arbitrary $\mathbf{m}_1 \in \mathcal{R}$, we denote the projection of $$\mathbf{y}_1 = \mathbf{A} \mathbf{m}_1$$ onto $$\mathcal{B}_{\{m+1, \dots, k\}}$$ as $\hat{\mathbf{y}}_1$, which is given by:

$$
\begin{equation}
\label{eq:proj-y-1}
\begin{aligned}
\hat{\mathbf{y}}_1 
  &= \mathbf{P}_{\mathbf{M}} \mathbf{y}_1 \\
  &= \left[  \mathbb{I}_{k \times k} - (\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot})^\top [\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot} (\mathbf{A}^{-1}_{\mathbf{M}^c, \cdot})^\top]^{-1} \mathbf{A}^{-1}_{\mathbf{M}^c, \cdot} \right] \mathbf{y}_1
\end{aligned}
\end{equation}
$$

This leads us to the following result.

<div class="theorem">
<strong>Properties. </strong>
{% tabs kudo-3-1 %}
{% tab kudo-3-1 statement %}
A necessary and sufficient condition for $$\mathbf{m}_1 \in \mathcal{R}_{\{m+1, \dots, k\}}$$ is that:

$$
\begin{align}
t_i &= \langle \mathbf{a}^{(i)} , \mathbf{y}_1 \rangle \leq 0 & i = 1, \dots, m \label{eq:t-i-cdtn}\\
s_i &= \langle \mathbf{a}_i, \hat{\mathbf{y}}_1 \rangle > 0 & i = m + 1, \dots, k \label{eq:s-i-cdtn}
\end{align}
$$
{% endtab %}
{% tab kudo-3-1 proof %}
To do.
{% endtab %}
{% endtabs %}
</div>


### Derivation

Since the parts of the partitions are disjoint, we can compute the tail probability over each part separately and take the sum. We first consider the case when $\emptyset \subsetneq mathbf{M} \subsetneq \mathbf{K}$. 

Recall the definitions of $\mathbf{I}$ and $\mathbf{C}$ as in Eq. \eqref{eq:lemma-3-1-defns} and that $$(\mathbf{A}^{-1})^\top = \mathbf{I} \mathbf{C}$$. We will define:

$$
\mathbf{z} = \mathbf{I}^\top \mathbf{y} \iff \mathbf{I} \mathbf{z} = \mathbf{y}
$$

<aside><p>$\mathbf{I}^\top = \mathbf{I}^{-1}$ for orthonormal $\mathbf{I}$.</p></aside>

We can also show that:

$$
\hat{\mathbf{z}} = \mathbf{I}^\top \hat{\mathbf{y}} 
\hspace{5mm} \text{where} \hspace{5mm} 
\hat{\mathbf{z}} = \begin{bmatrix}
0 \\
\vdots \\
0 \\
\mathbf{z}_{m+1} \\
\vdots \\
\mathbf{z}_k
\end{bmatrix}
$$

<aside><p>HOW DO YOU SHOW THIS OMG?</p></aside>

Since the columns of $\mathbf{I}$ are orthonormal, we have that the coordinates of $\mathbf{z}$ are independent and have unit variance (as linear combinations of these orthonormal vectors). We can then write:

<aside><p>CHECK INDEPENDENCE CLAIM!</p></aside>

$$
\mathbf{m} = \mathbf{A}^{-1} \mathbf{y} = (\mathbf{I} \mathbf{C})^\top \mathbf{I}\mathbf{z} = \mathbf{C}^\top \mathbf{z}
\hspace{5mm}
\implies
\hspace{5mm}
\frac{1}{n} \boldsymbol{\Sigma} = \mathbf{C}^\top \mathbf{C}
$$

<details>
<summary>Proof.</summary>
$$
\text{var}\left( \mathbf{m} \right)
  = \text{var}\left(\mathbf{C}^\top  \mathbf{z} \right) 
  = \mathbf{C}^\top  \text{var}(\mathbf{z}) \mathbf{C}
  = \mathbf{C}^\top \mathbf{C}
$$
However, we also know that:
$$
\text{var}\left(\mathbf{m}\right)
  = \text{var}\left(\frac{1}{n} \sum_{i = 1}^n \mathbf{x}_i \right) 
  = \frac{1}{n^2} \sum_{i = 1}^n \text{var}(\mathbf{x}_i) 
  = \frac{1}{n} \boldsymbol{\Sigma}
$$
And therefore:
$$
\frac{1}{n} \boldsymbol{\Sigma} = \mathbf{C}^\top  \mathbf{C}
$$
</details>

Let us define:

$$
\begin{equation}
\label{eq:c-i-sup}
\mathbf{C}^{\cdot, i}  = \mathbf{I}^\top \mathbf{a}_{i, \cdot}^\top = \begin{bmatrix}
\mathbf{i}_1^\top \mathbf{a}_{i, \cdot} \\
\vdots \\
\mathbf{i}_k^\top \mathbf{a}_{i, \cdot}
\end{bmatrix}
\end{equation}
$$

It is fairly simple to show that the coordinates of $\mathbf{C}^{\cdot, i}$ are linear combinations of the columns of $\mathbf{C}$ and that:

$$
\begin{equation}
\label{eq:c-sup-c-sub-orthog}
\begin{aligned}
\langle \mathbf{C}^{\cdot, i}, \mathbf{C}_{\cdot, j} \rangle &= \begin{cases}
0 & i \neq j \\
1 & i = j
\end{cases}
&(\mathbf{C}^{\cdot, \mathbf{M}^c})^\top \mathbf{C}_{\cdot, \mathbf{M}^c}
  &= \mathbb{I}_{m, m} \\
\mathbf{a}_{i, \cdot}^\top &= \mathbf{I} \mathbf{C}^{\cdot, i}
&\mathbf{C} &= \mathbf{I}^\top (\mathbf{A}^{-1})^\top
\end{aligned}
\end{equation}
$$

<details>
<summary>Proof.</summary>
Since $(\mathbf{a}_{i, \cdot})^\top = (\mathbf{A}^{-1})^\top \boldsymbol{\alpha}_i$ for some vector of appropriate coefficients $\boldsymbol{\alpha}_i$, and the columns of $\mathbf{I}$ are orthonormal, we see that:
$$
\begin{aligned}
  \mathbf{C}^{\cdot, i} 
    = \mathbf{I}^\top (\mathbf{a}_{i, \cdot})^\top
    = \mathbf{I}^\top (\mathbf{A}^{-1})^\top \boldsymbol{\alpha}_i
    = \mathbf{I}^\top \mathbf{I} \mathbf{C} \boldsymbol{\alpha}_i 
    = \mathbf{C} \boldsymbol{\alpha}_i
\end{aligned}
$$
where the second to last equality follows from the fact that $\mathbf{I}\mathbf{I}^\top = \mathbf{I}^\top \mathbf{I} = \mathbf{I}^{-1}$ for an orthonormal matrix. We also see that:
$$
\begin{aligned}
(\mathbf{C}^{\cdot, i})^\top \mathbf{C}_{\cdot, j}
  =  (\mathbf{I}^\top (\mathbf{a}_{i, \cdot})^\top)^\top \mathbf{I}^{-1} (\mathbf{A}^{-1}_{j, \cdot})^\top
  = \mathbf{a}_{i, \cdot} \mathbf{I} \mathbf{I}^{-1} (\mathbf{A}^{-1}_{j, \cdot})^\top
  = \begin{cases}
  0 & i \neq j \\
  1 & i = j
  \end{cases}
\end{aligned}
$$
We then can see that, for any set of indices $\mathbf{B}$ where $\mathbf{B}_i$ denotes its $i$-th element and $\rvert \mathbf{B} \rvert = b$, we get:
$$
(\mathbf{C}^{\cdot, \mathbf{B}})^\top \mathbf{C}_{\cdot, \mathbf{B}}
  = \begin{bmatrix}
  \langle \mathbf{C}^{\cdot, \mathbf{B}_1}, \mathbf{C}_{\cdot, \mathbf{B}_1} \rangle & \dots & \langle \mathbf{C}^{\cdot, \mathbf{B}_1}, \mathbf{C}_{\cdot, \mathbf{B}_b}  \rangle \\
  \vdots & \ddots & \vdots \\
  \langle \mathbf{C}^{\cdot, \mathbf{B}_b}, \mathbf{C}_{\cdot, \mathbf{B}_1} \rangle & \dots & \langle \mathbf{C}^{\cdot, \mathbf{B}_b}, \mathbf{C}_{\cdot, \mathbf{B}_b}  \rangle
  \end{bmatrix} = \mathbb{I}_{b \times b}
$$
</details>

For $i \in \mathbf{M}$, we have:

$$
\begin{aligned}
s_i 
  &= \langle (\mathbf{A}^{-1}_{i, \cdot})^\top, \hat{\mathbf{y}} \rangle
  = \mathbf{A}^{-1}_{i, \cdot} \hat{\mathbf{y}} 
  = (\mathbf{I} \mathbf{C}_{\cdot, i})^\top \hat{\mathbf{y}}
  = \mathbf{C}_{\cdot, i}^\top \mathbf{I}^\top \hat{\mathbf{y}}
  = \mathbf{C}_{\cdot, i}^\top \hat{\mathbf{z}} 
  = \langle \mathbf{C}_{\cdot, i}, \hat{\mathbf{z}} \rangle \\
t_i 
  &= \langle (\mathbf{a}_{i, \cdot})^\top, \mathbf{y} \rangle 
  = \mathbf{a}_{i, \cdot} \mathbf{y} 
  = (\mathbf{I} \mathbf{C}^{\cdot, i})^\top \mathbf{y} 
  = (\mathbf{C}^{\cdot, i})^\top \mathbf{I}^\top \mathbf{y} 
  = \langle \mathbf{C}^{\cdot, i}, \mathbf{z} \rangle
\end{aligned}
$$ 



---

## Maximum Likelihood Estimation
If $\mathbf{y} \notin \mathcal{C}$, then the second term does not vanish. Because $\mathcal{C}$ is closed, at least one minimum will exist, and because it is also convex, there will only be a unique solution. 

<aside><p>We have a quadratic program over a closed, convex set.</p></aside>

Denote this solution by $$\hat{\boldsymbol{\mu}}^*$$. It is fairly easy to see that $$\hat{\boldsymbol{\mu}}^*$$ is the MLE of $\boldsymbol{\mu}$ under the alternative hypothesis. 

<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\hat{\boldsymbol{\mu}}^* 
&= \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\arg \min} \left\{ (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\} \\
&= \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\arg \max} \left\{ - \frac{n}{2} (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\} \\
&= \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\arg \max} \left\{ - \frac{n}{2} \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1} \bar{\mathbf{x}} - \frac{n}{2} (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right\} \\ 
&= \underset{\substack{\mu_i \geq 0 \\ i = 1, \dots, k}}{\arg \max} \left\{ \frac{1}{(\sqrt{2 \pi})^{kn} \rvert \boldsymbol{\Sigma} \rvert^n} \exp \left( - \frac{n}{2} \bar{\mathbf{x}}^\top \boldsymbol{\Sigma}^{-1} \bar{\mathbf{x}} - \frac{n}{2} (\bar{\mathbf{x}} - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \boldsymbol{\mu}) \right) \right\} \\
\end{aligned}
$$
</details>

We can then derive the following result regarding $$\hat{\boldsymbol{\mu}}^*$$.

<div class="theorem">
<strong>Theorem 2.1.</strong><d-cite key=kudo1963></d-cite>
{% tabs kudo21 %}
{% tab kudo21 statement %}
A necessary and sufficient condition for $$\hat{\boldsymbol{\mu}}^*$$ to be the MLE is: 

For any $i \in [k]$, one of $$\hat{\boldsymbol{\mu}}^*_i$$ and $$\left[ - \boldsymbol{\Sigma}^{-1} (\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*) \right]_i$$ is equal to zero and the other is non-negative.
{% endtab %}
{% tab kudo21 proof %}
Consider some other $\hat{\boldsymbol{\mu}} \in \mathcal{C}$ (i.e. $\hat{\boldsymbol{\mu}}$ such that $\hat{\mu}_i \geq 0$ for all $i$). Consider the difference in the negative log-likelihoods evaluated at $\hat{\boldsymbol{\mu}}$ and $\hat{\boldsymbol{\mu}}^*$:

$$
(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}})^\top\boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}) - (\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*) \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*)
$$

where we have dropped the terms that do not depend on the parameter vector. With some rearranging, this is equivalent to:

$$
\underbrace{(\hat{\boldsymbol{\mu}} -\hat{\boldsymbol{\mu}}^*)^\top \boldsymbol{\Sigma}^{-1} (\hat{\boldsymbol{\mu}} - \hat{\boldsymbol{\mu}}^*)}_{(a)} + \underbrace{-2 (\hat{\boldsymbol{\mu}} - \hat{\boldsymbol{\mu}}^*)^\top \boldsymbol{\Sigma}^{-1} (\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*)}_{(b)}
$$

Since we assumed $\boldsymbol{\Sigma}$ to be invertible, $\boldsymbol{\Sigma}^{-1}$ is also (trivially) invertible. Thus, both are positive definite, which implies that $(a)$ is strictly positive. We now need to show that the second term is non-negative to get to our desired result. First, we will show sufficiency, so assume the condition is true. 

Notice that we can define $$\mathbf{v}= -2 \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*)$$ and $$\mathbf{u} = \hat{\boldsymbol{\mu}} - \hat{\boldsymbol{\mu}}^*$$, so that $$(b) = \mathbf{v}^\top \mathbf{u}$$. 

Suppose $\mathbf{v}_i \geq 0$. This implies that $$\left[-\boldsymbol{\Sigma} (\bar{\mathbf{x}} - \bar{\boldsymbol{\mu}}^*)\right]_i \geq 0$$. By assumption, $$\hat{\boldsymbol{\mu}}^*_i = 0$$, which implies that:

$$
\mathbf{u}_i = \hat{\boldsymbol{\mu}}_i - \underbrace{\hat{\boldsymbol{\mu}}^*_i}_{=0} \geq 0
$$

because $\hat{\boldsymbol{\mu}}_i \geq 0$ by construction. Thus, $\mathbf{v}_i \mathbf{u}_i \geq 0$. Any coordinate $j$ of $\mathbf{v}$ that is not strictly positive must be equal to zero (by assumption), so that $\mathbf{v}_j \mathbf{u}_j = 0$. 

Thus, $(b) \geq 0$, and the entire quantity is positive. Since this implies that the likelihood evaluated at $\hat{\boldsymbol{\mu}}^*$ is greater than the likelihood evaluated at any other feasible $\hat{\boldsymbol{\mu}}$, we conclude it is the MLE. 

We now show necessity. Assume that $$\hat{\boldsymbol{\mu}}^*$$ is the MLE. Note that:

$$
\left. \frac{\partial}{\partial \mathbf{z}} \left[ (\bar{\mathbf{x}} - \mathbf{z})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \mathbf{z}) \right] \right\rvert_{\mathbf{z} = \hat{\boldsymbol{\mu}}^*} = -2 \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*) 
$$

I do not really understand the remainder of Kudô's proof, but he states that if the condition is violated, then we can construct $\hat{\boldsymbol{\mu}}$ close to $$\boldsymbol{\mu}^*$$ such that each $\hat{\boldsymbol{\mu}}_i \geq 0$ such that $$(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}})^\top \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}})$$ is smaller than that of $$\hat{\boldsymbol{\mu}}^*$$. This contradicts the assumption that $$\hat{\boldsymbol{\mu}}^*$$ is the MLE. 

<!-- 
Suppose the condition is not true, which means that either both quantities are strictly positive or one is strictly negative while the other one is permitted to be whatever. Consider the first case in which $$\left[ -2 \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*)\right]_i > 0$$. This implies that $\hat{\boldsymbol{\mu}}^*$ is not a stationary point of the likelihood and thus cannot be the MLE, which is a contradiction. 

Now consider the second case. If $$\left[ -2 \boldsymbol{\Sigma}^{-1}(\bar{\mathbf{x}} - \hat{\boldsymbol{\mu}}^*)\right]_i < 0$$, then we arrive at a contradiction by the same argument as the first case. If $$\hat{\boldsymbol{\mu}}^* < 0$$, then  -->
{% endtab %}
{% endtabs %}
</div>
















<!-- 





<strong>Necessity.</strong>
Suppose $$\mathbf{m}_1 \in \mathcal{R}_{\{m+1, \dots, k \}}$$. We want to show that it must be the case that Conditions \eqref{eq:t-i-cdtn} and \eqref{eq:s-i-cdtn} hold. 

By definition of $$\mathcal{R}_{\{m+1, \dots, k\}}$$, the first $m$ entries of $\mathbf{m}_1$ are zero, and the rest are positive. Thus:

$$
\begin{equation}
\label{eq:a-inv-m}
\mathbf{y}_1 = \mathbf{A} \mathbf{m}_1
\implies
\mathbf{A}^{-1} \mathbf{y}_1 = \mathbf{m}_1 = \begin{bmatrix}
0 \\ \vdots \\ 0 \\ \mathbf{m}_{1,m+1} \\ \vdots \\ \mathbf{m}_{1, k}
\end{bmatrix}
\end{equation}
$$

This means $$\mathbf{y}_1 \in \mathcal{C}_{\{m+1, \dots, k\}} \subseteq \mathcal{C}$$. This implies that $\hat{\mathbf{y}}_1 = \mathbf{y}_1$ and therefore that $\langle \mathbf{A}^{-1}_i, \hat{\mathbf{y}}_1 \rangle > 0$ for any $i \in \{ m+1, \dots, k \}$. This shows Condition \eqref{eq:s-i-cdtn}.

To show Condition \eqref{eq:t-i-cdtn}, suppose there exists some $$l \in \{1, \dots, m\}$$ such that $\langle \mathbf{a}^{(l)}, \mathbf{y}_1 \rangle > 0$. 
















































Consider projecting $\mathbf{y}_1$ onto $$\mathcal{B}_{\{l, m+1, \dots, k \}}$$ (for some $$l \in \\{ 1, \dots, m\\}$$), which we denote with $\hat{\mathbf{y}}_1^{l}$. We can then write $$\hat{\mathbf{m}}_1^l = \mathbf{A}^{-1} \hat{\mathbf{y}}_1^l$$. 

By the fact that $$\hat{\mathbf{y}}_1^l \in \mathcal{B}_{\{l, m+1, \dots, k \}}$$, we know that $$\mathbf{A}^{-1}_i \hat{\mathbf{y}}_{1}^l = 0$$ for $$i = \{ 1, \dots m \} \setminus \{ l \}$$, which implies that $$\hat{\mathbf{m}}_{1,i}^l = 0$$. 


<i>Case 1.</i>
Suppose first that $$\hat{\mathbf{m}}_{1,i}^l > 0$$ for $$i \in \\{ l, m+1, \dots, k \\}$$. In this case, $$\hat{\mathbf{y}}_{1,i}^l = 0$$ for $i \notin \\{l, m+1, \dots, k\\}$ and $> 0$ otherwise. Thus, $$\hat{\mathbf{y}}_1^l \in \mathcal{C}_{\{l, m+1, \dots, k\}}$$. 

By the definition of the linear subspaces, any $$\mathbf{b} \in \mathcal{B}_{\{ m+1, \dots, k \}}$$ is also in $$\mathcal{B}_{\{ l, m+1, \dots, k \}}$$ (since we can set the appropriate coefficients in the linear combination to zero). Since we can also choose the appropriate coefficient values to get $\mathbf{b}$ such that $$\mathbf{b} \in \mathcal{B}_{\{ m+1, \dots, k \}}$$ and $$\mathbf{b} \notin \mathcal{B}_{\{ m+1, \dots, k \}}$$, we have that $$\mathcal{B}_{\{ m+1, \dots, k \}} \subset \mathcal{B}_{\{ l, m+1, \dots, k \}}$$. 

<!-- For some $$l \in \{ 1, \dots, m\}$$, let $$\hat{\mathbf{y}}_1^{(l)}$$ denote the orthogonal projection of $\mathbf{y}_1$ onto $$\mathcal{B}_{\{l, m+1, \dots, k\}}$$, let $$\hat{\mathbf{m}}_1^{(l)} = \mathbf{A}^{-1} \hat{\mathbf{y}}_1^{(l)}$$, and let $$\hat{\mathbf{m}}_{1,i}^{(l)}$$ denote the $i$-th coordinate of $$\hat{\mathbf{m}}^{(l)}_1$$. Since $$\hat{\mathbf{y}}_1^{(l)} \in \mathcal{B}_{\{l, m+1, \dots, k\}}$$, we know that:

$$
\hat{\mathbf{m}}_{1,j}^{(l)} = [\mathbf{A}^{-1} \hat{\mathbf{y}}^{(l)}_1]_j = 0 \hspace{4mm} \text{ for } j \notin \{ l, m+1, \dots, k\}
$$

All that remains is to show that $$\hat{\mathbf{m}}^{(l)}_1 = \langle \mathbf{a}^{(l)}, \mathbf{y}_1 \rangle \leq 0$$. We will proceed by contradiction in cases. Assume that $$\hat{\mathbf{x}}_{1,l}^{(l)} > 0$$.  -->

This fact implies that:

$$
\rvert \rvert \hat{\mathbf{y}}^{l}_1 - \mathbf{y} \rvert \rvert_2^2 \leq \rvert \rvert \hat{\mathbf{y}}_1 - \mathbf{y} \rvert \rvert_2^2
$$

since the orthogonal projection onto a subset cannot be farther than the orthogonal projection onto the superset.

Also recall that we assumed $$\mathbf{m}_1 \in \mathcal{R}_{\{m+1, \dots, k\}}$$, which implied that $$\mathbf{A}_i^{-1} \hat{\mathbf{y}}_1 = 0$$ for all $$i \notin \{m+1, \dots, k\}$$. Clearly, $$\hat{\mathbf{y}}_1^{l} \neq \hat{\mathbf{y}}_1$$ since the $l$-th element is strictly zero for $\hat{\mathbf{y}}_1$ but is positive for $$\hat{\mathbf{y}}_1^{l}$$.

This fact implies strict inequality:

$$
\begin{equation}
\label{eq:contradiction-ineq}
\rvert \rvert \hat{\mathbf{y}}^{(l)}_1 - \mathbf{y} \rvert \rvert_2^2 < \rvert \rvert \hat{\mathbf{y}}_1 - \mathbf{y} \rvert \rvert_2^2 
\end{equation}
$$

since the projections are not the same. By Eq. \eqref{eq:a-inv-m}, $\hat{\mathbf{y}}_1$ is the closest point in $\mathcal{C}$ to $\mathbf{y}_1$. However, by Eq. \eqref{eq:contradiction-ineq}, $$\hat{\mathbf{y}}^{l}_1$$ is a closer point, which contradicts the assumption that $$\mathbf{m}_1 \in \mathcal{R}_{\{m+1, \dots, k\}}$$. 

<i>Case 2.</i> Again assume $$\hat{\mathbf{x}}_{1,l}^{l} > 0$$ but also that for some $i = m+1, \dots, k$:

$$
\begin{equation}
\label{eq:case-2-assumption}
\hat{\mathbf{m}}^{l}_i \leq 0
\end{equation}
$$






Suppose we interpolate between the projection of $\mathbf{y}_1$ on $\mathcal{B}_{\{m+1, \dots, k\}}$ and on $\mathcal{B}_{\{l, \dots, k\}}$. This is just the line:

$$
\mathbf{y}(p) = p \hat{\mathbf{y}}_1 + (1 - p) \hat{\mathbf{y}}^{(l)}_1
\hspace{4mm}
\text{for } 
p \in [0, 1]
$$

We can also construct a corresponding line between $\mathbf{A}^{-1} \hat{\mathbf{y}}_1 = \hat{\mathbf{m}}_1$ and $\mathbf{A}^{-1}\hat{\mathbf{y}}_1^{(l)} = \hat{\mathbf{m}}_1^{(l)}$:

$$
\mathbf{m}(p) = p \hat{\mathbf{m}}_1 + (1 - p) \hat{\mathbf{m}}^{(l)}_1
\hspace{4mm}
\text{for } 
p \in [0, 1]
$$

Consider the $i$-th coordinate of $\mathbf{m}(p)$, which we denote with $\mathbf{m}_i(p)$. For $i = 1, \dots, m$ such that $i \neq l$, we have:

$$
\begin{aligned}
\mathbf{m}_i(p) 
  &= p \underbrace{[\mathbf{A}^{-1} \hat{\mathbf{y}}_1]_i}_{=0} + (1 - p) \underbrace{[ \mathbf{A}^{-1} \hat{\mathbf{y}}_1^{(l)}]_i}_{=0} \\
  &= 0
\end{aligned}
$$

For $i = l$, we have:

$$
\begin{aligned}
\mathbf{m}_l(p) 
  &= p \underbrace{[\mathbf{A}^{-1} \hat{\mathbf{y}}_1]_l}_{=0} + (1 - p) \underbrace{[ \mathbf{A}^{-1} \hat{\mathbf{y}}_1^{(l)}]_l}_{\geq 0} \\
  &= (1 - p) \hat{\mathbf{m}}_l^{(l)} \\
  &\geq 0 
\end{aligned}
$$

For $i = m+1, \dots, k$, we have:

$$
\begin{aligned}
\mathbf{m}_i(p) 
  &= p \underbrace{[\mathbf{A}^{-1} \hat{\mathbf{y}}_1]_i}_{>0} + (1 - p) [ \mathbf{A}^{-1} \hat{\mathbf{y}}_1^{(l)}]_i \\
  &= p \hat{\mathbf{m}}_i + (1 - p) \hat{\mathbf{m}}_i^{(l)}
\end{aligned}
$$

If $p = 1$, then $\mathbf{m}_i(p)$ is positive for all $i = m+1, \dots, k$. If $p = 0$, then $\mathbf{m}_i(p) \leq 0$ for some $i = m+1, \dots, k$ by the assumption in Eq. \eqref{eq:case-2-assumption}. For $p_0$ such that $0 < p_0 < 1$, 





























<strong>Sufficiency.</strong>
Assume that Conditions \eqref{eq:t-i-cdtn} and \eqref{eq:s-i-cdtn} hold. The orthogonal decomposition of $\mathbf{y}_1$ with respect to $$\mathcal{B}_{\mathbf{M}}$$ is given by:

$$
\mathbf{y}_1 = \hat{\mathbf{y}}_1 + \hat{\mathbf{y}}_1^\perp
$$

where $\hat{\mathbf{y}}_1$ lies on $$\mathcal{B}_{\mathbf{M}}$$, and $$\hat{\mathbf{y}}_1^\perp$$ lies in the orthogonal complement of $$\mathcal{B}_{\mathbf{M}}$$ (which is the same as the row space of $$\mathbf{A}^{-1}_{\mathbf{M}}$$). 

Using Eq. \eqref{eq:proj-y-1} and <a href="#kudo-lemma-3-1">Lemma 3.1</a>, we see:

$$
\begin{aligned}
\hat{\mathbf{y}}_1
  &= \mathbf{y}_1 -  (\mathbf{A}^{-1}_{\mathbf{M}})^\top [\mathbf{A}^{-1}_{\mathbf{M}} (\mathbf{A}^{-1}_\mathbf{M})^\top]^{-1} \mathbf{A}^{-1}_{\mathbf{M}}\mathbf{y}_1 \\
    &= \mathbf{y}_1 -  (\mathbf{A}^{-1}_{\mathbf{M}})^\top 
    \begin{bmatrix}
    - & (\mathbf{a}^{(1)})^\top - \\
    & \vdots & \\
    - & (\mathbf{a}^{(m)})^\top - 
    \end{bmatrix}
    \begin{bmatrix}
    | & & | \\
    \mathbf{a}^{(1)} & \dots & \mathbf{a}^{(m)} \\
    | & & |
    \end{bmatrix}
    \mathbf{A}^{-1}_{\mathbf{M}}\mathbf{y}_1 \\
  &= \mathbf{y}_1 - (\mathbf{A}^{-1}_\mathbf{M})^\top (\mathbf{a}^{(1:m)})^\top \mathbf{a}^{(1:m)} \mathbf{A}^{-1}_\mathbf{M} \mathbf{y}_1
\end{aligned}
$$

where we let $\mathbf{a}^{(1:m)}$ denote the matrix whose columns are the vectors $$\{ \mathbf{a}^{(i)} \}_{i = 1}^m$$.

Let us consider the equation of the plane that passes through $\hat{\mathbf{y}}_1$ and is orthogonal to $-\hat{\mathbf{y}}_1^\perp$. This is the plane that divides the space into two parts where $$\mathcal{B}_{\mathbf{M}}$$ lies in the halfspace corresponding to the positive side of the plane. For $\mathbf{v} \in \mathbb{R}^k$, we have:

$$
\begin{aligned}
0 &= \underbrace{(\hat{\mathbf{y}}_1 - \mathbf{y}_1)}_{= - \hat{\mathbf{y}}_1^\perp}^\top (\mathbf{v} - \hat{\mathbf{y}}_1) \\
\implies 0 &= (\hat{\mathbf{y}}_1 - \mathbf{y}_1)^\top \mathbf{v} - \underbrace{(\hat{\mathbf{y}}_1^\perp)^\top \hat{\mathbf{y}}_1}_{\hat{\mathbf{y}}_1^\perp \perp \hat{\mathbf{y}}_1} \\
\implies 0 &= (\hat{\mathbf{y}}_1 - \mathbf{y}_1)^\top \mathbf{v}
\end{aligned}
$$

And thus:

$$
\begin{aligned}
(\hat{\mathbf{y}}_1 - \mathbf{y}_1)^\top \mathbf{v}
  &= (\mathbf{y}_1  - (\mathbf{A}^{-1}_\mathbf{M})^\top (\mathbf{a}^{(1:m)})^\top \mathbf{a}^{(1:m)} \mathbf{A}^{-1}_\mathbf{M} \mathbf{y}_1 - \mathbf{y}_1)^\top \mathbf{v} \\
  &= - \mathbf{y}_1^\top (\mathbf{A}^{-1}_\mathbf{M})^\top (\mathbf{a}^{(1:m)})^\top \mathbf{a}^{(1:m)} \mathbf{A}^{-1}_\mathbf{M} \mathbf{v}\\
  &= -  \left[ \mathbf{a}^{(1:m)} \mathbf{A}^{-1}_\mathbf{M} \mathbf{y}_1 \right]^\top \left[ \mathbf{a}^{(1:m)} \mathbf{A}^{-1}_\mathbf{M} \mathbf{v}\right] 
\end{aligned}
$$

If Condition \eqref{eq:s-i-cdtn} holds, then:

$$
\begin{aligned}
\mathbf{0}_k &< \mathbf{A}^{-1}_{\mathbf{M}} \hat{\mathbf{y}}_1 \\
\implies \mathbf{0}_k &< \mathbf{A}^{-1}_{\mathbf{M}} \left[ \mathbf{y}_1 - (\mathbf{A}^{-1}_\mathbf{M})^\top (\mathbf{a}^{(1:m)})^\top \mathbf{a}^{(1:m)} \mathbf{A}^{-1}_\mathbf{M} \mathbf{y}_1\right] \\
\implies \mathbf{0}_k &< \left[ \mathbb{I}_{k \times k} - (\mathbf{A}^{-1}_\mathbf{M})^\top (\mathbf{a}^{(1:m)})^\top \mathbf{a}^{(1:m)} \right] \mathbf{A}^{-1}_{\mathbf{M}} \mathbf{y}_1
\end{aligned}
$$


Now consider any $\mathbf{v} \in \mathcal{C}$. Since we assumed Conditions \eqref{eq:t-i-cdtn} and \eqref{eq:s-i-cdtn} hold, we have:
































For $i = 1, \dots, m$:

$$
\begin{aligned}
\langle \mathbf{a}^{(i)}, \mathbf{y}_1 \rangle &\leq 0 \\
\implies \langle \mathbf{a}^{(i)}, \mathbf{A} \mathbf{x}_1 \rangle &\leq 0 \\
\implies \langle \mathbf{A}^\top \mathbf{a}^{(i)}, \mathbf{x}_1 \rangle &\leq 0
\end{aligned}
$$






For a fixed $l = 1, \dots, m$, we can write $\mathbf{A}^{-1}_j$ for $j = 1, \dots, k$ as:

$$
\begin{aligned}
\mathbf{A}^{-1}_j &= \lambda_{j,1} \mathbf{A}^{-1}_{j, 1} + \dots + \lambda_{j, l-1} \mathbf{A}^{-1}_{l - 1} + \lambda_{j, l} \mathbf{a}^{(l)} + \lambda_{j, l + 1} \mathbf{A}^{-1}_{l + 1} + \dots + \lambda_{j,m} \mathbf{A}^{-1}_m + \mathbf{b}_j \\
\end{aligned}
$$

where $\mathbf{b}_j$ is the difference between $\mathbf{A}^{-1}_j$ and its orthogonal projection onto $$\text{span}\left(\mathbf{A}^{-1}_1, \dots, \mathbf{A}^{-1}_{i - 1}, \mathbf{a}^{(i)}, \mathbf{A}_{i + 1}^{-1}, \dots, \mathbf{A}_{m}^{-1} \right)$$. Note that $\mathbf{b}_j = \mathbf{0}_k$ for $j = 1, \dots, m$ since the vector is within the span. 

Consider some $l = 1, \dots, m$ and let $j = l$. We have that:

$$
\begin{aligned}
\mathbf{A}^{-1}_l
  &= \lambda_{l,1} \mathbf{A}^{-1}_{l, 1} + \dots + \lambda_{l, l-1} \mathbf{A}^{-1}_{l - 1} + \lambda_{l, l} \mathbf{a}^{(l)} + \lambda_{l, l + 1} \mathbf{A}^{-1}_{l + 1} + \dots + \lambda_{l,m} \mathbf{A}^{-1}_m + \mathbf{b}_l \\
\implies (\mathbf{a}^{(l)})^\top \mathbf{A}^{-1}_l
  &= \lambda_{l, 1} \underbrace{(\mathbf{a}^{(l)})^\top \mathbf{A}^{-1}_{l, 1}}_{=0} + \dots + \underbrace{\lambda_{l, l-1} (\mathbf{a}^{(l)})^\top \mathbf{A}^{-1}_{l - 1}}_{=0} + \lambda_{l, l} (\mathbf{a}^{(l)})^\top \mathbf{a}^{(l)} + \lambda_{l, l + 1} \underbrace{(\mathbf{a}^{(l)})^\top \mathbf{A}^{-1}_{l + 1}}_{=0} + \dots + \lambda_{l,m} \underbrace{(\mathbf{a}^{(l)})^\top \mathbf{A}^{-1}_m}_{=0} \\
\implies \lambda_{l, l} 
  &= \frac{1}{\rvert \rvert \mathbf{a}^{(l)}\rvert\rvert_2^2} > 0
\end{aligned}
$$




















































































 -->
