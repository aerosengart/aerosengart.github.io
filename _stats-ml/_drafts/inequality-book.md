---
layout: distill
title: Constrained Statistical Inference Notes
description: 
date: 2026-07-13
tabs: true
tags: likelihood theory testing constraints
toc:
  - name: Background
  - name: A Couple of Illustrative Examples
    subsections:
      - name: Example 1
      - name: Example 2
  - name: The Chi-Bar-Square Distribution
  - name: General Tests in Parametric Settings
    subsections:
      - name: Set-Up
bibliography: stats-ml.bib
---

In this post, I will be going through some of the content in Silvapulle and Sen's <i>Constrained Statistical Inference: Inequality, Order, and Shape Restrictions</i>.<d-cite key=silvapulle2005></d-cite>.

---

## Background 


---

## A Couple Illustrative Examples
### Example 1
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

### Example 2
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

## The Chi-Bar-Square Distribution
We begin with the definition of the <i>chi-bar-square</i> distribution.

<div id="chi-bar-square"></div>
<div class="definition">
  <strong>Definition (Chi-Bar-Square Distribution).</strong>
  <br>
  Let $\mathcal{C} \subset \mathbb{R}^p$ be a closed, convex cone, and let $\mathbf{z} \sim \mathcal{N}(\mathbf{0}_p, \mathbf{V})$ for positive definite $\mathbf{V}$. 
  The <i>chi-bar-square distribution</i> is the distribution of the quantity:
  $$
  \mathbf{z}^\top \mathbf{V}^{-1} \mathbf{z} - \underset{\mathbf{c} \in \mathcal{C}}{\min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\}
  $$
  We let $\bar{\chi}^2(\mathbf{V}, \mathcal{C})$ denote such a random variable.
</div>

Since $\mathbf{V}$ is positive definite, so is $\mathbf{V}^{-1}$. We can then define an inner product $\langle \mathbf{a}, \mathbf{b} \rangle_{\mathbf{V}} = \mathbf{a}^\top \mathbf{V}^{-1} \mathbf{b}$ for any $\mathbf{a}, \mathbf{b} \in \mathbb{R}^p$, which induces the norm $\rvert \rvert \mathbf{d} \rvert \rvert_\mathbf{V} = \sqrt{\mathbf{d}^\top \mathbf{V}^{-1} \mathbf{d}}$. 

Let $\Pi_{\mathbf{V}}(\mathbf{z}, \mathcal{C})$ denote the orthogonal projection of $\mathbf{z}$ onto $\mathcal{C}$ in terms fo $\langle \cdot, \cdot \rangle_{\mathbf{V}}$. We call this the <strong>$\mathbf{V}$-projection</strong> of $\mathbf{z}$ on $\mathcal{C}$. Let $\tilde{\mathbf{z}}$ denote this projection. It is clear that:

$$
\tilde{\mathbf{z}} = \underset{\mathbf{c} \in \mathcal{C}}{\arg \min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\}
$$

By the orthogonal decomposition of $\mathbf{z}$ and the Pythagorean Theorem, we see that:

$$
\begin{equation}
\label{eq:cbs-dist-1}
\begin{aligned}
\rvert \rvert \mathbf{z} \rvert \rvert_\mathbf{V}
  &= \rvert \rvert \mathbf{z} - \tilde{\mathbf{z}} \rvert \rvert_{\mathbf{V}} + \rvert \rvert \tilde{\mathbf{z}} \rvert \rvert_{\mathbf{V}} \\
\implies
\mathbf{z}^\top \mathbf{V}^{-1} \mathbf{z} &= (\mathbf{z} - \tilde{\mathbf{z}})^\top \mathbf{V}^{-1} (\mathbf{z} - \tilde{\mathbf{z}}) + \underset{\mathbf{c} \in \mathcal{C}}{\min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\} \\
\implies
  \bar{\chi}^2 &= \rvert \rvert \tilde{\mathbf{z}} \rvert \rvert_\mathbf{V}
\end{aligned}
\end{equation}
$$

In words, the chi-bar-square random variable is the $\mathbf{V}$-length of the $\mathbf{V}$-projection of $\mathbf{z}$ onto $\mathcal{C}$. 

Let $$\mathcal{C}_0$$ be the polar cone to $\mathcal{C}$ with respect to $$\langle \cdot, \cdot \rangle_{\mathbf{V}}$$. Since $\mathcal{C}$ is closed and convex, so is $\mathcal{C}_0$. Now, let $$\tilde{\mathbf{z}}_0$$ denote the orthogonal projection of $\mathbf{z}$ onto $$\mathcal{C}_0$$ via $$\langle \cdot, \cdot \rangle_\mathbf{V}$$. By the same argument as above, we have:

$$
\begin{equation}
\label{eq:cbs-dist-2}
\begin{aligned}
\rvert \rvert \mathbf{z} \rvert \rvert_\mathbf{V}
  &= \rvert \rvert \mathbf{z} - \tilde{\mathbf{z}}_0 \rvert \rvert_{\mathbf{V}} + \rvert \rvert \tilde{\mathbf{z}}_0 \rvert \rvert_{\mathbf{V}} \\
\implies
\mathbf{z}^\top \mathbf{V}^{-1} \mathbf{z} &= (\mathbf{z} - \tilde{\mathbf{z}})_0^\top \mathbf{V}^{-1} (\mathbf{z} - \tilde{\mathbf{z}}_0) + \underset{\mathbf{c} \in \mathcal{C}_0}{\min} \left\{ (\mathbf{z} - \mathbf{c})^\top \mathbf{V}^{-1} (\mathbf{z} - \mathbf{c}) \right\} \\
\implies
  \bar{\chi}^2 &= \rvert \rvert \tilde{\mathbf{z}}_0 \rvert \rvert_\mathbf{V}
\end{aligned}
\end{equation}
$$

Note that, by the definition of the polar cone, we know that $$\langle \tilde{\mathbf{z}}_0, \mathbf{z} \rangle_\mathbf{V} = 0$$, $$\langle \mathbf{z} - \tilde{\mathbf{z}} , \tilde{\mathbf{z}} \rangle_{\mathbf{V}} = 0$$, and $$\langle \mathbf{z} - \tilde{\mathbf{z}}_0 , \tilde{\mathbf{z}}_0 \rangle_{\mathbf{V}} = 0$$. These vectors then define three sides of a rectangle, so we must have that $$\langle \mathbf{z} - \tilde{\mathbf{z}}_0, \mathbf{z} - \tilde{\mathbf{z}}_0 \rangle = 0$$. In addition, it must be the case that $$\rvert \rvert \tilde{\mathbf{z}} \rvert \rvert_\mathbf{V} = \rvert \rvert \mathbf{z} - \tilde{\mathbf{z}}_0 \rvert \rvert_{\mathbf{V}}$$ and $$\rvert \rvert \tilde{\mathbf{z}}_0 \rvert \rvert_\mathbf{V} = \rvert \rvert \mathbf{z} - \tilde{\mathbf{z}} \rvert \rvert_{\mathbf{V}}$$.

<aside><p>See Proposition 3.4.1 in Silvapulle and Sen.<d-cite key=silvapulle2005></d-cite></p></aside>

The tail probabilities of a $\bar{\chi}^2$ random variable are a mixture of the tail probabilities of $\chi^2$ distributions. 

<div id="chi-bar-square-distribution"></div>
<div class="theorem">
  <strong>Theorem 3.4.2.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  Let $\mathcal{C} \in \mathbb{R}^p$ be a closed, convex cone, and let $\mathbf{V}$ be a $p \times p$ positive definite matrix. Then for $c > 0$:
  $$
  \mathbb{P}(\bar{\chi}^2(\mathbf{V}, \mathcal{C}) \leq c \rvert H_0) = \sum_{i = 0}^p \omega_i(p, \mathbf{V}, \mathcal{C}) \mathbb{P}(\chi_i^2 \leq c)
  $$
  where $\omega_i(p, \mathbf{V}, \mathcal{C}) \geq 0$ for $i = 0, \dots, p$ and $\sum_{i = 0}^p \omega_i(p, \mathbf{V}, \mathcal{C}) = 1$.
</div>

These weights can sometimes be found analytically or via simulation. In some cases they can be very difficult to find or estimate. A similar result holds for replacing the null parameter space with a linear space contained in $\mathcal{C}

<div id="test-lin-space"></div>
<div class="theorem">
  <strong>Theorem 3.7.1.<d-cite key=silvapulle2005></d-cite></strong>
  <br>
  Let $\mathcal{M} \subset \mathcal{C}$ be a linear space where $\mathcal{C} \in \mathbb{R}^p$ is a closed, convex cone, and let $\mathbf{V}$ be a $p \times p$ positive definite matrix. Let $\mathbf{x} \sim \mathcal{N}(\boldsymbol{\mu}, \mathbf{V})$, and consider testing:
  $$
  H_0: \boldsymbol{\mu} \in \mathcal{M}
  \hspace{5mm} vs. \hspace{5mm}
  H_1: \boldsymbol{\mu} \in \mathcal{C} \setminus \mathcal{M}
  $$
  Then for $c > 0$:
  $$
  \mathbb{P}(\text{LRT} \leq c \rvert H_0)  = \sum_{i = 0}^p \omega_i(p, \mathbf{V}, \mathcal{C} \cap \mathcal{M}^\perp) \mathbb{P}(\chi_i^2 \leq c)
  $$
  where $\mathcal{M}^\perp$ is the orthogonal complement of $\mathcal{M}$ with respect to $\langle \cdot, \cdot \rangle_\mathbf{V}$. 
</div>

---

## General Tests in Parametric Settings
In this section, we will look at parts of Chapter 4, which covers tests in general parametric models. The score test for variance components falls under this umbrella. 

### Set-Up
We will assume to have independent and identically distributed $$\mathbf{x}_1, \dots, \mathbf{x}_n$$ from distribution $$P_{\boldsymbol{\theta}}$$ with density function $$f(\mathbf{x}; \boldsymbol{\theta})$$ for $$\boldsymbol{\theta} \in \boldsymbol{\Theta} \subset \mathbb{R}^p$$. We will let $\boldsymbol{\theta}_0$ denote the true value of $\boldsymbol{\theta}$, and $\hat{\boldsymbol{\theta}}$ the global maximum likelihood estimate over $\boldsymbol{\Theta}$. When relevant, $\tilde{\boldsymbol{\theta}}$ and $\bar{\boldsymbol{\theta}}$ will denote the estimators over $\boldsymbol{\Theta}$ subject to the constraints imposed by $H_0$ and $H_1$, respectively. We will also let $$\hat{\boldsymbol{\theta}}_\mathbf{A}$$ denote the MLE over some subset $$\mathbf{A} \subset \boldsymbol{\Theta}$$, which satisfies:

<aside><p>Though we write $\mathbf{x}$ in bold, it can be scalar-valued.</p></aside>

$$
\ell(\hat{\boldsymbol{\theta}}_\mathbf{A}; \mathbf{x}) = \underset{\boldsymbol{\theta} \in \mathbf{A}}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{x}) \right\} + o_p(1)
$$

In what follows, we assume that $\hat{\boldsymbol{\theta}}_\mathbf{A}$ is consistent when $\boldsymbol{\theta}_0 \in \mathbf{A}$. We also assume that the following regularity conditions hold:

<ul>
<li>If $P_{\boldsymbol{\theta}_1} = P_{\boldsymbol{\theta}_2}$ for $\boldsymbol{\theta}_1, \boldsymbol{\theta}_2 \in \boldsymbol{\Theta}$, then $\boldsymbol{\theta}_1 = \boldsymbol{\theta}_2$.</li>
<li>$\frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}$, $\frac{\partial^2 \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}$, and $\frac{\partial^3 \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_k}$ exist almost everywhere for all $i, j, k = 1, \dots, p$.</li>
<li>These exists $G(\mathbf{x})$ such that $\int G(\mathbf{x}) d\mathbf{x} < \infty$ and the absolute values of the partial derivatives of $\log f(\mathbf{x}; \boldsymbol{\theta})$ with respect to $\boldsymbol{\theta}$ up to order three are bounded by $G(\mathbf{x})$ in a neighborhood of the true parameter value, $\boldsymbol{\theta}_0$.</li>
<li>$\mathcal{I}(\boldsymbol{\theta})$ has finite entries and is positive definite.</li>
</ul>

We write the log-likelihood, score function, and information matrix for one observation as:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{x}) &= \log f(\mathbf{x}; \boldsymbol{\theta}) \\
\mathbf{S}(\boldsymbol{\theta}) &= \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x})}{\partial \boldsymbol{\theta}} \\
\mathcal{I}_{\boldsymbol{\theta}} &= \mathcal{I}(\boldsymbol{\theta}) = \mathbb{E}_{\boldsymbol{\theta}} \left[ \frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}^\top}\right]
\end{aligned}
$$

### Simple Null Score Test
Consider testing:

$$
H_0: \boldsymbol{\theta} = \boldsymbol{\theta}_0
\hspace{5mm}
\text{vs.}
\hspace{5mm}
H_n: \boldsymbol{\theta} \neq \boldsymbol{\theta}_0
$$

The score test statistic is given by:

$$
S_L = \frac{1}{n} \mathbf{S}^\top(\boldsymbol{\theta}_0) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) \mathbf{S}(\boldsymbol{\theta}_0)
$$

Under $H_0$, $$S_L \rightsquigarrow \chi_p^2$$ because $$\frac{1}{n} \mathbf{S}(\boldsymbol{\theta}_0) \overset{p}{\rightarrow} \mathbf{0}_p$$ and $$\frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_0) \rightsquigarrow \mathcal{N}(\mathbf{0}_p, \mathcal{I}(\boldsymbol{\theta}_0))$$. 

### Adding Nuisance Parameters
Suppose we can partition the parameter vector as $$\boldsymbol{\theta} = (\boldsymbol{\lambda}^\top, \boldsymbol{\psi}^\top)^\top = (\boldsymbol{\lambda} : \boldsymbol{\psi})$$ where $\boldsymbol{\psi} \in \mathbb{R}^r$. We will consider testing:

$$
H_0: \boldsymbol{\psi} = \boldsymbol{\psi}_0
\hspace{5mm}
\text{vs.}
\hspace{5mm}
H_1: \boldsymbol{\psi}  \neq \boldsymbol{\psi}_0
$$

In this case, we will partition $$\mathbf{S}(\boldsymbol{\theta})$$ as $$(\mathbf{S}_{\boldsymbol{\lambda}} : \mathbf{S}_{\boldsymbol{\psi}})$$, $$\mathcal{I}(\boldsymbol{\theta})$$ as $$[ \mathcal{I}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}}, \mathcal{I}_{\boldsymbol{\lambda}, \boldsymbol{\psi}} | \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\lambda}}, \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\psi}}]$$, and $\mathcal{I}^{-1}$ as $$\mathcal{I}(\boldsymbol{\theta})$$ as $$[ \mathcal{I}^{\boldsymbol{\lambda}, \boldsymbol{\lambda}}, \mathcal{I}^{\boldsymbol{\lambda}, \boldsymbol{\psi}} | \mathcal{I}^{\boldsymbol{\psi}, \boldsymbol{\lambda}}, \mathcal{I}^{\boldsymbol{\psi}, \boldsymbol{\psi}}]$$ to conform with the partitioning of $\boldsymbol{\theta}$. 

Using some facts about linear algebra, we have that:

$$
\mathcal{I}^{\boldsymbol{\psi}, \boldsymbol{\psi}} = \left[ \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\psi}} - \mathcal{I}_{\boldsymbol{\psi}, \boldsymbol{\lambda}} \mathcal{I}^{-1}_{\boldsymbol{\lambda}, \boldsymbol{\lambda}} \mathcal{I}_{\boldsymbol{\lambda}, \boldsymbol{\psi}} \right]^{-1}
$$

Furthermore, if $H_0$ is true, then the true value of $\boldsymbol{\theta}$ will have the form $\boldsymbol{\theta}_0 = (\boldsymbol{\lambda}_0 : \boldsymbol{\psi}_0)$. We let $\bar{\boldsymbol{\lambda}}$ and $\bar{\boldsymbol{\theta}} = (\bar{\boldsymbol{\lambda}} : \boldsymbol{\psi}_0)$ denote the MLE of $\boldsymbol{\lambda}$ and $\boldsymbol{\theta}$, respectively, under $H_0$.

The score test statistic is given by:

$$
S_L = \frac{1}{n}\mathbf{S}^\top_{\boldsymbol{\psi}}(\bar{\boldsymbol{\theta}}) (\mathcal{I}^{\boldsymbol{\psi}, \boldsymbol{\psi}}(\bar{\boldsymbol{\theta}}))^{-1} \mathbf{S}_{\boldsymbol{\psi}}(\bar{\boldsymbol{\theta}})
$$

### Local Score Test
Consider testing:

$$
H_0: \boldsymbol{\theta} = \mathbf{0}_p
\hspace{5mm}
\text{vs.}
\hspace{5mm}
H_1: \boldsymbol{\theta} \in \mathcal{C}
$$

for some closed, convex cone $\mathcal{C}$. Assume that the regularity conditions above hold and that $\mathbf{0}_p$ is on the interior of $\boldsymbol{\Theta}$. We define:

$$
\mathbf{U} = \frac{1}{\sqrt{n}} \mathcal{I}^{-1}(\mathbf{0}_p) \mathbf{S}(\mathbf{0}_0)
$$

For a fixed point $\boldsymbol{\delta} \in \mathbb{R}^p$, we can let $\boldsymbol{\theta}_n = \frac{1}{\sqrt{n}} \boldsymbol{\delta}$ and define a sequence of alternatives $H_n: \boldsymbol{\theta} = \boldsymbol{\theta}_n$. It can be shown that, under $H_n$:

$$
\mathbf{U} \rightsquigarrow \mathcal{N}(\boldsymbol{\delta}, \mathcal{I}^{-1}(\mathbf{0}_p))
$$

<aside><p>How?</p></aside>

Suppose the true value of $\boldsymbol{\theta}$ is $\frac{1}{\sqrt{n}} \boldsymbol{\delta}$. Then $\boldsymbol{\theta} \in \mathcal{C}$ if $\boldsymbol{\delta} \in \mathcal{C}$, and $\boldsymbol{\theta} = \mathbf{0}_p$ if $\boldsymbol{\delta} = \mathbf{0}_p$. Thus, we can alternatively test $\boldsymbol{\delta} = \mathbf{0}_p$ against $\boldsymbol{\delta} \in \mathcal{C}$ instead of our original hypotheses.





---



$$
H_0: \boldsymbol{\theta} = \boldsymbol{\theta}_0
\hspace{5mm}
\text{vs.}
\hspace{5mm}
H_n: \boldsymbol{\theta} = \boldsymbol{\theta}_n := \boldsymbol{\theta}_0 + \frac{1}{\sqrt{n}} \boldsymbol{\delta}
$$

for fixed $\boldsymbol{\delta}$. Here, $H_n$ is a sequence of alternatives. We can expand $\mathbf{S}(\boldsymbol{\theta}_0)$ about $\boldsymbol{\theta}_n$ under $H_n$ as:

$$
\begin{aligned}
\frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_0)
  &= \frac{1}{\sqrt{n}} \mathbf{S}(\boldsymbol{\theta}_n) + \left[ \frac{\partial \mathbf{S}(\boldsymbol{\theta})}{\partial \boldsymbol{\theta}^\top}\right]_{\boldsymbol{\theta} = \boldsymbol{\theta}_n} (\boldsymbol{\theta}_n - \boldsymbol{\theta}_0) + o_p(1) \\
  &= \frac{1}{\sqrt{n}}\mathbf{S}(\boldsymbol{\theta}_n) + \mathcal{I}(\boldsymbol{\theta}_n) \boldsymbol{\delta} + o_p(1)
\end{aligned}
$$

<aside><p>Why?</p></aside>





We have:

$$
\begin{aligned}
\mathbb{E}_{\boldsymbol{\theta}} \left[ \log\left(\frac{f(\mathbf{x}; \boldsymbol{\theta}^*)}{f(\mathbf{x}; \boldsymbol{\theta})}\right)\right] 
  &\leq \mathbb{E}_{\boldsymbol{\theta}} \left[ \frac{f(\mathbf{x}; \boldsymbol{\theta}^*)}{f(\mathbf{x}; \boldsymbol{\theta})} \right] - 1 & \left( \log(x) \leq x - 1 \right) \\
  &= 0  & \left(\mathbb{E}_{\boldsymbol{\theta}} \left[ \frac{f(\mathbf{x}; \boldsymbol{\theta}^*)}{f(\mathbf{x}; \boldsymbol{\theta})} \right] = 1 \text{  because } \boldsymbol{\theta}^*$ \text{ is consistent?}\right)
\end{aligned}
$$
