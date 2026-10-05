---
layout: distill
title: Definitions
description: 
date: 1998-10-26
tabs: true
tags: 
toc: 
    - name: Geometry and Linear Algebra
bibliography: stats-ml.bib
---

Out of convenience, I've put some very useful definitions for my own reference here.

---

## Statistics

### Consistency
There are several, slightly different, definitions of <i>consistency</i> for an estimator.

#### Weak Consistency
An estimator $\hat{\boldsymbol{\theta}}_n$ of a parameter $\boldsymbol{\theta}$ is called <i>weakly consistent</i> if it converges in probability to $\boldsymbol{\theta}$. That is, for all $\epsilon > 0$:

$$
\underset{n \rightarrow \infty}{\lim} \left\{\rvert \rvert \hat{\boldsymbol{\theta}}_n - \boldsymbol{\theta} \rvert\rvert > \epsilon \right\} = 0 
$$

This can be written alternatively as $$\hat{\boldsymbol{\theta}}_n = \boldsymbol{\theta} + o_p(1)$$.

#### Strong Consistency
An estimator $\hat{\boldsymbol{\theta}}_n$ of a parameter $\boldsymbol{\theta}$ is called <i>strongly consistent</i> if it converges almost surely to $\boldsymbol{\theta}$. That is:

$$
\mathbb{P}\left(\underset{n \rightarrow \infty}{\lim} \left\{  \hat{\boldsymbol{\theta}}_n \right\} = \boldsymbol{\theta} \right) = 1
$$

#### Fisher Consistency
Suppose we have independent and identically distributed random variables $x_1, \dots, x_n$ from some distribution $P_{\boldsymbol{\theta}}$ parametrized by $\boldsymbol{\theta}$. Define the empirical distribution function as:

$$
\tilde{\mathcal{F}}_n(x) \frac{1}{n} \sum_{i = 1}^n \text{hv}(x - x_{(i)}); \hspace{5mm}
\text{hv}(x - x_{(i)}) = \begin{cases} 1 & x - x_{(i)} \geq 0 \\
0 & x - x_{(i)} < 0
\end{cases}
$$

Let $\hat{\boldsymbol{\theta}}_n$ be an estimator that is some function of the empirical distribution function. That is:

$$
\hat{\boldsymbol{\theta}}_n = T\left(\tilde{\mathcal{F}}_n(\cdot)\right)
$$

We call $\hat{\boldsymbol{\theta}}_n$ <i>Fisher consistent</i> if $T(\mathcal{F}(\cdot, \boldsymbol{\theta})) = \boldsymbol{\theta}$. A Fisher consistent estimator will also be weakly consistent.


---

## Geometry and Linear Algebra

### Cone 
Let $\mathbf{V}$ be a vector space. A subset $\mathbf{A} \subset \mathbf{V}$ is called a <strong>cone</strong> if:

$$
\lambda \mathbf{x} \in \mathbf{A};
\hspace{5mm}
\text{ for } \mathbf{x} \in \mathbf{A}
\text{ and } \lambda > 0
$$


### Convex Set
Let $\mathbf{V}$ be a vector space over the reals. A subset $\mathbf{A} \subset \mathbf{V}$ is called <strong>convex</strong> if:

<aside><p>This definition extends to vector spaces over ordered fields.</p></aside>

$$
\lambda \mathbf{x} + (1 - \lambda) \mathbf{y} \in \mathbf{A}; \hspace{5mm}
\text{ for } \mathbf{x}, \mathbf{y} \in \mathbf{A} 
\text{ and } \lambda \in [0, 1]
$$


### Field
A set $\mathbf{F}$ is called a <strong>field</strong> if it is equipped with the two binary relations (i.e. mappings of the form $\mathbf{F} \times \mathbf{F} \rightarrow \mathbf{F}$), <i>addition</i> ($+$) and <i>multiplication</i> ($\times$) that satisfy for any $a, b, c \in \mathbf{F}$:

<ol>
<li><strong>Associativity</strong>: $a + (b + c) = (a + b) + c$ and $a \cdot (b \cdot c) = (a \cdot b) \cdot c$</li>
<li><strong>Commutativity</strong>: $a + b = b + a$ and $a \cdot b = b \cdot a$</li>
<li><strong>Identity</strong>: $\exists 0, 1 \in \mathbf{F}$ such that $a + 0 = a$ and $a \cdot 1 = a$</li>
<li><strong>Inverses</strong>: for any $a \in F$, $\exists -a \in F$ such that $a + (-a) = 0$, and if $a \neq 0$, then $\exists a^{-1} \in F$ such that $a \cdot a^{-1} = 1$</li>
<li><strong>Distributivty</strong>:$a \cdot (b + c) = (a \cdot b) + (a \cdot c)$</li>
</ol>

### Hilbert Space
An inner product space $\mathbf{V}$ is also a <strong>Hilbert space</strong> if it is complete with respect to the distance function induced by the inner product. 

In other words, let $\mathbf{x}_1, \mathbf{x}_2, \dots$ be a sequence of elements in $\mathbf{V}$. Suppose for every $r > 0$, there exists $n \in \mathbb{N}$ such that for all $m, n > N$:

$$
d(\mathbf{x}_m, \mathbf{x}_n) < r
$$

If every such sequence in $\mathbf{V}$ converges to some $\mathbf{x} \in \mathbf{V}$, then $\mathbf{V}$ is <strong>complete</strong>.


### Inner Product Space
An <strong>inner product space</strong> is a vector space $\mathbf{V}$ over some field $\mathbf{F}$ endowed with a mapping (called an <i>inner product</i>) $\langle \cdot, \cdot \rangle : \mathbf{V} \times \mathbf{V} \rightarrow \mathbf{F}$ that satisfies for any $\mathbf{x}, \mathbf{y}, \mathbf{z} \in \mathbf{V}$ and $a, b \in \mathbf{F}$:

<ol>
<li><strong>Conjugate Symmetry</strong>: $\langle \mathbf{x}, \mathbf{y} \rangle = \overbar{\langle \mathbf{y}, \mathbf{x} \rangle}$ (if $\mathbf{F} = \mathbb{R}$, then this is just regular symmetry)</li>
<li><strong>Linearity</strong>: $\langle a \mathbf{x} + b \mathbf{y}, \mathbf{z} \rangle = a \langle \mathbf{x}, \mathbf{z} \rangle + b \langle \mathbf{y}, \mathbf{z} \rangle$</li>
<li>Positive-Definiteness</strong>: $\langle \mathbf{x}, \mathbf{x} \rangle = 0$ if $\mathbf{x} \neq \mathbf{0}$</li>
</ol>

### Orthogonal Complement
For a vector subspace $\mathbf{W}$ of an inner product space $\mathbf{H}$, the <strong>orthogonal complement</strong> of $\mathbf{W}$ is the vector subspace of vectors in $\mathbf{H}$ that are orthogonal to all vectors in $\mathbf{W}$:

$$
\mathbf{W}^\perp = \left\{ \mathbf{x} \in \mathbf{H} : \langle \mathbf{x}, \mathbf{w} \rangle = 0 \hspace{2mm} \forall \mathbf{w} \in \mathbf{W} \right\}
$$

### Orthogonal Decomposition
Let $\mathbf{W} \subset \mathbb{R}^p$ and let $\mathbf{y} \in \mathbb{R}^p$. The <strong>orthogonal decomposition</strong> of $\mathbf{y}$ is the unique sum:

$$
\mathbf{y} = \mathbf{y}_{\mathbf{W}} + \mathbf{y}_{\mathbf{W}^\perp}
$$

where $\mathbf{y}_{\mathbf{W}} \in \mathbf{W}$ and $\mathbf{y}_{\mathbf{W}^\perp} \in \mathbf{W}^\perp$, the orthogonal complement of $\mathbf{W}$.

### Polar Cone
Let $\mathcal{C}$ denote a cone in vector space $\mathbf{V}$ with inner product denoted by $\langle \cdot, \cdot \rangle$. The <strong>polar cone</strong> of $\mathcal{C}$ is the set:

$$
\mathcal{C}^0 = \left\{ \mathbf{c} \in \mathbf{V} : \langle \mathbf{y}, \mathbf{x} \rangle \leq 0 \hspace{3mm} \forall \mathbf{x} \in \mathcal{C} \right\}
$$

Geometrically speaking, the polar cone of $\mathcal{C}$ is the set of vectors that form non-acute angles with all vectors in $\mathcal{C}$. Thus, vectors along the boundary of the polar cone will be orthogonal to the those on the boundary of $\mathcal{C}$.

### Projection
Let $\mathbf{V}$ be a vector space. A <strong>projection</strong> is a linear operator $\mathbf{P}: \mathbf{V} \rightarrow \mathbf{V}$ satisfying $\mathbf{P}^2 = \mathbf{P}$. A square matrix $\mathbf{P}$ is called a <strong>projection matrix</strong> if $\mathbf{P}^2 = \mathbf{P}$. It will have eigenvalues equal to $0$ or $1$. 

If $\mathbf{V}$ is a Hilbert space, then $\mathbf{P}$ is an <strong>orthogonal projection</strong> if:

$$
\langle \mathbf{P}\mathbf{x}, \mathbf{y} \rangle = \langle \mathbf{x}, \mathbf{P}\mathbf{y} \rangle \hspace{5mm} \forall \mathbf{x}, \mathbf{y} \in \mathbf{V}
$$

We denote the <strong>orthogonal projection of $\mathbf{x}$ onto $\mathbf{V}$</strong> with respect to the inner product induced by positive definite matrix $\mathbf{M}$ (i.e. $\langle \mathbf{x}, \mathbf{y} \angle_{\mathbf{M}} = \mathbf{x} \mathbf{M}^{-1} \matbf{y}$) as:

$$
\Pi_{\mathbf{M}}(\mathbf{x} \rvert \mathbf{V}) = \underset{\mathbf{v} \mathbf{V}}{\arg \min} \left\{ (\mathbf{x} - \mathbf{v})^\top \mathbf{M}^{-1}(\mathbf{x} - \mathbf{v}) \right\}
$$



### Vector Space
A non-empty set $\mathbf{V}$ is called a <strong>vector space</strong> over a field $\mathbf{F}$ if it is equipped with the binary operation <strong>vector addition</strong> and the binary function <strong>scalar multiplication</strong> which satisfy for any $\mathbf{u}, \mathbf{v}, \mathbf{w} \in \mathbf{V}$ and $a, b \in \mathbf{F}$:

<ol>
<li><strong>Associativity</strong>: $\mathbf{u} + (\mathbf{v} + \mathbf{w}) = (\mathbf{u} + \mathbf{v}) + \mathbf{w}$</li>
<li><strong>Commutativity</strong>: $\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$</li>
<li><strong>Identity</strong>: $\exists \mathbf{0} \in \mathbf{V}$ such that $\mathbf{v} + \mathbf{0} = \mathbf{v}$ for all $\mathbf{v} \in \mathbf{V}$ and $1 \mathbf{v} = \mathbf{v}$ for $1 \in \mathbf{F}$</li>
<li><strong>Inverse</strong>: for every $\mathbf{v} \in \mathbf{V}$, $\exists -\mathbf{v} \in \mathbf{V}$ such that $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$</li>
<li>$a(b \mathbf{v}) = (a b) \mathbf{v}$</li>
<li><strong>Distributivity</strong>: $a(\mathbf{u} + \mathbf{v}) = a \mathbf{u} + a \mathbf{v}$ and $(a + b) \mathbf{v} = a \mathbf{v} + b \mathbf{v}$</li>
</ol>