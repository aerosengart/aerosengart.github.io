---
layout: distill
title:  Useful Definitions
description:
date: 2026-01-16
tabs: true
tags: reference
toc:
  - name: The Very Basics
    subsections:
            - name: Sets
            - name: Metrics
            - name: Sequences
            - name: Functions
    
bibliography: measure.bib
---

For easy reference, I've decided to put some of the more "building block" definitions in a separate post. Most of the definitions are from Wikipedia and Axler. In addition, most concepts are introduced in the context of metric spaces with further discussion in Euclidean space. 

---

## Sets

### Types

We'll begin with the <i>field</i>, which is just a set with some binary operations. Binary operations are basically just functions taking in two elements from a set, $S$, and mapping them to some (unique) element also in $S$ (i.e. $f: S \times S \rightarrow S$). 

<div id="field"></div>
<div class="definition">
  <strong>Definition (Field).<d-cite key=field2026></d-cite></strong>
  <br>
  A <i>field</i> is a set, $F$, and the two binary operations of <i>addition</i> (which we denote with $+$) and <i>multiplication</i> (which we denote with $\cdot$) that satisfy the following axioms. For any $a, b, c$ in $F$:

  <ul>
    <li><strong>Associativity</strong>: $(a + (b + c)) = (a + b) + c$ and $a \cdot (b \cdot c) = (a \cdot b) \cdot c$</li>
    <li><strong>Commutativity</strong>: $a+b = b+a$ and $a \cdot b = b \cdot a$</li>
    <li><strong>Identity</strong>: $0$ and $1$ are elements in $F$ such that $a + 0 = a$ and $a \cdot 1 = a$</li>
    <li><strong>Distributivity</strong>: $a \cdot (b + c) = (a \cdot b) + (a \cdot c)$</li>
    <li><strong>Additive Inverse</strong>: For any $a$ in $F$, there exists its <i>additive inverse</i> in $F$, denoted by $-a$, which satisfies $a + (-a) = 0$</li>
    <li><strong>Multiplicative Inverse</strong>: For any $a \neq 0$ in $F$, there exists its <i>multiplicative inverse</i> in $F$, denote by $a^{-1}$, which satisfies $a \cdot a^{-1} = 1$</li>
  </ul>
</div>

A set, $M$, with a metric, $d$, together as an ordered pair, $(M, d)$, forms a <strong>metric space</strong>. 


### Characteristics

Sets (i.e. spaces for our purposes) can be described in a lot of different ways. First, let's look at what it means to be <i>complete</i>.

<div id="complete-space"></div>
<div class="definition">
  <strong>Definition (Complete).<d-cite key=complete2025></d-cite></strong>
  <br>
  A metric space, $(M, d)$, is called <i>complete</i> if every <a href="#cauchy">Cauchy sequence</a> of points in $M$ has its limit in $M$. 

  The metric, $d$, is also called complete. 
</div>

This implies a <strong>complete field</strong> is a field equipped with complete metric.

We can also describe a set as <strong>open</strong> or <strong>closed</strong>.

<div id="open-closed"></div>
<div class="definition">
  <strong>Definition (Open/Closed).</strong><d-cite key=open2025></d-cite>
  <br>
  Let $U$ be a subset of metric space $(M, d)$. We call $U$ <i>open</i> if, for any $x \in U$, there exists a real number $\epsilon > 0$ such that for any $y \in M$ where $d(x, y) < \epsilon$, we have that $y \in M$. The complement of an open set is called <i>closed</i>.
  <br>
  Equivalently, $U$ is called open if for every $u \in U$, there exists some $\delta > 0$ such that the <i>open ball</i> centered at $u$ with radius $\delta$ is contained in $U$. That is:
  $$
  B(x, \delta) = \{ y \in M: d(y - x) < \delta \} \subseteq U
  $$
</div>

<aside><p>One way to think about what it means to be open is that there is some distance, $\epsilon$, associated with each element, $x$, in the set, and all of the points that are less than $\epsilon$ away from $x$ are also in the set.</p></aside>

Though a bit confusing, sets can also be both closed and open (also called <strong>clopen</strong>). The empty set is an example.

Open sets have some nice properties including the fact that the (finite or infinite) union of open sets is open, and the finite intersection of open sets is open. In contrast, the (finite or infinite) intersection of closed sets is closed, and the finite union of closed sets is closed.

### Ordering 

Elements can be compared using <i>orders</i>, which are binary relations denoted by $\leq$. 

<div id="partial-order"></div>
<div class="definition">
  <strong>Definition (Partial Order).</strong><d-cite key=partial_order2025></d-cite>
  <br>
  A <i>partial order</i> is a binary relation, $\leq$, between a set, $X$, and itself satisfying the following for any $a, b, c \in X$:

  <ol>
    <li> Reflexivity: $a \leq a$ </li>
    <li> Antisymmetry: $a \leq b$ and $b \leq a \implies a = b$ </li>
    <li> Transitivity: $a \leq b$ and $b \leq c \implies a \leq c$ </li>
  </ol>

  Partial orders as defined above are sometimes called <i>reflexive</i>, <i>weak</i>, or <i>non-strict</i>. A <i>strict partial order</i> is a binary relation, $<$, between a set, $X$, and itself satisfying the following for all $a, b, c \in X$:

  <ol>
    <li> Irreflexivity: $\neg(a < a)$ </li>
    <li> Asymmetry: $a < b \implies \neg (b < a)$ </li>
    <li> Transitivity: $a < b$ and $b < c \implies a < c$ </li>
  </ol>
</div>

A <i>total order</i> is a partial order satisfying one additional property. 

<div id="total-order"></div>
<div class="definition">
  <strong>Definition (Total Order).</strong><d-cite key=total_order2025></d-cite>
  <br>
  A <i>total order</i>, also called a <i>linear order</i>, is a partial order satisfying the additional property for all $a, b, c \in X$:

  <ul>
    <li> Totality: $a \leq b$ or $b \leq a$ </li>
  </ul>

 
  Total orders as defined above are sometimes called <i>non-strict</i>. A <i>strict total order</i> is a strict partial order that satisfies the following additional proerpty for all $a, b \in X$:

  <ul>
    <li> Connectivity: $a \neq b \implies a < b$ or $b < a$ </li>
  </ul>
</div>

Joining a field, $F$, with a total order, $\leq$, defines an <strong>ordered field</strong> if the following are satisfied for all $a,b,c \in F$:

<ul>
    <li>If $a \leq b$, then $a + c \leq b + c$</li>
    <li>If $0 \leq a$ and $0 \leq b$ then $0 \leq a \cdot b$</li>
</ul>



---

## Sequences

Sequences will also come up a lot, and it's important to have a good handle on how these can be characterized.

<div id="cauchy"></div>
<div class="definition">
  <strong>Definition (Cauchy).<d-cite key=complete_field2025></d-cite></strong>
  <br>
  Let $x_1, x_2, x_3, \dots$ denote a sequence defined on some metric space, $(X, d)$. The sequence is called <i>Cauchy</i> if, for any positive $r \in \mathbb{R}$, there exists a positive $N \in \mathbb{N}$ such that for all $m, n > N$:

  $$
    d(x_m, x_n) < r
  $$
</div>

<aside><p> Intuitively, a Cauchy sequence is one in which its elements eventually get arbitrarily "close" together (with closeness defined by the metric $d$). </p></aside>

Some sequences <i>converge</i>, which means they have a <i>limit</i>.

<div id="limit"></div>
<div class="definition">
  <strong>Definition (Limit).<d-cite key=limit2026></d-cite></strong>
  <br>
  Let $\{ x_n \}_{n \geq 0}$ be a sequence of elements in a metric space $(M, d)$, and let $L \in M$. If, for any $\epsilon > 0$, there exists $N$ such that for all $n > N$:

  $$
  d(L, x_n) < \epsilon
  $$

  then we call $L$ the <i>limit</i> of the sequence $\{ x_n \}_{n \geq 0}$ and denote it by:

  $$
  \underset{n \rightarrow \infty}{\lim} x_n = L
  $$
</div>

The above definition is for general metric spaces. Often we are dealing with $\mathbb{R}^n$, in which case we will choose either the $\ell^2$ norm or the infinity norm. Limits can also be defined coordinate-wise for vector-valued elements. 




---


## Norms

First, the <strong>$\ell^2$ norm</strong>, which is one of the most commonly encountered norms (at least for me). It is also called the <i>Euclidean norm</i>.

<div id="2-norm"></div>
<div class="definition">
  <strong>Definition ($\rvert \rvert \cdot \rvert \rvert_2$).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $x = (x_1, \dots, x_n) \in \mathbb{R}^n$. The <i>$\ell^2$ norm</i> is defined as:

  $$
  \rvert \rvert x \rvert \rvert_2 = \sqrt{ \rvert x_1 \rvert ^2 + \dots + \rvert x_n \rvert^2}
  $$
</div>


Now, the <strong>infinity norm</strong>.

<div id="infinity-norm"></div>
<div class="definition">
  <strong>Definition ($\rvert \rvert \cdot \rvert \rvert_{\infty}$).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $x = (x_1, \dots, x_n) \in \mathbb{R}^n$. The <i>infinity norm</i> is defined as:

  $$
  \rvert \rvert x \rvert \rvert_{\infty} = \max \left\{ \rvert x_1 \rvert,   \dots, \rvert x_n \rvert \right\}
  $$
</div>

---


## Functions
First thing's first: what's a function?

<div id="function"></div>
<div class="definition">
  <strong>Definition (Function).<d-cite key=function2025></d-cite></strong>
  <br>
  For sets $X$ and $Y$, a function, $f: X \rightarrow Y$, assigns a value in $Y$ to each value in $X$. 
</div>

The set $X$ is called its <strong>domain</strong>, and the set $Y$ is called its <strong>codomain</strong>. The <strong>image</strong> of a function (also called its <strong>range</strong>) is the subset of elements in $Y$ that are mapped to by elements in $x$. In notation:

$$
\text{Im}(f) = \{ f(x): x \in X \}
$$


### Types
A <i>distance metric</i> or <i>function</i> (or just <i>metric</i> for short) assigns a value to represent how "far apart" two elements are in a set.

<div id="distance-metric"></div>
<div class="definition">
  <strong>Definition (Distance Metric).<d-cite key=metric2025></d-cite></strong>
  <br>
  Let $X$ be a set. A <i>distance metric</i> is any function $d: X \times X \rightarrow \mathbb{R}$ satisfying the following for all $x, y, z, \in X$:

  <ul>
    <li>$d(x, x) = 0$</li>
    <li><strong>Positivity</strong>: $x \neq y \implies d(x,y) > 0$ </li>
    <li><strong>Symmetry</strong>: $d(x,y) = d(y,x)$</li>
    <li><strong>Triangle Inequality</strong>: $d(x,z) \leq d(x,y) + d(y, z)$</li>
  </ul>
</div>


### Characteristics
One of the most important characteristics of functions is <strong>continuity</strong>.

<div id="continuity"></div>
<div class="definition">
  <strong>Definition (Continuity).<d-cite key=continuous2025></d-cite></strong>
  <br>
  Let $f: X \rightarrow Y$ be a function between metric spaces $(X, d_X)$ and $(Y, d_Y)$. The function $f$ is called <i>continuous</i> at a point $c \in X$ if, for any positive real $\epsilon > 0$, there exists a positive real $\delta > 0$ such that $d_Y(f(x), f(c)) < \epsilon$ for all $x \in X$ such that $d_X(x, c) < \delta$. 
</div>

<aside><p>The basic idea is that we can find a value, $\delta$, that is governed by $\epsilon$, such that all points that are $\delta$-close to $c$ yield a change in the output (from $f(c)$) no greater than $\epsilon$. That is, the function does not change too much too fast.</p></aside>

If $f$ is a real-valued function, then the definition is equivalent to saying that $\underset{n \rightarrow \infty}{\lim} f(x_n) = f(c)$ for every $\{ x_n \}_{n \geq 0}$ in $X$ with $\underset{n \rightarrow \infty}{\lim} x_n = c$.

A stronger kind of continuity is <strong>uniform continuity</strong>.

<div id="unif-continuity"></div>
<div class="definition">
  <strong>Definition (Uniform Continuity).<d-cite key=ucontinuous2025></d-cite></strong>
  <br>
  Let $f: X \rightarrow Y$ be a function between metric spaces $(X, d_X)$ and $(Y, d_Y)$. The function $f$ is called <i>uniformly continuous</i> if, for every real $\epsilon > 0$, there exists a real $\delta > 0$ such that $d_x(f(x), f(y)) < \epsilon$ for all $x, y \in X$ such that $d_1(x, y) < \delta$. 
</div>

<aside><p>Uniform continuity is basically continuity but with a <i>single</i> $\delta$ for the entire domain of $f$. This means that <i>any</i> points that are close enough together will yield a difference in outputs that is smaller than the specified $\epsilon$.</p></aside>

Continuous, real-valued functions defined on a closed, bounded subset of the reals will be uniformly continuous. In addition, they will achieve their maxima and minima.