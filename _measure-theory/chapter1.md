---
layout: distill
title: Chapter 1
description: Riemann Integrals
date: 2026-01-17
tabs: true
tags: riemann
toc:
  - name: Riemann Sums
  - name: Riemann Integrals
  - name: Limitations
bibliography: measure.bib
---

In this post, I'll be going through Chapter 1 in Axler<d-cite key=axler2025></d-cite>. It covers Riemann integration.

In what follows, we'll let $\mathbb{R}$ denote the <a href="../../measure/useful-definitions/#open-closed">closed</a>, <a href="../../measure/useful-definitions/#total-order">ordered field</a> of real numbers.

---

## Riemann Sums

First, let's recall the building block of the Riemann integral: the <i>Riemann sum</i>.

<div id="riemann-sum"></div>
<div class="definition">
  <strong>Definition (Lower/Upper Riemann Sum).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $f: [a, b] \rightarrow \mathbb{R}$ be a bounded function, and let $P = \{ x_0, x_1, \dots, x_n \}$ denote a partition of $[a, b]$. The <i>lower Riemann sum</i> and <i>upper Riemann sum</i>, denoted respectively with $L(f, P, [a, b])$ and $U(f, P, [a, b])$, are given by:

  $$
  \begin{equation}
  \label{eq:lu-riemann-sum}
  \begin{aligned}
    L(f, P, [a, b]) &= \sum_{j = 1}^n (x_j - x_{j-1}) \underset{x \in [x_{j-1}, x_j]}{\inf} \left\{ f(x) \right\} \\
    U(f, P, [a, b]) &= \sum_{j = 1}^n (x_j - x_{j-1}) \underset{x \in [x_{j-1}, x_j]}{\sup} \left\{ f(x) \right\} \\
  \end{aligned}
  \end{equation}
  $$
</div>

<aside><p>Put into words, a Riemann sum of a function over some interval with respect to a given partition is the sum of areas of the rectangles (defined by our partition) that we can draw under the function. For a lower Riemann sum, we restrict the rectangles' heights to be the lowest value of the function on each subinterval of the partition. For an upper Riemann sum, we use the greatest value.</p></aside>


### Properties

One property of Riemann sums is that finer partitions (ones will more points) will yield greater lower Riemann sums and smaller upper Riemann sums. That is, if we let $P$ and $P'$ be partitions of $[a, b]$ such that $x \in P'$ for all $x \in P$, then the following holds:

$$
L(f, P, [a, b]) \leq L(f, P', [a, b]) \leq U(f, P', [a,b]) \leq U(f, P, [a, b])
$$

Another property is that lower Riemann sums are no greater than an upper Riemann sum (for the same function and interval). That is, for partitions $P$ and $S$ on $[a, b]$:

$$
L(f, P, [a, b]) \leq U(f, P', [a, b])
$$

Proofs of these properties are fairly simple and can be found in Axler. 

---

## Riemann Integrals
We now come to two stars of the chapter: the lower and upper Riemann integrals. 

<div id="riemann-integrals"></div>
<div class="definition">
  <strong>Definition (Lower/Upper Riemann Integral).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $f: [a, b] \rightarrow \mathbb{R}$ be a bounded function. The <i>lower Riemann integral</i> and <i>upper Riemann integral</i>, denoted respectively by $L(f, [a, b])$ and $U(f, [a, b])$, are given by:
  
  $$
  \begin{equation}
  \label{eq:lu-riemann-integral}
  \begin{aligned}
  L(f, [a, b]) &= \underset{P}{\sup} \left\{ L(f, P, [a, b]) \right\} \\
  U(f, [a, b]) &= \underset{P}{\sup} \left\{ U(f, P, [a, b]) \right\}
  \end{aligned}
  \end{equation}
  $$
</div>

We have a similar relationship between the lower and super Riemann integrals as Riemann sums:

$$
\begin{equation}
\label{eq:integral-prop-1}
L(f, [a, b]) \leq U(f, [a, b])
\end{equation}
$$

Now, the main event: the Riemann integral. First, we note that a bounded function $f: [a, b] \rightarrow \mathbb{R}$ is called <strong>Riemann integral</strong> if $L(f, [a, b]) = U(f, [a, b])$.

<div id="riemann-integral"></div>
<div class="definition">
  <strong>Definition (Riemann Integral).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $f: [a, b] \rightarrow \mathbb{R}$ be Riemann integrable. Its <i>Riemann integral</i> is defined by:

  $$
  \begin{equation}
  \label{eq:riemann-integral}
  \int_a^b f(x) dx = \int_a^b f = L(f, [a, b]) = U(f, [a, b])
  \end{equation}
  $$
</div>

<aside><p>Technically, these are <a href="https://en.wikipedia.org/wiki/Darboux_integral">Darboux integrals</a>, but they are equivalent to Riemann integrals.</p></aside>


### Properties
Below we have some important results about Riemann integrals.

<!-- #region axler-1-11 -->
<div class="theorem">
<strong>Claim.<d-cite key="axler2025"></d-cite></strong>
{% tabs axler-1-11 %}
{% tab axler-1-11 statement %}
All continuous real-valued functions on all closed, bounded intervals are Riemann integrable. 
{% endtab %}
{% tab axler-1-11 proof %}
Suppose that $a, b \in \mathbb{R}$ with $a < b$ and $f: [a, b] \rightarrow \mathbb{R}$ is continuous. Since $[a, b]$ is a closed interval on the reals, it is compact. Thus, by the <a href="https://en.wikipedia.org/wiki/Heine–Cantor_theorem">Heine-Cantor Theorem</a>, $f$ is <a href="/measure-theory/useful-definitions#unif-continuity">uniformly continuous</a>.
<br>
Fix $\epsilon > 0$. By the uniform continuity of $f$ on the reals, we have that there exists $\delta > 0$ such that:

$$
\begin{equation}
\label{eq:assumption}
\rvert f(s) - f(t) \rvert < \epsilon \text{ for all } s,t \in [a, b] \text{ with } \rvert s - t \rvert < \delta
\end{equation}
$$

Choose some positive integer $n$ such that $\frac{b - a}{n} < \delta$, and construct an evenly spaced partition of $[a, b]$ as:

$$
P = \{ x_0, x_1, \dots, x_n \} \text{ such that } x_j - x_{j - 1} = \frac{b - a}{n}
$$

We then have that:

$$
\begin{aligned}
U(f, [a, b]) - L(f, [a, b]) &\leq U(f, P, [a, b]) - L(f, P, [a, b]) & \left(\text{Eq. } \eqref{eq:lu-riemann-integral}\right) \\
&= \sum_{j = 1}^n (x_j - x_{j-1}) \underset{x \in [x_{j - 1}, x_j]}{\sup} \left\{ f(x) \right\} - \sum_{j = 1}^n (x_j - x_{j-1}) \underset{x \in [x_{j - 1}, x_j]}{\inf} \left\{ f(x) \right\} & \left(\text{Eq. } \eqref{eq:lu-riemann-sum} \right) \\
&= \frac{b - a}{n} \sum_{j = 1}^n \left( \underset{x \in [x_{j - 1}, x_j]}{\sup} \left\{ f(x) \right\}  - \underset{x \in [x_{j - 1}, x_j]}{\inf} \left\{ f(x) \right\}  \right) \\
&\leq \frac{b - a}{n} \sum_{j = 1}^n \epsilon & \left(\text{Eq. } \eqref{eq:assumption}\right) \\
&= (b - a)\epsilon
\end{aligned}
$$

Note that Eq. \eqref{eq:integral-prop-1} states that:

$$
U(f, [a, b]) \geq  L(f, [a, b]) 
$$

However, we have just shown that, for all $\epsilon > 0$:

$$
\begin{aligned}
U(f, [a, b]) - L(f, [a, b]) &\leq (b - a) \epsilon  \\
U(f, [a, b]) &\leq (b - a) \epsilon + L(f, [a, b]) \\
\overset{\epsilon \rightarrow 0}{\implies} U(f, [a, b]) &\leq L(f, [a, b])
\end{aligned}
$$

Together, these two points imply:

$$
U(f, [a, b]) = L(f, [a, b])
$$

which implies $f$ is Riemann integrable by Eq. \eqref{eq:riemann-integral}. 
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->



<div class="theorem">
<strong>Claim.<d-cite key="axler2025"></d-cite></strong>
{% tabs axler-1-18 %}
{% tab axler-1-18 statement %}
Let $a, b, M \in \mathbb{R}$ with $a < b$. Let $f_1, f_2, \dots$ be a sequence of Riemann integrable functions on $[a, b]$ such that, for all positive integers $k$ and all $x \in [a, b]$:

$$
\rvert f_k(x) \rvert \leq M
$$

Assume:

$$
f(x) = \underset{k \rightarrow \infty}{\lim} f_k(x)
$$

exists for all $x \in [a, b]$. If $f$ is Riemann integrable on $[a, b]$ then:

$$
\int_a^b f(x) dx = \underset{k \rightarrow \infty}{\lim} \int_a^b f_k(x) dx
$$

{% endtab %}
{% endtabs %}
</div>

---

## Limitations
Unfortunately, for a lot of the things we want to do, the Riemann integral is insufficient. The first limitation that Axler introduces is that many functions that are "not so nice" are not Riemann integrable. 

<div class="example">
    <strong>Example.<d-cite key="axler2025"></d-cite></strong>
    Let $f: [0, 1] \rightarrow \mathbb{R}$ be the function defined by:
    $$
    f(x) = \begin{cases}
    1 & \text{if } x \text{ is rational} \\
    0 & \text{if } x \text{ is irrational}
    \end{cases}
    $$
    This function is not Riemann integrable because its lower and upper Riemann integrals are not equal. For any $[a, b] \subseteq [0, 1]$ with $a < b$:
    $$
    \underset{x \in [a,b]}{\inf} \left\{ f(x) \right\} = 0 \hspace{10mm}
    \underset{x \in [a,b]}{\inf} \left\{ f(x) \right\} = 1
    $$
    since all intervals contain at least one rational and irrational number. This implies:
    $$
    L(f, P, [0, 1]) = 0 \hspace{10mm} U(f, P, [0, 1]) = 1
    $$
    for all partitions of $[0, 1]$. Thus, $L(f, [a, b]) = 0$ and $U(f, [a, b]) = 1$.
</div>

The example above is a function with many discontinuities (infinitely many, in fact). Another problem that arises is related to the interchanging of the Riemann integral and the limit function.

<div class="example">
    <strong>Example.<d-cite key="axler2025"></d-cite></strong>
    Let $r_1, r_2, \dots$ be the sequence of all rational numbers in the interval $[0, 1]$ (once and in order). For positive integer $k$, let $f_k: [0, 1] \rightarrow \mathbb{R}$ be the function defined by:
    $$
    f_k(x) = \begin{cases}
    1 & \text{if } x \in \{ r_1, \dots, r_k \} \\
    0 & \text{else}
    \end{cases}
    $$
    Each $f_k$ is Riemann integrable, and we also have that:
    $$
    \underset{k \rightarrow \infty}{\lim} f_k(x) = f(x)
    $$
    However, as we showed in the previous example, $f(x)$ is not Riemann integrable.
</div>

Unbounded functions are also not Riemann integrable and reveal some deficiencies in Riemann integration.

<div class="example">
    <strong>Example.<d-cite key="axler2025"></d-cite></strong>
    Let $f: [0, 1] \rightarrow \mathbb{R}$ be the function defined by:
    $$
    f(x) = \begin{cases}
    \frac{1}{\sqrt{x}} & \text{if } 0 < x \leq 1 \\
    0 & \text{if } x = 0
    \end{cases}
    $$
    For any partition $P = \{ x_0, x_1, \dots, x_n \}$ of $[0, 1]$, $\underset{x \in [x_0, x_1]}{\sup}\left\{ f(x) \right\} = \infty$ and $\underset{x \in [x_0, x_1]}{\inf} \left\{ f(x) \right\} = 0$. These results imply the upper Riemann integral is infinity while the lower one is $0$. 
    <br>
    Another concern is the fact that the area under $f$ on $[0, 1]$ to be equal to $2$. 
</div>



