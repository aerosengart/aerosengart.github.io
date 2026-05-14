---
layout: distill
title: Chapter 2
description: Measures
date: 2026-05-08
tabs: true
tags: measure
toc:
  - name: Outer Measure
    subsections:
        - name: Open Intervals
        - name: Closed Intervals
  - name: Measurable Sets
    subsections:
        - name: Sigma-Algebra
        - name: Borel Sets
        - name: Functions
  - name: Measures
    subsections:
        - name: Characteristics
        - name: Lebesgue Measure
  - name: Convergence
bibliography: measure.bib
---

In this post, I'll be going through Chapter 2 in Axler<d-cite key=axler2025></d-cite>. It covers the idea of <i>measures</i>.

As before, we'll let $\mathbb{R}$ denote the <a href="../../measure/useful-definitions/#open-closed">closed</a>, <a href="../../measure/useful-definitions/#total-order">ordered field</a> of real numbers.


## Outer Measure

### Open Intervals
First, we begin with the idea of <i>length</i> for open intervals in $\mathbb{R}$. 

<div id="length"></div>
<div class="definition">
  <strong>Definition (Length).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $I$ be an open interval in $\mathbb{R}$. The <i>length</i> of $I$, denoted by $\ell(I)$, is defined as:

  $$
  \ell(I) =
  \begin{cases}
  b - a & \text{ if } I = (a, b) \text{ for some } a, b \in \mathbb{R} \text{ with } a < b \\
  0  & \text{ if } I = \emptyset \\
  \infty & \text{ if } I = (- \infty, a) \text{ or } I = (a, \infty) \text{ for some } a \in \mathbb{R} \\
  \infty & \text{ if } I = (- \infty, \infty) \\
  \end{cases}
  $$
</div>

We can use this (pretty intuitive!) definition of length to define the <i>outer measure</i>, which gives us a way to describe how big a set is.

<div id="outer-measure"></div>
<div class="definition">
  <strong>Definition (Outer Measure).<d-cite key=axler2025></d-cite></strong>
  <br>
  The <i>outer measure</i> of a set, $A \subseteq \mathbb{R}$, denoted by $\rvert A \rvert$, is defined as:

  $$
  \begin{aligned}
  \rvert A \rvert 
  &= \underset{I_1, I_2, \dots}{\inf} \left\{ \sum_{k = 1}^{\infty} \ell(I_k) \right\}
  \end{aligned}
  $$

  where $I_1, I_2, \dots$ are open intervals in $\mathbf{R}$ such that $A \subseteq \cup_{k = 1}^\infty I_K$. 
</div>

The outer measure satisfies some properties that make it nice to work with. 

<!-- #region outer-measure-properties -->
<div id="outer-measure-properties">
{% tabs outer-measure-properties %}
{% tab outer-measure-properties statement %}
<ol>
<li>Every countable subset of $\mathbb{R}$ has outer measure $0$.</li>
<li><strong>Order Preservation:</strong> $\rvert A \rvert \leq \rvert B \rvert$ for any $A, B \subseteq \mathbb{R}$ with $A \subseteq B$</li>
<li><strong>Translation Invariance:</strong> $\rvert t + A \rvert = \rvert A \rvert$ for any $t \in \mathbb{R}$ and $A \subseteq \mathbb{R}$</li>
<li><strong>Countable Subadditivity:</strong> $\rvert \cup_{k = 1}^\infty A_k \rvert \leq \sum_{k = 1}^\infty \rvert A_k \rvert$ for any sequence of subsets $A_1, A_2, \dots \subseteq \mathbb{R}$.</li>
</ol>
{% endtab %}
{% tab outer-measure-properties proof %}
<strong>Every countable subset of $\mathbb{R}$ has outer measure $0$</strong>

Let $A \subseteq \mathbb{R}$ be a countable subset with $A = \{ a_1, a_2, \dots \}$. Let $\epsilon > 0$ be an arbitrary, positive, real number. For any positive integer $$k \in \mathbb{Z}_+$$, we define:

$$
I_k = \left(a_k - \frac{\epsilon}{2^k}, a_k + \frac{\epsilon}_{2^k} \right)
$$

which is the open interval centered at $a_k$ with length $\frac{\epsilon}_{2^{k-1}}$. Clearly, the union of $I_1, I_2, \dots$ contains $A$. Furthermore, notice that:

$$
\sum_{k = 1}^\infty \ell(I_k) = \sum_{k = 1}^\infty \frac{\epsilon}{2^{k-1}} = 2\epsilon \sum_{k = 1}^\infty \frac{1}{2^k} = 2 \epsilon
$$

We then take the infimum, which means we shrink $\epsilon$ as small as possible, which implies $\rvert A \rvert = 0$. 
<br><br>
<strong>Order Preservation</strong>

Let $I_1, I_2, \dots$ be a sequence of open intervals in $\mathbb{R}$ whose union contains $B$. Since $A \subseteq B$, this union also contains $A$. Thus:

$$
\rvert A \rvert \leq \sum_{k = 1}^\infty \ell(I_k) 
$$

Since the above holds for all such $I_1, I_2, \dots$, it will also hold for the infimum. Thus:

$$
\rvert A \rvert \leq \underset{I_1, I_2, \dots}{\inf} \left\{ \sum_{k = 1}^{\infty} \ell(I_k) \right\} = \rvert B \rvert
$$

<strong>Translation Invariance</strong>

Let $I_1, I_2, \dots$ be a sequence of open intervals whose union contains $A$. Clearly, the union of the sequence of translated open intervals $t + I_1, t + I_2, \dots$ will contain $t + A$. Thus:

$$
\rvert t + A \rvert \sum_{k = 1}^\infty \ell(t + I_k) = \sum_{k = 1}^\infty \ell(I_k)
$$

since translation will not affect the length of the intervals. We can then take the infimum over the right-hand side to see that:

$$
\rvert t + A \rvert \leq \underset{I_1, I_2, \dots}{\inf} \left\{ \sum_{k = 1}^{\infty} \ell(I_k) \right\} = \rvert A \rvert
$$

Now, we notice that $A = t + (-t + A)$. We can repeat the above proof with $(-t + A)$ as $A$ and $-t$ as $t$ to see that:

$$
\rvert t + (-t + A) \rvert \leq \underset{I_1, I_2, \dots}{\inf} \left\{ \sum_{k = 1}^{\infty} \ell(I_k) \right\} = \rvert -t + A \rvert 
$$

This shows the other direction, and so we conclude that $\rvert t + A \rvert = \rvert A \rvert$.
<br><br>
<strong>Countable Additivity</strong>

Assume $$\rvert A_k \rvert < \infty$ for all $k \in \mathbb{Z}_+$$. For each $k \in \mathbb{Z}_+$, define a sequence of open intervals $I_{1,k}, I_{2,k}, \dots$ whose union contains $A_k$. In the event that some fo the $A_k$ overlap, we have the inequality:

$$
\cup_{k = 1}^\infty A_k \subseteq \cup_{k = 1}^\infty \cup_{j = 1}^\infty I_{j, k}
$$

By the order preservation property of the outer measure, we have

$$
\rvert \cup_{k = 1}^\infty A_k \rvert \leq \rvert \cup_{k = 1}^\infty \cup_{j = 1}^\infty I_{j, k}\rvert 
$$

Taking the infimum over the sequences of open intervals yields the result. Note that if the outer measure of any $A_k$ is not finite, then the claim clearly holds.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Note that, since finite sets are also countable, the countable subadditivity property implies finite subadditivity of the outer measure.

### Closed Intervals
Before we can discuss the outer measure of closed (bounded) intervals, we need to introduce a few more tools.

<div id="open-cover"></div>
<div class="definition">
  <strong>Definition (Open Cover/Finite Subcover).<d-cite key=axler2025></d-cite></strong>
  <br>
  Let $A \subseteq \mathbb{R}$, and let $\mathcal{C}$ be a collection of open subsets, $c_1, c_2, \dots$ of $\mathbb{R}$. $\mathcal{C}$ is called an <i>open cover</i> of $A$ if it is contained in the union of all sets in $\mathcal{C}$. That is, if:

  $$
  A \subseteq \cup_{k = 1}^\infty c_k
  $$

  $\mathcal{C}$ is called a <i>finite subcover</i> of $A$ if $A$ is contained in the union of finitely many sets in $\mathcal{C}$. That is, if:

  $$
  A \subseteq \cup_{k = 1}^n c_k
  $$
</div>

As an example, consider $A = [2, \infty)$. The collection $\mathcal{C} = \{ (k, k + 2) \rvert k \in \mathbb{Z}_+ \}$ is an open cover of $A$; however, $A$ does not have a finite subcover. 

We now come to an important theorem in measure theory: the <i>Heine-Borel Theorem</i>.

<div id="heine-borel"></div>
<div class="theorem">
  <strong>Theorem (Heine-Borel).<d-cite key=axler2025></d-cite></strong>
  <br>
{% tabs heine-borel %}
{% tab heine-borel statement %}
Every open cover of a closed, bounded subset of $\mathbb{R}$ has a finite subcover. In other words, any closed and bounded subset of $\mathbb{R}$ is <i>compact</i>. 
{% endtab %}
{% tab heine-borel proof %}
This proof follows <a href="https://en.wikipedia.org/wiki/Heine–Borel_theorem#Proof">the wikipedia article</a>.
Let $F \subseteq \mathbb{R}$ be a closed, bounded subset of $\mathbb{R}$, and let $\mathcal{C}$ be an open cover of $F$. 
<br><br>
<strong>Compact sets are closed.</strong>

Let $S \subseteq \mathbb{R}$ be a compact set, and let $h$ be a limit point of $S$ (i.e. every neighborhood of $h$ contains at least one $x \in S$ such that $x \neq h$). 

Suppose $S$ is not closed. This implies that there exists a limit point $h \notin S$. We will consider an open neighborhood about each $x \in S$, denoted by $N_x$, all of which are disjoint from some neighborhood of $h$, $V_x$ (not necessarily the same one!). Let $$\mathcal{C}'$$ be the collection of such $N_x$. Furthermore, for any $x \in S$, we can find $C \in \mathcal{C}'$ such that $x \in C$. Since each $N_x$ is open, this implies that $\mathcal{C}'$ is an open cover of $S$. 

Consider an arbitrary, finite subcollection of $\mathcal{C}'$, denoted by $\mathcal{K}$. Let $$X_\mathcal{K}$$ denote the set of all $x$ corresponding to $N_x \in \mathcal{K}$. For every $N_x \in \mathcal{K}$, there exists a $V_x$ that $N_x$ does not intersect with. Thus, $$\cup_{x \in X_{\mathcal{K}}} N_x$$ will be disjoint from $$\cap_{x \in X_{\mathcal{K}}} V_x$$. 

Since each $V_x$ is a neighborhood of $h$, there exist open sets $U_x$ such that $h \in U_x \subseteq V_x$ (by definition of neighborhood). It follows that:

$$
h \in \cap_{x \in X_{\mathcal{K}}} U_x \subseteq \cap_{x \in X_{\mathcal{K}}}
$$

which implies that the intersection of such $V_x$ is a neighborhood of $h$. Note, the intersection will not be the set containing only $h$ because we can construct an open ball contained in each $U_x$ (and therefore each $V_k$) that is centered at $h$. The intersection will contain the smallest open ball out of all of the $U_x$.

Since the intersection is a neighborhood of $h$, there exists $x \in \cap_{x \in X_{\mathcal{K}}}$ such that $x \in S$. However, because of how we constructed $\mathcal{C}'$, $x \notin \cup_{x \in X_{\mathcal{K}}} N_x$. Thus, $\mathcal{K}$ is not an open subcover of $S$. This implies $S$ is not compact, which is a contradiction. Thus, all limit points must be in $S$, so $S$ is closed.
<br><br>
<strong>Compact sets are bounded.</strong>

Let $S \subseteq \mathbb{R}$ be a compact set, and let $U_x$ be an open ball centered at some $x \in \mathbb{R}$ with radius $1$. Clearly, $S \subseteq \cup_{x \in S} U_x$, so $U_x$ is an open cover of $S$. 

Consider some finite subcover of $S$ of this union. Let $M$ be the maximum distance between the centers of any two balls (denote these centers as $C_p$ and $C_q$). Take any arbitrary points $p, q \in S$ such that $p$ is in the ball centered at $C_p$, and $q$ is in the ball centered at $C_q$. Then, by the triangle inquality:

$$
d(p, q) \leq d(p, C_p) + d(C_p, C_q) + d(C_q, q) \leq 1 + M + 1 \leq 2 + M
$$

So the diameter of $S$ is bounded by $2 + M$.
<br><br>
<strong>Closed subsets of compact sets are compact.</strong>

Let $K$ be a closed subset of $T \subseteq \mathbb{R}$ where $T$ is compact. Let $\mathcal{C}_K$ be an open cover of $K$. Since $K$ is closed, $U = \mathbb{R} \setminus K$ is an open set, and $\mathcal{C}_T = \mathcal{C}_K \cup \left\{ U \right\}$ is an open cover of $T$. 

Since $T$ is compact, there exists a finite subcover of $T$ of $\mathcal{C}_T$, which we denote by $$\mathcal{C}_T'$$, which also covers $K$. Moreover, $K$ is also covered by $$\mathcal{C}_{K}' = \mathcal{C}_{T}' \setminus \{ U \}$$ since there are no points in $U$ that are in $K$. 

Since $$\mathcal{C}_T'$$ is a finite subcover, $$\mathcal{C}_{K}'$$ is also finite, and every $$C \in \mathcal{C}_{K}'$$ is in $$\mathcal{C}_T'$$, which implies $$C \in \mathcal{C}_T$$. By the definition of $T$ and the fact that $C \notin U$, we have that $C \in \mathcal{C}_K$. Thus, $$\mathcal{C}_K'$$ is a finite subcollection of $C_K$ and therefore a finite subcover of $K$ (i.e. $K$ is compact).
<br><br>
<strong>Closed, bounded sets are compact.</strong>
Let $S \in \mathbb{R}$ be bounded. Let $T_0$ be a closed and bounded interval enclosing $S$. WLOG, we will say that $T_0 = [-a, a]$ for some $a > 0$. 

Suppose $T_0$ is not compact. Let $\mathcal{C}$ be the (infinite) open cover of $T_0$ with no finite subcover of $T_0$. We can divide $T_0$ into two pieces, one of which <i>must</i> have an infinite subcover of $\mathcal{C}$ (if neither did, then the union of the two subcovers would be a subcover of $\mathcal{C}$). Denote this half $T_1$, and note that $T_1 \subset T_0$.

We can repeat the same argument for $T_1$ to get $T_2$ where $T_2 \subset T_1$. We can continue this to construct a sequence of nested intervals:

$$
T_0 \supset T_1 \supset T_2 \supset \dots \supset T_k \supset \dots
$$

The length of interval $T_k$ is the length of $T_0$ divided by $2^k$. Thus:

$$
\underset{k \rightarrow \infty}{\lim} \frac{2a}{2^k} = 0
$$

Let $$(x_k)_k$$ such that $x_k \in T_k$. Clearly, this sequence is Cauchy, since $x_k$ and $x_{k+1}$ can be, at most, $\frac{2a}{2^k} + \frac{2a}{2^{k+1}}$ apart. Since the sequence is Cauchy, it must converge (let $L$ denote this limit). Since the boxes are nested, for any $k$, the sequence will eventually always be inside $T_k$, and each $T_k$ is closed, so $L \in T_k$. 

Recall that $\mathcal{C}$ covers $T_0$, so since $L \in T_k$ and $T_k \subset T_0$, there must exist some open set $U \in \mathcal{C}$ such that $L \in U$. As an open set, one can define an open ball centered at $L$, denoted by $B_L$, that is contained in $U$. 

If we pick a big enough value $k$, then $T_k \subseteq B_L \subseteq U$. However, this implies that the single subset $U$ is an open cover of $T_k$, which means it is a finite subcover of $T_k$. This contradicts the fact that all of the $T_k$s do not have finite subcovers. Thus, $T_0$ must be compact.

By the previous result, since $S$ is closed and a subset of $T_0$, $S$ is also compact.
{% endtab %}
{% endtabs %}
</div>

With the Heine-Borel Theorem, we can show that the outer measure of a closed interval in $\mathbb{R}$ is simply the difference betwen the endpoints.

<!-- #region axler-2-14 -->
<div class="theorem">
  <strong>Theorem 2.14.<d-cite key=axler2025></d-cite></strong>
  <br>
{% tabs axler-2-14 %}
{% tab axler-2-14 statement %}
Let $a, b \in \mathbb{R}$ with $a < b$. The outer measure of $[a, b]$ is given by:

$$
\rvert [a, b] \rvert = b - a
$$
{% endtab %}
{% tab axler-2-14 proof %}
First, consider $\epsilon > 0$. We can define $(a-\epsilon, b + \epsilon), \emptyset, \emptyset, \dots$, which is a sequence of open intervals whose union contains $[a, b]$. Let $u_i$ denote the $i$-th interval in this sequence. Then:

$$
\begin{aligned}
\rvert \cup_{i = 1}^\infty u_i \rvert &\geq \rvert [a, b] \rvert \\
\implies \rvert (a - \epsilon, b + \epsilon) \rvert &\geq \rvert [a, b] \rvert \\
\implies (b - a) + 2 \epsilon &\geq \rvert [a, b] \rvert
\end{aligned}
$$

Taking $\epsilon \rightarrow 0$ yields one direction of the result:

$$
\rvert [a, b] \rvert \leq b - a
$$
<br><br>
Now, let $I_1, I_2, \dots$ be a sequence of open intervals satisfying $[a, b] \subseteq \cup_{k = 1}^\infty I_k$. By the Heine-Borel theorem, there exists a finite subcover of $[a, b]$ from this open cover. That is, there exists some $n \in \mathbb{Z}_+$ such that:

$$
[a, b] \subseteq I_1 \cup \dots \cup I_n
$$

<strong>Base case.</strong>
Suppose that $n = 1$. Clearly:

$$
b - a \leq \ell(I_1) \leq \sum_{k = 1}^\infty \ell(I_k)
$$

<strong>Inductive hypothesis.</strong>
Now, suppose $n > 1$ and that the following holds for all $a, b \in \mathbb{R}$ with $a < b$:

$$
[a, b] \subseteq I_1 \cup \dots \cup I_n \implies \sum_{k = 1}^n \ell(I_k) \geq b - a
$$

Consider appened one more open interval such that:

$$
[a, b] \subseteq I_1 \cup \dots \cup I_n \cup I_{n+1}
$$

WLOG, assume $b \in I_{n+1} = (c, d)$ (we can relabel the intervals to make this true). If $c \leq a$, then the length of $I_{n+1}$ must be at least as large as $b - a$, and it follows that:

$$
b - a \leq \sum_{k = 1}^n \ell(I_k)
$$

The other case is that $c > a$, so we have $a < c < b < d$, which implies:

$$
[a, c] \subseteq I_1 \cup \dots \cup I_n
$$

By assumption, we have that:

$$
c - a \leq \sum_{k = 1}^n \ell(I_k)
$$

We then have:

$$
\begin{aligned}
b - a &\leq d - a \\
&= (c - a) + (d - c) \\
&= (c - a) + \ell(I_{n+1}) \\
&\leq \sum_{k = 1}^{n+1} \ell(I_k)
\end{aligned}
$$

which completes the proof.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Unfortunately, the outer measure does not satisfy <i>additivity</i>, which is a property that would make it behave more intuitively and lend itself naturally to more advanced usage. More specifically, there exist $A, B \subset \mathbb{R}$ such that $A \cap B = \emptyset$ and:

$$
\rvert A \cup B \rvert \neq \rvert A \rvert + \rvert B \rvert
$$ 

Perhaps somewhat surprising, but there does not exist a way to describe the size of a subset of $\mathbb{R}$ that satisfies <i>all</i> of the nice properties we desire. 


<!-- #region measure-properties -->
<div id="theorem-2-22"></div>
<div class="theorem">
  <strong>Theorem 2.22.<d-cite key=axler2025></d-cite></strong>
  <br>
{% tabs axler-2-22 %}
{% tab axler-2-22 statement %}
There does not exist a function, $\mu$, that satisfies all of the following:

<ol>
<li>$\mu$ is a function from the set of all subsets of $\mathbb{R}$ to $[0, \infty]$</li>
<li>$\mu(I) = \ell(I)$ for every open interval $I \subseteq \mathbb{R}$</li>
<li>$\mu\left(\cup_{k = 1}^\infty A_k \right) = \sum_{k = 1}^\infty \mu(A_k)$ for any disjoint sequence of subsets of $\mathbb{R}$</li>
<li>$\mu(t + A) = \mu(A)$ for every $A \subseteq \mathbb{R}$ and every $t \in \mathbb{R}$</li>
</ol>
{% endtab %}
{% tab axler-2-22 proof %}
Let $\mu$ be a function that satisfies all four properties. 

<strong>Some consequences.</strong>
Notice that $\mu(\emptyset) = 0$ by Property 2. Let $A \subseteq B \subseteq \mathbb{R}$. Notice that $B = \cup_{i = 1}^\infty C_i$ where $C_1 = A$, $C_2 = B \setminus A$, and $C_i = \emptyset$ for $i \geq 3$ (i.e. $B$ can be written as a sequence of disjoint sets). By Property 3, we see:

$$
\begin{equation}
\label{eq:order-preserving}
\mu(B) = \mu(A) + \mu(B \setminus A) + \mu(\empty) + \dots = \mu(A) + \mu(B \setminus A) \geq \mu(A)
\end{equation}
$$

Let $a, b \in \mathbb{R}$ with $a < b$. For any choice of $\epsilon > 0$, we have $(a, b) \subseteq [a, b] \subseteq (a - \epsilon, b + \epsilon)$. This implies:

$$
\begin{equation}
\label{eq:closed-size}
b - a \leq \mu([a, b]) \leq b - a + 2 \epsilon \implies \mu([a, b]) = b - a
\end{equation}
$$

Let $A_1, A_2, \dots \subseteq \mathbb{R}$. We consider the sequence of disjoint subsets constructed as $A_1, A_2 \setminus A_1, A_3 \setminus (A_2 \cup A_1), \dots$. We then have:

$$
\begin{aligned}
\mu\left(\cup_{k = 1}^\infty A_k \right) &= \mu\left( A_1 \cup (A_2 \setminus A_1) \cup (A_3 \setminus (A_1 \cup A_2)) \cup \dots\right) \\
&= \mu(A_1) + \mu(A_2 \setminus A_1) + \mu(A_3 \setminus (A_1 \cup A_2)) + \dots & \left(\text{Property } 2 \right)\\
&\leq \sum_{k = 1}^\infty \mu(A_k)
\end{aligned}
$$

<strong>A Contradiction.</strong>
We now show that there exist disjoint $A, B \subseteq \mathbb{R}$ such that $\rvert A \cup B \rvert \neq \rvert A \rvert + \rvert B \rvert$, which contradicts Property 3. 

Let $a \in [-1, 1]$, and define the set of numbers in $[-1, 1]$ that differ from an element in $a$ by a rational number:

$$
\tilde{a} = \left\{ c \in [-1, 1] \rvert a - c \in \mathbb{Q} \right\}
$$

where $\mathbb{Q}$ is the set of rational numbers. 

Suppose $a, b \in [-1, 1]$ and $\tilde{a} \cap \tilde{b} \neq \emptyset$. Then it must be the case that $\tilde{a} = \tilde{b}$. This follows from the fact that, if there exists a $d \in \tilde{a} \cap \tilde{b}$, then $a - d, b - d \in \mathbb{Q}$. If we subtract the two, we see that $(a - d) - (b - d) = a - b \in \mathbb{Q}$ since the difference of rational numbers is rational. Consider now the equation $a - c = (a - b) + (b - c)$, which we obtain by adding and subtracting $b$. If $c \in [-1, 1]$, then $a - c$ is rational if, and only if, $b - c$ is also rational. Thus, $\tilde{a} = \tilde{b}$. 

For every $a \in [-1, 1]$, $a \in \tilde{a}$ (since $a$ differs from $a$ by $0$, which is rational). Thus:

$$
[-1, 1] = \bigcup_{a \in [-1, 1]} \tilde{a}
$$

We define the set, $V$, which is constructed by taking a single element from each of the distinct sets in $\{ \tilde{a} \rvert a \in [-1, 1] \}$. Clearly, $a - v \in \mathbb{Q}$ because $v \in \tilde{a}$, and all elements in $\tilde{a}$ differ from $a$ by a rational number. 

We then define $r_1, r_2, \dots$, which is a sequence of distinct rational number that satisfy:

$$
[-2, 2] \cap \mathbb{Q} = \{ r_1, r_2, \dots \}
$$

Since $a - v \in \mathbb{Q}$, it follows that $a = r_k + v \in r_k + V$ for some $k \in \mathbb{Z}_+$. Thus:

$$
[-1, 1] \subseteq \bigcup_{k = 1}^\infty (r_k + V)
$$

By Property 3 and Eq. \eqref{eq:order-preserving}, we have that:

$$
\mu([-1, 1]) \leq \sum_{k = 1}^\infty \mu(r_k + V)
$$

And by Eq. \eqref{eq:closed-size}:

$$
\mu([-1, 1]) = 2
$$

Then, by Property 4, we can translate the subsets on the right-hand side to get:

$$
2 \leq \sum_{k = 1}^\infty \mu(V) \implies \mu(V) > 0
$$

Now, if there exists $t \in (r_j + V) \cap (r_k + V)$, then it must be the case that $t = r_j + v_1 = r_k + v_2$ for some $v_1, v_2 \in V$. Subtracting $v_2$ and $r_j$ from both sides yields the equality $v_1 - v_2 = r_k - r_j$. Clearly, $v_1 - v_2 = r_k - r_j \in \mathbb{Q}$, since $r_k, r_j \in \mathbb{Q}$. However, by the construction of $V$, $v_1 = v_2$, so $r_k = r_j$, and thus $j = k$. This implies $(r_1 + V), (r_2 + V), \dots$ are disjoint sets.

For arbitrary $n \in \mathbb{Z}_+$, we have:

$$
\bigcup_{k = 1}^n (r_k + V) \subseteq [-3, 3]
$$

since $V \subseteq [-1, 1]$ and $r_k \in [-2, 2]$ for all $k$. By Property 2 and Property 3, we have:

$$
\mu\left( \bigcup_{k = 1}^n (r_k + V) \right) \leq 6
$$

Let's pick $n \in \mathbb{Z}_+$ such that $n \mu(V) > 6$. By Property 3 and Property 4:

$$
\begin{equation}
\label{eq:2-22}
\mu\left( \bigcup_{k = 1}^n (r_k + V)  \right) \leq 6 < \sum_{k = 1}^n \mu(r_k + V) = \sum_{k = 1}^n \mu(V) = n \mu(V)
\end{equation}
$$

This is basically the end of the proof, since we have constructed an example of such disjoint subsets. We can do a little more explanation, though.

If it were the case that $\mu(A \cup B) = \mu(A) + \mu(B)$ for any disjoint subsets $A, B \subseteq \mathbb{R}$, then (by induction), we would have that:

$$
\mu\left(\bigcup_{k = 1}^n A_k \right) = \sum_{k = 1}^n \mu(A_k)
$$

for any disjoint subsets $A_1, \dots, A_n \subseteq \mathbb{R}$. However, Eq. \eqref{eq:2-22} shows this is not the case. Thus, Property 3 is violated.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

---

## Measurable Spaces and Functions
The failure to find a function that satisfies all of the properties listed in <a href="#theorem-2-22">Theorem 2.22</a> means we should probably drop one so that we can continue with analysis... Because Properties 2-4 are so important to our intuition and the utility of this idea of size, we must be satisfied with only these three. We now come to some fundamental building blocks of measure theory.


### Sigma-Algebra
We'll begin with a $\sigma$-algebra. This is simply a collection of subsets of some other set. 

<div id="sigma-field"></div>
<div class="definition">
  <strong>Definition ($\sigma$-Algebra).</strong>
  <br>
  Let $X$ be a set, and let $\mathcal{P}(X)$ denote its power set (the set of all possible subsets of $X$). A $\sigma$<i>-algebra</i> is any subset $\mathcal{S} \subseteq \mathcal{P}(X)$ satisfying the following:
  <ol>
    <li> $X \in \mathcal{S}$ </li>
    <li> <strong>Closed Under Complementation:</strong> $A \in \mathcal{S} \implies A^c \in \mathcal{S}$ </li>
    <li> <strong>Closed Under Countable Unions:</strong> $A_1, A_2, \dots \in \mathcal{S} \implies A = \bigcup_{i = 1}^\infty A_i \in \mathcal{S}$ </li>
  </ol>
</div>

An example may make it a bit more concrete in one's mind...

<div class="example">
  <strong>Example ($\sigma$-Algebra).</strong>
  <br>
  Let $X = \{ 1, 2, 3\}$. One $\sigma$-algebra on $X$ is $\mathcal{S} = \{ \emptyset, \{ 1 \}, \{2, 3 \}, \{1, 2, 3 \} \}$. 
  <br>
  It is easy to see that the first two properties are satisfied. The entirety of $X$ is in $\mathcal{S}$ by construction. Any element in $\mathcal{S}$, its complement is also in $\mathcal{S}$ (the complement of $\emptyset$ is $X$, and the complement of $\{ 1\}$ is $\{ 2, 3 \}$). 
  <br>
  Taking any countable union of elements of $\mathcal{S}$ also yields an element of $\mathcal{S}$ (the union of any $A \in \mathcal{S}$ and the empty set is the set itself; the complement of the entirety of $X$ and any $A \in \mathcal{S}$ will just be $X$; and the union of $\{ 1\}$ and $\{ 2, 3\}$ is $X$). This satisfies the third property, which completes the proof that $\mathcal{S}$ is a $\sigma$-algebra on $X$.
</div>

Notice that the first and second properties in the above definition imply that $\emptyset \in \mathcal{S}$ as well. The properties also imply that a $\sigma$-algebra must be closed under countable intersection. That is, $\cap_{i = 1}^\infty A_i \in \mathcal{S}$ for some sequence of $A_1, A_2, \dots \in \mathcal{S}$. 

A $\sigma$-field is a generalization of the concept of an <i>algebra</i> (also called a <i>field</i>). 

<div id="algebra"></div>
<div class="definition">
  <strong>Definition (Algebra).</strong>
  <br>
  Let $X$ be a set. A collection, $\mathcal{A}$, of subsets of $X$ is an <i>algebra</i> (or <i>field</i>) if the following are satisfied:
  <ol>
  <li>$X \in \mathcal{A}$</li>
  <li><strong>Closed Under Complementation:</strong> $A \in \mathcal{A} \implies A^c \in \mathcal{A}$</li>
  <li><strong>Closed Under Finite Unions:</strong> $A, B \in \mathcal{A} \implies A \cup B \in \mathcal{A}$</li>
  </ol>
</div>

Now we can define <i>measurable spaces</i>!

<div id="measurable-space"></div>
<div class="definition">
  <strong>Definition (Measurable Space).</strong>
  <br>
  Let $X$ be some set, and let $\mathcal{S}$ be a $\sigma$-algebra on $X$. The tuple $(X, \mathcal{S})$ is called a <i>measurable space</i>, and any element of $\mathcal{S}$ is called an <i>$\mathcal{S}$-measurable set</i>.
</div>

We'll come back to this definition later when we discuss measures, but a measurable space is just a space that <i>could</i> be assigned a measure. 

Let's finish up this sub-section by introducing <i>topological spaces</i> and <i>Borel sets</i>.

<div id="topology"></div>
<div class="definition">
  <strong>Definition (Topology).</strong>
  <br>
  Let $X$ be a non-empty space. A <i>topology</i>, $\tau$, on $X$ is any collection of subsets of $X$ that satisfy:
  <ol>
    <li>The empty set, $\emptyset$, and the entirety of $X$ are in $\tau$.</li>
    <li>The union (finite or infinite) of any subset in $\tau$ is also in $\tau$.</li>
    <li>The intersection of a finite number of subsets in $\tau$ is also in $\tau$.</li>
  </ol>
  We call the tuple $(X, \tau)$ a <i>topological space</i>.
</div>


### Borel Sets
A lot of proofs I've come across implicitly define random variables on probability spaces (to be covered later!) where the $\sigma$-algebra is the Borel $\sigma$-algebra on $\mathbb{R}$.

<div id="borel-set"></div>
<div class="definition">
  <strong>Definition (Borel Set).</strong>
  <br>
  A <i>Borel set</i> on a topological space, $X$, denoted by $\mathcal{B}(X)$, is any subset of $X$ that can be constructed from open sets on that space in $X$ via countable unions, countable intersections, and set differences. 
  <br>
  The smallest $\sigma$-algebra containing all open subsets of $X$ is called the <i>Borel $\sigma$-algebra</i> (or <i>Borel field</i>, or the collection of <i>Borel subsets</i>).
</div>

This definition is a bit tricky to develop intuition for. The Borel $\sigma$-algebra is just the collection of <i>all possible open sets</i> in a given space, $X$ (usually $\mathbb{R}$). 

An important Borel $\sigma$-algebra that will come up again when we discuss measures and probability is the Borel $\sigma$-algebra on the real line. Several examples follow from our definition:

<ol>
  <li>Any <i>closed</i> subset of $\mathbb{R}$ is a Borel set because $\sigma$-algebras are closed under complementation.</li>
  <li>Any <i>countable</i> subset of $\mathbb{R}$ is a Borel set because $\sigma$-algebras are closed under countable unions, and a single point is a closed subset of $\mathbb{R}$.</li>
  <li>Any half-open interval is a Borel set because $\sigma$-algebras are closed under countable intersections.</li>
</ol>

Borel sets on $\mathbb{R}$ can also be extended to $[-\infty, \infty]$. However there are many sequences that are <i>not</i> Borel sets. For example, the set of all countable unions of countable intersections of open subsets of $\mathbb{R}$ is not the set of Borel sets. 

Along with the Borel set and the $\sigma$-algebra is the <i>semialgebra</i>.

<div id="semialgebra"></div>
<div class="definition">
  <strong>Definition (Semialgebra).</strong>
  <br>
  Let $\mathcal{S}$ be a collection of sets. $\mathcal{S}$ is called a <i>semialgebra</i> if it satisfies the following properties:
  <ol>
  <li><strong>Closed Under Intersection:</strong> $S,T \in \mathcal{S} \implies S \cap T \in \mathcal{S}$</li>
  <li><strong>(Sort Of) Closed Under Complementation:</strong> $S \in \mathcal{S} \implies S^c$ is a finite disjoint union of $T \in \mathcal{S}$</li>
  </ol>
</div>

This concept will not be as useful in later discussions, but we include it for completeness. An example of a semialgebra is the union of $\{ \emptyset \}$ and the collection of sets that can be written as:

$$
(a_1, b_1] \times \dots \times (a_d, b_d] \subset \mathbb{R}^d \hspace{5mm} \text{for } -\infty \leq a_i < b_i \leq \infty
$$

Given a semialgebra, $\mathcal{S}$, the collection of finite disjoint unions of sets in $\mathcal{S}$ forms an algebra called the <i>algebra generated by $\mathcal{S}$</i>.

### Functions
We now need to define a concept that is at the crux of our discussions of mappings: the <i>inverse image</i>.

<div id="inverse-image"></div>
<div class="definition">
  <strong>Definition (Inverse Image/Pre-Image).</strong>
  <br>
  Let $f: X \rightarrow Y$ be some function, and let $A \subseteq Y$. The <i>inverse image</i>, also called the <i>pre-image</i>, of subset $A$ is defined as the set:

  $$
  f^{-1}(A) = \left\{x \in X \rvert f(x) \in A \right\}
  \nonumber
  $$

  The inverse image satisfies nice properties:

  <ol>
    <li>For any $A \subseteq Y$, $f^{-1}(Y \setminus A) = X \setminus f^{-1}(A)$</li>
    <li>For any set $\mathcal{A}$ of subsets of $Y$: $f^{-1}(\cup_{A \in \mathcal{A}} A) = \cup_{A \in \mathcal{A}}f^{-1}(A)$</li>
    <li>For any set $\mathcal{A}$ of subsets of $Y$: $f^{-1}(\cap_{A \in \mathcal{A}} A) = \cap_{A \in \mathcal{A}}f^{-1}(A)$</li>
    <li>For function $g: Y \rightarrow W$: $(g \circ f)^{-1}(A) = f^{-1}(g^{-1}(A))$s for any $A \subseteq W$</li>
  </ol>
</div>

In words, the inverse image of a subset $A$ of $Y$ under function $f$ is the subset of elements in the domain $X$ that map to elements in $A$. It's important to note that the inverse image of the whole of $Y$ does not necessarily have to be the whole of $X$!

We now introduce a definition that describes what it means for functions of a certain type to be "nice" with respect to a $\sigma$-field.

<div id="measurable-function"></div>
<div class="definition">
  <strong>Definition (Measurable Function).</strong>
  <br>
  Let $(X, \mathcal{S})$ be a measurable space, and let $f: X \rightarrow [-\infty, \infty]$ be a function mapping to the extended real line. We say $f$ is <i>$\mathcal{S}$-measurable</i> if $f^{-1}(B) \in \mathcal{S}$ for every Borel set $B \subseteq [-\infty, \infty]$. 
  Any function from $X$ to $\mathbb{R}$ is $\mathcal{S}$-measurable if $\mathcal{S} = \mathcal{P}(X)$, the power set of $X$.
  <br> 
  Furthermore, for $\mathcal{S}$-measurable functions $f,g: X \rightarrow \mathbb{R}$:
  <ol>
    <li>$f+g$, $f-g$, and $fg$ are $\mathcal{S}$-measurable.</li>
    <li>$f/g$ is $\mathcal{S}$-measurable if $g(x) \neq 0$ for all $x \in X$.</li>
  </ol>
  More generally, for measurable spaces $(X, \mathcal{S})$ and $(Y, \mathcal{S}')$, $f: X \rightarrow Y$ is $(\mathcal{S}, \mathcal{S}')$-measurable if, for all $E \in \mathcal{S}'$, we have $f^{-1}(E) \in \mathcal{S}$. 
</div>

The basic idea behind an $\mathcal{S}$-measurable function is that we should be able to achieve any Borel set as output for <i>some</i> part of $\mathcal{S}$, which is in its domain (since $\text{dom}(f) = X$). It is important to remember that measurability is with respect to the $\sigma$-fields of the two measure spaces of interest.

To put it intuitively, a measurable function $f$ needs to take on values that "make sense" with respect to the $\sigma$-field of interest. For example, only constant functions are measurable with respect to the trivial $\sigma$-field $\{ \emptyset, \Omega \}$ for some $\Omega$. In addition, we have the following claim:

<div class="theorem">
<strong>Claim.</strong>
{% tabs const-funcs %}
{% tab const-funcs statement %}
Constant functions are measurable with respect to <i>any</i> $\sigma$-field.
{% endtab %}
{% tab const-funcs proof %}
Suppose we have measurable spaces $(X, \mathcal{S})$ and $(Y, \mathcal{S}')$. Let $\mathcal{S} = \{ \emptyset, \Omega \}$, and suppose we have non-constant function $f: X \rightarrow Y$. That is, there exist $a, b \in \Omega$ such that $f(a), f(b) \in \mathcal{S}'$ and $f(a) \neq f(b)$. 

Consider the pre-image of one of these points. We know that $f^{-1}(f(a)) = a \notin \mathcal{S}$ since $a$ is neither the null set nor the entirety of $\Omega$ (since we also have $b$ and, necessarily, $a \neq b$). 

To prove the second claim, consider $\mathcal{S} = \{ \emptyset, \Omega \}$ and arbitrary $\mathcal{S}'$ in the previous set-up. Since $f$ is constant, it must be the case that $f(x) = a$ for all $x \in X$ and some $a$. Pick any $s \in \mathcal{S}'$. If $a \in s$, then $f^{-1}(s) = \Omega$, since any input value maps to $a$ ($f$ is constant). If $a \notin s$, then $f^{-1}(s) = \Omega^c = \emptyset$ by the same argument. 

Thus, for any $s \in \mathcal{S}'$, $f^{-1}(s) \in \mathcal{S}$, implying that $f$ is $(\mathcal{S}, \mathcal{S}')$-measurable for <i>any</i> $\mathcal{S}'$. 
{% endtab %}
{% endtabs %}
</div>

To check whether a function is $\mathcal{S}$-measurable, it is sufficient to check whether $$f^{-1}((a, \infty]) = \{ x \in X \rvert f(x) > a \} \in \mathcal{S}$$ for all $a \in \mathbb{R}$. 

Furthermore, in the special case that $X \subseteq \mathbb{R}$ and $\mathcal{S}$ is the set of Borel subsets of $\mathbb{R}$ that are contained in $X$, then a function $f: X \rightarrow \mathbb{R}$ is called <strong>Borel measurable</strong> if $f^{-1}(B)$ is a Borel set for all Borel sets $B \subseteq \mathbb{R}$. It can be shown that any <i>continuous</i> or <i>increasing</i> function $f: X \rightarrow \mathbb{R}$ where $X$ is a Borel subset of $\mathbb{R}$ is Borel measurable (see pg. 33 for a proof).

Later, it will be useful to have a sense of limits of measurable functions. 

<div id="limit-function"></div>
<div class="definition">
  <strong>Definition (Limit of Measurable Functions).</strong>
  <br>
  Let $(X, \mathcal{S})$ be a measure space, and let $f_1, f_2, \dots$ be a sequence of $\mathcal{S}$-measurable functions with $f_i: X \rightarrow \mathbb{R}$. Assume $\underset{k \rightarrow \infty}{\lim} f_k(x)$ exists for all $x \in X$. The function $f: X \rightarrow \mathbb{R}$, defined as:

  $$
  f(x) = \underset{k \rightarrow \infty}{\lim} f_k(x)
  $$

  is $\mathcal{S}$-measurable.
</div>

A similar result holds if $f_1, f_2, \dots$ are $\mathcal{S}$-measurable functions from $X$ to $[-\infty, \infty]$. Then $g, h: X \rightarrow [-\infty, \infty]$ defined as:

$$
g(x) = \inf\left\{ f_k(x) \rvert k \in \mathbb{Z}_+ \right\} \hspace{5mm} h(x) = \sup \left\{ f_k(x) \rvert k \in \mathbb{Z}_+ \right\}
$$

are also $\mathcal{S}$-measurable. 

<strong>Simple functions</strong> are ones that take on only finitely many values. For example, a constant function is simple because it takes on only one value. We have the following feature of approximating measurable function with simple functions.

In other words, for a measurable space $(X, \mathcal{S})$, a simple function $f: X \rightarrow \mathbb{R}$ can be written as:

$$
f = c_1 \chi_{E_1} + \dots + c_n \chi_{E_n}
$$

where $c_1, \dots, c_n \in \mathbb{R}$ are the distinct values that $f$ can take on, and $E_k = f^{-1}(\{ c_k \})$. We also have the following result that will help us prove other things later.

<div id="simple-properties"></div>
<div class="theorem">
  <strong>Theorem.</strong>
  <br>
  Let $(X, \mathcal{S})$ be a measurable space, and let $f: X \rightarrow [-\infty, \infty]$ be an $\mathcal{S}$-measurable function. There exists a sequence $f_1, f_2, \dots$ of functions with $f_i: X \rightarrow \mathbb{R}$ that satisfies:

  <ol>
  <li>$f_i$ is a simple and $\mathcal{S}$-measurable for all $i$</li>
  <li>for all $i \in \mathbb{Z}_+$ and for all $x \in X$, $\rvert f_i(x) \rvert \leq \rvert f_{i + 1}(x) \rvert \leq \rvert f(x) \rvert$</li>
  <li>for every $x \in X$, $\underset{i \rightarrow \infty}{\lim} f_i(x) = f(x)$</li>
  <li>if $f$ is bounded, then $f_1, f_2, \dots$ converges uniformly on $X$ to $f$</li>
  </ol>
</div>



---

## Measures
We have finally come to the star of our discussion: the <i>measure</i>. A measure is a function that assigns a "size" to sets (it is similar to the idea of length for intervals or area for two dimensional regions).

<div id="measure"></div>
<div class="definition">
  <strong>Definition (Measure).</strong>
  <br>
  Let $(X, \mathcal{S})$ be a measure space. A function $\mu: \mathcal{S} \rightarrow [0, \infty]$ is called a <i>measure</i> on $(X, \mathcal{S})$ if:
  <ol>
    <li>$\mu(\emptyset) = 0$</li>
    <li>$\mu\left(\bigcup_{i = 1}^\infty A_i \right) = \sum_{i = 1}^\infty \mu(A_i)$ for every disjoint (i.e. $A_i \cap A_j = \emptyset$ for all $i \neq j$) sequence $A_1, A_2, \dots$ of sets in $\mathcal{S}$</li>
  </ol>
</div>

With this definition, we define a <strong>measure space</strong>, which is the tuple $(X, \mathcal{S}, \mu)$. 

### Properties
For measure space $(X, \mathcal{S}, \mu)$ and $A, B \in \mathcal{S}$ such that $A \subseteq B$, we have that $\mu(A) \leq \mu(B)$ (<strong>order preserving</strong>) and $\mu(B \setminus A) = \mu(B) - \mu(A)$ (assuming that $\mu(A)$ is finite). 

We also have the additional property of <strong>countable subadditivity</strong>, which is basically a generalization of Boole's inequality:

$$
\mu\left(\bigcup_{i = 1}^\infty A_i \right) \leq \sum_{i = 1}^\infty \mu(A_i)
$$

for any sequence of sets $A_1, A_2, \dots \in \mathcal{S}$. Measures also satisfy:

$$
\mu(A \cup B) = \mu(A) + \mu(B) - \mu(A \cap B)$$

when $\mu(A \cap B)$ is finite.

We also have that the increasing union of countably many sets in a measurable space has measure equal to the limit of the measures of the sets in the sequence. The same holds for a decreasing intersection. More formally, for measure space, $(X, \mathcal{S}, \mu)$, an increasing sequence of sets in $\mathcal{S}$, $E_1 \subseteq E_2 \subseteq \dots$, and a decreasing sequence of sets in $\mathcal{S}$, $F_1 \supseteq F_2 \supseteq \dots$, we have:

$$
\begin{equation}
\label{eq:increase-decrease}
\mu\left(\bigcup_{k = 1}^\infty E_k \right) = \underset{k \rightarrow \infty}{\lim} \mu(E_k)
\hspace{5mm}
\mu\left(\bigcap_{k = 1}^\infty F_k \right) = \underset{k \rightarrow \infty}{\lim} \mu(F_k)
\end{equation}
$$

<aside><p>We also need $\mu(F_1) < \infty$.</p></aside>




If we have two $\sigma$-finite (see <a href="#sigma-finite">below</a>) measure spaces, $(X, \mathcal{S}, \mu_1)$ and $(Y, \mathcal{T}, \mu_2)$, we can define two addition sets:

$$
\begin{aligned}
\Omega &= X \times Y = \{ (x, y): x \in X, y \in Y\} \\
\mathcal{U} &= \{ S \times T: S \in \mathcal{S}, T \in \mathcal{T}\}
\end{aligned} 
$$

Sets $U \in \mathcal{U}$ are <i>rectangles</i>. Let $\mathcal{F} = \mathcal{S} \times \mathcal{T}$ be the $\sigma$-algebra generated by $\mathcal{U}$. The unique measure $\mu = \mu_1 \times \mu_2$ on $\mathcal{F}$ defined as $\mu(S \times T) = \mu_1(S) \mu_2(T)$ is called a <i>product measure</i>. This result can be extended to finitely many $\sigma$-finite measurable spaces. 

### Characteristics
Measures can be characterized in a variety of ways. First, consider the $\sigma$-finite measure.

<div id="sigma-finite"></div>
<div class="definition">
  <strong>Definition ($\sigma$-Finite).</strong>
  <br>
  Let $(X, \mathcal{S})$ be a measure space, and let $\mu$ be a measure defined on it. We call $\mu$ a <i>$\sigma$-finite measure</i> if any of the following are satisfied:
  <ul>
    <li>There exist countably many $A_1, A_2, \dots \in \mathcal{S}$ with $\mu(A_n) < \infty$ for all $n \in \mathbb{N}$ such that $\bigcap_{n \in \mathbb{N}} = X$. That is, $X$ can be covered with the intersection of countably many measurable sets in $\mathcal{S}$.</li>
    <li>There exist disjoint and countably many $B_1, B_2, \dots \in \mathcal{S}$ with $\mu(B_n) < \infty$ for all $n \in \mathbb{N}$ such that $\bigcup_{n \in \mathbb{N}} = X$. That is, $X$ can be covered by the union of countably many disjoint measurable sets in $\mathcal{S}$.</li>
    <li>There exist countably many $C_1, C_2, \dots \in \mathcal{S}$ with $C_1 \subset C_2 \subset \dots$ with $\mu(C_n) < \infty$ for all $n \in \mathbb{N}$ such that $\bigcup_{n \in \mathbb{N}} C_n = X$. That is, $X$ can be covered with the union of a countable monotone sequence of measurable sets in $\mathcal{S}$.</li>
    <li>There exists a function $f$ such that $f(x) > 0$ for all $x \in X$ and $\int f(x) \mu(dx) < \infty$. That is, there exists a strictly positive function with finite integral that is measurable with respect to $\mu$.</li>
  </ul>
</div>

We can also define a sense of continuity to measures.

<div id="absolute-continuity"></div>
<div class="definition">
  <strong>Definition (Absolute Continuity).</strong>
  <br>
  Let $\mu$ be a measure on the Borel subsets of $\mathbb{R}$. We call $\mu$ <i>absolutely continuous</i> with respect to the Lebesgue measure, $\lambda$, (see below for definition), if, for every $\lambda$-measurable set $A$, $\lambda(A) = 0$ implies $\mu(A) = 0$. This condition is denoted by $\mu << \lambda$, and we say that $\mu$ is <i>dominated</i> by $\lambda$. 
</div>

Measures can also be "coarsened" by restricting the $\sigma$-field on which they operate.

<div id="restricted-measure"></div>
<div class="definition">
  <strong>Definition (Restricted Measure).</strong>
  <br>
  Let $(\Omega, \mathcal{F}, \mu)$ be a measure space. Let $\mathcal{F}'$ be a sub-$\sigma$-field of $\mathcal{F}$. The <i>restricted measure</i> of $\mu$ to $\mathcal{F}'$ is the mapping $\nu: \mathcal{F}' \rightarrow \mathbb{R} \cup \{ -\infty, +\infty\}$ such that $\nu(E') = \mu(E')$ for all $E' \in \mathcal{F}'$. 
</div>

A restricted measure is basically the original measure but its domain is shrunken to whatever sub-$\sigma$-field it is restricted to. Measures also satisfy several properties.

<div id="theorem-1-1-1"></div>
<div class="theorem">
  <strong>Theorem 1.1.1.<d-cite key=durett2019></d-cite></strong>
  {% tabs theorem-1-1-1 %}
  {% tab theorem-1-1-1 statement %}
  Let $\mu$ be a measure on $(\Omega, \mathcal{F})$, and let $A_i \uparrow A$ denote $A_1 \subset A_2 \subset \dots$ with $\cup_i A_i = A$. The measure $\mu$ satisfies the following:
  <ul>
  <li>Monotonicity: $A \subset B \implies \mu(A) \leq \mu(B)$</li>
  <li>Subadditivity: $A \subset \cup_{m = 1}^\infty A_m \implies \mu(A) \leq \sum_{m= 1}^\infty \mu(A_m)$</li>
  <li>Continuity From Below: $A_i \uparrow A \implies \mu(A_i) \uparrow \mu(A)$</li>
  <li>Continuity From Above: $A_i \downarrow A \implies \mu(A_i) \downarrow \mu(A)$</li>
  </ul>
  {% endtab %}
  {% tab theorem-1-1-1 proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

A sense of "convergence" with respect to a measure can be defined for measurable functions.

<div id="in-measure"></div>
<div class="definition">
  <strong>Definition (In Measure).</strong>
  <br>
  Let $\mu$ be a $\sigma$-finite probability measure, and let $f$ be a measurable function and let $\{ f_n \}_{n = 1}^\infty$ be a sequence of measurable functions. We say that $f_n \rightarrow f$ <i>in measure</i> if, for any $e > 0$, we have that:
  $$
  \mu\left(\{ x : \rvert f_n(x) - f(x) \rvert > e \} \right) \rightarrow 0 \hspace{5mm} \text{ as } n \rightarrow \infty
  $$
</div>

Before we can move on to some of the core concepts in probability theory, we need one more definition. 

<div id="measurable-map"></div>
<div class="definition">
  <strong>Definition (Measurable Map).</strong>
  <br>
  Let $(\Omega, \mathcal{F})$ and $(S, \mathcal{S})$ be two measurable spaces. The function $X: \Omega \rightarrow S$ is a <i>measurable map</i> (from $(\Omega, \mathcal{F})$ to $(S, \mathcal{S})$) if, for all $B \in \mathcal{S}$:
  $$
  X^{-1}(B) = \{ \omega: X(\omega) \in B \} \in \mathcal{F} 
  \nonumber
  $$
</div>

### Lebesgue Measure
An example of a measure is the <strong>Lebesgue measure</strong>, which is the outer measure on $(\mathbb{R}, \mathcal{B})$, where $\mathcal{B}$ is the $\sigma$-algebra of Borel subsets of $\mathbb{R}$. The outer measure is <i>not</i> a measure in general, but it is when it is applied to Borel sets. 

A set $A \subseteq \mathbb{R}$ is called a <strong>Lebesgue measurable set</strong> if there exists a Borel set $B \subset A$ satisfying $\rvert A \setminus B \rvert = 0$. All Borel sets are Lebesgue measurable (set $B = A$). Once can also define a Lebesgue measurable set in many, equivalent ways.

<aside><p>Sometimes <i>Lebesgue measure</i> refers to the outer measure restricted to the $\sigma$-algebra of Lebesgue measurable subsets of $\mathbb{R}$.</p></aside>

<div id="lebesgue-measurable-set"></div>
<div class="theorem">
  <strong>Theorem 2.71 (Lebesgue Measurable Equivalences).</strong>
  <br>
  Let $A \subseteq \mathbb{R}$. The following statements are equivalent:

  <ol>
  <li>$A$ is Lebesgue measurable</li>
  <li>For each $\epsilon > 0$, there exists a closed set $F \subseteq A$ such that $\rvert A \setminus F \rvert < \epsilon$</li>
  <li>These exist closed sets $F_1, F_2, \dots$ contained in $A$ such that $\left\rvert A \setminus \cup_{k = 1}^\infty F_k \right\rvert = 0$</li>
  <li>There exists a Borel set $B \subseteq A$ such that $\rvert A \setminus B \rvert = 0$</li>
  <li>For each $\epsilon > 0$, there exists an open set $G \supseteq A$ such that $\rvert G \setminus A \rvert < \epsilon$</li>
  <li>These exist open sets $G_1, G_2, \dots$ containing $A$ such that $\left\rvert \left( \cap_{k = 1}^\infty G_k \right) \setminus A \right\rvert = 0$</li>
  <li>These exists a Borel set $B \supseteq A$ such that $\rvert B \setminus A \rvert = 0$</li>
  </ol>
</div>

In words, a Lebesgue measurable set is an element of the union of the Borel set and the set of all subsets of $\mathbb{R}$ with outer measure zero. 

We can use the above definition to define <i>Lebesgue measurable functions</i>.

<div id="lebesgue-measurable-function"></div>
<div class="definition">
  <strong>Definition (Lebesgue Measurable Function).</strong>
  <br>
  We call a function $f: A \rightarrow \mathbb{R}$ with $A \subseteq \mathbb{R}$ <i>Lebesgue measurable</i> if $f^{-1}(B)$ is a Lebesgue measurable set for every Borel set $B \subseteq \mathbb{R}$. 
</div>

The above definition means that, if $A$ is a Lebesgue measurable subset of $\mathbb{R}$, then $f$ is $\mathcal{L}$-measurable, where $\mathcal{L}$ is the $\sigma$-algebra of all Lebesgue measurable subsets of $A$. 

We have the following result that will allow us to just consider Borel measurable functions if we think it's okay to ignore the sets of outer measure zero. First, we'll define the <strong>characteristic function</strong>. Let $E \subseteq X$. The characteristic function of $E$, denoted by $\chi_E: X \rightarrow \mathbb{R}$, is defined as:

$$
\chi_E(x) = \begin{cases}
1 & \text{if } x \in E \\
0 & \text{if } x \notin E
\end{cases}
$$

<div class="theorem">
  <strong>Theorem 2.95.</strong>
  <br>
  {% tabs theorem-2-95 %}
  {% tab theorem-2-95 statement %}
  Let $f: \mathbb{R} \rightarrow \mathbb{R}$ be a Lebesgue measurable function. There exists a Borel measurable function $g: \mathbb{R} \rightarrow \mathbb{R}$ satisfying:

  $$
  \rvert \left\{ x \in \mathbb{R} \rvert g(x) \neq f(x) \right\} \rvert = 0
  $$
  {% endtab %}
  {% tab theorem-2-95 proof %}
  TO DO.
  <!-- By the equivalences given <a href="#simple-properties">here</a>, there exists a sequence, $f_1, f_2, \dots$, of Lebesgue measurable, simple functions $f_i: \mathbb{R} \rightarrow \mathbb{R}$ that converges pointwise on $\mathbb{R}$ to $f$. 

  Let $k \in \mathbb{Z}_+$. Since the $f_i$ functions are simple, there exist $c_1, \dots, c_n \in \mathbb{R}$ and disjoint, Lebesgue measurable sets $A_1, \dots, A_n \subseteq \mathbb{R}$ such that:

  $$
  f_k = c_1 \chi_{A_1} + \dots + c_n \chi_{A_n}
  $$

  By the first and fourth equivalences given <a href="#lebesgue-measurable-set">here</a> for Lebesgue measurable sets, we have that, for each $j \in \{ 1, \dots, n \}$, there exists a Borel set $B_j \subseteq A_j$ satisfying $\rvert A_j \setminus B_j \rvert = 0$. We define:

  $$
  g_k = c_1 \chi_{B_1} + \dots + c_n \chi_{B_n}
  $$

  Clearly, $\rvert \{ x \in \mathbb{R} \rvert g_k(x) \neq f_k(x) \} \rvert = 0$ because $g_k$ takes on the exact same values as $f_k$ on subsets of the corresponding sets in the domain of $f$. Furthermore, $g_k$ is Borel measurable (proof?).

  Suppose that $x \notin \bigcup_{i = 1}^\infty \left\{ x \in \mathbb{R} \rvert g_i(x) \neq f_i(x)  \right\}$. This implies that  -->
  {% endtab %}
  {% endtabs %}
</div>

---

## Convergence
Convergence gives us a way to describe the behavior of a sequence of functions. In what follows, we will let $(X, \mathcal{S})$ be our measurable space with set $X$ and $\sigma$-algebra $\mathcal{S}$ on $X$. 

<div id="pointwise"></div>
<div class="definition">
  <strong>Definition (Pointwise Convergence).</strong>
  <br>
  Let $X$ be a set, and let $f_1, f_2, \dots$ be a sequence of functions with $f_i: X \rightarrow \mathbb{R}$. Let $f: X \rightarrow \mathbb{R}$ as well. We say that the sequence <i>converges pointwise</i> to $f$ if, for each $x \in X$:

  $$
  \underset{k \rightarrow \infty}{\lim} f_k(x) = f(x)
  $$ 
</div>

Another way to state the definition of pointwise convergence is that for each $x \in X$ and every $\epsilon > 0$, there exists $n \in \mathbb{Z}_+$ such that:

$$
\rvert f_k(x) - f(x) \rvert < \epsilon \hspace{5mm} \forall k \geq n
$$

Similarly, we have <i>uniform convergence</i>.

<div id="uniform"></div>
<div class="definition">
  <strong>Definition (Uniform Convergence).</strong>
  <br>
  Let $X$ be a set, and let $f_1, f_2, \dots$ be a sequence of functions with $f_i: X \rightarrow \mathbb{R}$. Let $f: X \rightarrow \mathbb{R}$ as well. We say that the sequence <i>converges uniformly</i> to $f$ if, for each $\epsilon > 0$, there exists $n \in \mathbb{Z}_+$ such that:

  $$
  \rvert f_k(x) - f(x) \rvert < \epsilon  \hspace{5mm} \forall k \geq n \text{ and } \forall x \in X
  $$ 
</div>

The difference between the definitions means that uniform convergence implies pointwise convergence (on the same set). The converse is not always true; however, it's <i>kind of</i> true for sequences of functions on measure spaces with finite total measure. 

<!-- #region egorov -->
<div id="egorov"></div>
<div class="theorem">
  <strong>Egorov's Theorem.</strong>
  <br>
  {% tabs egorov %}
  {% tab egorov statement %}
  Let $(X, \mathcal{S}, \mu)$ be a measure space with finite total measure (i.e $\mu(X) < \infty$). Let $f_1, f_2, \dots$ be a sequence of $\mathcal{S}$-measurable functions where $f_i: X \rightarrow \mathbb{R}$ that converges pointwise on $X$ to $f: X \rightarrow \mathbb{R}$. Then, for every $\epsilon > 0$, there exists $E \in \mathcal{S}$ such that $\mu(X \setminus E) < \epsilon$, and so $f_1, f_2, \dots$ converges uniformly to $f$ on $E$.
  {% endtab %}
  {% tab egorov proof %}
  Let $\epsilon > 0$, and fix some $$n \in \mathbb{Z}_+$$. We define, for $$m \in \mathbb{Z}_+$$:
  
  $$
  A_{m,n} = \bigcap_{k = m}^\infty \left\{ x \in X \rvert \frac{1}{n} > \rvert f_k(x) - f(x) \rvert \right\} 
  $$
  
  Note that each $A_{m,n} \in \mathcal{S}$ since $f_k - f$ is an $\mathcal{S}$-measurable function by the properties of measurable functions covered <a href="#limit-function">here</a> and <a href="#measurable-function">here</a>. Furthermore, we see that $A_{1,n} \subseteq A_{2,n} \subseteq \dots$ is an increasing sequence of sets (since we are doing fewer and fewer intersections as $m$ increases). Thus, by the <a href="#pointwise">definition</a> of pointwise convergence, it must be the case that:

  $$
  \bigcup_{m = 1}^\infty A_{m,n} = X
  $$

  To see why, set $\epsilon = \frac{1}{n}$. By the definition, there must be some $m \in \mathbb{Z}_+$ such that $\rvert f_k(x) - f(x) \rvert < \epsilon$ for all $k \geq m$ for all $x \in X$. Since the above holds for $x \in X$, taking the union of the sets of $x \in X$ satisfying the condition over all such choices of $m$ will yield the original set $X$.

  By Eq. \eqref{eq:increase-decrease}, we have that:

  $$
  \underset{m \rightarrow \infty}{\lim} \mu(A_{m,n}) = \mu(X)
  $$

  This implies that there must be an index $$m_n \in \mathbb{Z}_+$$ such that the lefthand and righthand sides of the prior equation are arbitrarily close. That is:

  $$
  \mu(X) - \mu(A_{m_n, n}) < \frac{\epsilon}{2^n}
  $$

  We will define:

  $$
  E = \bigcap_{n=1}^\infty A_{m_n, n}
  $$

  It follows that:

  $$
  \begin{aligned}
  \mu(X \setminus E)
  &= \mu\left(X \setminus \bigcap_{n = 1}^\infty A_{m_n, n} \right) \\
  &= \mu \left(\bigcup_{n = 1}^\infty \left(X \setminus A_{m_n, n} \right) \right)  & \left(\text{DeMorgan's Laws}\right) \\
  &\leq \sum_{n = 1}^\infty \mu\left(X \setminus A_{m_n, n}\right) & \left(\text{countable subadditivity}\right) \\
  &< \epsilon  & \left(  \mu(X) - \mu(A_{m_n, n}) < \frac{\epsilon}{2^n} \right)
  \end{aligned}
  $$

  This shows that the measure of the set $X \setminus E$ has measure going to zero. Now we need to show that $f_1, f_2, \dots$ converges uniformly to $f$ on $E$.

  Let $\epsilon' > 0$, and choose $$n \in \mathbb{Z}_+$$ such that $\frac{n} < \epsilon'$. By construction, $E \subseteq A_{m_n, n}$, and thus, for any $x \in E$:

  $$
  \rvert f_k(x) - f(x) \rvert < \frac{1}{n} < \epsilon'
  $$

  for all $k \geq m_n$, which concludes the proof.
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endregion -->

By <i>kind of true</i>, we mean that $f_1, f_2, \dots$ converges uniformly to $f$ except on sets of arbitrarily small measure. 

