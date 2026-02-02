---
layout: distill
title:  Measure Theory
description: A Primer
date: 2025-04-30
tabs: true
tags: theory likelihood
toc:
  - name: Building Blocks
    subsections:
        - name: Sets
        - name: Functions
  - name: Measures
    subsections:
        - name: Characteristics
  - name: Probability
    subsections:
        - name: Random Variables
        - name: Stochastic Processes
  - name: Integration
    subsections:
        - name: Non-Negative Functions
        - name: Real-Valued Functions
        - name: Radon-Nikodym
        - name: Integral Properties
  - name: Expectation
  - name: Miscellaneous
    subsections:
        - name: Vector Spaces
        - name: Helpful Definitions
        - name: Assorted Results
bibliography: measure.bib
---

My work has become much more technical that I am used to, so I thought it would be good to take some notes on basic measure and probability theory in anticipation of working through several theoretical papers. A lot of the definitions below come from Wikipedia, Durrett<d-cite key=durrett2019></d-cite>, and Axler<d-cite key=axler2025></d-cite>, but I've added some more intuitive ways of thinking about these concepts that I've come up with or collected from some really helpful sites I've found.

Note: Not all of the proofs are finished/included. I am hoping to find the time to return to this post and complete them.

---

## Building Blocks
Measure theory, in my mind, is just about sets, mappings, and ways to describe them. To formalize these ideas, however, we need to define some basic building blocks.

### Sets

We'll begin with a $\sigma$-field. This is simply a collection of subsets of some other set. 

<div id="sigma-field"></div>
<div class="definition">
  <strong>Definition ($\sigma$-Field).</strong>
  <br>
  Let $X$ be a set, and let $\mathcal{P}(X)$ denote its power set (the set of all possible subsets of $X$). A $\sigma$<i>-field</i> is any subset $\mathcal{S} \subseteq \mathcal{P}(X)$ satisfying the following:
  <ol>
    <li> $X \in \mathcal{S}$ </li>
    <li> Closed Under Complementation: $A \in \mathcal{S} \implies A^c \in \mathcal{S}$ </li>
    <li> Closed Under Countable Unions: $A_1, A_2, \dots \in \mathcal{S} \implies A = \bigcup_{i = 1}^\infty A_i \in \mathcal{S}$ </li>
  </ol>
</div>

<aside><p>For the technical reader, a $\sigma$-field, also called a $\sigma$-<i>algebra</i>, is a generalization of the <i>algebra</i>, which has the same definition except must be closed under <i>finite</i> unions.</p></aside>

An example may make it a bit more concrete in one's mind...

<div class="example">
  <strong>Example ($\sigma$-Field).</strong>
  <br>
  Let $X = \{ 1, 2, 3\}$. One $\sigma$-field on $X$ is $\mathcal{S} = \{ \emptyset, \{ 1 \}, \{2, 3 \}, \{1, 2, 3 \} \}$. 
  <br>
  It is easy to see that the first two properties are satisfied. The entirety of $X$ is in $\mathcal{S}$ by construction. Any element in $\mathcal{S}$, its complement is also in $\mathcal{S}$ (the complement of $\emptyset$ is $X$, and the complement of $\{ 1\}$ is $\{ 2, 3 \}$). 
  <br>
  Taking any countable union of elements of $\mathcal{S}$ also yields an element of $\mathcal{S}$ (the union of any $A \in \mathcal{S}$ and the empty set is the set itself; the complement of the entirety of $X$ and any $A \in \mathcal{S}$ will just be $X$; and the union of $\{ 1\}$ and $\{ 2, 3\}$ is $X$). This satisfies the third property, which completes the proof that $\mathcal{S}$ is a $\sigma$-field on $X$.
</div>

Notice that the first and second properties in the above definition imply that $\emptyset \in \mathcal{S}$ as well. The properties also imply that a $\sigma$-field must be closed under countable intersection. That is, $\cap_{i = 1}^\infty A_i \in \mathcal{S}$ for some sequence of $A_1, A_2, \dots \in \mathcal{S}$. 

A $\sigma$-field is a generalization of the concept of an <i>algebra</i> (also called a <i>field</i>). 

<div id="algebra"></div>
<div class="definition">
  <strong>Definition (Algebra).</strong>
  <br>
  Let $X$ be a set. A collection, $\mathcal{A}$, of subsets of $X$ is an <i>algebra</i> (or <i>field</i>) if the following are satisfied:
  <ol>
  <li>$X \in \mathcal{A}$</li>
  <li>Closed Under Complementation: $A \in \mathcal{A} \implies A^c \in \mathcal{A}$</li>
  <li>Closed Under Finite Unions: $A, B \in \mathcal{A} \implies A \cup B \in \mathcal{A}$</li>
  </ol>
</div>

Now we can define <i>measurable spaces</i>!

<div id="measurable-space"></div>
<div class="definition">
  <strong>Definition (Measurable Space).</strong>
  <br>
  Let $X$ be some set, and let $\mathcal{S}$ be a $\sigma$-field on $X$. The tuple $(X, \mathcal{S})$ is called a <i>measurable space</i>, and any element of $\mathcal{S}$ is called an <i>$\mathcal{S}$-measurable set</i>.
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

<div id="borel-set"></div>
<div class="definition">
  <strong>Definition (Borel Set).</strong>
  <br>
  A <i>Borel set</i> on a topological space, $X$, denoted by $\mathcal{B}(X)$, is any subset of $X$ that can be constructed from open sets on that space in $X$ via countable unions, countable intersections, and set differences. 
  <br>
  The <i>Borel $\sigma$-field</i> (or <i>Borel algebra</i>) is the collection of all Borel sets on a space. 
</div>

This definition is a bit tricky to develop intuition for. The Borel $\sigma$-field is just the collection of <i>all possible open sets</i> in a given space, $X$. 

An important Borel $\sigma$-field that will come up again when we discuss measures and probability is the Borel $\sigma$-field on the real line. Several examples follow from our definition:

<ol>
  <li>Any <i>closed</i> subset of $\mathbb{R}$ is a Borel set because $\sigma$-fields are closed under complementation.</li>
  <li>Any <i>countable</i> subset of $\mathbb{R}$ is a Borel set because $\sigma$-fields are closed under countable unions, and a single point is a closed subset of $\mathbb{R}$.</li>
  <li>Any half-open interval is a Borel set because $\sigma$-fields are closed under countable intersections.</li>
</ol>

Borel sets on $\mathbb{R}$ can also be extended to $[-\infty, \infty]$. 

Along with the Borel set and the $\sigma$-field is the <i>semialgebra</i>.

<div id="semialgebra"></div>
<div class="definition">
  <strong>Definition (Semialgebra).</strong>
  <br>
  Let $\mathcal{S}$ be a collection of sets. $\mathcal{S}$ is called a <i>semialgebra</i> if it satisfies the following properties:
  <ol>
  <li>Closed Under Intersection: $S,T \in \mathcal{S} \implies S \cap T \in \mathcal{S}$</li>
  <li>(Sort Of) Closed Under Complementation: $S \in \mathcal{S} \implies S^c$ is a finite disjoint union of $T \in \mathcal{S}$</li>
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

Furthermore, in the special case that $X \subseteq \mathbb{R}$ and $\mathcal{S}$ is the set of Borel subsets of $\mathbb{R}$ that are contained in $X$, then a function $f: X \rightarrow \mathbb{R}$ is called <i>Borel measurable</i> if $f^{-1}(B)$ is a Borel set for all Borel sets $B \subseteq \mathbb{R}$. It can be shown that any <i>continuous</i> or <i>increasing</i> function $f: X \rightarrow \mathbb{R}$ where $X$ is a Borel subset of $\mathbb{R}$ is Borel measurable.

---

## Measures
We have finally come to the star of our discussion: the measure. A measure is a function that assigns a "size" to sets (it is similar to the idea of length for intervals or area for two dimensional regions).

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

With this definition, we define a <i>measure space</i>, which is the tuple $(X, \mathcal{S}, \mu)$. For measure space $(X, \mathcal{S}, \mu)$ and $A, B \in \mathcal{S}$ such that $A \subseteq B$, we have that $\mu(A) \leq \mu(B)$ and $\mu(B \setminus A) = \mu(B) - \mu(A)$ (assuming that $\mu(A)$ is finite). We also have the additional property of <i>countable subadditivity</i>, which is basically a generalization of Boole's inequality:

$$
\mu\left(\bigcup_{i = 1}^\infty A_i \right) \leq \sum_{i = 1}^\infty \mu(A_i)
$$

for any sequence of sets $A_1, A_2, \dots \in \mathcal{S}$. Measures also satisfy $\mu(A \cup B) = \mu(A) + \mu(B) - \mu(A \cap B)$ (assuming that $\mu(A \cap B)$ is finite).

If we have two $\sigma$-finite (see <a href="#sigma-finite">below</a>) measure spaces, $(X, \mathcal{S}, \mu_1)$ and $(Y, \mathcal{T}, \mu_2)$, we can define two addition sets:

$$
\begin{aligned}
\Omega &= X \times Y = \{ (x, y): x \in X, y \in Y\} \\
\mathcal{U} &= \{ S \times T: S \in \mathcal{S}, T \in \mathcal{T}\}
\end{aligned} 
$$

Sets $U \in \mathcal{U}$ are <i>rectangles</i>. Let $\mathcal{F} = \mathcal{S} \times \mathcal{T}$ be the $\sigma$-filed generated by $\mathcal{U}$. The unique measure $\mu = \mu_1 \times \mu_2$ on $\mathcal{F}$ defined as $\mu(S \times T) = \mu_1(S) \mu_2(T)$ is called a <i>product measure</i>. This result can be extended to finitely many $\sigma$-finite measurable spaces. 

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
  \mu\left(\{ x : \rvert f_n(x) - f(x) \rvert > e \right) \rightarrow 0 \hspace{5mm} \text{ as } n \rightarrow \infty
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


---

## Probability
With our building blocks in place, we can move on to probability theory. We'll start with a fundamental definition: the <i>probability space</i>, which is just a special measure space!

<div id="probability-space"></div>
<div class="definition">
  <strong>Definition (Probability Space).</strong>
  <br>
  A <i>probability space</i> is a measure space $(\Omega, \mathcal{F}, P)$ such that $P(\Omega) = 1$. We have the following conventions:

  <ol>
    <li>The set $\Omega$ is called the <i>sample space</i>. It is the set of all possible outcomes.</li>
    <li>The $\sigma$-field $\mathcal{F}$ over $\Omega$ is called the <i>event space</i>. It is a measurable set of subsets of the sample space.</li>
    <li>The measure $P$ is called the <i>probability measure</i>. It assigns a value (in $[0, 1]$) to give an event to give a sense of that event's likelihood of occurring.</li>
  </ol>
</div>

### Random Variables
Using the above, we can define random variables and vectors in a rigorous way. Note that the following can be generalized to the extended real line (i.e. $\mathbb{R} \cup \{-\infty, \infty \}$). 

<div id="random-variable"></div>
<div class="definition">
  <strong>Definition (Random Variable/Vector).</strong>
  <br>
  Let $(\Omega, \mathcal{F}, P)$ be a probability space, and let $(\mathbb{R}^d, \mathcal{B}^d)$ be our measurable space, where $\mathcal{B}^d(\mathbb{R}^d)$ is the $\sigma$-field on $\mathbb{R}^d$. A measurable map $X: \Omega \rightarrow \mathbb{R}^d$ is called a <i>random vector</i> if $d > 1$ and a <i>random variable</i> otherwise.
</div>

Random variables map each element in the sample space to an element in $H$, which is the set of all possible values the variable can take on. Naturally, we need the pre-image of all elements in $\mathcal{H}$ to be in $\mathcal{F}$. When we refer to a random variable being measurable with respect to some $\mathcal{F}'$ (a sub-$\sigma$-field of $\mathcal{F}$), we mean that it is $(\mathcal{F}', \mathcal{B}(\mathbb{R}))$-measurable. 

<div id="theorem-1-3-567"></div>
<div class="theorem">
<strong>Theorems 1.3.5, 1.3.6, 1.3.7.<d-cite key=durett2019></d-cite></strong>
{% tabs theorems-5-6-7 %}
{% tab theorems-5-6-7 statement %}
Let $X_1, X_2, \dots$ be random variables, and let $f: (\mathbb{R}^n, \mathcal{B}^n) \rightarrow (\mathbb{R}, \mathcal{B})$ be a measurable function. Then the following are also random variables:
<ul>
<li>$f(X_1, \dots, X_n)$</li>
<li>$X_1 + \dots + X_n$</li>
<li>$\underset{n}{\inf} X_n$</li>
<li>$\underset{n}{\sup} X_n$</li>
<li>$\underset{n}{\lim \inf} X_n$</li>
<li>$\underset{n}{\lim \sup} X_n$</li>
</ul>
{% endtab %}
{% tab theorems-5-6-7 proof %}
Proof to be completed.
{% endtab %}
{% endtabs %}
</div>

The <i>distribution</i> of a random variable can also be defined from a measure theoretic perspective.

<div id="distribution"></div>
<div class="definition">
  <strong>Definition (Distribution).</strong>
  <br>
  Let $(\Omega, \mathcal{F}, P)$ be a probability space. The probability measure induced by a random variable, $X$, defined on this space is called its <i>distribution</i> and is defined as: $\mu(A) = \mathbb{P}_P(X \in A)$ for all Borel sets $A \in \mathcal{B}$. 
  <br>
  The <i>distribution function</i>, $F(x) = \mathbb{P}_P(X \leq x)$, describes the distribution of $X$ and satisfies the following properties:
  <ol>
  <li>$F$ is non-decreasing</li>
  <li>$\underset{x \rightarrow \infty}{\lim} F(x) = 1$ and $\underset{x \rightarrow -\infty}{\lim} F(x) = 0$</li>
  <li>$\underset{y \downarrow x}{\lim} F(x) = 0$ (right continuous)</li>
  <li>$F(x-) = \mathbb{P}_P(X < x)$ </li>
  <li>$\mathbb{P}_P(X = x) = F(x) - F(x-)$</li>
  </ol>
  where $F(x-) = \underset{y \uparrow x}{\lim} F(y)$.
</div>

Durrett provides the best intuition for the distribution of a random variable: "In words, we pull $A \in \mathbb{B}$ back to $X^{-1}(A) \in \mathcal{F}$ and then take $P$ of that set"[^fn-durrett]. 

It's important to remember that two different random variables can induce the same distribution. In this case, we say that the random variables (denote them by $X$ and $Y$) are <i>equal in distribution</i>, which we denote with $X \overset{d}{=} Y$. 

A distribution function with the form $F(x) = \int_{-\infty}^x f(y) dy$ can also be described by its <i>density function</i>, $f$, satisfying:

$$
\mathbb{P}(X = x) = \underset{e \rightarrow 0}{\lim} \int_{x - e}^{x + e} f(y) dy = 0
$$

In this case, we say that $F$ is <i>absolutely continuous</i>. Integrating the density function over the entire sample space/support will equal $1$, and the density function will always be non-negative.

Similarly, we can define a <i>discrete</i> distribution function (i.e. an induced probability measure) as one in which there exists a countable set $S$ such that $P(S^c) = 0$. 

<div class="example">
  <strong>Example (Discrete Distribution).</strong>
  <br>
  Suppose a random variable on $(\mathbb{R}, \mathcal{B})$ induces distribution $F$ such that $F(x) = 1$ for $x \geq 0$ and $F(x) = 0$ for $x < 0$. This measure is discrete, and we call it a <i>point mass</i> at $0$.
</div>

Now, since random variables are just measurable functions, we can use its mapping to define special $\sigma$-fields. 

<div id="generated-sigma-field"></div>
<div class="definition">
  <strong>Definition (Generated $\sigma$-Field - Function).</strong>
  <br>
  Let $(\Omega, \mathcal{A})$ and $(S, \mathcal{B})$ be measurable spaces and let $f: \Omega \rightarrow S$ be a measurable function, then the $\sigma$-field $f^{-1}(\mathcal{A})$ is called the <i>$\sigma$-field generated by $f$</i>. We denote it by $\sigma(f)$. 
  <br>
  More intuitively, $\sigma(f)$ is the smallest $\sigma$-field on which $f$ is measurable. 
</div>

Put intuitively, the $\sigma$-field generated by random variable $X$ is the collection of all possible subsets of the set of values $X$ can take on such that the probability of the event that $X$ takes on that value can be determined (i.e. is measurable). 

We can also define $\sigma$-fields generated by arbitary subsets. This is the smallest $\sigma$-field containing a given collection of subsets.

<div id="generated-sigma-field2"></div>
<div class="definition">
  <strong>Definition (Generated $\sigma$-Field - Family).</strong>
  <br>
  Let $X$ be a set, and let $A$ be a collection of subsets of $X$. The $\sigma$-field generated by $A$, notated as $\sigma(A)$, is the collection of all subsets of $X$ that can be constructed from elements of $A$ under countable unions, intersections, and complementations. 
  <br>
  If $A = \emptyset$, then $\sigma(A) = \{ \emptyset, X\}$. If $A$ is a single event, then $\sigma(A) = \{ \emptyset, A, A^c, X \}$. 
</div>

<div class="example">
  <strong>Example (Generated $\sigma$-Field).</strong>
  <br>
  Consider the Borel $\sigma$-field. This $\sigma$-field can be generated by any of the following sets:
  <ul>
    <li>$\{ (a, b) \rvert a, b \in \mathbb{R}, a < b \}$ or $\{ (a, \infty) \rvert a \in \mathbb{R} \}$</li>
    <li>$\{ (a, b] \rvert a, b \in \mathbb{R}, a < b \}$ or $\{ (-\infty, a) \rvert a \in \mathbb{R} \}$</li>
    <li>$\{ [a, b) \rvert a, b \in \mathbb{R}, a < b \}$ or $\{ [a, \infty) \rvert a \in \mathbb{R} \}$</li>
    <li>$\{ [a, b] \rvert a, b \in \mathbb{R}, a < b \}$ or $\{ (-\infty, a] \rvert a \in \mathbb{R} \}$</li>
  </ul>
</div>


### Stochastic Processes 
We can also think of having many random variables, each associated with some step in a sequence (perhaps time or space). We call this a <i>stochastic process</i>.

<div id="stochastic-process"></div>
<div class="definition">
  <strong>Definition (Stochastic Process).</strong>
  <br>
  Let $(\Omega, \mathcal{F}, P)$ be a probability space, and let $(S, \Sigma)$ be a measurable space. Let $T$ be an index set. We call a collection of random variables, $\{ X(t, \omega): t \in T, \omega \in \Omega \}$, taking on values in $S$ a <i>stochastic process</i> on $(\Omega, \mathcal{F}, P)$ with state space $(S, \Sigma)$. In other words, a stochastic process is a random function $X: T \times \Omega \rightarrow S$. 
  <br>
  We will sometimes alternatively denote a stochastic process with $\{ X(t): t \in T\}$ and $(X_t)_{t \in T}$. 
</div>

Stochastic processes can be characterized by their <i>continuity</i> (or lack thereof).

<div id="sp-continuity"></div>
<div class="definition">
  <strong>Definition (Right-Continuous).</strong>
  <br>
  Let $X: T \times \Omega \rightarrow S$ be a stochastic process on probability space $(\Omega, \mathcal{F}, P)$ with (measurable) state space $(S, \Sigma)$. If, for all $\omega \in \Omega$, there exists $\epsilon > 0$ such that $X_s(\omega) = X_t(\omega)$ for all $s, t$ such that $t \leq s \leq t + \epsilon$, then we call $(X_t)_{t \in T}$ a <i>right-continuous</i> stochastic process.
  <br>
  Note that this definition requires the index set to be over the non-negative reals.
</div>

---

## Integration
Before we can look at random variables any further, we need to discuss a <i>very</i> important concept in mathematics. In the following, we will restrict our discussion to $\mathbb{R}$, but the definitions can easily be generalized to higher dimensions by exchanging lengths for volumes via Cartesian products. 

First, we define a special indicator function that got a fancy name (not sure why).

<div id="characteristic-function"></div>
<div class="definition">
  <strong>Definition (Characteristic Function).</strong>
  <br>
  Let $E \subseteq X$. We define the <i>characteristic function</i> of $E$ as the function $\chi_E: X \rightarrow \mathbb{R}$ defined by:

  $$
  \chi_E(x) = 
  \begin{cases}
    1 & \text{if } x \in E \\
    0 & \text{if } x \notin E
  \end{cases}
  $$
</div>

Though not very useful for our discussion, we'll define the <i>outer measure</i> of a set $A \subseteq \mathbb{R}$. The outer measure formalizes the size of a set by using the lengths of open intervals. 

<div id="outer-measure"></div>
<div class="definition">
  <strong>Definition (Outer Measure).</strong>
  <br>
  Let $I$ be some open interval on the real line, and let $\ell(I)$ denote the <i>length</i> of $I$ defined as:

  $$
  \ell(I) = 
  \begin{cases}
    b - a & \text{if } I = (a, b) \text{ for some } a, b \in \mathbb{R} \text{ such that } a < b \\
    0 & \text{if } I = \emptyset \\
    \infty & \text{otherwise }
  \end{cases}
  \nonumber
  $$

  We define the <i>outer measure</i> of $A \subseteq \mathbb{R}$ as:

  $$
  \rvert A \rvert = \inf\left\{ \sum_{i = 1}^\infty \ell(I_i) \bigg\rvert I_1, I_2, \dots \text{ are open intervals such that } A \subseteq \bigcup_{i = 1}^\infty I_i \right\}
  \nonumber
  $$

  The outer measure satisfies:

  <ol>
    <li>$\rvert A \rvert \leq \rvert B \rvert$ for $A \subseteq B \subseteq \mathbb{R}$</li>
    <li>$\rvert \{ t + a \rvert a \in A \}\rvert = \rvert A \rvert$ for $t \in \mathbb{R}$ and $A \subseteq \mathbb{R}$</li>
    <li>$\rvert \cup_{i = 1]^\infty} A_i \rvert \leq \sum_{i = 1}^\infty \rvert A_i \rvert$ for $A_1, A_2, \dots \subseteq \mathbb{R}$</li>
  </ol>
</div>

In words, the outer measure of a set is the smallest total length of some sequence of open intervals of $\mathbb{R}$ that, together, contain $A$. Finite sets have outer measure $0$ because we can make our open intervals arbitrarily "short" (i.e. force them to have length approaching $0$). By similar reasoning, any countable subset of $\mathbb{R}$ also has outer measure $0$.

It's important to remember that the outer measure is not a true measure in the sense that we defined. However, the outer measure allows us to define a special (and true) measure called the <i>Lebesgue measure</i>. 

<div id="lebesgue-measure"></div>
<div class="definition">
  <strong>Definition (Lebesgue Measure).</strong>
  <br>
  Let $\mathcal{B}$ be the $\sigma$-field of Borel subsets of $\mathbb{R}$, and let $(\mathbb{R}, \mathcal{B})$ be our measurable space. The <i>Lebesgue measure</i> on $(\mathbb{R}, \mathcal{B})$ is the measure such that $\mu(B) = \rvert B \rvert$ for any Borel set $B \in \mathcal{B}$. 
</div>

In words, the outer measure becomes a true measure if we restrict ourselves to only Borel sets. The Lebesgue measure leads to a refinement of the idea of a measurable set. A set $A \subseteq \mathbb{R}$ is called <i>Lebesgue measurable</i> if it is really "close" to being a Borel set. Put formally, $A$ is Lebesgue measurable if there exists a Borel set $B \subseteq A$ such that $\rvert A \setminus B \rvert = 0$. There are also many equivalent definitions (see pg. 52 of Axler (2025)). 

Note that sometimes the definition of the Lebesgue measure is <span class="popup" onclick="PopupFunc('pop2')">altered<span class="popuptext" id="pop2">The change is limited to the function's domain (Borel vs. Lebesgue measurable sets).</span></span> to mean the measure on $(\mathbb{R}, \mathcal{L})$ where $\mathcal{L}$ is the $\sigma$-field of Lebesgue measurable subsets of $\mathbb{R}$. 

A function $f: A \rightarrow \mathbb{R}$ for $A \subseteq \mathbb{R}$ is <i>Lebesgue measurable</i> if $f^{-1}(B)$ is a Lebesgue measurable set for every Borel set $B \subseteq \mathbb{R}$. 

A lot of things in probability depend upon integration. For example, expectation, variance, cumulative probability, and many more things can all be stated as some type of integral. Thus, it's important we have a solid understanding of the integral.


### Non-Negative Functions
We start with the integral of the characteristic function:

$$
\int \chi_E d\mu = \mu(E) \hspace{5mm} \forall E \in \mathcal{S}
$$

Recall that a simple function is any function that takes on finitely many values. Any piecewise function with finitely many pieces is simple. We can use the integral of the characteristic function to that of simple functions by taking a linear combination. 

Let $(X, \mathcal{S}, \mu)$ be a measure space, let $A_1, \dots, A_n$ be disjoint set in $\mathcal{S}$, and let $c_1, \dots, c_n \in [0, \infty]$. Then:

$$
\int \left(\sum_{i = 1}^n c_i \chi_{A_i} \right) d\mu = \sum_{i = 1}^n c_i \mu(A_i)
$$

With these definitions in mind, we can define the integral of any non-negative function. 

<div id="integral-nonnegative"></div>
<div class="definition">
  <strong>Definition (Integral of a Non-Negative Function).</strong>
  <br>
  Let $(X, \mathcal{S}, \mu)$ be a measure space, and let $f: X \rightarrow [0, \infty]$ be an $\mathcal{S}$-measurable function. Its integral with respect to $\mu$ is defined as:

  $$
  \int f d\mu = \sup \left\{ \sum_{i = 1}^n c_i \mu(A_i) \bigg\rvert \text{ disjoint } A_1, \dots, A_n \in \mathcal{S}; \hspace{2mm} c_1, \dots, c_n \in [0, \infty); \hspace{2mm} f(x) \geq \sum_{i = 1}^n c_i \chi_{A_{i}}(x) \text{ for all } x \in X \right\}
  $$
</div>

### Real-Valued Functions
We begin with a definition.

<div id="f-plus-minus"></div>
<div class="definition">
  <strong>Definition ($f^+$ and $f^-$).</strong>
  <br>
  Let $f: X \rightarrow [-\infty, \infty]$ be a function. We have that $f = f^+ - f^-$ and $\rvert f \rvert = f^+ + f^-$ for piecewise functions:

  $$
  f^+ = \begin{cases}
  f(x) & \text{ if } f(x) \geq 0 \\
  0    & \text{ if } f(x) < 0
  \end{cases}
  \hspace{10mm} \text{and} \hspace{10mm}
  f^- = \begin{cases}
  0 & \text{ if } f(x) \geq 0 \\
  f(x) & \text{ if } -f(x) < 0
  \end{cases}
  \nonumber
  $$
</div>

Notice that if $f(x) \geq 0$, then $f^+(x) \geq 0$ and $f^-(x) = 0$. Alternatively, if $f(x) < 0$, then $f^+(x) = 0$ and $f^-(x) = -f(x) > 0$. Thus, $f^+$ and $f^-$ are both non-negative functions. This allows us extend the definition of the integral to real-valued functions.

<div id="integral"></div>
<div class="definition">
  <strong>Definition (Integral of a Real-Valued Function).</strong>
  <br>
  Let $(X, \mathcal{S}, \mu)$ be a measure space, and let $f: X \rightarrow [-\infty, \infty]$ be an $\mathcal{S}$-measurable function such that $\int f^+ d\mu < \infty$, $\int f^- d\mu < \infty$, or both. The <i>integral</i> of $f$ with respect to $\mu$ is defined as:

  $$
  \int f d\mu = \int f^+ d\mu - \int f^- d\mu
  \nonumber
  $$

  The integral is homogeneous (i.e. $\int c f d\mu = c \int f d\mu$ for any $c \in \mathbb{R}$) and additive (i.e. $\int (f + g) d\mu = \int f d\mu + \int g d\mu$ for $\mathcal{S}$-measurable $f$ and $g$ satisfying $\int \rvert f \rvert d\mu < \infty$ and $\int \rvert g \rvert d \mu < \infty$).
</div>

If we have $(\Omega, \mathcal{F}, \mu) = (\mathbb{R}^d, \mathcal{B}^d, \lambda)$, then we denote $\int f d\lambda$ with $\int f(x) dx$, and if $(\Omega, \mathcal{F}, \mu) = (\mathbb{R}, \mathcal{B}, \lambda)$ and we have some interval $E = [a, b]$, we write $\int_a^b f(x) dx$ instead of $\int_E f d\lambda$.

Integration can be restricted to a subset of the domain of a function. That is, for $E \in \mathcal{S}$:

$$
\int_E f d\mu = \int f \chi_E d \mu
$$

It can also be restricted to an interval of the extended real line. First, we call a bounded function $f: [a, b] \rightarrow \mathbb{R}$ <i>Riemann integrable</i> if the set of points in $[a, b]$ at which $f$ is not continuous has length $0$. If we have Lebesgue measure on $\mathbb{R}$, $\lambda$, and $f: (a, b) \rightarrow \mathbb{R}$ is a Lebesgue measurable function, then for $-\infty \leq a < b \leq \infty$ we let $\int_a^b f(x) dx = \int_{(a,b)} f d\lambda$. 

### Radon-Nikodym
Two different measures can be related via the <a href="https://en.wikipedia.org/w/index.php?title=Radon%E2%80%93Nikodym_theorem&oldid=1288156682">Radon-Nikodym Theorem</a>, which states that (under certain conditions), there exists a function such that one measure is equivalent to the integral of the function with respect to a second measure.

<div id="rn-theorem"></div>
<div class="theorem">
  <strong>Radon-Nikodym Theorem.</strong>
  {% tabs radon-nikodym %}
  {% tab radon-nikodym statement %}
  Let $(X, \mathcal{S})$ be a measurable space, and let $\mu$ and $\nu$ denote two $\sigma$-finite measures on this space such that $\nu << \mu$ ($\nu$ is absolutely continuous with respect to $\mu$). Then there existgs a $\mathcal{S}$-measurable function, $f: X \rightarrow [0, \infty)$ such that, for any measurable $A \subset \mathcal{S}$:

  $$
  \nu(A) = \int_A f d\mu
  $$
  {% endtab %}
  {% tab radon-nikodym proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

A fun fact is that $f$ is unique up to some set of measure $0$ with respect to $0$. That is, for any other $g$ that satisfies the definition, $f(x) = g(x)$ for all $x \in X$ except some $x \in X' \subset X$ such that $\mu(X') = 0$. Such a function, $f$, is called the <i>Radon-Nikodym derivative</i> and can be denoted by $\frac{d \nu}{d \mu}$.

### Integral Properties
Here we list and prove several properties of integrals that are ubiquitous in theoretical statistics.

<div id="jensen"></div>
<div class="theorem">
  <strong>Jensen's Inequality.</strong>
  {% tabs jensen %}
  {% tab jensen statement %}
  Let $\phi$ be a convex function (i.e. $\lambda \phi(x) + (1- \lambda)\phi(y) \geq \phi(\lambda x + (1-\lambda)y)$ for all $\lambda \in (0, 1)$, $x,y \in \mathbb{R}$). Let $\mu$ be a probability measure, and let $f$ and $\phi(f)$ be integrable. Jensen's inequality states:

  $$
  \phi\left(\int f d\mu \right) \leq \int \phi(f)d\mu
  $$
  {% endtab %}
  {% tab jensen proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

<div id="holder"></div>
<div class="theorem">
  <strong>Hölder's Inequality.</strong>
  {% tabs holder %}
  {% tab holder statement %}
  Let $\mu$ be a probability measure, and let $p, q \in (1, \infty)$ such that $\frac{1}{p} + \frac{1}{q} = 1$. Hölder's inequality states:

  $$
  \int \rvert fg \rvert d\mu \leq \rvert \rvert f \rvert \rvert_p \rvert \rvert g \rvert \rvert_1
  $$

  where $\rvert \rvert f \rvert \rvert_p = (\int \rvert f \rvert^p d\mu)^{\frac{1}{p}}$ for $1 \leq p < \infty$. 
  {% endtab %}
  {% tab holder proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

If $p = q = 2$, the above is called the <i>Cauchy-Schwarz inequality</i>.

<div id="bounded-convergence"></div>
<div class="theorem">
  <strong>Bounded Convergence Theorem.</strong>
  {% tabs bounded-conv %}
  {% tab bounded-conv statement %}
  Let $E$ be a set of finite measure (i.e. $\mu(E) < \infty$), and let $\{ f_n \}$ be a sequence of functions that vanish on $E^c$, are uniformly pointwise bounded (i.e. $\rvert f_n(x) \rvert \leq M$), and $f_n \rightarrow f$ in measure. Then:

  $$
  \int f d\mu = \underset{n \rightarrow \infty}{\lim} \int f_n d\mu
  $$
  {% endtab %}
  {% tab bounded-conv proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

<div id="monotone-convergence"></div>
<div class="theorem">
  <strong>Monotone Convergence Theorem.</strong>
  {% tabs monotone-conv %}
  {% tab monotone-conv statement %}
  Let $(\Omega, \mathcal{F}, \mu)$ be a measure space, and let $X \in \mathcal{F}$ be a measurable set. Let $$\{ f_k \}_{k = 0}^\infty$$ be a pointwise non-decreasing sequence of $$(\mathcal{F}, \mathbb{B}(\bar{\mathbb{R}}_{\geq 0})$$-measurable, non-negative functions (i.e. $$0 \leq \dots \leq f_k(x) \leq f_{k+1}(x) \leq \dots \leq \infty$$ for every $k \geq 1$ and $x \in X$). Then the pointwise supremum, defined as the function:

  $$
  \underset{k}{\sup} f_k: x \rightarrow \underset{k}{\sup} f_k(x)
  $$

  is $(\mathcal{F}, \mathbb{B}(\bar{\mathbb{R}}_{\geq 0}))$-measurable and satisfies:

  $$
  \underset{k}{\sup} \int_X f_k d\mu = \int_X \underset{k}{\sup} f_k d\mu
  $$
  {% endtab %}
  {% tab monotone-conv proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

<div id="dominated-convergence"></div>
<div class="theorem">
  <strong>Dominated Convergence Theorem.</strong>
  {% tabs domin-conv %}
  {% tab domin-conv statement %}
  Let $(\Omega, \mathcal{F}, \mu)$ be a measure space, and let $$\{ f_k \}_{k \in T}$$ be a sequence of measurable functions (with index set $T$) on this space such that $$\underset{n \rightarrow \infty}{\lim} f_n(x) = f(x)$$ for some function $f$ for all $x \in \Omega$ (i.e. $$\{ f_k \}_{k \in T}$$ converges pointwise to $f$). Suppose that our sequence is <i>dominated</i> by some other integrable function, $g$; that is:

  $$
  \rvert f_n(x) \rvert \leq g(x) \hspace{5mm} \forall x \in \Omega, \hspace{2mm} \forall n \in T
  $$

  The Dominated Convergence Theorem states that $f_n$ and $f$ are both (Lebesgue) integrable and:

  $$
  \underset{n \rightarrow \infty}{\lim} \int_\Omega f_n d\mu = \int_\Omega \underset{n \rightarrow \infty}{\lim} f_n d\mu = \int_\Omega f d\mu
  $$
  {% endtab %}
  {% tab domin-conv proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

<div id="fatou"></div>
<div class="theorem">
  <strong>Fatou's Lemma.</strong>
  {% tabs fatou %}
  {% tab fatou statement %}
  Let $(\Omega, \mathcal{F}, \mu)$ be a measure space, and let $$\{ f_n: \Omega \rightarrow [0, \infty]\}$$ be a sequence of non-negative measurable functions. Then:

  $$
  \int_X \underset{n \rightarrow \infty}{\lim} \underset{m \geq n}{\inf} f_n d\mu \leq \underset{n \rightarrow \infty}{\lim} \underset{m \geq n}{\inf} \int_X f_n d\mu
  $$
  {% endtab %}
  {% tab fatou proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

<div id="fubini"></div>
<div class="theorem">
  <strong>Fubini's Theorem.</strong>
  {% tabs fubini %}
  {% tab fubini statement %}
  Let $(X, \mathcal{S}, \mu_1)$ and $(Y, \mathcal{T}, \mu_2)$ be $\sigma$-finite measurable spaces, and let $\mu = \mu_1 \times \mu_2$ (the product meeasure). If we have a function $f$ such that $f \geq 0$ or $\int \rvert f \rvert d \mu$, then:

  $$
  \int_X \int_Y f(x,y) \mu_2(dy) \mu_1(dx) = \int_{X \times Y} f d \mu = \int_Y \int_X f(x,y) \mu_1(dx) \mu_2(dy)
  $$
  {% endtab %}
  {% tab fubini proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

Fubini's Theorem tells us when it is okay to exchange the order of a double integral and to compute a double integral as an interated integral. 

---

## Expectation
For a random variable, $X$ on probability space $(\Omega, \mathcal{F}, P)$, how can we describe its central tendency (i.e. what values $X$ usually takes on)? We answer this question with the following definitions, which use ideas from integration (see <a href="#integration">later in this post</a>). 

<div id="expectation"></div>
<div class="definition">
  <strong>Definition (Expectation).</strong>
  <br>
  Let $X$ be a real-valued random variable on probability space $(\Omega, \mathcal{F}, P)$. Its <i>expectation</i> or <i>expected value</i> is defined as the following Lebesgue integral:

  $$
  \mathbb{E}[X] = \int_\Omega X dP
  $$
</div>

The expected value or expectation of a random variable is basically just integration with respect to the probability measure of the space the variable is defined on. It can be any real number or even $\infty$. Since it is just an integral, we can extend all of the results in the previous section to the expectation. The results are the same, just rewritten with $\mathbb{E}[X]$ instead of $\int_\Omega X dP$. 

We can also define the <i>conditional expectation</i> of a random variable with respect to a particular sub-$\sigma$-field. 

<div id="conditional-expectation"></div>
<div class="definition">
  <strong>Definition (Conditional Expectation).</strong>
  <br>
  Let $(\Omega, \mathcal{F}, P)$ be a probability space, let $X: \Omega \rightarrow \mathbb{R}^n$ be a real-valued random variable with finite expectation, and let $\mathcal{H} \subseteq \mathcal{F}$ be a sub-$\sigma$-field of $\mathcal{F}$. A <i>conditional expectation of $X$ given $\mathcal{H}$</i> is any $\mathcal{H}$-measurable function, $\mathbb{E}(X \rvert \mathcal{H}): \Omega \rightarrow \mathbb{R}^n$, satisfying:

  $$
  \int_H \mathbb{E}[X \rvert \mathcal{H}] dP = \int_H X dP \hspace{5mm} \forall H \in \mathcal{H}
  $$

  This function exists and is unique. 
</div>

<div id="conditional-expectation-ex"></div>
<div class="example">
  <strong>Definition (Conditional Expectation).</strong>
  <br>
  Suppose $X \in \mathcal{F}$. Then $\mathbb{E}_P[X \rvert \mathcal{F}] = X$ itself. In this scenario, we have perfect information. Since $X \in \mathcal{F}$, we have complete knowledge of whether it occurred or not when we are given all of $\mathcal{F}$. 
</div>

These definitions are a bit confusing, so let's parse them by coming at the topic from a different angle (see <a href="https://math.stackexchange.com/questions/375994/intuition-behind-conditional-expectation-when-sigma-algebra-isnt-generated-by-a">this post</a>).

Let's say we have a random variable $X$ on some probability space $(\Omega, \mathcal{F}, P)$. We don't know anything about it, so our best guess at its value would be some sort of weighted average over all of the possible values it could take on. These weights are determined by the probability measure, $P$, since a good guess should be closer to the more likely outcomes. 

Now, suppose we know some information about $X$'s outcome (i.e. we can answer some set of questions about $X$). We could formulate this as a collection of subsets of $\Omega$. For example, if we were rolling dice, the question "Is $X$ odd?" could be contained in the set $$\{1, 3, 5\}$$ or $$\{2, 4, 6\}$$. We could imagine outputting a different best guess depending upon what set of information we are given, which is basically what the conditional expectation does. 

In one way of thinking, $\mathbb{E}[X \rvert \mathcal{H}]$ is a random variable mapping from the possible values of $X$ to the best guesses. The condition $\int_H \mathbb{E}[X \rvert \mathcal{H}] dP = \int_H X dP$ for all $H \in \mathcal{H}$ can be thought of as enforcing the idea that, if we only are guessing values that are consistent with $H$, then our best guess using <i>only</i> the information in $H$ should be the same as the weighted average of $X$ itself over $H$. More concretely, if $$H = \{ 1, 3, 5\}$$ in our dice rolling example, having $\int_H \mathbb{E}[X \rvert \mathcal{H} dP = \int_H X dP$ implies that, given $H$, we can guess the average of $X$ perfectly. 

We can relate this measure theoretic definition with the more common ones learned in statistics courses. First, partition the sample space, $\Omega$, into disjoint sets $\Omega_1, \Omega_2, \dots$ such that $\mu(\Omega_i) > 0$ for all $i$. Let $\mathcal{F} = \sigma(\Omega_1, \Omega_2, \dots)$ be the $\sigma$-field generated by this collection of sets. For random variable $X$ defined on $(\Omega, \mathcal{F}, \mu)$, we have:

$$
\mathbb{E}_\mu[X \rvert \mathcal{F}] = \frac{\mathbb{E}_\mu[X \rvert \Omega_i]}{\mu(\Omega_i)} \hspace{5mm} \text{on } \Omega_i
$$

When we are given some information, $\Omega_i$, about which set in our partition $X$ can be found in, our best guess at $X$ becomes the average of $X$ over that set. 

In most probability courses, we also learn about conditional expectations with respect to some other random variable. In this case, we write $\mathbb{E}[X \rvert Y]$ to mean $\mathbb{E}[X \rvert \sigma(Y)]$. 


---

## Miscellaneous

### Vector Spaces

We now introduce the idea of <i>vector spaces</i>. We begin with the definition of a <i>field</i>.


<div id="field"></div>
<div class="definition">
  <strong>Definition (Field).</strong>
  <br>
  A <i>field</i> is a set $F$ with the operations of addition ($+$) and multiplication ($\cdot$), which satisfy for any $a, b, c \in F$:

  <ol>
    <li>Associativity: $a + (b + c) = (a + b) + c$ and $a \cdot (b \cdot c) = (a \cdot b) \cdot c$</li>
    <li>Communtatitvity: $a + b = b + a$ and $a \cdot b = b \cdot a$</li>
    <li>Identity: $\exists 0, 1 \in F$ such taht $a + 0 = a$ and $a \cdot 1 = a$</li>
    <li>Additive Inverses: $\forall a \in F$, $\exists -a \in F$ such taht $a + (-a) = 0$</li>
    <li>Multiplicative Inverses: $\forall a \in F$ such that $a \neq 0$, $\exists a^{-1} \in F$ such that $a \cdot a^{-1} = 1$</li>
    <li>Distributivity: $a \cdot (b + c) = (a \cdot b) + (a \cdot c)$</li>
  </ol>
</div>

A vector space is defined with respect to a field. In generality, it is a set of elements that satisfy some special properties in relation to some field. 

<div id="vector-space"></div>
<div class="definition">
  <strong>Definition (Vector Space).</strong>
  <br>
  Let $F$ be a field. A <i>vector space</i>, $V$, is some (non-empty) set with the operation of <i>vector addition</i> ($+$) and the function of <i>scalar multiplication</i>. Vector addition takes two vectors in $V$ and assigns them a sum, which is just a third vector in $V$. Scalar multiplication takes any a vector in $V$ and any scalar $a$ in $F$ and assigns it a product, which is another vector in $V$. This operation and function satisfy the following for any $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and any $a, b \in F$:

  <ol>
    <li>Associativity: $\mathbf{u} + (\mathbf{v} + \mathbf{w}) = (\mathbf{u} + \mathbf{v}) + \mathbf{w}$</li>
    <li>Commutativity: $\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$</li>
    <li>Vector Addition Identity: $\exists \mathbf{0} \in V$ such that $\mathbf{v} + \mathbf{0} = \mathbf{v}$ for all $\mathbf{v} \in V$</li>
    <li>Scalar Multiplication Identity: $1 \mathbf{v} = \mathbf{v}$ where $1$ is the multiplicative identity in $F$</li>
    <li>Inverses: $\exists -\mathbf{v} \in V$ for every $ \mathbf{v} \in V$ such that $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$</li>
    <li>Compatibility: $a(b \mathbf{v}) = (ab) \mathbf{v}$</li>
    <li>Distributivity: $a(\mathbf{u} + \mathbf{v}) = a \mathbf{u} + a \mathbf{v}$ and $(a + b) \mathbf{v} = a \mathbf{v} + b \mathbf{v}$</li>
  </ol>
</div>

Many concepts in linear algebra and general mathematics are derived from the vector space, including linear combinations, subspaces, and bases. It's important to note that, though we usually think of vectors as tuples, they don't need to be. You could define a vector to be different cheeses, and as long as the definition is satisfied, it will be a valid vector space.

If we equip a vector space with a special type of map, then we get an <i>inner product space</i>. 

<div id="inner-product-space"></div>
<div class="definition">
  <strong>Definition (Inner Product Space).</strong>
  <br>
  An <i>inner product space</i> is a vector space, $V$, over the field, $F$, of real numbers or complex numbers with the map $\langle, \cdot, \cdot, \rangle: V \times V \rightarrow F$, called an <i>inner product</i>, which satisfies the following for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and all $a, b \in F$:

  <ol>
    <li>Conjugate Symmetry: $\langle \mathbf{u}, \mathbf{v} \rangle = \overline{\langle \mathbf{v}, \mathbf{u} \rangle}$</li>
    <li>Linearity (in the first argument): $\langle a \mathbf{u} + b \mathbf{v}, \mathbf{w} \rangle = a \langle \mathbf{u}, \mathbf{w} \rangle + b \langle \mathbf{v}, \mathbf{w} \rangle$</li>
    <li>Positive-Definiteness: $\langle \mathbf{u}, \mathbf{u} \rangle > 0$ for any $\mathbf{u} \neq \mathbf{0}$</li>
  </ol>
</div>

Something that will be very useful is a map from a vector space to the real numbers that can be thought of as assigning a "size" to vectors in the space. We call this a <i>norm</i>, and if we equip a vector space with a norm, then we have a <i>normed vector space</i>. 

<div id="norm"></div>
<div class="definition">
  <strong>Definition (Norm).</strong>
  <br>
  Let $V$ be a vector space over a scalar field $K$. A <i>norm</i>, $\rvert \rvert \cdot \rvert \rvert: V \rightarrow \mathbb{R}$, is a map satisfying:

  <ol>
    <li>Non-negativity: $\rvert \rvert x \rvert \rvert \geq 0$ for all $x \in V$ </li>
    <li>Positive definiteness: $\rvert \rvert x \rvert \rvert = 0$ if and only if $x$ is the zero vector, for all $x \in V$ </li>
    <li>Absolute homogeneity: $\rvert \rvert \lambda x \rvert \rvert = \rvert \lambda \rvert \rvert \rvert x \rvert \rvert$ for all $\lambda \in K$ and $x \in V$</li>
    <li>Triangle inequality: $\rvert \rvert x + y \rvert \rvert \leq \rvert \rvert x \rvert \rvert + \rvert \rvert y \rvert \rvert$ for all $x, y \in V$</li>
  </ol>
</div>

We can define the <i>canonical norm</i> of an inner product space as $\rvert \rvert x \rvert \rvert \sqrt{\langle x, x \rangle}$. Thus, any inner product space is a normed vector space. A special type of normed vector space is the <i>Banach space</i>.


<div id="banach-space"></div>
<div class="definition">
  <strong>Definition (Banach Space).</strong>
  <br>
  Let $X$ be a vector space over a scalar field (perhaps $\mathbb{R}$ or $\mathbb{C}$), and let $\rvert \rvert \cdot \rvert \rvert: X \rightarrow \mathcal{R}$ be a norm. Together, $(X, \rvert \rvert \cdot \rvert \rvert)$ form a <i>normed space</i>. If this space is also complete, then $(X, \rvert \rvert \cdot \rvert \rvert)$ is a <i>Banach space</i>.
  <br>
  <i>Note: Any finite-dimensional normed vector space is a Banach space. This includes finite-dimensional Euclidean spaces.</i>
</div>

By "complete", we mean that the space does not have any "holes" in it. Formally put, any Cauchy sequence taking values in $X$ converges to a point in $X$ as well.


The <i>induced metric</i> (i.e. the distance metric induced by the norm of a vector space) is the function $d: V \times V \rightarrow \mathbb{R}$ satisfying $d(x, y) = \rvert \rvert x - y \rvert \rvert$ for all $x,y \in V$. If we combine a metric with a set, then we get a <i>metric space</i>, which is just a set on which we have a particular sense of distance between its elements.

<div id="metric-space"></div>
<div class="definition">
  <strong>Definition (Metric Space).</strong>
  <br>
  Let $X$ be a set, and let $d$ be a distance metric. A <i>metric space</i> is the ordered pair $(X, d)$. A metric space is called <i>complete</i> if every Cauchy sequence in $X$ converges to a point in $X$.
</div>


Using our definitions of inner product and complete metric spaces, we can define what is known as a <i>Hilbert space</i>.

<div id="hilbert-space"></div>
<div class="definition">
  <strong>Definition (Hilbert Space).</strong>
  <br>
  A <i>Hilbert space</i> is a real (or complex) inner product space that is also a complete metric space where the distance metric is that induced by its inner product. 
</div>

### Helpful Definitions

<div id="open-cover"></div>
<div class="definition">
  <strong>Definition (Open Cover).</strong>
  <br>
  An <i>open cover</i> of subset $A \subseteq \mathbb{R}$ is any collection $\mathcal{C}$ of open subsets of $\mathbb{R}$ such that $A \subseteq \bigcup_{C \in \mathcal{C}} C$. 
  <br>
  A <i>finite subcover</i> of an open cover $\mathcal{C}$ of $A$ is any finite subset of sets in $\mathcal{C}$. 
</div>


<div id="almost-every"></div>
<div class="definition">
  <strong>Definition (Almost Every).</strong>
  <br>
  Let $(X, \mathcal{S}, \mu)$ be a measure space, and let $A \in \mathcal{S}$. We say that $A$ contains $\mu$-<i>almost every</i> element of $X$ if $\mu(X \setminus A) = 0$ (in words, if $A$ contains all of $X$ except a subset of measure $0$). 
</div>


<div id="almost-everywhere"></div>
<div class="definition">
  <strong>Definition (Almost Everywhere).</strong>
  <br>
  Let $\mu$ be a $\sigma$-finite measure on $(\Omega, \mathcal{F})$, and let $\phi$ and $\psi$ be functions on $(\Omega, \mathcal{F}, \mu)$. We say that $\phi \geq \psi$ <i>$\mu$-almost everywhere</i> if $\mu(\{ \omega: \phi(\omega) < \psi(\omega) \}) = 0$. 
</div>


<div id="almost-surely"></div>
<div class="definition">
  <strong>Definition (Almost Surely).</strong>
  <br>
  Let $(\Omega, \mathcal{F}, P)$ be a probability space. An event $A$ happens <i>almost surely</i> if $P(A) = 1$. 
</div>




### Assorted Results

<div id="doob"></div>
<div class="theorem">
  <strong>Doob's (First) Convergence Theorem.</strong>
  {% tabs doob %}
  {% tab doob statement %}
  Let $(\Omega, \mathcal{F}, P)$ be a probability space, and let $$\mathbb{F} = (\mathcal{F}_t)_{t \geq 0}$$ be a filtration such that $\mathcal{F}_t$ is a sub-$\sigma$-field of $\mathcal{F}$ for all $t$. (That is, $(\Omega, \mathcal{F}, \mathbb{F}, P)$ is a <i>filtered probability space</i>). Suppose we also have $$X: [0, \infty) \times \Omega \rightarrow \mathbb{R}$$, a right-continuous supermartingale with respect to $\mathbb{F}$. 
  <br>
  For $t \geq 0$, define $$X^-_t = \max\{-X_t, 0 \}$$. Assume $$\underset{t > 0}{\sup} \mathbb{E}[X_t^-] < +\infty$$. Then (the point-wise limit) $$X(\omega) = \underset{t \rightarrow + \infty}{\lim} X_t(\omega)$$ exists and is finite for all $\omega \in \Omega$ except a $P$-null set.
  {% endtab %}
  {% tab doob proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>
