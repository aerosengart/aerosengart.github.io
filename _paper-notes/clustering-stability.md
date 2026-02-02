---
layout: distill
title: Clustering Stability
description: What Does It Mean To Be Stable?
date: 2025-01-31
tabs: true
tags: clustering philosophy paper-review
toc:
    - name: Stability (Conceptually)
    - name: Background And Notation
    - name: Results
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2

bibliography: 2025-01-31-clustering-stability.bib
---

My last post related to clustering discussed how to describe a "good" clustering algorithm. One way to measure this is by <i>stability</i>, which I'll define more rigorously later. 

The main paper I'll be summarizing, <i>Sober Look at Clustering Stability</i>, by Shai Ben-David, Ulrike von Luxberg, and Dávid Pál<d-cite key=bendavid2006></d-cite> takes a decision theoretic perspective. This allows the authors to derive some results about the behavior of (un)stable algorithms; however it limits their discussion to clustering algorithms that optimize some objective function. I'll restrict myself to these as well. 

---

## Stability (Conceptually)
Similar to the general idea of a robust estimation procedure, a <i>stable</i> clustering algorithm should produce similar results for small perturbations of the data. As Ben-David, von Luxburg, and Pál point out, stability is often used for tuning $k$, the number of clusters the algorithm should yield. Intuitively, choosing $k$ too large will result in arbitrary splitting of true clusters. On the other hand, a $k$ that is too small will do the opposite and arbitrarily merge true clusters.

For optimization-based clustering algorithms, there are two sources of instability:

- <strong>Sampling</strong>: Variability in our observations will result in different values of our objective function. This, in turn, will influence the clustering. 

- <strong>Data Distribution</strong>: If the data distribution is such that there are multiple optimizers of the objective function, then there is ambiguity in which clustering is the best. 

## Background And Notation
I'll be following the notation in Ben-David, von Luxburg, and Pál fairly closely throughout this post with a few slight tweaks to make it align better with my habits. Let's do a quick run-down of their main definitions.

Our observations $S = \{x_1, \dots, x_n\}$ will be assumed to be i.i.d. from some sample space, $X$, with probability measure $P$. If $X$ is a metric space, its metric will be denoted by $\ell$. $P_S$ denotes the uniform distribution over the finite sample $S$.

<div id="clustering"></div>
<div class="definition">
  <strong>Definition (Clustering).</strong>
  <br>
  Formally, a <i>clustering</i>, $\mathcal{C}: X \rightarrow \mathbb{N}$, of a set, $X$, is a finite partition. <i>Clusters</i> are the sets of data points in $X$ that are in the same group, and we denote these by $C_i := \{ x \in X; \mathcal{C}(x) = i \}$. A <i>clustering algorithm</i> is a function $A$ that outputs a clustering of $X$ given a finite sample $S \subset X$.
</div>

We'll denote the cluster of point $x$ with $C(x)$. If $C(x) = C(y)$, then we write $x \underset{\mathcal{C}}{\sim} y$. 


<div id="clustering-distance"></div>
<div class="definition">
  <strong>Definition (Clustering Distance).</strong>
  <br>
  For family of probability distributions $\mathcal{P}$ over some data space, $X$, and family of clusterings of $X$, $\mathcal{S}$, a <i>clustering distance</i> will be a function $d: \mathcal{P} \times \mathcal{S} \times \mathcal{S} \rightarrow [0, 1]$ that satisfies the following properties for any $P \in \mathcal{P}$ and $\mathcal{C}_1, \mathcal{C}_2, \mathcal{C}_3 \in \mathcal{S}$:
  <ul>
    <li>$d_P(\mathcal{C}_1, \mathcal{C}_1) = 0$</li>
    <li>$d_P(\mathcal{C}_1, \mathcal{C}_2) = d_P(\mathcal{C}_2, \mathcal{C}_1)$</li>
    <li>$d_P(\mathcal{C}_1, \mathcal{C}_3) \leq d_P(\mathcal{C}_1, \mathcal{C}_2) + d_P(\mathcal{C}_2, \mathcal{C}_3)$</li>
  </ul>
</div>

Note that the above definition differs from the traditional definition of a distance metric in that it is not required that if $d_P(\mathcal{C}_1, \mathcal{C}_2) = 0$ then $\mathcal{C}_1 = \mathcal{C}_2$.


With these definitions out of the way, we can finally define stability formally.

<div id="stability"></div>
<div class="definition">
  <strong>Definition (Stability).</strong>
  <br>
  For probability distribution $P$ over data space $X$, clustering distance $d$, and clustering algorithm $A$, the <i>stability of $A$ for a sample of size $n$ with respect to $P$</i> is:

  $$
  \text{stab}(A, P, n) = \underset{S_1 \sim P^n; \\ S_2 \sim P^n}{\mathbb{E}} \left[ d_P(A(S_1), A(S_2)) \right]\nonumber
  $$

  The <i>stability of $A$ with respect to $P$</i> is:

  $$
  \text{stab}(A, P) = \underset{n \rightarrow \infty}{\lim \sup} \left[ \text{stab}(A, P, n) \right]
  $$

  where $\underset{n \rightarrow \infty}{\lim \sup} := \underset{n \rightarrow \infty}{\lim} \left[ \underset{m \geq n}{\sup} \text{stab}(A, P, n) \right]$.
</div>

In words, $\text{stab}(A, P, n)$ is the expected distance between the clusterings output by $A$ whe given two i.i.d. samples of size $n$ drawn from $P$. 

<div id="r-minimizing"></div>
<div class="definition">
  <strong>Definition ($R$-Minimizing).</strong>
  <br>
  Let $opt(P) := \underset{\mathcal{C} \in \mathcal{S}}{\inf} \left[ R(P, \mathcal{C}) \right]$, the optimal risk achieved by clusterings in $\mathcal{S}$ for distribution $P$. For a sample $S \subseteq X$, the <i>empirical risk</i> is the risk over the empirical data distribution, $R(P_S, \mathcal{C})$. A clustering algorithm, $A$, is then called <i>$R$-minimizing</i> if $R(P_S, \mathcal{C}) = opt(P_S)$ for any $S \subseteq X$.
</div>

An $R$-minimizing algorithm will achieve empirical risk equal to the optimal risk (for that distribution).

<div id="risk-optimizing-converging"></div>
<div class="definition">
  <strong>Definition (Risk Optimizing and Risk Converging).</strong>
  <br>
  For some domain set, $X$, some set of clusterings of $X$, $\mathcal{S}$, some set of probability distributions on $X$, $\mathcal{P}$, a clustering algorithm is said to be <i>risk optimizing</i> if its output clustering is chosen to minimize (or maximize if you switch the signs) an objective function (or risk function), $R: \mathcal{P} \times \mathcal{S} \rightarrow \mathbb{R}^+_0$. 
  <br>
  An $R$-minimizing algorithm $A$ is called <i>risk converging</i> if, for every $\epsilon > 0$ and every $\delta \in (0, 1)$, $\exists n_0$ such that $\forall n > n_0$:

  $$
  \underset{S \sim P^n}{\mathbb{P}}\left( R(P, A(S)) < opt(P) + \epsilon \right) > 1 - \delta
  \nonumber
  $$

  Essentially, a risk converging clustering algorithm will achieve risk within $\epsilon$ of the optimal risk with high probability for a large enough sample size.
</div>

For clustering distance, $d$, a probability distribution, $P$, has a <i>unique minimizer</i>, $$\mathcal{C}^*$$ if for any $\eta > 0$, there exists $$\epsilon > 0$$ such that if $$R(P, \mathcal{C}) < opt(P) + \epsilon$$, then $$d_P(\mathcal{C}^*, \mathcal{C}) < \eta$$. In words, this means that any clustering with risk that is within some $\epsilon$ of the risk achieved by the unique minimizer should be very close (by $\eta$) to the unique minimizer. 

$P$ is said to have $m$ distinct minimizers if there are $$\mathcal{C}^*_1, \dots, \mathcal{C}^*_m$$ such that $$d_P(\mathcal{C}^*_i, \mathcal{C}^*_j) > 0$$ for all $$i \neq j$$ such that for all $$\eta > 0$$, there exists $$\epsilon > 0$$ such that if $$R(P, \mathcal{C}) < opt(P) + \epsilon$$, then there exists $$1 \leq i \leq m$$ such that $$d_P(\mathcal{C}^*_i, \mathcal{C}) < \eta$$. 


<div id="measure-preserving-symmetry"></div>
<div class="definition">
  <strong>Definition (Measure-Preserving Symmetry).</strong>
  <br>
  For probability distribution, $P$, over $(X, \ell)$, a function $g: X \rightarrow X$ is a <i>$P$-preserving symmetry</i> of $(X, \ell)$ if:
  <ul>
    <li>$\mathbb{P}(A) = \mathbb{P}(g(A))$ for any $P$-measurable set $A \subseteq X$</li>
    <li>$\underset{x,y \sim P}{\mathbb{P}}(\ell(x,y) = \ell(g(x), g(y))) = 1$</li>
  </ul>
</div>

The above distribution gives us a way to characterize risk functions and clustering distances.

<div id="odd"></div>
<div class="definition">
  <strong>Definition (Distance-Distribution Dependent Risk and Clustering Distance).</strong>
  <br>
   A risk function, $R$ is called <i>ODD</i> if, for every distribution $P$, every $P$-preserving symmetry $g$, and every clustering $\mathcal{C}$, we have that $R(P, \mathcal{C}) = R(P, g(\mathcal{C}))$.
   <br>
   A clustering distance, $d$, is called <i>ODD</i> if, for every distribution $P$, every $P$-preserving symmetry $g$, and any clusterings $\mathcal{C}_1, \mathcal{C}_2$, we have that $d_P(\mathcal{C}_1, \mathcal{C}_2) = d_P(g(\mathcal{C}_1), g(\mathcal{C}_2))$. 
   <br>
  Intuitively, the above state that ODD risks and distances only depend on distances and distributions. That is, transforming the points in ways that do not affect their probabilities or distances will not affect the risk or distance between the points. All common distances are ODD. 
</div>


----

## Results
Ben-David et al. present their first theorem as follows:

<div class="theorem">
<strong>Theorem 10.<d-cite key=bendavid2006></d-cite></strong>
{% tabs theorem-10-bendavid %}
{% tab theorem-10-bendavid statement %}
Any risk converging, $R$-minimizing clustering algorithm will be stable on distribution $P$ if $P$ has a unique minimizer.
{% endtab %}
{% tab theorem-10-bendavid proof %}
The basic idea is to show that, for any value $\xi > 0$, the stability of an algorithm $A$ satisfying the stated properties is less than $\xi$ for sufficiently large sample size $m$. This implies that $A$ is stable on $P$.

Define $A$ as a risk converging, $R$-minimizing algorithm, and fix some $\xi > 0$. Let $$\mathcal{C}^*$$ be the unique minimizer of $P$. First, set some $\delta \in (0, 1)$ and $\eta > 0$ such that $2(\eta + \delta) < \xi$.

Since $$\mathcal{C}^*$$ is the unique minimizer of $P$, we have that $\exists \epsilon > 0$ such that:

$$
R(P, \mathcal{C}) < opt(P) + \epsilon \implies d_P(\mathcal{C}^*, \mathcal{C}) < \eta
\nonumber
$$

by the definition of a unique minimizer for a distribution. Furthermore, by the definition of risk convergence, $\exists m_0$ such that for all $m > m_0$ we have:

$$ 
\underset{S \sim P^m}{\mathbb{P}} \left( R(P, A(S)) \geq opt(P) + \epsilon \right) < \delta
\nonumber
$$

Consider $m > m_0$ and $S \sim P^m$. If $R(P, A(S)) < opt(P) + \epsilon$, then $$d_P(A(S), \mathcal{C}^*) < \eta$$. Since the implication is only one way, the probability of the $R(P, A(S)) < opt(P) + \epsilon$ is less than or equal to the probability of $$d_P(A(S), \mathcal{C}^*) < \eta$$ (since there may be some cases where the latter holds without the former). Thus:

$$
\underset{S \sim P^m}{\mathbb{P}}\left( d_P(A(S), \mathcal{C}^*) \geq \eta \right) \leq \underset{S \sim P^m}{\mathbb{P}}\left( R(P, A(S)) \geq opt(P) + \epsilon\right) < \delta
\nonumber
$$

The bound then follows:

$$
\begin{aligned}
\text{stab}(A, P, m) &= \underset{S_1, S_2 \sim P^m}{\mathbb{E}} \left[ d_P(A(S_1), A(S_2)) \right] \\
&\leq \underset{S_1, S_2 \sim P^m}{\mathbb{E}} \left[ d_P(A(S_1), \mathcal{C}^*) + d_P(\mathcal{C}^*, A(S_2)) \right] & \left(\text{triangle ineq.}\right)\\
&= 2 \underset{S \sim P^m}{\mathbb{E}}\left[ d_P(A(S), \mathcal{C}^*) \right] & \left(\text{i.i.d. samples}\right) \\
&\overset{(i)}{\leq} 2\left( \eta \cdot \underset{S \sim P^m}{\mathbb{P}}\left( d_P(A(S), \mathcal{C}^*) < \eta \right) + 1 \cdot \underset{S \sim P^m}{\mathbb{P}} \left( d_P(A(S), \mathcal{C}^*) \geq \eta \right) \right) &\left(\text{defn. of expectation}\right) \\ 
&\leq 2\left( \eta + \underset{S \sim P^m}{\mathbb{P}} \left( R(P, A(S)) \geq opt(P) + \epsilon \right) \right) \\
&\leq 2(\eta + \delta) \\
&< \xi
\end{aligned}
$$

In $(i)$, I am not entirely sure of the choice of $1$. However, I think we can take $\xi$ to be arbitrarily small as $m \rightarrow \infty$?
{% endtab %}
{% endtabs %}
</div>

The authors emphasize a couple points (see the paper for illustrations):

<ul>
  <li> Algorithms (such as $k$-means and $k$-medians) can be stable for <i>any</i> choice of $k$ </li>
  <li> A stable algorithm does not mean the choice of number of clusters is correct </li>
</ul>

With these points in mind, the authors prove a related theorem.

<div class="theorem">
<strong>Theorem 15.<d-cite key=bendavid2006></d-cite></strong>
{% tabs theorem-15-bd %}
{% tab theorem-15-bd statement %}
Suppose we have an ODD risk function, $R$, an ODD clustering distance $d$, a probability distribution, $P$, with some number $n$ distinct minimizers, and $P$-symmetry, $g$, such that $$d_P(\mathcal{C}^*, g(\mathcal{C}^*)) > 0$$ for every $R$-minimizer $$\mathcal{C}^*$$. Any risk convergent, $R$-minimizing clustering algorithm is <i>unstable</i> on $P$.
{% endtab %}
{% tab theorem-15-bd proof %}
For the $R$-minimizers, denote the clustering they result in as $$\{ \mathcal{C}_1^*, \dots, \mathcal{C}_n^*\}$$. Furthermore, let $$r := \underset{1 \leq i \leq n}{\min} \left( d_P(\mathcal{C}^*_i, g(\mathcal{C}^*_i)) \right)$$. 

Choose $\epsilon > 0$ such that a clustering achieving risk within $\epsilon$ of the best possible risk implies the existence of some $R$-minimizer that is close to $\mathcal{C}$ (less than $r/4$ distance between them). That is, $R(P, \mathcal{C}) < opt(P) + \epsilon$ implies $\exists 1 \leq i \leq n$ such that $d_P(\mathcal{C}^*_i, \mathcal{C}) < r/4$. 

Define $$T = \{ S \in X^m \rvert R(P, A(S)) < opt(P) + \epsilon \}$$, the set of samples from $X$ of size $m$ whose clustering using $A$ achieves risk within $\epsilon$ of $opt(P)$. $A$ is assumed to be risk convergent, so there exists $m_0$ such that for $m > m_0$, $P(T) > 0.9$ (by choosing $\delta = 0.1$).

Now, for $1 \leq i \leq n$, define $$T_i = \{ S \in T \rvert d_p(\mathcal{C}^*_i, A(S)) \leq r/4 \}$$, the samples in $T$ whose clusterings with $A$ are within $r/4$ distance of the $i$-th $R$-minimizing clustering, $$\mathcal{C}^*_i$$. This step is just defining subsets of $T$ that we know exist since $P$ has $n$ distinct minimizers (which implies this $\epsilon$ exists).

Since $P(T) > 0.9$ for sufficiently large $m$, there must be some $i_0$ such that $$P(T_{i_0}) \geq 0.9/n$$ (<i>why?</i>). Since $g$ is a $P$-preserving symmetry and $R$ is assumed to be ODD, we know that for any $S \in T$, $g(S) \in T$ (by definition of ODD). The fact that $g$ is a $P$-preserving symmetry also implies that $$P(g[T_{i_0}]) \geq 0.9/n$$ since $$P(T_{i_0}) \geq 0.9/n$$.

By the definition of $r$, we have that $$d_P(\mathcal{C}^*_{i_0}, g(\mathcal{C}^*_{i_0})) \geq r$$, and by the construction of $T_{i_0}$, we also have that for all $S \in T_{i_0}$, $$d_P(\mathcal{C}^*_{i_0}, A(S)) \leq r/4$$. Recall that $d$ is ODD as well, so for all $$S \in T_{i_0}$$, $$d_P(g(\mathcal{C}^*_{i_0}), A(g(S))) \leq r/4$$. 

By our definition of clustering distances, $d$ must satisfy the triangle inequality. Thus, for any $S \in T_{i_0}$ and $S' \in g[T_{i_0}]$:

$$
d_P(A(S), A(S')) \leq d_P(\mathcal{C}^*_i, A(S)) + d_P(\mathcal{C}^*_i, A(S')) \leq \frac{r}{2}
$$

All of the above leads to the following. For any $m \geq m_0$:

$$
\begin{aligned}
stab(A, P, m) &= \underset{S, S' \sim P^m}{\mathbb{E}} \left[ d_P(A(S), A(S')) \right] \\
&\geq \frac{r}{2} \underset{S, S' \sim P^m}{\mathbb{P}} \left( d_P(A(S), A(S')) \geq \frac{r}{2} \right) &\left(\text{ ignore other case (which must be non-negative) }\right) \\
&\geq \frac{r}{2}\underset{S, S' \sim P^m}{\mathbb{P}} \left( A \in T_{i_0} \cap S' \in g[T_{i_0}] \right) \\
&= \frac{r}{2} \underset{S \sim P^m}{\mathbb{P}} \left( S \in T_{i_0} \right) \underset{S' \sim P^m}{\mathbb{P}}\left( S' \in g[T_{i_0}] \right) &\left(S, S' \text{ independent }\right)\\
&\geq \frac{r0.9^2}{2 n^2}
\end{aligned} 
$$

The stability is therefore always positive (regardless of $m$), which implies that $stab(A,P)$ is as well (since $r > 0$). Thus, $A$ is unstable on $P$.
{% endtab %}
{% endtabs %}
</div>

