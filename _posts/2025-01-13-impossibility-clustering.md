---
layout: post
title:  "Clustering: An Axiomatic Approach"
date: 14 January 2025
categories: posts
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

Though the journey to this point is a bit confusing, I have recently become interesting in clustering metrics and evaluation. In this post, I'll work through Jon Kleinberg's "An Impossibility Theorem for Clustering".

## Background

First, let's more rigorously define what we mean by _clustering_ and our problem setting as defined by Kleinberg.

Suppose we have a set $S$ of $n \geq 2$ observations, which we notate as $S = \{ 1, 2, \dots, n\}$. We'll define a _distance function_ $d: S \times S \rightarrow \mathbb{R}$ as a function satisfying, for all $i,j,k \in S$:

- Non-Negativity: $d(i,j) \geq 0$
- Symmetry: $d(i,j) = d(j,i)$ 
- Identity of Indiscernibles: $d(i,j) = 0$ if, and only if, $i = j$
- Triangle Inequality: $d(i,k) \leq d(i,j) + d(j,k)$

The last condition (Triangle Inequality) is not necessary for the following discussion, but it is a nice property and requiring it makes a choice of $d$ a metric. 

A _clustering function_, $f$, is a function operating on distance functions $d$ on $S$ that outputs a partition $\Gamma$ of $S$. 

Finally, we'll define some helpful terms and quantities that will be used later.

Let $\alpha \cdot d$ denote scaling the distance function $d$ by $\alpha$. That is $\alpha \cdot d$ is the distance function who assigns distance $\alpha d(i,j)$ between points $i$ and $j$. We'll denote all possible output partitions of $f$ as $\text{Range}(f)$. For a partition $\Gamma$, we'll denote the cluster that point $i$ is put into as $\Gamma^i$. 


#### $\Gamma$-transformation
Let $d$ and $d'$ be distance functions on $S$, and let $\Gamma$ be some partition of $S$. We call $d'$ a $\Gamma$-_transformation_ of $d$ if:

- $d'(i,j) \leq d(i,j)$ for all $i,j \in S$ such that $\Gamma^i = \Gamma^j$
- $d'(i,j) \geq d(i,j)$ for all $i,j \in S$ such that $\Gamma^i \neq \Gamma^j$

In words, a $\Gamma$-transformation of a distance function $d$ will assigns smaller distances to points in the same cluster and larger distances to points in different clusters. 

#### Refinement
Let $\Gamma$ be some partition of $S$. We call another partition $\Gamma'$ of $S$ a _refinement_ of $\Gamma$ if, for every cluster $C' \in \Gamma'$, there is a cluster $C \in \Gamma$ such that $C' \subseteq C$. 

More intuitively, a refinement is just a finer partition, which means each cluster in $\Gamma$ is either also in $\Gamma'$ or is split into multiple smaller clusters in $\Gamma'$.

We can define the partial order $\Gamma' \preccurlyeq \Gamma$ if $\Gamma'$ is a refinement of $\Gamma$.

#### Antichain
An _antichain_ is a collection $\mathcal{A}$ of partitions of $S$ such that for all $\Gamma, \Gamma' \in \mathcal{A}$ where neither $\Gamma$ is a refinement of $\Gamma'$ nor is $\Gamma'$ a refinement of $\Gamma$. We also require that $\Gamma$ and $\Gamma'$ be different partitions. 

#### $(a,b)$-conforming
Let $\Gamma$ be a partition of $S$. A distance function $d$ is said to $(a,b)$-_conform_ to $\Gamma$ if:

- $d(i,j) \leq a$ for all $i,j \in S$ such that $\Gamma^i = \Gamma^j$
- $d(i,j) \geq b$ for all $i,j \in S$ such that $\Gamma^i \neq \Gamma^j$

#### $\Gamma$-forcing
For a partition $\Gamma$ of $S$ and clustering function $f$, a pair of real numbers $(a, b)$ is called $\Gamma$-_forcing_ with respect to $f$ if $f(d) = \Gamma$ for all distance functions that are $(a,b)$-conforming to $\Gamma$.

---

## The Impossibility Theorem

Kleinberg's impossibility theorem is based upon a set of three axioms that characterize a "good" clustering function, $f$.

#### Scale-Invariance
For any distance function $d$ and positive real $\alpha > 0$, $f(d) = f(\alpha \cdot d)$. 

A scale-invariant clustering function should not change its output when all distances between points are changed by some factor. 

#### Richness
$\text{Range}(f)$ should be all possible paritions of $S$.

Richness implies that our clustering function is flexible enough that any partition of $S$ can be achieved if we find the right distance function.

#### Consistency
For any two distance functions $d$ and $d'$ such that $f(d) = \Gamma$ and $d'$ is a $\Gamma$-transformation of $d$, $f(d') = \Gamma$.

A consistent clustering function will output the same partition if we make points within the same cluster closer together and make points in different clusters farther apart.


<div class="theorem">
  <strong>Theorem 2.1 (Impossibility).</strong>
  <br>
  For $n \geq 2$, there does not exists a clustering function $f$ that is scale-invariant, rich, and consistent.

  <details>
  <summary>Proof.</summary> 
  The proof of Theorem 2.1 rests upon the claim that a scale-invariant and consistent clustering fucntion $f$ cannot possibly be rich as $\text{Range(f)}$ forms an antichain (Theorem 3.1 in paper). 
  <br>
  First, suppose we have a consistent clustering function $f$, and let $\Gamma$ be some partition in $\text{Range}(f)$. Because $\Gamma \in \text{Range}(f)$, $\exists d$ such that $f(d) = \Gamma$ (by definition). Define $a'$ as the minimum distance between any two points in the same cluster over all clusters in $\Gamma$, and define $b'$ as the maximum distance between any two points in different clusters over all clusters in $\Gamma$.
  <br>
  Select positive real numbers $a, b$ such that $a < b$, $a \leq a'$, and $b \geq b'$. Notice that, for any distance function $d'$ that $(a,b)$-conforms to $\Gamma$ is a $\Gamma$-transformation of $d$. This is due to the fact that $d'(i,j) \leq a \leq a' \leq d(i,j)$ for $i,j \in S$ such that $\Gamma^i = \Gamma^j$ and $d'(i,j) \geq b \geq b' \geq d(i,j)$ for $i,j \in S$ such that $\Gamma^i \neq \Gamma^j$, which is precisely the definition of a $\Gamma$-transformation.
  <br>
  Since $d'$ is a $\Gamma$-transformation of $d$ and $f$ is consistent by assumption, $(a,b)$ is $\Gamma$-forcing.  
  <br>
  Now we further assume that $f$ is scale-invariant and that $\exists \Gamma_0, \Gamma_1 \in \text{Range}(f)$ such that $\Gamma_0$ is a refinement of $\Gamma_1$. Define $(a_0, b_0)$ and $(a_1, b_1)$ be pairs of reals that are $\Gamma_0$-forcing and $\Gamma_1$-forcing, respectively, such that $a_0 < b_0$ and $a_1 < b_1$. Such pairs can be found by taking any two partitions in $\text{Range}(f)$ and setting $\Gamma_0$ and $\Gamma_1$ according to the argument above.
  <br>
  Let $a_2$ be some number such that $a_2 \leq a_1$ and set $\epsilon$ such that $0 < \epsilon < \frac{a_0 a_2}{b_0}$. We now construct the distance function $d$ that satisfies the following:
  <ul>
  <li><body>$d(i,j) \leq \epsilon$ for $i,j$ such that $\Gamma_0^i = \Gamma_0^j$</body></li>
  <li>$a_2 \leq d(i,j) \leq a_1$ for $i,j$ such that $\Gamma_1^i = \Gamma_1^j$ and $\Gamma_0^i \neq \Gamma_0^j$</li>
  <li>$d(i,j) \geq b_1$ for $i,j$ such that $\Gamma_1^i \neq \Gamma_1^j$</li>
  </ul>
  According to the above conditions, $d$ $(a_1, b_1)$-conforms to $\Gamma_1$. Since $f$ is still assumed to be consistent, this implies $f(d) = \Gamma_1$. 
  <br>
  Let $\alpha = \frac{b_0}{a_2}$ and $d' = \alpha \cdot d$. Since $f$ is scale-invariant, $f(d') = f(d) = \Gamma_1$. However, for $i,j$ such that $\Gamma_0^i = \Gamma_0^j$, $d'(i,j) = \alpha d(i,j) \leq \alpha \epsilon < \frac{a_0 a_2 b_0}{a_2 b_0} = a_0$. In addition, for $i,j$ such that $\Gamma_0^i \neq \Gamma_0^j$, $d'(i,j) = \alpha d(i,j) \geq \alpha a_2 = \frac{b_0 a_2}{a_2} = b_0$. This implies that $d'$ $(a_0, b_0)$-conforms to $\Gamma_0$, which implies that $f(d') = \Gamma_0$. 
  <br>
  Since $\Gamma_0$ is a refinement of $\Gamma_1$, they are distinct partitions, so $\Gamma_0 \neq \Gamma_1$, and we arrive at a contradition.
  </details>
</div>

---

## Characterization Theorem

Kleinberg goes on to prove an additional theorem that describes the partitions achievable by scale-invariant, consistent clustering functions $f$.

<div class="theorem">
  <strong>Theorem 3.2 Characterization.</strong>
  <br>
  Let $\mathcal{A}$ be any antichain of partitions. There exists a scale-invariant, consistent clustering function $f$ such that $\text{Range}(f) = \mathcal{A}$. 

  <details>
  <summary>Proof.</summary>
  Kleinberg's proof uses the <i>sum-of-pairs</i> clustering method which outputs the partition $\Gamma \in \mathcal{A}$ that minimizes $\Phi_d(\Gamma) = \sum_{(i,j) \sim \Gamma} d(i,j)$ where the notation $(i,j) \sim \Gamma$ indicates $\Gamma^i = \Gamma^j$. 
  <br>
  Notice that for any $\alpha > 0$, $\Phi_{\alpha \cdot d}(\Gamma) = \sum_{(i,j) \sim \Gamma} \alpha d(i,j) = \alpha \sum_{(i,j) \sim \Gamma} d(i,j) = \alpha \Phi_{d}(\Gamma)$. Since $\alpha$ is positive, the argmin of $\Phi_{\alpha \cdot d}(\Gamma)$ is equivalent to the argmin of $\Phi_d(\Gamma)$, which implies $f(d) = f(\alpha \cdot d)$ so $f$ is scale-invariant.
  <br>
  <br>
  Fix some $\Gamma \in \mathcal{A}$. Let $d$ be a distance function satisfying:
  <ul>
  <li>$d(i,j) < \frac{1}{n^3}$ for $i,j$ such that $\Gamma^i = \Gamma^j$</li>
  <li>$d(i,j) \geq 1$ for $(i,j)$ such that $\Gamma^i \neq \Gamma^j$</li>
  </ul>
  Notice that $\Phi_d(\Gamma) = \sum_{(i,j) \sim \Gamma} d(i,j) < 1$ since the summation is only over $i,j$ in the same clusters of which there can be, at most, $n^2$ pairs (if there is only one cluster). 
  <br> 
  Also notice that $\Phi_d(\Gamma') < 1$ only if $\Gamma'$ is a refinement of $\Gamma$. To see why, consider $\Gamma'$ that is <i>not</i> a refinement of $\Gamma$. This implies that there is some cluster $C' \in \Gamma'$ such that $C' \not\subseteq C$ for all $C \in \Gamma$. This implies that there exist points $i,j$ such that $\Gamma'^i = \Gamma'^j$, so they are included in the summation in $\Phi_d(\Gamma')$, but $\Gamma^i \neq \Gamma^j$, so $d(i,j) > 1$. If $\Gamma'$ is a refinement of $\Gamma$, then $\Gamma^i = \Gamma^j$, so $d(i,j) < \frac{1}{n^3}$. 
  <br>
  Because $\mathcal{A}$ is an antichain, there is no refinement $\Gamma' \in \mathcal{A}$ of $\Gamma$. Thus, $\Gamma = \underset{\Gamma^* \in \mathcal{A}}{\arg\min} \Phi_d(\Gamma^*)$, which implies $f(d) = \Gamma$. Thus, $f$ is rich, since there exists a $d$ such that $f(d) = \Gamma$ for any $\Gamma \in \mathcal{A}$. 
  <br>
  <br>
  Let $d$ be a distance function such that $f(d) = \Gamma$, and let $d'$ be a $\Gamma$-transformation of $d$. Furthermore, let $\Gamma'$ be some other partition and define $\Delta(\Gamma') := \Phi_d(\Gamma') - \Phi_{d'}(\Gamma')$. We have that $\Delta(\Gamma) = \sum_{(i,j) \sim \Gamma}d(i,j) - \sum_{(i,j) \sim \Gamma} d'(i,j)$ and $\Delta(\Gamma) = \sum_{(i,j) \sim \Gamma'} d(i,j)  - \sum_{(i,j) \sim \Gamma'} d'(i,j)$. Then:
  $$
  \begin{aligned}
  \Delta(\Gamma) &= \sum_{(i,j) \sim \Gamma} d(i,j)  - \sum_{(i,j) \sim \Gamma} d'(i,j) \\
  &= \sum_{(i,j) \sim \Gamma} d(i,j) - d'(i,j) \\
  &\overset{(i)}{\geq} \sum_{(i,j) \sim \Gamma, (i,j) \sim \Gamma'} d(i,j) - d'(i,j) \\
  &\overset{(ii)}{\geq} \sum_{(i,j) \sim \Gamma'} d(i,j) - d'(i,j) \\
  &= \Delta(\Gamma')
  \end{aligned}
  \nonumber
  $$  
  We know $d'$ is a $\Gamma$-transformation of $d$. Thus, $d'(i,j) \leq d(i,j)$ for $(i,j) \sim \Gamma$, implying that $d(i,j) - d'(i,j) \geq 0$ for all $(i,j) \sim \Gamma$. $(i)$ follows from the fact that $\rvert (i,j) \sim \Gamma \cap (i,j) \sim \Gamma' \rvert \leq \rvert (i,j) \sim \Gamma \rvert$.
  <br>
  Similarly, $(ii)$ is equivalent to the summation in $(i)$ plus $(i,j) \sim \Gamma'$ such that $(i,j) \not\sim \Gamma$. For these pairs, $d'(i,j) \geq d(i,j)$, since $d'$ is a $\Gamma$-transformation of $d$. This means $d(i,j) - d'(i,j) \leq 0$ for these pairs, which implies $(ii)$. 
  <br>
  The above argument implies $\Delta(\Gamma) \geq \Delta(\Gamma')$ for any $\Gamma' \in \mathcal{A}$. This implies $\Phi_d(\Gamma) - \Phi_{d'}(\Gamma) \geq \Phi_d(\Gamma') - \Phi_{d'}(\Gamma')$, which implies $\Phi_d(\Gamma) - \Phi_d(\Gamma') \geq \Phi_{d'}(\Gamma) - \Phi_{d'}(\Gamma')$. 
  <br>
  Since $\Gamma'$ minimizes $\Phi_{d'}$ and $\Gamma$ minimizs $\Phi_{d}$, we have that $\Phi_d(\Gamma) - \Phi_{d'}(\Gamma) \leq 0$. Chaining this together with the previous yields $0 \geq \Phi_d'(\Gamma) - \Phi_{d'}(\Gamma')$, which implies $\Phi_{d'}(\Gamma') \geq \Phi_{d'}(\Gamma)$. If $\Gamma'$ minimizes, $\Phi_{d'}$, then it must be the case that $\Gamma' = \Gamma$. Thus, $f(d') = \Gamma$, and therefore $f$ is consistent.
  </details>
</div>

---

## Relaxations

Due to the impossibility theorem, it might be worthwhile to look into easing up the conditions in the axioms. Kleinberg provides a few examples. 

#### Relaxing Richness
Theorem 3.2 (Characterization) is an example of a relaxation of the richness property. If we were satisfied with a clustering function that is scale-invariant, consistent, but only achieves an antichain as its range, then the sum-of-pairs method will work. 

#### Relaxing Consistency
Kleinberg proposes _Refinement-Consistency_, in which a $\Gamma$-transformation should output a refinement of $\Gamma$. Unfortunately, this is not yet enough; a scale-invariant, rich, and refinement-consistent clustering function does not exist. However, if one also relaxes richness to say that all but one (trivial) partition can be achieved by $f$ — termed _Near-Richness_ —, then Theorem 2.1 does not hold. 

An alternative is to relax consistency to something I'll call _Weak Consistency_, which is where if $d'$ is a $f(d)$-transformation of $d$, then either $f(d')$ is a refinement of $f(d)$ or $f(d)$ is a refinement of $f(d')$. There do exist clustering functions that satisfy all three of scale-invariance, richness, and weak consistency.

#### Relaxing Scale-Invariance
Kleinberg does not discuss a relaxation of scale-invariance at length. He mentions that 
<span class="popup" onclick="PopupFunc('pop3')">
  single-linkage clustering
  <span class="popuptext" id="pop3">
    Single-linkage clustering is an agglomerative method in which we start by letting each observation be its own cluster. Then we iteratively combine clusters based upon distances until some stopping criterion is met. 
  </span>
</span>
where we stop combining clusters when their distances exceed some value $r$ satisfies consistency and richness. It satisfies a weaker scale-invariance where we let $f(\alpha \cdot d)$ be a refinement of $f(d)$ when $\alpha > 1$. 

Besides this, I don't see why relaxing scale-invariance any more than this would be an attractive property of a clustering function.


---
# References

Jon Kleinberg. 2002. An impossibility theorem for clustering. In Proceedings of the 16th International Conference on Neural Information Processing Systems (NIPS'02). MIT Press, Cambridge, MA, USA, 463–470.