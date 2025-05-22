---
layout: post
title:  "Clustering: An Axiomatic Approach"
date: 14 January 2025
categories: posts
tags: ["clustering"]
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

Though the journey to this point is a bit confusing, I have recently become interesting in clustering metrics and evaluation. In this post, I'll work through a couple papers on describing how good a clustering function is based upon a set of axioms. These include Jon Kleinberg's "An Impossibility Theorem for Clustering"[^fn-kleinberg] and Margareta Ackerman and Shai Ben-David's "Measures of Clustering Quality: A Working Set of Axioms for Clustering"[^fn-ackerman].

---

## Background

First, let's more rigorously define what we mean by _clustering_ and our problem setting as defined by Kleinberg.

<div id="clustering"></div>
<div class="definition">
  <body>
  <strong>Definition (Clustering).</strong>
  <br>
  Suppose we have a set $S$ of $n \geq 2$ observations, which we notate as $S = \{ 1, 2, \dots, n\}$. We'll define a <i>distance function</i> $d: S \times S \rightarrow \mathbb{R}$ as a function satisfying, for all $i,j,k \in S$:

  <ul>
    <li>Non-Negativity: $d(i,j) \geq 0$</li>
    <li>Symmetry: $d(i,j) = d(j,i)$</li>
    <li>Identity of Indiscernibles: $d(i,j) = 0$ if, and only if, $i = j$</li>
    <li>Triangle Inequality: $d(i,k) \leq d(i,j) + d(j,k)$</li>
  </ul>

  The last condition (Triangle Inequality) is not necessary for the following discussion, but it is a nice property and requiring it makes a choice of $d$ a metric. 
  <br>
  A <i>clustering function</i>, $f$, is a function operating on a distance function $d$ on $S$ that outputs a partition $\Gamma$ of $S$. A partition is called <i>trivial</i> if each cluster contains one point or if it has only one cluster.
  </body>
</div>

We can evaluate a clustering function with what the author terms a <i>clustering quality metric</i>.

<div id="cqm"></div>
<div class="definition">
  <body>
  <strong>Definition (Clustering Quality Metric).</strong>
  <br>
  A <i>clustering-quality measure (CQM)</i> is a function operating on a partition $\Gamma$ of a set $S$ with respect to distance $d$ that outputs a non-negative real number which represents the "goodness" of the clustering $\Gamma$. That is, a CQM $m$ is the function $m: (\Gamma, S, d) \rightarrow \mathbb{R}^+_0$. 
  </body>
</div>


Finally, we'll define some helpful terms and quantities that will be used later. Let $\alpha \cdot d$ denote scaling the distance function $d$ by $\alpha$. That is $\alpha \cdot d$ is the distance function who assigns distance $\alpha d(i,j)$ between points $i$ and $j$. We'll denote all possible output partitions of $f$ as $\text{Range}(f)$. We denote the fact points $i$ and $j$ are in the same cluster in partition $\Gamma$ with $i \underset{\Gamma}{\sim} j$, and we use the notation $i \underset{\Gamma}{\not \sim} j$ if they are not in the same cluster. 


<div id="rep-set"></div>
<div class="definition">
  <body>
  <strong>Definition (Representative Set).</strong>
  <br>
  We call a set $G$ a <i>representative set</i> of $\Gamma$ if it contains a single observation from each cluster in $\Gamma$. That is, for $\Gamma = \{ \gamma_1, \dots, \gamma_g\}$, $G$ is a representative set of $\Gamma$ if $\rvert G \rvert = g$ and $G\cap \gamma_i \neq \emptyset$ \forall i$. 
  </body>
</div>

<div id="isomorphism"></div>
<div class="definition">
  <body>
  <strong>Definition (Isomorphism).</strong>
  <br>
  Let $\Gamma$ and $\Gamma'$ be partitions of $S$ with $d$. We call $\Gamma$ and $\Gamma'$ <i>isomorphic</i>, denoted by $\Gamma \underset{d}{\approx} \Gamma'$ if there is a
  <span class="popup" onclick="PopupFunc('pop4')">
    <i>distance-preserving isomorphism</i>
    <span class="popuptext" id="pop4">
      We can think of this as a mapping between our source set, $S$, and some target set $\phi(S)$ such that for any $i \in S$, there exists $\phi(i) \in \phi(S)$; for any $i, j \in S$, $d(i, j) = d(\phi(i), \phi(j))$; and there exists an inverse mapping, $\phi'$, such that $\phi'(\phi(i)) = i$. 
    </span>
  </span>
 $\phi: S \rightarrow S$ such that $\forall i,j \in S$, we have $i \underset{\Gamma}{\sim} j$ if, and only if, $\phi(i) \underset{\Gamma'}{\sim} \phi(j)$. 
  </body>
</div>

<div id="gamma-transformation"></div>
<div class="definition">
  <body>
  <strong>Definition ($\Gamma$-Transformation).</strong>
  <br>
  Let $d$ and $d'$ be distance functions on $S$, and let $\Gamma$ be some partition of $S$. We call $d'$ a <i>$\Gamma$-transformation</i> of $d$ if:
  <ul>
    <li>$d'(i,j) \leq d(i,j)$ for all $i,j \in S$ such that $i \underset{\Gamma}{\sim} j$</li>
    <li>$d'(i,j) \geq d(i,j)$ for all $i,j \in S$ such that $i \underset{\Gamma}{\not \sim} j$</li>
  </ul>
  In words, a $\Gamma$-transformation of a distance function $d$ will assigns smaller distances to points in the same cluster and larger distances to points in different clusters. 
  </body>
</div>

<div id="refinement"></div>
<div class="definition">
  <body>
  <strong>Definition (Refinement).</strong>
  <br>
  Let $\Gamma$ be some partition of $S$. We call another partition $\Gamma'$ of $S$ a <i>refinement</i> of $\Gamma$ if, for every cluster $C' \in \Gamma'$, there is a cluster $C \in \Gamma$ such that $C' \subseteq C$. 
  <br>
  More intuitively, a refinement is just a finer partition, which means each cluster in $\Gamma$ is either also in $\Gamma'$ or is split into multiple smaller clusters in $\Gamma'$.
  We can define the partial order $\Gamma' \preccurlyeq \Gamma$ if $\Gamma'$ is a refinement of $\Gamma$.
  </body>
</div>

<div id="antichain"></div>
<div class="definition">
  <body>
  <strong>Definition (Antichain).</strong>
  <br>
  An <i>antichain</i> is a collection $\mathcal{A}$ of partitions of $S$ such that for all $\Gamma, \Gamma' \in \mathcal{A}$ where neither $\Gamma$ is a refinement of $\Gamma'$ nor is $\Gamma'$ a refinement of $\Gamma$. We also require that $\Gamma$ and $\Gamma'$ be different partitions. 
  </body>
</div>

<div id="ab-conforming"></div>
<div class="definition">
  <body>
  <strong>Definition ($(a,b)$-Conforming).</strong>
  <br>
  Let $\Gamma$ be a partition of $S$. A distance function $d$ is said to <i>$(a,b)$-conform</i> to $\Gamma$ if:
  <ul>
    <li>$d(i,j) \leq a$ for all $i,j \in S$ such that $i \underset{\Gamma}{\sim} j$</li>
    <li>$d(i,j) \geq b$ for all $i,j \in S$ such that $i \underset{\Gamma}{\not \sim} j$</li>
  </ul>
  </body>
</div>


<div id="gamma-forcing"></div>
<div class="definition">
  <body>
  <strong>Definition ($\Gamma$-Forcing).</strong>
  <br>
  For a partition $\Gamma$ of $S$ and clustering function $f$, a pair of real numbers $(a, b)$ is called <i>$\Gamma$-forcing</i> with respect to $f$ if $f(d) = \Gamma$ for all distance functions that are $(a,b)$-conforming to $\Gamma$.
  </body>
</div>

---

## The Impossibility Theorem

Kleinberg's impossibility theorem is based upon a set of three axioms that characterize a "good" clustering function, $f$.

<div id="axioms"></div>
<div class="definition">
  <body>
  <strong>Definition (Kleinberg's Axioms<span markdown="1">[^fn-kleinberg]</span>).</strong>
  <ul>
    <li><i>Scale Invariance:</i> For any distance function $d$ and positive real $\alpha > 0$, $f(d) = f(\alpha \cdot d)$. A scale-invariant clustering function should not change its output when all distances between points are changed by some factor. </li>
    <li><i>Richness:</i> $\text{Range}(f)$ should be all possible paritions of $S$. Richness implies that our clustering function is flexible enough that any partition of $S$ can be achieved if we find the right distance function.</li>
    <li><i>Consistency:</i> For any two distance functions $d$ and $d'$ such that $f(d) = \Gamma$ and $d'$ is a $\Gamma$-transformation of $d$, $f(d') = \Gamma$. A consistent clustering function will output the same partition if we make points within the same cluster closer together and make points in different clusters farther apart.</li>
  </ul>
  </body>
</div>

Kleinberg goes on to show that these axioms imply a semi-surprising result.

<div id="theorem2.1"></div>
<div class="theorem">
  <strong>Theorem 2.1 (Kleinberg<span markdown="1">[^fn-kleinberg]</span>).</strong>
  <br>
  For $n \geq 2$, there does not exists a clustering function $f$ that is scale-invariant, rich, and consistent.

  <details>
  <summary>Proof.</summary> 
  The proof of Theorem 2.1 rests upon the claim that a scale-invariant and consistent clustering fucntion $f$ cannot possibly be rich as $\text{Range(f)}$ forms an antichain (Theorem 3.1 in paper). 
  <br>
  First, suppose we have a consistent clustering function $f$, and let $\Gamma$ be some partition in $\text{Range}(f)$. Because $\Gamma \in \text{Range}(f)$, $\exists d$ such that $f(d) = \Gamma$ (by definition). Define $a'$ as the minimum distance between any two points in the same cluster over all clusters in $\Gamma$, and define $b'$ as the maximum distance between any two points in different clusters over all clusters in $\Gamma$.
  <br>
  Select positive real numbers $a, b$ such that $a < b$, $a \leq a'$, and $b \geq b'$. Notice that, for any distance function $d'$ that $(a,b)$-conforms to $\Gamma$ is a $\Gamma$-transformation of $d$. This is due to the fact that $d'(i,j) \leq a \leq a' \leq d(i,j)$ for $i,j \in S$ such that $i \underset{\Gamma}{\sim} j$ and $d'(i,j) \geq b \geq b' \geq d(i,j)$ for $i,j \in S$ such that $i \underset{\Gamma}{\not \sim} j$, which is precisely the definition of a $\Gamma$-transformation.
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

Kleinberg proves an additional theorem that describes the partitions achievable by scale-invariant, consistent clustering functions $f$.

<div id="theorem3.2"></div>
<div class="theorem">
  <strong>Theorem 3.2 (Kleinberg<span markdown="1">[^fn-kleinberg]</span>).</strong>
  <br>
  Let $\mathcal{A}$ be any antichain of partitions. There exists a scale-invariant, consistent clustering function $f$ such that $\text{Range}(f) = \mathcal{A}$. 

  <details>
  <summary>Proof.</summary>
  Kleinberg's proof uses the <i>sum-of-pairs</i> clustering method which outputs the partition $\Gamma \in \mathcal{A}$ that minimizes $\Phi_d(\Gamma) = \sum_{(i,j) \sim \Gamma} d(i,j)$ where the notation $(i,j) \sim \Gamma$ indicates $i \underset{\Gamma}{\sim} j$. 
  <br>
  Notice that for any $\alpha > 0$, $\Phi_{\alpha \cdot d}(\Gamma) = \sum_{(i,j) \sim \Gamma} \alpha d(i,j) = \alpha \sum_{(i,j) \sim \Gamma} d(i,j) = \alpha \Phi_{d}(\Gamma)$. Since $\alpha$ is positive, the argmin of $\Phi_{\alpha \cdot d}(\Gamma)$ is equivalent to the argmin of $\Phi_d(\Gamma)$, which implies $f(d) = f(\alpha \cdot d)$ so $f$ is scale-invariant.
  <br>
  <br>
  Fix some $\Gamma \in \mathcal{A}$. Let $d$ be a distance function satisfying:
  <ul>
  <li>$d(i,j) < \frac{1}{n^3}$ for $i,j$ such that $i \underset{\Gamma}{\sim} j$</li>
  <li>$d(i,j) \geq 1$ for $(i,j)$ such that $i \underset{\Gamma}{\not \sim} j$</li>
  </ul>
  Notice that $\Phi_d(\Gamma) = \sum_{(i,j) \sim \Gamma} d(i,j) < 1$ since the summation is only over $i,j$ in the same clusters of which there can be, at most, $n^2$ pairs (if there is only one cluster). 
  <br> 
  Also notice that $\Phi_d(\Gamma') < 1$ only if $\Gamma'$ is a refinement of $\Gamma$. To see why, consider $\Gamma'$ that is <i>not</i> a refinement of $\Gamma$. This implies that there is some cluster $C' \in \Gamma'$ such that $C' \not\subseteq C$ for all $C \in \Gamma$. This implies that there exist points $i,j$ such that $i \underset{\Gamma'}{\sim} j$, so they are included in the summation in $\Phi_d(\Gamma')$, but $i \underset{\Gamma}{\not \sim} j$, so $d(i,j) > 1$. If $\Gamma'$ is a refinement of $\Gamma$, then $i \underset{\Gamma}{\sim} j$, so $d(i,j) < \frac{1}{n^3}$. 
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
<a href="#theorem3.2">Theorem 3.2</a> is an example of a relaxation of the richness property. If we were satisfied with a clustering function that is scale-invariant, consistent, but only achieves an antichain as its range, then the sum-of-pairs method will work. 


#### Relaxing Consistency
Kleinberg proposes <i>Refinement-Consistency</i>, in which a $\Gamma$-transformation should output a refinement of $\Gamma$. Unfortunately, this is not yet enough; a scale-invariant, rich, and refinement-consistent clustering function does not exist. However, if one also relaxes richness to say that all but one (trivial) partition can be achieved by $f$ — termed _Near-Richness_ —, then Theorem 2.1 does not hold. 

An alternative is to relax consistency to something I'll call <i>Weak Consistency</i>, which is where if $d'$ is a $f(d)$-transformation of $d$, then either $f(d')$ is a refinement of $f(d)$ or $f(d)$ is a refinement of $f(d')$. There do exist clustering functions that satisfy all three of scale-invariance, richness, and weak consistency. 

In some ways, this relaxation may be more reasonable. For example, consider some partition that results in four clusters arranged in a square. Now, construct a distance function that just puts more distance between the left clusters and the right clusters. Although this new distance function is a $\Gamma$-transformation, it might be better to combine the left clusters and right clusters to have a partition with two groups. Ackerman and Ben-David provide this example as an illustration in Figure 1 of their paper.  

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

## The Consistency Theorem
Kleinberg's main result is a theorem (and its proof) stating that there does not exist a clustering function that satisfies three simple and desirable properties. This seems quite disappointing as clustering is a popular subtopic in unsupervised learning, and it seems to imply that clustering is at least impossibly difficult to define precisely if not simply impossible to actually do. 

Ackerman and Ben-David push back against this interpretation and state that the way in which Kleinberg axiomatized clustering was part of the reason for the impossibility result. They provide a slightly different perspective and provide a set of axioms for the function used to assess clustering quality rather than the clustering function itself.

Ackerman and Ben-David adjust Steinberg's work to apply to CQMs rather than clustering functions, which results in a consistent set of axioms. They also add one additional property in order to make their set of axioms satisfy two properties (<i>soundness</i> and <i>completeness</i>, which I won't go into here) that make it more useful for defining what methods should be used for clustering.


<div id="axioms2"></div>
<div class="definition">
  <body>
  <strong>Definition (Ackerman and Ben-David Axioms<span markdown="1">[^fn-ackerman]</span>).</strong>
  <ul>
    <li><i>Scale Invariance:</i>A CQM $m$ is called <i>scale-invariant</i> if for every partition $\Gamma$ of set $S$ with respect to distance $d$ and ever $\alpha > 0$, we have $m(\Gamma, S, d) = m(\Gamma, S, \alpha \cdot d)$. </li>
    <li><i>Richness:</i>A CQM $m$ is called _rich_ if for each non-trivial partitions $\Gamma$ of $S$, there exists distance function $d$ such that $\Gamma = \underset{\Gamma'}{\arg \max} \left[ m(\Gamma', S, d) \right]$. A CQM will satisfy the richness property if, for each non-trivial $\Gamma$, we have $m(\Gamma, S, d) \geq m(\Gamma', S, d)$ for all possible partitions $\Gamma'$ of $S$ with $d$ (which may be chosen for each $\Gamma$). The $\max$ becomes a $\min$, and the inequality is reversed for CQMs that assign lower values to better clusterings.</li>
    <li><i>Consistency:</i>A CQM $m$ is called <i>consistent</i> if for every $\Gamma$ of $S$ with $d$, $m(\Gamma, X, d') \geq m(\Gamma, S, d)$ for any $\Gamma$-transformation, $d'$, of $d$. This condition is weaker than the consistency defined by Kleinberg as it does not penalize clustering functions that are only weakly consistent (recall that this means $f(d)$ and $f(d')$ are permitted to be refinements of each other).</li>
    <li><i>Isomorphism Invariance:</i>A CQM $m$ is called <i>isomorphism-invariant</i> if, for all $\Gamma$, $\Gamma'$ of $S$ with $d$ such that $\Gamma$ and $\Gamma'$ are isomorphic, we have $m(\Gamma, S, d) = m(\Gamma', S, d)$. This condition basically states that if we have two clusterings that would be the same if we swapped the points around (in a special way...), they should have the same score according to the CQM.</li>
  </ul>
  </body>
</div>

Ackerman and Ben-David then prove that consistent CQMs exist.

<div class="theorem">
  <strong>Theorem 2 & 3 (Ackerman and Ben-David<span markdown="1">[^fn-ackerman]</span>).</strong>
  <br>
  There exists a clustering quality measure that satisfies all four of scale-invariance, richness, consistency. and isomorphism invariance. That is, the four properties comprise a consistent set of axioms.

  <details>
  <summary>Proof.</summary>
  It suffices to construct a CQM and prove it satisfies the three axioms. 
  <br>
  First, we define the <i>Relative Point Margin</i>. For distance $d$ and clustering $\Gamma$ containing $k$ clusters, the <i>$G$-Relative Point Margin</i> is defined as $G$-$RM_{S, d}(x) = \frac{d(x, g_x)}{d(x, g'_{x})}$ where $g_x \in G$ is the closest cluster center to $x$, $g'_{x} \in G$ is the second closest cluster center to $x$, and $G \subseteq S$. 
  <br>
  Next, define the <i>Relative Margin</i>. Let $\mathcal{G}$ denote the set of all possible representative sets of $\Gamma$. Then the relative margin of $\Gamma$ over $(S, d)$ is $RM_{S, d}(\Gamma) = \underset{G \in \mathcal{G}}{\min} \left\{ \underset{x \in S \setminus G}{\textit{avg}} G\text{-}RM_{S, d}(x) \right\}$. This is the representative set that achieves the minimum average relative point margin where the average is taken over all points not in the representative set. The relative margin assigns lower values to better clusterings (so the inequalities in the richness and consistency definitions will be reversed).
  <br>
  Let $\Gamma$ be an arbitrary clustering of the set $S$ with distance function $d$ in the following.
  <br>
  <br>
  <i>Scale-Invariance.</i> 
  Let $d'$ be a distance function satisfying $d'(i,j) = \alpha d(i, j)$ for all $i,j \in S$ and all $\alpha > 0$. For any $i,j,k \in S$, we have:
  $$
  \frac{d'(i,j)}{d'(i,k)} = \frac{\alpha d(i,j)}{\alpha d(i,k)} = \frac{d(i,j)}{d(i,k)} \implies \frac{d'(i, g_i)}{d'(i, g'_i)} = \frac{d(i, g_i)}{d(i, g'_i)}
  \nonumber  
  $$
  since scaling all distances will result in the same centers. This implies that $RM_{S, d'}(\Gamma) = RM_{S, d}(\Gamma)$.
  <br>
  <br>
  <i>Consistency.</i> 
  Let $d'$ be a $\Gamma$-transformation of $d$. We have:
  $$
  \begin{aligned}
  \begin{cases}
  d'(i,j) \leq d(i,j) & \text{for } i \underset{\Gamma}{\sim} j \\
  d'(i,j) \geq d(i,j) & \text{for } i \underset{\Gamma}{\not \sim} j
  \end{cases}
  &\implies 
  \frac{d'(i,j)}{d'(i,k)} \leq \frac{d(i,j)}{d(i,k)} \text{ for } i \underset{\Gamma}{\sim} j; i \underset{\Gamma}{\not \sim} k  \\
  &\implies
  G\text{-}RM_{S,d'}(i) \leq G\text{-}RM_{S,d}(i) \text{ for any } G \in \mathcal{G}
  \end{aligned}
  \nonumber
  $$
  The first implication follows from the fact that the closest cluster center to $i$ will be a point in the same cluster, and the second closest cluster center to $i$ will be a point in a different cluster (by the definition of a representative set).
  <br>
  <br>
  <i>Richness.</i> 
  Let $\Gamma$ be an arbitrary non-trivial clustering of $S$ with $d$. Define the distance function $d$ as:
  $$
  \begin{cases}
  d(i,j) = 1 & \text{for } i \underset{\Gamma}{\sim} j \\
  d(i,j) = 10 & \text{for } i \underset{\Gamma}{\not \sim} j
  \end{cases}
  \nonumber
  $$
  It follows that:
  $$
  \begin{aligned}
  \underset{\Gamma'}{\arg \min} \left\{ RM_{S, d}(\Gamma') \right\} &= \underset{\Gamma'}{\arg \min} \left\{ \underset{G \in \mathcal{G}}{\min} \left[ \underset{x \in S \setminus G}{\textit{avg}} G\text{-}RM_{S,d}(x) \right]  \right\} \\
  &= \underset{\Gamma'}{\arg \min} \left\{ \underset{G \in \mathcal{G}}{\min} \left[ \frac{1}{\rvert S \setminus G \rvert} \sum_{x \in S \setminus G} \frac{d(x, g_x)}{d(x, g'_x)} \right]  \right\} \\ 
  &= \underset{\Gamma'}{\arg \min} \left\{ \underset{G \in \mathcal{G}}{\min} \left[ \frac{1}{\rvert S \setminus G \rvert} \sum_{x \in S \setminus G} \frac{1}{10} \right]  \right\} \\
  &= \underset{\Gamma'}{\arg \min} \left\{ \underset{G \in \mathcal{G}}{\min} \left[ \frac{1}{10} \right] \right\}
  \end{aligned}
  \nonumber
  $$
  We arrive at the fact that, for this choice of $d$, the minimum $RM_{S,d}(\Gamma')$ is achieved by every non-trivial partition of $S$. Thus, $\Gamma = \underset{\Gamma'}{\arg \min} \left\{ RM_{S, d}(\Gamma') \right\}$.
  <br>
  <br>
  <i>Isomorphism Invariance.</i>
  Let $\Gamma'$ be a partition such that $\Gamma \underset{d}{\approx} \Gamma'$. Since they are isomorphic, there exists a distance-preserving isomorphism $\phi$. Let $G' := \{ \phi(x): x \in G \}$, and let $\mathcal{G}'$ be the set of all $G'$. Thus:
  $$
  \begin{aligned}
  \frac{d(i, g_i)}{d(i, g'_i)} = \frac{d(\phi, g_{\phi(i)})}{d(\phi, g'_{\phi(i)})} &\implies G\text{-}RM_{S,d}(i) = G'\text{-}RM_{S,d}(i) \\
  &\implies \underset{i \in S \setminus G}{avg} \left\{ G\text{-}RM_{S,d}(i)\right\} = \underset{i \in S \setminus G'}{avg} \left\{ G'\text{-}RM_{S,d}(i) \right\} \\
  &\implies \underset{G \in \mathcal{G}}{\min} \left[ \underset{i \in S \setminus G}{avg} \left\{ G\text{-}RM_{S,d}(i)\right\}\right] = \underset{G \in \mathcal{G}}{\min} \left[ \underset{i \in S \setminus G'}{avg} \left\{ G'\text{-}RM_{S,d}(i)\right\}\right] \\
  &\implies RM_{S,d}(\Gamma) = RM_{S,d}(\Gamma')
  \end{aligned}
  \nonumber
  $$

  </details>
</div>

## Examples and Extensions

Ackerman and Ben-David provide a few different examples of CQMs that satisfy their axioms.

#### Weakest Link
Suppose we are in a linkage-based regime. Define the _Weakest Link Between Points_ as the following. Let $\gamma_k$ be the $k$-th cluster in $\Gamma$. For $\Gamma$ over $S$ with $d$:
$$
\Gamma\text{-}WL_{S, d}(i, j) = \underset{x \in \gamma_k \forall k}{\min} \left\{ \max \left[ d(i, x), d(i, j) \right] \right\}
\nonumber
$$
The _Weakest Link_ of $\Gamma$ over $S$ with $d$ is:
$$
WL(\Gamma) = \frac{\underset{i \underset{\Gamma}{\sim} j}{\max} \Gamma\text{-}WL_{S,d}(i,j)}{\underset{i \underset{\Gamma}{\not \sim}}{\min} d(i,j)}
\nonumber
$$
This is calculable in $O(n^3)$ time.

#### Additive Margin
Suppose we are in a center-based clustering (like $k$-means). Define the _Additive Point Margin_ as:
$$
G\text{-}AM_{S, d}(i) = d(i, \gamma'_i) - d(i, \gamma_i)
\nonumber
$$
The _Additive Margin_ of $\Gamma$ over $S$ with $d$ is:
$$
AM_{S,d}(\Gamma) = \underset{G \in \mathcal{G}}{\min} \left\{ 
\frac{
  \frac{1}{\rvert S \rvert} \sum_{i \in S} G\text{-}AM_{S, d}(i)
}{
  \frac{1}{\rvert \left\{ \{ i, j\} \subseteq S \rvert i \underset{\Gamma}{\sim} j \right\}} \sum_{i \underset{\Gamma}{\sim} j} d(i, j)
}
\right\}
\nonumber
$$
This is calculable in $O(n^{k+1})$ time.

#### Functions of CQMs
One can also use functions of clustering-quality measures to create a new one. If one had a CQM defined for partitions of $k$ clusters, one could consider taking the minimum, maximum, or average over all subsets of size $k$ for a clustering of arbitrary size greater than $k$.

#### Cluster Number Dependence
The above CQMs do not depend on the number of clusters in a partition, which makes it easy to compare clusterings that are of different sizes. Ackerman and Ben-David extend their framework to CQMs that _do_ depend on cluster number, such as those that are based upon the objective functions of some clustering methods (such as $k$-means). 

In order to make these CQMs compliant with the scale-invariance property, the quality scores must be normalized in some way. An example is _$\mathcal{L}$-normalization_, which scales the loss of a clustering by the loss of a trivial clustering that has a single cluster of all observations. 

Loss-based clustering functions also tend to either reward or punish more clusters. Ackerman and Ben-David term CQMs based on a loss function that "prefers" more clusters as _refinement-preferring_ and those based on losses that prefer fewer clusters as _coarsening-preferring_. More explicitly, a refinement-preferring CQM will assign a better quality score to refinements of $\Gamma$ than to $\Gamma$, and a coarsening-preferring CQM will assign better quality scores to $\Gamma$ than to its refinements. These CQMs do not satisfy the richness property.

---

## References
[^fn-ackerman]: Ackerman, M. and Ben-David, S. (2008). Measures of clustering quality: a working set of axioms for clustering. In Proceedings of the 22nd International Conference on Neural Information Processing Systems (NIPS'08). Curran Associates Inc., Red Hook, NY, USA, 121–128.

[^fn-kleinberg]: Kleinberg, J. (n.d.). An Impossibility Theorem for Clustering. In Proceedings of the 16th International Conference on Neural Information Processing Systems (NIPS'02). MIT Press, Cambridge, MA, USA, 463–470.