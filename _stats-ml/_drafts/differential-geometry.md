---
layout: distill
title: Differential Geometry
description: A Primer
date: 2026-04-01
tabs: true
tags: primer geometry calculus
toc:
    - name: Background
bibliography: 2026-04-02-diff-geo.bib
---

In this post, I will be covering some of the bones of differential geometry using <i>Introduction To Differential Geometry</i> by Robbin and Salamon.<d-cite key=robbin2024></d-cite> At its core, <a href="https://en.wikipedia.org/wiki/Differential_geometry">differential geometry</a> is the study of the properties of smooth (i.e. differentiable) shapes and spaces. It's useful for some of the papers related to my work on the score test, and I am having trouble conceptualizing some of the definitions, so I thought it would be best to try to formalize some of the basics for myself.

---

## Background

### Topology 
We begin with some foundational definitions starting with a <i>topological space</i>.

<div class="definition">
<strong>Definition (<a href="https://en.wikipedia.org/wiki/Topological_space#Definition_via_open_sets">Topological Space</a>).</strong>
<br>
Let $M$ be a set, and let $\tau$ be a collection of subsets of $M$ (called <strong>open sets</strong>) that satisfy:
<ul>
<li>Both $\emptyset$ and $M$ belong to $\tau$.</li>
<li>Any finite or infinite union of elements in $\tau$ belongs to $\tau$.</li>
<li>Any finite intersection of elements in $\tau$ belongs to $\tau$.</li>
</ul>
We call $\tau$ a <strong>topology</strong> on $M$. The pair $(M, \tau)$ is called a <i>topological space</i>.
</div>

Intutitively, a topological space is a space in which we have an idea of what it means to be "close" but we cannot, in general, measure closeness as in a metric space. We can describe a topological space as <i>homeomorphic</i> to some other topological space if there exists a <i>homeomorphism</i> between the two.

<div class="definition">
<strong>Definition (<a href="https://en.wikipedia.org/wiki/Homeomorphism">Homeomorphism</a>).</strong>
<br>
Let $M$ and $Y$ be topological spaces. A function $f: M \rightarrow Y$ is called a <i>homeomorphism</i> if it is:
<ul>
<li>bijective</li>
<li>continuous</li>
<li>has a continuous inverse</li>
</ul>
</div>

We now come to the <i>topological manifold</i>, which will be the main subject of this primer.

<div class="definition">
<strong>Definition (Topological Manifold).</strong>
<br>
A topological space, $M$, is called a <i>topological manifold</i> if there exists an <a href="https://en.wikipedia.org/wiki/Neighbourhood_(mathematics)">open neighborhood</a>, $U$, about each point $p \in M$ that is homeomorphic to an open subset of a Euclidean space. 
</div>

From a more physical perspective, a topological manifold is a topological space in which we can find some open set containing $p$ that can be stretched, shrunk, or bent into to "look" like an open subset of a Euclidean space. If each neighborhood is homeomorphic to $\mathbb{R}^n$, then we call $M$ an $n$-dimensional manifold. We can now define <i>charts</i> and <i>atlases</i>.

<div class="definition">
<strong>Definition (Chart and Atlas).</strong>
<br>
Let $M$ be a set, and let $U \subset M$. A <i>chart</i> on $M$ is a pair, $(\phi, U)$, where $\phi: U \rightarrow \phi(U)$ is a bijective mapping from $U$ to an open set $\phi(U) \subset \mathbb{R}^m$. 
<br>
An <i>atlas</i> on $M$ is a collection of charts, $\mathcal{A} = \{ (\phi_\alpha, U_{\alpha}) \rvert \alpha \in A \}$ such that $M = \cup_{\alpha \in A} U_\alpha$. 
</div>

We can relate two different charts, $(\phi_1, U_1)$ and $(\phi_2, U_2)$, as being <strong>topologically compatible</strong> if, and only if, $\phi_1(U_1 \cap U_2)$ and $\phi_2(U_1 \cap U_2)$ are open subsets of $\mathbb{R}^m$ and the map:

$$
\phi_2 \circ \phi_1^{-1}: \phi_1(U_1 \cap U_2) \rightarrow \phi_2(U_1 \cap U_2)
$$

is a homeomorphism. $\phi_2 \circ \phi_1^{-1}$ is called a <strong>transition map</strong>, and an atlas, $\mathcal{A}$, is called a <strong>topological atlas</strong> if $\phi_1$ and $\phi_2$ are topologically compatible for any $\phi_1, \phi_2 \in \mathcal{A}$. If $\mathcal{A}$ contains every chart that is smoothly compatible with the others in it, then we call $\mathcal{A}$ a <strong>maximal smooth atlas</strong>.


### Derivatives
Differential geometry requires derivatives, which themselves help us define what it means for a space to be <i>smooth</i>. 

<div class="definition">
<strong>Definition (Smooth).</strong>
<br>
For open sets $U \in \mathbb{R}^n$ and $V \in \mathbb{R}^m$, a map $f: U \rightarrow V$ is called <i>smooth</i> if, and only if, it is infinitely differentiable (i.e. all of its partial derivatives exist and are continuous). 
</div>

We define the derivative as follows.

<div class="definition">
<strong>Definition (Derivative).</strong>
Let $U \subset \mathbb{R}^m$ and $L \subset \mathbb{R}^n$. Let $f = (f_1, \dots, f_m): U \rightarrow V$ be infinitely differentiable, and fix a point $x \in U$. The <i>derivative of $f$ at $x$</i> is the linear map $df(x): \mathbb{R}^n \rightarrow \mathbb{R}^m$ satisfying

$$
df(x) \xi := \left. \frac{d}{dt} \left[ f(x + t \xi) \right] \right\rvert_{t = 0} = \underset{t \rightarrow 0}{\lim} \left[ \frac{f(x + t \xi) - f(x)}{t}\right]
$$

for $\xi \in \mathbb{R}^n$. We can write this map as the <strong>Jacobian matrix</strong> of $f$ at $x$ (the matrix of partial derivatives):

$$
df(x) := \begin{bmatrix}
\frac{\partial f_1}{\partial x_1}(x) & \dots & \frac{\partial f_1}{\partial x_n} (x) \\
\vdots & \ddots & \vdots \\
\frac{\partial f_m}{\partial x_1}(x) & \dots & \frac{\partial f_m}{\partial x_n}(x) 
\end{bmatrix}
$$
</div>

Compositions of smooth maps are also smooth. That is, for open sets $U \subset \mathbb{R}^n$, $V \subset \mathbb{R}^m$, and $W \subset \mathbb{R}^l$ and smooth maps $f: U \rightarrow V$ and $g: V \rightarrow W$, the composition $g \circ f: U \rightarrow W$ is also smooth. $g\circ f$ also satisfies:

$$
d(g \circ f)(x) = dg(f(x)) \circ df(x) : \mathbb{R}^n \rightarrow \mathbb{R}^l; \hspace{5mm} \forall x \in U
$$

<aside><p>The chain rule!</p></aside>

<div class="definition">
<strong>Definition (Diffeomorphism).</strong>
<br>
A map $f: U \rightarrow V$ is called a <i>diffeomorphism</i> if it satisfies:
<ul>
<li>$f$ is bijective</li>
<li>$f$ is smooth</li>
<li>$f^{-1}$ is smooth</li>
</ul>
</div>

### Manifolds
We can now define one of the main concepts in differential geometry: the smooth manifold.

<div class="definition">
<strong>Definition (Smooth Manifold).</strong>
<br>
Let $M$ be a set and let $\mathcalP{A}$ be a maximal smooth atlas on $M$. We call the pair $(M, \mathcal{A})$ a <i>smooth manifold</i>.
</div>





Maps are generally used to figure out how far locations are from each other, so a good $\phi$ should, in theory, preserve the lengths of curves in $\mathbf{S}^2$. That is, for a curve, $\gamma$, between $p, q \in \mathbf{S}^2$, we want $L(\gamma) = L(\phi \circ \gamma)$ for all $\gamma$ where:

$$
L(\gamma) = \int_0^1 \rvert \dot{\gamma}(t) \rvert dt, \hspace{5mm} \gamma(0) = p, \hspace{2mm} \gamma(1) = q
$$

A slightly less strict "ideal" $\phi$ would map the shortest path between $p$ and $q$, which is called a <strong>geodesic</strong>. That is, $\phi$ satisfies $

Unfortunately, no such chart exists. 


