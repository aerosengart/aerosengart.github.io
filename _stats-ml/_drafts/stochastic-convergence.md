---
layout: distill
title: Stochastic Convergence
description: A Primer
date: 2026-05-04
tabs: true
tags: theory statistics
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Convergence
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: stats-ml.bib
---

This post covers the main points of Chapter 2 of van der Vaart's <i>Asymptotic Statistics</i>.<d-cite key=vaart1998></d-cite>


## Convergence
We will assume to have a random vector $\mathbf{X} = (X_1, \dots, X_k) \in \mathbb{R}^k$, which is a vector of real random variables. 

<details>
<summary>Rigorous Definition.</summary>
Suppose we have some probability space $(\Omega, \mathcal{F}, P)$, where $\Omega$ is our sample space, $\mathcal{F}$ is our event space (a $\sigma$-algebra), and $P$ is some probability measure. Let $\mathcal{B}(\mathbb{R}^k)$ be the Borel $\sigma$-algebra on $\mathbb{R}^k$. 
<br> 
We call a function $\mathbf{X}: \Omega \rightarrow \mathbb{R}^k$ a <strong>random vector</strong> on $\Omega$ if:
$$
\left\{ \omega \in \Omega \rvert \mathbf{X}(\omega) \in B \right\} \in \mathcal{F}, \hspace{5mm} \forall B \in \mathcal{B}(\mathbb{R}^k)
$$
That is, $\mathbf{X}$ is a Borel measurable function from $(\Omega, \mathcal{F}, P)$ to $\mathbb{R}^k$. 
<br>
We can thus define the probability that $\mathbf{X} \in B$ as:
$$
\mathbb{P}(\mathbf{X} \in B) := P(\left\{ \omega \in \Omega \rvert \mathbf{X}(\omega) \in B \right\})
$$
</details>

We denote the <strong>distribution function</strong> of $\mathbf{X}$ as $F: \mathbf{x} \rightarrow \mathbb{P}(\mathbf{X} \leq \mathbf{x})$.

<div class=definition>
<strong>Definition. (Convergence in Distribution)</strong>
<br>
We say that a sequence of random vectors $\mathbf{X}_n$ <i>converges in distribution</i> to a random vector $\mathbf{X}$ if:
$$
\mathbb{P}(\mathbf{X}_n \leq \mathbf{x}) \rightarrow \mathbb{P}(\mathbf{X} \leq \mathbf{x})
$$
for all $\mathbf{x}$ such that $\mathbf{F}$ is continuous.
Alternative notation includes $\mathbf{X}_n \rightsquigarrow \mathbf{X}$ or $\mathbf{X}_n \rightsquigarrow L$ for some distribution, $L$. This is also called <i>weak convergence</i> or <i>convergence in law</i>.
</div>

Let $d: \mathbb{R}^k \times \mathbb{R}^k \rightarrow \rightarrow \mathbb{R}_{\geq 0}$ denote a distance function. 

<div class=definition>
<strong>Definition. (Convergence in Probability)</strong>
<br>
We say that a sequence of random vectors $\mathbf{X}_n$ <i>converges in probability</i> to a random vector $\mathbf{X}$ if, for all $\epsilon > 0$:
$$
\mathbb{P}(d(\mathbf{X}_n, \mathbf{X}) > \epsilon) \rightarrow 0
$$
Alternative notation includes $\mathbf{X}_n \overset{p}{\rightarrow} \mathbf{X}$. 
</div>

We can also define a stronger sense of convergence.

<div class=definition>
<strong>Definition. (Almost Sure Convergence)</strong>
<br>
We say that a sequence of random vectors $\mathbf{X}_n$ <i>converges almost surely</i> to a random vector $\mathbf{X}$ if:
$$
\mathbb{P}(\underset{n \rightarrow \infty}{\lim} \left\{ d(\mathbf{X}_n, \mathbf{X}) \right\} = 0) = 1
$$
This is denoted by $\mathbf{X}_n \overset{as}{\rightarrow} \mathbf{X}$.
</div>

Convergence in distribution is the weakest sense, followed by convergence in probability, and then almost sure convergence. The latter two also require that each $\mathbf{X}_n$ and $\mathbf{X}$ are defined on the same probability space.

We have the follow relationship between the three senses of convergence.

<div class=theorem>
<strong>Theorem 2.7.</strong>
{% tabs theorem-2-7-vaart %}
{% tab theorem-2-7-vaart statement %}
Let $\mathbf{X}_n$, $\mathbf{Y}_n$, $\mathbf{X}$, and $\mathbf{Y}$ be random vectors, and let $c$ be some constant. The following hold:
<ol>
<li>$\mathbf{X}_n \overset{as}{\rightarrow} \mathbf{X}$ implies $\mathbf{X}_n \overset{p}{\rightarrow} \mathbf{X}$</li>
<li>$\mathbf{X}_n \overset{p}{\rightarrow} \mathbf{X}$ implies $\mathbf{X}_n \rightsquigarrow \mathbf{X}$</li>
<li>$\mathbf{X}_n \overset{p}{\rightarrow} c$ if and only if $\mathbf{X}_n \rightsquigarrow c$</li>
<li>if $\mathbf{X}_n \rightsquigarrow \mathbf{X}$ and $d(\mathbf{X}_n, \mathbf{Y}_n) \overset{p}{\rightarrow} 0$, then $\mathbf{Y}_n \rightsquigarrow \mathbf{X}$</li>
<li>if $\mathbf{X}_n \rightsquigarrow \mathbf{X}$ and $\mathbf{Y}_n \overset{p}{\rightarrow} c$, then $(\mathbf{X}_n, \mathbf{Y}_n) \rightsquigarrow (\mathbf{X}, c)$</li>
<li>if $\mathbf{X}_n \overset{p}{\rightarrow} \mathbf{X}$ and $\mathbf{Y}_n \overset{p}{\rightarrow} \mathbf{Y}$, then $(\mathbf{X}_n, \mathbf{Y}_n) \overset{p}{\rightarrow} (\mathbf{X}, \mathbf{Y})$</li>
</ol>
{% endtab %}
{% tab theorem-2-7-vaart proof %}
TO DO. 
{% endtab %}
{% endtabs %}
</div>

---

## Helpful Theorems

We begin with the <i>Portmanteau Lemma</i>, which provides many different ways to define convergence in distribution. 

<div class=theorem>
<strong>Portmanteau Lemma.</strong>
{% tabs pormanteau-vaart %}
{% tab pormanteau-vaart statement %}
Let $\mathbf{X}_n$ and $\mathbf{X}$ be random vectors. The following are equivalent:
<ol>
<li></li>
</ol>
{% endtab %}
{% tab pormanteau-vaart proof %}
TO DO. 
{% endtab %}
{% endtabs %}
</div>


