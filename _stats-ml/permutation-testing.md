---
layout: distill
title: Permutation Testing
description: A Primer
date: 2026-01-29
tabs: true
tags: theory testing primer
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Background
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: stats-ml.bib
---

Recently, I've spent a bunch of time invested in parametric hypothesis testing settings. However, in real life, it isn't always reasonable to assume that we know the distribution of our data or to rely upon assumptotic results. In this post, I'll provide some background and theory behind <i>permutation testing</i>, which is a way to conduct distribution-free hypothesis tests.

I'll be using Good's <i>Permutation Testing</i><d-cite key=good2010></d-cite> and Persarin and Salmas' <i>Permutations Tests for Complex Data: Theory, Applications and Software</i> <d-cite key=pesarin2010></d-cite> for the most part. 

---

## Background
We will have the same kind of set-up as for general testing procedures. We will assume that we have some observation of random variable(s), $X$, and we will assume that $X$ has distribution $P_{\theta} \in \mathcal{P} = \{ P_\theta \rvert \theta \in \Theta \}$. We also have a set of possible decisions, $\mathcal{D} = \{ d \}$, that can result from a decision rule, $\delta$, and some chosen loss function, $\mathcal{L}(\theta, \delta)$, that is a measure of the performance of $\delta$ when $\theta$ is true. 

A test, denoted by $\phi$, is a decision rule with values in $[0, 1]$. If $\phi(x) = 1$, we reject our null hypothesis in favor of our alternative hypothesis. If $\phi(x) = 0$, we fail to reject the null hypothesis. If $\phi(x) \in (0, 1)$, then reject the null with probability $\phi(x)$. 

---










We can define a permutation test mathematically as follows.

<div class="definition">
<strong>Definition (Permuation Test).<d-cite key=good2010></d-cite></strong>
<br>
Let $\mathbf{z}$ be an $N$-dimensional vector of observations, let $T(\mathbf{z}) \in \mathbb{R}$ be some statistic computed for $\mathbf{z}$, and let $A: \mathbb{R} \times \mathbb{R} \rightarrow [0, 1]$ be an acceptance criterion. An <i>$\alpha$-level permutation test</i>, denoted by $\phi$, is a decision rule such that, for all $\mathbf{z}$, $\phi(\mathbf{z}) = 1$ if, and only if:

$$
W(\mathbf{z}) = \sum_{\pi \in \Pi} A(T(\mathbf{z}), \mathbf{T}(\pi \mathbf{z})) \leq \alpha N!
$$

where $\Pi$ is the set of all possible rearrangements of observations, and $\pi \mathbf{z}$ produces on such permutation of the vector $\mathbf{z}$.
</div>
