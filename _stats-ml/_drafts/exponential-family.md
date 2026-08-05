---
layout: distill
title: Exponential Families
description: A Primer
date: 2026-01-27
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

This post is a brief reference for exponential family distributions. 

---

## Background
We'll begin with a few crucial definitions.

<div class="definition">
<strong>Definition (Exponential Family).</strong><d-cite key=schervish1995></d-cite>
Let $f_{X \rvert \Theta}(x \rvert \theta)$ be the density function with respect to a measure, $\nu$, on $(\mathcal{X}, \mathcal{B})$ for some family of parametric distributions with parameter space, $\Omega$. We call this family an <i>exponential family</i> if it can be written as:

$$
\begin{equation}
\label{eq:ef}
f_{X \rvert \Theta}(x \rvert \theta) = c(\theta) h(x) \exp\left( \sum_{i = 1}^k \eta_i(\theta) T_i(x) \right)
\end{equation}
$$

for some measurable functions, $\eta_1, \dots, \eta_k, T_1, \dots, T_k$, and some $k \in \mathbb{Z}$.
</div>

In the above definition, the functions $T_1, \dots, T_k$ are sufficient statistics for $\theta$ and $\eta = (\eta_1, \dots, \eta_k)$ is a function of $\theta$ called the <strong>natural parameter</strong>. 

The function $h(x)$ is non-negative and called the <strong>base measure</strong>. It can be absorbed into the measure $\mu$. 

The function $c(\theta)$ is a normalization constant to ensure that we get a valid density function. It has the form:

$$
c(\theta) = \frac{1}{\int h(x) \exp\left( \sum_{i = 1}^k \eta_i(\theta) T_i(x) \right) \nu(dx) }
$$

With this, we can rewrite Eq. \eqref{eq:ef} as:

$$
\begin{equation}
\label{eq:ef2}
f_{X \rvert \Theta}(x \rvert \theta) = h(x) \exp\left( \sum_{i = 1}^k \eta_i(\theta) T_i(x) - A(\theta) \right)
\end{equation}
$$

where $A(\theta) = \log(c(\theta))$ is the <strong>log partition function</strong> (also called the <strong>cumulant function</strong>). 

If the natural parameter equals $\theta$, then we say that the family is in <strong>canonical form</strong>, and if the dimension of $\theta$ is less than that of $\eta$, then we say the family is <strong>curved</strong>.

We can define the <strong>link function</strong> as the function, $g(\cdot)$, that maps between the mean of the distribution and the natural (or canonical) parameters:

$$
g(\mu) = \eta \iff \mu = g^{-1}(\eta)
$$


---

## Examples

### Poisson
Suppose we have a Poisson random variable, $X$, with parameter, $\lambda$. Its probability mass function is given by:

$$
\begin{aligned}
f(x \rvert \lambda) &= \frac{\lambda^x \exp(- \lambda)}{x!}
&= \frac{1}{x!} \exp\left( \log(\lambda^x) - \lambda \right) \\
&= \frac{1}{x!} \exp\left( x \log(\lambda) - \lambda \right)
\end{aligned}
$$

The Poisson distribution is an exponential family distribution with parameter $\lambda$ and the following:

$$
\begin{aligned}
h(x) &= \frac{1}{x!} \\
T(x) &= x \\
\eta(\theta) &= \log(\theta) \\
A(\eta) &= \theta = \exp(\eta) \\
\log(\mu) &= \eta
\end{aligned}
$$

