---
layout: distill
title: Score Tests
description: A Primer
date: 2026-07-20
tabs: true
tags: likelihood theory testing score
toc:
  - name: Background
    subsections:
        - name: Notation
        - name: Assumptions
bibliography: stats-ml.bib
---

In this post, I cover the basics of score-type tests. 

For a more in-depth discussion of one-sided tests/constrained alternative parameter spaces, see Kudô (1967)<d-cite key=kudo1967></d-cite> and Silvapulle and Silvapulle (1997)<d-cite key=silvapulle1997></d-cite>.


## Background

### Notation
We will concern ourselves with the random variable $X$, which we assume follows a distribution from some class of distributions, $$\mathcal{P}_{\boldsymbol{\theta}} = \left\{ P_{\boldsymbol{\theta}} : \boldsymbol{\theta} \in \boldsymbol{\theta} \right\}$$, parametrized by the $k$-dimensional vector $\boldsymbol{\theta}$ with parameter space $\boldsymbol{\theta}$. We let $$\mathbf{X} = (X_1, \dots, X_n)$$ denote the collection of $n$ such random variables, and $$\mathbf{x} = (x_1, \dots, x_n)$$ to denote the collection of their realizations.

We write the log-likelihood, score function, and information matrix for one observation as:

$$
\begin{aligned}
\ell(\boldsymbol{\theta}; \mathbf{x}) &= \log f(\mathbf{x}; \boldsymbol{\theta}) \\
\mathbf{S}(\boldsymbol{\theta}) &= \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x})}{\partial \boldsymbol{\theta}} \\
\mathcal{I}(\boldsymbol{\theta}) &= \mathbb{E}_{\boldsymbol{\theta}} \left[ \frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \frac{\partial \log f(\mathbf{x}; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}^\top}\right]
\end{aligned}
$$

We will let $$\boldsymbol{\theta}_0$$ denote the true value of the parameter vector, and we will let $$\hat{\boldsymbol{\theta}}$$ denote the global maximum likelihood estimate of the parameter vector over all of $\boldsymbol{\Theta}$. For some $$\mathbf{A} \in \boldsymbol{\Theta}$$ such that $$\boldsymbol{\theta}_0 \in \mathbf{A}$$, we define the MLE over $\mathbf{A}$ as the value of the parameter vector, $$\hat{\boldsymbol{\theta}}_\mathbf{A}$$, satisfying:

$$
\ell(\hat{\boldsymbol{\theta}}_\mathbf{A}; \mathbf{x})
  = \underset{\boldsymbol{\theta} \in \mathbf{A}}{\sup} \left\{ \ell(\boldsymbol{\theta}; \mathbf{x}) \right\} + o_p(1)
$$

For some sequence $\{ \mathbf{x}_n \}_n$, we will denote its probability limit when $\boldsymbol{\theta}$ is true value with $$\underset{\boldsymbol{\theta}}{plim} [\mathbf{x}_n]$$, so:

$$
\underset{n \rightarrow \infty}{\lim} \left\{ \rvert \mathbf{x}_n - \underset{\boldsymbol{\theta}}{plim} [\mathbf{x}_n] \rvert > \epsilon \right\} = 0; \hspace{2mm} \forall \epsilon > 0 
$$


---

## Two-Sided Score-Type Tests
We will assume to have $n$ observations which depend on some $p$-dimensional parameter vector, $\boldsymbol{\theta}$. We assume we can partition $\boldsymbol{\theta}$ into the $(p-k)$-dimensional incidental/nuisance parameters, $\boldsymbol{\lambda}$, and the $k$-dimensional structural parameters, $\boldsymbol{\psi}$, as $\boldsymbol{\theta} = (\boldsymbol{\lambda}^\top, \boldsymbol{\psi}^\top)^\top$. We considering testing


$$
\begin{equation}
\label{eq:hypothesis-two-sided}
H_0: \boldsymbol{\psi} = \mathbf{0}_k 
\hspace{5mm} \text{vs.} \hspace{5mm}
H_1: \boldsymbol{\psi} \neq \mathbf{0}_k
\end{equation}
$$

We will let $\mathbf{S}_n(\boldsymbol{\theta}) = \mathbf{0}_k$ be a $k$-dimensional estimating equation for $\boldsymbol{\theta}$. As an example later, we will let $\mathbf{S}_n(\boldsymbol{\theta})$ be the (classical) score function. We assume the following condition holds:

<div class="theorem">
<strong>Condition A (Silvapulle and Silvapulle, pg. 343).</strong>
<br>
There exist non-singular matrices, $\mathbf{G}(\boldsymbol{\theta})$ and $\mathbf{V}(\boldsymbol{\theta})$, such that:

$$
\begin{align}
&\frac{1}{\sqrt{n}} \mathbf{S}_n(\boldsymbol{\theta}) \rightsquigarrow \mathcal{N}(\boldsymbol{0}_p, \mathbf{V}(\boldsymbol{\theta})) \label{eq:cdtn-a-1} \\
&\underset{\rvert \rvert \mathbf{h} \rvert \rvert \leq a}{\sup} \left\{ \frac{1}{\sqrt{n}} \left[ \mathbf{S}_n\left(\boldsymbol{\theta} + \frac{1}{\sqrt{n}} \mathbf{h} \right) - \mathbf{S}_n(\boldsymbol{\theta}) \right] + \mathbf{G}(\boldsymbol{\theta}) \mathbf{h} \right\} = o_p(1) \label{eq:cdtn-a-2}
\end{align}
$$

as $n \rightarrow \infty$ for any $a > 0$.
</div>

Eq. \eqref{eq:cdtn-a-1} states that $\mathbf{S}_n(\boldsymbol{\theta})$ is $\mathbf{G}(\boldsymbol{\theta})$ is asymptotically normal with mean zero when suitably scaled. Eq. \eqref{eq:cdtn-a-2} essentially implies that $\mathbf{G}(\boldsymbol{\theta}) = \frac{\partial}{\partial \boldsymbol{\theta}} \left[ - \frac{1}{n} \mathbf{S}_n(\boldsymbol{\theta}) \right]$ in the limit. 

<aside><p>See <a href="/paper-notes/score-test-one-sided/">this post</a> for more discussion on the interpretation of the above conditions.</p></aside>

### Example --- Score Function
Suppose we pick $\mathbf{S}_n(\boldsymbol{\theta}) = \frac{\partial \ell(\boldsymbol{\theta}; \mathbf{x})}{\partial \boldsymbol{\theta}}$, the score function. Suppose the 

<details>
<summary>Conditions.</summary>
We will assume the following regularity conditions hold (see <a href="/stats-ml/maximum-likelihood#regularity-conditions">this post</a> for details).

<ol>
<li>If $P_{\boldsymbol{\theta}_1} = P_{\boldsymbol{\theta}_2}$ for $\boldsymbol{\theta}_1, \boldsymbol{\theta}_2 \in \boldsymbol{\theta}$, then $\boldsymbol{\theta}_1 = \boldsymbol{\theta}_2 \label{cdtn:reg-1}$.</li>
<li>There exists a neighborhood, $N_{\boldsymbol{\theta}_0}$, of the true parameter value, $\boldsymbol{\theta}_0$, satisfying $N_{\boldsymbol{\theta}_0} \subseteq \boldsymbol{\theta}\label{cdtn:reg-2}$.</li>
<li>The $X_1, \dots, X_n$ are i.i.d. with continuous probability density or mass function $f(x; \boldsymbol{\theta}) := f_{X | \boldsymbol{\theta}}(x; \boldsymbol{\theta}) \label{cdtn:reg-3}$.</li>
<li>The support of $f(x; \boldsymbol{\theta})$, denoted by $\mathbf{A}$, is independent of $\boldsymbol{\theta} \label{cdtn:reg-4}$. (AND IS THE SAME FOR ALL $X_i$?) </li>
<li>The first-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to the coordinates of $\boldsymbol{\theta}$, $\frac{\partial f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i}$, exist for all $x \in \mathbf{A}$, $\boldsymbol{\theta} \in \boldsymbol{\Theta}$, and $i = 1, \dots, k\label{cdtn:reg-5}$.</li>
<ol>
<li>The second-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to the coordinates of $\boldsymbol{\theta}$, $\frac{\partial^2 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j}$, exist for all $x \in \mathbf{A}$, $\boldsymbol{\theta} \in \boldsymbol{\Theta}$, and $i,j = 1, \dots, k\label{cdtn:reg-5a}$.</li>
<li>The third-order partial derivatives of $f(x; \boldsymbol{\theta})$ with respect to the coordinates of $\boldsymbol{\theta}$, $\frac{\partial^3 f(x; \boldsymbol{\theta})}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_l}$, exist for all $x \in \mathbf{A}$, $\boldsymbol{\theta} \in \boldsymbol{\Theta}$, and $i,j,l = 1, \dots, k\label{cdtn:reg-5b}$.</li>
</ol>
<li>The order of differentiation with respect to $\boldsymbol{\theta}$ and integration can be exchanged when taking first-order partial derivatives of $\int f(x; \boldsymbol{\theta}) dx \label{cdtn:reg-6}$.</li>
<ol>
<li>The order of differentiation with respect to $\boldsymbol{\theta}$ and integration can be exchanged when taking second-order partial derivatives of $\int f(x; \boldsymbol{\theta}) dx \label{cdtn:reg-6a}$. </li>
<li>The order of differentiation with respect to $\boldsymbol{\theta}$ and integration can be exchanged when taking third-order partial derivatives of $\int f(x; \boldsymbol{\theta}) dx \label{cdtn:reg-6b}$.</li>
</ol>
<li>There exist $c(\boldsymbol{\theta}_0) \in \mathbb{R}^+$ and $M_{i,j,l}(x)$ such that:
    \begin{eqn}
        M_{i,j,l}(x) 
            &\geq \left\rvert \frac{\partial^3  \log( f(x; \boldsymbol{\theta})) }{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j \partial \boldsymbol{\theta}_l } \right\rvert \hspace{5mm} \forall \boldsymbol{\theta} \in \mathbf{C}_{\boldsymbol{\theta}_0} \\
        \text{where } \mathbf{C}_{\boldsymbol{\theta}_0} &= \left\{ \rvert \rvert \boldsymbol{\theta} - \boldsymbol{\theta}_0 \rvert \rvert_2^2 < c(\boldsymbol{\theta}_0) \hspace{2mm} \rvert \hspace{2mm} \boldsymbol{\theta} \in \boldsymbol{\theta} \right\}
    \end{eqn}
and $\mathbb{E}_{\boldsymbol{\theta}_0} \left[ M_{i,j,l}(x) \right] < \infty \label{cdtn:reg-7}$.</li>
<li>For all $i,j \in [k]$, $\mathcal{I}_{i,j}(\boldsymbol{\theta}) < \infty$, and $\mathcal{I}(\boldsymbol{\theta}) \label{cdtn:reg-8}$ is positive definite.
</li>
</ol>
</details>






---

## One-Sided Score-Type Tests

