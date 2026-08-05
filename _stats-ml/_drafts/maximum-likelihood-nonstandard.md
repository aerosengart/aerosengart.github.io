---
layout: distill
title: Maximum Likelihood Estimation (Non-Standard)
description: A Primer
date: 2026-07-28
tabs: true
tags: theory likelihood primer
toc:  
  - name: Background
    subsections:
        - name: Notation and Set-Up
        - name: Definitions
        - name: Properties
bibliography: stats-ml.bib
---

This post mostly covers <i>Asymptotic Properties of Maximum Likelihood Estimators and Likelihood Ratio Tests Under Nonstandard Conditions</i> by Self and Liang.<d-cite key=self1987></d-cite>

---

## Background 

### Notation and Set-Up
In what follows, we will basically use the same definitions and set-up in my <a href="/stats-ml/maximum-likelihood">previous maximum likelihood post</a>. 

We will let $X_1, \dots, X_n$ be i.i.d., and $X$ will denote an arbitrary observation. We'll use lowercase $x$ to denote realizations of $X$. A parameter will be denoted with $\boldsymbol{\theta}$ with particular values denoted by superscripts and subscripts and its true value by $$\boldsymbol{\theta}^*$$. We will also use the notation $\mathbf{X} = (X_1, \dots, X_n)$ and $\mathbf{x} = (x_1, \dots, x_n)$ to denote a collection of $n$ random variables or vectors $X$ and observations of each, respectively. 

Let $$\mathcal{P}_\boldsymbol{\theta}$$ be a family of distributions for $X$ parametrized by $\boldsymbol{\theta} \in \boldsymbol{\Theta} \subseteq \mathbb{R}^p$. We will assume $X_1, \dots, X_n \sim P_{\boldsymbol{\theta}^*}$. The probability density or mass function of $X$ given parameter $\boldsymbol{\theta}$ evaluated at a particular $x$ and $$\boldsymbol{\theta}_0$$ will be denoted by $$f_{X \rvert \boldsymbol{\theta}}(x; \boldsymbol{\theta}_0)$$ or, more compactly, $$f(x; \boldsymbol{\theta}_0)$$. 

---

## Regularity Conditions
Below are the conditions needed for stating some of the results discussed later in this post. The statements come from pg. 605 in Self and Liang.<d-cite key=self1987></d-cite>

#### Condition 1
If $P_{\boldsymbol{\theta}_1} = P_{\boldsymbol{\theta}_2}$ for $\boldsymbol{\theta}_1, \boldsymbol{\theta}_2 \in \boldsymbol{\Theta}$, then $\boldsymbol{\theta}_1 = \boldsymbol{\theta}_2$.

<aside><p>Same as in standard settings.</p></aside>

This implies that the distributions in the class $\mathcal{P}_{\boldsymbol{\theta}}$ are distinct.

#### Condition 2


#### Condition 3