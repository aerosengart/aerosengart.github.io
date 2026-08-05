---
layout: distill
title: Hypothesis Testing
description: A Primer
date: 2026-07-01
tabs: true
tags: theory likelihood primer
toc:  
  - name: Background  
  - name: Pure Significance Tests
  - name: Types of Tests
    subsections:
        - name: Simple Null and Simple Alternative
---

This post mostly covers Chapter 4 <i>Theory of Statistics</i> by Mark Schervish.<d-cite key=schervish1995></d-cite> I omit the discussion of the Bayesian perspective for simplicity.

---

## Background
We will let $(S, \mathcal{A}, \mu)$ be our probability space and $V: S \rightarrow \mathcal{V}$ some function that tells us the value of a yet to be known quantity. For example, $V$ could be $\Theta$, which parametrizes the distribution of our observations, or a measurable fucntion of $\Theta$. We will assume to have an <i>action space</i>, denoted by $\aleph$, that represents the space of all possible actions (decisions) we could make. We measure how poor our action $a \in \aleph$ is when $V = v \in \mathcal{V}$ with a <i>loss function</i>, $L: \mathcal{V} \times \aleph \rightarrow \mathbb{R}$. 

We come to the main concept of this section: <i>the hypothesis</i>.

<div id="hypothesis"></div>
<div class="definition">
<strong>Definition (Hypothesis; Schervish, pg. 214).</strong>
<br>
Assume that $\mathcal{V}$ can be partitioned into two portions, $\mathcal{V}_0$ and $\mathcal{V}_1$, such that $\mathcal{V}_0 \cup \mathcal{V}_1 = \mathcal{V}$ and $\mathcal{V}_0 \cap \mathcal{V}_1 = \emptyset$. We call the statement $V \in \mathcal{V}_0$ the <strong>null hypothesis</strong>, and we denote it with $H_0$. Similarly, we call the statement $V \in \mathcal{V}_1$ the <strong>alternative hypothesis</strong> and denote it with $H_1$. 
<br>
If $\aleph = \{ 0, 1\}$ and $L(v, a)$ satisfies both $L(v, 1) > L(v, 0)$ for all $v \in \mathcal{V}_0$ and $L(v, 0) > L(v, 1)$ for all $v \in \mathcal{V}_1$, then we call the problem <strong>hypothesis testing</strong>. We say that we <i>reject the null hypothesis</i> if $a = 1$, and we say we <i>fail to reject</i> if $a = 0$. 
</div>

<aside><p>Schervish uses the notation $H$ for the <strong>hypothesis</strong> and $A$ for the <strong>alternative</strong>, but I prefer to use null/alternative.</p></aside>

Suppose we are in the situation where $V = \Theta$. Then we instead partition the parameter space into $\Omega_0$ and $\Omega_1$ and use $H_0: \Theta \in \Omega_0$ and $H_1: \Theta \in \Omega_1$. We say that a hypothesis is <strong>simple</strong> if its corresponding partition has only one element. 

We can define a <strong>randomized decision rule</strong>, $\delta$, where we reject $H_0$ with probability given by the measurable function $\phi: \mathcal{X} \rightarrow [0, 1]$ (see <a href="/stats-ml/decision-theory">my post on decision theory</a> for more details). This function is called the <strong>test function</strong> and defined as:

$$
\phi(x) = \delta(x, 1) = \mathbb{P}(\text{reject } H_0 \rvert X = x)
$$

Thus, it is possible for us to reject $H_0$ when it is true and fail to reject when it is not true. These are called <strong>Type I</strong> and <strong>Type II</strong> errors, respectively. We now define characteristics and concepts related to a test, $\phi$.

<div id="power-hypothesis"></div>
<div class="definition">
<strong>Definition (Schervish, pg. 215).</strong>
<br>
Assume $V = \Theta$. We define the <strong>power function</strong> of $\phi$ as:

$$
\beta_\phi(\theta) = \mathbb{E}_{\theta}\left[ \phi(X) \right]
$$

The <strong>operating characteristic curve</strong> is:

$$
\rho_{\phi}(\theta) = 1 - \beta_{\phi}(\theta)
$$

The <strong>size</strong> and <strong>base</strong> of the test $\phi$ are, respectively:

$$
\underset{\theta \in \Omega_0}{\sup} \left\{ \beta_\phi(\theta) \right\}
\hspace{5mm} \text{and} \hspace{5mm}
\underset{\theta \in \Omega_1}{\inf} \left\{ \beta_\phi (\theta) \right\}
$$

If $\phi$ has a size of at most $\alpha$ for $0 \leq \alpha \leq 1$, then it is said to have <strong>level</strong> $\alpha$. Similarly, if $\phi$ has base at least $\gamma$, then it is said to have <strong>floor</strong> $\gamma$. 
</div>

How do we determine which tests are better than others? One way is consider each test's power. 

<div id="mp-test"></div>
<div class="definition">
<strong>Definition (Most Powerful Test; Schervish, pg. 230).</strong>
<br>
Consider $\Omega = \Omega_0 \cup \{ \theta_1 \}$ where $\theta_1 \notin \Omega_0$, and let $\phi$ be a level $\alpha$ test for $H_0: \Theta \in \Omega_0$ vs. $H_1: \Theta = \theta_1$. $\phi$ is called a <i>most powerful (MP) level $\alpha$ test</i> if, for every other level $\alpha$ test $\psi$, $\beta_{\psi}(\theta_1) \leq \beta_\phi(\theta_1)$. 
<br>
If, more generally, we have $\Omega = \Omega_0 \cup \Omega_1$, we call $\phi$ a <strong>uniformly most powerful (UMP) level $\alpha$ test</strong> if $\beta_\psi(\theta) \leq \beta_\phi(\theta)$ for all other $\theta \in \Omega_1$. 
</div>

Somewhat less important is the <strong>most cautious (MC) floor $\alpha$ test</strong>. If $\Omega = \Omega_1 \cup \{ \theta_0 \}$ with $\theta_0 \notin \Omega_1$, and $\phi$ is a level $\alpha$ test for $H_0: \Theta = \theta_1$ vs. $H_1: \Theta \in \Omega_1$, we call $\phi$ a most cautious floor $\alpha$ test if, for every other level $\alpha$ test $\psi$, $\beta_\psi(\theta_0) \geq \beta_\phi(\theta_0)$. A <strong>uniformly most cautious (UMC) level $\alpha$ test</strong> has an analagous definition to UMP tests. 

For these definitions, notice that we fix the test level at $\alpha$ and then consider test with the best power to be "best". This basically fixes our tolerance for Type I error, which is often considered more "costly" than Type II errors, and minimizes the Type II error probability subject to this constraint. 


---

## Pure Significance Tests
Suppose we define a <a href="https://en.wikipedia.org/wiki/Weak_ordering"><i>weak order</i></a> on the sample space, $\mathcal{X}$, denoted by $\preceq$, and we <i>only</i> state a null hypothesis, $H_0$. For $x, y \in \mathcal{X}$, we interpret $x \preceq y$ to mean that y is "more at odds" with $H_0$ than $x$.<d-cite key=schervish1995></d-cite> Suppose, now, that we observe $x \in \mathcal{X}$. We can define a <strong>pure significance test</strong> based upon the probability of observing $y \in \mathcal{X}$ such that $x \preceq y$. 

<div id="pure-significance-test"></div>
<div class="definition">
<strong>Definition (Pure Significance Test; Schervish, pg. 217).</strong>
<br>
Let $M = (\mathcal{X}, \mathcal{B})$ be a measurable space (Borel space), and let $H_0$ be a hypothesis defined on $M$. Assume that $p_{Q}(x) = Q(\{ y: x \preceq \})$ is (approximately) the same for all $Q \in H_0$. We call $p_Q(x)$ the <strong>significance probability of the data $x$ relative to the weak order $\preceq$</strong>. A <i>pure significance test</i> is one in which we reject $H_0$ when $p_{H_0}(x)$ is small enough. 
</div>

One way to make a pure significance test is to define a test statistic, $T: \mathcal{X} \rightarrow \mathbb{R}$. Our weak order can then be defined as satisfying $x \preceq y$ if, and only if, $T(x) \leq T(y)$. 

---

## Types of Tests
In this section, we will be comparing the <i>risk functions</i> of different tests. Recall that the risk function is the expected loss over the sample space (see <a href="/stats-ml/decision-theory">my post on decision theory</a>):

$$
R(\theta, \delta) = \int_{\mathcal{X}} \int_{\mathcal{V}} L(v, \delta(x)) dP_{\theta, V}(v) dP_{\theta}(x) = \int_{\mathcal{X}} L(\theta, \delta(x)) dP_{\theta}(x)
$$

Also recall that a decision rule $\delta$ is called <strong>inadmissible</strong> if there exists another $\delta'$ such that $R(\theta, \delta) \leq R(\theta, \delta')$ for all $\theta$ and strict inequality for at least one $\theta$. 

### Simple Null and Simple Alternative
Tests of this kind have null and alternative hypotheses that specify single points in the null and alternative parameter spaces. That is, the tests are of the form:

$$
H_0: \Theta = \theta_0 
\hspace{5mm} \text{vs.} \hspace{5mm} 
H_1: \Theta = \theta_1
$$




### Simple Null and Simple Alternative


### Simple Null and Simple Alternative







<div id="ump-test"></div>
<div class="definition">
<strong>Definition (Schervish, pg. 215).</strong>
<br>
</div>