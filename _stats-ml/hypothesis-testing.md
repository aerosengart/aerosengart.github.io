---
layout: distill
title: Hypothesis Testing
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
    - name: Testing
      subsections:
            - name: Non-Randomized and Randomized Tests
            - name: Simple and Composite Hypotheses
            - name: Errors
bibliography: stats-ml.bib
---

A large part of statistical inference is <i>hypothesis testing</i> in which we wish to learn about some aspect of how a set of random variables behave via a sample. For example, we might be interested in learning about the mean of a random variable or about whether it follows a Gaussian distribution. 

In this post, I'll cover some of the basics of hypothesis testing using Lehmann (2008).<d-cite key=lehmann2008></d-cite> Just because I find it easier to conceptualize this way, I'll be considering the case of hypothesis tests concerning a parameter value.

---

## Background
We'll focus on parametric settings, so we'll assume to have a random variable, $X$, taking on values $x \in \mathcal{X}$ with distribution $P_{\theta}$. We assume $P_{\theta}$ falls within some class of distributions, $\mathcal{P} = \{ P_{\theta} \rvert \theta \in \Theta \}$, parametrized by $\theta$ with parameter space, $\Theta$. 

Hypothesis testing involves constructing a <i>decision rule</i>, which is a function that takes in data and outputs a decision that relates to the inferential goals at hand. Letting $\mathcal{D}$ be the set of all possible decisions, we can denote a decision rule with:

$$
\begin{equation}
\label{eq:decision-rule}
\delta: X \rightarrow \mathcal{D}
\end{equation}
$$

We quantify how good a decision rule is with a <i>loss function</i>, which is a function of the choice of parameter (and, therefore, $P_{\theta}$ as $\theta$ uniquely labels each distribution within $\mathcal{P}$) and decision and defined as:

$$
\begin{equation}
\label{eq:loss-function}
\mathcal{L}: \Theta \times \mathcal{D} \rightarrow \mathbb{R}_{\geq 0}
\end{equation}
$$

We define the <i>risk</i> of $\delta$ as the average loss under the assumption that $P_{\theta}$ is the <i>true</i> distribution of $X$:

$$
\begin{equation}
\label{eq:risk}
R(\theta, \delta) = \mathbb{E}_{X \sim P_{\theta}}\left[  \mathcal{L}(\theta, \delta(X)) \right]
\end{equation}
$$

<aside><p>The <strong>risk function</strong> of a decision rule, $\delta$, is its risk as a function of $\theta$.</p></aside>

The question then becomes identifying the decision rule that performs best with respect to the choice of loss function (e.g. minimizes the risk). It's important to keep in mind that the solution to the problem we've described above is dependent upon the assumed distribution class, $\mathcal{P}$, the choice of loss function, $\mathcal{L}$, and the decision space, $\mathcal{D}$. There is a bit of art that goes into making assumptions on these three aspects that are restrictive enough to allow us to identify a solution but not too restrictive so as to make the problem trivial. 

Sometimes, it is too restrictive to use risk minimization as our selection criterion. For one, the true value of $\theta$ is often unknown, making $R(\theta, \delta)$ impossible to compute. One could also imagine that for many values of $\theta$ besides the true one, the decision has very bad performance. 

To combat this, we define <strong>decision procedures</strong>, which are methods that tell us which decision rule to prefer, even in cases when one rule does not always have smaller risk than another. A decision procedure is defined by how we want to judge the risk functions of our decision rules. 

Suppose we do not know the true value of $\theta$, but we do know a priori that it follows some distribution with density function $\rho(\theta)$. We can then use the <i>Bayes risk</i>. 

<div class="definition">
<strong>Definition (Bayes Risk).</strong>
<br>
Assume that $\theta \sim \rho(\theta)$ a priori. For a decision rule, $\delta$, the <i>Bayes risk</i> is defined as the average risk with respect to the parameter's density:

$$
\begin{equation}
\label{eq:bayes-risk}
r(\rho, \delta) = \int \mathbb{E}_{\theta} \left[ \mathcal{L}(\theta, \delta(X)) \right] \rho(\theta) d\theta
\end{equation}
$$
</div>

A decision rule that minimizes a Bayes risk is called a <strong>Bayes solution</strong>. One can also imagine identifying the Bayes solution subject to some constraint on its Bayes risk. This is called a <strong>restricted Bayes solution</strong>.

If we don't know anything about the distribution of $\theta$, then we could instead decide that we want to use the <i>maximum</i> risk as our criterion. 

<div class="definition">
<strong>Definition (Minimax).</strong>
<br>
For a decision rule, $\delta^*$, is called <i>minimax</i> if it minimizes the maximum risk over all parameter values. That is:

$$
\begin{equation}
\label{eq:minimax}
\underset{\theta \in \Theta}{\sup} \left\{ R(\theta, \delta^*) \right\} = \underset{\delta}{\inf} \underset{\theta \in \Theta}{\sup} \left\{ R(\theta, \delta) \right\}
\end{equation}
$$
</div>


---

## Testing
A <strong>hypothesis test</strong> is just a type of decision procedure: we want to decide whether some preconceived notion (our hypothesis) is true. We'll denote our hypothesis with $H$.

We make this actionable by assuming that, if we knew the true value of $\theta$, then we would know whether to accept or reject our hypothesis. This induces a partition of the distribution class, $\mathcal{P}$: $\mathcal{H} \subseteq \mathcal{P}$ is the subset of distributions labelled by values of $\theta$ for which we accept our hypothesis, and $\mathcal{K} \subseteq \mathcal{P}$ (the <strong>class of alternatives</strong>) is the subset for which we reject it. 

If we let $\Theta_H, \Theta_K \subseteq \Theta$ be the corresponding subsets of the parameter space, then we have that:

$$
\mathcal{H} \cup \mathcal{K} = \mathcal{P} \hspace{15mm} \text{and} \hspace{15mm} \Theta_H \cup \Theta_K = \Theta
$$

We denote the decision of accepting $H$ with $d_0$ and the decision of rejecting it with $d_1$. 

### Non-Randomized and Randomized Tests
A <strong>non-randomized test</strong> will assign either $d_0$ or $d_1$ to each value $x \in \mathcal{X}$ (with probability $1$). Thus, we can define $S_0$ and $S_1$ as the subsets of $\mathcal{X}$ that contain the values for which the test assigns $d_0$ and $d_1$, respectively. $S_0$ is the <strong>acceptance region</strong>, and $S_1$ is the <strong>critical region</strong> or the <strong>rejection region</strong>.

A <strong>randomized test</strong> will assigned either $d_1$ or $d_0$ to each value $x \in \mathcal{X}$ with probabilities $\phi(x)$ and $1 - \phi(x)$, respectively. We call $0 \leq \phi(x) \leq 1$ the <strong>critical function</strong> of the test, and it characterizes it completely.

In contrast to a non-randomized test, there is some uncertainty associated with the assignment of decisions; a non-randomized test basically "knows" which decision to assign to all possible realizations of $X$

### Simple and Composite Hypotheses
A class of distributions is called <strong>simple</strong> if it contains only a single distribution. Otherwise, it is called <strong>composite</strong>. The same concept extends to hypotheses. For example, if $\mathcal{H}$ contains only a single distribution, then it is called simple.

### Errors
One can imagine conducting a test and coming to the incorrect conclusion. There are two different ways this can happen: we decide to reject $H$ in favor of $K$ when $H$ is true, or we decide to accept $H$ when it is not true. The former is called a <strong>Type I error</strong>, and the latter is a <strong>Type II error</strong>. 

In a perfect world, we could choose a test procedure that minimizes the probability of both, but this is not possible. Instead, we design our test to be such that the probability of a Type I error is no greater than some value, called the <i>significance level</i>.

<div class="definition">
<strong>Definition (Significance Level).</strong>
<br>
The <i>significance level</i>, $\alpha \in (0, 1)$, of a hypothesis test is a selected value such that:

$$
\begin{equation}
\label{eq:sig-level}
\mathbb{P}_{X \sim P_{\theta}}\left(\delta(X) = d_1 \right) = \mathbb{P}_{X \sim P_{\theta}}\left(X \in S_1 \right) \leq \alpha \hspace{5mm} \forall \theta \in \Theta_H
\end{equation}
$$
</div>

A test that has a signficance level of $\alpha$ will (usually) satisfy:

$$
\underset{\theta \in \Theta_H}{\sup} \left\{ \mathbb{P}_{X \sim P_{\theta}}\left(X \in S_1 \right) \right\} = \alpha
$$

<aside><p>The lefthand side of the above is called the <strong>size</strong> of the hypothesis test.</p></aside>

We then minimize the probability of a Type II error subject to this constraint. That is, we minimize $P_{X \sim P_\theta}\left(\delta(X) = d_0 \right)$ for all $\theta \in \Theta_K$. This is equivalent to maximizing the probability of rejection of the hypothesis for all $\theta \in \Theta_K$.

<div class="definition">
<strong>Definition (Power Function).</strong>
<br>
The <i>power function</i> of a test, $\beta(\theta)$, is a function of $\theta \in \Theta$ and defined as:

$$
\begin{equation}
\label{eq:power-function}
\beta(\theta) = \mathbb{P}_{X \sim P_\theta}\left(\delta(X) = d_1 \right) = \mathbb{P}_{X \sim P_{\theta}}\left(X \in S_1\right)
\end{equation}
$$ 
</div>

Evaluating the power function at a particular value of $\theta \in \Theta_K$ yields the <strong>power of the test against the alternative, $\theta$</strong>. 

### Basic Characteristics
We can characterize tests in a number of ways. Here, we cover a few of the simpler ones. Let $\phi$ be a level $\alpha$ (randomized) test. A test for which the probability of rejecting a false hypothesis is greater than the probability of rejecting a true hypothesis is called <i>unbiased</i>.

<div class="definition">
<strong>Definition (Unbiased).</strong>
The test, $\phi$, is called an <i>unbiased</i> test if it satisfies:

$$
\begin{aligned}
& \beta_\phi(\theta) \leq \alpha & \text{ if }  \theta \in \Theta_H \\
& \beta_{\phi}(\theta) \geq \alpha & \text{ if } \theta \in \Theta_K
\end{aligned}
$$
</div>

A related concept is the <i>exact test</i>, which will control the Type I error.

<div class="definition">
<strong>Definition (Exact).</strong>
The test, $\phi$, is called <i>exact</i> if it satisfies:

$$
\begin{aligned}
\mathbb{E}_{X \sim P_{\theta}} \left[ \phi(x) \right] = \alpha
\end{aligned}
$$

for all $\theta \in \Theta_H$.
</div>

<aside><p>If, instead, $\mathbb{E}_{X \sim P_{\theta}} \left[ \phi(x) \right] \leq \alpha$, then we call $\phi$ <strong>conservative</strong>.</p></aside>



<!-- ## Characteristics

### Uniformly Most Powerful
Consider a randomized test with critical function, $\phi$, and suppose $X \sim P_{\theta}$. The probability of rejection is:

$$
\mathbb{E}_{X \sim P_{\theta}}\left[ \phi(X) \right] = \int \phi(x) d P_{\theta}(x)
$$

A test that maximizes $\beta_(\theta)$ for all $\theta \in \Theta_{K}$ is called <strong>uniformly most powerful (UMP)</strong>. Moreover, if $\mathcal{K}$ is simple, then we can completely specify an ideal test as one satisfying:

$$
\begin{aligned}
&\underset{\phi}{\min} \left\{ \beta_{\phi} (\theta) \right\} = \underset{\phi}{\min} \left\{ \mathbb{E}_{X \sim P_{\theta}}\left[ \phi(X) \right] \right\} & \forall \theta \in \Omega_K \\
&\text{subject to } \mathbb{E}_{X \sim P_{\theta}} \left[ \phi(X) \right] \leq \alpha & \forall \theta \in \Theta_H 
\end{aligned}
$$

This leads us to a fundamental theorem. 

<div class="theorem">
<strong>Neyman-Pearson Lemma.<d-cite key=lehmann2008></d-cite></strong>
{% tabs neyman-pearson %}
{% tab neyman-pearson statement %}
Let $P_0$ and $P_1$ be distributions with respective density functions, $p_0$ and $p_1$, with respect to measure $\mu$. When testing $H: p_0$ against $K: p_1$, there exists a test, $\phi$, with constant $k$ such that, for significance level $\alpha$:

$$
\begin{equation}
\label{eq:n-p}
\mathbb{E}_{X \sim p_0}\left[ \phi(X) \right] = \alpha
\end{equation}
$$

and:

$$
\begin{equation}
\label{eq:n-p2}
\phi(x) = \begin{cases}
1 & \text{when } p_1(x) > k p_0(x) \\
0 & \text{when } p_1(x) < k p_0(x)
\end{cases}
\end{equation}
$$

Furthermore, if $\phi$ satisfies Eq. \eqref{eq:n-p} for some $k$, then it is <strong>most powerful</strong> for testing $p_0$ against $p_1$ at level $\alpha$. 
<br>
In addition, if $\phi$ is most powerful for testing $p_0$ against $p_1$ at level $\alpha$, then it satisfies Eq. \eqref{eq:n-p2} for some $k$ almost everywhere (with respect to $\mu$). It will also satisfy Eq. \eqref{eq:n-p} unless there exists a test of size small than $\alpha$ with power $1$.
{% endtab %}
{% tab neyman-pearson proof %}
Proof to be completed. See 3.2 in Lehmann.<d-cite key=lehmann2008></d-cite>
{% endtab %}
{% endtabs %}
</div>

The Neyman-Pearson Lemma states that a most powerful test (in this particular simple setting) is uniquely (almost everywhere) determined by Eqs. \eqref{eq:n-p} and \eqref{eq:n-p2} when $\{ x \rvert p_1(x) = k p_0(x) \}$ has measure zero.  


## Characteristics
There are many ways to characterize decision rules. Here, we briefly introduce several.

<div class="definition">
<strong>Definition (Dominating).</strong>
<br>
A decision rule, $\delta^*$, is said to <i>dominate</i> (for some particular choice of loss function) another rule, $\delta$, if:

$$
\begin{cases}
R(\theta, \delta^*) \leq R(\theta, \delta) & \forall \theta \\
R(\theta, \delta^*) < R(\theta, \delta) & \text{for some } \theta
\end{cases}
$$

In this case, we call $\delta$ <i>inadmissible</i>.
</div>

Let $\mathcal{C}$ be a class of decision procedures. We call $\mathcal{C}$ <strong>complete</strong> if all of the rules in it have no dominating rule also in $\mathcal{C}$. The class is also called <strong>minimal</strong> if it does not contain a complete subclass, and is called <strong>essentially complete</strong> if, for any $\delta \in \mathcal{C}$, there exists $\delta^* \in \mathcal{C}$ such that $R(\theta, \delta^*) \leq R(\theta, \delta)$ for all $\theta$.  -->
