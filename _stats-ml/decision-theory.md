---
layout: distill
title: Decision Theory
description: A Primer
date: 2026-02-03
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
    - name: Decision Rules
    - name: Decision Principles
      subsections:
            - name: Bayes Principle
            - name: Minimax Principle
    - name: Assorted Results
      subsections:
            - name: An Aside On Sufficiency
            - name: Convex Loss and Rao-Blackwell
bibliography: stats-ml.bib
---

A large part of statistical inference is <i>hypothesis testing</i> in which we wish to learn about some aspect of how a set of random variables behave via a sample. For example, we might be interested in learning about the mean of a random variable or about whether it follows a Gaussian distribution. We can frame testing as a <i>decision</i> we must make, and there is a whole branch of statistical/machine learning literature about <i>decision theory</i>.

In this post, I'm going to cover the very basics (mostly first principles and definitions) to provide a base for some of my other posts. Hopefully this will be a self-contained reference for my future self. I'll mostly use Berger's <i>Statistical Decision Theory</i>.<d-cite key=berger1980></d-cite>

---

## Background
The basic set-up is as follows. We have some unknown quantity, denoted by $\theta$, and we want to reach some conclusion about its state (e.g. if it is positive or not, if it is exactly equal to zero or not, etc.). This quantity can take on many different possible states, and we denote the set of these states with $\Theta$. 

We make <strong>decisions</strong>, or take <strong>actions</strong>, denoted by $d$, which lead us to incur some sort of penalty that measures how bad our decision was. We denote the set of all possible decisions with $\mathcal{D}$. A <strong>loss function</strong> specifies the penalty incurred when we make decision $d$ and $\theta$ is the true state of nature:

$$
\mathcal{L}: \Theta \times \mathcal{D} \rightarrow \mathbb{R}
$$

In general, $\theta$ is referred to as the <strong>state of nature</strong>, but when $\theta$ labels a probability distribution, then we call it a <strong>parameter</strong>. In most statistical settings, we will be dealing with the latter case, and the decisions we make are usually estimates or conclusions about the value of $\theta$.

We will make our decisions about the parameter using a random sample of the random variable, $X$, that follows probability distribution, $P_\theta$, which depends on $\theta$. We observe, say, $n$ realizations of $X$ as $x_1, \dots, x_n$. We denote the <strong>sample space</strong> (set of values $X$ can take) with $\mathcal{X}$. 


### Sufficiency
Statistical inference and testing theory relies upon the idea that we can learn things about the state of nature by looking at data. In a way, we can consider a sample as containing information about the state of nature. By obtaining a sample, we gather evidence to help guide us in our decision making.

<aside><p>One way to quantify this is via the <strong>likelihood function</strong>, which is just the probability density/mass function of the random variable being observed but considered as a function of the parameter.</p></aside>

It is also nice to consider function of a sample that "renders down" all of the information about $\theta$ into a simple value. These are called <i>sufficient statistics</i>.

<div class="definition">
<strong>Definition (Sufficient Statistic).</strong>
<br>
Let $X \sim P_\theta$ where $\theta$ is unknown. A function $T$ of $X$ is called a <i>sufficient statistic</i> for $\theta$ if the conditional distribution of $X$ given $T(X) = t$ is independent of $\theta$ with probability $1$.
</div>

Given a statistic, $T$, we can partition the sample space into subsets that yield the same values of the statistic. Define the range of $T$ as $\mathcal{T} = \{ T(x): x \in \mathcal{X} \}$. This partition is defined as the set of subsets defined by:

$$
\mathcal{X}_t = \{ x \in \mathcal{X}: T(x) = t \} \hspace{5mm} \text{for } t \in \mathcal{T}
$$

<aside><p>If $T$ is sufficient for $\theta$, then the partition is also called sufficient.</p></aside>

Rather than the standard way of thinking about data generation, we can instead consider getting a sufficient statistic value of $t$, then selecting a value of $x \in \mathcal{X}_t$ according to the probability density/mass function of $X$ over this subset. 

---

## Decision Rules
How do we make our decisions? We use a <i>decision rule</i>, which we define now.

<div class="definition">
<strong>Definition (Non-Randomized Decision Rule).</strong>
<br>
For true state of nature, $\theta$, and $X \sim P_\theta$, a <i>non-randomized decision rule</i> is a (measurable) function specifying the decision to make when observing $X = x$:

$$
\delta: \mathcal{X} \rightarrow \mathcal{D}
$$
</div>

Note that the loss function is a function of both the decision taken <i>and</i> the true state of nature. It is therefore useful to define the <i>risk</i> of a decision rule, which averaged over the values of $\theta$.

<div class="definition">
<strong>Definition (Risk Function).</strong>
<br>
For a given loss function and decision rule, the <i>risk function</i> of the decision rule is its expected loss computed with respect to $X$:

$$
R(\theta, \delta) = \mathbb{E}_{X \sim P_\theta}\left[ \mathcal{L}(\theta, \delta(X)) \right]
$$
</div>

The definition of a decision rule we have above is one for a <i>non-randomized</i> rule. That is, when observing a sample, the rule outputs the appropriate decision with probability $1$. In contrast, we can define a <i>randomized</i> decision rule.

<div class="definition">
<strong>Definition (Randomized Decision Rule).</strong>
<br>
For true state of nature, $\theta$, and given an observation of $X \sim P_\theta$, a <i>randomized decision rule</i> is a (measurable) function specifying the probability that each decision, $d$, in $D \subseteq \mathcal{D}$ will be selected:

$$
\delta^*: \mathcal{X} \times D \rightarrow [0, 1]
$$

Thus, when considered as a function of the sample, a randomized decision rule is a probability distribution over $\mathcal{D}$.
</div>

We use the notation $$\delta^*(x, d)$$ to denote the probability that we choose decision $d$ based upon observing $X = x$, and we'll use $$\delta^*(x)$$ to denote the probability distribution generated when observing $X = x$. As such, a non-randomized decision rule can be thought of as a randomized decision rule that is just an indicator function for $$\delta(x) \in D$$. 

Because a randomized decision rule is a function of both the random sample and the subset of the decision space, its loss and risk must be defined in slightly different ways.

<div class="definition">
<strong>Definition (Randomized Decision Rule).</strong>
<br>
Let $\delta^*$ be a randomized decision rule, and let $\delta^*(x)$ be the distribution over $\mathcal{D}$ generated by this decision rule when observing $X = x$. The <i>loss function</i> of $\delta^*$ is the average loss over decisions made with $\delta^*(x)$:

$$
\mathcal{L}(\theta, \delta^*(x)) = \mathbb{E}_{d \sim \delta^*(x)}\left[ \mathcal{L}(\theta, d) \right]
$$

The <i>risk function</i> of $\delta^*$ is then:

$$
R(\theta, \delta^*) = \mathbb{E}_{X \sim P_\theta} \left[ \mathcal{L}(\theta, \delta^*(X))\right]
$$

which is the average loss computed over $X$. 
</div>

<aside><p>It's important to remember that the loss function of a randomized decision rule is dependent upon the sample, while the risk is not as it averages over all possible samples.</p></aside>

Going back to the idea of sufficiency, it is enough to consider only those decision rules that are based on sufficient statistics. Intuitively, this makes sense because a sufficient statistic should contain all of the information about $\theta$ that we could get from a sample.

<!-- #region berger-32 -->
<div class="theorem">
<strong>Theorem.</strong><d-cite key=berger1980></d-cite>
<br>
{% tabs berger-32 %}
{% tab berger-32 statement %}
Let $T$ be a sufficient statistic for $\theta$, and let $$\delta^*_0(x)$$ be a randomized decision rule. Then, under certain conditions, there exists a randomized rule $$\delta_1^*(t)$$ that only depends upon $T(x)$ such that:

$$
R(\theta, \delta_1^*) = R(\theta, \delta_0^*)
$$
{% endtab %}
{% tab berger-32 proof %}
Proof to be completed. See Berger pg. 32. 
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

---

## Decision Principles
Though we've defined two types of decision rules and their risks, we still have not explained how to judge decision rules and decisions. A <strong>decision principle</strong> is sort of like a philosophy by which we choose what makes a decision rule good or bad. 

### Bayes Principle
Suppose we have some idea about what the true state of nature is before we begin our experiment. We formalize this knowledge with a <strong>prior distribution</strong>, denoted by $\pi(\theta)$, which is a distribution over $\Theta$. The <strong>Bayes Principle</strong> states that a decision rule, $\delta_1$, is better than another decision rule, $\delta_2$, if its average risk with respect to $\pi(\theta)$ is smaller. That is:

$$
\mathbb{E}_{\theta \sim \pi(\theta)}\left[ R(\theta, \delta_1) \right] < \mathbb{E}_{\theta \sim \pi(\theta)}\left[ R(\theta, \delta_2)\right]
$$

We can then define the best rule according to this principle.

<div class="definition">
<strong>Definition (Bayes Rule).</strong>
<br>
The <i>Bayes risk</i> of a decision rule, $\delta$, with respect to $\pi$, is defined as:

$$
r(\pi, \delta) = \mathbb{E}_{\theta \sim \pi}\left[ R(\theta, \delta)\right]
$$

A <i>Bayes rule</i> is a decision rule, $\delta^\pi$, that minimizes the Bayes risk for a given prior:

$$
\delta^\pi = \underset{\delta}{\arg \min} \left\{ r(\pi, \delta) \right\}
$$

The <i>Bayes risk</i> of $\pi$, denoted by $r(\pi)$, is the Bayes risk of the Bayes rule (i.e. $r(\pi) = r(\pi, \delta^\pi)$). 
</div>

### Minimax Principle
Another way to rank decision rules is by the penalty they incur in the worst case scenario. The <strong>Minimax Principle</strong> states that a decision rule, $\delta_1$, is better than another decision rule, $\delta_2$, if its risk in the worst possible case is smaller. That is:

$$
\underset{\theta \in \Theta}{\sup} \left\{  R(\theta, \delta_1) \right\} < \underset{\theta \in \Theta}{\sup} \left\{  R(\theta, \delta_2) \right\} 
$$

<aside><p>We usually consider randomized rules here.</p></aside>

We can then define the best rule according to this principle.

<div class="definition">
<strong>Definition (Minimax Rule).</strong>
<br>
A decision rule, $\delta^M$, is called a <i>minimax decision rule</i> if it minimizes the maximum risk over states of nature:

$$
\underset{\theta \in \Theta}{\sup} \left\{ R(\theta, \delta^M) \right\} = \underbrace{\underset{\delta \in \Delta}{\inf} \underset{\theta \in \Theta}{\sup} \left\{ R(\theta, \delta) \right\}}_{(\star)}
$$

where $\Delta$ is the set of all (randomized) rules. We call $(\star)$ the <strong>minimax value</strong>
</div>


---

## Assorted Results
In this section, we'll just go through several important results concerning decision rules and related concepts. 

### An Aside On Sufficiency

First, let $T$ be a sufficient statistic for $\theta$. The conditional distribution of $X$ given $T(X) = t$ gives probability $1$ to $\mathcal{X}_t$ (the subset of $\mathcal{X}$ that yields a value of $t$ for the sufficient statistic) since $X \notin \mathcal{X}_t$ with probability zero if $T(X) = t$. Thus, we can define $f_t(x)$ to be the probability density/mass function of $X$ on $\mathcal{X}_t$ specifically. These functions are independent of $\theta$ by the definition of a sufficient statistic.

Furthermore, using these functions, we can define the conditional expectation of a quantity given $T(X) = t$ as:

$$
\mathbb{E}_{X \rvert t}\left[ h(X) \right] =
\begin{cases}
\int_{\mathcal{X}_t} h(x) f_t(x) dx & \text{continuous}\\
\sum_{x \in \mathcal{X}_t} h(x) f_t(x) & \text{discrete}
\end{cases}
$$


### Convex Loss and Rao-Blackwell
The conditional expectation we defined above is helpful in the following theorems. The first implies that only non-randomized rules are viable when the loss is convex.

<!-- #region berger-35 -->
<div class="theorem">
<strong>Theorem.</strong><d-cite key=berger1980></d-cite>
<br>
{% tabs berger-35 %}
{% tab berger-35 statement %}
Let $\mathcal{D} \subseteq \mathbb{R}^m$ be convex, and assume that $\mathcal{L}(\theta, \mathbf{d})$ be convex for each $\theta \in \Theta$ and $\mathbf{d} \in \mathcal{D}$. Let $\delta^*$ be a randomized decision rule satisfying:

$$
\mathbb{E}_{\mathbf{d} \sim \delta^*(x)}\left[ \rvert \mathbf{d} \rvert \right] < \infty \hspace{5mm} \forall x \in \mathcal{X}
$$

Then, under certain condtions, the non-randomized rule defined by:

$$
\delta(x) = \mathbb{E}_{\mathbf{d} \sim \delta^*(x)} \left[ \mathbf{d} \right]
$$

satisfies:

$$
\mathcal{L}(\theta, \delta(x)) \leq \mathcal{L}(\theta, \delta^*(x)) \hspace{5mm} \forall x, \theta
$$
{% endtab %}
{% tab berger-35 proof %}
Proof to be completed. See Berger pg. 35.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->


We also have the famous <i>Rao-Blackwell Theorem</i>.

<!-- #region berger-36 -->
<div class="theorem">
<strong>Theorem (Rao-Blackwell).</strong><d-cite key=berger1980></d-cite>
<br>
{% tabs berger-36 %}
{% tab berger-36 statement %}
Let $\mathcal{D} \subseteq \mathbb{R}^m$ be convex, and assume that $\mathcal{L}(\theta, \mathbf{d})$ is convex in $\mathbf{d}$ for all $\theta \in \Theta$. Let $T$ be a sufficient statistic for $\theta$, and let $\delta_0(x)$ be a non-randomized decision rule in $\mathbf{D}$. Then, if the following rule exists:

$$
\delta_1(t) = \mathbb{E}_{X \sim P_\theta \rvert t}\left[ \delta_0(X) \right]
$$

it will satisfy:

$$
R(\theta, \delta_0) \geq R(\theta, \delta_1)
$$
{% endtab %}
{% tab berger-36 proof %}
Proof to be completed. See Berger pg. 36. 
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->



