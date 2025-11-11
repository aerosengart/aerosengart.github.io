---
layout: post
title:  "Concentration Inequalities"
date: 20 February 2025
categories: posts
tags: ["theory", "probability"]
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---

In this post, I'm going to take myself through a review of standard concentration inequalities in probability theory with an ultimate goal of exploring empirical Bernstein bounds.

Note: Not all of the proofs are finished/included. I am hoping to find the time to return to thist post and complete them.

---
## Concentration Inequalities

In this section, we'll cover four main concentration inequalities that will aid in understanding later sections of this post. These are covered in Vershynin[^fn-vershynin] and Wasserman[^fn-wasserman] Two very standard concentration inequalities are Markov's inequality and Chebyshev's inequality.

<div class="theorem">
  <body>
  <strong>Markov's Inequality.</strong>
  <br>
  Let $X$ be a non-negative random variable with finite mean (i.e. $\mathbb{E}[X]$ exists). For any choice of $t > 0$, we have:

  $$
  \mathbb{P}(X \geq t) \leq \frac{\mathbb{E}[X]}{t}
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  Let $f(x)$ denote the probability density function of $X$. Since $X \geq 0$, we have:

  $$
  \begin{aligned}
  \mathbb{E}[X] = \int_0^\infty x f(x) dx = \underbrace{\int_0^t x f(x) dx}_{\geq 0} + \int_t^\infty xf(x) dx \geq \int_t^\infty x f(x) \overset{(i)}{\geq} t \int_t^\infty f(x) dx = t \mathbb{P}(X \geq t)
  \end{aligned}
  \nonumber
  $$ 

  In $(i)$, we use the fact that on the interval $[t, \infty]$, $x \geq t$. 
  </details>
  </body>
</div>

<div class="theorem">
  <body>
  <strong>Chebyshev's Inequality.</strong>
  <br>
  Let $X$ be a random variable with finite mean $\mathbb{E}\left[X\right] = \mu$ and finite variance $V(X) = \sigma^2$. Define $Z = \frac{X - \mu}{\sigma}$. We have:

  $$
  \mathbb{P}(\rvert X - \mu \rvert \geq t) \leq \frac{\sigma^2}{t^2} \hspace{10mm} \text{ and } \hspace{10mm} \mathbb{P}(\rvert Z \rvert \geq k) \leq \frac{1}{k^2}
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  Since $\rvert X - \mu \rvert \geq 0$, we can use Markov's Inequality:

  $$
  \mathbb{P}(\rvert X - \mu \rvert \geq t) = \mathbb{P}(\rvert X - \mu \rvert^2 \geq t^2) \leq \frac{\mathbb{E}\left[ (X - \mu)^2 \right]}{t^2} \overset{(i)}{=} \frac{\sigma^2}{t^2}
  \nonumber
  $$

  In $(i)$, we choose $t = k \sigma$.
  </details>
  </body>
</div>

A very important tail bound that will be helpful for the remainder of this post is Hoeffding's inequality, which bounds the deviation of the sum of independent, bounded random variables from their mean. 

<div class="theorem">
  <body>
  <strong>Hoeffding's Inequality.</strong>
  <br>
  Let $X_1, \dots, X_n$ be independent random variables such that $a_i \leq X_i \leq b_i$ almost surely. Let $\mathbb{E}[X_i] =: \mu_i$. For any $t > 0$:

  $$
  \mathbb{P}\left(\sum_{i = 1}^n (X_i - \mu_i) \geq t \right) \leq \exp \left( - \frac{2 t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \right)
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  First, we show the reformulation. Choose some $t > 0$ and let $\delta = \exp \left( - \frac{2t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \right)$. Then:
  $$
  \begin{aligned}
  -\log(\delta) &= \frac{2t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \\
  \implies t^2 &= \frac{\sum_{i = 1}^n (b_i - a_i)^2 \log\left(\frac{1}{\delta} \right)}{2} \\
  \implies t &= \sqrt{\frac{1}{2}\sum_{i = 1}^n (b_i - a_i)^2 \log\left(\frac{1}{\delta}\right)}
  \end{aligned}
  \nonumber
  $$
  Putting the above together with the converse of the first probability expression yields the reformulation.
    <p style="color:red;">TODO: FINISH PROOF</p>
  </details>
  </body>
</div>


An important note made in Vershynin[^fn-vershynin] is the fact that this bound is non-asymptotic; that is, it holds for any choice of $n$. The bound just becomes tighter as the sample size grows.

Unfortunately, Hoeffding's inequality is independent of the variances of the random variables in question. We can imagine that a better bound could be obtained when the variances are small since we know they should remain close to the means. Bernstein inequalities aim to do this by trading sharpness in cases with large variance for improved bounds in cases with small variance. 

There are several similar Berstein inequalities, one of which I've written below:

<div class="theorem">
  <body>
  <strong>A Bernstein Inequality.</strong>
  <br>
  Let $X_1, \dots, X_n$ be i.i.d. bounded random variables (i.e. their support is $[a, b]$ for finite $a, b$). Let $\mathbb{E}[X_i] =: \mu$ and $\mathbb{E}\left[ (X - \mu)^2 \right] =: \sigma^2$. For any $t \geq $, We have:

  $$
  \mathbb{P}\left( \bigg\rvert \frac{1}{n} \sum_{i = 1}^n X_i - \mu \bigg\rvert \geq t \right) \leq 2 \exp\left( -\frac{n t^2}{2(\sigma^2 + (b-a) t)} \right)
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
    <p style="color:red;">TODO: FINISH PROOF</p>
  </details>
  </body>
</div>

Bennett's inequality is similar to the Bernstein inequalities in that it bounds the sum of independent random variables in a variance-dependent fashion. 

<div id="theorem-bennett"></div>
<div class="theorem">
  <body>
  <strong>Bennett's Inequality.<span markdown="1">[^fn-bennett]</span></strong>
  <br>
  Let $X_1, \dots, X_n$ be independent random variables with finite mean $\mathbb{E}[X_i] =: \mu$. Let $\sigma^2 := \sum_{i = 1}^n \mathbb{E}\left[ (X_i - \mathbb{E}[X_i])^2\right]$. For all $i$, suppose $\rvert X_i - \mathbb{E}[X_i] \rvert \leq a$ almost surely. For any $t \geq $, We have:

  $$
  \mathbb{P}\left( \sum_{i = 1}^n (X_i - \mathbb{E}[X_i]) \geq t \right) \leq \exp\left( -\frac{\sigma^2}{a^2} \left[ \left(1 + \frac{at}{\sigma^2}\right) \log\left(1 + \frac{at}{\sigma^2}\right) - \frac{at}{\sigma^2} \right] \right)
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
    <p style="color:red;">TODO: FINISH PROOF</p>
  </details>
  </body>
</div>

We also have the following variation on Bennett's inequality from Audibert et al. (2009).

<div class="theorem-bennett2"></div>
<div class="theorem">
  <body>
  <strong>Lemma 5 (Audibert et al.<span markdown="1">[^fn-audibert]</span>).</strong>
  <br>
  Let $X \in \mathbb{R}$ be a random variable such that $X \leq s$ for some $s \in \mathbb{R}$ almost surely. Denote $\mu := \mathbb{E}[X]$ and $s_+ := s \vee 0$. Define $X_1, \dots, X_n$ as independent and identically distributed copies of $X$. Let $\bar{U}_m := \frac{1}{m} \sum_{i = 1}^m U_i$. Then, for any $\delta > 0$ and simultaneously for all $1 \leq m \leq n$, with probability at least $ 1 - \exp(-\delta)$:

  $$
  \begin{aligned}
  &m(\bar{U}_m - \mu) \leq \sqrt{2 n \mathbb{E}[U^2]\delta} + \frac{s_+ \delta}{3}\\
  &m(\bar{U}_m - \mu) \leq \sqrt{2 n \text{Var}(U)\delta} + \frac{(s - \mu)\delta}{3}
  \end{aligned}
  \nonumber
  $$

  A proof can be found in the publication.
  </body>
</div>

---

## Empirical Bernstein bounds

Equipped with a solid understanding of the above inequalities, we can dive into "Empirical Bernstein Bounds and Sample Variance Penalization" by Maurer and Pontil.[^fn-maurer-b]

### Background

We'll assume we have some random variable, $X$, under study that comes from some sample space, $\mathcal{X}$. Let $\mu$ denote the law of $X$, and let $\ell$ be some loss function. We'll use $\ell_h(X)$ to denote the loss of hypothesis $h$ for observed value $X$. 

<div id="risk"></div>
<div class="definition">
  <body>
  <strong>Definition (Risk).</strong>
  <br>
  The <i>risk</i> is defined by the expected loss with respect to $\mu$, $R(\ell_h, \mu) := \mathbb{E}_{X \sim \mu}\left[ \ell_h(X) \right]$. For an i.i.d. sample, $\mathbf{X} := (X_1, \dots, X_n)$, we can estimate the risk of hypothesis $h$ with the <i>empirical risk</i>, which is the empirical mean of the loss, $R_n(\ell_h, \mathbf{X}) := \frac{1}{n} \sum_{i = 1}^n \ell_h(X_i)$. 
  </body>
</div>

Assuming that the boundedness conditions are satisfied by $\ell_h$, we can use Hoeffding's inequality to get, with probability $1 - \delta$ (for any $\delta > 0$):

$$
\big\rvert R_n(\ell_h, \mathbf{X}) - R(\ell_h, \mu) \big \rvert \leq \sqrt{\frac{1}{2}\sum_{i = 1}^n (b_i - a_i)^2 \log\left(\frac{2}{\delta} \right)}
\nonumber
$$

Recall that machine learning is often concerned with hypothesis selection (e.g. model fitting, parameter estimation, etc.) via <i>empirical risk minimization</i>, which is a procedure that picks the hypothesis that minimizes $R_n(f, \mu)$ over all $f \in \mathcal{F}$. The basic idea is that if the empirical risk is similar to the true risk, then we will end up picking the best hypothesis. If we have some guarantee on the difference between the empirical and true risk, then learning will be possible with large enough sample sizes!

In later proofs, we will need to use some measure of the complexity of a function class. For a fixed $\epsilon > 0$, function class $\mathcal{F}$, and sample size $n$, we'll define the <i>growth function</i>, which is an extension of covering numbers. However, we'll need some intermediate definitions first.

<div id="epsilon-net"></div>
<div class="definition">
  <body>
  <strong>Definition ($\epsilon$-Net).</strong>
  <br>
  Let $T$ be some set, let $d$ be some distance metric, and let $(T, d)$ be the metric space created by their pairing. Pick $\epsilon > 0$ and some subset $K \subset T$. An <i>$\epsilon$-net of $K$</i> is defined as any subset $\mathcal{N} \subseteq K$ such that:

  $$
  \forall x \in K, \hspace{3mm} \exists x_0 \in \mathcal{N} \hspace{3mm} \text{such that} \hspace{3mm} d(x, x_0) \leq \epsilon
  \nonumber
  $$
  </body>
</div>

Intuitively, a subset $\mathcal{N}$ is an $\epsilon$-net of $K$ if all the points in $K$ can be contained in a ball of radius $\epsilon$ centered at some points in $\mathcal{N}$ (i.e. $K$ can be covered).

<div id="covering-number"></div>
<div class="definition">
  <body>
  <strong>Definition (Covering Number).</strong>
  <br>
  The <i>covering number of subset $K$</i> is defined as the smallest cardinality of any $\epsilon$-net of $K$. It will depend on the distance metric, $K$, and the choice of $\epsilon$, so we denote it by $\mathcal{N}(\epsilon, K, d)$. 
  </body>
</div>

Vershynin provides a very good alternative explanation of the covering number: the smallest number of balls with centers in $K$ and radii of $\epsilon$ that completely cover $K$. 

<div id="growth-function"></div>
<div class="definition">
  <body>
  <strong>Definition (Growth Function).</strong>
  <br>
  For function class $\mathcal{F}$, the <i>growth function</i> is defined by:
  $$
  \mathcal{N}_{\infty}(\epsilon, \mathcal{F}, n) = \underset{\mathbf{x} \in X^n}{\sup} \mathcal{N}\left( \epsilon, \mathcal{F}(\mathbf{x}), \rvert \rvert \cdot \rvert \rvert_{\infty} \right)
  \nonumber
  $$
  where $\mathcal{F}(\mathbf{X}) = \{ (f(x_1), \dots, f(x_n)) \rvert f \in \mathcal{F} \} \subseteq \mathbb{R}^n$ and $\rvert \rvert f \rvert \rvert_\infty = \underset{X \in \mathcal{X}}{\sup} \rvert f(X) \rvert$ is the <i>uniform norm</i>.
  </body>
</div>

The growth function is the least upper bound (over all samples of size $n$) on the covering number of the set of all hypotheses in $\mathcal{F}$ applied to the sample $\mathbf{x}$ with respect to the distance metric associated with $\rvert \rvert \cdot \rvert \rvert_{\infty}$. 

The growth function is also called the <i>shattering number</i> or the <i>shattering coefficient</i> of a set, but growth function is the terminology used in the original work of Vapnik and Chervonenkis. 

---

### Results

This section will rely on common notation. Let $[0, 1]^n$ be the $n$-dimensional unit cube, and let $$\mathbf{x} := (x_1, \dots, x_n) \in [0, 1]^n$$ be a point in the cube. Denote and define the sample mean and variance with $$P_n(\mathbf{x}) := \frac{1}{n}\sum_{i = 1}^n x_i$$ and $$V_n(\mathbf{x}) := \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j = 1}^n \frac{(x_i - x_j)^2}{2}$$, respectively. 

For a function $f$, we will use $$f(\mathbf{x}) := (f(x_1), \dots, f(x_n))$$ to denote the function applied element-wise to the sample vector $\mathbf{x}$. Similarly, we will use $$P_n(f, \mathbf{x}) := P_n(f(\mathbf{x}))$$ and $$V_n(f, \mathbf{x}) := V_n(f(\mathbf{x}))$$.

For distribution $\mu$, the $n$-fold product is $\mu^n$. General product measures are denoted with $\times$ or $\prod$. For $X \sim \mu$ and function $f$, $$P(f, \mu) := \mathbb{E}_{X \sim \mu}[f(X)]$$ and $$V(f, \mu) := V_{X \sim \mu} f(X)$$ which will be notated more simply as $\mathbb{E}[f(X)]$ and $V(f(X))$, respectively.

#### Inequalities 

Notice that, although the Bernstein inequalities and Bennett's inequality can provide better bounds, they require knowledge of the unknown variance. Maurer and Pontil first derive an inequality similar to those of Bernstein and Bennett that only requires calculation of the sample variance.

<div id="theorem-11"></div>
<div class="theorem">
  <body>
  <strong>Empirical Bernstein Bound (Theorem 11 Maurer and Pontil<span markdown="1">[^fn-maurer-b]</span>).</strong>
  <br>
  Let $\mathbf{X} = (X_1, \dots, X_n)$ be a vector with independent random variables taking values in $[0, 1]$. Let $\delta > 0$. With probability at least $1 - \delta$ in $\mathbf{X}$:
  $$
  \mathbb{E}\left[ P_n(\mathbf{X}) \right] \leq P_n(\mathbf{X}) + \sqrt{\frac{2V_n(\mathbf{X}) \log(2/\delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)}
  \nonumber
  $$
  where we recall that $P_n(\mathbf{X}) := \frac{1}{n} \sum_{i = 1}^n X_i$, the sample average.

  <details>
  <summary>Proof.</summary>
    We'll first write $W := \frac{1}{n} \sum_{i = 1}^n V(X_i)$, the average variance of each $X_i$ (since they aren't necessarily identically distributed). Note that $W \leq \mathbb{E}[V_n(\mathbf{X})]$.
    <details>
      <summary>Proof.</summary>
      $$
      \begin{aligned}
      W &= \frac{1}{n} \sum_{i = 1}^n V(X_i) \\
        &= \frac{1}{n} \sum_{i = 1}^n \mathbb{E}\left[ (X_i - \mathbb{E}[X_i])^2 \right] \\
        &\leq \frac{1}{n} \sum_{i = 1}^n \mathbb{E}\left[ (X_i - \mathbb{E}[X_i])^2 \right] + \frac{1}{2n(n-1)}\sum_{i = 1}^n \sum_{i \neq j} \left( \mathbb{E}[X_i] - \mathbb{E}[X_j] \right)^2 \hspace{10mm} \text{(Adding a positive quantity)} \\
        &= \frac{1}{n}\sum_{i = 1}^n \mathbb{E}\left[ X_i^2 + (\mathbb{E}[X_i])^2 - 2 X_i \mathbb{E}[X_i] \right]
          + \frac{1}{n(n-1)}\sum_{i = 1}^n \sum_{j < i}\left( (\mathbb{E}[X_i])^2 + (\mathbb{E}[X_j])^2 - 2 \mathbb{E}[X_i]\mathbb{E}[X_j] \right) \\
        &= \frac{1}{n} \sum_{i = 1}^n \left( \mathbb{E}\left[ X_i^2 \right] - (\mathbb{E}[X_i])^2 \right) + \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j < i} \left[(\mathbb{E}[X_i])^2 + (\mathbb{E}[X_j])^2 - 2\mathbb{E}[X_i]\mathbb{E}[X_j] \right]  \\
        &= \frac{1}{n} \sum_{i = 1}^n \mathbb{E}\left[ X_i^2 \right]  - \frac{1}{n} \sum_{i = 1}^n (\mathbb{E}[X_i])^2  + \frac{n-1}{n(n-1)} \sum_{i = 1}^n  (\mathbb{E}[X_i])^2 + \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j < i} \left[ (\mathbb{E}[X_j])^2 - 2\mathbb{E}[X_i]\mathbb{E}[X_j] \right] \\
        &= \frac{1}{n} \sum_{i = 1}^n  \mathbb{E}\left[ X_i^2 \right] + \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j < i} \left[ (\mathbb{E}[X_j])^2 - 2\mathbb{E}[X_i]\mathbb{E}[X_j] \right] \\
        &= \sum_{i = 1}^n \left( \frac{1}{n} \mathbb{E}\left[ X_i^2 \right] + \frac{1}{n(n-1)} \sum_{j < i} \left[ (\mathbb{E}[X_j])^2 - 2\mathbb{E}[X_i]\mathbb{E}[X_j] \right]\right) \\
        &= \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j < i} \left( \mathbb{E}\left[ X_i^2 \right] + (\mathbb{E}[X_j])^2 - 2\mathbb{E}[X_i]\mathbb{E}[X_j]\right) \\
        &= \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j < i}\left( \mathbb{E}[(X_i - X_j)^2] \right) \\
        &= \frac{1}{2n(n-1)} \sum_{i = 1}^n \sum_{j \neq i}\left( \mathbb{E}[(X_i - X_j)^2] \right) \\
        &= \frac{1}{2n(n-1)} \sum_{i = 1}^n \sum_{j = 1}^n \left( \mathbb{E}[(X_i - X_j)^2] \right) \hspace{10mm} X_i - X_j = 0 \text{ for } i = j \\
        &= \mathbb{E}[V_n(\mathbf{X})]
        \end{aligned}
        \nonumber
        $$
    </details>
      Since $X_i \in [0, 1]$ for all $i$, we can use a derivative of Bennett's inequality (see <a href="#bennett-derivative">Corollary 2</a>) to get:
      $$
      \mathbb{E}[P_n(\mathbf{X})] \leq P_n(\mathbf{X}) + \sqrt{\frac{2 W \log \left(\frac{1}{\delta}\right)}{n}} + \frac{\log \left( \frac{1}{\delta} \right)}{3n} \leq P_n(\mathbf{X}) + \sqrt{\frac{2 \mathbb{E}[V_n(\mathbf{X})] \log \left(\frac{1}{\delta}\right)}{n}} + \frac{\log \left( \frac{1}{\delta} \right)}{3n}
      $$
    Looking at <a href="#theorem-10">Theorem 10</a>, we also have that $\sqrt{V_n(\mathbf{X})} - \sqrt{\frac{2\log\left(\frac{1}{\delta}\right)}{n-1}} > \sqrt{\mathbb{E}[V_n(\mathbf{X})]}$ with probability $1 - \delta$ for any $\delta > 0$. 
      <br> 
      Since the probability of having both statements above satisfied is less than or equal to the sum of their individual probabilities (via the union bound), we can pick $\delta = 2\epsilon$ for any $\epsilon > 0$ to ensure that, with probability at least $1 - 2\epsilon = 1 - \delta$:
      $$
      \begin{aligned}
      P_n(\mathbf{X}) + 2\sqrt{\mathbb{E}[V_n(\mathbf{X})]} \sqrt{\frac{\log\left(\frac{1}{\epsilon}\right)}{2n}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} 
      &\leq P_n(\mathbf{X}) + 2\left(\sqrt{V_n(\mathbf{X})} - \sqrt{\frac{2 \log\left(\frac{1}{\epsilon}\right)}{n-1}} \right)\sqrt{\frac{\log\left(\frac{1}{\epsilon}\right)}{2n}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
      &= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - \sqrt{\frac{4\log^2\left(\frac{1}{\epsilon}\right)}{n(n-1)}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
      &= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - \frac{2\log\left(\frac{1}{\epsilon}\right)}{\sqrt{n(n-1)}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
      &\leq P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - \frac{2\log\left(\frac{1}{\epsilon}\right)}{\sqrt{n^2}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} & \left( n^2 \geq n(n-1) \implies -\frac{x}{n^{2}} \geq -\frac{x}{n(n-1)} \right)\\
      &= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - \frac{6\log\left(\frac{1}{\epsilon}\right)}{3n} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
      &= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - \frac{5\log\left(\frac{1}{\epsilon}\right)}{3n} \\
      &= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{2}{\delta}\right)}{n}} - \frac{5\log\left(\frac{2}{\delta}\right)}{3n} \hspace{10mm} \epsilon = \frac{\delta}{2}
      \end{aligned}
      \nonumber
      $$
    </details>

  If we want a probability inequality that matches the form of Markov's, Chebyshev's, and Hoeffding's, we can rearrange some things to get:
  $$
  \mathbb{P}\left(\mathbb{E}[P_n(\mathbf{X})]- P_n(\mathbf{X}) \geq t \right) \leq 2
  \exp\left(- 2\left[ \sqrt{\frac{3(n-1)}{7} \left( t + \sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)} - \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}} \right]\right)
  $$

  <details>
  <summary>Proof.</summary>
  We have that:
  $$
  \begin{aligned}
  &\mathbb{P}\left(  \mathbb{E}\left[ \frac{1}{n}\sum_{i = 1}^n X_i \right]  -  \frac{1}{n} \sum_{i = 1}^n X_i  \leq \sqrt{\frac{2 V_n(\mathbf{X})\log(2/\delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)} \right) \geq 1 - \delta \\
  \implies
  &1 - \mathbb{P}\left(  \mathbb{E}\left[ \frac{1}{n}\sum_{i = 1}^n X_i \right]  -  \frac{1}{n} \sum_{i = 1}^n X_i  \leq \sqrt{\frac{2 V_n(\mathbf{X})\log(2/\delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)} \right) \leq \delta \\
  \implies
  &\mathbb{P}\left(  \mathbb{E}\left[ \frac{1}{n}\sum_{i = 1}^n X_i \right]  -  \frac{1}{n} \sum_{i = 1}^n X_i  \geq \sqrt{\frac{2 V_n(\mathbf{X})\log(2/\delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)} \right) \leq \delta
  \end{aligned}
  \nonumber
  $$
  We then set the righthand side of the first inequality equal to $t$ and solve for $\delta$:
  $$
  \begin{aligned}
    t &= \sqrt{\frac{2 V_n(\mathbf{X})\log(2/\delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)} \\
    t &= 2 \sqrt{\log(2/\delta)} \sqrt{\frac{V_n(\mathbf{X})}{2n}} + \frac{7}{3(n-1)}\log(2/\delta) \\
    \frac{3(n-1)}{7} t &= 2 \frac{3(n-1)}{7} \sqrt{\log(2/\delta)} \sqrt{\frac{V_n(\mathbf{X})}{2n}} + \log(2/\delta)  \\
    \frac{3(n-1)}{7} t + \frac{3(n - 1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}} &= \frac{3(n - 1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}} + 2\frac{3(n-1)}{7} \sqrt{\log(2/\delta)} \sqrt{\frac{V_n(\mathbf{X})}{2n}} + \log(2/\delta)  \\
    \frac{3(n-1)}{7} t + \frac{3(n - 1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}} &= \left(\sqrt{\log(2/\delta)} + \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)^2 \\
    \sqrt{\log(2/\delta)} + \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}}  &= \sqrt{\frac{3(n-1)}{7} \left( t + \sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)} \\
    \sqrt{\log(2/\delta)} &= \sqrt{\frac{3(n-1)}{7} \left( t + \sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)} - \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}} \\
    \log(2/\delta) &= \left(\sqrt{\frac{3(n-1)}{7} \left( t + \sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)} - \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)^2 \\
    \frac{2}{\delta} &= \exp\left(2\left(\sqrt{\frac{3(n-1)}{7} \left( t + \sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)} - \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)\right) \\
    \delta &= 2\exp\left(- 2\left[ \sqrt{\frac{3(n-1)}{7} \left( t + \sqrt{\frac{V_n(\mathbf{X})}{2n}}\right)} - \frac{3(n-1)}{7}\sqrt{\frac{V_n(\mathbf{X})}{2n}} \right]\right)
  \end{aligned}
  \nonumber
  $$
  </details>
</body>
</div>

  If the $X_i \in [a, b]$ instead, then with probability at least $1 - \delta$:

  $$
 \mathbb{E}[P_n(\mathbf{X})]- P_n(\mathbf{X}) \leq \sqrt{\frac{2V_n(\mathbf{X}) \log(2 / \delta)}{n}} + \frac{7(b-a)\log(2/\delta)}{3(n-1)} 
  \nonumber
  $$

  <details>
    <summary>Proof.</summary>
      Suppose that the $X_i \in [a, b]$. We can repeat the prior argument with $Z_i := \frac{X_i - a}{b - a}$ so that $Z_i \in [0,1]$ and then rearrange the result to get a bound for $X_i$.
      $$
      \begin{aligned}
      V_n(\mathbf{Z}) &= \frac{1}{2n(n-1)}\sum_{i = 1}^n \sum_{j = 1}^n (Z_i - Z_j)^2 \\
      &= \frac{1}{2n(n-1)}\sum_{i = 1}^n \sum_{j = 1}^n \left( \frac{X_i - a}{b-a} - \frac{X_j - a}{b-a}\right)^2 \\
      &= \frac{1}{2n(n-1)(b-a)^2}\sum_{i = 1}^n \sum_{j = 1}^n \left( X_i - X_j \right)^2 \\
      &= \frac{1}{(b-a)^2} V_n(\mathbf{X})
      \end{aligned}
      \nonumber
      $$
      Applying Theorem 11 to the $Z_i$'s, we get that, with probability at least $1 - \delta$:
      $$
      \begin{aligned}
      &\mathbb{E}\left[ \frac{1}{n}\sum_{i = 1}^n \frac{X_i - a}{b-a} \right] \leq \frac{1}{n} \sum_{i = 1}^n \frac{X_i - a}{b-a} + \sqrt{\frac{2V_n(\mathbf{Z}) \log(2 / \delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)} \\
      \implies &\frac{1}{n(b-a)}\sum_{i = 1}^n \left[ \mathbb{E}\left[ X_i \right] - a \right] \leq \frac{1}{n(b-a)}\sum_{i = 1}^n (X_i - a) + \sqrt{\frac{2V_n(\mathbf{Z}) \log(2 / \delta)}{n}} + \frac{7 \log(2/\delta)}{3(n-1)} \\
      \implies &\frac{1}{n(b-a)}\sum_{i = 1}^n \left[  \mathbb{E}\left[ X_i \right] - a - X_i + a\right] \leq \sqrt{\frac{2V_n(\mathbf{Z}) \log(2 / \delta)}{n}} + \frac{7\log(2/\delta)}{3(n-1)} \\
      \implies &\frac{1}{n}\sum_{i = 1}^n \left[ \mathbb{E}\left[ X_i \right] - X_i \right] \leq (b-a)\sqrt{\frac{2V_n(\mathbf{Z}) \log(2 / \delta)}{n}} + \frac{7 (b-a)\log(2/\delta)}{3(n-1)} \\
      \implies  & \frac{1}{n}\sum_{i = 1}^n \left[ \mathbb{E}\left[ X_i \right] - X_i \right] \leq (b-a)\sqrt{\frac{2V_n(\mathbf{X}) \log(2 / \delta)}{(b-a)^2n}} + \frac{7 (b-a)\log(2/\delta)}{3(n-1)} \\
      \implies  &\frac{1}{n}\sum_{i = 1}^n \left[ \mathbb{E}\left[ X_i \right] - X_i \right] \leq \sqrt{\frac{2V_n(\mathbf{X}) \log(2 / \delta)}{n}} + \frac{7 (b-a)\log(2/\delta)}{3(n-1)} \\
      \end{aligned}
      \nonumber
      $$
  </details>

  We can also write Theorem 11 in a slightly different way. Again assume $X_i \in [a, b]$. Then, for any $t > 0$:

  $$
  \mathbb{P}\left( \bigg\rvert \sum_{i = 1}^n \left[X_i - \mathbb{E}[X_i] \right]\bigg\rvert \geq t\right) \leq 4\exp\left(- 2\left(\sqrt{\frac{3(n-1)t}{7 (b-a) n} + \left(\frac{3(n-1)\sqrt{\frac{n V_n(\mathbf{X})}{2}}}{7 (b-a) n}\right)^2}  - \frac{3(n-1)\sqrt{\frac{nV_n(\mathbf{X})}{2}}}{7(b-a)n} \right)\right)
  \nonumber
  $$
  
<details>
  <summary>Details Of Restatement.</summary>
  Theorem 11 states that, with probability at least $1 - \delta$:
  $$
  \mathbb{E}[P_n(\mathbf{X})]- P_n(\mathbf{X}) \leq \sqrt{\frac{2V_n(\mathbf{X}) \log(2 / \delta)}{n}} + \frac{7 (b-a)\log(2/\delta)}{3(n-1)}     
  \nonumber
  $$
  We set $t$ equal to righthand side of the above and solve for $\delta$:
  $$
  \begin{aligned}
  &t =  \sqrt{\frac{2 V_n(\mathbf{X}) \log(2/\delta)}{n}} + \frac{7 (b-a) \log(2/\delta)}{3(n-1)} \\
  \implies &t = 2\sqrt{\log(2/\delta)}\sqrt{\frac{V_n(\mathbf{X})}{2n}} + \frac{7 (b-a) \log(2/\delta)}{3(n-1)} \\
  \implies &\frac{3(n-1)t}{7 (b-a)} = 2\sqrt{\log(2/\delta)}\left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right) + \log(2/\delta) \\
  \implies &\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2= 2\sqrt{\log(2/\delta)}\left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right) + \log(2/\delta) + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2 \\
  \implies &\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2 = \left(\sqrt{\log(2/\delta)} + \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)}\right)^2 \\
  \implies &\sqrt{\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2} = \sqrt{\log(2/\delta)} + \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)} \\
  \implies &\sqrt{\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2}  - \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)} = \sqrt{\log(2/\delta)} \\
  \implies &\left(\sqrt{\frac{3(n-1)t}{7(b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2}  - \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)} \right)^2 = \log(2/\delta) \\
  \implies &2 / \delta = \exp \left( \left(\sqrt{\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2}  - \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)} \right)^2  \right) \\
  \implies &\delta = 2\exp\left(- 2\left(\sqrt{\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{n_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2}  - \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)} \right)\right) 
  \end{aligned}
  \nonumber
  $$
  Thus, for any $t > 0$:
  $$
  \mathbb{P}\left( \mathbb{E}[P_n(\mathbf{X})] - P_n(\mathbf{X}) \geq t\right) \leq 2\exp\left(- 2\left(\sqrt{\frac{3(n-1)t}{7 (b-a)} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a)}\right)^2}  - \frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7(b-a)} \right)\right)
  \nonumber
  $$
</details>


To make the statement two-sided, we just need to bound:

$$
\mathbb{P}\left(\mathbb{E}[P_n(\mathbf{X})] - P_n(\mathbf{X}) \leq -t \right)
\nonumber
$$

From the proof of <a href="#theorem-11">Theorem 11</a>:

$$
\mathbb{E}[P_n(\mathbf{X})] \leq P_n(\mathbf{X}) + \sqrt{\frac{2 W \log \left(\frac{1}{\delta}\right)}{n}} + \frac{\log \left( \frac{1}{\delta} \right)}{3n} \leq P_n(\mathbf{X}) + \sqrt{\frac{2 \mathbb{E}[V_n(\mathbf{X})] \log \left(\frac{1}{\delta}\right)}{n}} + \frac{\log \left( \frac{1}{\delta} \right)}{3n}
\nonumber
$$

Looking at $(b)$ in <a href="#theorem-10">Theorem 10</a>, we also have that:

$$
\sqrt{V_n(\mathbf{X})} > \sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{2\log\left(\frac{1}{\delta}\right)}{n-1}}
\implies
\sqrt{\mathbb{E}[V_n(\mathbf{X})]} < \sqrt{V_n(\mathbf{X})}  - \sqrt{\frac{2\log\left(\frac{1}{\delta}\right)}{n-1}}
\nonumber
$$

with probability $1 - \delta$ for any $\delta > 0$. Since the probability of having both above statements satisfied is less than or equal to the sum of their individual probabilities (via the union bound), we can pick $\delta = 2\epsilon$ for any $\epsilon > 0$ to ensure that, with probability at least $1 - 2\epsilon = 1 - \delta$:
  
$$
\begin{aligned}
P_n(\mathbf{X}) + 2\sqrt{\mathbb{E}[V_n(\mathbf{X})]} \sqrt{\frac{\log\left(\frac{1}{\epsilon}\right)}{2n}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} 
&\leq P_n(\mathbf{X}) + 2\left(\sqrt{V_n(\mathbf{X})} - \sqrt{\frac{2 \log\left(\frac{1}{\epsilon}\right)}{n-1}} \right)\sqrt{\frac{\log\left(\frac{1}{\epsilon}\right)}{2n}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
&= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - 2\sqrt{\frac{\log^2\left(\frac{1}{\epsilon}\right)}{n(n-1)}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
&= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} -\frac{2\log\left(\frac{1}{\epsilon}\right)}{\sqrt{n(n-1)}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
&\leq P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} - \frac{2\log\left(\frac{1}{\epsilon}\right)}{\sqrt{n^2}} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \hspace{10mm} n^2 \geq n(n-1) \implies n^{-1} \leq (n(n-1))^{-1/2} \\
&= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} +\frac{6\log\left(\frac{1}{\epsilon}\right)}{3n} + \frac{\log\left(\frac{1}{\epsilon}\right)}{3n} \\
&= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{1}{\epsilon}\right)}{n}} +\frac{7\log\left(\frac{1}{\epsilon}\right)}{3n} \\
&= P_n(\mathbf{X}) + \sqrt{\frac{2 V_n(\mathbf{X}) \log\left(\frac{2}{\delta}\right)}{n}} +\frac{7\log\left(\frac{2}{\delta}\right)}{3n} \hspace{10mm} \epsilon = \frac{\delta}{2}
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
&\mathbb{P}\left( \mathbb{E}[P_n(\mathbf{X})] \leq P_n(\mathbf{X}) + \sqrt{\frac{2V_n(\mathbf{X}) \log(2/\delta)}{n}} + \frac{7 \log(2/\delta)}{3n} \right) \geq 1 - \delta \\
\implies 
\end{aligned}
\nonumber
$$



---

### Additional Results

The authors then extend this result to finite function classes and present a uniform bound for particular function classes. See the original paper for proofs.

<div class="theorem">
  <body>
  <strong>Empirical Bernstein Bound For Finite Function Class (Corollary 5 Maurer and Pontil<span markdown="1">[^fn-maurer-b]</span>).</strong>
  <br>
  Let $X$ be a random variable taking values in $\mathcal{X}$ and with distribution $\mu$. Let $\mathbf{X} = (X_1, \dots, X_n)$ be an i.i.d. sample with distribution $\mu$. Let $\mathcal{F}$ be a class of hypotheses $f: \mathcal{X} \rightarrow [0, 1]$ such that $\rvert \mathcal{F} \rvert < \infty$. For any choice of $\delta > 0$ and $n \geq 2$, with probability at least $1 - \delta$ for all $f \in \mathcal{F}$:

  $$
  R(f, \mu) - R_n(f, \mathbf{X}) \leq \sqrt{\frac{2 V_n(f, \mathbf{X}) \log(2 \rvert \mathcal{F} \rvert/ \delta)}{n}} + \frac{7 \log(2 \rvert \mathcal{F} \rvert / \delta)}{3(n - 1)}
  \nonumber
  $$

  with $V_n(f, \mathbf{X}) := V_n(f(X_1), \dots, f(X_n))$.
  </body>
</div>

<div class="theorem">
  <body>
  <strong>Uniform Empirical Bernstein Bound (Theorem 6 Maurer and Pontil<span markdown="1">[^fn-maurer-b]</span>).</strong>
  <br>
  Let $X$ be a random variable taking values in $\mathcal{X}$ and with distribution $\mu$. Let $\mathbf{X} = (X_1, \dots, X_n)$ be an i.i.d. sample with distribution $\mu$. Let $\mathcal{F}$ be the class of hypotheses $f: \mathcal{X} \rightarrow [0, 1]$. Pick $\delta \in (0, 1)$, $n \geq 16$ and let $\mathcal{M}(n) = 10 \times \mathcal{N}_{\infty}\left( \frac{1}{n}, \mathcal{F}, 2n \right)$. With probability at least $1 - \delta$, for all $f \in \mathcal{F}$:

  $$
  R(f, \mu) - R_n(f, \mathbf{X}) \leq \sqrt{\frac{18 V_n(f, \mathbf{X}) \log(\mathcal{M}(n)/\delta)}{n}} + \frac{15 \log(\mathcal{M}(n) / \delta)}{n - 1}
  \nonumber
  $$
  </body>
</div>

The previous theorem essentially bounds (with high probability) the difference between the empirical and true risk of all hypotheses ins the defined function class as a function of the sample variance and the function class's complexity (as measured by the growth function).

---

### Helper Proofs

Some additional theorems/lemmata are necessary for the proofs of the main results. 

<div id="theorem-7"></div>
<div class="theorem">
  <body>
  <strong>Concentration Of Self-Bounding Random Variables (Theorem 7 Maurer and Pontil<span markdown="1">[^fn-maurer-b]</span>).</strong>
  Define $\mathbf{X} = (X_1, \dots, X_n)$ to be a vector of independent random variables taking on values in $\mathcal{X}$. For any $1 \leq k \leq n$ and any $y \in \mathcal{X}$, let $\mathbf{X}_{y,k} := (X_1, \dots, X_{k-1}, y, X_{k+1}, \dots, X_n)$ denote the replacement of the $k$-th cooordinate of $\mathbf{X}$ with $y$.

  Suppose we have some $a \geq 1$ and that $Z(\mathbf{X})$ satisfies, almost surely:

  $$
  \begin{aligned}
  (a) &\hspace{5mm} 1 \geq Z(\mathbf{X}) - \underset{y \in \mathcal{X}}{\inf} \left\{ Z(\mathbf{X}_{y,k}) \right\} \hspace{10mm} \forall k \\
  (b) &\hspace{5mm} aZ(\mathbf{X}) \geq \sum_{k = 1}^n \left( Z(\mathbf{X}) - \underset{y \in \mathcal{X}}{\inf} \left\{ Z(\mathbf{X}_{y,k}) \right\} \right)^2 
  \end{aligned}
  \nonumber
  $$

  For any $t > 0$:

  $$
  \mathbb{P}\left( \mathbb{E}\left[ Z(\mathbf{X}) \right] -  Z(\mathbf{X}) > t \right) \leq \exp \left( -\frac{t^2}{2a \mathbb{E}[Z(\mathbf{X})]} \right)
  \nonumber
  $$

  And if $Z(\mathbf{X})$ only satisfies the condition $(b)$ above, then:

  $$
  \mathbb{P}\left( Z(\mathbf{X}) - \mathbb{E}[Z(\mathbf{X})] \right) \leq \exp \left( -\frac{t^2}{2a\mathbb{E}[Z(\mathbf{X})]+ at}  \right)
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
    <p style="color:red;">TODO: FINISH PROOF</p>
    Proof is in <i>Concentration Inequalities for Functions of Independent Variables</i> by Andreas Maurer (2006). 
  </details>
  </body>
</div>

<div id="lemma-8"></div>
<div class="theorem">
  <body>
  <strong>Lemma 8 (Maurer and Pontil<span markdown="1">[^fn-maurer-b]</span>).</strong>
  <br>
  Suppose $X$ and $Y$ are i.i.d. random variables that are bounded in $[a, a + 1]$. We have that:

  $$
  \mathbb{E}_X\left[ \left( \mathbb{E}_Y \left[ (X - Y)^2 \right] \right)^2 \right] \leq \frac{1}{2}\mathbb{E}\left[ (X - Y)^2 \right]
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  First, the righthand side is equal to the variance.

  $$
  \begin{aligned}
  \frac{1}{2}\mathbb{E}\left[ (X - Y)^2 \right] &= \frac{1}{2} \mathbb{E}\left[ X^2 + Y^2 - 2XY \right] \\
  &= \frac{1}{2}\mathbb{E}\left[ X^2 \right] + \frac{1}{2}\mathbb{E}\left[ Y^2 \right] - \mathbb{E}\left[ X\right] \mathbb{E}\left[ Y \right] \hspace{10mm} (X, Y \text{ are independent})\\
  &= \mathbb{E}\left[ X^2 \right] - \left(\mathbb{E}\left[ X \right]\right)^2 \hspace{10mm} (X, Y \text{ are identically distributed}) \\
  &= \mathbb{E}\left[ X^2 - XY \right] \hspace{10mm} (X, Y \text{ are identically distributed})
  \end{aligned}
  \nonumber
  $$

  Expanding out the lefthand side yields:
  $$
  \begin{aligned}
  \mathbb{E}_X\left[ \left( \mathbb{E}_Y \left[ (X - Y)^2 \right] \right)^2 \right] &= \mathbb{E}_X\left[ \left( \mathbb{E}_Y \left[ X^2 + \mathbb{E}[Y^2] - 2X\mathbb{E}[Y] \right] \right)^2 \right] \\
  &= \mathbb{E} \left[ \left( X^2 + Y^2 - 2XY \right)^2 \right] \\
  &= \mathbb{E}_X \left[ X^4 + X^2 \mathbb{E}[Y^2] - 2X^3\mathbb{E}[Y] + X^2 \mathbb{E}[Y^2] + \left( \mathbb{E}[Y^2] \right)^2 - 2X \mathbb{E}[Y] \mathbb{E}[Y^2] - 2X^3\mathbb{E}[Y] - 2X \mathbb{E}[Y] \mathbb{E}[Y^2] + 4 X^2\left(\mathbb{E}[Y]\right)^2\right] \\
  &= \mathbb{E}_X\left[ X^4 + 2X^2\mathbb{E}[Y^2] - 4X^3 \mathbb{E}[Y] + (\mathbb{E}[Y^2])^2 -4X\mathbb{E}[Y]\mathbb{E}[Y^2] + 4X^2 (\mathbb{E}[Y])^2\right] \\
  &= \mathbb{E}[X^4] + 2 \mathbb{E}[X^2]\mathbb{E}[Y^2] - 4 \mathbb{E}[X^3]\mathbb{E}[Y] + (\mathbb{E}[Y^2])^2 - 4 \mathbb{E}[X] \mathbb{E}[Y] \mathbb{E}[Y^2] + 4 \mathbb{E}[X^2] (\mathbb{E}[Y])^2 \\
  &= \mathbb{E}[X^4] + 3\mathbb{E}[X^2]\mathbb{E}[Y^2] - 4\mathbb{E}[X^3]\mathbb{E}[Y] \hspace{10mm} (X,Y \text{ are identically distributed})
  \end{aligned}
  \nonumber
  $$

  Using the above facts, it suffices to show that the difference between the righthand side and lefthand side is greater than or equal to zero. That is, to show:

  $$
  \mathbb{E}[g(X,Y)] \geq 0 \hspace{5mm} \text{ where } \hspace{5mm} g(X,Y) := X^2 - XY - X^4 - 3X^2Y^2 + 4X^3Y
  \nonumber
  $$

  Consider the following:
  
  $$
  \begin{aligned}
  (X - Y + 1)(Y - X + 1)(Y-X)^2 &= (Y-X)^2(XY - X^2 + X - Y^2 + XY - Y + Y - X + 1) \\
  &= (Y - X)^2(2XY - X^2 - Y^2 + 1) \\
  &= (Y^2 + X^2 - 2XY)(2XY - X^2 - Y^2 + 1) \\
  &= 2Y^3X - X^2 Y^2 - Y^4 + Y^2 + 2X^3Y - X^4 - X^2 Y^2 + X^2 - 4X^2Y^2 + 2X^3Y + 2Y^3X - 2XY \\
  &= 4Y^3X + 4X^3Y -3X^2Y^2 - 3X^2Y^2 - XY - XY + Y^2 + X^2 - Y^4 - X^4 \\
  &=  (4X^3Y - X^4+ 4X^3Y - 3X^2Y^2- XY + X^2) + (4Y^3X -3X^2Y^2 - XY + Y^2 - Y^4) \\
  &= g(X,Y) + g(Y,X)
  \end{aligned}
  \nonumber
  $$

  Since $\rvert X - Y \rvert \leq 1$ due to their boundedness, $(X - Y + 1) \geq 0$ and $(Y - X + 1) \geq 0$ as well. Thus, $g(X,Y) + g(Y,X) \geq 0$.

  $$
  2\mathbb{E}[g(X,Y)] = \mathbb{E}[g(X,Y) + g(Y,X)] \geq 0 \implies \mathbb{E}[g(X,Y)] \geq 0
  \nonumber
  $$
  </details>
  </body>
</div>

<a href="#lemma-8">Lemma 8</a> yields the following corollary if $X$ and $Y$ are uniformly distributed over some set of values $\{ x_1, \dots, x_n \}$ with finite cardinality. If we assume $\{ x_1, \dots, x_n \} \subset [0, 1]$, then:

$$
\frac{1}{n} \sum_{k = 1}^n \left( \frac{1}{n}\sum_{j = 1}^n (x_k - x_j)^2 \right)^2 \leq \frac{1}{2n^2} \sum_{k = 1}^n \sum_{j = 1}^j (x_k - x_j)^2
\nonumber
$$

We use the above corollary in the proof of the following theorem.


<div id="theorem-10"></div>
<div class="theorem">
  <body>
  <strong>Theorem 10 (Maurer and Pontil<span markdown="1">[^fn-maurer-b]</span>).</strong>
  <br>
  Let $n \geq 2$, and let $\mathbf{X} = (X_1, \dots, X_n)$ be a vector of independent random variables such that $X_i \in [0, 1]$ for all $i$. Let $\mathbb{E}[V_n] := \mathbb{E}_\mathbf{X}[V_n(\mathbf{X})]$ where $V_n := V_n(\mathbf{x}) = \frac{1}{n(n-1)} \sum_{i = 1}^n \sum_{j = 1}^n \frac{(x_i - x_j)^2}{2}$ is the sample variance. For any $\delta > 0$, we have:

  $$
  \begin{aligned}
  (a) &\hspace{10mm} \mathbb{P}\left( \sqrt{\mathbb{E}[V_n]} > \sqrt{V_n(\mathbf{X})} + \sqrt{\frac{2 \log(1/\delta)}{n - 1}}\right) \leq \delta \\
  (b) &\hspace{10mm} \mathbb{P} \left( \sqrt{V_n(\mathbf{X})} > \sqrt{\mathbb{E}[V_n]} + \sqrt{\frac{2\log(1/\delta)}{n - 1}} \right) \leq \delta
  \end{aligned}
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  Let $Z(\mathbf{X}) := n V_n(\mathbf{X})$. Choose some $k$ and any $y \in [0, 1]$. It follows that:

  $$
  \begin{aligned}
  Z(\mathbf{X}) - Z(\mathbf{X}_{y,k}) &= \frac{n}{n(n-1)}\sum_{i = 1}^n \sum_{j = 1}^n \frac{(X_i - X_j)^2}{2} - \frac{1}{n(n-1)} \left(\frac{(y - y)^2}{2} + \sum_{j \neq k} \frac{(y - X_j)^2}{2} + \sum_{i \neq k} \left( \frac{(X_i - y)^2}{2} +  \sum_{j \neq k} \frac{(X_i - X_j)^2}{2}\right) \right) \\
  &= \frac{2}{n-1}\sum_{j = 1}^n \left[ \frac{(X_k - X_j)^2}{2} - \frac{(y - X_j)^2}{2} \right] \\
  &= \frac{1}{n-1} \sum_{j = 1}^n \left[ (X_k - X_j)^2 - \underbrace{(y - X_j)^2}_{\geq 0} \right] \\
  &\leq \frac{1}{n-1} \sum_{j = 1}^n (X_k - X_j)^2  \\
  \end{aligned}
  \nonumber
  $$
  
  We should note that the summation is over $n$ indices, but we could easily reduce it to only $n-1$ (restricting $j \neq k$). Thus, the last line above is the sample average of the $(X_k - X_j)^2$'s.

  Since each $X_i \in [0, 1]$ and $y \in [0, 1]$ it is clear that $Z(\mathbf{X}) - \underset{y \in [0, 1]}{\inf} Z(\mathbf{X}_{y,k}) \leq 1$. Then:

  $$
  \begin{aligned}
  \sum_{k = 1}^n \left( Z(\mathbf{X}) - \underset{y \in [0, 1]}{\inf} Z(\mathbf{X}_{y,k}) \right)^2 &\leq \sum_{k = 1}^n \left( \frac{1}{n-1} \sum_{j = 1}^n (X_k - X_j)^2 \right)^2 \hspace{10mm} \text{(above result)} \\
  &= \sum_{k = 1}^n \left( \frac{1}{n-1}\sum_{j = 1}^n (X_k - X_j)^2 \right)^2 \\
  &\leq \frac{n}{2(n-1)^2} \sum_{k = 1}^n \sum_{j = 1}^n (X_k - X_j)^2 \hspace{10mm} \text{(corollary of Lemma 8)}\\
  &= \frac{n}{n-1}\left[ \frac{1}{n-1} \sum_{k = 1}^n \sum_{j = 1}^n \frac{(X_k - X_j)^2}{2} \right] \\
  &= \frac{n}{n-1}Z(\mathbf{X})
  \end{aligned}
  \nonumber
  $$

  The above steps show that $Z(\mathbf{X})$ satisfy conditions $(a)$ and $(b)$ in <a href="#theorem-7">Theorem 7</a> where $a = \frac{n}{n-1}$. Now, notice that:
  
  $$
  \begin{aligned}
  &\mathbb{P}\left(\mathbb{E}\left[ V_n(\mathbf{X}) \right] - V_n(\mathbf{X}) > t\right) = \mathbb{P}\left(\mathbb{E}\left[ \frac{1}{n}V_n(\mathbf{X})\right] - \frac{1}{n}Z(\mathbf{X})> t \right) = \mathbb{P}\left( \mathbb{E}\left[ Z(\mathbf{X}) \right] - Z(\mathbf{X}) > nt \right) \\
  &\mathbb{P}\left(-\mathbb{E}\left[ V_n(\mathbf{X}) \right] + V_n(\mathbf{X}) > t\right) = \mathbb{P}\left(-\mathbb{E}\left[ \frac{1}{n}V_n(\mathbf{X})\right] + \frac{1}{n}Z(\mathbf{X})> t \right) = \mathbb{P}\left( -\mathbb{E}\left[ Z(\mathbf{X}) \right] + Z(\mathbf{X}) > nt \right) \\
  \end{aligned}
  \nonumber
  $$

  Applying Theorem 7, we have, for any $t > 0$:

  $$
  \begin{aligned}
  &\mathbb{P}\left(\mathbb{E}[Z(\mathbf{X})] - Z(\mathbf{X}) > nt \right) \leq \exp\left( -\frac{t^2n^2(n-1)}{2n\mathbb{E}[Z(\mathbf{X})]}\right) = \exp \left(- \frac{t^2n(n-1)}{2\mathbb{E}\left[n V_n(\mathbf{X}) \right]} \right) = \exp \left(-\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right]} \right)\\
  &\mathbb{P}\left(Z(\mathbf{X}) - \mathbb{E}[Z(\mathbf{X})] > nt\right) \leq \exp \left( - \frac{t^2n^2(n-1)}{2n\mathbb{E}[Z(\mathbf{X})] + n^2t} \right) = \exp \left(-\frac{t^2n(n-1)}{2\mathbb{E}\left[nV_n(\mathbf{X})\right] + nt} \right) = \exp \left( - \frac{t^2(n-1)}{2\mathbb{E}[V_n(\mathbf{X})] + t} \right)
  \end{aligned}
  \nonumber
  $$

  This implies that:

  $$
  \begin{aligned}
  &\mathbb{P}\left(\mathbb{E}\left[ V_n(\mathbf{X}) \right] - V_n(\mathbf{X}) > t\right) \leq \exp \left(-\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right]} \right) \\
  &\mathbb{P}\left(-\mathbb{E}\left[ V_n(\mathbf{X}) \right] + V_n(\mathbf{X}) > t\right) \leq \exp \left( - \frac{t^2(n-1)}{2\mathbb{E}[V_n(\mathbf{X})] + t} \right)
  \end{aligned}
  \label{eq:intermed-probs}
  $$

  The first line in Eq. \eqref{eq:intermed-probs} implies that with probability at least $1 - \delta$:

  $$
  \mathbb{E}\left[ V_n(\mathbf{X}) \right] - V_n(\mathbf{X}) \leq \sqrt{\frac{2 \mathbb{E}\left[V_n(\mathbf{X}) \right] \log \left( \frac{1}{\delta} \right)}{n-1}} \implies \mathbb{E}\left[ V_n(\mathbf{X}) \right] - \sqrt{\frac{4 \mathbb{E}\left[V_n(\mathbf{X}) \right] \log \left( \frac{1}{\delta} \right)}{2(n-1)}}  \leq V_n(\mathbf{X}) \implies  \mathbb{E}\left[ V_n(\mathbf{X}) \right] - 2\sqrt{\mathbb{E}\left[V_n(\mathbf{X})\right]}\sqrt{ \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}}  \leq V_n(\mathbf{X})
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  We simply solve the following for $t$:
  $$
  \mathbb{P}\left(\mathbb{E}\left[ V_n(\mathbf{X}) \right] - V_n(\mathbf{X}) > t\right) \leq \delta = \exp \left(-\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right]} \right)
  \nonumber
  $$
  Thus:
  $$
  \begin{aligned}
  \delta &= \exp \left(-\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right]} \right) \\
  &\implies \log (\delta) = -\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right]}  \\
  &\implies \log \left( \frac{1}{\delta} \right) = \frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right]}  \\
  &\implies 2 \mathbb{E}\left[V_n(\mathbf{X}) \right] \log \left( \frac{1}{\delta} \right) =t^2(n-1) \\
  &\implies \sqrt{\frac{2 \mathbb{E}\left[V_n(\mathbf{X}) \right] \log \left( \frac{1}{\delta} \right)}{n-1}} = t
  \end{aligned}
  \nonumber
  $$
  </details>
  
  With some manipulations, this implies that, with probability at least $1 - \delta$:

  $$
  \sqrt{\mathbb{E}\left[ V_n(\mathbf{X}) \right]} \leq \sqrt{V_n(\mathbf{X})} + \sqrt{\frac{2\log \left( \frac{1}{\delta} \right)}{n-1}}
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  $$
  \begin{aligned}
  \mathbb{E}\left[ V_n(\mathbf{X}) \right] - 2\sqrt{\mathbb{E}\left[V_n(\mathbf{X})\right]}\sqrt{ \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}}  \leq V_n(\mathbf{X}) &\implies \mathbb{E}\left[ V_n(\mathbf{X}) \right] - 2\sqrt{\mathbb{E}\left[V_n(\mathbf{X})\right]}\sqrt{ \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}} + \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)} \leq V_n(\mathbf{X}) + \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)} \\
  &\implies  \left(\sqrt{\mathbb{E}\left[ V_n(\mathbf{X}) \right]} - \sqrt{\frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}}\right) \leq V_n(\mathbf{X}) + \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)} \\
  &\implies \sqrt{\mathbb{E}\left[ V_n(\mathbf{X}) \right]} - \sqrt{\frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}} \leq \sqrt{V_n(\mathbf{X}) + \frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}} \overset{(i)}{\leq} \sqrt{V_n(\mathbf{X})} + \sqrt{\frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}} \\
  &\implies \sqrt{\mathbb{E}\left[ V_n(\mathbf{X}) \right]} \leq \sqrt{V_n(\mathbf{X})} + 2\sqrt{\frac{\log \left( \frac{1}{\delta} \right)}{2(n-1)}} = \sqrt{V_n(\mathbf{X})} + \sqrt{\frac{4\log \left( \frac{1}{\delta} \right)}{2(n-1)}} = \sqrt{V_n(\mathbf{X})} + \sqrt{\frac{2\log \left( \frac{1}{\delta} \right)}{n-1}}
  \end{aligned}
  \nonumber
  $$

  where we use the fact that $\sqrt{a + b} \leq \sqrt{a} + \sqrt{b}$ in $(i)$.
  </details>

  which is Claim $(a)$ in the Theorem statement. 
  <br>
  We can do similar manipulations of the second line in Eq. \eqref{eq:intermed-probs} to get, with probability $1 - \delta$:

  $$
  \begin{aligned}
  \sqrt{V_n(\mathbf{X})} &\leq  \sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{2\log\left(\frac{1}{\delta}\right)}{n-1}}
  \end{aligned}
  \nonumber
  $$

  <details>
  <summary>Proof.</summary>
  We simply solve the following for $t$:
  $$
  \mathbb{P}\left(V_n(\mathbf{X}) - \mathbb{E}\left[ V_n(\mathbf{X}) \right] > t\right) \leq \delta = \exp \left(-\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) + t\right]} \right)
  \nonumber
  $$

  Thus:
  $$
  \begin{aligned}
  \delta &= \exp \left(-\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) + t\right]} \right) \\
  &\implies \log (\delta) = -\frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right] + t}  \\
  &\implies \log \left( \frac{1}{\delta} \right) = \frac{t^2(n-1)}{2\mathbb{E}\left[V_n(\mathbf{X}) \right] + t}  \\
  &\implies \left(2 \mathbb{E}\left[V_n(\mathbf{X}) \right] + t \right) \log \left( \frac{1}{\delta} \right) =t^2 \\
  &\implies \frac{\left(2 \mathbb{E}\left[V_n(\mathbf{X}) \right] + t \right) \log \left( \frac{1}{\delta} \right)}{n-1} = t^2 \\
  &\implies \frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta} \right)}{n-1} = t^2 - \frac{2\log\left(\frac{1}{\delta}\right)t}{2(n-1)} \\
  &\implies \frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta} \right)}{n-1}  + \left(\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2 = t^2 - \frac{2\log\left(\frac{1}{\delta}\right)t}{2(n-1)} + \left(\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2 \\
  &\implies \frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta} \right)}{n-1}  + \left(\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2 = \left(t - \frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2 \\
  &\implies \sqrt{\frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta} \right)}{n-1}  + \left(\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2} = t - \frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)} \\
  &\implies \sqrt{\frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta} \right)}{n-1} + \left(\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2} + \frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)} = t 
  \end{aligned}
  \nonumber
  $$

  Thus, with probability at least $1 - \delta$:

  $$
  \begin{aligned}
  V_n(\mathbf{X}) - \mathbb{E}\left[ V_n(\mathbf{X}) \right] &\leq \sqrt{\frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta} \right)}{n-1} + \left(\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}\right)^2} + \frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)} \\
  &\overset{(i)}{\leq} \sqrt{\frac{2 \mathbb{E}\left[ V_n(\mathbf{X}) \right]\log\left(\frac{1}{\delta}\right)}{n-1}} + \frac{\log\left(\frac{1}{\delta} \right)}{2(n - 1)} + \frac{\log\left(\frac{1}{\delta} \right)}{2(n - 1)} \\
  \implies V_n(\mathbf{X}) &\leq \mathbb{E}\left[ V_n(\mathbf{X}) \right] + 2 \sqrt{\mathbb{E}\left[ V_n(\mathbf{X}) \right]}\sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}} + \frac{2\log\left(\frac{1}{\delta} \right)}{2(n - 1)} \\
  \implies V_n(\mathbf{X}) + \frac{\log\left(\frac{1}{\delta} \right)}{2(n-1)} &\leq \mathbb{E}\left[ V_n(\mathbf{X}) \right] + 2 \sqrt{\mathbb{E}\left[ V_n(\mathbf{X}) \right]}\sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}} +\frac{\log\left(\frac{1}{\delta} \right)}{2(n-1)} + \frac{\log\left(\frac{1}{\delta} \right)}{n - 1}  \\
  &= \left(\sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}}\right)^2 + \frac{\log\left(\frac{1}{\delta}\right)}{n-1} \\
  \implies V_n(\mathbf{X}) &\leq \left(\sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}}\right)^2 + \frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)} \\
  \implies \sqrt{V_n(\mathbf{X})} \leq \sqrt{\left(\sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}}\right)^2 + \frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}} \\
  &\overset{(ii)}{\leq} \sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}} + \sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}} \\
  &= \sqrt{\mathbb{E}[V_n(\mathbf{X})]} + 2\sqrt{\frac{\log\left(\frac{1}{\delta}\right)}{2(n-1)}} \\
  &=  \sqrt{\mathbb{E}[V_n(\mathbf{X})]} + \sqrt{\frac{2\log\left(\frac{1}{\delta}\right)}{n-1}}
  \end{aligned}
  \nonumber
  $$

  In $(i)$ and $(ii)$, we use the fact that $\sqrt{a + b} \leq \sqrt{a} + \sqrt{b}$.
  </details>

  And this is Claim $(b)$ in the Theorem.
  </details>
  </body>
</div>

Finally, we just need an extension of Bennett's inequality that is integral to the proof of <a href="#theorem-11">Theorem 11</a>.

<div id="bennett-derivative"></div>
<div class="theorem">
  <body>
  <strong>Theorem 3 (A Derivative Of Bennett's Inequality).</strong>
  <br>
  Let $\mathbf{Z} = (Z_1, \dots, Z_n)$ be a vector of i.i.d. random variables with $Z_i \in [0,1]$. With probability at least $1 - \delta$, for any $\delta > 0$, we have that:

  $$
  \begin{aligned}
  \mathbb{E}[Z] - \frac{1}{n} \sum_{i = 1}^n Z_i \leq \sqrt{\frac{2 V(Z) \log\left(\frac{1}{\delta}\right)}{n}} + \frac{\log\left(\frac{1}{\delta}\right)}{3n}
  \end{aligned}
  \nonumber
  $$

  where $Z$ is any one of the i.i.d. coordinates of $\mathbf{Z}$, and $V(Z) := \mathbb{E}[(Z - \mathbb{E}[Z])^2]$. See Audibert et al.<span markdown="1">[^fn-audibert]</span>for a proof.
  </body>
</div>



<!-- 
Since the $Z_i$'s are in $[0, 1]$ and are independent, we have by the [Chernoff bound](https://en.wikipedia.org/wiki/Chernoff_bound#Sums_of_independent_random_variables) for independent random variables:

$$
\mathbb{P}\left( \sum_{i = 1}^n Z_i \geq t\right) \leq \underset{s > 0}{\inf} \left\{ \exp(-s) \prod_{i = 1}^n  \mathbb{E}\left[ \exp(s Z_i) \right]\right\}
\label{eq:chernoff}
$$

Recall that we assume the $Z_i$'s are in $[0,1]$ almost surely. That is, $\mathbb{P}(\rvert Z_i \rvert \leq 1) = 1$ for all $i$. This implies that $\mathbb{P}\left( \rvert Z_i - \mathbb{E}[Z_i] \rvert \leq 2\right) = 1$. Almost surely bounded random variables have the [following property](https://www.planetmath.org/relationbetweenalmostsurelyabsolutelyboundedrandomvariablesandtheirabsolutemoments). For any integer $k \geq 1$:

$$
\mathbb{E}\left[ \rvert Z_i - \mathbb{E}[Z_i] \rvert^k \right] \leq 2^k
\nonumber
$$

Let $a = 2$. Then:

$$
\begin{aligned}
\mathbb{E}\left[ \exp\left(s(Z_i - \mathbb{E}[Z_i]) \right)\right] &= \mathbb{E}\left[ \sum_{k = 0}^\infty \frac{(s(Z_i -\mathbb{E}[Z_i]))^k}{k!}\right] \\
&= \sum_{k = 0}^\infty \frac{s^k}{k!}\mathbb{E}\left[ (Z_i - \mathbb{E}[Z_i])^k \right] \\
&= 1 + s\mathbb{E}\left[ Z_i - \mathbb{E}[Z_i] \right] + \sum_{k = 2}^\infty \frac{s^k}{k!}\mathbb{E}\left[ (Z_i - \mathbb{E}[Z_i])^{k-2} (Z_i - \mathbb{E}[Z_i])^{2} \right] \\
&\leq 1 + \sum_{k = 2}^\infty \frac{s^k a^(k-2)}{k!}\mathbb{E}\left[(Z_i - \mathbb{E}[Z_i])^{2} \right] \\
&= 1 + \sum_{k = 2}^\infty \frac{s^k a^k V(Z_i)}{a^2 k!}\\ 
&= 1 + \sum_{k = 0}^\infty \frac{s^k a^k V(Z_i)}{a^2 k!} - \frac{V(Z_i)}{a^2} - \frac{V(Z_i)}{a^2}as\\
&= 1 +\frac{V(Z_i)}{a^2} \left(\exp(sa) - as - 1 \right) \\
&\overset{(i)}{\leq} \exp \left( \frac{V(Z_i)}{a^2} \left(\exp(sa) - as - 1 \right) \right)
\end{aligned}
\nonumber
$$

where $(i)$ follows from the fact that $1 + x \leq e^x$. 
<br>
Centering and putting this together with Eq. \eqref{eq:chernoff}, we get:

$$
\begin{aligned}
\mathbb{P}\left( \sum_{i = 1}^n (Z_i - \mathbb{E}[Z_i]) \geq t \right) &\leq \underset{s > 0}{\inf} \left\{ \exp(-s) \prod_{i = 1}^n  \mathbb{E}\left[ \exp(s (Z_i - \mathbb{E}[Z_i])) \right]\right\} \\
&= \underset{s > 0}{\inf} \left\{ \exp(-s) \mathbb{E}\left[ \exp\left(s \sum_{i = 1}^n (Z_i - \mathbb{E}[Z_i])\right) \right]\right\} \\
&= \underset{s > 0}{\inf} \left\{ \exp(-s)\prod_{i = 1}^n \mathbb{E}\left[ \exp\left(s (Z_i - \mathbb{E}[Z_i])\right) \right]\right\} \\
&\leq \underset{s > 0}{\inf} \left\{ \exp(-s)\prod_{i = 1}^n \exp \left( \frac{V(Z_i)}{a^2} \left(\exp(sa) - as - 1 \right) \right) \right\} \\
&= \underset{s > 0}{\inf} \left\{ \exp \left(-s + \sum_{i = 1}^n \frac{V(Z_i)}{a^2} \left(\exp(sa) - as - 1 \right) \right) \right\} \\
\end{aligned}
\nonumber
$$

We can minimize the above to find the best choice of $s$. Since argument minima/maxima are invariant to strictly monotonic transformations, we can set the derivative of the part inside the exponential function equal to 0 and solve:

$$
\begin{aligned}
&0 = \frac{d}{ds}\left[ -s + \sum_{i = 1}^n \frac{V(Z_i)}{a^2} \left(\exp(sa) - as - 1 \right) \right] \\
\implies &0 = -1 + \sum_{i = 1}^n \frac{V(Z_i)}{a^2}\left(s\exp(as) -s  \right) \\
\implies &1 = (a \exp(as) - a)\sum_{i = 1}^n \frac{V(Z_i)}{a^2} \\
\implies &\frac{a^2}{\sum_{i = 1}^n V(Z_i)} = a\exp(as) - a \\
\implies &\frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 = \exp(as)\\
\implies &\frac{1}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) = s
\end{aligned}
\nonumber
$$

Plugging in yields:

$$
\begin{aligned}
\mathbb{P}\left( \sum_{i = 1}^n (Z_i - \mathbb{E}[Z_i]) \geq t \right) &\leq \exp \left(-\frac{1}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) + \sum_{i = 1}^n \frac{V(Z_i)}{a^2} \left(\exp\left(\frac{a}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right)\right) - \frac{a}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) - 1 \right) \right)  \\
&=  \exp \left(-\frac{1}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) + \sum_{i = 1}^n \frac{V(Z_i)}{a^2} \left(\frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 - \log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) - 1 \right) \right)  \\
&= \exp \left(-\frac{1}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) + \sum_{i = 1}^n \frac{V(Z_i)}{a^2} \left(\frac{a}{\sum_{i = 1}^n V(Z_i)} - \log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) \right) \right)  \\
&= \exp \left(-\frac{1}{a}\log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) + \frac{1}{a} - \sum_{i = 1}^n \frac{V(Z_i)}{a^2} \log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right) \right)   \\
&= \exp \left(\frac{1}{a} - \log\left( \frac{a}{\sum_{i = 1}^n V(Z_i)} + 1 \right)\left(\frac{1}{a}+ \sum_{i = 1}^n \frac{V(Z_i)}{a^2}\right) \right)   \\
\end{aligned}
\nonumber
$$ -->



---
## References

Additional content from Maurer[^fn-maurer-a] and Anthony and Bartlett[^fn-anthony].

[^fn-audibert]: Audibert, J.-Y., Munos, R., & Szepesvári, C. (2009). Exploration–exploitation tradeoff using variance estimates in multi-armed bandits. Theoretical Computer Science, 410(19), 1876–1902. https://doi.org/10.1016/j.tcs.2009.01.016

[^fn-anthony]: Anthony, M., & Bartlett, P. L. (1999). Neural Network Learning: Theoretical Foundations (1st ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511624216

[^fn-maurer-a]: Maurer, A. (2006). Concentration inequalities for functions of independent variables. Random Structures & Algorithms, 29(2), 121–138. https://doi.org/10.1002/rsa.20105

[^fn-maurer-b]: Maurer, A., & Pontil, M. (2009). Empirical Bernstein Bounds and Sample Variance Penalization (arXiv:0907.3740). arXiv. https://doi.org/10.48550/arXiv.0907.3740

[^fn-vershynin]: Vershynin, R. (2018). High-Dimensional Probability: An Introduction with Applications in Data Science. Cambridge: Cambridge University Press.

[^fn-wasserman]: Wasserman, Larry. All of Statistics: a concise course in statistical inference. New York: Springer, 2010.

[^fn-bennett]: Wikimedia Foundation. (2024, May 2). Bennett’s Inequality. Wikipedia. https://en.wikipedia.org/wiki/Bennett%27s_inequality 
