---
layout: distill
title:  "Concentration Inequalities"
description: A Primer
date: 2025-02-20
tabs: true
tags: theory primer
bibliography: stats-ml.bib

---

In this post, I'm going to take myself through a review of standard concentration inequalities in probability theory with an ultimate goal of exploring empirical Bernstein bounds (<a href="/paper-notes/empir-bern-bound">in this post</a>).

Note: Not all of the proofs are finished/included. I am hoping to find the time to return to this post and complete them.

---

In this section, we'll cover four main concentration inequalities that will aid in understanding later sections of this post. These are covered in Vershynin <d-cite key=vershynin2018></d-cite> and Wasserman<d-cite key=wasserman2004></d-cite> Two very standard concentration inequalities are Markov's inequality and Chebyshev's inequality.

<!-- #region markov-ineq -->
<div class="theorem">
  <strong>Markov's Inequality.</strong>
  {% tabs markov-ineq %}
  {% tab markov-ineq statement %}
  Let $X$ be a non-negative random variable with finite mean (i.e. $\mathbb{E}[X]$ exists). For any choice of $t > 0$, we have:

  $$
  \mathbb{P}(X \geq t) \leq \frac{\mathbb{E}[X]}{t}
  $$

  {% endtab %}
  {% tab markov-ineq proof %}
  Let $f(x)$ denote the probability density function of $X$. Since $X \geq 0$, we have:

  $$
  \begin{aligned}
  \mathbb{E}[X] = \int_0^\infty x f(x) dx = \underbrace{\int_0^t x f(x) dx}_{\geq 0} + \int_t^\infty xf(x) dx \geq \int_t^\infty x f(x) \overset{(i)}{\geq} t \int_t^\infty f(x) dx = t \mathbb{P}(X \geq t)
  \end{aligned}
  $$ 

  In $(i)$, we use the fact that on the interval $[t, \infty]$, $x \geq t$. 
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endregion -->

<!-- #region cheby-ineq -->
<div class="theorem">
  <strong>Chebyshev's Inequality.</strong>
  {% tabs cheby-ineq %}
  {% tab cheby-ineq statement %}
  Let $X$ be a random variable with finite mean $\mathbb{E}\left[X\right] = \mu$ and finite variance $V(X) = \sigma^2$. Define $Z = \frac{X - \mu}{\sigma}$. We have:

  $$
  \mathbb{P}(\rvert X - \mu \rvert \geq t) \leq \frac{\sigma^2}{t^2} \hspace{10mm} \text{ and } \hspace{10mm} \mathbb{P}(\rvert Z \rvert \geq k) \leq \frac{1}{k^2}
  $$
  {% endtab %}
  {% tab cheby-ineq proof %}
  Since $\rvert X - \mu \rvert \geq 0$, we can use Markov's Inequality:

  $$
  \mathbb{P}(\rvert X - \mu \rvert \geq t) = \mathbb{P}(\rvert X - \mu \rvert^2 \geq t^2) \leq \frac{\mathbb{E}\left[ (X - \mu)^2 \right]}{t^2} \overset{(i)}{=} \frac{\sigma^2}{t^2}
  \nonumber
  $$

  In $(i)$, we choose $t = k \sigma$.
  {% endtab %}
  {% endtabs %}
</div>
<!--#endregion -->

A very important tail bound that will be helpful for the remainder of this post is Hoeffding's inequality, which bounds the deviation of the sum of independent, bounded random variables from their mean. 

<!-- #region hoeff-ineq -->
<div class="theorem">
  <strong>Hoeffding's Inequality.</strong>
  {% tabs hoeff-ineq %}
  {% tab hoeff-ineq statement %}
  Let $X_1, \dots, X_n$ be independent random variables such that $a_i \leq X_i \leq b_i$ almost surely. Let $\mathbb{E}[X_i] =: \mu_i$. For any $t > 0$:

  $$
  \mathbb{P}\left(\sum_{i = 1}^n (X_i - \mu_i) \geq t \right) \leq \exp \left( - \frac{2 t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \right)
  $$
  {% endtab %}
  {% tab hoeff-ineq proof %}
  First, we show the reformulation. Choose some $t > 0$ and let $\delta = \exp \left( - \frac{2t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \right)$. Then:
  
  $$
  \begin{aligned}
  -\log(\delta) &= \frac{2t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \\
  \implies t^2 &= \frac{\sum_{i = 1}^n (b_i - a_i)^2 \log\left(\frac{1}{\delta} \right)}{2} \\
  \implies t &= \sqrt{\frac{1}{2}\sum_{i = 1}^n (b_i - a_i)^2 \log\left(\frac{1}{\delta}\right)}
  \end{aligned}
  $$

  Putting the above together with the converse of the first probability expression yields the reformulation.
    <p style="color:red;">TODO: FINISH PROOF</p>
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endregion -->

An important note made in Vershynin<d-cite key=vershynin2018></d-cite> is the fact that this bound is non-asymptotic; that is, it holds for any choice of $n$. The bound just becomes tighter as the sample size grows.

Unfortunately, Hoeffding's inequality is independent of the variances of the random variables in question. We can imagine that a better bound could be obtained when the variances are small since we know they should remain close to the means. Bernstein inequalities aim to do this by trading sharpness in cases with large variance for improved bounds in cases with small variance. 

There are several similar Berstein inequalities, one of which I've written below:

<!-- #region bern-ineq -->
<div class="theorem">
  <strong>A Bernstein Inequality.</strong>
  {% tabs bern-ineq %}
  {% tab bern-ineq statement %}
  Let $X_1, \dots, X_n$ be i.i.d. bounded random variables (i.e. their support is $[a, b]$ for finite $a, b$). Let $\mathbb{E}[X_i] =: \mu$ and $\mathbb{E}\left[ (X - \mu)^2 \right] =: \sigma^2$. For any $t \geq $, We have:

  $$
  \mathbb{P}\left( \bigg\rvert \frac{1}{n} \sum_{i = 1}^n X_i - \mu \bigg\rvert \geq t \right) \leq 2 \exp\left( -\frac{n t^2}{2(\sigma^2 + (b-a) t)} \right)
  $$
  {% endtab %}
  {% tab bern-ineq proof %}
  <p style="color:red;">TODO: FINISH PROOF</p>
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endregion -->

Bennett's inequality is similar to the Bernstein inequalities in that it bounds the sum of independent random variables in a variance-dependent fashion. 

<!-- #region benn-ineq -->
<div id="theorem-bennett"></div>
<div class="theorem">
  <strong>Bennett's Inequality.<d-cite key=bennett2025></d-cite></strong>
  {% tabs benn-ineq %}
  {% tab benn-ineq statement %}
  Let $X_1, \dots, X_n$ be independent random variables with finite mean $\mathbb{E}[X_i] =: \mu$. Let $\sigma^2 := \sum_{i = 1}^n \mathbb{E}\left[ (X_i - \mathbb{E}[X_i])^2\right]$. For all $i$, suppose $\rvert X_i - \mathbb{E}[X_i] \rvert \leq a$ almost surely. For any $t \geq $, We have:

  $$
  \mathbb{P}\left( \sum_{i = 1}^n (X_i - \mathbb{E}[X_i]) \geq t \right) \leq \exp\left( -\frac{\sigma^2}{a^2} \left[ \left(1 + \frac{at}{\sigma^2}\right) \log\left(1 + \frac{at}{\sigma^2}\right) - \frac{at}{\sigma^2} \right] \right)
  $$
  {% endtab %}
  {% tab benn-ineq proof %}
  <p style="color:red;">TODO: FINISH PROOF</p>
  {% endtab %}
  {% endtabs %}
</div>
<!-- #endreigon -->

We also have the following variation on Bennett's inequality from Audibert et al. (2009).<d-cite key=audibert2009></d-cite>

<div class="theorem-bennett2"></div>
<div class="theorem">
  <strong>Lemma 5 (Audibert et al.<d-cite key=audibert2009></d-cite>).</strong>
  {% tabs audi-ineq %}
  {% tab audi-ineq statement %}
  Let $X \in \mathbb{R}$ be a random variable such that $X \leq s$ for some $s \in \mathbb{R}$ almost surely. Denote $$\mu := \mathbb{E}[X]$$ and $$s_+ := s \vee 0$$. Define $X_1, \dots, X_n$ as independent and identically distributed copies of $X$. Let $$\bar{U}_m := \frac{1}{m} \sum_{i = 1}^m U_i$$. Then, for any $\delta > 0$ and simultaneously for all $1 \leq m \leq n$, with probability at least $ 1 - \exp(-\delta)$:

  $$
  \begin{aligned}
  &m(\bar{U}_m - \mu) \leq \sqrt{2 n \mathbb{E}[U^2]\delta} + \frac{s_+ \delta}{3}\\
  &m(\bar{U}_m - \mu) \leq \sqrt{2 n \text{Var}(U)\delta} + \frac{(s - \mu)\delta}{3}
  \end{aligned}
  $$
  {% endtab %}
  {% tab audi-ineq proof %}
  A proof can be found in the publication.
  {% endtab %}
  {% endtabs %}
</div>
