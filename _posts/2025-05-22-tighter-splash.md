---
layout: post
title:  "Improved Inference For SPLASH"
date: 22 May 2025
categories: posts
tags: ["sequencing", "genomics", "probability"]
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
    "/assets/js/modal.js",
]
---


This post builds off my <a href="/posts/2025/02/02/concentration.html">earlier post</a> on concentration inequalities and empirical Bernstein bounds. Here, I'm going to try to apply those ideas to get a better bound on the $p$-values from SPLASH[^fn-chaung] and OASIS[^fn-baharav].

The $p$-value bound in SPLASH/OASIS are based on an application of Hoeffding's inequality. However, the authors mention that it may be possible to do better through a Bernstein inequality. The issue is that we do not know the variances of the summands in the test statistic. This is where empirical Bernstein bounds come into play. Here I will try to apply the ideas explained above in Maurer and Pontil (2009)[^fn-maurer] to the SPLASH setting.

---

## SPLASH

Recall that the SPLASH setting is one in which we have a count matrix $X$ that is $I \times J$ where $I$ is the number of targets and $J$ is the number of samples. We use $n_j$ and $r_i$ to denote the sums of column $j$ and row $i$, respectively, and $M := \sum_{j = 1}^J n_J = \sum_{i = 1}^I r_I$ is the total count of the entire table. We additionally have fixed input vectors $\mathbf{f} \in \mathbb{R}^I$ and $\mathbf{c} \in \mathbb{R}^J$, which we may later restrict to some subset of the reals. The test statistic is $S(\mathbf{f}, \mathbf{c}) := \mathbf{f}^\top \tilde{X} \mathbf{c}$ where $\tilde{X} := (X - E)\text{diag}((X^\top \mathbf{1})^{-1/2}),$ and $E := \frac{1}{M} X\mathbf{1}\mathbf{1}^\top X$. For more details, see my <a href="/posts/2025/02/27/splash.html">post on SPLASH</a>.

As before, we need to rewrite the test statistic as the sum of independent bounded random variables. Define $Z_{j,k}$ as the random variable that denotes the target index for the $k$-th target observed in the $j$-th sample. Also define $$\bar{\mu} := \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}}$$, $$\hat{\mu}_j := \frac{1}{n_j} \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}}$$, and $$\mu := \mathbb{E}_{Z \sim \mathbf{p}}[\mathbf{f}_Z]$$. It can be shown that we can rewrite the test statistic as:

$$
S(\mathbf{f}, \mathbf{c}) = \sum_{j=1}^J \left(  \mathbf{c}_j\sqrt{n_j} (\hat{\mu}_j - \bar{\mu}) \right) = \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}  \sum_{k = 1}^{n_j} \left(\frac{1}{n_j} \mathbf{f}_{Z_{j,k}} - \frac{1}{n_j} \bar{\mu} \right) = \sum_{j = 1}^J \sum_{k =1 }^{n_j}\frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})
$$

Define the following:

$$
Y_{j,k} := \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})
$$

It follows that $S(\mathbf{f}, \mathbf{c}) = \sum_{j = 1}^J \sum_{k = 1}^{n_j} Y_{j,k}$. Furthermore, $\mathbb{E}[Y_{j,k}] = 0$ for all $j,k$.

<details>
<summary>Proof Of $\mathbb{E}[Y_{j,k}] = 0$.</summary>
$$
\begin{aligned}
\mathbb{E}[Y_{j,k}] &= \mathbb{E}\left[ \frac{\mathbf{c}_j}{\sqrt{n_j}} (\mathbf{f}_{Z_{j,k}} - \bar{\mu})\right] \\
&= \frac{\mathbf{c}_j}{\sqrt{n_j}} \left( \mathbb{E}[\mathbf{f}_{Z_{j,k}}] -  \frac{1}{M} \sum_{l = 1}^J \sum_{h = 1}^{n_l} \mathbb{E}[\mathbf{f}_{Z_{l,h}}] \right) \\
&= \frac{\mathbf{c}_j}{\sqrt{n_j}} \left( \mu -  \frac{1}{M} \sum_{l = 1}^J \sum_{h = 1}^{n_l} \mu \right) \\ 
&= \frac{\mathbf{c}_j}{\sqrt{n_j}} \left( \mu -  \frac{M}{M} \mu \right) \\ 
&= 0
\end{aligned} 
\nonumber
$$
</details>


Let's also assume that $\mathbf{f} \in [0, 1]^I$. This implies that:

$$
Y_{j,k} \in \left[-\bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert \bar{\mu}, \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert (1 - \bar{\mu}) \right]
$$


<details>
<summary>Proof.</summary>
$$
\begin{aligned}
\min Y_{j,k} &= \min \left\{ \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})\right\} \\
&\geq\min\left\{ \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert  (\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right\} \\
&=  \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert \min\left\{  \mathbf{f}_{Z_{j,k}} - \bar{\mu} \right\} \\
&= -\bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert \bar{\mu} \\
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
\max Y_{j,k} &= \max \left\{ \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})\right\} \\
&\leq\max\left\{ \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert  (\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right\} \\
&= \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert (\max\{ \mathbf{f}_{Z_{j,k}} \} - \bar{\mu}) \\
&= \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} \bigg\rvert (1 - \bar{\mu}) \\
\end{aligned}
\nonumber
$$
</details>

Since $\bar{\mu} \in [0, 1]$, $1-\bar{\mu} \in [0, 1]$ as well. Thus, for all $(j,k)$:

$$
\rvert Y_{j,k} \rvert \leq \max_h \bigg\rvert \frac{\mathbf{c}_h}{\sqrt{n_h}} \bigg\rvert := \alpha
\label{eq:y-bounds}
$$

### OASIS Bound

Our interest is in bounding the probability that the absolute value of $S$ exceeds some value. In OASIS, this is done with an application of Hoeffding's inequality. If $\rvert \rvert \mathbf{c} \rvert \rvert_2^2 \leq 1$ and $\gamma < 1$, then, for any $t > 0$:

$$
\mathbb{P}\left(\rvert S(\mathbf{f}, \mathbf{c}) \rvert \geq t \right) \leq  2 \exp\left(-\frac{2 t^2}{1 - \gamma}\right)
\label{eq:oasis-bound}
$$

---

## Empirical Bernstein Bounds

The SPLASH and OASIS $p$-value bounds use Hoeffding's inequality, which only relies upon the <i>range</i> of the bounded random variables we have. In our case, we would not expect to be in a "high variance regime" because it is, at least at face value, unlikely that an anchor will see all possible target sequences in all samples. What would be ideal would be a bound that incorporates the <i>variance</i> of the random variables, since we think that most take on values in a smaller range. Unfortunately, bounds of this type (e.g. Bernstein and Bennett) require knowledge of the <i>true variance</i> of the data, which is <span class="popup" onclick="PopupFunc('pop1')">much too complicated to even attempt to derive<span class="popuptext" id="pop1">Recall that we are dealing with categorical random variables that are indexing vectors and doing summations and scaling 🥲...</span></span>. Fortunately, methods have been developed to use the <i>empirical variance</i> from the sample at hand, making it possible to compute tighter bounds (in some cases) than ones that arise from Hoeffding's inequality. 

Maurer and Pontil[^fn-maurer] developed one of the first <i>empirical Bernstein bounds</i>. For independent random variables $X_1, \dots, X_n \in [a, b]$, for any $t > 0$:

$$
\mathbb{P}\left( \bigg\rvert \sum_{i = 1}^n \left[X_i - \mathbb{E}[X_i] \right]\bigg\rvert \geq t\right) \leq 4\exp\left(- 2\left(\sqrt{\frac{3(n-1)t}{7 (b-a) n} + \left(\frac{3(n-1)\sqrt{\frac{V_n(\mathbf{X})}{2n}}}{7 (b-a) }\right)^2}  - \frac{3(n-1)\sqrt{V_n(\mathbf{X})}}{2n} \right)\right)
$$

See my <a href="/posts/2025/02/02/concentration.html">concentration inequalities post</a> for details.


### Application To SPLASH
Applying this theorem to our $Y_{j,k}$, we see that, for any $t > 0$:

$$
\mathbb{P}(\rvert S(\mathbf{f}, \mathbf{c}) \rvert \geq t) = \mathbb{P}\left(\bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} Y_{j,k} \bigg\rvert \geq t \right) \leq 4 \exp\left( -2\left( \sqrt{\frac{3(M-1)t}{7(2\alpha)M} + \left( \frac{3(M-1)\sqrt{V_n(\mathbf{Y})}}{7(2\alpha)\sqrt{2M}} \right)^2} - \frac{3(M - 1)\sqrt{V_n(\mathbf{Y})}}{2M} \right)\right)
\label{eq:emp-bern-bound}
$$

where $V_n(\mathbf{Y}) = \frac{1}{2M(M-1)} \sum_{j, k}\sum_{j', k'}(Y_{j,k} - Y_{j',k'})^2$.

<details>
<summary>$V_n(\mathbf{Y})$ Another Way.</summary>
$$
\begin{aligned}
V_n(\mathbf{Y}) &= \frac{1}{2M(M-1)} \sum_{j, k}\sum_{j', k'}(Y_{j,k} - Y_{j',k'})^2 \\
&= \frac{1}{2M(M-1)} \sum_{j, k}\sum_{j', k'}\left(\frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu}) - \frac{\mathbf{c}_{j'}}{\sqrt{n_{j'}}}(\mathbf{f}_{Z_{j',k'}} - \bar{\mu})\right)^2 \\
&= \frac{1}{2M(M-1)} \sum_{j,k} \sum_{j', k'} \left(\frac{\mathbf{c}_j^2}{n_j}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})^2 + \frac{\mathbf{c}^2_{j'}}{n_{j'}}(\mathbf{f}_{Z_{j',k'}} - \bar{\mu})^2 - \frac{2\mathbf{c}_j\mathbf{c}_{j'}}{\sqrt{n_j}\sqrt{n_{j'}}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})(\mathbf{f}_{Z_{j',k'}} - \bar{\mu}) \right) \\
&= \frac{1}{2M(M-1)}\left[ 2\sum_{j,k}\sum_{j',k'}\frac{\mathbf{c}_j^2}{n_j}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})^2  - 2 \sum_{j,k} \sum_{j',k'} \frac{\mathbf{c}_j\mathbf{c}_{j'}}{\sqrt{n_j}\sqrt{n_{j'}}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})(\mathbf{f}_{Z_{j',k'}} - \bar{\mu}) \right] \\
&= \frac{2}{2M(M-1)}\left[ M\sum_{j,k} \frac{\mathbf{c}_j^2}{n_j}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})^2 - \sum_{j,k}\frac{\mathbf{c}_j}{\sqrt{n_j}} (\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \sum_{j',k'} \frac{\mathbf{c}_{j'}}{\sqrt{n_{j'}}} (\mathbf{f}_{Z_{j',k'}} - \bar{\mu}) \right] \\
&= \frac{1}{M-1} \left[ \sum_{j,k}\left(  \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right)^2 - \frac{1}{M}\left( \sum_{j,k} \frac{\mathbf{c}_j}{\sqrt{n_{j'}}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right)^2\right]
\end{aligned}
\nonumber
$$
</details>


A drawback of this method is that, in practice, it often requires very favorable samples. Simulations[^fn-hult] have shown that only in the case of very small variance and very large sample sizes do empirical Bernstein bounds beat Hoeffding. In my own explorations, I've found that, even in some seemingly favorable cases, the bound in Eq. \ref{eq:emp-bern-bound} decays faster than Hoeffding for a bit but eventually loses out as $t$ grows large. 

<div class="row">
    <div class="column">
        <img src="/img/2025-05-22-tighter-splash/bounds_binom_n500_p0.99.pdf" alt="Empirical Bernstein Simulation (p = 0.99)"/>
        <p>Figure A</p>
    </div>
    <div class="column">
        <img src="/img/2025-05-22-tighter-splash/bounds_binom_n500_p0.9.pdf" alt="Empirical Bernstein Simulation (p = 0.90)"/>
        <p>Figure B</p>
    </div>
    <div class="column">
        <img src="/img/2025-05-22-tighter-splash/bounds_binom_n500_p0.5.pdf" alt="Empirical Bernstein Simulation (p = 0.50)"/>
        <p>Figure C</p>
    </div>
</div>

The above figures show simulations of $500$ binomial random variables with different success probabilities. The red lines show the $p$-value bounds that result from using Hoeffding's inequality. In blue and green, we show the empirical Bernstein bounds and an estimate of the true $p$-value made from resampling. In yellow, we show the true probabilities.

As we go from Figure A to C, the variance increases, and we clearly see good performance in the low variance case (when $p = 0.99$), though it is difficult to tell exactly what is happening in the tail. In contrast, even when the success probability decreases just to $0.9$, the Hoeffding-derived $p$-value bound beats the empirical Bernstein-derived bound for moderate values of $t$. This effect is at an extreme when we have the highest variance ($p = 0.5$).

---

## Betting Bounds
In all three figures, we see that the Hoeffding-derived $p$-value bound is still quite far from the true $p$-value. This begs the question: can we do any better? Luckily, results from game theoretic statistics say yes (in a way).

A recent paper from Waudby-Smith and Ramdas[^fn-waudby-smith] introduces methods for tight confidence intervals and sequences for estimating the means of bounded random variables. I have <a href="/posts/2025/05/21/betting-bounds.html">another post</a> about the nitty gritty details of this paper, but I'll briefly introduce the pertinent result here. Before we dive in, it might be useful to review some of the concepts in my <a href="/posts/2025/04/30/measure-theory.html">post on measure theory</a> and my <a href="/posts/2025/05/20/martingales.html">post on martingales</a>.

Suppose we have a sequence of random variables $$(X_t)_{t=1}^\infty$$ from some distribution $P \in \mathcal{P}^\mu$ where $\mathcal{P}^\mu$ is the set of distributions on $[0, 1]^\infty$ such that $$\mathbb{E}_P[X_t \rvert \mathcal{F}_{t-1}] = \mu$$ for each $t$ where $$\mathcal{F} = (\mathcal{F}_t)_{t=0}^\infty$$ is the canonical filtration of $$(X_t)_{t=1}^\infty$$.

Let $m, \theta \in [0, 1]$, $c = \frac{1}{2}$ or $\frac{3}{4}$, and define:

$$
\mathcal{K}_t^\pm = \max\left\{ \theta \underbrace{\prod_{i = 1}^t (1 + \lambda_i^+(m) \cdot(X_i - m))}_{\mathcal{K}_t^+(m)}, \hspace{2mm} (1-\theta) \underbrace{\prod_{i = 1}^t(1 + \lambda_i^-(m)\cdot (X_i - m))}_{\mathcal{K}_t^-(m)} \right\}
\label{eq:hedged-capital-process}
$$

$$
\begin{aligned}
\lambda_i^+(m) &= \min\left\{ \bigg\rvert \sqrt{\frac{2\log(2/\delta)}{n \hat{\sigma}_{i-1}^2}}\bigg\rvert, \frac{c}{m}\right\} 
&\hspace{5mm}
\lambda_i^-(m) &= \min\left\{ \bigg\rvert \sqrt{\frac{2\log(2/\delta)}{n \hat{\sigma}_{i-1}^2}}\bigg\rvert, \frac{c}{1-m}\right\} \\
\hat{\mu}_i &= \frac{\frac{1}{2} + \sum_{k = 1}^i X_k}{i + 1} 
&\hspace{5mm}
\hat{\sigma}_i^2 &= \frac{\frac{1}{4} + \sum_{k = 1}^i(X_k - \hat{\mu}_k)^2}{i + 1}
\end{aligned}
$$

The following set forms a valid $(1 - \delta)$ confidence interval for $\mu$ for any fixed $n$:

$$
C_n = \bigcap_{i = 1}^n \mathcal{B}_i^\pm
\hspace{5mm} \text{ where }
\mathcal{B}_i^\pm = \left\{ m \in [0, 1] \bigg\rvert \mathcal{K}_i^\pm(m) < \frac{1}{\delta} \right\}
\label{eq:hedged-capital-ci}
$$

By a $(1 - \delta)$ confidence interval, we mean that $\mathbb{P}_P\left( \mu \not\in C_n \right) \leq \delta$. 

### Application To SPLASH
Using Eq. \ref{eq:y-bounds}, we can transform the $Y_{j,k}$ variables to be $[0, 1]$-valued so that we can use the above results by substituting $\tilde{Y}_{j,k}$ in place of $X_i$. We'll denote these shifted and scaled values:

$$
\tilde{Y}_{j,k} = \frac{Y_{j,k} - \alpha}{2\alpha}
$$

But what do we do with a confidence interval? Thus far, we have been concerned with bounding the probability $\mathbb{P}\left( \rvert S(\mathbf{f}, \mathbf{c}) \rvert \geq t \right)$ for some $t > 0$. If the confidence interval were less complicated, we could may invert it to obtain an expression for $\alpha$ and use that as a bound. However, the above confidence interval is very involved; there is many maximums, minimums, and intersections taken amongst many different sets that are, themselves, not entirely simple to compute. 

The authors of SPLASH state that <i>"...SPLASH calculates a p value bound for the null hypothesis that the observed target frequencies in samples all come from the same distribution, i.e., that there is no underlying variation of targets between samples"</i>.[^fn-chaung] Let $$\mathbf{p} \in [0, 1]^I$$ such that $$\sum_{i = 1}^I \mathbf{p}_i = 1$$. In statistical terms, this null hypothesis is that:

$$
H_0: \hspace{2mm} \forall j \in [J], \hspace{2mm} X_{\cdot, j} \rvert n_j \overset{iid}{\sim} \text{Multinom}\left(\mathbf{p}\right) 
\label{eq:oasis-null}
$$

This can be alternatively stated using the $Z_{j,k}$ random variables:

$$
H_0: \hspace{2mm} \forall j \in [J], k \in [n_j], \hspace{2mm} Z_{j,k} \overset{iid}{\sim} \text{Categorical}\left( \mathbf{p} \right)
\label{eq:oasis-null-z}
$$

For a fixed $\mathbf{f}$, $$\mathbb{E}_{Z \sim \mathbf{p}}\left[ \mathbf{f}_{Z_{j, k}} \right]$$ should be the same for all $j \in [J]$ and $k \in [n_j]$. That implies that, conditional on the column sums and fixing $\mathbf{f}$ and $\mathbf{c}$, all $Y_{j,k}$'s have mean zero. We can therefore restate the hypothesis as:

$$
H_0: \hspace{2mm} \forall j \in [J], k \in [n_j], \hspace{2mm} \mathbb{E}\left[ Y_{j,k} \right] = 0
\label{eq:oasis-null-y}
$$

Under the null, we have that:

$$
\begin{aligned}
\mathbb{P}\left( \rvert S(\mathbf{f}, \mathbf{c}) \rvert \geq t \right) 
= \mathbb{P}\left( \rvert S(\mathbf{f}, \mathbf{c}) - \mathbb{E}\left[ S(\mathbf{f}, \mathbf{c}) \right] \rvert \geq t \right) 
= \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} (Y_{j,k} - \mathbb{E}[Y_{j,k}]) \bigg\rvert \geq t\right) 
= \mathbb{P}\left( \bigg\rvert \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} (Y_{j,k} - \mathbb{E}[Y_{j,k}]) \bigg\rvert \geq \frac{t}{M}\right) 
= \mathbb{P}\left( \rvert \bar{Y} - \mathbb{E}[Y] \rvert \geq \frac{t}{M} \right)
\end{aligned}
\nonumber
$$

For any $s = \frac{t}{M} \geq 0$, the Hoeffding probability bound (see Eq. \ref{eq:oasis-bound}) yields:

$$
\mathbb{P}\left( \rvert \bar{Y} - \mathbb{E}[Y] \rvert \geq s \right) \leq 2 \exp\left( - \frac{2 s^2}{1 - \gamma} \right)
$$

Setting the righthandside of the above equal to $\delta$ and solving for $s$ gives:

$$
s = \sqrt{\frac{(1 - \gamma) \log(2/\delta)}{2}} 
\hspace{5mm} \implies \hspace{5mm}
\mathbb{P}\left(\rvert \bar{Y} - \mathbb{E}[Y] \rvert \geq  \sqrt{\frac{(1- \gamma)\log(2/\delta)}{2}} \right) \leq \delta
$$

Inverting the above probability statement gives us the following $(1 - \delta)$ confidence interval for $\mathbb{E}[Y]$:

$$
C_M = \bar{Y} \pm \sqrt{\frac{(1-\gamma) \log(2/\delta)}{2}} 
$$

If $0 \not \in C_M$, then we would reject $H_0$!

Since the problem thus reduces to one of testing the mean of bounded random variables, we can use the results of Waudby-Smith and Ramdas[^fn-waudby-smith]! We can test the null hypothesis using the confidence interval in Eq. \ref{eq:hedged-capital-ci} instead of $C_M$ described above. Since the OASIS bound improves on a strictly Hoeffding-based bound by using information we assume about the data, we don't necessarily have the same convergence and optimality results for the betting intervals that Waudby-Smith and Ramdas prove. 

However, from some elementary simulations, it looks like the betting-based confidence intervals are always less wide than the OASIS ones! Instead of looking at anchors that have $p$-values less than some designated $\delta$, we could look at the anchors whose associated confidence intervals exclude $0$. And instead of ranking the statistically significant anchors by their $p$-values, we could rank them by confidence interval width. 

Unfortunately, in practice, this method does not improve the power that significantly. In addition, it appears to loose the robustness property against biologically uninteresting alternatives that is a large strength of OASIS.


---
## References
[^fn-baharav]: Baharav, T. Z., Tse, D., & Salzman, J. (2024). OASIS: An interpretable, finite-sample valid alternative to Pearson’s X2 for scientific discovery. Proceedings of the National Academy of Sciences, 121(15), e2304671121. doi:10.1073/pnas.2304671121

[^fn-chaung]: Chaung, K., Baharav, T. Z., Henderson, G., Zheludev, I. N., Wang, P. L., & Salzman, J. (2023). SPLASH: A statistical, reference-free genomic algorithm unifies biological discovery. Cell, 186(25), 5440-5456.e26. https://doi.org/10.1016/j.cell.2023.10.028

[^fn-hult]: Hult, L. (2022, March 18). Empirical Bernstein bounds. Ludvig Hult. https://el-hult.github.io/2022/03/18/empirical-bernstein-bounds.html 

[^fn-maurer]: Maurer, A., & Pontil, M. (2009). Empirical Bernstein Bounds and Sample Variance Penalization (arXiv:0907.3740). arXiv. https://doi.org/10.48550/arXiv.0907.3740

[^fn-ravikumar]: Ravikumar, P. (2019). Carnegie Mellon 10-716 Lecture 6 Notes. https://www.cs.cmu.edu/~pradeepr/courses/716/2019-spring/notes/lec6.pdf

[^fn-vershynin]: Vershynin, R. (2018). High-Dimensional Probability: An Introduction with Applications in Data Science. Cambridge: Cambridge University Press.

[^fn-waudby-smith]: Waudby-Smith, I., & Ramdas, A. (2023). Estimating means of bounded random variables by betting. Journal of the Royal Statistical Society Series B: Statistical Methodology, 86(1), 1–27. doi:10.1093/jrsssb/qkad009.








<!-- 
### Optimization

Minimizing the bound is the same as maximizing the test statistic since the stuff inside the exponential function is a function of $t$. Since that function has a negative coefficient, we want to make the part in the innermost parentheses the biggest we can.

Assume $\mathbf{c}$ is known. Our goal is to find:

$$
\underset{\mathbf{f} \in [0, 1]^J}{\arg \min} \left\{ 4 \exp\left( -2 \left[ \sqrt{\frac{3(M-1)t}{7(2\alpha)M} + \left(\frac{3(M-1) \sqrt{V_n(\mathbf{Y})}}{7(2\alpha) \sqrt{2M}}\right)^2} - \frac{3(M-1)\sqrt{V_n(\mathbf{Y})}}{2M} \right] \right) \right\}
\nonumber
$$

where we take $t = S(\mathbf{f}, \mathbf{c})$. For convenience, we restate the variables below:

$$
\begin{aligned}
M &= \sum_{j = 1}^J n_j \\
\alpha &= \underset{h}{\max} \bigg\rvert \frac{\mathbf{c}_h}{\sqrt{n_j}} \bigg\rvert \\
Y_{j,k} &= \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \\
S(\mathbf{f}, \mathbf{c}) &= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu})
\end{aligned}
\nonumber
$$

Define the set $$\Lambda_i = \left\{ (j,k) \rvert Z_{j,k} = i \right\}$$ as the set of index pairs $(j,k)$ such that $Z_{j,k} = i$ for $i \in [I]$. We then have:

$$
\begin{aligned}
\frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} \right] &= 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{\mathbf{c}_j}{\sqrt{n_j}} (\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right] = 
\begin{cases}
\frac{\mathbf{c}_j}{\sqrt{n_j}} & \text{ if } (j,k) \in \Lambda_i \\
0 & \text{ otherwise }
\end{cases}
\\ \\
\frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right] &= 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{\mathbf{c}^2_j}{n_j} (\mathbf{f}_{Z_{j,k}} - \bar{\mu})^2 \right] = 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{\mathbf{c}^2_j}{n_j} (\mathbf{f}^2_{Z_{j,k}} - 2\mathbf{f}_{Z_{j,k}}\bar{\mu} + \bar{\mu}^2) \right] = 
\begin{cases}
\frac{2\mathbf{c}^2_j}{n_j} \left(\mathbf{f}_{i} - \bar{\mu}\right) & \text{ if } (j,k) \in \Lambda_i \\
0 & \text{ otherwise }
\end{cases}
\\ \\
\frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] &= 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{\mathbf{c}^2_j}{n_j} (\mathbf{f}_{Z_{j,k}} - \bar{\mu})(\mathbf{f}_{Z_{j, k'}} - \bar{\mu}) \right] = 
\begin{cases}
\frac{2\mathbf{c}^2_j}{n_j} \left(\mathbf{f}_{i}- \bar{\mu}\right) & \text{ if } (j,k), (j, k') \in \Lambda_i \\
\frac{\mathbf{c}^2_j}{n_j} \left(\mathbf{f}_{Z_{j,k'}} - \bar{\mu}\right) & \text{ if } (j,k) \in \Lambda_i, (j, k') \not\in \Lambda_i \\
\frac{\mathbf{c}^2_j}{n_j} \left(\mathbf{f}_{Z_{j,k}} - \bar{\mu}\right) & \text{ if } (j,k) \not\in \Lambda_i, (j, k') \in \Lambda_i \\
0 & \text{ otherwise }
\end{cases}
\\ \\
\frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right] &= 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} (\mathbf{f}_{Z_{j,k}} - \bar{\mu})(\mathbf{f}_{Z_{j', k}} - \bar{\mu}) \right] = 
\begin{cases}
\frac{2 \mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{i} - \bar{\mu} \right) & \text{ if } (j,k), (j', k) \in \Lambda_i \\
\frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j,k'}} - \bar{\mu}\right) & \text{ if } (j,k) \in \Lambda_i, (j', k) \not\in \Lambda_i \\
\frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j,k}} - \bar{\mu}\right) & \text{ if } (j,k) \not\in \Lambda_i, (j', k) \in \Lambda_i \\
0 & \text{ otherwise }
\end{cases}
\\ \\
\frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] &= 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} (\mathbf{f}_{Z_{j,k}} - \bar{\mu})(\mathbf{f}_{Z_{j', k'}} - \bar{\mu}) \right] = 
\begin{cases}
\frac{2 \mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{i} - \bar{\mu}\right) & \text{ if } (j,k), (j', k') \in \Lambda_i \\
\frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j',k'}} - \bar{\mu}\right) & \text{ if } (j,k) \in \Lambda_i, (j', k') \not\in \Lambda_i \\
\frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j,k}} - \bar{\mu}\right) & \text{ if } (j,k) \not\in \Lambda_i, (j', k') \in \Lambda_i \\
0 & \text{ otherwise }
\end{cases}
\end{aligned} 
\nonumber
$$

We then see that:

$$
\begin{aligned}
\frac{\partial}{\partial \mathbf{f}_i} \left[ S(\mathbf{f}, \mathbf{c}) \right] &= 
\frac{\partial}{\partial \mathbf{f}_i} \left[ \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\mathbf{c}_j}{\sqrt{n_j}} (\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right] 
= \sum_{(j,k) \in \Lambda_i} \frac{\mathbf{c}_j}{\sqrt{n_j}}
\end{aligned}
\nonumber
$$



And also:

$$
\begin{aligned}
\frac{\partial}{\partial \mathbf{f}_i} \left[ V_n(\mathbf{Y}) \right] &= \frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{1}{2M(M-1)} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{j' = 1}^J \sum_{k' = 1}^{n_{k'}} (Y_{j,k} - Y_{j', k'})^2 \right] \\
&= \frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{1}{2M(M-1)} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{j' = 1}^J \sum_{k' = 1}^{n_{k'}} \left( Y^2_{j,k} + Y^2_{j', k'} - 2 Y_{j, k} Y_{j', k'}  \right) \right] \\
&= \frac{\partial}{\partial \mathbf{f}_i} \left[ \frac{1}{2M(M-1)} \left( M\sum_{j = 1}^J \sum_{k = 1}^{n_j} Y^2_{j,k} + M \sum_{j' = 1}^J \sum_{k' = 1}^{n_{j'}} Y^2_{j', k'} -2 \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{j' = 1}^J \sum_{k' = 1}^{n_{j'}} Y_{j, k} Y_{j', k'}  \right) \right] \\
&= \frac{1}{2(M-1)} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right] + \frac{1}{2(M-1)} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right] - \frac{1}{M(M-1)} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{j' = 1}^J \sum_{k' = 1}^{n_{j'}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j, k} Y_{j', k'}\right] \\
&= \frac{1}{M-1} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right] - \frac{1}{M(M-1)} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{j' = 1}^J \sum_{k' = 1}^{n_{j'}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j, k} Y_{j', k'}\right] \\
&= \frac{2}{M-1} \sum_{(j,k) \in \Lambda_i} \frac{\mathbf{c}^2_j}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right) + \frac{1}{M(M-1)} \underbrace{\sum_{j = 1}^J  \sum_{k = 1}^{n_j} \sum_{k' = 1}^{n_{j}}  \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right]}_{(i)} + \frac{1}{M(M-1)} \underbrace{\sum_{j = 1}^J \sum_{j' \neq j}  \sum_{k = 1}^{n_j} \sum_{k' = 1}^{n_{j'}}  \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right]}_{(ii)} \\
&= \frac{2}{M-1} \sum_{(j,k) \in \Lambda_i} \frac{\mathbf{c}^2_j}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right)  + \frac{1}{M(M-1)} \left[ 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{k' : (j,k') \in \Lambda_i}}{\sum}\frac{\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right)\right)
+ 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{k' : (j,k') \not\in \Lambda_i \\ k' \neq k}} \frac{\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_{Z_{j,k'}} - \bar{\mu} \right) \right) \right. \\
&\hspace{10mm} + \left. 2 \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_i - \bar{\mu} \right)\right) 
+ 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left( \mathbf{f}_{Z_{j',k'}} - \bar{\mu} \right) \right) \right] \\
&= \frac{2}{M-1} \sum_{(j,k) \in \Lambda_i} \frac{\mathbf{c}^2_j}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right)  + \frac{1}{M(M-1)} \left[ 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{k' : (j,k') \not\in \Lambda_i \\ k' \neq k}} \frac{\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_{Z_{j,k'}} - \bar{\mu} \right) \right) \right. \\
&\hspace{10mm} + \left. 2 \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_i - \bar{\mu} \right)\right) 
+ 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left( \mathbf{f}_{Z_{j',k'}} - \bar{\mu} \right) \right) \right] \\
\end{aligned}
\nonumber
$$




<details>
<summary>Details Of $(i), (ii)$.</summary>
$$
\begin{aligned}
(i) &= \sum_{j = 1}^J  \sum_{k = 1}^{n_j} \sum_{k' = 1}^{n_{j'}}  \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right] + \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \\
&=  \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right]
    + \underbrace{\underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right]}_{=0}
    + \sum_{j = 1}^J \sum_{k = 1}^{n_j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \\
&=  \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right]
    + \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right]
    + \underset{(j,k) \not \in \Lambda_i}{\sum \sum}\sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \\
&=  \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right]
    + \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left[ \underset{\substack{k' : (j,k') \in \Lambda_i \\ k' \neq k}}{\sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] + \underset{\substack{k': (j,k') \not\in \Lambda_i \\ k' \neq k}}{\sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \right] \\
&\hspace{10mm} + \underset{j,k: (j,k) \not \in \Lambda_i}{\sum \sum} \left[ \underset{\substack{k': (j,k') \in \Lambda_i \\ k' \neq k}}{\sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] + \underbrace{\underset{\substack{k': (j,k') \not\in \Lambda_i \\ k' \neq k}}{\sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right]}_{=0} \right] \\
&= \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum}\left( \frac{\partial}{\partial \mathbf{f}_i} \left[ Y^2_{j,k} \right] 
    + \underset{\substack{k' : (j,k') \in \Lambda_i \\ k' \neq k}}{\sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \right)
    + 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \sum_{\substack{k' : (j,k') \not\in \Lambda_i \\ k' \neq k}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right]  \\
&= \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \frac{2\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right)
    + \underset{\substack{k' : (j,k') \in \Lambda_i \\ k' \neq k}}{\sum} \frac{2\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right) \right) 
    + 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \sum_{\substack{k' : (j,k') \not\in \Lambda_i \\ k' \neq k}} \frac{\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_{Z_{j,k'}} - \bar{\mu} \right) \\
&=  2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{k' : (j,k') \in \Lambda_i}}{\sum}\frac{\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_i - \bar{\mu}\right)\right)
     + 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{k' : (j,k') \not\in \Lambda_i \\ k' \neq k}} \frac{\mathbf{c}_j^2}{n_j}\left(\mathbf{f}_{Z_{j,k'}} - \bar{\mu} \right) \right)
\end{aligned}
\nonumber
$$

$$
\begin{aligned}
    (ii) &=  \sum_{j = 1}^J  \sum_{j' \neq j} \sum_{k = 1}^{n_j} \sum_{k' = 1}^{n_{j'}}  \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j, k'} \right] \\
&= \sum_{j = 1}^J \sum_{j' \neq j} \sum_{k = 1}^{n_j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right] + \sum_{j = 1}^J \sum_{j' \neq j} \sum_{k = 1}^{n_j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] \\
&=  \underbrace{\sum_{(j,k) \in \Lambda_i} \sum_{j' \neq j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right]}_{(a)} + 
    \underbrace{\sum_{(j,k) \not\in \Lambda_i} \sum_{j' \neq j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right]}_{(b)} + 
    \underbrace{\sum_{(j,k) \in \Lambda_i} \sum_{j' \neq j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right]}_{(c)} +
    \underbrace{\sum_{(j,k) \not\in \Lambda_i} \sum_{j' \neq j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right]}_{(d)} \\
&= 2 \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{ \mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}}\left(\mathbf{f}_i - \bar{\mu}\right)  \right) 
    + \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \not\in \Lambda_i \\ j' \neq j}} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j', k}} - \bar{\mu} \right) \right)\\
&\hspace{10mm} + \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}}\left(\mathbf{f}_{Z_{j,k}} - \bar{\mu} \right) \right)\\
&\hspace{10mm} +2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_i - \bar{\mu} \right)\right) +
    \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left( \mathbf{f}_{Z_{j',k'}} - \bar{\mu} \right) \right) \\
&\hspace{10mm} + \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left(\underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j,k}} - \bar{\mu} \right)\right) \\
&= 2 \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{ \mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}}\left(\mathbf{f}_i - \bar{\mu}\right)  \right) 
    + 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \not\in \Lambda_i \\ j' \neq j}} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j', k}} - \bar{\mu} \right) \right)\\
&\hspace{10mm} +2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_i - \bar{\mu} \right)\right) + 2
    \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left( \mathbf{f}_{Z_{j',k'}} - \bar{\mu} \right) \right) \\
&= 2 \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_i - \bar{\mu} \right)\right) + 2
    \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left( \mathbf{f}_{Z_{j',k'}} - \bar{\mu} \right) \right)
\end{aligned}
\nonumber
$$

<details>
<summary>Details Of $(a), (b), (c), (d)$.</summary>
$$
\begin{aligned}
(a) &= \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum}\left(\sum_{j' \neq j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right]\right) \\
    &= \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left[ \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right] + \sum_{\substack{j': (j', k) \not\in \Lambda_i \\ j' \neq j}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right] \right]\\
    &= 2 \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{ \mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}}\left(\mathbf{f}_i - \bar{\mu}\right)  \right) 
    + \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \not\in \Lambda_i \\ j' \neq j}} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j', k}} - \bar{\mu} \right) \right) \\ \\
(b) &= \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left( \sum_{j' \neq j} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right] \right) \\
    &= \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum}\left( \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right] \right) + \underbrace{\underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \not\in \Lambda_i \\ j' \neq j}} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k} \right]\right)}_{=0} \\
    &= \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left( \sum_{\substack{j': (j', k) \in \Lambda_i \\ j' \neq j}} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}}\left(\mathbf{f}_{Z_{j,k}} - \bar{\mu} \right) \right) \\ \\
(c) &= \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \sum_{j' \neq j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right]\right) \\
    &=  \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] \right) +
    \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] \right) \\
    &= 2\underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_i - \bar{\mu} \right)\right) +
    \underset{j,k : (j,k) \in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left( \mathbf{f}_{Z_{j',k'}} - \bar{\mu} \right) \right) \\ \\
(d) &= \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left( \sum_{j' \neq j} \sum_{k' \neq k} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] \right) \\
    &= \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left( \underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k}}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] \right) +
    \underbrace{\underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum}\left(\underset{\substack{j', k': (j',k') \not\in \Lambda_i \\ j' \neq j \\ k' \neq k }}{\sum \sum} \frac{\partial}{\partial \mathbf{f}_i} \left[ Y_{j,k} Y_{j', k'} \right] \right) }_{=0} \\
    &= \underset{j,k : (j,k) \not\in \Lambda_i}{\sum \sum} \left(\underset{\substack{j', k': (j',k') \in \Lambda_i \\ j' \neq j \\ k' \neq k}}{\sum \sum} \frac{\mathbf{c}_j \mathbf{c}_{j'}}{\sqrt{n_j n_{j'}}} \left(\mathbf{f}_{Z_{j,k}} - \bar{\mu} \right) \right)
\end{aligned}
\nonumber
$$
</details>
</details>

If we let $$\mathbf{e}_{j,k} = (0, 0, \dots, 1, \dots, 0, 0)$$ be the standard basis vector (with zeroes at all coordinates except the $Z_{j,k}$-th), then we can write the partial derivative of $Y_{j,k}$ with respect to $\mathbf{f}$ as:

$$
\frac{\partial }{\partial \mathbf{f}} \left[ Y_{j,k} \right] = \frac{\partial }{\partial \mathbf{f}} \left[ \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right]  = \frac{\mathbf{c}_{j}}{\sqrt{n_j}} \mathbf{e}_{j,k}
\nonumber
$$

$$
\begin{aligned}
\frac{\partial}{\partial \mathbf{f}} \left[ t \right] &= \frac{\partial}{\partial \mathbf{f}} \left[ \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \bar{\mu}) \right]
\end{aligned}
\nonumber
$$

This result implies that:

$$
\frac{\partial}{\partial \mathbf{f}} \left[ \mathbf{Y} \right] = 
\begin{bmatrix}
    \rule[.5ex]{2.5ex}{0.5pt} & \frac{\partial}{\partial \mathbf{f}} \left[ Y_{1, 1} \right] & \rule[.5ex]{2.5ex}{0.5pt} \\
    & \vdots & \\
    \rule[.5ex]{2.5ex}{0.5pt} & \frac{\partial}{\partial \mathbf{f}} \left[ Y_{J, n_J} \right] & \rule[.5ex]{2.5ex}{0.5pt}
\end{bmatrix}
\nonumber
$$ -->
