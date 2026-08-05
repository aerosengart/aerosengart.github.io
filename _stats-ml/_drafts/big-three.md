---
layout: distill
title: The Big Three
description: The Likelihood Ratio, Wald, and Score Tests
date: 2026-06-30
tabs: true
tags: theory likelihood primer
toc:  
  - name: Background  
  - name: The Big Three
    subsections:
      - name: The Wald Test
      - name: The Likelihood Ratio Test
      - name: The Rao Score Test
bibliography: stats-ml.bib
---

This post mostly covers the end of Chapter 7 of <i>Elements of Large-Sample Theory</i> by Erich Lehmann.<d-cite key=lehmann2004></d-cite> For context and preliminary results and definitions, see <a href="/stats-mle/maximum-likelihood">my post</a> on maximum likelihood estimation.

---

## Background
Again assume that $X_1, \dots, X_n$ are i.i.d. with $X_i \sim P_{\boldsymbol{\theta}}$, and let $\hat{\boldsymbol{\theta}}$ be an estimator for $\boldsymbol{\theta}$ that satisfies:

$$
\sqrt{n} \left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}\right) \rightsquigarrow \mathcal{N} (0, \mathcal{I}^{-1}(\boldsymbol{\theta}))
$$

<aside><p>That is, it is efficient.</p></aside>

Let $\hat{\mathcal{I}}(\boldsymbol{\theta})$ be any <i>consistent</i> estimate of $\mathcal{I}(\boldsymbol{\theta})$. Then:

$$
\sqrt{n} \left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta} \right) \hat{\mathcal{I}}^{\frac{1}{2}}(\boldsymbol{\theta}) \rightsquigarrow \mathcal{N}(0, 1)
$$

Letting $z_{\frac{\alpha}{2}}$ be the $\frac{\alpha}{2}$ quantile of the standard normal distribution, we can rearrange the above to get:

$$
\hat{\boldsymbol{\theta}} - \frac{z_{\boldsymbol{\alpha}{2}}}{\sqrt{n}}\hat{\mathcal{I}(\boldsymbol{\theta})} < \boldsymbol{\theta} < \hat{\boldsymbol{\theta}} + \frac{z_{\frac{\alpha}{2}}}{\sqrt{n \hat{\mathcal{I}(\boldsymbol{\theta})}}}
$$

which is a level $\alpha$ asymptotic confidence interval for $\boldsymbol{\theta}$. 

In the event that $\mathcal{I}(\boldsymbol{\theta})$ is continuous in $\boldsymbol{\theta}$, then one choice of estimator is $\hat{\mathcal{I}}(\boldsymbol{\theta}) = \mathcal{I}(\hat{\boldsymbol{\theta}})$. Another choice is to use $\hat{\mathcal{I}}(\boldsymbol{\theta}) = - \frac{1}{n} \left. \frac{\partial^2 \ell(\boldsymbol{\theta}; \mathbf{x})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right \rvert_{\boldsymbol{\theta} = \hat{\boldsymbol{\theta}}}$ under the conditions of Theorem 7.5.2.

In what follows, we will consider testing:

$$
\begin{equation}
\label{eq:test-1}
H_0: \boldsymbol{\theta} = \boldsymbol{\theta}_0
\hspace{5mm} vs. \hspace{5mm}
H_1: \boldsymbol{\theta} \neq \boldsymbol{\theta}_0
\end{equation}
$$

Let $\nu_{\alpha}$ be the $1-\alpha$ quantile for the $\chi^2_k$ distribution. 

---

## The Big Three

### The Wald Test
The <strong>Wald test</strong> statistic is constructed as:

$$
W_n = n \left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0 \right)^\top \hat{\mathcal{I}}(\boldsymbol{\theta}) \left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0\right)
$$

By Theorem 7.5.2., $$W_n \rightsquigarrow \chi^2_k$$, so we reject $H_0$ if $ W_n \geq \nu_{\alpha}$.


### The Likelihood Ratio Test
The <strong>likelihood ratio test</strong> rejects $H_0$ when the likelihood function (or log-likelihood function) is must greater at $\hat{\boldsymbol{\theta}}$ (a consistent root of $\ell(\boldsymbol{\theta}; \mathbf{x})$), than at the hypothesized value. This test statistic is given by:

$$
\Delta_n = 2 \left[ \ell(\hat{\boldsymbol{\theta}}; \mathbf{x}) - \ell(\boldsymbol{\theta}_0; \mathbf{x}) \right]
$$


### The Rao Score Test
The <strong>score test</strong> rejects $H_0$ when the curvature of the log-likelihood function is large at $\boldsymbol{\theta}_0$. The statistic is given by:

$$
R_n = \frac{1}{n} U_{\boldsymbol{\theta}}^\top(\boldsymbol{\theta}_0; \mathbf{x}) \mathcal{I}^{-1}(\boldsymbol{\theta}_0) U_{\boldsymbol{\theta}}(\boldsymbol{\theta}_0; \mathbf{x})
$$

If Conditions 1-5 and 6 hold, then:

$$
\begin{aligned}
\frac{1}{\sqrt{n}} U_{\boldsymbol{\theta}}(\boldsymbol{\theta}_0; \mathbf{x}) 
  &= \sqrt{n} \frac{1}{n} \sum_{i = 1}^n \frac{\partial \log (p_{\boldsymbol{\theta}_0}(x_i))}{\partial \boldsymbol{\theta}}
  \rightsquigarrow \mathcal{N}\left( \mathbf{0}_k, \mathcal{I}(\boldsymbol{\theta}_0) \right)
\end{aligned}
$$

where the convergence in distribution follows from an application of the central limit theorem. Thus, $R_n \rightsquigarrow \chi^2_k$, and we reject if $R_n \geq \nu_\alpha$. 

We come to the following result concerning the Big Three tests.

<div class="theorem">
  <strong>Theorem 7.7.3.</strong>
  {% tabs thrm-7-7-3 %}
  {% tab thrm-7-7-3 statement %}
  Assume that Conditions 1-8 hold. Then $W_n$, $R_n$, and $\Delta_n$ are asymptotically equivalent under $H_0$ and have asymptotic level $\alpha$. 
  {% endtab %}
  {% tab thrm-7-7-3 proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>

---

## Composite Hypotheses
Suppose we are instead interested in testing, for $1 \leq r < k$:

$$
\begin{equation}
\label{eq:comp-hypo}
H_0: \boldsymbol{\theta}_1 = \boldsymbol{\theta}_{0, 1}, \dots, \boldsymbol{\theta}_r = \boldsymbol{\theta}_{0, r}
\hspace{5mm} vs. \hspace{5mm}
H_1: \boldsymbol{\theta}_i \neq \boldsymbol{\theat}_{0, i} \hspace{2mm} \text{for at least one } i = 1, \dots, r
\end{equation}
$$

As before, let $\hat{\boldsymbol{\theta}}$ be a solution to $U_{\boldsymbol{\theta}}(\hat{\boldsymbol{\theta}}; \mathbf{x}) = \mathbf{0}_k$ that is consistent for $\boldsymbol{\theta}$. 

### Likelihood Ratio Test
Let $\bar{\boldsymbol{\theta}}_{r+1}, \dots, \bar{\boldsymbol{\theta}}_k$ be the respective solutions to $U_{\boldsymbol{\theta}_{r+1}}(\bar{\boldsymbol{\theta}}; \mathbf{x}) = \mathbf{0}_k$ where $\bar{\boldsymbol{\theta}} = (\boldsymbol{\theta}_{0, 1}, \dots, \boldsymbol{\theta}_{0, r}, \bar{\boldsymbol{\theta}}_{r+1}, \dots, \bar{\boldsymbol{\theta}}_k)^\top$. The likelihood ratio test statistic for Eq. \eqref{eq:comp-hypo} is:

$$
\Delta_n = 2\left[ \ell(\hat{\boldsymbol{\theta}}; \mathbf{x}) - \ell(\bar{\boldsymbol{\theta}}; \mathbf{x}) \right]
$$


### Wald Test
Let $\boldsymbol{\Sigma}^{(r)}(\boldsymbol{\theta})$ be the $r \times r$ sub-matrix of $\boldsymbol{\Sigma}(\boldsymbol{\theta})$ corresponding to the asymptotic covariances of the first $r$ coordinates of $\boldsymbol{\theta}$. The Wald test statistic then extends to:

$$
W_n = n  \left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}\right)^\top \left[ \boldsymbol{\Sigma}^{(r)}(\boldsymbol{\theta}) \right]^{-1} \left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}\right)
$$


Let $\mathbf{v}_{1:r}$ denote the first $r$ coordinates of vector $\mathbf{v}$. It can be shown that:

$$
\sqrt{n}\left(\hat{\boldsymbol{\theta}}_{1:r} - \boldsymbol{\theta}_{1:r} \right) \rightsquigarrow \mathcal{N}\left(\mathbf{0}_r, \boldsymbol{\Sigma}^{(r)}(\boldsymbol{\theta}) \right)
$$




### Score Test

It turns out that the three tests are still asymptotically equivalent! First, note that the LRT, which is defined in the same way as before, has a similar limiting distribution (under appropriate conditions).

<div class="theorem">
  <strong>Theorem 7.7.4.</strong>
  {% tabs thrm-7-7-4 %}
  {% tab thrm-7-7-4 statement %}
  Assume that Conditions 1-8 hold and that $\hat{\boldsymbol{\theta}}$ is a solution to $U_{\boldsymbol{\theta}}(\hat{\boldsymbol{\theta}}; \mathbf{x}) = \mathbf{0}_k$ that is consistent for $\boldsymbol{\theta}$. Further assume that Conditions 1-8 hold for the vector $\boldsymbol{\theta}' = (\boldsymbol{\theta}_{0,1}, \dots, \boldsymbol{\theta}_{0, r}, \boldsymbol{\theta}_{r+1}, \dots, \boldsymbol{\theta}_k)^\top$. Let $\bar{\boldsymbol{\theta}}_{r+1}, \dots, \bar{\boldsymbol{\theta}}_k$ be the respective solutions to $U_{\boldsymbol{\theta}_{r+1}}(\bar{\boldsymbol{\theta}}; \mathbf{x}) = \mathbf{0}_k$ where $\bar{\boldsymbol{\theta}} = (\boldsymbol{\theta}_{0, 1}, \dots, \boldsymbol{\theta}_{0, r}, \bar{\boldsymbol{\theta}}_{r+1}, \dots, \bar{\boldsymbol{\theta}}_k)^\top$. 

  Then, under $H_0$, $\Delta_n \rightsquigarrow \chi^2_r$. 
  {% endtab %}
  {% tab thrm-7-7-4 proof %}
  Proof to be completed.
  {% endtab %}
  {% endtabs %}
</div>
