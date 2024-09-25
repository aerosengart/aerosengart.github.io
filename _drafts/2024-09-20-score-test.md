---
layout: post
title:  "Rao Score Test"
date: 20 September 2024
categories: posts
---

I've been reading about the score test, but I've found the [Wikipedia page](https://en.wikipedia.org/wiki/Score_test) to be a bit lacking. 

What better way to clarify my understanding than to go back to the beginning and read C. R. Rao's original work from 1948?

## Problem Set-Up

Let $x_1, \dots, x_p; y_1, \dots, y_q; \dots$ be observations with probability densities $f_1(x; \theta), f_2(y; \theta), \dots$ where $\theta \in \mathbb{R}^k$ is some vector-valued parameter. Further assume that each density function $f_j(\cdot; \theta)$ depends on at least one component of $\theta$. D

Due to independence, the likelihood is given by $\mathcal{L}(\theta; x, y, \dots) = f_1(x; \theta) \times f_2(y; \theta) \times \dots$. Taking the natural logarithm yields the log-likelihood, which is the logarithm of the joint density: $\ell(\theta) = \log \left( \mathcal{L}(\theta) \right) = \log \left( f(x, y, \dots; \theta) \right)$. 

We define the _score_ as the first-order gradient of the log-likelihood: $\phi_i = \frac{\partial \log \mathcal{L}(\theta)}{\partial \theta_i}$ for $i = 1, \dots, k$.


Let $A := \text{Cov}(\phi\phi^\top)$ denote the covariance matrix of the scores, and let $\alpha^{i,j}$ denote the $i,j$-th entry of its inverse. It can also be shown that $\mathbb{E}[\phi_i] = 0$ for all $i$. 

---
<details>
  <summary>Proof.</summary>
  There's a little bit of abuse of notation...sorry.
  $$
  \begin{align}
  \mathbb{E}[\phi_i] &= \mathbb{E}\left[ \frac{\partial \log \mathcal{L}}{\partial \theta_i} \right] \\
  				     &= \int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm} 
                  \frac{1}{f(x, y, \dots; \theta)} \frac{\partial}{\partial \theta_i} \log \mathcal{L}\theta \hspace{2mm} 
                  dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &= \int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm}
                  \frac{1}{f(x, y, \dots; \theta)} \frac{\partial f(x, y, \dots; \theta)}{\partial \theta_i} 
                  f(x, y, \dots; \theta) \hspace{2mm} 
                  dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &= int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm}
                  \frac{\partial f(x, y, \dots; \theta)}{\partial \theta_i} \hspace{2mm} dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &\overset{(i)}{=} \frac{\partial}{\partial \theta_i} \int_{\mathbb{R}^p}\int_{\mathbb{R}^q}\cdots \hspace{2mm} f(x, y, \dots; \theta) \hspace{2mm} dx \hspace{2mm} dy \hspace{2mm}d\dots \\
               &\overset{(ii)}{=} \frac{\partial}{\partial \theta_i} 1\\
               &= 0
  \end{align}
  $$
  In $(i)$ we use the fact that, under [certain regularity condition](https://en.wikipedia.org/wiki/Leibniz_integral_rule), we can interchange the order of integration and differentiation. In $(ii)$ we use the fact that the integral of a probability density function evaluates to $1$.
</details>
---
Under certain conditions, we can also show that the score vector at the true value of the parameter is multivariate Normal with mean zero and covariance matrix $A$.

---
<details>
  <summary>Proof.</summary>
  As a function of independent random variables, the score is also a random variable. In fact, it is the sum of independent random variables. Let $f_j(\theta)$ denote the density of random variable $j$ (i.e. for $j = 1$, we have $f_j(\theta) = f_1(x; \theta)$). Let $\ell_j(\theta)$ denote the log-likelihood of $\theta$ using only observations of random variable $j$. We then have:
  $$
  \begin{align}
  \phi_i &= \frac{\partial \log \mathcal{L}(\theta)}{\partial \theta_i} \\
         &= \frac{\partial \log \prod_j f_j(\theta)}{\partial \theta_i} \\
         &= \frac{\partial \sum_j \log f_j(\theta)}{\partial \theta_i} \\
         &= \sum_j \frac{\partial \log f_j(\theta)}{\partial \theta_i} \\
         &= \sum_j \frac{\partial \ell_j(\theta)}{\partial \theta_i}
  \end{align}
  $$
  As shown above, $\mathbb{E}[\phi_i] = 0 < \infty$. Assume that $A$ has finite diagonal elements.
   
  If we further assume that, for some $\delta > 0$, 
  $$
  \underset{n \rightarrow \infty}{\lim} \frac{1}{\alpha_{i,i}} \sum_j \mathbb{E}[\rvert \phi_i - \mathbb{E}[\phi_i] \rvert^{2+\delta}] = 0
  $$
  then the Lyapunov CLT states that, as $n \rightarrow \infty$, $\sum_{j} (\phi_i - \mathbb{E}[\phi_i])$ converges in distribution to a multivariate Normal distribution with mean zero and covariance $A$.
</details>
---

This implies that the statistic $\chi^2 = \sum_i \sum_j \alpha^{i,j} \phi_i \phi_j$ is asymptotically $\chi^2$ with $k$ degrees of freedom (since $\theta$ is $k$-dimensional). 

We now assume that we have $s$ restrictions on our parameters: $\psi_i(\theta_1, \dots, \theta_k) = 0$ for $i = 1, 2, \dots, s$. Using the method of Lagrange multipliers, the maximum likelihood estimates of the parameter, denoted by $\hat{\theta}_1, \dots, \hat{\theta}_k$ is:
$$
\begin{align}
&\phi_i + \sum_j \lambda_j \frac{\partial \psi_j}{\partial \theta_i} = 0
&\psi_l = 0
\end{align}
$$
for $i = 1, 2, \dots, k$ and $l = 1, 2, \dots, s$.








---

C. R. Rao, “Large Sample Tests of Statistical Hypotheses Concerning Several Parameters with Applications to Problems of Estimation,” _Proceedings of the Cambridge Philosophical Society_, Vol. 44, No. 1, 1948, pp. 50-57. [doi:10.1017/S0305004100023987](https://doi.org/10.1017/S0305004100023987).