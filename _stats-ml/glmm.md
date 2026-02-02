---
layout: distill
title: Generalized Linear Mixed Models
description: A Primer
date: 2025-06-04
tabs: true
tags: regression likelihood models primer
toc:
  - name: Background
  - name: Set-Up
  - name: Likelihood
  - name: Quasi-Likelihood
    subsections:
        - name: Penalized Quasi-Likelihood
        - name: Marginal Quasi-Likelihood
        - name: An Aside

bibliography: 2025-06-04-glmm.bib
---

This post is a primer on generalized linear mixed models (GLMM) and their associated estimation procedures. Throughout this post, we perform derivatives of matrices, vectors, and scalars. We denote all of them with $\partial$ even though they may not <i>technically</i> be partial derivatives. The meaning should be clear from the context if you keep track of the dimensions of the variables. 

We'll use the notation $\text{vec}(\mathbf{M})$ denote the vectorization of matrix $\mathbf{M}$, where we flatten the matrix row-by-row. Also, any function $f: \mathbb{R} \rightarrow \mathbb{R}$ will be applied element-wise to vectors and matrices.

---

## Background
Let's first explain <i>why</i> generalized linear mixed models are of interest and helpful in data analysis. Suppose we have a dataset in which the observations can be grouped in some way (i.e. we have repeated observations from the same unit). One example is longitudinal studies where we obtain measurements from the same individuals over time. Another example is in education where we may observe students from the same classes or schools. In these settings, it seems reasonable to assume that observations from the same unit would be more similar that those from different ones. This can be realized in our model by including a way for the effects of covariates to differ across units. 

As in generalized linear models, we assume the mean response is some function of a linear combination of covariates. However, GLMMs are called <i>mixed</i> because the effects of the covariates are of two types: <i>fixed</i> and <i>random</i>. The fixed effects represent population-level relationships, while the random effects represent unit-specific ones. In this way, we can use GLMMs to analyze between-subject variation (via the fixed effects) and within-subject variation (via the random effects): the fixed effects determine the relationship between the covariates and the mean response for the population (i.e. overall), and the random effects describe how each group's mean response deviates from that. Later on, we introduce measurement errors, which account for deviation of each <i>observation</i> from its group's mean. 

---

## Set-Up
We have $n$ observations coming from $k$ different clusters, each of size $n_t$ for $t \in [k]$. The full data will be denoted by $\mathbf{y}$. Though $\mathbf{y}$ is a vector, we'll denote the $j$-th observation from cluster $i$ with $\mathbf{y}_{i,j}$. For example, $$\mathbf{y}_{i,j}$$ denotes element $$\sum_{l = 1}^{i - 1} n_l + j$$ of $\mathbf{y}$. We'll denote the $n_i$-dimensional vector of responses for cluster $i$ with $\mathbf{y}_i$. 

For each observation, we will have $p$ fixed effect covariates arranged in a $p$-dimensional vector, $$\mathbf{x}_{i, j}$$, and $q$ random effects covariates in a $q$-dimensional vector, $$\mathbf{z}_{i,j}$$. We'll assume that the observations within the same cluster are independent.

We'll let $\alpha$ denote a $p$-dimensional fixed effect coefficient vector, $\alpha$, and $\beta_t$ will denote a $q$-dimensional random effect coefficient vector corresponding to group $t$. We assume $\beta_t \overset{iid}{\sim} \mathcal{F}$ for all $t \in [k]$ for an exponential family distribution, $\mathcal{F}$, and $q \times q$ covariance matrix, $D(\tau^2)$, depending on an $m$-dimensional variance component vector, $\tau^2$. These random effects are assumed to be independent of the covariates. When they are not, we can run into issues with bias when estimating the fixed effects.

<details>
<summary><strong>Additional Assumptions.</strong></summary>
A few auxiliary assumptions must also be made for the analysis later, which we list here:
<ul>
  <li>The third moment and higher moments of $\beta$ are of order $o(\rvert \rvert \tau^2 \rvert \rvert)$.</li>
  <li>The entries of $D(\tau^2)$ are linear in $\tau^2$.</li>
  <li>$D(\tau^2) = \text{vec}(0)$ if $\tau^2 = \text{vec}(0)$</li>
</ul>
</details>

Our model comes in the form of a specification of the conditional mean, $\mu_{i,j} = \mathbb{E}[\mathbf{y}_{i,j} \rvert \beta_i]$ (where we suppress the addition conditioning on the covariates themselves). For a monotonic and differentiable link function (e.g. $\log(\cdot)$ or $\text{logit}(\cdot)$), the conditional mean of the $j$-th observation in group $i$ is assumed to be given by:

$$
\begin{equation}
\label{eq:glmm}
\mu_{i,j} = g^{-1}\left(\mathbf{x}_{i,j}^\top \alpha + \mathbf{z}_{i,j}^\top \beta_i \right)
\end{equation}
$$

We then assume that the observations themselves follow some exponential family distribution with measurement errors, $\epsilon_{i,j}$, which is the deviation of the response from its (unit-specific) conditional mean. These errors are assumed to have mean zero and be independent of each other and of the random effects. We further assume the responses, $\mathbf{y}_{i,j}$, conditional on the random effects (and the covariates), are independent with variances equal to some function of the conditional mean. 

To write Eq. \eqref{eq:glmm} in matrix notation, we assume that the observations are ordered by group, so the first $n_1$ observations are all from group $1$. We can then define the following vectors/matrices:

$$
\mu = 
\begin{bmatrix}
\mu_{1, 1} \\
\vdots \\
\mu_{k, n_k}
\end{bmatrix},
\hspace{2mm}
\mathbf{X} = 
\begin{bmatrix}
 — & \mathbf{x}_{1,1} & — \\
 & \dots &  \\
 — & \mathbf{x}_{k, n_k} & —  \\
\end{bmatrix},
\hspace{2mm}
\tilde{\mathbf{Z}}_i = 
\begin{bmatrix}
 — & \mathbf{z}_{i, 1} & — \\
 & \dots &  \\
 — & \mathbf{z}_{i, n_i} & —  \\
\end{bmatrix},
\hspace{2mm}
\mathbf{Z} = 
\begin{bmatrix}
\tilde{\mathbf{Z}}_1 & \dots & 0 \\
\vdots & \ddots & \vdots \\
0 & \dots & \tilde{\mathbf{Z}}_k
\end{bmatrix},
\hspace{2mm}
\beta = 
\begin{bmatrix}
\beta_1  \\
\vdots \\
\beta_k 
\end{bmatrix}
$$

where $\tilde{\mathbf{Z}}_i$ is constructed by stacking the $$\mathbf{z}_{i,j}$$ vectors for group $i$ into a matrix, and where $\beta$ is the vector created by stacking together all of the $\beta_i$ vectors vertically. As a reminder, $$\mu \in \mathbb{R}^{n \times 1}$$, $$\mathbf{X} \in \mathbb{R}^{n \times p}$$, $$\tilde{\mathbf{Z}}_i \in \mathbb{R}^{n_i \times q}$$, $$\mathbf{Z} \in \mathbb{R}^{n \times (q \times k)}$$, $\alpha \in \mathbb{R}^p$, and $$\beta \in \mathbb{R}^{(q \times k) \times 1}$$.

The model is then:

$$
\begin{equation}
\label{eq:glmm-matrix}
\mu = g^{-1}(\mathbf{X} \alpha + \mathbf{Z} \beta)
\end{equation}
$$

We'll use $[\cdot] \rvert_{\beta = \beta_0}$ to denote evaluation of the function in brackets when setting $\beta$ equal to $\beta_0$. Similarly, we'll use $[\cdot] \rvert_{H_0}$ to denote evaluation under the null hypothesis. We'll also use a superscript $0$ (e.g. $\mu^0$, $\eta^0$, etc.) to denote the quantity under the null hypothesis (i.e. $\tau^2 = \mathbf{0} \implies \beta = \mathbf{0}$).


<div class="example">
<strong>Example.</strong>
<br>
To keep things simple, we'll assume to a fixed intercept (that is cluster-specific) and a random slope. Thus, $\alpha, \beta \in \mathbb{R}^{k}$, and $\mathbf{z}_{i, j} \in \mathbb{R}$ as well. 
<br>
In this simple case, $\tilde{\mathbf{Z}}_i \in \mathbb{R}^{n_t}$, so $\mathbf{Z}$ is $n \times k$. Our model is then:
$$
\mathbf{y}_{i,j} \rvert \beta_i \sim \text{Poi}\left(\mu_{i,j}\right), 
\hspace{5mm}
\mu_{i,j} = \exp\left(\alpha_i + \beta_i \mathbf{z}_{i,j} \right)
$$
We'll have a scalar-valued variance component that we call $\tau^2$. The random effects will have distribution:
$$
\beta_i \overset{iid}{\sim} \mathcal{N}(0, \tau^2), \hspace{5mm} \forall i \in [k]
$$
If we let $\mathbf{A}$ denote the $n$-dimensional vector of intercepts where the $i$-th entry of $\alpha$ is repeated $n_i$ times, then we can write our model in vector form as:
$$
\mu = \exp(\mathbf{A} + \mathbf{Z} \beta)
$$
</div>

---

## Likelihood
We can write the conditional log-likelihood using the exponential family form (see my <a href="/stats-ml/glm">generalized linear models post</a>). 

$$
\begin{equation}
\label{eq:condition-log-lik}
\begin{aligned}
\ell(\mathbf{y}; \alpha, \tau^2 \rvert \beta) &= \sum_{i = 1}^{k} \sum_{j = 1}^{n_i} \left[ \frac{\zeta_{i,j} \mathbf{y}_{i,j} - A(\zeta_{i,j})}{d(\phi, \omega_{i,j})} + \log(h(\mathbf{y}_{i,j}, \phi, \omega_{i,j}))\right] \\
\mathcal{L}(\mathbf{y}; \alpha, \tau^2 \rvert \beta) &= \exp \left( \sum_{i = 1}^n \left[ \frac{\zeta_i \mathbf{y}_i - A(\zeta_i)}{d(\phi, \omega_i)} + \log(h(\mathbf{y}_i, \phi, \omega_i))\right] \right)
\end{aligned}
\end{equation}
$$

where $\phi > 0$ is a dispersion/scale parameter; $\omega_{i,j}$ is a (prior) dispersion weights; $\zeta_{i,j}$ is a distribution parameter; and $h(\cdot)$ and $A(\cdot)$ are known functions. We assume that $d(\phi, \omega_{i,j}) = \phi \omega_i^{-1}$. 

The conditional variance of the responses is given by:

$$
\text{Cov}(\mathbf{y} \rvert \beta) = \text{diag}(d(\phi, \omega)) \underbrace{\frac{\partial^2 A(\tau^2)}{\partial \tau^2 \partial (\tau^2)^\top}}_{=V(\mu)} = \text{diag}^{-1}\left( \frac{\omega}{\phi}\right) V(\mu)
$$

---

## Quasi-Likelihood
We don't actually need to specify the true likelihood, even though (in practice) we will usually be able to do so. Using quasi-likelihood methods can often be less computationally expensive than maximum likelihood for count data. 

In the quasi-likelihood scheme, we keep the conditional mean assumption that we discussed above and also make the assumption that the conditional variance can be expressed as a function of the mean in the form $\text{Cov}(\mathbf{y}) = \text{diag}\left(d(\phi, \omega)\right)V(\mu)$ for some function $V(\cdot)$. Later, we'll denote diagonal elements of $\text{Cov}(\mathbf{y})$ with $v(\mu_{i,j})$. 

The conditional log quasi-likelihood and quasi-likelihood are given by:

$$
\begin{equation}
\label{eq:quasi-lik}
\begin{aligned}
\ell_q(\mathbf{y}_{i,j}; \alpha, \tau^2 \rvert \beta_i) &= \int_{\mathbf{y}_{i,j}}^{\mu_{i,j}} \frac{\omega_{i,j}(\mathbf{y}_{i,j} - u)}{\phi V(u)} du \\
\mathcal{L}_q(\mathbf{y}_{i,j}; \alpha, \tau^2 \rvert \beta_i) &= \exp\left(\int_{\mathbf{y}_{i,j}}^{\mu_{i,j}} \frac{\omega_{i,j}(\mathbf{y}_{i,j} - u)}{\phi V(u)} du \right)
\end{aligned}
\end{equation}
$$

When $\mathcal{F}$ <i>is</i> an exponential family distribution, the log-likelihood and log quasi-likelihood are equal. See my <a href="/stats-ml/quasi-likelihood">post on quasi-likelihood</a> for a more in-depth discussion.

Let $f$ denote the density associated with $\mathcal{F}$, the distribution of the random effects. The unconditional quasi-likelihood is then given by integrating out the random effects from the joint quasi-likelihood:

$$
\begin{equation}
\label{eq:uncondition-log-lik}
\begin{aligned}
\mathcal{L}_q(\mathbf{y}; \alpha, \tau^2) &= \prod_{i = 1}^k \mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2) =  \prod_{i = 1}^k \int\prod_{j = 1}^{n_i} \mathcal{L}_q(\mathbf{y}_{i,j}; \alpha, \tau^2 \rvert \beta_i) f(\beta_i)  d\beta_i
\end{aligned}
\end{equation}
$$

Remember that the above integral is multi-dimensional if $\beta_i$ has dimension greater than $1$!

We usually do maximum likelihood estimation for the parameter vector $\beta$ by taking the derivative of the log of the above, setting it equal to zero, and solving for $\beta$. However, depending upon the form of $f(\cdot)$ and the number of random effects we have, Eq. \eqref{eq:uncondition-log-lik} may be a multi-dimensional integral that is difficult to evaluate, let alone differentiate. One way is to use quadrature to approximate the integral as a summation. However, there are two other, computationally more feasible, methods that provide a work-around. 

There are two common methods of approximate inference for GLMMs: penalized quasi-likelihood (PQL) and marginal quasi-likelihood (MQL). They are very similar in that they both perform a Taylor approximation of the conditional log quasi-likelihood to evaluate the integral in Eq. \eqref{eq:uncondition-log-lik}. However, PQL uses the maximum a priori estimates of the random effects as the operating point for the expansion, while MQL uses the random effects mean. See below for a discussion of the approaches.

### Penalized Quasi-Likelihood
Penalized quasi-likelihood (PQL) essentially uses Laplace's method to approximate the integral above. The general idea behind <a href="https://en.wikipedia.org/wiki/Laplace%27s_method">Laplace's method</a> is to approximate integrals of a certain form as:

$$
\int \exp(M f(x)) dx \approx \exp(M f(x_0)) \int \exp\left( - \frac{1}{2} M \rvert f''(x_0) \rvert (x - x_0)^2 \right) dx
$$

where $M$ is a big scalar, $f$ is a twice-differentiable function, and $x_0$ is a global maximum of $f$. If $$f''(x_0) < 0$$ and $M \rightarrow \infty$, then the integrand above is basically a Gaussian kernel and the approximation becomes:

$$
\begin{equation}
\label{eq:ql-gauss}
\int \exp(M f(x)) dx \approx \left(\frac{2 \pi}{M \rvert f''(x_0) \rvert} \right)^{\frac{1}{2}} \exp(M f(x_0))
\end{equation}
$$

<div class="theorem">
<strong>Claim (Rewriting the Integral).</strong>
{% tabs claim-3 %}
{% tab claim-3 statement %}
Letting $c = (2 \pi)^{-\frac{m}{2}}$ and $d_{i,j}(y_{i,j}, \mu_{i,j}) = -2 \int_{y_{i,j}}^{\mu_{i,j}} \frac{y_{i,j} - z}{a_{i,j} V(z)} dz$, we can rewrite the integral as:

$$
\begin{aligned}
\mathcal{L}_q(\mathbf{y}; \alpha, \tau^2) &=  c \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \int \exp\left( -\kappa(\beta) \right) d \beta; \\
\kappa(\beta) &= \frac{1}{2\phi}\sum_{i = 1}^n \sum_{j = 1}^{n_i} d_{i,j}(y_{i,j}, \mu_{i,j}) + \frac{1}{2}\beta^\top D^{-1}(\tau^2) \beta
\end{aligned}
$$
{% endtab %}
{% tab claim-3 proof %}
$$
\begin{aligned}
\mathcal{L}_q(\mathbf{y}; \alpha, \tau^2) &= 
\int \mathcal{L}_q(\mathbf{y}; \alpha, \tau^2 \rvert \beta) \mathcal{L}(\beta) d \beta \\
&= \int \left[ \prod_{i = 1}^n \prod_{j = 1}^{n_i} \exp\left( \ell_q(y_{i,j}; \mu_{i,j} \rvert \beta) \right) \right] (2 \pi)^{-\frac{m}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}}  \exp\left(- \frac{1}{2}\beta^\top D^{-1}(\tau^2) \beta \right) d\beta \\
&= (2 \pi)^{-\frac{m}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \int \exp\left(\sum_{i = 1}^n  \sum_{j = 1}^{n_i} \left( \int_{y_{i,j}}^{\mu_{i,j}} \frac{y_{i,j} - z}{\phi a_{i,j} V(z)} dz \right) - \frac{1}{2}\beta^\top D^{-1}(\tau^2) \beta \right) d\beta \\
&= (2 \pi)^{-\frac{m}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \int \exp\left(- \frac{1}{2 \phi} \sum_{i = 1}^n \sum_{j = 1}^{n_i} \left( -2 \int_{y_{i,j}}^{\mu_{i,j}} \frac{y_{i,j} - z}{a_{i,j} V(z)} dz \right) - \frac{1}{2}\beta^\top D^{-1}(\tau^2) \beta \right) d\beta \\
&= (2 \pi)^{-\frac{m}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \int \exp\left(- \frac{1}{2 \phi} \sum_{i = 1}^n \sum_{j = 1}^{n_i} d_{i,j}(y_{i,j}, \mu_{i,j}) - \frac{1}{2}\beta^\top D^{-1}(\tau^2) \beta \right) d\beta & \left( d_{i,j}(y_{i,j}, \mu_{i,j}) = -2 \int_{y_{i,j}}^{\mu_{i,j}} \frac{y_{i,j} - z}{a_{i,j} V(z)} dz \right)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Note that the first term in $\kappa(\beta)$ <i>does</i> involve $\beta$ since $$\mu_{i,j} = g^{-1}(\mathbf{x}_{i,j}^\top \alpha + \mathbf{z}_{i,j}^\top \beta_j)$$. Let $\kappa'$ and $\kappa''$ denote the vector of first-order partial derivatives and the matrix of second-order partial derivatives, respectively, of $\kappa$ with respect to $\beta$. 

<div class="theorem">
<strong>Claim (Partial Derivatives).</strong>
{% tabs claim-4 %}
{% tab claim-4 statement %}
These have the form:

$$
\begin{aligned}
\kappa'_k(\beta) &= -\sum_{j = 1}^{n_k} \left(\frac{y_{k,j} - \mu_{k,j}}{\phi a_{k,j} V(\mu_{k,j})\frac{\partial g(\mu_{i,j})}{\partial \mu_{i,j}}}\right)\mathbf{z}_{i,j} + D^{-1}(\tau^2) \beta_k  \\
\kappa''(\beta) &\approx \tilde{\mathbf{Z}}^\top \mathbf{W} \tilde{\mathbf{Z}} + D^{-1}(\tau^2)
\end{aligned}
$$

where $\mathbf{W}$ is the $n \times n$ diagonal matrix with elements $\frac{1}{\phi a_{i,j} V(\mu_{i,j}) \frac{\partial g(\mu_{i,j})}{\partial \mu_{i,j}}}$, and $\tilde{\mathbf{Z}}$ is the $n \times q$ matrix formed by concatenating the $\tilde{\mathbf{Z}}^t$ matrices.
{% endtab %}
{% tab claim-4 proof %}
$$
\begin{aligned}
\kappa'_k(\beta) &= \frac{\partial}{\partial \beta_k} \left[\frac{1}{2 \phi} \sum_{i = 1}^k \sum_{j = 1}^{n_i} d_{i,j}(y_{i,j}, \mu_{i,j}) + \frac{1}{2} \beta^\top D^{-1}(\tau^2) \beta \right] \\
&= \frac{1}{2 \phi} \sum_{i = 1}^k \sum_{j = 1}^{n_i} \frac{\partial}{\partial \beta_k} \left[d_{i,j}(y_{i,j}, \mu_{i,j}) \right] + D^{-1}(\tau^2) \beta  \\
&= \frac{1}{2 \phi} \sum_{j = 1}^{n_k} \frac{\partial \mu_{k,j}}{\partial \beta_k} \frac{\partial}{\partial \mu_{k,j}}\left[ -2 \int_{y_{k,j}}^{\mu_{k,j}} \frac{y_{k,j} - z}{a_{k,j} V(z)}dz \right] + D^{-1}(\tau^2) \beta_k  \\
&= \frac{1}{2 \phi} \sum_{j = 1}^{n_k} \frac{\partial \mu_{k,j}}{\partial \beta_k} \left(-2  \frac{y_{k,j} - \mu_{k,j}}{a_{k,j} V(\mu_{k,j})}\right) + D^{-1}(\tau^2) \beta_k  \\
&= -\sum_{j = 1}^{n_k} \frac{\partial \mu_{k,j}}{\partial \beta_k} \left(\frac{y_{k,j} - \mu_{k,j}}{\phi a_{k,j} V(\mu_{k,j})}\right) + D^{-1}(\tau^2) \beta_k \\
\end{aligned}
$$

Assuming $g(\cdot)$ is <i>strictly monotone and continuous</i>, then we have:

$$
\frac{\partial g^{-1}(\eta_{i,j})}{\partial \eta_{i,j}} = \frac{\partial \mu_{i,j}}{\partial \eta_{i,j}} = \left(\frac{\partial \eta_{i,j}}{\partial \mu_{i,j}}\right)^{-1} = \left(\frac{\partial g(\mu_{i,j})}{\partial \mu_{i,j}}\right)^{-1}
$$

This implies:

$$
\frac{\partial \mu_{k,j}}{\partial \beta_k} = \frac{\partial \eta_{k,j}}{\partial \beta_k} \frac{\partial \mu_{k, j}}{\partial \eta_{i,j}} = \mathbf{z}_{i,j} \left(\frac{\partial \eta_{i,j}}{\partial \mu_{i,j}}\right)^{-1} = \frac{\mathbf{z}_{i,j}}{\frac{\partial g(\mu_{i,j})}{\partial \mu_{i,j}}}
$$

Thus, $\kappa'_k(\beta)$ is:

$$
\kappa'_k(\beta) = -\sum_{j = 1}^{n_k} \left(\frac{y_{k,j} - \mu_{k,j}}{\phi a_{k,j} V(\mu_{k,j})\frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}}\right)\mathbf{z}_{k,j} + D^{-1}(\tau^2) \beta_k 
$$

Looking at the second-order partial derivatives:

$$
\begin{aligned}
\kappa''_{k,k}(\beta) &= 
\frac{\partial}{\partial \beta_k^\top} \left[ -\sum_{j = 1}^{n_k} \left(\frac{y_{k,j} - \mu_{k,j}}{\phi a_{k,j} V(\mu_{k,j})\frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}}\right)\mathbf{z}_{k,j} + D^{-1}(\tau^2) \beta_k  \right] \\
&= - \sum_{j = 1}^{n_k} \left( \frac{\mathbf{z}_{k,j}}{\phi a_{k,j} V(\mu_{k,j}) \frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}} \frac{\partial}{\partial \beta_k^\top} \left[ y_{k,j} - \mu_{k,j} \right] + \mathbf{z}_{k,j}(y_{k,j} - \mu_{k,j})\frac{\partial}{\partial \beta_k^\top} \left[ \frac{1}{\phi a_{k,j} V(\mu_{k,j}) \frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}} \right] \right) + D^{-1}(\tau^2) \\
&\overset{(i)}{=} - \sum_{j = 1}^{n_k} \left( \frac{-\mathbf{z}_{k,j}\mathbf{z}_{k, j}^\top}{\phi a_{k,j} V(\mu_{k,j}) \left(\frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}\right)^2} + \mathbf{z}_{k,j}(y_{k,j} - \mu_{k,j})\frac{\partial}{\partial \beta_k^\top} \left[ \frac{1}{\phi a_{k,j} V(\mu_{k,j}) \frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}} \right] \right) + D^{-1}(\tau^2) \\
&\overset{(ii)}{=} \sum_{j = 1}^{n_k} \frac{\mathbf{z}_{k,j}\mathbf{z}_{k, j}^\top}{\phi a_{k,j} V(\mu_{k,j}) \left(\frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}\right)^2} + D^{-1}(\tau^2) + R_k \\
&=(\tilde{\mathbf{Z}}^k)^\top \mathbf{W}_k \tilde{\mathbf{Z}}^k + D^{-1}(\tau^2) + R_k
\end{aligned}
$$

where in $(i)$, we set:

$$
R_k = -\sum_{j = 1}^{n_k} (y_{k,j} - \mu_{k,j}) \mathbf{z}_{k,j} \frac{\partial}{\partial \beta_k^\top} \left[ \frac{1}{\phi a_{k,j} V(\mu_{k,j}) \frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}} \right]
$$

and in $(ii)$, $\mathbf{W}_k$ is the diagonal matrix with elements:

$$
\frac{1}{\phi a_{k,j} V(\mu_{k,j}) \frac{\partial g(\mu_{k,j})}{\partial \mu_{k,j}}}
$$

Since $$\kappa'_k(\beta)$$ only involves $\beta_k$, $$\kappa''_{k, k'}(\beta) = \mathbf{0}$$. Stacking all of the $\mathbf{W}_k$ and $R_k$ together, we get:

$$
\kappa''(\beta) = \tilde{\mathbf{Z}}^\top \mathbf{W} \tilde{\mathbf{Z}} + D^{-1}(\tau^2) + R
\approx  \tilde{\mathbf{Z}}^\top \mathbf{W} \tilde{\mathbf{Z}} + D^{-1}(\tau^2)
$$

where we drop $R$ for the approximation. 
{% endtab %}
{% endtabs %}
</div>


Letting $\hat{\beta}$ be a solution to the equation $\kappa'(\beta) = \mathbf{0}$, we can ignore constants and remainders (by assuming they are negligible in the long run...) to obtain the approximation in Eq. \eqref{eq:ql-gauss}:

$$
\begin{aligned}
\ell_q(\mathbf{y}; \alpha, \tau^2) 
&\approx \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \log\left( \left( \frac{2 \pi}{\rvert \kappa''(\beta_0) \rvert}\right)^{\frac{1}{2}} \exp(-\kappa(\beta_0)) \right) \\
&= \frac{1}{2}\log(2 \pi) - \frac{1}{2} \log(\kappa''(\beta_0)) + \kappa(\beta_0) \\
&\approx -\frac{1}{2} \log\left( \rvert D(\tau^2) \rvert \right) -\frac{1}{2} \log\left( \rvert \tilde{\mathbf{Z}}^\top \mathbf{W} \mathbf{\tilde{Z}} + D^{-1}(\tau^2) \rvert \right) + \frac{1}{2\phi} \sum_{i = 1}^k \sum_{j = 1}^{n_i} d_{i,j}(y_{i,j}, \mu_{i,j}) + \frac{1}{2}\beta^\top D^{-1}(\tau^2) \beta
\end{aligned}
$$


### Marginal Quasi-Likelihood
An alternative way to think about a generalized linear mixed model is from a marginal perspective, where we instead focus on the <i>unconditional</i> expectation of the response:

$$
\mathbb{E}[y_{i,j}] = \mu_{i,j} = g(\mathbf{x}_{i,j}^\top \alpha)
$$

In contrast to PQL, marginal quasi-likelihood (MQL) performs a Taylor approximation of the log quasi-likelihood about the random effects (prior) mean $\beta = \mathbf{0}$ to get a closed form solution to the multi-dimensional integral. For the $i$-th cluster, this approximation is:

$$
\begin{aligned}
\ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \approx \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)\rvert_{\beta_i = \mathbf{0}} + \beta_i^\top \left[ \frac{\partial \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i} \right] \bigg\rvert_{\beta_i = \mathbf{0}} + \frac{1}{2} \beta_i^\top  \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}} \beta_i
\end{aligned}
$$

Thus:

$$
\begin{aligned}
\mathcal{L}_q(\mathbf{y}; \alpha, \tau^2) 
&= \prod_{i = 1}^k (2 \pi)^{-\frac{q}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \int \exp\left( \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) - \frac{1}{2}\beta_i^\top D^{-1}(\tau^2) \beta_i \right) d \beta_i \\
&\approx \prod_{i = 1}^k (2 \pi)^{-\frac{q}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \int \exp\left( [\ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)]_{\beta_i = \mathbf{0}} + \beta_i^\top \left[ \frac{\partial \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i} \right] \bigg\rvert_{\beta_i = \mathbf{0}} + \frac{1}{2} \beta_i^\top \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}} \beta_i - \frac{1}{2}\beta_i^\top D^{-1}(\tau^2) \beta_i \right) d \beta_i \\
&= \prod_{i = 1}^k (2 \pi)^{-\frac{q}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right] \bigg\rvert_{\beta_i = \mathbf{0}} \int  \exp\left( \beta_i^\top \left[ \frac{\partial \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i} \right] \bigg\rvert_{\beta_i = \mathbf{0}} - \frac{1}{2} \beta_i^\top \left[ D^{-1}(\tau^2) - \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}} \right] \beta_i\right) d \beta_i
\end{aligned}
$$

We can rewrite the above integral in <i>canonical form</i> by defining:

$$
\begin{aligned}
\mathbf{K}_i &= D^{-1}(\tau^2) - \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}} \\
\mathbf{h}_i &= \left[ \frac{\partial \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i} \right] \bigg\rvert_{\beta_i = \mathbf{0}}
\end{aligned}
$$

And:

$$
\bar{\Sigma}_i^{-1} = \mathbf{K}_i;
\hspace{5mm}
\bar{\mu}_i = \bar{\Sigma}_i \mathbf{h}_i;
\hspace{5mm}
g_i = - \frac{1}{2} \bar{\mu}_i^\top \bar{\Sigma}_i^{-1} \bar{\mu}_i - \log\left((2 \pi)^{\frac{q}{2}} \rvert \bar{\Sigma}_i \rvert^{\frac{1}{2}} \right)
$$

<div class="theorem">
<strong>Claim (Likelihood Approximation).</strong> 
{% tabs claim-5 %}
{% tab claim-5 statement %}
It follows that:

$$
\begin{aligned}
\mathcal{L}_q(\mathbf{y}; \alpha, \tau^2)
&\approx
\prod_{i = 1}^k \left\rvert \mathbb{I}_{q \times q} - D(\tau^2) \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}}\right\rvert^{-\frac{1}{2}}  \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right]\bigg\rvert_{\beta_i = \mathbf{0}}  \exp\left(\frac{1}{2} \left[\mathbf{h}_i \mathbf{K}_i^{-1} \mathbf{h}_i \right] \right) \\
\implies
\ell_q(\mathbf{y}; \alpha, \tau^2) 
&\approx
\sum_{i = 1}^k \left[ -\frac{1}{2} \log\left( \left\rvert \mathbb{I}_{q \times q} - D(\tau^2) \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}}\right\rvert \right) + \left[ \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right] \bigg\rvert_{\beta_i = \mathbf{0}} + \frac{1}{2} \mathbf{h}_i^\top \mathbf{K}_i^{-1} \mathbf{h}_i \right]
\end{aligned}
$$
{% endtab %}
{% tab claim-5 proof %}
$$
\begin{aligned}
 \mathcal{L}_q(\mathbf{y}; \alpha, \tau^2)
 &\approx 
 \prod_{i = 1}^k (2 \pi)^{-\frac{q}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right] \bigg\rvert_{\beta_i = \mathbf{0}} \underbrace{\int  \exp\left( \beta_i^\top \mathbf{h}_i - \frac{1}{2} \beta_i^\top \mathbf{K}_i \beta_i\right) d \beta_i}_{= \exp(-g_i)} \\
 &= \prod_{i = 1}^k (2 \pi)^{-\frac{q}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right] \bigg\rvert_{\beta_i = \mathbf{0}} \exp\left(\frac{1}{2}\bar{\mu}_i^\top \bar{\Sigma}_i^{-1}\bar{\mu}_i + \log\left((2 \pi)^{\frac{q}{2}} \rvert \bar{\Sigma}_i \rvert^{\frac{1}{2}}\right) \right) \\
 &= \prod_{i = 1}^k (2 \pi)^{-\frac{q}{2}} \rvert D(\tau^2) \rvert^{-\frac{1}{2}} \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right]\bigg\rvert_{\beta_i = \mathbf{0}}  \exp\left(\frac{1}{2}\bar{\mu}_i^\top \bar{\Sigma}_i^{-1}\bar{\mu}_i\right) (2 \pi)^{\frac{q}{2}} \rvert \bar{\Sigma}_i \rvert^{\frac{1}{2}} \\
 &= \prod_{i = 1}^k \rvert D(\tau^2)\rvert^{-\frac{1}{2}} \rvert  \bar{\Sigma}_i^{-1} \rvert^{-\frac{1}{2}}\left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right]\bigg\rvert_{\beta_i = \mathbf{0}}  \exp\left(\frac{1}{2}\bar{\mu}_i^\top \bar{\Sigma}_i^{-1}\bar{\mu}_i\right) \\
 &= \prod_{i = 1}^k \rvert D(\tau^2) \mathbf{K}_i \rvert^{-\frac{1}{2}}  \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right]\bigg\rvert_{\beta_i = \mathbf{0}}  \exp\left(\frac{1}{2} \left[\mathbf{h}_i (\mathbf{K}_i^{-1})^\top \mathbf{K}_i \mathbf{K}_i^{-1} \mathbf{h}_i \right] \right) \\
 &= \prod_{i = 1}^k \left\rvert D(\tau^2) \left[ D^{-1}(\tau^2) - \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}}\right] \right\rvert^{-\frac{1}{2}}  \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right]\bigg\rvert_{\beta_i = \mathbf{0}}  \exp\left(\frac{1}{2} \left[\mathbf{h}_i \mathbf{K}_i^{-1} \mathbf{h}_i \right] \right) \\
 &= \prod_{i = 1}^k \left\rvert \mathbb{I}_{q \times q} - D(\tau^2) \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}} \right\rvert^{-\frac{1}{2}}  \left[\mathcal{L}_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i) \right]\bigg\rvert_{\beta_i = \mathbf{0}}  \exp\left(\frac{1}{2}  \left[ \frac{\partial \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i^\top } \right] \bigg\rvert_{\beta_i = \mathbf{0}}  \left[ D^{-1}(\tau^2) - \left[ \frac{\partial^2 \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i \partial \beta_i^\top} \right] \bigg\rvert_{\beta_i = \mathbf{0}} \right]^{-1} \left[ \frac{\partial \ell_q(\mathbf{y}_i; \alpha, \tau^2 \rvert \beta_i)}{\partial \beta_i} \right] \bigg\rvert_{\beta_i = \mathbf{0}}\right)
 \end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

### An Aside
Penalized quasi-likelihood and marginal quasi-likelihood are quite similar, but they do have notable differences. For one, MQL will not result in predictions for the random effects. Thus, it is not suitable for inference at the group level. In addition, the estimates of the fixed effects are essentially for the marginal model. Fitzmaurice et al.<d-cite key=fitzmaurice2011></d-cite> note that MQL estimates are highly biased unless the random effects variance is near zero, regardless of the number of repeated measurements.

In contrast, though PQL can be used to estimate the fixed effects and predict the random effects, there are some complications. For small counts or low numbers of repeated measurements, PQL underestimates the fixed effects as well as the variance components. PQL and MQL can also be reframed as approximations of a GLMM with an adjusted LMM (see below). 

We can reframe penalized and marginal quasi-likelihood as arising from a Gaussian approximation of the GLMM at hand. We'll follow Chapter 15 in Fitzmaurice et al.<d-cite key=fitzmaurice2011></d-cite>  We assume to have the following generalized linear mixed model:

$$
\begin{equation}
\label{eq:glmm-y}
\begin{aligned}
\mathbf{y}_{i,j} \rvert \beta_i &\sim \mathcal{F}(\mu_{i,j}); \\
\mu_{i,j} &= g^{-1}\left(\eta_{i,j}\right) = g^{-1}\left(\alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j}\right)
\end{aligned}
\end{equation}
$$

We first approximate the response by assuming it can be represented as its conditional mean plus some Gaussian error:

$$
\mathbf{y}_{i,j} \approx \mu_{i,j} + \epsilon_{i,j}
$$

where, for some variance function $V(\cdot)$ of the conditional mean (e.g. $V(\mu_{i,j}) = \mu_{i,j}$ in the Poisson case), the errors satisfy:

$$
\mathbb{E}\left[ \epsilon_{i,j} \right] = 0;
\hspace{5mm}
\text{Var}(\epsilon_{i,j} ) = V(\mu_{i,j});
\hspace{5mm}
\text{Cov}(\epsilon_{i,j}, \epsilon_{i', j'} ) = 0
$$

<div class="theorem">
<strong>Claim (Approximating the Response).</strong>
{% tabs claim-6 %}
{% tab claim-6 statement %}
Since the link function is non-linear, we will linearize the mean with a first-order Taylor approximation about $\hat{\eta}_{i,j}$, the linear predictor estimated under $H_0$:

$$
\mathbf{y}_{i,j} \approx g^{-1}(\hat{\eta}_{i,j}) + \delta(\hat{\eta}_{i,j}) (\eta_{i,j} - \hat{\eta}_{i,j}) + \epsilon_{i,j} 
$$

where $$\delta(\hat{\eta}_{i,j}) = \frac{\partial g^{-1}(\eta_{i,j})}{\partial \eta_{i,j}}\bigg\rvert_{\eta_{i,j} = \hat{\eta}_{i,j}}$$, the gradient of the inverse link function w.r.t $\eta_{i,j}$ evaluated at the estimate under $H_0$.
{% endtab %}
{% tab claim-6 proof %}
$$
\begin{aligned}
\mathbf{y}_{i,j} &\approx \left[g^{-1}(\eta_{i,j})\right]\bigg\rvert_{\eta_{i,j} = \hat{\eta}_{i,j}} + \frac{\partial g^{-1}(\eta_{i,j})}{\partial \eta_{i,j}}\bigg\rvert_{\eta_{i,j} = \hat{\eta}^{(c)}_{i,j}} (\eta_{i,j} - \hat{\eta}_{i,j}) + \epsilon_{i,j} \\
&= g^{-1}(\hat{\eta}_{i,j}) + \delta(\hat{\eta}_{i,j}) (\eta_{i,j} - \hat{\eta}_{i,j}) + \epsilon_{i,j} & \left(\delta(\hat{\eta}_{i,j}) = \frac{\partial g^{-1}(\eta_{i,j})}{\partial \eta_{i,j}}\bigg\rvert_{\eta_{i,j} = \hat{\eta}^{(c)}_{i,j}}\right)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

<div class="theorem">
<strong>Claim (Working Response).</strong>
{% tabs claim-7 %}
{% tab claim-7 claim %}
We can rearrange the above as:

$$
\mathbf{y}_{i,j}^\star \approx \alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j} + \epsilon^\star_{i,j}
$$

where $$\mathbf{y}_{i,j}^\star = \frac{1}{\delta(\hat{\eta}_{i,j})}\left(\mathbf{y}_{i,j}  - \hat{\mu}_{i,j}\right) + \hat{\eta}_{i,j}$$ and $$\epsilon^\star_{i,j} = \frac{1}{\delta(\hat{\eta}_{i,j})}\epsilon_{i,j}$$.

{% endtab %}
{% tab claim-7 proof %}
$$
\begin{aligned}
&\mathbf{y}_{i,j} \approx g^{-1}(\hat{\eta}_{i,j}) + \delta(\hat{\eta}_{i,j}) (\eta_{i,j} - \hat{\eta}_{i,j}) + \epsilon_{i,j} \\
\implies &\mathbf{y}_{i,j} \approx \hat{\mu}_{i,j} + \delta(\hat{\eta}_{i,j})(\alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j}) -  \delta(\hat{\eta}_{i,j})(\hat{\alpha}^\top \mathbf{x}_{i,j} + \mathbf{0}^\top \mathbf{z}_{i,j}) + \epsilon_{i,j} \\
\implies &\mathbf{y}_{i,j} + \delta(\hat{\eta}_{i,j})(\hat{\alpha}^\top \mathbf{x}_{i,j}) \approx \hat{\mu}_{i,j} + \delta(\hat{\eta}_{i,j})(\alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j}) + \epsilon_{i,j} \\
\implies &\mathbf{y}_{i,j}  - \hat{\mu}_{i,j} + \delta(\hat{\eta}_{i,j})(\hat{\alpha}^\top \mathbf{x}_{i,j}) \approx \delta(\hat{\eta}_{i,j})(\alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j}) + \epsilon_{i,j} \\
\implies &\frac{1}{\delta(\hat{\eta}_{i,j})}\left(\mathbf{y}_{i,j}  - \hat{\mu}_{i,j}\right) + \hat{\alpha}^\top \mathbf{x}_{i,j} \approx \alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j} + \frac{1}{\delta(\hat{\eta}_{i,j})}\epsilon_{i,j} \\
\implies &\frac{1}{\delta(\hat{\eta}_{i,j})}\left(\mathbf{y}_{i,j}  - \hat{\mu}_{i,j}\right) + \hat{\eta}_{i,j} \approx \alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j} + \frac{1}{\delta(\hat{\eta}_{i,j})}\epsilon_{i,j}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Notice that since $$\epsilon^*_{i,j}$$ is just a scaled version of $\epsilon_{i,j}$, it is also Gaussian with:

$$
\begin{aligned}
\mathbb{E}\left[ \epsilon^\star_{i,j} \right] &= \mathbb{E}\left[ \frac{1}{\delta(\hat{\eta}_{i,j})} \epsilon_{i,j} \right] = 0; \\
\text{Var}\left(\epsilon^\star_{i,j} \right) &= \text{Var}\left( \frac{1}{\delta(\hat{\eta}_{i,j})} \epsilon_{i,j}  \right) = \frac{1}{\delta^2(\hat{\eta}_{i,j})} \text{Var}(\epsilon_{i,j}) = \frac{V(\hat{\mu}_{i,j})}{\delta^2(\hat{\eta}_{i,j})}
\end{aligned}
$$

The above essentially specifies a linear mixed model:

$$
\mathbf{y}^\star_{i,j} \sim \mathcal{N}\left(\alpha^\top \mathbf{x}_{i,j} + \beta_i^\top \mathbf{z}_{i,j}, \frac{V(\hat{\mu}_{i,j})}{\delta^2(\hat{\eta}_{i,j})}\right)
$$

To use this in practice, we can first estimate the parameters with iteratively reweighted least squares (or something like that), then use the estimates to compute the working response and errors. Then we can proceed how we would with a linear mixed model. 

