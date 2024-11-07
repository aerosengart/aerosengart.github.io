---
layout: post
title:  "Deriving A Variance Component Score Test - Poisson"
date: 06 November 2024
categories: posts
use_math: true
include_scripts: [
    "/assets/js/snackbar.js",
    "/assets/js/popup.js",
]
---


## Set-Up

The set-up is the standard generalized linear mixed model setting. Suppose we have $n$ samples and $k$ targets (or phenotypes, if you're dealing with GWAS). 

Our data will include a count table, $Y$, that is $n \times k$, some fixed covariates, $X$, that is $n \times p$, and some random covariates, $Z$, that is $n \times q$. 

We'll assume that the counts are Poisson, so our link function is $\log(\cdot)$. Let $\mu_{j, i}$ denote the mean of $Y_{j, i}$ (the count for target $j$ of sample $i$). We further assume:

$$
\log(\mu_{j,i}) = X \alpha + Z\beta
$$

for $p \times k$ fixed effects matrix $\alpha$ and $q \times k$ random effects matrix $\beta$. Conditional on the random effects, we assume the counts are independent and identically distributed.

Let $\text{vec}(\cdot)$ denote vectorization of the input where we take the matrix of interest and concanate its columns into one long vector. As is standard with mixed models, we assume:

$$
\text{vec}(\beta) \sim \mathcal{N}(\vec{0}, \tau^2 \left( \Sigma_Z \otimes \Sigma_T \right))
$$

where $\tau^2$ is a variance parameter, $\Sigma_Z$ is an $q \times q$ covariance matrix associated with the random effects, and $\Sigma_T$ is a $k \times k$ covariance matrix associated with the targets.

Our interest lies in testing whether the columns of $\beta$ are $\vec{0}$ (i.e. $\beta_1 = \beta_2 = \dots = \beta_k = \vec{0}$). Equivalently, we can just test whether $\tau = 0$. If $\tau = 0$, then $\text{vec}(\beta)$ has $0$ variance, and all of the random effects are $0$. Thus, our null hypothesis will be $H_0: \tau = 0$. 

---

## Score Test

We will need three ingredients to write out the score test:

1. The log likelihood function
2. The score, $U(\tau)$, which is the gradient of the log-likelihood with respect to $\tau$
3. The Fisher Information, $\mathcal{I}(\tau)$, which is the variance of the score

Once we have the above pieces, the test statistic for a null hypothesis $H_0: \tau = \tau_0$ is given by:

$$
S(\tau_0) = U^\top(\hat{\tau}_0) \mathcal{I}^{-1}(\hat{\tau}_0) U(\hat{\tau}_0)
$$

where $\hat{\tau}_0$ is the maximum likelihood estimate of $\tau$ under $H_0$.

---

## Case 1

For simplicity, let's assume that $\Sigma_T = \mathbb{I}_{k \times k}$, the $k$-dimensional identity matrix and $q = 1$. This implies that the random effects are independent and that $\beta$ and $Z$ are $k$-dimensional vectors. Let's assume they are also transposed to they are column vectors to make notation more consistent.

#### Likelihood

Let $y_i$ ($1 \times k$), $x_i$ ($1 \times p$), and $z_i$ ($1 \times q$) denote the $i$-th rows of $Y$, $X$, and $Z$, respectively. Let $\alpha_j$ ($p \times 1$) and $\beta_j$ ($q \times 1$) denote the $j$-th columns of $\alpha$ and $\beta$, respectively. Recall that since we assumed $q = 1$, $\beta_j$ and $z_i$ are scalar-valued.

The likelihood is given by:

$$
\begin{aligned}
\mathcal{L}(\alpha, \tau; Y, X, Z) &= \int \mathcal{L}(\alpha, \tau, \beta; Y, X, Z) d\beta \\
&= \int \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \mathcal{L}(\beta) d\beta \\
&= \int \exp \left(  \ell(\alpha, \tau; Y, X, Z\rvert \beta\right) \left( \frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert }} \exp \left(- \frac{1}{2}\beta^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta \right) \right) d\beta\\
&= \frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert }} \int \exp \left( \ell(\alpha, \tau; Y, X, Z \rvert \beta) - \frac{1}{2} \beta^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta \right) d\beta\\
\end{aligned}
\label{eq:lik}
$$

Since Eq. \eqref{eq:lik} is of the form $c \rvert \tau^2 \Sigma_T \rvert^{-1/2} \int \exp(f(\beta)) d\beta$ for some constant $c$ and function $f$, we can use Laplace's method to approximate the integral part. 

We first expand $f(\beta) = \ell(\alpha, \tau; Y, X, Z \rvert \beta) - \frac{1}{2}\beta^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta$ about $\tilde{\beta}$, which requires the first- and second-order partial derivatives of $f$ with respect to $\beta$:

$$
\begin{aligned}
\frac{\partial}{\partial \beta} f(\beta) &= \begin{bmatrix}
\sum_{i = 1}^n [y_{i, 1}z_i - \exp(x_i \alpha_1 + z_i \beta_1) z_i] \\
\vdots \\
\sum_{i = 1}^n [y_{i, k}z_i - \exp(x_i \alpha_k + z_i \beta_k) z_i]
\end{bmatrix} 
- (\tau^2 \Sigma_T)^{-1}\beta \\
\frac{\partial^2}{\partial \beta \partial \beta^\top} f(\beta) &= \begin{bmatrix}
- \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \beta_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \sum_{i = 1}^n \exp(x_i \alpha_k + z_i \beta_k)z_i^2
\end{bmatrix} 
- (\tau^2 \Sigma_T)^{-1}
\end{aligned}
$$

---
<details>
<summary>Derivation Of First-Order Partial Derivatives.</summary>

The conditional log-likelihood for the complete data is the sum of the conditional log-likelihoods for each observation. Furthermore, since the elements of $\beta$ are also independent, we can just sum over the targets as well:

$$
\begin{aligned}
\ell(\alpha, \tau; Y, X, Z \rvert \beta) &=  \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j + z_i \beta_j \right) - \exp\left( x_i \alpha_j + z_i \beta_j \right) - \log(y_{i,j}!) \right]
\end{aligned}
$$

The $m$-th element of the gradient vector with respect to $\beta$ is given by:

$$
\begin{aligned}
\frac{\partial}{\partial \beta_m} \ell(\alpha, \tau; Y, X, Z \rvert \beta)&= \frac{\partial}{\partial \beta_m} \left[ \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta)\right]  \\
&=  \sum_{i = 1}^n \sum_{j = 1}^k \frac{\partial}{\partial \beta_m} \left[ y_{i,j} \left( x_i \alpha_j + z_i \beta_j \right) - \exp\left( x_i \alpha_j + z_i \beta_j \right) - \log(y_{i,j}!) \right]  \\
&= \sum_{i = 1}^n [y_{i, m}z_i - \exp(x_i \alpha_m + z_i \beta_m) z_i]
\end{aligned}
$$

The gradient of $\frac{1}{2}\beta^\top (\tau^2 \Sigma_T)^{-1} \beta$ is:

$$
\frac{\partial}{\partial \beta} \left[ \frac{1}{2}\beta^\top (\tau^2 \Sigma_T)^{-1} \beta \right] = (\tau^2 \Sigma_T)^{-1} \beta 
$$

Thus, the $k \times 1$ first-order partial derivative vector is:

$$
\begin{bmatrix}
\sum_{i = 1}^n [y_{i, 1}z_i - \exp(x_i \alpha_1 + z_i \beta_1) z_i] \\
\vdots \\
\sum_{i = 1}^n [y_{i, k}z_i - \exp(x_i \alpha_k + z_i \beta_k) z_i]
\end{bmatrix} 
- (\tau^2 \Sigma_T)^{-1}\beta
$$
</details>
---
<details>
<summary>Derivation Of Second-Order Partial Derivatives.</summary>
To calculate the second-order partial derivatives, we just find the partial derivatives of the first-order partial derivative vector. 

We have:

$$
\begin{aligned}
\frac{\partial^2}{\partial \beta_m^2} \ell(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial}{\partial \beta_m} \left[ \sum_{i = 1}^n [y_{i, 1}z_i - \exp(x_i \alpha_m + z_i \beta_m) z_i]
\right] \\
&= - \sum_{i = 1}^n \exp(x_i \alpha_m + z_i \beta_m)z_i^2 \\
\frac{\partial^2}{\partial \beta_l \partial \beta_m} \ell(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial}{\partial \beta_l} \left[ \sum_{i = 1}^n [y_{i, 1}z_i - \exp(x_i \alpha_m + z_i \beta_m) z_i]
\right] \\
&= 0
\end{aligned}
\label{eq:second-a}
$$

Furthermore, we have that:
$$
\frac{\partial^2}{\partial \beta \partial \beta^\top} \beta^\top (\tau^2\Sigma_T)^{-1} \beta = (\tau^2 \Sigma_T)^{-1}
\label{eq:second-b}
$$

Putting Eqs. \eqref{eq:second-a} and \eqref{eq:second-b} yields the matrix of partial derivatives:

$$
\begin{bmatrix}
- \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \beta_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \sum_{i = 1}^n \exp(x_i \alpha_k + z_i \beta_k)z_i^2
\end{bmatrix} 
- (\tau^2 \Sigma_T)^{-1}
$$
</details>
---

We choose $\tilde{\beta}$ to be the solution to $\frac{\partial}{\partial \beta} f(\beta) = 0$. Our expansion is then:

$$
\begin{aligned}
f(\beta) &= f(\tilde{\beta}) + (\beta - \tilde{\beta})^\top \frac{\partial}{\partial \beta} f(\tilde{\beta}) + \frac{1}{2}(\beta - \tilde{\beta})^\top \frac{\partial^2}{\partial \beta \partial \beta^\top} f(\tilde{\beta}) (\beta - \tilde{\beta}) + R \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] \\
&\hspace{6mm} -\frac{1}{2}
\sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right]
- \frac{1}{2} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right]
- \frac{1}{2}\beta^\top  (\tau^2 \Sigma_T)^{-1} \beta + R
\end{aligned}
$$

where $R$ contains the higher order terms of the expansion.

---
<details>
<summary>Taylor Expansion Simplification.</summary>
$$
\begin{aligned}
f(\beta) &= f(\tilde{\beta}) + (\beta - \tilde{\beta})^\top \frac{\partial}{\partial \beta} f(\tilde{\beta}) + \frac{1}{2}(\beta - \tilde{\beta})^\top \frac{\partial^2}{\partial \beta \partial \beta^\top} f(\tilde{\beta}) (\beta - \tilde{\beta}) + R \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] - \frac{1}{2}\tilde{\beta}^\top \left( \tau^2 \Sigma_T \right)^{-1} \tilde{\beta} \\
&\hspace{6mm} + \frac{1}{2}(\beta - \tilde{\beta})^\top 
\left(
\begin{bmatrix}
- \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \tilde{\beta}_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \sum_{i = 1}^n \exp(x_i \alpha_k + z_i \tilde{\beta}_k)z_i^2
\end{bmatrix} 
- (\tau^2 \Sigma_T)^{-1}\right)
(\beta - \tilde{\beta}) + R \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] - \frac{1}{2}\tilde{\beta}^\top \left( \tau^2 \Sigma_T \right)^{-1} \tilde{\beta} \\
&\hspace{6mm} + 
\frac{1}{2}(\beta - \tilde{\beta})^\top  \begin{bmatrix}
- \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \tilde{\beta}_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \sum_{i = 1}^n \exp(x_i \alpha_k + z_i \tilde{\beta}_k)z_i^2
\end{bmatrix} (\beta - \tilde{\beta})
- \frac{1}{2}(\beta - \tilde{\beta})^\top  (\tau^2 \Sigma_T)^{-1} (\beta - \tilde{\beta}) 
 + R \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] \\
&\hspace{6mm} + 
\frac{1}{2}\beta^\top  \begin{bmatrix}
- \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \tilde{\beta}_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \sum_{i = 1}^n \exp(x_i \alpha_k + z_i \tilde{\beta}_k)z_i^2
\end{bmatrix} \beta \\
&\hspace{6mm} + 
\frac{1}{2} (-\tilde{\beta})^\top  \begin{bmatrix}
- \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \tilde{\beta}_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \sum_{i = 1}^n \exp(x_i \alpha_k + z_i \tilde{\beta}_k)z_i^2
\end{bmatrix} (-\tilde{\beta})
- \frac{1}{2}\beta^\top  (\tau^2 \Sigma_T)^{-1} \beta 
+ R \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] \\
&\hspace{6mm} + 
\begin{bmatrix}
- \frac{1}{2}\beta_1 \sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \tilde{\beta}_1)z_i^2 \\
0 & \ddots & 0 \\
0 & \dots & - \frac{1}{2}\beta_k\sum_{i = 1}^n \exp(x_i \alpha_k + z_i \tilde{\beta}_k)z_i^2
\end{bmatrix} \beta \\
&\hspace{6mm} + 
\begin{bmatrix}
\frac{1}{2} \tilde{\beta}_1\sum_{i = 1}^n \exp(x_i \alpha_1 + z_i \tilde{\beta}_1)z_i^2 & \dots & 0 \\
0 & \ddots & 0 \\
0 & \dots & - \frac{1}{2} \tilde{\beta}_k\sum_{i = 1}^n \exp(x_i \alpha_k + z_i \tilde{\beta}_k)z_i^2
\end{bmatrix} (-\tilde{\beta})
- \frac{1}{2}\beta^\top  (\tau^2 \Sigma_T)^{-1} \beta + R \\
&= \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] \\
&\hspace{6mm} -\frac{1}{2}
\sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right]
- \frac{1}{2} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right]
- \frac{1}{2}\beta^\top  (\tau^2 \Sigma_T)^{-1} \beta + R
\end{aligned}
$$
</details>
---

Supposing that $R$ is negligible, we can <button class="bland-button" onclick="SnackbarFunc()">approximate the integral</button> in Eq. \eqref{eq:lik} as:

<div id="snackbar"> 😮‍💨 Phew, that was a lot of work for just the first step! 😮‍💨</div>

$$
\begin{aligned}
\mathcal{L}(\alpha, \tau; Y, X, Z) &= \frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert }} \int \exp \left( \ell(\alpha, \tau; Y, X, Z \rvert \beta) - \frac{1}{2} \beta^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta \right) d\beta \\
&\approx 
\frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}} 
\sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] \\
&\hspace{6mm} - \frac{1}{2\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}} 
\sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right] \\
&\hspace{6mm} -\frac{1}{2\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right] \\
&\hspace{6mm} -\frac{1}{2\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}}\beta^\top  (\tau^2 \Sigma_T)^{-1} \beta 
\end{aligned}
$$


#### Score

Now that we have the (approximate) likelihood $\mathcal{L}(\alpha, \tau; Y, X, Z)$, we can take the first derivative with respect to $\tau$ to get the score. Let $A(\tau) = \sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}$.

$$
\begin{aligned}
U(\tau) &= \frac{d}{d \tau} \left[ \mathcal{L}(\alpha, \tau; Y, X, Z) \right] \\
&\approx \frac{d}{d \tau} \left[ 
\frac{1}{A(\tau)} 
\sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right] \right] \\
&\hspace{6mm} - \frac{d}{d \tau} \left[\frac{1}{A(\tau)}  
\sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right] \right] \\
&\hspace{6mm} - \frac{d}{d \tau} \left[ \frac{1}{A(\tau)} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right] \right]\\
&\hspace{6mm} - \frac{d}{d \tau} \left[ \frac{1}{A(\tau)} \beta^\top  (\tau^2 \Sigma_T)^{-1} \beta \right] \\

&= -\frac{k}{\tau}
\sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right]  \\
&\hspace{6mm} + \frac{k}{\tau} \sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right]  \\
&\hspace{6mm} + \frac{k}{\tau} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right] \\
&\hspace{6mm} - \frac{d}{d \tau} \left[ \frac{1}{A(\tau)} \right] \beta^\top  (\tau^2 \Sigma_T)^{-1} \beta - \frac{1}{A(\tau)} \frac{d}{d \tau} \left[ \beta^\top  (\tau^2 \Sigma_T)^{-1} \beta \right] \\

&= -\frac{k}{\tau}
\sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right]  \\
&\hspace{6mm} + \frac{k}{\tau} \sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right]  \\
&\hspace{6mm} + \frac{k}{\tau} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right] \\
&\hspace{6mm} + \frac{k}{\tau} \beta^\top  (\tau^2 \Sigma_T)^{-1} \beta - \frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}} \times \frac{1}{\tau^3} \beta^\top (\Sigma_T)^{-1} \beta \\ 

&\overset{(i)}{=} \frac{k}{\tau}
\sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!)\right]  \\
&\hspace{6mm} + \frac{k}{\tau} \sum_{j = 1}^k
\beta_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2\right]  \\
&\hspace{6mm} + \frac{k}{\tau} \sum_{j = 1}^k \tilde{\beta}_j^2 \left[ \sum_{i = 1}^n \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2 \right] \\
&\hspace{6mm} + \frac{k}{\tau^3} \beta^\top  (\Sigma_T)^{-1} \beta - \frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}} \times \frac{1}{\tau^3} \beta^\top (\Sigma_T)^{-1} \beta \\ 

&\overset{(i)}{=} \frac{k}{\tau}
\sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j}(x_i \alpha_j + z_i \tilde{\beta_j}) - \exp(x_i \alpha_j + z_i \tilde{\beta}_j) - \log(y_{i,j}!) +  (\beta_j^2 + \tilde{\beta}_j^2) \exp(x_i \alpha_j + z_i \tilde{\beta}_j)z_i^2  \right] + \frac{k}{\tau^3} \beta^\top  (\Sigma_T)^{-1} \beta \left( 1 - \frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}} \right) \\ 
\end{aligned}
$$

where $(i)$ follows from the fact that for a constant $c$ and invertible, $k \times k$ matrix $M$:

$$
(cM)^{-1} = (c \mathbb{I}_{k \times k} M)^{-1} = M^{-1}(c \mathbb{I}_{k \times k})^{-1} = \frac{1}{c}M^{-1}
$$

---
<details>
<summary>Derivation of $\frac{d}{d\tau} (A(\tau))^{-1}$.</summary>
Recall the following fact for matrices. Let $c$ be some constant, and let $M$ be a $k \times k$ matrix.

$$
\begin{aligned}
\rvert cM \rvert &= \rvert c \mathbb{I}_{k \times k} M \rvert \\
&= \rvert c \mathbb{I}_{k \times k} \rvert \rvert M \rvert \\
&= c^k \rvert M \rvert
\end{aligned}
$$

Thus, the derivative of $A(\tau)$ is given by:

$$
\begin{aligned}
\frac{d}{d\tau} A(\tau) &= \frac{d}{d\tau} \sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert} \\
&= \frac{d}{d\tau} \sqrt{(2\pi)^k \tau^{2k}\rvert  \Sigma_T \rvert} \\
&= \frac{1}{2} \left( (2\pi)^k \tau^{2k}\rvert \Sigma_T \rvert) \right)^{-1/2} \frac{d}{d\tau} \left[ (2\pi)^k \tau^{2k}\rvert  \Sigma_T \rvert \right] \\
&= \frac{1}{2} \left( (2\pi)^k \tau^{2k}\rvert \Sigma_T \rvert) \right)^{-1/2}\left( 2k (2\pi)^k \tau^{2k - 1} \rvert \Sigma_T \rvert\right) \\
&=  ((2\pi)^{-k/2} \tau^{-k}\rvert \Sigma_T \rvert^{-1/2})\left(k (2\pi)^k \tau^{2k - 1} \rvert \Sigma_T \rvert\right) \\ 
&= k (2 \pi)^{k/2}\tau^{k - 1}\rvert \Sigma_T \rvert^{1/2}
\end{aligned}
$$

Using standard derivative rules yields:

$$
\begin{aligned}
\frac{d}{d \tau} (A(\tau))^{-1} &= -\frac{1}{A(\tau)^2}\frac{d}{d \tau} \left[ A(\tau) \right] \\
&= -\frac{1}{\sqrt{(2\pi)^k \rvert \tau^2 \Sigma_T \rvert}}\left( k (2 \pi)^{k/2}\tau^{k - 1}\rvert \Sigma_T \rvert^{1/2} \right) \\
&= -\frac{k (2 \pi)^{k/2}\tau^{k - 1}\rvert \Sigma_T \rvert^{1/2}}{(2\pi)^{k/2} \tau^{k} \rvert  \Sigma_T \rvert^{1/2}} \\
&= -\frac{k}{\tau}
\end{aligned}
$$
</details>
---
<details>
<summary>Derivation of $\frac{d}{d\tau} [\beta^\top (\tau^2 \Sigma_T)^{-1} \beta]$.</summary>
Note that $\tau^2 \Sigma_T = \tau^2 \mathbb{I}_{k \times k} \Sigma_T$. Assuming that $\Sigma_T$ is invertible, we have that:

$$
(\tau^2 \mathbb{I}_{k \times k} \Sigma_T)^{-1} = (\Sigma_T)^{-1} (\tau^2 \mathbb{I}_{k \times k})^{-1} = (\Sigma_T)^{-1}\left(\frac{1}{\tau^2} \mathbb{I}_{k \times k}\right)
$$

The derivative is then:

$$
\begin{aligned}
\frac{d}{d\tau} [\beta^\top (\tau^2 \Sigma_T)^{-1} \beta] &= \beta^\top \frac{d}{d\tau} \left[ (\Sigma_T)^{-1}\left(\frac{1}{\tau^2} \mathbb{I}_{k \times k}\right) \right]  \beta \\
&= \beta^\top (\Sigma_T)^{-1} \left(-2 \frac{1}{\tau^3}\right) \mathbb{I}_{k \times k} \beta \\
&= \frac{1}{\tau^3} \beta^\top (\Sigma_T)^{-1} \beta
\end{aligned}
$$
</details>


---

## References

Lin, X. "Variance component testing in generalised linear models with random effects". Biometrika, Volume 84, Issue 2, June 1997. Pages 309–326. [doi:10.1093/biomet/84.2.309](https://doi.org/10.1093/biomet/84.2.309).






<!-- q > 1 -->
<!-- #### Log-Likelihood

The log-likelihood is given by:

$$
\begin{aligned}
\mathcal{L}(\alpha, \tau; Y, X, Z) &= \int \cdots \int \mathcal{L}(\alpha, \tau, \beta; Y, X, Z) d\beta_1 \dots d\beta_k \\
&= \int \cdots \int \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \mathcal{L}(\beta) d\beta_1 \dots d\beta_k \\
&= \int \cdots \int \exp \left(  \ell(\alpha, \tau; Y, X, Z\rvert \beta\right) \prod_{j = 1}^k \left( \frac{1}{\sqrt{(2\pi)^q \rvert \tau^2 \Sigma_T \rvert }} \exp \left(- \frac{1}{2}\beta_j^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta_j \right) \right) d\beta_1 \dots d\beta_k\\
&= \frac{1}{\sqrt{(2\pi)^q \rvert \tau^2 \Sigma_T \rvert }} \int \cdots \int \exp \left( \ell(\alpha, \tau; Y, X, Z \rvert \beta) - \frac{1}{2} \sum_{j = 1}^k \beta_j^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta_j \right) d\beta_1 \dots d\beta_k \\
\end{aligned}
\label{eq:lik}
$$

Since Eq. \eqref{eq:lik} is of the form $c \rvert \tau^2 \Sigma_T \rvert^{-1/2} \int \exp(f(\beta)) d\beta$ for some constant $c$ and function $f$, we can use Laplace's method to approximate the integral part. 

We first expand $f(\beta) = \ell(\alpha, \tau; Y, X, Z \rvert \beta) - \frac{1}{2}\beta^\top \left( \tau^2 \Sigma_T \right)^{-1} \beta$ about $\tilde{\beta}$, which requires the first- and second-order partial derivatives of $f$ with respect to $\beta$:

$$
\begin{aligned}
\frac{\partial}{\partial \beta} f(\beta) &= ?? \\
\frac{\partial^2}{\partial \beta \partial \beta^\top} f(\beta) &= ??
\end{aligned}
$$ 


---
<details>
<summary>First-Order Partial Derivatives.</summary>

Let $y_i$ ($1 \times k$), $x_i^\top$ ($1 \times p$), and $z_i^\top$ ($1 \times q$) denote the $i$-th rows of $Y$, $X^\top$, and $Z^\top$, respectively. Let $\alpha_j$ ($p \times 1$) and $\beta_j$ ($q \times 1$) denote the $j$-th columns of $\alpha$ and $\beta$, respectively. 

Since we only have a single fixed effect and a single random effect, $x_i$, $z_i$, $\alpha_j$, and $\beta_j$ are all scalars, so we can drop the transpositions. Let $\mathbf{1}_k$ denote a $k$-dimensional vector of ones, and let $\exp(v)$ for vector $v$ be the element-wise exponentiation. Note that:

$$
\begin{aligned}
\ell(\alpha, \tau; Y, X, Z \rvert \beta) &=  \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i^\top, z_i^\top \rvert \beta) \\
&=  \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i^\top \alpha_j + z_i^\top \beta_j \right) - \exp\left( x_i^\top \alpha_j + z_i^\top \beta_j \right) - \log(y_{i,j}!) \right]
\end{aligned}
$$

The $q$-dimensional gradient vector with respect to $\beta_m$ is given by:

$$
\begin{aligned}
\frac{\partial}{\partial \beta_m} \ell(\alpha, \tau; Y, X, Z \rvert \beta)&= \frac{\partial}{\partial \beta_m} \left[ \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta)\right]  \\
&= \frac{\partial}{\partial \beta_m} \left[ \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i^\top \alpha_j + z_i^\top \beta_j \right) - \exp\left( x_i^\top \alpha_j + z_i^\top \beta_j \right) - \log(y_{i,j}!) \right] \right] \\
&= \sum_{i = 1}^n \left[ y_{i,m} z_i^\top - \exp(x_i^\top \alpha_m + z_i^\top \beta_m)z_i^\top \right] \\
&= \sum_{i = 1}^n \left(y_{i,m} - \exp(x_i^\top \alpha_m + z_i^\top \beta_m) \right)z_i^\top
\end{aligned}
$$

The gradient of $\frac{1}{2}\beta^\top (\tau^2 \Sigma_T)^{-1} \beta$ is $(\tau^2 \Sigma_T)^{-1} \beta_m$, so the first-order partial derivative vector is:

$$
\begin{bmatrix}
\left( \sum_{i = 1}^n \left(y_{i,1} - \exp(x_i^\top \alpha_1 + z_i^\top \beta_1) \right)z_i^\top \right)^\top & \dots &
\left( \sum_{i = 1}^n \left(y_{i,k} - \exp(x_i^\top \alpha_k + z_i^\top \beta_k) \right)z_i^\top \right)^\top
\end{bmatrix} - (\tau^2 \Sigma_T)^{-1}\beta
$$
</details>
---


<details>
<summary>Second-Order Partial Derivatives.</summary>
The second-order partial derivatives will be a $k \times k$ matrix. First we will do the conditional log-likelihood to get its diagonal components:

$$
\begin{aligned}
\frac{\partial^2}{\partial \beta_m^2} \ell(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial}{\partial \beta_m} \left[ \sum_{i = 1}^n [y_{i, m}z_i - \exp(x_i \alpha_m + z_i \beta_m)z_i]\right] \\
&= -\sum_{i = 1}^m \exp(x_i \alpha_m + z_i \beta_m)z_i^2
\end{aligned}
$$


And the off-diagonal elements will be given by:

$$
\begin{aligned}
\frac{\partial^2}{\partial \beta_m \partial \beta_l} \ell(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial}{\partial \beta_l} \left[ \sum_{i = 1}^n [y_{i, m}z_i - \exp(x_i \alpha_m + z_i \beta_m)z_i] \right] \\
&= - 
\end{aligned}
$$

</details>
---


#### Score

Taking the derivative (or gradient in multi-dimensional settings) of Eq. \eqref{eq:lik} with respect to $\tau$ can be difficult because the integral is cumbersome. As in (Lin, 1997), we will do a Taylor expansion of $\mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta)$ about $\beta = \vec{0}$. 

---
<details>
  <summary>Taylor Expansion Details.</summary>
  A Taylor expansion about $\beta = \vec{0}$ is given by:

  $$
  \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) = \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} + \beta^\top \left[ \frac{\partial}{\partial \beta} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \vec{0}) \bigg \rvert_{\beta = \vec{0}} \right]  + \frac{1}{2} \beta^\top \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \beta + e
  $$

  where $e$ contains all higher order partial derivative terms. To write the above equation, we need to find the first- and second-order partial derivatives of $\mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta)$ with respect to $\beta$. Let $y_i$, $x_i^\top$, and $z_i^\top$ denote the $i$-th rows of $Y$, $X$, and $Z$, respectively. Let $\alpha_j$ and $\beta_j$ denote the $j$-th columns of $\alpha$ and $\beta$, respectively. Since we only have a single fixed effect and a single random effect, $x_i$, $z_i$, $\alpha_j$, and $\beta_j$ are all scalars, so we can drop the transpositions. We can write our likelihood as:

  $$
  \begin{aligned}
  \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \\
  &= \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j + z_i \beta_j \right) - \exp\left( x_i \alpha_j + z_i \beta_j \right) - \log(y_{i,j}!) \right] \right)
  \end{aligned}
  \label{eq:split-lik}
  $$

  Evaluating Eq. \eqref{eq:split-lik} at $\beta = \vec{0}$ yields:

  $$
  \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} = \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right)
  $$

  The elements of the gradient are given by:

  $$
  \begin{aligned}
  \frac{\partial}{\partial \beta_j} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial}{\partial \beta_j} \left[ \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \right] \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \frac{\partial}{\partial \beta_j} \left[ \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right] \right) \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n \frac{\partial}{\partial \beta_j} \left[ \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j + z_i \beta_j \right) - \exp\left( x_i \alpha_j + z_i \beta_j \right) - \log(y_{i,j}!) \right]  \right] \right) \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n \frac{\partial}{\partial \beta_j} \left[ y_{i,j} \left( x_i \alpha_j + z_i \beta_j \right) - \exp\left( x_i \alpha_j + z_i \beta_j \right) - \log(y_{i,j}!) \right] \right) \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right) \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n (y_{i,j} - \mu_{i,j}) z_i\right) \\
  \implies  \frac{\partial}{\partial \beta_j} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} &= \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right) \left( \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp(x_i \alpha_j) z_i \right]\right)
  \end{aligned}
  \label{eq:first-order}
  $$

  The diagonal elements of the Hessian are given by:

  $$
  \begin{aligned}
  \frac{\partial^2}{\partial \beta_j^2} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial^2}{\partial \beta_j^2} \left[ \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \right] \\
  &=  \frac{\partial}{\partial \beta_j} \left[ \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right) \right] \\
  &=  \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \frac{\partial}{\partial \beta_j} \left[ \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right] + \frac{\partial}{\partial \beta_j} \left[ \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \right] \left( \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right)  \\
  &=  \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right)\left( \sum_{i = 1}^n - \exp(x_i \alpha_j + z_i \beta_j) z_i^2 \right) + \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left(\sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right)^2 \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left[ \left(\sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right)^2 -\sum_{i = 1}^n \exp(x_i \alpha_j + z_i \beta_j) z_i^2 \right]  \\
  \implies  \frac{\partial^2}{\partial \beta_j^2} \mathcal{L}(\alpha, \tau; Y, X, Z, \rvert \beta) \bigg \rvert_{\beta = \vec{0}} &= \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right) \left[ \left( y_{i,j} z_i - \exp(x_i \alpha_j) z_i \right)^2  - \sum_{i = 1}^n \exp(x_i \alpha_j)z_i^2 \right] 
  \end{aligned}
  \label{eq:second-order-diag}
  $$

  The off-diagonal elements of the Hessian are:

  $$
  \begin{aligned}
  \frac{\partial^2}{\partial \beta_j \partial \beta_m} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) &= \frac{\partial^2}{\partial \beta_j \partial \beta_m} \left[ \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \right] \\
  &= \frac{\partial}{\partial \beta_m} \left[ \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right) \right] \\
  &= \exp \left( \sum_{i = 1}^n \ell(\alpha, \tau; y_i, x_i, z_i \rvert \beta) \right) \left( \sum_{i = 1}^n \left[ y_{i,m} z_i - \exp\left( x_i \alpha_m + z_i \beta_m \right)z_i \right] \right)\left( \sum_{i = 1}^n \left[ y_{i,j} z_i - \exp\left( x_i \alpha_j + z_i \beta_j \right)z_i \right] \right) \\
  \implies \frac{\partial^2}{\partial \beta_j \partial \beta_m} \mathcal{L}(\alpha, \tau; Y, X, Z, \rvert \beta) \bigg \rvert_{\beta = \vec{0}} &=  \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right) \left( \sum_{i = 1}^n[y_{i,m}z_i - \exp(x_i \alpha_m) z_i]\right) \left( \sum_{i = 1}^n [y_{i,j}z_i - \exp(x_i \alpha_j)z_i ]\right)
  \end{aligned}
  \label{eq:second-order-off-diag}
  $$
</details>
---

If we assume that terms of higher order than two are small enough, then we can use a Taylor approximation by truncating the expansion to:

$$
\mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \approx \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} + \beta^\top \left[ \frac{\partial}{\partial \beta} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right]  + \frac{1}{2} \beta^\top \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \beta
\label{eq:taylor}
$$

Denote the right-hand side of Eq. \eqref{eq:taylor} with $\mathcal{L}_T(\alpha, \tau; Y, X, Z \rvert \beta)$. Notice that Eq. \eqref{eq:lik} is the expectation of the conditional likelihood of $(\alpha, \tau)$ given $\beta$ (subject to the probability function used for $\beta$). We can then get around the complicated integral by taking the expectation of our Taylor approximation.

$$
\begin{aligned}
\mathcal{L}(\alpha, \theta; Y, X, Z) &= \int \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \mathcal{L}(\beta) d \beta \\
&\approx \int \mathcal{L}_T(\alpha, \tau; Y, X, Z \rvert \beta) \mathcal{L}(\beta) d\beta \\
&= \mathbb{E}\left[ \mathcal{L}_T(\alpha, \tau; Y, X, Z \rvert \beta)  \right] \\
&= \underbrace{\mathbb{E}\left[ \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right]}_{(a)} + 
\underbrace{\mathbb{E} \left[ \beta^\top \left[ \frac{\partial}{\partial \beta} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \vec{0}) \bigg \rvert_{\beta = \vec{0}} \right] \right]}_{(b)} + 
\frac{1}{2} \underbrace{\mathbb{E} \left[ \beta^\top \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \beta \right]}_{(c)} \\
&= \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}}
+ 0 + \frac{1}{2} \text{tr}\left[ \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \tau^2 \Sigma_T \right] \\
&\overset{(i)}{=} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}}
+ \frac{\tau^2}{2} \text{tr}\left[ \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \right] \\
&= \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}}
+ \frac{\tau^2}{2} \sum_{j = 1}^k \left[ \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \right]_{(j, j)} \\
&= \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right) + \frac{1}{2} \sum_{m = 1}^k  \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right) \left[ \left( y_{i,m} z_i - \exp(x_i \alpha_m) z_i \right)^2  - \sum_{i = 1}^n \exp(x_i \alpha_m)z_i^2 \right]  \\
&= \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right)  \left(1 + \frac{\tau^2}{2} \sum_{m = 1}^k  \left( y_{i,m} z_i - \exp(x_i \alpha_m) z_i \right)^2  - \sum_{i = 1}^n \exp(x_i \alpha_m)z_i^2  \right) 
\end{aligned}
\label{eq:approx-marg-lik}
$$

where $(i)$ follows from the fact that we take $\Sigma_T = \mathbb{I}_{k \times k}$. 

---
<details>
  <summary>Expectation Calculations - $(a)$</summary>
  $$
  \begin{aligned}
  \mathbb{E}\left[ \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] &= \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}}
  \end{aligned}
  $$
</details>
---
<details>
  <summary>Expectation Calculations - $(b)$</summary>
  $$
  \begin{aligned}
  \mathbb{E} \left[ \beta^\top \left[ \frac{\partial}{\partial \beta} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \vec{0}) \bigg \rvert_{\beta = \vec{0}} \right] \right] = \vec{0}^\top \frac{\partial}{\partial \beta} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \vec{0}) \bigg \rvert_{\beta = \vec{0}} = 0
  \end{aligned}
  $$
</details>
---
<details>
  <summary>Expectation Calculations - $(c)$</summary>
  $$
  \begin{aligned}
  \mathbb{E} \left[ \beta^\top \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \beta \right] 
  &= \mathbb{E}\left[\sum_{i = 1}^k \sum_{j = 1}^k  \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right]_{i,j} \beta_i \beta_j \right] \\
  &= \sum_{i =1}^k \sum_{j = 1}^k \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right]_{i,j} \mathbb{E}\left[ \beta_i \beta_j \right] \\
  &\overset{(i)}{=} \sum_{i =1}^k \sum_{j = 1}^k \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right]_{i,j} \text{Cov}(\beta_i, \beta_j) \\
  &= \sum_{i = 1}^k \sum_{j = 1}^k \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right]_{i,j} (\tau^2 \Sigma_T)_{i,j} \\
  &= \text{tr}\left[ \left[ \frac{\partial^2}{\partial \beta \beta^\top} \mathcal{L}(\alpha, \tau; Y, X, Z \rvert \beta) \bigg \rvert_{\beta = \vec{0}} \right] \tau^2 \Sigma_T \right]
  \end{aligned}
  $$

  $(i)$ is due to the fact that $\mathbb{E}[\beta] = \vec{0}$. 
</details>
---

We can then approximate the score with the derivative of Eq. \eqref{eq:approx-marg-lik} with respect to $\tau$:

$$
\begin{aligned}
U(\tau) &= \frac{d}{d\tau} \left[ \mathcal{L}(\alpha, \theta; Y, X, Z) \right] \\
&\approx \frac{d}{d\tau} \left[ \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right)  \left(1 + \frac{\tau^2}{2} \sum_{m = 1}^k  \left( y_{i,m} z_i - \exp(x_i \alpha_m) z_i \right)^2  - \sum_{i = 1}^n \exp(x_i \alpha_m)z_i^2  \right) \right]  \\
&= \tau \exp \left( \sum_{i = 1}^n \sum_{j = 1}^k \left[ y_{i,j} \left( x_i \alpha_j\right) - \exp\left( x_i \alpha_j \right) - \log(y_{i,j}!) \right] \right) \left( \sum_{m = 1}^k  \left( y_{i,m} z_i - \exp(x_i \alpha_m) z_i \right)^2  - \sum_{i = 1}^n \exp(x_i \alpha_m)z_i^2  \right) 
\end{aligned}
$$

-->