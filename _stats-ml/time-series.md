---
layout: distill
title: Time Series
description: A Primer
date: 2026-01-22
tabs: true
tags: theory time-series primer
toc:  
  - name: Background
  - name: Characteristics
    subsections:
        - name: Properties
        - name: Sample Functions
  - name: Examples
  - name: Models
  - name: Estimation
    subsections:
        - name: Non-Seasonal Model With Trend
        - name: Classical Decomposition Model
bibliography: 2026-01-22-time-series.bib
---

In this post, I'll work through some of the elementary concepts in time series analysis. I'll mostly focus on features of time series and popular (basic) modeling techniques. I'm doing this review to prepare for an interview, so I'll reserve more advanced methods for a later post. I'll also mostly use <i>Introduction to Time Series and Forecasting</i> by Brockwell and Davis.<d-cite key=brockwell2016></d-cite> 

---

## Background
A good place to start is with an explanation of what a time series is. In simple terms, a time series is a dataset where each observation is associated with some timepoint. Timepoints can be discrete or continuous, and they yield <i>discrete time series</i> and <i>continuous time series</i>, respectively. 

We will use $$\{ x_t \}_{t = 0}^{n}$$ to denote the observed values of a time series (where $n$ is permitted to be $\infty$). We assume that each $x_t$ is a realization of some random variable, $X_t$, and so the observed time series is a realization of the sequence $$\{ X_t \}_{t = 0}^m$$ (where $m \geq n$ and is permitted to be infinity).

In what follows, we'll omit the indexing of time for sequences of random variables (i.e. $$\{ X_t \}$$), but we will usually include the indexing for samples. The time series should be assumed to be either discrete or continuous unless context dictates otherwise. We'll also do most of our work for univariate time series with the extension to multiple dimensions being fairly straightforward.

---

## Characteristics
We now move to characteristics of time series, namely <i>stationarity</i> and <i>autocorrelation</i>.

<div class="definition">
<strong>Definition (Mean and Covariance).</strong><d-cite key=brockwell2016></d-cite>
<br>
Let $\{ X_t \}$ be a time series with finite second moment (i.e. $\mathbb{E}[X_t^2] < \infty$). We define the <i>mean</i> and <i>covariance</i> functions, respectively, of $\{ X_t \}$ as:

$$
\begin{aligned}
\mu_X(t) &= \mathbb{E}[X_t] \\
\gamma_X(r, s) &= \text{Cov}(X_r, X_s) = \mathbb{E}\left[ (X_r - \mu_X(r))(X_s - \mu_X(s)) \right]
\end{aligned}
$$

for all $r, s \in \mathbb{N}$. 
</div>

With these definitions, we come to stationarity. In words, a stationary time series will have a mean function that doesn't change over time.

<div class="definition">
<strong>Definition (Weak Stationarity).</strong><d-cite key=brockwell2016></d-cite>
<br>
A time series, $\{ X_t \}$, is called <i>weakly stationary</i> (or, often, just <i>stationary</i>) if:
<ul>
<li>$\mu_X(t)$ is independent of $t$</li>
<li>$\gamma_X(t+h, t)$ is independent of $t$ for each choice of $h$</li>
</ul>
</div>

Similarly, a time series, $$\{ X_t \}$$, is called <strong>strictly stationary</strong> if, for any $h, n \in \mathbb{Z}^+$, $(X_1, \dots, X_n)$ and $(X_{1 + h}, \dots, X_{n + h})$ have the same joint distribution.

<div class="definition">
<strong>Definition (Strict Stationarity).</strong><d-cite key=brockwell2016></d-cite>
<br>
A time series, $\{ X_t \}$, is called <i>strictly stationary</i> if:

$$
(X_1, \dots, X_n)^\top \overset{d}{=} (X_{1 + h}, \dots, X_{n + h})^\top
\hspace{10mm}
\text{for all } h \in \mathbb{Z}; \text{2mm} n \geq 1
$$
</div>

A strictly stationary time series will also be weakly stationary if $\mathbb{E}[X_t^2] < \infty$ for all $t$. Furthermore, a strictly stationary time series will be a sequence of identically distributed random variables. 

For a stationary time series, we can define the <i>autocovariance</i> and <i>autocorrelation functions</i>. These functions describe the dependence (or lack thereof) between the random variables within the sequence.

<div class="definition">
<strong>Definition (Autocovariance and Autocorrelation).</strong><d-cite key=brockwell2016></d-cite>
<br>
Let $\{ X_t \}$ be a (weakly) stationary time series, and let $h \in \mathbb{Z}^+$ (called a <i>lag</i>). The <i>autocovariance function (ACVF)</i> and <i>autocorrelation function (ACF)</i> of $\{ X_t \}$ at lag $h$ are defined as, respectively:

$$
\begin{aligned}
\gamma_X(h) &= \text{Cov}(X_{t + h}, X_t) \\
\rho_X(h) &= \frac{\gamma_X(h)}{\gamma_X(0)} = \text{Cor}(X_{t+h}, X_t)
\end{aligned}
$$
</div>

### Properties
The above functions have some nice properties for stationary processes, which we explore below.

<!-- #region autocovar-prop1 -->
<div class="theorem">
<strong>Claim.</strong>
{% tabs autocovar-prop1 %}
{% tab autocovar-prop1 statement %}
Let $\\{ X_t \\}$ be a stationary process. Then:
<ul>
<li>$\gamma_X(0) \geq 0$</li>
<li>$\rvert \gamma_X(h) \rvert \leq \gamma(0)$ for all $h \in \mathbb{Z}$</li>
<li>$\gamma_X(h) = \gamma(-h)$ for all $h \in \mathbb{Z}$ (i.e. $\gamma_X(\cdot)$ is even)</li>
</ul>
{% endtab %}
{% tab autocovar-prop1 proof %}
To show the first claim, note that:

$$
\gamma_X(0) = \text{Cov}(X_{t}, X_{t}) = \text{Var}(X_t)
$$

which necessarily must be non-negative. To show the second, recall that correlations have, at most, an absolute value of $1$. Then see that:

$$
\begin{aligned}
\rho_X(h) &= \frac{\gamma_X(h)}{\gamma_X(0)} \\
\implies \gamma_X(0) \rho_X(h) &= \gamma_X(h) \\
\implies \gamma_X(0) &\geq \gamma_X(h)
\end{aligned}
$$

The third claim follows from the following argument:

$$
\begin{aligned}
\gamma_X(h) &= \text{Cov}(X_{t + h}, X_t) \\
&= \text{Cov}(X_t, X_{t+h}) \\
&= \gamma_X(-h)
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

We also have that the autocovariance function must be positive semi-definite. A real-valued function, $f: \mathbb{Z} \rightarrow \mathbb{R}$, is <strong>positive semi-definite</strong> if, for any vector $\mathbf{a} \in \mathbb{R}^n$ and any $n \in \mathbb{Z}^+$:

$$
\begin{equation}
\label{eq:non-negative-definite}
\sum_{i = 1}^n \sum_{j = 1}^n \mathbf{a}_i f(i - j) \mathbf{a}_j \geq 0
\end{equation}
$$

<!-- #region autocovar-prop2 -->
<div class="theorem">
{% tabs autocovar-prop2 %}
{% tab autocovar-prop2 statement %}
A real-valued function, $f: \mathbb{Z} \rightarrow \mathbb{R}$, is the autocovariance function for some stationary time series if, and only if, it is even (i.e $f(x) = f(-x)$) and positive semi-definite.
{% endtab %}
{% tab autocovar-prop2 proof %}
Fix $\mathbf{a} \in \mathbb{R}^n$, and let $\mathbf{X}_n = (X_1, \dots, X_n)^\top$. Define the covariance matrix of $\mathbf{X}_n$ as:

$$
\Gamma_n = \mathbb{E}\left[ (\mathbf{X}_n - \mathbb{E}[\mathbf{X}_n])(\mathbf{X}_n - \mathbb{E}[\mathbf{X}_n])^\top \right]
$$

We then have:

$$
\begin{aligned}
\text{Var}(\mathbf{a}^\top \mathbf{X}_n) 
&= \mathbf{a}^\top \text{Var}(\mathbf{X}_n) \mathbf{a} \\
&= \mathbf{a}^\top \mathbb{E}\left[ (\mathbf{X}_n - \mathbb{E}[\mathbf{X}_n])(\mathbf{X}_n - \mathbb{E}[\mathbf{X}_n])^\top \right] \mathbf{a} \\
&= \mathbf{a}^\top \Gamma_n \mathbf{a}
&\geq 0
\end{aligned}
$$

where the inequality follows from the fact that $\Gamma_n$ is a covariance matrix (so it is positive semi-definite), and the equation is a quadratic form so it must be non-negative (see <a href="https://en.wikipedia.org/wiki/Definite_matrix#Definitions">here</a>).

The converse is more involved and omitted.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

The above properties extend to the autocorrelation function with the addition that $\rho_X(\cdot)$ also satisfies $\rho_X(0) = 1$. 

### Sample Functions
The above functions were defined for sequences of random variables. We can define analogous functions for observations.

<div class="definition">
<strong>Definition (Sample Mean, Autocovariance, and Autocorrelation).</strong><d-cite key=brockwell2016></d-cite>
<br>
Let $\{ x_t \}_{t = 0}^n$ be observations of a time series. The <i>sample mean</i>, <i>sample autocovariance</i>, and <i>sample autocorrelation</i> functions are defined as:

$$
\begin{aligned}
\bar{x} &= \frac{1}{n} \sum_{t = 0}^n x_t \\
\hat{\gamma}(h) &= \frac{1}{n} \sum_{t = 0}^{n - \rvert h \rvert} (x_{t + \rvert h \rvert} - \bar{x})(x_t - \bar{x}) \\
\hat{\rho}(h) &= \frac{\hat{\gamma}(h)}{\hat{\gamma}(0)}
\end{aligned}
$$

for $-n < h < n$. 
</div>

---

## Examples
We'll introduce a few important examples from the book to help illustrate some of the concepts.

<div class="definition">
<strong>Definition (White Noise).</strong><d-cite key=brockwell2016></d-cite>
<br>
<i>White noise</i> is a sequence of random variables, $\{ X_t \}$,  where the $X_t$ are i.i.d. with mean zero and variance $\sigma^2 < \infty$. We denote this as:

$$
\{ X_t \} \sim \text{WN}(0, \sigma^2)
$$
</div>

Since $\mathbb{E}\left[ X_t \right] = 0$ for all $t$ and $\gamma_X(t+h, t) = 0$ for all $h > 0$ (and equal to $\sigma^2$ for $h = 0$) by independence, white noise is a stationary time series. 

<div class="definition">
<strong>Definition (Random Walk).</strong><d-cite key=brockwell2016></d-cite>
<br>
A <i>random walk</i>, which we denote with $\{ S_t \}_{t = 0}^\infty$, is defined as a sequence of independent random variables satisfying $S_0 = 0$ and $S_t = \sum_{i = 0}^{t} X_t$, where $\{ X_t \}$ is a sequence of i.i.d. random variables. 
<br>
If $\{ X_t \}$ is a <i>binary process</i> (i.e. $\mathbb{P}(X_t = 1) = \frac{1}{2}$ and $\mathbb{P}(X_t = -1) = \frac{1}{2}$), then $\{ S_t \}$ is a <i>simple symmetric random walk</i>. 
</div>

Let $\{ X_t \}$ be i.i.d. noise with mean zero and variance $\sigma^2 < \infty$. The corresponding random walk $\{ S_t \}$ also has mean zero by the linearity of expectation. We also have:

$$
\begin{aligned}
\mathbb{E}\left[ S_t^2 \right] &= \mathbb{E}\left[ \left( \sum_{i = 0}^t X_i \right)^2 \right] \\
&= \sum_{i = 0}^t \sum_{j = 0}^t \mathbb{E}\left[ X_i X_j \right] \\
&= \sum_{i = 0}^t \mathbb{E}\left[ X_i^2 \right] & \left( X_i, X_j \text{ ind.}\right) \\
&= t \sigma^2
\end{aligned}
$$

However: 

$$
\begin{aligned}
\gamma_S(t + h, t) &= \text{Cov}(S_{t+h}, S_t) \\
&= \text{Cov}(S_t + \sum_{i = t+1}^{t+h} X_i, S_t) \\
&= \text{Cov}(S_t, S_t) \\
&= t \sigma^2
\end{aligned}
$$

Thus, this random walk is not stationary.

---

## Models
We now introduce several stationary processes that are very important to time series analysis.

<div class="definition">
<strong>Definition ($MA(q)$ Process).</strong><d-cite key=brockwell2016></d-cite>
<br>
A <i>moving average process of order $q$</i> is a time series $\{ X_t \}$ with the form:

$$
X_t = Z_t + \theta_1 Z_{t - 1} + \dots + \theta_q Z_{t - q}
$$

where $\{ Z_t \} \sim WN(0, \sigma^2)$ and $\theta_1, \dots, \theta_q \in \mathbb{R}$ and constant.
</div>

An $MA(q)$ process will be stationary, and it will also be strictly stationary because $\{ Z_t \}$ is i.i.d. noise. 

<div class="definition">
<strong>Definition (Linear Process).</strong><d-cite key=brockwell2016></d-cite>
<br>
A <i>linear process</i> is a time series $\{ X_t \}$ with the form:

$$
X_t = \sum_{j = -\infty}^{\infty} \psi_j Z_{t - j} \hspace{10mm} \text{for all } t
$$

where $\{ Z_t \} \sim WN(0, \sigma^2)$ and $\{ \psi_j \}$ is a sequence of constants satisfying $\sum_{j = -\infty}^\infty \rvert \psi_j \rvert < \infty$.
</div>

A linear process can be written with a backward shift operator as:

$$
\begin{aligned}
\psi(B) &= \sum_{j = -\infty}^\infty \psi_j B^j \\
\X_t &= \psi(B) Z_t
\end{aligned}
$$

If $\psi_j = 0$ for all $j < 0$, the we call the linear process a <strong>moving average</strong> (or an <strong>MA($\infty$)</strong>). 

<div class="definition">
<strong>Definition (Autoregressive Process).</strong><d-cite key=brockwell2016></d-cite>
<br>
An <i>autoregressive process of order $p$</i> is a time series $\{ X_t \}$ with the form:

$$
X_t = \sum_{i = 1}^p \phi_i X_{t - i} + Z_t
\hspace{10mm} \text{for } t \in \mathbb{Z}
$$

where $\{ Z_t \} \sim WN(0, \sigma^2)$, $Z_t$ is uncorrelated with $X_s$ for $s < t$, and $\rvert \phi_i \rvert < 1$ for $i = 1, \dots, p$ and constant.
</div>

Combining a moving average model with an autoregressive model yields an <i>autoregressive moving average</i> model.

<div class="definition">
<strong>Definition (Autoregressive Moving Average Process).</strong><d-cite key=brockwell2016></d-cite>
<br>
An <i>autoregressive moving average process (ARMA($p, q$))</i> is a stationary time series $\{ X_t \}$ with the form:

$$
\begin{aligned}
X_t - \sum_{j = 1}^p \phi_j X_{t - j} &= Z_t + \sum_{i = 1}^q \theta_i Z_{t - i} \\
\implies X_t &= Z_t + \sum_{i = 1}^q \theta_i Z_{t - i} + \sum_{j = 1}^p \phi_j X_{t - j} 
\hspace{10mm} \text{for } t \in \mathbb{Z}
\end{aligned}
$$

where $\{ Z_t \} \sim WN(0, \sigma^2)$ and the polynomials $\phi(y) = 1 - \sum_{i = 1}^p \phi_i y^i$ and $\theta(y) = 1 + \sum_{j = 1}^q \theta_j y^j$ have no common factors.
</div>

An ARMA process can also be written with a backward shift operator. If we define $B^jX_t = X_{t - j}$ and $B^j Z_t = Z_{t - j}$ for $j \in \mathbb{Z}$, then:

$$
\phi(B) X_t = \theta(B) Z_t
$$

We can also define a special function for ARMA processes. The <strong>partial autocorrelation function (PACF)</strong> of an ARMA process is the function $\alpha(h) = \phi_{h, h}$ (with $\alpha(0) = 1$) for $h \geq 1$ where:

$$
\phi_{h,h} = \left[
\begin{bmatrix}
\gamma_X(1-1) & \gamma_X(1-2) & \dots & \gamma_X(1-h) \\
\gamma_X(2-1) & \gamma_X(2-2) & \dots & \gamma_X(2-h) \\
\vdots & \vdots & \ddots & \vdots \\
\gamma_X(h-1) & \gamma_X(h-2) & \dots & \gamma_X(h-h)
\end{bmatrix}
\begin{bmatrix}
\gamma_X(1) \\
\vdots \\
\gamma_X(h)
\end{bmatrix}
\right]_{h}
$$

A sample version can also be computed using the sample autocovariance function.

---

## Estimation
There are many different models that can be used for time series analysis. Here, we'll cover some of the most foundational and common.

### Non-Seasonal Model With Trend

The <i>non-seasonal model with trend</i> has the form:

$$
\begin{equation}
\label{eq:non-seasonal-with-trend}
X_t = m_t + Y_t
\end{equation}
$$

Here, $m_t$ is the <strong>trend component</strong>. This part is meant to capture long-term/slow moving change in $X_t$. We can think of the trend component as representing the overall, more general, change in the sequence. $Y_t$ is weakly stationary random noise and captures anything in the sequence that cannot be captured by $m_t$.

We are often interested estimating $m_t$ so that we can remove its effect, identify a reasonable (stationary!) noise distribution, and then use our knowledge of stationary processes for prediction and inference. There are several popular ways to do the initial estimation. We'll assume we have some process, $$\{ X_t \}$$, as described in Eq. \eqref{eq:non-seasonal-with-trend}.

#### Moving Average
One way to estimate $m_t$ is to use a <strong>moving average filter</strong>. For finite $q \in \mathbb{N}$, the <strong>two-sided moving average </strong> of $$\{ X_t \}$$ is defined as:

$$
W_t = \frac{1}{2q + 1} \sum_{j = - q}^q X_{t - j}
$$

Under the assumption that $m_t$ is approximately linear over $[t - q, t + q]$, and that the mean of the noise component over $[t - 1, t + q]$ is approximately zero, we can estimate $m_t$ as $W_t$:

$$
\begin{aligned}
W_t &= \frac{1}{2q + 1} \sum_{j = - q}^q (m_{t - j} + Y_{t - j}) \approx m_t \\
\implies
\hat{m}_t &= \sum_{j = -\infty}^\infty a_j X_{t - j} \hspace{10mm} \text{for } a_j = \frac{1}{2q + 1}; \hspace{2mm} -q \leq j \leq q
\end{aligned}
$$

One drawback is that this filter is two-sided, so it uses information from before time $t$ when estimating $m_t$. Additionally, our samples are usually finite (e.g. over time $[0, n]$), so the moving average estimate given above cannot be used for $t \leq q$ or $t > n - q$.

#### Exponential Smoothing
<strong>Exponential smoothing</strong> provides a one-sided alternative to the moving average filter. One can fix $\alpha \in [0, 1]$ and estimate $m_t$ with:

$$
\hat{m}_t = \alpha X_t + (1 - \alpha) \hat{m}_{t - 1} \hspace{10mm} \text{for } t = 2, \dots, n
$$

where we define $\hat{m}_1 = X_1$. The name comes from the fact that the recursion compounds the multiplication by $\alpha$ for $t \geq 2$:

$$
\hat{m}_t =  (1 - \alpha)^{t-1} X_1 + \sum_{j = 0}^{t - 2} \alpha(1 -\alpha)^j X_{t - j} 
$$

#### Polynomial Fitting
A simple alternative, though still useful, is <strong>polynomial fitting</strong> in which we simply fit a polynomial of degree $k$ that has the form:

$$
m_t = a_0 + a_1 t + a_2 t^2 + \dots + a_k t^k
$$

These polynomials are usually fit via <a href="/stats-ml/regression/#estimation">least squares</a>. 

#### Differencing
To explain <strong>differencing</strong>, we first define the <strong>lag-$1$ difference operator</strong> and <strong>backward shift operator</strong>, respectively:

$$
\begin{aligned}
\nabla X_t &= X_t - X_{t-1} = (1 - B) X_t \\
B X_t &= X_{t-1}
\end{aligned}
$$

We can compound both to the $k$-th degree as:

$$
\begin{aligned}
\nabla^k X_t &= \nabla(\nabla^{k - 1} X_t ) = (1 - B)^k X_t\\
B^k X_t &= X_{t - k}
\end{aligned}
$$

Under the assumption that the trend component is a $k$-degree polynomial, we can apply the difference operator $k$ to yield:

$$
\begin{equation}
\label{eq:differencing}
\begin{aligned}
\nabla^k X_t &= \nabla^k m_t + \nabla^k Y_t \\
&= \nabla^{k}\left[ \sum_{j = 0}^k c_j t^j \right] + \nabla^k Y_t \\
&= k! c_k + \nabla^k Y_t
\end{aligned}
\end{equation}
$$

If $X_t = m_t + Y_t$ and $Y_t$ is stationary with mean zero, then Eq. \eqref{eq:differencing} is stationary with mean $k! c_k$. 


### Classical Decomposition Model
The <i>classical decomposition model</i> has the form:

$$
\begin{equation}
\label{eq:classical-decomposition-model}
X_t = m_t + s_t + Y_t
\end{equation}
$$

The $m_t$ and $Y_t$ components are the same as in the <a href="#non-seasonal-model-with-trend">non-seasonal model with trend</a>, and $s_t$ is the <strong>seasonal component</strong> and has known <strong>period</strong>, $d$. The seasonal component captures patterns that might occur repeatedly in the sequence and should satisfy $\sum_{j = 1}^d s_j = 0$. The period represents the amount of time between repetitions, and thus $s_{t + d} = s_t$. 

#### Moving Average Filter
Because there is an additional component in the classical decomposition model, we must do several estimation steps to remove both the trend and seasonality to form a noise process. 

For example, one way is to first apply a filter to estimate the trend. For an odd period, we can use a two-sided moving average centered at $t$ and of window equal to the period. For an even period (e.g. $d = 2q$), we can instead do a slightly modified moving average:

$$
\hat{m}_t = \frac{1}{d} \left(\frac{1}{2} x_{t - q} + x_{t - q + 1} + \dots + x_{t + q - 1} + \frac{1}{2} x_{t+q} \right)
\hspace{10mm} \text{for } q < t \leq n - q
$$

We can compute the sequence of deviations $$\{ x_t - \hat{m}_t \}_{t = q + 1}^{n - q}$$. We can then compute the average of the deviations at each "position" in the period, which we denote with $w_k$, as the average of the set $$\{ x_t - \hat{m}_t \}$$ for $q < k + jd \leq n - q$ for $k = 1, \dots, d$. We next estimate the seasonal component as:

$$
\begin{aligned}
\hat{s}_k &= w_k - \frac{1}{d} \sum_{i = 1}^d w_i
\hspace{10mm} \text{for } k = 1, \dots, d \\
\hat{s}_k &= \hat{s}_{k - d} \hspace{10mm} \text{for } k > d
\end{aligned}
$$

Removing the estimate for the seasonal component yields the <strong>deseasonalized</strong> data:

$$
d_t = x_t - \hat{s}_t
$$

and we can use one of the methods from <a href="#non-seasonal-model-with-trend">the previous section</a> to estimate the trend for this non-seasonal data. The estimated noise series is:

$$
\hat{Y}_t = x_t - \hat{m}_t - \hat{s}_t
$$

where $\hat{s}_t$ is the estimated seasonal component, and $\hat{m}_t$ is the estimated trend component for the non-seasonal data.


#### Differencing
Similar to the non-seasonal model, we can use differencing for model with a seasonality component. The <strong>lag-$d$ differencing operator</strong> is defined as:

$$
\nabla_d X_t = X_t - X_{t - d} = X_t - B^d X_t
$$

Using this operator on the classical decomposition model where $s_t$ has a period of $d$ yields:

$$
\nabla_d X_t = m_t - m_{t - d} + Y_t - Y_{t - d}
$$

The above is simply a non-seasonal model with trend component $m_t - m_{t-d}$ and noise $Y_t - Y_{t - d}$, and the methods from <a href="#non-seasonal-model-with-trend">the previous section</a> can be used. 
