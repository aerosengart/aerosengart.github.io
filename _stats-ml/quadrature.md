---
layout: distill
title: Quadrature
description: A Primer
date: 2026-02-05
tabs: true
tags: methods approximation theory primer
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Step Function Quadrature
      subsections:
        - name: Left and Right Rules
        - name: Midpoint Rule
        - name: Trapezoidal Rule
    - name: Interpolating Function Quadrature
      subsections:
        - name: Newton-Cotes Rules
        - name: Simpson's Rules
    - name: Gaussian Quadrature
      subsections:
        - name: Gauss-Legendre Quadrature
        - name: Gauss-Jacobi Quadrature
        - name: Gauss-Laguerre Quadrature
        - name: Gauss-Hermite Quadrature
    - name: Adaptive Quadrature
    - name: Example
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: stats-ml.bib
---

Since I've been working on integral approximations via marginal quasi-likelihood (and some penalized quasi-likelihood), it's probably a good idea to review <i>numerical integration</i> techniques to see how good the approximations that the more analytical techniques provide are. I'll focus on one-dimensional cases for now.

In essence, numerical integration<d-cite key=numericalintegration2025></d-cite> is a way to compute an integral without evaluating it analytically. In one-dimensional cases, it is fairly simple to approximate definite integrals (under some niceness conditions) by evaluating the function at a finite number of points and adding together the results in a weighted sum. The way this sum is defined is called a <strong>quadrature rule</strong>, and using these rules to approximate an integral is called <strong>quadrature</strong>.

<aside><p>That is, a Riemann sum. As such, we will assume that all integrands are Riemann integrable.</p></aside>

In what follows, we will be considering approximations of the following one-dimensional integral with $f(x)$ continuous on $[a, b]$:

$$
\int_a^b f(x) dx
$$

We will devide the domain of integration into $n \in \mathbb{N}$ sub-intervals each of length $\delta_x = \frac{b - a}{n}$. We will let $P = \{ x_0, x_1, \dots, x_n \}$ be the set of endpoints of these sub-intervals. 

<aside><p>If $n > 2$, then these are called <strong>composite</strong> rules.</p></aside>

---

## Step Function Quadrature

### Left and Right Rules
The <strong>left rule</strong> and <strong>right rule</strong> are quadratures rules that uses the function's output at the left and right, respectively, endpoints of sub-intervals defined on the domain of integration as the weights for the summation.<d-cite key=riemann2025></d-cite>

The left rule is given by:

$$
\begin{equation}
\label{eq:left}
L_n = \sum_{i = 1}^n \delta_x f(x_{i - 1})
\end{equation}
$$

Similarly, the right rule is:

$$
\begin{equation}
\label{eq:right}
R_n = \sum_{i = 1}^n \delta_x f(x_{i})
\end{equation}
$$

In words, these rules approximate the integral using the sum of the areas of the rectangles that have widths equal to $\delta_x$ and heights equal to the function's output at the left or right endpoints of the sub-intervals.

The left rule and right rule can both result in overestimation or underestimation in certain cases. If the integrand is monotonically decreasing on $[a, b]$, then the left rule will overestimate while the right rule will underestimation. If it is monotonically increasing, then the left rule will underestimate while the right rule will over estimate.

If we let $S_n$ be either $L_n$ or $R_n$, the error of the approximation is bounded by:

$$
\left\rvert \int_a^b f(x) dx - S_n \right\rvert \leq \frac{(b - a)^2}{2n} \underset{x \in [a, b]}{\max} \left\{ \rvert f'(x) \rvert \right\}
$$

### Midpoint Rule
The <strong>midpoint rule</strong><d-cite key=riemann2025></d-cite> is essentially splits the difference between the left and right rules by using the points halfway between the left and right endpoints of the sub-intervals as the weights.

The midpoint rule is given by:

$$
\begin{equation}
\label{eq:midpoint}
M_n = \sum_{i = 1}^n f(m_i) \delta_x
\end{equation}
$$

The midpoint rule's error can be bounded as:

$$
\left\rvert \int_a^b f(x) dx - M_n \right\rvert \leq \frac{(b - a)^3}{24n^2} \underset{x \in [a, b]}{\max} \left\{ \rvert f''(x) \rvert \right\}
$$

<aside><p>The midpoint rule can be generalized to obtain a better approximation error, but it is kind of complicated so we omit the details here.</p></aside>


---

## Interpolating Function Quadrature

### Trapezoidal Rule
The previous three rules all defined <i>step functions</i> on the sub-intervals over the domain of integration. We can also extend our approximations to more complex functions. 

The <strong>trapezoidal rule</strong> is similar to the midpoint rule but uses trapezoids rather than rectangles. This is a fairly simple interpolating function; rather than a step function, we have extended to general affine functions (i.e. linear functions, or polynomials of degree $1$)

Recall that the area of a trapezoid with height, $h$, and base lengths $b_1$ and $b_2$ is given by $\frac{1}{2} h(b_1 + b_2)$. The trapezoids we define will have heights equal to the length of the sub-intervals, and the base lengths will be the function's outputs at the two endpoints of each sub-interval.

We have:

$$
\begin{equation}
\label{eq:trapezoid}
T_n = \sum_{i = 1}^{n} \frac{1}{2} \delta_x (f(x_{i - 1}) + f(x_{i}))
\end{equation}
$$

By construction, the trapezoidal rule's approximation can also be obtained as the average of the left and right rules'. The trapezoidal rule will also tend to overestimate the integral on sub-intervals where the function is concave up and underestimate on sub-intervals where the function is concave down. Its error can be bounded as:

$$
\left\rvert \int_a^b f(x) dx - T_n \right\rvert \leq \frac{(b - a)^3}{12n^2} \underset{x \in [a, b]}{\max} \left\{ \rvert f''(x) \rvert \right\}
$$

This is exactly double the error bound for the midpoint rule's approximation. 

### Newton-Cotes Rules
To get a little more complicated, we can use a <strong>Newton-Cotes rule</strong><d-cite key=newtoncotes2026></d-cite>, which uses the integral of Langrange basis polynomials as the weights. 

<aside><p>Some of the other rules are specific types of Newton-Cotes rules!</p></aside>

Assuming the endpoints in $P$ are distinct (which should be true), the <strong>Lagrange basis</strong> for polynomials of degree $\leq n$ (for $P$) is the set of polynomials of degree $n$, $\{ l_0(x), l_1(x), \dots, l_n(x) \}$, where each polynomial is defined as:

$$
l_i(x) = \prod_{0 \leq m \leq n}{m \neq i} \frac{x - x_m}{x_i - x_m}
$$

The <strong>Lagrange interpolating polynomial</strong>, also called the <strong>Lagrange basis polynomial</strong> is the unique polynomial that interpolated the points in $P$. It is defined as:

$$
L(x) = \sum_{i = 0}^n x_i l_i(x)
$$

A Newton-Cotes rule is then defined as:

$$
\begin{equation}
\label{eq:newton-cotes}
C_n = \sum_{i = 0}^n f(x_i) \int_a^b l_i(x) dx
\end{equation}
$$

which arises from the approximation of the function with some Lagrange interpolating polynomial. 


### Simpson's Rules
Simpson's rules<d-cite key=simpson2026></d-cite> are a few different types of interpolating quadrature rules. If $n > 2$ is even, then <strong>Simpson's 1/3 rule</strong> is given by:

$$
\begin{equation}
\label{eq:simpson-1-3}
S^{1/3}_n = \frac{1}{3} \delta_x \left[ f(x_) + 4 \sum_{i = 1}^{\frac{n}{2}} f(x_{2i - 1}) + 2 \sum_{i = 1}^{\frac{n}{2} - 1} f(x_{2i}) + f(x_n) \right]
\end{equation}
$$

We can bound the approximation error as:

$$
\begin{aligned}
\left\rvert \int_a^b f(x) dx - S^{1/3}_n \right\rvert \leq \frac{\delta_x^3(b - a)}{180} \underset{x \in [a, b]}{\max} \left\{ \rvert f^(4)(x) \rvert \right\}
\end{aligned}
$$

If $n$ is a multiple of $3$, we can use <strong>Simpson's 3/8 rule</strong>, which is defined as:

$$
\begin{equation}
\label{eq:simpson-3-8}
S^{3/8}_n = \frac{3}{8} \delta_x \left[ f(x_) + 3 \sum_{i = 1;  i \nmid 3}^{n - 1} f(x_i) + 2 \sum_{i = 1}^{\frac{n}{3} - 1} f(x_{3i}) + f(x_n) \right]
\end{equation}
$$

<aside><p>$i \nmid 3$ means $i$ is not divisible by $3$.</p></aside>

Simpson's rules are examples of <i>closed</i> Newton-Cotes rules. 

<aside><p>A <strong>closed</strong> Newton-Cotes rule is such that $x_0 = a$ and $x_n = b$. If $x_0 > a$ and $x_n < b$, then the rule is <strong>open</strong>.</p></aside>


---

## Gaussian Quadrature
A <strong>Gaussian quadrature rule</strong><d-cite key=gaussquad2025></d-cite> using interpolating functions as well. In general, a Gaussian quadrature rule is an approximation of the form:

$$
\int_{a}^b \omega(x) f(x) dx \approx \sum_{i = 1}^n w_i f(x_i)
$$

where $\{ x_1, \dots, x_n \}$ are the roots of some specially chosen polynomial of degree $n$ from a family of orthogonal polynomials.

### Gauss-Legendre Quadrature
In the simplest setting, we have $[a, b] = [-1, 1]$ and $\omega(x) = 1$. The <a href="https://en.wikipedia.org/wiki/Legendre_polynomials"><strong>Legendre polynomials</strong></a> are the system of polynomials satisfying, for $k = 0, 1, 2, \dots, n$:

$$
\begin{aligned}
P_k(1) &= 1 \\
\int_{-1}^1 P_m(x) P_k(x) dx = 0 \text{ if } k \neq m
\end{aligned}
$$

This implies $P_0(x) = 1$, $P_1(x)$, etc. <strong>Gauss-Legendre quadrature</strong><d-cite key=gausslegendre2025></d-cite> is then:

$$
\begin{equation}
\label{eq:legendre}
G_n = \sum_{i = 1}^n \frac{2}{(1 - x_i)^2 \left[ P'_n(x_i) \right]^2} f(x_i)
\end{equation}
$$

where $\{ x_1, \dots, x_n \} are the roots of the $n$-th degree Legendre polynomial. 

<aside><p>For polynomials of degree $2n - 1$, this approximation will be exact.</p></aside>

### Gauss-Jacobi Quadrature
Slightly more complicated is <strong>Gauss-Jacobi quadrature</strong>,<d-cite key=gaussjacobi2025></d-cite> which supposes $[a, b] = [-1, 1]$ but that $\omega(x) = (1 - x)^\alpha(1 + x)^\beta$ for $\alpha, \beta > -1$. It is useful over Gauss-Legendre quadrature because it can be uses when $f$ is not defined at the endpoints, $a$ and/or $b$.

It uses <a href="https://en.wikipedia.org/wiki/Jacobi_polynomials"><strong>Jacobi polynomials</strong></a>, which are polynomials of the form:

$$
\begin{aligned}
P_n^{(\alpha, \beta)}(z) &= \frac{(\alpha+1)^{(n)}}{n!}  _2F_1 \left(-n, 1 + \alpha + \beta + n; \alpha + 1; \frac{1 - z}{2}\right) \\
&= \frac{\Gamma(\alpha + n + 1)}{n! \Gamma(\alpha + \beta + n + 1)}  \sum_{m = 0}^n {n \choose m} \frac{\Gamma(\alpha + \beta + n + m + 1)}{\Gamma(\alpha + m + 1)} \left(\frac{z - 1}{2} \right)^m
\end{aligned}
$$

where:

$$
\begin{aligned}
(q)^{(n)} &= \prod_{k = 1}^n (q + k - 1) & \left(\text{rising factorial}\right) \\
_2 F_1 (a, b; c; z) &= \sum_{n = 0}^\infty \frac{(a)^{(n)} (b)^{(n)}}{(c)^{(n)}} \frac{z^n}{n!} \text{ for } \rvert z \rvert < 1 & \left(\text{hypergeometric function}\right)
\end{aligned}
$$

and $\Gamma(\cdot)$ denotes the <a href="https://en.wikipedia.org/wiki/Gamma_function">Gamma function</a>.

Letting $\{ x_1, \dots x_n \}$ be the set of roots of the $n$-th Jacobi polynomial, the Gauss-Jacobi quadrature approximation on $n$ points is:

$$
\begin{equation}
\label{eq:jacobi}
\begin{aligned}
GJ_n &= \sum_{i = 1}^n \lambda_i f(x_i) \\
\lambda_i &= - \left(\frac{2n + \alpha + \beta + 2}{n + \alpha + \beta + 1}\right)\left(\frac{\Gamma(n + \alpha + 1) \Gamma(n + \beta + 1)}{\Gamma(n + \alpha + \beta + 1)(n + 1)!} \right)\left(\frac{2^{\alpha + \beta}}{\left(P_n^{(\alpha,\beta)}(x_i)\right)' P_{n + 1}^{(\alpha, \beta)}(x_i)}\right)
\end{aligned}
\end{equation}
$$

Its error can be bounded as:

$$
\begin{aligned}
\left\rvert \int_{-1}^1 f(x) (1 - x)^\alpha (1 + x)^\beta dx - GJ_n \right\rvert \leq \left(\frac{\Gamma(n + \alpha + 1)\Gamma(n + \beta + 1) \Gamma(n + \alpha + \beta + 1)}{(2n + \alpha + \beta + 1)\left[ \Gamma(2n + \alpha + \beta + 1) \right]^2}\right) \left(\frac{2^{2 + \alpha + \beta + 1}}{(2n)!}\right) \underset{x \in (-1, 1)}{\max} \left\{ \rvert f^(2n)(x) \rvert \right\}
\end{aligned}
$$

### Gauss-Laguerre Quadrature
If we instead have $[a, b) = [0, +\infty)$ and $\omega(x) = \exp(-x)$, we can use <strong>Gauss-Laguerre quadrature</strong><d-cite key=gausslaguerre2025></d-cite>. It uses the roots of the $n$-th <a href="https://en.wikipedia.org/wiki/Laguerre_polynomials"><strong>Laguerre polynomial</strong></a>. These polynomials are defined to satisfy:

$$
\begin{aligned}
L_0(x) &= 1 \\
L_1(x) &= 1 - x \\
L_{k + 1}(x) = \frac{(2k + 1 - x)L_k(x) - kL_{k - 1}(x)}{k + !} \text{ for } k \geq 1
\end{aligned}
$$

Thus, the $n$-th Laguerre polynomial is given by:

$$
L_n(x) = \sum_{k = 0}^n {n \choose k} \frac{(-1)^k}{k!} x^k \text{ for } n = 0, 1, \dots
$$

The Gauss-Laguerre approximation is then:

$$
\begin{equation}
\label{eq:laguerre}
GL_n(x) = \sum_{i = 1}^n \left( \frac{x_i}{(n + 1)^2} \left[ L_{n + 1}(x_i) \right]^2 \right) f(x_i)
\end{equation}
$$

Note that this quadrature rule can be used for general function, $g(x)$, by defining $f(x) =  \exp(x) g(x)$, yielding:

$$
g(x) = g(x) \exp(x) \exp(-x) = \exp(-x) f(x)
$$

### Gauss-Hermite Quadrature
Finally, suppose we have an indefinite integral with $\omega(x) = \exp\left(-x^2\right)$. We can then choose to use the (physicist's) <a href="https://en.wikipedia.org/wiki/Hermite_polynomials"><strong>Hermite polynomials</strong></a>, which are defined as:

$$
H_n(x) = (-1)^n \exp\left( x^2\right) \frac{d^n}{dx^n} \left[ \exp\left(- x^2 \right) \right]
\text{ for } n = 0, 1, \dots
$$

The <strong>Gauss-Hermite quadrature</strong><d-cite key=gausshermite2025></d-cite> approximation is then:

$$
GH_n(x) = \sum_{i = 1}^n \frac{2^{n-1} n! \sqrt{\pi}}{n^2 \left[ H_{n-1}(x_i)\right]^2} f(x_i)
$$

where $\{ x_1, \dots, x_n \}$ are the roots of the $n$-th Hermite polynomial.

This approximation has the following error bound:

$$
\left\rvert \int_{-\infty}^{+\infty} f(x) \exp(-x^2) dx - GH_n \right\rvert \leq \frac{n! \sqrt{\pi}}{2^n (2n)!}\underset{x \in (-\infty, +\infty)}{\max} \left\{ \rvert f^(2n)(x) \rvert \right\}
$$


---

## Adaptive Quadrature
The basic idea behind <strong>adaptive quadrature</strong><d-cite key=adaptquad2025></d-cite> is to repeatedly approximate an integral using some other quadrature rule by defining smaller and smaller sub-intervals.

<ol>
<li>Approximate integral with quadrature.</li>
<li>Compute absolute error of approximation.</li>
<li>If error is too large, split domain of integration in half.</li>
<li>Repeat for both halves of the domain.</li>
</ol>


---

## Example
To illustrate its use, let's use Gauss-Hermite quadrature to approximate the derivation of the marginal log-likelihood in a <a href="/stats-ml/glmm">generalized linear mixed model</a> with just a random intercept. We'll assume our $N$ datapoints come from $K$ clusters, and the $i$-th response in the $k$-th cluster will be denoted by $y_i^k$. For simplicity, we'll assume the clusters are all of size $n$.

We'll further assume an explicit distribution for our responses (Poisson) and the random effects (Gaussian), so our set-up is:

$$
\begin{aligned}
\mathbb{E}\left[ y_i^k \rvert \beta_k \right] &= \mu_i^k \\
\log(\mu_i^k) &= \alpha + \beta_k \\
\beta_k &\overset{iid}{\sim} \mathcal{N}(0, \tau^2)
\end{aligned}
$$

We'll let $\mathbf{y} = (y_1^1, \dots, y_n^1, \dots, y_1^K, \dots, y_n^K)^\top$ and $\beta = (\beta_1, \dots, \beta_K)^\top$. To obtain the marginal likelihood for each cluster, we multiply by the likelihood of its random intercept and then integrate that parameter out:

$$
\mathcal{L}(\mathbf{y}^k; \alpha, \tau^2) = \int \mathcal{L}(\mathbf{y}^k; \alpha, \tau^2 \rvert \beta_k) \mathcal{L}(\beta_k) d \beta_k
$$

Since we assumed Poisson responses and Gaussian random intercepts, we have:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}^k; \alpha, \tau^2) 
&= \int \left[ \mathcal{L}(\mathbf{y}^k; \alpha, \tau^2 \rvert \beta_k) \mathcal{L}(\beta_k) \right] d \beta_k \\
&= \int \left[\mathcal{L}(\mathbf{y}^k; \alpha, \tau^2 \rvert \beta_k) \left( \frac{1}{\sqrt{2 \pi \tau^2}} \exp\left(- \frac{\beta_k^2}{2 \tau^2} \right) \right) \right] d \beta_k \\
&= \frac{1}{\sqrt{2 \pi \tau^2}} \int \left[  \mathcal{L}(\mathbf{y}^k; \alpha, \tau^2 \rvert \beta_k) \exp\left(- \frac{\beta_k^2}{2 \tau^2} \right) \right] d \beta_k 
\end{aligned}
$$

Let $x_k^2 = \frac{\beta_k^2}{2 \tau^2}$, which implies $\beta_k = \sqrt{2 \tau^2 x_k^2}$. We can then define:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}^k; \alpha, \tau^2 \rvert \beta_k)
&= \prod_{i = 1}^n \frac{(\mu_i^k)^{y_i^k} \exp(- \mu_i^k)}{y_i^k!} \\
&= \prod_{i = 1}^n \frac{(\exp(\alpha + \beta_k))^{y_i^k} \exp\left( - \exp(\alpha + \beta_k)\right)}{y_i^k!} \\
&= \prod_{i = 1}^n \frac{(\exp(\alpha + \sqrt{2 \tau^2 x_k^2}))^{y_i^k} \exp\left( - \exp(\alpha + \sqrt{2 \tau^2 x_k^2})\right)}{y_i^k!} \\
&= \frac{1}{\prod_{i = 1}^n y_i^k!} \left[ \exp(\alpha + \sqrt{2 \tau^2 x_k^2})\right]^{\sum_{i = 1}^n y_i^k} \left[ \exp\left( - \exp( \alpha + \sqrt{2 \tau^2 x_k^2}) \right) \right]^n
\end{aligned}
$$

which we can call $f(x_k)$. Then, via $u$-substitution, we have:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}^k; \alpha, \tau^2) 
&= \frac{1}{\sqrt{2 \pi \tau^2}} \int f(x_k) \exp(- x_k^2) \left(\sqrt{2 \tau^2} \right) dx_k \\
&= \frac{1}{\sqrt{\pi}} \int f(x_k) \exp(- x_k^2) dx_k
\end{aligned}
$$

This is the form we need for Gauss-Hermite quadrature, and the order $m$ approximation is:

$$
\begin{aligned}
\mathcal{L}(\mathbf{y}; \alpha, \tau^2) 
&= \prod_{k = 1}^K \mathcal{L}(\mathbf{y}^k; \alpha, \tau^2) \\
&= \prod_{k = 1}^K \left(\frac{1}{\sqrt{\pi}} \right) \int f(x_k) \exp(- x_k^2) dx_k \\
&\approx \frac{1}{\pi^{\frac{K}{2}}} \prod_{k = 1}^K \left( \sum_{i = 1}^m \frac{2^{m - 1} m! \sqrt{\pi}}{m^2 \left[ H_{m - 1} (\bar{x}_i) \right]^2} f(\bar{x}_i) \right) \\
\implies
\ell(\mathbf{y}; \alpha, \tau^2) 
&= \sum_{k = 1}^K \left[ -\frac{1}{K} \log(\pi) + \log\left( \int f(x_k) \exp(- x_k^2) dx_k\right) \right] \\
&\approx - \frac{K}{2} \log(\pi) + \sum_{k = 1}^K \log\left( \sum_{i = 1}^m \frac{2^{m - 1} m! \sqrt{\pi}}{m^2 \left[ H_{m - 1} (\bar{x}_i) \right]^2} f(\bar{x}_i) \right)
\end{aligned}
$$

where $m$ is the desired polynomial degree for the approximation and $\bar{x}_i$ are the roots of the $m$-th Hermite polynomial.



<!-- We can also approximate the derivatives of the log-likelihood with respect to any of the model parameters (which we denote with $\gamma$). Let $f_k(\bar{x}_i)$ denote the function, $f$, defined above for the $k$-th random intercept evaluated at the $i$-th root of the $m$-th Hermite polynomial:

$$
f_k(\bar{x}_i) = \frac{1}{\prod_{j = 1}^n y_j^k!} \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k} \left[ \exp\left(- \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \right) \right]^n
$$

Noticing that the weights are all fixed, we have:

$$
\begin{aligned}
\ell(\mathbf{y}; \alpha, \tau^2) 
&= \sum_{k = 1}^K \left[ -\frac{1}{K} \log(\pi) + \log\left( \int f(x_k) \exp(- x_k^2) dx_k\right) \right] \\
\implies
\frac{\partial \ell(\mathbf{y}; \alpha, \tau^2)}{\partial \gamma}
&= \sum_{k = 1}^K \left[\frac{\frac{\partial}{\partial \gamma} \left[ \int f(x_k) \exp(- x_k^2) dx_k \right]}{\int f(x_k) \exp(- x_k^2) dx_k} \right] \\ 
&\approx \sum_{k = 1}^K \left[ \frac{\frac{\partial}{\partial \gamma} \left[ \sum_{i = 1}^m \omega(\bar{x}_i) f_k(\bar{x}_i) \right]}{\sum_{i = 1}^m \omega(\bar{x}_i) f_k(\bar{x}_i)} \right] \\
&= \sum_{k = 1}^K \left[ \frac{ \sum_{i = 1}^m \omega(\bar{x}_i) \frac{\partial}{\partial \gamma} \left[ f_k(\bar{x}_i) \right]}{\sum_{i = 1}^m \omega(\bar{x}_i) f_k(\bar{x}_i)} \right]
\end{aligned}
$$

Taking the derivative of $f_k(\bar{x}_i)$ with respect to $\alpha$ and $\tau^2$ yields:

$$
\begin{aligned}
\frac{\partial f_k(\bar{x}_i)}{\partial \tau^2} 
&= \frac{\partial}{\partial \tau^2} \frac{1}{\prod_{j = 1}^n y_j^k!} \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k} \left[ \exp\left(- \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \right) \right]^n \\
&= \frac{1}{\prod_{j = 1}^n y_j^k!} \left( \left[ \exp\left(- \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \right) \right]^n \frac{\partial}{\partial \tau^2} \left[  \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k} \right] +  \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k} \frac{\partial}{\partial \tau^2} \left[ \left[ \exp\left(- \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \right) \right]^n\right] \right) \\
&= \frac{1}{\prod_{j = 1}^n y_j^k!} \left( \frac{\left(\sum_{j = 1}^n y_j^k\right)\left[\exp\left( \alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k - 1}\exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \left(\frac{1}{2}(2 \bar{x}_i^2)(2 \tau^2 \bar{x}_i^2)^{- \frac{1}{2}} \right)}{\left[ \exp\left(\exp\left( \alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right)\right)\right]^n}+ \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k} \left[ - \frac{n \bar{x}_i^2 \exp\left( \alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \left(\exp\left(- \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right)\right)^n}{\sqrt{2 \tau^2 \bar{x}_i^2}} \right] \right) \\
&= \frac{1}{\prod_{j = 1}^n y_j^k!} \left( \frac{\bar{x}_i^2 \left(\sum_{j = 1}^n y_j^k\right)\left[\exp\left( \alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k}}{\sqrt{2 \tau^2 \bar{x}_i^2} \left[ \exp\left(\exp\left( \alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right)\right)\right]^n}+ \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k} \left[ - \frac{n \bar{x}_i^2 \exp\left( \alpha + \sqrt{2 \tau^2 \bar{x}_i^2}\right) \left(\exp\left(- \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right)\right)^n}{\sqrt{2 \tau^2 \bar{x}_i^2}} \right] \right) \\
&= \frac{\bar{x}_i^2 \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k}}{\sqrt{2 \tau^2 \bar{x}_i^2} \left(\prod_{j = 1}^n y_j^k!\right)\left[ \exp\left(\exp \left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right)\right)\right]^n} \left[ \sum_{j = 1}^n y_j^k - n \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right)\right] \\
&= \frac{\bar{x}_i^2 \left[ \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right]^{\sum_{j = 1}^n y_j^k}}{\sqrt{2 \tau^2 \bar{x}_i^2} \left(\prod_{j = 1}^n y_j^k!\right)\left[ \exp\left(\exp \left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right)\right)\right]^n} \left[ \sum_{j = 1}^n\left(y_j^k - \exp\left(\alpha + \sqrt{2 \tau^2 \bar{x}_i^2} \right) \right)\right]
\end{aligned}
$$ -->


