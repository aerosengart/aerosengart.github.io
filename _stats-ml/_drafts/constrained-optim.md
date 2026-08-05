---
layout: distill
title: Constrained Optimization
description: A Primer
date: 2026-03-05
tabs: true
tags: theory optimization primer
# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
    - name: Background
    - name: Equality Constraints
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
bibliography: stats-ml.bib
---

In many machine learning methods, we have to optimize some objective function subject to constraints. In this post, I'll go over some of the basics concepts behind solving such problems. 

---

## Background
We assume to have some scalar- or vector-valued variable (parameter), $\theta$, and a loss function, $\mathcal{L}(\theta)$. Our goal is to find:

$$
\begin{equation}
\label{eq:const-optim}
\begin{aligned}
\theta^* &= \underset{\theta \in \mathcal{C}}{\arg \min} \left\{ \mathcal{L}(\theta) \right\} \\
\mathcal{C} &= \left\{ \theta \in \mathbb{R}^p \hspace{1mm} \rvert \hspace{1mm} h_i(\theta) = 0, i \in \mathcal{E} \text{; } g_j(\theta) \leq 0, j \in \mathcal{I} \right\}
\end{aligned}
\end{equation}
$$

<aside><p>We write $h_i(\theta) = 0$ and $g_j(\theta) \leq 0$, but with a simple shifting, we could generalize to $h_i(\theta) = c_i$ and $g_j(\theta) \leq d_i$.</p></aside>

The set, $\mathcal{C}$, is our <strong>feasible</strong> or <strong>constraint</strong> set, and it dictates the possible valid values of $\theta$. The set, $\mathcal{E}$, are our <strong>equality constraints</strong>, and the set, $\mathcal{I}$, is the set of <strong>inequality constraints</strong>.



---

## Equality Constraints
Suppose we only have a single equality constraint:

$$
\begin{equation}
\label{eq:const-optim-one-eq}
\begin{aligned}
\theta^* &= \underset{\theta \in \mathcal{C}}{\arg \min} \left\{ \mathcal{L}(\theta) \right\} \\
\mathcal{C} &= \left\{ \theta \in \mathbb{R}^p \hspace{1mm} \rvert \hspace{1mm} h(\theta) = 0 \right\}
\end{aligned}
\end{equation}
$$

We have the following result that will provide some foundation for later theory.

<!-- #region constraint-surface-->
<div class="theorem">
<strong>Claim.<d-cite key=murphy2025></d-cite></strong>
{% tabs constraint-surface %}
{% tab constraint-surface statement %}
In the setting described in Eq. \eqref{eq:const-optim} with a single equality constraint, $\frac{\partial h(\theta)}{\partial \theta}$, the gradient of the equality constraint with respect to $\theta$, at any point on the constraint surface will be orthogonal to the surface.
{% endtab %}
{% tab constraint-surface proof %}
Consider a first-order Taylor expansion of the equality constraint about $\theta$ on the constraint surface. Let $\theta + \epsilon$ be some point on the constraint surface that is close to $\theta$. We have:

$$
\begin{aligned}
h(\theta) &\approx h(\theta + \epsilon) + (\theta - (\theta + \epsilon))^\top \frac{\partial h(\theta)}{\partial \theta} \\
\implies h(\theta + \epsilon) &\approx h(\theta) + \epsilon^\top \frac{\partial h(\theta)}{\partial \theta}
\end{aligned}
$$

Since both $\theta$ and $\theta + \epsilon$ lie on the constraint surface, it must be that $h(\theta) = h(\theta + \epsilon)$. Thus:

$$
\epsilon^\top \frac{\partial h(\theta)}{\partial \theta} \approx 0
$$

which implies that $\epsilon$ and $\frac{\partial h(\theta)}{\partial \theta}$ are orthogonal. For $\theta + \epsilon$ <i>and</i> $\theta$ to lie on the constraint surface, $\epsilon$ must be parallel to the surface. Thus, $\frac{\partial h(\theta)}{\partial \theta}$ is perpendicular to the surface.
{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Since it lies on the constraint surface, the solution, $$\theta^*$$, in Eq. \eqref{eq:const-optim-one-eq} will also satisfy $$\left. \frac{\partial h(\theta)}{\partial \theta} \right\rvert_{\theta = \theta^*} = 0$$. Furthermore, since $\theta^*$ minimizes the loss function, the following must hold.

<!-- #region constraint-surface-2 -->
<div class="theorem">
<strong>Claim.<d-cite key=murphy2025></d-cite></strong>
{% tabs constraint-surface-2 %}
{% tab constraint-surface-2 statement %}
In the setting described in Eq. \eqref{eq:const-optim} with a single equality constraint, $\frac{\partial \mathcal{L}(\theta)}{\partial \theta}$, the gradient of the loss with respect to $\theta$, at any point on the constraint surface will be orthogonal to the surface.
{% endtab %}
{% tab constraint-surface-2 proof %}
Consider a first-order Taylor expansion of the loss about $\theta$ on the constraint surface. Let $\theta + \epsilon$ be some point on the constraint surface that is close to $\theta$. We have:

$$
\begin{aligned}
\mathcal{L}(\theta) &\approx \mathcal{L}(\theta + \epsilon) + (\theta - (\theta + \epsilon))^\top \frac{\partial \mathcal{L}(\theta)}{\partial \theta} \\
\implies \mathcal{L}(\theta + \epsilon) &\approx \mathcal{L}(\theta) + \epsilon^\top \frac{\partial \mathcal{L}(\theta)}{\partial \theta}
\end{aligned}
$$

Now, consider $\theta^*$, the minimizer of Eq. \eqref{eq:const-optim-one-eq}:

$$
\mathcal{L}(\theta^* + \epsilon) \approx \mathcal{L}(\theta^*) + \epsilon^\top \left.\frac{\partial \mathcal{L}(\theta)}{\partial \theta} \right\rvert_{\theta = \theta^*}
$$

Since $\theta^*$ is the minimizer, it must satisfy:

$$
\mathcal{L}(\theta^*) \leq \mathcal{L}(\theta)
$$

for any other $\theta$ that lies on the constraint surface. 


{% endtab %}
{% endtabs %}
</div>
<!-- #endregion -->

Since the two gradients must be orthogonal to the constraint surface at $$\theta^*$$, they must be parallel or anti-parallel to each other. This implies that there exists some constant, $$\lambda^* \in \mathbb{R}$$, such that:

$$
\left.\frac{\partial \mathcal{L}(\theta)}{\partial \theta} \right\rvert_{\theta=\theta^*} = \lambda^* \left. \frac{\partial h(\theta)}{\partial \theta}\right\rvert_{\theta = \theta^*}
$$

This constant, $$\lambda^*$$, is called the <strong>Lagrange multiplier</strong>. We can construct an alternative objective function, called the <strong>Lagrangian</strong>:

$$
\begin{equation}
\label{eq:lagrangian-1}
\mathcal{L}(\theta, \lambda) = \mathcal{L}(\theta) + \lambda h(\theta)
\end{equation}
$$

We can then find solutions by identifying stationary points, $$(\theta^*, \lambda^*)$$, of the Lagrangian because:

$$
\begin{aligned}
\left. \frac{\partial \mathcal{L}(\theta, \lambda)}{\partial \theta} \right\rvert_{\theta = \theta^*, \lambda = \lambda^*} = 0
&\text{;     } \left. \frac{\partial \mathcal{L}(\theta, \lambda)}{\partial \theta} \right\rvert_{\theta = \theta^*, \lambda = \lambda^*} = 0 \\
&\iff \\
\lambda^* \left.\frac{\partial h(\theta)}{\partial \theta} \right\rvert_{\theta = \theta^* } = \left.\frac{\partial \mathcal{L}(\theta)}{\partial \theta}\right\rvert_{\theta = \theta^*}
&\text{;     } 
h(\theta^*) = 0
\end{aligned}
$$

<aside><p>We call $(\theta^*, \lambda^*)$ a <strong>critical point</strong>.</p></aside>

### Multiple Constraints
Multiple constrains can be accommodated by letting $\lambda \in \mathhbb{R}^m$ (one for each constraint):

$$
\mathcal{L}(\theta, \lambda) = \mathcal{L}(\theta) + \sum_{j = 1}^m \lambda_j h_j(\theta)=
$$

---

## Inequality Constraints


