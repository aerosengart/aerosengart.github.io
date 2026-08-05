---
layout: distill
title: Relative Curvature Measures of Nonlinearity
description:
date: 2026-03-25
tabs: true
tags: glmm inference models
toc:
    - name: Background
      subsections:
        - name: Set-Up
        - name: Solution Locus
        - name: Measures of Non-Linearity
    - name: Geometry
      subsections:
        - name: Differential Geometry
        - name: Physical Interpretation
bibliography: 2026-03-25-rel-curvature.bib
---

## Background

### Set-Up
For each of $n$ observations, we have the response, $y_i$, and some covariate vector, $\mathbf{x}_i$. We assume that the responses follow a standard non-linear model parametrized by $\boldsymbol{\theta} = (\theta_1, \dots, \theta_p)^\top$. That is:

$$
\begin{aligned}
y_i &= f(\mathbf{x}_i, \boldsymbol{\theta}) + \epsilon_i \\
\mathbb{E}\left[ \epsilon_i \right] &= 0 \\
\text{Cov}(\epsilon_i, \epsilon_{i'}) &= \sigma^2 \delta_{i,i'}
\end{aligned}
$$

The function $f(\mathbf{x}_i, \boldsymbol{\theta})$ is called the <strong>model function</strong> and represents the mean of $y$ conditional on the covariates, $\mathbf{x}$, and the parameter $\boldsymbol{\theta}$. For fixed $\mathbf{x}_i$, if we define:

$$
\begin{aligned}
\eta_i(\boldsymbol{\theta}) &= \mathbb{E}\left[ y_i \rvert \boldsymbol{\theta} \right] = f(\mathbf{x}_i, \boldsymbol{\theta}) \\
\boldsymbol{\eta}(\boldsymbol{\theta}) &= (\eta_1(\boldsymbol{\theta}), \dots, \eta_n(\boldsymbol{\theta}))^\top
\end{aligned}
$$

then we can write the <strong>least squares estimate</strong> of $\boldsymbol{\theta}$ as:

$$
\hat{\boldsymbol{\theta}} = \underset{\boldsymbol{\theta}}{\arg \min} \left\{ S(\boldsymbol{\theta})\right\} = \underset{\boldsymbol{\theta}}{\arg \min} \left\{ \rvert \rvert \mathbf{y} - \boldsymbol{\eta}(\boldsymbol{\theta})\rvert \rvert^2 \right\}
$$

### Solution Locus
Given a point, $\boldsymbol{\theta}_0$, in the parameter space, $\Theta$, and fixing the covariates, we can map $\boldsymbol{\theta}_0$ to some point $P(\boldsymbol{\theta}_0)$ in the $n$-dimensional sample space, $\mathcal{D}$, whose $i$-th coordinate is determined by $\eta_i(\boldsymbol{\theta}) = f(\mathbf{x}_i, \boldsymbol{\theta}_0)$. Consider the sets:

$$
\begin{aligned}
\Theta_S &= \left\{ \boldsymbol{\theta} \in \Theta \rvert S(\boldsymbol{\theta}) = 0 \right\} \\
\mathbf{P}_S &= \left\{ P(\boldsymbol{\theta}) \rvert \boldsymbol{\theta} \in \Theta_S \right\}
\end{aligned}
$$

The set $\mathbf{P}_S$ is termed the <strong>solution locus</strong>, and $S(\boldsymbol{\theta})$ is the squared distance between the responses and the conditional expected responses for a particular value of $\boldsymbol{\theta}$. It is unlikely that a perfect fit can be achieved; in other words, $\mathbf{y}$ will usually not lie in the solution locus. Instead, we find $\hat{\boldsymbol{\theta}}$ such that our conditional expected responses, $\boldsymbol{\eta}(\hat{\boldsymbol{\theta}})$, is closest to $\mathbf{y}$ (in terms of this squared distance).

Because $f(\mathbf{x}, \boldsymbol{\theta})$ may be non-linear in $\boldsymbol{\theta}$, we often perform a local linear approximation of the model function about some chosen value $\boldsymbol{\theta}_0$:

$$
\begin{aligned}
f(\mathbf{x}, \boldsymbol{\theta}) 
&\approx f(\mathbf{x}, \boldsymbol{\theta}_0) + (\boldsymbol{\theta} - \boldsymbol{\theta}_0)^\top \left. \left[ \frac{\partial f(\mathbf{x}, \boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \right]\right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \\
&= \boldsymbol{\eta}(\boldsymbol{\theta}_0) + (\boldsymbol{\theta} - \boldsymbol{\theta}_0)^\top \left. \left[ \frac{\partial \boldsymbol{\eta}(\boldsymbol{\theta})}{\partial \boldsymbol{\theta}} \right]\right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \\
\end{aligned}
$$

This replaces the solution locus with a plane that is tangent to the locus at the point $\boldsymbol{\eta}(\boldsymbol{\theta}_0)$ (called the <strong>planar assumption</strong>). This also places a uniform coordinate system onto that tangent plane (called the <strong>uniform coordinate assumption</strong>). 

<aside><p>I am not entirely sure, but I think they mean Euclidean coordinate system?</p></aside>

### Measures of Non-Linearity
<strong>Measures of non-linearity</strong> are quantities that convey how close the above linear approximation is to the model and how that relationship affects downstream inferences. Several were proposed in Beale (1960)<d-cite key=beale1960></d-cite>, and Bates and Watts provide a general overview of a few (which we'll omit here). New measures of non-linearity in the setting specified above are the subject of this paper.


---

## Geometry
The new measures of non-linearity proposed by Bates and Watts are developed from a geometric perspective. We re-iterate the relevant concepts here.

Recall that the $p$-dimensional parameter space is denoted by $\Theta$, and the $n$-dimensional sample space is $\mathcal{D}$. For fixed $\mathbf{x}$, the model function, $f(\mathbf{x}, \boldsymbol{\theta})$, is just a function of the parameter vector and maps between $\Theta$ and $\mathcal{D}$. 

### Differential Geometry
Consider some straight line in $\Theta$ through some arbitrary point $\boldsymbol{\theta}_0 \in \Theta$:

$$
\boldsymbol{\theta}(b) = \boldsymbol{\theta}_0 + b \mathbf{h}
$$

where $\mathbf{h}$ is any non-zero, $p$-dimensional vector. This line maps to the curve (also called a <strong>lifted line</strong> since we bring it up the line in $p$ dimensions to $n$ dimensions):

$$
\boldsymbol{\eta}_\mathbf{h}(b) = \boldsymbol{\eta}(\boldsymbol{\theta}_0 + b \mathbf{h}) \in \mathcal{D}
$$

The tangent vector to $\boldsymbol{\eta}_\mathbf{h}(b)$ at $b = 0$ is given by:

$$
\begin{equation}
\label{eq:velocity}
\dot{\boldsymbol{\eta}}_\mathbf{h} = \left. \left[ \frac{d \boldsymbol{\eta}_\mathbf{h}(b)}{d b} \right] \right\rvert_{b = 0} = \sum_{i = 1}^p \left. \left[ \frac{\partial \boldsymbol{\eta}(\boldsymbol{\theta})}{\partial \theta_i} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} \left. \left[ \frac{d \theta_i}{d b} \right] \right\rvert_{b = 0} = \sum_{i = 1}^p \left. \left[ \frac{\partial \boldsymbol{\eta}(\boldsymbol{\theta})}{\partial \theta_i} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} h_i
\end{equation}
$$

Thus, this tangent vector is just a linear combination of the first-order partial derivatives of the conditional expected response vectors with respect to the components of $\boldsymbol{\theta}$ evaluated at our choice of $\boldsymbol{\theta}_0$. Note that the <strong>tangent plane</strong> at $\boldsymbol{\eta}(\boldsymbol{\theta}_0)$ is given by all possible such linear combinations. 

Consider, now, the second derivative of $\boldsymbol{\eta}_\mathbf{h}(b)$ at $b = 0$. 

$$
\begin{equation}
\label{eq:acceleration}
\begin{aligned}
\ddot{\boldsymbol{\eta}}_\mathbf{h}
&= \left. \left[ \frac{d^2 \boldsymbol{\eta}_\mathbf{h}(b)}{d b^2} \right] \right\rvert_{b = 0} \\
&= \sum_{j = 1}^p \left. \left[ \frac{\partial}{\partial \theta_j} \left[ \frac{d \boldsymbol{\eta}_{\mathbf{h}}(b)}{d b} \right] \right] \right\rvert_{b = 0} \left. \left[ \frac{d \theta_j}{d b} \right] \right\rvert_{b = 0} \\
&= \sum_{j = 1}^p \left( \sum_{i = 1}^p \left. \left[ \frac{\partial \boldsymbol{\eta}(\boldsymbol{\theta})}{\partial \theta_i} \right]\right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0}  \left. \left[ \frac{d \theta_i}{d b} \right] \right\rvert_{b = 0} \right) \left. \left[ \frac{d \theta_j}{d b} \right] \right\rvert_{b = 0} \\
&= \sum_{j = 1}^p  \sum_{i = 1}^p \left. \left[ \frac{\partial \boldsymbol{\eta}(\boldsymbol{\theta})}{\partial \theta_i} \right] \right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}_0} h_i h_j \\ 
&= \mathbf{h}^\top  \left. \left[ \frac{\partial \boldsymbol{\eta}(\boldsymbol{\theta})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right]\right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}} \mathbf{h}
\end{aligned}
\end{equation}
$$

where $$\ddot{\boldsymbol{\eta}}_\mathbf{h}$$ is a $p \times p \times n$ array (because $\boldsymbol{\theta}$ is $p$-dimensional and $\boldsymbol{\eta}(\boldsymbol{\theta})$ is $n$-dimensional). The quadratic form in the last line above evaluates to an $n$-dimensional vector whose $i$-th element is given by $$\mathbf{h}^\top  \left. \left[ \frac{\partial \eta_i(\boldsymbol{\theta})}{\partial \boldsymbol{\theta} \partial \boldsymbol{\theta}^\top} \right]\right\rvert_{\boldsymbol{\theta} = \boldsymbol{\theta}} \mathbf{h}$$. 

### Physical Interpretation
Imagine, now, that $b$ represents time (or something similar), and we follow a particle moving through the sample space such that at time $b$ the particle is at $$\boldsymbol{\eta}_\mathbf{h}(b)$$ (this is our position). The first and second derivatives of our position with respect to time give us our <i>velocity</i> and <i>acceleration</i>. Thus, $$\dot{\boldsymbol{\eta}}_\mathbf{h}$$ and $$\ddot{\boldsymbol{\eta}}_\mathbf{h}$$ in Eq. \eqref{eq:velocity} and \eqref{eq:acceleration} give our instantaneous velocity and acceleration, respectively, when $b = 0$. 

Since acceleration is itself a vector, we can decompose it into two components: the <strong>tangential</strong> and <strong>normal</strong> acceleration.

<!-- #region accel -->
{% tabs accel %}
{% tab accel statement %}
Let $$\bar{\boldsymbol{\eta}}_\mathbf{h}(s)$$ be the arc length parametrization of $$\boldsymbol{\eta}_\mathbf{h}(b)$$; that is, for a given $s$, we let $b = b(s)$, the inverse of the arc length function.  

$$
\ddot{\boldsymbol{\eta}}_{\mathbf{h}} = \ddot{\boldsymbol{\eta}}_\mathbf{h}^T + \ddot{\boldsymbol{\eta}}_\mathbf{h}^N 
$$


Recall that the velocity of our particle is given by:

$$
\begin{aligned}
\frac{\mathbf{T}'(t)}{\rvert \rvert \mathbf{T}'(t) \rvert \rvert}
\end{aligned}
$$

Thus:

$$
\begin{aligned}
\ddot{\boldsymbol{\eta}}_{\mathbf{h}}(b)
&= \frac{d}{d b} \left[ \ddot{\boldsymbol{\eta}}_\mathbf{h}(b) \right] \\
&= \frac{d}{d b} \left[ \mathbf{T}(s) \frac{ds}{db} \right] \\ 
&= \frac{d}{d b} \left[ \mathbf{T}(s) \right] \frac{ds}{db} + \mathbf{T}(s) \frac{d^2 s}{db^2} \\ 
&= \frac{d \mathbf{T}(s)}{ds} \left(\frac{ds}{db}\right)^2 + \left( \frac{d^2 s}{db^2} \right) \mathbf{T}(s)
\end{aligned}
$$

The arc length function is defined as:

$$
s(b) = \int_a^b \rvert \rvert \dot{\boldsymbol{\eta}}_\mathbf{h}(x) \rvert \rvert dx 
$$

By the fundamental theorem of calculus:

$$
\begin{aligned}
\frac{ds}{db} &=  \rvert \rvert \dot{\boldsymbol{\eta}}_\mathbf{h}(b) \rvert \rvert \\
\frac{d^2s}{db^2} &= \frac{d}{db} \left[ \rvert \rvert \dot{\boldsymbol{\eta}}_{\mathbf{h}}(b) \rvert \rvert \right]
\end{aligned}
$$

Recall also that:

$$
\mathbf{N}(s) = \frac{d \mathbf{T}(s)}{ds}
$$

Using what we found in our tangent on vector calculus:

$$
\begin{aligned}
\ddot{\boldsymbol{\eta}}_{\mathbf{h}}(b)
&= \frac{d \mathbf{T}(s)}{ds} \left(\frac{ds}{db}\right)^2 + \left( \frac{d^2 s}{db^2} \right) \mathbf{T}(s) \\
&= \mathbf{T}'(s) \left(\frac{ds}{db}\right)^2 +  \left( \frac{d^2 s}{db^2} \right)  \frac{\mathbf{T}(t)}{\rvert \rvert \mathbf{T}(t) \rvert \rvert}\\ 
&= \kappa(s) \left(\frac{ds}{db}\right)^2 \mathbf{N}(s) \left( \frac{d^2 s}{db^2} \right)  \frac{\mathbf{T}(t)}{\rvert \rvert \mathbf{T}(t) \rvert \rvert} \\
&= a_\mathbf{T}()
\end{aligned}
$$
{% endtab %}
{% tab accel proof %}
Proof?
{% endtab %}
{% endtabs %}
<!-- #endregion -->

These components have special interpretations because the acceleration is the derivative of the velocity vector. The tangential acceleration is the component that is parallel to the velocity vector at the given point, and it can be interpreted as the rate of change in the <i>magnitude</i> of the velocity vector. The normal acceleration is the component that is orthogonal to the velocity vector at the given point, and it can be interpreted as the rate of change in the <i>direction</i> of the velocity vector. 

As we noted above, the velocity vector is the same as the tangent vector, which implies that the normal acceleration is orthogonal to the tangent plane at $\boldsymbol{\eta}(\boldsymbol{\theta}_0)$. We can therefore interpret the normal acceleration as describing the change in <i>direction</i> in the velocity that occurs normal to the tangent plane. 

We can further decompose the tangent acceleration with respect to the tangent plane as:

$$
\ddot{\boldsymbol{\eta}}_\mathbf{h}^T = \ddot{\boldsymbol{\eta}}_\mathbf{h}^P + \ddot{\boldsymbol{\eta}}_\mathbf{h}^G
$$

<aside><p>Proof?</p></aside>


These vectors are the tangent and normal components of the tangent acceleration and represent its change in direction that are parallel and normal to the tangent plane at $\boldsymbol{\eta}(\boldsymbol{\theta}_0)$, respectively.

As such, $$\ddot{\boldsymbol{\eta}}_\mathbf{h}^P$$, which we call the <strong>parallel tangent component</strong>, describes the change in <i>speed</i> of the particle. It is the component that lies parallel to $$\dot{\boldsymbol{\eta}}_\mathbf{h}$$; that is, it is parallel to the tangent plane. The component $$\ddot{\boldsymbol{\eta}}_\mathbf{h}^G$$, which we call the <strong>geodesic acceleration</strong>, is normal to $$\dot{\boldsymbol{\eta}}_\mathbf{h}$$ and represents the change in direction of the acceleration that is parallel to the tangent plane lying normal to $$\dot{\boldsymbol{\eta}}(\boldsymbol{\theta}_0)$$.

We can also define the <strong>curvature</strong>, which is the magnitude of the acceleration. Each of the components of our decomposed acceleration vector can be converted into a curvature:

$$
\begin{aligned}
\kappa_\mathbf{h} 
&= \rvert \rvert \ddot{\boldsymbol{\eta}}_\mathbf{h} \rvert \rvert \\
&= \rvert \rvert \ddot{\boldsymbol{\eta}}_\mathbf{h}^N + \ddot{\boldsymbol{\eta}}_\mathbf{h}^T \rvert \rvert \\
&= 
\end{aligned}
$$



---






---







