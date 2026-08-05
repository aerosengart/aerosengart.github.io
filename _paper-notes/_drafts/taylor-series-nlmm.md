---
layout: distill
title: Two Taylor-Series Approximation Methods for Non-Linear Mixed Models
description:
date: 2026-03-20
tabs: true
tags: glmm inference models
toc:
    - name: Non-Linear Mixed Model
    - name: Zero-Expansion Method
bibliography: 2026-03-20-taylor-series-nlmm.bib
---

Wolfinger and Lin propose two methods for approximate inference in non-linear mixed models in their work <i>Two Taylor series approximation methods for nonlinear mixed models</i>.<d-cite key=wolfinger1997></d-cite> These methods are based upon Taylor series approximations of the integrated likelihood and differ in the point of expansion. One expands about the random effects' mean and the other expands about the <i>empirical best linear unbiased predictors (EBLUPs)</i> for the random effects. 

---

## Non-Linear Mixed Model
We let $\boldsymbol{\alpha}$ be a $p$-dimensional vector of fixed effects and $\boldsymbol{\beta}$ a $q$-dimensional vector of random effects. We assume to have the $n \times r$ and $n \times s$ matrices of covariates, $\mathbf{X}$ and $\mathbf{Z}$, associated with the fixed and random effects, respectively. For an $n$-dimensional response vector, $\mathbf{y}$, a <i>non-linear mixed model</i> is given by:

$$
\begin{equation}
\label{eq:gen-model}
\begin{aligned}
\mathbf{y} &= f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) + \mathbf{e} \\
\mathbf{e} &\sim \mathcal{N}(\mathbf{0}, \mathbf{R}(\boldsymbol{\theta})) \\
\boldsymbol{\beta} &\sim \mathcal{N}(\mathbf{0}, \mathbf{D}(\boldsymbol{\theta}))
\end{aligned}
\end{equation}
$$

for some covariance parameter vector, $\boldsymbol{\theta}$. We assume that the random effects, $\boldsymbol{\beta}$, and the errors, $\mathbf{e}$, are independent. 

The model in Eq. \eqref{eq:gen-model} encapsulates many more common forms such as an additive random effects model:

$$
\mathbf{y} = f(\mathbf{X}, \boldsymbol{\alpha}) + \mathbf{Z} \boldsymbol{\beta} + \mathbf{e}
$$

and a special case of <a href="/stats-ml/glmm/">generalized linear model</a> (if $f(\cdot)$ is invertible). 

The <i>integrated likelihood</i> for Eq. \eqref{eq:gen-model} removes the random effects by integrating the joint likelihood of the responses and random effects:

$$
\begin{equation}
\label{eq:integ-lik}
\begin{aligned}
\mathcal{L}(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta}) &= \int \exp\left(\ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right) \exp\left(- \ell(\boldsymbol{\beta}) \right) d \boldsymbol{\beta} \\
\ell(\boldsymbol{\beta}) &= - \frac{q}{2} \log(2 \pi) - \frac{1}{2} \log\left( \rvert \mathbf{D}(\boldsymbol{\theta}) \rvert\right) - \frac{1}{2} \boldsymbol{\beta}^\top \mathbf{D}^{-1}(\boldsymbol{\theta}) \boldsymbol{\beta} \\
\ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) &= - \frac{n}{2} \log(2 \pi) - \frac{1}{2} \log\left( \rvert \mathbf{R}(\boldsymbol{\theta}) \rvert\right) - \frac{1}{2} \left(\mathbf{y} - f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})\right)^\top \mathbf{R}^{-1}(\boldsymbol{\theta}) \left(\mathbf{y} - f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})\right) 
\end{aligned}
\end{equation}
$$

where we have suppressed the covariates from the notation.

Because the first line in Eq. \eqref{eq:integ-lik} may be a high-dimensional integral, it may be intractable and require approximation.

---

## Zero-Expansion Method
The <i>zero-expansion method</i> of Wolfinger and Lin is a second-order Taylor approximation of the conditional log-likelihood of the responses on the random effects about the random effects mean:

$$
\begin{aligned}
\left. \left[ \frac{\partial}{\partial \boldsymbol{\beta}} \left[ \ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right] \right]  \right\rvert_{\boldsymbol{\beta} = \mathbf{0}}  
&= - \frac{1}{2} \left(- \left(\mathbf{y} - f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \mathbf{0}) \right)^\top \mathbf{R}^{-1}(\boldsymbol{\theta})  \left. \left[ \frac{\partial f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}}  - \left(\mathbf{y} - f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \mathbf{0}) \right)^\top \mathbf{R}^{-1}(\boldsymbol{\theta})  \left. \left[ \frac{\partial f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}}\right) \\
&= \left(\mathbf{y} - f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \mathbf{0}) \right)^\top \mathbf{R}^{-1}(\boldsymbol{\theta}) \left. \left[ \frac{\partial f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}} 
\end{aligned}
$$

$$
\begin{aligned}
\frac{\partial^2}{\partial \beta_i \partial \beta_{i'}} \left[ \ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right]
&= \frac{\partial^2}{\partial \beta_i \partial \beta_{i'}} \left[ - \frac{1}{2} \sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} (\mathbf{y}_j - f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}))  (\mathbf{y}_{j'} - f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})) \right] \\
&= - \frac{1}{2} \sum_{j = 1}^n \sum_{j' = 1}^n  \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'}  \frac{\partial}{\partial \beta_i} \left[ (\mathbf{y}_j - f_j(\mathbf{X}, \boldsymbol{\alpha},\mathbf{Z}, \boldsymbol{\beta}))\left( \frac{\partial}{\partial \beta_{i'}} \left[ -f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) + (\mathbf{y}_{j'} - f_{j'}(\mathbf{X}, \boldsymbol{\alpha},\mathbf{Z}, \boldsymbol{\beta}))\left( \frac{\partial}{\partial \beta_{i'}} \left[ -f_{j}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) \right] \\
&= - \frac{1}{2} \sum_{j = 1}^n \sum_{j' = 1}^n  \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} \left( (\mathbf{y}_j - f_j(\mathbf{X}, \boldsymbol{\alpha},\mathbf{Z}, \boldsymbol{\beta}))\left( \frac{\partial^2}{\partial \beta_i \partial \beta_{i'}} \left[ -f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) 
    + \left( \frac{\partial}{\partial \beta_{i}} \left[ -f_{j}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) \left( \frac{\partial}{\partial \beta_{i'}} \left[ -f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) 
    + (\mathbf{y}_{j'} - f_{j'}(\mathbf{X}, \boldsymbol{\alpha},\mathbf{Z}, \boldsymbol{\beta}))\left( \frac{\partial^2}{\partial \beta_i \partial \beta_{i'}} \left[ -f_{j}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) 
    + \left( \frac{\partial}{\partial \beta_i} \left[ -f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right)\left( \frac{\partial}{\partial \beta_{i'}} \left[ -f_{j}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}) \right] \right) \right) \\
&= - \frac{1}{2} \left[ -\sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} (\mathbf{y}_j - f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})) \frac{\partial^2 f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i \partial \beta_{i'}} 
+ \sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} \frac{\partial f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i} \frac{\partial f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_{i'}}
\right]
- \sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} (\mathbf{y}_{j'} - f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})) \frac{\partial^2 f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i \partial \beta_{i'}} 
+ \sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} \frac{\partial f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i} \frac{\partial f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_{i'}} \\
&= \sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} (\mathbf{y}_j - f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})) \frac{\partial^2 f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i \partial \beta_{i'}} - \sum_{j = 1}^n \sum_{j' = 1}^n \left[ \mathbf{R}^{-1}(\boldsymbol{\theta}) \right]_{j,j'} \frac{\partial f_j(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i} \frac{\partial f_{j'}(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_{i'}} & \left(\mathbf{R}(\boldsymbol{\theta}) \text{ sym.} \right) \\
&= (\mathbf{y} - f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta}))^\top \mathbf{R}^{-1}(\boldsymbol{\theta}) \left[ \frac{\partial^2 f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i \partial \beta_{i'}}\right] - \left[\frac{\partial f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_i}\right]^\top \mathbf{R}^{-1}(\boldsymbol{\theta})  \left[ \frac{\partial f(\mathbf{X}, \boldsymbol{\alpha}, \mathbf{Z}, \boldsymbol{\beta})}{\partial \beta_{i'}} \right]
\end{aligned}
$$


$$
\begin{aligned}
\left. \left[ \frac{\partial}{\partial \boldsymbol{\beta} \partial \boldsymbol{\beta}^\top} \left[ \ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right] \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}}
&= \left. \left[  \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}} \\
\end{aligned}
$$

$$
\begin{aligned}
\ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta})
&\approx 
\left. \left[ \ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}} 
+ \left. \boldsymbol{\beta}^\top \left[ \frac{\partial}{\partial \boldsymbol{\beta}} \left[ \ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right] \right]  \right\rvert_{\boldsymbol{\beta} = \mathbf{0}} 
+ \frac{1}{2} \boldsymbol{\beta}^\top \left. \left[ \frac{\partial^2}{\partial \boldsymbol{\beta} \partial \boldsymbol{\beta}^\top} \left[ \ell(\mathbf{y}; \boldsymbol{\alpha}, \boldsymbol{\theta} \rvert \boldsymbol{\beta}) \right] \right] \right\rvert_{\boldsymbol{\beta} = \mathbf{0}}  \\
&= 
\end{aligned}
$$