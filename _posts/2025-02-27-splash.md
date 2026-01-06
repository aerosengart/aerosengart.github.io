---
layout: distill
title: SPLASH and OASIS
description: 
date: 2025-02-27
tabs: true
tags: sequencing genomics paper-review
toc:
  - name: Background
    subsections:
        - name: Notation
        - name: Rewriting $S$
  - name: Bounds
  - name: Optimization
bibliography: 2025-02-27-splash.bib
---

There is an exciting new framework for reference-free genomic discovery called <i>SPLASH</i> (Statistically Primary aLignment Agnostic Sequence Homing) from the Salzman Lab at Stanford<d-cite key=chaung2023></d-cite><d-cite key=baharav2024></d-cite>. It's super cool, and there are a lot of new questions that arise from their work. In this post, I'm going to be working through some of the proofs of their $p$-value bounds and asymptotic results (just for my own understanding).

---

## Background
We can imagine that the genetic sequences (DNA, RNA) can be $k$-merized into little chunks. We'll define an <i>anchor</i> to be $k$-mer and a <i>target</i> to be another one that is downstream (by some specified number of base pairs) from an anchor. We will have many, many anchors when analyzing samples, and we will also have many, many targets. 

Since there is variability in genetic sequences, an anchor will be associated with many different targets. Our concern will be with testing for structured variation in the targets across multiple samples. Our data come to us as count tables, and, in this case, we will be concerned with a single table that corresponds to one anchor.

This table, $X$, will be $I \times J$ where $I$ is the number of targets and $J$ is the number of samples. The null hypothesis is that (conditional on the column sums, denoted by $n_j$ for $j \in [n]$), each column of $X$ is independent and multinomial with shared row probability vector, $\mathbf{p}$:

$$
\begin{equation}
\label{eq:oasis-null}
X^{(j)} \rvert n_j \overset{ind}{\sim} \text{Multinomial}(n_j, \mathbf{p}) \hspace{5mm} \forall j \in [n]
\end{equation}
$$


### Notation 
Denote the total counts in column $j$ with $n_j := \sum_{i = 1}^I X_{i,j}$ and the $J$-vector of column sums with $\mathbf{n}:= X^\top \mathbf{1}$. Denote the total counts in row $i$ with $r_i := \sum_{j = 1}^J X_{i, j}$ and the $I$-vector of row sums with $\mathbf{r} := X \mathbf{1}$. Let the grand total be $M := \sum_{j = 1}^J n_j = \sum_{i = 1}^I r_i$. 

Further define:

$$
\gamma := \gamma(\mathbf{n}, \mathbf{c}) = \bigg\langle \frac{\mathbf{c}}{\rvert \rvert \mathbf{c} \rvert \rvert}, \sqrt{\frac{\mathbf{n}}{M}} \bigg\rangle^2
$$

We have the following property for $\gamma$:

<div class="theorem">
<strong>Claim (Bounds of $\gamma$).</strong>
{% tabs gamma-defn %}
{% tab gamma-defn statement %}
$$\gamma \in [0, 1]$$
{% endtab %}
{% tab gamma-defn proof %}
$$
\begin{aligned}
\gamma &= \left\langle \frac{\mathbf{c}}{\rvert\rvert \mathbf{c} \rvert\rvert}, \sqrt{\frac{n}{M}} \right\rangle^2 \\
&\overset{(i)}{\leq} \left\langle \frac{\mathbf{c}}{\rvert\rvert \mathbf{c} \rvert\rvert}, \frac{\mathbf{c}}{\rvert\rvert \mathbf{c} \rvert\rvert} \right\rangle \cdot \left\langle \sqrt{\frac{n}{M}}, \sqrt{\frac{n}{M}} \right\rangle \\
&= \bigg\rvert\bigg\rvert \frac{\mathbf{c}}{\rvert\rvert \mathbf{c} \rvert\rvert} \bigg\rvert\bigg\rvert^2 \cdot \bigg\rvert\bigg\rvert \sqrt{\frac{n}{M}} \bigg\rvert\bigg\rvert^2 \\
&\overset{(ii)}{=} \bigg\rvert\bigg\rvert \sqrt{\frac{n}{M}} \bigg\rvert\bigg\rvert^2 \\
&\overset{(iii)}{=} \sum_{j = 1}^J \left( \sqrt{\frac{n_j}{M}} \right)^2 \\
&= \frac{1}{M} \sum_{j = 1}^J n_j \\
&\overset{(iv)}{=} 1
\end{aligned}
$$

$(i)$ is due to Cauchy-Schwarz. $(ii)$ follows from the fact that $\frac{\mathbf{c}}{\rvert\rvert \mathbf{c} \rvert\rvert}$ has unit norm. $(iii)$ is from the fact that we are dealing with the Euclidean norm. $(iv)$ is due to the fact that $M = \sum_{j = 1}^J n_j$. Since $\gamma$ is a squared value, it cannot be negative, thus $\gamma \in [0, 1]$ as desired. 
{% endtab %}
{% endtabs %}
</div>

We also know that:

<div class="theorem">
<strong>Claim (Expected Counts).</strong>
{% tabs exp-count %}
{% tab exp-count statement %}
Under the null, the expected count table is given by:

$$
E := \frac{1}{M} X\mathbf{1}\mathbf{1}^\top X
$$
{% endtab %}
{% tab exp-count proof %}
Under the null (see Eq. \eqref{eq:oasis-null}), we have i.i.d. multinomial observations $$$X_{\cdot, 1}, \dots, X_{\cdot, J}$$, each of dimension $$I \times 1$$. The log-likelihood (conditional on the column sums) can be written as:

$$
\begin{aligned}
\ell(\mathbf{p}; X_{\cdot, 1}, \dots, X_{\cdot, J}, \mathbf{n}) &= \log\left(\prod_{j = 1}^J \frac{n_j!}{X_{1, j}! \times \dots \times X_{I,j}!} \times \mathbf{p}_1^{X_{1, j}} \times \dots \times \mathbf{p}_I^{X_{I, j}} \right) \\
&= \sum_{j = 1}^J \left( \log(n_j!) + \sum_{i = 1}^I \log\left(\mathbf{p}_i^{X_{i,j}}\right) - \sum_{i = 1}^I \log(X_{i,j}!)\right)
\end{aligned}
$$

We want to maximize the above over $$\mathbf{p}$$, but we have to enforce the constraint that $$\sum_{i = 1}^I \mathbf{p}_i = 1$$ (equivalently, $$1 - \sum_{i = 1}^I \mathbf{p}_i = 0$$). We use the method of Lagrangian multipliers. The Lagrangian for this problem is:

$$
\mathcal{L}(\mathbf{p}) = \ell(\mathbf{p}; X_{\cdot, 1}, \dots, X_{\cdot, J}, \mathbf{n}) + \lambda\left(1 - \sum_{i = 1}^I \mathbf{p}_i\right)
$$

The partial derivative of the Lagrangian with respect to the $k$-th element of $\mathbf{p}$ is:

$$
\begin{aligned}
\frac{\partial}{\partial \mathbf{p}_k}\left[ \mathcal{L}(\mathbf{p}) \right] &= \frac{\partial}{\partial \mathbf{p}_k} \left[ \ell(\mathbf{p}; X_{\cdot, 1}, \dots, X_{\cdot, J}, \mathbf{n}) \right] + \frac{\partial}{\partial \mathbf{p}_k} \left[ \lambda\left(1 - \sum_{i = 1}^I \mathbf{p}_i\right)\right] \\
&= \frac{\partial}{\partial \mathbf{p}_k} \left[ \sum_{j = 1}^J \left(\log(n_j!) + \sum_{i = 1}^I \log \left( \mathbf{p}_i^{X_{i,j}} \right) - \sum_{i = 1}^I \log(X_{i,j}!) \right) \right] - \lambda \\
&= \sum_{j = 1}^J \sum_{i = 1}^I \frac{\partial}{\partial \mathbf{p}_k} \left[ \log\left(\mathbf{p}_i^{X_{i,j}} \right) - \log(X_{i,j}!)\right] - \lambda \\
&= \sum_{j = 1}^J \frac{X_{k,j}}{\mathbf{p}_k} - \lambda
\end{aligned}
$$

We set the above equal to $0$ and solve for $\mathbf{p}_k$:

$$
\begin{aligned}
&0 =  \sum_{j = 1}^J \frac{X_{k,j}}{\mathbf{p}_k} - \lambda  \\
\implies
&\mathbf{p}_k = \sum_{j = 1}^J \frac{X_{k,j}}{\lambda}
\end{aligned}
$$

Next we take the above and solve for $\lambda$ using our constraint:

$$
\begin{aligned}
&\mathbf{p}_k = \sum_{j = 1}^J \frac{X_{k,j}}{\lambda} \\
\implies 
&\sum_{k = 1}^I \mathbf{p}_k = \sum_{k = 1}^I \sum_{j = 1}^J \frac{X_{k,j}}{\lambda} \\
\implies
&1 = \sum_{k = 1}^I \sum_{j = 1}^J \frac{X_{k,j}}{\lambda} \\
\implies
&\lambda = M
\end{aligned}
$$

where $M$ is the table sum (i.e. sum over all $X_{i,j}$). Thus, the MLE for $\mathbf{p}_k$ is $$\frac{1}{M} \sum_{j = 1}^J X_{k,j}$$, and the MLE for $$\mathbf{p}$$ is $$\frac{1}{M} X \mathbf{1}$$ where $\mathbf{1}$ is a $J \times 1$ vector of ones. We then see that the $(i,j)$-th entry of $E$ is: 

$$
\begin{aligned}
&E_{i,j} = n_j \mathbf{p}_i = \frac{n_j}{M} \sum_{k = 1}^J X_{i,k} = \frac{1}{M} n_j r_i \\
\implies
&E = 
\begin{bmatrix}
\frac{1}{M} n_1 r_1 & \dots &  \frac{1}{M} n_J r_1 \\
\vdots & \ddots & \vdots \\
\frac{1}{M} n_1 r_I & \dots &  \frac{1}{M} n_J r_I 
\end{bmatrix} =
\frac{1}{M} \mathbf{r} \mathbf{n}^\top = 
\frac{1}{M} X \mathbf{1} \mathbf{1}^\top X
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Denote the centered and normalized count table with $\tilde{X} := (X - E) \text{diag}\left[ (X^\top \mathbf{1})^{-1/2}\right]$. Define input vectors $\mathbf{f} \in \mathbb{R}^I$ and $\mathbf{c} \in \mathbb{R}^J$ which are meant to partition the count table according to the targets and samples, respectively. The test statistic for SPLASH is given by:

$$
S = S(\mathbf{f}, \mathbf{c}) = \mathbf{f}^\top \tilde{X} \mathbf{c}
$$

#### Rewriting $S$
Later proofs will require rewriting the above statistic in the following way. Define $Z_{j,k}$ as a random variable that denotes the target index for the $k$-th target observed in the $j$-th sample. Consider sample $j$, which is associated with $n_j$ total counts. $Z_{j,k}$ will be the target index (that is, $k$ will be between $1$ and $n_j$) that the $k$-th count for sample $j$ is associated with. For example, suppose $X_{i,j} = 3$, then we will have $Z_{j, 1}$, $Z_{j, 2}$, and $Z_{j, 3}$, and all of them will equal $i$. 

Since $\mathbf{f}$ is an $I$-vector, we can use the $Z$ variable to index $\mathbf{f}$. Also notice that, using the above as an example, if $X_{i,j} = 3$, then $f_i X_{i,j} = 3 f_i$, which is equivalent to $3f_{Z_{j, k}}$. Define the following:

$$
\begin{aligned}
\bar{\mu} &:= \frac{1}{M} \sum_{i = 1}^I \sum_{j = 1}^J \mathbf{f}_i X_{i,j} = \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
\hat{\mu}_j &:= \frac{1}{n_j} \sum_{i = 1}^I \mathbf{f}_i X_{i,j} = \frac{1}{n_j} \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
S_j &:= \sqrt{n_j}(\hat{\mu}_j - \bar{\mu})
\end{aligned}
$$

<div class="theorem">
<strong>Claim (Rewriting $S$).</strong>
{% tabs s-rewrite %}
{% tab s-rewrite statement %}
With the above, we can rewrite the test statistic as the weighted sum of the $S_j$ values where the weights are determined by $\mathbf{c}$:

$$
\begin{equation}
\label{eq:s-rewrite}
S = \sum_{j = 1}^J \mathbf{c}_j  S_j = \sum_{j = 1}^J \mathbf{c}_j \left(\sqrt{n_j} (\hat{\mu}_j - \bar{\mu}) \right) 
\end{equation}
$$
{% endtab %}
{% tab s-rewrite proof %}
First notice that:

$$
\begin{aligned}
\bar{\mu} &= \frac{1}{M} \sum_{i = 1}^I \sum_{j = 1}^J \mathbf{f}_i X_{i,j} = \frac{1}{M} \sum_{i = 1}^I \mathbf{f}_{i} \underbrace{\sum_{j = 1}^J X_{i,j}}_{\text{row sums}} = \frac{1}{M} \mathbf{f}^\top X \mathbf{1}
\end{aligned}
$$

We can also write the vector of $\hat{\mu}_j$'s as:

$$
\begin{aligned}
\hat{\mu}^\top &:= \begin{bmatrix} \hat{\mu}_1 \\ \vdots \\ \hat{\mu}_I \end{bmatrix}^\top \\
&= \begin{bmatrix}
\frac{1}{n_1} \sum_{i = 1}^I \mathbf{f}_i X_{i,1} \\
\vdots \\
\frac{1}{n_J} \sum_{i = 1}^I \mathbf{f}_i X_{i,J}
\end{bmatrix}^\top \\
&= 
\mathbf{f}^\top X
\begin{bmatrix}
\frac{1}{n_1} & \dots & 0 \\
\vdots & \ddots & \vdots \\
0 & \dots & \frac{1}{n_J}
\end{bmatrix}
= \mathbf{f}^\top X \text{diag}\left( \mathbf{n}\right)^{-1}
\end{aligned}
$$

We can then define the vector $$S^*$$, which is the vector of the $S_j$'s:

$$
\begin{aligned}
S^* &:= (\hat{\mu}^\top - \bar{\mu}\mathbf{1}^\top)\text{diag}(\mathbf{n})^{1/2}\\
&= \left(\mathbf{f}^\top X \text{diag}(\mathbf{n})^{-1} - \frac{1}{M} \mathbf{f}^\top X \mathbf{1} \mathbf{1}^\top \right)\text{diag}(\mathbf{n})^{1/2} \\
&= \mathbf{f}^\top \left(X \text{diag}(\mathbf{n})^{-1} - \frac{1}{M} X \mathbf{1} \mathbf{1}^\top\right)\text{diag}(\mathbf{n})^{1/2} \\
&= \mathbf{f}^\top \left(X  - \frac{1}{M} X \mathbf{1} \mathbf{1}^\top\text{diag}\left(\mathbf{n}\right) \right)\text{diag}(\mathbf{n})^{-1} \text{diag} (\mathbf{n})^{1/2} \\
&= \mathbf{f}^\top\left(X - \frac{1}{M} X \mathbf{1}\mathbf{1}^\top \text{diag}(\mathbf{n}) \right) \text{diag}(\mathbf{n})^{-1/2}
\end{aligned}
$$

Notice that:

$$
\mathbf{1}^\top \text{diag}(\mathbf{n}) = \mathbf{n}^\top =  \mathbf{1}^\top X
$$

Plugging this result into the equation for $$S^*$$ yields:

$$
\begin{aligned}
S^* &= \mathbf{f}^\top\left(X - \frac{1}{M} X \mathbf{1}\mathbf{1}^\top X\right)\text{diag}(\mathbf{n})^{-1/2} \\
&= \mathbf{f}^\top\left(X - E \right)\text{diag}(\mathbf{n})^{-1/2} \\
&= \mathbf{f}^\top \tilde{X}
\end{aligned}
$$

Thus, the test statistic is:

$$
S := S^* \mathbf{c} = \mathbf{f}^\top (X - E) \text{diag}(\mathbf{n})^{-1/2} \mathbf{c}
$$
{% endtab %}
{% endtabs %}
</div>

---

## Bounds
Though an exact $p$-value is difficult (maybe impossible?) to derive, one can derive a finite-sample value $p$-value <i>bound</i>. 

<div class="theorem">
<strong>Proposition 1.<d-cite key=chaung2023></d-cite></strong>
{% tabs chaung-1 %}
{% tab chaung-1 statement %}
For fixed $\mathbf{f} \in [0, 1]^I$ and $\mathbf{c} \in \mathbb{R}^J$ satisfying $\rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1$, if $\gamma < 1$, then, under the null hypothesis:

$$
\mathbb{P}(\rvert S \rvert \geq s) \leq 2 \exp\left(-\frac{2(1-a)^2s^2}{\sum_{j: n_j > 0} \mathbf{c}_j^2}\right) + 2 \exp \left(-\frac{2 M(as)^2}{\left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}\right)^2} \right)
$$

where:

$$
a = \left(1 + \frac{\sqrt{M \sum_{j: n_j < 0} \mathbf{c}_j^2}}{\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}}\right)^{-1}
$$
{% endtab %}
{% tab chaung-1 proof %}
Under the null hypothesis, the $Z$'s are just categorical random variables with event probabilities, $\mathbf{p}$. We'll define the following $\mu := \mathbb{E}_{Z \sim \mathbf{p}}[\mathbf{f}_Z]$. We can manipulate the rewritten version of $S$ in Eq. \eqref{eq:s-rewrite} as:

$$
S = \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \left(\frac{1}{n_j}\sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \bar{\mu} \right) = \sum_{j = 1}^J \mathbf{c}_j \sum_{k = 1}^{n_j} \sqrt{n_j} \left( \frac{1}{n_j} \mathbf{f}_{Z_{j,k}} - \frac{1}{n_j} \bar{\mu} \right) = \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \bar{\mu}}{\sqrt{n_j}} \right)
$$

Adding and subtracting $\mu$ in the numerator of the righthand side above yields:

$$
\begin{aligned}
\sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \bar{\mu}}{\sqrt{n_j}} \right) &= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu + \mu - \bar{\mu}}{\sqrt{n_j}} \right) \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right) +  \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mu - \bar{\mu}}{\sqrt{n_j}}\right) \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right) +  (\mu - \bar{\mu}) \sum_{j = 1}^J  \frac{n_j \mathbf{c}_j}{\sqrt{n_j}} \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right) +  (\mu - \bar{\mu}) \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j} \\
\end{aligned}
$$

We can use the above in the probability statement:

$$
\begin{aligned}
\mathbb{P}(\rvert S \rvert \geq s) &= \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \bar{\mu}}{\sqrt{n_j}} \right) \bigg\rvert \geq s\right) \\
&= \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right) +  (\mu - \bar{\mu}) \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j}  \bigg\rvert \geq s \right) \\
&\overset{(i)}{\leq} \underset{a \in (0, 1)}{\min} \left\{ \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right)  \bigg\rvert \geq  (1 - a) s\right)  + \mathbb{P}\left( \bigg\rvert (\mu - \bar{\mu}) \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j} \bigg\rvert \geq a s \right) \right\} \\
&\overset{(ii)}{=} \underset{a \in (0, 1)}{\min} \left\{ \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right)  \bigg\rvert \geq  (1 - a) s\right)  + \mathbb{P}\left(  \bigg\rvert \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} (\mathbf{f}_{Z_{j,k}} - \mu) \bigg\rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big \rvert} \right) \right\} \\
&\overset{(iii)}{\leq} \underset{a \in (0, 1)}{\min} \left\{ 2 \exp\left(- \frac{2((1-a)s)^2}{\sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{ \mathbf{c}_j}{\sqrt{n_j}}\right)^2} \right) + 2 \exp \left(-\frac{2 M(as)^2}{\left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}\right)^2} \right) \right\} \\
&\overset{(iv)}{\leq} \underset{a \in (0, 1)}{\min} \left\{ 2 \exp\left(-\frac{2(1-a)^2s^2}{\sum_{j: n_j > 0} \mathbf{c}_j^2}\right) + 2 \exp \left(-\frac{2 M(as)^2}{\left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}\right)^2} \right) \right\} 
\end{aligned}
$$

<details>
  <summary>Details Of $(i)$.</summary>
  In $(i)$, we use the triangle inequality and probability laws. Recall that the triangle inequality states $\rvert A + B \rvert \leq \rvert A \rvert + \rvert B \rvert$. It follows that:
  $$
  \begin{aligned}
  \mathbb{P}(\rvert A + B \rvert \geq \epsilon) &\leq \mathbb{P}(\rvert A \rvert + \rvert B \rvert \geq \epsilon) \\
  &= \mathbb{P}(\rvert A \rvert \geq (1 - a)\epsilon \cap \rvert B \rvert \geq a \epsilon) \hspace{15mm} \text{ for any } a \in [0, 1] \\
  &\leq \mathbb{P}(\rvert A \rvert \geq (1 - a) \epsilon \cup \rvert B \rvert \geq a \epsilon) \\
  &\leq \mathbb{P}(\rvert A \rvert \geq (1 - a) \epsilon) + \mathbb{P}(\rvert B \rvert \geq a \epsilon) \\
  \implies \mathbb{P}(\rvert A + B \rvert \geq \epsilon) &\leq \underset{a \in (0, 1)}{\min} \left\{\mathbb{P}(\rvert A \rvert \geq (1 - a) \epsilon) + \mathbb{P}(\rvert B \rvert \geq a \epsilon)  \right\}
  \end{aligned}
  \nonumber
  $$
</details>

<details>
  <summary>Details Of $(ii)$.</summary>
  In $(ii)$, we divide both sides of inequality in the second term by $\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert$, which we assume is not $0$:
  $$
  \begin{aligned}
  \bigg\rvert (\mu - \bar{\mu}) \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j} \bigg\rvert \geq a s \hspace{10mm} &\implies \hspace{10mm} \frac{\bigg\rvert (\mu - \bar{\mu}) \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j} \bigg\rvert}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} \\
  &\implies \hspace{10mm} \rvert \mu - \bar{\mu} \rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} \\
  &\implies \hspace{10mm} \rvert \bar{\mu} - \mu \rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} \\
  &\implies \hspace{10mm} \bigg\rvert \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \mu \bigg\rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} \hspace{15mm} \left(\text{ plug-in for } \hat{\mu} \right)\\
  &\implies \hspace{10mm}  \bigg\rvert \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_J} \mu \bigg\rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} \hspace{15mm} \left(\sum_{j = 1}^J \sum_{k = 1}^{n_j} 1 = M\right) \\
  &\implies \hspace{10mm} \bigg\rvert \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} (\mathbf{f}_{Z_{j,k}} - \mu) \bigg\rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big\rvert} 
  \end{aligned}
  \nonumber
  $$
</details>

<details>
  <summary>Details Of $(iii)$.</summary>
  Recall Hoeffding's inequality for bounded random variables. 
  <div class="theorem">
    <body>
    <strong>Theorem. (Hoeffding's Inequality)</strong>
    <br>
    Let $X_1, \dots, X_n$ be independent and bounded random variables (i.e. $a_i \leq X_i \leq b_i$ for all $i = 1, \dots, n$ almost surely). For any choice of $t > 0$:
    $$
    \begin{aligned}
    &\mathbb{P}\left( \sum_{i = 1}^n X_i - \mathbb{E}\left[ \sum_{i = 1}^n X_i \right] \geq t \right)  \leq  \exp\left(-\frac{2t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \right) \\
    \implies &\mathbb{P}\left( \bigg\rvert \sum_{i = 1}^n X_i - \mathbb{E}\left[ \sum_{i = 1}^n X_i \right] \bigg\rvert \geq t \right)  \leq  2\exp\left(-\frac{2t^2}{\sum_{i = 1}^n (b_i - a_i)^2} \right)
    \end{aligned}
    $$
    </body>
  </div>
  Consider the first term, $X_{j,k} := \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \mu)$. This has zero mean:
  $$
  \mathbb{E}\left[ X_{j,k} \right] =\mathbb{E}\left[ \frac{\mathbf{c}_j}{\sqrt{n_j}}(\mathbf{f}_{Z_{j,k}} - \mu) \right] = \frac{\mathbf{c}_j}{\sqrt{n_j}} \mathbb{E}\left[ \mathbf{f}_{Z_{j,k}} - \mu\right] = \frac{\mathbf{c}_j}{\sqrt{n_j}} (\mathbb{E}[\mathbf{f}_{Z_{j,k}}] - \mu) = 0
  \nonumber
  $$
  Further note that the $X_{j,k}$'s are bounded and independent by construction. This is due to the fact that $\mathbf{f} \in [0, 1]^I$, thus $\mu \in [0, 1]$ as well. Thus:
  $$
  X_{j,k} \in \left[ -\frac{\mathbf{c}_j \mu}{\sqrt{n_j}}, \frac{ \mathbf{c}_j(1-\mu)}{\sqrt{n_j}} \right] \nonumber
  $$
  Letting $a_{j,k} := - \frac{\mathbf{c}_j \mu }{\sqrt{n_j}}$ and $b_{j,k} := \frac{\mathbf{c}_j (1- \mu)}{\sqrt{n_j}}$, we can use the triangle inequality to bound their difference.
  $$
  b_{j,k} - a_{j,k} =  \frac{\mathbf{c}_j (1- \mu) }{\sqrt{n_j}} - \left( - \frac{ \mathbf{c}_j \mu }{\sqrt{n_j}} \right) = \frac{ \mathbf{c}_j + (- \mathbf{c}_j \mu)  +  \mathbf{c}_j \mu }{\sqrt{n_j}} = \frac{ \mathbf{c}_j }{\sqrt{n_j}}
  \nonumber
  $$
  Choosing $t = (1- a)s$, we can apply Hoeffding's inequality to the probability term:
  $$
  \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{c}_j \left( \frac{\mathbf{f}_{Z_{j,k}} - \mu }{\sqrt{n_j}}\right)  \bigg\rvert \geq  (1 - a) s\right) = \mathbb{P}\left( \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} X_{j,k} \bigg\rvert \geq (1-a)s \right) \leq 2 \exp\left(- \frac{2((1-a)s)^2}{ \sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j }{\sqrt{n_j}}\right)^2} \right)
  \nonumber
  $$
  For the second probability term, we can let $Y_{j,k} := \mathbf{f}_{Z_{j,k}} - \mu$, which clearly has mean $0$ and is bounded between $- \mu$ and $1 - \mu$. Letting $b_{j,k} = 1 - \mu$ and $a_{j,k} = -\mu$, we can see that $b_{j,k} - a_{j,k} = 1$ for all $j,k$. Thus:
  $$
  \mathbb{P}\left(  \bigg\rvert \frac{1}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} (\mathbf{f}_{Z_{j,k}} - \mu) \bigg\rvert \geq \frac{a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big \rvert} \right) = \mathbb{P}\left(  \bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} (\mathbf{f}_{Z_{j,k}} - \mu) \bigg\rvert \geq \frac{M a s}{\big\rvert \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \big \rvert} \right) \leq 2 \exp \left(-\frac{\frac{2 (Mas)^2}{\left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}\right)^2}}{\sum_{j = 1}^J \sum_{k = 1}^{n_j} 1^2} \right) = 2 \exp \left(-\frac{2 M(as)^2}{\left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}\right)^2} \right)
  \nonumber
  $$
</details>

<details>
<summary>Details Of $(iv)$.</summary>
  In $(iv)$, we do some manipulation and simplification of the terms. If we assume all $n_j > 0$, then:
  $$
  \sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} \right)^2 = \sum_{j = 1}^J \mathbf{c}_j^2 \sum_{k = 1}^{n_j} \frac{1}{n_j} = \sum_{j = 1}^J \mathbf{c}_j^2 \geq \sum_{j: n_j > 0} \mathbf{c}_j^2 \implies \exp\left(- \frac{2(1-a)^2 s^2}{\sum_{j = 1}^J \mathbf{c}_j^2}\right) \leq \exp\left(-\frac{2(1-a)^2s^2}{\sum_{j: n_j > 0} \mathbf{c}_j^2}\right)
  \nonumber
  $$
</details>

The authors claim that the bound can be optimized within a factor of $2$ by setting the two terms equal to each other and solving for $a$ (I am not sure why, and I cannot easily find a reference). Thus, the $p$-value bound is:

$$
\mathbb{P}(\rvert S \rvert \geq s) \leq 2 \exp\left(-\frac{2(1-a)^2s^2}{\sum_{j: n_j > 0} \mathbf{c}_j^2}\right) + 2 \exp \left(-\frac{2 M(as)^2}{\left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}\right)^2} \right)
$$

where:

$$
a = \left(1 + \frac{\sqrt{M \sum_{j: n_j < 0} \mathbf{c}_j^2}}{\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}}\right)^{-1}
$$

<details>
  <summary>Details Of $a$.</summary>
  $$
  \begin{aligned}
  2 \exp\left( -\frac{2(1-a)^2 s^2}{\sum_{j: n_j < 0} \mathbf{c}_j^2} \right) = 2 \exp\left(-\frac{2M(as)^2}{\left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2}\right) \implies  &\frac{(1-a)^2}{\sum_{j: n_j < 0} \mathbf{c}_j^2} = \frac{Ma^2}{\left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2} \\
  \implies &\frac{(1-a)^2}{M a^2} = \frac{\sum_{j: n_j < 0} \mathbf{c}_J^2}{\left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2} \\
  \implies &\frac{1-a}{a\sqrt{M}} = \sqrt{\frac{\sum_{j: n_j < 0}\mathbf{c}_J^2}{\left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2}} \\
  \implies &\frac{1}{a} - 1 = \frac{\sqrt{M \sum_{j: n_j < 0} \mathbf{c}_j^2}}{\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}} \\
  \implies &a = \left(1 + \frac{\sqrt{M \sum_{j: n_j < 0} \mathbf{c}_j^2}}{\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j}}\right)^{-1}
  \end{aligned}
  \nonumber
  $$
</details>

{% endtab %}
{% endtabs %}
</div>

The previous proposition asserts a bound on the $p$-value for the SPLASH test statistic that can be improved. This improvement is shown in Baharav et al.<d-cite key=baharav2024></d-cite>

<div class="theorem">
<strong>Proposition 3.<d-cite key=baharav2024></d-cite></strong>
{% tabs baharav-3 %}
{% tab baharav-3 statement %}
For fixed $\mathbf{f} \in \mathbb{R}^{I}$ and $\mathbf{c} \in \mathbb{R}^J$, if $\gamma = 1$, then $S = 0$ with probability $1$. If $\gamma < 1$, then, under the null hypothesis:

$$
\mathbb{P}(\rvert S(\mathbf{f}, \mathbf{c}) \rvert \geq s) \leq 2 \exp \left( - \frac{s^2}{2 \rvert \rvert \mathbf{f} \rvert \rvert_{\infty}^2 \rvert \rvert \mathbf{c} \rvert \rvert_2^2  ( 1- \gamma)} \right)
$$
{% endtab %}
{% tab baharav-3 proof %}
The first claim is fairly simple to prove, so I've hidden it.

<details>
  <summary>Proof Of First Claim.</summary>
  For the first claim, note that, if $\gamma = 1$, then:
  $$
  \begin{aligned}
  1 = \gamma &\implies 1 = \left( \sum_{j =1 }^J \frac{\mathbf{c}_j}{\rvert \rvert \mathbf{c} \rvert \rvert_2}\sqrt{\frac{n_j}{M}}\right)^2 \implies \rvert \rvert \mathbf{c} \rvert \rvert_2^2 = \left(\sum_{j = 1}^J \mathbf{c}_j \sqrt{\frac{n_j}{M}}\right)^2 \implies \rvert \rvert \mathbf{c} \rvert \rvert_2^2 = \bigg\langle \mathbf{c}, \sqrt{\frac{\mathbf{n}}{M}} \bigg\rangle^2 
  \end{aligned}
  \nonumber
  $$
  By inspection, it is clear that the last equality holds only if $\mathbf{c} = \pm \sqrt{\frac{\mathbf{n}}{M}}$. Thus:
  $$
  \begin{aligned}
  S(\mathbf{f}, \mathbf{c}) &= \sum_{j = 1}^J \mathbf{c}_j(\sqrt{n_j}(\hat{\mu}_j - \bar{\mu})) \\
  &= \sum_{j = 1}^J \frac{\mathbf{c}_j \sqrt{n_j}}{n_j} \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \sum_{l = 1}^J \frac{\mathbf{c}_l \sqrt{n_l}}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
  &= \sum_{j = 1}^J \frac{\pm\sqrt{\frac{n_j}{M}} \sqrt{n_j}}{n_j} \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j, k}} - \sum_{l = 1}^J \frac{\pm \sqrt{\frac{n_l}{M}}\sqrt{n_l}}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \hspace{20mm} (\text{ plugging in } \mathbf{c}) \\
  &= \sum_{j = 1}^J \pm\sqrt{\frac{1}{M}} \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \sum_{l = 1}^J \pm \frac{n_l}{M \sqrt{M}} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
  &= \pm\sqrt{\frac{1}{M}} \sum_{j = 1}^J\sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \pm \frac{1}{\sqrt{M}} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
  &= 0
  \end{aligned}
  \nonumber
  $$
</details>

As before, we denote the expectation of $\mathbf{f}_Z$ with $\mu$. The improvement begins with a manipulation of $S$:

$$
\begin{equation}
\label{eq:s-rewrite-2}
\begin{aligned}
S(\mathbf{f}, \mathbf{c}) &= \sum_{j = 1}^J \mathbf{c}_j(\sqrt{n_j} (\hat{\mu}_j - \bar{\mu})) \\
&= \sum_{j = 1}^J \frac{\mathbf{c}_j \sqrt{n_j}}{n_j} \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} - \sum_{l = 1}^J \frac{\mathbf{c}_l \sqrt{n_l}}{M} \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \frac{\mathbf{c}_j}{\sqrt{n_j}} \mathbf{f}_{Z_{j,k}} - \frac{\sum_{l = 1}^J \mathbf{c}_j \sqrt{n_l}}{M}\sum_{j = 1}^J \sum_{k = 1}^{n_j} \mathbf{f}_{Z_{j,k}} \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right) \mathbf{f}_{Z_{j,k}} \\
&\overset{(i)}{=}\sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right) \mathbf{f}_{Z_{j,k}} - \sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right) \mu \\
&= \sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right) (\mathbf{f}_{Z_{j,k}} - \mu)
\end{aligned}
\end{equation}
$$

<details>
  <summary>Details Of $(i)$.</summary>
  In $(i)$, we note that:
  $$
  \begin{aligned}
  \sum_{j = 1}^J \sum_{k = 1}^{n_j} \mu \left( \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right) &= \mu \left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} - \sum_{j = 1}^J n_j \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}\right) \\
  &= \mu \left( \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}  \sum_{j = 1}^J n_j\right) \\
  &= \mu \left( \sum_{j = 1}^J  \mathbf{c}_j \sqrt{n_j} - \sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}\right) \\
  &= 0
  \end{aligned}
  \nonumber
  $$
</details>

Now define $$Y_{j,k} := \left( \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right)(\mathbf{f}_{Z_{j,k}} - \mu)$$. Recall that we define $$\mu$$ as the expectation of $$\mathbf{f}_Z$$ for $$Z \sim \mathbf{p}$$. This implies that the $$Y_{j,k}$$'s have mean zero. Furthermore, for fixed $\mathbf{c}$ and $\mathbf{f}$, the largest and smallest that $Y_{j,k}$ can be only depends on the biggest and smallest values of $\mathbf{f}$. It follows that:

$$
\begin{aligned}
\max Y_{j,k} - \min Y_{j,k} &= \max \left\{ \left( \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right)(\mathbf{f}_{Z_{j,k}} - \mu) \right\} - \min\left \{ \left( \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right)(\mathbf{f}_{Z_{j,k}} - \mu) \right\} \\
&\leq \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \bigg\rvert (\max \left\{ \mathbf{f}_{Z_{j,k}} \right\} - \mu) - \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \bigg\rvert \left( \min \left\{ \mathbf{f}_{Z_{j,k}} \right\} - \mu \right) \\
&= \bigg\rvert \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \bigg\rvert \left( \max \left\{ \mathbf{f}_{Z_{j,k}} \right\} - \min \left\{ \mathbf{f}_{Z_{j,k}} \right\} \right)
\end{aligned}
$$

We then apply Hoeffding's inequality as we did in the original proof, this time using the form of $S$ in Eq. \eqref{eq:s-rewrite-2}:

$$
\begin{aligned}
\mathbb{P}\left( \rvert S \rvert \geq s\right) &= \mathbb{P}\left(\bigg\rvert \sum_{j = 1}^J \sum_{k = 1}^{n_j} \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right) (\mathbf{f}_{Z_{j,k}} - \mu) \bigg\rvert \geq s \right) \\
&\leq 2 \exp \left(- \frac{2s^2}{\sum_{j = 1}^J \sum_{k = 1}^{n_j} \left( \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right)\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)\right)^2 }\right) \\
&\overset{(i)}{=} 2 \exp\left(- \frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 \left(\rvert \rvert \mathbf{c} \rvert \rvert_2^2 + \sum_{j = 1}^J \left(n_j\frac{\left(\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l} \right)^2}{M^2} - 2 \mathbf{c}_j \sqrt{n_j} \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}  \right) \right)}\right) \\
&\overset{(ii)}{=} 2 \exp\left(- \frac{2s^2}{ \left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 \left(\rvert \rvert \mathbf{c} \rvert \rvert_2^2 - \frac{1}{M} \left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2\right)}\right) \\
\end{aligned}
\nonumber
$$

<details>
  <summary>Details Of $(i)$.</summary>
  This follows from straightforward manipulation of the denominator:
  $$
  \begin{aligned}
   2 \exp \left(- \frac{2s^2}{\sum_{j = 1}^J \sum_{k = 1}^{n_j} \left( \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \right)\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)\right)^2 }\right) &= 2 \exp \left(-\frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2\sum_{j = 1}^J \sum_{k = 1}^{n_j}\left( \frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}\right)^2}\right) \\
  &= 2\exp \left(- \frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2\sum_{j = 1}^J n_j \left(\frac{\mathbf{c}_j}{\sqrt{n_j}} - \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}\right)^2}\right) \\
  &= 2 \exp \left( - \frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2\sum_{j = 1}^J \left(\mathbf{c}_j - \sqrt{n_j}\frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}\right)^2} \right) \\
  &= 2 \exp \left(-\frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2\sum_{j = 1}^J \left( \mathbf{c}_j^2 - 2 \mathbf{c}_j \sqrt{n_j} \frac{\sum_{l = 1}^J \mathbf{c}_l\sqrt{n_l}}{M} + n_j \frac{\left(\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l} \right)^2}{M^2} \right)} \right) \\
  &= 2 \exp\left(- \frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 \left(\rvert \rvert \mathbf{c} \rvert \rvert_2^2 + \sum_{j = 1}^J \left(n_j\frac{\left(\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l} \right)^2}{M^2} - 2 \mathbf{c}_j \sqrt{n_j} \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}  \right) \right)}\right) 
  \end{aligned}
  \nonumber
  $$
  In this step, we have assumed that $\sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \neq 0$. If this is not the case, then $a$ can be taken to be $0$, and the $p$-value bound follows from there.
</details>

<details>
  <summary>Details Of $(ii)$.</summary>
  In $(i)$, we simplify the summation in the denominator:
  $$
  \begin{aligned}
  \sum_{j = 1}^J \left(n_j\frac{\left(\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l} \right)^2}{M^2} - 2 \mathbf{c}_j \sqrt{n_j} \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M}  \right) &= \frac{\left(\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}\right)^2}{M^2} \sum_{j = 1}^J n_j - 2 \frac{\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l}}{M} \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \\
  &= \frac{\left(\sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l} \right)^2}{M} - \frac{2\left( \sum_{l = 1}^J \mathbf{c}_l \sqrt{n_l} \right)^2}{M} \\
  &= - \frac{1}{M} \left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2
  \end{aligned}
  \nonumber
  $$
</details>

Next, we rewrite $\gamma$:

$$
\gamma = \bigg\langle \frac{\mathbf{c}}{\rvert \rvert \mathbf{c} \rvert \rvert_2}, \sqrt{\frac{\mathbf{n}}{M}} \bigg\rangle^2 = \left(\sum_{j = 1}^J \frac{\mathbf{c}_j \sqrt{n_j}}{\sqrt{M} \rvert \rvert \mathbf{c} \rvert \rvert_2} \right)^2 = \frac{1}{M \rvert \rvert \mathbf{c} \rvert \rvert_2^2}\left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2 \implies \rvert \rvert \mathbf{c}\rvert\rvert_2^2 = \frac{1}{M} \left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2
\nonumber
$$

We can plug this into the bound above to get:

$$
\begin{equation}
\label{eq:improved-bound}
\mathbb{P}\left(\rvert S \rvert \geq s\right) \leq  2\exp\left( - \frac{2s^2}{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 \left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2(1 - \gamma)}\right)  \leq 2 \exp\left(- \frac{s^2}{2 \rvert \rvert \mathbf{f}\rvert\rvert_{\infty}^2 \rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma)}\right) 
\end{equation}
$$

<details>
<summary>Details Of Bound.</summary>
$$
\begin{aligned}
\mathbb{P}\left(\rvert S \rvert \geq s\right) &\leq 2 \exp\left( - \frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2\left(\rvert \rvert \mathbf{c} \rvert \rvert_2^2 - \frac{1}{M}\left( \sum_{j = 1}^J \mathbf{c}_j \sqrt{n_j} \right)^2\right)}  \right) \\
&= 2 \exp\left(- \frac{2s^2}{\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2\left(\rvert \rvert \mathbf{c} \rvert \rvert_2^2 - \gamma\rvert\rvert \mathbf{c} \rvert \rvert^2 \right)}\right) \\
&= 2\exp\left( - \frac{2s^2}{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 \left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2(1 - \gamma)}\right) \\
&\overset{(i)}{\leq} 2 \exp\left(- \frac{2s^2}{4 \rvert \rvert \mathbf{f}\rvert\rvert_{\infty}^2 \rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma)}\right) \\
&= 2 \exp\left(- \frac{s^2}{2 \rvert \rvert \mathbf{f}\rvert\rvert_{\infty}^2 \rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma)}\right) 
\end{aligned}
\nonumber
$$

where $(i)$ follows from:

$$
\begin{aligned}
\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 &\leq \left(\underset{i \in [I]}{\max} \rvert\mathbf{f}_i\rvert - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 \\
&\leq \left(\underset{i \in [I]}{\max} \rvert\mathbf{f}_i\rvert + \underset{i \in [I]}{\min} \rvert\mathbf{f}_i \rvert \right)^2 \\
&\leq \left( \underset{i \in [I]}{\max} \rvert\mathbf{f}_i\rvert+ \underset{i \in [I]}{\max} \rvert\mathbf{f}_i \rvert  \right)^2 \\
&= \left(2 \underset{i \in [I]}{\max} \rvert \mathbf{f}_i \rvert \right)^2 \\
&= 4 \rvert \rvert \mathbf{f} \rvert \rvert_{\infty}^2
\end{aligned}
\nonumber
$$
</details>
{% endtab %}
{% endtabs %}
</div>

Baharav et al. also provide a simplified bound that can be useful in practice.

<div class="theorem">
<strong>Proposition 1.<d-cite key=baharav2024></d-cite></strong>
{% tabs baharav-1 %}
{% tab baharav-1 statement %}
For fixed $\mathbf{f} \in [0, 1]^I$ and $\mathbf{c} \in \mathbb{R}^J$ such that $\rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1$, if $\gamma = 1$, then $S = 0$ with probability $1$. If $\gamma < 1$, then, under the null hypothesis:

$$
\mathbb{P}(\rvert S(\mathbf{f}, \mathbf{c}) \rvert \geq s) \leq 2 \exp\left(- \frac{2s^2}{1 - \gamma}\right)
\nonumber
$$
{% endtab %}
{% tab baharav-1 proof %}
This follows from the previous proof. Using the intermediate bound in Eq. \eqref{eq:improved-bound}, we know that $\mathbf{f} \in [0, 1]^I$ and $\rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1$, which implies:

$$
\begin{aligned}
&\left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 \leq  1 \\
\implies 
&\rvert \rvert \mathbf{c} \rvert \rvert_2^2 \left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 \leq  1 \\
\implies &\exp\left(-\frac{2s^2}{\rvert \rvert \mathbf{c} \rvert \rvert^2_2 \left(\underset{i \in [I]}{\max} \mathbf{f}_i - \underset{i \in [I]}{\min} \mathbf{f}_i \right)^2 (1 - \gamma)}\right) \leq \exp\left(-\frac{2s^2}{1-\gamma}\right)
\end{aligned}
$$

since shrinking the denominator will increase the value inside the exponential, which is negative. 
{% endtab %}
{% endtabs %}
</div>

---

## Optimization

The optimization procedure for choosing $\mathbf{f}$ and $\mathbf{c}$ is based upon a lemma concerning the containment of solutions to the optimization problem at hand within the set of solutions to a problem that is simpler to solve. In the following I use $\left[ \cdot \right]$ to denote a set. 

<div class="theorem">
<strong>Lemma 1.<d-cite key=baharav2024></d-cite></strong>
{% tabs lemma-1-baharav %}
{% tab lemma-1-baharav statement %}
The optimal $\mathbf{f}$ and $\mathbf{c}$ for the $p$-value bound in Proposition 1 are contained in the set of optimal solution to maximizing $\mathbf{f}^\top \tilde{X} \mathbf{c}$. Written differently:

$$
\left[ \underset{0 \leq \mathbf{f} \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c} \right\} \right] \subseteq \left[ \underset{0 \leq \mathbf{f} \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \min} \left\{ 2 \exp\left( - \frac{2(\mathbf{f}^\top \tilde{X} \mathbf{c})^2}{1 - \gamma}\right) \right\} \right]
$$
{% endtab %}
{% tab lemma-1-baharav proof %}
To begin, the authors shift and rescale $\mathbf{f}$ to be $\rvert \mathbf{f}_i \rvert \leq 1$ for all $i$. This doesn't really matter since one can recover the argument maximum/minimum of the original problem by un-shifting and un-scaling the argument maximum/minimum of the alternative problem. Furthermore, note that the restriction that $\rvert \mathbf{f}_i \rvert \leq 1$ is equivalent to $\max_i \mathbf{f}_i \leq 1$, which is the same as $$\rvert \rvert \mathbf{f} \rvert \rvert_\infty \leq 1$$. Rewriting the constraints on $\mathbf{f}$, our goal is to show that:

$$
\left[\underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c} \right\} \right] \subseteq \left[ \underset{\rvert \rvert \mathbf{f} \rvert \rvert_\infty \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \min} \left\{ \frac{(\mathbf{f}^\top \tilde{X}\mathbf{c})^2}{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma)} \right\} \right] 
$$

Recall the following properties of $\arg \max$ and $\arg \min$. For any function, $g$, any monotonic function, $f$, and positive scalar, $\alpha$:

$$
\begin{array}{l}
\arg\max_x g(x) = \arg\min_x -g(x) \\
\arg\max_x g(x) = \arg\max_x f(g(x)) \\
\arg\max_x g(x) = \arg\max_x \alpha g(x)
\end{array}
$$

Taking $g(\cdot) = \frac{(\mathbf{f}^\top \tilde{X} \mathbf{c})^2}{\rvert \rvert\mathbf{c} \rvert \rvert_2^2 (1 - \gamma)}$, $f(\cdot) = 2\exp(\cdot)$ and $\alpha = 2$, these properties imply we can instead show:

$$
\left[\underset{\rvert \rvert \mathbf{f}\rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c} \right\} \right] \subseteq \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \frac{(\mathbf{f}^\top \tilde{X}\mathbf{c})^2}{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma)} \right\} \right] 
$$

Recall that $\gamma = \big\langle \frac{\mathbf{c}}{\rvert\rvert\mathbf{c}\rvert\rvert_2}, \sqrt{\frac{\mathbf{n}}{M}} \big \rangle^2$. The denominator of the righthand side of the above can be rewritten as:

$$
\begin{aligned}
\rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma) &= \rvert \rvert \mathbf{c} \rvert \rvert_2^2 \left(1 - \bigg\langle \frac{\mathbf{c}}{\rvert \rvert \mathbf{c} \rvert \rvert_2}, \sqrt{\frac{\mathbf{n}}{M}} \bigg \rangle^2\right) \\
&= \rvert \rvert \mathbf{c} \rvert \rvert_2^2  - \frac{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 }{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 M}\big\langle \mathbf{c}, \sqrt{\mathbf{n}} \big\rangle^2 \\
&= \rvert \rvert \mathbf{c} \rvert \rvert_2^2 - \frac{1}{M} \big\langle \mathbf{c}, \sqrt{\mathbf{n}} \big\rangle^2
\end{aligned}
$$

Next, we decompose $\mathbf{c}$ into the part that is orthogonal to $\sqrt{\mathbf{n}}$ and the part parallel to $\sqrt{\mathbf{n}}$. That is, for some $\alpha \in \mathbb{R}$, we can write $\mathbf{c} = \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}}$. Plugging this into the objective function in the righthand side yields:

$$
\begin{aligned}
\frac{(\mathbf{f}^\top \tilde{X} \mathbf{c})^2}{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 - \frac{1}{M} \big\langle \mathbf{c}, \sqrt{\mathbf{n}} \big\rangle^2} 
&= \frac{(\mathbf{f}^\top \tilde{X}(\mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}}))^2 }{\rvert \rvert \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}} \rvert \rvert_2^2 - \frac{1}{M} \langle (\mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}}), \sqrt{\mathbf{n}} \rangle^2} \\
&= \frac{(\mathbf{f}^\top \tilde{X} \mathbf{c}^\perp + \alpha \mathbf{f}^\top \tilde{X} \sqrt{\mathbf{n}})^2}{\rvert \rvert \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}} \rvert \rvert_2^2 - \frac{1}{M} (\langle \mathbf{c}^\perp, \sqrt{\mathbf{n}} \rangle + \alpha \langle \sqrt{\mathbf{n}}, \sqrt{\mathbf{n}} \rangle)^2} \\
&= \frac{(\mathbf{f}^\top \tilde{X} \mathbf{c}^\perp + \alpha \mathbf{f}^\top \tilde{X} \sqrt{\mathbf{n}})^2}{\rvert \rvert \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}} \rvert \rvert_2^2 - \frac{\alpha^2}{M} \rvert \rvert \sqrt{\mathbf{n}} \rvert \rvert_2^2}  \hspace{20mm} (\mathbf{c}^\perp \perp \sqrt{\mathbf{n}} \implies \langle \mathbf{c}^\perp, \sqrt{\mathbf{n}} \rangle = 0 ) \\
&\overset{(i)}{=} \frac{(\mathbf{f}^\top \tilde{X} \mathbf{c}^\perp + \alpha \mathbf{f}^\top \tilde{X} \sqrt{\mathbf{n}})^2}{\rvert \rvert \mathbf{c}^\perp \rvert\rvert_2^2 + \alpha^2 \rvert \rvert \sqrt{\mathbf{n}} \rvert \rvert_2^2 - \frac{\alpha^2}{M} \rvert \rvert \sqrt{\mathbf{n}} \rvert \rvert_2^2} \\
&= \frac{(\mathbf{f}^\top \tilde{X} \mathbf{c}^\perp + \alpha \mathbf{f}^\top \tilde{X} \sqrt{\mathbf{n}})^2}{\rvert \rvert \mathbf{c}^\perp \rvert\rvert_2^2 + \alpha^2 M - \frac{\alpha^2M^2}{M}}  \hspace{20mm} (\rvert \rvert \sqrt{\mathbf{n}} \rvert \rvert_2^2 = \sum_{j = 1}^J \sqrt{n_j}\sqrt{n_j} = \sum_{j = 1}^J n_j = M) \\
&\overset{(ii)}{=} \frac{(\mathbf{f}^\top \tilde{X} \mathbf{c}^\perp)^2}{\rvert \rvert \mathbf{c}^\perp \rvert\rvert_2^2}
\end{aligned}
$$

<details>
  <summary>Details Of $(i)$.</summary>
  In $(i)$, we use the fact that:
  $$
  \begin{aligned}
  \rvert \rvert \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}} \rvert \rvert_2^2 &= \big\langle \mathbf{c}^\perp+ \alpha \sqrt{\mathbf{n}}, \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}}\big\rangle \\
  &= \langle \mathbf{c}^\perp, \mathbf{c}^\perp \rangle + 2 \alpha \langle \mathbf{c}^\perp, \sqrt{\mathbf{n}} \rangle + \alpha^2 \langle \sqrt{\mathbf{n}}, \sqrt{\mathbf{n}} \rangle \\
  &= \langle \mathbf{c}^\perp, \mathbf{c}^\perp \rangle + \alpha^2 \langle \sqrt{\mathbf{n}}, \sqrt{\mathbf{n}} \rangle \hspace{30mm} (\mathbf{c}^\perp \perp \sqrt{\mathbf{n}} \implies \langle \mathbf{c}^\perp, \sqrt{\mathbf{n}} \rangle = 0 ) \\
  &= \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2^2 + \alpha^2 \rvert \rvert \sqrt{\mathbf{n}} \rvert \rvert_2^2 
  \end{aligned}
  $$
</details>

<details>
  <summary>Details Of $(ii)$.</summary>
  $$
  \begin{aligned}
  \tilde{X}\sqrt{\mathbf{n}} &= \left(X - \frac{1}{M} X \mathbf{1} \mathbf{1}^\top X \right) \text{diag}\left(\sqrt{\mathbf{n}}\right)^{-1} \sqrt{\mathbf{n}} \\
  &= \left(X - \frac{1}{M}X \mathbf{1} \mathbf{1}^\top X \right) \mathbf{1} \\
  &= X \mathbf{1} - \frac{1}{M}X \mathbf{1}\mathbf{1}^\top X \mathbf{1} \\
  &= X \mathbf{1} - \frac{1}{M}(M) X\mathbf{1}  \hspace{20mm} (\mathbf{1}^\top X \mathbf{1} = \sum_{i = 1}^I \sum_{j = 1}^J X_{i,j} = M) \\
  &= 0
  \end{aligned}
  $$
</details>

The above derivations imply:

$$
\left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \frac{(\mathbf{f}^\top \tilde{X}\mathbf{c})^2}{\rvert \rvert \mathbf{c} \rvert \rvert_2^2 (1 - \gamma)} \right\} \right] = \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \frac{(\mathbf{f}^\top \tilde{X}\mathbf{c}^\perp)^2}{\rvert \rvert \mathbf{c}^\perp \rvert \rvert^2_2} \right\} \right] 
$$

It is clear that $\max_x g(x) \leq \max_x (g(x))^2$ since $g(x) \leq (g(x))^2$ for any $x$. Thus, any $x$ (and $-x$) achieving a maximum of $g(x)$ will also achieve a maximum of $(g(x))^2$. This implies:

$$
\left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \frac{(\mathbf{f}^\top \tilde{X}\mathbf{c}^\perp)^2}{\rvert \rvert \mathbf{c}^\perp \rvert \rvert^2_2} \right\} \right] \supseteq \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \frac{\mathbf{f}^\top \tilde{X} \mathbf{c}^\perp}{\rvert \rvert \mathbf{c}^\perp \rvert \rvert_2} \right\} \right] = \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 = 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c}^\perp\right\} \right]
$$

Looking back at the fact that $\tilde{X}\sqrt{\mathbf{n}} = 0$, we can then add back in the part of $\mathbf{c}$ that is parallel to $\sqrt{n}$:

$$
\begin{aligned}
\left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 = 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c}^\perp\right\} \right] &= \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 = 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c}^\perp + \mathbf{f}^\top \tilde{X}(\alpha \sqrt{\mathbf{n}}) \right\} \right] \\
&= \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 = 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} (\mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}}) \right\} \right] 
\end{aligned}
$$

Restricting $\rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 = 1$ implies that, by the triangle inequality:

$$
\begin{aligned}
\rvert \rvert \mathbf{c} \rvert \rvert_2 &= \rvert \rvert \mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}} \rvert \rvert_2 \\
&\leq \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 + \rvert \rvert \alpha\sqrt{\mathbf{n}} \rvert \rvert_2 \\
\implies \rvert \rvert \mathbf{c} \rvert \rvert_2 &\leq 1 + \rvert \rvert \alpha \sqrt{\mathbf{n}} \rvert \rvert_2
\end{aligned}
$$

Now, notice that the set of $\mathbf{c}$ vectors satisfying the last line above is a subset of the vectors satisfying $\rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1$, since $\rvert \rvert \alpha \sqrt{\mathbf{n}} \rvert \rvert_2 \geq 0$. Thus:

$$
\left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c}^\perp \rvert \rvert_2 = 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} (\mathbf{c}^\perp + \alpha \sqrt{\mathbf{n}}) \right\} \right] \supseteq \left[ \underset{\rvert \rvert \mathbf{f} \rvert\rvert_\infty \leq 1, \rvert \rvert \mathbf{c} \rvert \rvert_2 \leq 1}{\arg \max} \left\{ \mathbf{f}^\top \tilde{X} \mathbf{c} \right\} \right] 
$$
{% endtab %}
{% endtabs %}
</div>