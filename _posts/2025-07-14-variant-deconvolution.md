---
layout: distill
title: Variant Deconvolution
description: A Primer
date: 2025-07-14
tabs: true
tags: virus deconvolution epidemiology primer
toc:
  - name: VirPool, Model-Based Estimation of SARS-CoV-2 Variant Proportions in Wastewater Samples
    subsections:
        - name: Variant Profiles
        - name: Model
  - name: Estimating the Relative Proportions of SARS-CoV-2 Haplotypes from Wastewater Samples
    subsections:
        - name: Model
        - name: Imputation
        - name: Estimation
  - name: Unsupervised Detection and Fitness Estimation of Emerging SARS-CoV-2 Variants
    subsections:
        - name: Model
  - name: Tracking SARS-CoV-2 Genomic Variants in Wastewater Sequencing Data with LolliPop
    subsections:
        - name: Model
        - name: Temporal Regularization
        - name: Uncertianty Quantification
bibliography: 2025-07-14-variant-deconvolution.bib
---

In the case of a rapidly mutating virus such as SARS-CoV-2, the genetic material in a wastewater sample will be a mixture from different variants. More specifically, the sample will contain small segments of viral RNA, and these small segments will be from some unknown mixture of genomes where the mixture proportions are determined by the prevalence of the different variants in the population served by that sewer system. The goal of variant deconvolution is to estimate these proportions in order to detect the presence of a new or rising variant in an area or to aid in assessing variant fitness. 

One idea of doing this estimation is to use the number of times one sees different variant-defining mutations in the sample. However, an added layer of difficult is that there may be new variants (called "cryptic" variants), and their defining mutations may be unknown as well. Several, relatively new, methods have been proposed to tackle the problem of variant deconvolution in an unsupervised manner by using mutation frequencies. Though this does not allow for the identification of the variants, it does allow for cryptic variants. One could also simply match the variants to the frequent mutations to "annotate" the results. 

One problem that arises with wastewater samples is that the sequencing is performed on the amplified genetic material. Since the amount of viral RNA is quite small compared to the other things in the sample, tiling amplicon PCR is used to increase the quantity. However, this procedure could be biased towards or against particular sequences, which could lead to overrepresentation of certain regions of the genome or missing coverage of others. Another issue is that the genetic material is fragmented, which means that it is no longer possible to identify which mutations occured on the same molecule. Moreover, the sample preparation process can introduce additional bias and errors. Variant deconvolution methods are often adjusted to account for one or more of these effects.<d-cite key=dreifuss2022></d-cite>

In this post, I'm going to go through a brief literature review of variant deconvolution.

---

## VirPool, Model-Based Estimation of SARS-CoV-2 Variant Proportions in Wastewater Samples
In the work of Gafurov et al.<d-cite key=gafurov2022></d-cite> the variant proportions are estimated using a mixture model. This method is supervised in that the number and "identity" of the variants must be specified in the model. 

### Variant Profiles
The identity is determined by a <i>variant profile</i> where the variant profile for variant $k$ is denoted by $P_k(i, a)$ and represents the probability of it having nucleotide $a$ at position $i$ (in the reference genome). 

These profiles are constructed for the variants of interest using data from GISAID. Each variant is associated with a Pango lineage identifier, and each sample in the GISAID database can be assigned to the closest ancestor clade of the lineages of our variants of interest (or assigned to an "other" category).

Due to minimal coverage at certain positions caused by deletions, the relative frequency needed to calculate $P_k(i, a)$ must be adjusted, leading to the form:

$$
P_k(i, a) = \frac{C_k(i, a)}{\max \left\{ \gamma_k, \sum_{b} C_k(i, b) \right\}}
$$

where $\gamma_k$ is a coverage threshold for variant $k$, and $C_k(i, a)$ is the number of times we observe base $a$ at position $i$ in the genomes (i.e. GISAID samples) associated with variant $k$. 


### Model
Let $R = (r_1, r_2, \dots, r_{\rvert R \rvert})$ be a sequencing read. Suppose it starts at position $s$ and comes from variant $k$. We define the probability of this read as:

$$
\mathbb{P}(R \rvert k) = \prod{i = 1}^{\rvert R \rvert} P_k(s + i - 1, r_i)
$$

Let $\rho = (R_1, R_2, \dots, R_m)$ denote a sequence of $m$ reads. Assuming independence of reads and given variant weights $w_1, \dots, w_K$, the model assigns a likelihood to the probability of observing $\rho$:

$$
\mathcal{L}(w_1, \dots, w_K; \rho) = \mathbb{P}\left( \rho \rvert w_1, \dots, w_K \right) = \prod_{i = 1}^m \left[ \sum_{k = 1}^k w_k \mathbb{P}(R_i \rvert k) \right]
$$

where $w_i$ reflect the proportion of reads from variant $i$ in the sample. The weights $W = (w_1, \dots, w_K)$ for reads $R = (R_1, \dots, R_m)$ are then estimated as those values that maximize this likelihood as:

$$
W^* = \text{softmax} \left\{ \underset{(\xi_1, \dots, \xi_K) \in \mathbb{R}^K}{\arg \min} \left\{ - \sum_{i = 1}^m \log\left( \sum_{k = 1}^K \frac{\exp(\xi_k)}{\sum_{j = 1}^K \exp(\xi_j)} \mathbb{P}(R_i \rvert k) \right) \right\} \right\},
\hspace{5mm}
w_i = \frac{\exp(\xi_i)}{\sum_{j = 1}^K \exp(\xi_j)}
$$

The model is also extend to include uniformly at random sequencing error (though it does not address insertions). Let $\epsilon \geq 0$ be the error rate. The model becomes:

$$
\mathbb{P}\left( R  \rvert k, \epsilon \right) = \prod_{i = 1}^{\rvert R \rvert} \left[ (1 - \epsilon) P_k(s + i - 1, r_i) + \frac{\epsilon}{3} \sum_{ a \neq r_i} P_k(s + i - 1, a) \right]
$$

From what I can tell, they do not justify the factor $\frac{1}{3}$ in the second term...

---

## Estimating the Relative Proportions of SARS-CoV-2 Haplotypes from Wastewater Samples
Pipes et al.<d-cite key=pipes2022></d-cite> take a similar approach and construct a likelihood function of haplotype proportions that they maximize using an expectation-maximization algorithm. Since a haplotype is just a collection of mutations that are usually inherited together, we can think of them as representative of a variant (since the genomes associated with a given variant are highly similar). 

### Model
Let $r$ denote the number of sequencing reads in the experiment, and let $g$ denote the number of different genomes comprising the constructed reference database (i.e. a collection of samples from GISAID). The authors first define an $r \times g$ matrix, $D$, where the $(i,j)$-th element, $D_{i,j}$, denotes the number of sequences mismatches between sequencing read $i$ and genome $j$. (These mismatches can be found because the reads should be aligned to the same reference genome as the sequences in the MSA used to form the reference database.)

Let $q_{i,j}$ denote the probability of seeing read $i$ given it is from haplotype $j$. This is assumed to have the form:

$$
\begin{equation}
\label{eq:q-defn}
q_{i,j} = \alpha^{D_{i,j}} \times (1 - \alpha)^{n_j - D_{i,j}}
\end{equation}
$$

where $\alpha \geq 0$ is some (chosen) error rate, and $n_j$ is the length of read $j$. Assuming independence of the reads, the log-likelihood can be written as:

$$
\ell(p_1, \dots, p_k) = \sum_{j = 1}^r \log\left( \sum_{i = 1}^k q_{i,j} p_{i} \right)
$$

where $p_i$ denotes the relative proportion of haplotype $i$. Haplotypes with identical values for $q_{i,j}$ over all $j \in [n]$ are collapsed into a single group to keep the model identifiable.

### Imputation
Pipes et al. also introduce a tree-based method of imputing missing nucleotides via phylogenetic data. Their imputation method is based on the Jukes-Cantor model<d-cite key=jukes1969></d-cite> (also called the JC69 model) for DNA evolution, which itself relies on some Markov chain theory (see below).

<details>
<summary>Markov Chain Review.</summary>
We'll restrict our discussion to continuous-time Markov chains.
<div id="markov-chain"></div>
<div class="definition">
<strong>Definition (Continuous-Time Markov Chain).<d-cite key=anderson1991></d-cite></strong>
<br>
Let $(\Omega, \mathcal{A}, P)$ be a probability space, and let $S$ be a countable, non-empty set (called the <i>state space</i>). A stochastic process, $\{ X(t), t \in [0, \infty) \}$, defined on $(\Omega, \mathcal{A}, P)$ is called a <i>continuous-time parameter Markov chain</i> if, for any set $\{ t_i \}_{i = 1}^{n+1}$ such that $0 \leq t_1 < t_2 < \dots < t_n < t_{n+1}$ with corresponding set of states $\{ s_i \}_{i = 1}^{n+1}$ in $S$ such that:
$$
\mathbb{P}\left(X(t_n) = s_n, X(t_{n-1}) = X_{s_{n-1}}, \dots, X(t_1) = s_1 \right) > 0
\nonumber
$$
satisfies the Markov property. That is:
$$
\mathbb{P}\left(X(t_{n+1}) = s_{n+1} \rvert X(t_n) = s_n, X(t_{n-1}) = s_{n-1}, \dots, X(t_1) = s_1\right)  \mathbb{P}\left(X(t_{n+1}) = s_{n+1} \rvert X(t_n) = s_n \right)
\nonumber
$$
</div>
A continuous-time Markov chain can be described by its <i>transition function</i>, which is essentially the probability of transitioning between states in $S$.
<div id="markov-chain"></div>
<div class="definition">
<strong>Definition (Transition Function).<d-cite key=anderson1991</d-cite></strong>
<br>
Let $S$ be a state space. A <i>transition function</i>, denoted by $Q_{i,j}(t)$ for $i,j \in S$ and $t \geq 0$, is a function defined on $S$ that satisfies:
<ol>
    <li>
        $Q_{i,j}(t) \geq 0$ for all $t \geq 0$ and $i,j \in S$.
    </li>
    <li>
        $Q_{i,j}(0) = \delta_{i,j}$ for any $i,j \in S$ where $\delta_{i,j} = 1$ if $i = j$ and $0$ otherwise.
    </li>
    <li>
        $\sum_{j \in S} Q_{i,j}(t) \leq 1$ for all $t \geq 0$ and for all $i \in S$.
    </li>
    <li>
        $P_{i,j}(t + r) = \sum_{k \in S} P_{i,k}(t) P_{k,j}(r)$ for all $t, r \geq 0$ and for all $i,j \in S$.
    </li>
</ol>
In Property (3) above, if $\sum_{j \in S} Q_{i,j}(t) = 1$, then the transition function is called <i>honest</i>. Property $(4)$ is called the <i>Chapman-Kolmogorov equation</i>.
<br>
If the transition function additionally satisfies:
$$
\underset{t \rightarrow 0}{\lim} Q_{i,i}(t) = 1
$$
for all $i \in S$, then it is called <i>standard</i>. This property implies that $Q_{i,j}(t) \rightarrow \delta_{i,j}$ for all $i,j \in S$.
</div>
<br>
Since $Q_{i,j}(t)$ is just a function mapping from $S^2$ to $\mathbb{R}$, we can use it to define a (possibly infinitely dimensional) matrix whose $(i,j)$-th entry corresponds to the value $Q_{i,j}(t)$ and represents the probability of the process transitioning from state $i$ to state $j$ at time $t$. 
<br>
If we define $\pi^(0)$ to be the vector of initial state probabilities (i.e. the probability of the Markov chain starting in each state $s \in S$), we can define $\pi^(n)$ (the $n$-step transition probabilities) which is the vector of state probabilities at time $n$. $\pi^{(n)}_{i,j}$ represents the probability of the Markov chain being in state $j$ after $n$ transitions when starting in state $i$. 
<br>
It follows that:
$$
\pi^{(n)} = \pi^{(0)} Q^n, \hspace{5mm} \pi^{(n+1)} = \pi^{(n)} Q
\nonumber
$$
We now cover several important properties of Markov chain states:
<ul>
    <li>
        <i>Communicable</i>: A state $j \in S$ is called <i>accessible</i> (denoted $i \rightarrow j$) from state $j \in S$ if $\pi_{i,j}^{(n)}(t) > 0$. States $i,j \in S$ can <i>communicate</i> (denoted $i \leftrightarrow j$) if they are accessible from each other. 
    </li>
    <li>
        <i>Periodicity</i>: Let $k = \text{GCD}\left\{ n > 0 \rvert \mathbb{P}\left(X(n) = s \rvert X(0) = s \right) > 0 \right\}$ for some state $s \in S$. That is, $k$ is the greatest common divisor of the number of transitions needed for the Markov chain to start at state $s$ and return to $s$. $k$ is called the <i>period</i> of state $s$. If $k > 1$, then $s$ is called <i>periodic</i> (and <i>aperiodic</i> if $k = 1$).
    </li>
    <li>
        <i>Absorbing</i>: A state $s \in S$ is called <i>absorbing</i> if the probability of transitioning out of $s$ is zero.
    </li>
    <li>
        <i>Recurrence</i>: A state $s \in S$ is called <i>recurrent</i> if the probability of the process never returning to $s$ is zero. Otherwise, it is called <i>transient</i>.
        <ul>
            <li>
                The <i>hitting time</i>, $H_A$, of a set $A \subset S$ is the random variable representing the first time at which the Markov chain transitions to a state in $A$:
                    $$
                        H_A = \min\left\{ t \geq 0 \rvert X(t) \in A \right\}
                        \nonumber
                    $$
            </li>
            <li>
                The <i>hitting probability</i>, $h_{s, A}$, of a set $A \subset S$ is the probability from starting state $s$ is the probability of the Markov chain ever transitioning to a state in $A$ when its starting state is $s$:
                    $$
                        h_{s, A} = \mathbb{P}\left( X(t) \in A \hspace{2mm} \text{ for some } t \geq 0 \rvert X(0) = s \right) = \mathbb{P}\left( H_A < \infty \rvert X(0) = s \right)
                        \nonumber
                    $$
            </li>
            <li>
                The <i>expected hitting time</i>, $\eta_{s, A}$, is the expected hitting time of set $A \subset S$ from state $s$:
                    $$
                        \eta_{s, A} = \mathbb{E}\left[ H_A \rvert X(0) = s \right]
                        \nonumber
                    $$
                $\eta_{i, A} < \infty$ only if $h_{s, A} = 1$.
            </li>
        </ul>
        A state $s$ is called <i>positive recurrent</i> if $\eta_{s, A} < 1$. That is, if the expected time to return to state $s$ is finite. Otherwise, $s$ is called <i>null recurrent</i>.
    </li>
    <li>
        A state $s \in S$ is called <i>ergodic</i> if it is aperiodic and positive recurrent.
    </li>
</ul>
Similar properties can be defined for the Markov chain itself:
<ul>
    <li>
        <i>Irreducibility</i>: A Markov chain is called <i>irreducible</i> if, and only if, $i \leftrightarrow j$ for all $i,j \in S$.
    </li>
    <li>
        <i>Ergodic</i>: A Markov chain is called <i>ergodic</i> if it is irreducible and all of its states are ergodic.
    </li>
    <li>
        <i>Regular</i>: A Markov chain is called <i>regular</i> if, for some power $n$, $Q^n$ has all positive entries. A regular Markov chain can be represented by a fully connected graph.
    </li>
    <li>
        <i>Homogeneous</i>: A Markov chain is called <i>homogeneous</i> if, and only if, the transition probabilities are independent of time. In this case, $Q^(n) = Q^(n-1)$ for all $n$.
    </li>
    <li>
        <i>Stationary</i>: A Markov chain is called <i>stationary</i> if it has a stationary distribution, defined as a probability vector, $pi$, such that $\pi Q = \pi$. If the Markov chain is irreducible and aperiodic, then $\pi$ is unique. 
    </li>
</ul>
</details>

<details>
<summary>Jukes-Cantor Model.</summary>
The basic idea is to model the substitution of base pairs with a continuous-time Markov model. Suppose we have an ancestor sequence and a descendant sequence of this ancestor; in other words, we have observe a Markov process at time $0$ and some later time $n$. If we know the starting state at a position $i$ in the ancestor sequence (i.e. the initial state), the transition rate matrix ($R$ below), and the expected number of changes (<i>branch length</i>) between ancestor and descendant, then we can recover the probability of observing the descendant (i.e. the transition probability matrix).
<br>
Since this model is for DNA evolution, there are four states, each corresponding to one of the nucelotide bases: $S = \{A, G, C, T\}$. We define the following transition matrix:
$$
Q(t) = 
\begin{bmatrix}
    Q_{A,A}(t) & Q_{A,G}(t) & Q_{A,C}(t) & Q_{A,T}(t) \\
    Q_{G,A}(t) & Q_{G,G}(t) & Q_{G,C}(t) & Q_{G,T}(t) \\
    Q_{C,A}(t) & Q_{C,G}(t) & Q_{C,C}(t) & Q_{C,T}(t) \\
    Q_{T,A}(t) & Q_{T,G}(t) & Q_{T,C}(t) & Q_{T,T}(t)
\end{bmatrix}
\nonumber
$$
We define the <i>transition rate</i> matrix as:
$$
R = 
\begin{bmatrix}
    - \mu_A & \mu_{A, G} & \mu_{A, C} & \mu_{A, T} \\
    \mu_{G, A} & -\mu_{G} & \mu_{G, C} & \mu_{G, T} \\
    \mu_{C, A} & \mu_{C, G} & -\mu_{C} & \mu_{C, T} \\
    \mu_{T, A} & \mu_{T, G} & \mu_{T, C} & -\mu_{T}
\end{bmatrix}
\nonumber
$$
where $\mu_{i, j}$ is the instantaneous rate of change from state $i$ to state $j$. The diagonal elements, $\mu_i$, represent the total rate of change to state $i$. 
<br>
The Jukes-Cantor model makes a few simplifying assumptions. First, the state probabilities are assumed to be equal and unchanging over time, so $\pi(t) = \left(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4} \right)$. The model is parametrized by $\mu$, which represents the substitution rate in the genome. Additionally, the transition probability matrix is assumed to be multiplicative, so the transition probability at time $t + t'$ is equal to the product of the probabilities at time $t$ and time $t'$. Thus, the transition rate and transition probability matrices are given by:
$$
R = 
\begin{bmatrix}
    -3\mu & \mu & \mu & \mu \\
    \mu & - 3\mu & \mu & \mu \\
    \mu & \mu & -3\mu & \mu \\
    \mu & \mu & \mu & -3\mu
\end{bmatrix},
\hspace{5mm}
Q(t) =
\begin{bmatrix}
    r_t & q_t & q_t & q_t \\
    q_t & r_t & q_t & q_t \\
    q_t & q_t & r_t & q_t \\
    q_t & q_t & q_t & r_t
\end{bmatrix}
\nonumber
$$
Finding $Q(t)$ is our goals. Notice that, for some very short period of time, $\epsilon > 0$:
$$
Q(\epsilon) =
\mathbb{I}_{4 \times 4} + R\epsilon =
\begin{bmatrix}
    1 - 3 \mu \epsilon & \mu \epsilon & \mu \epsilon & \mu_\epsilon \\
    \mu_\epsilon & 1 - 3 \mu \epsilon & \mu_\epsilon & \mu_\epsilon \\
    \mu_\epsilon & \mu_\epsilon & 1 - 3 \mu \epsilon & \mu_\epsilon \\
    \mu_\epsilon & \mu_\epsilon & \mu_\epsilon & 1 - 3 \mu \epsilon
\end{bmatrix}
\nonumber
$$
We can use the multiplicativity property to see:
$$
\begin{aligned}
&Q(t + \epsilon) = Q(t) Q(\epsilon) = Q(t) (\mathbb{I}_{4 \times 4} + R \epsilon) \\
\implies
&Q(t + \epsilon) - Q(t) = Q(t) (\mathbb{I}_{4 \times 4} + R \epsilon) - Q(t)  \\
\implies 
&Q(t + \epsilon) - Q(t) = Q(t) R \epsilon  \\
\implies
&\frac{1}{\epsilon} \left[ Q(t + \epsilon) - Q(t) \right] = Q(t) R
\end{aligned}
\nonumber
$$
Taking $\epsilon \rightarrow 0$, we get the differential equations:
$$
\begin{aligned}
\frac{\partial}{\partial t} Q(t) 
&= Q(t) R  \\
&= \begin{bmatrix}
    r_t & q_t & q_t & q_t \\
    q_t & r_t & q_t & q_t \\
    q_t & q_t & r_t & q_t \\
    q_t & q_t & q_t & r_t
\end{bmatrix}
\begin{bmatrix}
    -3\mu & \mu & \mu & \mu \\
    \mu & - 3\mu & \mu & \mu \\
    \mu & \mu & -3\mu & \mu \\
    \mu & \mu & \mu & -3\mu
\end{bmatrix} \\
&= \begin{bmatrix}
-3 \mu r_t + 3 \mu q_t & \mu r_t - \mu q_t & \mu r_t - \mu q_t & \mu r_t - \mu q_t \\
\mu r_t - \mu q_t & -3 \mu r_t + 3 \mu q_t & \mu r_t - \mu q_t & \mu r_t - \mu q_t  \\
\mu r_t - \mu q_t & \mu r_t - \mu q_t & -3 \mu r_t + 3 \mu q_t & \mu r_t - \mu q_t \\
\mu r_t - \mu q_t & \mu r_t - \mu q_t & \mu r_t - \mu q_t & -3 \mu r_t + 3 \mu q_t
\end{bmatrix}
\end{aligned}
\nonumber
$$
Put simply, we have the system:
$$
\frac{d}{dt}[r_t] = -3 \mu r_t + 3 \mu q_t \hspace{4mm} \text{and} \hspace{4mm} \frac{d}{dt}[q_t] = \mu r_t + \mu q_t
\nonumber
$$
If we further assume that $r_t = q_t = \frac{1}{4}$ as $t \rightarrow \infty$, then we get the unique solutions:
$$
r_t = \frac{1}{4}(1 + 3 \exp(-4 \mu t)) \hspace{4mm} \text{and} \hspace{4mm} q_t = \frac{1}{4}(1 - \exp(-4 \mu t))
\nonumber
$$
One last concept we need to introduce is <i>branch length</i>. In this context, it is a quantification of the difference between a parent and child node in a phylogenetic tree (i.e. an ancestor and descendant) based upon the number of substitutions between the two sequences. Usually, this is normalized to the sequence lengths and also a value representing the <i>expected</i> number or proportion of substitutions. 
</details>

The general idea is that the missing nucleotides are imputed as the value maximizing the posterior probability of the descendant sequence using the Jukes-Cantor model as the probabilistic framework. There are some more details involving computational efficiency and details. Briefly, they do some preprocessing of the phylogenetic tree, and they split the entire tree into three parts to improve computational feasibility. The result of this step is a multiple sequence alignment table that has no gaps (missing bases in non-reference sequences are imputed, and positions corresponding to gaps in the reference strain are removed) and can serve as the <i>reference database</i>.

### Estimation
An expectation-maximization (EM) algorithm is used to estimatem $q_{i,j}$ in Eq. \eqref{eq:q-defn}. Given $q_{i,j}$ for all $i, j$, we first initialize the relative haplotype proportions as $$p^0_i \sim \text{Unif}(0, 1)$$ for all $i \in [k]$ and divide by $$\sum_{i = 1}^k p^0_i$$ so the vector sums to $1$. Then for $t = 1, 2, \dots$ until convergence:

<ol>
    <li>
        Compute the relative haplotype proportions as $p^t_i = \frac{1}{n} \sum_{j = 1}^n \frac{p_i^{t-1}q_{i,j}}{\sum_{l = 1}^k p_l^{t-1}q_{i,j}}$.
    </li>
    <li>
        Compute the log-likelihood at step $t$ as $\ell^t(p_1, \dots, p_k) = \sum_{j = 1}^r \log\left( \sum_{i = 1}^k q_{i,j} p^t_{t} \right)$
    </li>
</ol>

There is an unidentifiability issue in the case that there exist haplotypes $i$ and $i'$ such that $q_{i,j} = q_{i', j}$ for all reads $j \in [n]$. To overcome this, the authors group the haplotypes with identical $q$ values across all reads into a single "haplotype group". 

---

## Unsupervised Detection and Fitness Estimation of Emerging SARS-CoV-2 Variants
Lefebvre et al.<d-cite key=lefebvre2025></d-cite> use a clustering-based approach to solve the variant deconvolution problem. By leveraging the fact that observations from the same sewer system are correlated over time, they learn which mutations belong to the same latent group based upon their frequencies. 

Since their method simultaneously estimates variant fitness, they use some terminology that is not really relevant to our discussion (e.g. variants under positive/negative selection). I'll just use it to keep this post consistent with the original manuscript.

### Model
Let $K \geq 0$ be the number of latent groups (i.e. variants) that are under some selective pressure (either positive or negative), which the authors call "non-neutral groups". The model always includes a neutral group that is under no selective pressure, and we'll use a group number of $0$ for this one. Let $G_k$ for $k \in \{ 0, \dots, K \}$ denote group number.

Let $n$ be the total (over all samples) number of mutations found in the dataset (compared to some reference sequence). We use $Z = \{ Z_1, \dots, Z_n \} \in \{ 0, \dots, K \}^n$ to be the group assignments for all mutations, where $Z_i$ for $i \in \{ 1, \dots, n\}$ is the group assignment of the $i$-th mutation. 

Let $t_0, \dots, t_m$ denote the timepoints at which the data are collected, and let $\mathcal{T} = (t_0 - t_0 = 0, t_1 - t_0, \dots, t_m - t_0 = T)$ be the time differences between the observations. We will use $X_{i, t}$ to denote the count for mutation $i$ at time $t \in \mathcal{T}$, and we'll use $d_{i, t}$ to denote the read depth for the genome position of mutation $i$ at time $t \in \mathcal{T}$. 

Let $\pi = \{ \pi_0, \dots, \pi_K \}$ such that $\sum_{i = 0}^K \pi_i = 1$ where $\pi_k$ is the probability of a mutation belonging to group $k$ for $k \in \{ 0, \dots, K \}$. We assume a multinomial distribution for the elements of $Z$. 

The authors choose a generalized linear model for the mutation counts conditional on non-neutral group:

$$
X_{i,t} \rvert Z_i = k, k \neq 0 \sim \text{Binom}\left(d_{i,t}, \text{logit}^{-1}(\exp(\mu_k + s_k t))\right)
$$

where $\mu_k$ and $s_k$ are the intercept and slope for group $k$ when $k$ is not the neutral group. The authors call $s_k$ a "selection coefficient" because a negative/positive value represents the decrease/increase in the mutation frequency. 

For the neutral group, they use a Beta-binomial distribution:

$$
X_{i,t} \rvert Z_i = 0 \sim \text{Binom}(d_{i,t}, u), 
\hspace{5mm}
u \sim \text{Beta}(\alpha, \beta)
$$

for some learned hyperparameters, $\alpha$ and $\beta$. There is no selective pressure, so it seems reasonable to allow $s_0 = 0$ in this way. However, I'm not convinced one needs to learn $\alpha$ and $\beta$ over just $\mu_0$...

The parameter vector is denoted by $$\theta = (\pi, \mu, s, \alpha, \beta)$$, where $$\pi = (\pi_0, \dots, \pi_K) \in [0, 1]^{K + 1}$$, $$\mu = (\mu_1, \dots, \mu_K) \in \mathbb{R}^K$$, $$s = (s_1, \dots, s_K) \in \mathbb{R}^K_{\geq 0}$$, and $$\alpha, \beta \in \mathbb{R}_{> 0}$$. One can then write the joint distribution of $$X = \{ X_1, \dots, X_n \}$$ with $$X_i = \{ X_{i,t} \}_{t = 0, \dots, T}$$ and $Z$, conditional on $\theta$, as:

$$
\mathbb{P}(X, Z \rvert \theta) = \prod_{i = 1}^n \prod_{t \in \mathcal{T}}  \mathbb{P}(Z_i \rvert \pi)\mathbb{P}\left(X_{i,t} \rvert Z_i; \mu, s, \alpha, \beta \right) 
$$

Lefebvre et al. use an expectation-maximization algorithm to estimate the parameter values. 

---

## Tracking SARS-CoV-2 Genomic Variants in Wastewater Sequencing Data with LolliPop

Dreifuss et al.<d-cite key=dreifuss2022></d-cite> use a regression-based framework to perform variant deconvolution. Their novelty comes from an additional regularization step that imposes a kind of "smoothness over time", making their method particularly suited to time series analysis. 

### Model
We assume to have $V$ possible variants, each of which is characterized by some subset of the mutation set, $S_M = \{ 1, \dots, M\}$, defined with respect to a single, fixed reference strain. The characteristic mutations are collected into a <i>profile</i> matrix, $\mathbf{X} \in \{0, 1\}^{M \times V}$, where $\mathbf{X}_{i,j} = 1$ if mutation $i \in S_M$ is characteristic of variant $j \in [V]$ and $0$ otherwise. 

We let $$\mathbf{y}_t \in [0, 1]^M$$ denote the <i>mutation frequency vector</i> at time $$t \in [T]$$ where $$\mathbf{y}_{t, i}$$ is the observed frequency of mutation $i$ at time $t$ (computed as the proportion of reads with mutation $i$ from a sample taken at time $t$). The goal is to estimate the <i>relative variant abundance vector</i>, $$\mathbf{b}_t \in [0, 1]^V$$, where $$\mathbf{b}_{t, j}$$ denotes the relative proportion of variant $j \in [V]$ in the sample. We enforce that $$\rvert \rvert \mathbf{b}_t \rvert \rvert_1 = 1$$. 

The model is a linear regression framework. Under the assumption that:

$$
\begin{equation}
\label{eq:model-lollipop}
\mathbb{E}\left[ \mathbf{y}_t \rvert \mathbf{b}_t \right] = \mathbf{X} \mathbf{b}_t; \hspace{5mm} \forall t \in [T]
\end{equation}
$$

we aim to estimate $\mathbf{b}_1, \dots, \mathbf{b}_T$ such that $\mathbf{y}_t = \mathbf{X} \mathbf{b}_t$ for all $t \in [T]$. The solution is taken as:

$$
\begin{equation}
\label{eq:model-soln}
\hat{\mathbf{b}}_t = \underset{\mathbf{b}_t \in [0, 1]^V}{\arg \min}  \left\{ \mathcal{L}(\mathbf{y}_t - \mathbf{X} \mathbf{b}_t) \right\}
\end{equation}
$$

where $\mathcal{L}$ is some loss function.

### Temporal Regularization
The temporal regularization is implemented with a <a href="https://web.stanford.edu/group/SOL/papers/fused-lasso-JRSSB.pdf">fused</a> ridge penalty (i.e. an $L_2$ penalty on the differences rather than the vector themselves) that enforces "temporal continuity" between the $\mathbf{b}_t$ over values of $t$. Intuitively, applying the penalty to the differences in abundance vectors enforces sparsity in the entries of the difference, which keeps the two vectors "close" from timepoint to timepoint.

Letting $$\lambda_{i,j} > 0$$ denote the penalty parameter for abundance vectors, $\mathbf{b}_i$ and $\mathbf{b}_j$, at times $i$ and $j$, the loss function is:

$$
\begin{equation}
\label{eq:model-loss}
\mathcal{L}(\mathbf{b}_1, \dots, \mathbf{b}_T; \lambda) = \sum_{i = 1}^T \rvert \rvert \mathbf{y}_i - \mathbf{X} \mathbf{b}_i \rvert \rvert_2^2 + \sum_{i = 1}^T \sum_{j = 1}^T \lambda_{i,j} \rvert \rvert \mathbf{X} (\mathbf{b}_i - \mathbf{b}_j) \rvert \rvert_2^2
\end{equation}
$$

Defining $\mu_i = 1 + \sum_{j = 1}^T \lambda_{i,j}$ and the <a href="https://en.wikipedia.org/wiki/Doubly_stochastic_matrix">doubly stochastic matrix</a>:

$$
\Lambda = \begin{bmatrix}
\mu_1 - \lambda_{1, 1} & \dots & -\lambda_{T, 1} \\
\vdots & \ddots & \vdots \\
-\lambda_{1, T} & \dots & \mu_T - \lambda_{T, T}
\end{bmatrix}
$$

<div class="theorem">
<strong>Claim (Solution).</strong>
{% tabs lambda-soln %}
{% tab lambda-soln statement %}
If $\Lambda$ is has full rank ($T$), then the solution can be found:

$$
\begin{equation}
\label{eq:model-soln-full-rank}
\hat{\mathbf{B}} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y} \Lambda^{-1}
\end{equation}
$$

where $\mathbf{Y}$ is a matrix whose columns correspond to $\mathbf{y}_1, \dots, \mathbf{y}_T$, and $\mathbf{B}$ is a matrix whose columns correspond to $\mathbf{b}_1, \dots, \mathbf{b}_T$. 
{% endtab %}
{% tab lambda-soln proof %}
Taking the gradient of Eq. \eqref{eq:model-loss} with respect to $\mathbf{b}_i$ gives:

$$
\begin{aligned}
\frac{\partial \mathcal{L}(\mathbf{b}_1, \dots, \mathbf{b}_T; \lambda)}{\partial \mathbf{b}_i} 
&= \frac{\partial}{\partial \mathbf{b}_i} \left[ \sum_{i = 1}^T \rvert \rvert \mathbf{y}_i - \mathbf{X} \mathbf{b}_i \rvert \rvert_2^2 + \sum_{i = 1}^T \sum_{j = 1}^T \lambda_{i,j} \rvert \rvert \mathbf{X} (\mathbf{b}_i - \mathbf{b}_j) \rvert \rvert_2^2 \right] \\
&= \frac{\partial}{\partial \mathbf{b}_i} \left[ \sum_{i = 1}^T (\mathbf{y}_i - \mathbf{X} \mathbf{b}_i)^\top (\mathbf{y}_i - \mathbf{X} \mathbf{b}_i) + \sum_{i = 1}^T \sum_{j = 1}^T \lambda_{i,j} (\mathbf{X} (\mathbf{b}_i - \mathbf{b}_j))^\top(\mathbf{X}(\mathbf{b}_i - \mathbf{b}_j)) \right] \\
&= \frac{\partial}{\partial \mathbf{b}_i} \left[ \sum_{i = 1}^T (\mathbf{y}_i^\top \mathbf{y}_i - \mathbf{y}_i^\top \mathbf{X} \mathbf{b}_i - \mathbf{b}_i^\top \mathbf{X}^\top \mathbf{y}_i + \mathbf{b}_i^\top \mathbf{X}^\top \mathbf{X} \mathbf{b}_i) + \sum_{i = 1}^T \sum_{j = 1}^T \lambda_{i,j} (\mathbf{b}_i - \mathbf{b}_j)^\top \mathbf{X}^\top \mathbf{X}(\mathbf{b}_i - \mathbf{b}_j) \right] \\
&= - 2\mathbf{X}^\top \mathbf{y}_i  + 2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i + \sum_{j = 1}^T 2\lambda_{i,j} \mathbf{X}^\top \mathbf{X}(\mathbf{b}_i - \mathbf{b}_j) \\
&= 2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i - 2\mathbf{X}^\top \mathbf{y}_i + 2 \mathbf{X}^\top \mathbf{X} \sum_{j = 1}^T \lambda_{i,j} (\mathbf{b}_i - \mathbf{b}_j) \\
&= 2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i - 2\mathbf{X}^\top \mathbf{y}_i + 2 \mathbf{X}^\top \mathbf{X} \left( \mathbf{b}_i \sum_{j = 1}^T \lambda_{i,j}  - \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j  \right) \\
&= 2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i - 2\mathbf{X}^\top \mathbf{y}_i + 2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i \sum_{j = 1}^T \lambda_{i,j}  - 2 \mathbf{X}^\top \mathbf{X} \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j
\end{aligned}
$$

Setting the above equal to $0$ and solving for $\mathbf{b}_i$:

$$
\begin{aligned}
&\frac{\partial \mathcal{L}(\mathbf{b}_1, \dots, \mathbf{b}_T; \lambda)}{\partial \mathbf{b}_i} = 0 \\
\implies &2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i -  2\mathbf{X}^\top \mathbf{y}_i + 2 \mathbf{X}^\top \mathbf{X} \mathbf{b}_i \sum_{j = 1}^T \lambda_{i,j}  - 2 \mathbf{X}^\top \mathbf{X} \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j = 0 \\
\implies & \mathbf{X}^\top \mathbf{y}_i+ \mathbf{X}^\top \mathbf{X} \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j = \mathbf{X}^\top \mathbf{X} \mathbf{b}_i + \mathbf{X}^\top \mathbf{X} \mathbf{b}_i \sum_{j = 1}^T \lambda_{i,j} \\
\implies & (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_i + \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j = \mathbf{b}_i \left(1 + \sum_{j = 1}^T \lambda_{i,j}\right) \\
\implies & (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_i + \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j = \mu_i \mathbf{b}_i & \left(\mu_i = 1 + \sum_{j = 1}^T \lambda_{i,j}\right)\\
\implies & (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_i = \mu_i \mathbf{b}_i - \sum_{j = 1}^T \lambda_{i,j} \mathbf{b}_j \\
\implies & (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_i = (\mu_i - \lambda_{i,i}) \mathbf{b}_i - \sum_{j = 1 \\ j \neq i}^T \lambda_{i,j} \mathbf{b}_j \\
&=(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_i = (\mu_i - \lambda_{i,i}) \mathbf{b}_i - \sum_{j = 1 \\ j \neq i}^T \lambda_{i,j} \mathbf{b}_j 
\end{aligned}
$$

where the substitution of $\mu_i$ relies upon the assumption that $\lambda_{i,j} = \lambda_{j,i}$ for all $i,j$. Organizing $\mathbf{b}_1, \dots, \mathbf{b}_T$ as the columns of a matrix, $\mathbf{B}$, and $\mathbf{y}_1, \dots, \mathbf{y}_T$ into the columns of a matrix, $\mathbf{Y}$, we can rewrite the system of equations for $i \in [T]$ as:

$$
\begin{aligned}
\begin{Bmatrix}
(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_1 = (\mu_1 - \lambda_{1, 1}) \mathbf{b}_1 - \sum_{j = 1 \\ j \neq i}^T \lambda_{1,j} \mathbf{b}_j \\
(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_2 = (\mu_2 - \lambda_{2, 2}) \mathbf{b}_2 - \sum_{j = 1 \\ j \neq i}^T \lambda_{2,j} \mathbf{b}_j \\
\vdots \\
(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_T = (\mu_T - \lambda_{T, T}) \mathbf{b}_T - \sum_{j = 1 \\ j \neq i}^T \lambda_{T,j} \mathbf{b}_j
\end{Bmatrix}
&= 
\begin{Bmatrix}
(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_1 = (\mu_1 - \lambda_{1, 1}) \mathbf{b}_1 - \lambda_{1, 2} \mathbf{b}_2 - \dots - \lambda_{1, T} \mathbf{b}_T \\
(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_2 = -\lambda_{2, 1} \mathbf{b}_1 + (\mu_2 - \lambda_{2, 2}) \mathbf{b}_2 - \dots - \lambda_{2, T} \mathbf{b}_T \\
\vdots \\
(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{y}_T = -\lambda_{T, 1} \mathbf{b}_1 - \lambda_{T,2} \mathbf{b}_{2} - \dots + (\mu_T - \lambda_{T, T}) \mathbf{b}_T 
\end{Bmatrix}
\end{aligned}
$$

And the righthand side of the above is equivalent to:

$$
\begin{aligned}
&(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \begin{bmatrix} \rvert & & \rvert \\ \mathbf{y}_1 & \dots & \mathbf{y}_T \\ \rvert & & \rvert \end{bmatrix}
= \begin{bmatrix} 
\rvert & & \rvert \\
\mathbf{b}_1 & \dots & \mathbf{b}_T\\
\rvert &  & \rvert
\end{bmatrix}
\begin{bmatrix}
\mu_1 - \lambda_{1, 1} & \dots & -\lambda_{T, 1} \\
\vdots & \ddots & \vdots \\
-\lambda_{1, T} & \dots & \mu_T - \lambda_{T, T}
\end{bmatrix}  \\
\implies &(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y}
= \mathbf{B} \Lambda \\
\implies &(\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y} \Lambda^{-1} = \mathbf{B}
\end{aligned}
$$
{% endtab %}
{% endtabs %}
</div>

Note that the above can easily be generalized by instead specifying a (non-negative, non-increasing) kernel smoother function of $i - j$, $k(i, j)$, instead of scalar penalty terms $\lambda_{i,j}$. Let $\bar{\mathbf{K}}$ denote the matrix where:

$$
\bar{\mathbf{K}}_{i,j} = \frac{k(i,j)}{\sum_{l = 1}^T k(i, l)}
$$

For example, since the rows and columns of $$\Lambda^{-1}$$ sum to $1$ (because $\Lambda$ is doubly stochastic with non-negative entries), we can define $$\bar{\mathbf{K}}_{i,j} = (\Lambda^{-1})_{i,j}$$ as one such kernel matrix. We can then see that Eq. \eqref{eq:model-soln-full-rank} is equivalent to:

$$
\hat{\mathbf{B}} = (\mathbf{X}^\top \mathbf{X})^{-1} \mathbf{X}^\top \mathbf{Y} \bar{\mathbf{K}}
$$

which itself is equivalent to the least squares solution for:

$$
\mathbb{E}\left[ \mathbf{Y} \bar{\mathbf{K}} \rvert \mathbf{B} \right] = \mathbf{X} \mathbf{B}
\hspace{5mm} \iff \hspace{5mm}
\frac{1}{\sum_{l = 1}^T k(i, l)} \mathbb{E}\left[\mathbf{Y} \rvert \mathbf{b}_i \right]\begin{bmatrix} k(i, 1) \\ \vdots \\ k(i, T) \end{bmatrix} = \mathbf{X} \mathbf{b}_i \hspace{2mm} \forall i \in [T]
$$

In <i>Lollipop</i>, they choose the Gaussian kernel:

$$
k_{G}(i, j; \kappa) = \exp\left( \frac{-(i - j)^2}{2 \kappa}\right); \hspace{8mm} \kappa > 0
$$

where $\kappa$ must be tuned.

### Uncertainty Quantification
The authors provide two methods of getting confidence intervals for $\mathbf{B}$: a bootstrap method and an asymptotic method. The bootstrap method is essentially just resampling of the row indices of $\mathbf{Y}$ with replacement $B$ times. The deconvolution procedure is performed for all of the bootstrap samples, which results in $B$ estimates of $\mathbf{b}_t$ for some time point $t$. The empirical quantiles of the components of estimated variance abundance vector $\mathbf{b}_t$ can be used to make confidence intervals. 

Asymptotic confidence intervals are computed by assuming some parametric distribution for the responses. A simple version is to assume that for timepoint $t$, the proportion of reads for mutation $j$:

$$
n \mathbf{y}_{t, j} \sim \text{Binom}(n, \pi_j)
$$

where $n$ is the total number of reads in the sample. One can then easily derive an $\alpha$-level confidence intervals (assuming independence of reads) from the fact that:

$$
\frac{\hat{\mathbf{b}}_i - \mathbf{b}_i}{\sqrt{\text{Var}(\hat{\mathbf{b}}_i)}} \overset{asymp.}{\sim} \mathcal{N}(\mathbf{0}_V, \mathbb{I}_{V \times V}) 
$$

