---
layout: distill
title: Cellular Deconvolution
description: A Primer
date: 2025-11-06
tabs: true
tags: biology deconvolution sequencing primer
toc:
  - name: Background
    subsections:
        - name: Reference-Based Methods
        - name: Reference-Free Methods
  - name: Regression Framework
  - name: Probabilistic Framework
  - name: Dimension Reduction Framework
    subsections: 
        - name: Linseed
  - name: Future Work
    subsections:
        - name: Mutually Exclusive Marker Genes
        - name: Fixed Number of Cell Types
        - name: Confounding
        - name: Evaluation
bibliography: 2025-11-06-cellular-deconvolution.bib
---

Cellular deconvolution is an important analytical step for experiments that rely on samples that could comprise several different cell types (e.g. bulk RNA-seq, spatial genomics, etc.). Though these experiments are more cost effective compared to single cell RNA-seq (scRNA-seq), the gene expression data that result are an average over all of the cells in the sample, and the expression profiles of cells of interest that are scarce may be masked by those are that more abundant. At its core, cellular deconvolution is the problem of teasing apart the gene expression patterns of each cell type from bulk data. 

In this post, I'm going to reviewing the problem of cellular deconvolution and some of the work done in that area. I'll use the case of bulk RNA-seq for my discussions, but keep in mind that the problem is much more general than this particular application.

---
 
## Background
Deconvolution methods can be divided into two main types: <i>reference-based</i> (also called <i>partial</i> deconvolution methods) and <i>reference-free</i> (also called <i>complete</i> deconvolution methods).

### Reference-Based Methods
A majority of the more popular/recent cellular deconvolution methods are reference-based, which makes sense because of the rise of single cell sequencing. In reference-based methods, single cell data (the <i>reference</i> dataset) are required as inputs, and are used to help define the expression signature of each cell type or group. 

This often requires defining a set of (sometimes mutually exclusive) marker genes for each cell type, which can be done using prior knowledge like a literature review. However, if we know the cell types for each cell in our reference dataset, we can use it to identify the marker genes with some sort of statistical comparison. For example, the average expression of gene $i$ in cell type $j$ can be compared to the average of gene $i$ across all other cells via a $t$-test or likelihood ratio test. If the expression is statistically significantly different, then gene $i$ can be used as a marker for cell type $j$. 

Beyond picking marker genes, the reference dataset is used to compute the expression signatures by averaging the gene expression levels within each cell type. If the cell types are unknown in the reference dataset, then an unsupervised clustering step can be done first to group by cell type.

### Reference-Free Methods
In contrast to reference-based procedures, reference-free methods do not require any additional data. From the bulk data itself, marker genes are identified and the expression signatures are estimated, usually with unsupervised machine learning techniques.

One way to select the marker genes is to project the bulk gene expression data into an alternative (often lower-dimensional) space and choose the genes that are close enough to the corners of the simplex in which all of the data fall (<i>LinSeed</i><d-cite key=zaitsev2019></d-cite>). Another way is to treat the samples as documents and use existing topic modelling algorithms (<i>CellDistinguisher</i><d-cite key=newberg2018></d-cite>). 

---

## Regression Framework
A common way that cellular deconvolution methods approach the problem is to frame the average gene expression obtained via sequencing as a linear combination of the gene expression of each cell type. 

We'll assume to measure the expression levels of $n$ marker genes from $m$ tissue samples containing cells of $k$ types. The measured gene expression for $n$ genes will be collected in the matrix $$\mathbf{G} \in \mathbb{R}^{n \times m}_{\geq 0}$$ where $\mathbf{G}_{i,l}$ denotes the average gene expression of gene $i$ over all cells in sample $l$. We'll let $$\mathbf{S} \in \mathbb{R}^{n \times k}_{\geq 0}$$ denote the <i>signature</i> matrix, where $$\mathbf{S}_{i,j}$$ denotes the gene expression of gene $i$ in cell type $j$. We'll also have the matrix $$\mathbf{P} \in \mathbb{R}^{k \times m}_{\geq 0}$$, where $$\mathbf{P}_{j, l}$$ denotes the proportion of cell type $j$ in sample $l$. We enforce that $$\sum_{j = 1}^{k} \mathbf{p}_{j,l} = 1$$ for all $l$.

The linear combination framework aims to estimate $\mathbf{P}$ under the assumption that the following model holds:

$$
\begin{equation}
\label{eq:lin-comb-framework}
\mathbf{G} = \mathbf{S} \mathbf{P}
\end{equation}
$$

We can then find the cell type proportions as the solution to:

$$
\begin{equation}
\label{eq:reference-based}
\underset{\mathbf{P} \in \mathbb{R}^{k,m}_{\geq 0} \text{: } \sum_{j = 1}^k \mathbf{P}_{j,\cdot} = 1}{\arg \min}
\left\{ \mathcal{L}(\mathbf{G}, \mathbf{S}\mathbf{P}) \right\}
\end{equation}
$$

where $\mathbf{L}$ is some loss or distance function (e.g. constrained least squares, weight constrained least squares, etc.<d-cite key=nguyen2024></d-cite>). 

---

## Probabilistic Framework
An alternative to the linear combination framework is to fit a Bayesian or hierarchical model to the data. 

<i>Section to be completed</i>.

---

## Dimension Reduction Framework
Dimension reduction frameworks are often reference-free and generally rely upon <a href="https://en.wikipedia.org/wiki/Matrix_decomposition">(non-negative) matrix factorization</a> of $\mathbf{G}$ or <a href="https://en.wikipedia.org/wiki/Independent_component_analysis">independent component analysis</a>. At a very high level, a lower dimensional space may have dimensions corresponding to cell types (compared to the original space that had genes as its dimensions). If the subspace is learned, then this could lead to the discovery of new cell types or cell states.

As is characteristic of this kind of method, the cellular deconvolution via dimensionality reduction is usually very computationally efficient, and the results are often easily interpretable (especially in the case of independent component analysis). However, there is no guarantee that the dimensions of the subspace correspond to cell types/states, and these methods may perform poorly in cases where the data are very noisy or the cell types in the sample are highly similar.<d-cite key=gaspard2025></d-cite>


### Linseed
Linseed (linear subspace identification for gene expression deconvolution)<d-cite key=zaitsev2019></d-cite> is a dimension reduction framework that is entirely reference-free. It hinges upon the idea of <i>mutual linearity</i>, which is a property of genes that are markers for the same cell type. 

Linseed also assumes a true model that is linear:

$$
\mathbf{G} = \mathbf{W} \mathbf{P}
$$

where $\mathbf{W} \in \mathbb{R}_{\geq 0}^{n \times k}$ is an signature matrix that specifies the gene expression in each cell type. 


#### Mutual Linearity
Suppose we have $n$ samples of some cell type $a$ that has three true marker genes, $x$, $y$, and $z$. Let $x_j$, $y_j$, and $z_j$ denote the expression levels of these genes. Mutual linearity assumes that there exist corresponding proportionality constants, $k_{x,y}$ and $k_{x, z}$, such that the following holds:

$$
\mathbf{G} = 
\begin{bmatrix}
x_1 & x_2 & \dots & x_n \\
y_1 & y_2 & \dots & y_n \\
z_1 & z_2 & \dots & z_n
\end{bmatrix} =
\begin{bmatrix}
x_1 & x_2 & \dots & x_n \\
k_{x,y} x_1 & k_{x,y} x_2 & \dots & k_{x,y} x_n \\
k_{x,z} x_1 & k_{x,z} x_2 & \dots & k_{x,z} x_n
\end{bmatrix}
$$

#### Row Normalization
The first step of <i>Linseed</i> is to normalize by dividing each entry by the sum of its row:

$$
\tilde{\mathbf{G}} = 
\begin{bmatrix}
\frac{x_1}{\sum_{j = 1}^n x_j} & \frac{x_2}{\sum_{j = 1}^n x_j} & \dots & \frac{x_n}{\sum_{j = 1}^n x_j} \\
\frac{y_1}{\sum_{j = 1}^n y_j} & \frac{y_2}{\sum_{j = 1}^n y_j} & \dots & \frac{y_n}{\sum_{j = 1}^n y_j} \\
\frac{z_1}{\sum_{j = 1}^n z_j} & \frac{z_2}{\sum_{j = 1}^n z_j} & \dots & \frac{z_n}{\sum_{j = 1}^n z_j} \\
\end{bmatrix} =
\begin{bmatrix}
\frac{x_1}{\sum_{j = 1}^n x_j} & \frac{x_2}{\sum_{j = 1}^n x_j} & \dots & \frac{x_n}{\sum_{j = 1}^n x_j} \\
\frac{k_{x,y} x_1}{\sum_{j = 1}^n k_{x,y} x_j } & \frac{k_{x,y}x_2}{\sum_{j = 1}^n k_{x,y} x_j } & \dots & \frac{k_{x,y}x_n}{\sum_{j = 1}^n k_{x,y} x_j } \\
\frac{k_{x,z}x_1}{\sum_{j = 1}^n k_{x,z} x_j} & \frac{k_{x,z}x_2}{\sum_{j = 1}^n k_{x,z} x_j} & \dots & \frac{k_{x,z}x_n}{\sum_{j = 1}^n k_{x,z} x_j}
\end{bmatrix} =
\begin{bmatrix}
\frac{x_1}{\sum_{j = 1}^n x_j} & \frac{x_2}{\sum_{j = 1}^n x_j} & \dots & \frac{x_n}{\sum_{j = 1}^n x_j} \\
\frac{x_1}{\sum_{j = 1}^n x_j} & \frac{x_2}{\sum_{j = 1}^n x_j} & \dots & \frac{x_n}{\sum_{j = 1}^n x_j} \\
\frac{x_1}{\sum_{j = 1}^n x_j} & \frac{x_2}{\sum_{j = 1}^n x_j} & \dots & \frac{x_n}{\sum_{j = 1}^n x_j}
\end{bmatrix}
$$

This implies that all genes have a proportionality constant of $\tilde{k} = 1$ when we consider $\tilde{\mathbf{G}}$. 

#### Simplex

Transposing both of these matrices yields $\mathbf{G}^\top$ and $\mathbf{P}^\top$, which can be thought of as collections of $g$ and $c$ vectors, respectively, in an $m$-dimensional space defined by the samples. These matrices are then column-normalized (which is equivalent to row-normalizing the original matrices and then transposing) to get $\tilde{\mathbf{G}}^\top$ and $\tilde{\mathbf{P}}^\top$.

<div class="theorem">
<strong>Claim (Simplex).</strong>
{% tabs simplex-1 %}
{% tab simplex-1 statement %}
It can be shown that the columns of $\tilde{\mathbf{G}}$ lie in the $c - 1$ simplex generated by the columns of $\tilde{\mathbf{P}}$. 
{% endtab %}
{% tab simplex-1 proof %}
Under the assumption of a true linear model:

$$
\begin{aligned}
\mathbf{G}_{i,\cdot} &= \mathbf{W}_{i, \cdot} \mathbf{P} = \sum_{j = 1}^c \mathbf{W}_{i,j} \mathbf{P}_{j, \cdot} \\
\implies \tilde{\mathbf{G}}_{i,j} &= \frac{\sum_{j = 1}^c \mathbf{W}_{i,j} \mathbf{P}_{j, \cdot}}{\sum_{l = 1}^m \mathbf{G}_{i,l}}  & \left(\text{normalize } \mathbf{G} \right) \\
&=  \frac{\sum_{j = 1}^c \mathbf{W}_{i,j} \left( \sum_{l = 1}^m \mathbf{P}_{j, l}\right) \tilde{\mathbf{P}}_{j, \cdot}}{\sum_{l = 1}^m \mathbf{G}_{i,l}} & \left( \text{normalize } \mathbf{P} \right) \\
&=  \sum_{j = 1}^c  \alpha_{i,j} \tilde{\mathbf{P}}_{j, \cdot} & \left(\frac{\mathbf{W}_{i,j} \left( \sum_{l = 1}^m \mathbf{P}_{j, l}\right)}{\sum_{l = 1}^m \mathbf{G}_{i,l}} = \alpha_{i,j} \right)
\end{aligned}
$$

Now notice that $$\sum_{h = 1}^m \tilde{\mathbf{G}}_{i, h} = 1$$ and $$\sum_{h = 1}^m \tilde{\mathbf{P}}_{l, h} = 1$$ because the matrices are row normalized. Thus:

$$
\begin{aligned}
1 &= \sum_{h = 1}^m \tilde{\mathbf{G}}_{i, h} \\
&= \sum_{h = 1}^m \sum_{j = 1}^c  \alpha_{i,j} \tilde{\mathbf{P}}_{j, h} \\
&=  \sum_{j = 1}^c  \alpha_{i,j} \underbrace{\sum_{h = 1}^m \tilde{\mathbf{P}}_{j, h}}_{=1} \\ 
&= \sum_{j = 1}^c \alpha_{i,j}
\end{aligned}
$$

Since all of the matrix entries are non-negative, and the weights given by $\alpha$ sum to $1$, we have that the rows of $\tilde{\mathbf{G}}$ can be written as convex combinations of the rows of $\tilde{\mathbf{P}}$. This itself implies that the columns of $\tilde{\mathbf{G}}^\top$ can be written as convex combinations of the columns of $\tilde{\mathbf{P}}^\top$, and therefore that the columns of $\tilde{\mathbf{G}}$ lie in the $c - 1$ simplex generated by the columns of $\tilde{\mathbf{P}}$. 
{% endtab %}
{% endtabs %}
</div>

Furthermore, we have the followign result:

<div class="theorem">
<strong>Claim (Simplex Corners).</strong>
{% tabs simplex-2 %}
{% tab simplex-2 statement %}
The corners of the simplex correspond to the signature genes for each cell type. In other words, the normalized expression vector for gene $i$ (i.e. the $i$-th column of $\tilde{\mathbf{G}}$) will exactly equal one of the columns of $\tilde{\mathbf{P}}$ (specifically the $j$-th one).
{% endtab %}
{% tab simplex-2 proof %}
The proof is similar to the prior one. Let gene $i$ be a signature gene for cell type $j$. This means that only the part of each sample that corresponds to that cell type (i.e. the $j$-th row of $\mathbf{H}$) contributes to the expression of gene $i$ (i.e. the $i$-th row of $\mathbf{G}$), and only the $i$-th entry of $\mathbf{W}_{\cdot, j}$ is non-zero. That is:

$$
\begin{aligned}
\mathbf{G}_{i, \cdot} &= \mathbf{W}_{i,j} \mathbf{P}_{j, \cdot} \\
\implies \tilde{\mathbf{G}}_{i,\cdot} &= \frac{\mathbf{W}_{i,j} \mathbf{P}_{j, \cdot}}{\sum_{l = 1}^m \mathbf{G}_{i, l}} \\
&= \frac{\mathbf{W}_{i,j} \mathbf{P}_{j, \cdot}}{\sum_{l = 1}^m \mathbf{W}_{i,j} \mathbf{P}_{j, l}} \\
&= \frac{\mathbf{W}_{i,j} \mathbf{P}_{j, \cdot}}{\mathbf{W}_{i,j} \sum_{l = 1}^m \mathbf{P}_{j, l}} \\
&=  \frac{ \mathbf{P}_{j, \cdot}}{\sum_{l = 1}^m \mathbf{P}_{j, l}} \\
&= \tilde{\mathbf{P}}_{j, \cdot}
\end{aligned}
$$

The above implies that $\tilde{\mathbf{G}}^\top_{\cdot, i} = \tilde{\mathbf{P}}^\top_{\cdot, j}$. 
{% endtab %}
{% endtabs %}
</div>

---

## Future Work
Though cellular deconvolution has developed rapidly recently, there are still open problems.<d-cite key=nguyen2024></d-cite>

#### Mutually Exclusive Marker Genes
Most methods assume that the defining genes for each cell types are not expressed in other cell types (mutually exclusive), which does not always hold. For example, subtypes of a given cell type can have overlapping markers due to their shared lineage. Additionally, if the average expression levels in the bulk sample are very similar to a particular cell type, then that cell type's proportion may be overestimated. 

#### Fixed Number Of Cell Types
All methods assume that the number of cell types is known. Thus, estimates may be biased when a cell type is present in the sample but not represented in $\mathbf{S}$. 

#### Confounding
Most methods also fail to account for confounding effects during estimation. For example, laboratory protocols may introduce bias (e.g. sequencing errors, batch effects). In addition, cells are neither independent nor static. In the case of the immune system, the expression levels of certain cells is affected by the presence and signaling of other cells. Moreover, cells undergo state changes and differentiation that are often defined by changes in gene expression. 

#### Evaluation
Cellular deconvolution methods are quite difficult to validate due to the lack of a clear ground truth. However, this is not a problem unique to cellular deconvolution, and there are plenty proposed solutions. One example is taking two samples: one for bulk RNA-seq and one for scRNA-seq.

