# Differentially regulated gene network identification

## TODO/ideas

- Make subnetworks from B more realistic than band
- How are the candidate subnetworks actually identified???
  - "We considered $p = 1000$ and $μ_{10} = (0.3, 0.3, . . ., 0.3)^T$, and the gene regulatory networks were estimated using the lasso method, where the response and predictor variables are the expression levels of the target gene and regulator genes, respectively."
  - "We consider the gene regulatory network structure based on not only expression levels of genes but also edge structures and regulatory effects of regulator
genes to their target genes and incorporate the comprehensive information into a dissimilarity measure to compare the gene regulatory networks."
  - In the monte carlo setting they may just be looking at the selected subnetworks themselves??


## Background

- Task: Differentially regulated gene network identification
- Undirected graphs modeling conditional independence (TODO: review the proper name for this)
- Gaussian graphical models
- Bulk cell line gene networks
  - aggregates data across all cells in a cell line
  - identifies universal processes for the cell line
  - Differential analysis focuses on differential regulation but not variation in the graph structure
- Cell line characteristic-specific gene networks
  - separate cells with specific characteristics (ex: bulk RNA-seq data collected separately for sensitive vs resistant)
- Precision matrix: inverse of covariance matrix
  - if $\Omega_{i,j}$ is non-zero, these genes are dependent conditioned on all other genes
- SiGN-L1: infers gene networks via sparse precision matrix estimation (variation of Lasso)
- Huge.generator: Generates multivariate Gaussian data with various graph structures

## Monte Carlo simulations

As describedin CIdrgn

- Two phenotypes A and B
  - 50 samples of each
- 10 subnetworks, each consisting of 10 genes
  - Also did a set with 50 genes each
- 4 common subnetworks
- 6 A-specific
- 6 B-specific
- 4 scenarios
  - 1: "random" graph structure
  - 2: "cluster" graph structure
  - 3: "scale-free" graph structure
  - 4: "hub" graph structure
- $p = 100, \mu_{10} = [0.3, 0.3, ... 0.3]^T$
- Then added in cell-line specific characteristics (variance across samples) -- TODO: understand this better
  - sample-specific gene regulatory strength represents biological diversity in activation/suppression effect
  - modulator adjusts gene network structure
- Idea: do things other than just band for $\Omega_{B_i}$

Generating gene expression for 10 genes in each of 4 common subnetworks:

- 4 precision matrices $\Omega_{A,i} = \Omega_{B,i}$ for common subnetworks from "random" graph structure using huge.generator in R package Huge
  - off-diagaol elements set to 0.5 (controls magnitude of partial correlations)
  - 0.2 added to diagonal elements of precision matrix (controls magnitude of partial correlations)
  - Gene expression for the 100 cell lines sampled from $N(0_{10}, \Omega_{C,i}^{-1})$
  
Generating gene expression for 10 genes in each of 6 differential subnetworks:

- 6 precision matrices $\Omega_{A,i}$ generated from "random" graph structure
  - Gene expression for 50 cell lines sampled from $N(0_{10}, \Omega_{A,i}^{-1})$
- 6 precision matrices $\Omega_{B,i}$ generated from "band" graph structures
  - Gene expression for 50 cell lines sampled from $N(\mu_{10}, \Omega_{B,i}^{-1})$

Generating gene expression for the remaining $p-100$ genes

- Gene expression for $p-100$ genes sampled from $N(0_{p-100}, I_p)$

## Gene regulatory network estimation

- Lasso method
  - predictor variables $x^i$ are the expression levels of regulator genes
  - response variables $y_i$ are expression levels of each target gene
  - TODO: How are the regulator and response variables separated in the Monte-Carlo simulation??

- Cell line characteristic-specific gene network estimated using SiGN-L1

## Differential subnetwork identification

- Statistics compared
  - CIdrgn
  - SAM-GS
  - GSCA
- permutation $p$-values computed and responsive subnetworks were those with $p < 0.05$
- Metrics used for comparison
  - recall
  - precision
  - TNR
  - F-measure
  - accuracy
  - TODO: compare with their results in column “subnetwork size:10” of Table 1
- Findings
  - They find CIdrgn improves TNR of four common subnetworks
  - CIdrgn effective for gene network identification in bulk cell line gene networks and cell line characteristic-specific gene networks TODO: understand the difference

## Statistics

- SAM-GS statistic
  - $D_{SAM-GS}$ measures the difference in expression levels of gene $j$ in phenotypes $A$ and $B$
  - $D_{SAM-GS} = \sum_{j \in V} \frac{(\bar{x}_{Aj} - \bar{x}_{Bj})^2}{s_j + s_0}$
  - $s_0$: tuning parameter (Parameter 3.3 comes from SAM paper. They chose to minimize the coefficient of variation)
  - $s_j = \sqrt{a \cdot \left( \sum_{i=1}^{n_A} (x_{ij} - \bar{x}_{Aj})^2 \cdot \sum_{k=1}^{n_B} (x_{kj} - \bar{x}_{Bj})^2 \right) }$
  - $a = \frac{1/n_A + 1/n_B}{n_A + n_B -2}$

- Gene set co-expression analysis
  - statistic to identify differentially co-expressed genes
  - computes pairwise correlations for all gene pairs
  - $D_{GSCA}$ measures the dispersion of correlations between phenotypes A and B
  - $D_{GSCA} = \sqrt{\frac{1}{|V|(|V|-1)/2} \sum_{k=2}^{|V|}\sum_{j=1}^{k-1}(C_{kj}^A - C_{kj}^B)^2}$
  - $c_{kj}^A$ is the Pearson correlation between the $k$th and $j$th genes

## References

- [CIdrgn paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10446197/pdf/pone.0286044.pdf)
- [CIdrgn simulation results](https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0286044.t001)
- [SiGN-L1 code](https://sign.hgc.jp/signl1/)
- [SiGN-L1 paper 1](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-3-41)
- [SiGN-L1 paper 2](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0020804&type=printable)
- [LASSO](https://webdoc.agsci.colostate.edu/koontz/arec-econ535/papers/Tibshirani%20(JRSS-B%201996).pdf)
- [Huge R package](https://cran.r-project.org/web/packages/huge/vignettes/vignette.pdf)
- [Statistical test methods paper](https://www.nature.com/articles/s41598-019-47362-7)
- [SAM statistic](https://www.pnas.org/doi/full/10.1073/pnas.091062498)
- [SAM-GS statistic](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-242)
- [SAM-GS excel add-in](https://sites.ualberta.ca/~yyasui/SAM-GS/SAM-GS%20Documentation.pdf)
- [GSCA statistic](https://academic.oup.com/bioinformatics/article/25/21/2780/226874)
- [GSCA RData](https://www.biostat.wisc.edu/~kendzior/GSCA/GSCA.RData)