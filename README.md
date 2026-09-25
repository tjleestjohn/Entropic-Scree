# The Entropic Scree:<br>An Information-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems

###### Initial Methods & Function Release: August 16, 2026 (Happy Birthday, Dad)

*[Terrence J. Lee-St. John, PhD](mailto:terry@enli.com.au)*  
*[Enli: Predictive systems that remain stable under change](https://www.enli.com.au)*

[![Read Preprint](https://img.shields.io/badge/Read_Preprint-Zenodo-blue?style=for-the-badge)](https://doi.org/10.5281/zenodo.22028086).
[![CRAN](https://img.shields.io/badge/CRAN-v1.0.1-blue?style=for-the-badge)](https://CRAN.R-project.org/package=Entropic.Scree)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg?style=for-the-badge)](https://github.com/tjleestjohn/Entropic-Scree/blob/main/LICENSE)

**[Install the R package](#r-package-installation) · [Understand the information map](#an-information-map-of-measurements) · [Run the simulations](#simulation-code)**

> **At a glance**
>
> A small number of generating processes can produce many nonlinear, thresholded and categorical measurements. Those measurements may require many directions for linear representation, even when they share a much smaller generating system.
>
> **Entropic Scree maps relationships among measurements using normalized mutual information.** Its spectrum supports selection of a primary informational rank and an Extended Signal Tail; its loadings describe bipolar groups of observed measurements; and its gravity metrics summarize the allocation of informational volume.
>
> In the paper's main synthetic experiment, the selected primary boundary matched the 20 generating coordinates, whereas the evaluated PCA, rank-PCA and RBF Kernel PCA diagnostics did not give the same selection. The framework combines this dimensional diagnostic with an interpretable map of the measurements' informational organization.

## R package installation

The published implementation is **`Entropic.Scree` 1.0.1**, available through [CRAN](https://CRAN.R-project.org/package=Entropic.Scree). It provides a compiled Rcpp/C++ backend with OpenMP parallelism where supported by the build and platform.

```r
install.packages("Entropic.Scree")
library(Entropic.Scree)
library(data.table)

# Supply observations in rows and measurements in columns.
# The input must be a data.table, not a correlation matrix.
dt <- as.data.table(your_dataset)

# Encode nominal measurements as factors or character columns.
# For example: dt[, diagnosis := factor(diagnosis)]

results <- Entropic.Scree(
  dt,
  interactive_mode = TRUE,
  extract_bipolar_modules = TRUE
)

# Review the selected boundaries, allocations and pole coefficients.
results$K_roots
results$K_extended
results$AIG
results$FSIG_final
results$bipolar_modules

# Optional: revise boundaries after inspecting the spectrum.
# Choose values appropriate to the spectrum of your own dataset.
# updated_results <- Update.Entropic.Scree(
#   results, new_K_roots = 3, new_K_extended = 12
# )
```

The interactive interface permits review and revision of the automatic recommendations. Set `interactive_mode = FALSE` for unattended execution. `Update.Entropic.Scree()` recalculates rank-dependent summaries without repeating pairwise MI estimation or eigendecomposition.

See the [package reference manual](https://stat.ethz.ch/CRAN/web/packages/Entropic.Scree/refman/Entropic.Scree.html) for arguments and returned fields.

### Preparing mixed data

Continuous measurements are discretized into equal-frequency categories. The default target is

$$B=\max\{2,\lceil cN^{1/3}\rceil\},$$

where `bin_multiplier` supplies $c$ and defaults to one; `num_bins` can specify the target directly. Numeric columns with no more than the target number of distinct values are encoded by state. Nominal numeric labels should be explicitly represented as categorical columns.

Version 1.0.1 requires complete prepared records and a bin target at least as large as every categorical column's number of encoded states. If necessary, increase `num_bins`; this also changes the resolution target for continuous measurements. Resolve missing values deliberately before analysis.

The package checks for constant and duplicate columns and, where applicable, linear collinearity. Identical discretized columns and columns below the marginal-entropy threshold are also removed. These later checks remain active when the initial purge and collinearity options are disabled. In the formulas below, $m$ is the **retained** measurement count.

## Python package — coming soon

A native Python package is in development. Installation instructions and a PyPI link will be added here when it is released. In the meantime, the published R package is available through CRAN.

## Why distinguish linear dimension from generating dimension?

PCA describes linear representational structure. Its limitations become important when its selected dimension is interpreted as a count of generating mechanisms:

1. **Mixed measurement forms:** A continuous measurement and its deterministic threshold need not have unit Pearson correlation. Their linear association depends on the threshold and marginal distribution.
2. **Nonlinear dependence:** Zero covariance does not imply independence. For symmetric $X$ with suitable moments, $X$ and $X^2$ have zero covariance despite their deterministic relationship.
3. **The sample-covariance rank ceiling:** A centered $N\times m$ data matrix has covariance rank at most $\min(m,N-1)$.
4. **Linear representational expansion:** Different functions of one generating variable can span several linear directions. For standard-normal $X$, the functions $X$, $X^2-1$ and $X^3-3X$ are mutually uncorrelated, although all are generated by the same $X$.

These are substantive reasons to examine an information-based representation. PCA can correctly describe a large linear span without that span representing equally many independent causes.

Distance-based neighborhood methods address a different object: the geometry of observations. Under distance-concentration conditions, nearest and farthest distances lose relative contrast, weakening neighborhood-based inference. Entropic Scree instead compares measurements through their information relationships across observations, avoiding the need to place nominal category codes and continuous magnitudes on one sample-distance scale.

## An information map of measurements

**Each point in the map is a measurement, not an observation.** The framework asks which measurements retain overlapping information and how their relationships are organized across the dataset.

For discrete or discretized measurements, the similarity is mutual information normalized by joint entropy:

$$\mathcal M_{ij}=\frac{I(X_i;X_j)}{H(X_i,X_j)}
=\frac{I(X_i;X_j)}{H(X_i)+H(X_j)-I(X_i;X_j)}.$$

This measures shared information relative to joint information. Zero denotes independence of the represented discrete variables; one denotes informational equivalence, with each determining the other almost surely. Its complement is normalized variation of information, a metric with informationally equivalent representations identified.

The normalization retains differences in informational resolution. For example, if a four-state measurement is reduced to two equally sized groups, the coarser measurement is completely determined by the finer one but retains only half its information; their normalized MI is $1/2$.

Population MI is invariant to invertible marginal transformations. Continuous empirical measurements enter this implementation through discretization, so binning resolution, ties and finite-sample estimation remain relevant. The representation is also pairwise: purely joint dependencies with zero pairwise MI, such as a fair-bit XOR construction, are not identified by this matrix.

### Centering and Euclidean coordinates

The map is centered at the measurements' unweighted centroid:

$$\mathcal M_c=\mathbf H\mathcal M\mathbf H,
\qquad \mathbf H=\mathbf I-\frac{\mathbf1\mathbf1^T}{m}.$$

For $D_{ij}=1-\mathcal M_{ij}$, this is also the classical-scaling construction

$$\mathcal M_c=-\tfrac12\mathbf H(2\mathbf D)\mathbf H.$$

If the centered matrix is positive semidefinite, its coordinates reproduce distances $\sqrt{2(1-\mathcal M_{ij})}$. Otherwise, retaining its positive spectrum provides a Euclidean approximation. The negative spectral mass diagnostic describes the extent of discarded negative eigenvalue mass.

Centering establishes a relative reference origin; it does not remove all finite-sample MI bias. The nonlinear similarity matrix is not generally subject to the covariance-specific $N-1$ ceiling, although exact centering limits its rank to $m-1$ and accurate estimation still requires sufficient observations.

### Bipolar modules: interpreting the poles

For positive eigenvalues, variable coordinates are $L_{ik}=\sqrt{\lambda_k^+}v_{ik}$.

- **Large same-sign loadings** contribute positively to a pair's represented inner product on that axis and identify a candidate descriptive anchor group.
- **Large opposite-sign loadings** contribute negatively, expressing relative informational estrangement along that contrast. Loading magnitude determines how strongly each measurement participates.
- **Near-zero loadings** indicate weak participation in that particular contrast, not necessarily low relevance elsewhere in the map.

Opposite poles do not mean negative correlation. Other axes can reinforce or offset a pair's relationship: vectors $(1,2)$ and $(-1,2)$ occupy opposite poles on the first axis but have a positive total inner product of $3$. Full-space orthogonality requires the summed inner product to be zero.

Pole groups describe observed informational configurations; they are not automatically isolated generating roots. Groups and individual measurements can recur across axes because several contrasts involve them.

With `extract_bipolar_modules = TRUE`, version 1.0.1 exports **eigenvector coefficients** $v_{ik}$. By default it selects up to $\max\{1,\lfloor0.20m\rfloor\}$ variables separately from each sign, ordered by absolute magnitude. `bipolar_top_n` adjusts that selection. These coefficients preserve within-axis signs and ordering but differ in scale from the geometric coordinates $L_{ik}$.

## Reading the scree and selecting boundaries

The logarithmic spectrum makes proportional changes visible across a wide range of eigenvalues. Two boundaries serve different purposes:

- **Observed Generative Rank ($K_{roots}$):** the selected endpoint of the primary spectral region, proposed as a diagnostic of resolvable generating structure.
- **Extended Signal Tail ($K_{extended}$):** the endpoint of the broader region included in shared-volume accounting. Extra directions in this region are not each counted as another root.

Both indices identify the **last retained axis**. A boundary of 20 retains axes 1 through 20.

The automated Dual-Diagnostic Ensemble first proposes a macroscopic search region. Its Log-Gap rule recommends the primary boundary using relative eigenvalue drops; its Triple-Tap rule recommends the extended boundary using departure from a local tail trajectory. Automatically,

$$K_{roots}=K_{log},\qquad K_{extended}=\max(K_{log},K_{tap}).$$

These are heuristic recommendations to inspect against the plotted spectrum. A smooth decline need not resolve a distinct boundary. The nominal threshold used by Triple-Tap is not an established whole-procedure error guarantee.

The relationship between $K_{roots}$ and generating count is evaluated empirically. They need not coincide in every construction: $r$ independent information-equivalent blocks, for example, produce a centered map of rank $r-1$.

## Informational volume and gravity

The published implementation uses working spectral weights $\eta_i=\max(10^{-9},\lambda_i)$. Define

$$p_i=\frac{\eta_i}{\sum_j\eta_j},\qquad
R_{eff}=\exp\left(-\sum_i p_i\ln p_i\right),\qquad
P_{sig}=\frac{\sum_{i=1}^{K_{extended}}\eta_i}{\sum_j\eta_j}.$$

**Total Unique Probabilistic Volume ($R_{eff}$)** is an entropy effective rank: the effective number of equally weighted spectral directions. It describes the distribution of spectral weight, rather than the number of generating mechanisms.

The framework expresses its allocations in **variable equivalents**, relative to the retained measurement count:

| Quantity | Definition | Interpretation |
|---|---|---|
| Redundant Signal Volume | $m-R_{eff}$ | Allocation associated with concentration of spectral weight. |
| Unique Signal Volume | $P_{sig}R_{eff}$ | Effective spectral volume allocated through the extended boundary. |
| Total Shared Signal Volume ($TSV$) | $(m-R_{eff})+P_{sig}R_{eff}$ | Sum of the redundant and unique shared allocations. |
| Idiosyncratic Informational Volume | $(1-P_{sig})R_{eff}$ | Remaining allocation outside that shared-volume accounting. |

**Total shared plus idiosyncratic volume equals $m$.** Unique and redundant shared volumes are components of the total shared volume, not additional amounts to add on top of it.

For example, $m=100$, $R_{eff}=80$ and $P_{sig}=0.25$ give 20 redundant, 20 unique shared and 60 idiosyncratic variable equivalents. Total shared volume is 40.

These are spectral allocations, not direct measurements of random-error fractions or joint Shannon entropy. Idiosyncratic volume can include measurement error, estimation effects and unresolved structure. Centering also contributes to the accounting: exactly independent population measurements have $R_{eff}=m-1$ under zero clipping, producing one redundant unit by convention.

### AIG, FSIG and the Structural Topology Profile

**Average Informational Gravity** allocates total shared volume per primary axis:

$$AIG=\frac{TSV}{K_{roots}}.$$

**Factor-Specific Informational Gravity**, the paper's formal name for FSIG, distributes that volume across individual primary informational axes:

$$FSIG_i=\frac{\eta_i}{\sum_{j=1}^{K_{roots}}\eta_j}TSV.$$

Consequently, $\sum_i FSIG_i=TSV$ and the mean FSIG is AIG. These summarize axis allocations; an axis need not correspond to one generating root.

The **Structural Topology Profile** is

$$T_i=\frac{FSIG_i}{FSIG_1}=\frac{\eta_i}{\eta_1}.$$

It shows whether the primary map is dominated by one direction, several leading directions, or a more even distribution. Loadings and domain knowledge help interpret that shape: eigenvalues alone do not determine cluster membership, network centralization or independent mechanisms.

## Comparing with PCA: the Dimensional Inflation Index

$$\Delta_K=K_{PCA}-K_{roots}.$$

This compares the dimension selected from a specified PCA representation with the primary Entropic Scree boundary. Report the preparation and selection rule used for each.

- **Similar selected dimensions** are compatible with a compact linear description but do not establish linear dependence or identical captured structure.
- **Positive divergence** can reveal linear representational expansion relative to the informational boundary, motivating a more compact nonlinear representation. Coding, estimation and selection rules also affect the comparison.
- **Negative differences** are possible and warrant examining the representations and boundary choices.

The main experiment reports a PCA Kaiser count of 5,762 and $K_{roots}=20$, giving $\Delta_K=5,742$. The Kaiser count is an explicitly identified alternative selection rule, not a claimed visual elbow.

## From diagnostic mapping to downstream models

Entropic Scree maps **measurements**. It does not supply an observation-level reconstruction rule. Multiplication of a numeric data matrix by its eigenvectors is algebraically possible, but those scores depend on coding and units and do not automatically recover latent record coordinates.

The diagnostic can guide a downstream representation: primary rank suggests its size, pole loadings describe observed groupings, and gravity summarizes spectral allocation. A nonlinear autoencoder can then learn from the original records. Reconstruction quality, stability and interpretation of that representation require their own assessment.

## Simulation code

The paper's experiments use the self-contained **1.0.0 beta** function embedded in the historical simulation script. Ordinary package use should start with CRAN; the following GitHub scripts identify the reviewed experimental implementation:

- [Main simulation](https://github.com/tjleestjohn/Entropic-Scree/blob/064dcbd1401633c16669fe3d9826c6b37a9fe2d8/Entropic.Scree.R.Simulation%20-%20ENLI.R)
- [Binning-ablation companion](https://github.com/tjleestjohn/Entropic-Scree/blob/064dcbd1401633c16669fe3d9826c6b37a9fe2d8/Appendix.B.Binning.Ablation.R)

### What the paper reports

The main construction uses 10,000 observations, 20,000 generated measurements and 20 Gaussian generating coordinates. Measurements combine nonlinear terms, continuous and binary manifestations, proxy-score perturbations and measurement error. A pool of 21,759 candidate generating terms supplies the nonlinear expansion.

The manuscript reports $K_{roots}=20$ and $K_{extended}=211$ in the main run. The primary boundary also remained 20 at the three tested bin targets—22, 100 and 208—while effective volume and gravity changed substantially. The 208-bin target was derived from Freedman–Diaconis recommendations but supplied to the same equal-frequency discretizer; it was not a switch to ordinary equal-width FD binning.

Additional correlated-root runs retained a transition near the generating count. These are findings for the configurations studied, not a general invariance theorem. The evaluated comparators were standardized PCA, rank-PCA and one specified RBF Kernel PCA construction.

### Historical implementation versus CRAN 1.0.1

The MI estimator, normalization, centering, eigenvalue floor, spectral entropy and gravity formulas are shared. Automatic recommendations differ:

- The historical macro-gap reference uses a 5% offset and unfiltered reference gaps; CRAN 1.0.1 uses 10%, filters large reference gaps and restricts candidate positions.
- The historical tail scan defaults to a fixed 20-value window and scans to index 1; the published version adapts the window, exposes `fwer_alpha` and ends at the primary recommendation.
- Historical optional pole output defaults to ten coefficients per sign; the package defaults to 20% of the retained measurement count per sign.

These differences can change boundaries and therefore gravity summaries. Use the archived implementation when reproducing the reported experiments, and label any comparison with the CRAN package separately.

### Download and inspect the scripts

```r
revision <- "064dcbd1401633c16669fe3d9826c6b37a9fe2d8"
base_url <- paste0(
  "https://raw.githubusercontent.com/tjleestjohn/Entropic-Scree/",
  revision, "/"
)

file_name <- "Entropic.Scree.R.Simulation - ENLI.R"
download.file(
  paste0(base_url, "Entropic.Scree.R.Simulation%20-%20ENLI.R"),
  destfile = file_name,
  mode = "wb"
)
file.edit(file_name)

download.file(
  paste0(base_url, "Appendix.B.Binning.Ablation.R"),
  destfile = "Appendix.B.Binning.Ablation.R",
  mode = "wb"
)
```

Inspect the scripts before running them in a dedicated R session: the main script clears the workspace, installs dependencies and compiles its embedded C++ code. Source compilation requires an appropriate compiler toolchain, such as Rtools on Windows. Also ensure `MASS` and `stringr` are installed; correlated-root configurations require `clusterGeneration`.

Run the ablation **after the main simulation in the same session**. It reuses the generated data, baseline result and embedded function. The main simulation disables interactive rank revision for unattended execution.

**Resource planning:** This is a large dense-matrix experiment. A single $20{,}000\times20{,}000$ double-precision matrix occupies approximately 3.2 GB before temporary copies, eigenvectors, generated data and comparator models. Dense eigendecomposition generally scales as $O(m^3)$ time with $O(m^2)$ matrix storage. Test a smaller configuration before attempting the full experiment; runtime and peak memory depend on hardware, numerical libraries and settings.

For reproducibility, retain the script revision, package or embedded-function version, seed, preparation settings, retained measurement count, both selected boundaries and numerical outputs with each run.

## Citation and related resources

```bibtex
@misc{leestjohn2026entropic-scree,
  title = {The Entropic Scree: An Information-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems},
  author = {Lee-St. John, Terrence J.},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.22028087},
  url = {https://doi.org/10.5281/zenodo.22028087},
  year = {2026}
}
```

For the software citation, run `citation("Entropic.Scree")` in R. The package identifier is [10.32614/CRAN.package.Entropic.Scree](https://doi.org/10.32614/CRAN.package.Entropic.Scree).

| Resource | Link |
|---|---|
| Entropic Scree preprint | [Zenodo](https://doi.org/10.5281/zenodo.22028086)|
| Published R package | [CRAN](https://CRAN.R-project.org/package=Entropic.Scree) |
| From Garbage to Gold (G2G) preprint | [arXiv](https://arxiv.org/abs/2603.12288) |
| G2G simulation repository | [GitHub](https://github.com/tjleestjohn/from-garbage-to-gold) |
| Author contact | [terry@enli.com.au](mailto:terry@enli.com.au) |
| Enli | [www.enli.com.au](https://www.enli.com.au) |

Distributed under the [Apache License 2.0](https://github.com/tjleestjohn/Entropic-Scree/blob/main/LICENSE).
