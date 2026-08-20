# The Entropic Scree:<br>An Information-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems

###### Initial Methods & Function Release: August 16, 2026 (Happy Birthday, Dad)

*[Terrence J. Lee-St. John, PhD](mailto:terry@enli.com.au)*

*[Enli: Predictive systems that remain stable under change](https://www.enli.com.au)*

**Links**

[![Read arXiv Preprint (Coming Soon)](https://img.shields.io/badge/arXiv_Preprint-Coming_Soon-lightgrey?style=for-the-badge)](#)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg?style=for-the-badge)](https://opensource.org/licenses/Apache-2.0)

<p align="center">
  <a href="#-usage-r-script"><strong>R Simulation Script 💻</strong></a><br>
  <a href="#-installation"><strong>Python and R Packages 📦 (Coming Soon)</strong></a>
</p>

> **TL;DR**
>
> If you are working with high-dimensional, mixed-type, noisy tabular data, standard PCA fundamentally fractures non-linear dependencies into "Spurious Orthogonal Dimensions," drastically overestimating the true rank of the system. Meanwhile, non-linear alternatives like Kernel PCA and Euclidean nearest-neighbor estimators suffer structural collapse when generative roots are entangled or sparse.
> 
> Offered as an upgrade over these baseline methods that rely on strict assumptions about the underlying generative process, the Entropic Scree framework uses Normalized Mutual Information to collapse those spurious expansions back towards their true generative roots. By more faithfully identifying the true generative rank, the results can be used to explicitly size neural bottlenecks for downstream non-parametric manifold extractors (like autoencoders).
> 
> The Entropic Scree also:
> * Quantifies the underlying "informational gravity" of the roots, offering insight into the system's overall average stability, as well as the relative informational dominance of individual roots;
> * Estimates the data's overall ratio of shared signal to unshared idiosyncratic informational variance (noise);
> * Serves as a powerful exploratory map that separates unrelated clusters of variables, allowing you to easily identify decoupled sub-networks.
>
> Ultimately, it provides an non-parametric, model-agnostic diagnostic whose comparative advantage aggressively compounds at scale and with system complexity.

## Quick Start - Entropic Scree Function v1.0.0 (R)

While the official CRAN and PyPI packages are under active development, you can use the Entropic Scree Function v1.0.0 in R immediately by sourcing the standalone function file.

Just copy and paste the following block into your R console (and press Enter) to automatically download and load the Entropic.Scree() function into your R environment:

```R
# 1. Define the direct URL to the raw function script on GitHub
url <- "https://raw.githubusercontent.com/tjleestjohn/entropic-scree/main/Entropic.Scree.v1.0.0%20-%20ENLI.R"

# 2. Define what you want to name the file on your computer
file_name <- "Entropic.Scree.v1.0.0 - ENLI.R"

# 3. Download the script to your current working directory
download.file(url, destfile = file_name)

# 4. Source the core function into your R environment
source(file_name)

# 5. Ex. To estimate the intrinsic rank and informational gravity of your tabular dataset (df)
# results <- Entropic.Scree(df)
```
---

## Structural Constraints of Linear Estimators

For over a century, the universal standard for evaluating a dataset's representational rank has been PCA and its variance-based scree plot. However, when deployed in modern, complex data environments, standard linear matrices systematically degrade across four dimensions:

1. **Mixed-Data Penalty:** Linear correlation deflates when continuous waves are evaluated against discrete categorical step-functions.
2. **Non-Linear Blindness:** Pure linear estimators ignore synergistic, thresholded, or polynomial dependencies.
3. **The Algebraic Rank Constraint:** If you have more variables than observations ($m > N$), PCA hits a hard algebraic wall, permanently capping extractable rank at $N-1$.
4. **Spurious Orthogonalization (Dimensional Inflation):** Because linear matrices cannot map non-linear states, they fragment continuous generative drivers into hundreds of spurious, independent linear dimensions. Because this fragmentation multiplies combinatorially as predictor breadth ($m$) grows, classical matrices catastrophically fracture at scale.

**The Result:** PCA tells you your data is driven by hundreds of weak linear components, when it is actually driven by a smaller set of highly non-linear, robust macro-structures.

## The Solution: Information-Theoretic Geometry

The Entropic Scree methodology resolves this by shifting the math from linear Euclidean space into topological information space.

To guarantee global geometric coherence and enforce a strict metric space, the framework constructs a pairwise Normalized Mutual Information (NMI) matrix utilizing Information-Theoretic Jaccard Similarity:

$$ \mathcal{M}_{i,j} = \frac{I(X_i; X_j)}{H(X_i) + H(X_j) - I(X_i; X_j)} $$

### Double-Centering Bias Correction (cMDS)
Because eigendecomposition cannot operate on raw similarities, and empirical mutual information estimators suffer from a strictly positive finite-sample bias, the framework executes a
**double-centering transformation** ($\mathcal{M}_c = \mathbf{H} \mathcal{M} \mathbf{H}$) prior to decomposition. This single operation serves a dual mathematical purpose:
1. It safely converts the distance manifold into a coordinate-ready inner-product (Gram) space, natively embedding the square root of twice the Normalized Variation of Information ($\sqrt{2 \cdot NVI}$) to ensure Positive Semi-Definiteness.
2. It algebraically mitigates positive estimation bias, perfectly centering the macroscopic noise bulk at zero and leaving the matrix to map pure **Topological Information Variance**.

By utilizing a highly optimized C++ backend to evaluate this matrix, the Entropic Scree:
* Evaluates pure shared dependency via Copula Theory (Sklar's Theorem), completely immune to marginal shape mismatches.
* Subsumes non-linear and discrete relationships back into their root generative source.
* Easily computes an $m \times m$ pairwise matrix regardless of sample size, utterly breaking the $N-1$ algebraic ceiling enforced by standard PCA.
* Survives Root Entanglement: While standard non-linear baselines (like Kernel PCA) suffer structural collapse under even mild generative root entanglement, this geometry maintains a rigid, near-invariant boundary at the true generative rank.
* Maps Structural Estrangement: Because it evaluates shared information rather than linear direction, the extracted axes physically push unrelated clusters of variables to opposite geometric poles. This allows practitioners to cleanly identify and untangle decoupled sub-networks directly from the primary factor loadings.

### The Diagnostic Framework and Automated Scanners
Exactly like Cattell's classical variance-based scree test, the Entropic Scree is fundamentally designed as a **visual diagnostic framework**. Visual inspection of the log-linear spectral decay remains the gold standard for identifying the structural elbow that separates the generative signal from the idiosyncratic noise baseline.

However, to provide an optional baseline convenience utility for rapid exploratory analysis, the script employs a Dual-Diagnostic Ensemble to automate extraction. First, it estimates a strict boundary for the macroscopic cliff (the top of the Idiosyncratic Informational Variance bulk) to ensure the search never wanders into the unstructured continuous floor. Then, operating exclusively within this bounded signal space, two complementary engines map distinct boundaries of the non-linear manifold:

* **Engine A (Log-Gap) isolates the Observed Generative Rank ($K_{roots}$):** Identifies the primary structural elbow by maximizing the logarithmic percentage drop between successive eigenvalues, successfully separating the core generative drivers from their own combinatorial expansions.

* **Engine B (Triple-Tap) maps the Extended Signal Tail ($K_{extended}$):** Applies a "Topological Stitch" to mathematically close the macro gap, then scans backward using a dynamically scaled 20-point linear regression to identify the exact index where the residual combinatorial signal significantly breaks out of the expected Idiosyncratic Informational Variance trajectory.

⚠️ Automated Elbow Detection Heuristic Warning

The Entropic Scree does not claim to possess the exact analytical Random Matrix Theory (RMT) bounds (such as the Marchenko-Pastur law) that strictly govern linear sample covariance matrices. Consequently, automated extraction in this space inherently relies on empirical heuristics.

However, from a pragmatic engineering perspective: **applying an approximate empirical heuristic to a structurally valid, non-linear metric space represents a strict, objective upgrade over applying any "standard" threshold (like the Kaiser criterion) to a fundamentally distorted linear space.**

Because real-world systems frequently exhibit complex internal hierarchies among correlated drivers — which can produce large internal informational variance drops independent of the idiosyncratic baseline — the automated scanner is provided strictly as an analytical baseline, not a universal algorithmic law. Practitioners should always visually inspect the generated entropic scree plot to formally confirm (or manually override) the detected elbow and verify the true macroscopic generative boundary.

---

## Actionable Engineering Metrics (AIG & FSIG)

The Entropic Scree doesn't just count dimensions; it calculates their exact probabilistic weight, translating abstract eigenvalues into physical **Variable Equivalents**.

* **Total Unique Probabilistic Volume:** The dataset's total continuous probability volume, containing both the unique signal volume and the system's Idiosyncratic Informational Variance (Structural Uncertainty, Independent Measurement Error, and Unshared Signal Geometry).
* **Unique Signal Volume:** The specific proportion of the Total Unique Probabilistic Volume strictly controlled by the signal axes.
* **Redundant Signal Volume:** The overlapping topological redundancy ($m - R_{eff}$) representing the signature of repeating signal axes.
* **Total Shared Signal Volume:** The combined volume of the signal axes (the sum of the Unique and Redundant Signal Volumes).
* **Idiosyncratic Informational Variance:** The remaining probability volume consisting of Structural Uncertainty, Independent Measurement Error, and Unshared Signal Geometry.
* **AIG (Average Informational Gravity):** How much physical data (in variable equivalents) the average extracted signal factor accounts for.
* **FSIG (Factor-Specific Informational Gravity):** The specific structural weight of individual signal axes, allowing you to assess the ability to disentangle dominant signals from weak, secondary signals.
* **Structural Topology Profile:** The normalized mass distribution of the complete extracted signal against the primary topological axis ($\text{FSIG}_{1-K} / \text{FSIG}_1$). This acts as a direct diagnostic of the system's macroscopic network topology, determining whether the variables form a highly centralized, entangled web or a decentralized, modular environment.

---

## Testing Linear Sufficiency ($\Delta_K$)

The Entropic Scree can also be utilized as a formal diagnostic bounding box for PCA itself. By comparing the rank extracted by classical PCA ($K_{PCA}$) against the structural rank mapped by the Entropic Scree ($K_{roots}$), practitioners can calculate the
**Dimensional Inflation Index ($\Delta_K$)**:

$$ \Delta_K = K_{PCA} - K_{roots} $$

* **Convergence ($\Delta_K \approx 0$):** Linear sufficiency confirmed. The data is well-approximated by a simple linear factor model, meaning spurious orthogonalization is negligible and classical PCA is likely sufficient.
* **Divergence ($\Delta_K \gg 0$):** Severe dimensional inflation detected. Standard linear estimators are fragmenting non-linear synergies or mixed-data shapes, strongly motivating the use of non-linear manifold learning architectures.

---

## ⚠️ Important Note on Factor Extraction / Dimensionality Reduction

As detailed in the formal paper, the Entropic Scree is a **diagnostic oracle, not a linear projection matrix**.

Do not attempt to project your raw data onto the extracted eigenvectors via a standard linear dot product ($X \cdot V$). The eigenvectors map shared probability mass (topological geometry), not continuous physical magnitude. You should utilize the Entropic Scree to accurately identify your true generative rank ($K_{roots}$), and then, assuming $\Delta_K$ is not negligible, pass that rank parameter into a non-linear manifold learner (e.g., Autoencoders, UMAP) to execute the actual physical data reduction.

---

## <a id="-installation"></a>📦 Python and R Package Installation (Coming Soon)

*Native packages for Python (via PyPI) and R (via CRAN) are currently in active development and will be released shortly.*

**For Python (Coming Soon):**
```bash
pip install Entropic-Scree
```

**For R (Coming Soon):**
```R
install.packages("Entropic.Scree")
```

---

## <a id="-usage-r-script"></a>💻 R Simulation (Utilizes Entropic Scree Function v1.0.0)

This repository includes a fully-annotated simulation in R that is available to run now. The script generates a hostile, high-dimensional synthetic environment ($m=20,000$, $N=10,000$, highly centralized and entangled network topology, $\sim 98.5\%$ Idiosyncratic Informational Variance, non-linear distortion), demonstrates the systematic degradation of standard PCA and non-linear baselines (which suffer total structural collapse under even mild generative root entanglement), and utilizes the Entropic Scree to extract the true generative rank ($r=20$).

> ⏳ Hardware & Runtime Warning: This simulation is extremely heavy on both RAM and CPU. Generating the synthetic environment (expanding 20 generative roots into 20,000 proxies via a 21,759-term non-linear design matrix) is just as memory- and time-intensive as computing the 200 million pairwise dependencies for the final NMI matrix. A minimum of 16GB of RAM (32GB+ recommended) is strongly advised to prevent out-of-memory crashes. Depending on your hardware, generating the data and executing all four baseline models may take anywhere from 1 to 4+ hours. Step away - the script will automatically generate the final comparison plots when it finishes.

**Notes:**
* **Automatic Setup:** The script is self-contained. It will automatically detect and install missing dependencies (e.g., `Rcpp`, `data.table`, `ggplot2`) upon the first run.
* **C++ Backend:** The pairwise mutual information engine is written in C++ via `Rcpp` and utilizes `OpenMP` for rapid multi-threading in RAM.
* **Automated End-to-End Execution:** The interactive prompt has been disabled for this simulation so you can seamlessly "select all and run" the entire file from beginning to end. The engine relies on the dual Log-Gap and Triple-Tap convergence to extract the rank and complete the baseline comparisons without requiring manual console input.

### Quick Start
Just copy and paste the following block into your R console (and press Enter) to automatically download and open the simulation script directly:

```R
# 1. Define the direct URL to the raw script on GitHub
url <- "https://raw.githubusercontent.com/tjleestjohn/entropic-scree/main/Entropic.Scree.R.Simulation%20-%20ENLI.R"

# 2. Define what you want to name the file on your computer
file_name <- "Entropic.Scree.R.Simulation - ENLI.R"

# 3. Download just the script
download.file(url, destfile = file_name)

# 4. Open the script in your editor (like RStudio)
file.edit(file_name)
```

### Advanced Validation: Discretization Ablation (Appendix B)
The repository also includes `Appendix.B.Binning.Ablation.R`. This script reproduces the stress-test from the paper's appendix, proving that while extreme binning heuristics (like Freedman-Diaconis) artificially compress probabilistic volume ($R_{eff}$), the extraction of the underlying generative rank ($K_{roots} = 20$) remains mathematically invariant. It also automatically generates the 1x3 panel graph showing the generative signal remaining cleanly separated despite the artificial compression of the topological space (driven by finite-sample estimation error).

**Instructions:** Run this script in the exact same R workspace **immediately following** the successful execution of the main simulation script. It relies on the synthetic `observed_data` matrix already generated in your computer's memory by the primary simulation, ensuring you don't have to wait for the hostile environment to be generated twice.

---

## Citation & Contact

The full methodology is formally introduced in an upcoming preprint. Once published, the arXiv link and full BibTeX citation will be updated here.

```bibtex
% Coming Soon
@article{leestjohn2026entropic-scree,
  title={The Entropic Scree: An Information-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems},
  author={Lee-St. John, Terrence J.},
  journal={arXiv preprint (*Coming Soon*)},
  year={2026}
}
```

| **Related Resources** | **Link** |
| --- | --- |
| **The Entropic Scree: An Information-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems (arXiv)** | *Coming Soon* |
| **From Garbage to Gold (G2G): A Data Architectural Theory of Predictive Robustness (arXiv)** | [G2G arXiv](https://arxiv.org/abs/2603.12288) |
| **G2G arXiv Simulation Repository** | [G2G GitHub](https://github.com/tjleestjohn/from-garbage-to-gold) |
| **Contact First Author** | [Email Me](mailto:terry@enli.com.au) |
| **Enli Official Website** | [Enli: Predictive systems that remain stable under change](https://www.enli.com.au) |
