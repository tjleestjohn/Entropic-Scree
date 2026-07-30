# The Entropic Scree:<br>An Informational-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems

*[Terrence J. Lee-St. John, PhD](mailto:terry@enli.com.au)*

*[Enli: Predictive systems that remain stable under change](https://www.enli.com.au)*

**Links**

###### Initial Methods & Function Release: July 2026

[![Read arXiv Preprint (Coming Soon)](https://img.shields.io/badge/arXiv_Preprint-Coming_Soon-lightgrey?style=for-the-badge)](#)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg?style=for-the-badge)](https://opensource.org/licenses/Apache-2.0)

<p align="center">
  <a href="#-installation"><strong>Python and R Packages 📦 (Coming Soon)</strong></a> &nbsp;|&nbsp;
  <a href="#-usage-r-script"><strong>Jump to R Simulation Script 💻 (Available Now)</strong></a>
</p>

> **TL;DR**
>
> Offered as an upgrade over standard PCA and other methods that rely on strict assumptions about the nature of the underlying generative process, the Entropic Scree allows practitioners to estimate true generative rank of their tabular data by evaluating a valid non-linear, information-theoretic space.
>
> If you have high-dimensional, mixed-type, noisy tabular data, standard PCA will systematically misrepresent its true dimensionality, while other topological estimators (e.g., TWO-NN, MLE) also fail in these regimes due to distance concentration.
>
> The **Entropic Scree** replaces variance with **Normalized Mutual Information** to bypass algebraic sample-size limits ($m > N$), natively handle non-linear interactions, and collapse spurious linear expansions back to their generative roots.
>
> While all statistical estimators are universally bound by the resolution limits of finite sample sizes and predictor breadth, this framework ensures that complex non-linearities are permanently collapsed to their root rather than fragmented into spurious dimensions.
>
> Ultimately, what it extracts is not an unobservable physical parts list, but an operational map of what your environment possesses the statistical and structural power to resolve.

---

## Structural Constraints of Linear Estimators

For over a century, the universal standard for evaluating a dataset's representational rank has been PCA and its variance-based scree plot. However, when deployed in modern, complex data environments, standard linear matrices systematically degrade across four dimensions:

1. **Mixed-Data Penalty:** Linear correlation deflates when continuous waves are evaluated against discrete categorical step-functions.
2. **Non-Linear Blindness:** Pure linear estimators ignore synergistic, thresholded, or polynomial dependencies.
3. **The Algebraic Rank Constraint:** If you have more variables than observations ($m > N$), PCA hits a hard algebraic wall, permanently capping extractable rank at $N-1$.
4. **Spurious Orthogonalization (Dimensional Inflation):** Because linear matrices cannot map non-linear states, they fragment continuous generative drivers into hundreds of spurious, independent linear dimensions.

**The Result:** PCA tells you your data is driven by 600 weak linear components, when it is actually driven by 10 highly non-linear, robust macro-structures.

## The Solution: Information-Theoretic Geometry

The Entropic Scree methodology resolves this by shifting the math from linear Euclidean space into topological information space.

To guarantee global geometric coherence and enforce a strict metric space, the framework constructs a pairwise Normalized Mutual Information (NMI) matrix utilizing Information-Theoretic Jaccard Similarity:

$$ \mathcal{M}_{i,j} = \frac{I(X_i; X_j)}{H(X_i) + H(X_j) - I(X_i; X_j)} $$

### Double-Centering Bias Correction (cMDS)
Because eigendecomposition cannot operate on raw similarities, and empirical mutual information estimators suffer from a strictly positive finite-sample bias, the framework executes a **double-centering transformation** ($\mathcal{M}_c = H \mathcal{M} H$) prior to decomposition. This single operation serves a dual mathematical purpose:
1. It safely converts the distance manifold into a coordinate-ready inner-product (Gram) space, natively embedding the square root of twice the Normalized Variation of Information ($\sqrt{2 \cdot NVI}$) to ensure Positive Semi-Definiteness. 
2. It algebraically mitigates positive estimation bias, perfectly centering the macroscopic noise bulk at zero and leaving the matrix to map pure **Topological Information Variance**.

By utilizing a highly optimized C++ backend to evaluate this matrix, the Entropic Scree:
* Evaluates pure shared dependency via Copula Theory (Sklar's Theorem), completely immune to marginal shape mismatches.
* Subsumes non-linear and discrete relationships back into their root generative source.
* Easily computes an $m \times m$ pairwise matrix regardless of sample size, utterly breaking the $N-1$ algebraic ceiling enforced by standard PCA.

### The Diagnostic Framework and Automated Scanner
Exactly like Cattell's classical variance-based scree test, the Entropic Scree is fundamentally designed as a **visual diagnostic framework**. Visual inspection of the log-linear spectral decay remains the gold standard for identifying the structural elbow that separates the generative signal from the idiosyncratic noise baseline.

However, to provide an optional baseline convenience utility for rapid exploratory analysis, the script employs a Dual-Diagnostic Ensemble to automate extraction. First, it estimates a strict boundary for the macroscopic noise cliff (the top of the unstructured noise floor) to ensure the search never wanders into the idiosyncratic baseline. Then, operating exclusively within this bounded signal space, two independent engines evaluate the phase transition:

* **Engine A (Log-Gap):** Identifies the structural elbow by maximizing the logarithmic percentage drop between successive eigenvalues, resolving the scale imbalance between massive and subtle generative factors.

* **Engine B (Triple-Tap):** Applies a "Topological Stitch" to mathematically close the macro gap, then scans backward using a dynamically scaled 10-point quadratic regression to identify the exact index where the signal violently breaks out of the expected noise trajectory.

**⚠️ Heuristic Warning & The "Upgrade" Proposition ⚠️** 
The Entropic Scree does not claim to possess the exact analytical Random Matrix Theory (RMT) bounds (such as the Marchenko-Pastur law) that strictly govern linear sample covariance matrices. Consequently, automated extraction in this space inherently relies on empirical heuristics. 

However, from a pragmatic engineering perspective: **applying an approximate empirical heuristic to a structurally valid, non-linear metric space represents a strict, objective upgrade over applying any "standard" threshold (like the Kaiser criterion) to a fundamentally distorted linear space.** 

Because real-world systems frequently exhibit complex internal hierarchies among correlated drivers — which can produce large internal informational variance drops independent of the noise floor — the automated scanner is provided strictly as an analytical baseline, not a universal algorithmic law. Practitioners should always visually inspect the generated entropic scree plot to formally confirm (or manually override) the detected elbow and verify the true macroscopic generative boundary.

---

## Actionable Engineering Metrics (AIG & FSIG)

The Entropic Scree doesn't just count dimensions; it calculates their exact probabilistic weight, translating abstract eigenvalues into physical **Variable Equivalents**.

* **Total Unique Probabilistic Volume:** The dataset's total continuous probability volume, containing both unique signal volume and idiosyncratic noise (Structural Uncertainty and independent measurement error) volume.
* **Unique Signal Volume:** The specific proportion of the Total Unique Probabilistic Volume strictly controlled by the signal axes.
* **Redundant Signal Volume:** The overlapping topological redundancy ($m - R_{eff}$) representing the signature of repeating signal axes.
* **Total Signal Volume:** The combined volume of the signal axes (the sum of the Unique and Redundant Signal Volumes).
* **Idiosyncratic Noise Volume:** The remaining unshared probability volume consisting of Structural Uncertainty and independent Measurement Error.
* **AIG (Average Informational Gravity):** How much physical data (in column equivalents) the average extracted signal factor accounts for.
* **FSIG (Factor-Specific Informational Gravity):** The specific structural weight of individual signal axes, allowing you to assess the ability to disentangle dominant signals from weak, secondary signals.
* **Structural Topology Ratio:** The ratio of the secondary factors to the primary topological axis ($FSIG_{2-K} / FSIG_1$). This acts as a direct diagnostic of the system's macroscopic network topology, determining whether the variables form a highly centralized, entangled web or a decentralized, modular environment.

---

## Testing Linear Sufficiency ($\Delta_K$)

The Entropic Scree can also be utilized as a formal diagnostic bounding box for PCA itself. By comparing the rank extracted by classical PCA ($K_{rlzd}$) against the structural rank mapped by the Entropic Scree ($K_{elbow}$), practitioners can calculate the **Dimensional Inflation Index ($\Delta_K$)**:

$$ \Delta_K = K_{rlzd} - K_{elbow} $$

* **Convergence ($\Delta_K \approx 0$):** Linear sufficiency confirmed. The data is well-approximated by a simple linear factor model, meaning spurious orthogonalization is negligible and classical PCA is likely sufficient.
* **Divergence ($\Delta_K \gg 0$):** Severe dimensional inflation detected. Standard linear estimators are fragmenting non-linear synergies or mixed-data shapes, strongly motivating the use of non-linear manifold learning architectures.

---

## ⚠️ Important Note on Factor Extraction / Dimensionality Reduction ⚠️

As detailed in the formal paper, the Entropic Scree is a **diagnostic oracle, not a linear projection matrix**.

Do not attempt to project your raw data onto the extracted eigenvectors via a standard linear dot product ($X \cdot V$). The eigenvectors map shared probability mass (topological geometry), not continuous physical magnitude. You should utilize the Entropic Scree to accurately identify your true generative rank ($K_{elbow}$), and then, assuming $\Delta_K$ is not negligible, pass that rank parameter into a non-linear manifold learner (e.g., Autoencoders, UMAP) to execute the actual physical data reduction.

---

## <a id="-installation"></a>📦 Python and R Package Installation (Coming Soon)

*Native packages for Python (via PyPI) and R (via CRAN) are currently in active development and will be released shortly. In the meantime, please utilize the standalone R simulation script below.*

**For Python (Future Release):**
```bash
# pip install entropic-scree (Coming Soon)
```

**For R (Future Release):**
```R
# install.packages("entropicscree") (Coming Soon)
```

---

## <a id="-usage-r-script"></a>💻 Simulation (R Script)

This repository includes a fully-annotated simulation in R that is available to run now. The script generates a hostile, high-dimensional synthetic environment ($m=10,000$, $N=5,000$, highly centralized and entangled network topology, $\sim 99.5\%$ idiosyncratic noise, non-linear distortion), demonstrates the systematic degradation of standard PCA, and utilizes the Entropic Scree to extract the true generative rank ($r=10$).

**Notes:**
* **Automatic Setup:** The script is self-contained. It will automatically detect and install missing dependencies (e.g., `Rcpp`, `data.table`, `ggplot2`) upon the first run.
* **C++ Backend:** The pairwise mutual information engine is written in C++ via `Rcpp` and utilizes `OpenMP` for rapid multi-threading natively in RAM.
* **Automated End-to-End Execution:** The interactive prompt has been disabled for this simulation so you can seamlessly "select all and run" the entire file from beginning to end. The engine relies on the dual Log-Gap and Triple-Tap convergence to extract the rank and complete the baseline comparisons without requiring manual console input.

### Quick Start
Copy and paste the following code block into your R console or RStudio to download and open the script directly:

```R
# 1. Define the direct URL to the raw script on GitHub
url <- "https://raw.githubusercontent.com/tjleestjohn/entropic-scree/main/Entropic_Scree_R_Simulation%20-%20ENLI.R"

# 2. Define what you want to name the file on your computer
file_name <- "Entropic_Scree_R_Simulation - ENLI.R"

# 3. Download just the script
download.file(url, destfile = file_name)

# 4. Open the script in your editor (like RStudio)
file.edit(file_name)
```

### Advanced Validation: Discretization Ablation (Appendix B)
The repository also includes `Appendix_B_Binning_Ablation.R`. This script reproduces the stress-test from the paper's appendix, proving that while extreme binning heuristics (like Freedman-Diaconis) artificially compress probabilistic volume ($R_{eff}$), the extraction of the underlying generative rank ($K_{elbow} = 10$) remains mathematically invariant. It also automatically generates the 1x3 panel graph showing the generative signal remaining cleanly separated despite the artificial inflation of the continuous noise tail.

**Instructions:** Run this script in the exact same R workspace **immediately following** the successful execution of the main simulation script. It relies on the synthetic `observed_data` matrix already generated in your computer's memory by the primary simulation, ensuring you don't have to wait for the hostile environment to be generated twice.

---

## Citation & Contact

The full methodology is formally introduced in an upcoming preprint. Once published, the arXiv link and full BibTeX citation will be updated here.

```bibtex
% Coming Soon
@article{leestjohn2026entropic-scree,
  title={The Entropic Scree: An Informational-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems},
  author={Lee-St. John, Terrence J.},
  journal={arXiv preprint (*Coming Soon*)},
  year={2026}
}
```

| **Related Resources** | **Link** |
| --- | --- |
| **The Entropic Scree: An Informational-Theoretic Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems (Preprint)** | *Coming Soon* |
| **From Garbage to Gold: A Data Architectural Theory of Predictive Robustness (Preprint)** | [arXiv cs.LG](https://arxiv.org/abs/2603.12288) |
| **G2G Preprint Simulation Repository** | [From Garbage to Gold GitHub](https://github.com/tjleestjohn/from-garbage-to-gold) |
| **Contact First Author** | [Email Me](mailto:terry@enli.com.au) |
| **Enli Official Website** | [Enli: Predictive systems that remain stable under change](https://www.enli.com.au) |
