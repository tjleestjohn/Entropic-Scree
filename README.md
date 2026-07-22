# The Entropic Scree:<br>A Universal Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems
### Initial Methods & Function Release: July 2026

*[Terrence J. Lee-St. John, PhD](mailto:terry@enli.com.au)*

*[Enli: Predictive systems that remain stable under change](https://www.enli.com.au)*

**Links**

[![Read arXiv Preprint (Coming Soon)](https://img.shields.io/badge/arXiv_Preprint-Coming_Soon-lightgrey?style=for-the-badge)](#)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg?style=for-the-badge)](https://opensource.org/licenses/Apache-2.0)

<p align="center">
  <a href="#-installation"><strong>Python and R Packages 📦 (Coming Soon)</strong></a> &nbsp;|&nbsp;
  <a href="#-usage-r-script"><strong>Jump to R Simulation Script 💻 (Available Now)</strong></a>
</p>

> **TL;DR**
>
> If you have high-dimensional, mixed-type, noisy tabular data, standard Principal Component Analysis (PCA) will lie to you about its true dimensionality.
>
> The **Entropic Scree** replaces variance with **Normalized Mutual Information** to reveal the true underlying generative rank of your data. It bypasses algebraic sample-size limits ($m > N$), handles non-linear interactions natively, and provides actionable engineering metrics (Informational Gravity) to map the true geometry of your data.

---

## Standard PCA Limitations

For over a century, the universal standard for evaluating a dataset's representational rank has been PCA and its variance-based scree plot. However, when deployed in modern, complex data environments, standard linear matrices suffer a **Structural Collapse** across four dimensions:

1. **Mixed-Data Penalty:** Linear correlation deflates when continuous waves are evaluated against discrete categorical step-functions.
2. **Non-Linear Blindness:** Pure linear estimators ignore synergistic, thresholded, or polynomial dependencies.
3. **Sample-Size Prison:** If you have more variables than observations ($m > N$), PCA hits a hard algebraic wall, permanently capping extractable rank at $N-1$.
4. **Orthogonal Splintering:** Because linear matrices cannot map non-linear states, they shatter continuous generative drivers into hundreds of fragmented, spurious linear dimensions (Dimensional Inflation).

**The Result:** PCA tells you your data is driven by 600 weak linear components, when it is actually driven by 10 highly non-linear, robust macro-structures.

## The Solution: Information-Theoretic Geometry

The Entropic Scree methodology resolves this by shifting the math from linear Euclidean space into topological information space.

To guarantee global geometric coherence and enforce a strict metric space, the framework constructs a pairwise Normalized Mutual Information (NMI) matrix utilizing Information-Theoretic Jaccard Similarity:

$$ \mathcal{M}_{i,j} = \frac{I(X_i; X_j)}{H(X_i) + H(X_j) - I(X_i; X_j)} $$

By utilizing a highly optimized, C++ backend to evaluate this matrix, the Entropic Scree:
* Evaluates pure shared dependency via Copula Theory (Sklar's Theorem), completely immune to marginal shape mismatches.
* Subsumes non-linear and discrete relationships back into their root generative source.
* Easily computes an $m \times m$ pairwise matrix regardless of sample size, utterly breaking the $N-1$ algebraic ceiling enforced by standard PCA.

### Automated Elbow Detection
To identify the boundary between true structural signal and finite-sample noise (the Marchenko-Pastur bulk), the script employs a log-space algorithmic assessment. It scans backward from the deep noise tail, estimating the true structural elbow only when it detects a massive, sustained phase transition, safely ignoring localized noise ripples.

**🔍 Heuristic Warning:** The current form of the automatic elbow detector is provided strictly as a convenience heuristic. Because real-world noise distributions can vary unpredictably, the user should visually inspect the generated entropic scree plot and formally confirm (or manually override) the detected elbow to ensure the correct generative rank is selected.

---

## ⚖️ Actionable Engineering Metrics (AIG & FSIG)

The Entropic Scree doesn't just count dimensions; it calculates their exact probabilistic weight, translating abstract eigenvalues into physical **Variable Equivalents**.

* **Total Unique Probabilistic Volume:** The dataset's total continuous probability volume, containing both unique signal volume and idiosyncratic noise (Structural Uncertainty and independent measurement error) volume.
* **Unique Signal Volume:** The specific proportion of the Total Unique Probabilistic Volume strictly controlled by the signal axes.
* **Redundant Signal Volume:** The overlapping topological redundancy ($m - R_{eff}$) representing the signature of repeating signal axes.
* **Total Signal Volume:** The combined volume of the signal axes (the sum of the Unique and Redundant Signal Volumes).
* **Idiosyncratic Noise Volume:** The remaining unshared probability volume consisting of Structural Uncertainty and independent Measurement Error.
* **AIG (Average Informational Gravity):** How much physical data (in column equivalents) the average extracted signal factor accounts for.
* **FSIG (Factor-Specific Informational Gravity):** The specific structural weight of individual signal axes, allowing you to assess the ability to disentangle dominant signals from weak, secondary signals.

---

## 📐 Testing Linear Sufficiency ($\Delta_K$)

The Entropic Scree can also be utilized as a formal diagnostic bounding box for PCA itself. By comparing the rank extracted by classical PCA ($K_{rlzd}$) against the structural rank mapped by the Entropic Scree ($K_{elbow}$), practitioners can calculate the **Dimensional Inflation Index ($\Delta_K$)**:

$$ \Delta_K = K_{rlzd} - K_{elbow} $$

* **Convergence ($\Delta_K \approx 0$):** Linear sufficiency confirmed. The data is well-approximated by a simple linear factor model, meaning orthogonal splintering is negligible and classical PCA is safe to use.
* **Divergence ($\Delta_K \gg 0$):** Severe dimensional inflation detected. Standard linear estimators are shattering non-linear synergies or mixed-data shapes, strongly motivating the use of non-linear manifold learning architectures.

---

## ⚠️ Important Note on Factor Extraction / Dimensionality Reduction

As detailed in the formal paper, the Entropic Scree is a **diagnostic oracle, not a linear projection matrix**.

Do not attempt to project your raw data onto the extracted eigenvectors via a standard linear dot product ($X \cdot V$). The eigenvectors map shared probability mass (topological geometry), not continuous physical magnitude. You should utilize the Entropic Scree to accurately identify your true generative rank ($K_{elbow}$), and then, assuming $\Delta_K$ is not negligible, pass that rank parameter into a non-linear manifold learner (e.g., Autoencoders, UMAP) to execute the actual physical data reduction.

---

## <a id="-installation"></a>📦 Python and R Installation (Coming Soon)

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

This repository includes a fully-annotated simulation in R that is **available to run now**. The script generates a hostile, high-dimensional synthetic environment ($m=10,000$, $N=5,000$, 97% pure noise, non-linear distortion), demonstrates the structural collapse of standard PCA, and utilizes the Entropic Scree to flawlessly extract the true generative rank ($r=10$).

**Notes:**
* **Automatic Setup:** The script is self-contained. It will automatically detect and install missing dependencies (e.g., `Rcpp`, `data.table`, `ggplot2`) upon the first run.
* **C++ Backend:** The pairwise mutual information engine is written in C++ via `Rcpp` and utilizes `OpenMP` for rapid multi-threading natively in RAM.
* **Interactive Mode (User Input Required):** The script pauses at the diagnostic elbow to provide a visual plot preview, requiring the user to confirm or manually override the rank in the R console. **While you cannot "select all and run" the entire file at once, you can safely select and run everything from the top down to (and including) the `calculate_entropic_scree()` function call in Section 4. RUN ENTROPIC SCREE.** Once you confirm the rank in the console, you can run the rest of the script.

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

---

## 📬 Citation & Contact

This methodology is formally introduced in an upcoming preprint. Once published, the arXiv link and full BibTeX citation will be updated here.

```bibtex
% Coming Soon
@article{leestjohn2026entropic-scree,
  title={The Entropic Scree: A Universal Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems},
  author={Lee-St. John, Terrence J.},
  journal={arXiv preprint (*Coming Soon*)},
  year={*Coming Soon*}
}
```

| **Related Resources** | **Link** |
| --- | --- |
| **The Entropic Scree: A Universal Diagnostic Framework for Intrinsic Rank and Informational Gravity in Tabular Systems (Preprint)** | *Coming Soon* |
| **From Garbage to Gold: A Data Architectural Theory of Predictive Robustness (Preprint)** | [arXiv cs.LG](https://arxiv.org/abs/2603.12288) |
| **G2G Preprint Simulation Repository** | [From Garbage to Gold GitHub](https://github.com/tjleestjohn/from-garbage-to-gold) |
| **Contact First Author** | [Email Me](mailto:terry@enli.com.au) |
| **Enli Official Website** | [Enli: Predictive systems that remain stable under change](https://www.enli.com.au) |
