# Decomposition-Driven DNA methylation Biomarker Identification

This repository provides the implementation of a decomposition-driven framework for identifying cancer-specific and tissue-of-origin (TOO) DNA methylation biomarkers from bulk tissue data using Non-negative Matrix Factorization (NMF).
A detailed description of the methodology and its applications can be found in our manuscript:
"Decomposition-Driven Biomarker Identification Enhances non-Invasive Early Cancer Detection".

A Python 3 environment is required, with commonly used packages including numpy, scipy, scikit-learn, and pandas.

The framework consists of several steps, each implemented as an independent script. The workflow and required data formats for each step are described below.

## 1. NMF Decomposition of Cancer and Normal Tissue Samples
Code: 01.CancerMarker.NMF.py
In this step, NMF is applied to bulk DNA methylation data from cancer and matched normal (e.g., para-carcinoma) tissue samples to decompose composite methylation signals into multiple latent components. These components are subsequently used to identify cancer-associated components and enhance cancer-specific marker selection.

### 1.1 Input Data and Format
Two data files are used in the step: (1) Bulk DNA methylation matrix; (2) Sample annotation table.

#### (1) Raw DNA methylation matrix
* Row: CpGs
* Column: Samples
* Value: DNA methylation levels (range: 0-1) 

| ID | Sample:1 | Sample:2 | Sample:3 | ... | Sample:N |
| --- | --- | --- | --- | --- | --- |
| CpG:1 | 0.0 | 0.23 | 0.86 | ... | 1.0 |
| CpG:2 | 0.08 | 0.07 | 0.34 | ... | 0.98 |
| CpG:3 | 0.10 | 0.32 | 0.19 | ... | 0.65 |
| ... | ... | ... | ... | ... | ... |
| CpG:M | 0.55 | 0.0 | 0.62 | ... | 0.04 |

#### (2) Sample annotation table 
* Must include a sample ID column and a pathological status column (e.g., Cancer or Normal).
* Used for CpG pre-filtering prior to NMF.

| Sample_ID | Pathological |
| --- | --- |
| Sample:1 | Cancer |
| Sample:2 | Cancer |
| Sample:3 | Cancer |
| ... | ... |
| Sample:N-1 | Normal |
| Sample:N | Normal |

### 1.2 Workflow
* CpG pre-selection based on cancer–normal differences.
    * CpGs with negligible differences between cancer and normal samples are filtered out.
    * ROC-AUC is calculated for each CpG; CpGs with AUC > 0.65 or < 0.35 are retained.
* Coverage filtering
    * CpGs with ≥10× sequencing depth in ≥70% of samples are retained.
* Missing value imputation
    * KNN-based imputation is applied.
* Internal control
    * An artificial CpG with 100% methylation is added as an internal control.
* NMF decomposition
    * The bulk methylation matrix is factorized into component-specific methylation profiles and sample-wise component fractions.
    * The number of components (K) is typically set to ~10, reflecting the cellular and sub-tissue heterogeneity of cancer tissues (e.g., malignant cells, immune infiltration, cancer-associated fibroblasts, stromal cells, endothelial cells, and residual normal cells).
    * Alternative values of K (e.g., 3, 5, 20, 25, 30) are encouraged for sensitivity exploration.
* Post-processing
    * The artificial CpG is used to rescale the H matrix to approximately 0–1. Correspondingly, component fractions in the W matrix are normalized such that fractions per sample sum to ~1.
    * Values exceeding 1 in the H matrix are truncated.

### 1.3 Output data and the format

#### (1) H matrix: DNA methylation profiles for components
* Rows: CpGs
* Columns: Components
* Values: Methylation levels (0-1) 

| ID | Component:1 | ... | Component:K |
| --- | --- | --- | --- |
| CpG:1 | 0.0 | ... | 1.0 |
| CpG:2 | 0.25 | ... | 0.84 |
| CpG:3 | 0.16 | ... | 0.65 |
| ... | ... | ... | ... |
| CpG:M | 0.77 | ... | 0.02 |

#### (2) W matrix: Fraction of components for each sample
* Rows: Components
* Columns: Samples
* Values: Relative contribution of each component per sample

| ID | Sample:1 | Sample:2 | Sample:3 | ... | Sample:N |
| --- | --- | --- | --- | --- | --- |
| Component:1 | 0.21 | 0.43 | 0.16 | ... | 0.17 |
| ... | ... | ... | ... | ... | ... |
| Component:K | 0.34 | 0.06 | 0.32 | ... | 0.25 |



## 2. Organize and Visualize the W matrix
Code: 02.WmatrixAnalysis.py
The W matrix is reformatted and visualized as a heatmap to illustrate the distribution of components across samples.

* Input:
    * W matrix
    * Sample annotation table
* Samples can be sorted by phenotype (e.g., cancer vs. normal).
* This step enables intuitive identification of cancer- and normal-associated components. 



## 3. Comparison of Top CpG Markers Between NMF-Based and Conventional Selection Methods 
Code: 03.cfDNAsignalDistribution.py
This step evaluates the effectiveness of CpG marker selection by comparing cfDNA methylation changes between cancer patients and healthy individuals for CpGs selected using the NMF framework versus conventional bulk-level comparisons.

### 3.1 Workflow
* Preparation
    * cfDNA methylation data for each CpG were obtained from the CFEA database. Cancer–healthy cfDNA methylation changes for each CpG were quantified using multiple metrics, including:
        * Mean distance between cancer and healthy samples.
        * Log fold change between the mean methylation levels of cancer and healthy samples.
        * Median distance between cancer and healthy samples.
        * Log fold change between the median methylation levels of cancer and healthy samples.
        * Difference in average methylation fraction between cancer and healthy samples.
        * Log fold change of the average methylation fraction between cancer and healthy samples.
    * NMF-based CpG marker selection. The methylation profile of the cancer-specific component was compared with that of the healthy blood cell background by subtracting the corresponding CpG methylation signals. CpGs were then ranked by their differential methylation relative to the blood background for NMF-based marker selection.
    * Conventional CpG marker selection. Bulk tissue methylation signals were compared between cancer and normal samples for each CpG using multiple metrics, including center distance, log fold change, ROC AUC, and p-values from the Wilcoxon rank-sum (U) test.
* Evaluation
    * For each selection strategy, the top 1,000 hyper-methylated CpGs and the top 1,000 hypo-methylated CpGs were selected based on:
        * NMF-derived cancer-specific component vs. healthy blood cells.
        * Bulk cancer vs. normal tissue, CpGs ranked by center distance.
        * Bulk cancer vs. normal tissue, CpGs ranked by log fold change.
        * Bulk cancer vs. normal tissue, CpGs ranked by ROC AUC.
        * Bulk cancer vs. normal tissue, CpGs ranked by p-value from the U-test.
    * For hyper- and hypo-methylated CpGs separately, CpGs selected by the NMF framework were compared with those selected by each conventional approach based on their cancer–healthy cfDNA methylation changes.

### 3.2 Input data and format

#### (1) CpG table for cancer-specific component vs. healthy blood cells
This table is used to select top CpGs with hyper- or hypo-methylation signals based on the cancer-specific component in comparison with the healthy blood cell background. Each row corresponds to a single CpG. The table should include the following columns:
* CpG genomic coordinate
* Methylation level in the cancer-specific component
* Average methylation level in healthy blood cells
* Methylation difference between the cancer-specific component and healthy blood cells
* Cancer–healthy cfDNA methylation change quantified by mean distance
* Cancer–healthy cfDNA methylation change quantified by log fold change of the mean
* Cancer–healthy cfDNA methylation change quantified by median distance
* Cancer–healthy cfDNA methylation change quantified by log fold change of the median
* Cancer–healthy cfDNA methylation change quantified by the difference in average methylation fraction across samples
* Cancer–healthy cfDNA methylation change quantified by log fold change of the average methylation fraction across samples

#### (2) CpG table for bulk cancer vs. normal tissue
This table is used to select top CpGs with hyper- or hypo-methylation signals in cancer tissues compared with normal tissue samples. Each row corresponds to a single CpG. The table should include the following columns:
* CpG genomic coordinate
* Mean methylation level in cancer samples
* Mean methylation level in normal samples
* Median methylation level in cancer samples
* Median methylation level in normal samples
* Center distance between cancer and normal samples
* Log fold change between cancer and normal samples
* ROC AUC distinguishing cancer and normal samples
* P-value from the Wilcoxon rank-sum (U) test between cancer and normal samples
* Cancer–healthy cfDNA methylation change quantified by mean distance
* Cancer–healthy cfDNA methylation change quantified by log fold change of the mean
* Cancer–healthy cfDNA methylation change quantified by median distance
* Cancer–healthy cfDNA methylation change quantified by log fold change of the median
* Cancer–healthy cfDNA methylation change quantified by the difference in average methylation fraction across samples
* Cancer–healthy cfDNA methylation change quantified by log fold change of the average methylation fraction across samples



## 4. NMF Across Multiple Cancer Types
Code: 04.TOOMarker.NMF.py
In this step, NMF is performed jointly on bulk DNA methylation profiles from cancer tissues and their corresponding normal tissues (e.g., para-carcinoma samples) across multiple cancer types. The objective is to decompose bulk methylation signals into a set of latent components, which are subsequently used to identify cancer type–specific components and improve the effectiveness of cancer type–specific marker discovery.

It is strongly recommended to include healthy blood cell samples in the NMF decomposition, so that blood-derived background signals can be separated from tissue-derived methylation signals.

The cancer type–specific markers identified in this step will be used in the tissue-of-origin localization module of a multi-cancer early detection (MCED) assay, enabling accurate cancer type determination in positive cases.

### 4.1 Input Data and Format
Three input data files are required for this step:
(1) Raw bulk DNA methylation matrix
(2) Sample annotation table
(3) CpG-level statistical summary table.

#### (1) Raw bulk DNA methylation value matrix
* Rows: CpGs
* Columns: Samples
* Values: Bulk DNA methylation levels ranging from 0 to 1

| ID | Sample:1 | Sample:2 | Sample:3 | ... | Sample:N |
| --- | --- | --- | --- | --- | --- |
| CpG:1 | 0.07 | 0.54 | 0.66 | ... | 1.0 |
| CpG:2 | 0.21 | 0.17 | 0.31 | ... | 0.83 |
| CpG:3 | 0.10 | 0.0 | 0.18 | ... | 0.57 |
| ... | ... | ... | ... | ... | ... |
| CpG:M | 0.35 | 0.05 | 0.92 | ... | 1.0 |

#### (2) Sample annotation table 
This table must contain:
* A sample ID column
* A pathological status column (e.g., Cancer or Normal)
* A tissue type column
This table is used for CpG pre-filtering prior to NMF decomposition.

| Sample_ID | Pathological | Tissue_type |
| --- | --- | --- |
| Sample:1 | Cancer | Esophageal |
| Sample:2 | Cancer | Esophageal |
| ... | ... | ... |
| Sample:# | Normal | Esophageal |
| Sample:# | Normal | Esophageal |
| ... | ... | ... |
| Sample:# | Cancer | Gastric |
| Sample:# | Cancer | Gastric |
| ... | ... | ... |
| Sample:# | Normal | Gastric |
| Sample:# | Normal | Gastric |
| ... | ... | ... |
| Sample:# | Cancer | Colorectal |
| Sample:# | Cancer | Colorectal |
| ... | ... | ... |
| Sample:# | Normal | Colorectal |
| Sample:# | Normal | Colorectal |
| ... | ... | ... |
| Sample:# | Cancer | Liver |
| Sample:# | Cancer | Liver |
| ... | ... | ... |
| Sample:# | Normal | Liver |
| Sample:# | Normal | Liver |
| ... | ... | ... |
| Sample:# | Cancer | Pancreatic |
| Sample:# | Cancer | Pancreatic |
| ... | ... | ... |
| Sample:# | Normal | Pancreatic |
| Sample:# | Normal | Pancreatic |
| ... | ... | ... |
| Sample:# | Normal | Blood_cell |
| Sample:# | Normal | Blood_cell |
| ... | ... | ... |
| Sample:N-1 | Normal | Blood_cell |
| Sample:N | Normal | Blood_cell |

#### (3) CpG-level statistical summary table.
This table contains precomputed statistical metrics for each CpG site and is used in the CpG pre-selection step. Each row corresponds to a single CpG site. The table should include the following columns:
* ANOVA p-values
    * The p-value obtained from ANOVA performed across cancer samples from multiple cancer types.
    * The p-value obtained from ANOVA performed across normal (para-carcinoma) samples from multiple tissue types.
* Sequencing coverage statistics for cancer samples
    * For each cancer type, the percentage of samples in which the corresponding CpG site has a sequencing coverage of ≥10×.
* Sequencing coverage statistics for normal samples
    * For each normal tissue type, the percentage of samples in which the corresponding CpG site has a sequencing coverage of ≥10×.
These metrics collectively enable the identification of CpG sites that exhibit tissue- or cancer-type–associated methylation patterns while maintaining sufficient sequencing coverage for robust downstream analysis.

### 4.2 Workflow
* Preparation
    * For each CpG site, perform:
        * ANOVA across cancer samples from different cancer types
        * ANOVA across normal (para-carcinoma) samples from different tissue types
    * This evaluates whether the CpG harbors tissue-type– or cancer-type–associated methylation signals.
* CpG Pre-selection
    * CpGs are retained if they satisfy either of the following:
        * ANOVA p < 0.05 among cancer samples across cancer types
        * ANOVA p < 0.05 among normal samples across tissue types
    * Additionally, CpGs must have:
        * ≥10× sequencing coverage in ≥50% of cancer samples and ≥50% of normal samples
* Missing Value Imputation
    * Missing methylation values are imputed using the K-nearest neighbors (KNN) method.
* Internal Control CpG Addition
    * An artificial CpG with 100% methylation is added as an internal control for downstream normalization.
* NMF Decomposition
    * NMF is applied to the filtered bulk DNA methylation matrix.
    * Parameter K (number of components) should be selected based on the number of cancer types. If there are X cancer types, it is recommended to set K ≥ 2X + 3, accounting for:
        * X cancer type–specific components
        * X corresponding normal tissue–specific components
        * One shared cancer component
        * One shared normal tissue component
        * One immune-related component
    * Testing alternative K values is encouraged to gain intuitive understanding of the decomposition.
* Initial NMF Outputs
    * H matrix (methylation profile matrix):
        * Rows: CpGs
        * Columns: Components
    * W matrix (fraction matrix):
        * Rows: Components
        * Columns: Samples
* Normalization
    * The artificial CpG is used to rescale: Methylation values in the H matrix to approximately the range 0–1
    * Component fractions in the W matrix so that the sum per sample is approximately 1
* Final Adjustment
    * Values in the H matrix exceeding 1 are truncated to 1.
    * The resulting H and W matrices are treated as the final NMF outputs.

### 4.3 Output data and the format

#### (1) H Matrix: Component-specific DNA Methylation Profiles
* Rows: CpGs
* Columns: Components
* Values: Methylation levels (0–1)

| ID | Component:1 | ... | Component:K |
| --- | --- | --- | --- |
| CpG:1 | 1.0 | ... | 0.0 |
| CpG:2 | 0.55 | ... | 0.89 |
| CpG:3 | 0.24 | ... | 0.35 |
| ... | ... | ... | ... |
| CpG:M | 0.37 | ... | 0.08 |

#### (2) W matrix: Component Fractions per Sample
* Rows: Components
* Columns: Samples
* Values: Component fractions

| ID | Sample:1 | Sample:2 | Sample:3 | ... | Sample:N |
| --- | --- | --- | --- | --- | --- |
| Component:1 | 0.15 | 0.23 | 0.11 | ... | 0.16 |
| ... | ... | ... | ... | ... | ... |
| Component:K | 0.14 | 0.03 | 0.28 | ... | 0.05 |



## 5. Identification of Cancer Type–Specific Components
Code: 05.WmatrixStat.py
In this step, statistical tests are performed on the component fractions in the W matrix to determine whether individual components are significantly enriched in specific sample groups, such as a particular cancer type. This analysis enables the identification of cancer type–specific components, which will be used in subsequent downstream analyses.

### 5.1 Input Data and Required Information
* W matrix generated from multi-cancer NMF
    * Please refer to Section 4.3 for a detailed description of the W matrix format.
* Sample table
    * Please refer to Section 4.1 for the description of the sample table.
    * Group annotations in the sample table (e.g., cancer type and healthy blood cells) are used to define sample groups for statistical testing.

### 5.2 Workflow
* For each NMF component, pairwise t-tests are performed on the component fraction values between each cancer type and all other cancer types, as well as healthy blood cell samples.
* A component is considered a candidate cancer type–specific component if it is significantly enriched in one cancer type compared with all other cancer types and healthy blood cells.



## 6. Organizing the H Matrix for Methylation Profile Visualization
Code: 06.HmatrixAnalysis.py
In this step, the H matrix is reformatted and organized to facilitate downstream visualization of methylation profiles and the identification of cancer type–specific CpG markers.

* The H matrix generated from the multi-cancer NMF decomposition is required as the input for this step.



