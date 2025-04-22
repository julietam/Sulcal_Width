# 🧠 Sulcal Width Measurement Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ## 📄 Description

This repository contains scripts for a two-stage pipeline to measure cortical sulcal width from T1-weighted MRI images processed with FreeSurfer. The method calculates mean distance metrics within dilated sulcal label regions based on local maxima identified on a distance map.

This implementation is based on the method described in:
> Mateos MJ, Gastelum-Strozzi A, Barrios FA, Bribiesca E, Alcauter S, Marquez-Flores JA. A novel voxel-based method to estimate cortical sulci width and its application to compare patients with Alzheimer’s disease to controls. *NeuroImage*. 2020 Feb 15;207:116343. doi: 10.1016/j.neuroimage.2019.116343.

---

## 🗺️ Overview

The pipeline consists of two main stages:

1.  **Stage A: Generate Intermediate Volumes (`StageA_GenerateVolumes2.m`)** 
    * Takes FreeSurfer outputs as input (`mri/brain.nii`, `labels/*.label.nii`).
    * Generates intermediate files (distance map, brain closure, local maxima, combined/dilated sulci labels) required for metric calculation.

2.  **Stage B: Calculate Mean Sulci Metrics (`StageB_CalculateMetrics.m`)** 
    * Loads intermediate volumes from Stage A.
    * Calculates mean sulcal width metrics per label.
    * Appends results to a central summary CSV file.
    * Saves a per-subject backup file with means and counts.

---

## ✅ Requirements

* **MATLAB:** Including the Image Processing Toolbox.
* **FreeSurfer:** Requires `recon-all` output directories and the `mri_convert` tool available in the system PATH.
* **RegionGrowing24 Executable:** Provided in the `bin/` directory of this repository.
    * **Note:** The provided executable might be specific to a certain OS (e.g., Linux). If it doesn't run on your system, you may need to compile it from source (source not provided here).
* **`bwdistsc` MATLAB function:** This custom function for anisotropic distance transform must be available to MATLAB.
* **NIfTI MATLAB Toolbox:** Required for reading/writing NIfTI files. 
* **(Optional) SGE Cluster:** 🖥️ The provided `submit_*.sh` scripts are configured for Sun Grid Engine. Adaptation is needed for other schedulers (Slurm, LSF, PBS) or local execution.

---

## ⚙️ Setup

1.  **Clone Repository:**
    ```bash
    git clone <your-repository-url>
    cd <your-repository-directory>
    ```
2.  **Configure Scripts:** 🔧
    * **Crucially, edit paths** within the `config = struct();` sections of `StageA_GenerateVolumes2.m` and `StageB_CalculateMetrics.m` to match your environment (e.g., `config.rootFolder`, tool paths, output file). See script comments for details.
    * Adjust processing parameters (`config.imageSize`, `config.maxExpectedLabels`, etc.) if needed.
3.  **Prepare Subject List:** 
    * Create `subject_list.txt` listing one FreeSurfer subject ID per line. Place it where `SCRIPT_DIR` points in the submission scripts.
4.  **Configure Submission Scripts (if using a cluster):**
    * Adapt SGE options (`#$`) for your cluster/scheduler.
    * Set paths (`SCRIPT_DIR`, `LOG_DIR`, `MATLAB_CMD`).
    * **Update Job Dependency:** In `submit_stage_B.sh`, replace `<StageA_JobID>` with the actual job ID from Stage A submission (syntax may vary by scheduler).

---

## Usage

### Input Data Structure

1.  **Root Folder:** Specified as `config.rootFolder` in the scripts.
2.  **Subject Folders:** Inside the root folder, one directory per subject, named with their `<SubjectID>` (should be the FreeSurfer `recon-all` output directory).
3.  **Required Files:** Within each subject folder, the pipeline needs:
    * `mri/brain.nii`
    * `mri/ribbon.nii`
    * `labels/*.label.nii` (Sulcal label files from FreeSurfer).

### Running the Pipeline (Cluster Example using SGE)

1.  **Submit Stage A:**
    ```bash
    qsub submit_stage_A.sh
    # Note the returned Job ID (e.g., 12345)
    ```
2.  **Update Stage B Dependency:** Edit `submit_stage_B.sh`, setting the dependency flag (e.g., `-hold_jid 12345`).
3.  **Submit Stage B:**
    ```bash
    qsub submit_stage_B.sh
    ```

### Running Manually

Adapt the scripts by hardcoding the `datasetID` variable and removing `exit()` calls for interactive/single-subject execution. Run Stage A first, then Stage B.

To run the pipeline for a single subject directly in MATLAB without using the cluster submission scripts:

1.  **Adapt `StageA_GenerateVolumes/-git.m`**:
    * Find the line `datasetID = getenv('JOB_SUBJECT_ID');`.
    * Comment out this line by adding a `%` at the beginning: `% datasetID = getenv('JOB_SUBJECT_ID');`.
    * Add a new line below it to define the subject ID directly, replacing `'YourSubjectIDHere'` with the actual folder name of your subject: `datasetID = 'YourSubjectIDHere';`.
    * (Optional) For interactive debugging, you might also want to comment out the `exit(0);` and `exit(1);` lines near the end of the main function and its `catch` block.
2.  **Run `StageA_GenerateVolumes_git.m`** in MATLAB.
3.  **Adapt `StageB_CalculateMetrics_git.m`**:
    * Apply the same modifications as in Stage A: comment out the `getenv` line, add a line to hardcode the same `datasetID`, and optionally comment out the `exit()` calls.
4.  **Run `StageB_CalculateMetrics_git.m`** in MATLAB after Stage A has finished successfully for that subject.

### Output

* **Intermediate files:** Stored within each subject's directory.
* **Per-subject backup:** `<SubjectID>_sulci_metrics_counts_backup.txt` in each subject's directory.
* **Main Result:** A CSV file (default `summary_sulci_metrics_MEAN_v3.csv`) in `config.rootFolder` containing mean metrics per sulcus for all subjects.

---

## 📚 File Descriptions

* `StageA_GenerateVolumes2.m`: Stage 1 script - generates intermediate volumes.
* `StageB_CalculateMetrics.m`: Stage 2 script - calculates metrics.
* `submit_stage_A.sh`: Example SGE submission script for Stage A.
* `submit_stage_B.sh`: Example SGE submission script for Stage B.
* `bin/RegionGrowing24`: Pre-compiled executable for local maxima detection (used by Stage A). *(Adjust path if placed elsewhere)*
* `subject_list_TEMPLATE.txt`: Template file showing the format for listing subject IDs (one per line). Rename to `subject_list.txt` and populate with your subject IDs.
* `helper_functions/bwdistsc.m`: Custom anisotropic distance transform function.
* `LICENSE.md`: *(Add this line after creating the file)* Defines the permissions and limitations for using this code.
* `.gitignore`: *(Add this line after creating the file)* Specifies files intentionally untracked by Git (e.g., logs, results).
---
## 🌳 Repository Structure

```text
.
├── bin/
│   └── RegionGrowing24         # Executable dependency
├── helper_functions/
│   └── bwdistsc.m              # Custom distance function
├── StageA_GenerateVolumes2.m     # Stage 1 script
├── StageB_CalculateMetrics.m     # Stage 2 script
├── submit_stage_A.sh           # Stage 1 submission script template
├── submit_stage_B.sh           # Stage 2 submission script template
├── subject_list_TEMPLATE.txt   # Example subject list
├── LICENSE.md                  # License file (You need to create this)
├── README.md                   # This file
└── .gitignore                  # Git ignore file (You need to create this)

## 📜 License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. ## ✍️ Citation

If you use this pipeline in your research, please cite:

> Mateos MJ, Gastelum-Strozzi A, Barrios FA, Bribiesca E, Alcauter S, Marquez-Flores JA. A novel voxel-based method to estimate cortical sulci width and its application to compare patients with Alzheimer’s disease to controls. *NeuroImage*. 2020 Feb 15;207:116343. doi: 10.1016/j.neuroimage.2019.116343.
