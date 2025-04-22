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
* **NIfTI MATLAB Toolbox:** For reading/writing NIfTI files. Path configuration required in scripts.
* **`bwdistsc` MATLAB function:** Custom anisotropic distance transform. Path configuration required.
* **RegionGrowing24 Executable:** External tool for local maxima detection. Path configuration required.
* **(Optional) SGE Cluster:** The provided `submit_*.sh` scripts are configured for Sun Grid Engine. Adaptation is needed for other schedulers (Slurm, LSF, PBS) or local execution.

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

### Output

* **Intermediate files:** Stored within each subject's directory.
* **Per-subject backup:** `<SubjectID>_sulci_metrics_counts_backup.txt` in each subject's directory.
* **Main Result:** A CSV file (default `summary_sulci_metrics_MEAN_v3.csv`) in `config.rootFolder` containing mean metrics per sulcus for all subjects.

---

## File Descriptions

* `StageA_GenerateVolumes2.m`: Stage 1 script - generates intermediate volumes.
* `StageB_CalculateMetrics.m`: Stage 2 script - calculates metrics.
* `submit_stage_A.sh`: Example SGE submission script for Stage A.
* `submit_stage_B.sh`: Example SGE submission script for Stage B.
* `subject_list.txt` (User-created): List of subject IDs.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. ## ✍️ Citation

If you use this pipeline in your research, please cite:

> Mateos MJ, Gastelum-Strozzi A, Barrios FA, Bribiesca E, Alcauter S, Marquez-Flores JA. A novel voxel-based method to estimate cortical sulci width and its application to compare patients with Alzheimer’s disease to controls. *NeuroImage*. 2020 Feb 15;207:116343. doi: 10.1016/j.neuroimage.2019.116343.
