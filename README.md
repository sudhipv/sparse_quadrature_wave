# Sparse Quadrature Wave

This repository contains two closely related workflows for stochastic acoustic wave propagation:

- **`acoustic_MCS/`** implements Monte Carlo simulations (serial and MPI) together with scripts for generating Karhunen–Loève/PCE data and meshes.
- **`NISP_wave/`** implements non-intrusive spectral projection (NISP) with sparse quadrature, sharing the same preprocessing utilities.
- **`matlab_scripts/`** provides MATLAB post-processing utilities for inspecting Monte Carlo outputs (`MCS_pdf`, `MCS_variance`) and for comparing solver results against reference data (`MCS_process`).

The repo now tracks only source code and reproducible input decks. Large solver outputs such as `u_nisp_mean.mat` and `u_MCS_mean.mat` have been removed, and `.gitignore` prevents them from being re-added inadvertently.

---

## Repository Layout

- `acoustic_MCS/`
  - `preprocess.sh` – builds KLE/PCE data (`klePceData/`) and meshes (`meshData/`) via `gfortran` and `gmsh`.
  - `acoustic_parallel_MCS.py`, `acoustic_parallel_3RV.py`, `acoustic_pdf.py` – Monte Carlo drivers (MPI-enabled).
  - `klePceData/`, `meshData/`, `ssfem_solver/` – generated input data, mesh converters, and reference CSV outputs from the SSFEM solver.
- `NISP_wave/`
  - `preprocess.sh`, `quad.sh` – same preprocessing pipeline plus UQTk quadrature generation.
  - `nisp_serial.py`, `nisp_parallel*.py`, `nisp_SparseQuad.py` – NISP solvers matching different execution modes.
  - `serial/`, `quadData/`, `klePceData/`, `meshData/` – solver outputs and reusable data.
- `matlab_scripts/`
  - `MCS_pdf.m` – aggregates PDF samples `u_pdf_*.mat` from `outputs/<case>/...` to a single `u_pdf.mat`.
  - `MCS_variance.m` – computes confidence bounds across experiments stored under `outputs/3rv_sigma1/3rv_pdf_<idx>/`.
  - `MCS_process.m` – loads Monte Carlo results plus SSFEM reference CSV files to compare mean/variance histories.
- `.gitignore` – ignores heavy generated data (`outputs/`, solver `.mat` files, etc.).

---

## Requirements

- Python 3 with `dolfin`/FEniCS, `numpy`, `scipy`, `mpi4py`.
- GNU Fortran (`gfortran`) for building KLE/PCE data generators.
- `gmsh` for mesh generation.
- [UQTk](https://github.com/sandialabs/UQTk) binaries for quadrature generation (referenced in `NISP_wave/quad.sh`).
- MATLAB (or Octave) for the scripts inside `matlab_scripts/`.
- MPI runtime (e.g., OpenMPI) for parallel simulations.

---

## Typical Workflow

1. **Preprocess shared data**
   ```bash
   cd acoustic_MCS
   ./preprocess.sh        # builds klePceData/ and meshData/
   cd ../NISP_wave
   ./preprocess.sh        # builds equivalent data for NISP
   ./quad.sh              # generates quadrature points via UQTk
   ```
   Adjust mesh density (`lc`) and UQTk parameters inside the scripts as needed.

2. **Run simulations**
   - Monte Carlo (MPI):
     ```bash
     mpirun -np <ranks> python3 acoustic_parallel_MCS.py
     ```
     This produces `u_MCS_mean.mat`, `u_MCS_sd.mat`, `nS.mat`, etc., which remain untracked because of `.gitignore`.
   - NISP:
     ```bash
     python3 nisp_serial.py
     # or mpirun -np <ranks> python3 nisp_parallel1.py
     ```
     Outputs such as `serial/d3l*/u_nisp_{mean,sd}.mat` are likewise ignored.

3. **Post-process in MATLAB**
   - `matlab_scripts/MCS_pdf.m` expects Monte Carlo pdf chunks under `outputs/3rv_sigma3/3rv_sigma3_pdf_20/`.
   - `matlab_scripts/MCS_variance.m` reads aggregated pdf files under `outputs/3rv_sigma1/3rv_pdf_<idx>/`.
   - `matlab_scripts/MCS_process.m` compares Monte Carlo results against reference CSVs under `acoustic_MCS/ssfem_solver/13k_3rv_3s/`.
   Each script now infers the repository root automatically:
   ```matlab
   script_dir = fileparts(mfilename('fullpath'));
   repo_root = fileparts(script_dir);
   ```
   Update the case-specific subdirectories at the top of each script to match your experiment layout.

---

## Regenerating Removed Files

All deleted `.mat` artifacts are solver outputs and can be regenerated:

- **NISP results** (`NISP_wave/serial/d3l*/u_nisp_mean.mat`, `u_nisp_sd.mat`): rerun `nisp_serial.py` or the desired parallel variant.
- **Monte Carlo summaries** (`acoustic_MCS/u_MCS_mean.mat`, `u_MCS_sd.mat`, `nS.mat`): rerun `acoustic_parallel_MCS.py` (or related scripts).
- **PDF chunk consolidation**: run `matlab_scripts/MCS_pdf.m` after Monte Carlo chunk generation.

If additional heavy outputs should remain local only, add their patterns to `.gitignore`.

---

## Notes and Tips

- Keep large experiment folders (e.g., `outputs/`) outside version control; the scripts now point to `repo_root/outputs/...` by default, so create those directories locally.
- `NISP_wave/quad.sh` assumes a local UQTk installation path – adjust `cd /Users/.../UQTk-install/bin` to your environment.
- To keep the repository lightweight before pushing, periodically remove generated `.mat` files via the provided clean-up scripts or `rm` commands; `.gitignore` prevents accidental staging.
- The MATLAB scripts assume the SSFEM reference CSVs in `acoustic_MCS/ssfem_solver/` stay in place; keep those files tracked to reproduce figures.

With the repo trimmed to sources and configuration files only, it is well below GitHub’s size thresholds and safe to push without Git LFS. Refer to this README whenever you regenerate data to ensure the structure stays organized.
