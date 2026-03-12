# Dependencies

## MATLAB Toolboxes

| Toolbox | Required by | Notes |
|---------|-------------|-------|
| Image Processing Toolbox | `find_bound_particles_2.m` (`bwboundaries`) | Required |
| Statistics & Machine Learning Toolbox | `classify_tracks.m` (`datasample`, `ecdf`) | Required |
| Optimization Toolbox | `RL_HILO.m` (`lsqnonlin`) | Required for HILO analysis |
| Parallel Computing Toolbox | `TrajectoryCollection` bootstrap CI | Optional — `parfor` falls back to `for` |

---

## matlab-toml

Bundled as a git submodule in `third_party/matlab-toml` (MIT license).
No additional action required after cloning with `--recurse-submodules`.

- Author: George Kaplan
- Repository: https://github.com/g-s-k/matlab-toml

---

## pEMv2

**Not bundled** (no redistribution license found).

pEMv2 is required for Bayesian diffusivity classification
(`tc.getBayesianDiffusivity()`). Install it manually:

1. Download from https://github.com/mcculloughlab/pEMv2 (or the lab's
   distribution page).
2. Unzip so that the following exist on disk:
   ```
   pEMv2-master/pEMv2/          % contains pEMv2.m, pEMv2_core.m, etc.
   pEMv2-master/pEMv2_bootstrap_CI.m
   ```
3. Add both directories to the MATLAB path **before** calling any
   `getBayesianDiffusivity` method:
   ```matlab
   addpath('path/to/pEMv2-master/pEMv2');
   addpath('path/to/pEMv2-master');
   ```
   You may add these lines to `setup_path.m` once you know the local install path.
