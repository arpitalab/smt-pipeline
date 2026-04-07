# smt-pipeline

MATLAB toolbox for single-molecule tracking (SMT) analysis built around
`TrajectoryWrapper` and `TrajectoryCollection`. Trajectories can be loaded from
two sources:

- **CSV files** produced by [quot](https://github.com/alecheckert/quot), a
  Python-based single-molecule localizer and tracker
- **SMD objects** produced by the
  [SMD tracking pipeline](https://github.com/arpitalab/sr_tracking), a
  MATLAB-based localization and tracking pipeline for TIFF and ND2 movies

Both paths produce identical `TrajectoryWrapper` objects that flow into the same
downstream analyses: track culling by nucleus mask and quality criteria, global
cell-motion removal, ensemble MSD, Richardson-Lucy MSD decomposition, Bayesian
diffusivity classification (pEM), van Hove correlation, and particle lifetime
analysis.

---

## Requirements

- MATLAB R2021b or later
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox *(optional — used for bootstrap CI; falls back to serial)*

See [DEPENDENCIES.md](DEPENDENCIES.md) for third-party package requirements.

---

## Installation

```bash
git clone --recurse-submodules https://github.com/<your-org>/smt-pipeline.git
```

Then, at the start of each MATLAB session:

```matlab
run('path/to/smt-pipeline/setup_path.m')
```

> **pEMv2** must be installed separately — see [DEPENDENCIES.md](DEPENDENCIES.md).

---

## Quick start

**From quot CSV files:**

```matlab
run('setup_path.m')

tc = TrajectoryCollection();
tc.ReadParams('experiment.toml');   % sets PixelSize, FrameInterval, etc.

files = dir('data/*.csv');
for k = 1:numel(files)
    tc.addFromFile(fullfile(files(k).folder, files(k).name), ...
                   'FileID', num2str(k), 'Condition', 'Control');
end

tc.summary();
tc.getBayesianDiffusivity('Condition', 'Control');
```

**From SMD objects** (after running `smd.localize()` → `smd.track()` → `smd.cull_tracks()`):

```matlab
run('setup_path.m')

tc = TrajectoryCollection();
for k = 1:numel(smd_list)
    tc.addFromSMD(smd_list{k}, 'FileID', num2str(k), 'Condition', 'Control');
end

tc.summary();
tc.getBayesianDiffusivity('Condition', 'Control');
```

See [examples/example_workflow.m](examples/example_workflow.m) and
[examples/example_smd_to_rl_msd.m](examples/example_smd_to_rl_msd.m) for
fully-commented end-to-end scripts.

---

## SMTDatabase

`SMTDatabase` is a SQLite-backed registry that tracks experiments, condition groups, and analysis results across sessions. It requires **MATLAB R2022b+** (uses the built-in `sqlite` function — no Database Toolbox needed).

### Concepts

| Level | What it represents |
|---|---|
| **Experiment** | One FOV / file (`TrajectoryWrapper`). Acquisition date is auto-extracted from the file path. |
| **Collection** | All FOVs sharing the same molecule / treatment / substrate / imaging rate. Spans multiple dates. |
| **Analysis run** | One completed MSD, pEM, or Lifetime analysis. Key numbers stored as JSON; full results live in the `.mat`. |

### Typical session

```matlab
db = SMTDatabase('smt_pipeline.db');   % creates the file if it doesn't exist

% --- first time: register and analyse ---
tc = TrajectoryCollection();
tc.ReadParams('experiment.toml');
tc.addFromFile('cell1.csv', 'Condition', 'Control');

db.registerCollection(tc, 'Name', 'H2B_Control_1kPa_slow', 'Treatment', 'Control');

tc.getRLDecomposition('LagTime', 4, 'String', 'H2B');
tc.getMSD('ExposureTime', 0.01);             % ensemble MSD + bootstrap fits
tc.getBayesianDiffusivity('Condition', 'Control');
tc.save('results/H2B_Control_1kPa_slow.mat');

db.updateCollection('H2B_Control_1kPa_slow', 'MatPath', 'results/H2B_Control_1kPa_slow.mat');
db.logAnalysis('H2B_Control_1kPa_slow', tc, 'MSD');
db.logAnalysis('H2B_Control_1kPa_slow', tc, 'pEM');

% --- later sessions: browse and reload ---
db.listCollections('Molecule', 'H2B')
db.listAnalyses('AnalysisType', 'pEM')

info = db.getCollection('H2B_Control_1kPa_slow');   % metadata + decoded summaries
tc2  = db.loadCollection('H2B_Control_1kPa_slow');  % loads the saved .mat

% Raw SQL for anything else
T = db.query('SELECT * FROM analysis_runs WHERE optimal_k = 3');
```

Molecule, Substrate, ImagingHz, and ExposureTime are **auto-parsed from the file path** using the HILO directory convention (`HILO/{molecule}/{eXms-iYms}/{substrate}/{YYYYMMDD}/...`). Preview what will be inferred before registering:

```matlab
SMTDatabase.parsePathMetadata('/path/to/MCF7_H2B_1kPa_e10ms_i200ms_trajs.csv')
```

See [examples/example_smtdatabase.m](examples/example_smtdatabase.m) for a fully-commented walkthrough.

---

## Repo layout

```
smt-pipeline/
├── TrajectoryWrapper.m       % Single-file trajectory loader + filter
├── TrajectoryCollection.m    % Multi-file aggregator + analysis runner
├── TrajectoryAdapter.m       % SMD ↔ CSV bridge
├── TrackUtils.m              % Static track manipulation helpers
├── SMTDatabase.m             % SQLite experiment registry
├── setup_path.m              % Run once to configure MATLAB path
├── utils/
│   ├── maxSpanningForest.m
│   ├── relativeTracksFromForest.m
│   ├── find_bound_particles_2.m
│   ├── fun_v.m
│   └── RL_analysis/          % MSD, van Hove, pEM helpers
├── third_party/
│   └── matlab-toml/          % git submodule (MIT license)
└── examples/
    ├── example_workflow.m
    ├── example_smd_to_rl_msd.m
    └── example_smtdatabase.m
```

---

## License

Source code in this repository (excluding `third_party/`) is released under the
MIT License. See individual files in `third_party/` for their respective
licenses.
