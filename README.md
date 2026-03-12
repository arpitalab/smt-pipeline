# smt-pipeline

MATLAB toolbox for single-molecule tracking (SMT) analysis built around
`TrajectoryWrapper` and `TrajectoryCollection`. It loads per-cell CSV trajectory
files, culls tracks by nucleus mask and quality criteria, removes global cell
motion, and runs MSD, Bayesian diffusivity (pEM), van Hove, and lifetime
analyses.

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

```matlab
run('setup_path.m')

files = dir('data/*.csv');
tc = TrajectoryCollection();
tc.ReadParams('experiment.toml');   % sets PixelSize, FrameInterval, etc.

for k = 1:numel(files)
    tc.addFromFile(fullfile(files(k).folder, files(k).name), ...
                   'FileID',    num2str(k), ...
                   'Condition', 'Control');
end

tc.summary();
tc.getBayesianDiffusivity('Condition', 'Control');
```

See [examples/example_workflow.m](examples/example_workflow.m) for a
fully-commented end-to-end script.

---

## Repo layout

```
smt-pipeline/
├── TrajectoryWrapper.m       % Single-file trajectory loader + filter
├── TrajectoryCollection.m    % Multi-file aggregator + analysis runner
├── TrajectoryAdapter.m       % SMD ↔ CSV bridge
├── TrackUtils.m              % Static track manipulation helpers
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
    └── example_workflow.m
```

---

## License

Source code in this repository (excluding `third_party/`) is released under the
MIT License. See individual files in `third_party/` for their respective
licenses.
