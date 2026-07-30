# Milestone 0 baseline manifest

## Status

The code-side Milestone 0 gate passed on 2026-07-30. The long-duration
conservation and experiment-validation gate remains pending and must be run by
the case owner before Milestone 0 is declared fully closed.

## Reproducible environment

```text
Project baseline commit: e340ee0313410a7739203ae862435be3277f1839
Project development branch: development/condensation-model
OpenFOAM version: Foundation v14
OpenFOAM commit: 0db507470b
OpenFOAM build: 14-0db507470b5e
Compiler: g++ 8.4.0
WM_OPTIONS: linux64GccDPInt32Opt
```

The baseline case pair is:

```text
Standard: CCSJ01_multiphaseEuler_wallCondensation_CHT
Modular:  CCSJ01_cigmaMPE_wallCondensation_CHT
```

## Clean-build result

`./Allwmake` completed successfully. The rebuilt library checksums were:

```text
e924bb678a4ffb837a5240fda44602b2be5f05afd1bdf06ffb6264c3f19dffa6  libcigmaMPECondensationModels.so
fb4c3a2a1e674f186697f9847c87b6fb3188e6e8a0a05c60cebd41608b7b74c3  libcigmaMPEFluidSolver.so
4f990bd8298d14fdb1242bff1f4129db27fa561b8de01b459fe88a3cab478a2f  libcigmaMPESolidSolver.so
```

The condensation-library hash changed when the unchanged OpenFOAM v14
`massDiffusionLimitedPhaseChange` source was added for the CONAN validation
case. The fluid- and solid-solver hashes remained unchanged.

## One-step regression result

Both serial cases ran from `0` to `1e-6 s` with function objects disabled.
The test used temporary case copies and did not change either source case.

```text
Standard exit status: 0
Modular exit status: 0
Written files per case: 32
Raw file differences: 2
Normalised field differences: 0
Normalised aggregate hash: ef993c3745520c30d86bd60ab80fc60a3c68dd8acf3563cc582c41d672c46ed4
Residual and mDot trace hash: 1f4a9313138f804af1ae603acf123230183b338ab0e985c398416f1620618dd9
Standard warnings: 12
Modular warnings: 12
Standard peak RSS: 193048 kB
Modular peak RSS: 191152 kB
```

The two raw differences were `fluid/alphat.gas` and
`fluid/alphat.liquid`. They contained only the expected library-name metadata:

```text
libmultiphaseEulerFvModels.so
libcigmaMPECondensationModels.so
```

After normalising that metadata, all 32 output files were identical. The full
solver logs were also identical after normalising solver names, case paths,
process identifiers, and timing metadata. The standard and modular runs
reported the same pre-existing 12 wedge-planarity warnings.

The observed wall times were `0.42 s` for the standard case and `0.34 s` for
the modular case. These single-run startup values are recorded only as smoke
metrics and must not be used as performance evidence.

## Source audit

The local copies of the following model directories were byte-identical to the
OpenFOAM v14 source tree:

```text
phaseTurbulenceStabilisation
massDiffusionLimitedPhaseChange
wallPhaseChange
wallCondensation
```

Generated `Make/linux*` files are excluded from version control. Only the
portable `Make/files` and `Make/options` inputs belong to the source baseline.

## Pending validation gate

The one-step initial state has zero wall-condensation rate and is insufficient
for physics validation. The case owner must run the selected long validation
case and provide:

1. total gas and liquid mass histories;
2. integrated condensation rate for every active wall;
3. total fluid and solid energy histories;
4. pressure, water-species, and energy residual histories;
5. wall-clock time, peak RSS, MPI process count, and write settings;
6. the experimental comparison quantities selected for the paper.

The baseline must be archived before any new condensation or transport model
is enabled.
