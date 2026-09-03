# cigmaMPEFoam

`cigmaMPEFoam` is a modular OpenFOAM Foundation v14 project for multi-region
conjugate heat-transfer simulations with multiphase wall condensation.

## Project layout

```text
cigmaMPEFoam/
├── applications/
│   ├── solvers/
│   │   ├── fluid/
│   │   │   ├── cigmaMPEFluid/
│   │   │   └── cigmaSinglePhaseFluid/
│   │   └── solid/
│   │       └── cigmaMPESolid/
│   └── utilities/
│       └── cigmaRepairSinglePhaseCheckpoint/
├── src/
│   ├── condensationModels/
│   ├── multiphaseMomentumTransportModels/
│   └── singlePhaseCondensation/
├── Allwmake
├── Allclean
└── README.md
```

## Solver modules

`cigmaMPEFluid` is a runtime-selectable subclass of the standard v14
`Foam::solvers::multiphaseEuler` module. It retains the standard equations and
solution stages and adds an optional phase-specific Courant-number limiter.

`cigmaMPESolid` is a runtime-selectable subclass of the standard v14
`Foam::solvers::solid` module. It initially overrides no equations or
solution-stage methods.

`cigmaSinglePhaseFluid` is a runtime-selectable subclass of the standard v14
`Foam::solvers::multicomponentFluid` module. It solves one gas momentum,
pressure and energy system with multicomponent species transport. When no
`constant/condensationProperties` dictionary is present, its solved fields are
identical to the native module in the one-step regression test.

The modules are installed as `libcigmaMPEFluidSolver.so`,
`libcigmaSinglePhaseFluidSolver.so` and `libcigmaMPESolidSolver.so`. Their names follow the
`lib<solverName>Solver.so` convention used automatically by `foamMultiRun`.

## Phase-specific Courant-number limits

`cigmaMPEFluid` supports optional per-phase Courant-number limits for
transient calculations. The standard global `maxCo` remains active and is
used as the fallback for cases without `phaseMaxCo`.

```text
adjustTimeStep  yes;
maxCo           5;
maxDeltaT       0.01;

phaseMaxCo
{
    gas         5;
    liquid      1;
}
```

The applied time-step is the minimum allowed by the global limit, all listed
moving-phase limits, `maxDeltaT`, function objects, and fvModels. Cases that
do not define `phaseMaxCo` retain the standard OpenFOAM v14 behaviour.

## Condensation models

`src/condensationModels` contains the minimum OpenFOAM v14 source set required
by the CIGMA and CONAN wall-condensation cases:

- `wallCondensation` and its phase-change-rate patch field;
- `wallPhaseChange` and the `alphatPhaseChangeWallFunction` patch field;
- `phaseTurbulenceStabilisation`, which is also selected by the case;
- `massDiffusionLimitedPhaseChange`, required by the OpenFOAM v14 CONAN
  validation tutorial.

These four model groups are unchanged copies of the OpenFOAM Foundation v14
sources. Their runtime type names, equations, correlations, coefficients, and
dictionary interfaces are therefore unchanged. The selected `ArdenBuck`
saturation-pressure model remains supplied by the standard OpenFOAM v14
saturation-model library and is not duplicated here.

The library also contains
`timeVaryingExternalWallLayersHeatTransferCoefficient`, an optional boundary
utility that places conductive wall layers in series with a time-varying
external heat-transfer coefficient supplied as a `Function1`. This utility
does not modify the baseline wall-condensation closure.

The library is installed as `libcigmaMPECondensationModels.so`. A case using
this library should not also load `libmultiphaseEulerFvModels.so` for these same
types, because duplicate runtime-selection registrations are unnecessary.

## Multiphase buoyancy-aware turbulence

`src/multiphaseMomentumTransportModels` provides the
`cigmaBuoyantKOmegaSST` RAS model for the `phaseCompressible` runtime-selection
table used by `multiphaseEuler`. It is installed as the additive library
`libcigmaMPEMomentumTransportModels.so` and does not replace the corresponding
single-phase cigma momentum-transport library.

Load it from `system/controlDict`:

```text
libs
(
    "libcigmaMPEMomentumTransportModels.so"
);
```

For an HPC installation with active cigma jobs, build only this new library:

```sh
wmake libso /path/to/cigmaMPEFoam/src/multiphaseMomentumTransportModels
```

Do not run `Allwmake` in that situation because it intentionally starts with
`Allclean` and rebuilds all project libraries.

## Single-phase diffusion-layer wall condensation

`src/singlePhaseCondensation` adapts the containmentFOAM single-phase
diffusion-layer formulation and its `saturatedSteam`,
`nonCondensableMassFraction`, and `condensingWallVelocity` boundary
conditions to OpenFOAM Foundation v14.
The implementation uses the native v14 effective species flux returned by the
selected multicomponent thermophysical-transport model. This includes the
molecular and turbulent contributions configured by that model.

Select the solver in `system/controlDict`:

```text
solver  cigmaSinglePhaseFluid;
```

Enable the model in `constant/condensationProperties`:

```text
solveWallCondensation   yes;
wallCondensationModel   diffusionLayer;
condensingSpecie        H2O;
maxCondWallMassFlux     100;

saturationPressure
{
    type    ArdenBuck;
}

latentHeatModel
{
    type    polynomial;
}
```

On each condensing wall, use:

```text
// 0/H2O
type    saturatedSteam;
value   uniform 0;

// 0/<solved non-condensable specie>, if more than one NCG is present
type    nonCondensableMassFraction;
value   uniform 0;

// 0/U
type    condensingWallVelocity;
value   uniform (0 0 0);

// 0/T, for an external-temperature wall balance
type                externalCondensationTemperature;
Ta                  constant 304;
h                   uniform 4000;
qcond               qcond;
qcondRelaxation     0.2;
relaxation          0.2;
value               uniform 304;
```

The model writes `massTransferRate` and `qcond`. Their wall-patch values are
the condensation mass flux `[kg/m2/s]` and latent heat flux `[W/m2]`
respectively. `externalCondensationTemperature` adds the positive latent heat
flux to the native OpenFOAM v14 external-temperature wall balance. The default
latent-heat polynomial coefficients are the containmentFOAM values; they may
be overridden with `A`, `B`, `C` and `D` in `latentHeatModel`.

For mixtures with two or more non-condensable species,
`nonCondensableMassFraction` balances the outward condensation-suction flux
with the native v14 effective species diffusivity. Apply it to each solved
non-condensable mass fraction on a condensing wall; the default specie remains
the algebraic mass-fraction closure. Its condensing-wall boundary condition
must be `calculated`, allowing `normaliseY()` to maintain
`Ydefault = 1 - sum(Ysolved)` at the wall and in restart checkpoints.

The legacy containmentFOAM `DtWallFunction`, `velocityScale` and top-level
`Sct` entries are not read by this v14 port. Equivalent diffusivity settings
belong to the selected v14 thermophysical-transport model. This distinction
must be accounted for when constructing a numerically matched validation case.

### Restarting single-phase condensation

The condensing specie is kept active so that `H2O` is written at every new
checkpoint. Before each write, solved-species wall conditions are finalised and
the algebraic default specie is rebuilt, so a new simulation started from time
zero can restart directly without a repair step. This operation does not alter
the solved internal fields of `H2O` or `HE`.

Checkpoints created by older builds may contain either no active `H2O` or a
default-specie `AIR` field that is inconsistent with the wall values of `H2O`
and `HE`. The repair utility below is only for these legacy checkpoints.

Build only the required components while other jobs are active:

```sh
wmake libso src/singlePhaseCondensation
wmake applications/utilities/cigmaRepairSinglePhaseCheckpoint
```

Repair a decomposed checkpoint before restarting it:

```sh
mpirun -np 16 cigmaRepairSinglePhaseCheckpoint \
    -case /path/to/case -parallel -region fluid -time 2000
```

For a serial checkpoint, omit `mpirun` and `-parallel`. If `H2O` exists, the
utility preserves it and only rebuilds `AIR`. If `H2O` is missing, the original
`0/fluid/H2O` boundary-condition template and `H2O_0` are used to recover it.
In both modes, the old default-specie field is preserved as
`AIR_beforeCheckpointRepair`, and a second repair of the same checkpoint is
refused.

## Build

`main` is the single supported integration branch. Update it before building:

```sh
git switch main
git pull --ff-only origin main
```

Source OpenFOAM Foundation v14 and run:

```sh
./Allwmake
```

The build script first runs `Allclean`, then builds the condensation library,
the multiphase momentum-transport library, the single-phase condensation
library, both fluid solver modules, the solid solver module, and the checkpoint
repair utility into the user OpenFOAM platform directories.

To clean the project, run:

```sh
./Allclean
```

## Multi-region selection

The modular CIGMA case uses:

```text
regionSolvers
{
    fluid  cigmaMPEFluid;
    plate1 cigmaMPESolid;
    plate2 cigmaMPESolid;
    wall1  cigmaMPESolid;
    wall2  cigmaMPESolid;
    wall3  cigmaMPESolid;
}
```

This baseline intentionally changes only module ownership and library loading;
it does not change condensation physics or CHT configuration.

## Baseline regression

The model equations, sign conventions, dimensions, and conservation targets
are defined in [`MODEL_DESIGN.md`](MODEL_DESIGN.md).

After sourcing OpenFOAM Foundation v14, run the standard-versus-modular
one-step regression test with:

```sh
./tests/run_baseline_regression.sh
```

The test copies only the initial and configuration directories to a temporary
location. It does not modify simulation results in either source case. Set
`KEEP_TEST_CASES=1` to retain the temporary cases and logs for inspection.

The reference environment and the remaining long-run validation requirements
are recorded in [`tests/BASELINE_MANIFEST.md`](tests/BASELINE_MANIFEST.md).

The single-phase baseline and active-condensation smoke test is run with:

```sh
./tests/run_single_phase_smoke.sh
```

It compares the no-condensation result against native `multicomponentFluid`
and then exercises `diffusionLayer`, `saturatedSteam`,
`condensingWallVelocity`, and `externalCondensationTemperature` for one time
step on an isolated temporary case.
