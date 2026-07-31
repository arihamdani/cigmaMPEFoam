# cigmaMPEFoam

`cigmaMPEFoam` is a modular OpenFOAM Foundation v14 project for multi-region
conjugate heat-transfer simulations with multiphase wall condensation.

## Project layout

```text
cigmaMPEFoam/
├── applications/
│   └── solvers/
│       ├── fluid/
│       │   └── cigmaMPEFluid/
│       └── solid/
│           └── cigmaMPESolid/
├── src/
│   └── condensationModels/
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

The modules are installed as `libcigmaMPEFluidSolver.so` and
`libcigmaMPESolidSolver.so`. Their names follow the
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

## Build

Source OpenFOAM Foundation v14 and run:

```sh
./Allwmake
```

The build script first runs `Allclean`, then builds the condensation library,
the fluid solver module, and the solid solver module into `$FOAM_USER_LIBBIN`.

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
