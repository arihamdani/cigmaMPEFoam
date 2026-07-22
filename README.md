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
`Foam::solvers::multiphaseEuler` module. It initially overrides no equations or
solution-stage methods.

`cigmaMPESolid` is a runtime-selectable subclass of the standard v14
`Foam::solvers::solid` module. It initially overrides no equations or
solution-stage methods.

The modules are installed as `libcigmaMPEFluidSolver.so` and
`libcigmaMPESolidSolver.so`. Their names follow the
`lib<solverName>Solver.so` convention used automatically by `foamMultiRun`.

## Condensation models

`src/condensationModels` contains the minimum OpenFOAM v14 source set required
by the CIGMA wall-condensation case:

- `wallCondensation` and its phase-change-rate patch field;
- `wallPhaseChange` and the `alphatPhaseChangeWallFunction` patch field;
- `phaseTurbulenceStabilisation`, which is also selected by the case.

These files are unchanged copies of the OpenFOAM Foundation v14 sources. The
runtime type names, equations, correlations, coefficients, and dictionary
interfaces are therefore unchanged. The selected `ArdenBuck` saturation
pressure model remains supplied by the standard OpenFOAM v14 saturation-model
library and is not duplicated here.

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
