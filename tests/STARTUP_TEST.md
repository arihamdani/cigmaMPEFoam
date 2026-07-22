# Modular startup test

Date: 2026-07-22  
OpenFOAM: Foundation v14, build `14-0db507470b5e`

## Method

The original standard case and the modular `cigmaMPEFoam` case were cloned to
a temporary directory. The following temporary controls were applied equally
to both clones:

```text
startFrom       startTime;
endTime         1e-6;
adjustTimeStep  no;
writeControl    timeStep;
writeInterval   1;
```

Each clone was run for one time step with:

```sh
foamMultiRun -case <temporary-case> -noFunctionObjects
```

Disabling function objects affected only diagnostic output. The region
solvers, CHT coupling, phase models, `phaseTurbulenceStabilisation`, and
`wallCondensation` remained active.

## Result

Both runs exited successfully at `1e-6 s`. The modular run selected
`cigmaMPEFluid` for `fluid` and `cigmaMPESolid` for all five solid regions. It
also selected `wallCondensation` and the standard `ArdenBuck` saturation model
from the user condensation library without duplicate-registration warnings.

Courant numbers, equation residuals, iteration counts, temperatures, phase
fractions, and wall-condensation output were identical to the standard run and
to the previous working `cigmaMPE` baseline.

Each run wrote 32 files. Thirty files were byte-identical. The only raw-file
differences were the expected `libs` metadata in `alphat.gas` and
`alphat.liquid`:

```text
libmultiphaseEulerFvModels.so
libcigmaMPECondensationModels.so
```

After normalising that library-name metadata, both `alphat` files were also
byte-identical. The standard-run aggregate hash remained the same as the
previous `cigmaMPE` baseline:

```text
ba5b54f0955b505ea23b39e1438183dd8a6fc640b56dffb92951280f05a9e280
```

Both runs emitted the same pre-existing wedge-patch planarity warnings. No new
warning, fatal error, or numerical difference was introduced by the modular
project.
