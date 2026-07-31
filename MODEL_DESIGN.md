# cigmaMPEFoam model design

## 1. Purpose and baseline

`cigmaMPEFoam` is an OpenFOAM Foundation v14 research solver for multiphase,
multi-region conjugate heat transfer with wall condensation. The baseline is
defined by project commit `e340ee0313410a7739203ae862435be3277f1839` and
OpenFOAM commit `0db507470b`.

The fluid and solid solver modules are deliberately thin subclasses. They
inherit the equations and solution sequence of `multiphaseEuler` and `solid`,
respectively. The copied `phaseTurbulenceStabilisation`,
`massDiffusionLimitedPhaseChange`, `wallPhaseChange`, and `wallCondensation`
sources are identical to the corresponding OpenFOAM v14 sources. The
time-varying external heat-transfer-coefficient model is an additional
boundary utility and does not alter the condensation closure.

Future model development must preserve runtime selection of the baseline
implementation so that every extension can be compared with an unchanged
reference model.

## 2. Phase and sign conventions

The phase order in the wall-condensation dictionary is significant:

```foam
phases (gas liquid);
specie H2O;
```

The first phase is the vapour-containing donor phase and the second phase is
the liquid receiver phase. The wall-condensation model is one-way:

- `mDot > 0` means condensation from gas to liquid;
- the gas phase receives a mass sink;
- the liquid phase receives an equal mass source;
- evaporation is suppressed by clipping `mDot` to zero.

The volumetric condensation source has dimensions

```text
[mass length^-3 time^-1] = [1 -3 -1 0 0 0 0].
```

The integrated condensation rate is

```text
MdotCond = integral_V(mDot dV) [kg/s].
```

## 3. Baseline wall-condensation closure

For each active wall face, the OpenFOAM v14 model evaluates the effective
vapour diffusivity

```text
DEff = Dm + (Prt/Sct) alphatConv,V.
```

The adjacent-cell and wall vapour mole fractions are

```text
xc = Yc (Wi/Wmix)^-1,
xw = psat(Tw)/p.
```

Using `AbyV = Awall/Vcell` and the wall-normal inverse distance `delta`, the
volumetric condensation rate is

```text
mDot = -alphaV AbyV (Wi/Wmix) DEff delta
       log(max(1 - xc, 0.001)/max(1 - xw, 0.001)).
```

The final source is limited to condensation:

```text
mDot = max(mDot, 0).
```

The latent wall heat flux used by the thermal wall function is

```text
qCond = mDot L/AbyV [W/m2].
```

The baseline defaults are `Prt = 0.85`, `Sct = 0.7`, and
`specieSemiImplicit = false`, unless overridden by the case dictionary.

## 4. Source coupling and conservation requirements

All new implementations must satisfy the following discrete requirements:

```text
Sgas,mass + Sliquid,mass = 0
Sgas,H2O  + Sliquid,H2O  = 0
sum_i Ji = 0
```

The latent-energy source must use the same phase-change rate as the mass and
species equations. No model may independently recompute a different `mDot`
for the energy equation.

For a closed validation domain, the monitored conservation errors are

```text
epsilonMass = abs(M(t) - M(0))/max(abs(M(0)), small)

epsilonEnergy = abs(E(t) - E(0) - Qboundary(t))
                /max(abs(E(0)), abs(Qboundary(t)), small).
```

The initial development targets are:

- phase-transfer imbalance below `1e-10` of the integrated transfer rate;
- total mass drift below `1e-8` in the one-step regression test;
- no degradation of the validated long-run mass and energy errors relative to
  the baseline;
- identical results when every new model is disabled.

The final tolerances for the paper must be fixed before the validation sweep
and must not be tuned after examining the new-model results.

## 5. Numerical design rules

1. Every new closure is runtime-selectable and has an explicit baseline mode.
2. Dictionary defaults reproduce the OpenFOAM v14 baseline.
3. A limiter must preserve sign, dimensions, and donor availability.
4. Linearised source terms must provide a non-negative diagonal contribution.
5. Expensive multicomponent coefficients are cached and updated only when a
   documented temperature, pressure, or composition threshold is exceeded.
6. Diagnostic fields are optional and default to `NO_WRITE`.
7. Parallel reductions are performed only for user-requested diagnostics.
8. Model changes are accepted only after build, startup, regression,
   conservation, mesh, time-step, and MPI checks.

### 5.1 Phase-specific Courant-number limiter

For each configured moving phase `i`, `cigmaMPEFluid` evaluates

```text
Co_i = 0.5 deltaT max_cells(sum_faces(abs(phi_i))/V).
```

The phase-specific time-step candidate is

```text
deltaT_i = maxCo_i/Co_i deltaT.
```

The applied time-step is the minimum of the standard OpenFOAM global
time-step limit and every configured `deltaT_i`. The optional control is

```text
phaseMaxCo
{
    gas     5;
    liquid  1;
}
```

If `phaseMaxCo` is absent, the implementation returns the unmodified standard
OpenFOAM v14 `basicFluidSolver::maxDeltaT()` result. This disabled-state
behaviour is protected by the baseline regression test.

## 6. Regression invariants

The one-step standard-versus-modular regression test must satisfy all of the
following:

- both simulations exit successfully;
- `multiphaseEuler` and `cigmaMPEFluid` execute the same fluid equations;
- `solid` and `cigmaMPESolid` execute the same solid equations;
- the number and names of written fields are identical;
- all written fields are byte-identical after normalising only the expected
  condensation-library name in the `alphat` boundary metadata;
- solver iterations, residuals, and wall-condensation diagnostics are
  identical after normalising solver names and runtime metadata;
- no additional warning or fatal error is introduced.

The automated implementation is `tests/run_baseline_regression.sh`.

## 7. Validation data contract

Each long validation run supplied for model comparison must archive:

- case Git commit and OpenFOAM build identifier;
- mesh checksum and cell count;
- dictionary checksum;
- total gas mass and total liquid mass versus time;
- domain-integrated wall-condensation rate versus time;
- total fluid and solid energy versus time;
- pressure, water-species, and energy residual histories;
- wall-clock time, peak resident memory, process count, and write settings;
- experimental observables used by the paper.

Comparisons are invalid if mesh, boundary conditions, numerical controls, or
output intervals differ between baseline and candidate cases, unless the
specific sensitivity study declares that difference in advance.

## 8. Planned extension boundary

The first new transport model will be implemented without changing the public
dictionary contract of the baseline `wallCondensation` model. New source code,
class names, dictionary keys, tests, comments, logs, and technical
documentation must be written in English.

The baseline remains available throughout development for an A/B comparison:

```text
OpenFOAM v14 wallCondensation  -> reference
cigmaMPEFoam extended model    -> candidate
```
