/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2026 cigmaMPEFoam contributors
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "cigmaSinglePhaseFluid.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(cigmaSinglePhaseFluid, 0);
    addToRunTimeSelectionTable(solver, cigmaSinglePhaseFluid, fvMesh);
}
}

Foam::solvers::cigmaSinglePhaseFluid::cigmaSinglePhaseFluid(fvMesh& mesh)
:
    multicomponentFluid(mesh),
    condensation_(thermo_, thermophysicalTransport(), rho)
{}

Foam::solvers::cigmaSinglePhaseFluid::~cigmaSinglePhaseFluid()
{}


void Foam::solvers::cigmaSinglePhaseFluid::boundAndNormaliseMassFractions
(
    const bool reportCorrection
)
{
    tmp<volScalarField> tSolvedYSum
    (
        volScalarField::New
        (
            "cigmaSolvedYSum",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField& solvedYSum = tSolvedYSum.ref();

    forAll(Y_, i)
    {
        if (thermo_.solveSpecie(i))
        {
            Y_[i].max(scalar(0));
            solvedYSum += Y_[i];
        }
    }

    const dimensionedScalar maxSolvedY(max(solvedYSum));

    // OpenFOAM normaliseY() clips the default specie to zero but does not
    // reduce solved species whose sum exceeds one.  Scale only those cells,
    // preserving the relative composition of all solved species.
    if (maxSolvedY.value() > 1)
    {
        const tmp<volScalarField> tScale(max(solvedYSum, scalar(1)));

        forAll(Y_, i)
        {
            if (thermo_.solveSpecie(i))
            {
                Y_[i] /= tScale();
            }
        }
    }

    thermo_.normaliseY();

    if (reportCorrection && maxSolvedY.value() > 1 + small)
    {
        Info<< "Normalised solved-species maximum sum from "
            << maxSolvedY.value() << " to 1 before checkpoint write" << endl;
    }
}


void Foam::solvers::cigmaSinglePhaseFluid::prePredictor()
{
    multicomponentFluid::prePredictor();
    condensation_.predict();
}


void Foam::solvers::cigmaSinglePhaseFluid::thermophysicalPredictor()
{
    multicomponentFluid::thermophysicalPredictor();

    if (condensation_.active())
    {
        boundAndNormaliseMassFractions(false);
    }
}


void Foam::solvers::cigmaSinglePhaseFluid::postSolve()
{
    multicomponentFluid::postSolve();

    if (condensation_.active())
    {
        // Finish every checkpoint with evaluated solved-species boundaries.
        // This also resets their update state consistently for continuous and
        // restarted execution.
        forAll(Y_, i)
        {
            if (thermo_.solveSpecie(i))
            {
                Y_[i].correctBoundaryConditions();
            }
        }

        // Enforce the same bounded closure in continuous and restarted runs.
        boundAndNormaliseMassFractions(runTime.writeTime());
    }
}

// ************************************************************************* //
