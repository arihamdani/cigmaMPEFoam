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


void Foam::solvers::cigmaSinglePhaseFluid::prePredictor()
{
    multicomponentFluid::prePredictor();
    condensation_.predict();
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

        // Rebuild only the algebraic default specie before runTime.write().
        thermo_.normaliseY();
    }
}

// ************************************************************************* //
