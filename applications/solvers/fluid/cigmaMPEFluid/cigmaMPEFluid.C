/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2026 cigmaMPEFoam contributors
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cigmaMPEFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSurfaceIntegrate.H"

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(cigmaMPEFluid, 0);
    addToRunTimeSelectionTable(solver, cigmaMPEFluid, fvMesh);
}
}

Foam::solvers::cigmaMPEFluid::cigmaMPEFluid(fvMesh& mesh)
:
    multiphaseEuler(mesh)
{}

Foam::solvers::cigmaMPEFluid::~cigmaMPEFluid()
{}


Foam::scalar Foam::solvers::cigmaMPEFluid::maxDeltaT() const
{
    scalar deltaT = basicFluidSolver::maxDeltaT();

    const dictionary& controlDict = runTime.controlDict();

    if (!controlDict.found("phaseMaxCo"))
    {
        return deltaT;
    }

    const dictionary& phaseMaxCoDict = controlDict.subDict("phaseMaxCo");

    Info<< "Phase Courant numbers and limits:";

    forAll(movingPhases_, phasei)
    {
        const phaseModel& phase = movingPhases_[phasei];

        if (!phaseMaxCoDict.found(phase.name()))
        {
            continue;
        }

        const scalar maxPhaseCo =
            phaseMaxCoDict.lookup<scalar>(phase.name());

        if (maxPhaseCo <= 0)
        {
            FatalIOErrorInFunction(phaseMaxCoDict)
                << "The Courant-number limit for phase " << phase.name()
                << " must be greater than zero, but is " << maxPhaseCo
                << exit(FatalIOError);
        }

        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(phase.phi()))().primitiveField()
        );

        const scalar phaseCo =
            0.5
           *gMax(sumPhi/mesh.V().primitiveField())
           *runTime.deltaTValue();

        Info<< ' ' << phase.name() << '=' << phaseCo
            << " (limit " << maxPhaseCo << ')';

        if (phaseCo > small)
        {
            deltaT = min
            (
                deltaT,
                maxPhaseCo/phaseCo*runTime.deltaTValue()
            );
        }
    }

    Info<< endl;

    return deltaT;
}

// ************************************************************************* //
