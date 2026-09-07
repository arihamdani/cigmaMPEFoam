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
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "multivariateScheme.H"

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
    multiphaseEuler(mesh),
    singlePhaseSpeciesConvection_
    (
        pimple.dict().lookupOrDefault<Switch>
        (
            "singlePhaseSpeciesConvection",
            false
        )
    ),
    singlePhaseSpeciesPhase_
    (
        pimple.dict().lookupOrDefault<word>
        (
            "singlePhaseSpeciesPhase",
            "gas"
        )
    )
{
    if (singlePhaseSpeciesConvection_)
    {
        bool phaseFound = false;

        forAll(fluid.multicomponentPhases(), phasei)
        {
            phaseFound =
                phaseFound
             || fluid.multicomponentPhases()[phasei].name()
             == singlePhaseSpeciesPhase_;
        }

        if (!phaseFound)
        {
            FatalIOErrorInFunction(pimple.dict())
                << "singlePhaseSpeciesPhase " << singlePhaseSpeciesPhase_
                << " is not a multicomponent phase" << exit(FatalIOError);
        }

        Info<< "Single-phase-style multivariate species convection is ON"
            << nl << "    phase = " << singlePhaseSpeciesPhase_ << endl;
    }
}

Foam::solvers::cigmaMPEFluid::~cigmaMPEFluid()
{}


void Foam::solvers::cigmaMPEFluid::compositionPredictor()
{
    autoPtr<HashPtrTable<fvScalarMatrix>> popBalSpecieTransferPtr =
        populationBalanceSystem_.specieTransfer();
    HashPtrTable<fvScalarMatrix>& popBalSpecieTransfer =
        popBalSpecieTransferPtr();

    fluid_.correctReactions();

    forAll(fluid.multicomponentPhases(), multicomponentPhasei)
    {
        phaseModel& phase =
            fluid_.multicomponentPhases()[multicomponentPhasei];

        UPtrList<volScalarField>& Y = phase.YRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        const bool useMultivariateConvection =
            phase.name() == singlePhaseSpeciesPhase_;

        tmp<surfaceScalarField> talphaRhoPhi;
        tmp<fv::convectionScheme<scalar>> tmvConvection;
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

        if (useMultivariateConvection)
        {
            talphaRhoPhi = phase.alphaRhoPhi();
            const surfaceScalarField& alphaRhoPhi = talphaRhoPhi();

            forAll(Y, i)
            {
                fields.add(Y[i]);
            }
            fields.add(phase.thermo().he());

            const word schemeName
            (
                "div(" + alphaRhoPhi.name() + ",Yi_h)"
            );

            tmvConvection = fv::convectionScheme<scalar>::New
            (
                mesh_,
                fields,
                alphaRhoPhi,
                mesh_.schemes().div(schemeName)
            );
        }

        forAll(Y, i)
        {
            if (phase.solveSpecie(i))
            {
                fvScalarMatrix YiEqn
                (
                    phase.YiEqn(Y[i])
                 ==
                    popBalSpecieTransfer[Y[i].name()]
                  + fvModels().source(alpha, rho, Y[i])
                );

                if (useMultivariateConvection)
                {
                    const surfaceScalarField& alphaRhoPhi = talphaRhoPhi();

                    // Replace the native per-species convection matrix with
                    // the shared multivariate limiter used by
                    // multicomponentFluid.
                    YiEqn -= fvm::div
                    (
                        alphaRhoPhi,
                        Y[i],
                        "div(" + alphaRhoPhi.name() + ",Yi)"
                    );
                    YiEqn += tmvConvection->fvmDiv(alphaRhoPhi, Y[i]);
                }

                YiEqn.relax();
                fvConstraints().constrain(YiEqn);
                YiEqn.solve("Yi");
                fvConstraints().constrain(Y[i]);
            }
            else
            {
                Y[i].correctBoundaryConditions();
            }
        }
    }

    fluid_.correctSpecies();
}


void Foam::solvers::cigmaMPEFluid::energyPredictor()
{
    autoPtr<HashPtrTable<fvScalarMatrix>> heatTransferPtr =
        heatTransferSystem_.heatTransfer();
    HashPtrTable<fvScalarMatrix>& heatTransfer = heatTransferPtr();

    autoPtr<HashPtrTable<fvScalarMatrix>> popBalHeatTransferPtr =
        populationBalanceSystem_.heatTransfer();
    HashPtrTable<fvScalarMatrix>& popBalHeatTransfer =
        popBalHeatTransferPtr();

    forAll(fluid.thermalPhases(), thermalPhasei)
    {
        phaseModel& phase = fluid_.thermalPhases()[thermalPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        const bool useMultivariateConvection =
            phase.name() == singlePhaseSpeciesPhase_;

        tmp<surfaceScalarField> talphaRhoPhi;
        tmp<fv::convectionScheme<scalar>> tmvConvection;
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

        if (useMultivariateConvection)
        {
            talphaRhoPhi = phase.alphaRhoPhi();
            const surfaceScalarField& alphaRhoPhi = talphaRhoPhi();
            const UPtrList<volScalarField>& Y = phase.Y();

            forAll(Y, i)
            {
                fields.add(Y[i]);
            }
            fields.add(phase.thermo().he());

            const word schemeName
            (
                "div(" + alphaRhoPhi.name() + ",Yi_h)"
            );

            tmvConvection = fv::convectionScheme<scalar>::New
            (
                mesh_,
                fields,
                alphaRhoPhi,
                mesh_.schemes().div(schemeName)
            );
        }

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
            heatTransfer[phase.name()]
          + popBalHeatTransfer[phase.name()]
          + fvModels().source(alpha, rho, phase.thermo().he())
        );

        if (useMultivariateConvection)
        {
            const surfaceScalarField& alphaRhoPhi = talphaRhoPhi();
            volScalarField& he = phase.thermo().he();

            EEqn -= fvm::div(alphaRhoPhi, he);
            EEqn += tmvConvection->fvmDiv(alphaRhoPhi, he);
        }

        EEqn.relax();
        fvConstraints().constrain(EEqn);
        EEqn.solve();
        fvConstraints().constrain(phase.thermo().he());
    }

    fluid_.correctThermo();
    fluid_.correctContinuityError(populationBalanceSystem_.dmdts());
}


void Foam::solvers::cigmaMPEFluid::thermophysicalPredictor()
{
    if (!singlePhaseSpeciesConvection_)
    {
        multiphaseEuler::thermophysicalPredictor();
        return;
    }

    for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
    {
        fluid_.predictThermophysicalTransport();
        compositionPredictor();
        energyPredictor();

        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];

            Info<< phase.name() << " min/max T "
                << min(phase.thermo().T()).value()
                << " - " << max(phase.thermo().T()).value() << endl;
        }
    }
}


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
