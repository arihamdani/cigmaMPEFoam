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

Class
    Foam::solvers::cigmaMPEFluid

Description
    Runtime-selectable cigmaMPEFoam fluid module derived from the OpenFOAM
    Foundation v14 multiphaseEuler solver. It optionally applies per-phase
    Courant-number limits in addition to the standard global limit and can use
    single-phase-style multivariate convection for one multicomponent phase.

\*---------------------------------------------------------------------------*/

#ifndef cigmaMPEFluid_H
#define cigmaMPEFluid_H

#include "multiphaseEuler.H"

namespace Foam
{
namespace solvers
{

class cigmaMPEFluid
:
    public multiphaseEuler
{
    // Private data

        //- Use coupled multivariate convection for species and enthalpy
        Switch singlePhaseSpeciesConvection_;

        //- Multicomponent phase using single-phase-style species convection
        word singlePhaseSpeciesPhase_;


    // Private member functions

        //- Solve species with multivariate convection for the selected phase
        void compositionPredictor();

        //- Solve phase energy equations with optional multivariate convection
        void energyPredictor();


public:

    TypeName("cigmaMPEFluid");

    //- Construct from a region mesh
    explicit cigmaMPEFluid(fvMesh& mesh);

    //- Disallow copy construction
    cigmaMPEFluid(const cigmaMPEFluid&) = delete;

    //- Destructor
    virtual ~cigmaMPEFluid();

    //- Return the maximum time-step including optional per-phase Co limits
    virtual scalar maxDeltaT() const;

    //- Solve composition and energy with optional species-convection parity
    virtual void thermophysicalPredictor();

    //- Disallow assignment
    void operator=(const cigmaMPEFluid&) = delete;
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
