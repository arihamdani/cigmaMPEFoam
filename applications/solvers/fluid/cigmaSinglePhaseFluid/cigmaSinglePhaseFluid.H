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

Class
    Foam::solvers::cigmaSinglePhaseFluid

Description
    Runtime-selectable single-phase, multicomponent fluid module for
    cigmaMPEFoam. The baseline equations are inherited unchanged from the
    OpenFOAM Foundation v14 multicomponentFluid module.

\*---------------------------------------------------------------------------*/

#ifndef cigmaSinglePhaseFluid_H
#define cigmaSinglePhaseFluid_H

#include "multicomponentFluid.H"
#include "cigmaDiffusionLayerCondensation.H"

namespace Foam
{
namespace solvers
{

class cigmaSinglePhaseFluid
:
    public multicomponentFluid
{
    // Private data

        cigmaDiffusionLayerCondensation condensation_;


    // Private member functions

        //- Bound solved species and rebuild the algebraic default specie
        void boundAndNormaliseMassFractions(const bool reportCorrection);


public:

    TypeName("cigmaSinglePhaseFluid");

    //- Construct from a region mesh
    explicit cigmaSinglePhaseFluid(fvMesh& mesh);

    //- Disallow copy construction
    cigmaSinglePhaseFluid(const cigmaSinglePhaseFluid&) = delete;

    //- Destructor
    virtual ~cigmaSinglePhaseFluid();

    //- Update the wall-condensation flux before solving the flow equations
    virtual void prePredictor();

    //- Solve thermophysics and enforce a conservative mass-fraction closure
    virtual void thermophysicalPredictor();

    //- Finalise species boundaries and checkpoint mass-fraction closure
    virtual void postSolve();

    //- Disallow assignment
    void operator=(const cigmaSinglePhaseFluid&) = delete;
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
