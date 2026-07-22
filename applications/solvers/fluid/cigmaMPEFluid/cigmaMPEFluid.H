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
    Runtime-selectable cigmaMPEFoam fluid module. The initial implementation
    inherits the OpenFOAM Foundation v14 multiphaseEuler solver without
    overriding any equation or solution-stage method.

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
public:

    TypeName("cigmaMPEFluid");

    //- Construct from a region mesh
    explicit cigmaMPEFluid(fvMesh& mesh);

    //- Disallow copy construction
    cigmaMPEFluid(const cigmaMPEFluid&) = delete;

    //- Destructor
    virtual ~cigmaMPEFluid();

    //- Disallow assignment
    void operator=(const cigmaMPEFluid&) = delete;
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
