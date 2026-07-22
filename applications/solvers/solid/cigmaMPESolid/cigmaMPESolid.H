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
    Foam::solvers::cigmaMPESolid

Description
    Runtime-selectable cigmaMPEFoam solid module. The initial implementation
    inherits the OpenFOAM Foundation v14 solid solver without overriding any
    equation or solution-stage method.

\*---------------------------------------------------------------------------*/

#ifndef cigmaMPESolid_H
#define cigmaMPESolid_H

#include "solid.H"

namespace Foam
{
namespace solvers
{

class cigmaMPESolid
:
    public solid
{
public:

    TypeName("cigmaMPESolid");

    //- Construct from a region mesh
    explicit cigmaMPESolid(fvMesh& mesh);

    //- Disallow copy construction
    cigmaMPESolid(const cigmaMPESolid&) = delete;

    //- Destructor
    virtual ~cigmaMPESolid();

    //- Disallow assignment
    void operator=(const cigmaMPESolid&) = delete;
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
