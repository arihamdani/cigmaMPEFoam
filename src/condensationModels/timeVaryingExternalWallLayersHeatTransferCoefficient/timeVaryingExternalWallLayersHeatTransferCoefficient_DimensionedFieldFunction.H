/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::DimensionedFieldFunctions::
    timeVaryingExternalWallLayersHeatTransferCoefficient

Description
    Heat-transfer coefficient for conductive wall layers in series with a
    time-varying external coefficient supplied as a Function1.

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingExternalWallLayersHeatTransferCoefficient_DimensionedFieldFunction_H
#define timeVaryingExternalWallLayersHeatTransferCoefficient_DimensionedFieldFunction_H

#include "wallLayersHeatTransferCoefficient_DimensionedFieldFunction.H"
#include "Function1.H"

namespace Foam
{
namespace DimensionedFieldFunctions
{

class timeVaryingExternalWallLayersHeatTransferCoefficient
:
    public wallLayersHeatTransferCoefficient
{
    //- External heat-transfer coefficient [W/m2/K]
    autoPtr<Function1<scalar>> h_;

public:

    TypeName("timeVaryingExternalWallLayersHeatTransferCoefficient");

    timeVaryingExternalWallLayersHeatTransferCoefficient
    (
        const dictionary& dict,
        DimensionedField<scalar, fvPatch>& field
    );

    timeVaryingExternalWallLayersHeatTransferCoefficient
    (
        const timeVaryingExternalWallLayersHeatTransferCoefficient& dff,
        DimensionedField<scalar, fvPatch>& field
    );

    virtual autoPtr
    <
        DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>
    > clone(DimensionedField<scalar, fvPatch>& field) const;

    virtual ~timeVaryingExternalWallLayersHeatTransferCoefficient()
    {}

    virtual void evaluate();

    virtual bool update();

    virtual void write(Ostream& os) const;
};

} // End namespace DimensionedFieldFunctions
} // End namespace Foam

#endif

// ************************************************************************* //
