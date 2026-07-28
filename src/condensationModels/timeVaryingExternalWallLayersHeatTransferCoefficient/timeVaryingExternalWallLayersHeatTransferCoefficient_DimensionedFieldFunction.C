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

\*---------------------------------------------------------------------------*/

#include "timeVaryingExternalWallLayersHeatTransferCoefficient_DimensionedFieldFunction.H"
#include "DimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace DimensionedFieldFunctions
{
    defineTypeNameAndDebug
    (
        timeVaryingExternalWallLayersHeatTransferCoefficient,
        0
    );

    typedef DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>
        scalarFvPatchDimensionedFieldFunction;

    addToRunTimeSelectionTable
    (
        scalarFvPatchDimensionedFieldFunction,
        timeVaryingExternalWallLayersHeatTransferCoefficient,
        dictionary
    );
}
}

Foam::DimensionedFieldFunctions::
timeVaryingExternalWallLayersHeatTransferCoefficient::
timeVaryingExternalWallLayersHeatTransferCoefficient
(
    const dictionary& dict,
    DimensionedField<scalar, fvPatch>& field
)
:
    wallLayersHeatTransferCoefficient(dict, field),
    h_
    (
        Function1<scalar>::New
        (
            "h",
            field.time().userUnits(),
            dimPower/dimArea/dimTemperature,
            dict
        )
    )
{}

Foam::DimensionedFieldFunctions::
timeVaryingExternalWallLayersHeatTransferCoefficient::
timeVaryingExternalWallLayersHeatTransferCoefficient
(
    const timeVaryingExternalWallLayersHeatTransferCoefficient& dff,
    DimensionedField<scalar, fvPatch>& field
)
:
    wallLayersHeatTransferCoefficient(dff, field),
    h_(dff.h_, false)
{}

Foam::autoPtr
<
    Foam::DimensionedFieldFunction
    <
        Foam::DimensionedField<Foam::scalar, Foam::fvPatch>
    >
>
Foam::DimensionedFieldFunctions::
timeVaryingExternalWallLayersHeatTransferCoefficient::clone
(
    DimensionedField<scalar, fvPatch>& field
) const
{
    return autoPtr
    <
        DimensionedFieldFunction<DimensionedField<scalar, fvPatch>>
    >
    (
        new timeVaryingExternalWallLayersHeatTransferCoefficient(*this, field)
    );
}

void Foam::DimensionedFieldFunctions::
timeVaryingExternalWallLayersHeatTransferCoefficient::evaluate()
{
    wallLayersHeatTransferCoefficient::evaluate();

    const scalar externalH = h_->value(this->field_.time().value());

    // Product-over-sum is equivalent to the resistance formulation used by
    // the native model and remains finite for the supplied h(0) = 0 tables.
    this->field_.primitiveFieldRef() =
        externalH*this->field_.primitiveField()
       /(externalH + this->field_.primitiveField() + vSmall);
}

bool Foam::DimensionedFieldFunctions::
timeVaryingExternalWallLayersHeatTransferCoefficient::update()
{
    evaluate();
    return true;
}

void Foam::DimensionedFieldFunctions::
timeVaryingExternalWallLayersHeatTransferCoefficient::write
(
    Ostream& os
) const
{
    wallLayersHeatTransferCoefficient::write(os);
    writeEntry
    (
        os,
        this->field_.time().userUnits(),
        dimPower/dimArea/dimTemperature,
        h_()
    );
}

// ************************************************************************* //
