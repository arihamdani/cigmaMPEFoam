/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM, Copyright (C) Forschungszentrum Juelich.
\*---------------------------------------------------------------------------*/

#include "cigmaCondensingWallVelocityFvPatchVectorField.H"
#include "cigmaDiffusionLayerCondensation.H"
#include "addToRunTimeSelectionTable.H"

Foam::cigmaCondensingWallVelocityFvPatchVectorField::
cigmaCondensingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false)
{
    fvPatchField<vector>::operator=
    (
        dict.found("value")
      ? vectorField("value", dict, p.size())
      : vectorField(p.size(), Zero)
    );
}


Foam::cigmaCondensingWallVelocityFvPatchVectorField::
cigmaCondensingWallVelocityFvPatchVectorField
(
    const cigmaCondensingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::cigmaCondensingWallVelocityFvPatchVectorField::
cigmaCondensingWallVelocityFvPatchVectorField
(
    const cigmaCondensingWallVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}


void Foam::cigmaCondensingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const cigmaDiffusionLayerCondensation& condensation =
        db().lookupObject<cigmaDiffusionLayerCondensation>
        (
            cigmaDiffusionLayerCondensation::dictionaryName
        );

    operator==
    (
        patch().nf()*condensation.condSuctionVelocity(patch().index())
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::cigmaCondensingWallVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "value", *this);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cigmaCondensingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
