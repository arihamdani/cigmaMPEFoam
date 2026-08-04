/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM, Copyright (C) Forschungszentrum Juelich.
\*---------------------------------------------------------------------------*/

#include "cigmaNonCondensableMassFractionFvPatchScalarField.H"
#include "cigmaDiffusionLayerCondensation.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

namespace Foam
{
namespace compressible
{

void cigmaNonCondensableMassFractionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch " << patch().name() << " must be a wall, not "
            << patch().type() << exit(FatalError);
    }

    if (internalField().name() == "H2O")
    {
        FatalErrorInFunction
            << typeName << " must not be applied to H2O"
            << exit(FatalError);
    }
}


cigmaNonCondensableMassFractionFvPatchScalarField::
cigmaNonCondensableMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict, false)
{
    checkType();

    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = Zero;
    }
}


cigmaNonCondensableMassFractionFvPatchScalarField::
cigmaNonCondensableMassFractionFvPatchScalarField
(
    const cigmaNonCondensableMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


cigmaNonCondensableMassFractionFvPatchScalarField::
cigmaNonCondensableMassFractionFvPatchScalarField
(
    const cigmaNonCondensableMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


void cigmaNonCondensableMassFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label patchi = patch().index();

    const cigmaDiffusionLayerCondensation& condensation =
        mesh.lookupObject<cigmaDiffusionLayerCondensation>
        (
            cigmaDiffusionLayerCondensation::dictionaryName
        );

    const volScalarField& Yi =
        mesh.lookupObject<volScalarField>(internalField().name());

    const tmp<scalarField> tDEff =
        condensation.wallEffectiveMassDiffusivity(Yi, patchi);

    gradient() =
        condensation.condSuctionVelocity(patchi)
       *(*this)
       /max(tDEff(), scalar(SMALL));

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void cigmaNonCondensableMassFractionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


makePatchTypeField
(
    fvPatchScalarField,
    cigmaNonCondensableMassFractionFvPatchScalarField
);

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
