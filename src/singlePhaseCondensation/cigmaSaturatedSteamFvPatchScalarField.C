/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM, Copyright (C) Forschungszentrum Juelich.
\*---------------------------------------------------------------------------*/

#include "cigmaSaturatedSteamFvPatchScalarField.H"
#include "cigmaDiffusionLayerCondensation.H"
#include "addToRunTimeSelectionTable.H"
#include "fluidMulticomponentThermo.H"
#include "wallFvPatch.H"

namespace Foam
{
namespace compressible
{

void cigmaSaturatedSteamFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch " << patch().name() << " must be a wall, not "
            << patch().type() << exit(FatalError);
    }

    if (internalField().name() != "H2O")
    {
        FatalErrorInFunction
            << typeName << " can only be applied to the H2O field"
            << exit(FatalError);
    }
}


cigmaSaturatedSteamFvPatchScalarField::
cigmaSaturatedSteamFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false)
{
    checkType();

    if (dict.found("value"))
    {
        refValue() = scalarField("value", dict, p.size());
    }
    else
    {
        refValue() = patchInternalField();
    }
    refGrad() = Zero;
    valueFraction() = Zero;
    fvPatchScalarField::operator=(patchInternalField());
}


cigmaSaturatedSteamFvPatchScalarField::
cigmaSaturatedSteamFvPatchScalarField
(
    const cigmaSaturatedSteamFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    // Mapping can use temporary internal fields during reconstruction, so the
    // dictionary-only field-name validation must not be repeated here.
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


cigmaSaturatedSteamFvPatchScalarField::
cigmaSaturatedSteamFvPatchScalarField
(
    const cigmaSaturatedSteamFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{}


void cigmaSaturatedSteamFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label patchi = patch().index();

    const fluidMulticomponentThermo& thermo =
        mesh.lookupObject<fluidMulticomponentThermo>
        (
            physicalProperties::typeName
        );

    const cigmaDiffusionLayerCondensation& condensation =
        mesh.lookupObject<cigmaDiffusionLayerCondensation>
        (
            cigmaDiffusionLayerCondensation::dictionaryName
        );

    const label H2Oi = condensation.condensingSpeciei();
    const PtrList<volScalarField>& Y = thermo.Y();
    const scalar WH2O = thermo.Wi(H2Oi).value();

    const scalarField Tw(thermo.T().boundaryField()[patchi]);
    const scalarField pSat(condensation.pSat(Tw));
    const scalarField pCell
    (
        thermo.p().boundaryField()[patchi].patchInternalField()
    );
    const scalarField pWall(thermo.p().boundaryField()[patchi]);

    scalarField sumYByW(patch().size(), 0);
    scalarField sumXWnc(patch().size(), 0);
    scalarField sumXnc(patch().size(), 0);

    forAll(Y, i)
    {
        const scalar Wi = thermo.Wi(i).value();
        sumYByW +=
            Y[i].boundaryField()[patchi].patchInternalField()/Wi;
    }

    const scalarField xH2OCell
    (
        Y[H2Oi].boundaryField()[patchi].patchInternalField()
       /WH2O/max(sumYByW, scalar(SMALL))
    );

    forAll(Y, i)
    {
        if (i == H2Oi)
        {
            continue;
        }

        const scalar Wi = thermo.Wi(i).value();
        const scalarField xi
        (
            Y[i].boundaryField()[patchi].patchInternalField()
           /Wi/max(sumYByW, scalar(SMALL))
        );

        sumXnc += xi;
        sumXWnc += xi*Wi;
    }

    const scalarField Wnc(sumXWnc/max(sumXnc, scalar(SMALL)));
    const scalarField condensing(pos(xH2OCell*pCell - pSat));
    const scalarField xSat
    (
        min(max(pSat/max(pWall, scalar(SMALL)), scalar(0)), scalar(1))
    );
    const scalarField WSat(xSat*WH2O + (1 - xSat)*Wnc);
    const scalarField YSat(xSat*WH2O/max(WSat, scalar(SMALL)));

    refValue() = condensing*YSat + (1 - condensing)*patchInternalField();
    valueFraction() = condensing;
    refGrad() = Zero;

    mixedFvPatchScalarField::updateCoeffs();
}


void cigmaSaturatedSteamFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


makePatchTypeField
(
    fvPatchScalarField,
    cigmaSaturatedSteamFvPatchScalarField
);

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
