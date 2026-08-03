/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.
\*---------------------------------------------------------------------------*/

#include "cigmaExternalCondensationTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

void Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::getKappa
(
    scalarField& kappa,
    tmp<scalarField>& sumKappaTcByDelta,
    tmp<scalarField>& sumKappaByDelta,
    tmp<scalarField>& T,
    tmp<scalarField>& sumq
) const
{
    externalTemperatureFvPatchScalarField::getKappa
    (
        kappa,
        sumKappaTcByDelta,
        sumKappaByDelta,
        T,
        sumq
    );

    if
    (
        qcondName_ != word::null
     && db().foundObject<volScalarField>(qcondName_)
    )
    {
        const fvPatchScalarField& qcondCurrent =
            patch().lookupPatchField<volScalarField, scalar>(qcondName_);

        const scalarField qcond
        (
            qcondRelaxation_*qcondCurrent
          + (1 - qcondRelaxation_)*qcondPrevious_
        );

        qcondPrevious_ = qcond;
        plusEqOp(sumq, qcond);
    }
}


Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::
cigmaExternalCondensationTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    externalTemperatureFvPatchScalarField(p, iF, dict),
    qcondName_(dict.lookupOrDefault<word>("qcond", "qcond")),
    qcondRelaxation_
    (
        dict.lookupOrDefault<scalar>("qcondRelaxation", 1)
    ),
    qcondPrevious_
    (
        dict.found("qcondPrevious")
      ? scalarField
        (
            "qcondPrevious",
            dimPower/dimArea,
            dict,
            p.size()
        )
      : scalarField(p.size(), 0)
    )
{
    if (qcondRelaxation_ < 0 || qcondRelaxation_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "qcondRelaxation must be between zero and one"
            << exit(FatalIOError);
    }
}


Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::
cigmaExternalCondensationTemperatureFvPatchScalarField
(
    const cigmaExternalCondensationTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    externalTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    qcondName_(ptf.qcondName_),
    qcondRelaxation_(ptf.qcondRelaxation_),
    qcondPrevious_(mapper(ptf.qcondPrevious_)())
{}


Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::
cigmaExternalCondensationTemperatureFvPatchScalarField
(
    const cigmaExternalCondensationTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    externalTemperatureFvPatchScalarField(ptf, iF),
    qcondName_(ptf.qcondName_),
    qcondRelaxation_(ptf.qcondRelaxation_),
    qcondPrevious_(ptf.qcondPrevious_)
{}


void Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    externalTemperatureFvPatchScalarField::map(ptf, mapper);

    const cigmaExternalCondensationTemperatureFvPatchScalarField& cptf =
        refCast<const cigmaExternalCondensationTemperatureFvPatchScalarField>
        (
            ptf
        );

    mapper(qcondPrevious_, cptf.qcondPrevious_);
}


void Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    externalTemperatureFvPatchScalarField::reset(ptf);

    const cigmaExternalCondensationTemperatureFvPatchScalarField& cptf =
        refCast<const cigmaExternalCondensationTemperatureFvPatchScalarField>
        (
            ptf
        );

    qcondPrevious_.reset(cptf.qcondPrevious_);
}


void Foam::cigmaExternalCondensationTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    externalTemperatureFvPatchScalarField::write(os);
    writeEntry(os, "qcond", qcondName_);
    writeEntry(os, "qcondRelaxation", qcondRelaxation_);
    writeEntry(os, "qcondPrevious", qcondPrevious_);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        cigmaExternalCondensationTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
