/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 cigmaMPEFoam contributors
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    GNU General Public License, version 3 or later.

Class
    Foam::cigmaExternalCondensationTemperatureFvPatchScalarField

Description
    External-temperature boundary condition which adds the single-phase
    condensation heat-flux field to the external wall heat balance.
\*---------------------------------------------------------------------------*/

#ifndef cigmaExternalCondensationTemperatureFvPatchScalarField_H
#define cigmaExternalCondensationTemperatureFvPatchScalarField_H

#include "externalTemperatureFvPatchScalarField.H"

namespace Foam
{

class cigmaExternalCondensationTemperatureFvPatchScalarField
:
    public externalTemperatureFvPatchScalarField
{
    // Private data

        const word qcondName_;

        scalar qcondRelaxation_;

        mutable scalarField qcondPrevious_;


protected:

    virtual void getKappa
    (
        scalarField& kappa,
        tmp<scalarField>& sumKappaTcByDelta,
        tmp<scalarField>& sumKappaByDelta,
        tmp<scalarField>& T,
        tmp<scalarField>& sumq
    ) const;


public:

    TypeName("externalCondensationTemperature");

    cigmaExternalCondensationTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const dictionary&
    );

    cigmaExternalCondensationTemperatureFvPatchScalarField
    (
        const cigmaExternalCondensationTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const fieldMapper&
    );

    cigmaExternalCondensationTemperatureFvPatchScalarField
    (
        const cigmaExternalCondensationTemperatureFvPatchScalarField&
    ) = delete;

    cigmaExternalCondensationTemperatureFvPatchScalarField
    (
        const cigmaExternalCondensationTemperatureFvPatchScalarField&,
        const DimensionedField<scalar, fvMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, fvMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new cigmaExternalCondensationTemperatureFvPatchScalarField
            (
                *this,
                iF
            )
        );
    }

    virtual void map(const fvPatchScalarField&, const fieldMapper&);

    virtual void reset(const fvPatchScalarField&);

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
