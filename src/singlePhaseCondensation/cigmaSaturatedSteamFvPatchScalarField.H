/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM's saturatedSteam boundary condition,
    Copyright (C) Forschungszentrum Juelich.

Class
    Foam::compressible::cigmaSaturatedSteamFvPatchScalarField

Description
    Mixed H2O mass-fraction boundary condition. It applies the saturated
    steam mass fraction when the adjacent gas is supersaturated at the wall
    temperature and otherwise applies zero gradient.
\*---------------------------------------------------------------------------*/

#ifndef cigmaSaturatedSteamFvPatchScalarField_H
#define cigmaSaturatedSteamFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{
namespace compressible
{

class cigmaSaturatedSteamFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    void checkType() const;

public:

    TypeName("saturatedSteam");

    cigmaSaturatedSteamFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const dictionary&
    );

    cigmaSaturatedSteamFvPatchScalarField
    (
        const cigmaSaturatedSteamFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const fvPatchFieldMapper&
    );

    cigmaSaturatedSteamFvPatchScalarField
    (
        const cigmaSaturatedSteamFvPatchScalarField&
    ) = delete;

    cigmaSaturatedSteamFvPatchScalarField
    (
        const cigmaSaturatedSteamFvPatchScalarField&,
        const DimensionedField<scalar, fvMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, fvMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new cigmaSaturatedSteamFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace compressible
} // End namespace Foam

#endif

// ************************************************************************* //
