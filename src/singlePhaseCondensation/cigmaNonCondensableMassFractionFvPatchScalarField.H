/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM, Copyright (C) Forschungszentrum Juelich.
\*---------------------------------------------------------------------------*/

#ifndef cigmaNonCondensableMassFractionFvPatchScalarField_H
#define cigmaNonCondensableMassFractionFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{
namespace compressible
{

class cigmaNonCondensableMassFractionFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    void checkType() const;

public:

    TypeName("nonCondensableMassFraction");

    cigmaNonCondensableMassFractionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const dictionary&
    );

    cigmaNonCondensableMassFractionFvPatchScalarField
    (
        const cigmaNonCondensableMassFractionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const fieldMapper&
    );

    cigmaNonCondensableMassFractionFvPatchScalarField
    (
        const cigmaNonCondensableMassFractionFvPatchScalarField&,
        const DimensionedField<scalar, fvMesh>&
    );

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace compressible
} // End namespace Foam

#endif

// ************************************************************************* //
