/*---------------------------------------------------------------------------*\
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM's condensingWallVelocity boundary condition,
    Copyright (C) Forschungszentrum Juelich.
\*---------------------------------------------------------------------------*/

#ifndef cigmaCondensingWallVelocityFvPatchVectorField_H
#define cigmaCondensingWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class cigmaCondensingWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    TypeName("condensingWallVelocity");

    cigmaCondensingWallVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, fvMesh>&,
        const dictionary&
    );

    cigmaCondensingWallVelocityFvPatchVectorField
    (
        const cigmaCondensingWallVelocityFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, fvMesh>&,
        const fvPatchFieldMapper&
    );

    cigmaCondensingWallVelocityFvPatchVectorField
    (
        const cigmaCondensingWallVelocityFvPatchVectorField&
    ) = delete;

    cigmaCondensingWallVelocityFvPatchVectorField
    (
        const cigmaCondensingWallVelocityFvPatchVectorField&,
        const DimensionedField<vector, fvMesh>&
    );

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, fvMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new cigmaCondensingWallVelocityFvPatchVectorField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
