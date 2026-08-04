/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | cigmaFOAM
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef cigmaBuoyantKOmegaSST_H
#define cigmaBuoyantKOmegaSST_H

#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcGrad.H"
#include "fvcMeshPhi.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "kOmegaSST.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class cigmaBuoyantKOmegaSST
:
    public kOmegaSST<BasicMomentumTransportModel>
{
protected:

    dictionary buoyancyDict_;
    word buoyancyModel_;
    Switch buoyancyDissipation_;
    Switch directionalDissipation_;

    dimensionedScalar C3_;
    dimensionedScalar sigmaRho_;

    volScalarField C3e_;
    volScalarField PkBbyNut_;

    void validateBuoyancySettings() const;

    const uniformDimensionedVectorField& g() const;

    void updateBuoyancyProduction();

    tmp<volScalarField::Internal> PwB
    (
        const volScalarField::Internal& gamma
    );


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef eddyViscosity<RASModel<BasicMomentumTransportModel>>
        BaseMomentumTransportModel;

    TypeName("cigmaBuoyantKOmegaSST");

    cigmaBuoyantKOmegaSST
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~cigmaBuoyantKOmegaSST()
    {}

    virtual bool read();

    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "cigmaBuoyantKOmegaSST.C"
#endif

#endif

// ************************************************************************* //
