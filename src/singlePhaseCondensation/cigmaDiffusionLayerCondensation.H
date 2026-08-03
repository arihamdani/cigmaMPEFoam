/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 cigmaMPEFoam contributors
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM and is distributed under the GNU General
    Public License, version 3 or later.

    This implementation adapts the diffusion-layer wall-condensation model
    from containmentFOAM (Forschungszentrum Juelich) to the native OpenFOAM
    Foundation v14 multicomponent thermophysical-transport interface.

Class
    Foam::cigmaDiffusionLayerCondensation

Description
    Single-phase wall-condensation model based on the diffusion-layer
    formulation of Vijaya Kumar et al. (2021). The condensate mass flux is
    obtained from the native v14 effective H2O diffusion flux and the Stefan
    correction, mDot = -max(jH2O, 0)/(1 - YH2O).

\*---------------------------------------------------------------------------*/

#ifndef cigmaDiffusionLayerCondensation_H
#define cigmaDiffusionLayerCondensation_H

#include "IOdictionary.H"
#include "fluidMulticomponentThermo.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "saturationPressureModel.H"
#include "volFields.H"

namespace Foam
{

class cigmaDiffusionLayerCondensation
:
    public IOdictionary
{
    // Private data

        const fluidMulticomponentThermo& thermo_;

        const fluidMulticomponentThermophysicalTransportModel&
            thermophysicalTransport_;

        const volScalarField& rho_;

        Switch active_;

        word condensingSpecie_;

        label condensingSpeciei_;

        autoPtr<saturationPressureModel> saturationPressure_;

        scalar maxCondWallMassFlux_;

        scalar latentHeatA_;
        scalar latentHeatB_;
        scalar latentHeatC_;
        scalar latentHeatD_;

        volScalarField massTransferRate_;

        volScalarField qCondensation_;

        volScalarField::Boundary condSuctionVelocity_;


    // Private member functions

        IOobject createIOobject(const fvMesh& mesh) const;

        autoPtr<saturationPressureModel> createSaturationModel() const;

        void validate() const;


public:

    TypeName("cigmaDiffusionLayerCondensation");

    static const word dictionaryName;


    // Constructors

        cigmaDiffusionLayerCondensation
        (
            const fluidMulticomponentThermo& thermo,
            const fluidMulticomponentThermophysicalTransportModel&
                thermophysicalTransport,
            const volScalarField& rho
        );

        cigmaDiffusionLayerCondensation
        (
            const cigmaDiffusionLayerCondensation&
        ) = delete;


    //- Destructor
    virtual ~cigmaDiffusionLayerCondensation();


    // Member functions

        bool active() const
        {
            return active_;
        }

        const fluidMulticomponentThermo& thermo() const
        {
            return thermo_;
        }

        label condensingSpeciei() const
        {
            return condensingSpeciei_;
        }

        tmp<scalarField> pSat(const scalarField& T) const
        {
            return saturationPressure_->pSat(T);
        }

        tmp<scalarField> condSuctionVelocity(const label patchi) const
        {
            return condSuctionVelocity_[patchi];
        }

        const volScalarField& massTransferRate() const
        {
            return massTransferRate_;
        }

        const volScalarField& qCondensation() const
        {
            return qCondensation_;
        }

        scalarField latentHeat(const scalarField& T) const;

        void predict();


    // Member operators

        void operator=(const cigmaDiffusionLayerCondensation&) = delete;
};

} // End namespace Foam

#endif

// ************************************************************************* //
