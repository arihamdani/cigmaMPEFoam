/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 cigmaMPEFoam contributors
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    GNU General Public License, version 3 or later.

    Adapted from containmentFOAM's diffusionLayerWallCondensation model,
    Copyright (C) Forschungszentrum Juelich.
\*---------------------------------------------------------------------------*/

#include "cigmaDiffusionLayerCondensation.H"
#include "cigmaSaturatedSteamFvPatchScalarField.H"
#include "wallFvPatch.H"

namespace Foam
{
    defineTypeNameAndDebug(cigmaDiffusionLayerCondensation, 0);
}

const Foam::word Foam::cigmaDiffusionLayerCondensation::dictionaryName
(
    "condensationProperties"
);


Foam::IOobject Foam::cigmaDiffusionLayerCondensation::createIOobject
(
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        dictionaryName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
    }

    return io;
}


Foam::autoPtr<Foam::saturationPressureModel>
Foam::cigmaDiffusionLayerCondensation::createSaturationModel() const
{
    if (found("saturationPressure"))
    {
        return saturationPressureModel::New("saturationPressure", *this);
    }

    dictionary defaultModel;
    dictionary modelCoeffs;
    modelCoeffs.add("type", "ArdenBuck");
    defaultModel.add("saturationPressure", modelCoeffs);

    Info<< "    saturationPressure not specified; using ArdenBuck" << endl;

    return saturationPressureModel::New
    (
        "saturationPressure",
        defaultModel
    );
}


void Foam::cigmaDiffusionLayerCondensation::validate() const
{
    if (!active_)
    {
        return;
    }

    if
    (
        lookupOrDefault<word>("wallCondensationModel", "diffusionLayer")
     != "diffusionLayer"
    )
    {
        FatalIOErrorInFunction(*this)
            << "cigmaSinglePhaseFluid currently supports only the "
            << "diffusionLayer wallCondensationModel"
            << exit(FatalIOError);
    }

    if (!thermo_.containsSpecie(condensingSpecie_))
    {
        FatalIOErrorInFunction(*this)
            << "Condensing specie " << condensingSpecie_
            << " is not present in " << physicalProperties::typeName
            << exit(FatalIOError);
    }

    if (thermo_.species().size() < 2)
    {
        FatalIOErrorInFunction(*this)
            << "Wall condensation requires H2O and at least one "
            << "non-condensable specie"
            << exit(FatalIOError);
    }

    label nCondensingPatches = 0;
    const volScalarField& Yc = thermo_.Y()[condensingSpeciei_];

    forAll(Yc.boundaryField(), patchi)
    {
        if
        (
            Yc.boundaryField()[patchi].type()
         == compressible::cigmaSaturatedSteamFvPatchScalarField::typeName
        )
        {
            ++nCondensingPatches;
        }
    }

    if (!nCondensingPatches)
    {
        FatalIOErrorInFunction(*this)
            << "No saturatedSteam boundary condition was found for field "
            << Yc.name() << exit(FatalIOError);
    }
}


Foam::cigmaDiffusionLayerCondensation::cigmaDiffusionLayerCondensation
(
    const fluidMulticomponentThermo& thermo,
    const fluidMulticomponentThermophysicalTransportModel&
        thermophysicalTransport,
    const volScalarField& rho
)
:
    IOdictionary(createIOobject(thermo.mesh())),
    thermo_(thermo),
    thermophysicalTransport_(thermophysicalTransport),
    rho_(rho),
    active_
    (
        readOpt() == IOobject::NO_READ
      ? Switch(false)
      : lookupOrDefault<Switch>("solveWallCondensation", false)
    ),
    condensingSpecie_
    (
        readOpt() == IOobject::NO_READ
      ? word("H2O")
      : lookupOrDefault<word>("condensingSpecie", "H2O")
    ),
    condensingSpeciei_
    (
        thermo.containsSpecie(condensingSpecie_)
      ? thermo.species()[condensingSpecie_]
      : -1
    ),
    saturationPressure_(nullptr),
    maxCondWallMassFlux_
    (
        readOpt() == IOobject::NO_READ
      ? 100.0
      : lookupOrDefault<scalar>("maxCondWallMassFlux", 100.0)
    ),
    latentHeatA_(3557560.510766),
    latentHeatB_(-6509.752221),
    latentHeatC_(13.940664),
    latentHeatD_(-0.015642),
    massTransferRate_
    (
        IOobject
        (
            "massTransferRate",
            thermo.mesh().time().name(),
            thermo.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.mesh(),
        dimensionedScalar
        (
            "massTransferRate",
            dimMass/dimArea/dimTime,
            0
        )
    ),
    qCondensation_
    (
        IOobject
        (
            "qcond",
            thermo.mesh().time().name(),
            thermo.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.mesh(),
        dimensionedScalar("qcond", dimPower/dimArea, 0)
    ),
    condSuctionVelocity_
    (
        thermo.mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    )
{
    if (!active_)
    {
        Info<< "Single-phase wall condensation is OFF" << endl;
        return;
    }

    saturationPressure_ = createSaturationModel();

    if (isDict("latentHeatModel"))
    {
        const dictionary& latentHeatDict = subDict("latentHeatModel");

        if (latentHeatDict.lookupOrDefault<word>("type", "polynomial")
         != "polynomial")
        {
            FatalIOErrorInFunction(latentHeatDict)
                << "Only the polynomial latentHeatModel is currently "
                << "supported" << exit(FatalIOError);
        }

        latentHeatA_ = latentHeatDict.lookupOrDefault<scalar>
        (
            "A",
            latentHeatA_
        );
        latentHeatB_ = latentHeatDict.lookupOrDefault<scalar>
        (
            "B",
            latentHeatB_
        );
        latentHeatC_ = latentHeatDict.lookupOrDefault<scalar>
        (
            "C",
            latentHeatC_
        );
        latentHeatD_ = latentHeatDict.lookupOrDefault<scalar>
        (
            "D",
            latentHeatD_
        );
    }

    validate();

    Info<< "Single-phase wall condensation is ON" << nl
        << "    wallCondensationModel = diffusionLayer" << nl
        << "    condensingSpecie = " << condensingSpecie_ << nl
        << "    maxCondWallMassFlux = " << maxCondWallMassFlux_ << endl;
}


Foam::cigmaDiffusionLayerCondensation::~cigmaDiffusionLayerCondensation()
{}


Foam::scalarField Foam::cigmaDiffusionLayerCondensation::latentHeat
(
    const scalarField& T
) const
{
    return
        latentHeatA_
      + latentHeatB_*T
      + latentHeatC_*sqr(T)
      + latentHeatD_*pow(T, 3);
}


void Foam::cigmaDiffusionLayerCondensation::predict()
{
    if (!active_)
    {
        return;
    }

    Info<< "Updating single-phase diffusion-layer wall condensation" << endl;

    massTransferRate_.boundaryFieldRef() = Zero;
    qCondensation_.boundaryFieldRef() = Zero;
    condSuctionVelocity_ = Zero;

    volScalarField& Yc =
        const_cast<volScalarField&>(thermo_.Y()[condensingSpeciei_]);

    Yc.correctBoundaryConditions();

    const volScalarField::Boundary& rhoBf = rho_.boundaryField();
    const volScalarField::Boundary& TBf = thermo_.T().boundaryField();

    forAll(Yc.boundaryField(), patchi)
    {
        if
        (
            Yc.boundaryField()[patchi].type()
         != compressible::cigmaSaturatedSteamFvPatchScalarField::typeName
        )
        {
            continue;
        }

        const scalarField jH2O
        (
            thermophysicalTransport_.j(Yc, patchi)
        );

        const scalarField& Yw = Yc.boundaryField()[patchi];
        scalarField& mDot = massTransferRate_.boundaryFieldRef()[patchi];
        scalarField& qCond = qCondensation_.boundaryFieldRef()[patchi];
        scalarField& suction = condSuctionVelocity_[patchi];

        mDot =
            -min
            (
                posPart(jH2O)/max(1 - Yw, scalar(SMALL)),
                scalarField(jH2O.size(), maxCondWallMassFlux_)
            );

        suction = -mDot/max(rhoBf[patchi], scalar(SMALL));
        qCond = -mDot*latentHeat(TBf[patchi]);
    }
}

// ************************************************************************* //
