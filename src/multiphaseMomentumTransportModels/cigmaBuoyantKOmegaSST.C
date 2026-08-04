/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | cigmaFOAM
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "cigmaBuoyantKOmegaSST.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "fvc.H"
#include "fvm.H"

#include <cmath>

namespace Foam
{
namespace RASModels
{

namespace
{
    inline scalar finiteLimitedScalar
    (
        const scalar value,
        const scalar minValue,
        const scalar maxValue
    )
    {
        if (!std::isfinite(value))
        {
            return 0;
        }

        if (value < minValue)
        {
            return minValue;
        }

        if (value > maxValue)
        {
            return maxValue;
        }

        return value;
    }


    inline void limitScalarList
    (
        scalarField& values,
        const scalar minValue,
        const scalar maxValue
    )
    {
        forAll(values, i)
        {
            values[i] = finiteLimitedScalar(values[i], minValue, maxValue);
        }
    }
}


template<class BasicMomentumTransportModel>
void cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::
validateBuoyancySettings() const
{
    if
    (
        buoyancyModel_ != "SGDH"
     && buoyancyModel_ != "GGDH"
     && buoyancyModel_ != "none"
    )
    {
        FatalErrorInFunction
            << "Invalid buoyancyTurbulence model '" << buoyancyModel_
            << "'. Valid entries are SGDH, GGDH, or none."
            << exit(FatalError);
    }

    if (buoyancyModel_ == "none" && buoyancyDissipation_)
    {
        FatalErrorInFunction
            << "buoyancyTurbulence.dissipation must be false when "
            << "buoyancyTurbulence.model is none."
            << exit(FatalError);
    }

    if (!buoyancyDissipation_ && directionalDissipation_)
    {
        FatalErrorInFunction
            << "buoyancyTurbulence.directionalDissipation requires "
            << "buoyancyTurbulence.dissipation true."
            << exit(FatalError);
    }

    if (sigmaRho_.value() <= SMALL)
    {
        FatalErrorInFunction
            << "buoyancyTurbulence.sigmaRho must be positive."
            << exit(FatalError);
    }

}


template<class BasicMomentumTransportModel>
const uniformDimensionedVectorField&
cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::g() const
{
    if
    (
       !this->mesh_.objectRegistry::template
        foundObject<uniformDimensionedVectorField>("g")
    )
    {
        FatalErrorInFunction
            << "The cigmaBuoyantKOmegaSST model requires a uniform "
            << "dimensionedVectorField named g. Add constant/g or set "
            << "buoyancyTurbulence.model none."
            << exit(FatalError);
    }

    return
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
}


template<class BasicMomentumTransportModel>
void cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::
updateBuoyancyProduction()
{
    if (buoyancyModel_ == "none")
    {
        PkBbyNut_ = dimensionedScalar(PkBbyNut_.dimensions(), 0);
        return;
    }

    const uniformDimensionedVectorField& gField = g();
    if (buoyancyModel_ == "SGDH")
    {
        PkBbyNut_ = -1.0*(gField & fvc::grad(this->rho_))/sigmaRho_;
    }
    else if (buoyancyModel_ == "GGDH")
    {
        Info<< "cigmaBuoyantKOmegaSST: GGDH buoyancy production is active. "
            << "Use with care." << endl;

        tmp<volTensorField> tgradU = fvc::grad(this->U_);

        PkBbyNut_ =
            3.0/(2.0*sigmaRho_*max(this->k_, this->kMin_))
           *(
                gField
              & (
                    (
                        this->nut_*twoSymm(tgradU())
                      - (2.0/3.0)*this->k_*tensor::I
                    )
                  & fvc::grad(this->rho_)
                )
            );
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::PwB
(
    const volScalarField::Internal& gamma
)
{
    if (!buoyancyDissipation_)
    {
        return volScalarField::Internal::New
        (
            "PwB",
            this->mesh_,
            dimensionedScalar(PkBbyNut_.dimensions(), 0)
        );
    }

    C3e_ = dimensionedScalar(dimless, 0);

    const uniformDimensionedVectorField& gField = g();

    if (directionalDissipation_ && mag(gField.value()) > small)
    {
        const vector gHat(gField.value()/mag(gField.value()));

        const volScalarField v(gHat & this->U_);

        const volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("small", dimVelocity, small)
        );

        C3e_ = C3_*tanh(mag(v)/u);
    }
    else
    {
        C3e_ = C3_;
    }

    return
        (gamma + scalar(1))*C3e_()
       *max(PkBbyNut_(), dimensionedScalar(PkBbyNut_().dimensions(), 0))
      - PkBbyNut_();
}


template<class BasicMomentumTransportModel>
cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::cigmaBuoyantKOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    kOmegaSST<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),

    buoyancyDict_(this->RASDict().subOrEmptyDict("buoyancyTurbulence")),
    buoyancyModel_(buoyancyDict_.lookupOrDefault<word>("model", "SGDH")),
    buoyancyDissipation_
    (
        buoyancyDict_.lookupOrDefault<Switch>("dissipation", true)
    ),
    directionalDissipation_
    (
        buoyancyDict_.lookupOrDefault<Switch>
        (
            "directionalDissipation",
            false
        )
    ),
    C3_("C3", dimless, buoyancyDict_.lookupOrDefault<scalar>("C3", 1.0)),
    sigmaRho_
    (
        "sigmaRho",
        dimless,
        buoyancyDict_.lookupOrDefault<scalar>("sigmaRho", 1.0)
    ),
    C3e_
    (
        IOobject
        (
            this->groupName("C3e"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("C3e", dimless, 1.0)
    ),
    PkBbyNut_
    (
        IOobject
        (
            this->groupName("PkBbyNut"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "PkBbyNut",
            dimensionSet(1, -3, -2, 0, 0, 0, 0),
            0
        )
    )
{
    validateBuoyancySettings();

    Info<< "cigmaBuoyantKOmegaSST buoyancy turbulence model: "
        << buoyancyModel_ << nl
        << "    dissipation            " << buoyancyDissipation_ << nl
        << "    directionalDissipation " << directionalDissipation_ << nl
        << "    C3                     " << C3_.value() << nl
        << "    sigmaRho               " << sigmaRho_.value() << endl;
}


template<class BasicMomentumTransportModel>
bool cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST<BasicMomentumTransportModel>::read())
    {
        buoyancyDict_ <<= this->RASDict().subOrEmptyDict("buoyancyTurbulence");

        buoyancyModel_ = buoyancyDict_.lookupOrDefault<word>("model", "SGDH");
        buoyancyDissipation_ =
            buoyancyDict_.lookupOrDefault<Switch>("dissipation", true);
        directionalDissipation_ =
            buoyancyDict_.lookupOrDefault<Switch>
            (
                "directionalDissipation",
                false
            );

        C3_.readIfPresent(buoyancyDict_);
        sigmaRho_.readIfPresent(buoyancyDict_);

        validateBuoyancySettings();

        return true;
    }

    return false;
}


template<class BasicMomentumTransportModel>
void cigmaBuoyantKOmegaSST<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    BaseMomentumTransportModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    limitScalarList(S2.primitiveFieldRef(), 0, vGreat);
    limitScalarList(GbyNu.primitiveFieldRef(), -vGreat, vGreat);

    volScalarField::Internal G
    (
        IOobject
        (
            this->GName(),
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zeroG", nut.dimensions()*GbyNu.dimensions(), 0)
    );

    const scalarField& nutInternal = nut();
    const scalarField& GbyNuInternal = GbyNu;
    scalarField& GInternal = G.primitiveFieldRef();

    forAll(GInternal, celli)
    {
        GInternal[celli] =
            finiteLimitedScalar(nutInternal[celli], 0, vGreat)
           *finiteLimitedScalar(GbyNuInternal[celli], -vGreat, vGreat);
    }

    tgradU.clear();

    this->omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)
       *(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    updateBuoyancyProduction();

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        volScalarField::Internal omegaProduction
        (
            IOobject
            (
                "omegaProduction",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar
            (
                "zeroOmegaProduction",
                rho.dimensions()*GbyNu.dimensions(),
                0
            )
        );

        volScalarField::Internal omegaDivUSp
        (
            IOobject
            (
                "omegaDivUSp",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar
            (
                "zeroOmegaDivUSp",
                rho.dimensions()*divU.dimensions(),
                0
            )
        );

        volScalarField::Internal omegaDissipationSp
        (
            IOobject
            (
                "omegaDissipationSp",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar
            (
                "zeroOmegaDissipationSp",
                rho.dimensions()*this->omega_.dimensions(),
                0
            )
        );

        volScalarField::Internal omegaCrossDiffSp
        (
            IOobject
            (
                "omegaCrossDiffSp",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar
            (
                "zeroOmegaCrossDiffSp",
                rho.dimensions()*CDkOmega.dimensions()/this->omega_.dimensions(),
                0
            )
        );

        const scalarField& rhoI = rho();
        const scalarField& gammaI = gamma;
        const scalarField& betaI = beta;
        const scalarField& divUI = divU;
        const scalarField& omegaI = this->omega_();
        const scalarField& F1I = F1();
        const scalarField& F23I = F23();
        const scalarField& S2I = S2();
        const scalarField& CDkOmegaI = CDkOmega();

        scalarField& omegaProductionI = omegaProduction.primitiveFieldRef();
        scalarField& omegaDivUSpI = omegaDivUSp.primitiveFieldRef();
        scalarField& omegaDissipationSpI =
            omegaDissipationSp.primitiveFieldRef();
        scalarField& omegaCrossDiffSpI =
            omegaCrossDiffSp.primitiveFieldRef();

        forAll(omegaProductionI, celli)
        {
            const scalar alphaCell = 1;
            const scalar rhoCell =
                finiteLimitedScalar(rhoI[celli], small, 1e3);
            const scalar gammaCell =
                finiteLimitedScalar(gammaI[celli], -1e3, 1e3);
            const scalar betaCell =
                finiteLimitedScalar(betaI[celli], 0, 1e3);
            const scalar omegaCell =
                max(finiteLimitedScalar(omegaI[celli], small, 1e12), small);
            const scalar F1Cell =
                finiteLimitedScalar(F1I[celli], 0, 1);
            const scalar F23Cell =
                finiteLimitedScalar(F23I[celli], 0, 1);
            const scalar S2Cell =
                finiteLimitedScalar(S2I[celli], 0, 1e12);
            const scalar CDkOmegaCell =
                finiteLimitedScalar(CDkOmegaI[celli], -1e12, 1e12);
            const scalar GbyNuCell =
                finiteLimitedScalar(GbyNuInternal[celli], -1e12, 1e12);

            const scalar omegaLimiter =
                (this->c1_.value()/this->a1_.value())
               *this->betaStar_.value()*omegaCell
               *max
                (
                    this->a1_.value()*omegaCell,
                    this->b1_.value()*F23Cell*sqrt(S2Cell)
                );

            omegaProductionI[celli] =
                finiteLimitedScalar
                (
                    alphaCell*rhoCell*gammaCell
                   *min(GbyNuCell, omegaLimiter),
                    -1e12,
                    1e12
                );

            omegaDivUSpI[celli] =
                finiteLimitedScalar
                (
                    (2.0/3.0)*alphaCell*rhoCell*gammaCell
                   *finiteLimitedScalar(divUI[celli], -1e12, 1e12),
                    -1e12,
                    1e12
                );

            omegaDissipationSpI[celli] =
                finiteLimitedScalar
                (
                    alphaCell*rhoCell*betaCell*omegaCell,
                    0,
                    1e12
                );

            omegaCrossDiffSpI[celli] =
                finiteLimitedScalar
                (
                    alphaCell*rhoCell*(F1Cell - scalar(1))
                   *CDkOmegaCell/omegaCell,
                    -1e12,
                    1e12
                );
        }

        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            omegaProduction
          - fvm::SuSp(omegaDivUSp, this->omega_)
          + PwB(gamma)
          - fvm::Sp(omegaDissipationSp, this->omega_)
          - fvm::SuSp(omegaCrossDiffSp, this->omega_)
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvModels.source(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(this->omega_);
        this->boundOmega();
    }

    tmp<volScalarField::Internal> tepsilonByk =
        this->epsilonByk(F1, F23);
    const volScalarField::Internal& epsilonByk = tepsilonByk();

    volScalarField::Internal kProduction
    (
        IOobject
        (
            "kProduction",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "zeroKProduction",
            rho.dimensions()*G.dimensions(),
            0
        )
    );

    volScalarField::Internal kDivUSp
    (
        IOobject
        (
            "kDivUSp",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zeroKDivUSp", rho.dimensions()*divU.dimensions(), 0)
    );

    volScalarField::Internal kBuoyancyProduction
    (
        IOobject
        (
            "kBuoyancyProduction",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "zeroKBuoyancyProduction",
            PkBbyNut_.dimensions()*nut.dimensions(),
            0
        )
    );

    volScalarField::Internal kDissipationSp
    (
        IOobject
        (
            "kDissipationSp",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "zeroKDissipationSp",
            rho.dimensions()*epsilonByk.dimensions(),
            0
        )
    );

    const scalarField& kI = this->k_();
    const scalarField& rhoIk = rho();
    const scalarField& omegaIk = this->omega_();
    const scalarField& divUIk = divU;
    const scalarField& PkBbyNutI = PkBbyNut_();
    const scalarField& epsilonBykI = epsilonByk;

    scalarField& kProductionI = kProduction.primitiveFieldRef();
    scalarField& kDivUSpI = kDivUSp.primitiveFieldRef();
    scalarField& kBuoyancyProductionI =
        kBuoyancyProduction.primitiveFieldRef();
    scalarField& kDissipationSpI = kDissipationSp.primitiveFieldRef();

    forAll(kProductionI, celli)
    {
        const scalar rhoCell =
            finiteLimitedScalar(rhoIk[celli], small, 1e3);
        const scalar kCell =
            max(finiteLimitedScalar(kI[celli], small, 1e12), small);
        const scalar omegaCell =
            max(finiteLimitedScalar(omegaIk[celli], small, 1e12), small);
        const scalar GCell =
            finiteLimitedScalar(GInternal[celli], -1e12, 1e12);

        const scalar PkCell =
            min
            (
                GCell,
                this->c1_.value()*this->betaStar_.value()*kCell*omegaCell
            );

        kProductionI[celli] =
            finiteLimitedScalar(rhoCell*PkCell, -1e12, 1e12);

        kDivUSpI[celli] =
            finiteLimitedScalar
            (
                (2.0/3.0)*rhoCell
               *finiteLimitedScalar(divUIk[celli], -1e12, 1e12),
                -1e12,
                1e12
            );

        kBuoyancyProductionI[celli] =
            finiteLimitedScalar
            (
                finiteLimitedScalar(PkBbyNutI[celli], -1e12, 1e12)
               *finiteLimitedScalar(nutInternal[celli], 0, vGreat),
                -1e12,
                1e12
            );

        kDissipationSpI[celli] =
            finiteLimitedScalar
            (
                rhoCell*finiteLimitedScalar(epsilonBykI[celli], 0, 1e12),
                0,
                1e12
            );
    }

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
     ==
        kProduction
      - fvm::SuSp(kDivUSp, this->k_)
      + kBuoyancyProduction
      - fvm::Sp(kDissipationSp, this->k_)
      + this->kSource()
      + fvModels.source(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);
    this->boundOmega();

    this->correctNut(S2, F23);
}

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
