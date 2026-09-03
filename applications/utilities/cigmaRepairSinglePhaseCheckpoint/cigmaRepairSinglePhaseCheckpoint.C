/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | cigmaMPEFoam
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"

using namespace Foam;


int main(int argc, char* argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    timeSelector::select0(runTime, args);
    #include "createRegionMeshNoChangers.H"

    const word timeName(runTime.name());

    typeIOobject<volScalarField> existingH2O
    (
        "H2O",
        timeName,
        mesh,
        IOobject::NO_READ
    );

    if (existingH2O.headerOk())
    {
        FatalErrorInFunction
            << "H2O already exists at time " << timeName
            << ". Refusing to overwrite it." << exit(FatalError);
    }

    volScalarField H2O0
    (
        IOobject
        (
            "H2O_0",
            timeName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField H2O
    (
        IOobject
        (
            "H2O",
            "0",
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField HE
    (
        IOobject
        (
            "HE",
            timeName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField AIR
    (
        IOobject
        (
            "AIR",
            timeName,
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    H2O.primitiveFieldRef() = H2O0.primitiveField();
    forAll(H2O.boundaryField(), patchi)
    {
        if (H2O.boundaryField()[patchi].type() == "saturatedSteam")
        {
            // The thermo validates the species sum before saturatedSteam can
            // update.  Match the value reconstructed by its read constructor.
            H2O.boundaryFieldRef()[patchi] ==
                H2O.boundaryField()[patchi].patchInternalField();
        }
        else
        {
            H2O.boundaryFieldRef()[patchi] == H2O0.boundaryField()[patchi];
        }
    }
    H2O.instance() = timeName;

    // Reconstruct the wall values written as fixedGradient.  These values are
    // recomputed when the solver reads the checkpoint and must therefore be
    // included before rebuilding the default specie.
    forAll(HE.boundaryField(), patchi)
    {
        if (HE.boundaryField()[patchi].type() == "fixedGradient")
        {
            HE.boundaryFieldRef()[patchi].evaluate();
        }
    }

    const volScalarField activeSum(H2O + HE);
    const scalar minActive = min(activeSum).value();
    const scalar maxActive = max(activeSum).value();

    if (minActive < -small || maxActive > 1 + small)
    {
        FatalErrorInFunction
            << "H2O_0 + HE is outside [0, 1]: min=" << minActive
            << ", max=" << maxActive << exit(FatalError);
    }

    volScalarField AIRBackup
    (
        IOobject
        (
            "AIR_beforeCheckpointRepair",
            timeName,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        AIR
    );

    AIR == scalar(1) - activeSum;
    AIR.max(0);

    const volScalarField repairedSum(H2O + HE + AIR);

    Info<< "Repaired specie sum: min=" << min(repairedSum).value()
        << ", max=" << max(repairedSum).value() << nl << endl;

    AIRBackup.write();
    H2O.write();
    AIR.write();

    Info<< "Wrote H2O and AIR at time " << timeName << nl
        << "Preserved AIR as AIR_beforeCheckpointRepair" << nl
        << "End" << endl;

    return 0;
}

// ************************************************************************* //
