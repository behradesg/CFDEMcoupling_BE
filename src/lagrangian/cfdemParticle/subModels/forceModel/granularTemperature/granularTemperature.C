/*---------------------------------------------------------------------------*\
License
    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2024- Behrad Esgandari, JKU Linz, Austria

Description
    calculates the granular temperature using an adaptive filter

SourceFiles
    granularTemperature.C
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "granularTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(granularTemperature, 0);

addToRunTimeSelectionTable
(
    forceModel,
    granularTemperature,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
granularTemperature::granularTemperature
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    vflucRegName_(typeName + "vfluc_mag"),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    Np_(propsDict_.lookupOrDefault<scalar>("Np",10.0)),
    granularTemperature_
    (   IOobject
        (
            "granularTemperature",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-2,0,0), 0),
        voidfraction_.boundaryField().types()
    ),
    adaptiveFilterLength_
    (
        IOobject
        (
            "adaptiveFilterLength",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0),
        voidfraction_.boundaryField().types()
    ),
    vSmoothGranTempField_
    (
        IOobject
        (
            "vSmoothGranTempField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector::zero),
        UsField_.boundaryField().types()
    ),
    sSmoothGranTempField_
    (
        IOobject
        (
            "sSmoothGranTempField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,//READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0),
        voidfraction_.boundaryField().types()
    ),
    sSmoothGranTempFieldFinal_
    (
        IOobject
        (
            "sSmoothGranTempFieldFinal",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,//READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0),
        voidfraction_.boundaryField().types()
    ),
    DT_
    (   IOobject
        (
            "DT",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("DT", dimensionSet(0,2,-1,0,0,0,0), 0.0),
        voidfraction_.boundaryField().types()
    ),
    deltaX(adaptiveFilterLength_)


{
    particleCloud_.registerParticleProperty<double**>(vflucRegName_,1);
    granularTemperature_.write();


    // init force sub model
    setForceSubModels(propsDict_);
    deltaX.ref() = Foam::cbrt(particleCloud_.mesh().V());
    deltaX.correctBoundaryConditions();     
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

granularTemperature::~granularTemperature()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void granularTemperature::setForce() const
{
    double**& vfluc_ = particleCloud_.getParticlePropertyRef<double**>(vflucRegName_);


    volVectorField vSmoothField = vSmoothGranTempField_;
    volScalarField sSmoothField = sSmoothGranTempField_;

    volScalarField alphaParticle(scalar(1.0) - voidfraction_);

    vSmoothField.dimensions().reset(UsField_.dimensions()*alphaParticle.dimensions());
    vSmoothField.ref() = (UsField_.internalField())*(alphaParticle.internalField());
    vSmoothField.correctBoundaryConditions();
    vSmoothField.oldTime().dimensions().reset(UsField_.dimensions()*alphaParticle.dimensions());
    vSmoothField.oldTime() = UsField_*alphaParticle;
    vSmoothField.oldTime().correctBoundaryConditions();

    sSmoothField.dimensions().reset(alphaParticle.dimensions());
    sSmoothField.ref() = alphaParticle.internalField();
    sSmoothField.correctBoundaryConditions();
    sSmoothField.oldTime().dimensions().reset(alphaParticle.dimensions());
    sSmoothField.oldTime() = alphaParticle;
    sSmoothField.oldTime().correctBoundaryConditions();


    dimensionedScalar deltaT = vSmoothField.mesh().time().deltaT();
    
    dimensionedScalar d3("d3", dimensionSet(0,3,0,0,0,0,0), 0.0);
    
    
    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {  
       label cellI = particleCloud_.cellIDs()[index][0];  
       if (cellI > -1) // particle found
       {
           d3= dimensionedScalar("d3val", dimensionSet(0,3,0,0,0,0,0), Foam::pow(2.0*particleCloud_.radius(index), 3.0));
           break;
       }
    }
    
    adaptiveFilterLength_ = Foam::cbrt(Np_ * d3 / Foam::max(alphaParticle, scalar(0.0001)));

    dimensionedScalar zeroScalar
    (
        "zeroScalar",
        dimensionSet(0,2,0,0,0,0,0),
        0.0
    );
    // adaptive filter length corresponds to the full width at half height of the gaussian kernel.
    // to calculate DT_ refer to Page 34 of Capecelatro, J.S., 2014. A mesoscopic formalism for simulating
    // particle-laden flows  with applications in energy conversion processes. Cornell University.

    if (adaptiveFilterLength_ > deltaX)
    {
    	
    	DT_ = Foam::sqr(adaptiveFilterLength_) / ((16.0*log(2.0)) * deltaT);
    
    } else {
    
       DT_ = zeroScalar / ((16.0*log(2.0)) * deltaT);
    }

    // do the smoothing
    solve
    (
        fvm::ddt(vSmoothField)
       -fvm::laplacian( DT_, vSmoothField)
    );
    solve
    (
        fvm::ddt(sSmoothField)
       -fvm::laplacian( DT_, sSmoothField)
    );

    // get data from working vSmoothField
    volVectorField UsGranTemp(vSmoothField/max(sSmoothField,SMALL));
    UsGranTemp.correctBoundaryConditions();


    label cellI = 0;
    vector velfluc(0,0,0);
    vector position(0,0,0);
    vector UsGranTempInterp(0,0,0);

    interpolationCellPointWallModified<vector> UsGranTempInterpolator_(UsGranTemp);

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];

        if(cellI >= 0)
        {
            position     = particleCloud_.position(index);
            UsGranTempInterp = UsGranTempInterpolator_.interpolate(position,cellI);
            velfluc = particleCloud_.velocity(index) - UsGranTempInterp;
            vfluc_[index][0] = magSqr(velfluc);
        }
    }

    granularTemperature_.primitiveFieldRef() = 0.0;

    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        granularTemperature_,
        vfluc_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );

    granularTemperature_ *= 0.333;

    volScalarField sSmoothFieldFinal =  sSmoothGranTempFieldFinal_;

    sSmoothFieldFinal.dimensions().reset(granularTemperature_.dimensions()*sSmoothField.dimensions());
    sSmoothFieldFinal.ref() = granularTemperature_.internalField()*(sSmoothField.internalField());
    sSmoothFieldFinal.correctBoundaryConditions();
    sSmoothFieldFinal.oldTime().dimensions().reset(granularTemperature_.dimensions()*sSmoothField.dimensions());
    sSmoothFieldFinal.oldTime() = granularTemperature_*sSmoothField;
    sSmoothFieldFinal.oldTime().correctBoundaryConditions();

    solve
    (
        fvm::ddt(sSmoothFieldFinal)
       -fvm::laplacian( DT_, sSmoothFieldFinal)
    );

    granularTemperature_ = sSmoothFieldFinal/max(sSmoothField,SMALL);
    granularTemperature_.correctBoundaryConditions();

    // sSmoothField corresponds to the smoothed alpha.particles field from adaptive filter
    Info<< "Particle-phase averaged granular temperature: "
        << fvc::domainIntegrate(sSmoothField*granularTemperature_).value()
              /fvc::domainIntegrate(sSmoothField).value()
       << endl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
