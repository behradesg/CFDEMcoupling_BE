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

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

Description
    calculates the granular temperature with a fixed filter with the size of mesh
    cells

SourceFiles
    cellBasedGranularTemperature.C
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "cellBasedGranularTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cellBasedGranularTemperature, 0);

addToRunTimeSelectionTable
(
    forceModel,
    cellBasedGranularTemperature,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cellBasedGranularTemperature::cellBasedGranularTemperature
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    vflucRegName_(typeName + "vfluc_mag"),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    cellBasedGranularTemperature_
    (   IOobject
        (
            "cellBasedGranularTemperature",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-2,0,0), 0),
        "zeroGradient"
    )
{
    particleCloud_.registerParticleProperty<double**>(vflucRegName_,1);
    cellBasedGranularTemperature_.write();


    // init force sub model
    setForceSubModels(propsDict_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellBasedGranularTemperature::~cellBasedGranularTemperature()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void cellBasedGranularTemperature::setForce() const
{
    double**& vfluc_ = particleCloud_.getParticlePropertyRef<double**>(vflucRegName_);

    label cellI = 0;
    vector velfluc(0,0,0);
    vector position(0,0,0);
    vector UsGranTempInterp(0,0,0);

    interpolationCellPointWallModified<vector> UsGranTempInterpolator_(UsField_);


    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            position = particleCloud_.position(index);
            UsGranTempInterp = UsGranTempInterpolator_.interpolate(position,cellI);
            velfluc = particleCloud_.velocity(index) - UsGranTempInterp;
            vfluc_[index][0] = magSqr(velfluc);
        }
    }

    cellBasedGranularTemperature_.primitiveFieldRef() = 0.0;

    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        cellBasedGranularTemperature_,
        vfluc_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );

    cellBasedGranularTemperature_ *= 0.333;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
