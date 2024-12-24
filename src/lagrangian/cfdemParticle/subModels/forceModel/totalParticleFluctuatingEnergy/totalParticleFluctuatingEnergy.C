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
    calculates the total particle fluctuating energy

SourceFiles
    totalParticleFluctuatingEnergy.C
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "totalParticleFluctuatingEnergy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(totalParticleFluctuatingEnergy, 0);

addToRunTimeSelectionTable
(
    forceModel,
    totalParticleFluctuatingEnergy,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
totalParticleFluctuatingEnergy::totalParticleFluctuatingEnergy
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm)
{


}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

totalParticleFluctuatingEnergy::~totalParticleFluctuatingEnergy()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void totalParticleFluctuatingEnergy::setForce() const
{
    vector velfluc(0,0,0);
    scalar Vs       = scalar(0.0);
    scalar mp       = scalar(0.0);
    scalar rhos     = scalar(0.0);
    vector Up       = vector::zero;
    vector mpUp_proc    = vector::zero;
    scalar totalPartMass_proc = scalar(0.0);

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        Vs       = particleCloud_.particleVolume(index);
        rhos     = particleCloud_.particleDensity(index);
        mp       = Vs*rhos;
        Up       = particleCloud_.velocity(index);
        mpUp_proc     += mp*Up;
        totalPartMass_proc += mp;
    }

    vector mpUp          = returnReduce(mpUp_proc, sumOp<vector>());
    scalar totalPartMass = returnReduce(totalPartMass_proc, sumOp<scalar>());
    vector domainAveragedUp = mpUp/totalPartMass;

    scalar dotProductVelfluc = scalar(0.0);
    scalar sum_mpDotProductVelfluc_proc = scalar(0.0);
    scalar mpDotProductVelfluc_proc = scalar(0.0);

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        velfluc = particleCloud_.velocity(index) - domainAveragedUp;
        dotProductVelfluc = magSqr(velfluc);
        mpDotProductVelfluc_proc = mp*dotProductVelfluc;
        sum_mpDotProductVelfluc_proc += mpDotProductVelfluc_proc;
    }

    scalar totalMpDotProductVelfluc = returnReduce(sum_mpDotProductVelfluc_proc, sumOp<scalar>());
    Info << "Total particle fluctuating energy: " << 0.5 * totalMpDotProductVelfluc / totalPartMass << endl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
