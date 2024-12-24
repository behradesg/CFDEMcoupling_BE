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

    ref. Capecelatro, J. et al. "Numerical study of collisional particle dynamics
    in cluster-induced turbulence." Journal of Fluid Mechanics 747 (2014).

SourceFiles
    totalParticleFluctuatingEnergy.C
\*---------------------------------------------------------------------------*/

#ifndef totalParticleFluctuatingEnergy_H
#define totalParticleFluctuatingEnergy_H

#include "forceModel.H"
#include "averagingModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class totalParticleFluctuatingEnergy Declaration
\*---------------------------------------------------------------------------*/

class totalParticleFluctuatingEnergy
:
    public forceModel
{
private:

public:

    //- Runtime type information
    TypeName("totalParticleFluctuatingEnergy");

    // Constructors

        //- Construct from components
        totalParticleFluctuatingEnergy
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~totalParticleFluctuatingEnergy();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
