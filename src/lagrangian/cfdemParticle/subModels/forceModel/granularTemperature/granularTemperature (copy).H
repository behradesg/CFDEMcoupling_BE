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
    calculates the granular granular temperature using an adaptive filter
SourceFiles
    granularTemperature.C
\*---------------------------------------------------------------------------*/

#ifndef granularTemperature_H
#define granularTemperature_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint2.H"
#include "interpolationCellPointWallModified.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class granularTemperature Declaration
\*---------------------------------------------------------------------------*/

class granularTemperature
:
    public forceModel
{
private:

    dictionary propsDict_;

    const word vflucRegName_;

    const word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    const word UsFieldName_;

    const volVectorField& UsField_;

    const scalar Np_;

    mutable volScalarField granularTemperature_;

    mutable volScalarField adaptiveFilterLength_;

    mutable volVectorField vSmoothGranTempField_;

    mutable volScalarField sSmoothGranTempField_;

    mutable volScalarField sSmoothGranTempFieldFinal_;

    mutable volScalarField DT_;

    volScalarField deltaX;
    

public:

    //- Runtime type information
    TypeName("granularTemperature");

    // Constructors

        //- Construct from components
        granularTemperature
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~granularTemperature();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
