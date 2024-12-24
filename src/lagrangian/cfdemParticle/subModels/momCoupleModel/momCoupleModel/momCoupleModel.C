/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "momCoupleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(momCoupleModel, 0);

defineRunTimeSelectionTable(momCoupleModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> momCoupleModel::impMomSource()
{
    // FatalError<<"the solver calls for impMomSource()\n"
    //           <<", please set 'momCoupleModel' to type 'implicitCouple'\n"
    //           << abort(FatalError);

    return tmp<volScalarField>
    (
        new volScalarField 
        (
            IOobject
            (
                "tsource",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0),scalar(0.0))
        )
    );
}

tmp<volVectorField> momCoupleModel::expMomSource()
{
    // FatalError<<"the solver calls for expMomSource()\n"
    //           <<", please set 'momCoupleModel' to type 'explicitCouple'\n"
    //           << abort(FatalError);
    return tmp<volVectorField>
    (
        new volVectorField 
        (
            IOobject
            (
                "tsource",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedVector("zero",dimensionSet(1,-2,-2,0,0),vector::zero)
        )
    );
}

tmp<volScalarField> momCoupleModel::CdS()
{
    // FatalError<<"the solver calls for expMomSource()\n"
    //           <<", please set 'momCoupleModel' to type 'explicitCouple'\n"
    //           << abort(FatalError);
    return tmp<volScalarField>
    (
        new volScalarField 
        (
            IOobject
            (
                "tsource",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0),scalar(0.0))
        )
    );
}

void momCoupleModel::setSourceField(volVectorField & a) const
{
    //do nothing;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
momCoupleModel::momCoupleModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    maxAlpha_(1-SMALL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

momCoupleModel::~momCoupleModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
