/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright (C) 2013-     Graz University of
                                                        Technology, IPPT
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

Class
    constDiffSmoothingF

SourceFiles
    constDiffSmoothingF.C

\*---------------------------------------------------------------------------*/

#ifndef constDiffSmoothingF_H
#define constDiffSmoothingF_H

#include "smoothingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constDiffSmoothingF Declaration
\*---------------------------------------------------------------------------*/

class constDiffSmoothingF
:
    public smoothingModel
{

private:

    dictionary propsDict_;

    const scalar lowerLimit_;

    const scalar upperLimit_;

    const scalar filterLength_;

    const scalar smoothingLengthReference_;

    word smoothingLengthFieldName_;

    mutable volScalarField smoothingLengthField_;

    word smoothingLengthReferenceFieldName_;

    mutable volScalarField smoothingLengthReferenceField_;

    mutable volScalarField DT_;

    volScalarField deltaX;

    const bool verbose_;

public:

    //- Runtime type information
    TypeName("constDiffSmoothingF");


    // Constructors

        //- Construct from components
        constDiffSmoothingF
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~constDiffSmoothingF();


    // Member Functions
        bool doSmoothing() const;

        //void dSmoothing(volScalarField&) const;

        void smoothen(volScalarField&) const;

        void smoothen(volVectorField&) const;

        void smoothenReferenceField(volVectorField&, volScalarField&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
