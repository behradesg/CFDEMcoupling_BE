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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "laplaceFilterSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laplaceFilterSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    laplaceFilterSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
laplaceFilterSmoothing::laplaceFilterSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lowerLimit_(readScalar(propsDict_.lookup("lowerLimit"))),
    upperLimit_(readScalar(propsDict_.lookup("upperLimit"))),
    widthCoeff_(propsDict_.lookupOrDefault<scalar>("widthCoeff",6.0)),
    iter_(propsDict_.lookupOrDefault<label>("iterations",1)),
    verbose_(propsDict_.found("verbose"))
{
    checkFields(sSmoothField_);
    checkFields(vSmoothField_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

laplaceFilterSmoothing::~laplaceFilterSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool laplaceFilterSmoothing::doSmoothing() const
{
    return true;
}


void laplaceFilterSmoothing::smoothen(volScalarField& fieldSrc) const
{
    laplaceFilter filter(fieldSrc.mesh(),widthCoeff_);
    
    if(verbose_)
    {
        Info << "min/max(" 
             << fieldSrc.name() 
             << ") before smoothing " 
             << min(fieldSrc).value()
             << ", " 
             << max(fieldSrc).value() 
             << endl;
    }

    for (int i=0; i < iter_; i++) {
        fieldSrc = filter(fieldSrc);
        fieldSrc.correctBoundaryConditions();
    }
    fieldSrc.min(upperLimit_);
    fieldSrc.max(lowerLimit_);

    if(verbose_)
    {
        Info << "min/max(" 
             << fieldSrc.name() 
             << ") after " 
             << iter_ 
             << " smoothing iterations: " 
             << min(fieldSrc).value()
             << ", " 
             << max(fieldSrc).value() 
             << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void laplaceFilterSmoothing::smoothen(volVectorField& fieldSrc) const
{
    laplaceFilter filter(fieldSrc.mesh(),widthCoeff_);

    for (int i=0; i < iter_; i++) {
        fieldSrc = filter(fieldSrc);
        fieldSrc.correctBoundaryConditions();
    }

    if(verbose_)
    {
        Info << "min/max(" 
             << fieldSrc.name() 
             << ") after "
             << iter_ 
             << " smoothing iterations: " 
             << min(mag(fieldSrc)).value() 
             << ", " 
             << max(mag(fieldSrc)).value() 
             << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void laplaceFilterSmoothing::smoothenReferenceField(volVectorField& fieldSrc, volScalarField& sFieldSrc) const
{
    laplaceFilter filter(fieldSrc.mesh(),widthCoeff_);

    volScalarField sSmoothField(sFieldSrc);
    volVectorField vSmoothField(fieldSrc);

    for (int i=0; i < iter_; i++) {
        vSmoothField = filter(sSmoothField*vSmoothField);
        sSmoothField = filter(sSmoothField);
        sSmoothField.min(1.0);
        sSmoothField.max(0.0);
        sSmoothField.correctBoundaryConditions();
        vSmoothField /= (sSmoothField + SMALL);
        vSmoothField.correctBoundaryConditions();
    }  
    // ensure that velocity is 0 in empty cells
    // fieldSrc = pos(mag(fieldSrc))*vSmoothField;
    
    // get data from working vSmoothField
    fieldSrc = vSmoothField;
    fieldSrc.correctBoundaryConditions();

    if(verbose_)
    {
        Info << "min/max(" 
             << fieldSrc.name() 
             << ") after "
             << iter_ 
             << " smoothing iterations: " 
             << min(mag(fieldSrc)).value() 
             << ", " 
             << max(mag(fieldSrc)).value() 
             << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
