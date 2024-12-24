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

#include "simpleFilterSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleFilterSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    simpleFilterSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
simpleFilterSmoothing::simpleFilterSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lowerLimit_(readScalar(propsDict_.lookup("lowerLimit"))),
    upperLimit_(readScalar(propsDict_.lookup("upperLimit"))),
    iter_(propsDict_.lookupOrDefault<label>("iterations",1)),
    verbose_(propsDict_.found("verbose"))
{
    checkFields(sSmoothField_);
    checkFields(vSmoothField_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simpleFilterSmoothing::~simpleFilterSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool simpleFilterSmoothing::doSmoothing() const
{
    return true;
}


void simpleFilterSmoothing::smoothen(volScalarField& fieldSrc) const
{
    simpleFilter filter(fieldSrc.mesh());
    
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
             << " smoothing iteration(s): " 
             << min(fieldSrc).value()
             << ", " 
             << max(fieldSrc).value() 
             << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void simpleFilterSmoothing::smoothen(volVectorField& fieldSrc) const
{
    simpleFilter filter(fieldSrc.mesh());

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
             << " smoothing iteration(s): " 
             << min(mag(fieldSrc)).value() 
             << ", " 
             << max(mag(fieldSrc)).value() 
             << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void simpleFilterSmoothing::smoothenReferenceField(volVectorField& fieldSrc, volScalarField& sFieldSrc) const
{
    simpleFilter filter(fieldSrc.mesh());

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
             << ", ref) after "
             << iter_ 
             << " smoothing iteration(s): " 
             << min(mag(fieldSrc)).value() 
             << ", " 
             << max(mag(fieldSrc)).value() 
             << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
