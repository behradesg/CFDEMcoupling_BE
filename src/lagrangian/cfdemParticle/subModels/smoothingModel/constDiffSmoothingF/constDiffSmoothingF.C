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

#include "constDiffSmoothingF.H"
#include "addToRunTimeSelectionTable.H"
#include "interstitialInletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constDiffSmoothingF, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    constDiffSmoothingF,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
constDiffSmoothingF::constDiffSmoothingF
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lowerLimit_(readScalar(propsDict_.lookup("lowerLimit"))),
    upperLimit_(readScalar(propsDict_.lookup("upperLimit"))),
    filterLength_(propsDict_.lookupOrDefault<scalar>("filterLength", -1.0)),
    smoothingLengthReference_(propsDict_.lookupOrDefault<scalar>("smoothingLengthReference",filterLength_)),
    smoothingLengthFieldName_(propsDict_.lookupOrDefault<word>("smoothingLengthFieldName","smoothingLengthField")),
    smoothingLengthField_
    (   IOobject
        (
            smoothingLengthFieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("smoothingLength", dimensionSet(0,1,0,0,0,0,0), filterLength_),
        "zeroGradient"
    ),
    smoothingLengthReferenceFieldName_(propsDict_.lookupOrDefault<word>("smoothingLengthReferenceFieldName","smoothingLengthReferenceField")),
    smoothingLengthReferenceField_
    (   IOobject
        (
            smoothingLengthReferenceFieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("smoothingLengthReference", dimensionSet(0,1,0,0,0,0,0), smoothingLengthReference_),
        "zeroGradient"
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
        "zeroGradient"
    ),
    deltaX(smoothingLengthField_),
    verbose_(propsDict_.found("verbose"))
{
    // either use scalar or field parameters for smoothing
    if (filterLength_ > 0.0 || smoothingLengthReference_ > 0.0)
    {
        if (smoothingLengthField_.headerOk() || smoothingLengthReferenceField_.headerOk())
        {
            FatalError <<"constDiffSmoothingF: Either use scalar or field parameter for smoothing.\n" << abort(FatalError);
        }
    }

    if (filterLength_ < 0.0 && !smoothingLengthField_.headerOk())
    {
        FatalError <<"constDiffSmoothingF: Provide scalar or field parameter for smoothing.\n" << abort(FatalError);
    }
  
    // if no scalar length for smoothing wrt reference field is provided and no
    // such field, use smoothingLengthField
    if (smoothingLengthReference_ < 0.0 && !smoothingLengthReferenceField_.headerOk())
    {
        smoothingLengthReferenceField_ = smoothingLengthField_;
    }
    deltaX.ref() = Foam::cbrt(particleCloud_.mesh().V());
    deltaX.correctBoundaryConditions();
    checkFields(sSmoothField_);
    checkFields(vSmoothField_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constDiffSmoothingF::~constDiffSmoothingF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool constDiffSmoothingF::doSmoothing() const
{
    return true;
}


void constDiffSmoothingF::smoothen(volScalarField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    const fvPatchList& patches = fieldSrc.mesh().boundary(); 
    //Info << "field name is: " << fieldSrc.name() << endl;
    const volScalarField::Boundary fieldSrcPatchBC = fieldSrc.boundaryField();
    
    forAll(patches, iPatch) 
    {
	const fvPatch& curPatch = patches[iPatch];
	label idPatch = fieldSrc.mesh().boundaryMesh().findPatchID(curPatch.name()); 
	//Pout << "fieldSrc.boundaryField().types() is: " << fieldSrc.boundaryField()[idPatch].type() << endl;	
	autoPtr<fvPatchField<scalar>> clonedField(fieldSrcPatchBC[idPatch].clone().ptr());
	sSmoothField_.boundaryFieldRef().set(idPatch, clonedField);
	
	// For BCs without any dictionary values it is possible to use the following code to replace the BCs of smoothField,
	// but if a BC like interstitialInletVelocity is used which requires dictinary values using the following code will
	// initialize the boundary condition with the default values, e.g., alpha will look for volScalarField alpha which is not available!
        // vSmoothField_.boundaryFieldRef().set(idPatch, fvPatchField<scalar>::New(fieldSrc.boundaryField()[idPatch].type(), vSmoothField_.mesh().boundary()[idPatch], vSmoothField_));                                
    }
    
              
    volScalarField sSmoothField = sSmoothField_;

    sSmoothField.dimensions().reset(fieldSrc.dimensions());
    sSmoothField.ref()=fieldSrc.internalField();
    sSmoothField.correctBoundaryConditions();
    sSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    sSmoothField.oldTime()=fieldSrc;
    sSmoothField.oldTime().correctBoundaryConditions();

    dimensionedScalar zeroScalar
    (
        "zeroScalar",
        dimensionSet(0,2,0,0,0,0,0),
        0.0
    );

    dimensionedScalar deltaT = sSmoothField.mesh().time().deltaT();

    // DT_ is calculated based on the approach in Page 34 of Capecelatro, J.S., 2014. A mesoscopic formalism for simulating
    // particle-laden flows  with applications in energy conversion processes. Cornell University.
    if (smoothingLengthField_ > deltaX)
    {
    	
    	DT_ = Foam::sqr(smoothingLengthField_) / ((16.0*log(2.0)) * deltaT);
    
    } else {
    
        DT_ = zeroScalar/ ((16.0*log(2.0)) * deltaT);
    }

    // do smoothing
    solve
    (
        fvm::ddt(sSmoothField)
       -fvm::laplacian(DT_, sSmoothField)
    );

    // bound sSmoothField_
    forAll(sSmoothField,cellI)
    {
        sSmoothField[cellI]=max(lowerLimit_,min(upperLimit_,sSmoothField[cellI]));
    }

    // get data from working sSmoothField - will copy only values at new time
    fieldSrc=sSmoothField;
    fieldSrc.correctBoundaryConditions();

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(sSmoothField.oldTime()) << tab << max(sSmoothField.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void constDiffSmoothingF::smoothen(volVectorField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    const fvPatchList& patches = fieldSrc.mesh().boundary(); 
    //Info << "field name is: " << fieldSrc.name() << endl;
    const volVectorField::Boundary fieldSrcPatchBC = fieldSrc.boundaryField();
    
    forAll(patches, iPatch) 
    {
	const fvPatch& curPatch = patches[iPatch];
	label idPatch = fieldSrc.mesh().boundaryMesh().findPatchID(curPatch.name()); 
	//Pout << "fieldSrc.boundaryField().types() is: " << fieldSrc.boundaryField()[idPatch].type() << endl;	
	autoPtr<fvPatchField<vector>> clonedField(fieldSrcPatchBC[idPatch].clone().ptr());
	vSmoothField_.boundaryFieldRef().set(idPatch, clonedField);
	
	// For BCs without any dictionary values it is possible to use the following code to replace the BCs of smoothField,
	// but if a BC like interstitialInletVelocity is used which requires dictinary values using the following code will
	// initialize the boundary condition with the default values, e.g., alpha will look for volScalarField alpha which is not available!
        // vSmoothField_.boundaryFieldRef().set(idPatch, fvPatchField<vector>::New(fieldSrc.boundaryField()[idPatch].type(), vSmoothField_.mesh().boundary()[idPatch], vSmoothField_));                                
    }
    
    
    volVectorField vSmoothField = vSmoothField_;

    vSmoothField.dimensions().reset(fieldSrc.dimensions());
    vSmoothField.ref()=fieldSrc.internalField();
    vSmoothField.correctBoundaryConditions();
    vSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    vSmoothField.oldTime()=fieldSrc;
    vSmoothField.oldTime().correctBoundaryConditions();

    dimensionedScalar zeroScalar
    (
        "zeroScalar",
        dimensionSet(0,2,0,0,0,0,0),
        0.0
    );

    dimensionedScalar deltaT = vSmoothField.mesh().time().deltaT();

    if (smoothingLengthField_ > deltaX)
    {
    	
    	DT_ = Foam::sqr(smoothingLengthField_) / ((16.0*log(2.0)) * deltaT);
    
    } else {
    
        DT_ = zeroScalar/ ((16.0*log(2.0)) * deltaT);
    }

    // do smoothing
    solve
    (
        fvm::ddt(vSmoothField)
       -fvm::laplacian(DT_, vSmoothField)
    );

    // get data from working vSmoothField
    fieldSrc=vSmoothField;
    fieldSrc.correctBoundaryConditions();

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(vSmoothField.oldTime()) << tab << max(vSmoothField.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void constDiffSmoothingF::smoothenReferenceField(volVectorField& fieldSrc, volScalarField& sFieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template!
    const fvPatchList& patches = fieldSrc.mesh().boundary(); 
    //Info << "field name is: " << fieldSrc.name() << endl;
    const volVectorField::Boundary fieldSrcPatchBC = fieldSrc.boundaryField();
    
    forAll(patches, iPatch) 
    {
	const fvPatch& curPatch = patches[iPatch];
	label idPatch = fieldSrc.mesh().boundaryMesh().findPatchID(curPatch.name()); 
	//Pout << "fieldSrc.boundaryField().types() is: " << fieldSrc.boundaryField()[idPatch].type() << endl;	
	autoPtr<fvPatchField<vector>> clonedField(fieldSrcPatchBC[idPatch].clone().ptr());
	vSmoothField_.boundaryFieldRef().set(idPatch, clonedField);
	
	// For BCs without any dictionary values it is possible to use the following code to replace the BCs of smoothField,
	// but if a BC like interstitialInletVelocity is used which requires dictinary values using the following code will
	// initialize the boundary condition with the default values, e.g., alpha will look for volScalarField alpha which is not available!
        // vSmoothField_.boundaryFieldRef().set(idPatch, fvPatchField<vector>::New(fieldSrc.boundaryField()[idPatch].type(), vSmoothField_.mesh().boundary()[idPatch], vSmoothField_));                                
    }	
    // to print-out and check BCs of it is possible to use the following code,
    // const volVectorField::Boundary& vSmoothFieldBC = vSmoothField_.boundaryField();
    // Pout << "vSmoothFieldBC After change are: " << vSmoothFieldBC<< endl;
    
           
    volVectorField vSmoothField = vSmoothField_;
    volScalarField sSmoothField = sSmoothField_;

    vSmoothField.dimensions().reset(fieldSrc.dimensions()*sFieldSrc.dimensions());
    vSmoothField.ref() = (fieldSrc.internalField())*(sFieldSrc.internalField());
    vSmoothField.correctBoundaryConditions();
    vSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions()*sFieldSrc.dimensions());
    vSmoothField.oldTime() = fieldSrc*sFieldSrc;
    vSmoothField.oldTime().correctBoundaryConditions();

    sSmoothField.dimensions().reset(sFieldSrc.dimensions());
    sSmoothField.ref() = sFieldSrc.internalField();
    sSmoothField.correctBoundaryConditions();
    sSmoothField.oldTime().dimensions().reset(sFieldSrc.dimensions());
    sSmoothField.oldTime() = sFieldSrc;
    sSmoothField.oldTime().correctBoundaryConditions();

    dimensionedScalar zeroScalar
    (
        "zeroScalar",
        dimensionSet(0,2,0,0,0,0,0),
        0.0
    );

    dimensionedScalar deltaT = vSmoothField.mesh().time().deltaT();

    if (smoothingLengthField_ > deltaX)
    {
    	
    	DT_ = Foam::sqr(smoothingLengthField_) / ((16.0*log(2.0)) * deltaT);
    
    } else {
    
        DT_ = zeroScalar/ ((16.0*log(2.0)) * deltaT);
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
    fieldSrc = vSmoothField/max(sSmoothField, SMALL);
    fieldSrc.correctBoundaryConditions();

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(vSmoothField.oldTime()) << tab << max(vSmoothField.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
