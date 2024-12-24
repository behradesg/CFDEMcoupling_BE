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

#include "denseFace.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"
#include "syncTools.H"

//#include <mpi.h>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(denseFace, 0);

addToRunTimeSelectionTable
(
    averagingModel,
    denseFace,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
denseFace::denseFace
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    averagingModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

denseFace::~denseFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void denseFace::setScalarAverage
(
    volScalarField& field,
    double**& value,
    double**const& weight,
    volScalarField& weightField,
    double**const& mask
) const
{
    // get mesh
    const fvMesh& mesh(particleCloud_.mesh());

    // allocate face based averaging array
    List<scalar> averageField;
    List<scalar> particleWeights;
    List<scalar> faceWeights;
    averageField.resize(mesh.nFaces());
    particleWeights.resize(mesh.nFaces());
    faceWeights.resize(mesh.nFaces());

    forAll(averageField,i)
    {
        averageField[i]    = 0.;
        faceWeights[i]     = 0.5;
        particleWeights[i] = 0.;
    }
        // set correct faceWeights for facesums
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    forAll(boundaryMesh,patchi)
    {
        const label start = boundaryMesh[patchi].start();
        forAll(boundaryMesh[patchi],facei)
        {
            if (!boundaryMesh[patchi].coupled()) faceWeights[facei+start] = 1.0;
        }
    }

    scalar valueScal;
    scalar weightP;

    for(int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type
        // start at index 1, since 0 contains the particle cellID
        for(int faceI = 1; faceI < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; faceI++)
        {
            //Info << "subCell=" << subCell << endl;
            label globalFaceI = particleCloud_.cellIDs()[index][faceI];
            
            if (globalFaceI >= 0)
            {
                valueScal = value[index][0];
                weightP = weight[index][faceI];
                averageField[globalFaceI]    += valueScal*weightP;
                particleWeights[globalFaceI] += weightP;
            }
        }
    }
    // sync faceVoidFraction accross processor boundaries
    syncTools::syncFaceList
    (
        mesh,
        averageField,
        plusEqOp<scalar>()
    );
    syncTools::syncFaceList
    (
        mesh,
        particleWeights,
        plusEqOp<scalar>()
    );

    // reconstructing face based averaged field
    forAll(field,cellI)
    {
        scalar localAveragedValue(0.);
        scalar localWeightValue(SMALL);
        forAll(mesh.cells()[cellI], faceI)
        {
            label globalFaceI   = mesh.cells()[cellI][faceI];
            localAveragedValue += averageField[globalFaceI]*faceWeights[globalFaceI];
            localWeightValue   += particleWeights[globalFaceI]*faceWeights[globalFaceI];
        }   
        // Pout << "area: " << totalArea << " vof: " << localAlphaP << endl;
        field[cellI] = localAveragedValue/localWeightValue;
        weightField[cellI] = localWeightValue;
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void denseFace::setVectorAverage
(
    volVectorField& field,
    double**& value,
    double**const& weight,
    volScalarField& weightField,
    double**const& mask
) const
{
    // get mesh
    const fvMesh& mesh(particleCloud_.mesh());

    // allocate face based averaging array
    List<scalar> averageField1;
    List<scalar> averageField2;
    List<scalar> averageField3;
    List<scalar> particleWeights;
    List<scalar> faceWeights;
    averageField1.resize(mesh.nFaces());
    averageField2.resize(mesh.nFaces());
    averageField3.resize(mesh.nFaces());
    particleWeights.resize(mesh.nFaces());  
    faceWeights.resize(mesh.nFaces());

    forAll(averageField1,i)
    {
        averageField1[i]   = 0.;
        averageField2[i]   = 0.;
        averageField3[i]   = 0.;
        particleWeights[i] = 0.;
        faceWeights[i]     = 0.5;
    }
        // set correct faceWeights for facesums
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    forAll(boundaryMesh,patchi)
    {
        const label start = boundaryMesh[patchi].start();
        forAll(boundaryMesh[patchi],facei)
        {
            if (!boundaryMesh[patchi].coupled()) faceWeights[facei+start] = 1.0;
        }
    }

    vector valueVec;
    scalar weightP;

    for(int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type
        // start at index 1, since 0 contains the particle cellID
        for(int faceI = 1; faceI < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; faceI++)
        {
            //Info << "subCell=" << subCell << endl;
            label globalFaceI = particleCloud_.cellIDs()[index][faceI];
            label cellI = particleCloud_.cellIDs()[index][0];
            
            weightField[cellI] = 0.0;
            if (globalFaceI >= 0)
            {
                for(int i = 0; i < 3; i++)
                    valueVec[i] = value[index][i];
                weightP = weight[index][faceI];
                
                weightField[cellI] += weightP;
                averageField1[globalFaceI]   += valueVec[0]*weightP;
                averageField2[globalFaceI]   += valueVec[1]*weightP;
                averageField3[globalFaceI]   += valueVec[2]*weightP;
                particleWeights[globalFaceI] += weightP;
            }
        }
    }
    // sync faceVoidFraction accross processor boundaries
    syncTools::syncFaceList
    (
        mesh,
        averageField1,
        plusEqOp<scalar>()
    );
    syncTools::syncFaceList
    (
        mesh,
        averageField2,
        plusEqOp<scalar>()
    );
    syncTools::syncFaceList
    (
        mesh,
        averageField3,
        plusEqOp<scalar>()
    );
    syncTools::syncFaceList
    (
        mesh,
        particleWeights,
        plusEqOp<scalar>()
    );

    // reconstructing face based averaged field
    forAll(field,cellI)
    {
        vector localAveragedValue(vector::zero);
        scalar localWeightValue(SMALL);
        forAll(mesh.cells()[cellI], faceI)
        {
            label globalFaceI = mesh.cells()[cellI][faceI];
            localAveragedValue[0] += averageField1[globalFaceI]*faceWeights[globalFaceI];
            localAveragedValue[1] += averageField2[globalFaceI]*faceWeights[globalFaceI];
            localAveragedValue[2] += averageField3[globalFaceI]*faceWeights[globalFaceI];
            localWeightValue      += particleWeights[globalFaceI]*faceWeights[globalFaceI];
        }   
        // Pout << "area: " << totalArea << " vof: " << localAlphaP << endl;
        field[cellI] = localAveragedValue/localWeightValue;
        weightField[cellI] = localWeightValue;
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void denseFace::setScalarSum
(
    volScalarField& field,
    double**& value,
    double**const& weight,
    double**const& mask
) const
{
    // get mesh
    const fvMesh& mesh(particleCloud_.mesh());

    // allocate face based averaging array
    List<scalar> sumField;
    List<scalar> faceWeights;
    sumField.resize(mesh.nFaces());
    faceWeights.resize(mesh.nFaces());

    forAll(sumField,i)
    {
        sumField[i]    = 0.;
        faceWeights[i] = 0.5;
    }
        // set correct faceWeights for facesums
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    forAll(boundaryMesh,patchi)
    {
        const label start = boundaryMesh[patchi].start();
        forAll(boundaryMesh[patchi],facei)
        {
            if (!boundaryMesh[patchi].coupled()) faceWeights[facei+start] = 1.0;
        }
    }

    scalar valueScal;
    scalar weightP;

    for(int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type
        // start at index 1, since 0 contains the particle cellID
        for(int faceI = 1; faceI < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; faceI++)
        {
            //Info << "subCell=" << subCell << endl;
            label globalFaceI = particleCloud_.cellIDs()[index][faceI];
            
            if (globalFaceI >= 0)
            {
                valueScal = value[index][0];
                weightP = weight[index][faceI];
                sumField[globalFaceI] += valueScal*weightP;
            }
        }
    }
    // sync faceVoidFraction accross processor boundaries
    syncTools::syncFaceList
    (
        mesh,
        sumField,
        plusEqOp<scalar>()
    );

    // reconstructing face based averaged field
    forAll(field,cellI)
    {
        scalar localAveragedValue(0.);
        forAll(mesh.cells()[cellI], faceI)
        {
            label globalFaceI = mesh.cells()[cellI][faceI];
            localAveragedValue += sumField[globalFaceI]*faceWeights[globalFaceI];
        }   
        // Pout << "area: " << totalArea << " vof: " << localAlphaP << endl;
        field[cellI] = localAveragedValue;
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}

void denseFace::setVectorSum
(
    volVectorField& field,
    double**& value,
    double**const& weight,
    double**const& mask
) const
{
   // get mesh
    const fvMesh& mesh(particleCloud_.mesh());

    // allocate face based averaging array
    List<scalar> sumField1;
    List<scalar> sumField2;
    List<scalar> sumField3;
    List<scalar> faceWeights;
    sumField1.resize(mesh.nFaces());
    sumField2.resize(mesh.nFaces());
    sumField3.resize(mesh.nFaces());
    faceWeights.resize(mesh.nFaces());

    forAll(sumField1,i)
    {
        sumField1[i]     = 0.;
        sumField2[i]     = 0.;
        sumField3[i]     = 0.;
        faceWeights[i]   = 0.5;
    }
        // set correct faceWeights for facesums
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    forAll(boundaryMesh,patchi)
    {
        const label start = boundaryMesh[patchi].start();
        forAll(boundaryMesh[patchi],facei)
        {
            if (!boundaryMesh[patchi].coupled()) faceWeights[facei+start] = 1.0;
        }
    }

    vector valueVec;
    scalar weightP;

    for(int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        if(!checkParticleType(index)) continue; //skip this particle if not correct type
        // start at index 1, since 0 contains the particle cellID
        for(int faceI = 1; faceI < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; faceI++)
        {
            //Info << "subCell=" << subCell << endl;
            label globalFaceI = particleCloud_.cellIDs()[index][faceI];
            
            if (globalFaceI >= 0)
            {
                for(int i = 0; i < 3; i++)
                    valueVec[i] = value[index][i];
                weightP = weight[index][faceI];
                
                sumField1[globalFaceI] += valueVec[0]*weightP;
                sumField2[globalFaceI] += valueVec[1]*weightP;
                sumField3[globalFaceI] += valueVec[2]*weightP;
            }
        }
    }
    // sync faceVoidFraction accross processor boundaries
    syncTools::syncFaceList
    (
        mesh,
        sumField1,
        plusEqOp<scalar>()
    );
    syncTools::syncFaceList
    (
        mesh,
        sumField2,
        plusEqOp<scalar>()
    );
    syncTools::syncFaceList
    (
        mesh,
        sumField3,
        plusEqOp<scalar>()
    );

    // reconstructing face based averaged field
    forAll(field,cellI)
    {
        vector localAveragedValue(vector::zero);
        forAll(mesh.cells()[cellI], faceI)
        {
            label globalFaceI = mesh.cells()[cellI][faceI];
            localAveragedValue[0] += sumField1[globalFaceI]*faceWeights[globalFaceI];
            localAveragedValue[1] += sumField2[globalFaceI]*faceWeights[globalFaceI];
            localAveragedValue[2] += sumField3[globalFaceI]*faceWeights[globalFaceI];
        }   
        // Pout << "area: " << totalArea << " vof: " << localAlphaP << endl;
        field[cellI] = localAveragedValue;
    }

    // correct cell values to patches
    field.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
