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
#include "mathExtra.H"
#include "inverseDistanceVoidFractionFace.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(inverseDistanceVoidFractionFace, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    inverseDistanceVoidFractionFace,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inverseDistanceVoidFractionFace::inverseDistanceVoidFractionFace
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(propsDict_.lookupOrDefault<scalar>("alphaMin",0.01)),
    alphaLimited_(0)
{
    Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    FatalError << "\n\n This model does not yet work properly! \n\n" << endl;
    //reading maxCellsPerParticle from dictionary
    maxCellsPerParticle_ = propsDict_.lookupOrDefault<scalar>("maxCellsPerParticle",7);

    if(alphaMin_ > 1 || alphaMin_ < 0.01) { FatalError << "alphaMin must have a value between 0.01 and 1.0." << abort(FatalError); }

    checkWeightNporosity(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inverseDistanceVoidFractionFace::~inverseDistanceVoidFractionFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void inverseDistanceVoidFractionFace::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV)
{
    // get mesh
    const fvMesh& mesh(particleCloud_.mesh());

    scalar scaleVol = weight();
    List<scalar> particleVolumeFractionF;
    List<scalar> faceWeight;
    List<scalar> listMagSf;
    List<vector> listCf;
    List<vector> listSf;
    particleVolumeFractionF.resize(mesh.nFaces());
    faceWeight.resize(mesh.nFaces());
    listMagSf.resize(mesh.nFaces());
    listCf.resize(mesh.nFaces());
    listSf.resize(mesh.nFaces());

    forAll(particleVolumeFractionF,i)
    {
        particleVolumeFractionF[i]= 0.;
        faceWeight[i] = 0.5;
    }
    
    const surfaceVectorField& Cf = mesh.Cf();   // Face center coordinates
    const surfaceScalarField& magSf = mesh.magSf();           // face areas
    const surfaceVectorField& Sf = mesh.Sf();           // face vectors
    forAll(Cf, i)
    {
        listMagSf[i] = magSf[i];
        listCf[i] = Cf[i];
        listSf[i] = Sf[i]/magSf[i];
    }
    const surfaceScalarField::Boundary& magSfb = mesh.magSf().boundaryField();
    const surfaceVectorField::Boundary& Sfb = mesh.Sf().boundaryField();
    const surfaceVectorField::Boundary& Cfb = mesh.Cf().boundaryField();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    forAll(magSfb,patchi)
    {
        const label start = boundaryMesh[patchi].start();
        const scalarField& magSfbP = magSfb[patchi];
        const vectorField& SfbP = Sfb[patchi];
        const vectorField& CfbP = Cfb[patchi];
        forAll(boundaryMesh[patchi],facei)
        {
            if (!boundaryMesh[patchi].coupled()) faceWeight[facei+start] = 1.0;
            listMagSf[facei+start] = magSfbP[facei];
            listCf[facei+start] = CfbP[facei];
            listSf[facei+start] = SfbP[facei]/magSfbP[facei];
        }
    }

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            //reset
            for(int subcell = 0; subcell < maxCellsPerParticle_; subcell++)
            {
                particleWeights[index][subcell] = 0;
                particleVolumes[index][subcell] = 0;
            }
            cellsPerParticle()[index][0] = 1;
            particleV[index][0]          = 0;

            //collecting data
            label particleCenterCellID = particleCloud_.cellIDs()[index][0];

            if (particleCenterCellID >= 0)
            {
                scalar radius(particleCloud_.radius(index));
                scalar volume(constant::mathematical::fourPiByThree*radius*radius*radius*scaleVol);

                vector positionCenter(particleCloud_.position(index));
                scalar core;
                scalar dist2;

                if (multiWeights_) scaleVol = weight(index);
                scalar normalization(SMALL);
                label hits(0);

                forAll(mesh.cells()[particleCenterCellID], facei)
                {
                    const label newFaceI = mesh.cells()[particleCenterCellID][facei];
                    dist2 = Foam::magSqr((listCf[newFaceI] - positionCenter)&listSf[newFaceI]) + SMALL;
                    core = 
                        (1.0/dist2)
                       *listMagSf[newFaceI];
                    normalization += core;
                    hits++;
                }
                // increase hits by 1 to account for cellID at index 0
                hits++;
                cellsPerParticle()[index][0] = hits;

                if (hits > maxCellsPerParticle_)
                {
                    FatalError<< "particle algo found more cells ("<< hits
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                }
                else if (hits > 1)
                {

                    //==========================//
                    //setting the voidfractions
                    //==========================//
                    particleWeights[index][0] = 1;
                    particleVolumes[index][0] = volume;
                    particleV[index][0]       = volume;
                    forAll(mesh.cells()[particleCenterCellID], facei)
                    {
                        const label newFaceI = mesh.cells()[particleCenterCellID][facei];

                        dist2 = Foam::magSqr((listCf[newFaceI] - positionCenter)&listSf[newFaceI]) + SMALL;
                        core = 
                            (1.0/dist2)
                           *listMagSf[newFaceI];
                        particleVolumeFractionF[newFaceI] +=
                            faceWeight[newFaceI]
                           *(volume/mesh.V()[particleCenterCellID])
                           *(core/normalization);
                        particleCloud_.cellIDs()[index][facei + 1] = newFaceI; //adding faces representation
                        particleWeights[index][facei + 1]          = core/normalization;
                        particleVolumes[index][facei + 1]          = volume*core/normalization;
                    }
                }
            } // cellId > -1
        //}// end if masked
    }// end loop all particles

    // sync faceVoidFraction accross processor boundaries
    syncTools::syncFaceList
    (
        mesh,
        particleVolumeFractionF,
        plusEqOp<scalar>()
    );

    // reconstructing voidFraction
    forAll(voidfractionNext_,celli)
    {
        scalar localAlphaP(0.);
        forAll(mesh.cells()[celli], facei)
        {
            const label newFaceI = mesh.cells()[celli][facei];
            localAlphaP += particleVolumeFractionF[newFaceI];
        }   
        // Pout << "area: " << totalArea << " vof: " << localAlphaP << endl;
        voidfractionNext_[celli] = 1.0 - localAlphaP;
    }
    // limiting voidfraction
    voidfractionNext_.max(alphaMin_);
    voidfractionNext_.min(1.0);
    // correct Boundary Conditions
    voidfractionNext_.correctBoundaryConditions();

    //bringing eulerian field to particle array
    for(label index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        label cellID = particleCloud_.cellIDs()[index][0];

        if(cellID >= 0)
        {
            // set particle based voidfraction
            voidfractions[index][0] = voidfractionNext_[cellID];
            //Info<<"setting the voidfraction, index = "<<index<<endl;
        }
        else
        {
            voidfractions[index][0] = -1.;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
