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
#include "GaussVoidFractionUR.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(GaussVoidFractionUR, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    GaussVoidFractionUR,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
GaussVoidFractionUR::GaussVoidFractionUR
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(propsDict_.lookupOrDefault<scalar>("alphaMin",0.01)),
    stencil_(propsDict_.lookupOrDefault<scalar>("stencil",1)),
    alphaLimited_(0)
{
    Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    FatalError << "\n\n This model does not yet work properly! \n\n" << endl;
    //reading maxCellsPerParticle from dictionary
    maxCellsPerParticle_ = propsDict_.lookupOrDefault<scalar>("maxCellsPerParticle",64);

    if(alphaMin_ > 1 || alphaMin_ < 0.01) { FatalError << "alphaMin must have a value between 0.01 and 1.0." << abort(FatalError); }

    checkWeightNporosity(propsDict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

GaussVoidFractionUR::~GaussVoidFractionUR()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void GaussVoidFractionUR::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV)
{
    scalar cellSize(0);
    scalar radius(0);
    scalar volume(0);
    scalar scaleVol = weight();
    scalar stencil(0.);

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
            cellSize = Foam::cbrt(particleCloud_.mesh().V()[particleCenterCellID]);
            stencil = stencil_*cellSize;
            scalar stencil2 = Foam::sqr(stencil);
            scalar stencil3 = stencil2*stencil;
            scalar stencil5 = stencil2*stencil3;

            vector positionCenter = particleCloud_.position(index);
            scalar core;
            scalar dist2;

            if (particleCenterCellID >= 0)
            {
                if (multiWeights_) scaleVol = weight(index);
                radius = particleCloud_.radius(index);
                volume = constant::mathematical::fourPiByThree*radius*radius*radius*scaleVol;

                labelHashSet hashSet;

                //determining label and degree of coveredness of cells covered by the particle
                buildLabelHashSet(stencil, positionCenter, particleCenterCellID, hashSet);

                //generating list with cell and subcells
                label hashSetLength = hashSet.size();
                if (hashSetLength > maxCellsPerParticle_)
                {
                    FatalError<< "particle algo found more cells ("<< hashSetLength
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                }
                else if (hashSetLength > 0)
                {
                    cellsPerParticle()[index][0] = hashSetLength;

                    scalar normalization(SMALL);

                    // compute kernel normalization that sum of weights = 1
                    // use kernel function of Deen et al. 2004, Chemical Engineering Science 59 (2004) 1853–1861
                    for(label i = 0; i < hashSetLength; i++)
                    {
                        label cellI = hashSet.toc()[i];
                        dist2 = Foam::magSqr(particleCloud_.mesh().C()[cellI] - positionCenter);
                        core = 0.9375*(Foam::sqr(dist2)/stencil5 - 2.0*dist2/stencil3 + 1.0/stencil)*Foam::pos(stencil2 - dist2);
                        normalization += core;
                    }

                    //==========================//
                    //setting the voidfractions
                    //==========================//
                    // deleting the cell containing the center of the particle
                    // particleCenterCellID should get be first subCell in particleCloud_.cellIDs()
                    hashSet.erase(particleCenterCellID);
                    dist2 = Foam::magSqr(particleCloud_.mesh().C()[particleCenterCellID] - positionCenter);
                    core = 0.9375*(Foam::sqr(dist2)/stencil5 - 2.0*dist2/stencil3 + 1.0/stencil)*Foam::pos(stencil2 - dist2);
                    particleCloud_.cellIDs()[index][0] = particleCenterCellID; 
                    particleWeights[index][0] = core/normalization;
                    if (particleWeights[index][0] < SMALL) Warning << "too small particles weight: dist: " << sqrt(dist2) << endl;
                    particleVolumes[index][0] = volume*core/normalization;
                    voidfractionNext_[particleCenterCellID] -= particleVolumes[index][0]/particleCloud_.mesh().V()[particleCenterCellID];

                    // loop over neighbour cells
                    for(label i=0; i<hashSetLength-1;i++)
                    {
                        label cellI = hashSet.toc()[i];
                        dist2 = Foam::magSqr(particleCloud_.mesh().C()[cellI] - positionCenter);
                        core = 0.9375*(Foam::sqr(dist2)/stencil5 - 2.0*dist2/stencil3 + 1.0/stencil)*Foam::pos(stencil2 - dist2);
                        particleCloud_.cellIDs()[index][i+1] = cellI; //adding subcell represenation
                        particleWeights[index][i+1] = core/normalization;
                        if (particleWeights[index][i+1] < SMALL) Warning << "too small particles weight: dist: " << sqrt(dist2) << endl;
                        particleVolumes[index][i+1] = volume*core/normalization;
                        voidfractionNext_[cellI] -= particleVolumes[index][i+1]/particleCloud_.mesh().V()[cellI];
                    }
                    particleV[index][0] = volume;

                    // debug
                    if(index==0)
                    {
                        Info << "particle 0 is represented by " << hashSetLength << "cells; normalization: " << normalization << endl;
                    }
                    //==========================//
                }//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles

    // limiting voidfraction
    voidfractionNext_.max(alphaMin_);
    voidfractionNext_.min(1.0);
    // correct Boundary Conditions
    voidfractionNext_.correctBoundaryConditions();

    //bringing eulerian field to particle array
    for(label index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        for(label subcell = 0; subcell < cellsPerParticle()[index][0]; subcell++)
        {
            label cellID = particleCloud_.cellIDs()[index][subcell];

            if(cellID >= 0)
            {
                 // set particle based voidfraction
                 voidfractions[index][subcell] = voidfractionNext_[cellID];
                 //Info<<"setting the voidfraction, index = "<<index<<endl;
            }
            else
            {
                voidfractions[index][subcell] = -1.;
            }
        }
    }
}

void GaussVoidFractionUR::buildLabelHashSet
(
    const scalar radius,
    const vector position,
    const label cellID,
    labelHashSet& hashSet
)const
{
    hashSet.insert(cellID);
    //Info<<"cell inserted"<<cellID<<endl;
    const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
    forAll(nc,i)
    {
        label neighbor=nc[i];
        if(!hashSet.found(neighbor) && mag(position - particleCloud_.mesh().C()[neighbor]) < radius)
        {
            buildLabelHashSet(radius,position,neighbor,hashSet);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
