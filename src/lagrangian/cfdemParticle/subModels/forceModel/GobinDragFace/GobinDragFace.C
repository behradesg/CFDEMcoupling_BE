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

#include "GobinDragFace.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(GobinDragFace, 0);

addToRunTimeSelectionTable
(
    forceModel,
    GobinDragFace,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
GobinDragFace::GobinDragFace
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    phi_(readScalar(propsDict_.lookup("phi"))),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    minVoidfraction_(propsDict_.lookupOrDefault<scalar>("minVoidfraction",0.1)),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); // first entry must  be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("betaW");
    particleCloud_.probeM().scalarFields_.append("betaG");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().writeHeader();

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_IMPL_FORCE_DEM,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(SW_SCALAR_VISCOSITY,true); // activate scalarViscosity switch
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_ = scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_ = scalar(readScalar(propsDict_.lookup("scaleDrag")));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

GobinDragFace::~GobinDragFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void GobinDragFace::setForce() const
{
    // get mesh
    const fvMesh& mesh(voidfraction_.mesh());
    if (scaleDia_ > 1)
    {
        Info << "GobinDragFace using scale = " << scaleDia_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "GobinDragFace using scale from liggghts cg = " << scaleDia_ << endl;
    }

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);
    label cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar betaW(0);
    scalar betaG(0);
    scalar Vs(0);
    scalar localPhiP(0);

    scalar CdMagUrLag(0);       //Cd of the very particle
    
     // Fields for face-based interpolation
    surfaceScalarField voidFractionFace
    (
        fvc::interpolate(voidfraction_)
    );
    surfaceVectorField UfFace
    (
        fvc::interpolate(U_)
    );
    List<scalar> voidFractionList;
    List<vector> listUf;

    voidFractionList.resize(mesh.nFaces());
    listUf.resize(mesh.nFaces());

    // internalField
    forAll(voidFractionList,i)
    {
        voidFractionList[i] = voidFractionFace[i];
        listUf[i]           = UfFace[i];
    }
    // Boundary fields
    const surfaceScalarField::Boundary& voidFb = voidFractionFace.boundaryField();
    const surfaceVectorField::Boundary& Ufb = UfFace.boundaryField();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    forAll(voidFb,patchi)
    {
        const label start = boundaryMesh[patchi].start();
        const scalarField& voidFbP = voidFb[patchi];
        const vectorField& UfbP = Ufb[patchi];
        forAll(boundaryMesh[patchi],facei)
        {
            voidFractionList[facei+start] = voidFbP[facei];
            listUf[facei+start]           = UfbP[facei];
            // Pout << listUf[facei+start] << " p: " << patchi << endl;
        }
    } 

    #include "setupProbeModel.H"

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector::zero;
            dragExplicit = vector::zero;
            Vs = 0;
            Ufluid = vector::zero;
            voidfraction = 0;
            dragCoefficient = 0;

            if (cellI > -1) // particle Found
            {

                if(forceSubM(0).interpolation())
                {
                    scalar weightP;
                    // face based interpolation
                    for(int faceI = 1; faceI < particleCloud_.voidFractionM().cellsPerParticle()[index][0]; faceI++)
                    {
                        label globalFaceI = particleCloud_.cellIDs()[index][faceI];
                        weightP           = particleCloud_.particleWeights()[index][faceI];
                        Ufluid           += weightP*listUf[globalFaceI];
                        voidfraction     += weightP*voidFractionList[globalFaceI];
                    }
                    if (voidfraction > 1.00) voidfraction = 1.0;
                    if (voidfraction < minVoidfraction_) voidfraction = minVoidfraction_;
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                magUr = mag(Ur);
                ds = 2.0*particleCloud_.radius(index);
                rho = rhoField[cellI];
                nuf = nufField[cellI];

                Rep = 0.0;
                localPhiP = 1.0f-voidfraction+SMALL;
                Vs = ds*ds*ds*M_PI/6.0;

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and per m³ PARTICLE)
                Rep = ds/scaleDia_*voidfraction*magUr/nuf;
                CdMagUrLag = (24.0*nuf/(ds/scaleDia_*voidfraction)) //1/magUr missing here, but compensated in expression for betaP!
                            *(scalar(1.0)+0.15*Foam::pow(Rep, 0.687));

                betaW = 0.75
                       *(                                  //this is betaG = beta / localPhiP!
                            rho*voidfraction*CdMagUrLag
                           /(ds/scaleDia_*Foam::pow(voidfraction,2.65))
                        );

                betaG = (150.0 * localPhiP*nuf*rho)          //this is betaG = beta / localPhiP!
                            /  (voidfraction*ds/scaleDia_*phi_*ds/scaleDia_*phi_)
                        +
                            (1.75 * magUr * rho)
                            /((ds/scaleDia_*phi_));

                // calc particle's drag
                dragCoefficient = Vs*Foam::min(betaW,betaG)*scaleDrag_;
                if (modelType_ == "B")
                    dragCoefficient /= voidfraction;

                drag = dragCoefficient * Ur;

                // explicitCorr
                forceSubM(0).explicitCorr(drag,dragExplicit,dragCoefficient,Ufluid,U_[cellI],Us,UsField_[cellI],forceSubM(0).verbose());

                if (forceSubM(0).verbose() && index >= 0 && index < 2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds/scale = " << ds/scaleDia_ << endl;
                    Pout << "phi = " << phi_ << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "betaW = " << betaW << endl;
                    Pout << "betaG = " << betaG << endl;
                    Pout << "drag = " << drag << endl;
                }

                //Set value fields and write the probe
                if (probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(betaW);
                    sValues.append(betaG);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);

        //}// end if mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
