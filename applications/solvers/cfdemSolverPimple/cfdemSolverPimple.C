/*---------------------------------------------------------------------------*\
    Open Source CFD-DEM coupling

    Copyright (C) 2023  Behrad Esgandari, JKU Linz, Austria
-------------------------------------------------------------------------------
License

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

Application
    cfdemSolverPimple

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pimpleFoam in OpenFOAM(R) 6.0,
    where additional functionality for CFD-DEM coupling is added.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#include "cfdemCloud.H"
#include "implicitCouple.H"
#include "explicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloud particleCloud(mesh);
    #include "checkModelType.H"

    // switch for periodic box simulations
    Switch periodicBoxSwitch
    (
        pimple.dict().lookupOrDefault<Switch>("periodicBoxSwitch", false)
    );
    Switch dragCouplingNew
    (
        pimple.dict().lookupOrDefault<Switch>("dragCouplingNew", true)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);
        
        USmooth = U;
        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenReferenceField
            (
                USmooth,
                voidfraction
            );
            particleCloud.smoothingM().smoothen
            (
                particleCloud.forceM(0).impParticleForces()
            );
        }

        Info << "update Ksl and Fexp" << endl;
        if (!dragCouplingNew)
        {

            Ksl = particleCloud.momCoupleM(0).impMomSource();
            Fexp = -particleCloud.momCoupleM(0).expMomSource();
        } 
        else 
        {
            Ksl = particleCloud.momCoupleM(0).CdS();
            volVectorField& impParticleForces(particleCloud.forceM(0).impParticleForces());
            //volScalarField posKsl(pos(Ksl));
            //Fexp.ref() = -posKsl.ref()*impParticleForces.ref()/mesh.V();
	    Fexp.ref() = -impParticleForces.ref()/mesh.V();
            Fexp       -= particleCloud.momCoupleM(0).expMomSource();
            if (hasEvolved)
            {
                particleCloud.smoothingM().smoothen(Ksl);
            } 
            //Ksl *= posKsl;
        }
        Ksl.correctBoundaryConditions();
        Fexp.correctBoundaryConditions();

        //Force Checks
        dimensionedVector fImpTotal = fvc::domainIntegrate(Ksl*(Us-USmooth));
        Info << "TotalForceImp: " << fImpTotal.value() << endl;
        dimensionedVector fExpTotal = fvc::domainIntegrate(Fexp);
        Info << "TotalForceExp: " << fExpTotal.value() << endl;
        // Monitor domainIntegral of voidfraction
        dimensionedScalar meanAlphaP = 
            fvc::domainIntegrate(1.0-voidfraction)
           /fvc::domainIntegrate(unity);
        Info << "meanParticleVolumeFraction: " << meanAlphaP.value() << endl;

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PIMPLE corrector
            while (pimple.loop())
            {
                // Momentum predictor
                #include "UEqn.H"

                // --- Inner PIMPLE loop

                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }

            laminarTransport.correct();
            turbulence->correct();
        }
        else
        {
            Info << "skipping flow solution." << endl;
        }

        if (periodicBoxSwitch)
        {
           #include "periodicBoxProperties.H"
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
