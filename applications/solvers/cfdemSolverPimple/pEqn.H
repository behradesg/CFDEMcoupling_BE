{
    surfaceScalarField rAUfvoidfractionf("(voidfraction2|A(U))",rAUf*voidfractionf);

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(HbyA)
          + rAUfvoidfractionf*fvc::ddtCorr(U, phi/voidfractionf) //correction
        )
    );

    if (p.needReference())
    {
        adjustPhi(phiHbyA, U, p);
    }

    //explicit contribution
    phiHbyA += phicForces;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAU,MRF);

    if (modelType=="A")
        rAUfvoidfractionf *= voidfractionf;

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUfvoidfractionf, p) == fvc::div(voidfractionf*phiHbyA) + fvc::ddt(voidfraction)
        );
        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = voidfractionf*phiHbyA - pEqn.flux();

            #include "continuityErrorPhiPU.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            U = rAU*fvc::reconstruct((phi/voidfractionf)/rAUf);

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }
}
