particleCloud.otherForces(fOther);

{
    dimensionedVector ftotOther = fvc::domainIntegrate(fOther);
    Info << "TotalOtherForce: " << ftotOther.value() << endl;
}

fvVectorMatrix UEqn
(
    fvm::ddt(voidfraction,U)
  + fvm::div(phi,U)
  + fvm::SuSp(-fvc::ddt(voidfraction) - fvc::div(phi),U) // Continuity errors
    // in case of "periodic box" simulations the viscous term can be commented
    // during testing the effect of this term on the results was negligible (about 1-2%)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  + fvm::Sp(Ksl/rho,U)
  ==
    fvOptions(U)
);
UEqn.relax();

fvOptions.constrain(UEqn);

volScalarField rAU = 1.0/UEqn.A();
surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

surfaceScalarField phicForces
(
    fvc::flux(rAU*(fOther)/rho)
);

if (!dragCouplingNew)
{
    phicForces += fvc::flux(rAU*(Fexp + Ksl*Us)/rho);
}
else
{
    phicForces += fvc::flux(rAU*(Fexp + Ksl*USmooth)/rho);
}

if (pimple.momentumPredictor())
{
    if (modelType=="B" || modelType=="Bfull")
    {
        solve
        (
            UEqn 
         == 
            fvc::reconstruct
            (
                phicForces/rAUf - fvc::snGrad(p)*mesh.magSf()
            )
        );
    }
    else 
    {
        solve
        (
            UEqn 
         == 
            fvc::reconstruct
            (
                phicForces/rAUf - fvc::snGrad(p)*voidfractionf*mesh.magSf()
            )
        );
    }
}
fvOptions.correct(U);

