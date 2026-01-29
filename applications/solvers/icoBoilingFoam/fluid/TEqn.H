{

    kappa = alpha1*k1 + alpha2*k2;
    surfaceScalarField kappaf(fvc::interpolate(kappa));

    // Generate mixture phase Cp
    volScalarField Cpmix
    (
        IOobject
        (
            "Cpmix",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*cp1 + alpha2*cp2,
        alpha1.boundaryField().types()
    );
    Cpmix.oldTime();

    // The sensible heat source MUST BE calculated before "computeMassSource.H"
    // or the wrong T field will be obtaied!!!
    // Assumption: vapor appear in Tsat
    volScalarField impSenTSource(-continuitySource * Cpmix);
    volScalarField expSenTSource(continuitySource *Cpmix * Tsat);

    #include "computeMassSource.H"

    volScalarField impTSource(impLatTSource + impSenTSource);
    volScalarField expTSource(expLatTSource + expSenTSource);

    fvScalarMatrix TEqn
     (
         fvm::ddt(rhoCp,T)
       + fvm::div(rhoCpPhi,T)
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaf,T)
       ==
         fvm::Sp(impTSource, T)
       + expTSource
     );

     TEqn.relax();
     TEqn.solve();

     Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() <<endl;

     //update continuity source with new T
    #include "computeMassSource.H"

    if(mixture.found("writeQwall"))
    {
        qwFluid[i] = -kappa*fvc::grad(T);
    }
}
