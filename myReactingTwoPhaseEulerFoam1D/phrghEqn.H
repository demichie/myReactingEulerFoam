if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialization", false))
{
    volScalarField& ph_rgh = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    if (equal(runTime.value(), 0))
    {
        p = pMin;

        phase1.thermo().correct();
        phase2.thermo().correct();
        volScalarField rho("rho", fluid.rho());
        volScalarField rho_0("rho_0", rho);
        volScalarField delta_rho("delta_rho", rho);

        volScalarField p_0("p_0", p);
        volScalarField delta_p("delta_p", p);


        volVectorField U("U", fluid.U());
        // Limit static pressure
        p = max(ph_rgh + rho*gh, pMin);

        Info<< "min p " << min(p).value() <<
               " max p " << max(p).value() << endl;


        Info<< "min rho " << min(rho).value() <<
               " max rho " << max(rho).value() << endl;

        Info<< "min gh " << min(gh).value() <<
               " max gh " << max(gh).value() << endl;

        label nCorr0
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticIterations", 5)
        );


        for (label i=0; i<nCorr0; i++)
        {
 
            phase1.thermo().correct();
            phase2.thermo().correct();
            rho = fluid.rho();

            delta_rho = rho - rho_0;
            rho = rho_0 + 0.5 * delta_rho;

            rho_0 = rho;
            p_0 = p;

            // Limit static pressure
            p = max(ph_rgh + rho*gh, pMin);
            delta_p = p - p_0;

            // Limit ph_rgh
            ph_rgh = p - rho*gh;

            Info<< "min p " << min(p).value() <<
                   " max p " << max(p).value() << endl;

            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;

            Info<< "min delta_rho " << min(delta_rho).value() <<
                   " max delta_rho " << max(delta_rho).value() << endl;

            Info<< "min delta_p " << min(delta_p).value() <<
                   " max delta_p " << max(delta_p).value() << endl;

        }



        label nCorr
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticCorrectors", 5)
        );

        for (label i=0; i<nCorr; i++)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));
            surfaceScalarField pf("rhof", fvc::interpolate(p));

            surfaceScalarField phig
            (
                "phig",
               //-rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
               -ghf*fvc::snGrad(rho)*mesh.magSf()/rhof
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);

            fvScalarMatrix ph_rghEqn
            (
                //fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
                fvm::laplacian(1/rhof,ph_rgh) == fvc::div(phig)
            );
  
            //ph_rghEqn.relax();

            ph_rghEqn.solve();

            // Update and limit the static pressure
            p = max(ph_rgh + rho*gh, pMin);

            // Limit ph_rgh
            ph_rgh = p - rho*gh;

            phase1.thermo().correct();
            phase2.thermo().correct();
            
            // Correct p_rgh for consistency with p and the updated densities
            rho = fluid.rho();

            ph_rgh = p - rho*gh;

            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh) - min(ph_rgh)).value() << endl;

            Info<< "min p " << min(p).value() <<
                   " max p " << max(p).value() << endl;


            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;

        }


        volVectorField gradP
        (
            "gradP",
            fvc::grad(p)
        );
			

        // Limit ph_rgh
        ph_rgh = p - rho*gh;



        ph_rgh.write();
        p.write();
        rho.write();

        p_rgh = ph_rgh;
        p_rgh.write();
        gradP.write();

    }
}
