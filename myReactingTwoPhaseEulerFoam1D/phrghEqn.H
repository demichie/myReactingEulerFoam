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
    
    // phase 2 = liquid   
    PtrList<volScalarField>& Yl = phase2.Y();
    PtrList<volScalarField>& Yg = phase1.Y();
    phasePairKey keyAnd(phase2.name(), phase1.name());  // (liquid and gas)
    phasePairKey keyIn(phase2.name(), phase1.name(),true); // (liquid in gas)
    volScalarField& alphal = alpha2;
    volScalarField& alphag = alpha1;
    volScalarField& rhol = phase2.thermo().rho();
    volScalarField& rhog = phase1.thermo().rho();
    word inertSpecie(phase2.thermo().lookup("inertSpecie"));

    // phase 1 = liquid                     
    if (phase1.name() == "liquid")
    {       
        PtrList<volScalarField>& Yl = phase1.Y();
        PtrList<volScalarField>& Yg = phase2.Y();
        phasePairKey keyAnd(phase1.name(), phase2.name());  // (gas and liquid)
        phasePairKey keyIn(phase1.name(), phase2.name(),true); // (gas in liquid)
        volScalarField& alphal = alpha1;
        volScalarField& alphag = alpha2;
        volScalarField& rhol = phase1.thermo().rho();
        volScalarField& rhog = phase2.thermo().rho();
        word inertSpecie(phase1.thermo().lookup("inertSpecie"));        
    }      
    
    volScalarField alphal0("alphal0",alphal);
    volScalarField alphag0("alphag0",alphag);
    
    forAll(alphal.boundaryField(), patchI)//boundary field
    {
        Info<<alphal0.boundaryField()[patchI]<<endl;
    } 
      
    volScalarField one
    (
        IOobject
        (
	        "one",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.)
    );
    
    volScalarField zero
    (
        IOobject
        (
	        "zero",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0)
    );

    // total mass fractions of components
    PtrList<volScalarField> Y0(Yl.size());      
    
    forAll(Yl, i)
    {
        Y0.set
	    (
	        i,
	        new volScalarField
            (
	            IOobject
	            (
	                "Y0",
	                runTime.timeName(),
                    mesh,
	                IOobject::NO_READ,
	                IOobject::NO_WRITE
	            ),
	            Yl[i]
            )
        );
    }
            
    PtrList<volScalarField> Ylnew(Yl.size());      
    
    forAll(Yl, i)
    {
        Ylnew.set
	    (
	        i,
	        new volScalarField
            (
	            IOobject
	            (
	                "Ylnew",
	                runTime.timeName(),
                    mesh,
	                IOobject::NO_READ,
	                IOobject::NO_WRITE
	            ),
	            zero
            )
        );
    }
                 
       
    if (equal(runTime.value(), 0))
    {
        p = pMin;
        
        if (pimple.dict().lookupOrDefault<bool>("multicomponentPressureProfile", false))  
        {
              
            forAll(Yl,i)
            {    
                // component in the liquid phase      
                word component(Yl[i].member()); 
                 
                if (component == inertSpecie)
                {}
                else 
                {
                    // equilibrium mass fraction for the component  
                    volScalarField Yfd("Yfd", fluid.Yfd(keyAnd,keyIn,component));
                    Ylnew[i] == Yfd;
                    Info<<"Yfd "<<max(Yfd)<<endl;
                    Info<<"Yfd "<<min(Yfd)<<endl;
                }                      
            } 
            
            // update liquid components             
            volScalarField sumCompl("sumCompl",zero);
            
            volScalarField Yl_inert("Yl_inert",zero);
            
            forAll(Ylnew,i)
            {
                sumCompl += min( Y0[i], Ylnew[i] );
            }
            
            forAll(Yl,i)
            {
                Yl[i] == min(Y0[i], Ylnew[i]);
                
                if (Yl[i].member() == inertSpecie) 
                {
                    Yl[i] == one - sumCompl;
                    Yl_inert = Yl[i];
                }
                Info<<"Yl "<<max(Yl[i])<<endl;
                Info<<"Yl "<<min(Yl[i])<<endl;
            }
            
            // calculate the mass fraction of the liquid and gas phase
            volScalarField Ylphase("Ylphase",zero);
            volScalarField Ygphase("Ygphase",zero);
            
            forAll(Yl,i)
            {
                if (Yl[i].member() == inertSpecie)
                {}
                else
                {
                    Ygphase += Y0[i] - Yl[i]; 
                }
            }
            
            Ygphase = Ygphase/Yl_inert;
            
            Ylphase = one - Ygphase;
            
            Info<<"Ylphase "<<max(Ylphase)<<endl;
            Info<<"Ylphase "<<min(Ylphase)<<endl;
            Info<<"Ygphase "<<max(Ygphase)<<endl;
            Info<<"Ygphase "<<min(Ygphase)<<endl;
            
            // calculate the mass fractions of the components in the gas phase
            forAll(Yg,i)
            {
                forAll(Yl,j)
                {
                    if(Yg[i].member() == Yl[j].member()) 
                    {
                        Info<<"Y0[j]"<<max(Y0[j])<<endl; 
                        Info<<"Y0[j]"<<min(Y0[j])<<endl;
                        Info<<"Yl[j]"<<max(Yl[j])<<endl; 
                        Info<<"Yl[j]"<<min(Yl[j])<<endl;      
                        Yg[i] = one; //( Y0[j] - Yl[j]*Ylphase ) / max( scalar(1e-8), Ygphase );
                    }
                }
            }
            
            forAll(Yg,i)
            {
                Info<<"Yg[i] "<<max(Yg[i])<<endl;
                Info<<"Yg[i] "<<min(Yg[i])<<endl;
            }
            
            // calculate the gas volume fraction
            volScalarField rRho("rRho",zero);
            rRho = rhol/rhog;
            //volScalarField rRho = rhol/rhog;
            
            alphag = Ygphase*rRho / ( one-Ygphase+Ygphase*rRho ); 
            alphal = one - alphag; 
            
            Info<<"alphag "<<max(alphag)<<endl;
            Info<<"alphag "<<min(alphag)<<endl;
            Info<<"alphal "<<max(alphal)<<endl;
            Info<<"alphal "<<min(alphal)<<endl;
            
            // restore the originary boundaryFields of alpha
            /*Info<<"----------- "<<endl;
            forAll(alphag.boundaryField(), patchI)
            {
                Info<<alphag0.boundaryField()[patchI]<<endl;//
                Info<<alpha1.boundaryField()[patchI]<<endl;//
                alphag.boundaryField()[patchI] == alphag0.boundaryField()[patchI];
                //Info<<alphag0.internalField()[patchI]<<endl;//
                Info<<alpha1.boundaryField()[patchI]<<endl;//
            }
            Info<<"----------- "<<endl;*/
                        
            Info<<"alphal "<<max(alphal)<<endl;
            Info<<"alphal "<<min(alphal)<<endl;
            Info<<"alphag "<<max(alphag)<<endl;
            Info<<"alphag "<<min(alphag)<<endl;
            
            Info<<"alpha1 "<<max(alpha1)<<endl;
            Info<<"alpha1 "<<min(alpha1)<<endl;
            Info<<"alpha2 "<<max(alpha2)<<endl;
            Info<<"alpha2 "<<min(alpha2)<<endl;
           
        } 
        
        
        phase1.thermo().correct();
        phase2.thermo().correct();
        
        volScalarField rho("rho", fluid.rho());
        Info<< "min rho " << min(rho).value() <<
               " max rho " << max(rho).value() << endl;
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
     
        
        if (pimple.dict().lookupOrDefault<bool>("multicomponentPressureProfile", false))  
        {
              
            forAll(Yl,i)
            {    
                // component in the liquid phase      
                word component(Yl[i].member()); 
                 
                if (component == inertSpecie)
                {}
                else 
                {
                    // equilibrium mass fraction for the component  
                    volScalarField Yfd("Yfd", fluid.Yfd(keyAnd,keyIn,component));
                    Ylnew[i] == Yfd;
                    Info<<"Yfd "<<max(Yfd)<<endl;
                    Info<<"Yfd "<<min(Yfd)<<endl;
                }                      
            } 
            
            // update liquid components             
            volScalarField sumCompl("sumCompl",zero);
            
            volScalarField Yl_inert("Yl_inert",zero);
            
            forAll(Ylnew,i)
            {
                sumCompl += min( Y0[i], Ylnew[i] );
            }
            
            forAll(Yl,i)
            {
                Yl[i] == min(Y0[i], Ylnew[i]);
                
                if (Yl[i].member() == inertSpecie) 
                {
                    Yl[i] == one - sumCompl;
                    Yl_inert = Yl[i];
                }
                Info<<"Yl "<<max(Yl[i])<<endl;
                Info<<"Yl "<<min(Yl[i])<<endl;
            }
            
            // calculate the mass fraction of the liquid and gas phase
            volScalarField Ylphase("Ylphase",zero);
            volScalarField Ygphase("Ygphase",zero);
            
            forAll(Yl,i)
            {
                if (Yl[i].member() == inertSpecie)
                {}
                else
                {
                    Ygphase += Y0[i] - Yl[i]; 
                }
            }
            
            Ygphase = Ygphase/Yl_inert;
            
            Ylphase = one - Ygphase;
            
            Info<<"Ylphase "<<max(Ylphase)<<endl;
            Info<<"Ylphase "<<min(Ylphase)<<endl;
            Info<<"Ygphase "<<max(Ygphase)<<endl;
            Info<<"Ygphase "<<min(Ygphase)<<endl;
            
            // calculate the mass fractions of the components in the gas phase
            forAll(Yg,i)
            {
                forAll(Yl,j)
                {
                    if(Yg[i].member() == Yl[j].member()) 
                    {
                        Info<<"Y0[j]"<<max(Y0[j])<<endl; 
                        Info<<"Y0[j]"<<min(Y0[j])<<endl;
                        Info<<"Yl[j]"<<max(Yl[j])<<endl; 
                        Info<<"Yl[j]"<<min(Yl[j])<<endl;      
                        Yg[i] = one; //( Y0[j] - Yl[j]*Ylphase ) / max( scalar(1e-8), Ygphase );
                    }
                }
            }
            
            forAll(Yg,i)
            {
                Info<<"Yg[i] "<<max(Yg[i])<<endl;
                Info<<"Yg[i] "<<min(Yg[i])<<endl;
            }
            
            // calculate the gas volume fraction
	    volScalarField rRho("rRho",zero);
            rRho = rhol/rhog;
	    // volScalarField rRho = rhol/rhog;
            
            alphag = Ygphase*rRho / ( one-Ygphase+Ygphase*rRho ); 
            alphal = one - alphag; 
            
            Info<<"alphag "<<max(alphag)<<endl;
            Info<<"alphag "<<min(alphag)<<endl;
            Info<<"alphal "<<max(alphal)<<endl;
            Info<<"alphal "<<min(alphal)<<endl;
            
            // update boundary Fields
            /*Info<<"----------- "<<endl;
            forAll(alphag.boundaryField(), patchI)
            {
                Info<<alphag0.boundaryField()[patchI]<<endl;//
                Info<<alpha1.boundaryField()[patchI]<<endl;//
                alphag.boundaryField()[patchI] == alphag0.boundaryField()[patchI];
                //Info<<alphag0.internalField()[patchI]<<endl;//
                Info<<alpha1.boundaryField()[patchI]<<endl;//
            }
            Info<<"----------- "<<endl;*/
                        
            Info<<"alphal "<<max(alphal)<<endl;
            Info<<"alphal "<<min(alphal)<<endl;
            Info<<"alphag "<<max(alphag)<<endl;
            Info<<"alphag "<<min(alphag)<<endl;
            
            Info<<"alpha1 "<<max(alpha1)<<endl;
            Info<<"alpha1 "<<min(alpha1)<<endl;
            Info<<"alpha2 "<<max(alpha2)<<endl;
            Info<<"alpha2 "<<min(alpha2)<<endl;
           
        } 
        
        
            int integer;
            Info<<"----------------------"<<endl;
            cout << "insert an integer ";
            //cin >> integer;
                     
            phase1.thermo().correct();
            phase2.thermo().correct();
            rho = fluid.rho();

            delta_rho = rho - rho_0;
            rho = rho_0 + 0.5 * delta_rho;

            rho_0 = rho;
            p_0 = p; // only for Info

            // Limit static pressure
            p = max(ph_rgh + rho*gh, pMin);

            // Limit ph_rgh
            ph_rgh = p - rho*gh;

            Info<< "min p " << min(p).value() <<
                   " max p " << max(p).value() << endl;

            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;

            Info<< "min delta_rho " << min(delta_rho).value() <<
                   " max delta_rho " << max(delta_rho).value() << endl;

            delta_p = p - p_0;
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
            
        if (pimple.dict().lookupOrDefault<bool>("multicomponentPressureProfile", false))  
        {
              
            forAll(Yl,i)
            {    
                // component in the liquid phase      
                word component(Yl[i].member()); 
                 
                if (component == inertSpecie)
                {}
                else 
                {
                    // equilibrium mass fraction for the component  
                    volScalarField Yfd("Yfd", fluid.Yfd(keyAnd,keyIn,component));
                    Ylnew[i] == Yfd;
                    Info<<"Yfd "<<max(Yfd)<<endl;
                    Info<<"Yfd "<<min(Yfd)<<endl;
                }                      
            } 
            
            // update liquid components             
            volScalarField sumCompl("sumCompl",zero);
            
            volScalarField Yl_inert("Yl_inert",zero);
            
            forAll(Ylnew,i)
            {
                sumCompl += min( Y0[i], Ylnew[i] );
            }
            
            forAll(Yl,i)
            {
                Yl[i] == min(Y0[i], Ylnew[i]);
                
                if (Yl[i].member() == inertSpecie) 
                {
                    Yl[i] == one - sumCompl;
                    Yl_inert = Yl[i];
                }
                Info<<"Yl "<<max(Yl[i])<<endl;
                Info<<"Yl "<<min(Yl[i])<<endl;
            }
            
            // calculate the mass fraction of the liquid and gas phase
            volScalarField Ylphase("Ylphase",zero);
            volScalarField Ygphase("Ygphase",zero);
            
            forAll(Yl,i)
            {
                if (Yl[i].member() == inertSpecie)
                {}
                else
                {
                    Ygphase += Y0[i] - Yl[i]; 
                }
            }
            
            Ygphase = Ygphase/Yl_inert;
            
            Ylphase = one - Ygphase;
            
            Info<<"Ylphase "<<max(Ylphase)<<endl;
            Info<<"Ylphase "<<min(Ylphase)<<endl;
            Info<<"Ygphase "<<max(Ygphase)<<endl;
            Info<<"Ygphase "<<min(Ygphase)<<endl;
            
            // calculate the mass fractions of the components in the gas phase
            forAll(Yg,i)
            {
                forAll(Yl,j)
                {
                    if(Yg[i].member() == Yl[j].member()) 
                    {
                        Info<<"Y0[j]"<<max(Y0[j])<<endl; 
                        Info<<"Y0[j]"<<min(Y0[j])<<endl;
                        Info<<"Yl[j]"<<max(Yl[j])<<endl; 
                        Info<<"Yl[j]"<<min(Yl[j])<<endl;      
                        Yg[i] = one; //( Y0[j] - Yl[j]*Ylphase ) / max( scalar(1e-8), Ygphase );
                    }
                }
            }
            
            forAll(Yg,i)
            {
                Info<<"Yg[i] "<<max(Yg[i])<<endl;
                Info<<"Yg[i] "<<min(Yg[i])<<endl;
            }
            
            // calculate the gas volume fraction
            volScalarField rRho("rRho",zero);
            rRho = rhol/rhog;
	    // volScalarField rRho = rhol/rhog;
            
            alphag = Ygphase*rRho / ( one-Ygphase+Ygphase*rRho ); 
            alphal = one - alphag; 
            
            Info<<"alphag "<<max(alphag)<<endl;
            Info<<"alphag "<<min(alphag)<<endl;
            Info<<"alphal "<<max(alphal)<<endl;
            Info<<"alphal "<<min(alphal)<<endl;
            
            // update boundary Fields
            /*Info<<"----------- "<<endl;
            forAll(alphag.boundaryField(), patchI)
            {
                Info<<alphag0.boundaryField()[patchI]<<endl;//
                Info<<alpha1.boundaryField()[patchI]<<endl;//
                alphag.boundaryField()[patchI] == alphag0.boundaryField()[patchI];
                //Info<<alphag0.internalField()[patchI]<<endl;//
                Info<<alpha1.boundaryField()[patchI]<<endl;//
            }
            Info<<"----------- "<<endl;*/
                        
            Info<<"alphal "<<max(alphal)<<endl;
            Info<<"alphal "<<min(alphal)<<endl;
            Info<<"alphag "<<max(alphag)<<endl;
            Info<<"alphag "<<min(alphag)<<endl;
            
            Info<<"alpha1 "<<max(alpha1)<<endl;
            Info<<"alpha1 "<<min(alpha1)<<endl;
            Info<<"alpha2 "<<max(alpha2)<<endl;
            Info<<"alpha2 "<<min(alpha2)<<endl;
           
        } 
        
        
            int integer;
            Info<<"----------------------"<<endl;
            cout << "insert an integer ";
            //cin >> integer;
        
            phase1.thermo().correct();
            phase2.thermo().correct();
            
            // Correct p_rgh for consistency with p and the updated densities
            rho = fluid.rho();
            Info<<"rho"<<max(rho)<<min(rho)<<endl;

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
                
        forAll(Yl,i)
        {
            Yl[i].write();
        }
        
        forAll(Yg,i)
        {
            Yg[i].write();
        }
        
        alpha1.write();
        alpha2.write();
               
    }
}
