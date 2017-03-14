// Search for wall and patch (inlet/outlet) patches faces and store normals

label faceId(-1);
label patchId(-1);
label nPatchFaces(0);
vector patchNormal(vector::zero);
vector secondNormal(vector::zero);

const fvPatchList& patches = mesh.boundary();

forAll(patches, patchi)
{
    const fvPatch& currPatch = patches[patchi];
    
    Info<< " mesh type = " << mesh.boundary()[patchi].type() << endl; 

    

    // search for patches of type "patch" and count the
    // numebr of faces in the patches (if n=2 then 1D)
    if(mesh.boundary()[patchi].type() == "patch")  
    {
        const vectorField nf(currPatch.nf());
        
        forAll(nf, facei)
        {
            nPatchFaces++;
            
            if (nPatchFaces == 1)
            {
                patchNormal = -nf[facei];
                faceId = facei;
                patchId = patchi;
            }
            else if (nPatchFaces == 2)
            {
                const vector patchNormal2 = -nf[facei];
                
                //- Check that patch faces are parallel
                if
                    (
                        mag(patchNormal & patchNormal2) > 1.01
                        || mag(patchNormal & patchNormal2) < 0.99
                    )
                {
                    FatalErrorIn(args.executable())
                        << "patch faces are not parallel for patches "
                            << patches[patchId].name() << " and "
                            << currPatch.name() << nl
                            << exit(FatalError);
                }
            }
            
        }
    }
}

//Info<< "nPatchFaces" << nPatchFaces << endl;

if (nPatchFaces == 0)
{
    FatalErrorIn(args.executable()) << "No patch patches identified"
        << exit(FatalError);
}

if (nPatchFaces > 2)
{
    Info<< "2D/3D Simulation" << endl;
    OneD_Simulation = false;
}
else
{   
    OneD_Simulation = true;
    Info << "1D Simulation" << endl;

    const List<vector>& patchFound = mesh.boundaryMesh()[patchId].localPoints();
    //Info << "### localPoints" << patchFound << endl;

    secondNormal = patchFound[1] - patchFound[0];
    secondNormal = secondNormal /mag(secondNormal);

    vector thirdNormal = patchNormal^secondNormal;


    //Info<< " secondNormal = " << secondNormal << endl;
    //Info<< " thirdNormal = " << thirdNormal << endl;

    scalar diam1 = Foam::max(patchFound & secondNormal) - Foam::min(patchFound & secondNormal);
    //Info<< " diam1 " << diam1 << endl;

    scalar diam2 = Foam::max(patchFound & thirdNormal) - Foam::min(patchFound & thirdNormal);
    //Info<< " diam2 " << diam2 << endl;

    OneD_Radius.value() = 0.5*max(diam1,diam2);

    Info<< " Conduit radius 0 " << OneD_Radius << endl;
    
}
