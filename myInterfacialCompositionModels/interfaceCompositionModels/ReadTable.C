/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ReadTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::ReadTable<Thermo, OtherThermo>::ReadTable
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),

    N_(dict.lookup("N")),
    Nspecie_(dict.lookup("Nspecie")),
    speciemax_(dict.lookup("speciemax")),
    speciemin_(dict.lookup("speciemin")),
    Np_(dict.lookup("Np")),
    pmax_(dict.lookup("pmax")),
    pmin_(dict.lookup("pmin")),
    NT_(dict.lookup("NT")),
    Tmax_(dict.lookup("Tmax")),
    Tmin_(dict.lookup("Tmin")),
    fileName_(dict.lookup("fileName")), 

    YSolvent_
    (
        IOobject
        (
            IOobject::groupName("YSolvent", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar("one", dimless, 1)
    ),

    one_
    (
        IOobject
        (
            IOobject::groupName("one", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar("one", dimless, 1)
    )

{
    if (  Nspecie_.size()!=this->speciesNames_.size()
        ||speciemax_.size()!=this->speciesNames_.size()
        ||speciemin_.size()!=this->speciesNames_.size()
        ||fileName_.size()!=this->speciesNames_.size() )
    {
        FatalErrorInFunction
            << "Differing number of species and solubilities"
            << exit(FatalError);
    }

    //- Define indexes for multilinear interpolation
    if (  Nspecie_.size() == 1 )
    {
        // bilinear
        int count = 1;
        int tot_key = pow(2,Nspecie_.size()+2); 
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                // key [0,1,...,2^Nspecie )
                // values [0,1,...,Nspecie) 
                std::pair<int,std::vector<int>> key_values (tot_key-count,{i,j});
                yindex_.insert(key_values);
                count += 1;
            } 
        }
    } 
    else if (  Nspecie_.size() == 2 )
    {
        // quadrilinear
        int count = 1;
        int tot_key = pow(2,Nspecie_.size()+2); 
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<2;k++)
                {
                    for(int l=0;l<2;l++)
                    {
                        std::pair<int,std::vector<int>> key_values (tot_key-count,{i,j,k,l});
                        yindex_.insert(key_values);
                        count += 1;
                    }
                }
            } 
        }
    } 
    else
    {
        FatalErrorInFunction
            << "interpolation for more than 2 species to be defined yet!"
            << exit(FatalError);    
    }   
     
    //- Read from input file

    forAllConstIter(hashedWordList, this->speciesNames_, iter)
    {
        const label index = this->speciesNames_[*iter]; // specie number
       
        std::ifstream inFile(fileName_[index], std::ios::in | std::ios::binary);

        std::vector<double> dataTable(N_[0]);

        inFile.read(reinterpret_cast<char*>(&dataTable[0]), dataTable.size()*sizeof(dataTable[0]));

        int key = index;     
        std::pair<int,std::vector<double>> key_values (key,dataTable);
        saturation_.insert (key_values);

        inFile.close();  
    } 

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::ReadTable<Thermo, OtherThermo>::~ReadTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::interfaceCompositionModels::ReadTable<Thermo, OtherThermo>::update
(
    const volScalarField& Tf
)
{
    YSolvent_ = scalar(1);

    forAllConstIter(hashedWordList, this->speciesNames_, iter)
    {
        YSolvent_ -= Yf(*iter, Tf);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::ReadTable<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (this->speciesNames_.contains(speciesName))
    {

        const label index = this->speciesNames_[speciesName]; // specie number

        //- Define variables

        volScalarField p = this->thermo_.rhoThermo::p();
        volScalarField Tphase = this->thermo_.rhoThermo::T();
        volScalarField TotherPhase = this->otherThermo_.rhoThermo::T();
        volScalarField T = scalar(0.5)*(Tphase+TotherPhase); // try to use Tf
        PtrList<volScalarField> Yspecie(this->speciesNames_.size());

        //- Calculate total mass fractions of species

        // mass fraction of the phase (mass phase/total mass)
        volScalarField Yphase = ( this->pair_.phase1()/this->otherThermo_.rhoThermo::rho() )
                                 /(  this->pair_.phase2()/this->thermo_.rhoThermo::rho() 
                                   + this->pair_.phase1()/this->otherThermo_.rhoThermo::rho() );
        volScalarField YotherPhase = one_ - Yphase;

        forAllConstIter(hashedWordList, this->speciesNames_, iter)
        {
            const label indexList = this->speciesNames_[*iter];

            // mass fraction of specie in the phase (mass specie in the phase/mass phase)
            volScalarField YspeciePhase_phase = this->thermo_.composition().Y(*iter);
            volScalarField YspecieOtherPhase_otherPhase = this->otherThermo_.composition().Y(*iter);

            // total mass fraction of specie (total mass specie/total mass)
            Yspecie.set
	          (
	              indexList,
	              new volScalarField
                (
	                  IOobject
	                  (
	                      "Yspecie." + *iter,
	                      this->pair_.phase1().time().timeName(),
	                      this->pair_.phase1().mesh(),
	                      IOobject::NO_READ,
	                      IOobject::NO_WRITE
	                  ),
	                  YspeciePhase_phase*Yphase + YspecieOtherPhase_otherPhase*YotherPhase
                )
	          ); 
       
        }    

        //- Interpolate

        volScalarField one_p
        (
            IOobject
            (
                "one_p",
	              this->pair_.phase1().time().timeName(),
	              this->pair_.phase1().mesh(),
	              IOobject::NO_READ,
	              IOobject::NO_WRITE
            ),
            this->pair_.phase1().mesh(),
            dimensionedScalar("one_p", dimensionSet(1,-1,-2,0,0,0,0), 1)
        );

        volScalarField one_T
        (
            IOobject
            (
                "one_T",
	              this->pair_.phase1().time().timeName(),
	              this->pair_.phase1().mesh(),
	              IOobject::NO_READ,
	              IOobject::NO_WRITE
            ),
            this->pair_.phase1().mesh(),
            dimensionedScalar("one_T", dimensionSet(0,0,0,1,0,0,0), 1)
        );

        std::vector<int> yindex_0(yindex_.at(0));    

        // Lists of indexes
        PtrList<PtrList<volScalarField> > Index(yindex_0.size());
        // Lists of coefficients 
        PtrList<PtrList<volScalarField> > Alpha(yindex_0.size());

        // set pressure index [0]
        Index.set(0, new PtrList<volScalarField>(2));

        Index[0].set // set i   [0][0]
        (
            0,
            new volScalarField
            (
                IOobject
                (
                    "Index.pressure" ,
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                Np_[0]*(p - one_p*pmin_[0])/( one_p*pmax_[0] - one_p*pmin_[0] )
            )
        );

        // correct internalField
        forAll(Index[0][0],celli)
        {
            Index[0][0][celli] = floor(Index[0][0][celli]);
        }   

        Index[0].set // set (i+1)   [0][1]
        (
            1,
            new volScalarField
            (
                IOobject
                (
                    "Index.pressure",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                Index[0][0] + one_
            )
        );

        // set pressure coefficients [0]
        Alpha.set(0, new PtrList<volScalarField>(2));

        Alpha[0].set // set (1-alpha)   [0][0]
        (
            0,
            new volScalarField
            (
                IOobject
                (
                    "Alpha.pressure" ,
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                one_ - (Np_[0]*p - Np_[0]*one_p*pmin_[0] - one_p*(pmax_[0]-pmin_[0])*Index[0][0])/( one_p*(pmax_[0]-pmin_[0])) 
            )
        );

        Alpha[0].set // set alpha   [0][1]
        (
            1,
            new volScalarField
            (
                IOobject
                (
                    "Alpha.pressure",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                (Np_[0]*p - Np_[0]*one_p*pmin_[0] - one_p*(pmax_[0]-pmin_[0])*Index[0][0])/( one_p*(pmax_[0]-pmin_[0]))
            )
        );

        // set temperature index [1]
        Index.set(1, new PtrList<volScalarField>(2));

        Index[1].set // set i   [1][0]
        (
            0,
            new volScalarField
            (
                IOobject
                (
                    "Index.temperature" ,
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                NT_[0]*(T - one_T*Tmin_[0])/( one_T*Tmax_[0] - one_T*Tmin_[0] )
            )
        );

        // correct internalField
        forAll(Index[1][0],celli)
        {
            Index[1][0][celli] = floor(Index[1][0][celli]);
        }  

        Index[1].set // set (i+1)   [1][1]
        (
            1,
            new volScalarField
            (
                IOobject
                (
                    "Index.temperature",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                Index[1][0] + one_
            )
        );

        // set temperature coefficients [1]
        Alpha.set(1, new PtrList<volScalarField>(2));

        Alpha[1].set // set (1-alpha)   [1][0]
        (
            0,
            new volScalarField
            (
                IOobject
                (
                    "Alpha.temperature" ,
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                one_ - (NT_[0]*T - NT_[0]*one_T*Tmin_[0] - one_T*(Tmax_[0]-Tmin_[0])*Index[1][0])/( one_T*(Tmax_[0]-Tmin_[0])) 
            )
        );

        Alpha[1].set // set alpha   [1][1]
        (
            1,
            new volScalarField
            (
                IOobject
                (
                    "Alpha.temperature",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                (NT_[0]*T - NT_[0]*one_T*Tmin_[0] - one_T*(Tmax_[0]-Tmin_[0])*Index[1][0])/( one_T*(Tmax_[0]-Tmin_[0])) 
            )
        );

        if( this->speciesNames_.size() > 1 )
        {
            forAllConstIter(hashedWordList, this->speciesNames_, iter) // set coefficient for species
            {
                const label indexList = this->speciesNames_[*iter]; 

                // set species index [2:Nspecies-1]
                Index.set(indexList+2, new PtrList<volScalarField>(2));
 
                Index[indexList+2].set // set i   [indexList][0]
                (
                    0,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Index." + *iter,
	                          this->pair_.phase1().time().timeName(),
	                          this->pair_.phase1().mesh(),
	                          IOobject::NO_READ,
	                          IOobject::NO_WRITE
                        ),
                        Nspecie_[indexList]*( Yspecie[indexList] - one_*speciemin_[indexList])/( one_*speciemax_[indexList] - one_*speciemin_[indexList] )
                    )
                );

                // correct internalField
                forAll(Index[indexList+2][0],celli)
                {
                    Index[indexList+2][0][celli] = floor(Index[indexList+2][0][celli]);
                }  

                Index[indexList+2].set // set (i+1)   [indexList][1]
                (
                    1,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Index." + *iter,
	                          this->pair_.phase1().time().timeName(),
	                          this->pair_.phase1().mesh(),
	                          IOobject::NO_READ,
	                          IOobject::NO_WRITE
                        ),
                        Index[indexList+2][0] + one_
                    )
                );

                // set species coefficients [2:Nspecies-1]
                Alpha.set(indexList+2, new PtrList<volScalarField>(2));
 
                Alpha[indexList+2].set // set (1-alpha)   [indexList][0]
                (
                    0,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Alpha." + *iter,
	                          this->pair_.phase1().time().timeName(),
	                          this->pair_.phase1().mesh(),
	                          IOobject::NO_READ,
	                          IOobject::NO_WRITE
                        ),
                        one_ - (Nspecie_[indexList]*Yspecie[indexList] - Nspecie_[indexList]*one_*speciemin_[indexList] - (speciemax_[indexList]-speciemin_[indexList])*Index[indexList+2][0])/( (speciemax_[indexList]-speciemin_[indexList]))
                    )
                );

                Alpha[indexList+2].set // set alpha   [indexList][1]
                (
                    1,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Alpha." + *iter,
	                          this->pair_.phase1().time().timeName(),
	                          this->pair_.phase1().mesh(),
	                          IOobject::NO_READ,
	                          IOobject::NO_WRITE
                        ),
                        (Nspecie_[indexList]*Yspecie[indexList] - Nspecie_[indexList]*one_*speciemin_[indexList] - (speciemax_[indexList]-speciemin_[indexList])*Index[indexList+2][0])/( (speciemax_[indexList]-speciemin_[indexList]))
                    )
                );

            }
        }

        // interpolation

        volScalarField F_lookup
        (
            IOobject
            (
                "F_lookup",
	              this->pair_.phase1().time().timeName(),
	              this->pair_.phase1().mesh(),
	              IOobject::NO_READ,
	              IOobject::NO_WRITE
            ),
            this->pair_.phase1().mesh(),
            dimensionedScalar("F_lookup", dimensionSet(0,0,0,0,0,0,0), 0)
        );

        std::vector<int> Nvector(yindex_0.size()); // number of values for each variable
        Nvector[0] = Np_[0]+1;
        Nvector[1] = NT_[0]+1;
        for (int j=2; j< yindex_0.size(); j++)
        {
            Nvector[j] = Nspecie_[j-2]+1;
        }

        int key = index; // index for specie
        std::vector<double> saturation(saturation_.at(key)); // lookup table for specie     

        // loop over interpolation terms (from 0 to 2^Nspecie)
        for (int i=0;i<std::pow(2, yindex_0.size() );i++)
        {
             volScalarField Alpha_Pi
            (
                IOobject
                (
                    "Alpha_Pi",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                this->pair_.phase1().mesh(),
                dimensionedScalar("Alpha_Pi", dimensionSet(0,0,0,0,0,0,0), 1)
            );

            std::vector<int> yindex(yindex_.at(i)); 

            // loop over variables (p, T, species)  to calculate alpha Pi
            for (int j=0; j< yindex_0.size(); j++)
            {    
                Alpha_Pi = Alpha[j][yindex[j]]*Alpha_Pi;               
            } 

            // calculate single index for lookup table
            volScalarField K
            (
                IOobject
                (
                    "K",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                Index[yindex_0.size()-1][yindex[yindex_0.size()-1]] // initialize with value for last specie
            );

            for (int j=2; j<yindex_0.size()+1; j++)
            {
                int N_Pi = 1;
                for (int k=1; k< j; k++)
                { 
                    // calculate N Pi     
                    N_Pi = Nvector[yindex_0.size()-k]*N_Pi;
                } 
                K += Index[yindex_0.size()-j][yindex[yindex_0.size()-j]]*N_Pi;
            }
 
            // evaluate the ith interpolation term            
            volScalarField F_lookup_last 
            (
                IOobject
                (
                    "F_lookup_last",
	                  this->pair_.phase1().time().timeName(),
	                  this->pair_.phase1().mesh(),
	                  IOobject::NO_READ,
	                  IOobject::NO_WRITE
                ),
                this->pair_.phase1().mesh(),
                dimensionedScalar("F_lookup_last", dimensionSet(0,0,0,0,0,0,0), 1)
            );
            
            forAll(K,celli)
            { 
                int kindx = K[celli];
 
                F_lookup_last[celli] = Alpha_Pi[celli]*saturation[kindx];
            }

            // calculate iteratively the interpolation function for the specie               
            F_lookup += F_lookup_last;  
          
        }    

        return
            F_lookup*one_;

    }
    else
    {
        return
            YSolvent_
           *this->thermo_.composition().Y(speciesName);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::ReadTable<Thermo, OtherThermo>::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("YfPrime", this->pair_.name()),
                this->pair_.phase1().mesh().time().timeName(),
                this->pair_.phase1().mesh()
            ),
            this->pair_.phase1().mesh(),
            dimensionedScalar("zero", dimless/dimTemperature, 0)
        )
    );
}


// ************************************************************************* //
