/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::ReadTable

Description
    Reads from tables the saturation surface and interpolates on pressure,
    temperature and total mass fractions of volatile components. Tables have
    to be provided in binary files created with the following loop order: 
    count = 0
    for(i=0;i<Np;i++) // pressure
      for(j=0;j<Np;j++)  // temperature
        for(k=0;k<Np;k++)  // specie 1
          for(l=0;l<Np;l++) // specie 2
            ...
            table[count] = ...
            count += 1; 

SourceFiles
    ReadTable.C

\*---------------------------------------------------------------------------*/

#ifndef ReadTable_H
#define ReadTable_H

#include "InterfaceCompositionModel.H"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class ReadTable Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class ReadTable
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
private:

    // Private data

        // Table info
        const scalarList N_; // tot number of values 
        const scalarList Nspecie_; // number of intervals for specie
        const scalarList speciemax_,speciemin_;
        const scalarList Np_; // number of intervals for p
        const scalarList pmax_,pmin_;
        const scalarList NT_; // number of intervals for T
        const scalarList Tmax_,Tmin_;
        fileNameList fileName_;

        //- The remaining solvent species fraction
        volScalarField YSolvent_;
        
        volScalarField one_;

        // indexes for multilinear interpolation
        std::map<int,std::vector<int>> yindex_;

        //- Hash map for saturation surface
        std::map<int,std::vector<double>> saturation_;

public:

    //- Runtime type information
    TypeName("ReadTable");

    // Constructors

        //- Construct from components
        ReadTable
        (
            const dictionary& dict,
            const phasePair& pair
        );

    //- Destructor
    virtual ~ReadTable();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReadTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
