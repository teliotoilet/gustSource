/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::fv::gustSource

Description
    Provide a sinusoidal gust through perturbations to the momentum equation.
    Perturbations are applied uniformly along a specified line.

    Created from actuationDiskSource

SourceFiles
    gustSource.C
    gustSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef gustSource_H
#define gustSource_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class gustSource Declaration
\*---------------------------------------------------------------------------*/

class gustSource
:
    public option
{

protected:

    // Protected data
        vector      gustDirection_;
        scalarList  gustAmplitudes_;
        scalarList  gustFrequencies_;

        vector      r0_;
        vector      avec_;
        scalar      sourceRadius_;

private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addGustMomenta
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        gustSource(const gustSource&);

        //- Disallow default bitwise assignment
        void operator=(const gustSource&);


public:

    //- Runtime type information
    TypeName("gustSource");


    // Constructors

        //- Construct from components
        gustSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~gustSource()
    {}


    // Member Functions

        // Access


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gustSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
