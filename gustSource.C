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

\*---------------------------------------------------------------------------*/

#include "gustSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(gustSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        gustSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::gustSource::checkData() const
{
    if (gustAmplitudes_.size() != gustFrequencies_.size())
    {
        FatalErrorIn("Foam::fv::gustSource::checkData()")
           << "different numbers of amplitudes and frequencies specified: "
           << "w_g= " << gustAmplitudes_ << " "
           << "omega_g= " << gustFrequencies_ << " "
           << exit(FatalIOError);
    }
    if (gustAmplitudes_.size() == 0)
    {
        Info<< "Note: No amplitudes/frequencies specified" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::gustSource::gustSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    gustAmplitudes_(readList<scalar>(coeffs_.lookup("amplitude"))),
    gustFrequencies_(readList<scalar>(coeffs_.lookup("frequency")))
{
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating gusty zone: "
        << this->name() << endl;

    checkData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::gustSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addGustMomenta
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U
        );
    }
}


void Foam::fv::gustSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addGustMomenta
        (
            Usource,
            cells_,
            cellsV,
            rho,
            U
        );
    }
}


void Foam::fv::gustSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::gustSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("amplitude", gustAmplitudes_);
        coeffs_.readIfPresent("frequency", gustFrequencies_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
