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
    if (magSqr(diskArea_) <= VSMALL)
    {
        FatalErrorIn("Foam::fv::gustSource::checkData()")
           << "diskArea is approximately zero"
           << exit(FatalIOError);
    }
    if (Cp_ <= VSMALL || Ct_ <= VSMALL)
    {
        FatalErrorIn("Foam::fv::gustSource::checkData()")
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorIn("Foam::fv::gustSource::checkData()")
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }
    if (returnReduce(upstreamCellId_, maxOp<label>()) == -1)
    {
        FatalErrorIn("Foam::fv::gustSource::checkData()")
           << "upstream location " << upstreamPoint_  << " not found in mesh"
           << exit(FatalIOError);
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
    diskDir_(coeffs_.lookup("diskDir")),
    Cp_(readScalar(coeffs_.lookup("Cp"))),
    Ct_(readScalar(coeffs_.lookup("Ct"))),
    diskArea_(readScalar(coeffs_.lookup("diskArea"))),
    upstreamPoint_(coeffs_.lookup("upstreamPoint")),
    upstreamCellId_(-1)
{
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating gusty zone: "
        << this->name() << endl;

    upstreamCellId_ = mesh.findCell(upstreamPoint_);

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
        coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);
        coeffs_.readIfPresent("diskArea", diskArea_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //