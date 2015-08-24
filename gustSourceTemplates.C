/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "volFields.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::gustSource::addGustMomenta
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    scalar A(0.0);
    scalar omg(0.0);
    vector bvec(0.0,0.0,0.0);
    vector rvec(0.0,0.0,0.0);
    scalar expfn(0.0);

    scalar t = mesh().time().timeOutputValue();
    vectorField cc(mesh().C().internalField());

    forAll(gustFrequencies_, mode)
    {
        A = gustAmplitudes_[mode];
        omg = gustFrequencies_[mode];

        forAll(cells, i)
        {
            rvec = cc[cells[i]] - r0_;
            bvec = rvec - (rvec & avec_)*avec_;
            expfn = Foam::exp( -(bvec & bvec)/(2*sourceRadius_*sourceRadius_) );
            //Usource[cells[i]] += A*sin(omg*t) * gustDirection_;
            Usource[cells[i]] += A*sin(omg*t) * gustDirection_ * expfn;
        }
    }
}


// ************************************************************************* //
