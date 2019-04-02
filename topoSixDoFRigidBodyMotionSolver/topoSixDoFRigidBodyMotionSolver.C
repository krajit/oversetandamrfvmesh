/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "topoSixDoFRigidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSixDoFRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        topoSixDoFRigidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSixDoFRigidBodyMotionSolver::topoSixDoFRigidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    sixDoFRigidBodyMotionSolver(mesh,dict)
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




void Foam::topoSixDoFRigidBodyMotionSolver::solve()
{
    // reset points
    points0() = mesh().points();

    sixDoFRigidBodyMotionSolver::solve();

}


bool Foam::topoSixDoFRigidBodyMotionSolver::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    return sixDoFRigidBodyMotionSolver::writeObject(fmt, ver, cmp, valid);
}


bool Foam::topoSixDoFRigidBodyMotionSolver::read()
{
    return sixDoFRigidBodyMotionSolver::read();
}


// ************************************************************************* //