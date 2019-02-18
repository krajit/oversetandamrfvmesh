/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "oversetAndAMRFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetAndAMRFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, oversetAndAMRFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetAndAMRFvMesh::oversetAndAMRFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    meshOverset_(io),
    meshDynamicRefine_(io)
{}


Foam::oversetAndAMRFvMesh::oversetAndAMRFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    dynamicFvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
        meshOverset_(io),
        meshDynamicRefine_(io)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::oversetAndAMRFvMesh::update()
{

    //bool oversetMeshChange = meshOverset_.update();

    // TODO: add mesh rifine later
        bool dynamicRefineMeshChange = meshDynamicRefine_.update();

    return dynamicRefineMeshChange;
//    return oversetMeshChange + dynamicRefineMeshChange;
}


// ************************************************************************* //
