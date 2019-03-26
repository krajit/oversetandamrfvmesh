/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "refinementFieldFunction.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
defineTypeNameAndDebug(refinementFieldFunction, 0);

addToRunTimeSelectionTable(
    functionObject,
    refinementFieldFunction,
    dictionary);
} // namespace functionObjects
} // namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::refinementFieldFunction::writeFileHeader(Ostream &os) const
{
    writeHeader(os, "y+ ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::refinementFieldFunction::refinementFieldFunction(
    const word &name,
    const Time &runTime,
    const dictionary &dict)
    : fvMeshFunctionObject(name, runTime, dict),
      writeFile(obr_, name, typeName, dict),
      xBuffer_(readScalar(dict.lookup("xBuffer"))),
      yBuffer_(readScalar(dict.lookup("yBuffer"))),
      zBuffer_(readScalar(dict.lookup("zBuffer")))
{
    read(dict);

    writeFileHeader(file());

    volScalarField *refinementFieldPtr(
        new volScalarField(
            IOobject(
                "refinementField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            mesh_,
            dimensionedScalar(dimless, Zero)));

    mesh_.objectRegistry::store(refinementFieldPtr);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::refinementFieldFunction::~refinementFieldFunction()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::refinementFieldFunction::read(const dictionary &dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}

bool Foam::functionObjects::refinementFieldFunction::execute()
{

    volScalarField &refinementField =
        const_cast<volScalarField &>(
            lookupObject<volScalarField>("refinementField"));

    const labelIOList& zoneID = lookupObject<labelIOList>("zoneID");


    // get the bounding box containing the zoneID with labels 1
    const vectorField &cellCenters = mesh_.C();

    forAll(cellCenters, celli)
    {
        // reset
        refinementField[celli] = -1.0;

        // only loop over front mesh
        if (zoneID[celli] == 1)
        {
            refinementField[celli] = 2.0;

            scalar cellDistFrontToBack = GREAT;
            label backgroundNearestCell = 0;

            // now loop over background mesh to find the background cell closest to this cellj
            forAll(cellCenters, cellj)
            {
                // filter out the background mesh
                if (zoneID[cellj] == 0)
                {
                    if (magSqr(cellCenters[celli] - cellCenters[cellj]) <= cellDistFrontToBack)
                    {
                        backgroundNearestCell = cellj;
                        cellDistFrontToBack = magSqr(cellCenters[celli] - cellCenters[cellj]);
                    }
                }
            }


            refinementField[backgroundNearestCell] = 1.0;
        }
    }
    return true;
}

bool Foam::functionObjects::refinementFieldFunction::write()
{
    // const volScalarField& refinementField =
    //     obr_.lookupObject<volScalarField>(type());

    // Log << type() << " " << name() << " write:" << nl
    //     << "    writing field " << refinementField.name() << endl;

    // refinementField.write();

    // const volScalarField::Boundary& refinementFieldBf = refinementField.boundaryField();
    // const fvPatchList& patches = mesh_.boundary();

    // forAll(patches, patchi)
    // {
    //     const fvPatch& patch = patches[patchi];

    //     if (isA<wallFvPatch>(patch))
    //     {
    //         const scalarField& refinementFieldp = refinementFieldBf[patchi];

    //         const scalar minrefinementField = gMin(refinementFieldp);
    //         const scalar maxrefinementField = gMax(refinementFieldp);
    //         const scalar avgrefinementField = gAverage(refinementFieldp);

    //         if (Pstream::master())
    //         {
    //             Log << "    patch " << patch.name()
    //                 << " refinementField : min = " << minrefinementField << ", max = " << maxrefinementField
    //                 << ", average = " << avgrefinementField << nl;

    //             writeTime(file());
    //             file()
    //                 << token::TAB << patch.name()
    //                 << token::TAB << minrefinementField
    //                 << token::TAB << maxrefinementField
    //                 << token::TAB << avgrefinementField
    //                 << endl;
    //         }
    //     }
    // }

    return true;
}

// ************************************************************************* //
