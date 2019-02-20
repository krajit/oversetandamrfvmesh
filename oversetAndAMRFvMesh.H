/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef oversetAndAMRFvMesh_H
#define oversetAndAMRFvMesh_H

#include "dynamicFvMesh.H"
#include "dynamicOversetFvMesh.H"
#include "dynamicRefineFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class oversetAndAMRFvMesh Declaration
\*---------------------------------------------------------------------------*/

class oversetAndAMRFvMesh
:
    public dynamicOversetFvMesh
{
        //- No copy construct
        oversetAndAMRFvMesh(const oversetAndAMRFvMesh&) = delete;

        //- No copy assignment
        void operator=(const oversetAndAMRFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("oversetAndAMRFvMesh");


    // Constructors

        //- Construct from IOobject
        oversetAndAMRFvMesh(const IOobject& io);

    //- Destructor
    ~oversetAndAMRFvMesh() = default;

        //- Dummy update function which does not change the mesh
     virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
