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
    public dynamicFvMesh
{
    // Private Member Functions


    // overset mesh member
    dynamicOversetFvMesh meshOverset_;

    // dynamic refine member
    dynamicRefineFvMesh meshDynamicRefine_;



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

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        oversetAndAMRFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );


    //- Destructor
    ~oversetAndAMRFvMesh() = default;


    // Member Functions

        //- Is mesh dynamic
        virtual bool dynamic() const
        {
            return false;
        }

        //- Dummy update function which does not change the mesh
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
