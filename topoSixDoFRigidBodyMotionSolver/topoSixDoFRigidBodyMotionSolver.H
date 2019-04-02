/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#ifndef topoSixDoFRigidBodyMotionSolver_H
#define topoSixDoFRigidBodyMotionSolver_H

#include "sixDoFRigidBodyMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class topoSixDoFRigidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class topoSixDoFRigidBodyMotionSolver
:
    public sixDoFRigidBodyMotionSolver
{

    // Private Member Functions

        //- No copy construct
        topoSixDoFRigidBodyMotionSolver
        (
            const topoSixDoFRigidBodyMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const topoSixDoFRigidBodyMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("topoSixDoFRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        topoSixDoFRigidBodyMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~topoSixDoFRigidBodyMotionSolver() = default;


    // Member Functions

        //- Solve for motion
        virtual void solve();

        //- Write state using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool valid
        ) const;

        //- Read dynamicMeshDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
