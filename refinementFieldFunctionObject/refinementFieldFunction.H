/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
-------------------------------------------------------------------------------

Class
    Foam::functionObjects::refinementFieldFunction

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for turbulence models.

    The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

Usage
    Example of function object specification to calculate the y+ (LES):
    \verbatim
    refinementFieldFunction1
    {
        type        refinementFieldFunction;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: refinementFieldFunction        | yes         |
        result       | Name of y+ field        | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    refinementFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_refinementFieldFunction_H
#define functionObjects_refinementFieldFunction_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class refinementFieldFunction Declaration
\*---------------------------------------------------------------------------*/

class refinementFieldFunction
    : public fvMeshFunctionObject,
      public writeFile
{

    // private members
    scalar xBuffer_;
    scalar yBuffer_;
    scalar zBuffer_;

    // Private Member Functions

    //- File header information
    virtual void writeFileHeader(Ostream &os) const;

    //- No copy construct
    refinementFieldFunction(const refinementFieldFunction &) = delete;

    //- No copy assignment
    void operator=(const refinementFieldFunction &) = delete;

  public:
    //- Runtime type information
    TypeName("refinementFieldFunction");

    // Constructors

    //- Construct from Time and dictionary
    refinementFieldFunction(
        const word &name,
        const Time &runTime,
        const dictionary &dict);

    //- Destructor
    virtual ~refinementFieldFunction();

    // Member Functions

    //- Read the refinementFieldFunction data
    virtual bool read(const dictionary &);

    //- Calculate the refinementFieldFunction field
    virtual bool execute();

    //- Write the refinementFieldFunction field
    virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
