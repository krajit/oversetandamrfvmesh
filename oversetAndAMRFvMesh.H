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
//#include "dynamicRefine2DFvMesh.H"


#include "hexRef2D.H"
#include "bitSet.H"
#include "PackedBoolList.H"
#include "Switch.H"



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
    // protected members copied for dynamicRefine2DFvMesh
  protected:
    //- Mesh cutting engine
    hexRef2D meshCutter_;

    //- Dump cellLevel for postprocessing
    Switch dumpLevel_;

    //- Fluxes to map
    List<Pair<word>> correctFluxes_;

    //- Number of refinement/unrefinement steps done so far.
    label nRefinementIterations_;

    //- Protected cells (usually since not hexes)
    PackedBoolList protectedCell_;

    // Private Member Functions

    //- Count set/unset elements in packedlist.
    static label count(const PackedBoolList &, const unsigned int);

    //- Calculate cells that cannot be refined since would trigger
    //  refinement of protectedCell_ (since 2:1 refinement cascade)
    void calculateProtectedCells(PackedBoolList &unrefineableCell) const;

    //- Read the projection parameters from dictionary
    void readDict();

    //- Refine cells. Update mesh and fields.
    autoPtr<mapPolyMesh> refine(const labelList &);

    //- Unrefine cells. Gets passed in centre points of cells to combine.
    autoPtr<mapPolyMesh> unrefine(const labelList &);

    // Selection of cells to un/refine

    //- Calculates approximate value for refinement level so
    //  we don't go above maxCell
    scalar getRefineLevel(
        const label maxCells,
        const label maxRefinement,
        const scalar refineLevel,
        const scalarField &) const;

    //- Get per cell max of connected point
    scalarField maxPointField(const scalarField &) const;

    //- Get point min of connected cell
    scalarField minCellField(const volScalarField &) const;

    scalarField cellToPoint(const scalarField &vFld) const;

    scalarField error(
        const scalarField &fld,
        const scalar minLevel,
        const scalar maxLevel) const;

    //- Select candidate cells for refinement
    virtual void selectRefineCandidates(
        const scalar lowerRefineLevel,
        const scalar upperRefineLevel,
        const scalarField &vFld,
        PackedBoolList &candidateCell) const;

    //- Subset candidate cells for refinement
    virtual labelList selectRefineCells(
        const label maxCells,
        const label maxRefinement,
        const PackedBoolList &candidateCell) const;

    //- Select edges that can be unrefined (2D case).
    virtual labelList selectUnrefineEdges(
        const scalar unrefineLevel,
        const PackedBoolList &markedCell,
        const scalarField &pFld) const;

    //- Extend markedCell with cell-face-cell.
    void extendMarkedCells(PackedBoolList &markedCell) const;


    // end of copying of members from dynamicRefine2dFvMesh



private:
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

    //-------------------------------------------------
    // begin dynamicRefine2DFvMesh function declaration
    //-------------------------------------------------
   //- Direct access to the refinement engine
    const hexRef2D &meshCutter() const
    {
        return meshCutter_;
    }

    //- Cells which should not be refined/unrefined
    const PackedBoolList &protectedCell() const
    {
        return protectedCell_;
    }

    //- Cells which should not be refined/unrefined
    PackedBoolList &protectedCell()
    {
        return protectedCell_;
    }

    //- Update the mesh for both mesh motion and topology change
    // virtual bool update();
    // declared separately // AK: Mar 19, 2019

    // Writing

    //- Write using given format, version and compression
    // virtual bool writeObject(
    //     IOstream::streamFormat fmt,
    //     IOstream::versionNumber ver,
    //     IOstream::compressionType cmp,
    //     const bool valid) const;
    //-------------------------------------------------
    // end dynamicRefine2DFvMesh function declaration
    //-------------------------------------------------

    bool refineUpdate();

    virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
