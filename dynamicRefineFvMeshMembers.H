protected:
//public:
//- Mesh cutting engine
hexRef8 meshCutter_;

//- Dump cellLevel for post-processing
bool dumpLevel_;

//- Fluxes to map
HashTable<word> correctFluxes_;

//- Number of refinement/unrefinement steps done so far.
label nRefinementIterations_;

//- Protected cells (usually since not hexes)
bitSet protectedCell_;


// Protected Member Functions

    //- Count set/unset elements in packedlist.
    static label count(const bitSet&, const unsigned int);


    //- Calculate cells that cannot be refined since would trigger
//  refinement of protectedCell_ (since 2:1 refinement cascade)
void calculateProtectedCells(bitSet& unrefineableCell) const;

//- Read the projection parameters from dictionary
void readDict();


//- Refine cells. Update mesh and fields.
autoPtr<mapPolyMesh> refine(const labelList&);

//- Unrefine cells. Gets passed in centre points of cells to combine.
autoPtr<mapPolyMesh> unrefine(const labelList&);


// Selection of cells to un/refine

    //- Calculates approximate value for refinement level so
    //  we don't go above maxCell
    scalar getRefineLevel
    (
        const label maxCells,
        const label maxRefinement,
        const scalar refineLevel,
        const scalarField&
    ) const;

    //- Get per cell max of connected point
    scalarField maxPointField(const scalarField&) const;

    //- Get point max of connected cell
    scalarField maxCellField(const volScalarField&) const;

    scalarField cellToPoint(const scalarField& vFld) const;

    scalarField error
    (
        const scalarField& fld,
        const scalar minLevel,
        const scalar maxLevel
    ) const;


    //- Select candidate cells for refinement
    //virtual void selectRefineCandidates 
    //AK: was virtual in the dynamicRefineFvMesh
    void selectRefineCandidates
    (
        const scalar lowerRefineLevel,
        const scalar upperRefineLevel,
        const scalarField& vFld,
        bitSet& candidateCell
    ) const;

    //- Subset candidate cells for refinement
    //virtual labelList selectRefineCells
    //AK: was virtual in the dynamicRefineFvMesh
    labelList selectRefineCells
    (
        const label maxCells,
        const label maxRefinement,
        const bitSet& candidateCell
    ) const;

    //- Select points that can be unrefined.
    //virtual labelList selectUnrefinePoints
    //AK: was virtual in the dynamicRefineFvMesh
    labelList selectUnrefinePoints
    (
        const scalar unrefineLevel,
        const bitSet& markedCell,
        const scalarField& pFld
    ) const;

    //- Extend markedCell with cell-face-cell.
    void extendMarkedCells(bitSet& markedCell) const;

    //- Check all cells have 8 anchor points
    void checkEightAnchorPoints
    (
        bitSet& protectedCell,
        label& nProtected
    ) const;