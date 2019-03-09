// Member Functions

//- Direct access to the refinement engine
const hexRef8& meshCutter() const
{
    return meshCutter_;
}

//- Cells which should not be refined/unrefined
const bitSet& protectedCell() const
{
    return protectedCell_;
}


//- Cells which should not be refined/unrefined
bitSet& protectedCell()
{
    return protectedCell_;
}

//- Update the mesh for both mesh motion and topology change
// this will be merged with other
// virtual bool update();


// Writing

//- Write using given format, version and compression
// AK: was virtual earlier
//virtual bool writeObject

// AK: this method is creating problem. Work on it later
// bool writeObject
// (
//     IOstream::streamFormat fmt,
//     IOstream::versionNumber ver,
//     IOstream::compressionType cmp,
//     const bool valid
// ) const;
