/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR,AFFILIATION
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = b076c47cf5973be90ce21ceebfc3b51ab3b96e24
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void cricketBall_b076c47cf5973be90ce21ceebfc3b51ab3b96e24(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    cricketBallFixedValueFvPatchScalarField
);


const char* const cricketBallFixedValueFvPatchScalarField::SHA1sum =
    "b076c47cf5973be90ce21ceebfc3b51ab3b96e24";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cricketBallFixedValueFvPatchScalarField::
cricketBallFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24"
            " from patch/DimensionedField\n";
    }
}


cricketBallFixedValueFvPatchScalarField::
cricketBallFixedValueFvPatchScalarField
(
    const cricketBallFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24"
            " from patch/DimensionedField/mapper\n";
    }
}


cricketBallFixedValueFvPatchScalarField::
cricketBallFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24"
            " from patch/dictionary\n";
    }
}


cricketBallFixedValueFvPatchScalarField::
cricketBallFixedValueFvPatchScalarField
(
    const cricketBallFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24"
            " as copy\n";
    }
}


cricketBallFixedValueFvPatchScalarField::
cricketBallFixedValueFvPatchScalarField
(
    const cricketBallFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cricketBallFixedValueFvPatchScalarField::
~cricketBallFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cricketBallFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs cricketBall sha1: b076c47cf5973be90ce21ceebfc3b51ab3b96e24\n";
    }

//{{{ begin code
    #line 0 "/home/ajit/OpenFOAM/ajit-v1806/oversetAndAMRFvMesh/testCase/ballSwing/0/nut.boundaryField.ballIn"
const label topCell = 28;

            vectorField faceCenters(this->patch().Cf());

			vector currentCOM = sum(faceCenters)/faceCenters.size();

			scalarField nut(faceCenters.size(), 0.0);

			scalar smoothNut = 0.0;
			scalar seamNut   = 0.1;
			scalar roughNut  = 0.01;
			scalar PI = 3.14159265;


			scalar seamWidth = 20;   // in angle (deg) subtended at center

			
			forAll(faceCenters, faceI) {

				vector ot = faceCenters[topCell] - currentCOM; 
				vector op = faceCenters[faceI] - currentCOM; 

				if ((ot&op)/(mag(ot)*mag(op)) >= cos( (90.0 - seamWidth/2)*PI/180))
				{
					// side on top of seam
					// keep it smooth
					nut[faceI] = smoothNut;
				} else if ((ot&op)/(mag(ot)*mag(op)) >= cos( (90.0 + seamWidth/2)*PI/180))
				{
					// seam
					
					nut[faceI] = seamNut;
				} else
				{
					// remaining sides are rought
					nut[faceI] = roughNut;
				}


			}
            operator==(nut);
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

