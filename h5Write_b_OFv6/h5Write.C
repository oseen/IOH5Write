/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |		    2012-2014 Håkon Strandenes
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

#include "h5Write.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "scalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(h5Write, 0);
    addToRunTimeSelectionTable(functionObject, h5Write, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::h5Write::h5Write
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    name_(name)
{
    read(dict);

    // Calssify fields store in nFields_
    nFields_ = classifyFields();
    
    // Initialize file
    fileCreate(); 

    // Only do if some fields are to be written
    if (nFields_)
    {
        // Set length of cell numbers array of each proc
        nCells_.setSize(Pstream::nProcs());
        nPatchCellsPerProc_.setSize(Pstream::nProcs());

        // Write mesh and initial conditions
        meshWrite();
        write();
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::h5Write::~h5Write()
{
    // Close the HDF5 dataset
    fileClose();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::h5Write::read(const dictionary& dict)
{

    // Lookup in dictionary
    dict.lookup("objectNames") >> objectNames_;
    dict.lookup("cloudNames") >> cloudNames_;
    dict.lookup("cloudAttribs") >> cloudAttribs_;
    dict.lookup("writeInterval") >> writeInterval_;
    dict.lookup("patchNames") >> patchNames_;
    
    dict.lookup("suppressFieldDataWrite") >> suppressFieldDataWrite_;

    // Lookup chunk size if present
    chunkSize_ = dict.lookupOrDefault<label>("chunkSize", 0);

    // Set next write NOW
    nextWrite_ = 0;
    timeSteps_ = 0;

    // Print info to terminal
    int writeprec = sizeof(ioScalar);
    Info<< type() << " " << name() << ":" << endl
        << "  Compiled with " << writeprec << " bytes precision." << endl
        << "  Writing every " << writeInterval_ << " iterations:"
        << endl
        << "   ";

 
    // Do a basic check to see if the objectNames_ is accessible
    // print found objectNames_ (U, p, ...)
    forAll(objectNames_, i)
    {
        if (obr_.foundObject<regIOobject>(objectNames_[i]))
        {
            Info<< " " << objectNames_[i];
        }
        else
        {
            WarningIn
            (
                "Foam::writeRegisteredObject::read(const dictionary&)"
            )   << "Object " << objectNames_[i] << " not found in "
                << "database. Available objects:" << nl << obr_.sortedToc()
                << endl;
        }
    }


    // Do a basic check to see if patchNames_ is accessible
    // print found patchNames_ (inlet, outlet ...)
    Info<< "\n  Include patch/patches for writing: "
            << endl
            << "   ";

    forAll(patchNames_, i)
    {
            Info<< " " << patchNames_[i];
    }


    // Also print the cloud names
    forAll(cloudNames_, i)
    {
      Info<< " " << cloudNames_[i];
    }
    Info<< endl << endl;


    if (suppressFieldDataWrite_ == true){Info << "  Not included for writing: internal field data!\n" << endl;}


    // Check if writeInterval is a positive number
    if (writeInterval_ <= 0)
    {
        FatalIOErrorIn("h5Write::read(const dictionary&)", dict)
            << "Illegal value for writeInterval " << writeInterval_
            << ". It should be > 0."
            << exit(FatalIOError);
    }

        Info << endl;

    return true;
}


bool Foam::functionObjects::h5Write::execute()
{
    return true;
}


bool Foam::functionObjects::h5Write::end()
{
    return true;
}


bool Foam::functionObjects::h5Write::write()
{

    // Check if we are going to write
    if ( timeSteps_ == nextWrite_ )
    {
        // Write info to terminal
        Info << "Writing HDF5 data for time " << obr_.time().timeName() << endl;
        Info << endl;


	// Only write field data if fields are specified
        if (nFields_)
        {
            // Re-write mesh if dynamic
            if (mesh_.changing())
            {
                meshWrite();
            }

            // Write field data
            if (suppressFieldDataWrite_ == true){Info << "Suppress writing internal field data to HDF5!\n" << endl;}
	    else { fieldWrite(); }

            // Write patch data
            patchWrite();
        }


	// Only write cloud data if any clouds is specified
	if (cloudNames_.size())
        {
            //cloudWrite();
        }

	// Flush file cache (in case application crashes before it is finished)
        H5Fflush(fileID_, H5F_SCOPE_GLOBAL);
    
	// Calculate time of next write
        nextWrite_ = timeSteps_ + writeInterval_;
    }

    // Update counter
    timeSteps_++;



    return true;
}


// ************************************************************************* //
