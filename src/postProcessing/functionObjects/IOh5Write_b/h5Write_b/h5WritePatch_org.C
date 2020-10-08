/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2012-2014 Håkon Strandenes
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::h5Write::patchWrite()
{
    Info<< "  h5Write::patchWrite:"  << endl;
    patchWriteScalar();
    patchWriteVector();
}

//--------------------------------------------------------------------//


void Foam::h5Write::patchWriteScalar()
{
    Info<< "  h5Write::patchWriteScalar:"  << endl;

    forAll(scalarFields_, sfieldI) // ### FLO: Loop over all scalar fields, iterator = fieldI
    {
        // ### FLO: Info which scalar fields are to be written
        Info<< "    patchWriteScalar: " << scalarFields_[sfieldI] << endl;

	forAll(patchNames_, patchI) //loop over all patches with iterator patchI
        {
		Info <<"  Patch for loop: " << patchNames_[patchI] << endl;
		label patchID = mesh_.boundaryMesh().findPatchID( patchNames_[patchI]  ); // get patchID from each patch
		Info <<"  patchID: " << patchID << endl;
		
		// define scalarField and patch
    		const volScalarField& myScalarField = mesh_.lookupObject<volScalarField>( scalarFields_[sfieldI] );
    		const fvPatch& myPatch = mesh_.boundary()[patchID];

		// Initialize a plain continous array for the data
	        ioScalar* scalarPatchDataArray;
	        scalarPatchDataArray = new ioScalar[myPatch.size()]; // ### Flo: initialize array with field.size() components
        	Info << "Patch.size() = " << myPatch.size() << endl;

    		forAll(myPatch,faceI)
    		{
			Info <<"\t Patch[" << patchID << "]: " << patchNames_[patchI] << ", Skalarfeld: " << scalarFields_[sfieldI] <<
                        "   Face[" << faceI << "] = " << myScalarField.boundaryField()[patchID][faceI] << endl;
			
			// fill array scalarData with patch information
		        scalarPatchDataArray[faceI] = myScalarField.boundaryField()[patchID][faceI]; 
			Info <<"\t scalarPatchDataArray[" << faceI << "] = " << scalarPatchDataArray[faceI] << endl;

    		}


		Info << "Pstream::myProcNo() = " << Pstream::myProcNo() << endl;

/*

	        // Create the different datasets (needs to be done collectively)
	        char datasetName[80];
        	hsize_t dimsf[1];       // define dataset dimensions
	        hid_t fileSpace;        // fileSpace handle
       		hid_t dsetID;           // handle
      		hid_t plistID;          // handle
        	hid_t plistDCreate;     // handle

        	// loop through cells on patch of each processor with iterator = proc = processor number (0,1,2,..)
		forAll(nPatchCellsPerProc_, proc)
     	 	{
            		// Create the dataspace for the dataset with number of patch cells on the current processor & patch
            		dimsf[0] = nPatchCellsPerProc_[proc][patchI];

            		fileSpace = H5Screate_simple(1, dimsf, NULL); // H5Screate_simple(RANK, dimsf, NULL);

            		// Set property to create parent groups as neccesary
            		plistID = H5Pcreate(H5P_LINK_CREATE);
            		H5Pset_create_intermediate_group(plistID, 1);

            		// Set chunking, compression and other HDF5 dataset properties
            		plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
            		dsetSetProps(1, sizeof(ioScalar), nPatchCellsPerProc_[proc][patchI], plistDCreate);

            		// Create the dataset for points
            		// ### Flo: see h5dump_H --> GROUP "outlet" { GROUP "0.025" { GROUP "processor0" { DATASET "Lambda2" }}}}}}

            		sprintf
                	(
                    		datasetName,
                    		"%s/%s/processor%i/%s",
				patchNames_[patchI].c_str(),
                    		mesh_.time().timeName().c_str(),
                    		proc,
                    		scalarFields_[sfieldI].c_str()

                    		// %s --> mesh_.time().timeName().c_str(); %i --> proc; %s --> scalarFields_[fieldI].c_str()
                	);

            		dsetID = H5Dcreate2
                	(
                    		fileID_,
                    		datasetName,
                    		H5T_SCALAR,
                    		fileSpace,
                    		plistID,
                    		plistDCreate,
                    		H5P_DEFAULT
                	);
            		
			H5Dclose(dsetID);
            		H5Pclose(plistID);
            		H5Pclose(plistDCreate);
            		H5Sclose(fileSpace);

        	}// here does loop over all nPatchCellsPerProc_ end

		// Open correct dataset for this process
        	sprintf
        	(
                	datasetName,
                	"%s/%s/processor%i/%s",
			patchNames_[patchI].c_str(),
                	mesh_.time().timeName().c_str(),
                	Pstream::myProcNo(),
                	scalarFields_[sfieldI].c_str()
        	);
        
		dsetID = H5Dopen2(fileID_, datasetName, H5P_DEFAULT);

        	// Create property list for collective dataset write.
        	plistID = H5Pcreate(H5P_DATASET_XFER);
        	H5Pset_dxpl_mpio(plistID, H5_XFER_MODE);

        	// Do the actual write
        	H5Dwrite
            	(
                	dsetID,
                	H5T_SCALAR,
                	H5S_ALL,
                	H5S_ALL,
                	plistID,
                	scalarPatchDataArray
            	);

        	// Close/release resources.
        	H5Pclose(plistID);
        	H5Dclose(dsetID);

		// Release memory
        	delete [] scalarPatchDataArray;
*/


	} // end of forAll(patchNames_, patchI)

    } // end of forAll(scalarFields_, sfieldI)

} // end of void Foam::h5Write::patchWriteScalar()

//-------------------------------------------------------------------//
// ### Flo: The same for all vectors (velocity, ...)

void Foam::h5Write::patchWriteVector()
{
    Info<< "  h5Write::patchWriteVector:"  << endl;

    forAll(vectorFields_, vfieldI) // ### FLO: Loop over all scalar fields, iterator = fieldI
    {
        // ### FLO: Info which vector fields are to be written
        Info<< "    patchWriteVector: " << vectorFields_[vfieldI] << endl;

        forAll(patchNames_, patchI) //loop over all patches with iterator patchI
        {
                Info <<"  Patch for loop: " << patchNames_[patchI] << endl;
                label patchID = mesh_.boundaryMesh().findPatchID( patchNames_[patchI]  ); // get patchID from each patch
                Info <<"  patchID: " << patchID << endl;

                const volVectorField& myVectorField = mesh_.lookupObject<volVectorField>( vectorFields_[vfieldI] );
                const fvPatch& myPatch = mesh_.boundary()[patchID];

                forAll(myPatch,faceI)
                {
                        //Foam::label faceCellI = myPatch.faceCells()[faceI];
                        //Info << myScalarField[faceCellI] << endl;
                        Info <<"\t Patch[" << patchID << "]: " << patchNames_[patchI] << ", Vectorfeld: " << vectorFields_[vfieldI] <<
                        "   Face[" << faceI << "] = " << myVectorField.boundaryField()[patchID][faceI] << endl;
		}

	}

     }

}


// ************************************************************************* //
