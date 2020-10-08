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
//    Info<< "###  h5Write::patchWrite:"  << endl;
    patchWriteScalar();
    patchWriteVector();
}

//--------------------------------------------------------------------//


void Foam::h5Write::patchWriteScalar()
{
//    Info << "###  h5Write::patchWriteScalar:"  << endl;
//    Info << "  Master-Processor for writing: Processor[" << Pstream::myProcNo() << "] \n" << endl;

    forAll(patchNames_, patchI) //loop over all patches with iterator patchI
    {
//	Info <<"  Patch for loop: " << patchNames_[patchI] << endl;
	label patchID = mesh_.boundaryMesh().findPatchID( patchNames_[patchI]  ); // get patchID from each patch
//	Info <<"  patchID: " << patchID << endl;
    	const fvPatch& myPatch = mesh_.boundary()[patchID]; // define patch

	// Initialize a plain continous array for the data to write to hdf5 later
        ioScalar* scalarPatchDataArray;
        scalarPatchDataArray = new ioScalar[myPatch.size()]; //  initialize array with patch.size() components
//     	Info << "  Patch.size() on processor[" << Pstream::myProcNo() << "] = " << myPatch.size() << endl;

    	
	forAll(scalarFields_, sfieldI) // Loop over all scalar fields, iterator = fieldI
        {
		// define scalarField 
    		const volScalarField& myScalarField = mesh_.lookupObject<volScalarField>( scalarFields_[sfieldI] );

    		forAll(myPatch,faceI)
    		{
//			Info <<"\t Patch[" << patchID << "]: " << patchNames_[patchI] << ", Skalarfeld: " << scalarFields_[sfieldI] <<
//                        "   Face[" << faceI << "] = " << myScalarField.boundaryField()[patchID][faceI] << endl;
			
			// fill array scalarData with patch information
		        scalarPatchDataArray[faceI] = myScalarField.boundaryField()[patchID][faceI]; 
//			Info <<"\t scalarPatchDataArray[" << faceI << "] = " << scalarPatchDataArray[faceI] << endl;
    		}

//		Info << endl;
	        
		// Create the different datasets (needs to be done collectively)
	        char datasetName[80];
        	hsize_t dimsf[1];       // define dataset dimensions
	        hid_t fileSpace;        // fileSpace handle
       		hid_t dsetID;           // dataset handle ?
      		hid_t plistID;          // handle
        	hid_t plistDCreate;     // handle

        	// loop through cells on patch of each processor with iterator = proc = processor number (0,1,2,..)
		forAll(nPatchCellsPerProc_, proc)
     	 	{
            		// Create the dataspace for the dataset with number of patch cells on the current processor & patch
            		dimsf[0] = nPatchCellsPerProc_[proc][patchI];

            		fileSpace = H5Screate_simple(1, dimsf, NULL); // H5Screate_simple(RANK, dimsf, NULL); H5S -> dataSpace

            		// Set property to create parent groups as neccesary
            		plistID = H5Pcreate(H5P_LINK_CREATE); // H5P --> Properties
            		H5Pset_create_intermediate_group(plistID, 1);

            		// Set chunking, compression and other HDF5 dataset properties
            		plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
			// set properties for plistDCreate
            		dsetSetProps(1, sizeof(ioScalar), nPatchCellsPerProc_[proc][patchI], plistDCreate);

            		// Create the dataset for points
            		// h5dump_H --> GROUP "PATCHFIELDS" { GROUP "0.025" { GROUP "inlet" { GROUP "processor0" { DATASET "Lambda2" }}}}}}

            		sprintf
                	(
                    		datasetName,
                    		"PATCHFIELDS/%s/%s/processor%i/%s",
                    		mesh_.time().timeName().c_str(),
				patchNames_[patchI].c_str(),
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
            		
			H5Dclose(dsetID); // H5Dclose: discard data object
            		H5Pclose(plistID); // H5Pclose: discard property object
            		H5Pclose(plistDCreate); // H5Pclose:discard property object
            		H5Sclose(fileSpace); // H5Sclose: discard Space object

        	}// here does loop over all nPatchCellsPerProc_ end

    		if ( nPatchCellsPerProc_[Pstream::myProcNo()][patchI] != 0 )
    		{
			// Open correct dataset for this process
	        	sprintf
        		(
                		datasetName,
                		"PATCHFIELDS/%s/%s/processor%i/%s",
                		mesh_.time().timeName().c_str(),
				patchNames_[patchI].c_str(),
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

    		} // end of if loop

	} // end of forAll(scalarFields_, sfieldI) 

	// Release memory
        delete [] scalarPatchDataArray;

    } // end of forAll(patchNames_, patchI)

} // end of void Foam::h5Write::patchWriteScalar()

//-------------------------------------------------------------------//
// ### Flo: The same for all vectors (velocity, ...)

void Foam::h5Write::patchWriteVector()
{
//    Info<< "  ### h5Write::patchWriteVector:"  << endl;

    forAll(patchNames_, patchI) //loop over all patches with iterator patchI
    {
//        Info <<"  Patch for loop: " << patchNames_[patchI] << endl;
        label patchID = mesh_.boundaryMesh().findPatchID( patchNames_[patchI]  ); // get patchID from each patch
//        Info <<"  patchID: " << patchID << endl;
	const fvPatch& myPatch = mesh_.boundary()[patchID]; // define patch
//        Info << "  Patch.size() on processor[" << Pstream::myProcNo() << "] = " << myPatch.size() << endl;

        // Initialize a plain continous array for the data to write to hdf5 later
        ioScalar* vectorPatchDataArray;
        vectorPatchDataArray = new ioScalar[myPatch.size()*3]; //  initialize array with patch.size() components


    	forAll(vectorFields_, vfieldI) // Loop over all vector fields, iterator = fieldI
        {

		// Info which vector fields are to be written
//        	Info<< "    patchWriteVector: " << vectorFields_[vfieldI] << endl;
                const volVectorField& myVectorField = mesh_.lookupObject<volVectorField>( vectorFields_[vfieldI] );

                forAll(myPatch,faceI)
                {
		       // Loop through the field and construct the array
		        vectorPatchDataArray[3*faceI+0] = myVectorField.boundaryField()[patchID][faceI].x();
            		vectorPatchDataArray[3*faceI+1] = myVectorField.boundaryField()[patchID][faceI].y();
            		vectorPatchDataArray[3*faceI+2] = myVectorField.boundaryField()[patchID][faceI].z();
	               //vectorData[x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4......] this is how the array gets filled
	
//			Info <<"\t Patch[" << patchID << "]: " << patchNames_[patchI] << ", Vectorfeld: " << vectorFields_[vfieldI] <<
//                        "   Face[" << faceI << "] = " << myVectorField.boundaryField()[patchID][faceI] << endl;
	
//			Info <<"\t Array.x[" << faceI << "] = " << vectorPatchDataArray[3*faceI+0] << 
//				" , Array.y[" << faceI << "] = " << vectorPatchDataArray[3*faceI+1] <<
//				" , Array.z[" << faceI << "] = " << vectorPatchDataArray[3*faceI+2] << endl;

		} // end of loop that fills array

	        // Create the different datasets (needs to be done collectively)
        	char datasetName[80];
   		hsize_t dimsf[2]; // here 2 ... 2D array for vector
        	hid_t fileSpace;
        	hid_t dsetID;
        	hid_t plistID;
        	hid_t plistDCreate;

        	// loop through nCells_ with iterator proc
        	forAll(nPatchCellsPerProc_, proc)
		{

	           	// Create the dataspace for the dataset
	            	dimsf[0] = nPatchCellsPerProc_[proc][patchI]; // cells per patch per processor
	            	dimsf[1] = 3; // 3 velocity components
	            	fileSpace = H5Screate_simple(2, dimsf, NULL);

	            	// Set property to create parent groups as neccesary
	            	plistID = H5Pcreate(H5P_LINK_CREATE);
	            	H5Pset_create_intermediate_group(plistID, 1);

	            	// Set chunking, compression and other HDF5 dataset properties
	            	plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
	            	dsetSetProps(3, sizeof(ioScalar), nPatchCellsPerProc_[proc][patchI], plistDCreate);

	            	// Create the dataset for points
	            	sprintf
	                	(
	                    	datasetName,
	                    	"PATCHFIELDS/%s/%s/processor%i/%s",
	                    	mesh_.time().timeName().c_str(),
				patchNames_[patchI].c_str(),
	                    	proc,
	                    	vectorFields_[vfieldI].c_str()
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

        	} // end loop over nPatchCellsPerProc_ with iterator proc

		if ( nPatchCellsPerProc_[Pstream::myProcNo()][patchI] != 0 )
                {

        		// Open correct dataset for this process
        		sprintf
            			(
                		datasetName,
                       		"PATCHFIELDS/%s/%s/processor%i/%s",
                    		mesh_.time().timeName().c_str(),
				patchNames_[patchI].c_str(),
                		Pstream::myProcNo(),
                		vectorFields_[vfieldI].c_str()
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
                        	vectorPatchDataArray
                        	);

	        	// Close/release resources.
        		H5Pclose(plistID);
        		H5Dclose(dsetID);

		} // end of if loop

	} // end of loop over vectorFields

	// Release memory, vectorData array not needed anymore after hdf5 was written
        delete [] vectorPatchDataArray;


     }// end of loop over patches

}


// ************************************************************************* //
