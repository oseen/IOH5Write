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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::h5Write::fieldWrite()
{
    fieldWriteScalar();
    fieldWriteVector();
}

//--------------------------------------------------------------------//

void Foam::functionObjects::h5Write::fieldWriteScalar()
{
    forAll(scalarFields_, fieldI) // Loop over all scalar fields, iterator = fieldI
    {
        
        const volScalarField& field = obr_.lookupObject<volScalarField>
            (
                scalarFields_[fieldI]
            );
        
        
        // Initialize a plain continous array for the data
        ioScalar* scalarData;
        scalarData = new ioScalar[field.size()]; // initialize array with field.size() components
        
        // Loop through field (iterator = iter) and construct the array
        forAll(field, iter)
        {
            scalarData[iter] = field[iter]; // fill array scalarData with field information
        }
        
        
        // Create the different datasets (needs to be done collectively)
        char datasetName[80];
        hsize_t dimsf[1];	// define dataset dimensions
        hid_t fileSpace;	// fileSpace handle
        hid_t dsetID;		// handle
        hid_t plistID;		// handle
        hid_t plistDCreate;	// handle
        
	// loop through nCells_ with iterator = proc = processor number (0,1,2,..)
        forAll(nCells_, proc)
        {
            // Create the dataspace for the dataset with number of cells on the current processor
            dimsf[0] = nCells_[proc];
            
            fileSpace = H5Screate_simple(1, dimsf, NULL); // H5Screate_simple(RANK, dimsf, NULL);
            
            // Set property to create parent groups as neccesary
            plistID = H5Pcreate(H5P_LINK_CREATE);
            H5Pset_create_intermediate_group(plistID, 1);
            
            // Set chunking, compression and other HDF5 dataset properties
            plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
            dsetSetProps(1, sizeof(ioScalar), nCells_[proc], plistDCreate);
            
            // Create the dataset for points
	    // see h5dump_H --> GROUP "0.025" { GROUP "processor0" { DATASET "Lambda2" { } }}}

            sprintf
                (
                    datasetName,
                    "INTERNALFIELDS/%s/processor%i/%s",
                    mesh_.time().timeName().c_str(),
                    proc,
                    scalarFields_[fieldI].c_str()
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

        }// here does loop over all nCells_ end
        
        
        // Open correct dataset for this process
        sprintf
            (
                datasetName,
                "INTERNALFIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                scalarFields_[fieldI].c_str()
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
                scalarData
	    );
        
        // Close/release resources.
        H5Pclose(plistID);
        H5Dclose(dsetID);
        
        // Release memory
	// array scalarData obsolete after written to hdf5
        delete [] scalarData;

    } // here does loop over scalarFields_ end

}

//-------------------------------------------------------------------//
// Do the same for all vectors (velocity, ...)

void Foam::functionObjects::h5Write::fieldWriteVector()
{
    forAll(vectorFields_, fieldI) // Loop through all vectorFields_ with iterator fieldI
    {
        
        const volVectorField& field = obr_.lookupObject<volVectorField>
            (
                vectorFields_[fieldI]
            );
        
        //Initialize a plain continous array for the data
        ioScalar* vectorData;
        vectorData = new ioScalar[field.size()*3];
        
        
        // Loop through the field and construct the array
        forAll(field, iter)
        {
            vectorData[3*iter+0] = field[iter].x();
            vectorData[3*iter+1] = field[iter].y();
            vectorData[3*iter+2] = field[iter].z();
	    //vectorData[x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4......] this is how the array gets filled
        }
        
        
        // Create the different datasets (needs to be done collectively)
        char datasetName[80];
        hsize_t dimsf[2]; 
        hid_t fileSpace;
        hid_t dsetID;
        hid_t plistID;
        hid_t plistDCreate;
        
	// loop through nCells_ with iterator proc
        forAll(nCells_, proc)
        {
            
            // Create the dataspace for the dataset
            dimsf[0] = nCells_[proc];
            dimsf[1] = 3;
            fileSpace = H5Screate_simple(2, dimsf, NULL);
            
            // Set property to create parent groups as neccesary
            plistID = H5Pcreate(H5P_LINK_CREATE);
            H5Pset_create_intermediate_group(plistID, 1);
            
            // Set chunking, compression and other HDF5 dataset properties
            plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
            dsetSetProps(3, sizeof(ioScalar), nCells_[proc], plistDCreate);
            
            // Create the dataset for points
            sprintf
                (
                    datasetName,
                    "INTERNALFIELDS/%s/processor%i/%s",
                    mesh_.time().timeName().c_str(),
                    proc,
                    vectorFields_[fieldI].c_str()
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

        } // end loop over nCells_ with iterator proc
        
        
        // Open correct dataset for this process
        sprintf
            (
                datasetName,
                "INTERNALFIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                vectorFields_[fieldI].c_str()
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
                vectorData
		        );
        
        // Close/release resources.
        H5Pclose(plistID);
        H5Dclose(dsetID);
        
        // Release memory, vectorData array not needed anymore after hdf5 was written
        delete [] vectorData;

    } // end loop over vectorFields with iterator fieldI

}


// ************************************************************************* //
