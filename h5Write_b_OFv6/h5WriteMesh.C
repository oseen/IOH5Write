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
#include "cellModeller.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::h5Write::meshWrite()
{

    // Find over all (global) number of cells per process
    nCells_[Pstream::myProcNo()] = mesh_.cells().size();
    Pstream::gatherList(nCells_);
    Pstream::scatterList(nCells_);

    ////////////////////////////////////////////////////////////////////////////////

    // Find over all (global) number of patch faces per process

        nPatchCells_.setSize( patchNames_.size() ); // set size of List<label> nPatchCells_
	nPatchFacePoints_.setSize( patchNames_.size() );
	nPatchFacePointsPerProc_.setSize(Pstream::nProcs());

        forAll(patchNames_, patchI) //loop over all patches with iterator patchI
        {
                label patchID = mesh_.boundaryMesh().findPatchID( patchNames_[patchI]  ); // get patchID from each patch
	    if ( patchID >=0 )
	    {
		nPatchCells_[patchI] = mesh_.boundaryMesh()[patchID].size(); // fill list nPatchCells with nr. of patch cells for each patch

		const fvPatch& myPatch = mesh_.boundary()[patchID]; // define patch

		int myPatchFacePointsNumber = 0;

		forAll(myPatch, faceI)
		{
		    const label& faceID = myPatch.start() + faceI;

		    myPatchFacePointsNumber++; // shapeID of polygon, 3

		    myPatchFacePointsNumber++; // polygon face points count

		    myPatchFacePointsNumber += mesh_.faces()[faceID].size(); // for face points
		}
		nPatchFacePoints_[patchI] = myPatchFacePointsNumber;
	    }
	    else
	    {
		FatalErrorIn
		    (
		     "h5Write::meshWrite()"
		    )   << "Unknown patch "
		    << patchNames_[patchI] << endl
		    << exit(FatalError);
	    }
	}

        // fill List<List<label>> nPatchCellsPerProc_ with List<label> nPatchCells_
        nPatchCellsPerProc_[Pstream::myProcNo()] = nPatchCells_; // fill nPatchCellsPerProc_ with nr. of patch cells for each patch for each proc

        // share with all mpi processes
        Pstream::gatherList(nPatchCellsPerProc_);
        Pstream::scatterList(nPatchCellsPerProc_);

	// fill List<List<label>> nPatchFacePointsPerProc_
	// with List<label> nPatchFacePoints_;
	nPatchFacePointsPerProc_[Pstream::myProcNo()] = nPatchFacePoints_;

	// share with all mpi processes
	Pstream::gatherList(nPatchFacePointsPerProc_);
	Pstream::scatterList(nPatchFacePointsPerProc_);

    ////////////////////////////////////////////////////////////////////////////////


    // Write mesh
    meshWritePoints();
    meshWriteCells();

    // Write mesh patch faces
    meshWritePatchFaces();

}


void Foam::functionObjects::h5Write::meshWritePoints()
{


    const pointField& points = mesh_.points();


    // Find out how many points each process has
    List<label> nPoints(Pstream::nProcs());
    nPoints[Pstream::myProcNo()] = points.size();
    Pstream::gatherList(nPoints);
    Pstream::scatterList(nPoints);


    // Create the different datasets (needs to be done collectively)
    char datasetName[80];
    hsize_t dimsf[2];
    hid_t fileSpace;
    hid_t dsetID;
    hid_t plistID;
    hid_t plistDCreate;

    forAll(nPoints, proc)
    {

        // Create the dataspace for the dataset
        dimsf[0] = nPoints[proc];
        dimsf[1] = 3;
        fileSpace = H5Screate_simple(2, dimsf, NULL);

        // Set property to create parent groups as neccesary
        plistID = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(plistID, 1);

        // Set chunking, compression and other HDF5 dataset properties
        plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
        dsetSetProps(3, sizeof(ioScalar), nPoints[proc], plistDCreate);

        // Create the dataset for points
        sprintf
            (
                datasetName,
                "MESH/%s/processor%i/POINTS",
                mesh_.time().timeName().c_str(),
                proc
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
    }


    // Create a simple array of points (to pass on to H5Dwrite)
    ioScalar pointList[points.size()][3];
    forAll(points, ptI)
    {
        pointList[ptI][0] = points[ptI].x();
        pointList[ptI][1] = points[ptI].y();
        pointList[ptI][2] = points[ptI].z();
    }


    // Open correct dataset for this process
    // see with h5dump -H --> GROUP "MESH" { GROUP "0" { GROUP "processor0" { DATASET "POINTS" { ... } }}}

    sprintf
        (
            datasetName,
            "MESH/%s/processor%i/POINTS",
            mesh_.time().timeName().c_str(),
            Pstream::myProcNo()
            //  variables %s --> mesh_.time().timeName().c_str(); %i --> Pstream::myProcNo()
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
                        pointList
                    );

    // Close/release resources.
    H5Pclose(plistID);
    H5Dclose(dsetID);
}




void Foam::functionObjects::h5Write::meshWriteCells()
{


    // Map shapes OpenFOAM->XDMF
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(unknownModel->index(), 16);
    shapeLookupIndex.insert(hexModel->index(), 9);
    shapeLookupIndex.insert(prismModel->index(), 8);
    shapeLookupIndex.insert(pyrModel->index(), 7);
    shapeLookupIndex.insert(tetModel->index(), 6);


    const cellList& cells  = mesh_.cells();

    const faceList& faces = mesh_.faces();

    const cellShapeList& shapes = mesh_.cellShapes();

    int myDatasetSize = 0;

    forAll(cells, cellId)
    {
	myDatasetSize += 1; // Store the type/shape of the cell

	myDatasetSize += 1; // Store the number of faces for polyhedra cells,
		            // or just skipped for non-polyhedra cells

	const cell& c = cells[cellId];

	// Get all faces of the cell
        forAll(c, faceId)
        {
	    myDatasetSize += 1; // Store face points count

	    const face& f = faces[c[faceId]];

            myDatasetSize += f.size(); // Store all points
        }
    }

    // Find dataset length for this process and fill dataset in one operation
    // this will possible give a little overhead w.r.t. storage, but on a
    // hex-dominated mesh, this is OK.
    int j = 0;
    int myDataset[myDatasetSize];
    forAll(cells, cellId)
    {
	const cell& c = cells[cellId];
        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        // A registered primitive type
	if (shapeLookupIndex.found(mapIndex))
	{
	    label shapeId = shapeLookupIndex[mapIndex];

	    if (shapeId != 16) // shapeId !=16 for non-polyhedral cells
	    {
		myDataset[j] = shapeId; j++;

		const labelList& vrtList = shapes[cellId];

		forAll(vrtList, i)
		{

		    myDataset[j] = vrtList[i]; j++;

		}
	    }

	    // For polyhedral cells
	    else
	    {
		// Polyhedral cell type is 16
		// see https://www.xdmf.org/index.php/XDMF_Model_and_Format
		myDataset[j] = 16; j++;

		myDataset[j] = c.size(); j++; // write number of faces

		forAll(c, faceId)
		{

		    const face& f = faces[c[faceId]];

		    myDataset[j] = f.size(); j++; // write number of points of face

		    forAll(f, nodei)
		    {
			const label& nodeId = f[nodei];

			myDataset[j] = nodeId; j++; // write each point
		    }

		}
	    }

	}

        // If the cell is not a basic type, exit with an error
        else
        {
            FatalErrorIn
            (
                "h5Write::meshWriteCells()"
            )   << "Unsupported or unknown cell type for cell number "
                << cellId << endl
                << exit(FatalError);
        }
    }


    // Find out how many points each process has
    List<label> datasetSizes(Pstream::nProcs());
    datasetSizes[Pstream::myProcNo()] = j;
    Pstream::gatherList(datasetSizes);
    Pstream::scatterList(datasetSizes);

        // Create the different datasets (needs to be done collectively)
    char datasetName[80];
    hsize_t dimsf[1];
    hid_t fileSpace;
    hid_t dsetID;
    hid_t attrID;
    hid_t plistID;
    hid_t plistDCreate;

    forAll(datasetSizes, proc)
    {
        // Set property to create parent groups as neccesary
        plistID = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(plistID, 1);

        // Set chunking, compression and other HDF5 dataset properties
        plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
        dsetSetProps(1, sizeof(int), datasetSizes[proc], plistDCreate);

        // Create dataspace for cell list
        dimsf[0] = datasetSizes[proc];
        fileSpace = H5Screate_simple(1, dimsf, NULL);

        //  h5dump -H --> GROUP "MESH" { GROUP "0" { GROUP "processor0" { DATASET "CELLS" { ... } }}}
        sprintf
            (
                datasetName,
                "MESH/%s/processor%i/CELLS",
                mesh_.time().timeName().c_str(),
                proc

                //  variables %s --> mesh_.time().timeName().c_str(); %i --> proc
            );

        dsetID = H5Dcreate2
            (
                fileID_,
                datasetName,
                H5T_NATIVE_INT,
                fileSpace,
                plistID,
                plistDCreate,
                H5P_DEFAULT
            );
        H5Sclose(fileSpace);


        // Create and write attributte to store number of cells
        dimsf[0] = 1;
        fileSpace = H5Screate_simple(1, dimsf, NULL);

        attrID = H5Acreate2
            (
                dsetID,
                "nCells",
                H5T_NATIVE_INT,
                fileSpace,
                H5P_DEFAULT,
                H5P_DEFAULT
            );

        int nCells = nCells_[proc];
        H5Awrite
            (
                attrID,
                H5T_NATIVE_INT,
                      &nCells
                  );

        H5Aclose(attrID);
        H5Sclose(fileSpace);

        // Close last access pointers
        H5Dclose(dsetID);
        H5Pclose(plistID);
        H5Pclose(plistDCreate);
    }


    // Create property list for collective dataset write.
    plistID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistID, H5_XFER_MODE);


    // Open and write correct dataset for this process
    sprintf
        (
            datasetName,
            "MESH/%s/processor%i/CELLS",
            mesh_.time().timeName().c_str(),
            Pstream::myProcNo()
        );
    dsetID = H5Dopen2(fileID_, datasetName, H5P_DEFAULT);

    H5Dwrite
        (
            dsetID,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            plistID,
            myDataset
        );

    H5Pclose(plistID);
    H5Dclose(dsetID);
}// end of Foam::h5Write::meshWriteCells()

////////////////////////////////////////////////////////////////////


void Foam::functionObjects::h5Write::meshWritePatchFaces()
{


    // Map face shapes OpenFOAM -> XDMF (polygon=3,triangle=4, quadrilateral=5)
/*    Map<label> faceLookupIndex;
    faceLookupIndex.insert(hexModel->index(), 0);
    faceLookupIndex.insert(prismModel->index(), 1);
    faceLookupIndex.insert(tetModel->index(), 2);
    faceLookupIndex.insert(pyrModel->index(), 3);
*/

    forAll(patchNames_, patchI) //loop over all patches with iterator patchI
    {
        label patchID = mesh_.boundaryMesh().findPatchID( patchNames_[patchI]  ); // get patchID from each patch
        const fvPatch& myPatch = mesh_.boundary()[patchID]; // define patch
//      Info <<"  Patch for loop: " << patchNames_[patchI] <<" with patchID "<< patchID << endl;

        int j = 0;
        // initialize plain array with number of faces for each patch on each proc * 5 ( shapeID + 4 points --> quadrilateral faces)
        int myFaceDataset[nPatchFacePointsPerProc_[Pstream::myProcNo()][patchI]];
        label shapeId = 3; //treat every face as polygon

        forAll(myPatch,faceI) //loop over all faces of patch
        {
                const label& faceID = myPatch.start() + faceI;

                myFaceDataset[j] = shapeId; // first array component is shapeID, seperate nodes relating to a face with shapeID
                j++;

                myFaceDataset[j] = mesh_.faces()[faceID].size(); // second array component is number of points of the face
                j++;

                forAll(mesh_.faces()[faceID], nodei) //loop over all nodes of one face with faceID
                        {
                                const label& nodeID = mesh_.faces()[faceID][nodei];

                                // fill array with nodes that form a face , seperate each face with shape ID
                                myFaceDataset[j]=nodeID;
//                              Info << "myFaceDataset["<<j<<"]="<<myFaceDataset[j]<< endl;
                                j++;
                        }
        } //end of loop over faces for myPatch


        // Create the different datasets (needs to be done collectively)
        char datasetName[80];
        hsize_t dimsf[1];
        hid_t fileSpace;
        hid_t dsetID;
        hid_t attrID;
        hid_t plistID;
        hid_t plistDCreate;

        forAll(nPatchCellsPerProc_, proc)
        {

	   // write only whan proc has faces on that patch ... can also be commented since we write mesh only once
           if ( nPatchCellsPerProc_[proc][patchI] != 0 )
           {

                // Set property to create parent groups as neccesary
                plistID = H5Pcreate(H5P_LINK_CREATE);
                H5Pset_create_intermediate_group(plistID, 1);

                // Set chunking, compression and other HDF5 dataset properties
                plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
                dsetSetProps(1, sizeof(int), nPatchFacePointsPerProc_[proc][patchI], plistDCreate);

                // Create dataspace for cell list
                dimsf[0] = nPatchFacePointsPerProc_[proc][patchI];
                fileSpace = H5Screate_simple(1, dimsf, NULL);

                // h5dump -H --> GROUP "MESH" { GROUP "0" { GROUP "inlet" { GROUP "processor0" { DATASET "FACES" { ... } }}}}
                sprintf
                (
                        datasetName,
                        "MESHPATCHES/%s/%s/processor%i/FACES",
                        mesh_.time().timeName().c_str(),
                        patchNames_[patchI].c_str(),
                        proc
                );

                dsetID = H5Dcreate2
                (
                        fileID_,
                        datasetName,
                        H5T_NATIVE_INT,
                        fileSpace,
                        plistID,
                        plistDCreate,
                        H5P_DEFAULT
                );

                H5Sclose(fileSpace);


                // Create and write attributte to store number of faces of each processor
                dimsf[0] = 1;
                fileSpace = H5Screate_simple(1, dimsf, NULL);

                attrID = H5Acreate2
                (
                        dsetID,
                        "nFaces",
                        H5T_NATIVE_INT,
                        fileSpace,
                        H5P_DEFAULT,
                        H5P_DEFAULT
                );

                int nFaces = nPatchCellsPerProc_[proc][patchI];
                H5Awrite
                (
                        attrID,
                        H5T_NATIVE_INT,
                        &nFaces
                );

                H5Aclose(attrID);
                H5Sclose(fileSpace);

                // Close last access pointers
                H5Dclose(dsetID);
                H5Pclose(plistID);
                H5Pclose(plistDCreate);

           }// end "if ( nPatchCellsPerProc_[proc][patchI] != 0 )" 

        }// end loop over nPatchCellsPerProc_



        // Open and write correct dataset for this process and only write when processor has faces on that patch
        if ( nPatchCellsPerProc_[Pstream::myProcNo()][patchI] != 0 )
        {

                // Create property list for collective dataset write.
                plistID = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plistID, H5_XFER_MODE);

                sprintf
                (
                        datasetName,
                        "MESHPATCHES/%s/%s/processor%i/FACES",
                        mesh_.time().timeName().c_str(),
                        patchNames_[patchI].c_str(),
                        Pstream::myProcNo()
                );

                dsetID = H5Dopen2(fileID_, datasetName, H5P_DEFAULT);

                H5Dwrite
                (
                        dsetID,
                        H5T_NATIVE_INT,
                        H5S_ALL,
                        H5S_ALL,
                        plistID,
                        myFaceDataset
                );

                H5Pclose(plistID);
                H5Dclose(dsetID);

        } // end of if loop

    } // end of loop over patches

    // Release memory
    //delete [] myFaceDataset;

} // end of Foam::h5Write::meshWritePatchFaces()

////////////////////////////////////////////////////////////////////


const Foam::cellModel* Foam::functionObjects::h5Write::unknownModel = Foam::cellModel::ptr
(
    "unknown"
);


const Foam::cellModel* Foam::functionObjects::h5Write::tetModel = Foam::cellModel::
ptr
(
    "tet"
);


const Foam::cellModel* Foam::functionObjects::h5Write::pyrModel = Foam::cellModel::
ptr
(
    "pyr"
);


const Foam::cellModel* Foam::functionObjects::h5Write::prismModel = Foam::cellModel::ptr
(
    "prism"
);


const Foam::cellModel* Foam::functionObjects::h5Write::hexModel = Foam::cellModel::ptr
(
    "hex"
);


// ************************************************************************* //

