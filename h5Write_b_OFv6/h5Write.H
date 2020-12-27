/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

Class
    Foam::functionObjects::h5Write

Description
    Writes selected (vector and scalar) fields from the internalField or from patches or boundariesto a XDMF+HDF5 file archive.
    Require HDF5 version >= 1.8.7 (tested with HDF5-1.10.5)

Author
    Håkon Strandenes, hakon.strandenes@ntnu.no, 2012-2014
    Florian Wachter, flo.wachter@fau.de, 2019-2021

    Example of function object specification:
    \verbatim
    h5data1
    {
        type           h5Write;
        libs ("libIOh5Write.so");
	...
	suppressFieldDataWrite false;
	objectNames     ( U p k nut );
	patchNames 	( inlet outlet wall);
	cloudNames      (  );
	cloudAttribs    (  );
	chunkSize       0;
	writeInterval   20;
        ...
    }
    \endverbatim


SourceFiles
    h5Write.C

\*---------------------------------------------------------------------------*/

#ifndef h5Write_H
#define h5Write_H

#include "fvMeshFunctionObject.H"
#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "logFiles.H"

#include "mpi.h"
#include "hdf5.h"


// If OpenFOAM is compiled with single precision, the IO will be single
// precision even if WRITE_DP is set.

#if defined(WRITE_DP) && defined(WM_DP)
#   define H5T_SCALAR H5T_NATIVE_DOUBLE
#   define ioScalar   double
#elif defined(WRITE_SP)
#   define H5T_SCALAR H5T_NATIVE_FLOAT
#   define ioScalar   float
#endif

#define H5_XFER_MODE H5FD_MPIO_INDEPENDENT


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

 // Forward declaration of classes
 class objectRegistry;
 class dictionary;
 class fvMesh;


namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class h5Write Declaration
\*---------------------------------------------------------------------------*/

class h5Write
:
    public fvMeshFunctionObject
{

        //- Class template used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}

        };


    // Private data

        //- word
        word name_;


	// Read from dictionary

            //- Names of volume fields to write
            wordReList objectNames_;

            //- Names of clouds to write
            wordReList cloudNames_;

            //- Which attributes to write for particle data
            wordReList cloudAttribs_;

            //- HDF5 chunk size (in bytes)
            label chunkSize_;

            //- Write interval (in number of iterations)
            label writeInterval_;

            //- Names of patches to write
            wordReList patchNames_;

            //--> true/false
            bool suppressFieldDataWrite_ = false;


	// Categorized scalar/vector/tensor fields
        fieldGroup<scalar> scalarFields_;
        fieldGroup<vector> vectorFields_;


        // Used to track when to write
        label nextWrite_;
        label timeSteps_;


        // Number of fields (vector, scalar) to write
        label nFields_;


        // File and dataset identifier 
        hid_t fileID_;


        // Array containing number of cells per process (often used list)
        List<label> nCells_;
        // Array containing number of patch cells per process (often used list)
        List<label> nPatchCells_; // nPatchCells_[patch] = number of patch-cells
        List<List<label>> nPatchCellsPerProc_; // nPatchCells_[proc][patch] = number of patch-cells per processor

	List<label> nPatchFacePoints_; // sum of face points per patch
	List<List<label>> nPatchFacePointsPerProc_; // sum of face points per patch per processor


    // Private Member Functions
    
       // General functions

        // Create and open dataset for field data (vector and scalar data)
        void fileCreate();

        // Close dataset for field data
        void fileClose();

	// Set dataset creation properties such as chunking
        void dsetSetProps
          (
              label nCmps,
              label size,
              label maxRows,
              hid_t plistDCreate
          );



	  // Functions for handling field data

          // Field write wrapper
          void fieldWrite();
          void patchWrite();

          // Write scalar data
          void fieldWriteScalar();
          void patchWriteScalar();

          // Write vector data
          void fieldWriteVector();
          void patchWriteVector();




	  // Functions for handling particle cloud data

          // Cloud write wrapper
          void cloudWrite();

          // Write attributes for cloud particles
          void cloudWriteAttrib
          (
              label myParticles,
              label offset,
              label nTot,
              label nCmps,
              void* databuf,
              char* datasetName,
              hid_t H5type
          );

      // Functions for handling mesh

          // Mesh write wrapper
          void meshWrite();

          // Write points
          void meshWritePoints();

          // Write cells
          void meshWriteCells();

          // Write faces of patches
          void meshWritePatchFaces();



      // Miscallenous definitions

        //- Append fieldName to the appropriate group
        label appendFieldGroup(const word& fieldName, const word& fieldType);

        //- Classify field types, returns the number of fields
        label classifyFields();


        //- Disallow default bitwise copy construct
        h5Write(const h5Write&);

        //- Disallow default bitwise assignment
        void operator=(const h5Write&);


      //- Pointers to cell shape models
      static const cellModel* unknownModel;
      static const cellModel* tetModel;
      static const cellModel* pyrModel;
      static const cellModel* prismModel;
      static const cellModel* hexModel;



public:

    //- Runtime type information
    TypeName("h5Write");


    // Constructors

        //- Construct from name, time and dictionary
        h5Write
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~h5Write();


    // Member Functions


        //- Return name of the h5Write
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the h5Write data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the h5Write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
