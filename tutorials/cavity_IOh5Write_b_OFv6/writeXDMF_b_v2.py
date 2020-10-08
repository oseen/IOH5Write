#!/usr/bin/env python3


"""writeXDMF.py:
Script to parse an HDF5 file written by OpenFOAM and write
corresponding XDMF file(s), allowing the data to be postprocessed.
"""

__author__ = "Håkon Strandenes"
__email__ = "hakon.strandenes@ntnu.no"
__copyright__ = "Copyright 2012-2013"
__license__ = "GPL v3 or later"

# Extended for patch fields and further minor changes
__author2__ = "Flo Wachter"
__email2__ = "flo.wachter@fau.de"


import numpy
import h5py
import math
import re
import argparse
import os


# Definitions and Functions:
############################################################################

# Define a dictionary with the XDMF attribute types
xdmfAttrTypes = {1: 'Scalar', 3: 'Vector', 6: 'Tensor6', 9: 'Tensor'}

# Functions that checks if there are any clouds present and returns a list with their names
def detectClouds(f):
    # Check if the file contain any clouds
    if ('CLOUDS' in f['/']) is False:
        print('Note: No clouds present in file \'{}\''.format(args.input))
        return False

    clouds = list(f['CLOUDS'].keys())
    return clouds


# Function that checks if there are any fields present
def detectFields(f):
    # Check if the file contain any fields
    if ('MESH' in f) is False or ('INTERNALFIELDS' in f) is False:
        print('Note: No mesh and field data in file \'{}\''.format(args.input))
        return False

    return True


# Function that checks if there are any pacthes present
def detectPatches(f):
    # Check if the file contain any fields
    if ('MESHPATCHES' in f) is False or ('PATCHFIELDS' in f) is False:
        print('Note: No patch mesh and patch field data in file \'{}\''.format(args.input))
        return False

    return True


#-------------------------------------------------------------------------#

# Function to parse and write patch fields
def writePatchFields(f, fo, args):

    # Print header
    fo.write('<Xdmf>\n')
    fo.write('  <Domain>\n')
    fo.write('    <Grid Name="PatchFieldData" GridType="Collection" CollectionType="Temporal">\n\n')

    # Find name and path to HDF5-file relative to XDMF-file
    h5Path = os.path.relpath(f.filename, os.path.dirname(fo.name))

    # some useful handles - handles indicated with _handleName_
    _firstMeshPatchTime_ = list(f['MESHPATCHES'].keys())[0] # handle to first mesh time 
    _allPatchFieldTime_ = f['PATCHFIELDS'] # handle to simulation times
    _firstPatchFieldTime_ = list(f['PATCHFIELDS'].keys())[0] # handle to first simulation time
    _allPatches_ = f['PATCHFIELDS'][_firstPatchFieldTime_] # handle to patches
    _allMeshPatches_ = f['MESHPATCHES'][_firstMeshPatchTime_] # handle to patches
    _firstPatch_ = list(_allPatches_.keys())[0] # handle to first patch
    _allProcsFirstTime_ = f['PATCHFIELDS'][_firstPatchFieldTime_][_firstPatch_] # handle to first proc of first time


    nPatches = len(_allPatches_) # number of patches
    patchNameArray = [None]*nPatches #init. array that contains patches
    timeLoopCount = -1
    
    # LOOP OVER SIMULATION TIME NAMES
    for timeIter in _allPatchFieldTime_:
        
        timeLoopCount = timeLoopCount+1

        fo.write('\n')
        fo.write('      <Grid GridType="Collection" CollectionType="Spatial">\n')
        fo.write('       <Time Type="Single" Value="{}" />\n'.format(timeIter))

        # LOOP OVER NUMBER OF PATCHES
        for i in range(nPatches):


            # handle to the processors of each patch
            _procOfPatch_ = f['PATCHFIELDS'][list(_allPatchFieldTime_.keys())[timeLoopCount]][list(_allPatches_.keys())[i]]
            _procOfMeshPatch_ = f['MESHPATCHES'][_firstMeshPatchTime_][list(_allMeshPatches_.keys())[i]]

            nProcs = len(_procOfPatch_) # number of involved processors
            procArray = [None]*nProcs # init. array that contains number of processors
            
            patchName = list(_allPatches_.keys())[i]
            patchNameArray[i] = patchName # fill array
            nFacesArray = [None]*nProcs # init.array for number of faces of each patch
            faceDataLengthArray = [None]*nProcs # init.array for length of the data that describe the faces

            # LOOP OVER PROCESSOR NAMES  (procIter=processor0,processor1,...)
            count_procIter = 0
            for procIter in _procOfPatch_:
            
                attrs = h5py.AttributeManager(_procOfMeshPatch_[procIter]['FACES']) # handle to attribute
                nFacesArray[count_procIter]=attrs['nFaces'][0] # fill array
                faceDataLengthArray[count_procIter]=len(_procOfMeshPatch_[procIter]['FACES']) #fill array
                _meshPoints_ = f['MESH'][_firstMeshPatchTime_][procIter]['POINTS'] # handle to mesh point of each processor
                nPoints = len(_meshPoints_) # number of points per processor
                pointPrec = _meshPoints_.dtype.itemsize # point precision

                if nFacesArray[count_procIter] != 0: # if processor has no faces on the patch do not write anything

                    fo.write('\n')
                    fo.write('       <Grid Name="time{}-{}-{}" Type="Uniform">\n'.format(timeIter,patchName,procIter))
                    fo.write('         <Topology Type="Mixed" NumberOfElements="{}">\n'.format(nFacesArray[count_procIter]))
                    fo.write('           <DataStructure Dimensions="{}" NumberType="Int" Format="HDF">\n'.format(faceDataLengthArray[count_procIter]))

                    fo.write('             {}:/MESHPATCHES/{}/{}/{}/FACES\n'.format(h5Path,_firstMeshPatchTime_,patchNameArray[i],procIter))                
                    fo.write('           </DataStructure>\n')                
                    fo.write('         </Topology>\n')                
                    fo.write('         <Geometry GeometryType="XYZ">\n')                
                    fo.write('           <DataStructure Dimensions="{} 3" NumberType="Float" Precision="{}" Format="HDF" >\n'.format(nPoints,pointPrec))
                    fo.write('             {}:/MESH/{}/{}/POINTS\n'.format(h5Path,_firstMeshPatchTime_,procIter))
                    fo.write('           </DataStructure>\n')                
                    fo.write('         </Geometry>\n')  

                    _patchFields_ = f['PATCHFIELDS'][timeIter][patchName][procIter] # handle to fields of the patch (p,U,T,..)

                    #LOOP OVER FIELDS (p,U,...)
                    for fieldIter in _patchFields_:

                        patchFieldPrec = f['PATCHFIELDS'][timeIter][patchName][procIter][fieldIter].dtype.itemsize
                        si = _patchFields_[fieldIter].size
                        sh = _patchFields_[fieldIter].shape[0]
                        nComponents = int(si/sh)
                        tensorType = xdmfAttrTypes[nComponents]


                        fo.write('         <Attribute Name="{}" Center="Cell" AttributeType="{}">\n'.format(fieldIter,tensorType))
                        fo.write('           <DataStructure Format="HDF" DataType="Float" Precision="{}" Dimensions="{} {}">\n'.format(patchFieldPrec,nFacesArray[count_procIter],nComponents))                
                        fo.write('             {}:/PATCHFIELDS/{}/{}/{}/{}\n'.format(h5Path,timeIter,patchName,procIter,fieldIter))
                        fo.write('           </DataStructure>\n')                
                        fo.write('         </Attribute>\n')                
                    
                    fo.write('       </Grid>\n')  # </Grid> close for processors            

                count_procIter=+1
                
        fo.write('      </Grid>\n') # </Grid> close for time               

    # Print footer
    fo.write('    </Grid>\n')
    fo.write('  </Domain>\n')
    fo.write('</Xdmf>\n')

    # Return control
    return             


#-------------------------------------------------------------------------#

# Function to parse and write internal mesh/fields (fo = file open)
def writeInternalFields(f, fo, args):

    # Some useful handles
    meshTime = list(f['MESH'].keys())[0]
    mesh = f['MESH'][meshTime]
    fields = f['INTERNALFIELDS']

    # Find number of points, cells and dataset length for each process
    nProcs = len(mesh)
    nPoints = [None]*nProcs
    nCells = [None]*nProcs
    cellLength = [None]*nProcs
    procNo = [None]*nProcs

    i = 0
    for proc in mesh:
        attrs = h5py.AttributeManager(mesh[proc]['CELLS'])
        
        procNo[i] = int(re.findall(r'\d+', proc)[0])
        nCells[i] = attrs['nCells'][0]
        cellLength[i] = len(mesh[proc]['CELLS'])
        nPoints[i] = len(mesh[proc]['POINTS'])

        i += 1

    # Create a list of scalar time values
    timeNames = list(fields.keys())
    timeValues = [float(i) for i in timeNames] # transform time strings to floats
    timeIndex = numpy.argsort(timeValues) # create time step indices 1,2,3...lastTimeStep

    # Determine the precision of the HDF5 file
    someFieldName = list(fields[timeNames[0]]['processor0'].keys())[0]
    prec = fields[timeNames[0]]['processor0'][someFieldName].dtype.itemsize

    # Find name and path to HDF5-file relative to XDMF-file
    h5Path = os.path.relpath(f.filename, os.path.dirname(fo.name))

    # Print header
    fo.write('<Xdmf>\n')
    fo.write('  <Domain>\n')
    fo.write('    <Grid Name="InternalFieldData" GridType="Collection" '
             'CollectionType="Temporal">\n\n')

    # Loop over all times
    for index in timeIndex:

        # Skip time = 0 if required
        if (args.noZero and
                index == timeIndex[0] and
                math.fabs(timeValues[index] - 0) < 1e-8):
            continue

        # Skip if it only is to process latest timestep
        if args.latestTime and index != timeIndex[timeIndex.size-1]:
            continue

        timeName = timeNames[index]

        fo.write('      <Grid GridType="Collection" '
                 'CollectionType="Spatial">\n'
                 .format(timeValues[index]))

        fo.write('        <Time Type="Single" Value="{}" />\n'
                 .format(timeValues[index]))

        # Loop over all processes
        i = 0
        for proc in mesh:
            fo.write('        <Grid Name="time{}-{}" Type="Uniform">\n'
                     .format(timeValues[index], proc))

            # Geometry definition
            fo.write('          <Topology Type="Mixed" '
                     'NumberOfElements="{}">\n'
                     .format(nCells[i]))

            fo.write('            <DataStructure Dimensions="{}" '
                     'NumberType="Int" Format="HDF" >\n'.format(cellLength[i]))

            fo.write('              {}:/MESH/{}/{}/CELLS\n'
                     .format(h5Path, meshTime, proc))
            fo.write('            </DataStructure>\n')
            fo.write('          </Topology>\n')
            fo.write('          <Geometry GeometryType="XYZ">\n')
            fo.write('            <DataStructure Dimensions="{} 3" '
                     'NumberType="Float" Presicion="{}" Format="HDF" >\n'
                     .format(nPoints[i], prec))

            fo.write('              {}:/MESH/{}/{}/POINTS\n'
                     .format(h5Path, meshTime, proc))

            fo.write('            </DataStructure>\n')
            fo.write('          </Geometry>\n')

            # Loop over all fields
            for field in fields[timeName][proc]:
                h = fields[timeName][proc][field]
                nCmp = int(h.size/h.shape[0])

                fo.write('          <Attribute Name="{}" Center="Cell" '
                         'AttributeType="{}">\n'
                         .format(field, xdmfAttrTypes[nCmp]))

                fo.write('            <DataStructure Format="HDF" '
                         'DataType="Float" Precision="{}" '
                         'Dimensions="{} {}">\n'
                         .format(prec, nCells[i], nCmp))

                fo.write('              {}:/INTERNALFIELDS/{}/{}/{}\n'
                         .format(h5Path, timeName, proc, field))

                fo.write('            </DataStructure>\n')
                fo.write('          </Attribute>\n')

            # Write processor footer
            fo.write('        </Grid>\n')

            # Increment processor counter
            i += 1

        # Write footer for time
        fo.write('      </Grid>\n\n\n')


    # Print footer
    fo.write('    </Grid>\n')
    fo.write('  </Domain>\n')
    fo.write('</Xdmf>\n')

    # Return control
    return


############################################################################
""" Main part of program/script start below here
"""
# Set up and read command line arguments
parser = argparse.ArgumentParser(description='Script to parse an HDF5 file '
                                 'written by OpenFOAM  and write '
                                 'corresponding XDMF files')
parser.add_argument('-i', '--input',
                    default='h5Data/h5Data0.h5',
                    help='Input file name (default: \'h5Data/h5Data0.h5\')',
                    required=False)
parser.add_argument('-d', '--dir',
                    default='xdmf_b',
                    help='Output directory (default: \'xdmf_b\')',
                    required=False)
parser.add_argument('-l', '--noLagrangian',
                    action='store_true',
                    help='>>> Skip Lagrangian clouds',
                    required=False)
parser.add_argument('-f', '--noFields',
                    action='store_true',
                    help='>>> Skip internal fields data',
                    required=False)
parser.add_argument('-p', '--noPatches',
                    action='store_true',
                    help='>>> Skip patch fields data',
                    required=False)
parser.add_argument('-z', '--noZero',
                    action='store_true',
                    help='>>> Do not process time zero (t=0)',
                    required=False)
parser.add_argument('-t', '--latestTime',
                    action='store_true',
                    help='>>> Only process latest timestep present',
                    required=False)

args = parser.parse_args()

print('Reading from \'{}\''.format(args.input))
print('Printing XDMF files to folder \'{}\''.format(args.dir))

# Check if file exist
if os.path.isfile(args.input) is False:
    print('Error: File \'{}\' does not exist!'.format(args.input))
    exit(1)

# Open HDF5 file
f = h5py.File(args.input, 'r')

# Check if XDMF output directory exist, create if not
try:
    os.stat(args.dir)
except:
    os.mkdir(args.dir)

# Write patch fields if present
patchesPresent = detectPatches(f)
if args.noPatches is False and patchesPresent:
    fo = open('{}/patchFieldData.xdmf'.format(args.dir), 'w')
    print('Writing patch field data to file {}'.format(fo.name))
    writePatchFields(f, fo, args)
    fo.close()
else:
    print('>>> Skipping patch fields')


# Write internal fields if present
fieldsPresent = detectFields(f)
if args.noFields is False and fieldsPresent:
    fo = open('{}/internalFieldData.xdmf'.format(args.dir), 'w')
    print('Writing internal field data to file {}'.format(fo.name))
    writeInternalFields(f, fo, args)
    fo.close()
else:
    print('>>> Skipping internal fields')


# Detect, loop over clouds and write them
clouds = detectClouds(f)
if args.noLagrangian is False and clouds:
    for cloud in clouds:
        fo = open('{}/{}.xdmf'.format(args.dir, cloud), 'w')
        H = f['CLOUDS'][cloud]
        print('Writing cloud data for \'{}\' to {}'.format(cloud, fo.name))
        writeCloud(H, fo, args)
        fo.close()
else:
    print('>>> Skipping Lagrangian clouds')

# Close HDF5 file
f.close()
