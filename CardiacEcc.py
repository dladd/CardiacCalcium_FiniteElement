#!/usr/bin/env python
'''
!> \file
!> \author David Ladd
!> \brief Main program file to simulate cardiac ECC using opencmiss
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Vijay Rajagopal
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>
'''

# Add Python bindings directory to PATH
import numpy as np
import pandas as pd
from scipy.io import loadmat
import gc
import time
import requests
import shutil
import os
from opencmiss.iron import iron

# C l a s s e s   a n d    F u n c t i o n s 
# ----------------------------------------------------
class numberIncrementor:
    """Simple incrementor for unique integer user numbers"""
    def __init__(self):
        self.__value = 0
    def getInc(self):
        self.__value += 1
        return self.__value

def download_file(local_file, url):
    """Download the local_file from the URL"""
    r = requests.get(url, stream=True)
    with open(local_file, 'wb') as f:
        shutil.copyfileobj(r.raw, f)


"""
======================================================
C O N T R O L    P A N E L
======================================================
"""
# Solve time data (ms)
startTime = 0.0
endTime = 30.00001
pdeTimestep = 0.1
odeTimestep = 0.0001
outputFrequency = 10
# Mesh data
meshDir = './input/'
lowRyrDensity = True
mitochondria = False
# The CellML model file to use
cellmlFile = meshDir + 'ryrNtroponinNfluo3_wryrscaling_wtimelag_wtimecourse.xml'
iCa = 2.0e-15

# Initial conditions
initCa = 0.1
initF = 22.92
initFCa = 2.08
initCaM = 23.529
initCaMCa = 0.471
initATP = 454.682
initATPCa = 0.318
initCaTnC = 10.0
# Diffusion parameters
diffCa = [0.22, 0.22, 0.22]
diffF = [0.042, 0.042, 0.042]
diffFCa = [0.042, 0.042, 0.042]
diffCaM = [0.025, 0.025, 0.025]
diffCaMCa = [0.025, 0.025, 0.025]
diffATP = [0.14, 0.14, 0.14]
diffATPCa = [0.14, 0.14, 0.14]
"""
======================================================
"""

# G e t   m e s h   r e s o u r c e s
# ----------------------------------------------------
resourceLinks = {
    "Combined_8Sarc_1319kNodes_node.h5" : "https://www.dropbox.com/s/mid02d7u8anpll4/Combined_8Sarc_1319kNodes_node.h5?dl=1",
    "Combined_8Sarc_1319kNodes_elem.h5" : "https://www.dropbox.com/s/nkx1z7n3woh75a9/Combined_8Sarc_1319kNodes_elem.h5?dl=1",
    "Combined_8Sarc_1436kNodes_elem.h5" : "https://www.dropbox.com/s/7lmn3igzyzbk5bf/Combined_8Sarc_1436kNodes_elem.h5?dl=1",
    "Combined_8Sarc_1436kNodes_node.h5" : "https://www.dropbox.com/s/j3vldx5liqqgxy1/Combined_8Sarc_1436kNodes_node.h5?dl=1",
    "Cyto_8Sarc_1319kNodes_elem.h5" : "https://www.dropbox.com/s/redc6z8qsv35lwb/Cyto_8Sarc_1319kNodes_elem.h5?dl=1",
    "Cyto_8Sarc_1319kNodes_node.h5" : "https://www.dropbox.com/s/cstdrqo0mvf7cpz/Cyto_8Sarc_1319kNodes_node.h5?dl=1",
    "Cyto_8Sarc_1436kNodes_elem" : "https://www.dropbox.com/s/oqeu5y4iag3vye0/Cyto_8Sarc_1436kNodes_elem.h5?dl=1",
    "Cyto_8Sarc_1436kNodes_node.h5" : "https://www.dropbox.com/s/noo5u5512npvcbc/Cyto_8Sarc_1436kNodes_node.h5?dl=1",
    "Cyto_8Sarc_1319kNodes_spherical_ryr_kie_wh0.05.FixedNnd_offset05_N50_1umSpacing_tausimnum_1.mat" : "https://www.dropbox.com/s/xfuntofplizmixm/Cyto_8Sarc_1319kNodes_spherical_ryr_kie_wh0.05.FixedNnd_offset05_N50_1umSpacing_tausimnum_1.mat?dl=1",
    "Cyto_8Sarc_1436kNodes_spherical_ryr_kie_wh0.05.FixedNnd_offset05_N123_tausimnum_1.mat" : "https://www.dropbox.com/s/xc9p2lc6tdr87se/Cyto_8Sarc_1436kNodes_spherical_ryr_kie_wh0.05.FixedNnd_offset05_N123_tausimnum_1.mat?dl=1"
}

if mitochondria:
    meshName = 'Cyto_'
else:
    meshName = 'Combined_'
if lowRyrDensity:
    meshName = meshName + '8Sarc_1319kNodes_'
    ryrName = 'Cyto_8Sarc_1319kNodes_spherical_ryr_kie_wh0.05.FixedNnd_offset05_N50_1umSpacing_tausimnum_1.mat'
else:
    meshName = meshName + '8Sarc_1436kNodes_'
    ryrName = 'Cyto_8Sarc_1436kNodes_spherical_ryr_kie_wh0.05.FixedNnd_offset05_N123_tausimnum_1.mat'

# Download node file if it doesn't already exist
fileName = meshName + 'node.h5'
nodeFile = meshDir + fileName
if not os.path.exists(nodeFile):
    download_file(nodeFile, resourceLinks[fileName])

# Download elem file if it doesn't already exist
fileName = meshName + 'elem.h5'
elemFile = meshDir + fileName
if not os.path.exists(elemFile):
    download_file(elemFile, resourceLinks[fileName])

# Download ryr file if it doesn't already exist
ryrFile = meshDir + ryrName
if not os.path.exists(ryrFile):
    download_file(ryrFile, resourceLinks[ryrName])


# P r o b l e m    S e t u p
# ----------------------------------------------------
# Get the computational node information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()
userNumber = numberIncrementor()
if computationalNodeNumber == 0:
    print('Setting up problem...')

# Set up 3D RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(userNumber.getInc())
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create world region
region = iron.Region()
region.CreateStart(userNumber.getInc(), iron.WorldRegion)
region.label = "Region"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create basis
basis = iron.Basis()
basis.CreateStart(userNumber.getInc())
basis.numberOfXi = 3
basis.type = iron.BasisTypes.SIMPLEX
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*3
basis.CreateFinish()


# M e s h   S e t u p
# ----------------------------------------------------
# Open and close the node store to get the number of nodes
if computationalNodeNumber == 0:
    print('Setting up mesh...')
store = pd.HDFStore(nodeFile, 'r')
df = store['Node_Coordinates']
numberOfNodes = df.shape[0]
store.close()
nodes = iron.Nodes()
nodes.CreateStart(region, numberOfNodes)
nodes.CreateFinish()

# Open the elem store
store = pd.HDFStore(elemFile, 'r')
df = store['Element_Node_Map']
numberOfElements = df.shape[0]
mesh = iron.Mesh()
mesh.CreateStart(userNumber.getInc(), region, 3)
mesh.NumberOfElementsSet(numberOfElements)
mesh.NumberOfComponentsSet(1)

# Load in element dataframe
meshElements = iron.MeshElements()
meshElements.CreateStart(mesh, 1, basis)
elementNumber = 0
start = time.time()
for elementNodes in df.itertuples(index=False, name=None):
    elementNumber += 1
    meshElements.NodesSet(elementNumber, elementNodes[:4])

end = time.time()
if computationalNodeNumber == 0:
    print('Number of Nodes: ' + str(numberOfNodes))
    print('Number of Elements: ' + str(numberOfElements))
    print('Element read time: ' + str(end - start))
    print('Finalising mesh...')
meshElements.CreateFinish()
mesh.CreateFinish()
# Destroy element dataframe and collect garbage to free up memory
store.close()
gc.collect()

# Decompose mesh accross computational nodes
if computationalNodeNumber == 0:
    print('Decomposing mesh...')
decomposition = iron.Decomposition()
decomposition.CreateStart(userNumber.getInc(), mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CalculateFacesSet(False)
decomposition.CreateFinish()


# G e o m e t r i c   F i e l d
# ----------------------------------------------------
if computationalNodeNumber == 0:
    print('Setting up geometric field...')
geometricField = iron.Field()
geometricField.CreateStart(userNumber.getInc(), region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
geometricField.CreateFinish()

# Load in node df
store = pd.HDFStore(nodeFile, 'r')
df = store['Node_Coordinates']

start = time.time()
for node in range(numberOfNodes):
    nodeDomain = decomposition.NodeDomainGet(node+1, 1)
    if nodeDomain == computationalNodeNumber:
        for component in range(3):
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                    iron.FieldParameterSetTypes.VALUES,
                                                    1, 1, node+1, component+1, df.iat[node, component])

end = time.time()
if computationalNodeNumber == 0:
    print('Node read time: ' + str(end - start))
    print('Updating nodal fields')

geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)

# destroy node df and garbage collect
store.close()
del df
gc.collect()
df = pd.DataFrame()


# E q u a t i o n s   S e t s
# ----------------------------------------------------
if computationalNodeNumber == 0:
    print('Setting up equations sets...')

equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]

equationLabels = ['Ca', 'F', 'FCa', 'CaM', 'CaMCa', 'ATP', 'ATPCa']
dependentInit = [initCa, initF, initFCa, initCaM, initCaMCa, initATP, initATPCa]
materialsDiff = [diffCa, diffF, diffFCa, diffCaM, diffCaMCa, diffATP, diffATPCa]
equationsSets = []
dependentFields = []
materialsFields = []
sourceFields = []

#    M a i n   b u f f e r i n g    e q u a t i o n s
#    -------------------------------------------------
i = 0
for label in equationLabels:
    if computationalNodeNumber == 0:
        print('...' + label)
    # Equations set
    equationsSets.append(iron.EquationsSet())
    equationsSetField = iron.Field()
    equationsSets[i].CreateStart(userNumber.getInc(), region, geometricField, equationsSetSpecification,
                             userNumber.getInc(), equationsSetField)
    equationsSets[i].CreateFinish()

    # Dependent
    dependentFields.append(iron.Field())
    equationsSets[i].DependentCreateStart(userNumber.getInc(), dependentFields[i])
    dependentFields[i].LabelSet(label)
    equationsSets[i].DependentCreateFinish()
    dependentFields[i].ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                   dependentInit[i])

    # Materials
    materialsFields.append(iron.Field())
    equationsSets[i].MaterialsCreateStart(userNumber.getInc(), materialsFields[i])
    materialsFields[i].LabelSet(label + '_Materials')
    equationsSets[i].MaterialsCreateFinish()
    for c in range(3):
        materialsFields[i].ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                       c+1, materialsDiff[i][c])
    materialsFields[i].ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4, 1.0)

    # Source
    sourceFields.append(iron.Field())
    equationsSets[i].SourceCreateStart(userNumber.getInc(), sourceFields[i])
    sourceFields[i].VariableLabelSet(iron.FieldVariableTypes.U, 'i' + label)
    equationsSets[i].SourceCreateFinish()
    if label == 'Ca':
        sourceFields[i].ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, iCa)
    else:
        sourceFields[i].ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 0.0)

    i += 1


#    T r o p o n i n   C
#    -------------------
CaTnCField = iron.Field()
CaTnCField.CreateStart(userNumber.getInc(), region)
CaTnCField.TypeSet(iron.FieldTypes.GENERAL)
CaTnCField.MeshDecompositionSet(decomposition)
CaTnCField.GeometricFieldSet(geometricField)
CaTnCField.NumberOfVariablesSet(1)
CaTnCField.VariableTypesSet([iron.FieldVariableTypes.U])
CaTnCField.DataTypeSet(iron.FieldVariableTypes.U, iron.FieldDataTypes.DP)
CaTnCField.DimensionSet(iron.FieldVariableTypes.U, iron.FieldDimensionTypes.SCALAR)
CaTnCField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
CaTnCField.LabelSet('CaTnC')
CaTnCField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
CaTnCField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
CaTnCField.CreateFinish()
CaTnCField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, initCaTnC)


#    R y R    D e n s i t y   a n d    T i m e l a g
#    -----------------------------------------------
RyRDensityField = iron.Field()
RyRDensityField.CreateStart(userNumber.getInc(), region)
RyRDensityField.TypeSet(iron.FieldTypes.GENERAL)
RyRDensityField.MeshDecompositionSet(decomposition)
RyRDensityField.GeometricFieldSet(geometricField)
RyRDensityField.NumberOfVariablesSet(1)
RyRDensityField.VariableTypesSet([iron.FieldVariableTypes.U])
RyRDensityField.DataTypeSet(iron.FieldVariableTypes.U, iron.FieldDataTypes.DP)
RyRDensityField.DimensionSet(iron.FieldVariableTypes.U, iron.FieldDimensionTypes.SCALAR)
RyRDensityField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
RyRDensityField.LabelSet('RyRDensity')
RyRDensityField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
RyRDensityField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
RyRDensityField.CreateFinish()
RyRDensityField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 0.0)

RyRTimelagField = iron.Field()
RyRTimelagField.CreateStart(userNumber.getInc(), region)
RyRTimelagField.TypeSet(iron.FieldTypes.GENERAL)
RyRTimelagField.MeshDecompositionSet(decomposition)
RyRTimelagField.GeometricFieldSet(geometricField)
RyRTimelagField.NumberOfVariablesSet(1)
RyRTimelagField.VariableTypesSet([iron.FieldVariableTypes.U])
RyRTimelagField.DataTypeSet(iron.FieldVariableTypes.U, iron.FieldDataTypes.DP)
RyRTimelagField.DimensionSet(iron.FieldVariableTypes.U, iron.FieldDimensionTypes.SCALAR)
RyRTimelagField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
RyRTimelagField.LabelSet('RyRTimelag')
RyRTimelagField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
RyRTimelagField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
RyRTimelagField.CreateFinish()
RyRTimelagField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 0.0)

ryrData = loadmat(ryrFile)
for ryrNodeIdx in range(len(ryrData['nonzeroNodes'])):
    nodeNumber = np.asscalar(ryrData['nonzeroNodes'][ryrNodeIdx])
    nodeDomain = decomposition.NodeDomainGet(nodeNumber, 1)
    if nodeDomain == computationalNodeNumber:
        RyRDensityField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nodeNumber, 1,
                                                 np.asscalar(ryrData['nonzeroIntensities'][ryrNodeIdx]))
        RyRTimelagField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, nodeNumber, 1,
                                                 np.asscalar(ryrData['nonzeroTimelags'][ryrNodeIdx]))

RyRDensityField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
RyRDensityField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
RyRTimelagField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
RyRTimelagField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
del ryrData
gc.collect()



# C e l l M L
# -----------
if computationalNodeNumber == 0:
    print('Setting up CellML model and variable mapping...')
cellmlModelIndex = 1
cellml = iron.CellML()
cellml.CreateStart(userNumber.getInc(), region)
cellml.ModelImport(cellmlFile)

known = ["CRU/iCa",
         "CRU/ryrDensity",
         "CRU/timelag"]
for var in known:
    cellml.VariableSetAsKnown(cellmlModelIndex, var)

wanted = ["CRU/Jryr",
          "FluoBuffer/Jfluo",
          "TnCBuffer/Jtnc",
          "ATPBuffer/JATP",
          "CaMBuffer/JCaM"]
for var in wanted:
    cellml.VariableSetAsWanted(cellmlModelIndex, var)

cellml.CreateFinish()

# Field mapping
cellml.FieldMapsCreateStart()
fieldToCellmlMaps = [(dependentFields[equationLabels.index('Ca')], "CRU/Ca_free"),
                     (sourceFields[equationLabels.index('Ca')], "CRU/iCa"),
                     (RyRDensityField, "CRU/ryrDensity"),
                     (RyRTimelagField, "CRU/timelag"),
                     (dependentFields[equationLabels.index('F')], "FluoBuffer/Fluo_free"),
                     (dependentFields[equationLabels.index('FCa')], "FluoBuffer/FluoCa"),
                     (CaTnCField, "TnCBuffer/CaTnC"),
                     (dependentFields[equationLabels.index('CaM')], "CaMBuffer/CaM_free"),
                     (dependentFields[equationLabels.index('CaMCa')], "CaMBuffer/CaMCa"),
                     (dependentFields[equationLabels.index('ATP')], "ATPBuffer/ATP_free"),
                     (dependentFields[equationLabels.index('ATPCa')], "ATPBuffer/ATPCa")]

for pair in fieldToCellmlMaps:
    cellml.CreateFieldToCellMLMap(pair[0],
                                  iron.FieldVariableTypes.U, 1,
                                  iron.FieldParameterSetTypes.VALUES,
                                  cellmlModelIndex,
                                  pair[1],
                                  iron.FieldParameterSetTypes.VALUES)
    cellml.CreateCellMLToFieldMap(cellmlModelIndex,
                                  pair[1],
                                  iron.FieldParameterSetTypes.VALUES,
                                  pair[0],
                                  iron.FieldVariableTypes.U, 1,
                                  iron.FieldParameterSetTypes.VALUES)

cellml.FieldMapsCreateFinish()

# Models field
cellmlModelsField = iron.Field()
cellml.ModelsFieldCreateStart(userNumber.getInc(), cellmlModelsField)
cellml.ModelsFieldCreateFinish()
cellmlModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1)
cellmlModelsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                          iron.FieldParameterSetTypes.VALUES)
cellmlModelsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                           iron.FieldParameterSetTypes.VALUES)

# State field
cellmlStateField = iron.Field()
cellml.StateFieldCreateStart(userNumber.getInc(), cellmlStateField)
cellml.StateFieldCreateFinish()
cellmlStateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)
cellmlStateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                          iron.FieldParameterSetTypes.VALUES)

# Intermediate field
cellmlIntermediateField = iron.Field()
cellml.IntermediateFieldCreateStart(userNumber.getInc(), cellmlIntermediateField)
cellml.IntermediateFieldCreateFinish()
cellmlIntermediateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES)
cellmlIntermediateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                                 iron.FieldParameterSetTypes.VALUES)

# Parameters field
cellmlParametersField = iron.Field()
cellml.ParametersFieldCreateStart(userNumber.getInc(), cellmlParametersField)
cellml.ParametersFieldCreateFinish()
cellmlParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
cellmlParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

# E q u a t i o n s
# -----------------
if computationalNodeNumber == 0:
    print('Setting up equations for each buffer...')

for equationsSet in equationsSets:
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
    equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
    equationsSet.EquationsCreateFinish()

# P r o b l e m
# -------------
if computationalNodeNumber == 0:
    print('Setting up problem...')

problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                        iron.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                        iron.ProblemSubtypes.CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT]

problem = iron.Problem()
problem.CreateStart(userNumber.getInc(), problemSpecification)
problem.CreateFinish()

problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.TimesSet(startTime, endTime, pdeTimestep)
controlLoop.TimeOutputSet(outputFrequency)
controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.NONE)
problem.ControlLoopCreateFinish()

# S o l v e r s
# -------------
problem.SolversCreateStart()

# CellML solver
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
solver.DAETimeStepSet(odeTimestep)
solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

# PDE solver
solver = iron.Solver()
linearSolver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 2, solver)
solver.DynamicThetaSet([1.0])
solver.OutputTypeSet(iron.SolverOutputTypes.NONE)
solver.DynamicLinearSolverGet(linearSolver)
linearSolver.LinearIterativeMaximumIterations = 1000
linearSolver.linearIterativeAbsoluteTolerance = 1.0e-10
linearSolver.linearIterativeRelativeTolerance = 1.0e-8

# CellML solver
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 3, solver)
solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
solver.DAETimeStepSet(odeTimestep)
solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

problem.SolversCreateFinish()

# C e l l M L   e n v i r o n m e n t
# -----------------------------------
problem.CellMLEquationsCreateStart()
for solverIndex in (1, 3):
    solver = iron.Solver()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], solverIndex, solver)
    cellmlEquations = iron.CellMLEquations()
    solver.CellMLEquationsGet(cellmlEquations)
    cellmlEquations.CellMLAdd(cellml)
problem.CellMLEquationsCreateFinish()

# P D E   S o l v e r   e q u a t i o n s
# ---------------------------------------
problem.SolverEquationsCreateStart()
solver = iron.Solver()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 2, solver)
solverEquations = iron.SolverEquations()
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
for equationsSet in equationsSets:
    solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()


# B o u n d a r y   C o n d i t i o n s
# -------------------------------------
# Don't need to set BCs- source terms from the CellML model
# will be injecting Ca into the domain. Sarcolemma and mitochondrial
# boundaries will default to a no flux condition.
if computationalNodeNumber == 0:
    print('Setting up boundary conditions...')
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
solverEquations.BoundaryConditionsCreateFinish()


# S o l v e    p r o b le m
# -------------------------
gc.collect()
if computationalNodeNumber == 0:
    print('Solving problem...')

start = time.time()
problem.Solve()
end = time.time()

if computationalNodeNumber == 0:
    print('Success!')
    print('Solve time: ' + str(end - start))

iron.Finalise()
