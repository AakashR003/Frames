# FEM for Frame Elements

This repository contains Python files. You can either clone this repository and use the software locally, or:

1. Use the web interface: https://structuralanalysis.streamlit.app
2. Install it as a PyPI module: https://pypi.org/project/NStructAnaly

## 1.Introduction:

This is a Python Module for performing Linear, Nonlinear(2nd Order), Dynamic and Sensitivity Analysis of plane frames and Truss. You can use the package to calculate the following:

1. Member Bending Moment 
2. Member Shear Force
3. Member Normal Force
4. Nodal Reactions Fx, Fy, My
5. Nodal Displacement and Rotations

All of the above results can be obtained for both 1st order and 2nd order Analysis. In addtion you can also obtain following advance results:

1. Critical Buckling Eigen Value
2. Critical Buckling Eigen Mode
3. Frequency Eigen Value
4. Frequency Eigen Mode
5. Node Sensitivity
6. Axial Stiffness Sensitivity
7. Bending Sensitivity

You also have an options  

1. Dividing Strutures Individual Finite Elements (Reduces working time on creating Nodes and elements)
2. creating comparitive plot between two models (Presently its available only for Bending Moment)

## 2.How to Use:

1. Create and Assemble the Structural Model
2. Create Responses
3. Perform the Individual Response

The package contains following Modules:

1. StructuralElements - With this you can create Node, Member
2. Loads - With this you can create NeumanBC
3. main - with this you can assembale the model
4. FirstOrderResponse - With this you can create FirstOrderGLobalResponse and FirstOrderMemberResponse
5. SecondOrderResponse - With this you can create SecondOrderGLobalResponse and SecondOrderMemberResponse
6. DynamicResponse - With this you can create DynamicGlobalResponse
7. FiniteElementDivisor and Comparision




### 2.1 How to create Strucutural Model

Based on the Modules gives above use the following variables explained below to create strucutre

1. **Node:**  
   - **Node_Number** – Give the node number as an integer.  
   - **xcoordinate** – Specify the x-coordinate as a float.  
   - **ycoordinate** – Specify the y-coordinate as a float.  
   - **Support_Condition** – Define the type of support as a string from the given options-
   "Hinged Support", "Fixed Support", "Rigid Joint", "Roller in X-plane", "Roller in Y-plane", "Hinge Joint", "Glided Support", "Roller in X-plane-Hinge"

2. **Member:**  
   - **Beam_Number** – Assign a unique beam number as an integer.  
   - **Start_Node** – Specify the starting node number as an integer.  
   - **End_Node** – Specify the ending node number as an integer.  
   - **Area** – Define the cross-sectional area as a float.  
   - **Youngs_Modulus** – Specify Young’s modulus as a float.  
   - **Moment_of_Inertia** – Provide the moment of inertia as a float.  

3. **NeumanBC:**  
   - **type** – Define the Load type as a string.  
                "UDL" - Uniformly Distributed Load, "PL" - Point Load
   - **Magnitude** – Specify the magnitude of the Load as a float.  
   - **Distance1** – Provide the distance at which the Load is starts as a float.  
   - **Distance2** – Provide the distance at which the Load end as a float. ( Not needed if you use "PL")
   - **AssignedTo** – Indicate the member to which the load is applied as a string.  
   - **Members** – List the member numbers assigned to this boundary condition as integers.

Example:

Points = [
StructuralElements.Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Hinged Support"),
StructuralElements.Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
StructuralElements.Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Rigid Joint"),
StructuralElements.Node(Node_Number=4, xcoordinate=5, ycoordinate=0, Support_Condition="Hinged Support")
]


Members = [
StructuralElements.Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
StructuralElements.Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
StructuralElements.Member(Beam_Number=3, Start_Node=Points[2], End_Node=Points[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
] # square cross section - 0.3 x 0.3, units N, m


Loads = [
Loads.NeumanBC(type="PL", Magnitude=100000, Distance1= 2.5, AssignedTo="Member 2", Members = Members)
] 

This creates a Square frame with point Load on Middle


### 2.2 How to create responses:

Response is created based on the Structure. Example of creating each response is given

Model1 = main.Model(Points = Points, Members = Members, Loads = Loads)

GlobalRes1 = FirstOrderResponse.FirstOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)

NodalRes1 = NodalResponse(Points = Points, Members = Members, Loads = Loads)

MemberRes1 = FirstOrderResponse.FirstOrderMemberResponse(Points = Points, Members = Members, Loads = Loads)

SecondOrderResponse1 = SecondOrderResponse.SecondOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)

SecondOrderMemberResponse1 = SecondOrderResponse.SecondOrderMemberResponse(Points = Points, Members = Members, Loads = Loads)

Comparision1 = Comparision.Comparision(MainModel = MemberRes1, Model2 = SecondOrderMemberResponse1)

DynamicResponse1 = DynamicResponse.DynamicGlobalResponse(Points = Points, Members = Members, Loads = Loads)


Instead of Model1, GlobalResponse1 you can use anyother name. 
Based on this model Structural results can be obtained as shown below

Follwing Structural results are availabe for each response:


**1. FirstOrderGlobalResponse** - SupportForcesVector, DisplacementVectorDict

**2. FirstOrderNodalResponse** - NodeForce, NodeDisplacement

**3. FirstOrderMemberResponse** - MemberForceLocal, MemberBMD, MemberSFD PlotMemberBMD, PlotMemberSFD, PlotGlobalBMD, PlotGlobalSFD, PlotGlobalDeflection

**4. SecondOrderGlobalResponse** - SupportForcesVector, DisplacementVectorDict, BucklingEigenLoad, PlotEigenMode

**5. SecondOrderMemberResponse** - MemberForceLocal, MemberBMD, MemberSFD PlotMemberBMD, PlotMemberSFD, PlotGlobalBMD, PlotGlobalSFD, PlotGlobalDeflection

**6. Comparision** - PlotGlobalBMDComparison, PlotGlobalSFDComparison, PlotGlobalDeflectionComparison
    



