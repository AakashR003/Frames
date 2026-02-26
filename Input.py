from config import config 
from Model import Model
from StructuralElements import Node, Couple_Nodes, Member
from Loads import NeumanBC
from FirstOrderResponse import FirstOrderGlobalResponse, FirstOrderMemberResponse, FirstOrderNodalResponse
from SecondOrderResponse import  SecondOrderGlobalResponse, SecondOrderMemberResponse
from DynamicResponse import DynamicGlobalResponse
from Comparision import Comparision
from Sensitivity import Senstivity, SecondOrderSensitivity, FiniteDifferenceSensitivity, ComparisionSensitivity
from ApproximatedSecondOrderAnalysis import ApproximatedSecondOrderAnalysis, ApproximatedAnalysisDisplacement
from Strain_Energy import StrainEnergy
from FiniteElementDivisor import divide_into_finite_elements
from Functions import print_class_Objects

import numpy as np




#Model Parts - Basic essential for building a model
config.set_FEDivision(1000)

"""
Points = [
Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Hinged Support"),
#Node(Node_Number=2, xcoordinate=5, ycoordinate=0, Support_Condition="Hinge Joint"),
Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=4, xcoordinate=5, ycoordinate=0, Support_Condition="Hinged Support")
]

#Coupling = [
#Couple_Nodes(Main_Node=Points[1], Dependent_Node=Points[2], xDof=True, yDof=True, RotationDof=False),
#]


Members = [
Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=3, Start_Node=Points[2], End_Node=Points[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=4, Start_Node=Points[0], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
] # square cross section - 0.3 x 0.3, units N, m


Loads = [
#NeumanBC(type="UDL", Magnitude=10, Distance1= 2, Distance2= 6, AssignedTo="Member 1", Members = Members),
NeumanBC(type="PL", Magnitude=-79600, Distance1= 2.5, AssignedTo="Member 2", Members = Members),
#NeumanBC(type="NL", Magnitude=-10, AssignedTo="Node 2", Members = Members, Nodes = Points)
]

"""

Points = [
Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
#Node(Node_Number=2, xcoordinate=5, ycoordinate=0, Support_Condition="Hinge Joint"),
Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=4, xcoordinate=5, ycoordinate=0, Support_Condition="Hinged Support")
]

#Coupling = [
#Couple_Nodes(Main_Node=Points[1], Dependent_Node=Points[2], xDof=True, yDof=True, RotationDof=False),
#]


Members = [
Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=3, Start_Node=Points[2], End_Node=Points[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=4, Start_Node=Points[0], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
] # square cross section - 0.3 x 0.3, units N, m


Loads = [
#NeumanBC(type="UDL", Magnitude=10, Distance1= 2, Distance2= 6, AssignedTo="Member 1", Members = Members),
NeumanBC(type="PL", Magnitude=-60600, Distance1= 2.5, AssignedTo="Member 2", Members = Members),
#NeumanBC(type="NL", Magnitude=-10, AssignedTo="Node 2", Members = Members, Nodes = Points)
]

Points, Members, Loads = divide_into_finite_elements(Points, Members, Loads, 7)
#"""


"""
Points = [
Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
#Node(Node_Number=2, xcoordinate=5, ycoordinate=0, Support_Condition="Hinge Joint"),
Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Fixed Support"),
#Node(Node_Number=4, xcoordinate=5, ycoordinate=0, Support_Condition="Hinged Support")
]

#Coupling = [
#Couple_Nodes(Main_Node=Points[1], Dependent_Node=Points[2], xDof=True, yDof=True, RotationDof=False),
#]


Members = [
Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
#Member(Beam_Number=3, Start_Node=Points[2], End_Node=Points[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
] # square cross section - 0.3 x 0.3, units N, m


Loads = [
#NeumanBC(type="UDL", Magnitude=10, Distance1= 2, Distance2= 6, AssignedTo="Member 1", Members = Members),
NeumanBC(type="PL", Magnitude=-256000, Distance1= 2.5, AssignedTo="Member 2", Members = Members),
#NeumanBC(type="NL", Magnitude=-10, AssignedTo="Node 2", Members = Members, Nodes = Points)
]

Points, Members, Loads = divide_into_finite_elements(Points, Members, Loads, 7)
"""






#main Model part - Main mode part includes sub model part
Model1 = Model(Points = Points, Members = Members, Loads = Loads)
GlobalRes1 = FirstOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)
NodalRes1 = FirstOrderNodalResponse(Points = Points, Members = Members, Loads = Loads)
MemberRes1 = FirstOrderMemberResponse(Points = Points, Members = Members, Loads = Loads)
SecondOrderResponse1 = SecondOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)
SecondOrderMemberResponse1 = SecondOrderMemberResponse(Points = Points, Members = Members, Loads = Loads)
Comparision1 = Comparision(MainModel = MemberRes1, Model2 = SecondOrderMemberResponse1)
DynamicResponse1 = DynamicGlobalResponse(Points = Points, Members = Members, Loads = Loads)
SecondOrderSensitivity1 = SecondOrderSensitivity(Points = Points, Members = Members, Loads = Loads)
Senstivity1 = Senstivity(Points = Points, Members = Members, Loads = Loads)
ComparisionSensitivity1 = ComparisionSensitivity(Points = Points, Members = Members, Loads = Loads)
ApproximatedSecondOrderAnalysis1 = ApproximatedSecondOrderAnalysis(Points = Points, Members = Members, Loads = Loads)
ApproximatedAnalysisDisplacement1 = ApproximatedAnalysisDisplacement(Points = Points, Members = Members, Loads = Loads )
StrainEnergy1 = StrainEnergy(Points = Points, Members = Members, Loads = Loads)
FiniteDifferenceSensitivity1 = FiniteDifferenceSensitivity(Points = Points, Members = Members, Loads = Loads)


print("unconstrained dof",Model1.UnConstrainedDoF())
print("constrained dof",Model1.ConstrainedDoF())


Model1.PlotGlobalModel()

SecondOrderResponse1.PlotEigenMode(EigenModeNo = 1, Solver="eigsh", scale_factor = 1)
#FiniteDifferenceSensitivity1.Members_All_Sensitivity()
ComparisionSensitivity1.AnlSimpson_Vs_FDSimpson()
#ComparisionSensitivity1.AnlSimpson_Vs_FDNonLinear()


#print("FD_1_1",FiniteDifferenceSensitivity1.FiniteDifferenceFirstOrderLinearBendingSensitivity(MemberNumber=7, scale=0.000000000001))
#print("FD_2_1",FiniteDifferenceSensitivity1.FiniteDifferenceSecondOrderLinearBendingSensitivity(MemberNumber=7, scale=0.000000000001))
#print("FD_2_1.1_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",FiniteDifferenceSensitivity1.FiniteDifferenceApproximatedSecondOrderBendingSensitivity(MemberNumber=7, scale=0.000000000001))
#print("FD_2_1.1_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",FiniteDifferenceSensitivity1.FiniteDifferenceSimpsonsApproximatedSecondOrderBendingSensitivity(MemberNumber=7, scale=0.000000000001))
#print("FD_2_2_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",FiniteDifferenceSensitivity1.FiniteDifferenceSecondOrderNonLinearBendingSensitivity(MemberNumber=7, scale=0.000000000001))
#LinearSensitivity = Senstivity1.BendingMemberSensitivity(MemberNumber=7, scale=0.000000000001)
#print("Anl_1_1", Senstivity1.BendingMemberSensitivity(MemberNumber=3, scale=0.000000000001))
#print("Anl_2_1_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",SecondOrderSensitivity1.GlobalSecondOrderBendingSensitivity(MemberNumber=7, scale=0.000000000001))
#print("Anl_2_2.1_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",SecondOrderSensitivity1.try2SecondOrderPartSensitivity(MemberNumber=9, scale = 0.000000000001))
#print("Anl_2_2.1_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",SecondOrderSensitivity1.SimposonsApproximatedSecondOrderSensitivity(MemberNumber=7, scale = 0.000000000001))

#print("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB", SecondOrderSensitivity1.try3(7, LinearSensitivity, scale=0.000000000001))
#sum = (Senstivity1.BendingMemberSensitivity(MemberNumber=3, scale=0.000000000001)  + SecondOrderSensitivity1.SecondOrderPartBendingSensitivity(MemberNumber=3, scale=0.000000000001))

#print("FinalAnswer",sum)
#GlobalRes1.PlotEigenMode(EigenModeNo = 42, Solver="eigsh", scale_factor = 1)
#print(GlobalRes1.CalculateModalStiffness())
#GlobalRes1.OrthogonalForceVector()

#print(np.round(ApproximatedSecondOrderAnalysis1.CalculateModifiedOrthogonalStiffnessDifferenceMatrix().real, 2))
#ApproximatedSecondOrderAnalysis1.checkcorrectness()
#ApproximatedSecondOrderAnalysis1.PlotSecondOrderLoadDisplacementCurve(NodeNumber = 12, Direction = "y", division =10)
#ApproximatedAnalysisDisplacement1.PlotSecondOrderLoadDisplacementCurve(NodeNumber = 9, Direction = "y", division =10)

#Comparision1.PlotLoadDisplacementCurveComparison(NodeNumber = 15, Direction = "y", division =20)
#ApproximatedSecondOrderAnalysis1.SetModifiedValues()
#print(ApproximatedSecondOrderAnalysis1.ApproximatedSecondOrderDisplacementLocal())
#print(ApproximatedSecondOrderAnalysis1.CalculateApproximatedValueDisplacement())

#print("Linear Strain Energy",StrainEnergy1.CalculateLinearStrainEnergy())
#print("Linearized Strain Energy",StrainEnergy1.CalculateLinearizedStrainEnergy())
#print("Approximated NonLinear Strain Energy",StrainEnergy1.CalculateApproximatedNonlinearStrainEnergy())
#print("Simposons Approximated Strain Energy", StrainEnergy1.CalculateSimpsonsApproximatedNonLinearStrainEnergy())
#print("Finite DIfference NonLinear Strain Energy",StrainEnergy1.CalculateFiniteDifferenceNonLinearStrainEnergy())

#print("Completed Necessary")
#print("Starting Additional")



#print("dof number",Members[1].DoFNumber())
#print("dof number",Members[2].DoFNumber())
#print("Force Vector",Model1.ForceVector())

#MemberRes1.PlotMemberBMD(2)
#MemberRes1.PlotGlobalSFD()
#MemberRes1.PlotGlobalDeflection()
#print(Model1.ForceVector())
#print(GlobalRes1.DisplacementVectorDict())
#print(MemberRes1.MemberForceLocal(1,All = True))
#MemberRes1.PlotMemberBMD(1)
#MemberRes1.PlotMemberSFD(1)
#print(MemberRes1.MemberForceLocal(1, All = True))
#print(SecondOrderMemberResponse1.MemberForceLocal(1, All = True))
#MemberRes1.PlotGlobalBMD(show_structure=True, scale_factor=2)
#MemberRes1.PlotGlobalSFD(show_structure=True, scale_factor=2)


#print(SecondOrderResponse1.BucklingEigenLoad())
#mf = SecondOrderMemberResponse1.MemberForceLocal(1, All = True)
#print(SecondOrderResponse1.BucklingEigenLoad()[0])
#SecondOrderMemberResponse1.PlotGlobalDeflection()
#SecondOrderMemberResponse1.PlotMemberDeflection(1)
#SecondOrderMemberResponse1.PlotMemberBMD(1)
#SecondOrderMemberResponse1.PlotGlobalSFD()
#SecondOrderMemberResponse1.PlotMemberBMD(2)
#SecondOrderMemberResponse1.PlotMemberBMD(1)
#SecondOrderMemberResponse1.PlotGlobalBMD(show_structure=True, scale_factor=2)
#SecondOrderMemberResponse1.PlotGlobalSFD(show_structure=True)
#print(SecondOrderResponse1.BucklingEigenLoad())
#SecondOrderResponse1.PlotEigenMode(EigenModeNo = 3, Solver="eigsh", scale_factor = 1)
SecondOrderResponse1.PlotSecondOrderLoadDisplacementCurve(NodeNumber = 12, Direction = "y", LoadFactor = None, division =10)
#print(SecondOrderResponse1.SecondOrderDisplacementVector(iteration_steps=5))
#SecondOrderMemberResponse1.PlotMemberBMD(1)
#SecondOrderMemberResponse1.PlotGlobalBMD(show_structure=True)
#mf = SecondOrderMemberResponse1.MemberForceLocal(1, All = True)
#print("MemNo",MemNo,SecondOrderMemberResponse1.MemberBMD(MemNo, mf[i]))
#print(SecondOrderResponse1.MemberEigenMode(11, EigenModeNo = 1, scale_factor = 1000000))

#print("EigenFrequency", DynamicResponse1.EigenFrequency())
#DynamicResponse1.PlotDynamicEigenMode(1)


#Comparision1.PlotGlobalBMDComparison()
#Comparision1.PlotGlobalDeflectionComparison(scale_factor = 0.75)
Comparision1.PlotLoadDisplacementCurveComparison(NodeNumber = 12, Direction = "y", division =20)
