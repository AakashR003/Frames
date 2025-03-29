import main
from main import Node, Member, NeumanBC, Model, GlobalResponse, MemberResponse, NodalResponse, SecondOrderGlobalResponse, SecondOrderMemberResponse, Comparision, DynamicGlobalResponse
from FiniteElementDivisor import divide_into_finite_elements
from Functions import print_class_Objects



main.FEDivision = 20
#Model Parts - Basic essential for building a model
Points = [
Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Rigid Joint"),
Node(Node_Number=4, xcoordinate=5, ycoordinate=0, Support_Condition="Hinged Support")
]


Members = [
Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
Member(Beam_Number=3, Start_Node=Points[2], End_Node=Points[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
] # square cross section - 0.3 x 0.3, units N, m


Loads = [
NeumanBC(type="PL", Magnitude=10000, Distance1= 2.5, AssignedTo="Member 2", Members = Members)
] 





Points, Members, Loads = divide_into_finite_elements(Points, Members, Loads, 10)


#main Model part - Main mode part includes sub model part
Model1 = Model(Points = Points, Members = Members, Loads = Loads)
GlobalRes1 = GlobalResponse(Points = Points, Members = Members, Loads = Loads)
NodalRes1 = NodalResponse(Points = Points, Members = Members, Loads = Loads)
MemberRes1 = MemberResponse(Points = Points, Members = Members, Loads = Loads)
SecondOrderResponse1 = SecondOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)
SecondOrderMemberResponse1 = SecondOrderMemberResponse(Points = Points, Members = Members, Loads = Loads)
Comparision1 = Comparision(MainModel = MemberRes1, Model2 = SecondOrderMemberResponse1)
DynamicResponse1 = DynamicGlobalResponse(Points = Points, Members = Members, Loads = Loads)


Model1.PlotGlobalModel()

#MemberRes1.PlotMemberBMD(2)
#print(SecondOrderResponse1.BucklingEigenLoad()[0])

#SecondOrderMemberResponse1.PlotMemberBMD(1)
#mf = SecondOrderMemberResponse1.MemberForceLocal(1, All = True)



#MemberRes1.PlotMemberBMD(1)
#MemberRes1.PlotGlobalBMD(show_structure=True)
#print(SecondOrderResponse1.SecondOrderDisplacementVector(10))
#SecondOrderMemberResponse1.PlotMemberBMD(1)
#SecondOrderMemberResponse1.PlotGlobalBMD(show_structure=True)
#mf = SecondOrderMemberResponse1.MemberForceLocal(1, All = True)
#print("MemNo",MemNo,SecondOrderMemberResponse1.MemberBMD(MemNo, mf[i]))
#Comparision1.PlotGlobalBMDComparison()
#MemberRes1.PlotGlobalDeflection()
#print(SecondOrderResponse1.MemberEigenMode(11, EigenModeNo = 1, scale_factor = 1000000))
SecondOrderResponse1.PlotEigenMode(EigenModeNo = 1, scale_factor = 0.3)

#print("EigenFrequency", DynamicResponse1.EigenFrequency())
#DynamicResponse1.PlotDynamicEigenMode()
