
from main import Node, Member, NeumanBC, Model, GlobalResponse, MemberResponse, NodalResponse, SecondOrderGlobalResponse, Senstivity
from FiniteElementDivisor import divide_into_finite_elements
from Functions import print_class_Objects


"""
#Model Parts - Basic essential for building a model
Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Hinged Support"),
        Node(Node_Number=2,xcoordinate=10,ycoordinate=0,Support_Condition="Hinged Support"),
        Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Hinged Support"),
        ] 

Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
          Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
           ]

Loads = [NeumanBC(type="PL",Magnitude=5,Distance1=5,AssignedTo="Member 1", Members = Members),
        NeumanBC(type="UDL",Magnitude=5,Distance1=0, Distance2 = 10, AssignedTo="Member 2", Members = Members)
        ] 
"""

#Model Parts - Basic essential for building a model
Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Hinged Support"),
         Node(Node_Number=2,xcoordinate=0,ycoordinate=5,Support_Condition="Rigid Joint"),
         Node(Node_Number=3,xcoordinate=5,ycoordinate=10,Support_Condition="Rigid Joint"),
         Node(Node_Number=4,xcoordinate=5,ycoordinate=0,Support_Condition="Hinged Support"),
        ] 

Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=0.09,Youngs_Modulus=200000000,Moment_of_Inertia=0.000675),
           Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=0.09,Youngs_Modulus=200000000,Moment_of_Inertia=0.000675),
           Member(Beam_Number=3,Start_Node=Points[2],End_Node=Points[3],Area=0.09,Youngs_Modulus=200000000,Moment_of_Inertia=0.000675),
           ]

Loads = [#NeumanBC(type="PL",Magnitude=5,Distance1=5,AssignedTo="Member 1", Members = Members),
        NeumanBC(type="UDL",Magnitude=5,Distance1=0, Distance2 = 5, AssignedTo="Member 2", Members = Members),
        NeumanBC(type="UDL",Magnitude=20,Distance1=0, Distance2 = 3, AssignedTo="Member 1", Members = Members),
        NeumanBC(type="PL",Magnitude=10,Distance1=1, AssignedTo="Member 3", Members = Members)
        ]


Points, Members, Loads = divide_into_finite_elements(Points, Members, Loads, 20)


#main Model part - Main mode part includes sub model part
Model1 = Model(Points = Points, Members = Members, Loads = Loads)
GlobalRes1 = GlobalResponse(Points = Points, Members = Members, Loads = Loads)
NodalRes1 = NodalResponse(Points = Points, Members = Members, Loads = Loads)
MemberRes1 = MemberResponse(Points = Points, Members = Members, Loads = Loads)
Sensitivity1 = Senstivity(Points = Points, Members = Members, Loads = Loads)
SecondOrderResponse1 = SecondOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)


Model1.PlotGlobalModel()
print("mem1",MemberRes1.MemberForceGlobal(1))
print("mem2",MemberRes1.MemberForceGlobal(2))
print("mem1",MemberRes1.MemberForceLocal(1))
print("mem2",MemberRes1.MemberForceLocal(2))
print(SecondOrderResponse1.BucklingEigenLoad())
#MemberRes1.PlotMemberBMD(1)