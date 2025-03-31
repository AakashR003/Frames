
from main import Model
from StructuralElements import Node, Member
from Loads import NeumanBC
from Sensitivity import Senstivity

# Nodes for Howe Truss
Points = [
    Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Hinged Support"),
    Node(Node_Number=2, xcoordinate=5, ycoordinate=0, Support_Condition="Hinge Joint"),
    Node(Node_Number=3, xcoordinate=10, ycoordinate=0, Support_Condition="Hinge Joint"),
    Node(Node_Number=4, xcoordinate=15, ycoordinate=0, Support_Condition="Hinge Joint"),
    Node(Node_Number=5, xcoordinate=20, ycoordinate=0, Support_Condition="Hinged Support"),
    Node(Node_Number=6, xcoordinate=2.5, ycoordinate=5, Support_Condition="Hinge Joint"),
    Node(Node_Number=7, xcoordinate=7.5, ycoordinate=5, Support_Condition="Hinge Joint"),
    Node(Node_Number=8, xcoordinate=12.5, ycoordinate=5, Support_Condition="Hinge Joint"),
    Node(Node_Number=9, xcoordinate=17.5, ycoordinate=5, Support_Condition="Hinge Joint"),
]

# Members for Howe Truss
Members = [
    # Bottom Chords
    Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=3, Start_Node=Points[2], End_Node=Points[3], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=4, Start_Node=Points[3], End_Node=Points[4], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    
    # Top Chords
    Member(Beam_Number=5, Start_Node=Points[5], End_Node=Points[6], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=6, Start_Node=Points[6], End_Node=Points[7], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=7, Start_Node=Points[7], End_Node=Points[8], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    
    # Vertical Members
    Member(Beam_Number=8, Start_Node=Points[0], End_Node=Points[5], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=9, Start_Node=Points[1], End_Node=Points[6], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=10, Start_Node=Points[2], End_Node=Points[7], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=11, Start_Node=Points[3], End_Node=Points[8], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=12, Start_Node=Points[4], End_Node=Points[8], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    
    # Diagonal Members (Right Side)
    Member(Beam_Number=16, Start_Node=Points[5], End_Node=Points[1], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=17, Start_Node=Points[6], End_Node=Points[2], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=18, Start_Node=Points[7], End_Node=Points[3], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    Member(Beam_Number=19, Start_Node=Points[8], End_Node=Points[4], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    
    # Additional Members for Stability
    Member(Beam_Number=20, Start_Node=Points[5], End_Node=Points[7], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
]

# Loads for Howe Truss
Loads = [
    NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=5, AssignedTo="Member 1", Members=Members),
    NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=5, AssignedTo="Member 2", Members=Members),
    NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=5, AssignedTo="Member 3", Members=Members),
    NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=5, AssignedTo="Member 4", Members=Members),
]

Model1 = Model(Points = Points, Members = Members, Loads = Loads)
#GlobalRes1 = GlobalResponse(Points = Points, Members = Members, Loads = Loads)
#NodalRes1 = NodalResponse(Points = Points, Members = Members, Loads = Loads)
#MemberRes1 = MemberResponse(Points = Points, Members = Members, Loads = Loads)
Sensitivity1 = Senstivity(Points = Points, Members = Members, Loads = Loads)

Model1.PlotGlobalModel()
#print(MemberRes.PlotMemberBMD(2))
print(Sensitivity1.NodeYSensitivity(1,0.001))
print(Sensitivity1.GlobalSizeSensitivity("Axial"))

Sensitivity1.PlotSensitivity("Axial")
