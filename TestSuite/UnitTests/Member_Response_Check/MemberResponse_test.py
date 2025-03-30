
import pytest
import numpy as np

from config import config
from StructuralElements import Node, Member
from Loads import NeumanBC
from FirstOrderResponse import FirstOrderGlobalResponse, FirstOrderMemberResponse



def test_MemberForce():

    config.set_FEDivision(1000)
    PointsT = [
        Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Hinged Support"),
        Node(Node_Number=2, xcoordinate=10, ycoordinate=0, Support_Condition="Rigid Joint"),
        Node(Node_Number=3, xcoordinate=20, ycoordinate=0, Support_Condition="Hinged Support")
    ]
    MembersT = [
        Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
        Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    ]
    LoadsT = [
        NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=10, AssignedTo="Member 1", Members = MembersT),
        NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=10, AssignedTo="Member 2", Members = MembersT)
    ]
    #ModelT = Model(Points=PointsT, Members=MembersT, Loads=LoadsT)
    GlobalResponseT = FirstOrderGlobalResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    MemberResT = FirstOrderMemberResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)

    Member1Force = MemberResT.MemberForceLocal(1)
    print(GlobalResponseT.SupportForcesVector())
    print(GlobalResponseT.DisplacementVector())
    print(Member1Force)