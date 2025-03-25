
import pytest
import numpy as np
from main import Node, Member, NeumanBC, Model, GlobalResponse, MemberResponse, SecondOrderGlobalResponse
from FiniteElementDivisor import divide_into_finite_elements


@pytest.fixture
def setup_model():
    PointsT = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Hinged Support"),
            Node(Node_Number=2,xcoordinate=0,ycoordinate=5,Support_Condition="Rigid Joint"),
            Node(Node_Number=3,xcoordinate=5,ycoordinate=10,Support_Condition="Rigid Joint"),
            Node(Node_Number=4,xcoordinate=5,ycoordinate=0,Support_Condition="Hinged Support"),
            ] 

    MembersT = [Member(Beam_Number=1,Start_Node=PointsT[0],End_Node=PointsT[1],Area=0.09,Youngs_Modulus=200000000,Moment_of_Inertia=0.000675),
            Member(Beam_Number=2,Start_Node=PointsT[1],End_Node=PointsT[2],Area=0.09,Youngs_Modulus=200000000,Moment_of_Inertia=0.000675),
            Member(Beam_Number=3,Start_Node=PointsT[2],End_Node=PointsT[3],Area=0.09,Youngs_Modulus=200000000,Moment_of_Inertia=0.000675),
            ]

    LoadsT = [#NeumanBC(type="PL",Magnitude=5,Distance1=5,AssignedTo="Member 1", Members = Members),
            NeumanBC(type="UDL",Magnitude=5,Distance1=0, Distance2 = 5, AssignedTo="Member 2", Members = MembersT),
            NeumanBC(type="UDL",Magnitude=20,Distance1=0, Distance2 = 3, AssignedTo="Member 1", Members = MembersT),
            NeumanBC(type="PL",Magnitude=10,Distance1=1, AssignedTo="Member 3", Members = MembersT)
            ]

    PointsT, MembersT, LoadsT = divide_into_finite_elements(PointsT, MembersT, LoadsT, 20)
    SecondOrderResponse = SecondOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)
    return SecondOrderResponse


def test_FirstEigenLoad(setup_model):
    
    SecondOrderResponseT = setup_model

    EigenValueT = SecondOrderResponseT.BucklingEigenLoad()[0]
    EigenValueR = 174.45

    assert np.allclose(EigenValueT, EigenValueR, atol=1), "Eigen value is wrong."