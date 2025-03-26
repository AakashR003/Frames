import main
import pytest
import numpy as np
from main import Node, Member, NeumanBC, Model, GlobalResponse, MemberResponse, SecondOrderGlobalResponse
from FiniteElementDivisor import divide_into_finite_elements


@pytest.fixture
def setup_model():
    main.FEDivision = 20
    Points = [
            Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Hinged Support"),
            Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
            Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Hinged Support")
            ]


    Members = [
            Member(Beam_Number=1, Start_Node=Points[0], End_Node=Points[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
            Member(Beam_Number=2, Start_Node=Points[1], End_Node=Points[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
            ] # square cross section - 0.3 x 0.3, units N, m


    Loads = [
            NeumanBC(type="PL", Magnitude=100000, Distance1=2.5, AssignedTo="Member 2", Members = Members)
            ] 

    PointsT, MembersT, LoadsT = divide_into_finite_elements(PointsT, MembersT, LoadsT, 20)
    SecondOrderResponse = SecondOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)
    return SecondOrderResponse


def BendingMomentAll(setup_model): #check the bending moment of all members
    
    SecondOrderResponseT = setup_model

    BendingMomentall = SecondOrderResponseT.BucklingEigenLoad()[0]
    BendingMomentallT = 174.45

    assert np.allclose(BendingMomentall, BendingMomentallT, atol=0.05), "Eigen value is wrong."