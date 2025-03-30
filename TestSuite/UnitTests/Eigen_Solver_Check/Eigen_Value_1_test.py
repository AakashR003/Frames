
import pytest
import numpy as np

from config import config
from main import Model
from StructuralElements import Node, Member
from Loads import NeumanBC
from SecondOrderResponse import  SecondOrderGlobalResponse

@pytest.fixture
def setup_model1():

    # Cantilivered L-frame with Point load which can be also teated as cantilivered beam Eigen value problem
    config.set_FEDivision(1000)
    PointsT = [
        Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
        Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=3, xcoordinate=1, ycoordinate=5, Support_Condition="Rigid Joint")
    ]
    MembersT = [
        Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
    ] # square cross section - 0.3 x 0.3, units N, m
    LoadsT = [
        NeumanBC(type="PL", Magnitude=100, Distance1=1, AssignedTo="Member 2", Members = MembersT)
    ]

    ModelT = Model(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #ResT = GlobalResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #MemberResponseT = MemberResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    SecondOrderResponseT = SecondOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)

    return SecondOrderResponseT

def test_FirstEigenLoad_1(setup_model1):
    SecondOrderResponseT = setup_model1

    EigenValue = SecondOrderResponseT.BucklingEigenLoad()[0]
    EigenValueR = 133.0

    assert np.allclose(EigenValue, EigenValueR, atol=2), "Eigen value is wrong."


@pytest.fixture
def setup_model2():
    #Model Parts - Basic essential for building a model
    config.set_FEDivision(1000)
    PointsT = [
        Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
        Node(Node_Number=2, xcoordinate=0, ycoordinate=1, Support_Condition="Rigid Joint"),
        Node(Node_Number=3, xcoordinate=0, ycoordinate=2, Support_Condition="Rigid Joint"),
        Node(Node_Number=4, xcoordinate=0, ycoordinate=3, Support_Condition="Rigid Joint"),
        Node(Node_Number=5, xcoordinate=0, ycoordinate=4, Support_Condition="Rigid Joint"),
        Node(Node_Number=6, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=7, xcoordinate=1, ycoordinate=5, Support_Condition="Rigid Joint")
    ]
    MembersT = [
        Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=3, Start_Node=PointsT[2], End_Node=PointsT[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=4, Start_Node=PointsT[3], End_Node=PointsT[4], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=5, Start_Node=PointsT[4], End_Node=PointsT[5], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=6, Start_Node=PointsT[5], End_Node=PointsT[6], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        ] # square cross section - 0.3 x 0.3, units N, m
    LoadsT = [
        NeumanBC(type="PL", Magnitude=100, Distance1=0.001, AssignedTo="Member 6", Members = MembersT)
    ]

    ModelT = Model(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #ResT = GlobalResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #MemberResponseT = MemberResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    SecondOrderResponseT = SecondOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)

    return SecondOrderResponseT

def test_FirstEigenLoad_2(setup_model2):
    SecondOrderResponseT = setup_model2

    EigenValue = SecondOrderResponseT.BucklingEigenLoad()[0]
    EigenValueR = 133.0

    assert np.allclose(EigenValue, EigenValueR, atol=0.5), "Eigen value is wrong."