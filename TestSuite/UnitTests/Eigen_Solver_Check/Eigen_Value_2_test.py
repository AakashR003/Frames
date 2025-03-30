
import pytest
import numpy as np

from config import config
from main import Model
from StructuralElements import Node, Member
from Loads import NeumanBC
from SecondOrderResponse import  SecondOrderGlobalResponse


@pytest.fixture
def setup_model():
    config.set_FEDivision(1000)
    PointsT = [
        Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
        Node(Node_Number=2, xcoordinate=0, ycoordinate=1, Support_Condition="Rigid Joint"),
        Node(Node_Number=3, xcoordinate=0, ycoordinate=2, Support_Condition="Rigid Joint"),
        Node(Node_Number=4, xcoordinate=0, ycoordinate=3, Support_Condition="Rigid Joint"),
        Node(Node_Number=5, xcoordinate=0, ycoordinate=4, Support_Condition="Rigid Joint"),
        Node(Node_Number=6, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=7, xcoordinate=1, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=8, xcoordinate=2, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=9, xcoordinate=3, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=10, xcoordinate=4, ycoordinate=5, Support_Condition="Rigid Joint"),
        Node(Node_Number=11, xcoordinate=5, ycoordinate=5, Support_Condition="Fixed Support")
    ]
    MembersT = [
        Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=3, Start_Node=PointsT[2], End_Node=PointsT[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=4, Start_Node=PointsT[3], End_Node=PointsT[4], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=5, Start_Node=PointsT[4], End_Node=PointsT[5], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=6, Start_Node=PointsT[5], End_Node=PointsT[6], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=7, Start_Node=PointsT[6], End_Node=PointsT[7], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=8, Start_Node=PointsT[7], End_Node=PointsT[8], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=9, Start_Node=PointsT[8], End_Node=PointsT[9], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
        Member(Beam_Number=10, Start_Node=PointsT[9], End_Node=PointsT[10], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675)
    ] # square cross section - 0.3 x 0.3, units N, m
    LoadsT = [
        NeumanBC(type="UDL", Magnitude=100, Distance1=0, Distance2=1, AssignedTo="Member 6", Members = MembersT),
        NeumanBC(type="UDL", Magnitude=100, Distance1=0, Distance2=1, AssignedTo="Member 7", Members = MembersT),
        NeumanBC(type="UDL", Magnitude=100, Distance1=0, Distance2=1, AssignedTo="Member 8", Members = MembersT),
        NeumanBC(type="UDL", Magnitude=100, Distance1=0, Distance2=1, AssignedTo="Member 9", Members = MembersT),
        NeumanBC(type="UDL", Magnitude=100, Distance1=0, Distance2=1, AssignedTo="Member 10", Members = MembersT),
        #NeumanBC(type="PL", Magnitude=100, Distance1=1, AssignedTo="Member 4", Members = Members)
    ]

    ModelT = Model(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #ResT = GlobalResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #MemberResponseT = MemberResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    SecondOrderResponseT = SecondOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)

    return SecondOrderResponseT

def test_NormalForce(setup_model): #test the normal force which is the reason for reduction 2nd order matrix
    SecondOrderResponseT = setup_model

    NormalForceT = SecondOrderResponseT.NormalForce()
    NormalForceR = [218.3,218.3,218.3,218.3,218.3,30.8,30.8,30.8,30.8,30.8]

    assert np.allclose(NormalForceT, NormalForceR, atol=1e-1), "Normal force is wrong."


def test_FirstEigenLoad(setup_model):
    SecondOrderResponseT = setup_model

    EigenValueT = SecondOrderResponseT.BucklingEigenLoad()[0]
    EigenValueR = 670

    assert np.allclose(EigenValueT, EigenValueR, atol=15), "Eigen value is wrong."