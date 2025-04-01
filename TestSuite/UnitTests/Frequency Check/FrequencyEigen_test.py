
import pytest
import numpy as np

from config import config
from main import Model
from StructuralElements import Node, Member
from Loads import NeumanBC
from DynamicResponse import DynamicGlobalResponse
from FiniteElementDivisor import divide_into_finite_elements


@pytest.fixture
def setup_model():

    # Cantilivered L-frame with Point load which can be also teated as cantilivered beam Eigen value problem
    config.set_FEDivision(20)
    #Model Parts - Basic essential for building a model
    PointsT = [
    Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
    Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Fixed Support"),
    ]


    MembersT = [
    Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=0.016, Youngs_Modulus=200000000, Moment_of_Inertia=2.13333E-07),
    ] # square cross section - 0.3 x 0.3, units N, m


    LoadsT = [
    NeumanBC(type="PL", Magnitude=-10000, Distance1= 2.5, AssignedTo="Member 1", Members = MembersT)
    ] 

    PointsT, MembersT, LoadsT = divide_into_finite_elements(PointsT, MembersT, LoadsT, 10)

    #ModelT = Model(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #ResT = GlobalResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    #MemberResponseT = MemberResponse(Points=PointsT, Members=MembersT, Loads=LoadsT)
    DynamicResponseT = DynamicGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)

    return DynamicResponseT

def test_EigenFrequency(setup_model):
    DynamicResponseT = setup_model

    EigenValueT = DynamicResponseT.EigenFrequency()[1]
    EigenValueR = [0.08, 0.23, 0.45, 0.74, 1.11, 1.56, 2.1, 2.71, 3.38, 4.5, 5.47, 6.65, 8.06, 9.73, 11.67, 13.86, 16.03, 16.09, 17.9, 32.45, 49.67, 68.08, 88.0, 109.5, 131.96, 153.38, 169.73]

    assert np.allclose(EigenValueT, EigenValueR, atol=2), "Eigen Frequency is wrong."


