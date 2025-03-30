import pytest
import numpy as np

import main
from main import Model
from StructuralElements import Node, Member
import Loads
from Loads import NeumanBC
from SecondOrderResponse import  SecondOrderGlobalResponse, SecondOrderMemberResponse
from FiniteElementDivisor import divide_into_finite_elements


@pytest.fixture
def setup_model():

    # Cantilivered L-frame with Point load which can be also teated as cantilivered beam Eigen value problem
    Loads.FEDivision = 20
    #Model Parts - Basic essential for building a model
    PointsT = [
    Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
    Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
    Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Rigid Joint"),
    Node(Node_Number=4, xcoordinate=5, ycoordinate=0, Support_Condition="Hinged Support")
    ]


    MembersT = [
    Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
    Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
    Member(Beam_Number=3, Start_Node=PointsT[2], End_Node=PointsT[3], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
    ] # square cross section - 0.3 x 0.3, units N, m


    LoadsT = [
    NeumanBC(type="PL", Magnitude=10000, Distance1= 2.5, AssignedTo="Member 2", Members = MembersT)
    ] 





    PointsT, MembersT, LoadsT = divide_into_finite_elements(PointsT, MembersT, LoadsT, 10)
    ModelT = Model(Points=PointsT, Members=MembersT, Loads=LoadsT)
    SecondOrderResponseT = SecondOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)

    return SecondOrderResponseT

def test_BucklingEigenValue(setup_model):
    SecondOrderResponseT = setup_model

    EigenValueT = SecondOrderResponseT.BucklingEigenLoad()[1]
    EigenValueR = [4.74, 15.37, 29.83, 46.43, 71.93, 95.78, 131.98, 162.91, 207.2, 241.59, 280.72, 339.62, 380.98, 448.06, 501.49, 556.3, 646.93, 699.28, 820.02, 875.9, 975.31, 1069.85, 1235.49, 1337.12, 1461.76, 1571.79, 1738.44, 1853.11, 2061.18, 2208.83, 2383.18, 2581.07, 2681.71, 2789.69, 3056.47, 3275.76, 3481.55, 3481.55, 3481.55, 3481.55, 3481.55, 3481.55, 3481.55, 3481.55, 3481.55, 3493.73, 
                   3613.4, 3726.79, 3726.79, 3726.79, 3726.79, 3726.79, 3726.79, 3726.79, 3726.79, 3726.79, 3790.07, 4047.53, 4331.9, 4624.29, 4918.97, 5197.25, 5372.03, 5502.57, 5810.15, 5954.84, 6211.48, 6454.11, 6859.58, 7462.18, 9344.47, 11489.69, 14036.5, 17177.27, 20818.42, 21128.99, 21128.99, 21128.99, 21128.99, 21128.99, 21128.99, 21128.99, 21128.99, 21128.99, 25051.32, 29505.05, 33748.16, 36863.26]

    assert np.allclose(EigenValueT, EigenValueR, atol=2), "Buckling Eigen Mode is wrong."




def Test_BucklingEigenMode(setup_model):
    SecondOrderResponseT = setup_model
    "Not checked implement in future with text document to eigen modes"
    EigenMode1T = SecondOrderResponseT.BucklingEigenLoad(EigenModeNo = True)[2][:,0]
    EigenMode2T = SecondOrderResponseT.BucklingEigenLoad(EigenModeNo = True)[2][:,1]
    EigenMode3T = SecondOrderResponseT.BucklingEigenLoad(EigenModeNo = True)[2][:,2]
    EigenMode4T = SecondOrderResponseT.BucklingEigenLoad(EigenModeNo = True)[2][:,3]
    EigenMode5T = SecondOrderResponseT.BucklingEigenLoad(EigenModeNo = True)[2][:,4]
    EigenMode6T = SecondOrderResponseT.BucklingEigenLoad(EigenModeNo = True)[2][:,5]
    EigenValueR = [0.08, 0.23, 0.45, 0.74, 1.11, 1.56, 2.1, 2.71, 3.38, 4.5, 5.47, 6.65, 8.06, 9.73, 11.67, 13.86, 16.03, 16.09, 17.9, 32.45, 49.67, 68.08, 88.0, 109.5, 131.96, 153.38, 169.73]

    assert np.allclose(EigenValueR, EigenValueR, atol=2), "Buckling Eigen Mode is wrong."


