""" NPTEL Week - 7 example"""

import pytest
import numpy as np

from config import config
from main import Model
from StructuralElements import Node, Member
from Loads import NeumanBC
from FirstOrderResponse import FirstOrderGlobalResponse, FirstOrderMemberResponse, FirstOrderNodalResponse

"PointsT - T stands for test"

config.set_FEDivision(1000)
PointsT = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Fixed Support"),
              Node(Node_Number=2,xcoordinate=240,ycoordinate=180,Support_Condition="Rigid Joint"),
              Node(Node_Number=3,xcoordinate=480,ycoordinate=180,Support_Condition="Fixed Support")] 

MembersT = [Member(Beam_Number=1,Start_Node=PointsT[0],End_Node=PointsT[1],Area=12,Youngs_Modulus=29000,Moment_of_Inertia=600),
            Member(Beam_Number=2,Start_Node=PointsT[1],End_Node=PointsT[2],Area=12,Youngs_Modulus=29000,Moment_of_Inertia=600),]

LoadsT = [NeumanBC(type="UDL",Magnitude=-0.25,Distance1=0,Distance2=240,AssignedTo="Member 2", Members = MembersT), # converted 3ksi to 0.25 k/inch
        ] 


ModelT = Model(Points = PointsT, Members = MembersT, Loads = LoadsT)
GlobalResponseT = FirstOrderGlobalResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)
MemberResT = FirstOrderMemberResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)


def test_MSML(): #Member Stiffness Matrix Local Coordinate

    Element1MSM = MembersT[0].First_Order_Local_Stiffness_Matrix_1()

    Element2MSM = MembersT[1].First_Order_Local_Stiffness_Matrix_1()
    # Expected matrix

    Element1MSMR = [
    [1160, 0, 0, -1160, 0, 0],
    [0, 7.733333333, 1160, 0, -7.733333333, 1160],
    [0, 1160, 232000, 0, -1160, 116000],
    [-1160, 0, 0, 1160, 0, 0],
    [0, -7.733333333, -1160, 0, 7.733333333, -1160],
    [0, 1160, 116000, 0, -1160, 232000]
    ]

    Element2MSMR = [
    [1450, 0, 0, -1450, 0, 0],
    [0, 15.10416667, 1812.5, 0, -15.10416667, 1812.5],
    [0, 1812.5, 290000, 0, -1812.5, 145000],
    [-1450, 0, 0, 1450, 0, 0],
    [0, -15.10416667, -1812.5, 0, 15.10416667, -1812.5],
    [0, 1812.5, 145000, 0, -1812.5, 290000]
    ]
    
    assert len(Element1MSM) == 6, "Stiffness Matrix does not have 6 rows."
    assert all(len(row) == 6 for row in Element1MSM), "Stiffness Matrix does not have 6 columns."

    # Assert matrix content
    assert np.allclose(Element1MSM, Element1MSMR, atol=0), "Stiffness matrix of Member 1 is wrong."
    assert np.allclose(Element2MSM, Element2MSMR, atol=0), "Stiffness matrix of Member 2 is wrong."



def test_MSMG(): #Member Stiffness Matrix Global Coordinate

    Element1MSM = MembersT[0].First_Order_Global_Stiffness_Matrix_1()

    Element2MSM = MembersT[1].First_Order_Global_Stiffness_Matrix_1()
    # Expected matrix

    Element1MSM_T = [
    [745.184, 553.088, -696, -745.184, -553.088, -696],
    [553.088, 422.5493333, 928, -553.088, -422.5493333, 928],
    [-696, 928, 232000, 696, -928, 116000],
    [-745.184, -553.088, 696, 745.184, 553.088, 696],
    [-553.088, -422.5493333, -928, 553.088, 422.5493333, -928],
    [-696, 928, 116000, 696, -928, 232000]
    ]

    Element2MSM_T = [
    [1450, 0, 0, -1450, 0, 0],
    [0, 15.10416667, 1812.5, 0, -15.10416667, 1812.5],
    [0, 1812.5, 290000, 0, -1812.5, 145000],
    [-1450, 0, 0, 1450, 0, 0],
    [0, -15.10416667, -1812.5, 0, 15.10416667, -1812.5],
    [0, 1812.5, 145000, 0, -1812.5, 290000]
    ]
    
    assert len(Element1MSM) == 6, "Stiffness Matrix does not have 6 rows."
    assert all(len(row) == 6 for row in Element1MSM), "Stiffness Matrix does not have 6 columns."

    # Assert matrix content
    assert np.allclose(Element1MSM, Element1MSM_T, atol=0), "Stiffness matrix of Member 1 is wrong."
    assert np.allclose(Element2MSM, Element2MSM_T, atol=0), "Stiffness matrix of Member 2 is wrong."


def test_GSM(): #Global Stiffness Matrix A11 
    GSM = ModelT.GlobalStiffnessMatrixCondensed()

    GSMT = [
            [2195.18, 553.09, 696],
            [553.09, 437.65, 884.5],
            [696, 884.5, 522000]
        ]
    assert np.allclose(GSM, GSMT, atol=0.01), "Stiffness matrix A11  is wrong."

def test_GSMA21():  # Condensed Global Stiffness Matrix A21 testing

    GSMA21 = ModelT.GlobalStiffnessMatrixCondensedA21()

    GSMA21T = [
                [-745.18, -553.09, -696],
                [-553.09, -422.55, 928],
                [696, -928, 116000],
                [-1450, 0, 0],
                [0, -15.1, -1812.5],
                [0, 1812.5, 145000]
            ]
    assert np.allclose(GSMA21, GSMA21T, atol=0.01), "Stiffness matrix A21 is wrong." 


def test_Displacement(): # test the displacement response

    DisplacementList = GlobalResponseT.DisplacementVector()


    DisplacementListR = [
                             0.0247,
                            -0.0954,
                            -0.00217,
                            ]
    assert np.allclose(DisplacementList, DisplacementListR, atol=0.0001), "Displacement is wrong." 

def test_supportforce():

    ForceList = GlobalResponseT.SupportForcesVector()

    ForceListR = [
        35.86,
        24.63,
        -145.99,
        -35.85,
        5.37,
        -487.6,
    ]
    assert np.allclose(ForceList, ForceListR, atol=0.01), "Force Vector 1 is wrong." 


def test_MemberForce():

    Member1Force = MemberResT.MemberForceLocal(2)
    print(Member1Force)
