# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 04:23:42 2025

@author: aakas
"""

import pytest
import numpy as np
from main import Node, Member, NeumanBC, Model, GlobalResponse

@pytest.fixture
def setup_model():
    PointsT = [
        Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Fixed Support"),
        Node(Node_Number=2, xcoordinate=10, ycoordinate=0, Support_Condition="Rigid Joint"),
        Node(Node_Number=3, xcoordinate=20, ycoordinate=0, Support_Condition="Fixed Support")
    ]
    MembersT = [
        Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
        Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=1, Youngs_Modulus=1, Moment_of_Inertia=1),
    ]
    Loads = [
        NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=10, AssignedTo="Member 1", Members = MembersT),
        NeumanBC(type="UDL", Magnitude=5, Distance1=0, Distance2=10, AssignedTo="Member 2", Members = MembersT)
    ]
    Model1 = Model(Points=PointsT, Members=MembersT, Loads=Loads)
    Res = GlobalResponse(Points=PointsT, Members=MembersT, Loads=Loads)
    return Model1, Res, MembersT

def test_FixedBeamUDL(setup_model):
    Model1, Res, MembersT = setup_model

    # Member Stiffness Matrix of Member 1
    Element1MSM = MembersT[0].First_Order_Global_Stiffness_Matrix_1()
    Element2MSM = MembersT[1].First_Order_Global_Stiffness_Matrix_1()

    Element1MSM_R = [
        [0.1, 0, 0, -0.1, 0, 0],
        [0, 0.012, 0.06, 0, -0.012, 0.06],
        [0, 0.06, 0.4, 0, -0.06, 0.2],
        [-0.1, 0, 0, 0.1, 0, 0],
        [0, -0.012, -0.06, 0, 0.012, -0.06],
        [0, 0.06, 0.2, 0, -0.06, 0.4]
    ]

    Element2MSM_R = [
        [0.1, 0, 0, -0.1, 0, 0],
        [0, 0.012, 0.06, 0, -0.012, 0.06],
        [0, 0.06, 0.4, 0, -0.06, 0.2],
        [-0.1, 0, 0, 0.1, 0, 0],
        [0, -0.012, -0.06, 0, 0.012, -0.06],
        [0, 0.06, 0.2, 0, -0.06, 0.4]
    ]

    # Assert matrix dimensions
    assert len(Element1MSM) == 6, "Stiffness Matrix does not have 6 rows."
    assert all(len(row) == 6 for row in Element1MSM), "Stiffness Matrix does not have 6 columns."

    # Assert member stiffness matrix
    assert np.allclose(Element1MSM, Element1MSM_R, atol=1e-3), "Stiffness matrix of Member 1 is wrong."
    assert np.allclose(Element2MSM, Element2MSM_R, atol=1e-3), "Stiffness matrix of Member 2 is wrong."

    GSMS = Model1.GlobalStiffnessMatrixCondensed()
    GSMS_R = [
        [0.2, 0, 0],
        [0, 0.024, 0],
        [0, 0, 0.8]
    ]

    # Assert Structures Global Stiffness matrix - condensed a11
    assert np.allclose(GSMS, GSMS_R, atol=1e-3), "Structure Global Stiffness matrix A11 is wrong."

    GSMSa21 = Model1.GlobalStiffnessMatrixCondensedA21()
    GSMSa21_R = [
        [-0.1, 0, 0],
        [0, -0.012, 0.06],
        [0, -0.06, 0.2],
        [-0.1, 0, 0],
        [0, -0.012, -0.06],
        [0, 0.06, 0.2]
    ]

    # Assert Structures Global Stiffness matrix - condensed a21
    assert np.allclose(GSMSa21, GSMSa21_R, atol=1e-3), "Structure Global Stiffness matrix A21 is wrong"

    FV = Model1.ForceVector()
    FV_R = [0, -50, 0]
    assert np.allclose(FV, FV_R, atol=1e-3), "Force Vector is wrong"

    DV = Res.DisplacementVector()
    DV_R = [0, -2083.333333, -1.33227E-13]
    assert np.allclose(DV, DV_R, atol=1e-3), "Displacement Vector is wrong"
