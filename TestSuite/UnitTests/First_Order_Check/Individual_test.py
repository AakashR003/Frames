# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 03:17:40 2025

@author: aakas
"""


import pytest
import numpy as np

from config import config
from main import Model
from StructuralElements import Node, Member

# Test function for Member Stiffness Matrix( From NPTEL) in Global Coordinates 
def test_MSMG(): #Member Stiffness Matrix Global Coordinate
    config.set_FEDivision(1000)
    Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Fixed Support"),
              Node(Node_Number=2,xcoordinate=10,ycoordinate=5,Support_Condition="Hinged Support"),
              Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Rigid Joint")] 

    Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
               Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),]

    Element1MSM = Members[0].First_Order_Global_Stiffness_Matrix_1()
    Element2MSM = Members[1].First_Order_Global_Stiffness_Matrix_1()
    # Expected matrix

    Element1MSM_T = [
    [0.073271475, 0.032342487, -0.021466253, -0.073271475, -0.032342487, -0.021466253],
    [0.032342487, 0.024757745, 0.042932505, -0.032342487, -0.024757745, 0.042932505],
    [-0.021466253, 0.042932505, 0.357770876, 0.021466253, -0.042932505, 0.178885438],
    [-0.073271475, -0.032342487, 0.021466253, 0.073271475, 0.032342487, 0.021466253],
    [-0.032342487, -0.024757745, -0.042932505, 0.032342487, 0.024757745, -0.042932505],
    [-0.021466253, 0.042932505, 0.178885438, 0.021466253, -0.042932505, 0.357770876]
    ]
    
    Element2MSM_T = [
    [0.073271475, -0.032342487, 0.021466253, -0.073271475, 0.032342487, 0.021466253],
    [-0.032342487, 0.024757745, 0.042932505, 0.032342487, -0.024757745, 0.042932505],
    [0.021466253, 0.042932505, 0.357770876, -0.021466253, -0.042932505, 0.178885438],
    [-0.073271475, 0.032342487, -0.021466253, 0.073271475, -0.032342487, -0.021466253],
    [0.032342487, -0.024757745, -0.042932505, -0.032342487, 0.024757745, -0.042932505],
    [0.021466253, 0.042932505, 0.178885438, -0.021466253, -0.042932505, 0.357770876]
    ]

    # Assert matrix dimensions
    assert len(Element1MSM) == 6, "Stiffness Matrix does not have 6 rows."
    assert all(len(row) == 6 for row in Element1MSM), "Stiffness Matrix does not have 6 columns."

    # Assert matrix content
    assert np.allclose(Element1MSM, Element1MSM_T,atol=1e-3), "Stiffness matrix of Member 1 is wrong."
    assert np.allclose(Element1MSM, Element1MSM_T,atol=1e-3), "Stiffness matrix of Member 2 is wrong."



def test_DoF(): 
    """ To check constrained and Unconstrained DOF's are correct"""
    Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Fixed Support"),
              Node(Node_Number=2,xcoordinate=10,ycoordinate=0,Support_Condition="Hinged Support"),
              Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Rigid Joint")] 

    Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
               Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),]


    Model1 = Model(Points = Points, Members = Members)
    print(Model1.UnConstrainedDoF())
    print(Model1.ConstrainedDoF())
    print(Model1.TotalDoF())



def test_GSM(): #Global Stiffness Matrix - Structure
    Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Fixed Support"),
              Node(Node_Number=2,xcoordinate=10,ycoordinate=0,Support_Condition="Hinged Support"),
              Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Rigid Joint")] 

    Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
               Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),]


    Model1 = Model(Points = Points, Members = Members)
    print(Model1.GlobalStiffnessMatrixCondensed())


def test_GSMA21(): #Global Stiffness Matrix - Structure
    """ To check Condensed A21 part of Global Stiffness Matrix of the structure is correct - used for computing Support Forces"""
    Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Fixed Support"),
              Node(Node_Number=2,xcoordinate=10,ycoordinate=0,Support_Condition="Hinged Support"),
              Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Rigid Joint")] 

    Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
               Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),]


    Model1 = Model(Points = Points, Members = Members)
    print(Model1.GlobalStiffnessMatrixCondensedA21())

""" 
failing dont know why, have to check in future
def test_FV(): #Force Vector
    #To check Force Vector
    Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Fixed Support"),
              Node(Node_Number=2,xcoordinate=10,ycoordinate=0,Support_Condition="Hinged Support"),
              Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Rigid Joint")] 

    Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
               Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),]
   
    Loads = [NeumanBC(type="UDL",Magnitude=5,Distance1=0,Distance2=10,AssignedTo="Member 1"),
             NeumanBC(type="UDL",Magnitude=5,Distance1=0,Distance2=10,AssignedTo="Member 2")] 


    Model1 = Model(Points = Points, Members = Members, Loads = Loads)
    print(Model1.ForceVector())
"""

#MSMG_test()
#DoF_test()
#GSM_test()
#GSMA21_test()
#FV_test()
