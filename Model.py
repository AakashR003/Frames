# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 03:54:23 2024

@author: aakas
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
#import importlib
import math
import scipy.sparse as sp
from scipy.sparse.linalg import eigs, cg
from scipy.sparse import csc_matrix


try:
    from .Computer import Computer
    from .Functions import max_nested
    from .StructuralElements import Node, Member
    from .Loads import NeumanBC
except:
    from Computer import Computer
    from Functions import max_nested
    from StructuralElements import Node, Member
    from Loads import NeumanBC

#import time
#import FiniteElementDivisor


class Model():
    
    def __init__(self,**kwargs):
        
        self.Points = kwargs.get("Points", None)
        self.Members = kwargs.get("Members", None)
        self.Loads = kwargs.get("Loads", None)
        self.NoMembers = len(self.Members)
    
    def UnConstrainedDoF(self):
        """
        Returns a list of unconstrained degrees of freedom (DoF) for the model. 
        The unconstrained DoF are determined based on the support conditions of the nodes.
        The function iterates through each node in the model and appends the corresponding DoF to the UnConstrainedDoFList.
        The function asumes an additional DOF for each node for all dof for excess dof consideration which will be helpful in 
        hinge joint and glided support creation above 100000"""

        UnConstrainedDoFList=[]
        ConstrainedDoFList=[]
        for node in self.Points:
            if node.support_condition=="Hinged Support" :
                UnConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Fixed Support" :
                pass
            
            if node.support_condition=="Roller in X-plane" :
                UnConstrainedDoFList.append(node.dof_x)
                UnConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Roller in Y-plane" :
                UnConstrainedDoFList.append(node.dof_y)
                UnConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Glided Support" :
                ConstrainedDoFList.append(node.dof_x)
                ConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Hinge Joint" :
                UnConstrainedDoFList.append(node.dof_x)
                UnConstrainedDoFList.append(node.dof_y)
                UnConstrainedDoFList = UnConstrainedDoFList + node.additional_dof_tita
                
            if node.support_condition=="Hinged Joint Support" :
                UnConstrainedDoFList = UnConstrainedDoFList + node.additional_dof_tita
                
            if node.support_condition=="Roller in X-plane-Hinge" :
                UnConstrainedDoFList = UnConstrainedDoFList + node.additional_dof_tita
                UnConstrainedDoFList.append(node.dof_x)
                
            if(node.support_condition=="Rigid Joint"):
                UnConstrainedDoFList.append(node.dof_x)
                UnConstrainedDoFList.append(node.dof_y)
                UnConstrainedDoFList.append(node.dof_tita)
                
            else:
                pass
        #Additional check to remove the duplicate dof from the list
        UnConstrainedDoFList= list(set(UnConstrainedDoFList))
        return UnConstrainedDoFList
        
    def ConstrainedDoF(self):
        ConstrainedDoFList=[]
        for node in self.Points:
            if node.support_condition=="Hinged Support" :
                ConstrainedDoFList.append(node.dof_x)
                ConstrainedDoFList.append(node.dof_y)
                
            if node.support_condition=="Fixed Support" :
                ConstrainedDoFList.append(node.dof_x)
                ConstrainedDoFList.append(node.dof_y)
                ConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Roller in X-plane" :
                ConstrainedDoFList.append(node.dof_y)
                
            if node.support_condition=="Roller in Y-plane" :
                ConstrainedDoFList.append(node.dof_x)
                
            if node.support_condition=="Glided Support" :
                ConstrainedDoFList.append(node.dof_x)
                ConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Hinge Joint" :
                pass
            
            if node.support_condition=="Hinged Joint Support" :
                ConstrainedDoFList.append(node.dof_x)
                ConstrainedDoFList.append(node.dof_y)
                
            if node.support_condition=="Roller in X-plane-Hinge" :
                ConstrainedDoFList.append(node.dof_y)
                
            if node.support_condition=="Rigid Joint" :
                pass
            else:
                pass
        return ConstrainedDoFList
    
    def TotalDoF(self):
        return self.UnConstrainedDoF() + self.ConstrainedDoF()
    
    def UnConstrainedDoFDict(self):
        return {num: 0 for num in self.UnConstrainedDoF()}
    
    def TotalDoFDict(self):
        return {num: 0 for num in self.TotalDoF()}
    
    def GlobalStiffnessMatrix(self):
       
        C1 = Computer.StiffnessMatrixAssembler(self.TotalDoF(), self.Members, "First_Order_Global_Stiffness_Matrix_1")
        return C1
    
    def GlobalStiffnessMatrixCondensed(self):
        
        C1 = Computer.StiffnessMatrixAssembler(self.UnConstrainedDoF(), self.Members, "First_Order_Global_Stiffness_Matrix_1")
        return C1
    
    def GlobalStiffnessMatrixCondensedA21(self):
        C1=[]
        for Mc in self.ConstrainedDoF():
            R1=[]
            for Mr in self.UnConstrainedDoF():
                y=0
                for mn in range(0,self.NoMembers):
                    for mr in range(0,6):
                        if(self.Members[mn].DoFNumber()[mr]==Mr):
                            for mc in range(0,6):
                                if(self.Members[mn].DoFNumber()[mc]==Mc):
                                    x=self.Members[mn].First_Order_Global_Stiffness_Matrix_1()[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return C1
    
    def ForceVector(self):
        self.ForceVectorDict=self.TotalDoFDict()
        for var1 in self.Loads:
            if var1.type == "NL":
                self.ForceVectorDict[var1.NodalLoad()['Fx'][1]] = self.ForceVectorDict[var1.NodalLoad()['Fx'][1]] + var1.NodalLoad()['Fx'][0]
                self.ForceVectorDict[var1.NodalLoad()['Fy'][1]] = self.ForceVectorDict[var1.NodalLoad()['Fy'][1]] + var1.NodalLoad()['Fy'][0]
                self.ForceVectorDict[var1.NodalLoad()['Moment'][1]] = self.ForceVectorDict[var1.NodalLoad()['Moment'][1]] + var1.NodalLoad()['Moment'][0]
            
            else:
                self.ForceVectorDict[var1.EquivalentLoad()['Va'][1]] = self.ForceVectorDict[var1.EquivalentLoad()['Va'][1]] + var1.EquivalentLoad()['Va'][0]
                self.ForceVectorDict[var1.EquivalentLoad()['Vb'][1]] = self.ForceVectorDict[var1.EquivalentLoad()['Vb'][1]] + var1.EquivalentLoad()['Vb'][0]
                self.ForceVectorDict[var1.EquivalentLoad()['Ha'][1]] = self.ForceVectorDict[var1.EquivalentLoad()['Ha'][1]] + var1.EquivalentLoad()['Ha'][0]
                self.ForceVectorDict[var1.EquivalentLoad()['Hb'][1]] = self.ForceVectorDict[var1.EquivalentLoad()['Hb'][1]] + var1.EquivalentLoad()['Hb'][0]
                self.ForceVectorDict[var1.EquivalentLoad()['Ma'][1]] = self.ForceVectorDict[var1.EquivalentLoad()['Ma'][1]] + var1.EquivalentLoad()['Ma'][0]
                self.ForceVectorDict[var1.EquivalentLoad()['Mb'][1]] = self.ForceVectorDict[var1.EquivalentLoad()['Mb'][1]] + var1.EquivalentLoad()['Mb'][0]
        ForceVector = []
        for var2 in self.UnConstrainedDoF():
            ForceVector.append(self.ForceVectorDict[var2])
        return ForceVector
    
    def PlotGlobalModel(self, sensitivities=None):
        """
        Plots the structural model using matplotlib.
        """

        fig, ax = plt.subplots(figsize=(10, 6))
        plt.title("Structural Model")
        
        computer_instance = Computer()
        computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = True)
        
        # Find the maximum load magnitude across all loads (PL and UDL)
        max_load_magnitude = max(max(abs(load.Magnitude) for load in self.Loads), 1)  # Ensure at least 1 to avoid division by zero

        # Plot loads
        for load in self.Loads:
            if load.type == "UDL":
                self._plot_udl(load, max_load_magnitude)
            elif load.type == "PL":
                self._plot_point_load(load, max_load_magnitude)
            
        # Add labels and legend
        plt.xlabel("X-coordinate")
        plt.ylabel("Y-coordinate")
        plt.legend(loc="upper right")
        plt.grid(True)
        plt.axis('equal')  # Ensure equal scaling for x and y axes
        plt.show()

    def _plot_point_load(self, load, max_load_magnitude):
        """
        Plots a Point Load on the assigned member in its local coordinate system.
        Positive loads act upward in the local y’ direction.
        """
        member_number = int(load.AssignedTo.split()[1]) - 1
        if member_number < 0 or member_number >= len(self.Members):
            raise ValueError(f"Invalid member number {member_number + 1} for load: {load}")

        member = self.Members[member_number]
        start_node = member.Start_Node
        end_node = member.End_Node

        # Calculate member direction and orientation
        dx = end_node.xcoordinate - start_node.xcoordinate
        dy = end_node.ycoordinate - start_node.ycoordinate
        length = np.sqrt(dx**2 + dy**2)
        theta = np.arctan2(dy, dx)  # Angle of member relative to global x-axis

        # Position of the load along the member
        x = start_node.xcoordinate + (load.Distance1 / length) * dx
        y = start_node.ycoordinate + (load.Distance1 / length) * dy

        # Arrow direction in LOCAL y’ coordinates
        arrow_length = 0.5 * (abs(load.Magnitude) / max_load_magnitude)
        direction_sign = 1 if load.Magnitude > 0 else -1
        dx_arrow = -np.sin(theta) * direction_sign * arrow_length  # Local y’ x-component
        dy_arrow = np.cos(theta) * direction_sign * arrow_length  # Local y’ y-component

        # Plot arrow in local direction
        plt.arrow(x, y, dx_arrow, dy_arrow, head_width=0.2, head_length=0.2, fc='r', ec='r')

    def _plot_udl(self, load, max_load_magnitude):
        """
        Plots a UDL on the assigned member in its local coordinate system.
        Positive loads act upward in the local y’ direction.
        """
        member_number = int(load.AssignedTo.split()[1]) - 1
        if member_number < 0 or member_number >= len(self.Members):
            raise ValueError(f"Invalid member number {member_number + 1} for load: {load}")

        member = self.Members[member_number]
        start_node = member.Start_Node
        end_node = member.End_Node

        dx = end_node.xcoordinate - start_node.xcoordinate
        dy = end_node.ycoordinate - start_node.ycoordinate
        length = np.sqrt(dx**2 + dy**2)
        theta = np.arctan2(dy, dx)  # Angle of member relative to global x-axis

        # Start/end points of UDL segment
        x1 = start_node.xcoordinate + (load.Distance1 / length) * dx
        y1 = start_node.ycoordinate + (load.Distance1 / length) * dy
        x2 = start_node.xcoordinate + (load.Distance2 / length) * dx
        y2 = start_node.ycoordinate + (load.Distance2 / length) * dy

        # Arrow direction in LOCAL y’ coordinates
        arrow_length = 0.2 * (abs(load.Magnitude) / max_load_magnitude)
        direction_sign = 1 if load.Magnitude > 0 else -1
        dx_arrow = -np.sin(theta) * direction_sign * arrow_length  # Local y’ x-component
        dy_arrow = np.cos(theta) * direction_sign * arrow_length  # Local y’ y-component

        # Plot UDL as arrows along the member
        num_arrows = 15
        for i in range(num_arrows):
            xi = x1 + (x2 - x1) * (i / num_arrows)
            yi = y1 + (y2 - y1) * (i / num_arrows)
            plt.arrow(xi, yi, dx_arrow, dy_arrow, head_width=0.1, head_length=0.1, fc='g', ec='g')