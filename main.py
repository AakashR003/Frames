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
from sksparse.cholmod import cholesky
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
                UnConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Hinged Joint Support" :
                UnConstrainedDoFList.append(node.dof_tita)
                
            if node.support_condition=="Roller in X-plane-Hinge" :
                UnConstrainedDoFList.append(node.dof_tita)
                UnConstrainedDoFList.append(node.dof_x)
                
            if(node.support_condition=="Rigid Joint"):
                UnConstrainedDoFList.append(node.dof_x)
                UnConstrainedDoFList.append(node.dof_y)
                UnConstrainedDoFList.append(node.dof_tita)
                
            else:
                pass
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
        plt.figure(figsize=(10, 6))
        plt.title("Structural Model")

        # Normalize sensitivities if provided
        if sensitivities is not None:
            min_sensitivity = min(sensitivities)
            max_sensitivity = max(sensitivities)
            # Avoid division by zero if all sensitivities are the same
            if max_sensitivity == min_sensitivity:
                normalized_sensitivities = [0.5 for _ in sensitivities]
            else:
                normalized_sensitivities = [(s - min_sensitivity) / (max_sensitivity - min_sensitivity) for s in sensitivities]

        # Plot members
        for i, member in enumerate(self.Members):
            start_node = member.Start_Node
            end_node = member.End_Node
            if sensitivities is not None:
                # Map normalized sensitivity to a color gradient (red to yellow)
                color = plt.cm.OrRd(normalized_sensitivities[i])  # OrRd colormap: red to yellow
                plt.plot([start_node.xcoordinate, end_node.xcoordinate], 
                        [start_node.ycoordinate, end_node.ycoordinate], color=color, linewidth=2)
            else:
                plt.plot([start_node.xcoordinate, end_node.xcoordinate], 
                        [start_node.ycoordinate, end_node.ycoordinate], 'b-')

        # Plot nodes and support conditions
        for i, node in enumerate(self.Points):
            # Plot nodes
            plt.plot(node.xcoordinate, node.ycoordinate, 'ro')
            
            # Add node numbers
            plt.text(node.xcoordinate, node.ycoordinate + 0.2, f"{i+1}", fontsize=12, ha='center', va='bottom', color='black')

            # Plot support conditions
            if node.support_condition == 'Fixed Support':
                plt.plot(node.xcoordinate, node.ycoordinate, 'ks', markersize=10, label="Fixed Support" if i == 0 else "")
            elif node.support_condition == 'Hinged Support':
                plt.plot(node.xcoordinate, node.ycoordinate, 'g^', markersize=10, label="Hinged Support" if i == 0 else "")
            elif node.support_condition == 'Roller in X-plane':
                plt.plot(node.xcoordinate, node.ycoordinate, 'bv', markersize=10, label="Roller in X-plane" if i == 0 else "")
            elif node.support_condition == 'Roller in Y-plane':
                plt.plot(node.xcoordinate, node.ycoordinate, 'r>', markersize=10, label="Roller in Y-plane" if i == 0 else "")
            elif node.support_condition == 'Hinge Joint':
                plt.plot(node.xcoordinate, node.ycoordinate, 'go', markerfacecolor='none', markersize=10, label="Hinged Support" if i == 0 else "")

        # Find the maximum load magnitude across all loads (PL and UDL)
        max_load_magnitude = max(max(abs(load.Magnitude) for load in self.Loads), 1)  # Ensure at least 1 to avoid division by zero

        # Plot loads
        for load in self.Loads:
            if load.type == "UDL":
                self._plot_udl(load, max_load_magnitude)
            elif load.type == "PL":
                self._plot_point_load(load, max_load_magnitude)
            
        # Add a colorbar if sensitivities are provided
        if sensitivities is not None:
            sm = plt.cm.ScalarMappable(cmap=plt.cm.OrRd, norm=plt.Normalize(vmin=min_sensitivity, vmax=max_sensitivity))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=plt.gca(), label="Sensitivity")  # Explicitly associate with the current Axes
        
        # Add labels and legend
        plt.xlabel("X-coordinate")
        plt.ylabel("Y-coordinate")
        plt.legend(loc="upper right")
        plt.grid(True)
        plt.axis('equal')  # Ensure equal scaling for x and y axes
        plt.show()

    def _plot_point_load(self, load, max_load_magnitude):
        """
        Plots a Point Load on the assigned member.
        :param load: NeumanBC object representing the Point Load.
        :param max_load_magnitude: Maximum load magnitude for scaling.
        """
        # Extract the member number from the AssignedTo attribute (e.g., "Member 1" -> 1)
        member_number = int(load.AssignedTo.split()[1]) - 1  # Convert to zero-based index
        if member_number < 0 or member_number >= len(self.Members):
            raise ValueError(f"Invalid member number {member_number + 1} for load: {load}")

        # Get the assigned member
        member = self.Members[member_number]
        start_node = member.Start_Node
        end_node = member.End_Node

        # Calculate the position of the Point Load
        dx = end_node.xcoordinate - start_node.xcoordinate
        dy = end_node.ycoordinate - start_node.ycoordinate
        length = np.sqrt(dx**2 + dy**2)
        x = start_node.xcoordinate + (load.Distance1 / length) * dx
        y = start_node.ycoordinate + (load.Distance1 / length) * dy

        # Scale the arrow length based on the load magnitude
        arrow_length = 0.5 * (abs(load.Magnitude) / max_load_magnitude)  # Scale the arrow length

        # Plot the Point Load as an arrow on top of the beam
        arrow_dy = -arrow_length if load.Magnitude > 0 else arrow_length  # Arrow direction based on load sign
        plt.arrow(x, y, 0, arrow_dy, head_width=0.2, head_length=0.2, fc='r', ec='r')

    def _plot_udl(self, load, max_load_magnitude):
        """
        Plots a Uniformly Distributed Load (UDL) on the assigned member.
        :param load: NeumanBC object representing the UDL.
        :param max_load_magnitude: Maximum load magnitude for scaling.
        """
        # Extract the member number from the AssignedTo attribute (e.g., "Member 2" -> 2)
        member_number = int(load.AssignedTo.split()[1]) - 1  # Convert to zero-based index
        if member_number < 0 or member_number >= len(self.Members):
            raise ValueError(f"Invalid member number {member_number + 1} for load: {load}")

        # Get the assigned member
        member = self.Members[member_number]
        start_node = member.Start_Node
        end_node = member.End_Node

        # Calculate the direction of the member
        dx = end_node.xcoordinate - start_node.xcoordinate
        dy = end_node.ycoordinate - start_node.ycoordinate
        length = np.sqrt(dx**2 + dy**2)
        angle = np.arctan2(dy, dx)

        # Calculate the start and end points of the UDL
        x1 = start_node.xcoordinate + (load.Distance1 / length) * dx
        y1 = start_node.ycoordinate + (load.Distance1 / length) * dy
        x2 = start_node.xcoordinate + (load.Distance2 / length) * dx
        y2 = start_node.ycoordinate + (load.Distance2 / length) * dy

        # Scale the arrow length based on the load magnitude
        arrow_length = 0.2 * (abs(load.Magnitude) / max_load_magnitude)  # Scale the arrow length

        # Plot the UDL as a series of arrows on top of the beam
        num_arrows = 15  # Number of arrows to represent the UDL
        for i in range(num_arrows):
            xi = x1 + (x2 - x1) * (i / num_arrows)
            yi = y1 + (y2 - y1) * (i / num_arrows)
            # Adjust the arrow position to be on top of the beam
            arrow_dy = -arrow_length if load.Magnitude > 0 else arrow_length  # Arrow direction based on load sign
            plt.arrow(xi, yi, 0, arrow_dy, head_width=0.1, head_length=0.1, fc='g', ec='g')






