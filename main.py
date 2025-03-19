# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 03:54:23 2024

@author: aakas
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
import importlib


class Node():
    
    i=1
    titam = 1
    titay = 1
    def __init__ (self, Node_Number, xcoordinate, ycoordinate, Support_Condition):
        self.node_number = Node_Number
        self.xcoordinate = xcoordinate
        self.ycoordinate = ycoordinate
        self.support_condition = Support_Condition

        if self.support_condition in ["Hinged Support", "Fixed Support", "Rigid Joint", "Roller in X-plane", "Roller in Y-plane"]:
            self.dof_x=(self.node_number)*3-2
            self.dof_y=(self.node_number)*3-1
            self.dof_tita=(self.node_number)*3
        elif self.support_condition == "Hinge Joint" :
            self.dof_x=(self.node_number)*3-2
            self.dof_y=(self.node_number)*3-1
            self.dof_tita=3000+Node.titam
            Node.titam += 1
        elif self.support_condition == "Glided Support" :
            self.dof_x=(self.node_number)*3-2
            self.dof_y=2000+Node.titay
            self.dof_tita=(self.node_number)*3
            Node.titay += 1
        elif self.support_condition == "Hinged Joint Support" or self.support_condition == "Roller in X-plane-Hinge" :
            self.dof_x=(self.node_number)*3-2
            self.dof_y=(self.node_number)*3-1
            self.dof_tita=3000+Node.titam
            Node.titam += 1
        else:
            raise ValueError(f"Unsupported support condition: '{self.support_condition}'")
        
        
    def DoF(self):
        self.check1 = [self.dof_x, self.dof_y, self.dof_tita]    
        return [self.dof_x, self.dof_y, self.dof_tita]
        
    
class Member():
    
    def __init__ (self, Beam_Number, Start_Node, End_Node, Area, Youngs_Modulus, Moment_of_Inertia ):
        
        self.Beam_Number = Beam_Number
        self.Start_Node = Start_Node
        self.End_Node = End_Node
        
        self.area=Area
        self.youngs_modulus = Youngs_Modulus
        self.moment_of_inertia = Moment_of_Inertia
        
    def length(self):
        x=((self.End_Node.xcoordinate-self.Start_Node.xcoordinate)**2 +
           (self.End_Node.ycoordinate-self.Start_Node.ycoordinate)**2)**0.5
        return x
    
    def alpha(self):
        return (self.End_Node.xcoordinate-self.Start_Node.xcoordinate)/self.length()
   
    def beta(self):
        return (self.End_Node.ycoordinate-self.Start_Node.ycoordinate)/self.length()
    
    def DoFNumber(self):
        return(self.Start_Node.DoF() + self.End_Node.DoF())
    
    def Transformation_Matrix(self):
        Transformation_Matrix=[[self.alpha(),self.beta(),0,0,0,0],
                                [-self.beta(),self.alpha(),0,0,0,0],
                                [0,0,1,0,0,0],
                                [0,0,0,self.alpha(),self.beta(),0],
                                [0,0,0,-self.beta(),self.alpha(),0],
                                [0,0,0,0,0,1]]
     
        return Transformation_Matrix
    
    def First_Order_Local_Stiffness_Matrix_2(self, NormalForce = None):
        """ This Stiffness Matrix is from Lecture Notes of Stability of Structure 
        - typically used for 2nd order analysis"""
        length=self.length()
        ma11=self.area*self.youngs_modulus/length
        ma22=12*self.youngs_modulus*self.moment_of_inertia/length**3
        ma23=6*self.youngs_modulus*self.moment_of_inertia/length**2
        ma32=6*self.youngs_modulus*self.moment_of_inertia/length**2
        ma33=4*self.youngs_modulus*self.moment_of_inertia/length
        ma36=2*self.youngs_modulus*self.moment_of_inertia/length

        Stiffness_Matrix=[[ma11,0,0,-ma11,0,0],
                           [0,ma22,-ma23,0,-ma22,-ma23],
                           [0,-ma32,ma33,0,ma32,ma36],
                           [-ma11,0,0,ma11,0,0],
                           [0,-ma22,ma23,0,ma22,ma23],
                           [0,-ma32,ma36,0,ma32,ma33]]
         
        return Stiffness_Matrix
     
    def First_Order_Global_Stiffness_Matrix_2(self, NormalForce = None):    
        return np.transpose(self.Transformation_Matrix()) @ np.array(self.First_Order_Local_Stiffness_Matrix_2()) @ np.array(self.Transformation_Matrix())

    def First_Order_Local_Stiffness_Matrix_1(self, NormalForce = None):   
        """ This Stiffness Matrix is from NPTEL - Matrix Method """
        length=self.length()
        ma11=self.area*self.youngs_modulus/length
        ma22=12*self.youngs_modulus*self.moment_of_inertia/length**3
        ma23=6*self.youngs_modulus*self.moment_of_inertia/length**2
        ma32=6*self.youngs_modulus*self.moment_of_inertia/length**2
        ma33=4*self.youngs_modulus*self.moment_of_inertia/length
        ma36=2*self.youngs_modulus*self.moment_of_inertia/length

        Stiffness_Matrix=[[ma11,0,0,-ma11,0,0],
                           [0,ma22,ma23,0,-ma22,ma23],
                           [0,ma32,ma33,0,-ma32,ma36],
                           [-ma11,0,0,ma11,0,0],
                           [0,-ma22,-ma23,0,ma22,-ma23],
                           [0,ma32,ma36,0,-ma32,ma33]]
         
        return Stiffness_Matrix
     
    def First_Order_Global_Stiffness_Matrix_1(self, NormalForce = None):    

        return np.transpose(self.Transformation_Matrix()) @ np.array(self.First_Order_Local_Stiffness_Matrix_1()) @ np.array(self.Transformation_Matrix())

    def Second_Order_Reduction_Matrix(self, NormalForce):

        length=self.length()
        ReductionMatrix=[[0,0,0,0,0,0],
            [0,6/5/length,-1/10,0,-6/5/length,-1/10],
            [0,-1/10,2/15*length,0,1/10,-1/30*length],
            [0,0,0,0,0,0],
            [0,-6/5/length,1/10,0,6/5/length,1/10],
            [0,-1/10,-1/30*length,0,1/10,2/15*length]]
        
        StiffnessReductionMatrix = np.array(ReductionMatrix) * NormalForce
        
        return StiffnessReductionMatrix
    
    def Seconf_Order_Local_Stiffness_Matrix(self,NormalForce):

        SecondOrderLocalStiffnessMatrix = np.array(self.First_Order_Local_Stiffness_Matrix_2()) + self.Second_Order_Reduction_Matrix(NormalForce)
        return SecondOrderLocalStiffnessMatrix
    
    def Second_Order_Global_Stiffness_Matrix(self, NormalForce):

        return np.transpose(self.Transformation_Matrix()) @ np.array(self.Seconf_Order_Local_Stiffness_Matrix(NormalForce)) @ np.array(self.Transformation_Matrix())

class Stiffness_Matrix():
    
    def __init__ (self,member):
        self.member=member
        
        
    def First_Order_Local_Stiffness_Matrix_2(self):   
        length=self.member.length()
        ma11=self.member.area*self.member.youngs_modulus/length
        ma22=12*self.member.youngs_modulus*self.member.moment_of_inertia/length**3
        ma23=6*self.member.youngs_modulus*self.member.moment_of_inertia/length**2
        ma32=6*self.member.youngs_modulus*self.member.moment_of_inertia/length**2
        ma33=4*self.member.youngs_modulus*self.member.moment_of_inertia/length
        ma36=2*self.member.youngs_modulus*self.member.moment_of_inertia/length

        Stiffness_Matrix=[[ma11,0,0,-ma11,0,0],
                           [0,ma22,-ma23,0,-ma22,-ma23],
                           [0,-ma32,ma33,0,ma32,ma36],
                           [-ma11,0,0,ma11,0,0],
                           [0,-ma22,ma23,0,ma22,ma23],
                           [0,-ma32,ma36,0,ma32,ma33]]
         
        return Stiffness_Matrix
     
    def First_Order_Global_Stiffness_Matrix_2(self):    
        return np.transpose(self.member.Transformation_Matrix()) @ np.array(self.First_Order_Local_Stiffness_Matrix()) @ np.array(self.member.Transformation_Matrix())


class NeumanBC():
    
    def __init__(self,**kwargs):
        self.type = kwargs.get("type", None)
        self.Magnitude = kwargs.get("Magnitude", None)
        self.Distance1 = kwargs.get("Distance1", None)
        self.Distance2 = kwargs.get("Distance2", None)
        self.AssignedTo = kwargs.get("AssignedTo", None)
        self.Members = kwargs.get("Members", None)
        
        self.MemberNo = int(self.AssignedTo[-1])-1
    
    def EquivalentLoad(self):
        
        self.frml=[] # Free moment Distribution(Simply supported) along beam 
        tarea=0
        tyda=0
        mp=0
        if self.type == "PL" :
            va=-self.Magnitude*(self.Members[self.MemberNo].length()-self.Distance1)/(self.Members[self.MemberNo].length())
            vb=-self.Magnitude*self.Distance1/(self.Members[self.MemberNo].length())
            while mp<=self.Members[self.MemberNo].length():
                if mp > self.Distance1 :
                    mppl=(self.Magnitude)*(mp-self.Distance1)
                else:
                    mppl=0
                m=va*mp+mppl
                area=m*(self.Members[self.MemberNo].length()/1000)
                yda=area*mp
                tarea=area+tarea
                tyda=yda+tyda
                self.frml.append(m)
                mp=mp+(self.Members[self.MemberNo].length()/1000)
            if(tarea==0):
                centroid=0
            else:
                centroid=tyda/tarea
            
        elif self.type == "UDL" :
            self.Range = abs(self.Distance2 - self.Distance1)
            va=-self.Magnitude*self.Range*(self.Members[self.MemberNo].length()-self.Distance1-self.Range*0.5)/(self.Members[self.MemberNo].length())
            vb=-self.Magnitude*self.Range*(self.Distance1+self.Range*0.5)/(self.Members[self.MemberNo].length())
            while mp<=self.Members[self.MemberNo].length():
                if(mp>self.Distance1 and mp<=(self.Distance1+self.Range)):
                    mpu=self.Magnitude*0.5*(mp-self.Distance1)**2
                elif(mp>(self.Distance1+self.Range) and mp<self.Members[self.MemberNo].length()):
                    mpu=self.Magnitude*self.Range*(self.Range*0.5+(mp-(self.Distance1+self.Range)))
                else:
                    mpu=0
                m=va*mp+mpu
                area=m*(self.Members[self.MemberNo].length()/1000)
                yda=area*mp
                tarea=area+tarea
                tyda=yda+tyda
                self.frml.append(m)
                mp=mp+(self.Members[self.MemberNo].length()/1000)
            if(tarea==0):
                centroid=0
            else:
                centroid=tyda/tarea
        
        self.mfab=((2*(tarea*(self.Members[self.MemberNo].length()-centroid)*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length())-(tarea*centroid*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length()))/3)
        self.mfba=-(2*(tarea*centroid*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length())-(tarea*(self.Members[self.MemberNo].length()-centroid)*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length()))/3

        
        if(self.Members[self.MemberNo].alpha()>=0):
            self.V_b=(-self.mfab-self.mfba+vb*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
            self.V_a=(self.mfab+self.mfba+va*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
        else:
            self.V_b=-(-self.mfab-self.mfba+vb*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
            self.V_a=-(self.mfab+self.mfba+va*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
        
        
        return {"Ha":(0,self.Members[self.MemberNo].DoFNumber()[0]),
                "Va":(self.V_a,self.Members[self.MemberNo].DoFNumber()[1]),
                "Ma":(self.mfab,self.Members[self.MemberNo].DoFNumber()[2]),
                "Hb":(0,self.Members[self.MemberNo].DoFNumber()[3]),
                "Vb":(self.V_b,self.Members[self.MemberNo].DoFNumber()[4]),
                "Mb":(self.mfba,self.Members[self.MemberNo].DoFNumber()[5]),
                "FreeMoment": self.frml}#[self.V_a,self.mfab,self.V_b,self.mfba]

    
"""  
dont use error is there - because while computing 


    def Va(self):
        return self.V_a
    def Vb(self):
        return self.V_b
    def Ma(self):
        return self.mfab
    def Mb(self):
        return self.mfba
    
"""


class Computer():
    """
    This class shall be used in future for combining common computers on different class into gloabl computer
    """

    def GlobalStiffnessMatrix(TotalDoF,NoMembers,Members,StiffnessMatrixType):

        C1=[]
        for Mc in TotalDoF():
            R1=[]
            for Mr in TotalDoF():
                y=0
                for mn in range(0,NoMembers):
                    for mr in range(0,6):
                        if(Members[mn].DoFNumber()[mr]==Mr):
                            for mc in range(0,6):
                                if(Members[mn].DoFNumber()[mc]==Mc):
                                    x = getattr(Members[mn],StiffnessMatrixType)[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return None
    
    def GLobalStifnessMatrixCondensedA11(UnConstrainedDoF,NoMembers,Members,StiffnessMatrixType):
        C1=[]
        for Mc in UnConstrainedDoF():
            R1=[]
            for Mr in UnConstrainedDoF():
                y=0
                for mn in range(0,NoMembers):
                    for mr in range(0,6):
                        if(Members[mn].DoFNumber()[mr]==Mr):
                            for mc in range(0,6):
                                if(Members[mn].DoFNumber()[mc]==Mc):
                                    x = getattr(Members[mn],StiffnessMatrixType)[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return C1
    
    def GlobalStifnessMatrixA21():
        return None
    
    def DisplacementVector():
        return None
    
    def ModelDisplacementList_To_Dict(Displacement,UnConstrainedDoF,TotalDoF):

        DisplacementDict={}
        for i in range(len(TotalDoF())):
            if(i<(len(UnConstrainedDoF()))):
                DisplacementDict[str(TotalDoF()[i])] = Displacement[i]
            else:
                DisplacementDict[str(TotalDoF()[i])]=0
        return DisplacementDict

    def SupportForceVector():
        return None

    def ModelDisplacement_To_MemberDisplacement(MemberNumber,DisplacementDict,Members):
        MemberNo = int(MemberNumber)
        MemberDisplacement = [DisplacementDict[str(Members[MemberNo-1].DoFNumber()[0])],
                             DisplacementDict[str(Members[MemberNo-1].DoFNumber()[1])],
                             DisplacementDict[str(Members[MemberNo-1].DoFNumber()[2])],
                             DisplacementDict[str(Members[MemberNo-1].DoFNumber()[3])],
                             DisplacementDict[str(Members[MemberNo-1].DoFNumber()[4])],
                             DisplacementDict[str(Members[MemberNo-1].DoFNumber()[5])]]
        return MemberDisplacement
    
    
    def MemberDisplacement_To_ForceLocal(StiffnessMatrixType, MemberNumber, Members, MemberDisplacement, Loads, NormalForce = None):
        
        MemberNo = int(MemberNumber)
        MemberForce = np.dot(
                    np.dot(
                    getattr(Members[MemberNo-1],StiffnessMatrixType)(NormalForce),Members[MemberNo-1].Transformation_Matrix()),
                    MemberDisplacement)
        FixedendForce = [0, 0, 0, 0, 0, 0]
        for a in range(len(Loads)):
            if(int(Loads[a].AssignedTo[-1]) == MemberNo):
                FixedendForcei = list(Loads[a].EquivalentLoad().values())[:-1]
                FixedendForcei= [x[0] for x in FixedendForcei]
                FixedendForce = [x + y for x, y in zip(FixedendForce, FixedendForcei)]
        MemberForce = np.round(MemberForce - FixedendForce,2)

        return MemberForce


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
        C1=[]
        for Mc in self.TotalDoF():
            R1=[]
            for Mr in self.TotalDoF():
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
    
    def GlobalStiffnessMatrixCondensed(self):
        C1=[]
        for Mc in self.UnConstrainedDoF():
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
        member_number = int(load.AssignedTo[-1]) - 1  # Convert to zero-based index
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
        member_number = int(load.AssignedTo[-1]) - 1  # Convert to zero-based index
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


class GlobalResponse(Model):
    
    def DisplacementVector(self):
        self.Displacement=np.dot((np.linalg.inv(np.array(self.GlobalStiffnessMatrixCondensed()))),self.ForceVector())
        
        #DisplacementDict formation
        self.DisplacementDict={}
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.UnConstrainedDoF()))):
                self.DisplacementDict[str(self.TotalDoF()[i])] = self.Displacement[i]
            else:
                self.DisplacementDict[str(self.TotalDoF()[i])] = 0
        return self.Displacement
    
    def DisplacementVectorDict(self):
        self.DisplacementDict={}
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.UnConstrainedDoF()))):
                self.DisplacementDict[str(self.TotalDoF()[i])] = self.Displacement[i]
            else:
                self.DisplacementDict[str(self.TotalDoF()[i])]=0
        return self.DisplacementDict
    
    def SupportForcesVector(self):

        SupportForces = np.dot(np.array(self.GlobalStiffnessMatrixCondensedA21()),self.DisplacementVector())
        
        self.ForceVectorDict={}
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.ConstrainedDoF()))):
                self.ForceVectorDict[str(self.TotalDoF()[i])]=SupportForces[i]
            else:
                self.ForceVectorDict[str(self.TotalDoF()[i])] = 0
        
        #force dict formation
        return SupportForces
   

class NodalResponse(GlobalResponse):
    
    def NodeDisplacement(self,NodeNumber):
        self.NodeNo = int(NodeNumber)
        
        NodeDisplacement = []
        for i in range(3):
            if self.Points[self.NodeNo-1].DoF()[i] in self.UnConstrainedDoF():
                NodeDisplacement.append(self.DisplacementVectorDict()[str(self.Points[self.NodeNo-1].DoF()[i])])
            else:
                NodeDisplacement.append(0)
        
        return NodeDisplacement
    
    def NodeForce(self,NodeNumber):
        self.SupportForcesVector()
        self.NodeNo = int(NodeNumber)
        
        NodeForce =[]
        for i in range(3):
            if self.Points[self.NodeNo-1].DoF()[i] in self.ConstrainedDoF():
                NodeForce.append(self.ForceVectorDict[str(self.Points[self.NodeNo-1].DoF()[i])])
            else:
                NodeForce.append(0)
        
        return NodeForce


class MemberResponse(GlobalResponse):
    
    """
    init is fomred to have MemberNumber called single time all thorought the class, but class 
    is childrean of another, which variable to call will become a issue, hence not activated now
    
    def __init__(self, MemberNumber):
        
        self.MemberNo = int(MemberNumber)
        if self.MemberNo == "" or float(self.MemberNo) > self.NoMembers:
            self.MemberNo = 1
        else:
            self.MemberNo = int(self.MemberNo)
    """
    
    def MemberDisplacement(self, MemberNumber):
        self.MemberNo = int(MemberNumber)
        self.DisplacementVector()
        MemberDisplacement = [self.DisplacementDict[str(self.Members[self.MemberNo-1].DoFNumber()[0])],
                             self.DisplacementDict[str(self.Members[self.MemberNo-1].DoFNumber()[1])],
                             self.DisplacementDict[str(self.Members[self.MemberNo-1].DoFNumber()[2])],
                             self.DisplacementDict[str(self.Members[self.MemberNo-1].DoFNumber()[3])],
                             self.DisplacementDict[str(self.Members[self.MemberNo-1].DoFNumber()[4])],
                             self.DisplacementDict[str(self.Members[self.MemberNo-1].DoFNumber()[5])]]
        return MemberDisplacement
    
    
    def MemberForceLocal(self, MemberNumber):
        self.MemberNo = int(MemberNumber)
        self.ForceVector()
        MemberFixedEndForce = [self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[0]],
                             self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[1]],
                             self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[2]],
                             self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[3]],
                             self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[4]],
                             self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[5]]]
    
        MemberForce = np.dot(np.dot(self.Members[self.MemberNo-1].First_Order_Local_Stiffness_Matrix_1(),self.Members[self.MemberNo-1].Transformation_Matrix()),self.MemberDisplacement(MemberNumber))

        FixedendForce = [0, 0, 0, 0, 0, 0]
        for a in range(len(self.Loads)):
            if(int(self.Loads[a].AssignedTo[-1]) == self.MemberNo):
                FixedendForcei = list(self.Loads[a].EquivalentLoad().values())[:-1]
                FixedendForcei= [x[0] for x in FixedendForcei]
                FixedendForce = [x + y for x, y in zip(FixedendForce, FixedendForcei)]
        MemberForce = np.round(MemberForce - FixedendForce,2)

        return MemberForce
    
    def MemberForceGlobal(self,MemberNumber):
        
        MemberForce = self.MemberForceLocal(MemberNumber)
        MemberForceGlobal = np.dot(np.transpose(self.Members[self.MemberNo-1].Transformation_Matrix()),MemberForce)

        return MemberForceGlobal
    
    
    def MemberBMD(self, MemberNumber):
        
        self.MemberNo = int(MemberNumber)
        #MOMENT AND SHEAR FORCE
        if(self.Members[self.MemberNo-1].alpha()>=0):
            fem1=self.MemberForceLocal(self.MemberNo)[2] #Fixed End Moment
            fem2=self.MemberForceLocal(self.MemberNo)[5]
        else:
            fem1=self.MemberForceLocal(self.MemberNo)[5]
            fem2=self.MemberForceLocal(self.MemberNo)[2]
        print(fem1,fem2)

        amp=0
        abcd1=[]
        for l in range(1001):
            abcd1.append(0)
        abcd2=[]
        abcd3=[]
        abcd4=[]
        self.amplist=[]
        for a in range(len(self.Loads)):
            if(int(self.Loads[a].AssignedTo[-1]) == self.MemberNo):
                if(self.Members[self.MemberNo-1].alpha()>=0):
                    abcd1=[abcd1[m]+self.Loads[a].EquivalentLoad()['FreeMoment'][m] for m in range(1000)]
                else:
                    abcd1=[abcd1[m]-self.Loads[a].EquivalentLoad()['FreeMoment'][m] for m in range(1000)]
        while(amp<self.Members[self.MemberNo-1].length()):
            mapi=(amp/self.Members[self.MemberNo-1].length()*(-fem2-fem1))+fem1
            abcd2.append(mapi)
            self.amplist.append(amp)
            amp=amp+self.Members[self.MemberNo-1].length()/999 
        abcd3=[abcd1[n]+abcd2[n] for n in range(0,1000)]
        for i in range(0,999):
            ax=(abcd3[i+1]-abcd3[i])/(self.amplist[i+1]-self.amplist[i])
            abcd4.append(ax)
        self.MemberMoment = abcd3
        self.MemberShear = abcd4

        return self.MemberMoment
    
    def MemberSFD(self, MemberNumber):
        self.MemberBMD(MemberNumber)
        
        return self.MemberShear
    
    def MemberAmplitude(self, MemberNumber):
        
        return self.amplist
    
    def MemberNFD(self, MemberNumber):
        return None
    
    def PlotMemberBMD(self, MemberNumber):
        
        self.MemberNo = int(MemberNumber)
        x_max=int(self.Members[self.MemberNo-1].length())
        y_m_max = int(max(self.MemberBMD(self.MemberNo)) * 2)
        y_m_min = int(min(self.MemberBMD(self.MemberNo)) * 2)
        
        if y_m_max == 0:
            y_m_max = 5
        if y_m_min == 0:
            y_m_min = -5
        if y_m_max == y_m_min:
            y_m_max = abs(y_m_max)
            y_m_min = -abs(y_m_min)
        
        c = self.MemberAmplitude(self.MemberNo)
        d = self.MemberBMD(self.MemberNo)
        g = [0, self.Members[self.MemberNo-1].length()]
        h = [0, 0]
        
        plt.figure(figsize=(8, 5))
        plt.plot(c, d, label="Bending Moment", color='red', linewidth=1.5)
        plt.plot(g, h, label="Baseline", color='black', linewidth=1.5, linestyle='dashed')
        
        plt.xlabel('Distance (Meter)')
        plt.ylabel('Bending Moment (kNm)')
        plt.xticks(range(0, x_max + 1, max(1, round(self.Members[self.MemberNo-1].length() / 10))))
        plt.yticks(range(y_m_min, y_m_max + 1, max(1, round((abs(y_m_max) + abs(y_m_min)) / 10))))
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.title(f'Moment Diagram for Member {self.MemberNo}')
        plt.show()


class SecondOrderGlobalResponse(Model):

    def NormalForce(self):

        NormalForceList = []
        for i in range(0,self.Members):
            MemberNumber = i+1
            NormalForceList.append(self.MemberForceLocal(MemberNumber)()[0])

        return NormalForceList
    
    def SecondOrderGlobalStiffnessMatrix(self, NormalForceList):
        C1=[]
        for Mc in self.TotalDoF():
            R1=[]
            for Mr in self.TotalDoF():
                y=0
                for mn in range(0,self.NoMembers):
                    for mr in range(0,6):
                        if(self.Members[mn].DoFNumber()[mr]==Mr):
                            for mc in range(0,6):
                                if(self.Members[mn].DoFNumber()[mc]==Mc):
                                    x=self.Members[mn].Second_Order_Global_Stiffness_Matrix_1(NormalForceList[mn])[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return C1
    
    def SecondOrderGlobalStiffnessMatrixCondensed(self, NormalForceList):
        C1=[]
        for Mc in self.UnConstrainedDoF():
            R1=[]
            for Mr in self.UnConstrainedDoF():
                y=0
                for mn in range(0,self.NoMembers):
                    for mr in range(0,6):
                        if(self.Members[mn].DoFNumber()[mr]==Mr):
                            for mc in range(0,6):
                                if(self.Members[mn].DoFNumber()[mc]==Mc):
                                    x=self.Members[mn].Second_Order_Global_Stiffness_Matrix(NormalForceList[mn])[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return C1
    
    def SecondOrderGlobalStiffnessMatrixCondensedA21(self, NormalForceList):
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
                                    x=self.Members[mn].Second_Order_Global_Stiffness_Matrix_1(NormalForceList[mn])[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return C1
    
    def SecondOrderDisplacement(self, iteration_steps):

        NoMem = len(self.Members)

        #1st iteration
        FirstOderDisplacement = np.dot((np.linalg.inv(np.array(self.GlobalStiffnessMatrixCondensed()))),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Global_Stiffness_Matrix_2", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(MemberForceLocal[0])

        for j in range(0,iteration_steps):

            SecondOrderDisplacement=np.dot((np.linalg.inv(np.array(self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)))),self.ForceVector())
            DisplacementDict = Computer.ModelDisplacementList_To_Dict(SecondOrderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
            
            NorForList1 = NorForList
            NorForList=[]
            for i in range(NoMem):
                SecondOrderMemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
                SecondOrderMemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("Second_Order_Global_Stiffness_Matrix", i+1, self.Members, SecondOrderMemberDisplacement, self.Loads, NorForList1[i])
                NorForList.append(SecondOrderMemberForceLocal[0])
        
        return SecondOrderDisplacement

    def BucklingEigenLoad(self):

        #preparing variables for computation
        NormaldirectionDOFList = [i.dof_x for i in self.Points]
        NoMem = len(self.Members)

        new_list1 = [item for item in self.UnConstrainedDoF() if item not in NormaldirectionDOFList]
        new_list2=self.ConstrainedDoF() + NormaldirectionDOFList
        Detlist=im5.StaticCondensation(self.UnConstrainedDoF(),self.ConstrainedDoF(),ScOGSM)
        #print("Determinant",np.linalg.det(Detlist))
        AllScOGBMM=im7.LSMGlobal2ndOrder(NoMem,GeoMemProp,MatMemPropList,MemSuppNum,self.NormalForce())[2]
        AllScOGSMM_1st_Order=im7.LSMGlobal1stOrder(NoMem,GeoMemProp,MatMemPropList,MemSuppNum)

        gr_buck=new_list1+new_list2
        BGSMConden=im5.StaticCondensation(new_list1,new_list2,im2.GlobalStiffMatr(AllScOGBMM,gr_buck,NoMem))
        BGSMM_1st_Ord_condensed=im5.StaticCondensation(new_list1,new_list2,im2.GlobalStiffMatr(AllScOGSMM_1st_Order,gr_buck,NoMem))
        criticalload , mode=eig(BGSMM_1st_Ord_condensed,BGSMConden)

        return None


class Senstivity(GlobalResponse):

    def AxialMemberSensitivity(self,MemberNumber,scale):

        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].area += scale
                
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))/scale
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def BendingMemberSensitivity(self,MemberNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))/scale
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def MaterialSensitivity(self,MemberNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].youngs_modulus += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))/scale
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def NodeXSensitivity(self,NodeNumber,scale):

        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].xcoordinate += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))/scale
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def NodeYSensitivity(self,NodeNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].ycoordinate += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))/scale
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def GlobalShapeSensitivity(self,SensitivityType):
        return None
    
    def GlobalSizeSensitivity(self,SensitivityType):
        sensitivities = []
        for i in range(len(self.Members)):
            # Calculate the sensitivity for each member
            if SensitivityType == "Axial":
                sensitivity = self.AxialMemberSensitivity(i+1, 1e-6)  # Using a small scale factor 
            elif SensitivityType == "Bending":
                sensitivity = self.BendingMemberSensitivity(i+1, 1e-6)  # Using a small scale factor 
            elif SensitivityType == "Material":
                sensitivity = self.MaterialSensitivity(i+1, 1e-6)  # Using a small scale factor 
            else:
                raise ValueError("Unsupported SensitivityType. Currently, only 'Bending' is supported.")
            sensitivities.append(sensitivity)
        return sensitivities
    
    def PlotSensitivity(self,SensitivityType):
        sensitivities = self.GlobalSizeSensitivity(SensitivityType)
        self.PlotGlobalModel(sensitivities)




#Model Parts - Basic essential for building a model
Points = [Node(Node_Number=1,xcoordinate=0,ycoordinate=0,Support_Condition="Hinged Support"),
        Node(Node_Number=2,xcoordinate=10,ycoordinate=0,Support_Condition="Hinged Support"),
        Node(Node_Number=3,xcoordinate=20,ycoordinate=0,Support_Condition="Hinged Support"),
        ] 

Members = [Member(Beam_Number=1,Start_Node=Points[0],End_Node=Points[1],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
          Member(Beam_Number=2,Start_Node=Points[1],End_Node=Points[2],Area=1,Youngs_Modulus=1,Moment_of_Inertia=1),
           ]

Loads = [NeumanBC(type="PL",Magnitude=5,Distance1=5,AssignedTo="Member 1", Members = Members),
        NeumanBC(type="UDL",Magnitude=5,Distance1=0, Distance2 = 10, AssignedTo="Member 2", Members = Members)
        ] 



#main Model part - Main mode part includes sub model part
Model1 = Model(Points = Points, Members = Members, Loads = Loads)
GlobalRes1 = GlobalResponse(Points = Points, Members = Members, Loads = Loads)
NodalRes1 = NodalResponse(Points = Points, Members = Members, Loads = Loads)
MemberRes1 = MemberResponse(Points = Points, Members = Members, Loads = Loads)
Sensitivity1 = Senstivity(Points = Points, Members = Members, Loads = Loads)
SecondOrderResponse1 = SecondOrderGlobalResponse(Points = Points, Members = Members, Loads = Loads)


Model1.PlotGlobalModel()
print("Node1",NodalRes1.NodeForce(1))
print("Node3",NodalRes1.NodeForce(3))
print("NOde2",NodalRes1.NodeDisplacement(2))
print("mem1",MemberRes1.MemberForceGlobal(1))
print("mem2",MemberRes1.MemberForceGlobal(2))
print("mem1",MemberRes1.MemberForceLocal(1))
print("mem2",MemberRes1.MemberForceLocal(2))
print("mem1",MemberRes1.MemberBMD(1))
print(len(MemberRes1.MemberBMD(1)))
MemberRes1.PlotMemberBMD(1)






#LSM=Stiffness_Matrix(Members[0])
#fom=LSM.First_Order_Global_Stiffness_Matrix_1()
#print(Members[1].First_Order_Global_Stiffness_Matrix_1())
#print(Points[1].DoF())
#print(fom)
#print(Members[2].DoFNumber())
#print(Model1.UnConstrainedDoF())
#print(Model1.GlobalStiffnessMatrix())
#print(Model1.GlobalStiffnessMatrixCondensedA21())
#print(Loads[1].EquivalentLoad())
#print(Model1.ForceVector())
#print(GlobalRes.DisplacementVector())
#print(GlobalRes.SupportForcesVector())

#print("SupportForces", NodalRes.NodeForce(NodeNumber = 1))
#print("NodeDisplacement", NodalRes.NodeDisplacement(NodeNumber = 2))
#print(MemberRes.MemberDisplacement(MemberNumber = 1))
#print(MemberRes.MemberForce(MemberNumber = 1))
#print(MemberRes.PlotMemberBMD(MemberNumber = 1))
#print(MemberRes.PlotMemberBMD(2))
#print(Sensitivity1.NodeYSensitivity(1,0.001))
#print(Sensitivity1.GlobalSizeSensitivity("Bending"))



"""
Wrong

1. Hinged Joint Support is behaving as hinged support

"""
"""
Improvements

1. Include number of finit element forces as a variable

"""



