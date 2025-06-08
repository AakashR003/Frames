
import numpy as np

# Finite element Division


class Node():
    
    i=1
    titam = 1
    titay = 1
    def __init__ (self, Node_Number, xcoordinate, ycoordinate, Support_Condition):
        self.node_number = Node_Number
        self.xcoordinate = xcoordinate
        self.ycoordinate = ycoordinate
        self.support_condition = Support_Condition
        self.additional_dof_tita=[]

        self.LBNodeXSensitivity = None
        self.LBNodeYSensitivity = None

        if self.support_condition in ["Hinged Support", "Fixed Support", "Rigid Joint", "Roller in X-plane", "Roller in Y-plane"]:
            self.dof_x=(self.node_number)*3-2
            self.dof_y=(self.node_number)*3-1
            self.dof_tita=(self.node_number)*3

        elif self.support_condition == "Hinge Joint" :
            self.dof_x=(self.node_number)*3-2
            self.dof_y=(self.node_number)*3-1
            self.dof_tita=300000
            #self.additional_dof_tita=[]

        elif self.support_condition == "Glided Support" :
            self.dof_x=(self.node_number)*3-2
            self.dof_y=200000
            self.additional_dof_y = []
            self.dof_tita=(self.node_number)*3

        elif self.support_condition == "Hinged Joint Support" or self.support_condition == "Roller in X-plane-Hinge" :
            self.dof_x=(self.node_number)*3-2
            self.dof_y=(self.node_number)*3-1
            self.dof_tita=300000
            #self.additional_dof_tita=[]
        
        else:
            raise ValueError(f"Unsupported support condition: '{self.support_condition}'")
        
        
    def DoF(self):
        self.check1 = [self.dof_x, self.dof_y, self.dof_tita]    
        return [self.dof_x, self.dof_y, self.dof_tita]


class Couple_Nodes():
        
    def __init__ (self, Main_Node, Dependent_Node, xDof = True, yDof = True, RotationDof = True):
        self.Main_Node = Main_Node
        self.Dependent_Node = Dependent_Node
        self.xDof = xDof
        self.yDof = yDof
        self.titaDof = RotationDof

        if self.xDof == True:
            self.Dependent_Node.dof_x  = self.Main_Node.dof_x

        if self.yDof == True:
            self.Dependent_Node.dof_y  = self.Main_Node.dof_y

        if self.titaDof == True:
            self.Dependent_Node.dof_tita  = self.Main_Node.dof_tita
        
    def CoupledDoF(self):
        return [self.Main_Node_Number, self.Dependent_Node_Number, self.xDof, self.yDof, self.titaDof]   
            
    
class Member():
    
    def __init__ (self, Beam_Number, Start_Node, End_Node, Area, Youngs_Modulus, Moment_of_Inertia, Density = 7850 ):
        
        self.Beam_Number = Beam_Number
        self.Start_Node = Start_Node
        self.End_Node = End_Node
        
        self.area=Area
        self.youngs_modulus = Youngs_Modulus
        self.moment_of_inertia = Moment_of_Inertia
        self.Density = Density

        #calculated properties
        self.LBSensitivityBend = None
        self.LBSensitivityAxial = None

        #Additional for Hinge Joint
        if self.Start_Node.support_condition in ["Hinge Joint", "Hinged Joint Support", "Roller in X-plane-Hinge"]:
            self.Start_Node.additional_dof_tita.append(300000 + self.Beam_Number*2-1)

        if self.End_Node.support_condition in ["Hinge Joint", "Hinged Joint Support", "Roller in X-plane-Hinge"]:
            self.End_Node.additional_dof_tita.append(300000 + self.Beam_Number*2)
            #ends Here

    def length(self):
        x=((self.End_Node.xcoordinate-self.Start_Node.xcoordinate)**2 +
           (self.End_Node.ycoordinate-self.Start_Node.ycoordinate)**2)**0.5
        return x
    
    def alpha(self):
        return (self.End_Node.xcoordinate-self.Start_Node.xcoordinate)/self.length()
   
    def beta(self):
        return (self.End_Node.ycoordinate-self.Start_Node.ycoordinate)/self.length()
    
    def DoFNumber(self):

        #Additional for Hinge Joint
        if self.Start_Node.support_condition in ["Hinge Joint", "Hinged Joint Support", "Roller in X-plane-Hinge"]:
            self.Start_Node.dof_tita = 300000 + self.Beam_Number*2-1

        if self.End_Node.support_condition in ["Hinge Joint", "Hinged Joint Support", "Roller in X-plane-Hinge"]:
            self.End_Node.dof_tita = 300000 + self.Beam_Number*2
        #Ends Here
        
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

    def Second_Order_Reduction_Matrix_1(self, NormalForce):

        length=self.length()
        ReductionMatrix=[[1/length,  0,           0,           -1/length,  0,            0           ],
                        [0,          6/5/length,  1/10,         0,         -6/5/length,  1/10        ],
                        [0,          1/10,        2/15*length,  0,         -1/10,        -1/30*length],
                        [-1/length,  0,           0,            1/length,  0,            0           ],
                        [0,          -6/5/length,  -1/10,        0,         6/5/length,   -1/10        ],
                        [0,          1/10,        -1/30*length, 0,         -1/10,        2/15*length ]]
        # reference -2 Matrix structural Analysis and Stability by William McGuire and Richard H. Gallagher
        StiffnessReductionMatrix = np.array(ReductionMatrix) * NormalForce
        
        return StiffnessReductionMatrix

    def Second_Order_Global_Reduction_Matrix_1(self, NormalForce):
        return np.transpose(self.Transformation_Matrix()) @ np.array(self.Second_Order_Reduction_Matrix_1(NormalForce)) @ np.array(self.Transformation_Matrix())

    def Second_Order_Local_Stiffness_Matrix_1(self,NormalForce):

        SecondOrderLocalStiffnessMatrix = np.array(self.First_Order_Local_Stiffness_Matrix_1()) + self.Second_Order_Reduction_Matrix_1(NormalForce)
        return SecondOrderLocalStiffnessMatrix
    
    def Second_Order_Global_Stiffness_Matrix_1(self, NormalForce):

        return np.transpose(self.Transformation_Matrix()) @ np.array(self.Second_Order_Local_Stiffness_Matrix_1(NormalForce)) @ np.array(self.Transformation_Matrix())

    def Second_Order_Reduction_Matrix_2(self, NormalForce):

        length=self.length()
        ReductionMatrix=[[1/length,    0,           0,           -1/length,    0,            0           ],
                         [0,           6/5/length,  -1/10,       0,           -6/5/length,  -1/10       ],
                         [0,          -1/10,        2/15*length, 0,            1/10,        -1/30*length],
                         [0,           0,           0,           0,            0,            0           ],
                         [0,          -6/5/length,   1/10,       0,            6/5/length,   1/10        ],
                         [0,          -1/10,       -1/30*length, 0,            1/10,        2/15*length ]]
        # reference - 1.Stability Notes
        StiffnessReductionMatrix = np.array(ReductionMatrix) * NormalForce
        
        return StiffnessReductionMatrix
    
    def Second_Order_Global_Reduction_Matrix_2(self, NormalForce):
        return np.transpose(self.Transformation_Matrix()) @ np.array(self.Second_Order_Reduction_Matrix_2(NormalForce)) @ np.array(self.Transformation_Matrix())

    def Second_Order_Local_Stiffness_Matrix_2(self,NormalForce):

        SecondOrderLocalStiffnessMatrix = np.array(self.First_Order_Local_Stiffness_Matrix_2()) + self.Second_Order_Reduction_Matrix_2(NormalForce)
        return SecondOrderLocalStiffnessMatrix
    
    def Second_Order_Global_Stiffness_Matrix_2(self, NormalForce):

        return np.transpose(self.Transformation_Matrix()) @ np.array(self.Second_Order_Local_Stiffness_Matrix_2(NormalForce)) @ np.array(self.Transformation_Matrix())

    def Local_Mass_Matrix(self):

        L = self.length()
        mu = self.area * self.Density

        MassMatrix = [
                    [140 * mu * L / 420, 0, 0, 70 * mu * L / 420, 0, 0],
                    [0, 156 * mu * L / 420, 22 * mu * L**2 / 420, 0, 54 * mu * L / 420, -13 * mu * L**2 / 420],
                    [0, 22 * mu * L**2 / 420, 4 * mu * L**3 / 420, 0, 13 * mu * L**2 / 420, -3 * mu * L**3 / 420],
                    [70 * mu * L / 420, 0, 0, 140 * mu * L / 420, 0, 0],
                    [0, 54 * mu * L / 420, 13 * mu * L**2 / 420, 0, 156 * mu * L / 420, -22 * mu * L**2 / 420],
                    [0, -13 * mu * L**2 / 420, -3 * mu * L**3 / 420, 0, -22 * mu * L**2 / 420, 4 * mu * L**3 / 420]
                    ]


        return MassMatrix
    
    def Global_Mass_Matrix(self):
        
        return np.transpose(self.Transformation_Matrix()) @ np.array(self.Local_Mass_Matrix()) @ np.array(self.Transformation_Matrix())


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

