import numpy as np
import matplotlib.pyplot as plt


try:
    from .Model import Model
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Functions import max_nested
    from .FirstOrderResponse import FirstOrderGlobalResponse
    from .SecondOrderResponse import SecondOrderGlobalResponse
    from .Strain_Energy import StrainEnergy
except:
    from Model import Model
    from StructuralElements import Node, Member
    from Computer import Computer
    from Functions import max_nested
    from FirstOrderResponse import FirstOrderGlobalResponse
    from SecondOrderResponse import SecondOrderGlobalResponse
    from Strain_Energy import StrainEnergy

class Senstivity(FirstOrderGlobalResponse):

    def AxialMemberSensitivity(self,MemberNumber,scale):

        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].area += scale
                
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
        sensitivity = 0.5 * np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def BendingMemberSensitivity(self,MemberNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_BendingStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
        sensitivity = 0.5 * np.dot(np.dot(np.transpose(self.DisplacementVector()),d_BendingStiffness_ds),self.DisplacementVector())
        
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        return sensitivity
    
    def MaterialSensitivity(self,MemberNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].youngs_modulus += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def NodeXSensitivity(self,NodeNumber,scale):

        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].xcoordinate += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def NodeYSensitivity(self,NodeNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].ycoordinate += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
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

class EigenSensitivity(SecondOrderGlobalResponse):
    """
    This class is a placeholder for second-order sensitivity analysis.
    Currently, it does not implement any specific methods.
    """

    def compute_SM(self):
        NoMem = len(self.Members)
        SM = self.GlobalStiffnessMatrixCondensed()

        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        GeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        return SM, GeoSM
    
    def compute_sensitivity(self, UnModifiedSM, UnModifiedGeoSM, ModifiedSM, ModifiedGeoSM, Eigen, EigenModeNo):
        #Modified system
        ModifiedSM, ModifiedGeoSM = self.compute_SM()

        if Eigen is None:
            Eigen_mode = self.BucklingEigenLoad(Solver = "eigsh")[2][:,(EigenModeNo-1)]
            Eigen_load = self.BucklingEigenLoad(Solver = "eigsh")[1][EigenModeNo-1]
        else:
            Eigen_mode = Eigen[2][:,(EigenModeNo-1)]
            Eigen_load = Eigen[1][EigenModeNo-1]

        d_K_ds = (np.array(ModifiedSM) - np.array(UnModifiedSM))
        d_Kg_ds = (np.array(ModifiedGeoSM) - np.array(UnModifiedGeoSM))

        sensitivity = np.dot(np.dot(np.transpose(Eigen_mode),(d_K_ds - Eigen_load* d_Kg_ds)), Eigen_mode)
        return sensitivity        
    
    
    def NodeXSensitivity(self,NodeNumber,scale, Eigen = None, EigenModeNo = 1):

        #unModified system
        UnModifiedSM, UnModifiedGeoSM = self.compute_SM()

        #small perturbation
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].xcoordinate += scale
        
        #Modified system
        ModifiedSM, ModifiedGeoSM = self.compute_SM()

        # Compute the sensitivity
        sensitivity = self.compute_sensitivity(UnModifiedSM, UnModifiedGeoSM, ModifiedSM, ModifiedGeoSM, Eigen, EigenModeNo)
        #update sensitivity

        return sensitivity
    
    def NodeYSensitivity(self,NodeNumber,scale, Eigen = None, EigenModeNo = 1):
        #unModified system
        UnModifiedSM, UnModifiedGeoSM = self.compute_SM()

        #small perturbation
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].ycoordinate += scale
        
        #Modified system
        ModifiedSM, ModifiedGeoSM = self.compute_SM()

        # Compute the sensitivity
        sensitivity = self.compute_sensitivity(UnModifiedSM, UnModifiedGeoSM, ModifiedSM, ModifiedGeoSM, Eigen, EigenModeNo)
        #update sensitivity
        return sensitivity
    
    def BeamBendingSensitivity(self,MemberNumber,scale, Eigen = None, EigenModeNo = 1):

        UnModifiedSM, UnModifiedGeoSM = self.compute_SM()
        
        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
                
        #Modified system
        ModifiedSM, ModifiedGeoSM = self.compute_SM()

        # Compute the sensitivity
        sensitivity = self.compute_sensitivity(UnModifiedSM, UnModifiedGeoSM, ModifiedSM, ModifiedGeoSM, Eigen, EigenModeNo)
        #update sensitivity
        self.Members[MemberNumber-1].LBSensitivityBend = sensitivity

        return sensitivity

    def BeamAxialSensitivity(self,MemberNumber,scale, Eigen = None, EigenModeNo = 1):

        UnModifiedSM, UnModifiedGeoSM = self.compute_SM()
        
        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].area += scale
                
        #Modified system
        ModifiedSM, ModifiedGeoSM = self.compute_SM()

        # Compute the sensitivity
        sensitivity = self.compute_sensitivity(UnModifiedSM, UnModifiedGeoSM, ModifiedSM, ModifiedGeoSM, Eigen, EigenModeNo)
        #update sensitivity

        return sensitivity
    
    def GlobalSecondOrderBendingSensitivity(self, EigenModeNo = 1):
        Eigen = self.BucklingEigenLoad(Solver = "eigsh")
        for i in range(len(self.Members)):
            # Calculate the sensitivity for each member
            self.BeamBendingSensitivity(i+1, 1e-10, Eigen = Eigen, EigenModeNo = EigenModeNo)
        sensitivity_values = [member.LBSensitivityBend for member in self.Members]

        return sensitivity_values
    
    def GlobalSecondOrderAxialSensitivity(self, EigenModeNo = 1):
        Eigen = self.BucklingEigenLoad(Solver = "eigsh")
        for i in range(len(self.Members)):
            # Calculate the sensitivity for each member
            self.BeamBendingSensitivity(i+1, 1e-10, Eigen = Eigen, EigenModeNo = EigenModeNo)
        sensitivity_values = [member.LBSensitivityBend for member in self.Members]

        return sensitivity_values
    


    def PlotGlobalSecondOrderMemberSensitivity(self, Sensitivity_type = "Bending", EigenModeNo = 1):

        if Sensitivity_type == "Bending":
            sensitivity_values = self.GlobalSecondOrderBendingSensitivity(EigenModeNo = EigenModeNo)
        if Sensitivity_type == "Axial":
            sensitivity_values = self.GlobalSecondOrderAxialSensitivity(EigenModeNo = EigenModeNo)
        fig, ax = plt.subplots(figsize=(12, 8))
        computer_instance = Computer()
        computer_instance.PlotStructuralElements(ax, self.Members, self.Points, ShowNodeNumber=False, sensitivities= sensitivity_values)
                
        ax.set_title(f"Buckling Eigenmode Size Sensitivity", fontsize=16)
        ax.axis('equal')
        plt.show()
        return None
    
    def PlotGlobalSecondOrderNodeSensitivity(self, EigenModeNo = 1):
        
        return None

class SecondOrderSensitivity(StrainEnergy):
    """
    This class is a placeholder for second-order sensitivity analysis.
    Currently, it does not implement any specific methods.
    """

    def GlobalSecondOrderBendingSensitivity(self, MemberNumber, scale=1e-6):

        DisplacementVector = self.SecondOrderDisplacementVector(iteration_steps=5)

        StiffnessMatrix_Unmodified, self.SecondOderDisplacement = self.SecondOrderDisplacementVector(iteration_steps = 5, ReturnSM = True)

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
        
        StiffnessMatrix_Modified, self.SecondOderDisplacement = self.SecondOrderDisplacementVector(iteration_steps = 5, ReturnSM = True)

        dk_ds = (np.array(StiffnessMatrix_Modified) - np.array(StiffnessMatrix_Unmodified))        
        sensitivity = 0.5*np.dot(np.dot(np.transpose(DisplacementVector), dk_ds), DisplacementVector)

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale


        return sensitivity
    
    def FirstOrderPartBendingSensitivity(self, MemberNumber, scale=1e-6):
        """
        This method calculates the first part of the second-order bending sensitivity for a member.
        It perturbs the member's moment of inertia and computes the sensitivity.
        """
        DisplacementVector = self.SecondOrderDisplacementVector(iteration_steps=5)

        StiffnessMatrix_Unmodified, self.SecondOderDisplacement = self.SecondOrderDisplacementVector(iteration_steps = 5, ReturnSM = True)
        FirstOrderResponse1 = FirstOrderGlobalResponse(Points=self.Points, Members=self.Members, Loads=self.Loads)
        FirstOrder_stiffness_Unmodified = FirstOrderResponse1.GlobalStiffnessMatrixCondensed()

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
        
        StiffnessMatrix_Modified, self.SecondOderDisplacement = self.SecondOrderDisplacementVector(iteration_steps = 5, ReturnSM = True)

        dk_ds = (np.array(StiffnessMatrix_Modified) - np.array(StiffnessMatrix_Unmodified))
        sensitivity = 0.5*np.dot(np.dot(np.transpose(DisplacementVector), dk_ds), DisplacementVector)

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        return sensitivity
    
    def SecondOrderPartBendingSensitivity(self, MemberNumber, scale=1e-6):

        SecondOrder_stiffness = self.SecondOrderDisplacementVector(iteration_steps=8, ReturnSM=True)
        SecondOrder_DisplacementVector = self.SecondOrderDisplacementVector(iteration_steps=8)
        #SecondOrder_DisplacementVector_orthogonal = Computer.OrthogonalSolver(SecondOrder_DisplacementVector, SMatrix=SecondOrder_stiffness)
        

        FirstOrderResponse1 = FirstOrderGlobalResponse(Points=self.Points, Members=self.Members, Loads=self.Loads)

        FirstOrder_stiffness = FirstOrderResponse1.GlobalStiffnessMatrixCondensed()
        FirstOrder_DisplacementVector = FirstOrderResponse1.DisplacementVector()
        #FirstOrder_DisplacementVector_orthogonal = Computer.OrthogonalSolver(FirstOrder_DisplacementVector, SMatrix=FirstOrder_stiffness)

        FirstOrder_FlexibilityMatrix = Computer.FlexibilityMatrixSolver(FirstOrder_stiffness)
        SecondOrder_FlexibilityMatrix = Computer.FlexibilityMatrixSolver(SecondOrder_stiffness)

        self.SetModifiedValues()
        DifferenceStiffness_matrix1 = SecondOrder_stiffness - FirstOrder_stiffness
        #DifferenceStiffness_matrix1 = Computer.FlexibilityMatrixSolver(self.ModifiedFlexibilityMatrix)
        Difference_DisplacementVector1 = SecondOrder_DisplacementVector - FirstOrder_DisplacementVector
        Difference_FlexibilityMatrix1 = SecondOrder_FlexibilityMatrix - FirstOrder_FlexibilityMatrix
        Orthogonal_Difference_Flexibility_Matrix1 = Computer.OrthogonalSolver(Difference_FlexibilityMatrix1, SMatrix=Difference_FlexibilityMatrix1)
        ForceVector_Orthogonal1 = Computer.OrthogonalSolver(self.ForceVector(), SMatrix=Orthogonal_Difference_Flexibility_Matrix1)
        Orthogonal_Modified_Difference_Flexibility_Modified_matrix1 = Orthogonal_Difference_Flexibility_Matrix1 / ForceVector_Orthogonal1[:, np.newaxis] # U = Flexbility * F_square

        # WIth modified Flexibility, Flexibility *FSqaure = Nonlinear part of U. Hence Flexibility * F is not Nonlinear part of U. Hence computing it
        Orthogonal_Modified_Displacement_with_modified_Flexibility = Computer.DirectDisplacementSolver(Orthogonal_Modified_Difference_Flexibility_Modified_matrix1, self.ForceVector())
        Modified_Displacement_with_modified_Flexibility = Computer.OrthogonalSolver(Orthogonal_Modified_Displacement_with_modified_Flexibility, SMatrix=Difference_FlexibilityMatrix1, Back=True)
        Orthogonal_Difference_DisplacementVector1 = Computer.OrthogonalSolver(Difference_DisplacementVector1, SMatrix=Difference_FlexibilityMatrix1)
        
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale

        FirstOrderResponse2 = FirstOrderGlobalResponse(Points=self.Points, Members=self.Members, Loads=self.Loads)
        SecondOrder_DisplacementVector2 = self.SecondOrderDisplacementVector(iteration_steps=8)
        FirstOrder_DisplacementVector2 = FirstOrderResponse2.DisplacementVector()
        FirstOrder_stiffness2 = FirstOrderResponse2.GlobalStiffnessMatrixCondensed()
        SecondOrder_stiffness2 = self.SecondOrderDisplacementVector(iteration_steps=8, ReturnSM=True)

        FirstOrder_FlexibilityMatrix2 = Computer.FlexibilityMatrixSolver(FirstOrder_stiffness2)
        SecondOrder_FlexibilityMatrix2 = Computer.FlexibilityMatrixSolver(SecondOrder_stiffness2)
        
        self.SetModifiedValues() 
        DifferenceStiffness_matrix2 = SecondOrder_stiffness2 - FirstOrder_stiffness2
        #DifferenceStiffness_matrix2 = Computer.FlexibilityMatrixSolver(self.ModifiedFlexibilityMatrix)
        Difference_DisplacementVector2 = SecondOrder_DisplacementVector2 - FirstOrder_DisplacementVector2
        Difference_FlexibilityMatrix2 = SecondOrder_FlexibilityMatrix2 - FirstOrder_FlexibilityMatrix2
        Orthogonal_Difference_Flexibility_Matrix2 = Computer.OrthogonalSolver(Difference_FlexibilityMatrix2, SMatrix=Difference_FlexibilityMatrix2)
        ForceVector_Orthogonal2 = Computer.OrthogonalSolver(self.ForceVector(), SMatrix=Orthogonal_Difference_Flexibility_Matrix2)
        Orthogonal_Modified_Difference_Flexibility_Modified_matrix2 = Orthogonal_Difference_Flexibility_Matrix2 / ForceVector_Orthogonal2[:, np.newaxis] # U = Flexbility * F_square
        


        d_BendingStiffness_ds = DifferenceStiffness_matrix2 - DifferenceStiffness_matrix1
        print(d_BendingStiffness_ds)
        d_BendingStiffness_ds = np.linalg.inv(Orthogonal_Modified_Difference_Flexibility_Modified_matrix1) - np.linalg.inv(Orthogonal_Modified_Difference_Flexibility_Modified_matrix2)

        #sensitivity = 2/3 * np.dot(np.dot(np.transpose(Modified_Displacement_with_modified_Flexibility),d_BendingStiffness_ds),Difference_DisplacementVector1)
        
        sensitivity = 2/3 * np.dot(np.dot(np.transpose(Orthogonal_Modified_Displacement_with_modified_Flexibility),d_BendingStiffness_ds),Difference_DisplacementVector1)

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale
        
        
        return sensitivity
    
    def try2SecondOrderPartSensitivity(self, MemberNumber, scale = 1e-6):

        self.SetModifiedValues() 
        one = np.dot(self.ModifiedFlexibilityMatrix, self.ForceVector())
        one_1= self.CalculateApproximatedValueDisplacement(ReturnNonlinearDisplacement=True)
        one_2 = self.SecondOrderDisplacementVector(iteration_steps = 5)
        one_3 = self.ForceVector()

        three = self.CalculateApproximatedValueDisplacement(ReturnNonlinearDisplacement=True)

        k1 = np.linalg.inv(Computer.OrthogonalSolver(self.ModifiedFlexibilityMatrix, SMatrix=self.CalculateK2ndOrderMinusKLinear, Back=True))
        k1_1 = self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)[0] - self.GlobalStiffnessMatrixCondensed()
        k1_2 = self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)[0]
        k1_3 = self.GlobalStiffnessMatrixCondensed()
        k1_4 = np.linalg.inv(self.ModifiedFlexibilityMatrix)
        k1_5 = self.ModifiedFlexibilityMatrix
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
        
        self.SetModifiedValues()
        k2 = np.linalg.inv(Computer.OrthogonalSolver(self.ModifiedFlexibilityMatrix, SMatrix=self.CalculateK2ndOrderMinusKLinear, Back=True))
        k2_1 = self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)[0] - self.GlobalStiffnessMatrixCondensed()
        k2_2 = self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)[0]
        k2_3 = self.GlobalStiffnessMatrixCondensed()
        k2_4 = np.linalg.inv(self.ModifiedFlexibilityMatrix)
        k2_5 = self.ModifiedFlexibilityMatrix



        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        two = k2- k1
        two_1 = k2_1 - k1_1
        two_2 = k2_2 - k1_2
        two_3 = k2_3 - k1_3
        two_4 = k2_4 - k1_4
        two_5 = k2_5 - k1_5

        Senstivity = 1/3 * np.dot(np.dot(np.transpose(one_1), two), three)

        return Senstivity
    
    def try3(self, MemberNumber, linearsensitvity, scale = 1e-6 ):
        
        senoverall = self.GlobalSecondOrderBendingSensitivity( MemberNumber, scale) * 2
        print("Rectangle_AAAAAAAAAAAAAAAAAAAAAAAAA", senoverall)
        parabola = self.try2SecondOrderPartSensitivity( MemberNumber, scale)
        
        return senoverall - parabola - linearsensitvity


class FiniteDifferenceSensitivity(Model):
    
    def __init__(self, Points, Members, Loads):
        """
        Initialize the FiniteDifferenceSensitivity class with points, members, and loads.
        """
        super().__init__(Points=Points, Members=Members, Loads=Loads)
        self.FirstOrderResponse1 = FirstOrderGlobalResponse(Points=self.Points, Members=self.Members, Loads=self.Loads)
        self.SecondOrderResponse1 = SecondOrderGlobalResponse(Points=self.Points, Members=self.Members, Loads=self.Loads)

    
    def CalculateObjective(self, LinearElement = True, NonLinearElement = False):
        
        if LinearElement:
            linear_strain_energy = 0.5*np.dot(np.transpose(self.ForceVector()),self.FirstOrderResponse1.DisplacementVector())

        if NonLinearElement:
            linear_strain_energy = 0.5*np.dot(np.transpose(self.ForceVector()),self.SecondOrderResponse1.SecondOrderDisplacementVector(iteration_steps=5))
        
        return linear_strain_energy
    
    def FiniteDifferenceFirstOrderLinearBendingSensitivity(self,MemberNumber, scale=1e-6):
        """
        This method calculates the finite difference sensitivity of the first-order bending stiffness for a member.
        It perturbs the member's moment of inertia and computes the sensitivity.
        """

        linear_strain_energy = self.CalculateObjective(LinearElement = True, NonLinearElement = False)

        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale

        modified_linear_strain_energy = self.CalculateObjective(LinearElement = True, NonLinearElement = False)

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        #Calculate Sensitivity
        sensitivity = (modified_linear_strain_energy - linear_strain_energy) 
        
        return sensitivity
    
    def FiniteDifferenceSecondOrderLinearBendingSensitivity(self,MemberNumber, scale=1e-6):
        """
        This method calculates the finite difference sensitivity of the first-order bending stiffness for a member.
        It perturbs the member's moment of inertia and computes the sensitivity.
        """

        linear_strain_energy = self.CalculateObjective(LinearElement = False, NonLinearElement = True)

        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale

        modified_linear_strain_energy = self.CalculateObjective(LinearElement = False, NonLinearElement = True)

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        #Calculate Sensitivity
        sensitivity = (modified_linear_strain_energy - linear_strain_energy) 
        
        return sensitivity
    
    def FiniteDifferenceApproximatedSecondOrderBendingSensitivity(self, MemberNumber, scale =1e-6):
        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateApproximatedNonlinearStrainEnergy()
        nonlinearstrain_energy = objective
        
        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale

        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateApproximatedNonlinearStrainEnergy()
        modified_nonlinearstrain_energy = objective

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        #Calculate Sensitivity
        sensitivity = (modified_nonlinearstrain_energy - nonlinearstrain_energy) 
        
        return sensitivity
    
    def FiniteDifferenceSecondOrderNonLinearBendingSensitivity(self,MemberNumber, scale=1e-6):
        """
        This method calculates the finite difference sensitivity of the second-order bending stiffness for a member.
        It perturbs the member's moment of inertia and computes the sensitivity.
        """

        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateFiniteDifferenceNonLinearStrainEnergy()
        nonlinearstrain_energy = objective
        
        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale

        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateFiniteDifferenceNonLinearStrainEnergy()
        modified_nonlinearstrain_energy = objective

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        #Calculate Sensitivity
        sensitivity = (modified_nonlinearstrain_energy - nonlinearstrain_energy) 
        
        return sensitivity
    
    def Members_All_Sensitivity(self, all = True):

        Approximated_list = []
        Linearized_list = []
        Linear_list = []
        Nonlinear_list = []
        for i in range(len(self.Members)):
            Approximated = self.FiniteDifferenceApproximatedSecondOrderBendingSensitivity(MemberNumber=i+1, scale=0.000000000001)
            Linearized = self.FiniteDifferenceSecondOrderLinearBendingSensitivity(MemberNumber = i+1, scale=0.000000000001)
            Linear = self.FiniteDifferenceFirstOrderLinearBendingSensitivity(MemberNumber = i+1, scale=0.000000000001)
            Nonlinear = self.FiniteDifferenceSecondOrderNonLinearBendingSensitivity(MemberNumber = i+1, scale=0.000000000001)

            Approximated_list.append(Approximated)
            Linearized_list.append(Linearized)
            Linear_list.append(Linear)
            Nonlinear_list.append(Nonlinear)
        
        print(Approximated_list)
        print(Linearized_list)
        print(Linear_list)
        print(Nonlinear_list)
        
        # Calculate relative errors for each element
        print("Relative error (w.r.t Nonlinear) for each element:")
        methods = [
            ("Approximated", Approximated_list),
            ("Linearized", Linearized_list),
            ("Linear", Linear_list)
        ]
        n = len(self.Members)
        error_sums = {name: 0.0 for name, _ in methods}
        for i in range(n):
            print(f"Element {i+1}:")
            nonlinear = Nonlinear_list[i]
            for name, values in methods:
                if nonlinear != 0:
                    rel_error = abs((values[i] - nonlinear) / nonlinear) * 100
                else:
                    rel_error = float('inf')
                error_sums[name] += rel_error
                print(f"  {name} vs Nonlinear: {rel_error:.2f}%")
        print("\nOverall error percentage for each method:")
        for name in error_sums:
            overall_error = error_sums[name] / n
            print(f"  {name}: {overall_error:.2f}%")

        # Normalization by the first value of each list
        def normalize_list(lst):
            if lst[0] != 0:
                return [x / lst[0] for x in lst]
            else:
                return [float('nan') for _ in lst]

        Approximated_norm = normalize_list(Approximated_list)
        Linearized_norm = normalize_list(Linearized_list)
        Linear_norm = normalize_list(Linear_list)
        Nonlinear_norm = normalize_list(Nonlinear_list)

        print("\nNormalized lists (by first value of each):")
        print("Approximated_norm:", Approximated_norm)
        print("Linearized_norm:", Linearized_norm)
        print("Linear_norm:", Linear_norm)
        print("Nonlinear_norm:", Nonlinear_norm)

        # Calculate relative errors for normalized lists
        print("\nRelative error (w.r.t normalized Nonlinear) for each element:")
        norm_methods = [
            ("Approximated_norm", Approximated_norm),
            ("Linearized_norm", Linearized_norm),
            ("Linear_norm", Linear_norm)
        ]
        norm_error_sums = {name: 0.0 for name, _ in norm_methods}
        for i in range(n):
            print(f"Element {i+1}:")
            nonlinear = Nonlinear_norm[i]
            for name, values in norm_methods:
                if nonlinear != 0:
                    rel_error = abs((values[i] - nonlinear) / nonlinear) * 100
                else:
                    rel_error = float('inf')
                norm_error_sums[name] += rel_error
                print(f"  {name} vs Nonlinear_norm: {rel_error:.2f}%")
        print("\nOverall error percentage for each normalized method:")
        for name in norm_error_sums:
            overall_error = norm_error_sums[name] / n
            print(f"  {name}: {overall_error:.2f}%")



