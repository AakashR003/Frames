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
                sensitivity = self.AxialMemberSensitivity(i+1, 1e-12)  # Using a small scale factor 
            elif SensitivityType == "Bending":
                sensitivity = self.BendingMemberSensitivity(i+1, 1e-12)  # Using a small scale factor 
            elif SensitivityType == "Material":
                sensitivity = self.MaterialSensitivity(i+1, 1e-12)  # Using a small scale factor 
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
            self.BeamBendingSensitivity(i+1, 1e-12, Eigen = Eigen, EigenModeNo = EigenModeNo)
        sensitivity_values = [member.LBSensitivityBend for member in self.Members]

        return sensitivity_values
    
    def GlobalSecondOrderAxialSensitivity(self, EigenModeNo = 1):
        Eigen = self.BucklingEigenLoad(Solver = "eigsh")
        for i in range(len(self.Members)):
            # Calculate the sensitivity for each member
            self.BeamBendingSensitivity(i+1, 1e-12, Eigen = Eigen, EigenModeNo = EigenModeNo)
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
    This class calculates sensitivities of strain energy with respect to member properties
    """

    def GlobalSecondOrderBendingSensitivity(self, MemberNumber, scale=1e-12):

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
       
    def CalculateSimpsonsSensitivityAllMembers(self, scale =1e-12 ):

        #UnModified
        
        #SecondPart
        F2 = np.array(self.ForceVector())
        StiffnessMatrix_2_Unmodified, U2 = self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)

        #FirstPart
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude / 2
        
        StiffnessMatrix_1_Unmodified, U1 = self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)

        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude * 2

        
        #Modified
        SenstivityList = []
        for Member in self.Members:

            Member.moment_of_inertia += scale
            StiffnessMatrix_2_Modified= self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)[0]

            #FirstPart
            for i in range(len(self.Loads)):
                self.Loads[i].Magnitude = self.Loads[i].Magnitude / 2
            
            StiffnessMatrix_1_Modified= self.SecondOrderDisplacementVector(iteration_steps=5, ReturnSM=True)[0]
            
            for i in range(len(self.Loads)):
                self.Loads[i].Magnitude = self.Loads[i].Magnitude * 2

            difference_k1 = StiffnessMatrix_1_Modified - StiffnessMatrix_1_Unmodified
            difference_k2 = StiffnessMatrix_2_Modified - StiffnessMatrix_2_Unmodified

            FirstPart = 4/6*np.dot(np.dot(np.transpose(Computer.DirectInverseDisplacementSolver(StiffnessMatrix_1_Unmodified, F2)),difference_k1),U1)
            SecondPart = np.dot(np.dot(np.transpose(U2),difference_k2),U2) / 6

            complimentarySensitivity = FirstPart + SecondPart
            rectangleSensitivity = np.dot(np.dot(np.transpose(U2),difference_k2),U2)
            Sensitivity = rectangleSensitivity - complimentarySensitivity
            
            SenstivityList.append(-Sensitivity)
            Member.moment_of_inertia -= scale
        
        return SenstivityList

        


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
    
    def FiniteDifferenceFirstOrderLinearBendingSensitivity(self,MemberNumber, scale=1e-12):
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
    
    def FiniteDifferenceSecondOrderLinearBendingSensitivity(self,MemberNumber, scale=1e-12):
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
    
    def FiniteDifferenceApproximatedSecondOrderBendingSensitivity(self, MemberNumber, scale =1e-12):
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
    
    def FiniteDifferenceSimpsonsApproximatedSecondOrderBendingSensitivity(self, MemberNumber, scale = 1e-12):
        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateSimpsonsApproximatedNonLinearStrainEnergy()
        nonlinearstrain_energy = objective
        
        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale

        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateSimpsonsApproximatedNonLinearStrainEnergy()
        modified_nonlinearstrain_energy = objective

        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia -= scale

        #Calculate Sensitivity
        sensitivity = (modified_nonlinearstrain_energy - nonlinearstrain_energy) 
        
        return sensitivity
    
    def FiniteDifferenceSimpsonsApproximatedSecondOrderBendingSensitivityAllMember(self, scale = 1e-12):
        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateSimpsonsApproximatedNonLinearStrainEnergy()
        nonlinearstrain_energy = objective

        SensitivityList = []

        for Member in self.Members:
            Member.moment_of_inertia += scale

            StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
            objective = StrainEnergy1.CalculateSimpsonsApproximatedNonLinearStrainEnergy()
            modified_nonlinearstrain_energy = objective

            Member.moment_of_inertia -= scale
            #Calculate Sensitivity
            sensitivity = (modified_nonlinearstrain_energy - nonlinearstrain_energy) 
            SensitivityList.append(sensitivity)

        return SensitivityList            
    
    def FiniteDifferenceSecondOrderNonLinearBendingSensitivity(self,MemberNumber, scale=1e-12):
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
    
    def FiniteDifferenceSecondOrderNonLinearBendingSensitivityAllMember(self, scale = 1e-12):

        StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
        objective = StrainEnergy1.CalculateFiniteDifferenceNonLinearStrainEnergy()
        nonlinearstrain_energy = objective

        Senstivity_List = []
        for Member in self.Members:
            
            #small perturbation
            Member.moment_of_inertia += scale

            StrainEnergy1 = StrainEnergy(Points=self.Points, Members=self.Members, Loads=self.Loads)
            objective = StrainEnergy1.CalculateFiniteDifferenceNonLinearStrainEnergy()
            modified_nonlinearstrain_energy = objective

            Member.moment_of_inertia -= scale

            #Calculate Sensitivity
            sensitivity = (modified_nonlinearstrain_energy - nonlinearstrain_energy)
            Senstivity_List.append(sensitivity)
        
        return Senstivity_List
    


    
    def Members_All_Sensitivity(self, all = True):

        SimpsonsApproximated_List = []
        Approximated_list = []
        Linearized_list = []
        Linear_list = []
        Nonlinear_list = []
        for i in range(len(self.Members)):
            SimpsonsApproximated = self.FiniteDifferenceSimpsonsApproximatedSecondOrderBendingSensitivity(MemberNumber=i+1, scale=0.000000000001)
            Linearized = self.FiniteDifferenceSecondOrderLinearBendingSensitivity(MemberNumber = i+1, scale=0.000000000001)
            Linear = self.FiniteDifferenceFirstOrderLinearBendingSensitivity(MemberNumber = i+1, scale=0.000000000001)
            Nonlinear = self.FiniteDifferenceSecondOrderNonLinearBendingSensitivity(MemberNumber = i+1, scale=0.000000000001)

            SimpsonsApproximated_List.append(SimpsonsApproximated)
            Linearized_list.append(Linearized)
            Linear_list.append(Linear)
            Nonlinear_list.append(Nonlinear)
        
        print(Linearized_list)
        print(Linear_list)
        print(Nonlinear_list)
        print(SimpsonsApproximated_List)
        
        # Calculate relative errors for each element
        print("Relative error (w.r.t Nonlinear) for each element:")
        methods = [
            ("Linearized", Linearized_list),
            ("Linear", Linear_list),
            ("Simpsons Approxmated", SimpsonsApproximated_List)
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

        SimpsonsApproximated_norm = normalize_list(SimpsonsApproximated_List)
        Linearized_norm = normalize_list(Linearized_list)
        Linear_norm = normalize_list(Linear_list)
        Nonlinear_norm = normalize_list(Nonlinear_list)

        print("\nNormalized lists (by first value of each):")
        print("SimpsonsApproximated_norm", SimpsonsApproximated_norm)
        print("Linearized_norm:", Linearized_norm)
        print("Linear_norm:", Linear_norm)
        print("Nonlinear_norm:", Nonlinear_norm)

        # Calculate relative errors for normalized lists
        print("\nRelative error (w.r.t normalized Nonlinear) for each element:")
        norm_methods = [
            ("SimpsonsApproximated_norm", SimpsonsApproximated_norm),
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

class ComparisionSensitivity(Model):

    def __init__(self, **kwargs):
        self.Points = kwargs.get("Points", None)
        self.Members = kwargs.get("Members", None)
        self.Loads = kwargs.get("Loads", None)
        self.FiniteDifferenceSensitivity_object = FiniteDifferenceSensitivity(Points = self.Points, Members = self.Members, Loads = self.Loads)
        self.SecondOrderSensitivity_object = SecondOrderSensitivity(Points = self.Points, Members = self.Members, Loads = self.Loads)
    
    def AnlSimpson_Vs_FDNonLinear(self):

        SimpsonsApproximated_List = self.SecondOrderSensitivity_object.CalculateSimpsonsSensitivityAllMembers()
        Nonlinear_List = self.FiniteDifferenceSensitivity_object.FiniteDifferenceSecondOrderNonLinearBendingSensitivityAllMember()
        print("Simpsons", SimpsonsApproximated_List)
        print("NonLinear", Nonlinear_List)
        # Calculate relative errors for each element
        print("Relative error (w.r.t Nonlinear) for each element:")
        methods = [
            ("Simpsons Approxmated", SimpsonsApproximated_List)
        ]
        n = len(self.Members)
        error_sums = {name: 0.0 for name, _ in methods}
        for i in range(n):
            print(f"Element {i+1}:")
            nonlinear = Nonlinear_List[i]
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
    
    def AnlSimpson_Vs_FDSimpson(self):
        SimpsonsApproximated_List = self.SecondOrderSensitivity_object.CalculateSimpsonsSensitivityAllMembers()
        Numerical_SimpsonsApproximated = self.FiniteDifferenceSensitivity_object.FiniteDifferenceSimpsonsApproximatedSecondOrderBendingSensitivityAllMember()
        print("Analytical Simpsons", SimpsonsApproximated_List)
        print("Numerical Simpsons", Numerical_SimpsonsApproximated)
        # Calculate relative errors for each element
        print("Relative error (w.r.t Nonlinear) for each element:")
        methods = [
            ("Simpsons Approxmated", SimpsonsApproximated_List)
        ]
        n = len(self.Members)
        error_sums = {name: 0.0 for name, _ in methods}
        for i in range(n):
            print(f"Element {i+1}:")
            nonlinear = Numerical_SimpsonsApproximated[i]
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




