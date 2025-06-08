import numpy as np
import matplotlib.pyplot as plt


try:
    from .Model import Model
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Functions import max_nested
    from .FirstOrderResponse import FirstOrderGlobalResponse
    from .SecondOrderResponse import SecondOrderGlobalResponse
except:
    from Model import Model
    from StructuralElements import Node, Member
    from Computer import Computer
    from Functions import max_nested
    from FirstOrderResponse import FirstOrderGlobalResponse
    from SecondOrderResponse import SecondOrderGlobalResponse

class Senstivity(FirstOrderGlobalResponse):

    def AxialMemberSensitivity(self,MemberNumber,scale):

        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].area += scale
                
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

        return sensitivity
    
    def BendingMemberSensitivity(self,MemberNumber,scale):
        
        UnMOdifiedSM = self.GlobalStiffnessMatrixCondensed()
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                self.Members[i].moment_of_inertia += scale
        
        ModifiedSM= self.GlobalStiffnessMatrixCondensed()
        d_AxialStiffness_ds = (np.array(ModifiedSM) - np.array(UnMOdifiedSM))
        sensitivity = np.dot(np.dot(np.transpose(self.DisplacementVector()),d_AxialStiffness_ds),self.DisplacementVector())

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

class SecondOrderSensitivity(SecondOrderGlobalResponse):
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
