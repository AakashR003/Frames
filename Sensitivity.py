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
    
    def NodeXSensitivity(self,NodeNumber,scale, EigenModeNo = 1):

        #unchanged system
        NoMem = len(self.Members)
        UnModifiedSM = self.GlobalStiffnessMatrixCondensed()

        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        UnModifiedGeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        #small perturbation
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].xcoordinate += scale
        
        #changed system
        ModifiedSM = self.GlobalStiffnessMatrixCondensed()
        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        ModifiedGeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        Eigen_mode = self.BucklingEigenLoad(Solver = "eigsh")[2][:,(EigenModeNo-1)]
        Eigen_load = self.BucklingEigenLoad(Solver = "eigsh")[1]

        d_KnodeX_ds = (np.array(ModifiedSM) - np.array(UnModifiedSM))
        d_KgNodeX_ds = (np.array(ModifiedGeoSM) - np.array(UnModifiedGeoSM))

        sensitivity = np.dot(np.dot(np.transpose(Eigen_mode),(d_KnodeX_ds - Eigen_load* d_KgNodeX_ds)), Eigen_mode)

        return sensitivity
    
    def NodeYSensitivity(self,NodeNumber,scale, EigenModeNo = 1):
        #unchanged system
        NoMem = len(self.Members)
        UnModifiedSM = self.GlobalStiffnessMatrixCondensed()

        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        UnModifiedGeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        #small perturbation
        for i in range(len(self.Points)):
            if self.Points[i].node_number == NodeNumber:
                self.Points[i].ycoordinate += scale
        
        #changed system
        ModifiedSM = self.GlobalStiffnessMatrixCondensed()
        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        ModifiedGeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        Eigen_mode = self.BucklingEigenLoad(Solver = "eigsh")[2][:,(EigenModeNo-1)]
        Eigen_load = self.BucklingEigenLoad(Solver = "eigsh")[1]

        d_KnodeY_ds = (np.array(ModifiedSM) - np.array(UnModifiedSM))
        d_KgNodeY_ds = (np.array(ModifiedGeoSM) - np.array(UnModifiedGeoSM))
        sensitivity = np.dot(np.dot(np.transpose(Eigen_mode),(d_KnodeY_ds - Eigen_load* d_KgNodeY_ds)), Eigen_mode)
        return sensitivity
    
    def BeamBendingSensitivity(self,MemberNumber,scale, EigenModeNo = 1):

        #unchanged system
        NoMem = len(self.Members)
        UnModifiedSM = self.GlobalStiffnessMatrixCondensed()

        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        UnModifiedGeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        #small perturbation
        for i in range(len(self.Members)):
            if i == MemberNumber-1:
                print(self.Members[i].moment_of_inertia)
                self.Members[i].moment_of_inertia += scale
                print(self.Members[i].moment_of_inertia)
        
        #changed system
        ModifiedSM = self.GlobalStiffnessMatrixCondensed()
        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])
        ModifiedGeoSM = self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList)

        Eigen_mode = self.BucklingEigenLoad(Solver = "eigsh")[2][:,(EigenModeNo-1)]
        Eigen_load = self.BucklingEigenLoad(Solver = "eigsh")[1][EigenModeNo-1]

        d_KMemberBen_ds = (np.array(ModifiedSM) - np.array(UnModifiedSM))
        d_KgMemberBen_ds = (np.array(ModifiedGeoSM) - np.array(UnModifiedGeoSM))

        sensitivity = np.dot(np.dot(np.transpose(Eigen_mode),(d_KMemberBen_ds - Eigen_load* d_KgMemberBen_ds)), Eigen_mode)
        
        self.Members[MemberNumber-1].LBSensitivityBend = sensitivity

        return sensitivity


    def GlobalSecondOrderShapeSensitivity(self, EigenModeNo = 1):
        for i in range(len(self.Members)):
            # Calculate the sensitivity for each member
            self.BeamBendingSensitivity(i+1, 1e-10, EigenModeNo = EigenModeNo)

        fig, ax = plt.subplots(figsize=(12, 8))
        computer_instance = Computer()
        computer_instance.PlotStructuralElements(ax, self.Members, self.Points, ShowNodeNumber=False)

        # Color mapping code - add this after the above call
        # Extract sensitivity values from members
        sensitivity_values = [member.LBSensitivityBend for member in self.Members]
        print(sensitivity_values)
        
        # Get min and max values
        min_sensitivity = min(sensitivity_values)
        max_sensitivity = max(sensitivity_values)

        # Create color mapping (Red to Blue)
        for i, member in enumerate(self.Members):
            start_node = member.Start_Node
            end_node = member.End_Node
            
            # Normalize the sensitivity value (0 to 1)
            if max_sensitivity == min_sensitivity:
                normalized_sensitivity = 0.5  # Handle case where all values are the same
            else:
                normalized_sensitivity = (sensitivity_values[i] - min_sensitivity) / (max_sensitivity - min_sensitivity)
            
            print(normalized_sensitivity)
            
            # Create RGB color (Red to Blue mapping)
            # Red = (1, 0, 0) for min values, Blue = (0, 0, 1) for max values
            red_component = 1 - normalized_sensitivity  # More red for lower values
            blue_component = normalized_sensitivity     # More blue for higher values
            color = (red_component, 0, blue_component)

            color = plt.cm.RdBu_r(normalized_sensitivity)
            
            # Plot the colored member
            ax.plot([start_node.xcoordinate, end_node.xcoordinate], 
                    [start_node.ycoordinate, end_node.ycoordinate], 
                    color=color, linewidth=3)  # Made linewidth thicker to show colors better
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r, norm=plt.Normalize(vmin=min_sensitivity, vmax=max_sensitivity))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('LB Sensitivity Bend', rotation=270, labelpad=15)

        ax.set_title(f"Buckling Eigenmode Shape Sensitivity", fontsize=16)
        ax.axis('equal')
        plt.show()
        return None
    
    def GlobalSizeSensitivity(self):
        return None
    
    def PlotSensitivity(self):
        return None
