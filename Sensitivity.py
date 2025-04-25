import numpy as np


try:
    from .Model import Model
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Functions import max_nested
    from .FirstOrderResponse import FirstOrderGlobalResponse
except:
    from Model import Model
    from StructuralElements import Node, Member
    from Computer import Computer
    from Functions import max_nested
    from FirstOrderResponse import FirstOrderGlobalResponse

class Senstivity(FirstOrderGlobalResponse):

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
