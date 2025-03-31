import numpy as np
import matplotlib.pyplot as plt


try:
    from .main import Model
    from .StructuralElements import Node, Member
    from .config import config
    from .Computer import Computer
    from .Functions import max_nested
except:
    from main import Model
    from StructuralElements import Node, Member
    from config import config
    from Computer import Computer
    from Functions import max_nested

class FirstOrderGlobalResponse(Model):
    
    def DisplacementVector(self):
        self.Displacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        #DisplacementDict formation
        self.DisplacementDict={}
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.UnConstrainedDoF()))):
                self.DisplacementDict[str(self.TotalDoF()[i])] = self.Displacement[i]
            else:
                self.DisplacementDict[str(self.TotalDoF()[i])] = 0
        print("1st order displacement computed")
        return self.Displacement
    
    def DisplacementVectorDict(self):
        self.DisplacementVector()
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
   

class FirstOrderNodalResponse(FirstOrderGlobalResponse):
    
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


class FirstOrderMemberResponse(FirstOrderGlobalResponse):
    
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
        Displacement = self.DisplacementVector()
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(Displacement, self.UnConstrainedDoF, self.TotalDoF)
        MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(MemberNumber, DisplacementDict, self.Members)
        """
        MemberDisplacement = [self.DisplacementVectorDict()[str(self.Members[self.MemberNo-1].DoFNumber()[0])],
                             self.DisplacementVectorDict()[str(self.Members[self.MemberNo-1].DoFNumber()[1])],
                             self.DisplacementVectorDict()[str(self.Members[self.MemberNo-1].DoFNumber()[2])],
                             self.DisplacementVectorDict()[str(self.Members[self.MemberNo-1].DoFNumber()[3])],
                             self.DisplacementVectorDict()[str(self.Members[self.MemberNo-1].DoFNumber()[4])],
                             self.DisplacementVectorDict()[str(self.Members[self.MemberNo-1].DoFNumber()[5])]]
        """
        return MemberDisplacement
        
    def MemberForceLocal(self, MemberNumber, All = False):

        """ this function computes the local force in the member using the displacement vector.
        It uses the computer class to convert the model displacement to member displacement and 
        then computes the local force using the member displacement and member stiffness matrix.
        If All is True, it computes the local force for all members and returns a list of local forces"""

        self.MemberNo = int(MemberNumber)
        Displacement = self.DisplacementVector()
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(Displacement,self.UnConstrainedDoF,self.TotalDoF)
        MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(MemberNumber,DisplacementDict,self.Members)
        MemberForce = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", MemberNumber, self.Members, MemberDisplacement, self.Loads )

        if All == True:

            MemberForceLocalAll = []
            for i in range(self.NoMembers):
                MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
                MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
                MemberForceLocalAll.append(MemberForceLocal)

            return MemberForceLocalAll

        return MemberForce
    
    def MemberForceGlobal(self,MemberNumber):
        
        MemberForce = self.MemberForceLocal(MemberNumber)
        MemberForceGlobal = np.dot(np.transpose(self.Members[self.MemberNo-1].Transformation_Matrix()),MemberForce)

        return MemberForceGlobal

    def MemberBMD(self, MemberNumber, MemberForceLocal=None):

        """ This Function computes Bending moment diagram along the length of the beam by using MemberForceLocal.
        If MemberForceLocal is not provided, it computes MemberForceLocal using MemberForceLocal function.
        It divides the length of the beam into FEDivision parts and computes the BMD at each part. At each part, it computes
        the Fixed end moment and SS beam moment from Neuman class output and combines them to get the total moment distribution"""

        FEDivision = config.get_FEDivision()
        MemberNo = int(MemberNumber)
        member = self.Members[MemberNo - 1]
        alpha = member.alpha()
        length = member.length()

        if MemberForceLocal is None:
            local_forces = self.MemberForceLocal(MemberNo)
        else:
            local_forces = MemberForceLocal
        
        # Determine fem1 and fem2 based on alpha
        fem1, fem2 = (local_forces[2], local_forces[5]) if alpha >= 0 else (local_forces[5], local_forces[2])
        
        # Initialize abcd1 with NumPy for vector operations
        abcd1 = np.zeros(FEDivision)
        
        # Process loads using vectorization
        for load in self.Loads:
            if int(load.AssignedTo.split()[1]) == MemberNo:
                free_moment = np.array(load.EquivalentLoad()['FreeMoment'][:FEDivision])
                abcd1 += free_moment if alpha >= 0 else -free_moment
        
        # Vectorized calculation of abcd2 (linear moment component)
        amp_values = np.linspace(0, length, FEDivision)
        abcd2 = (amp_values / length) * (-fem2 - fem1) + fem1
        
        # Combine fixed end moments and linear component
        abcd3 = abcd1 + abcd2  # Total moment distribution
        
        # Calculate shear force using vectorized difference
        step = length / (FEDivision - 1)
        abcd4 = np.diff(abcd3) / step  # Shear distribution
        
        # Update instance variables (convert to lists if necessary)
        self.MemberMoment = abcd3.tolist()
        self.MemberShear = abcd4.tolist()
        self.amplist = amp_values.tolist()

        return self.MemberMoment
    
    def MemberSFD(self, MemberNumber, MemberForceLocal=None):
        self.MemberBMD(MemberNumber, MemberForceLocal)
        

        return self.MemberShear
    
    def MemberAmplitude(self, MemberNumber):
        
        return self.amplist
    
    def MemberNFD(self, MemberNumber):
        return None
    
    def MemberDeflection(self, MemberNumber, ScaleFactor = 1, DisplacementDict= None):
        
        FEDivision = config.get_FEDivision()
        MemberNo = int(MemberNumber)
        member = self.Members[MemberNo - 1]
        length = member.length()

        if DisplacementDict == None:
            Displacement = self.DisplacementVector()
            DisplacementDict = Computer.ModelDisplacementList_To_Dict(Displacement, self.UnConstrainedDoF, self.TotalDoF)

        MemberDisplacementGlobal = Computer.ModelDisplacement_To_MemberDisplacement(MemberNumber, DisplacementDict, self.Members)
        MemberDisplacementLocal = np.dot((self.Members[MemberNumber-1].Transformation_Matrix()),MemberDisplacementGlobal)
        self.DeflectionPosition, BeamDisplacement = Computer.Qudaratic_Interpolate_Displacements(MemberDisplacementLocal, length, FEDivision, ScaleFactor)
        

        return BeamDisplacement
    
    def PlotMemberBMD(self, MemberNumber):
        
        self.MemberNo = int(MemberNumber)
        x_max = int(self.Members[self.MemberNo-1].length())
        MemberBMD = self.MemberBMD(self.MemberNo)
        y_m_max = int(max(MemberBMD) * 2)
        y_m_min = int(min(MemberBMD) * 2)
        
        if y_m_max == 0:
            y_m_max = 5
        if y_m_min == 0:
            y_m_min = -5
        if y_m_max == y_m_min:
            y_m_max = abs(y_m_max)
            y_m_min = -abs(y_m_min)
        
        c = self.MemberAmplitude(self.MemberNo)
        d = MemberBMD
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
        plt.title(f'First Order Moment Diagram for Member {self.MemberNo}')
        plt.show()

    def PlotGlobalBMD(self, scale_factor=0.5, show_structure=True):

        """
        Plots bending moment diagram with optional structure visualization
        scale_factor: Controls the visual scaling of BMD magnitudes
        show_structure: If True, shows the structural elements
        """

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("First Order Bending Moment Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points)
        
        MemberForceLocalAll = self.MemberForceLocal(1, All=True)

        # Determine global maximum absolute moment for scaling
        max_abs_moment = max(max(abs(moment) for moment in member_forces) 
                         for member_forces in MemberForceLocalAll)
        
        # Plot BMD for each member as simple lines
        for member_idx, member in enumerate(self.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            moments = self.MemberBMD(member_idx+1, MemberForceLocal=MemberForceLocalAll[member_idx])
            positions = self.amplist
            
            # Calculate member orientation
            dx = end.xcoordinate - start.xcoordinate
            dy = end.ycoordinate - start.ycoordinate
            angle = np.arctan2(dy, dx)
            
            # Create perpendicular direction vector
            perp_dir = np.array([-np.sin(angle), np.cos(angle)])
            
            # Normalize moments and apply scaling
            scaled_moments = [m * scale_factor / max_abs_moment if max_abs_moment != 0 else 0 
                             for m in moments]
            
            # Create points for BMD visualization
            x_points = []
            y_points = []
            for pos, moment in zip(positions, scaled_moments):
                # Calculate position along member
                x_pos = start.xcoordinate + (dx * pos/L)
                y_pos = start.ycoordinate + (dy * pos/L)
                
                # Offset by moment value in perpendicular direction
                x_points.append(x_pos + perp_dir[0] * moment)
                y_points.append(y_pos + perp_dir[1] * moment)
            
            # Plot BMD as simple black line
            ax.plot(x_points, y_points, color='black', linewidth=1)

        ax.axis('equal')
        plt.show()

    def PlotMemberSFD(self, MemberNumber):
        
        self.MemberNo = int(MemberNumber)
        x_max = int(self.Members[self.MemberNo-1].length())
        MemberSFD = self.MemberSFD(self.MemberNo)
        y_m_max = int(max(MemberSFD) * 2)
        y_m_min = int(min(MemberSFD) * 2)
        
        if y_m_max == 0:
            y_m_max = 5
        if y_m_min == 0:
            y_m_min = -5
        if y_m_max == y_m_min:
            y_m_max = abs(y_m_max)
            y_m_min = -abs(y_m_min)
        
        c = self.MemberAmplitude(self.MemberNo)
        d = MemberSFD
        g = [0, self.Members[self.MemberNo-1].length()]
        h = [0, 0]
        
        plt.figure(figsize=(8, 5))
        plt.plot(c, d, label="Shear Force", color='red', linewidth=1.5)
        plt.plot(g, h, label="Baseline", color='black', linewidth=1.5, linestyle='dashed')
        
        plt.xlabel('Distance (Meter)')
        plt.ylabel('Shear Force (kN)')
        plt.xticks(range(0, x_max + 1, max(1, round(self.Members[self.MemberNo-1].length() / 10))))
        plt.yticks(range(y_m_min, y_m_max + 1, max(1, round((abs(y_m_max) + abs(y_m_min)) / 10))))
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.title(f'First Order Shear Force Diagram for Member {self.MemberNo}')
        plt.show()

    def PlotGlobalSFD(self, scale_factor=0.5, show_structure=True):

        """
        Plots Shear Force diagram with optional structure visualization
        scale_factor: Controls the visual scaling of BMD magnitudes
        show_structure: If True, shows the structural elements
        """

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("First Order Shear Force Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points)
        
        MemberForceLocalAll = self.MemberForceLocal(1, All=True)

        # Determine global maximum absolute shear for scaling
        max_abs_moment = max(max(abs(moment) for moment in member_forces) 
                         for member_forces in MemberForceLocalAll)
        
        # Plot SFD for each member as simple lines
        for member_idx, member in enumerate(self.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get SFD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            Shears = self.MemberSFD(member_idx+1, MemberForceLocal=MemberForceLocalAll[member_idx])
            positions = self.amplist
            
            # Calculate member orientation
            dx = end.xcoordinate - start.xcoordinate
            dy = end.ycoordinate - start.ycoordinate
            angle = np.arctan2(dy, dx)
            
            # Create perpendicular direction vector
            perp_dir = np.array([-np.sin(angle), np.cos(angle)])
            
            # Normalize moments and apply scaling
            scaled_moments = [m * scale_factor / max_abs_moment if max_abs_moment != 0 else 0 
                             for m in Shears]
            
            # Create points for SFD visualization
            x_points = []
            y_points = []
            for pos, moment in zip(positions, scaled_moments):
                # Calculate position along member
                x_pos = start.xcoordinate + (dx * pos/L)
                y_pos = start.ycoordinate + (dy * pos/L)
                
                # Offset by moment value in perpendicular direction
                x_points.append(x_pos + perp_dir[0] * moment)
                y_points.append(y_pos + perp_dir[1] * moment)
            
            # Plot SFD as simple black line
            ax.plot(x_points, y_points, color='black', linewidth=1)

        ax.axis('equal')
        plt.show()

    def PlotMemberDeflection(self, MemberNumber):
        
        self.MemberNo = int(MemberNumber)
        x_max = int(self.Members[self.MemberNo-1].length())
        MemberDeflection = self.MemberDeflection(self.MemberNo)
        y_m_max = int(max(MemberDeflection) * 2)
        y_m_min = int(min(MemberDeflection) * 2)
        
        if y_m_max == 0:
            y_m_max = 5
        if y_m_min == 0:
            y_m_min = -5
        if y_m_max == y_m_min:
            y_m_max = abs(y_m_max)
            y_m_min = -abs(y_m_min)
        
        c = self.DeflectionPosition
        d = MemberDeflection
        g = [0, self.Members[self.MemberNo-1].length()]
        h = [0, 0]
        
        plt.figure(figsize=(8, 5))
        plt.plot(c, d, label="Deflection", color='red', linewidth=1.5)
        plt.plot(g, h, label="Baseline", color='black', linewidth=1.5, linestyle='dashed')
        
        plt.xlabel('Distance (Meter)')
        plt.ylabel('Deflection (m)')
        plt.xticks(range(0, x_max + 1, max(1, round(self.Members[self.MemberNo-1].length() / 10))))
        plt.yticks(range(y_m_min, y_m_max + 1, max(1, round((abs(y_m_max) + abs(y_m_min)) / 10))))
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.title(f'First Order Deflection of Member {self.MemberNo}')
        plt.show()

    def PlotGlobalDeflection(self, scale_factor = 1, show_structure=True):

        """
        Plots Deflection with optional structure visualization
        scale_factor: Controls the visual scaling of BMD magnitudes
        show_structure: If True, shows the structural elements
        """

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("First Order Deflection Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = False)
        
        DisplacementList = self.DisplacementVector()
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(DisplacementList, self.UnConstrainedDoF, self.TotalDoF)


        # Determine global maximum absolute deflection for scaling
        max_abs_deflection = max(DisplacementList)
        scale_factor = scale_factor / max_abs_deflection
        
        # Plot Deflection for each member as simple lines
        for member_idx, member in enumerate(self.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get Deflection values and positions
            deflections = self.MemberDeflection(member_idx+1, scale_factor, DisplacementDict = DisplacementDict)
            positions = self.DeflectionPosition
            
            # Calculate member orientation
            dx = end.xcoordinate - start.xcoordinate
            dy = end.ycoordinate - start.ycoordinate
            angle = np.arctan2(dy, dx)
            
            # Create perpendicular direction vector
            perp_dir = np.array([-np.sin(angle), np.cos(angle)])
            
            # Create points for Deflection visualization
            x_points = []
            y_points = []
            for pos, deflection in zip(positions, deflections):
                # Calculate position along member
                x_pos = start.xcoordinate + (dx * pos/L)
                y_pos = start.ycoordinate + (dy * pos/L)
                
                # Offset by moment value in perpendicular direction
                x_points.append(x_pos + perp_dir[0] * deflection)
                y_points.append(y_pos + perp_dir[1] * deflection)
            
            # Plot BMD as simple black line
            ax.plot(x_points, y_points, color='black', linewidth=1)

        ax.axis('equal')
        plt.show()
