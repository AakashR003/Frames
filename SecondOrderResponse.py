import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
from scipy.sparse.linalg import eigsh
#import importlib
import math
import scipy.sparse as sp
from scipy.sparse.linalg import eigs, cg
from scipy.sparse import csc_matrix



try:
    from .Model import Model
    from .config import config
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Functions import max_nested
except:
    from Model import Model
    from config import config
    from StructuralElements import Node, Member
    from Computer import Computer
    from Functions import max_nested



class SecondOrderGlobalResponse(Model):

    def NormalForce(self):

        NoMem = len(self.Members)

        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(MemberForceLocal[0])

        return NorForList
    
    def SecondOrderGlobalStiffnessMatrix(self, NormalForceList):
        
        C1 = Computer.StiffnessMatrixAssembler(self.TotalDoF(), self.Members, "Second_Order_Global_Stiffness_Matrix_1", NormalForceList)
        return C1
    
    def SecondOrderGlobalStiffnessMatrixCondensed(self, NormalForceList):

        C1 = Computer.StiffnessMatrixAssembler(self.UnConstrainedDoF(), self.Members, "Second_Order_Global_Stiffness_Matrix_1", NormalForceList)
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
    
    def DisplacementVector(self, iteration_steps):

        NoMem = len(self.Members)

        #1st iteration
        FirstOderDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(FirstOderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        NorForList =[]
        for i in range(NoMem):
            MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
            MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("First_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads )
            NorForList.append(-MemberForceLocal[0])

        #2nd iteration
        for j in range(0,iteration_steps):

            SecondOrderDisplacement = Computer.DirectInverseDisplacementSolver(self.SecondOrderGlobalStiffnessMatrixCondensed(NorForList),self.ForceVector())
            DisplacementDict = Computer.ModelDisplacementList_To_Dict(SecondOrderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
            
            NorForList1 = NorForList
            NorForList=[]
            for i in range(NoMem):
                SecondOrderMemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
                SecondOrderMemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("Second_Order_Local_Stiffness_Matrix_1", i+1, self.Members, SecondOrderMemberDisplacement, self.Loads, NorForList1[i])
                NorForList.append(-SecondOrderMemberForceLocal[0])
        
        self.NormalForceList = NorForList
        print("2nd order displacement computed")
        
        return SecondOrderDisplacement
    
    def DisplacementVectorDict(self):
        self.DisplacementDict={}
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.UnConstrainedDoF()))):
                self.DisplacementDict[str(self.TotalDoF()[i])] = self.DisplacementVector(5)[i]
            else:
                self.DisplacementDict[str(self.TotalDoF()[i])]=0
        return self.DisplacementDict
    
    def SecondOrderSupportForcesVector(self):

        #self.DisplacementVector(5)
        SupportForces = np.dot(np.array(self.SecondOrderGlobalStiffnessMatrixCondensedA21(self.NormalForceList)),self.DisplacementVector(5))
        
        self.ForceVectorDict={}
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.ConstrainedDoF()))):
                self.ForceVectorDict[str(self.TotalDoF()[i])]=SupportForces[i]
            else:
                self.ForceVectorDict[str(self.TotalDoF()[i])] = 0
        
        return SupportForces
    
    def BucklingEigenLoad(self, Solver = False):

        gr_buck = self.UnConstrainedDoF()

        BGSMConden = Computer.StiffnessMatrixAssembler(gr_buck,self.Members,"Second_Order_Global_Reduction_Matrix_1", NormalForce = self.NormalForce())
        BGSMM_1st_Ord_condensed = Computer.StiffnessMatrixAssembler(gr_buck,self.Members,"First_Order_Global_Stiffness_Matrix_1")
        
        CriticalLoad , EigenMode = eig(BGSMM_1st_Ord_condensed,BGSMConden)

        if Solver == "eigs":
            x, EigenMode = eigs(
                                BGSMM_1st_Ord_condensed, 
                                M=BGSMConden, 
                                k=10, 
                                which='SM', 
                                maxiter=10000,    # Increase iterations
                                tol=1e-6,         # Loosen tolerance
                                ncv=50            # More Lanczos vectors
                                )
        
        if Solver == "eigsh":
            x, EigenMode = eigsh(
                                    BGSMM_1st_Ord_condensed, 
                                    M=BGSMConden, 
                                    k=10, 
                                    sigma=1e-6,       # Shift near zero (critical for stability)
                                    which='LM',       # Largest magnitude after shift-invert
                                    mode='buckling',  # For buckling problems (A x = λ M x)
                                    maxiter=10000,
                                    tol=1e-6,
                                    ncv=50
                                )
        
        CriticalLoad = sorted(CriticalLoad, key=lambda x: abs(x))
        CriticalLoad = [round(float(x.real), 2) for x in CriticalLoad]
        print("Stability Eigen Calculated")

        return min(filter(math.isfinite, [abs(z.real) for z in CriticalLoad])), CriticalLoad, EigenMode
    
    def MemberEigenMode(self, MemberNumber, scale_factor = 1, EigenModeNo = 1, EigenVectorDict = None):
        
        FEDivision = config.get_FEDivision()
        MemberNo = int(MemberNumber)
        member = self.Members[MemberNo - 1]
        length = member.length()

        if EigenVectorDict == None:
            EigenVector = self.BucklingEigenLoad(EigenModeNo = EigenModeNo)[2][:, (EigenModeNo-1)]
            EigenVectorDict = Computer.ModelDisplacementList_To_Dict(EigenVector, self.UnConstrainedDoF, self.TotalDoF)

        EigenVectorGlobal = Computer.ModelDisplacement_To_MemberDisplacement(MemberNo, EigenVectorDict, self.Members)
        EigenVectorLocal = np.dot((self.Members[MemberNo-1].Transformation_Matrix()), EigenVectorGlobal)
        self.DeflectionPosition, BeamEigenVector = Computer.Linear_Interpolate_Displacements(EigenVectorLocal, length, FEDivision, scale_factor)
        
        return BeamEigenVector

    def PlotEigenMode(self, EigenModeNo = 1, scale_factor = 1, Solver ="eigsh", show_structure = True):
        fig, ax = plt.subplots(figsize=(12, 8))

        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = False)

        Eigen = self.BucklingEigenLoad(Solver = Solver)
        print("Eigen Load", Eigen[0], Eigen[1])
        EigenVector = Eigen[2][:,(EigenModeNo-1)]
        EigenVectorDict = Computer.ModelDisplacementList_To_Dict(EigenVector, self.UnConstrainedDoF, self.TotalDoF)
        
        # Determine global maximum absolute EigenDisp for scaling
        max_abs_Eigendeflection = max_nested(EigenVector)
        scale_factor = scale_factor / max_abs_Eigendeflection
        
        # Plot BMD for each member as lines
        for member_idx, member in enumerate(self.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get values and positions
            EigenModeDeflections = self.MemberEigenMode(member_idx+1, scale_factor, EigenVectorDict = EigenVectorDict)
            positions = self.DeflectionPosition
            
            # Calculate member orientation
            dx = end.xcoordinate - start.xcoordinate
            dy = end.ycoordinate - start.ycoordinate
            angle = np.arctan2(dy, dx)
            
            # Create perpendicular direction vector
            perp_dir = np.array([-np.sin(angle), np.cos(angle)])
            
            # Create points for BMD visualization
            x_points = []
            y_points = []
            for pos, deflection in zip(positions, EigenModeDeflections):
                # Calculate position along member
                x_pos = start.xcoordinate + (dx * pos/L)
                y_pos = start.ycoordinate + (dy * pos/L)
                
                # Offset by value in perpendicular direction
                x_points.append(x_pos + perp_dir[0] * deflection)
                y_points.append(y_pos + perp_dir[1] * deflection)
            
            # Plot as simple black line
            ax.plot(x_points, y_points, color='red', linewidth = 2)

        ax.set_title(f"Buckling Eigen Mode {EigenModeNo} - {Eigen[1][EigenModeNo-1]}")
        ax.axis('equal')
        plt.show()


class SecondOrderMemberResponse(SecondOrderGlobalResponse):
    
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
        MemberNo = int(MemberNumber)
        self.DisplacementVector(5)
        MemberDisplacement = [self.DisplacementVectorDict()[str(self.Members[MemberNo-1].DoFNumber()[0])],
                             self.DisplacementVectorDict()[str(self.Members[MemberNo-1].DoFNumber()[1])],
                             self.DisplacementVectorDict()[str(self.Members[MemberNo-1].DoFNumber()[2])],
                             self.DisplacementVectorDict()[str(self.Members[MemberNo-1].DoFNumber()[3])],
                             self.DisplacementVectorDict()[str(self.Members[MemberNo-1].DoFNumber()[4])],
                             self.DisplacementVectorDict()[str(self.Members[MemberNo-1].DoFNumber()[5])]]
        return MemberDisplacement
       
    def MemberForceLocal(self, MemberNumber, All = False):
        
        MemberNo = int(MemberNumber)
        SecondOrderDisplacement = self.DisplacementVector(5)
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(SecondOrderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(MemberNumber,DisplacementDict,self.Members)
        MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("Second_Order_Local_Stiffness_Matrix_1", MemberNo, self.Members, MemberDisplacement, self.Loads, self.NormalForceList[MemberNo-1] )

        if All == True:

            MemberForceLocalAll = []
            for i in range(self.NoMembers):
                MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(i+1,DisplacementDict,self.Members)
                MemberForceLocal = Computer.MemberDisplacement_To_ForceLocal("Second_Order_Local_Stiffness_Matrix_1", i+1, self.Members, MemberDisplacement, self.Loads, self.NormalForceList[MemberNo-1] )
                MemberForceLocalAll.append(MemberForceLocal)

            return MemberForceLocalAll

        """
        MemberNo = int(MemberNumber)
        SecondOrderDisplacement = self.DisplacementVector(5)
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(SecondOrderDisplacement,self.UnConstrainedDoF,self.TotalDoF)
        MemberDisplacement = Computer.ModelDisplacement_To_MemberDisplacement(MemberNumber,DisplacementDict,self.Members)

        MemberForce = np.dot(np.dot(self.Members[MemberNo-1].Second_Order_Local_Stiffness_Matrix_1(self.NormalForceList[MemberNo-1]),self.Members[MemberNo-1].Transformation_Matrix()),MemberDisplacement)


        FixedendForce = [0, 0, 0, 0, 0, 0]
        for a in range(len(self.Loads)):
            if(int(self.Loads[a].AssignedTo[-1]) == MemberNo):
                FixedendForcei = list(self.Loads[a].EquivalentLoad().values())[:-1]
                FixedendForcei= [x[0] for x in FixedendForcei]
                FixedendForce = [x + y for x, y in zip(FixedendForce, FixedendForcei)]
        MemberForce = np.round(MemberForce - FixedendForce,2)
        #"""

        return MemberForceLocal
    
    def MemberForceGlobal(self,MemberNumber):
        
        MemberNo = int(MemberNumber)
        MemberForce = self.MemberForceLocal(MemberNo)
        MemberForceGlobal = np.dot(np.transpose(self.Members[MemberNo-1].Transformation_Matrix()),MemberForce)

        return MemberForceGlobal
        
    def MemberBMD(self, MemberNumber, MemberForceLocal=None):

        FEDivision = config.get_FEDivision()
        MemberNo = int(MemberNumber)
        member = self.Members[MemberNo-1]
        length = member.length()
        alpha = member.alpha()

        if MemberForceLocal is None:
            local_forces = self.MemberForceLocal(MemberNo)
        else:
            local_forces = MemberForceLocal
        
        # Get local forces once and reuse
        fem1, fem2 = (local_forces[2], local_forces[5]) if alpha >= 0 else (local_forces[5], local_forces[2])

        
        # Initialize arrays with NumPy
        abcd1 = np.zeros(FEDivision)
        amp_values = np.linspace(0, length, FEDivision)
        self.amplist = amp_values.tolist()

        # Process loads using vectorization
        for load in self.Loads:
            if int(load.AssignedTo.split()[1]) == MemberNo:
                
                free_moment = np.array(load.EquivalentLoad()['FreeMoment'][:FEDivision])
                abcd1 += free_moment if alpha >= 0 else -free_moment

        # Vectorized calculations for abcd2 and abcd3
        abcd2 = (amp_values / length) * (-fem2 - fem1) + fem1
        abcd3 = -abcd1 + abcd2

        # Calculate shear forces using vectorized difference
        step = length / (FEDivision - 1)
        abcd4 = np.diff(abcd3) / step

        # Convert to lists if needed for compatibility
        self.MemberMoment = abcd3.tolist()
        self.MemberShear = abcd4.tolist()

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
        iteration_steps = 5

        if DisplacementDict == None:
            Displacement = self.DisplacementVector(iteration_steps)
            DisplacementDict = Computer.ModelDisplacementList_To_Dict(Displacement, self.UnConstrainedDoF, self.TotalDoF)

        MemberDisplacementGlobal = Computer.ModelDisplacement_To_MemberDisplacement(MemberNumber, DisplacementDict, self.Members)
        MemberDisplacementLocal = np.dot((self.Members[MemberNumber-1].Transformation_Matrix()),MemberDisplacementGlobal)
        self.DeflectionPosition, BeamDisplacement = Computer.Qudaratic_Interpolate_Displacements(MemberDisplacementLocal, length, FEDivision, ScaleFactor)
        

        return BeamDisplacement
    
    def PlotMemberBMD(self, MemberNumber):
        
        MemberNo = int(MemberNumber)
        x_max = int(self.Members[MemberNo-1].length())
        MemberBMD = self.MemberBMD(MemberNo)
        y_m_max = int(max(MemberBMD) * 2)
        y_m_min = int(min(MemberBMD) * 2)
        
        if y_m_max == 0:
            y_m_max = 5
        if y_m_min == 0:
            y_m_min = -5
        if y_m_max == y_m_min:
            y_m_max = abs(y_m_max)
            y_m_min = -abs(y_m_min)
        
        d = MemberBMD
        c = self.MemberAmplitude(MemberNo)
        g = [0, self.Members[MemberNo-1].length()]
        h = [0, 0]
        
        plt.figure(figsize=(8, 5))
        plt.plot(c, d, label="Bending Moment", color='red', linewidth=1.5)
        plt.plot(g, h, label="Baseline", color='black', linewidth=1.5, linestyle='dashed')
        
        plt.xlabel('Distance (Meter)')
        plt.ylabel('Bending Moment (kNm)')
        plt.xticks(range(0, x_max + 1, max(1, round(self.Members[MemberNo-1].length() / 10))))
        plt.yticks(range(y_m_min, y_m_max + 1, max(1, round((abs(y_m_max) + abs(y_m_min)) / 10))))
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.title(f'Second Order Moment Diagram for Member {MemberNo}')
        plt.show()

    def PlotGlobalBMD(self, scale_factor=0.5, show_structure=True):

        """
        Plots bending moment diagram with optional structure visualization
        scale_factor: Controls the visual scaling of BMD magnitudes
        show_structure: If True, shows the structural elements
        """

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("Second Order Bending Moment Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = False)
        
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
            ax.plot(x_points, y_points, color='red', linewidth=1)

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
        
        self.MemberAmplitude(self.MemberNo).pop()
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
        plt.title(f'Second Order Shear Force Diagram for Member {self.MemberNo}')
        plt.show()

    def PlotGlobalSFD(self, scale_factor=0.5, show_structure=True):

        """
        Plots Shear Force diagram with optional structure visualization
        scale_factor: Controls the visual scaling of BMD magnitudes
        show_structure: If True, shows the structural elements
        """

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("Second Order Shear Force Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = False)
        
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
            ax.plot(x_points, y_points, color='red', linewidth=1)

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
        plt.title(f'Second Order Deflection of Member {self.MemberNo}')
        plt.show()

    def PlotGlobalDeflection(self, scale_factor = 1, show_structure=True):

        """
        Plots Deflection with optional structure visualization
        scale_factor: Controls the visual scaling of BMD magnitudes
        show_structure: If True, shows the structural elements
        """

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("Second Order Deflection Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = False)
        
        iteration_steps = 5
        DisplacementList = self.DisplacementVector(5)
        DisplacementDict = Computer.ModelDisplacementList_To_Dict(DisplacementList, self.UnConstrainedDoF, self.TotalDoF)


        # Determine global maximum absolute moment for scaling
        max_abs_deflection = max(DisplacementList)
        scale_factor = scale_factor / max_abs_deflection
        
        # Plot BMD for each member as simple lines
        for member_idx, member in enumerate(self.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            deflections = self.MemberDeflection(member_idx+1, scale_factor, DisplacementDict = DisplacementDict)
            positions = self.DeflectionPosition
            
            # Calculate member orientation
            dx = end.xcoordinate - start.xcoordinate
            dy = end.ycoordinate - start.ycoordinate
            angle = np.arctan2(dy, dx)
            
            # Create perpendicular direction vector
            perp_dir = np.array([-np.sin(angle), np.cos(angle)])
            
            # Create points for BMD visualization
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
            ax.plot(x_points, y_points, color='red', linewidth = 2)

        ax.axis('equal')
        plt.show()
