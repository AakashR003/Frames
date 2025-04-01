
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
    from .main import Model
    from .config import config
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Functions import max_nested
except:
    from main import Model
    from config import config
    from StructuralElements import Node, Member
    from Computer import Computer
    from Functions import max_nested


class DynamicGlobalResponse(Model):

    def EigenFrequency(self, EigenModeNo = False):

        dof = self.UnConstrainedDoF()

        MM_Conden = Computer.StiffnessMatrixAssembler(dof,self.Members,"Global_Mass_Matrix")
        _1st_OrdSM_condensed = Computer.StiffnessMatrixAssembler(dof,self.Members,"First_Order_Global_Stiffness_Matrix_1")
        
        EigenFreq , EigenMode = eig(_1st_OrdSM_condensed, MM_Conden)
        if EigenModeNo:
            x, EigenMode = eigsh(
                                    _1st_OrdSM_condensed, 
                                    M=MM_Conden, 
                                    k=10, 
                                    sigma=1e-6,       # Shift near zero (critical for stability)
                                    which='LM',       # Largest magnitude after shift-invert
                                    mode='buckling',  # For buckling problems (A x = λ M x)
                                    maxiter=10000,
                                    tol=1e-6,
                                    ncv=50
                                )
        
        EigenFreq = sorted(EigenFreq, key=lambda x: abs(x)) # answer will be in radians - converting to Hz 
        EigenFreq = [round(float(x.real)**0.5/(2*np.pi), 2) for x in EigenFreq]
        print("Dynamic Eigen Calculated in Hertz", EigenFreq)

        return min(filter(math.isfinite, [abs(z.real) for z in EigenFreq])), EigenFreq, EigenMode

    def MemberEigenMode(self, MemberNumber, scale_factor = 10000, EigenModeNo = 2, EigenVectorDict = None):
        
        FEDivision = config.get_FEDivision()
        MemberNo = int(MemberNumber)
        member = self.Members[MemberNo - 1]
        length = member.length()

        if EigenVectorDict == None:
            EigenVector = self.EigenFrequency(EigenModeNo = EigenModeNo)[2][:, (EigenModeNo-1)]
            EigenVectorDict = Computer.ModelDisplacementList_To_Dict(EigenVector, self.UnConstrainedDoF, self.TotalDoF)

        EigenVectorGlobal = Computer.ModelDisplacement_To_MemberDisplacement(MemberNo, EigenVectorDict, self.Members)
        EigenVectorLocal = np.dot((self.Members[MemberNo-1].Transformation_Matrix()), EigenVectorGlobal)
        self.DeflectionPosition, BeamEigenVector = Computer.Linear_Interpolate_Displacements(EigenVectorLocal, length, FEDivision, scale_factor)
        
        return BeamEigenVector


    def PlotDynamicEigenMode(self, EigenModeNo = 1, scale_factor = 1, show_structure = True):

        fig, ax = plt.subplots(figsize=(12, 8))
        
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.Members, self.Points, ShowNodeNumber = False)

        Eigen = self.EigenFrequency(EigenModeNo = EigenModeNo)
        print("Eigen Frequency", Eigen[0], Eigen[1])
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

        ax.set_title(f"Dynamic Eigen Mode {EigenModeNo} - {Eigen[1][EigenModeNo-1]}")
        ax.axis('equal')
        plt.show()

