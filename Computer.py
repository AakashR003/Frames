

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, cg
from scipy.linalg import eig
from sksparse.cholmod import cholesky

class Computer():
    """
    This class is used for combining common computers on different class into gloabl computer
    """

    def StiffnessMatrixAssembler(UnConstrainedDoF,Members,StiffnessMatrixType, NormalForce = None):
        
        unconstrained_dofs = UnConstrainedDoF
        num_dofs = len(unconstrained_dofs)
        NoMembers = len(Members)
        
        # Create a DoF mapping for quick lookup
        dof_index = {dof: i for i, dof in enumerate(unconstrained_dofs)}
        
        # Initialize stiffness matrix as a NumPy array
        C1 = np.zeros((num_dofs, num_dofs))
        
        # Precompute Second Order Global Stiffness Matrices for all members
        if(NormalForce == None):        
            member_matrices = [
            np.array(getattr(Members[mn],StiffnessMatrixType)())
            for mn in range(NoMembers)
            ]
        else:
            member_matrices = [
            np.array(getattr(Members[mn],StiffnessMatrixType)(NormalForce[mn]))
            for mn in range(NoMembers)
            ]

        # Loop efficiently over members and DoFs
        for mn in range(NoMembers):
            member = Members[mn]
            dof_numbers = member.DoFNumber()
            K_local = member_matrices[mn]

            for mc in range(6):
                if dof_numbers[mc] in dof_index:
                    row = dof_index[dof_numbers[mc]]
                    for mr in range(6):
                        if dof_numbers[mr] in dof_index:
                            col = dof_index[dof_numbers[mr]]
                            C1[row, col] += K_local[mc, mr]

        return C1
   
    def GlobalStifnessMatrixA21():
        return None
    
    def DirectInverseDisplacementSolver(StiffnessMatrix, ForceVector):
        
        Displacement = np.dot((np.linalg.inv(np.array(StiffnessMatrix))),ForceVector)

        return Displacement
    
    def CholeskyDisplacementSolver(StiffnessMatrix, ForceVector):

        K = sp.csc_matrix(StiffnessMatrix)
        # Perform Cholesky factorization
        factor = cholesky(K)
        # Solving for displacement without computing the inverse explicitly
        Displacement = factor.solve_A(ForceVector)

        return Displacement
    
    def ConjugateGradientDisplacementSolver():
        return None

    def SupportForceVector():
        return None

    def ModelDisplacementList_To_Dict(Displacement,UnConstrainedDoF,TotalDoF):

        DisplacementDict={}
        for i in range(len(TotalDoF())):
            if(i<(len(UnConstrainedDoF()))):
                DisplacementDict[str(TotalDoF()[i])] = Displacement[i]
            else:
                DisplacementDict[str(TotalDoF()[i])]=0
        return DisplacementDict

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
        
        if "global" in StiffnessMatrixType.lower():
            raise ValueError("Conversion to global is not allowed in this Function.")

        MemberNo = int(MemberNumber)
        MemberForce = np.dot(
                    np.dot(
                    getattr(Members[MemberNo-1],StiffnessMatrixType)(NormalForce),Members[MemberNo-1].Transformation_Matrix()),
                    MemberDisplacement)
        FixedendForce = [0, 0, 0, 0, 0, 0]
        for a in range(len(Loads)):
            if(int(Loads[a].AssignedTo.split()[1]) == MemberNo):
                FixedendForcei = list(Loads[a].EquivalentLoad().values())[:-1]
                FixedendForcei= [x[0] for x in FixedendForcei]
                FixedendForce = [x + y for x, y in zip(FixedendForce, FixedendForcei)]
        MemberForce = np.round(MemberForce - FixedendForce,2)

        return MemberForce

    def ForceLocal_To_ForceGlobal(StiffnessMatrixType, MemberNumber, Members, MemberDisplacement, Loads, NormalForce = None):
        return None
    
    def Interpolate_Displacements(MemberDisplacment, length, n_points, ScaleFcator = 1 ):
        """
        Compute displacements at `n_points` along a beam element using shape functions.

        Parameters:
            nodal_values (list): List of nodal values [v_i, θ_i, v_j, θ_j].
            length (float): Length of the beam element (must be > 0).
            n_points (int): Number of points to interpolate (including endpoints).

        Returns:
            tuple: (x_values, displacements)
                x_values (list): Positions along the beam from 0 to `length`.
                displacements (list): Interpolated displacements at each position.
        """
        MemberDisplacment = np.array(MemberDisplacment) * ScaleFcator
        u_i, v_i, theta_i, u_j, v_j, theta_j = MemberDisplacment
           
        # Generate x values from 0 to length
        if n_points <= 1:
            x_values = [0.0]
        else:
            x_values = []
            x_valuesOutput = []
            for i in range(n_points):
                x = i*length/(n_points - 1) 
                x_values.append(x)
                x_valuesOutput.append( x + (u_i * (1-x/length)) + (u_j * x/length) )
        
        displacements = []
        for x in x_values:
            L = length
            # Compute generalized shape functions for any beam length L

            N1 = 1 - 3 * (x**2/L**2) + 2 * (x**3/L**3)
            N2 = x - 2 * (x**2/L) +(x**3/L**2)
            N3 = 3 * (x**2/L**2) - 2 * (x**3/L**3)
            N4 = -(x**2/L) + (x**3/L**2)
            
            # Calculate displacement
            v = N1 * v_i + N2 * theta_i + N3 * v_j + N4 * theta_j
            displacements.append(v)
        
        return x_valuesOutput, displacements

    def PlotStructuralElements(self, ax, Members, Points, sensitivities=None):
        """
        Helper function to plot structural elements (members, nodes, supports)
        ax: matplotlib axes object to plot on
        sensitivities: optional list of sensitivity values for color coding
        """
        # Plot members
        for i, member in enumerate(Members):
            start_node = member.Start_Node
            end_node = member.End_Node
            if sensitivities is not None:
                # Normalize sensitivities
                min_sensitivity = min(sensitivities)
                max_sensitivity = max(sensitivities)
                if max_sensitivity == min_sensitivity:
                    normalized_sensitivity = 0.5
                else:
                    normalized_sensitivity = (sensitivities[i] - min_sensitivity) / (max_sensitivity - min_sensitivity)
                color = plt.cm.OrRd(normalized_sensitivity)
                ax.plot([start_node.xcoordinate, end_node.xcoordinate], 
                        [start_node.ycoordinate, end_node.ycoordinate], 
                        color=color, linewidth=2)
            else:
                ax.plot([start_node.xcoordinate, end_node.xcoordinate], 
                       [start_node.ycoordinate, end_node.ycoordinate], 'b-')

        # Plot nodes and support conditions
        for i, node in enumerate(Points):
            # Plot nodes
            ax.plot(node.xcoordinate, node.ycoordinate, 'ro')
            
            # Add node numbers
            ax.text(node.xcoordinate, node.ycoordinate + 0.2, f"{i+1}", 
                   fontsize=12, ha='center', va='bottom', color='black')

            # Plot support conditions
            if node.support_condition == 'Fixed Support':
                ax.plot(node.xcoordinate, node.ycoordinate, 'ks', 
                       markersize=10, label="Fixed Support" if i == 0 else "")
            elif node.support_condition == 'Hinged Support':
                ax.plot(node.xcoordinate, node.ycoordinate, 'g^', 
                       markersize=10, label="Hinged Support" if i == 0 else "")
            elif node.support_condition == 'Roller in X-plane':
                ax.plot(node.xcoordinate, node.ycoordinate, 'bv', 
                       markersize=10, label="Roller in X-plane" if i == 0 else "")
            elif node.support_condition == 'Roller in Y-plane':
                ax.plot(node.xcoordinate, node.ycoordinate, 'r>', 
                       markersize=10, label="Roller in Y-plane" if i == 0 else "")
            elif node.support_condition == 'Hinge Joint':
                ax.plot(node.xcoordinate, node.ycoordinate, 'go', 
                       markerfacecolor='none', markersize=10, 
                       label="Hinged Support" if i == 0 else "")
      
    def GLobalStifnessMatrixCondensedA11_old(UnConstrainedDoF,Members,StiffnessMatrixType, NormalForce = None): #Stiffness matrix type - name of definition of Stiffness matrix in Member class
        NoMembers = len(Members)
        C1=[]
        for Mc in UnConstrainedDoF:
            R1=[]
            for Mr in UnConstrainedDoF:
                y=0
                for mn in range(0,NoMembers):
                    for mr in range(0,6):
                        if(Members[mn].DoFNumber()[mr]==Mr):
                            for mc in range(0,6):
                                if(Members[mn].DoFNumber()[mc]==Mc):
                                    if(NormalForce == None):
                                        x = getattr(Members[mn],StiffnessMatrixType)()[mc][mr]
                                    else:
                                        x = getattr(Members[mn],StiffnessMatrixType)(NormalForce[mn])[mc][mr]
                                    y=y+x
                R1.append(y)
            C1.append(R1)
        return C1

    def GlobalStiffnessMatrixold(TotalDoF,NoMembers,Members,StiffnessMatrixType):

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

