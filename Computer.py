

import numpy as np
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, cg, eigs
from scipy.linalg import eig
#from sksparse.cholmod import cholesky # Uncomment if using sksparse for Cholesky factorization

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
    
    def FlexibilityMatrixSolver(StiffnessMatrix):
        """
        Computes the flexibility matrix as the inverse of the stiffness matrix.
        """
        # Ensure the stiffness matrix is a NumPy array
        stiffness_matrix = np.array(StiffnessMatrix)
        
        # Check if the matrix is square and invertible
        if stiffness_matrix.shape[0] != stiffness_matrix.shape[1]:
            raise ValueError("Stiffness matrix must be square.")
        
        # Compute the flexibility matrix
        flexibility_matrix = np.linalg.inv(stiffness_matrix)
        
        return flexibility_matrix
    
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
    
    def DirectDisplacementSolver(FlexibilityMatrix, ForceVector):
        """
        Computes the displacement vector by multiplying the flexibility matrix with the force vector.
        """
        # Ensure the flexibility matrix is a NumPy array
        flexibility_matrix = np.array(FlexibilityMatrix)
        
        # Check if the matrix is square and compatible with the force vector
        if flexibility_matrix.shape[0] != flexibility_matrix.shape[1]:
            raise ValueError("Flexibility matrix must be square.")
        
        # Compute the displacement vector
        Displacement = np.dot(flexibility_matrix, ForceVector)
        
        return Displacement

    def DeterminantSolver(StiffnessMatrix):

        """
        if CholeskyMethod == True:
            K = sp.csc_matrix(StiffnessMatrix)/Norm
            factor = cholesky(K)
            logdet = 2 * np.sum(np.log(factor.D()))
            Determinant = np.exp(logdet)

            return Determinant
        Determinant = np.linalg.det(StiffnessMatrix/Norm)

        
        Computes a numerically stable determinant of matrix K.
        Scales the matrix using median absolute value and uses SVD-based log-determinant.
        Returns the actual determinant and log-determinant.
        """
        K= np.array(StiffnessMatrix)
        precision = 5
        getcontext().prec = precision  # Set precision

        K = np.array(K, dtype=object)
        K = np.vectorize(lambda x: Decimal(str(x)))(K)  # Convert to Decimal

        n = K.shape[0]
        det = Decimal(1)

        # Basic LU decomposition with partial pivoting
        for i in range(n):
            pivot = i + np.argmax([abs(K[j, i]) for j in range(i, n)])
            if K[pivot, i] == 0:
                return Decimal(0)
            if pivot != i:
                K[[i, pivot]] = K[[pivot, i]]
                det *= -1  # row swap changes sign

            det *= K[i, i]

            for j in range(i + 1, n):
                factor = K[j, i] / K[i, i]
                K[j, i:] = [K[j, k] - factor * K[i, k] for k in range(i, n)]

        return +det  # unary plus rounds to current context precision
    
    def EigenSolver(Matrix, AdditionalMatrix = None, Solver = "eigs", k=10):
        if Solver == "eigs":
            x, EigenMode = eigs(
                                Matrix, 
                                M=AdditionalMatrix, 
                                k=k, 
                                which='SM', 
                                maxiter=10000,    # Increase iterations
                                tol=1e-6,         # Loosen tolerance
                                ncv=50            # More Lanczos vectors
                                )
        
        if Solver == "eigsh":
            x, EigenMode = eigsh(
                                    Matrix, 
                                    M=AdditionalMatrix, 
                                    k=k, 
                                    sigma=1e-6,       # Shift near zero (critical for stability)
                                    which='LM',       # Largest magnitude after shift-invert
                                    mode='buckling',  # For buckling problems (A x = λ M x)
                                    maxiter=10000,
                                    tol=1e-6,
                                    ncv=50
                                )
        return EigenMode.real
    
    def OrthogonalSolver (Matrix, SMatrix, Back = False):
        
        if Back == False:
            if Matrix.ndim == 2 and Matrix.shape[1] == 1:
                return np.dot(np.transpose(Computer.EigenSolver(SMatrix, k = len(SMatrix))), Matrix)
            
            elif Matrix.ndim == 1:
                return np.dot(np.transpose(Computer.EigenSolver(SMatrix, k = len(SMatrix))), np.transpose(Matrix))
            
            else:       
                return np.dot(np.dot(np.transpose(Computer.EigenSolver(SMatrix, k = len(SMatrix))), Matrix), Computer.EigenSolver(SMatrix, k = len(SMatrix)))
        
        else:
            print("Back is True")
            if Matrix.ndim == 2 and Matrix.shape[1] == 1:
                return np.dot(Computer.EigenSolver(SMatrix, k = len(SMatrix)), Matrix)
            
            elif Matrix.ndim == 1:
                return np.dot(Computer.EigenSolver(SMatrix, k = len(SMatrix)), np.transpose(Matrix))
            
            else:    # need to check this   
                return np.dot(Computer.EigenSolver(SMatrix, k = len(SMatrix)), np.dot(Matrix, Computer.EigenSolver(SMatrix, k = len(SMatrix))))
            
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
        MemberDisplacementLocal = np.dot((Members[MemberNo-1].Transformation_Matrix()), MemberDisplacement)
        MemberForce = np.dot(
                    getattr(Members[MemberNo-1],StiffnessMatrixType)(NormalForce),
                    MemberDisplacementLocal)
        FixedendForce = [0, 0, 0, 0, 0, 0]
        for a in range(len(Loads)):
            if(int(Loads[a].AssignedTo.split()[1]) == MemberNo):
                FixedendForcei = Loads[a].EquivalentLoad(ReturnLocal = True)
                FixedendForce = [x + y for x, y in zip(FixedendForce, FixedendForcei)]
        MemberForce = np.round(MemberForce - FixedendForce,2)

        return MemberForce

    def ForceLocal_To_ForceGlobal(StiffnessMatrixType, MemberNumber, Members, MemberDisplacement, Loads, NormalForce = None):
        return None
    
    def Linear_Interpolate_Displacements(MemberDisplacment, length, n_points, scale_factor = 1 ):
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
        MemberDisplacment = np.array(MemberDisplacment) * scale_factor
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

            N1 = (1-x/length)
            N3 = x/length
            N2 = 0
            N4 = 0
            
            # Calculate displacement
            v = N1 * v_i + N2 * theta_i + N3 * v_j + N4 * theta_j
            displacements.append(v)
        
        return x_valuesOutput, displacements
    
    def Qudaratic_Interpolate_Displacements(MemberDisplacment, length, n_points,scale_factor = 1 ):
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
        MemberDisplacment = np.array(MemberDisplacment) * scale_factor
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

    def PlotStructuralElements(self, ax, Members, Points, ShowNodeNumber = True, sensitivities=None):
        """
        Helper function to plot structural elements (members, nodes, supports)
        ax: matplotlib axes object to plot on
        sensitivities: optional list of sensitivity values for color coding
        """
        # Plot members
        if sensitivities is not None:
            min_sensitivity = min(sensitivities)
            max_sensitivity = max(sensitivities)
            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r, norm=plt.Normalize(vmin=min_sensitivity, vmax=max_sensitivity))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label('Sensitivity', rotation=270, labelpad=15)

            
        for i, member in enumerate(Members):
            start_node = member.Start_Node
            end_node = member.End_Node
            if sensitivities is not None:
                # Normalize sensitivities
                normalized_sensitivity = (sensitivities[i] - min_sensitivity) / (max_sensitivity - min_sensitivity)
                color = plt.cm.RdBu_r(normalized_sensitivity)
                ax.plot([start_node.xcoordinate, end_node.xcoordinate], 
                        [start_node.ycoordinate, end_node.ycoordinate], 
                        color=color, linewidth = 3)
            else:
                ax.plot([start_node.xcoordinate, end_node.xcoordinate], 
                       [start_node.ycoordinate, end_node.ycoordinate], 'b-',
                       linewidth = 2)

        # Plot nodes and support conditions
        for i, node in enumerate(Points):
            # Plot nodes
            ax.plot(node.xcoordinate, node.ycoordinate, 'o', color='violet',  markersize = 4)
            ax.set_facecolor('black')
            
            # Add node numbers
            if ShowNodeNumber == True:
                ax.text(node.xcoordinate, node.ycoordinate + 0.2, f"{i+1}", 
                   fontsize=12, ha='center', va='bottom', color='violet')

            # Plot support conditions
            if node.support_condition == 'Fixed Support':
                ax.plot(node.xcoordinate, node.ycoordinate, 'gs', 
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

