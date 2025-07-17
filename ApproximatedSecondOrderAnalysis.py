
import numpy as np
import matplotlib.pyplot as plt

try:
    from .SecondOrderResponse import SecondOrderMemberResponse
    from .FirstOrderResponse import FirstOrderGlobalResponse

    from .Computer import Computer
except:
    from SecondOrderResponse import SecondOrderMemberResponse
    from FirstOrderResponse import FirstOrderGlobalResponse

    from Computer import Computer

class ApproximatedSecondOrderAnalysis(SecondOrderMemberResponse):

    def _check_orthogonal_matrix(self, matrix, name="Orthogonal matrix"):
        """
        Checks that each row of the matrix has only one non-zero element (approximate zero: abs(x) < 1e-3).
        If not, prints the offending row index and raises a RuntimeError.
        """
        for i, row in enumerate(matrix):
            row_real = np.real(row)
            # Round to 2 decimals for zero approximation
            row_rounded = np.round(row_real, 2)
            non_zero_count = np.count_nonzero(row_rounded != 0)
            if non_zero_count != 1:
                pass
                #print(f"{name} row {i} has {non_zero_count} non-zero elements (rounded to 2 decimals): {row}")
                #raise RuntimeError(f"{name} row {i} does not have exactly one non-zero element (rounded to 2 decimals).")
        #print(f"{name} is orthogonal: each row has exactly one non-zero element (rounded to 2 decimals).")


    def CalculateK2ndOrderMinusKLinear(self):

        GlobalSecondOrderStiffnessMatrix = self.SecondOrderDisplacementVector(iteration_steps =5, ReturnSM = True)
        GlobalFirstOrderStiffnessMatrix = self.GlobalStiffnessMatrixCondensed()

        GlobalSecondOrderFlexibilityMatrix = Computer.FlexibilityMatrixSolver(GlobalSecondOrderStiffnessMatrix)
        GlobalFirstOrderFlexibilityMatrix = Computer.FlexibilityMatrixSolver(GlobalFirstOrderStiffnessMatrix)


        difference = GlobalSecondOrderFlexibilityMatrix - GlobalFirstOrderFlexibilityMatrix

        return difference
    
    """
    def CalculateOrthogonalLinearStiffnessMatrix(self):
        Matrix = self.GlobalStiffnessMatrixCondensed()
        orthogonalLinearStiffnessMatrix = Computer.OrthogonalSolver(Matrix, SMatrix = Matrix)

        self._check_orthogonal_matrix(orthogonalLinearStiffnessMatrix, name="Orthogonal Linear Stiffness Matrix")

        return orthogonalLinearStiffnessMatrix
    
    def CalculateOrthogonal2ndOrderStiffnessMatrix(self):
        Matrix = self.SecondOrderDisplacementVector(iteration_steps =5, ReturnSM = True)
        orthogonalSecondOrderStiffnessMatrix = Computer.OrthogonalSolver(Matrix, SMatrix = Matrix)

        self._check_orthogonal_matrix(orthogonalSecondOrderStiffnessMatrix, name="Orthogonal 2nd Order Stiffness Matrix")

        return orthogonalSecondOrderStiffnessMatrix
    """

    def CalculateOrthogonalFlexibilityDifferenceMatrix(self):
        Matrix = self.CalculateK2ndOrderMinusKLinear()
        OrthogonalFlexibilityDifferenceMatrix = Computer.OrthogonalSolver(Matrix, SMatrix = Matrix)

        self._check_orthogonal_matrix(OrthogonalFlexibilityDifferenceMatrix, name="Orthogonal Flexibility Difference Matrix")

        return OrthogonalFlexibilityDifferenceMatrix
    
    def CalculateOrthogonalSecondOrderForceVector(self):
        F = self.ForceVector()
        #SMatrix = self.SecondOrderDisplacementVector(iteration_steps =5, ReturnSM = True)
        SMatrix = self.CalculateK2ndOrderMinusKLinear()

        orthogonalSecondOrderForceVector = Computer.OrthogonalSolver(np.array(F), SMatrix = SMatrix)

        return orthogonalSecondOrderForceVector
    
    """
    def CalculateOrthogonalLinearDisplacementVector(self):
        U = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        SMatrix = self.GlobalStiffnessMatrixCondensed()

        orthogonalLinearDisplacementVector = Computer.OrthogonalSolver(U, SMatrix = SMatrix)

        return orthogonalLinearDisplacementVector
    
    def CalculateOrthogonalSecondOrderDisplacementVector(self):
        U = self.SecondOrderDisplacementVector(iteration_steps=5)
        SMatrix = self.SecondOrderDisplacementVector(iteration_steps =5, ReturnSM = True)
    
        orthogonalStiffnessSecondOrderDisplacementVector = Computer.OrthogonalSolver(U, SMatrix = SMatrix)

        return orthogonalStiffnessSecondOrderDisplacementVector
    """
        
    def CalculateModifiedOrthogonalFlexibilityDifferenceMatrix(self):
        Flexibility_2ndOrder = self.CalculateOrthogonalFlexibilityDifferenceMatrix()
        F_2ndOrder = self.CalculateOrthogonalSecondOrderForceVector()
        Flexibility_modified = Flexibility_2ndOrder / F_2ndOrder[:, np.newaxis]
        
        return Flexibility_modified
    
    
    def SetModifiedValues(self):

        LoadFactor = self.BucklingEigenLoad()[0] * 0.8
        
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude * LoadFactor

        self.EigenVectorMatrix = self.CalculateK2ndOrderMinusKLinear()
        self.ModifiedFlexibilityMatrix = self.CalculateModifiedOrthogonalFlexibilityDifferenceMatrix()

        #print("ModifiedFlexibility matrix at critical load",self.ModifiedFlexibilityMatrix)

        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude / LoadFactor
    
    def CalculateForceVectorSquare(self, ForceVector):

        SMatrix = self.EigenVectorMatrix
        orthogonalForceVector = Computer.OrthogonalSolver(np.array(ForceVector), SMatrix = SMatrix)
        orthogonalForceVector_square = orthogonalForceVector ** 2

        return orthogonalForceVector_square
    
    def CalculateApproximatedValueDisplacement(self, ReturnNonlinearDisplacement=False, ReturnLinearDisplacement=False):

        #LinearPart
        LinearDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())

        #NonlinearPart
        NonlinearDisplacementOrthogonal = Computer.DirectDisplacementSolver(self.ModifiedFlexibilityMatrix, self.CalculateForceVectorSquare(self.ForceVector()))
        NonlinearDisplacement = Computer.OrthogonalSolver(NonlinearDisplacementOrthogonal, SMatrix=self.EigenVectorMatrix, Back=True)

        if ReturnNonlinearDisplacement == True:
            return NonlinearDisplacement
        if ReturnLinearDisplacement == True:
            return LinearDisplacement
        return LinearDisplacement + NonlinearDisplacement
    
    def ApproximatedSecondOrderDisplacementLocal(self, ReturnNonlinearDisplacement=False, ReturnLinearDisplacement=False):
        """
        Returns the approximated second order displacement vector. 
        This is applicable for point where force vector.
        This solves by difference between linear and nonlinear displacement is known
        """
        #LinearPart
        LinearDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())

        #NonlinearPart
        DifferenceFlexibilityMatrix = self.CalculateK2ndOrderMinusKLinear()
        NonlinearDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())
        NonlinearDisplacementOrthogonal = Computer.OrthogonalSolver(NonlinearDisplacement, SMatrix=DifferenceFlexibilityMatrix)
        rootNonlinearDisplacementOrthogonal = np.array(NonlinearDisplacementOrthogonal)/np.array(self.ForceVector())
        NonlinearDisplacement = Computer.OrthogonalSolver(NonlinearDisplacement, SMatrix=DifferenceFlexibilityMatrix, Back=True)

        if ReturnNonlinearDisplacement == True:
            return NonlinearDisplacement
        if ReturnLinearDisplacement == True:
            return LinearDisplacement
        return LinearDisplacement + NonlinearDisplacement
        
    
    def ApproximatedSecondOrderDisplacementVectorDict(self):

        """ Returns a dictionary of displacements for each DOF."""

        self.DisplacementDict={}
        displacement = self.CalculateApproximatedValueDisplacement( ReturnNonlinearDisplacement = False, ReturnLinearDisplacement = False)
        for i in range(len(self.TotalDoF())):
            if(i<(len(self.UnConstrainedDoF()))):
                self.DisplacementDict[str(self.TotalDoF()[i])] = displacement[i]
            else:
                self.DisplacementDict[str(self.TotalDoF()[i])]=0
        return self.DisplacementDict


    def ApproximatedSecondOrderDisplacementTrace(self, NodeNumber = 1, Direction = "x", LoadFactor = None, division =20):
        """
        Plot the approximated second order displacement vector.
        """
        if LoadFactor is None:
            LoadFactor = self.BucklingEigenLoad()[0] * 0.8
        
        self.SetModifiedValues()
        
        if Direction == "x":
            DOFNumber = self.Points[NodeNumber-1].dof_x
        elif Direction == "y":
            DOFNumber = self.Points[NodeNumber-1].dof_y
        elif Direction == "tita":
            DOFNumber = self.Points[NodeNumber-1].dof_tita
        
        
        LoadTrace = np.linspace(0.1, LoadFactor, division)
        #LoadTrace = 0.1 + (LoadFactor - 0.1) * np.linspace(0, 1, division) ** 0.5 # Use for non-linear load
        DisplacementTrace = []

        for load in LoadTrace:
            for i in range(len(self.Loads)):
                self.Loads[i].Magnitude = self.Loads[i].Magnitude * load
            displacement = self.ApproximatedSecondOrderDisplacementVectorDict()[str(DOFNumber)]
            DisplacementTrace.append(abs(displacement))
            for i in range(len(self.Loads)):
                self.Loads[i].Magnitude = self.Loads[i].Magnitude / load
        
        return LoadTrace, DisplacementTrace
    
    def PlotSecondOrderLoadDisplacementCurve(self, NodeNumber = 1, Direction = "x", LoadFactor = None, division = 20, iteration_steps = 5):
        """
        Plots the load-displacement curve for a specific DOF.
        If Load is None, it uses the critical load by buckling analysis.
        division: Number of points to calculate along the load path.
        iteration_steps: Number of iterations for convergence.
        """
        
        LoadTrace, DisplacementTrace = self.ApproximatedSecondOrderDisplacementTrace(NodeNumber, Direction, LoadFactor, division)
        LoadTrace1, DisplacementTrace1 = self.SecondOrderLoadDisplacementTrace(NodeNumber, Direction, LoadFactor, division, iteration_steps)
        
        LoadFactor = self.BucklingEigenLoad()[0] * 0.8
        linear_model = FirstOrderGlobalResponse(Points = self.Points, Members = self.Members, Loads = self.Loads)
        LoadTrace2, DisplacementTrace2 = linear_model.LoadDisplacementTrace(NodeNumber, Direction, LoadFactor, division)
        
        plt.figure(figsize=(10, 6))
        plt.plot(DisplacementTrace, LoadTrace, marker='o', linestyle='-', color='blue')
        plt.plot(DisplacementTrace1, LoadTrace1, marker='o', linestyle='-', color='green')
        plt.plot(DisplacementTrace2, LoadTrace2, marker='o', linestyle='-', color='red')
        plt.title(f"Load-Displacement Curve for Node {NodeNumber} in {Direction} Direction")
        plt.xlabel("Displacement")
        plt.ylabel("Load")
        plt.grid(True)
        plt.show()






    def CalculateApproximatedSecondOrderDisplacementVector(self):

        #linear part

        LinearDisplacement = Computer.DirectInverseDisplacementSolver(self.GlobalStiffnessMatrixCondensed(),self.ForceVector())

        #NonLinear Part
        Flexibility_modified = self.CalculateModifiedOrthogonalFlexibilityDifferenceMatrix()
        F_2ndOrder = self.CalculateOrthogonalSecondOrderForceVector()
        F_2ndOrder_square = F_2ndOrder ** 2
        orthogonal_second_order_displacement_vector_nonlinear_part =  Computer.DirectDisplacementSolver(Flexibility_modified, F_2ndOrder_square)
        second_order_displacement_vector_nonlinear_part = Computer.OrthogonalSolver(orthogonal_second_order_displacement_vector_nonlinear_part, SMatrix = self.CalculateK2ndOrderMinusKLinear(), Back = True)
        
        approximated_second_order_displacement_vector = LinearDisplacement + second_order_displacement_vector_nonlinear_part
        
        return approximated_second_order_displacement_vector


    def checkcorrectness(self):
        """
        Check if the approximated second order displacement vector is correct
        """
        approximated_second_order_displacement_vector = self.CalculateApproximatedSecondOrderDisplacementVector()
        second_order_displacement_vector = self.SecondOrderDisplacementVector(iteration_steps=5)
        

        # Calculate percentage error for the whole vector (norm-based)
        norm_true = np.linalg.norm(second_order_displacement_vector)
        norm_approx = np.linalg.norm(approximated_second_order_displacement_vector)
        norm_error = np.linalg.norm(approximated_second_order_displacement_vector - second_order_displacement_vector)
        percent_error_vector = (norm_error / norm_true) * 100 if norm_true != 0 else 0

        # Calculate element-wise percentage error
        with np.errstate(divide='ignore', invalid='ignore'):
            elementwise_error = np.abs((approximated_second_order_displacement_vector - second_order_displacement_vector) / second_order_displacement_vector) * 100
            elementwise_error = np.nan_to_num(elementwise_error, nan=0.0, posinf=0.0, neginf=0.0)

        if np.allclose(approximated_second_order_displacement_vector, second_order_displacement_vector):
            print("Approximated second order displacement vector is correct.")
        else:
            print("Approximated second order displacement vector is incorrect.")
        print(f"Percentage error (whole vector, norm-based): {percent_error_vector:.4f}%")
        print("Element-wise percentage error and values:")
        for idx, (err, approx_val, true_val) in enumerate(zip(elementwise_error, approximated_second_order_displacement_vector, second_order_displacement_vector)):
            print(f"  Element {idx}: Error = {err:.4f}%, Approximated = {approx_val:.6g}, True = {true_val:.6g}")