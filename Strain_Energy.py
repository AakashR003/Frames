
import numpy as np

try:
    from .FirstOrderResponse import FirstOrderGlobalResponse
    from .ApproximatedSecondOrderAnalysis import ApproximatedSecondOrderAnalysis

    from .Computer import Computer
except:
    from FirstOrderResponse import FirstOrderGlobalResponse
    from ApproximatedSecondOrderAnalysis import ApproximatedSecondOrderAnalysis

    from Computer import Computer

class StrainEnergy(ApproximatedSecondOrderAnalysis):
    """
    Class to calculate strain energy.
    """
    
    def CalculateLinearStrainEnergy(self):

        #LoadFactor = self.BucklingEigenLoad()[0] * 0.8
        LoadFactor = 1
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude * LoadFactor
        
        linear_model = FirstOrderGlobalResponse(Points = self.Points, Members = self.Members, Loads = self.Loads)
        displacement_vector = linear_model.DisplacementVector()
        ForceVector = self.ForceVector()
        
        strain_energy = 0.5 * np.dot(displacement_vector, np.transpose(ForceVector))

        strain_energy_orthogonal = 0.5 * np.dot(Computer.OrthogonalSolver(displacement_vector, SMatrix = self.GlobalStiffnessMatrixCondensed()),
                                        np.transpose(Computer.OrthogonalSolver(ForceVector, SMatrix = self.GlobalStiffnessMatrixCondensed())))

        print("Orthogonal",strain_energy_orthogonal)
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude / LoadFactor
        
        return strain_energy
    
    def CalculateLinearizedStrainEnergy(self):

        #LoadFactor = self.BucklingEigenLoad()[0] * 0.8
        LoadFactor = 1
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude * LoadFactor
        
        displacement_vector = self.SecondOrderDisplacementVector(iteration_steps=5)
        ForceVector = self.ForceVector()
        
        strain_energy = 0.5 * np.dot(displacement_vector, np.transpose(ForceVector))

        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude / LoadFactor
        
        return strain_energy
    
    def CalculateSimpsonsApproximatedNonLinearStrainEnergy(self):

        F2 = np.array(self.ForceVector())
        U2 = self.SecondOrderDisplacementVector(iteration_steps=5)
    
        
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude / 2
        
        U1 = self.SecondOrderDisplacementVector(iteration_steps=5)
        
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude * 2
        
        

        ComplimentaryStrainEnergy = np.dot(F2/6, np.transpose((4*U1 + U2)))
        rectarea = np.dot(F2, U2)

        StrainEnergy = rectarea - ComplimentaryStrainEnergy 

        return StrainEnergy
    
    def CalculateFiniteDifferenceNonLinearStrainEnergy(self):
        """
        Calculate the finite difference nonlinear strain energy.
        """

        #LoadFactor = self.BucklingEigenLoad()[0] * 0.8 #Uncomment this line if you want to use a buckling load factor
        LoadFactor = 1.0
        division = 10
        LoadTrace = np.linspace(0.1, LoadFactor, division)

        displacement_Previous = np.zeros(len(self.UnConstrainedDoF()))
        Force_Previous = np.zeros(len(self.UnConstrainedDoF()))
        
        strain_energy = 0
        for load in LoadTrace:
            for i in range(len(self.Loads)):
                self.Loads[i].Magnitude = self.Loads[i].Magnitude * load

            displacement = self.SecondOrderDisplacementVector(iteration_steps=5)
            Force_vector = np.array(self.ForceVector())
            strainenergy_at_infitesimalArea = 0.5* np.dot((displacement - displacement_Previous), (Force_vector + Force_Previous))
            strain_energy  = strain_energy  + strainenergy_at_infitesimalArea

            Force_Previous = Force_vector
            displacement_Previous = displacement

            for i in range(len(self.Loads)):
                self.Loads[i].Magnitude = self.Loads[i].Magnitude / load
        
        return strain_energy