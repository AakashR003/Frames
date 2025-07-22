
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
    
    def CalculateApproximatedNonlinearStrainEnergy(self):

        """ calculate strain enery at critical point """
        #LoadFactor = self.BucklingEigenLoad()[0] * 0.8
        LoadFactor = 1
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude * LoadFactor
        
        self.SetModifiedValues()
        Force_vector = self.ForceVector()
        
        #Linear part
        linear_model = FirstOrderGlobalResponse(Points = self.Points, Members = self.Members, Loads = self.Loads)
        LinearPart_displacement_vector = linear_model.DisplacementVector()
        Linear_partStrainEnergy = 0.5 * np.dot(LinearPart_displacement_vector, np.transpose(Force_vector))

        #NonLinear Part 
        differencedisplacement = self.SecondOderDisplacement - LinearPart_displacement_vector
        
        #nonlinearpart_strain_energy = 1/3*np.dot(differencedisplacement, np.transpose(Force_vector))
        nonlinearpart_strain_energy = 1/3*np.dot(self.CalculateApproximatedValueDisplacement(ReturnNonlinearDisplacement=True), np.transpose(Force_vector))
        #nonlinearpart_strain_energy = 2/3*np.dot(self.CalculateApproximatedValueDisplacement(ReturnNonlinearDisplacement=True), np.transpose(Force_vector))
        
        rect_area = np.dot(self.SecondOderDisplacement, np.transpose(Force_vector)) # we are like computing complementaory strain energy, then subtracting with reactangle area
        strain_energy1 = rect_area - nonlinearpart_strain_energy - Linear_partStrainEnergy #use 13 for nonlinear strain energy
        #print("Rect Area_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa", rect_area, strain_energy1)
        for i in range(len(self.Loads)):
            self.Loads[i].Magnitude = self.Loads[i].Magnitude / LoadFactor
        
        return strain_energy1 #Linear_partStrainEnergy + nonlinearpart_strain_energy
    
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