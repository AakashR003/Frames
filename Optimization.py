
try:
    from .Model import Model
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Sensitivity import SecondOrderSensitivity
except:
    from Model import Model
    from StructuralElements import Node, Member
    from Computer import Computer
    from Sensitivity import SecondOrderSensitivity

class SizeOptimization(SecondOrderSensitivity):

    def SecondOderBendingOptimization(self):
        pass

    def SecondOderAxialOptimization(self, iterations = 8, stepsize = 1, type = "Bending"):
        
        #No constraint added for volume
        for i in range(iterations):
            print("Iteration", i, "Objective",self.BucklingEigenLoad(Solver = "eigsh")[0])
            #sensitivity analysis
            sensitivity_values = self.GlobalSecondOrderBendingSensitivity(EigenModeNo = 1)
            #Size update
            for i in range (len(self.Members)):
                self.Members[i].moment_of_inertia = self.Members[i].moment_of_inertia - stepsize * sensitivity_values[i]

        
        pass

    def ACSecondOderAxialOptimization(self, iterations = 8, stepsize = 3, type = "Bending"):
        
        #Aritificially volume constraint added

        volume = 0
        for member in self.Members:
            volume += member.moment_of_inertia
        
        for i in range(iterations):
            print("Iteration", i, "Volume", volume, "Objective",self.BucklingEigenLoad(Solver = "eigsh")[0])
            #sensitivity analysis
            sensitivity_values = self.GlobalSecondOrderBendingSensitivity(EigenModeNo = 1)
            #Size update
            for i in range (len(self.Members)):
                self.Members[i].moment_of_inertia = self.Members[i].moment_of_inertia - stepsize * sensitivity_values[i]
            
            #volume constraint
            volume1 = volume
            volume = 0
            for member in self.Members:
                volume += member.moment_of_inertia
            
            change = (volume - volume1)/(len(self.Members)-1)
            for member in self.Members:
                member.moment_of_inertia = member.moment_of_inertia- change
            
        
        pass