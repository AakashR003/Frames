import matplotlib.pyplot as plt

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

class ShapeOptimization(SecondOrderSensitivity):

    def SecondOderNodeOptimization(self, iterations = 8, stepsize = 0.001):

        for i in range(iterations):
            print("Iteration", i, "Objective",self.BucklingEigenLoad(Solver = "eigsh")[0])
            #sensitivity analysis
            sensitivity_values_x = self.GlobalSecondOrderNodeXSensitivity(EigenModeNo = 1)
            sensitivity_values_y = self.GlobalSecondOrderNodeYSensitivity(EigenModeNo = 1)
            #Size update
            for i in range (len(self.Points)):
                self.Points[i].xcoordinate = self.Points[i].xcoordinate - stepsize * sensitivity_values_x[i]
                self.Points[i].ycoordinate = self.Points[i].ycoordinate - stepsize * sensitivity_values_y[i]

        fig, ax = plt.subplots(figsize=(12, 8))
        computer_instance = Computer()
        computer_instance.PlotStructuralElements(ax, self.Members, self.Points, ShowNodeNumber=False)
                
        ax.set_title(f"Buckling Eigenmode Shape Optimized", fontsize=16)
        ax.axis('equal')
        plt.show()

        return None
        

class SizeOptimization(SecondOrderSensitivity):


    def SecondOrderBendingOptimization(self, iterations = 8, stepsize = 1, type = "Bending"):
        
        #No constraint added for volume
        for i in range(iterations):
            print("Iteration", i, "Objective",self.BucklingEigenLoad(Solver = "eigsh")[0])
            #sensitivity analysis
            sensitivity_values = self.GlobalSecondOrderBendingSensitivity(EigenModeNo = 1)
            #Size update
            for i in range (len(self.Members)):
                self.Members[i].moment_of_inertia = self.Members[i].moment_of_inertia - stepsize * sensitivity_values[i]

        
        return None

    def ACSecondOrderBendingOptimization(self, iterations = 8, stepsize = 3, type = "Bending"):
        
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
            
        
        return None