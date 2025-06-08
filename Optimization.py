import matplotlib.pyplot as plt

try:
    from .Model import Model
    from .StructuralElements import Node, Member
    from .Computer import Computer
    from .Sensitivity import SecondOrderSensitivity, Senstivity
except:
    from Model import Model
    from StructuralElements import Node, Member
    from Computer import Computer
    from Sensitivity import SecondOrderSensitivity, Senstivity

class SecondOrderShapeOptimization(SecondOrderSensitivity):

    def SecondOderNodeOptimization(self, iterations = 8, stepsize = 0.02):

        for i in range(iterations):
            print("Iteration", i, "Objective",self.BucklingEigenLoad(Solver = "eigsh")[0])
            #sensitivity analysis
            self.GlobalSecondOrderNodeXSensitivity(EigenModeNo = 1)
            self.GlobalSecondOrderNodeYSensitivity(EigenModeNo = 1)
            #Size update
            for i in range (len(self.Points)):
                if self.Points[i].RestrainedDOF()[0] == False:
                    self.Points[i].xcoordinate = self.Points[i].xcoordinate + stepsize * self.Points[i].LBNodeXSensitivity
                if self.Points[i].RestrainedDOF()[1] == False:
                    self.Points[i].ycoordinate = self.Points[i].ycoordinate + stepsize * self.Points[i].LBNodeYSensitivity

        fig, ax = plt.subplots(figsize=(12, 8))
        computer_instance = Computer()
        computer_instance.PlotStructuralElements(ax, self.Members, self.Points, ShowNodeNumber=False)
                
        ax.set_title(f"Buckling Eigenmode Shape Optimized", fontsize=16)
        ax.axis('equal')
        plt.show()

        return None
        

class SecondOrderSizeOptimization(SecondOrderSensitivity):


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

    def ACSecondOrderBendingOptimization(self, iterations = 8, stepsize = 3, type = "Bending"): # aritifical constriained
        
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
    
class ShapeOptimization(Senstivity):
    
    def NodeOptimization(self, iterations = 8, stepsize = 0.02):

        for i in range(iterations):
            print("Iteration", i, "Objective")
            #sensitivity analysis
            self.GlobalNodeXSensitivity()
            self.GlobalNodeYSensitivity()
            #Size update
            for i in range (len(self.Points)):
                if self.Points[i].RestrainedDOF()[0] == False:
                    self.Points[i].xcoordinate = self.Points[i].xcoordinate + stepsize * self.Points[i].NodeXSensitivity
                if self.Points[i].RestrainedDOF()[1] == False:
                    self.Points[i].ycoordinate = self.Points[i].ycoordinate + stepsize * self.Points[i].NodeYSensitivity

        fig, ax = plt.subplots(figsize=(12, 8))
        computer_instance = Computer()
        computer_instance.PlotStructuralElements(ax, self.Members, self.Points, ShowNodeNumber=False)
                
        ax.set_title(f"Shape Optimized", fontsize=16)
        ax.axis('equal')
        plt.show()

        return None