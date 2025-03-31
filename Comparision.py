
import matplotlib.pyplot as plt
import numpy as np


try:
    from .main import Model
    from .StructuralElements import Node, Member
    from .Computer import Computer
except:
    from main import Model
    from StructuralElements import Node, Member
    from Computer import Computer


class Comparision():
    
    def __init__(self,**kwargs):
        
        self.MainModel = kwargs.get("MainModel", None)
        self.Model2 = kwargs.get("Model2", None)

    def PlotGlobalBMDComparison(self, scale_factor=1.0, show_structure=True):
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("Comparision Bending Moment Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.MainModel.Members, self.MainModel.Points)
        
        MemberForceLocalAll1 = self.MainModel.MemberForceLocal(1, All=True)
        MemberForceLocalAll2 = self.Model2.MemberForceLocal(1, All=True)

        
        # Determine global maximum absolute moment for scaling
        max_abs_moment1 = max(max(abs(moment) for moment in member_forces) 
                         for member_forces in MemberForceLocalAll1)
        max_abs_moment2 = max(max(abs(moment) for moment in member_forces) 
                         for member_forces in MemberForceLocalAll1)
        
        max_abs_moment = max(max_abs_moment1, max_abs_moment2)
        
        # Plot BMD for each member as simple lines
        for member_idx, member in enumerate(self.MainModel.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            moments = self.MainModel.MemberBMD(member_idx+1, MemberForceLocal=MemberForceLocalAll1[member_idx])
            positions = self.MainModel.amplist
            
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
            ax.plot(x_points, y_points, color='green', linewidth=1)
        

        for member_idx, member in enumerate(self.Model2.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            moments = self.Model2.MemberBMD(member_idx+1, MemberForceLocal=MemberForceLocalAll2[member_idx])
            positions = self.Model2.amplist
            
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

    def PlotGlobalSFDComparison(self, scale_factor=1.0, show_structure=True):
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("Comparision Bending Moment Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.MainModel.Members, self.MainModel.Points)
        
        MemberForceLocalAll1 = self.MainModel.MemberForceLocal(1, All=True)
        MemberForceLocalAll2 = self.Model2.MemberForceLocal(1, All=True)

        
        # Determine global maximum absolute moment for scaling
        max_abs_moment1 = max(max(abs(moment) for moment in member_forces) 
                         for member_forces in MemberForceLocalAll1)
        max_abs_moment2 = max(max(abs(moment) for moment in member_forces) 
                         for member_forces in MemberForceLocalAll1)
        
        max_abs_moment = max(max_abs_moment1, max_abs_moment2)
        
        # Plot BMD for each member as simple lines
        for member_idx, member in enumerate(self.MainModel.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            moments = self.MainModel.MemberBMD(member_idx+1, MemberForceLocal=MemberForceLocalAll1[member_idx])
            positions = self.MainModel.amplist
            
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
            ax.plot(x_points, y_points, color='green', linewidth=1)
        

        for member_idx, member in enumerate(self.Model2.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            moments = self.Model2.MemberBMD(member_idx+1, MemberForceLocal=MemberForceLocalAll2[member_idx])
            positions = self.Model2.amplist
            
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

    def PlotGlobalDeflectionComparison(self, scale_factor=1.0, show_structure=True):
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("Comparision Deflection Diagram")
        
        if show_structure:
            computer_instance = Computer()
            computer_instance.PlotStructuralElements(ax,self.MainModel.Members, self.MainModel.Points)

        DisplacementList1 = self.MainModel.DisplacementVector()
        DisplacementDict1 = Computer.ModelDisplacementList_To_Dict(DisplacementList1, self.MainModel.UnConstrainedDoF, self.MainModel.TotalDoF)

        iteration_steps = 5
        DisplacementList2 = self.Model2.DisplacementVector(iteration_steps)
        DisplacementDict2 = Computer.ModelDisplacementList_To_Dict(DisplacementList2, self.Model2.UnConstrainedDoF, self.Model2.TotalDoF)


        
        # Determine global maximum absolute deflection for scaling
        max_abs_deflection1 = max(DisplacementList1)
        max_abs_deflection2 = max(DisplacementList2)
        
        max_abs_deflection = max(max_abs_deflection1, max_abs_deflection2)
        
        # Plot Deflection for each member as simple lines
        for member_idx, member in enumerate(self.MainModel.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get Deflection values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            deflections = self.MainModel.MemberDeflection(member_idx+1, scale_factor, DisplacementDict = DisplacementDict1)
            positions = self.MainModel.DeflectionPosition
            
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
            ax.plot(x_points, y_points, color='green', linewidth=1)
        

        for member_idx, member in enumerate(self.Model2.Members):
            # Get member properties
            start = member.Start_Node
            end = member.End_Node
            L = member.length()
            
            # Get BMD values and positions
            #positions = self.MemberAmplitude(member_idx+1)
            deflections = self.Model2.MemberDeflection(member_idx+1, scale_factor, DisplacementDict = DisplacementDict2)
            positions = self.Model2.DeflectionPosition
            
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
            ax.plot(x_points, y_points, color='red', linewidth=1)

        ax.axis('equal')
        plt.show()

