
try:
    from .StructuralElements import Node, Member
    from .main import Model
except:
    from StructuralElements import Node, Member
    from main import Model

def create_framed_structure(x_bays, y_bays, x_spacing, y_spacing):
    Points = []
    Members = []
    node_number = 1
    beam_number = 1

    # Generate Points (Nodes)
    for i in range(y_bays + 1):  # Rows of nodes
        for j in range(x_bays + 1):  # Columns of nodes
            x = j * x_spacing
            y = i * y_spacing
            support_condition = "Rigid Joint"  # Default support condition
            if i == 0:  # Bottom row
                support_condition = "Fixed Support"
            elif i == y_bays:  # Top row
                support_condition = "Hinged Support"
            Points.append(Node(Node_Number=node_number, xcoordinate=x, ycoordinate=y, Support_Condition=support_condition))
            node_number += 1

    # Generate Members (Beams)
    for i in range(y_bays + 1):  # Horizontal members
        for j in range(x_bays):
            start_node = Points[i * (x_bays + 1) + j]
            end_node = Points[i * (x_bays + 1) + j + 1]
            Members.append(Member(Beam_Number=beam_number, Start_Node=start_node, End_Node=end_node, Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675))
            beam_number += 1

    for i in range(y_bays):  # Vertical members
        for j in range(x_bays + 1):
            start_node = Points[i * (x_bays + 1) + j]
            end_node = Points[(i + 1) * (x_bays + 1) + j]
            Members.append(Member(Beam_Number=beam_number, Start_Node=start_node, End_Node=end_node, Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675))
            beam_number += 1

    return Points, Members

# Example usage
x_bays = 3  # Number of bays in the x direction
y_bays = 2  # Number of bays in the y direction
x_spacing = 5  # Spacing in the x direction (meters)
y_spacing = 4  # Spacing in the y direction (meters)

Points, Members = create_framed_structure(x_bays, y_bays, x_spacing, y_spacing)

# You can now use Points and Members in your finite element model
Model1 = Model(Points=Points, Members=Members, Loads=[])
Model1.PlotGlobalModel()