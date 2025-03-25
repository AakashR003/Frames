
from main import Node, Member, NeumanBC

# Finite Element maker
def divide_into_finite_elements(nodes, members, loads, num_elements):
    # Create new nodes and members
    new_nodes = []
    new_members = []
    
    # Copy fixed nodes first (they need to keep their support conditions)
    node_counter = 1
    original_node_map = {}  # Maps original node numbers to their indices in new_nodes
    
    # First add all original nodes to new_nodes (we'll add intermediate nodes later)
    for node in nodes:
        new_node = Node(Node_Number=node_counter, 
                       xcoordinate=node.xcoordinate, 
                       ycoordinate=node.ycoordinate, 
                       Support_Condition=node.Support_Condition)
        new_nodes.append(new_node)
        original_node_map[node.Node_Number] = node_counter - 1  # store index (0-based)
        node_counter += 1
    
    # Process each member to divide it into elements
    for member in members:
        start_node = member.Start_Node
        end_node = member.End_Node
        
        # Calculate increments
        dx = (end_node.xcoordinate - start_node.xcoordinate) / num_elements
        dy = (end_node.ycoordinate - start_node.ycoordinate) / num_elements
        
        # Get the start node index in new_nodes
        if member.Start_Node.Node_Number in original_node_map:
            current_start_idx = original_node_map[member.Start_Node.Node_Number]
        else:
            # Shouldn't happen if input is correct
            raise ValueError("Start node not found in original nodes")
        
        # Create intermediate nodes and members
        for i in range(1, num_elements):
            # Create new node
            x = start_node.xcoordinate + i * dx
            y = start_node.ycoordinate + i * dy
            new_node = Node(Node_Number=node_counter, 
                           xcoordinate=x, 
                           ycoordinate=y, 
                           Support_Condition="None")
            new_nodes.append(new_node)
            
            # Create member from current_start to this new node
            new_member = Member(Beam_Number=len(new_members)+1,
                              Start_Node=new_nodes[current_start_idx],
                              End_Node=new_node,
                              Area=member.Area,
                              Youngs_Modulus=member.Youngs_Modulus,
                              Moment_of_Inertia=member.Moment_of_Inertia)
            new_members.append(new_member)
            
            # Update current_start for next iteration
            current_start_idx = len(new_nodes) - 1
            node_counter += 1
        
        # Add final member segment
        if member.End_Node.Node_Number in original_node_map:
            end_idx = original_node_map[member.End_Node.Node_Number]
        else:
            # Shouldn't happen if input is correct
            raise ValueError("End node not found in original nodes")
        
        new_member = Member(Beam_Number=len(new_members)+1,
                          Start_Node=new_nodes[current_start_idx],
                          End_Node=new_nodes[end_idx],
                          Area=member.Area,
                          Youngs_Modulus=member.Youngs_Modulus,
                          Moment_of_Inertia=member.Moment_of_Inertia)
        new_members.append(new_member)
    
    # Process loads - distribute point loads to nearest node
    new_loads = []
    for load in loads:
        if load.type == "PL" and load.AssignedTo.startswith("Member"):
            member_num = int(load.AssignedTo.split()[1])
            original_member = members[member_num - 1]
            
            # Calculate the position of the load in the original member
            total_length = ((original_member.End_Node.xcoordinate - original_member.Start_Node.xcoordinate)**2 +
                          (original_member.End_Node.ycoordinate - original_member.Start_Node.ycoordinate)**2)**0.5
            load_position = load.Distance1 / total_length  # normalized position [0,1]
            
            # Find which new element this would be in
            element_num = int(load_position * num_elements)
            element_num = min(max(element_num, 0), num_elements-1)  # clamp to valid range
            
            # Find the corresponding new member
            first_new_member_for_original = (member_num - 1) * num_elements
            target_member_idx = first_new_member_for_original + element_num
            target_member = new_members[target_member_idx]
            
            # Create new load at the nearest node (we'll put it at the end node of the element)
            new_load = NeumanBC(type="PL",
                              Magnitude=load.Magnitude,
                              Distance1=0,  # point load at node
                              AssignedTo=f"Node {target_member.End_Node.Node_Number}",
                              Members=None)  # Not assigned to member directly
            new_loads.append(new_load)
        else:
            # For other load types, just copy them (may need more sophisticated handling)
            new_loads.append(load)
    
    return new_nodes, new_members, new_loads