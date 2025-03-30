
from StructuralElements import Node, Member

# Finite Element maker
def divide_into_finite_elements(nodes, members, loads, num_elements):
    new_nodes = nodes.copy()
    new_members = []
    new_loads = []
    
    node_counter = max(node.node_number for node in nodes) + 1
    member_mapping = {}

    # Divide members and create new nodes
    for member_idx, member in enumerate(members):
        start = member.Start_Node
        end = member.End_Node
        total_length = member.length()
        
        dx = (end.xcoordinate - start.xcoordinate) / num_elements
        dy = (end.ycoordinate - start.ycoordinate) / num_elements
        
        member_mapping[member_idx] = []
        prev_node = start
        
        for i in range(1, num_elements):
            new_node = Node(
                Node_Number=node_counter,
                xcoordinate=start.xcoordinate + i*dx,
                ycoordinate=start.ycoordinate + i*dy,
                Support_Condition="Rigid Joint"
            )
            new_nodes.append(new_node)
            
            new_member = Member(
                Beam_Number=len(new_members)+1,
                Start_Node=prev_node,
                End_Node=new_node,
                Area=member.area,
                Youngs_Modulus=member.youngs_modulus,
                Moment_of_Inertia=member.moment_of_inertia
            )
            new_members.append(new_member)
            member_mapping[member_idx].append(new_member)
            
            prev_node = new_node
            node_counter += 1
        
        # Add final segment
        final_member = Member(
            Beam_Number=len(new_members)+1,
            Start_Node=prev_node,
            End_Node=end,
            Area=member.area,
            Youngs_Modulus=member.youngs_modulus,
            Moment_of_Inertia=member.moment_of_inertia
        )
        new_members.append(final_member)
        member_mapping[member_idx].append(final_member)

    # Process loads
    for load in loads:
        if not hasattr(load, 'AssignedTo') or not load.AssignedTo.startswith('Member'):
            new_loads.append(load)
            continue
            
        member_num = int(load.AssignedTo.split()[1]) - 1
        original_member = members[member_num]
        total_length = original_member.length()
        sub_members = member_mapping[member_num]
        sub_length = total_length / num_elements
        
        if load.type == "PL":
            # Point Load Handling
            abs_position = load.Distance1
            sub_idx = min(int(abs_position // sub_length), num_elements-1)
            local_position = abs_position - (sub_idx * sub_length)
            
            new_load = type(load)(
                type="PL",
                Magnitude=load.Magnitude,
                Distance1=local_position,
                AssignedTo=f"Member {sub_members[sub_idx].Beam_Number}",
                Members=new_members
            )
            new_loads.append(new_load)
            
        elif load.type == "UDL":
            # Uniform Load Handling
            udl_start = getattr(load, 'Distance1', 0)
            udl_end = getattr(load, 'Distance2', total_length)
            
            for i, sub_member in enumerate(sub_members):
                sub_start = i * sub_length
                sub_end = (i+1) * sub_length
                
                # Calculate overlap with UDL region
                overlap_start = max(udl_start, sub_start)
                overlap_end = min(udl_end, sub_end)
                
                if overlap_start >= overlap_end:
                    continue  # No overlap
                
                # Convert to local coordinates
                local_start = overlap_start - sub_start
                local_end = overlap_end - sub_start
                
                new_load = type(load)(
                    type="UDL",
                    Magnitude=load.Magnitude,  # Same intensity (force/unit-length)
                    Distance1=local_start,
                    Distance2=local_end,
                    AssignedTo=f"Member {sub_member.Beam_Number}",
                    Members=new_members
                )
                new_loads.append(new_load)

    return new_nodes, new_members, new_loads