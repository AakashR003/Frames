


def print_class_Objects(objects, list_name=None):
    """
    Prints a list of objects in the same format as your input
    objects: List of Node/Member/Load objects
    list_name: Optional name for the list (defaults to object type name)
    """
    if not list_name:
        list_name = type(objects[0]).__name__ + "s"  # "Nodes", "Members", etc.
    
    print(f"{list_name} = [")
    for obj in objects:
        # Get all attributes of the object
        attrs = vars(obj)
        
        # Format the attributes
        attr_strs = []
        for name, value in attrs.items():
            if isinstance(value, str):
                attr_strs.append(f'{name}="{value}"')
            elif hasattr(value, 'node_number'):  # Handle Node references
                attr_strs.append(f'{name}=Points[{value.node_number-1}]')
            elif hasattr(value, 'Beam_Number'):  # Handle Member references
                attr_strs.append(f'{name}=Members[{value.Beam_Number-1}]')
            else:
                attr_strs.append(f'{name}={value}')
        
        # Print with proper indentation
        print(f"    {type(obj).__name__}({', '.join(attr_strs)}),")
    print("]")


def max_nested(lst):
    return max([max_nested(sub) if isinstance(sub, list) else sub for sub in lst])




"""
MemberFixedEndForce = [self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[0]],
                        self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[1]],
                        self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[2]],
                        self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[3]],
                        self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[4]],
                        self.ForceVectorDict[self.Members[self.MemberNo-1].DoFNumber()[5]]]
Use this if needed later
"""