# Create a variable named "name" with value equal to default, only if such a variable is not yet defined.

def declareDefault(name, default, glb):
    if name not in glb:
        glb[name]=default
