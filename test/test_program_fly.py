import types

# Define the contents of the module
module_code = """
img = "images/SingleCell_RNA_PIPseq.png"
xx = "xx"
"""

# Create a new module object
program = types.ModuleType('')

# Execute the module code within the module object's namespace
exec(module_code, program.__dict__)

# Now you can access the variables defined in the module
print(program.img)
print(program.xx)
