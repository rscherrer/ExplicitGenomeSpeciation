# Write an R function that

rm(list = ls())

# Take a template file (or blank file if none)
# For each parameter requested
# Replace or add the line for this parameter with the requested values
# Save the file to the requested location

set_param_file(pars = list(scaleA = "5 0 0"), template = "parameters.txt",
               saveto = "param2.txt")
