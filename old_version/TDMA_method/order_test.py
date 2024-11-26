#!/usr/bin/env python3

import subprocess
import sys

# Define the range of powers for N (from 2^3 to 2^12) and boundary conditions
Npower = list(range(3, 13))
boundary_conditions = ["PP", "NN", "ND"]
results = {bc: [] for bc in boundary_conditions}  # Dictionary to store results for each BC

# Path to configuration file and executable
config_file = "PP.input"
executable = "./build/poisson_1d.exe"

# Iterate over each boundary condition
for BC in boundary_conditions:
    print(f"Processing for boundary condition: {BC}")
    
    # For each boundary condition, vary N as 2^3 to 2^12
    for power in Npower:
        # Calculate N as 2^power
        N = 2 ** power

        # Update the configuration file with the current N and boundary condition BC
        with open(config_file, "w") as f:
            f.write("&config\n")
            f.write(f"  N = {N}\n")
            f.write(f"  BC = '{BC}'\n")
            f.write("/\n")

        # Run the executable with the updated configuration file
        result = subprocess.run([executable, config_file], capture_output=True, text=True)
        
        # Check if the execution was successful
        if result.returncode != 0:
            print(f"Error running {executable} with N = {N} and BC = {BC}")
            results[BC].append(None)  # Append None if there's an error
            continue
        
        # Assume the error value is the only output or the last line of the output
        output = result.stdout.strip()
        try:
            # Convert the output to a float (assuming it's the error value)
            error = float(output.splitlines()[-1])
            results[BC].append(error)  # Store the error for this N and BC
        except ValueError:
            print(f"Unexpected output format for N = {N} and BC = {BC}: {output}")
            results[BC].append(None)  # Append None if parsing fails
            continue

# Write results to the output file with four columns: Npower, error (PP), error (NN), error (ND)
output_filename = "Npower_vs_error.txt"
with open(output_filename, "w") as f:
    f.write("# Npower error_PP error_NN error_ND\n")
    for i, power in enumerate(Npower):
        error_pp = results["PP"][i] if results["PP"][i] is not None else "NA"
        error_nn = results["NN"][i] if results["NN"][i] is not None else "NA"
        error_nd = results["ND"][i] if results["ND"][i] is not None else "NA"
        f.write(f"{power} {error_pp} {error_nn} {error_nd}\n")

print(f"Data written to {output_filename}")
