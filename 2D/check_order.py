import subprocess
import os

# Define the cases with their corresponding integer values
cases = {
    "TAYLOR_GREEN": 0,
    "STEADY_CHANNEL": 1,
    "CAVITY": 2,
    "DEVELOPING_CHANNEL": 3
}

# Create the 'data' folder if it doesn't exist
output_dir = "data"
os.makedirs(output_dir, exist_ok=True)

# Exponents for grid sizes (N = 2^3 to 2^12)
exponents = list(range(3, 13))

# Function to run the executable and extract the last three lines for errors
def run_case_and_extract_errors(icase, nx, ny):
    # Create the input file content
    input_content = f"""
&case
    icase = {icase}
/
&gridsize
    nx = {nx}
    ny = {ny}
/
&geometry
    Lx = 1.0
    Ly = 1.0
/
"""
    # Write the input file
    with open("sample.input", "w") as input_file:
        input_file.write(input_content)

    # Run the executable and capture the output
    try:
        result = subprocess.run(
            ["./bin/PoissonSolver", "sample.input"],
            text=True,
            capture_output=True,
            check=True
        )
        # Extract the last three lines of the output for errors
        output_lines = result.stdout.strip().split("\n")[-3:]
        error_p = float(output_lines[0].split()[-1])
        error_dpd_x = float(output_lines[1].split()[-1])
        error_dpd_y = float(output_lines[2].split()[-1])
        return error_p, error_dpd_x, error_dpd_y
    except subprocess.CalledProcessError as e:
        print(f"Error running the executable for icase={icase}, nx={nx}: {e}")
        return None, None, None

# Loop through each case and grid size
for casename, icase in cases.items():
    # Prepare the output file
    output_file = os.path.join(output_dir, f"{casename}.txt")
    with open(output_file, "w") as f:
        # Write the header
        f.write("# N, Error in p, Error in dp/dx, Error in dp/dy\n")

        # Loop through each grid size (2^3 to 2^12)
        for exponent in exponents:
            nx = ny = 2 ** exponent
            error_p, error_dpd_x, error_dpd_y = run_case_and_extract_errors(icase, nx, ny)

            # Skip writing if there was an error in the execution
            if error_p is None:
                continue

            # Write the data line
            f.write(f"{exponent}, {error_p:.8e}, {error_dpd_x:.8e}, {error_dpd_y:.8e}\n")

print("Script completed successfully. Check the 'data' folder for the output files.")
