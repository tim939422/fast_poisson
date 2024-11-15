#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science", "nature"])

# Load data from Npower_vs_error.txt, skipping the header line with '#'
data = np.loadtxt("Npower_vs_error.txt", comments='#', dtype={'names': ('Npower', 'error_PP', 'error_NN', 'error_ND'), 'formats': (int, float, float, float)})

# Extract columns
Npower = data['Npower']
error_PP = data['error_PP']
error_NN = data['error_NN']
error_ND = data['error_ND']

# Create the figure and axis
fig, ax = plt.subplots()

# Plot each boundary condition error
ax.plot(Npower, error_PP, marker='o', label='PP')
ax.plot(Npower, error_NN, marker='s', label='NN')
ax.plot(Npower, error_ND, marker='^', label='ND')

# Set x-tick labels to display 2^N
ax.set_xticks(Npower)
ax.set_xticklabels([f"$2^{{{n}}}$" for n in Npower])

# Add a helper line for second-order convergence with a higher starting point and shorter length
# Using the first error value in error_PP as a reference, and starting it higher
higher_reference_error = error_PP[0] * 5  # Start at 5 times the initial error_PP[0] value
shortened_Npower = Npower[2:8]  # Use the first half of Npower for a shorter line
second_order_convergence = [higher_reference_error / (2**(2 * (n - Npower[0]))) for n in shortened_Npower]
ax.plot(shortened_Npower, second_order_convergence, 'k--', label=r'$\sim \Delta x^2$')

# Labeling the plot
ax.set_xlabel("Grid Size ($N$)")
ax.set_ylabel("Error")
ax.set_yscale("log")  # Use a logarithmic scale for error
ax.legend()

ax.set_yscale("log")
ax.grid(which="major", axis="x", linestyle="--", linewidth=0.5)  # Major grid lines on x-axis only
ax.grid(which="minor", axis="y", linestyle=":", linewidth=0.5)   # Minor grid lines on y-axis only
#ax.grid(linestyle="--", linewidth=0.5)

# Save and show the plot
fig.savefig("error_vs_npower.jpg", dpi=1000)
