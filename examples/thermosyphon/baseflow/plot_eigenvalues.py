#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import glob
import os

# Function to create a base plot for circle eigenvalue plots
def create_circle_plot():
    fig = plt.figure(figsize=(10, 10))
    ax = plt.gca()
    
    # Plot unit circle as reference
    unit_circle = Circle((0, 0), 1, fill=False, color='gray', linestyle='--', alpha=0.5)
    ax.add_patch(unit_circle)
    
    # Set equal aspect ratio
    plt.axis('equal')
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add labels
    plt.xlabel('Re(λ)', fontsize=14)
    plt.ylabel('Im(λ)', fontsize=14)
    
    # Add coordinate lines
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    return fig, ax

# NPY Plot - purely real eigenvalues (eigsNS.png)
plt.figure(figsize=(12, 6))

npy_files = glob.glob("*_eigenspectrum.npy")
npy_data_found = False

for npy_file in npy_files:
    try:
        # Load the NPY file
        data = np.load(npy_file)
        
        if data.ndim == 2:
            # Only use the real parts from column 0
            real_parts = data[:, 0]
            
            # Plot eigenvalues on a number line
            plt.stem(np.arange(len(real_parts)), real_parts, 
                    linefmt='b-', markerfmt='bo', basefmt='k-',
                    label=f'From {os.path.basename(npy_file)}')
            
            # Add value labels
            for i, v in enumerate(real_parts):
                plt.annotate(f'{v:.6f}', (i, v), 
                            xytext=(0, 5), textcoords='offset points',
                            ha='center', fontsize=10)
            
            npy_data_found = True
            print(f"Loaded {len(data)} eigenvalues from {npy_file}")
        else:
            print(f"Warning: {npy_file} has unexpected structure")
    except Exception as e:
        print(f"Error loading {npy_file}: {e}")

# Add labels and title
plt.xlabel('Index', fontsize=14)
plt.ylabel('Eigenvalue (Real)', fontsize=14)
plt.title('Real Eigenvalue Spectrum (NPY files)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.5)

# Add legend if we have data from multiple files
handles, labels = plt.gca().get_legend_handles_labels()
if len(handles) > 0:
    plt.legend(loc='best', fontsize=12)

# Save the NPY figure if data was found
if npy_data_found:
    plt.tight_layout()
    plt.savefig('eigsNS.png', dpi=300)
    print("NPY data plot saved as eigsNS.png")
else:
    print("No valid NPY data found, eigsNS.png not saved")
    plt.close()

# Text File Plot - eigs_output.txt (eigsH.png)
fig_txt, ax_txt = create_circle_plot()
plt.title('Eigenvalue Spectrum (Text file)', fontsize=16)

filename = 'eigs_output.txt'
if os.path.exists(filename):
    try:
        data = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 5:  # Ensure line has enough data
                        try:
                            re = float(parts[1])
                            im = float(parts[2])
                            residual = float(parts[4])
                            
                            data.append([re, im, residual])
                        except ValueError:
                            # Skip lines that can't be parsed
                            continue
        
        if data:
            data = np.array(data)
            real_parts = data[:, 0]
            imag_parts = data[:, 1]
            residuals = data[:, 2]
            
            # Plot with residual coloring
            sc = plt.scatter(
                real_parts, 
                imag_parts,
                c=residuals,
                cmap='plasma',
                s=100,
                marker='x',
                label=f'From {filename}'
            )
            
            plt.colorbar(sc, label='Residual')
            
            # Save the text file figure
            plt.tight_layout()
            plt.savefig('eigsH.png', dpi=300)
            print(f"Loaded {len(data)} eigenvalues from {filename}")
            print("Text file plot saved as eigsH.png")
        else:
            print(f"No data found in {filename}, eigsH.png not saved")
    except Exception as e:
        print(f"Error processing {filename}: {e}")
else:
    print(f"File {filename} not found, eigsH.png not saved")
