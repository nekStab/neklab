#!/usr/bin/env python3
import re
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from datetime import datetime
import matplotlib.patches as mpatches
from scipy import stats
from scipy.optimize import curve_fit

def parse_log_file(log_file='lightkrylov.log'):
    """
    Parse the lightkrylov.log file to extract iteration data for both Newton and GMRES steps.
    """
    # Data structures to store parsed information
    newton_steps = []
    newton_residuals = []
    current_newton_step = 0
    
    # For each Newton step, track the GMRES data
    gmres_data = {}  # Dictionary keyed by Newton step number
    
    # Regular expressions to match the different step lines
    newton_step_pattern = re.compile(r'.*: LK_NwtKryl % newton_rdp: Start step (\d+): rnorm=\s+([0-9.E+-]+), tol=\s+([0-9.E+-]+)')
    inner_step_pattern = re.compile(r'.*: LK_Solvers % gmres_rdp: INFO: GMRES\(k\)\s+inner step\s+(\d+): \|res\|=\s+([0-9.E+-]+), tol=\s+([0-9.E+-]+)')
    init_step_pattern = re.compile(r'.*: LK_Solvers % gmres_rdp: INFO: GMRES\(k\)\s+init step\s+: \|res\|=\s+([0-9.E+-]+), tol=\s+([0-9.E+-]+)')
    outer_step_pattern = re.compile(r'.*: LK_Solvers % gmres_rdp: INFO: GMRES\(k\) outer step\s+(\d+): \|res\|=\s+([0-9.E+-]+), tol=\s+([0-9.E+-]+)')
    newton_converged_pattern = re.compile(r'.*: LK_NwtKryl % newton_rdp: Newton iteration converged after (\d+) iterations')
    
    with open(log_file, 'r') as f:
        for line in f:
            # Check for Newton step
            newton_match = newton_step_pattern.search(line)
            if newton_match:
                step_num = int(newton_match.group(1))
                residual = float(newton_match.group(2))
                tolerance = float(newton_match.group(3))
                
                current_newton_step = step_num
                newton_steps.append(step_num)
                newton_residuals.append(residual)
                
                # Initialize GMRES data for this Newton step
                gmres_data[current_newton_step] = {
                    'steps': [],
                    'residuals': [],
                    'tolerance': tolerance
                }
                continue
            
            # Extract initial GMRES step info
            init_match = init_step_pattern.search(line)
            if init_match:
                if current_newton_step > 0:
                    residual = float(init_match.group(1))
                    gmres_data[current_newton_step]['steps'].append(0)
                    gmres_data[current_newton_step]['residuals'].append(residual)
                continue
            
            # Extract inner GMRES step info
            inner_match = inner_step_pattern.search(line)
            if inner_match:
                if current_newton_step > 0:
                    step_num = int(inner_match.group(1))
                    residual = float(inner_match.group(2))
                    
                    gmres_data[current_newton_step]['steps'].append(step_num)
                    gmres_data[current_newton_step]['residuals'].append(residual)
                continue
            
            # Extract outer GMRES step info (typically the final GMRES result for this Newton iteration)
            outer_match = outer_step_pattern.search(line)
            if outer_match:
                # We don't need to process this separately as we already have the inner steps
                pass
    
    return {
        'newton_steps': newton_steps,
        'newton_residuals': newton_residuals,
        'gmres_data': gmres_data
    }

def plot_basic_residuals(data, output_file='residual.png', logscale=True):
    """
    Plot the basic residual values:
    1. Overall Newton iteration convergence
    2. GMRES inner step convergence for each Newton iteration
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Colors for different Newton iterations
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    markers = ['o', 's', '^', 'D', 'v', '<', '>']
    
    #--------------------------------------------------------------
    # Plot 1: Newton iteration convergence (linear)
    #--------------------------------------------------------------
    newton_steps = data['newton_steps']
    newton_residuals = data['newton_residuals']
    
    ax1.plot(newton_steps, newton_residuals, 'ko-', markersize=10, linewidth=2)
    ax1.set_xlabel('Newton Iteration', fontsize=12)
    ax1.set_ylabel('Residual Norm', fontsize=12)
    ax1.set_title('Newton Iteration Convergence', fontsize=14)
    ax1.grid(True, which='both', linestyle='--', alpha=0.7)
    
    if logscale:
        ax1.set_yscale('log')
    
    #--------------------------------------------------------------
    # Plot 2: GMRES convergence for each Newton iteration
    #--------------------------------------------------------------
    gmres_data = data['gmres_data']
    legend_handles = []
    
    for newton_step, step_data in sorted(gmres_data.items()):
        if not step_data['steps'] or not step_data['residuals']:
            continue
            
        color_idx = (newton_step - 1) % len(colors)
        marker_idx = (newton_step - 1) % len(markers)
        
        ax2.plot(step_data['steps'], step_data['residuals'], 
                 marker=markers[marker_idx], color=colors[color_idx], 
                 linewidth=1.5, markersize=7)
        
        # Add tolerance line
        if 'tolerance' in step_data:
            ax2.axhline(y=step_data['tolerance'], color=colors[color_idx], 
                       linestyle='--', alpha=0.5)
        
        # Create custom legend handle
        patch = mpatches.Patch(color=colors[color_idx], 
                               label=f'Newton Step {newton_step}')
        legend_handles.append(patch)
    
    ax2.set_xlabel('GMRES Iteration', fontsize=12)
    ax2.set_ylabel('GMRES Residual', fontsize=12)
    ax2.set_title('GMRES Convergence per Newton Iteration', fontsize=14)
    ax2.grid(True, which='both', linestyle='--', alpha=0.7)
    ax2.legend(handles=legend_handles, loc='best')
    
    if logscale:
        ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Basic plot saved to {output_file}")

def power_law(x, a, b):
    """Power law function y = a*x^b"""
    return a * np.power(x, b)

def log_linear(x, a, b):
    """Linear function in log-log space: log(y) = log(a) + b*log(x)"""
    return a + b * x

def classify_convergence(rate):
    """
    Classify convergence rate based on measured value
    
    Convergence Classifications:
    - Sublinear/Algebraic (p < 1): Slower than linear, includes rates like O(1/n)
    - Linear (p ≈ 1): Errors reduce by constant factor each iteration
    - Superlinear/Subgeometric (1 < p < 2): Faster than linear but slower than quadratic
    - Quadratic (p ≈ 2): Errors are squared each iteration (typical of Newton's method)
    - Supergeometric (p > 2): Faster than quadratic convergence
    """
    if rate < 0.8:
        return "Sublinear/Algebraic (p < 1)"
    elif 0.8 <= rate < 1.2:
        return "Linear (p ≈ 1)"
    elif 1.2 <= rate < 1.8:
        return "Superlinear/Subgeometric (1 < p < 2)"
    elif 1.8 <= rate < 2.2:
        return "Quadratic (p ≈ 2)"
    else:
        return "Supergeometric (p > 2)"

def plot_quadratic_analysis(data, output_file='residual_quadratic.png', logscale=True):
    """
    Plot advanced convergence analysis to demonstrate quadratic convergence of Newton's method
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # Add a bit more space at the bottom for the legends
    plt.subplots_adjust(bottom=0.2)
    
    newton_steps = data['newton_steps']
    newton_residuals = data['newton_residuals']
    
    #--------------------------------------------------------------
    # Plot 1: Log-log plot of Newton residuals to detect convergence pattern
    #--------------------------------------------------------------
    if len(newton_steps) > 1:
        x_loglog = np.array(newton_steps) - min(newton_steps) + 1
        ax1.loglog(x_loglog, newton_residuals, 'ko-', markersize=8, linewidth=2, label='Newton Residuals')
        
        # Add reference lines for comparison
        x_ref = np.linspace(1, max(x_loglog), 100)
        
        # Sublinear convergence reference (slope = -0.5 on log-log plot)
        y_sublinear = newton_residuals[0] * 0.3**(x_ref-1)
        ax1.loglog(x_ref, y_sublinear, 'm--', linewidth=1, label='Sublinear (r ~ 0.3^n)')
        
        # Linear convergence reference (slope = -1 on log-log plot)
        y_linear = newton_residuals[0] * 0.1**(x_ref-1)
        ax1.loglog(x_ref, y_linear, 'b--', linewidth=1, label='Linear (r ~ 0.1^n)')
        
        # Superlinear convergence reference (between linear and quadratic)
        y_superlinear = newton_residuals[0] * 0.03**(x_ref-1)
        ax1.loglog(x_ref, y_superlinear, 'c--', linewidth=1, label='Superlinear (r ~ 0.03^n)')
        
        # Quadratic convergence reference (superlinear slope on log-log plot)
        if len(newton_residuals) > 2:
            # Try to estimate a reasonable quadratic model
            factor = min(1.0, (newton_residuals[1]/newton_residuals[0]**2) * 10)
            y_quadratic = newton_residuals[0] * (factor * np.array(newton_residuals[0])**(2**(x_ref-1) - 1))
            ax1.loglog(x_ref, y_quadratic, 'r--', linewidth=1, label='Quadratic (r ~ r_prev^2)')
            
            # Supergeometric reference (faster than quadratic)
            y_supergeometric = newton_residuals[0] * (factor * np.array(newton_residuals[0])**(3**(x_ref-1) - 1))
            ax1.loglog(x_ref, y_supergeometric, 'g--', linewidth=1, label='Supergeometric (r ~ r_prev^3)')
            
            # Fit a function to the data points to estimate convergence rate
            # We'll use an exponential decay model: r ~ C * a^n
            # In log space: log(r) ~ log(C) + n*log(a)
            log_residuals = np.log10(newton_residuals)
            x_fit = np.array(newton_steps) - min(newton_steps)
            x_fit[x_fit == 0] = 1  # Avoid log(0)
            
            # Only fit if we have enough points
            if len(x_fit) >= 3:
                try:
                    # Fit log-linear model
                    popt, pcov = curve_fit(lambda x, a, b: a + b*x, x_fit, log_residuals)
                    a_fit, b_fit = popt
                    fitted_rate = 10**b_fit  # Convert from log10 base
                    
                    # Plot the fitted line
                    x_dense = np.linspace(min(x_fit), max(x_fit), 100)
                    y_fitted = 10**(a_fit + b_fit * x_dense)
                    ax1.loglog(x_dense + min(newton_steps), y_fitted, 'g-', linewidth=1.5, 
                              label=f'Fitted: r ≈ {10**a_fit:.2e} × {fitted_rate:.4f}ⁿ')
                    
                    # Classify convergence type based on rate
                    convergence_class = "Unknown"
                    if fitted_rate > 0.3:
                        convergence_class = "Sublinear/Algebraic"
                    elif 0.05 < fitted_rate <= 0.3:
                        convergence_class = "Linear"
                    elif 0.01 < fitted_rate <= 0.05:
                        convergence_class = "Superlinear/Subgeometric"
                    elif 0.001 < fitted_rate <= 0.01:
                        convergence_class = "Quadratic"
                    else:
                        convergence_class = "Supergeometric"
                    
                    # Annotate with the fitted convergence rate
                    ax1.text(0.05, 0.05, 
                            f'Fitted rate: {fitted_rate:.4f}\n'
                            f'Classification: {convergence_class}\n'
                            f'Lin≈0.1, Quad≈0.01, Super≈0.001', 
                            transform=ax1.transAxes, fontsize=10, 
                            bbox=dict(facecolor='white', alpha=0.8))
                except:
                    # If fitting fails, just continue without the fit
                    pass
    
        ax1.set_xlabel('Iteration Number', fontsize=12)
        ax1.set_ylabel('Residual Norm (log scale)', fontsize=12)
        ax1.set_title('Newton Convergence Analysis (Log-Log)', fontsize=14)
        ax1.grid(True, which='both', linestyle='--', alpha=0.7)
        ax1.legend(loc='best')
        
        # Set y-limit to double precision accuracy (~1e-16)
        ax1.set_ylim(bottom=1e-16)
    
    #--------------------------------------------------------------
    # Plot 2: Residual_n+1 vs Residual_n to detect quadratic convergence
    #--------------------------------------------------------------
    if len(newton_residuals) > 1:
        res_current = newton_residuals[:-1]
        res_next = newton_residuals[1:]
        
        # For quadratic convergence, log(res_n+1) vs log(res_n) should have slope ~2
        ax2.loglog(res_current, res_next, 'ko-', markersize=8, linewidth=2, label='Residuals')
        
        # Add reference lines
        x_ref = np.logspace(np.log10(min(res_current)), np.log10(max(res_current)), 100)
        
        # Determine a common point for both reference lines to pass through
        # Use the first data point as reference
        common_x = res_current[0]
        common_y = res_next[0]
        
        # Reference lines for different convergence rates
        # Sublinear reference (slope = 0.5)
        y_sublin = (x_ref/common_x)**0.5 * common_y
        ax2.loglog(x_ref, y_sublin, 'm--', linewidth=1.5, label='p = 0.5 (Sublinear)')
        
        # Linear reference (slope = 1)
        y_lin = (x_ref/common_x) * common_y
        ax2.loglog(x_ref, y_lin, 'b--', linewidth=1.5, label='p = 1.0 (Linear)')
        
        # Superlinear reference (slope = 1.5)
        y_superlin = (x_ref/common_x)**1.5 * common_y
        ax2.loglog(x_ref, y_superlin, 'c--', linewidth=1.5, label='p = 1.5 (Superlinear)')
        
        # Quadratic reference (slope = 2)
        y_quad = (x_ref/common_x)**2 * common_y
        ax2.loglog(x_ref, y_quad, 'r--', linewidth=1.5, label='p = 2.0 (Quadratic)')
        
        # Supergeometric reference (slope = 3)
        y_supergeom = (x_ref/common_x)**3 * common_y
        ax2.loglog(x_ref, y_supergeom, 'g--', linewidth=1.5, label='p = 3.0 (Supergeometric)')
        
        # Fit power law model to estimate actual convergence order
        if len(res_current) >= 3:
            try:
                # Fit power law: y = a*x^b (in log-log space: log(y) = log(a) + b*log(x))
                log_x = np.log10(res_current)
                log_y = np.log10(res_next)
                
                # Linear regression in log-log space
                slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)
                
                # Plot fitted power law
                a_fit = 10**intercept
                b_fit = slope
                y_fitted = a_fit * x_ref**b_fit
                ax2.loglog(x_ref, y_fitted, 'g-', linewidth=1.5, 
                          label=f'Fitted: r₊₁ ≈ {a_fit:.2e} · rₙᵖ, p={b_fit:.4f}')
                
                # Classify convergence based on slope
                convergence_type = classify_convergence(b_fit)
                
                # Annotate with the fitted slope (convergence order)
                ax2.text(0.05, 0.05, 
                        f'Fitted order: p = {b_fit:.4f}\n'
                        f'R² = {r_value**2:.4f}\n'
                        f'Classification: {convergence_type}', 
                        transform=ax2.transAxes, fontsize=10, 
                        bbox=dict(facecolor='white', alpha=0.8))
            except:
                # If fitting fails, just continue without the fit
                pass
        
        ax2.set_xlabel('Residual at step n', fontsize=12)
        ax2.set_ylabel('Residual at step n+1', fontsize=12)
        ax2.set_title('Convergence Rate Analysis', fontsize=14)
        ax2.grid(True, which='both', linestyle='--', alpha=0.7)
        ax2.legend(loc='best')
    
    #--------------------------------------------------------------
    # Plot 3: Convergence rate analysis
    #--------------------------------------------------------------
    if len(newton_residuals) > 2:
        # Calculate empirical convergence rate
        # For order p convergence: |x_{n+1} - x*| ≈ C|x_n - x*|^p
        # Assume x* ≈ 0 (converged solution has ~0 residual)
        # Then r_{n+1} ≈ C * r_n^p
        # Taking logs: log(r_{n+1}) ≈ log(C) + p*log(r_n)
        # So the slope of log(r_{n+1}) vs log(r_n) is approximately p
        
        log_res_current = np.log10(newton_residuals[:-1])
        log_res_next = np.log10(newton_residuals[1:])
        
        # Calculate slope for each pair of points
        convergence_rates = []
        iterations = []
        
        for i in range(len(log_res_current)-1):
            rate = (log_res_next[i+1] - log_res_next[i]) / (log_res_current[i+1] - log_res_current[i])
            convergence_rates.append(rate)
            iterations.append(newton_steps[i+1])
        
        # Plot convergence rates
        ax3.plot(iterations, convergence_rates, 'ko-', markersize=8, linewidth=2, label='Point-wise rates')
        
        # Reference lines for different convergence classifications
        ax3.axhline(y=0.5, color='m', linestyle='--', label='p = 0.5 (Sublinear)')
        ax3.axhline(y=1.0, color='b', linestyle='--', label='p = 1.0 (Linear)')
        ax3.axhline(y=1.5, color='c', linestyle='--', label='p = 1.5 (Superlinear)')
        ax3.axhline(y=2.0, color='r', linestyle='--', label='p = 2.0 (Quadratic)')
        ax3.axhline(y=3.0, color='g', linestyle='--', label='p = 3.0 (Supergeometric)')
        
        # Calculate and show average convergence rate
        if convergence_rates:
            avg_rate = sum(convergence_rates) / len(convergence_rates)
            ax3.axhline(y=avg_rate, color='m', linestyle='-', linewidth=1.5, 
                       label=f'Average order: p = {avg_rate:.4f}')
            
            # Classify average convergence
            avg_convergence_type = classify_convergence(avg_rate)
            
            # Annotate with the statistics
            ax3.text(0.05, 0.05, 
                    f'Average order: p = {avg_rate:.4f}\n'
                    f'Range: [{min(convergence_rates):.2f}, {max(convergence_rates):.2f}]\n'
                    f'Classification: {avg_convergence_type}', 
                    transform=ax3.transAxes, fontsize=10, 
                    bbox=dict(facecolor='white', alpha=0.8))
        
        ax3.set_xlabel('Newton Iteration', fontsize=12)
        ax3.set_ylabel('Convergence Rate (p)', fontsize=12)
        ax3.set_title('Estimated Order of Convergence', fontsize=14)
        ax3.grid(True, which='both', linestyle='--', alpha=0.7)
        ax3.legend(loc='best')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Quadratic analysis plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Plot convergence data from LightKrylov log file')
    parser.add_argument('--output', '-o', default='residual.png', help='Output file prefix to save plots (default: residual.png and residual_quadratic.png)')
    parser.add_argument('--linear', '-l', action='store_true', help='Use linear scale instead of log scale')
    parser.add_argument('--file', '-f', default='lightkrylov.log', help='Path to the lightkrylov.log file (default: lightkrylov.log)')
    
    args = parser.parse_args()
    
    log_file = args.file
    if not os.path.isfile(log_file):
        print(f"Error: File {log_file} does not exist")
        return 1
    
    # Parse the log file
    data = parse_log_file(log_file)
    
    # Get output file names
    output_base = os.path.splitext(args.output)[0]
    basic_output = f"{output_base}.png"
    quadratic_output = f"{output_base}_quadratic.png"
    
    # Plot residuals
    plot_basic_residuals(data, basic_output, not args.linear)
    plot_quadratic_analysis(data, quadratic_output, not args.linear)
    
    return 0

if __name__ == "__main__":
    main()
