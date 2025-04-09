import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

dt = dtype=np.dtype([
   ("iter", int),
   ("Re", float),
   ("Im", float),
   ("modulus", float),
   ("residual", float),
   ("converged", "U12")
   ])

data = np.loadtxt('eigs_output.txt', dtype=dt, skiprows=1)

n = data.shape[0]

E = np.array([ d[1] + 1j*d[2] for row in data for d in [row] ])
conv = [ False for i in range(n) ]
for i in range(n):
   if data[i][5] == 'T':
      conv[i] = True

fig, axs = plt.subplots(1, 2, figsize=(20, 12))
ax = axs[0]
ax.axhline(0, color='black')
ax.axvline(0, color='black')
X, Y = np.meshgrid(np.linspace(-2, 2, 400), np.linspace(-2, 2, 400))
Z = X**2 + Y**2
ax.contourf(X, Y, Z, levels=[1, np.max(Z)], colors=['lightgrey'], alpha=0.5)
eit = np.exp(1j*np.linspace(0, 2*np.pi, 1000))
ax.plot(np.real(eit), np.imag(eit), color='black')
ax.scatter(np.real(E), np.imag(E), s=30, marker='+', color='black')
ax.scatter(np.real(E[conv]), np.imag(E[conv]), s=30, marker='o', color='red')
ax.axis('equal')
ax.set_title(r'Eigenvalues of $\exp(\tau L)$')
ax.set_xlim([-1.25, 1.25])
ax.set_ylim([-1.25, 1.25])

ax = axs[1]
tau = 1.0
Er = np.log(E)/tau
ax.axhline(0, color='black')
ax.axvline(0, color='black')
vertices = [(5, 0), (5, 1), (-5, 1), (-5, 0)]
rectangle = Polygon(vertices, closed=True, color='lightgrey', alpha=0.5, label='unstable region')
ax.add_patch(rectangle)
ax.scatter(np.imag(Er),       np.real(Er), s=30, marker='+', color='black', label='unconverged')
ax.scatter(np.imag(Er[conv]), np.real(Er[conv]), s=50, marker='o', color='red', label='converged')
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.5, 0.03])
ax.legend()
ax.set_title(r'Eigenvalues of $L$')
ax.set_ylabel('Re')
ax.set_xlabel('Im')

plt.show()
