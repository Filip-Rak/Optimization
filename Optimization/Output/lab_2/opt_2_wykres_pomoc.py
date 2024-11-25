import numpy as np
import matplotlib.pyplot as plt

# funkcji celu
def objective_function(x1, x2):
    return x1**2 + x2**2 - np.cos(2.5 * np.pi * x1) - np.cos(2.5 * np.pi * x2) + 2

# zakres dla x1 i x2
x1_vals = np.linspace(-1.5, 1.5, 500)
x2_vals = np.linspace(-1.5, 1.5, 500)
x1_vals_2 = np.linspace(-1.25, 0.25, 500)
x2_vals_2 = np.linspace(-0.5, 1.0, 500)

# siatka wartości funkcji celu
X1, X2 = np.meshgrid(x1_vals, x2_vals)
Z = objective_function(X1, X2)

X1_2, X2_2 = np.meshgrid(x1_vals_2, x2_vals_2)
Z_2 = objective_function(X1_2, X2_2)

# współrzędne punktów iteracyjnych z algorytmu Hooke’a-Jeevesa (dane z excela)
hooke_jeeves_points = np.array([
    [-0.933653, 0.582141], [-0.183653, -0.167859], [0.00384724, 0.0196406],
    [0.00384724, -0.00379693], [-0.00201214, 0.00206244], [0.000917552, -0.000867244],
    [-0.000547292, 0.0005976], [0.00018513, -0.000134822], [2.02476e-06, 4.82831e-05],
    [-8.36266e-07, -3.54292e-07], [5.94246e-07, -3.54292e-07]
])

# współrzędne punktów iteracyjnych z algorytmu Rosenbrocka
rosenbrock_points = np.array([
    [-0.933653, 0.582141], [-0.183653, 0.582141], [-0.183653, 0.207141],
    [-0.0159477, -0.064212], [0.0259786, 0.0196406], [0.015497, -0.00132257],
    [0.0147269, -7.64985e-05], [-0.0052102, -0.0123983], [-0.000153932, -0.00182651],
    [-0.00081467, -0.00151049], [-0.00100358, -0.000176136], [-0.000439404, 0.00029093],
    [-0.00029535, -4.57578e-05], [5.30734e-05, 5.94204e-05], [-3.4573e-05, 3.29627e-05],
    [-3.04291e-05, 1.04528e-05], [-1.06854e-05, 2.20309e-05], [-7.37814e-06, 1.10751e-05],
    [9.85086e-06, 1.5121e-06], [4.84782e-06, 4.28905e-06], [6.75948e-07, -1.12067e-06],
    [6.51838e-07, 3.09637e-07], [3.85296e-08, -5.83885e-08]
])

########## symetryczny wykres
# poziomice funkcji celu
plt.figure(figsize=(10, 8))
contour = plt.contour(X1, X2, Z, levels=50, cmap='viridis')
plt.colorbar(contour)
plt.title("Poziomice funkcji celu z punktami iteracji")
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")

# naniesienie punktów iteracyjnych z algorytmu Hooke’a-Jeevesa
plt.plot(hooke_jeeves_points[:, 0], hooke_jeeves_points[:, 1], 'o-', color='red', label='Hooke-Jeeves')

# naniesienie punktów iteracyjnych z algorytmu Rosenbrocka
plt.plot(rosenbrock_points[:, 0], rosenbrock_points[:, 1], 'o-', color='blue', label='Rosenbrock')

# początkowe i końcowe punkty
plt.plot(hooke_jeeves_points[0, 0], hooke_jeeves_points[0, 1], 'ro', markersize=8, label='Start Hooke-Jeeves')
plt.plot(rosenbrock_points[0, 0], rosenbrock_points[0, 1], 'bo', markersize=8, label='Start Rosenbrock')
plt.plot(hooke_jeeves_points[-1, 0], hooke_jeeves_points[-1, 1], 'r*', markersize=10, label='End Hooke-Jeeves')
plt.plot(rosenbrock_points[-1, 0], rosenbrock_points[-1, 1], 'b*', markersize=10, label='End Rosenbrock')

plt.legend()
plt.show()

########## wyraźniejszy wykres
plt.figure(figsize=(10, 8))
contour_2 = plt.contour(X1_2, X2_2, Z_2, levels=50, cmap='viridis')
plt.colorbar(contour_2)
plt.title("Poziomice funkcji celu z punktami iteracji")
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")

# naniesienie punktów iteracyjnych z algorytmu Hooke’a-Jeevesa
plt.plot(hooke_jeeves_points[:, 0], hooke_jeeves_points[:, 1], 'o-', color='red', label='Hooke-Jeeves')

# naniesienie punktów iteracyjnych z algorytmu Rosenbrocka
plt.plot(rosenbrock_points[:, 0], rosenbrock_points[:, 1], 'o-', color='blue', label='Rosenbrock')

# początkowe i końcowe punkty
plt.plot(hooke_jeeves_points[0, 0], hooke_jeeves_points[0, 1], 'ro', markersize=8, label='Start Hooke-Jeeves')
plt.plot(rosenbrock_points[0, 0], rosenbrock_points[0, 1], 'bo', markersize=8, label='Start Rosenbrock')
plt.plot(hooke_jeeves_points[-1, 0], hooke_jeeves_points[-1, 1], 'r*', markersize=10, label='End Hooke-Jeeves')
plt.plot(rosenbrock_points[-1, 0], rosenbrock_points[-1, 1], 'b*', markersize=10, label='End Rosenbrock')

plt.legend()
plt.show()
