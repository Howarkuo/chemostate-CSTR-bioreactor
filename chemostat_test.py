
import numpy as np
import matplotlib.pyplot as plt
import control as ct

D = 1.0
Y_XS = 0.6
Y_PS = 0.05
mu_max = 2.0
Km = 0.5
Cs0_initial = 3.0

# Positive-biomass steady state
Cs_ss = D * Km / (Y_XS * mu_max - D)
Cx_ss = Y_XS * (Cs0_initial - Cs_ss)
Cp_ss = Y_PS * (Cs0_initial - Cs_ss)
y_ss = np.array([Cs_ss, Cx_ss, Cp_ss])



# Linearized process model (A, B matrices), using correct Monod derivative
b = mu_max * Cs_ss / (Km + Cs_ss)
k = mu_max * Km * Cx_ss / (Km + Cs_ss)**2
A = np.array([[-D - k, -b, 0],
              [Y_XS*k, -D + Y_XS*b, 0],
              [Y_PS*k, Y_PS*b, -D]])
B = np.array([[D], [0], [0]])
C = np.eye(3)
E = np.zeros((3, 1))

# Use your existing A, B, C, E, and y_ss
system_ss = ct.ss(A, B, C, E)

# Convert state-space to transfer functions
G = ct.ss2tf(system_ss)

G_S = G[0, 0]
G_X = G[1, 0]
G_P = G[2, 0]

t = np.linspace(0, 30, 1501)
delta_u = 2.8 - 3.0  # -0.2 g/L
u = np.full_like(t, delta_u)

# Method 1: State-space forced response
response = ct.forced_response(
    system_ss, T=t, U=u, X0=np.zeros(3)
)
y_forced = np.asarray(response.outputs)

# Method 2: Transfer-function step response
y_step = []

for Gi in [G_S, G_X, G_P]:
    _, yi = ct.step_response(Gi, T=t)
    y_step.append(delta_u * np.asarray(yi))

y_step = np.array(y_step)

# Convert deviations to absolute concentrations
y_forced_abs = y_ss[:, None] + y_forced
y_step_abs = y_ss[:, None] + y_step

# Compare numerical differences
print("Maximum difference:",
      np.max(np.abs(y_forced_abs - y_step_abs)))

# Plot
labels = ["Substrate", "Biomass", "Product"]

fig, axes = plt.subplots(
    3, 1, figsize=(9, 8), sharex=True
)

for i, ax in enumerate(axes):
    ax.plot(t, y_forced_abs[i], label="forced_response")
    ax.plot(t, y_step_abs[i], "--", label="step_response")
    ax.set_ylabel(f"{labels[i]} (g/L)")
    ax.legend()
    ax.grid(True)

axes[-1].set_xlabel("Time (h)")
plt.tight_layout()
plt.show()
