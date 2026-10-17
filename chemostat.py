#Chemostat- microorganism-feed-product stirred-tank bioreactor (CSTR)

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import control as ct
HAVE_CONTROL = True

# try:
#     import control as ct
#     HAVE_CONTROL = True
# except ImportError:
#     # For environments without python-control, the exact same state-space
#     # response is computed with SciPy.
#     from scipy.signal import StateSpace, lsim
#     HAVE_CONTROL = False

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

# Nonlinear model (state order S, X, P)
# f(t,y)
def nonlinear_ode(t, y, Cs0):
    Cs, Cx, Cp = y
    r = mu_max * Cs * Cx / (Km + Cs)
    return [D * (Cs0 - Cs) - r,
            -D * Cx + Y_XS * r,
            -D * Cp + Y_PS * r]

# Linearized process model (A, B matrices), using correct Monod derivative
b = mu_max * Cs_ss / (Km + Cs_ss)
k = mu_max * Km * Cx_ss / (Km + Cs_ss)**2
A = np.array([[-D - k, -b, 0],
              [Y_XS*k, -D + Y_XS*b, 0],
              [Y_PS*k, Y_PS*b, -D]])
B = np.array([[D], [0], [0]])
C = np.eye(3)
E = np.zeros((3, 1))

print('Initial steady state (S, X, P):', y_ss)
print('Linear model A =\n', A)

# scipy.integrate.solve_ivp(
#     fun, t_span, y0, method='RK45', t_eval=None, 
#     dense_output=False, events=None, vectorized=False, args=None, **options
# 
def simulate(Cs0_after, t_end, basename):
    t = np.linspace(0, t_end, 1501)
    sol = solve_ivp(nonlinear_ode, (0, t_end), y_ss,
                    args=(Cs0_after,), t_eval=t,
                    rtol=1e-10, atol=1e-12)
    # Tolerances (rtol / atol)
    # System Default: rtol = 1e-3 (0.1% relative error) and atol = 1e-6.
    assert sol.success
    y_nonlinear = sol.y

    # Transfer-function/state-space representation of the coupled linear model
    # LLM suggest control.forced_response()
    # Use the control.step_response() from class instead
    # du = np.full_like(t, Cs0_after - Cs0_initial)
    # if HAVE_CONTROL:
    #     system = ct.ss(A, B, C, E)
    #     result = ct.forced_response(system, T=t, U=du, X0=np.zeros(3))
    #     y_deviation = np.asarray(result.outputs)
    # else:
    #     system = StateSpace(A, B, C, E)
    #     _, y_deviation, _ = lsim(system, U=du, T=t)
    #     y_deviation = y_deviation.T
    # y_linear = y_ss[:, None] + y_deviation

    # Linear model: step_response instead of forced_response

# Step change in inlet substrate concentration
    delta_Cs0 = Cs0_after - Cs0_initial

# Construct linear state-space model
    system = ct.ss(A, B, C, E)

# Calculate response to a UNIT step (amplitude = 1)
    response = ct.step_response(system, T=t)

# Result shape: (outputs, inputs, time)
# One input, three outputs
    y_unit_step = np.asarray(response.outputs)[:, 0, :]

# Scale unit-step response to actual input change
    y_deviation = delta_Cs0 * y_unit_step

# Convert deviations to absolute concentrations
    y_linear = y_ss[:, None] + y_deviation


    

    labels = ['Substrate $C_S$', 'Biomass $C_X$', 'Product $C_P$']
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True, constrained_layout=True)
    for i, ax in enumerate(axes):
        ax.plot(t, y_nonlinear[i], linewidth=2.2, label='Nonlinear ODE')
        ax.plot(t, y_linear[i], '--', linewidth=2, label='Linear transfer-function model')
        ax.set_ylabel(f'{labels[i]} (g/L)')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    axes[-1].set_xlabel('Time (h)')
    fig.suptitle(f'Feed substrate step: {Cs0_initial:g} → {Cs0_after:g} g/L,   0–{t_end:g} h', fontsize=13)
    output = (f'{basename}.png')
    fig.savefig(output, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'{basename}: final nonlinear S,X,P =', np.round(y_nonlinear[:,-1], 6))
    print(f'{basename}: final linear    S,X,P =', np.round(y_linear[:,-1], 6))
    return output

if __name__ == '__main__':
    simulate(2.8, 30, 'part_e_bioreactor')
    simulate(0.5, 10, 'part_f_bioreactor')
# (base) kuochingchia@kuochingchias-MacBook-Air-2966 CHEMENG162 % poetry run python chemostat.py
# Initial steady state (S, X, P): [2.5   0.3   0.025]
# Linear model A =
#  [[-1.03333333e+00 -1.66666667e+00  0.00000000e+00]
#  [ 2.00000000e-02  0.00000000e+00  0.00000000e+00]
#  [ 1.66666667e-03  8.33333333e-02 -1.00000000e+00]]
# part_e_bioreactor: final nonlinear S,X,P = [2.416201 0.23028  0.01919 ]
# part_e_bioreactor: final linear    S,X,P = [2.423887 0.225668 0.018806]
# part_f_bioreactor: final nonlinear S,X,P = [0.489588 0.006316 0.000526]
# part_f_bioreactor: final linear    S,X,P = [ 0.647019 -0.088143 -0.007345]