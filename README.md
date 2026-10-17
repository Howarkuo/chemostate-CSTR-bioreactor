# Chemostat Bioreactor: Nonlinear Dynamics and Linearization

**Substrate-limited microbial growth in a continuous stirred-tank reactor (CSTR)**
![Screenshot 2026-10-17 at 12.56.52 AM](https://hackmd.io/_uploads/ryNH5oe2Gl.png)

This project models substrate consumption, biomass growth, and product formation in a chemostat. It compares a **nonlinear Monod-kinetics model** integrated with SciPy against a **first-order linearized state-space model** simulated with `python-control`.

## Engineering problem statement
How does changing the inlet substrate concentration $\(C_{S0}\)$ affect the transient and steady-state concentrations of substrate $\(C_S\)$, biomass $\(C_X\)$, and product $\(C_P\)$ in a continuous bioreactor? How accurately can a linearized state-space model predict these changes compared with the original nonlinear Monod kinetics model?


## Concepts

Monod kinetics · Chemostat/CSTR · Substrate utilization · Microbial growth · Yield coefficients · Dynamic species balances · Nonlinear coupled ODEs · Steady states and washout · First-order Taylor linearization · Jacobian matrices · Deviation variables · State-space systems · Step responses · Eigenvalues and stability · Numerical ODE integration · Nonlinear versus linear model validation



## Governing equations

### Monod kinetics and substrate consumption
![Monods-function-of-saturation-curve](https://hackmd.io/_uploads/rkllRYil2Gx.webp)


The substrate consumption rate is modeled as

$$
\boxed{r=-r_S=\mu(C_S)C_X=\frac{\mu_{\max}C_S}{K_m+C_S}C_X}
$$

where:

- $r$: positive volumetric substrate consumption rate, g/(L·h).
- $C_S$: substrate concentration, g/L.
- $C_X$: biomass concentration, g/L.
- $\mu_{\max}$: maximum **substrate-specific utilization coefficient** in this model, h⁻¹.
- $K_m$: half-saturation constant, g/L, giving $\mu(K_m)=\mu_{\max}/2$.

> **Model convention:** Here, $\mu C_X$ denotes **substrate consumption**, and biomass production is $Y_{X/S}\mu C_X$. Accordingly, $\mu$ should not be confused with the usual biomass-specific growth rate $\mu_g=Y_{X/S}\mu$. This distinction ensures the rate and yield equations are consistent.

**Substrate-limited regime** ($C_S\ll K_m$):

$$
r\approx\frac{\mu_{\max}}{K_m}C_SC_X
$$

The rate is first-order in each concentration (second-order overall).

**Substrate-saturated regime** ($C_S\gg K_m$):

$$
r\approx\mu_{\max}C_X
$$

The rate becomes approximately independent of substrate concentration.

### Biomass and product formation

Assuming constant yield coefficients,

$$
\begin{aligned}
r_X &=Y_{X/S}r=Y_{X/S}\frac{\mu_{\max}C_S}{K_m+C_S}C_X\\
r_P &=Y_{P/S}r=Y_{P/S}\frac{\mu_{\max}C_S}{K_m+C_S}C_X
\end{aligned}
$$

with

$$
Y_{X/S}=\frac{r_X}{-r_S},\qquad Y_{P/S}=\frac{r_P}{-r_S}.
$$

### Dynamic material balances

For a constant-volume, perfectly mixed reactor, define the dilution rate as

$$
D=\frac{F}{V}\quad[\mathrm{h}^{-1}].
$$

Assuming no biomass or product in the inlet stream,

$$
\boxed{
\begin{aligned}
\frac{dC_S}{dt}&=D(C_{S0}-C_S)-r\\
\frac{dC_X}{dt}&=-DC_X+Y_{X/S}r\\
\frac{dC_P}{dt}&=-DC_P+Y_{P/S}r
\end{aligned}}
$$

Here $C_{S0}$ is the substrate concentration in the inlet stream.

## Process variables and parameters

| Category | Variable | Description |
|---|---|---|
| Manipulated input (MV) |$F$ | Volumetric Flow Rate |
| Controlled/observed outputs (CVs) | $C_S$, $C_X$, $C_P$ | Reactor substrate, biomass, and product concentrations |
| Potential disturbances (DVs) |  $C_{S0}$ | Inital Reactor substrate concentration  |

| Parameter | Value | Units |
|---|---:|---|
| Dilution rate, $D$ | 1.0 | h⁻¹ |
| Biomass yield, $Y_{X/S}$ | 0.6 | g biomass/g substrate |
| Product yield, $Y_{P/S}$ | 0.05 | g product/g substrate |
| Utilization coefficient, $\mu_{\max}$ | 2.0 | h⁻¹ |
| Half-saturation constant, $K_m$ | 0.5 | g/L |
| Initial feed concentration, $C_{S0}$ | 3.0 | g/L |
| Feed concentration after step | 2.8 | g/L |
| Simulation horizon | 30 | h |

## Steady-state calculation

At a non-washout steady state ($C_{X,ss}\ne 0$), the biomass balance gives

$$
0=-DC_{X,ss}+Y_{X/S}\frac{\mu_{\max}C_{S,ss}}{K_m+C_{S,ss}}C_{X,ss}.
$$

Divide by $C_{X,ss}$ and rearrange:

$$
D=Y_{X/S}\frac{\mu_{\max}C_{S,ss}}{K_m+C_{S,ss}}
$$

$$
\boxed{C_{S,ss}=\frac{DK_m}{Y_{X/S}\mu_{\max}-D}}
$$

With the specified parameters,

$$
C_{S,ss}=\frac{(1)(0.5)}{(0.6)(2)-1}=\boxed{2.5\ \mathrm{g/L}}.
$$

From the substrate steady-state balance,

$$
\boxed{C_{X,ss}=Y_{X/S}(C_{S0}-C_{S,ss})}
$$

$$
C_{X,ss}=0.6(3-2.5)=\boxed{0.3\ \mathrm{g/L}}.
$$

From the product steady-state balance,

$$
\boxed{C_{P,ss}=Y_{P/S}(C_{S0}-C_{S,ss})}
$$

$$
C_{P,ss}=0.05(3-2.5)=\boxed{0.025\ \mathrm{g/L}}.
$$

Therefore the pre-step equilibrium is

$$
\boxed{\mathbf y_{ss}=\begin{bmatrix}2.5&0.3&0.025\end{bmatrix}^{T}\ \mathrm{g/L}}.
$$

> The non-washout solution requires $Y_{X/S}\mu_{\max}>D$ and $C_{S0}>C_{S,ss}$. A washout steady state ($C_X=C_P=0$, $C_S=C_{S0}$) also exists.

## Nonlinear dynamic simulation

At $t=0$, the feed concentration decreases from **3.0 to 2.8 g/L**. The initial reactor concentrations remain at their pre-step steady-state values:

$$
\mathbf y(0)=\begin{bmatrix}2.5&0.3&0.025\end{bmatrix}^{T}\ \mathrm{g/L}.
$$

The initial consumption rate is

$$
r(0)=\frac{2(2.5)(0.3)}{0.5+2.5}=0.5\ \mathrm{g/(L\cdot h)}.
$$

Immediately after the step,

$$
\begin{aligned}
\left.\frac{dC_S}{dt}\right|_{0^+}&=1(2.8-2.5)-0.5=\boxed{-0.2}\\
\left.\frac{dC_X}{dt}\right|_{0^+}&=-0.3+0.6(0.5)=\boxed{0}\\
\left.\frac{dC_P}{dt}\right|_{0^+}&=-0.025+0.05(0.5)=\boxed{0}
\end{aligned}
$$

All derivatives have units g/(L·h).

### Numerical integration: Runge–Kutta

The system is written as

$$
\dot{\mathbf y}=f(t,\mathbf y,C_{S0}),\qquad
\mathbf y=\begin{bmatrix}C_S&C_X&C_P\end{bmatrix}^{T}.
$$

`scipy.integrate.solve_ivp` uses adaptive RK45 integration by default. For reference, one *fixed-step* classical fourth-order Runge–Kutta (RK4) update is

$$
\begin{aligned}
\mathbf k_1 &= f(t_n,\mathbf y_n)\\
\mathbf k_2 &= f(t_n+h/2,\mathbf y_n+h\mathbf k_1/2)\\
\mathbf k_3 &= f(t_n+h/2,\mathbf y_n+h\mathbf k_2/2)\\
\mathbf k_4 &= f(t_n+h,\mathbf y_n+h\mathbf k_3)\\
\mathbf y_{n+1}&=\mathbf y_n+\frac h6(\mathbf k_1+2\mathbf k_2+2\mathbf k_3+\mathbf k_4).
\end{aligned}
$$

> RK4 is shown conceptually; `solve_ivp(method="RK45")` uses an **adaptive embedded 5(4) pair**, not these exact fixed-step equations. `t_eval` specifies requested output times, not the internal integration step.

```python
import numpy as np
from scipy.integrate import solve_ivp

D = 1.0
Y_XS = 0.6
Y_PS = 0.05
mu_max = 2.0
Km = 0.5
Cs0_initial = 3.0

Cs_ss = D * Km / (Y_XS * mu_max - D)
Cx_ss = Y_XS * (Cs0_initial - Cs_ss)
Cp_ss = Y_PS * (Cs0_initial - Cs_ss)
y_ss = np.array([Cs_ss, Cx_ss, Cp_ss])


def nonlinear_ode(t, y, Cs0):
    Cs, Cx, Cp = y
    r = mu_max * Cs * Cx / (Km + Cs)
    return [
        D * (Cs0 - Cs) - r,
        -D * Cx + Y_XS * r,
        -D * Cp + Y_PS * r,
    ]


def simulate_nonlinear(Cs0_after=2.8, t_end=30):
    t = np.linspace(0, t_end, 1501)  # Output spacing = 0.02 h
    solution = solve_ivp(
        nonlinear_ode, (0, t_end), y_ss,
        args=(Cs0_after,), t_eval=t,
        method="RK45", rtol=1e-10, atol=1e-12,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return t, solution.y
```

`args=(Cs0_after,)` tells SciPy to evaluate `nonlinear_ode(t, y, Cs0_after)` at every solver call.

## Taylor linearization

Linearize the reaction rate $r(C_S,C_X)$ around the initial steady state. Introduce deviation variables

$$
\overline C_i=C_i-C_{i,ss},\qquad \overline C_{S0}=C_{S0}-C_{S0,initial}.
$$

The first-order expansion is

$$
\overline r\approx
\left.\frac{\partial r}{\partial C_S}\right|_{ss}\overline C_S+
\left.\frac{\partial r}{\partial C_X}\right|_{ss}\overline C_X
=k\overline C_S+b\overline C_X,
$$

where

$$
\boxed{
\begin{aligned}
k&=\frac{\mu_{\max}K_mC_{X,ss}}{(K_m+C_{S,ss})^2},\\
b&=\frac{\mu_{\max}C_{S,ss}}{K_m+C_{S,ss}}.
\end{aligned}}
$$

The linearized balances are

$$
\boxed{
\begin{aligned}
\frac{d\overline C_S}{dt}&=-(D+k)\overline C_S-b\overline C_X+D\overline C_{S0}\\
\frac{d\overline C_X}{dt}&=Y_{X/S}k\overline C_S+(-D+Y_{X/S}b)\overline C_X\\
\frac{d\overline C_P}{dt}&=Y_{P/S}k\overline C_S+Y_{P/S}b\overline C_X-D\overline C_P.
\end{aligned}}
$$

Only $C_{S0}$ changes in this scenario; hence $\overline D=0$.

## State-space model

Define the state vector and input:

$$
\overline{\mathbf x}=\begin{bmatrix}\overline C_S&\overline C_X&\overline C_P\end{bmatrix}^{T},
\qquad u=\overline C_{S0}.
$$

Then

$$
\boxed{\dot{\overline{\mathbf x}}=A\overline{\mathbf x}+Bu,\qquad
\mathbf y_{dev}=C\overline{\mathbf x}+Eu}
$$

where

$$
A=\begin{bmatrix}
-D-k&-b&0\\
Y_{X/S}k&-D+Y_{X/S}b&0\\
Y_{P/S}k&Y_{P/S}b&-D
\end{bmatrix},\qquad
B=\begin{bmatrix}D\\0\\0\end{bmatrix},
$$

$$
C=I_3,\qquad E=\mathbf 0_{3\times 1}.
$$

Numerically,

$$
\boxed{A=\begin{bmatrix}
-1.033333&-1.666667&0\\
0.020000&0&0\\
0.001667&0.083333&-1
\end{bmatrix},\quad
B=\begin{bmatrix}1\\0\\0\end{bmatrix}}
$$

| Matrix | Interpretation |
|---|---|
| $A$ | Coupling between concentration deviations and their time derivatives |
| $B$ | How the feed-substrate input enters the dynamics |
| $C=I_3$ | Outputs all three state deviations |
| $E=0$ | No instantaneous input-to-output feedthrough |

`E` is used for the feedthrough matrix, conventionally called `D` in control theory, to avoid confusion with the reactor dilution rate `D`.

```python
import control as ct

b = mu_max * Cs_ss / (Km + Cs_ss)
k = mu_max * Km * Cx_ss / (Km + Cs_ss)**2

A = np.array([
    [-D - k, -b, 0],
    [Y_XS * k, -D + Y_XS * b, 0],
    [Y_PS * k, Y_PS * b, -D],
])
B = np.array([[D], [0], [0]])
C = np.eye(3)
E = np.zeros((3, 1))

system = ct.ss(A, B, C, E)
```

## Linear step response

For a feed change from 3.0 to 2.8 g/L,

$$
\Delta C_{S0}=2.8-3.0=-0.2\ \mathrm{g/L}.
$$

`ct.step_response()` computes the response to a **unit step** input and zero initial state deviations. Scale that response by $\Delta C_{S0}$, then add the initial steady state back to obtain absolute concentrations.

```python
def simulate_linear(t, Cs0_after=2.8):
    delta_Cs0 = Cs0_after - Cs0_initial
    response = ct.step_response(system, T=t)

    # With one input and three outputs, squeeze=True gives (3, len(t)).
    y_unit_step = np.asarray(response.outputs)
    y_unit_step = y_unit_step.reshape(3, -1)

    y_deviation = delta_Cs0 * y_unit_step
    y_linear = y_ss[:, None] + y_deviation
    return y_linear
```

**Equivalent alternative:** Use `ct.forced_response(system, T=t, U=delta_Cs0 * np.ones_like(t))` to apply the input directly. Both functions simulate the **linear state-space model**, not `nonlinear_ode()`.

## Comparison and results
![part_e_bioreactor](https://hackmd.io/_uploads/BJZRKox2Mg.png)
![part_f_bioreactor](https://hackmd.io/_uploads/B1bRYoe3Gl.png)
# Chemostat Bioreactor: Nonlinear Dynamics and Linearization


### Values reported in the original notes

| Scenario | Model | $C_S$ (g/L) | $C_X$ (g/L) | $C_P$ (g/L) |
|---|---|---:|---:|---:|
| Part E (feed step to 2.8 g/L) | Nonlinear | 2.416201 | 0.230280 | 0.019190 |
| Part E (feed step to 2.8 g/L) | Linear | 2.423887 | 0.225668 | 0.018806 |
| Part F (input not specified in source notes) | Nonlinear | 0.489588 | 0.006316 | 0.000526 |
| Part F (input not specified in source notes) | Linear | 0.647019 | -0.088143 | -0.007345 |

The Part E results demonstrate reasonable agreement for a modest input change. The Part F results indicate that **linearization may become physically invalid for sufficiently large disturbances**, producing negative predicted biomass and product concentrations. Because the original notes do not specify the Part F input, that case is reported but not reproduced in the runnable example.

the input disturbance moves the reactor farther from that point.

## Physical Interpretation
We see that the linear model does a very good job of capturing the qualitative behavior in
part (e) as the tested perturbation was relatively small. It did a much poorer job in part
(f), and it even failed to account for the fact that biomass and product concentrations
cannot be negative. We can confidently use the linear transfer function models with small
process perturbations, but not large ones.
Physically, when we decrease the substrate concentration steady state, there are more cells
in the chemostat than the substrate concentration can support. This causes the cells to
die off. In part (e), when cells die off, the substrate is consumed at a lower rate, which
ends up increasing the concentration of substrate in the chemostat, which returns to a new
steady state. In part (f), there is not enough substrate in the feed to support any cells in
the bioreactor at all, so the cell biomass concentration approaches zero.
