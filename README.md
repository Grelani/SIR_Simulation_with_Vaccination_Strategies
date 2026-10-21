## Methodology

The project follows seven steps: choose the model, observe it, find the theory behind it, write the mathematical model, choose a numerical method, implement it, and verify it with convergence tests.

---

### 1. Choosing the model

We use the **SIR** model (Susceptible–Infected–Recovered), extended with a **Vaccinated (V)** group.

**Why this model:**
- It is the basic model for diseases where recovery gives lasting immunity.
- It has only three parameters ($\beta$, $\gamma$, $\nu$), so each one is easy to interpret and change.
- Adding V lets us ask a concrete question: **does the timing of vaccination matter?** We compare a preventive and a reactive strategy.

---

### 2. Observing the model

First, we look at how the plain SIR model behaves with no vaccination ($\nu = 0$), using $N = 1000$, $\beta = 0.3$, $\gamma = 0.1$ ($R_0 = 3$) and one initial infected person:

- $S(t)$ only goes down and $R(t)$ only goes up.
- $I(t)$ rises, reaches a single peak, and then falls to almost zero.
- The peak is about **301 infected people on day 38**, and about **941 out of 1000** people get infected at some point.
- The outbreak ends because there are no susceptible people left to infect, not because everyone got sick: about 59 people are never infected.

---

### 3. Theory

**Compartmental models.** The population is split into groups (compartments), and the model describes how people move from one group to another. The SIR model was introduced by Kermack and McKendrick (1927).

**Assumptions:**
- Closed population of constant size $N$ (no births, deaths or migration).
- Everyone mixes equally: any person can infect any other.
- New infections are proportional to $S \cdot I / N$.
- Recovered people stay immune.
- The vaccine works immediately and perfectly: a vaccinated person can no longer be infected.

**Basic reproduction number.**

$$R_0 = \frac{\beta}{\gamma}$$

It is the average number of people one infected person infects in a fully susceptible population. If $R_0 > 1$ there is an outbreak; if $R_0 < 1$ the infection dies out.

**Herd immunity.** The outbreak stops growing once the immune fraction of the population is above

$$p_c = 1 - \frac{1}{R_0}$$

With $R_0 = 3$, $p_c \approx 0.67$.

**Known results for SIR without vaccination** (we use them to check the simulation):
- The peak happens when $S = N / R_0$, and its height is

$$I_{\max} = N\left[1 - \frac{1 + \ln R_0}{R_0}\right]$$

- Final size: the number of people who never get infected, $S_\infty$, satisfies $\ln\left(S_\infty / S_0\right) = R_0\,(S_\infty - S_0)/N$ (approximately, since $S_0 \approx N$).

**Vaccination strategies:**
- **Preventive:** vaccinate from day 0, before the infection spreads.
- **Reactive:** start vaccinating only when the number of infected people passes a threshold.

**References**
- Kermack, W. O., McKendrick, A. G. (1927). *A contribution to the mathematical theory of epidemics*. Proc. R. Soc. A, 115, 700–721.
- Anderson, R. M., May, R. M. (1991). *Infectious Diseases of Humans: Dynamics and Control*. Oxford University Press.
- Hethcote, H. W. (2000). *The mathematics of infectious diseases*. SIAM Review, 42(4), 599–653.

---

### 4. Mathematical model

Four nonlinear ordinary differential equations:

$$\frac{dS}{dt} = -\beta \frac{S I}{N} - \nu(t)\, S$$

$$\frac{dI}{dt} = \beta \frac{S I}{N} - \gamma I$$

$$\frac{dR}{dt} = \gamma I$$

$$\frac{dV}{dt} = \nu(t)\, S$$

Each equation says how a group changes: what comes in minus what goes out.

**Variables:** $S, I, R, V$ are the number of people in each group, and $S + I + R + V = N$.

**Parameters:**

| Symbol | Meaning | Base value |
|---|---|---|
| $N$ | Total population | 1000 |
| $\beta$ | Transmission rate | 0.3 |
| $\gamma$ | Recovery rate ($1/\gamma$ = average length of infection) | 0.1 |
| $\nu_0$ | Vaccination rate when vaccination is on | to be defined |
| $I_{\text{threshold}}$ | Number of infected that triggers reactive vaccination | to be defined |

**Initial conditions:** $S(0) = N - 1$, $I(0) = 1$, $R(0) = V(0) = 0$.

**Strategies** (they define $\nu(t)$):
- Preventive: $\nu(t) = \nu_0$ for all $t \ge 0$.
- Reactive: $\nu(t) = 0$ while $I(t) < I_{\text{threshold}}$, and $\nu(t) = \nu_0$ once $I(t)$ reaches it.

**Conservation.** If we add the four equations, every term cancels out:

$$\frac{dS}{dt} + \frac{dI}{dt} + \frac{dR}{dt} + \frac{dV}{dt} = 0 \;\Rightarrow\; S + I + R + V = N$$

So nobody is created or lost, people only move between groups.

---

### 5. Numerical method

The model is an **initial value problem** with no closed-form solution (because of the nonlinear $SI$ term), so we solve it numerically.

**Chosen method: Runge–Kutta 4(5) by Dormand–Prince (RK45), with adaptive step size**, using `scipy.integrate.solve_ivp`.

**Why:**
- The system is **not stiff**: $\beta$ and $\gamma$ have similar time scales, so we do not need an implicit method.
- The adaptive step controls the error through the tolerances `rtol` and `atol`, so we do not pick the step by hand.
- It supports **event detection**, which we need for the reactive strategy.

**Comparison methods** (coded by hand for the convergence tests):
- Explicit Euler, order 1: $y_{n+1} = y_n + h\,f(t_n, y_n)$.
- Classical Runge–Kutta of order 4 (RK4), fixed step.

**Reactive strategy.** The switch of $\nu(t)$ at the threshold makes the equations discontinuous. To keep the solution accurate:
1. Integrate with $\nu = 0$ and an event $g(t, y) = I - I_{\text{threshold}}$, which stops the integration when $I$ crosses the threshold going up.
2. Use the state at that moment as the new initial condition.
3. Continue with $\nu = \nu_0$ until the final time $T$.

**Reference:** Dormand, J. R., Prince, P. J. (1980). *A family of embedded Runge–Kutta formulae*. J. Comput. Appl. Math., 6(1), 19–26.

---

### 6. Implementation

Python with NumPy, SciPy and Matplotlib. The code is in the notebook `sir_vacunacion.ipynb`.

Main parts:
- `sirv(t, y, beta, gamma, nu, N)`: the right-hand side of the equations.
- `simular(...)`: solves the system with `solve_ivp` (RK45, `rtol = 1e-8`, `atol = 1e-10`).
- Preventive strategy: constant `nu`.
- Reactive strategy: solve in two pieces using event detection (pending).
- `metricas(sol)`: peak of infected, day of the peak, total infected, total vaccinated.

---

### 7. Convergence tests and verification

Goal: check that the numerical solution respects the properties of the model and that it converges.

| Test | What we check | Status |
|---|---|---|
| Population conservation | $S + I + R + V = N$ at all times | Done |
| Reduction to classic SIR | With $\nu = 0$ it matches the 3-equation SIR and $V = 0$ | Done |
| Peak height | Simulated $I_{\max}$ matches the theoretical formula | Pending |
| Final size | Simulated $S_\infty$ satisfies the final size relation | Pending |
| Order of convergence | Error vs. step size $h$ (against a high-accuracy reference): slope $\approx 1$ for Euler and $\approx 4$ for RK4 | Pending |
| Tolerance sensitivity | RK45 results do not change when `rtol` is reduced | Pending |
| Reactive event | The time the threshold is crossed does not depend on the tolerance | Pending |
