# SIR Simulation with Vaccination Strategies

Simulation of an extended SIR epidemiological model with a vaccination compartment, comparing two intervention strategies: **preventive vaccination** (from the start of the outbreak) and **reactive vaccination** (triggered once infections cross a threshold).

## Description

This project extends the classic SIR (Susceptible–Infected–Recovered) model by adding a **Vaccinated (V)** compartment, and studies how the timing of vaccination affects the dynamics of an outbreak.

### Model

$$\frac{dS}{dt} = -\beta \frac{SI}{N} - \nu(t)\,S$$

$$\frac{dI}{dt} = \beta \frac{SI}{N} - \gamma I$$

$$\frac{dR}{dt} = \gamma I$$

$$\frac{dV}{dt} = \nu(t)\,S$$

Where:
- $\beta$ — transmission rate
- $\gamma$ — recovery rate
- $\nu(t)$ — vaccination rate, time-dependent according to the strategy

### Simulated strategies

| Strategy | Definition of $\nu(t)$ |
|---|---|
| **Preventive** | $\nu(t) = \nu_0$, constant from $t = 0$ |
| **Reactive** | $\nu(t) = 0$ until $I(t)$ exceeds a threshold $I_{\text{threshold}}$; then $\nu(t) = \nu_0$ |

