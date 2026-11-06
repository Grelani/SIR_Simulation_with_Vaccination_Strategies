# SIR Simulation with Vaccination Strategies

## Objective

Simulate the spread of an infectious disease and compare two vaccination strategies:

- **Preventive:** vaccination starts on day 0
- **Reactive:** vaccination starts when infections reach a certain threshold

The goal is to determine how vaccination timing affects the epidemic.

## Methodology

### 1. Model Selection

We use the **SIRV model**, which divides the population into four groups:

- **S:** Susceptible
- **I:** Infected
- **R:** Recovered
- **V:** Vaccinated

This model allows us to study how infections and vaccination change over time.

### 2. Model Observation

We begin with a population of 1,000 people, including one infected individual.

Without vaccination, infections increase, reach a peak, and eventually decrease.

### 3. Theoretical Background

The SIRV model uses differential equations to describe the spread of an infection and the effects of vaccination.

We assume a closed population, lasting immunity after recovery, and fully effective vaccination.

The model is based on the epidemiological approaches discussed by Marinov and Marinova (2022) and Lin et al. (2024).

### 4. Mathematical Model

The system is described by:

$$
\frac{dS}{dt}=-\beta\frac{SI}{N}-\nu(t)S
$$

$$
\frac{dI}{dt}=\beta\frac{SI}{N}-\gamma I
$$

$$
\frac{dR}{dt}=\gamma I
$$

$$
\frac{dV}{dt}=\nu(t)S
$$

Where:

- $\beta = 0.3$: transmission rate
- $\gamma = 0.1$: recovery rate
- $\nu(t)$: vaccination rate
- $N = 1000$: total population

The initial conditions are:

- $S(0) = 999$
- $I(0) = 1$
- $R(0) = V(0) = 0$

### 5. Numerical Method

We will use the **Runge-Kutta 4(5) method (RK45)** with **solve_ivp** from SciPy

RK45 automatically adjusts the step size to control numerical error, allowing us to approximate the solutions of the differential equations.
