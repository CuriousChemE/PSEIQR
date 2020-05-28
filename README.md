# PSEIQR
SEIR model interactive graphic with adjustable parameters, quarantine, and vaccination

P = Not Susceptible (e.g. protected by vaccine)
S = Susceptible (not yet exposed),
E = Exposed (infected, but not yet able to infect others),
I = Infective (able to infect others), 
Q = Quarantined (removed from the Infective general population),
R = Recovered (no longer infective, either survived or deceased).

N (population size) set at 10000 = P + S + E + I + Q + R

dS/dT = -βSI

dE/dt = βSI – αE

dI/dT = αE – [(φγI) + (λ(1 - φ)I)

dQ/dt = φγI - δQ

dR/dt = δQ + λ(1 - φ)I

Solved via semi-implicit Euler method with time step = 0.1 days (see Hobbs, 2020).
