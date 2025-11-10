# MS-NMPC_dFBA
Supporting codes of the paper: Multi-stage Economic Nonlinear Model Predictive Control of bioreactors using dynamic flux balance analysis models.

The paper describes a robust control approach based on the Multi-stage Economic Nonlinear Model Predictive Control, which allows handling with the degenerate solutions of the Flux Balance Analysis problem without being too conservative.

## Folders:
The files here are organized in folders by the different applications: 
1. Simulation
2. Dynamic Optimization
3. Closed-loop

## Scripts:
- main.jl: main file that defines all the parameters and initial conditions for the problem.
- pFBA_KKT_flux.jl: solves the optimization problem in JuMP.
- plotar.jl: creates figures in PDF from the output.

## Julia packages:
- JuMP 
- Ipopt
- HSL MA27
