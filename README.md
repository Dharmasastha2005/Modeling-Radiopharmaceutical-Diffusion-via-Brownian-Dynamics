# Modeling-Radiopharmaceutical-Diffusion-via-Brownian-Dynamics
Simulate how radiopharmaceuticals diffuse in normal vs. cancerous tissue environments.

Background
Drug molecules undergo Brownian motion (random walk).
Cancer tissue can be modeled with altered diffusion coefficient (D) or barriers.
Brownian motion is governed by stochastic differential equation:

dx= sqrt(2DΔt)⋅N(0,1)
Equations

Mean squared displacement (MSD):
(r^2(t)⟩=2nDt
(n = dimension: 2D or 3D).

Step Plan
Define environment (grid of normal + cancerous zones).
Assign diffusion coefficient D_normal, D_cancer.
Initialize N particles at origin (source).
Update positions with Brownian step each Δt.
Apply boundaries (absorbing = sink, reflecting = barrier).
Collect stats: concentration map, MSD, penetration depth.

Coding Outline (Python)
Inputs: Grid size, diffusion coefficients, particle count, time steps.
Loop: Move particles according to diffusion coefficients of local tissue.
Outputs: Heatmap of concentration vs. time, MSD curves, comparison normal vs. cancer.

Contents
brownian_diffusion.py → main code.
README: Biological relevance (drug penetration into tumors).
Figures: diffusion maps, MSD plots.
