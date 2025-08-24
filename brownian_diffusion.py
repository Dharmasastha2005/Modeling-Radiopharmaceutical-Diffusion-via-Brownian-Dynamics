# brownian_diffusion.py
# Minimal runnable Brownian diffusion simulator in heterogeneous tissue.
# Usage: python brownian_diffusion.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def brownian_diffusion_2d(grid_size=(100,100), cancer_region=((40,60),(40,60)), D_normal=1.0, D_cancer=0.3,
    n_particles=2000, n_steps=500, dt=1.0, rng_seed=42):
    rng = np.random.default_rng(rng_seed)
    nx, ny = grid_size
    D_map = np.full((nx, ny), D_normal)
    D_map[cancer_region[0][0]:cancer_region[0][1], cancer_region[1][0]:cancer_region[1][1]] = D_cancer
    positions = np.full((n_particles,2), [nx//2, ny//2], dtype=float)
    snapshots, msd = [], []
    for step in range(n_steps):
        if step % (n_steps//5) == 0 or step == n_steps-1:
            conc = np.zeros((nx,ny), dtype=int)
            for x,y in positions.astype(int):
                if 0 <= x < nx and 0 <= y < ny:
                    conc[x,y] += 1
            snapshots.append((step, conc))
        displacements = positions - np.array([nx//2, ny//2])
        msd.append(np.mean(np.sum(displacements**2, axis=1)))
        for i in range(n_particles):
            x,y = positions[i].astype(int)
            D = D_map[x,y] if (0 <= x < nx and 0 <= y < ny) else D_normal
            sigma = np.sqrt(2*D*dt)
            positions[i] += rng.normal(0, sigma, size=2)
            positions[i,0] = min(max(0, positions[i,0]), nx-1)
            positions[i,1] = min(max(0, positions[i,1]), ny-1)
    return snapshots, msd, D_map

if __name__ == '__main__':
    snapshots, msd, D_map = brownian_diffusion_2d()
    df_msd = pd.DataFrame({'step': np.arange(len(msd)), 'MSD': msd})
    df_msd.to_csv('diffusion_msd.csv', index=False)
    plt.figure(figsize=(12,3))
    for i,(step,conc) in enumerate(snapshots):
        plt.subplot(1,len(snapshots),i+1)
        plt.imshow(conc.T, origin='lower', cmap='hot', interpolation='nearest')
        plt.title(f'Step {step}')
        plt.axis('off')
    plt.tight_layout()
    plt.savefig('diffusion_snapshots.png', dpi=150)
    plt.close()
    plt.figure()
    plt.plot(np.arange(len(msd)), msd)
    plt.xlabel('Step')
    plt.ylabel('Mean Squared Displacement')
    plt.title('Brownian Diffusion MSD')
    plt.tight_layout()
    plt.savefig('diffusion_msd.png', dpi=150)
    print('Simulation complete: outputs saved.')
