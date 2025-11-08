# 2D Quantum Oscillator Viewer

Visual exploration of the isotropic two-dimensional quantum harmonic oscillator. A WebGL point cloud renders the Hermite–Gaussian eigenfunctions $\Psi_{n_x,n_y}(x,y)$ while a single slider enforces $n_x = n_y$ to give the familiar circularly degenerate ladder.

---

## Mathematical Background

### Hamiltonian and Eigenproblem

For the isotropic oscillator (mass $m$, frequency $\omega$) in two spatial dimensions,

$$
\hat{H} = \frac{\hat{p}_x^2 + \hat{p}_y^2}{2m} + \frac{1}{2}\ m\omega^2\left(\hat{x}^2 + \hat{y}^2\right).
$$

The separable eigenfunctions factor into 1D harmonic oscillator solutions,

$$
\Psi_{n_x,n_y}(x,y) = \psi_{n_x}(x)\ \psi_{n_y}(y),
$$

with

$$
\psi_n(x) = \frac{1}{\sqrt{2^n\ n!}}\ \left(\frac{m\omega}{\pi\hbar}\right)^{1/4}\ e^{-m\omega x^2 / 2\hbar}\ H_n\left(\sqrt{\frac{m\omega}{\hbar}}\ x\right),
$$

where $H_n(x)$ are the Hermite polynomials. The energy spectrum is

$$
E_{n_x,n_y} = \hbar\omega\left(n_x + n_y + 1\right),
$$

and each level with $N = n_x + n_y$ has degeneracy $g_N = N + 1$.

### What the Shader Evaluates

The Hermite polynomials are computed via the stable three-term recurrence relation

$$
\begin{aligned}
H_0(x) &= 1, \qquad H_1(x) = 2x, \\
H_{n+1}(x) &= 2xH_n(x) - 2nH_{n-1}(x)
\end{aligned}
$$

and combined with the Gaussian envelope to form the normalized eigenfunctions $\psi_n(x)$. The vertex shader evaluates this recurrence on the GPU to avoid factorial overflow. The rendered height is proportional to $\Psi_{n_x,n_y}(x,y)$; every frame we normalize by the maximum absolute amplitude so the surface always spans $z \in \left[-1,1\right]$.

### Color and Brightness

Each point is shaded with a combination of:

- Amplitude magnitude $|\psi|$ (controls intrinsic brightness);
- Lambert lighting with a fixed view-space light to reveal curvature;
- A cool $\to$ warm tint map driven by the normalized height $z$ so positive lobes glow warmer and negative lobes recede into blue shadows.

---
