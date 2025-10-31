# 2D Quantum Oscillator Viewer

Visual exploration of the isotropic two‑dimensional quantum harmonic oscillator. A WebGL point cloud renders the Hermite–Gaussian eigenfunctions $\psi_{n_x,n_y}(x,y)$ while a single slider enforces $n_x = n_y$ to give the familiar circularly degenerate ladder.

---

## Mathematical background

### Hamiltonian and eigenproblem
For the isotropic oscillator (mass $m$, frequency $\omega$) in two spatial dimensions,
$$
\hat{H} = \frac{\hat{p}_x^2 + \hat{p}_y^2}{2m} + \frac{1}{2}m\omega^2\left(\hat{x}^2 + \hat{y}^2\right).
$$
The separable eigenfunctions factor into 1D harmonic oscillator solutions,
$$
\psi_{n_x,n_y}(x,y) = \phi_{n_x}(x)\,\phi_{n_y}(y),
$$
with
$$
\phi_n(x) = \frac{1}{\pi^{1/4} \sqrt{2^n n!}}\,H_n(x)\,e^{-x^2/2},
$$
where $H_n$ denotes the Hermite polynomials in the dimensionless coordinate $x = \sqrt{m\omega/\hbar}\,x_{\text{phys}}$ (same for $y$). The energy spectrum is
$$
E_{n_x,n_y} = (n_x + n_y + 1)\hbar\omega,
$$
and each level with $N = n_x + n_y$ has degeneracy $g_N = N + 1$.

### What the shader evaluates
The vertex shader implements the stable three-term recurrence
$$
\phi_{0}(x)=\pi^{-1/4}e^{-x^2/2},\qquad
\phi_{1}(x)=\sqrt{2}\,x\,\phi_{0}(x),\qquad
\phi_{k+1}(x)=\frac{\sqrt{2}\,x\,\phi_{k}(x)-\sqrt{k}\,\phi_{k-1}(x)}{\sqrt{k+1}}
$$
to avoid factorial overflow on the GPU. The rendered height is proportional to $\psi_{n_x,n_y}(x,y)$; every frame we normalise by the maximum absolute amplitude so the surface always spans $z \in [-1,1]$.

### Colour and brightness
Each point is shaded with a combination of:
- amplitude magnitude $|\psi|$ (controls intrinsic brightness);
- Lambert lighting with a fixed view-space light to reveal curvature;
- a cool $\rightarrow$ warm tint map driven by the normalised height $z$ so positive lobes glow warmer and negative lobes recede into blue shadows.
This emphasises nodal structure without overwhelming the monochrome point-cloud aesthetic.

---
