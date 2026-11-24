# P3-FluidSimulator

<img width="426" height="240" alt="Video Project 2" src="https://github.com/user-attachments/assets/b49b884e-3796-4d89-99d1-abf0fccb11be" />



A browser-based 2D fluid simulation experiment, originally created as an extra project for CS 559 at the University of Wisconsin–Madison. JavaScript advances moving particles and renders them through the Canvas 2D API.

The project explores how local particle interactions can produce macroscopic fluid motion at low computational cost. Each particle queries a small local lattice, estimates pressure and velocity changes from nearby particles, and updates its own state. A shared spatial map makes those local queries practical.

## Run locally

Serve the repository with a local HTTP server and open `fluids.html`. For example, with Python available on your PATH, run the following from the project directory.

```powershell
python -m http.server 8000 --bind 127.0.0.1
```

Open [the simulator](http://127.0.0.1:8000/fluids.html). On installations with the Python launcher, `py` can replace `python`.

The animation starts automatically. Keep the server running, save edits, and refresh the page to load them. Use **Ctrl+F5** if the browser retains an earlier script. Stop the server with **Ctrl+C** in its terminal.

## Project layout

| File | Responsibility |
| --- | --- |
| `fluids.html` | Canvas, initial control values, and module entry point |
| `fluids.js` | Particle state, spatial indexing, physics, rendering, and profiling |
| `libs/custom.css` | Responsive canvas and control layout |
| `libs/runCanvas.js` | Reusable animation-control helper from the course library, currently unused by this simulator |

## Particle state and units

Each `Particle` stores the following quantities.

| Field | Symbol | Meaning |
| --- | --- | --- |
| `x`, `y` | $\mathbf{x}_i$ | Particle center in canvas coordinates |
| `ux`, `uy` | $\mathbf{u}_i$ | Velocity components |
| `d` | $\rho_i$ | Evolving density estimate |
| `mass` | $m_i$ | Particle mass |

Canvas coordinates increase to the right and downward, so positive `gy` accelerates particles toward the bottom of the screen. Lengths use pixels and time uses simulation units. Parameters currently form an experimental numerical model without a calibration to SI material properties.

Particles are rendered as translucent blue circles of radius 5. This visual radius is independent of the interaction radius. Neighbor selection uses particle centers.

## Local lattice and spatial lookup

Each particle sees a fixed **3 × 3 stencil of spatial cells**, including its central cell. The stencil is centered on the cell containing that particle and follows it as it moves across cell boundaries. Particle positions remain continuous.

The shared `lattice` map associates each occupied coordinate range with an array of particles. For cell width $s$, a center is assigned to

$$
c_x = \left\lfloor\frac{x}{s}\right\rfloor,
\qquad
c_y = \left\lfloor\frac{y}{s}\right\rfloor.
$$

`latticeKey()` joins those integer indices into a string. `localNeighbors()` queries the nine surrounding keys and yields every candidate except the particle itself. Several particles can occupy one cell.

The current cell width equals the kernel support radius, $s=h=20$. Therefore every particle within distance $h$ is contained in the nine queried cells. The derivative loop rejects candidates with squared distance at least $h^2$.

`updateMap()` clears the existing cell arrays and repopulates them from `particles` after every simulation step. It reuses map entries and their arrays. Empty entries remain allocated after particles leave them.

A fixed nine-cell stencil bounds the number of map lookups per particle. The number of particle interactions still grows with local occupancy. Dense clustering can make a step much more expensive.

## Kernel and reference density

The active interaction uses a compactly supported, two-dimensional Poly6 kernel,

$$
W(r,h)=\frac{4}{\pi h^8}(h^2-r^2)^3
$$

This expression applies for $0 \le r < h$. For $r \ge h$, the kernel is zero.

Its radial derivative inside support is

$$
W'(r,h)=-\frac{24}{\pi h^8}r(h^2-r^2)^2.
$$

For displacement $\mathbf{r}_{ij}=\mathbf{x}_i-\mathbf{x}_j$, the gradient with respect to particle $i$ is

$$
\nabla_i W_{ij}
=-\frac{24}{\pi h^8}(h^2-r_{ij}^2)^2\mathbf{r}_{ij}.
$$

`calculateLocalDerivatives()` uses this final expression directly. The factor of distance cancels analytically, allowing the loop to use squared distances throughout and avoid square roots.

`initSim()` estimates a reference density $\rho_0$ by summing mass-weighted kernel values over a regular lattice with spacing 10,

$$
\rho_0=m\sum_{a,b\in\mathbb{Z}}
W\!\left(10\sqrt{a^2+b^2},h\right).
$$

Only offsets inside kernel support contribute. The sum includes the central contribution. Every particle starts with this same density, including particles at the edges of the initial block. Subsequent density evolution uses velocity divergence rather than repeating this mass sum.

The file also contains Gaussian and double-cosine branches. The active solver uses Poly6. The Gaussian branches reference an undefined `pi` identifier and require correction before use.

## Pressure acceleration

Pressure is calculated from excess density,

$$
p_i=c^2\max(0,\rho_i-\rho_0),
$$

with the current speed-of-sound parameter $c=20$. Compression above reference density produces positive pressure. Below reference density, pressure is zero.

The pressure contribution to acceleration is

$$
\mathbf{a}_{i,p}
=-\sum_j m_j
\left(\frac{p_i}{\rho_i^2}+\frac{p_j}{\rho_j^2}\right)
\nabla_i W_{ij}.
$$

For positive pressure, the negative kernel derivative gives a separating response between particles. This supplies the model's interparticle repulsion. There is no separate geometric particle collision solver.

Exactly coincident particles receive velocity damping but skip pressure and divergence contributions. The Poly6 pressure gradient also approaches zero as separation approaches zero, limiting its ability to separate near-overlapping particles.

## Experimental viscosity

The model adds local damping based on differences in velocity from every neighbor within support,

$$
\mathbf{a}_{i,v}
=\eta\sum_j w_{ij}(\mathbf{u}_j-\mathbf{u}_i),
\qquad
w_{ij}=\left(1-\frac{r_{ij}^2}{h^2}\right)^2.
$$

The coefficient `this.viscosity` defaults to $\eta=0.1$. Setting it to zero disables this contribution. Equal velocities contribute zero damping. Positive weights exchange velocity between neighbors and oppose relative motion.

This is an experimental smoothing rule. Its strength increases with neighbor count, and its coefficient has no current calibration to physical dynamic or kinematic viscosity. It also damps rigid rotation. For equal masses, symmetric pair contributions conserve linear momentum at the acceleration level. Explicit integration and the velocity-increment wrapping described below affect the resulting dynamics.

## Density evolution

The divergence estimate uses relative velocities,

$$
D_i=\sum_j\frac{m_j}{\rho_j}
(\mathbf{u}_j-\mathbf{u}_i)\cdot\nabla_iW_{ij}.
$$

Uniform translation gives zero contribution. Density evolves through a discrete continuity update,

$$
\rho_i^{n+1}
=\max\left(0.01\rho_0,
\rho_i^n\left[1-\Delta t\,D_i\right]\right).
$$

Negative divergence increases density and positive divergence decreases it. The 1% reference-density floor prevents zero or negative denominators. Frequent activation indicates that the unconstrained update is leaving its useful operating range. A non-finite next density raises an error.

## Simulation step

Each animation callback advances one step with $\Delta t=0.1$ and draws the result.

1. Compute pressure acceleration, velocity damping, and divergence for every particle from the same current state.
2. Update each particle's velocity and density using the stored derivatives.
3. Advance position using the updated velocity.
4. Reflect particles at the canvas boundaries.
5. Rebuild spatial membership and update profiling counters.
6. Clear the canvas and render the particles.

The current velocity updates contain JavaScript remainder operations,

```js
p.ux += (ax + Number(this.g[0])) * delT % 10;
p.uy += (ay + Number(this.g[1])) * delT % 10;
```

These wrap each acceleration increment with remainder 10. They can produce discontinuities and do not impose a maximum accumulated speed. They are part of the current implementation and affect the interpretation of its physical behavior.

Position integration uses the updated velocity,

$$
\mathbf{x}_i^{n+1}
=\mathbf{x}_i^n+\Delta t\,\mathbf{u}_i^{n+1}.
$$

The step size is independent of the animation timestamp. Consequently, simulation time advances faster when more frames are rendered per real second.

## Boundaries

Particle centers reflect within a margin of 10 pixels from the canvas edges. `reflectBoundary()` folds an overshooting coordinate back into the domain, including multiple crossings during one update. It reverses the corresponding velocity component while preserving its magnitude.

For example, crossing the right wall reverses `ux` and preserves `uy`. This is an elastic wall response. There are no wall particles, wall friction, or explicit pressure corrections for missing neighbors near a wall.

Non-finite position or velocity raises an error before reflection. These checks stop the animation callback through an exception, leaving an error in the browser console.

## Current setup and controls

| Setting | Current implementation |
| --- | --- |
| Canvas coordinate dimensions | 1400 × 680 |
| Initial particles | 70 rows × 29 columns = 2,030 |
| Initial centers | $x=560\ldots840$, $y=10\ldots700$, in increments of 10 |
| Mass | HTML startup value `0.2` |
| Kernel support and cell width | 20 |
| Reference spacing | 10 |
| Speed-of-sound parameter | 20, assigned by the constructor |
| Experimental viscosity | 0.1 |
| Gravity | HTML startup values `gx = 0`, `gy = 1` |
| Step per animation frame | 0.1 |

CSS can scale the displayed canvas to fit the page. Its coordinate dimensions still determine simulation bounds.

The current initialization extends below the canvas. Those particles are folded inside by the first boundary update. The `N=100` constructor argument creates empty array slots, but the hard-coded initialization loops determine the actual population. Array iteration skips those initial empty slots. Changing `N` alone therefore does not change the number of particles.

The mass and gravity controls are sampled during startup. Changing them after startup has no live effect. For a different initial value, edit the HTML and reload. The constructor overrides the speed slider to 20. Change `this.c` in the constructor to alter that setting. The viscosity coefficient is also set in the constructor.

## Performance profiling

Open **Developer Tools → Console** to inspect a `Fluid profile` report approximately once per second. Reports contain averages and maximums over the reporting interval.

| Metric | Interpretation |
| --- | --- |
| Physics | Simulation duration excluding the measured map rebuild |
| Map rebuild | Clearing cell arrays, assigning membership, and updating occupancy metrics |
| Drawing calls | CPU time issuing canvas operations |
| Neighbor visits per step | Directed candidate visits before support filtering |
| Maximum cell occupancy | Largest cell population after a step |
| Occupied cells | Cells currently containing particles |
| Retained map entries | All allocated map keys, including empty cells |
| Frame interval | Time between animation callback timestamps |

The report also estimates frames per second. Its `elapsedSeconds` value uses the browser performance clock, which measures time relative to the page's time origin. Drawing timers can omit deferred browser rendering work. Instrumentation itself adds overhead.

Candidate pairs are currently evaluated in both directions. Increasing particle count or concentrating particles into fewer cells increases neighbor work. A rough cost model is $O(N+C+V)$ per step, where $C$ is the number of retained cells and $V$ is the number of directed candidate visits. A highly concentrated population can approach quadratic interaction cost.

A browser Performance recording can further distinguish arithmetic, generator traversal, allocation, garbage collection, and rendering costs.

## Physical interpretation and limitations

The design is inspired by particle discretizations of fluid equations. Its density evolution follows the material continuity equation,

$$
\frac{D\rho}{Dt}=-\rho\nabla\cdot\mathbf{u},
$$

and its momentum update combines pressure, experimental viscous damping, and gravity,

$$
\frac{D\mathbf{u}}{Dt}
=-\frac{1}{\rho}\nabla p
+\mathbf{a}_v+\mathbf{g}.
$$

Moving particles carry advection through their trajectories. The implementation's pressure discretization, density floor, remainder operations, fixed step per frame, and viscosity approximation determine how closely the discrete behavior follows those continuum relations.

There is no incompressibility projection, surface tension, or explicit cohesion model. Density can drift from a kernel-summed estimate. Boundary sample deficits and near-overlap behavior also influence the result. Research use would require validation of conservation, symmetry, compressibility, resolution dependence, and benchmark flows alongside runtime measurements.
