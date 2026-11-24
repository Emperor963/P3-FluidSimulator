# P3-FluidSimulator



https://github.com/user-attachments/assets/5d7c7ec4-d15e-4154-9888-68b46290bd40


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
| `fluids.html` | Page structure and experiment module entry point |
| `fluids.js` | Exported particle and solver classes, spatial indexing, physics, and solver metrics |
| `experiments.js` | Parameters, centered initialization, initial vortex, renderer selection, stirring, animation, and profiling reports |
| `field-renderer.js` | Reusable sampled density and velocity fields, opacity, and smooth image rendering |
| `experimental-fluid.js` | Compatibility entry point importing `experiments.js` |
| `libs/custom.css` | Full-viewport canvas layout with scrolling disabled |
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

## Rendering and color switches

Two switches above the canvas change the display immediately.

| Switch | Behavior |
| --- | --- |
| **Particles / Field** | Draw individual circles or a continuous sampled field |
| **Velocity / Density** | Color by speed or density |

The initial renderer comes from `experiment.rendering.mode`, currently `"field"`. Velocity coloring is selected initially. Drag stirring works in either renderer.

### Particle rendering

Particle mode draws one circle per particle, with radius set by `rendering.radius`, currently 10 simulation pixels. The radius controls appearance. Neighbor selection uses particle centers and the support radius `h`.

Velocity coloring maps speed to an HSL hue,

$$
v_i=\sqrt{u_{x,i}^2+u_{y,i}^2},
\qquad
H_i=H_0\left(1-\min\left(\frac{v_i}{v_{\mathrm{scale}}},1\right)\right).
$$

Here `slowHue` supplies $H_0=210$ and `speedScale` supplies $v_{\mathrm{scale}}=20$. Stationary particles appear blue. Higher speeds progress through cyan, green, and yellow to red. Speeds at or above 20 share the red endpoint.

Density coloring uses the evolved particle density `p.d`. Its normalized color value is $\min(\max(\rho_i/(\rho_0 S),0),1)$, where $S$ is `densityScale`. With the current value of 1, reference density reaches red. Increasing `densityScale` keeps reference-density regions toward the blue end of the palette. Saturation and lightness are set by the corresponding rendering parameters, currently 90% and 60%.

### Field rendering

`FieldRenderer` builds a low-resolution image from overlapping particle contributions. At each image sample $\mathbf{x}$, it accumulates the compact Poly6 field,

$$
\rho_{\mathrm{visual}}(\mathbf{x})
=\sum_i m_i W(\lVert\mathbf{x}-\mathbf{x}_i\rVert,h).
$$

Each particle updates only samples within its support radius. In velocity-color mode, the renderer also accumulates weighted velocity and divides by local field density,

$$
\mathbf{u}_{\mathrm{visual}}(\mathbf{x})
=\frac{\sum_i m_i\mathbf{u}_i W_i}{\rho_{\mathrm{visual}}(\mathbf{x})}.
$$

The speed of this interpolated vector selects the color. Opposing particle velocities can cancel locally. In density-color mode, the renderer uses $\rho_{\mathrm{visual}}$ relative to reference density. Particle mode displays the evolved `p.d`, while field mode displays density reconstructed from current positions and masses. These two estimates can differ.

Opacity fades toward the free surface. For `opacityDensity` $T$, the normalized opacity input and final alpha are

$$
z=\min\left(\frac{\rho_{\mathrm{visual}}}{T\rho_0},1\right),
\qquad
\alpha=z^2(3-2z).
$$

The current $T=0.35$ makes samples fully opaque at 35% of reference density. Lowering this setting makes sparse regions more opaque. Raising it makes them more transparent. Samples with zero field density are transparent.

The renderer reuses typed arrays and image storage, caches a 256-color palette, and draws its offscreen image onto the main canvas with image smoothing enabled. The visual field is reconstructed each draw. Its values feed the displayed image only.

### Field quality and cost

Adjust `rendering.fieldSpacing` in `experiments.js` and reload. Smaller values produce finer samples and smoother spatial detail. The current value is 5.

| `fieldSpacing` | Image samples for a 1400 by 680 domain | Relative sample count |
| --- | --- | --- |
| 6 | 234 by 114 = 26,676 | About 0.70 times the current count |
| 5 | 280 by 136 = 38,080 | Current setting |
| 4 | 350 by 170 = 59,500 | About 1.56 times |
| 3 | 467 by 227 = 106,009 | About 2.78 times |

Halving the spacing produces approximately four times as many samples and increases particle-to-sample accumulation work. Increasing `N` changes the physical particle sampling and physics cost. Reducing `fieldSpacing` improves image sampling at the existing particle count.

A field can join nearby particle contributions into a smooth surface. Regions outside every particle's support remain empty. Large particle separations still affect the result.

`rendering.radius` applies to particle mode. Field mode uses `h` for contribution support and `opacityDensity` for surface visibility. Both modes share `speedScale`, `densityScale`, `slowHue`, `saturation`, and `lightness`.

## Interactive stirring

Hold the left mouse button and drag through the fluid to add velocity along the drag direction. Pointer events also support primary touch and pen input. Releasing the pointer ends the interaction. A stationary press applies no impulse.

The interaction affects particles within 80 simulation pixels of each pointer-move location. For pointer displacement $\Delta\mathbf{q}$ and distance $r_i$ from the current pointer location, the handler adds

$$
\Delta\mathbf{u}_i
=0.3\left(1-\frac{r_i}{80}\right)\Delta\mathbf{q},
\qquad 0\le r_i<80.
$$

Particles outside that radius receive no contribution. The effect fades linearly toward the edge. Density continues to evolve through the solver after the velocity change.

`canvasPoint()` converts browser pointer coordinates into canvas coordinates, accounting for CSS scaling and the canvas border. Pointer capture keeps the drag active when the pointer leaves the canvas. Pointer release, cancellation, or loss of capture clears the stored drag state. Touch scrolling is disabled over the canvas so dragging can stir the fluid.

Tune `experiment.stirring.radius` and `experiment.stirring.strength` in `experiments.js`, currently 80 and 0.3. Increasing the radius affects a larger region. Increasing strength transfers more velocity for the same pointer displacement. The handler scans the particle array for each pointer-move event.

This feature replaces the earlier outward click impulse. It applies user-controlled external momentum through pointer displacement, with no time normalization or swept-path interpolation. Rapid dragging can leave gaps between affected regions, and the result can vary with pointer-event frequency.

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

`initializeExperiment()` estimates a reference density $\rho_0$ by summing mass-weighted kernel values over a regular lattice with configured spacing $\Delta$, currently 8,

$$
\rho_0=m\sum_{a,b\in\mathbb{Z}}
W\!\left(\Delta\sqrt{a^2+b^2},h\right).
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

The coefficient `experiment.viscosity` is currently $\eta=0.1$. Setting it to zero disables this contribution. Equal velocities contribute zero damping. Positive weights exchange velocity between neighbors and oppose relative motion.

This is an experimental smoothing rule. Its strength increases with neighbor count, and its coefficient has no current calibration to physical dynamic or kinematic viscosity. It also damps rigid rotation. For equal masses, symmetric pair contributions conserve linear momentum at the acceleration level. Explicit integration and neighbor occupancy affect the resulting dynamics.

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
6. Clear the canvas and draw the selected particle or field renderer.

Velocity is advanced using the combined pressure, damping, and gravity acceleration,

```js
p.ux += (ax + Number(this.g[0])) * delT;
p.uy += (ay + Number(this.g[1])) * delT;
```

All derivatives are evaluated before any particle state is updated, so each step uses a shared state for its interaction calculations.

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

## Experiment configuration

All adjustable settings are in the exported `experiment` object at the top of `experiments.js`. Edit the values, save the file, and reload the page. The canvas fills the available viewport beneath the title, with page scrolling disabled. Its display scales to the window while simulation coordinates use `experiment.canvas.width` and `height`.

| Setting | Current value | Effect |
| --- | --- | --- |
| `N` | 9,000 | Exact number of particles to create |
| `columns` | 125 | Maximum particles across each row |
| `spacing` | 8 | Initial center-to-center distance and reference-density lattice spacing |
| `h` | 20 | Kernel support radius and spatial cell width |
| `mass` | 0.2 | Mass assigned to each particle |
| `gravity` | `[0, 1]` | Horizontal and vertical acceleration |
| `speedOfSound` | 20 | Pressure stiffness through its squared value |
| `viscosity` | 0.1 | Neighbor velocity-damping strength |
| `timeStep` | 0.1 | Simulation time advanced per animation frame |
| `densityFloorFraction` | 0.01 | Minimum density as a fraction of reference density |
| `boundaryMargin` | 10 | Distance between the reflecting walls and canvas edges |
| `canvas` | `{ width: 1400, height: 680 }` | Simulation coordinate dimensions |
| `vortex` | `{ sigma: 120, angularSpeed: 0 }` | Spread and strength of the initial rotation |
| `stirring` | `{ radius: 80, strength: 0.3 }` | Drag interaction size and strength |
| `rendering.mode` | `"field"` | Initial renderer, switchable on the page |
| `rendering.fieldSpacing` | 5 | Approximate distance between field samples |
| `rendering.opacityDensity` | 0.35 | Reference-density fraction producing full field opacity |
| `rendering.densityScale` | 1 | Reference-density multiplier reaching the red endpoint |
| `rendering.radius` | 10 | Particle-mode circle radius |
| `rendering.speedScale` | 20 | Speed reaching the red endpoint |
| `rendering.slowHue` | 210 | Blue hue at zero speed or density |
| `rendering.saturation` | 90 | HSL saturation percentage |
| `rendering.lightness` | 60 | HSL lightness percentage |
| `profiling` | `{ enabled: false, reportIntervalMs: 1000 }` | Console report visibility and interval |

### Centered particle initialization

A single particle-creation loop generates exactly `N` particles. Rows are derived from the requested count and column count,

$$
\text{rows}=\left\lceil\frac{N}{\text{columns}}\right\rceil.
$$

For counts smaller than `columns`, initialization uses only the required columns. The vertical extent is centered in the canvas, and each row is centered horizontally, including a partially filled final row. With the current configuration, 125 columns and 72 rows occupy centers from $x=204$ to $1196$ and $y=56$ to $624$.

Initialization rejects layouts that exceed the available area inside the boundary margins. Reducing `spacing` produces a tighter block and changes the reference density. Increasing it spreads the same particles across a larger area. `N = 0` creates an empty state.

`fluids.js` receives its numeric parameters through `new Simulator(canvas, experiment)`. `initializeExperiment()` prepares the particles and reference density, then passes them to `sim.initSim(particles, restDensity)`.

### Gaussian-envelope initial rotation

The initial vortex is centered at the canvas midpoint. With displacement $(d_x,d_y)$ from that point, the initial velocity is

$$
E=\exp\left(-\frac{d_x^2+d_y^2}{2\sigma^2}\right),
\qquad
u_x=-\omega d_y E,
\qquad
u_y=\omega d_x E.
$$

Here $\sigma$ is `vortex.sigma` and $\omega$ is `vortex.angularSpeed`. The symbols $u_x$ and $u_y$ denote velocity components. Positive angular speed produces clockwise rotation in downward-positive canvas coordinates. A negative value reverses it. Setting angular speed to zero gives zero initial velocity. Larger sigma spreads rotation over a wider region.

This field is applied once during initialization. Pressure, damping, gravity, stirring, and boundaries subsequently evolve the velocities.

### Density floor

`densityFloorFraction` sets the minimum allowed density relative to the calibrated reference density. At 0.01, the lower bound is $0.01\rho_0$. This protects divisions by density. Raising the fraction increases the lower bound. Frequent activation means the safeguard is materially influencing the density evolution.

## Performance profiling

Set `experiment.profiling.enabled = true`, reload, and open the browser developer tools Console to inspect `Fluid profile` reports. `reportIntervalMs` sets the approximate time between reports. At 1000, reports appear about once per second. Each report contains averages and maximums over its reporting interval.

**Profiling is currently disabled for console reporting only.** With `enabled: false`, timers, per-neighbor counters, occupancy tracking, and metric aggregation still execute. The flag suppresses the console reports but retains instrumentation overhead. A complete instrumentation bypass would require a further code change.

Neighbor counters execute once per candidate visit, potentially hundreds of thousands of times per frame. Timing, aggregation, and console output also add overhead. Use these measurements to identify bottlenecks and interpret measured timings with that cost in mind.

| Metric | Interpretation |
| --- | --- |
| Physics | Simulation duration excluding the measured map rebuild |
| Map rebuild | Clearing cell arrays, assigning membership, and updating occupancy metrics |
| Drawing calls | CPU time for the selected renderer, including field accumulation and image generation in field mode |
| Neighbor visits per step | Directed candidate visits before support filtering |
| Maximum cell occupancy | Largest cell population after a step |
| Occupied cells | Cells currently containing particles |
| Retained map entries | All allocated map keys, including empty cells |
| Frame interval | Time between animation callback timestamps |

The report also estimates frames per second. Its `elapsedSeconds` value uses the browser performance clock, which measures time relative to the page's time origin. Drawing timers can omit deferred browser rendering work. Instrumentation itself adds overhead.

Candidate pairs are currently evaluated in both directions. Increasing particle count or concentrating particles into fewer cells increases neighbor work. A rough cost model is $O(N+C+V)$ per step, where $C$ is the number of retained cells and $V$ is the number of directed candidate visits. A highly concentrated population can approach quadratic interaction cost.

A browser Performance recording can further distinguish arithmetic, generator traversal, allocation, garbage collection, and rendering costs.

## Rolling back the field-rendering trial

For an immediate visual comparison, use the **Particles / Field** switch. Set `rendering.mode` to `"particles"` and reload to start in particle mode.

The pre-field source snapshot is in `rollback-field-20261124-100459`. To restore the active files from that snapshot, run these commands from the project directory, then refresh the page.

```powershell
Copy-Item -LiteralPath .\rollback-field-20261124-100459\experiments.js -Destination .\experiments.js -Force
Copy-Item -LiteralPath .\rollback-field-20261124-100459\fluids.html -Destination .\fluids.html -Force
```

These commands overwrite the current experiment settings and HTML with their snapshot versions, including settings changed since the snapshot. The added `field-renderer.js` can remain as an unused file after restoration.

The snapshot also contains `manifest.json` and `Restore.ps1` for restoring all saved project files and removing the added renderer. The local Windows execution policy blocked the unsigned script during verification. The direct copy commands above restore the active entry point and page through ordinary file operations.

A later snapshot in `rollback-render-switch-20261124-100947` contains the field-rendering version immediately before adding the on-page renderer switch.

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

Moving particles carry advection through their trajectories. The implementation's pressure discretization, density floor, fixed step per frame, and viscosity approximation determine how closely the discrete behavior follows those continuum relations.

There is no incompressibility projection, surface tension, or explicit cohesion model. Density can drift from a kernel-summed estimate. Boundary sample deficits and near-overlap behavior also influence the result. Research use would require validation of conservation, symmetry, compressibility, resolution dependence, and benchmark flows alongside runtime measurements.
