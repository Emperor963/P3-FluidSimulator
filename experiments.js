import { Particle, Simulator } from "./fluids.js";

// Edit experiment settings here. Reload the page after editing this object.
export const experiment = {
    N: 9000,
    columns: 125,
    spacing: 8,
    h: 20,
    mass: 0.2,
    gravity: [0, 0],
    speedOfSound: 20,
    viscosity: 0.1,
    timeStep: 0.1,
    densityFloorFraction: 0.01,
    boundaryMargin: 10,
    canvas: { width: 1400, height: 680 },
    vortex: { sigma: 120, angularSpeed: 0.0 },
    stirring: { radius: 80, strength: 0.3 },
    rendering: { radius: 10, speedScale: 20, slowHue: 210, saturation: 90, lightness: 60 },
    profiling: { enabled: false, reportIntervalMs: 1000 }
};

const canvas = document.getElementById("canvas");
canvas.width = experiment.canvas.width;
canvas.height = experiment.canvas.height;
const context = canvas.getContext("2d");
const sim = new Simulator(canvas, experiment);

function initializeExperiment(){
    const { N, columns, spacing, mass, boundaryMargin, vortex } = experiment;
    if(!Number.isInteger(N) || N < 0 || !Number.isInteger(columns) || columns < 1){
        throw new RangeError("N must be a nonnegative integer and columns a positive integer.");
    }
    const positive = [spacing, mass, experiment.h, experiment.speedOfSound,
        experiment.timeStep, experiment.densityFloorFraction, vortex.sigma];
    if(positive.some(value => !Number.isFinite(value) || value <= 0) ||
        !experiment.gravity.every(Number.isFinite) ||
        !Number.isFinite(experiment.viscosity) || experiment.viscosity < 0 ||
        !Number.isFinite(vortex.angularSpeed) ||
        !Number.isFinite(boundaryMargin) || boundaryMargin < 0){
        throw new RangeError("Experiment parameters must be finite and within their valid ranges.");
    }
    const actualColumns = Math.min(columns, N);
    const rows = N === 0 ? 0 : Math.ceil(N / actualColumns);
    const blockWidth = Math.max(0, actualColumns - 1)*spacing;
    const blockHeight = Math.max(0, rows - 1)*spacing;
    if(canvas.width <= 2*boundaryMargin || canvas.height <= 2*boundaryMargin ||
        blockWidth > canvas.width - 2*boundaryMargin ||
        blockHeight > canvas.height - 2*boundaryMargin){
        throw new RangeError("The particle layout exceeds the canvas. Adjust N, columns, spacing, or canvas size.");
    }

    // Calibrate reference density from a regular lattice at the chosen spacing.
    let restDensity = 0;
    const reach = Math.ceil(experiment.h / spacing);
    for(let row = -reach; row <= reach; row++){
        for(let column = -reach; column <= reach; column++){
            restDensity += mass*sim.kernel(Math.hypot(column, row)*spacing, 2);
        }
    }

    const particles = [];
    const startY = (canvas.height - blockHeight)/2;
    for(let index = 0; index < N; index++){
        const row = Math.floor(index / actualColumns);
        const column = index % actualColumns;
        const rowColumns = Math.min(actualColumns, N - row*actualColumns);
        const startX = (canvas.width - (rowColumns - 1)*spacing)/2;
        const x = startX + column*spacing;
        const y = startY + row*spacing;
        const dx = x - canvas.width/2;
        const dy = y - canvas.height/2;
        const envelope = Math.exp(-(dx*dx + dy*dy)/(2*vortex.sigma*vortex.sigma));
        const ux = -vortex.angularSpeed*dy*envelope;
        const uy = vortex.angularSpeed*dx*envelope;
        particles.push(new Particle(x, y, ux, uy, restDensity, mass));
    }
    sim.initSim(particles, restDensity);
}

initializeExperiment();

let previousPointer = null;

function canvasPoint(event) {
    const bounds = canvas.getBoundingClientRect();
    const scaleX = bounds.width / canvas.offsetWidth;
    const scaleY = bounds.height / canvas.offsetHeight;

    return {
        x: (event.clientX - bounds.left - canvas.clientLeft * scaleX) *
            canvas.width / (canvas.clientWidth * scaleX),
        y: (event.clientY - bounds.top - canvas.clientTop * scaleY) *
            canvas.height / (canvas.clientHeight * scaleY)
    };
}

canvas.style.touchAction = "none";

canvas.addEventListener("pointerdown", event => {
    if (!event.isPrimary || event.button !== 0) return;

    previousPointer = {
        ...canvasPoint(event),
        id: event.pointerId
    };
    canvas.setPointerCapture(event.pointerId);
});

canvas.addEventListener("pointermove", event => {
    if (!previousPointer || event.pointerId !== previousPointer.id) return;

    const point = canvasPoint(event);
    const dx = point.x - previousPointer.x;
    const dy = point.y - previousPointer.y;
    const radius = experiment.stirring.radius;
    const strength = experiment.stirring.strength;

    sim.particles.forEach(particle => {
        const distance = Math.hypot(
            particle.x - point.x,
            particle.y - point.y
        );
        if (distance >= radius) return;

        const weight = 1 - distance / radius;
        particle.ux += dx * strength * weight;
        particle.uy += dy * strength * weight;
    });

    previousPointer = { ...point, id: event.pointerId };
});

function stopStirring(event) {
    if (previousPointer?.id === event.pointerId) {
        previousPointer = null;
    }
}

canvas.addEventListener("pointerup", stopStirring);
canvas.addEventListener("pointercancel", stopStirring);
canvas.addEventListener("lostpointercapture", stopStirring);

//sim.simulation(experiment.timeStep);
function draw(){
    
    context.clearRect(0, 0, canvas.width, canvas.height);
    //sim.lattice.forEach(p => console.log(p.ux+","+p.uy));
    sim.particles.forEach(p => {
        const fraction = Math.min(Math.hypot(p.ux, p.uy) / experiment.rendering.speedScale, 1);
        const hue = experiment.rendering.slowHue * (1 - fraction);
        context.fillStyle = `hsl(${hue}, ${experiment.rendering.saturation}%, ${experiment.rendering.lightness}%)`;
        context.beginPath();
        context.arc(p.x, p.y, experiment.rendering.radius, 0, Math.PI * 2);
        context.fill();
    });

}
// Aggregate measurements to avoid logging on every animation frame.
const profile = {
    reportStart: performance.now(),
    previousTimestamp: undefined,
    frames: 0,
    totals: {},
    maximums: {},
    counts: {}
};

function recordMeasurement(name, value){
    profile.totals[name] = (profile.totals[name] || 0) + value;
    profile.maximums[name] = Math.max(profile.maximums[name] || 0, value);
    profile.counts[name] = (profile.counts[name] || 0) + 1;
}

function loop(timestamp){
    const start = performance.now();
    sim.simulation(experiment.timeStep);
    const afterSimulation = performance.now();
    draw();
    const afterDraw = performance.now();

    recordMeasurement("Physics (ms)", afterSimulation - start - sim.metrics.mapMs);
    recordMeasurement("Map rebuild (ms)", sim.metrics.mapMs);
    recordMeasurement("Drawing calls (ms)", afterDraw - afterSimulation);
    recordMeasurement("Neighbor visits per step", sim.metrics.neighborVisits);
    recordMeasurement("Maximum cell occupancy", sim.metrics.maxCellOccupancy);
    recordMeasurement("Occupied cells", sim.metrics.occupiedCells);
    recordMeasurement("Retained map entries", sim.lattice.size);
    if(profile.previousTimestamp !== undefined){
        recordMeasurement("Frame interval (ms)", timestamp - profile.previousTimestamp);
    }
    profile.previousTimestamp = timestamp;
    profile.frames++;

    if(experiment.profiling.enabled && afterDraw - profile.reportStart >= experiment.profiling.reportIntervalMs){
        const rows = {};
        for(const name of Object.keys(profile.totals)){
            rows[name] = {
                average: Number((profile.totals[name] / profile.counts[name]).toFixed(3)),
                maximum: Number(profile.maximums[name].toFixed(3))
            };
        }
        console.log("Fluid profile", {
            elapsedSeconds: Number((afterDraw / 1000).toFixed(1)),
            framesPerSecond: Number((profile.frames * 1000 /
                (afterDraw - profile.reportStart)).toFixed(1))
        });
        console.table(rows);
        profile.reportStart = afterDraw;
        profile.frames = 0;
        profile.totals = {};
        profile.maximums = {};
        profile.counts = {};
    }
    window.requestAnimationFrame(loop);
}

window.requestAnimationFrame(loop);