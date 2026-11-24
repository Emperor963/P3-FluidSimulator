/**
 * Primary Class for Fluid Simulation using Lattice based 
 * Smoothened Particle Hydrodynamics.
 */


//Particle class
export class Particle {
    /**
     * Constructor method for a particle.
     * @param {number} x   - x coordinate of position
     * @param {number} y   - y coordinate of position
     * @param {number} ux  - x coordinate of velocity
     * @param {number} uy  - y coordinate of velocity
     * @param {number} d   - local microscopic density
     * @param {number} mass- mass of the particle
     */
    constructor(x,y,ux,uy,d, mass){
        this.x = x;
        this.y = y;
        this.ux = ux;
        this.uy = uy;
        this.d = d;
        this.mass = mass;

    }

    get(key){
        if(key == 'x'){
            return this.x;
        }
        else if(key == 'y'){
            return this.y;
        }
        else if(key == 'ux'){
            return this.ux;
        }
        else if(key == 'uy'){
            return this.uy;
        }
        else if(key == 'd'){
            return this.d;
        }
    }

}

/**
 * Object oriented design of the simulator. This class should contain all relevant code pertaining to
 * Smoothening of microscopic variables, calculating gradients, time step solution to the modified NSE
 * Pseudo Collision Operator, etc.
 */
export class Simulator {

    /** Construct a solver from numeric parameters supplied by its caller. */
    constructor(canvas, options){
        this.canvas = canvas;
        this.h = options.h;
        this.cellSize = options.h;
        this.g = [...options.gravity];
        this.c = options.speedOfSound;
        this.viscosity = options.viscosity;
        this.densityFloorFraction = options.densityFloorFraction;
        this.boundaryMargin = options.boundaryMargin;
        this.particles = [];
        this.lattice = new Map();
        this.metrics = { mapMs: 0, neighborVisits: 0, maxCellOccupancy: 0, occupiedCells: 0 };
    }

    /** Install an initial state prepared by the caller. */
    initSim(particles, restDensity){
        if(!Number.isFinite(restDensity) || restDensity <= 0){
            throw new RangeError("Reference density must be positive and finite.");
        }
        this.particles = particles;
        this.N = particles.length;
        this.restDensity = restDensity;
        this.lattice.clear();
        this.updateMap();
    }
    /**
     * 
     * @param {number} qi     : coordinate value q of the ith particle
     * @param {number} qj     : coordinate value q of the jth particle
     * @param {number} choice : choice of kernel; 1 indicates gaussian, 2 indicates poly6
     * and 3 indicates double cosine
     */
    kernel(rj, choice){
        let r = Math.abs(rj);
        let h = this.h;

        switch(choice){

            //GAUSSIAN KERNEL
            //NOTE: DO NOT USE THE GAUSSIAN KERNEL! IT DOES NOT MAP THE NEIGHBORHOOD TIGHTLY ENOUGH
            //TO CAPTURE LOCAL PROPERTIES IN A ONE DEGREE LATTICE.
            case 1:
                let normalg = Math.sqrt(pi)*h;
                return normalg * Math.exp(- Math.pow(r/h, 2));

            //POLY6 KERNEL
            case 2:
                if(r >= h) return 0;
                let normalp6 = 4/(Math.PI*Math.pow(h,8));
                let distp6 = Math.pow(h,2) - Math.pow(r,2);

                return normalp6 * Math.pow(distp6, 3);

            //DOUBLE COSINE KERNEL
            case 3:
                let cos1 = Math.cos(Math.PI * r/h);
                let cos2 = Math.cos(2*Math.PI * r/h);

                return 4*cos1 + cos2 + 3;


            default:
                console.log("Select proper kernel!");
        }

    }

    /**
     * 
     * @param {number} rj :distance
     * @param {number} choice :choice of kernel 
     * @returns 
     */
    kernelGradient(rj, choice){
        let r = Math.abs(rj);
        let h = this.h;
        
        //We have to return the gradient along a specific direction. I think we can get away with a 
        //1D picture for this purpose.

        switch(choice){

            //GAUSSIAN KERNEL
            case 1:
                let normalg = 1/(Math.sqrt(pi)*h);
                return -2 * r * normalg * Math.exp(- Math.pow(r/h, 2)) * Math.pow(h,-2);
            
            //POLY6 KERNEL
            case 2:
                if(r >= h) return 0;
                let normalp6 = 4/(Math.PI*Math.pow(h,8));
                let distp6 = Math.pow(h,2) - Math.pow(r,2);

                return -6*normalp6*r*Math.pow(distp6,2);

            //DOUBLE COSINE KERNEL
            case 3:
                let sin1 = Math.sin(Math.PI* r/h);
                let sin2 = Math.sin(Math.PI* 2*r/h);

                return -2*(Math.PI/h) * (sin2 + 2*sin1);
        }

    }


    // Each particle queries its own fixed 3-by-3 local stencil.
    *localNeighbors(particle){
        const column = Math.floor(particle.x / this.cellSize);
        const row = Math.floor(particle.y / this.cellSize);
        for(let offsetY = -1; offsetY <= 1; offsetY++){
            for(let offsetX = -1; offsetX <= 1; offsetX++){
                const neighbors = this.lattice.get(`${column + offsetX},${row + offsetY}`);
                if(neighbors === undefined) continue;
                for(const neighbor of neighbors){
                    if(neighbor !== particle) yield neighbor;
                }
            }
        }
    }

    calculateLocalDerivatives(particle){
        let ax = 0;
        let ay = 0;
        let divergence = 0;
        const hSquared = this.h*this.h;
        const gradientScale = -24/(Math.PI*Math.pow(this.h, 8));
        const pressure = this.c*this.c*Math.max(0, particle.d - this.restDensity);
        const pressureWeight = pressure/(particle.d*particle.d);

        for(const neighbor of this.localNeighbors(particle)){
            this.metrics.neighborVisits++;
            const dx = particle.x - neighbor.x;
            const dy = particle.y - neighbor.y;
            const distanceSquared = dx*dx + dy*dy;
            if(distanceSquared >= hSquared) continue;

            const weight = (1 - distanceSquared / hSquared) ** 2;
            const damping = this.viscosity * weight;
            ax += damping * (neighbor.ux - particle.ux);
            ay += damping * (neighbor.uy - particle.uy);
            if(distanceSquared === 0) continue;

            // Poly6's radial derivative divided by distance cancels its r factor.
            const remaining = hSquared - distanceSquared;
            const gradientOverDistance = gradientScale*remaining*remaining;
            const neighborPressure = this.c*this.c*Math.max(0, neighbor.d - this.restDensity);
            const strength = -neighbor.mass*(
                pressureWeight + neighborPressure/(neighbor.d*neighbor.d)
            )*gradientOverDistance;
            ax += strength*dx;
            ay += strength*dy;
            divergence += (neighbor.mass/neighbor.d)*gradientOverDistance*
                ((neighbor.ux - particle.ux)*dx + (neighbor.uy - particle.uy)*dy);
        }
        return [ax, ay, divergence];
    }
    simulation(delT){
        this.metrics.neighborVisits = 0;
        
        // Evaluate all interactions before changing any particle's state.
        const derivatives = this.particles.map(p => this.calculateLocalDerivatives(p));
        this.particles.forEach((p, index) => {
            const [ax, ay, divU] = derivatives[index];

            
            //console.log(gradD[0] + "," + gradD[1]);
            p.ux += (ax + Number(this.g[0]))*delT;
            p.uy += (ay + Number(this.g[1]))*delT;

            const minimumDensity = this.restDensity * this.densityFloorFraction;
            const nextDensity = p.d * (1 - divU * delT);
            if (!Number.isFinite(nextDensity)) {
                throw new Error("Non-finite density during simulation");
            }
            p.d = Math.max(minimumDensity, nextDensity);

            p.x += p.ux * delT;
            p.y += p.uy * delT;
            if (![p.x, p.y, p.ux, p.uy].every(Number.isFinite)) {
                throw new Error("Non-finite position or velocity during simulation");
            }

            [p.x, p.ux] = this.reflectBoundary(p.x, p.ux, this.boundaryMargin, this.canvas.width - this.boundaryMargin);
            [p.y, p.uy] = this.reflectBoundary(p.y, p.uy, this.boundaryMargin, this.canvas.height - this.boundaryMargin);
        })

        const mapStart = performance.now();
        this.updateMap();
        this.metrics.mapMs = performance.now() - mapStart;
    }

    latticeKey(x, y){
        const column = Math.floor(x / this.cellSize);
        const row = Math.floor(y / this.cellSize);
        return `${column},${row}`;
    }

    updateMap(){
        // Reuse cell arrays and rebuild membership from all particles.
        for(const particles of this.lattice.values()){
            particles.length = 0;
        }
        this.metrics.maxCellOccupancy = 0;
        this.metrics.occupiedCells = 0;
        this.particles.forEach(particle => {
            const key = this.latticeKey(particle.x, particle.y);
            let particles = this.lattice.get(key);
            if(particles === undefined){
                particles = [];
                this.lattice.set(key, particles);
            }
            if(particles.length === 0) this.metrics.occupiedCells++;
            particles.push(particle);
            this.metrics.maxCellOccupancy = Math.max(this.metrics.maxCellOccupancy, particles.length);
        });
    }
    /**
     * Custom step function for clockwise iteration over neighbors
     * @param {Number} iter 
     */
    cwSteppper(iter){
        iter = ((iter % 8) + 8) % 8;

        if(iter === 0) return 0;
        else if(iter  - 4 > 0) return 1;
        else if(iter - 4 === 0) return 0;
        else if(iter - 4 < 0) return -1;
    }

    reflectBoundary(position, velocity, minimum, maximum){
        const width = maximum - minimum;
        const period = 2 * width;
        const offset = ((position - minimum) % period + period) % period;

        if (offset === 0) {
            return [minimum, Math.abs(velocity)];
        }
        if (offset === width) {
            return [maximum, -Math.abs(velocity)];
        }
        if (offset < width) {
            return [minimum + offset, velocity];
        }
        return [maximum - (offset - width), -velocity];
    }
}


