/**
 * Primary Class for Fluid Simulation using Lattice based 
 * Smoothened Particle Hydrodynamics.
 */


//Particle class
class Particle {
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

    /**
     * Draws the particle
     * @param {CanvasRenderingContext2D} context 
     * 
     * TODO: Add functionality for changing color based on 
     * local density/velocity profiles
     */
    draw(context){
        /*context.save();
            context.translate(this.x, this.y);
            
            context.beginPath();
            context.arc(0,0,20,0,2*Math.PI);
            context.fillStyle = '#0099aa';
            context.fill();
        context.restore();*/
        // Create a radial gradient for the blur effect
        var gradient = context.createRadialGradient(1000, 1000, 1000, 1000, 100, 1000);
        gradient.addColorStop(0, 'rgba(100, 200, 250, 0.9)');
        gradient.addColorStop(1, 'rgba(100, 200, 250, 0.9)');

        // Fill a circle with the gradient
        context.fillStyle = gradient;
        context.beginPath();
        context.arc(this.x, this.y, 10, 0, Math.PI * 2);
        context.fill();
    }

}


/**
 * Object oriented design of the simulator. This class should contain all relevant code pertaining to
 * Smoothening of microscopic variables, calculating gradients, time step solution to the modified NSE
 * Pseudo Collision Operator, etc.
 */
class Simulator {

    /**
     * Constructor method for the simulator
     * @param {number} N : number of particles in the simulation
     * Note: This number can not be changed once simulation is created!
     * Recommended value of N for fast perfomance : ####
     * @param {number} h : Kernel distance for smoothening
     * @param {number} gx : value of acceleration due to gravity
     * @param {number} gy :value of acceleration due to gravity in y direction
     */
    constructor(N, h,canvas){
        this.N = N;
        this.particles = Array(N);
        this.h = h;
        this.gx =Number(document.getElementById("gx").value);
        this.gy = Number(document.getElementById("gy").value);
        this.g = [gx.value, gy.value];
        this.particleSpacing = 10;
        this.cellSize = h;
        this.lattice = new Map();
        this.c = 20;
        document.getElementById("c").value = String(this.c);
        this.canvas = canvas;
    }
    /**
     * Initialize the map and draw it on the canvas
     */
    initSim(){

        let length = Math.floor(Math.sqrt(this.N));
        let width = Math.floor(this.N/length);
        let d0 = 0;
        const massValue = Number(document.getElementById("mass").value);
        const reach = Math.ceil(this.h / this.particleSpacing);
        for(let row = -reach; row <= reach; row++){
            for(let column = -reach; column <= reach; column++){
                d0 += massValue * this.kernel(
                    Math.hypot(column, row)*this.particleSpacing, 2);
            }
        }
        this.restDensity = d0;
        let mass = massValue;

        // Initialize particles and index their centers by cell.
        for(let i = 1; i <= 50; i++){
            for(let j = 1; j < 30; j++){
                let particle = new Particle(10*j+250,10*i,0,0,d0,mass);
                this.particles.push(particle);
            }
        }
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

    calculatePressureAcceleration(particle){
        let ax = 0;
        let ay = 0;
        const pressure = this.c*this.c * Math.max(0, particle.d - this.restDensity);
        for(const neighbor of this.localNeighbors(particle)){
            const dx = particle.x - neighbor.x;
            const dy = particle.y - neighbor.y;
            const distance = Math.hypot(dx, dy);
            if(distance === 0 || distance >= this.h) continue;
            const neighborPressure = this.c*this.c * Math.max(0, neighbor.d - this.restDensity);
            const strength = -neighbor.mass * (
                pressure/(particle.d*particle.d) +
                neighborPressure/(neighbor.d*neighbor.d)
            ) * this.kernelGradient(distance, 2);
            ax += strength * dx/distance;
            ay += strength * dy/distance;
        }
        return [ax, ay];
    }

    calculateDivU(particle){
        let divergence = 0;
        for(const neighbor of this.localNeighbors(particle)){
            const dx = particle.x - neighbor.x;
            const dy = particle.y - neighbor.y;
            const distance = Math.hypot(dx, dy);
            if(distance === 0 || distance >= this.h) continue;
            divergence += (neighbor.mass/neighbor.d) * this.kernelGradient(distance, 2) *
                ((neighbor.ux - particle.ux)*dx + (neighbor.uy - particle.uy)*dy)/distance;
        }
        return divergence;
    }
    simulation(delT){
        
        this.particles.forEach(p => {
            let acceleration = this.calculatePressureAcceleration(p);
            let divU = this.calculateDivU(p);
            
            //console.log(gradD[0] + "," + gradD[1]);
            p.ux += (acceleration[0] + Number(this.g[0]))*delT;
            p.uy += (acceleration[1] + Number(this.g[1]))*delT;

            p.d  -= p.d * (divU) *delT;
            //console.log(p.d);

            p.x += p.ux*delT;
            p.y += this.stepfn(p.y, 480)*p.uy*delT;
            if(p.y  >= 400 + 50 || p.y + delT*p.uy <= 10){
                //console.log(p.y);
                p.uy = -0.5*p.uy;
            }
            if(p.x + delT*p.ux >= canvas.width + 20 || p.x + delT* p.ux <= 10){
                p.ux *= -0.5;
            }
            
        })

        this.updateMap();
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
        this.particles.forEach(particle => {
            const key = this.latticeKey(particle.x, particle.y);
            let particles = this.lattice.get(key);
            if(particles === undefined){
                particles = [];
                this.lattice.set(key, particles);
            }
            particles.push(particle);
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

    stepfn(t,x){
        if(x > t) return 1;

        if(x < t) return 0;
    }
}


//THIS CODE GOES AT THE VERY BOTTOM

let canvas =/**  @type {HTMLCanvasElement} */ (document.getElementById("canvas"));
let context = canvas.getContext("2d");

let sim = new Simulator(100,20,canvas);
sim.initSim();
//sim.simulation(0.1);
function draw(){
    
    context.clearRect(0, 0, canvas.width, canvas.height);
    //sim.lattice.forEach(p => console.log(p.ux+","+p.uy));
    sim.particles.forEach(p => p.draw(context));
    console.log("FRAME OVER");
}
function loop(timestamp){
    sim.simulation(0.1);
    draw(context);
    window.requestAnimationFrame(loop);
}

window.requestAnimationFrame(loop);
