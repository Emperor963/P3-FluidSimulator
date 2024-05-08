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
        context.save();
            context.translate(this.x, this.y);
            
            context.beginPath();
            context.arc(0,0,2,0,2*Math.PI);
            context.fillStyle = '#0099aa';
            context.fill();
        context.restore();
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
    constructor(N, h, gx, gy){
        this.N = N;
        this.particles = Array(N);
        this.h = h;
        this.g = [0, -10];
        this.lattice = new Map();
        this.c = 2;
    }
    /**
     * Initialize the map and draw it on the canvas
     */
    initSim(){

        let length = Math.floor(Math.sqrt(this.N));
        let width = Math.floor(this.N/length);
        let d0 = 5;
        let mass = 1

        //Initiate the hashmap. We are assuming grid size to be one pixel but this can be tuned later
        for(let i = 1; i <= 39; i++){
            for(let j = 1; j < 59; j++){
                let particle = new Particle(10*j,10*i,0,0,d0,mass);
                let c = [particle.x , particle.y];
                let coords = c.join(",");
                this.lattice.set(coords,particle);
                this.particles.push(particle);
            }
        }

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
                let normalp6 = 315/(64*Math.PI*Math.pow(h,9));
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
                let normalp6 = 315/(64*Math.PI*Math.pow(h,9));
                let distp6 = Math.pow(h,2) - Math.pow(r,2);

                return -6*normalp6*r*Math.pow(distp6,2);

            //DOUBLE COSINE KERNEL
            case 3:
                let sin1 = sin(Math.PI* r/h);
                let sin2 = sin(Math.PI* 2*r/h);

                return -2*(Math.PI/h) * (sin2 + 2*sin1);
        }

    }


    /**
     * Property is calculated by using the smoothening kernel over the 8 nearest neighbors.
     * For some property A(x) = \sum_{neighbors j} A(x_j) * [m(x_j)/d(x_j)] * W(x_i, x_j);
     * We shall loop in clockwise direction, starting at 12 o' clock for reference when
     * looping over the neighbors.
     * @param {String} key : Property identifier 
     * @param {Particle} particle : an instance of a specific particle;
     */
    calculateProperty(key, particle){
        let kernel = 2; //Poly6 Kernel;
        let unitdist = 1;

        let property = 0;
        
        let k_Sides = kernel(unitdist, kernel);
        let k_Corners = kernel(Math.sqrt(2)*unitdist, kernel);

        for(let i = 0; i < 8; i++){
            let k = k_Corners;
            let nbrM = 0;
            let nbrD = particle.d;
            let nbr = this.lattice.get(JSON.stringify([particle.x + this.cwSteppper(i), particle.y+this.cwSteppper(i-2)]));
            if(i%2 == 0) k = k_Sides;
            if(nbr != null){
                nbrM = nbr.mass;
                nbrD = nbr.d;
            }
            property += k*nbr.get(key)*(nbrM/nbrD);
        }

        return property;
    }

    /**
     * Base formula: grad(A_i) = d_i * Sum_j [m_j * (A_i/d_i^2 + A_j/d_j^2) * grad(kernel)]; 
     * returns a 2-vector containing the X and Y components of the gradient of rho.
     * @param {Particle} particle : particle where gradrho needs to be calculated
     */
    calculateGradD(particle){
        let kernel = 2; //Poly6 kernel
        let unitdist = 1;

        let gradX = 0;
        let gradY = 0;

        let delK_Sides = this.kernelGradient(unitdist, kernel);
        let delK_Corners = this.kernelGradient(Math.sqrt(2)*unitdist, kernel);
        let theta = 0;

        for(let i = 0; i < 8; i++){
            let delK = delK_Corners;
            let nbrM = 0;
            let nbrD = particle.d;
            let nbr = this.lattice.get(JSON.stringify([particle.x + this.cwSteppper(i), particle.y+this.cwSteppper(i-2)]));
            if(i%2 == 0) delK = delK_Sides
            if(nbr != null){
                nbrM = nbr.mass;
                nbrD = nbr.d;
            }
            gradY += particle.d * (nbrM * delK * (1/nbrD + 1/particle.d)) * Math.cos(theta);
            gradX += particle.d * (nbrM * delK * (1/nbrD + 1/particle.d)) * Math.sin(theta);

            theta += Math.PI / 4;
        }
        return [gradX, gradY];
    }

    /**
     * Calculates the divergence of u at a specific point. Formula for divergence is as follows:
     * div(A) = Sum_j [m_j/d_j * dot(A_j, grad(kernel))];
     * @param {Particle} particle 
     */
    calculateDivU(particle){
        let kernel = 2 //Poly6 kernel
        let unitdist = 1;

        let div = 0;

        let delK_Sides = this.kernelGradient(unitdist, kernel);
        let delK_Corners = this.kernelGradient(Math.sqrt(2)*unitdist, kernel);
        let theta = 0;

        for(let i = 0; i < 8; i++){
            let delK = delK_Corners;
            let nbrM = 0;
            let nbrD = particle.d;
            let nbrUx = 0;
            let nbrUy = 0;
            let nbr = this.lattice.get(JSON.stringify([particle.x + this.cwSteppper(i), particle.y + this.cwSteppper(i-2)]));
            if(nbr != null){
                nbrM = nbr.mass;
                nbrD = nbr.d;
                nbrUx = nbr.ux;
                nbrUy = nbr.uy;
            }
            if(i%2 === 0) delK = delK_Sides;

            div += (nbrM/nbrD) * (nbrUx * Math.sin(theta) + nbrUy * Math.cos(theta));

            theta += Math.PI/4;
        }

        return div;

    }

    simulation(delT){
        
        this.particles.forEach(p => {
            let gradD = this.calculateGradD(p);
            let divU = this.calculateDivU(p);
            //console.log(this.c);
            p.ux += (-1*(Math.pow(this.c,2)/p.d)*gradD[0] + this.g[0])*delT;
            p.uy += (-1*(Math.pow(this.c,2)/p.d)*gradD[1] + this.g[1])*delT;

            p.d *= (1 + divU)*delT;

            p.x += p.ux*delT;
            p.y += p.uy*delT;
        })

        this.updateMap();
    }


    updateMap(){

    }

    /**
     * Custom step function for clockwise iteration over neighbors
     * @param {Number} iter 
     */
    cwSteppper(iter){
        if(iter < 0) iter = 7 + iter;

        if(iter === 0) return 0;
        else if(iter  - 4 > 0) return 1;
        else if(iter - 4 === 0) return 0;
        else if(iter - 4 < 0) return -1;
    }
}












//THIS CODE GOES AT THE VERY BOTTOM

let canvas =/**  @type {HTMLCanvasElement} */ (document.getElementById("canvas"));
let context = canvas.getContext("2d");

let sim = new Simulator(2400,3);
sim.initSim();
sim.simulation(0.1);
function draw(){
    context.clearRect(0, 0, canvas.width, canvas.height);
    //sim.particles.forEach(p => console.log(p.x+","+p.y));
    sim.lattice.forEach(p => p.draw(context));
    console.log("FRAME OVER");
}
function loop(timestamp){
    draw(context);
    window.requestAnimationFrame(loop);
}

window.requestAnimationFrame(loop);