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
     */
    constructor(x,y,ux,uy,d){
        this.x = x;
        this.y = y;
        this.ux = ux;
        this.uy = uy;
        this.d = d;

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
            context.arc(0,0,10,0,2*Math.PI);
            context.fillStyle = 'blue';
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
     */
    constructor(N, h){
        this.N = N;
        this.particles = Array(N);
        this.h = h;
    }

    /**
     * 
     * @param {number} qi     : coordinate value q of the ith particle
     * @param {number} qj     : coordinate value q of the jth particle
     * @param {number} choice : choice of kernel; 1 indicates gaussian, 2 indicates poly6
     * and 3 indicates double cosine
     */
    kernel(qi, qj, choice){
        let r = qi - qj;
        let h = this.h;

        switch(choice){

            //GAUSSIAN KERNEL
            case 1:
                let normalg = 1;
                return normalization * Math.exp(- Math.pow(r/h, 2));

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

}

//THIS CODE GOES AT THE VERY BOTTOM

let canvas =/**  @type {HTMLCanvasElement} */ (document.getElementById("canvas"));
let context = canvas.getContext("2d");

function loop(timestamp){


    window.requestAnimationFrame(loop);
}

window.requestAnimationFrame(loop);