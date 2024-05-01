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

//THIS CODE GOES AT THE VERY BOTTOM

let canvas =/**  @type {HTMLCanvasElement} */ (document.getElementById("canvas"));
let context = canvas.getContext("2d");

function loop(timestamp){


    window.requestAnimationFrame(loop);
}

window.requestAnimationFrame(loop);