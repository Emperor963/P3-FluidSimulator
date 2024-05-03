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
     * @param {number} g : value of acceleration due to gravity
     */
    constructor(N, h){
        this.N = N;
        this.particles = Array(N);
        this.h = h;
        this.g = g;
        this.lattice = new Map(String, Particle);
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
     * Property is calculated by using the smoothening kernel over the 8 nearest neighbors.
     * For some property A(x) = \sum_{neighbors j} A(x_j) * [m(x_j)/d(x_j)] * W(x_i, x_j);
     * We shall loop in clockwise direction, starting at 12 o' clock for reference when
     * looping over the neighbors.
     * @param {String} key : Property identifier 
     * @param {Particle} particle : an instance of a specific particle;
     */
    CalculateProperty(key, particle){
        let kernel = 2; //Poly6 Kernel;
        let unitdist = 1;

        let property = 0;
        
        let k_Sides = kernel(unitdist, kernel);
        let k_Corners = kernel(Math.sqrt(2)*unitdist, kernel);

        //YES I AM AWARE THIS IS HORRIBLE CODE ETHIC I'M SORRY D':
        //I CAN'T THINK OF A BETTER WAY TO LOOP TONIGHT, WILL REFACTOR LATER
        //It is really a parity issue that can be fixed quite easily with some sort of parity operator


        
        //NORTH
        nbr_N = this.lattice.get(JSON.stringify([particle.x, particle.y-1]));
        property += k_Sides*nbr_N.get(key)*(nbr_N.mass/nbr_N.d);

        //NORTH EAST
        nbr_NE = this.lattice.get(JSON.stringify([particle.x+1, particle.y-1]));
        property += k_Corners*nbr_NE.get(key)*(nbr_NE.mass/nbr_NE.d);

        //EAST
        nbr_E = this.lattice.get(JSON.stringify([particle.x+1, particle.y]));
        property += k_Sides*nbr_E.get(key)*(nbr_E.mass/nbr_N.d);

        //SOUTH EAST
        nbr_SE = this.lattice.get(JSON.stringify([particle.x+1, particle.y+1]));
        property += k_Corners*nbr_SE.get(key)*(nbr_SE.mass/nbr_SE.d);

        //SOUTH
        nbr_S = this.lattice.get(JSON.stringify([particle.x, particle.y+1]));
        property += k_Sides*nbr_S.get(key)*(nbr_S.mass/nbr_S.d);

        //SOUTH WEST
        nbr_SW = this.lattice.get(JSON.stringify([particle.x-1, particle.y+1]));
        property += k_Corners*nbr_SW.get(key)*(nbr_SW.mass/nbr_SW.d);

        //WEST
        nbr_W = this.lattice.get(JSON.stringify([particle.x-1, particle.y]));
        property += k_Sides*nbr_W.get(key)*(nbr_W.mass/nbr_W.d);

        //NORTH WEST
        nbr_NW = this.lattice.get(JSON.stringify([particle.x-1, particle.y-1]));
        property += k_Corners*nbr_NW.get(key)*(nbr_NW.mass/nbr_NW.d);


        return property;
    }

}












//THIS CODE GOES AT THE VERY BOTTOM

let canvas =/**  @type {HTMLCanvasElement} */ (document.getElementById("canvas"));
let context = canvas.getContext("2d");

function loop(timestamp){


    window.requestAnimationFrame(loop);
}

window.requestAnimationFrame(loop);