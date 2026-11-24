/** Low-resolution, mass-weighted visual field built from particle centers. */
export class FieldRenderer {
    constructor(canvas, options, supportRadius) {
        if (!Number.isFinite(options.fieldSpacing) || options.fieldSpacing <= 0 ||
            !Number.isFinite(options.opacityDensity) || options.opacityDensity <= 0) {
            throw new RangeError("Field spacing and opacity density must be positive.");
        }
        this.canvas = canvas;
        this.options = options;
        this.h = supportRadius;
        this.width = Math.max(1, Math.ceil(canvas.width / options.fieldSpacing));
        this.height = Math.max(1, Math.ceil(canvas.height / options.fieldSpacing));
        this.stepX = canvas.width / this.width;
        this.stepY = canvas.height / this.height;
        const size = this.width * this.height;
        this.density = new Float32Array(size);
        this.momentumX = new Float32Array(size);
        this.momentumY = new Float32Array(size);
        this.surface = document.createElement("canvas");
        this.surface.width = this.width;
        this.surface.height = this.height;
        this.context = this.surface.getContext("2d");
        this.image = this.context.createImageData(this.width, this.height);
        this.palette = new Uint8ClampedArray(256 * 3);
        this.paletteKey = "";
    }

    updatePalette() {
        const { slowHue, saturation, lightness } = this.options;
        const key = `${slowHue},${saturation},${lightness}`;
        if (key === this.paletteKey) return;
        const light = lightness / 100;
        const amplitude = saturation / 100 * Math.min(light, 1 - light);
        for (let index = 0; index < 256; index++) {
            const hue = ((slowHue * (1 - index / 255)) % 360 + 360) % 360;
            // Standard HSL-to-RGB conversion, cached for all field pixels.
            for (const [channel, phase] of [[0, 0], [1, 8], [2, 4]]) {
                const k = (phase + hue / 30) % 12;
                const component = light - amplitude * Math.max(-1, Math.min(k - 3, 9 - k, 1));
                this.palette[index * 3 + channel] = 255 * component;
            }
        }
        this.paletteKey = key;
    }

    draw(sim, colorByDensity) {
        this.density.fill(0);
        if (!colorByDensity) {
            this.momentumX.fill(0);
            this.momentumY.fill(0);
        }
        const hSquared = this.h * this.h;
        const normalization = 4 / (Math.PI * this.h ** 8);
        for (const particle of sim.particles) {
            const minX = Math.max(0, Math.ceil((particle.x - this.h) / this.stepX - 0.5));
            const maxX = Math.min(this.width - 1, Math.floor((particle.x + this.h) / this.stepX - 0.5));
            const minY = Math.max(0, Math.ceil((particle.y - this.h) / this.stepY - 0.5));
            const maxY = Math.min(this.height - 1, Math.floor((particle.y + this.h) / this.stepY - 0.5));
            const massScale = particle.mass * normalization;
            for (let row = minY; row <= maxY; row++) {
                const dy = (row + 0.5) * this.stepY - particle.y;
                const rowOffset = row * this.width;
                for (let column = minX; column <= maxX; column++) {
                    const dx = (column + 0.5) * this.stepX - particle.x;
                    const remaining = hSquared - dx * dx - dy * dy;
                    if (remaining <= 0) continue;
                    const contribution = massScale * remaining * remaining * remaining;
                    const index = rowOffset + column;
                    this.density[index] += contribution;
                    if (!colorByDensity) {
                        this.momentumX[index] += contribution * particle.ux;
                        this.momentumY[index] += contribution * particle.uy;
                    }
                }
            }
        }

        this.updatePalette();
        const pixels = this.image.data;
        const opacityDensity = sim.restDensity * this.options.opacityDensity;
        for (let index = 0; index < this.density.length; index++) {
            const density = this.density[index];
            const pixel = index * 4;
            if (density === 0) {
                pixels[pixel + 3] = 0;
                continue;
            }
            const fraction = colorByDensity
                ? density / sim.restDensity / this.options.densityScale
                : Math.hypot(this.momentumX[index], this.momentumY[index]) /
                    density / this.options.speedScale;
            const color = Math.round(Math.min(Math.max(fraction, 0), 1) * 255) * 3;
            pixels[pixel] = this.palette[color];
            pixels[pixel + 1] = this.palette[color + 1];
            pixels[pixel + 2] = this.palette[color + 2];
            const opacity = Math.min(density / opacityDensity, 1);
            pixels[pixel + 3] = 255 * opacity * opacity * (3 - 2 * opacity);
        }
        this.context.putImageData(this.image, 0, 0);
        const context = this.canvas.getContext("2d");
        context.save();
        context.imageSmoothingEnabled = true;
        context.drawImage(this.surface, 0, 0, this.canvas.width, this.canvas.height);
        context.restore();
    }
}
