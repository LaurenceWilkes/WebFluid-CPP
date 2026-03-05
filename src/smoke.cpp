#include <SDL3/SDL.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <string>

using namespace std;

const int N = 164; // simulation size is N x N voxels
const int iterGS = 10; // number of iterations of the Gauss--Seidel linear solvers
const float visc = 0.00005f;        // diffusion for velocities of fluid
const float diff = 0.000001f;        // diffusion for dye
const int particleCount = 10000;

struct ScalarField {
    int side; // in our case side = N + 2
    vector<float> data;
    vector<float> tmp;

    ScalarField(int side) : side(side), data(side * side, 0), tmp(side * side, 0) {}

    void swapTemp() {data.swap(tmp);}

    inline int IX(int i, int j) const {return i + side * j;}
    float& operator()(int i, int j) {return data[IX(i, j)];}
    const float& operator()(int i, int j) const {return data[IX(i, j)];}
    float& tmpAt(int i, int j) {return tmp[IX(i, j)];}
    const float& tmpAt(int i, int j) const {return tmp[IX(i, j)];}

    void setBound(int b) {
        for (int i = 1; i <= side-2; i++) {
            data[IX(0, i)]      = b == 1 ? -data[IX(1, i)] : data[IX(1, i)];
            data[IX(side-1, i)] = b == 1 ? -data[IX(side-2, i)] : data[IX(side-2, i)];
            data[IX(i, 0)]      = b == 2 ? -data[IX(i, 1)] : data[IX(i, 1)];
            data[IX(i, side-1)] = b == 2 ? -data[IX(i, side-2)] : data[IX(i, side-2)];
        }
        data[IX(0, 0)]           = 0.5 * (data[IX(1, 0)] + data[IX(0, 1)]);
        data[IX(0, side-1)]      = 0.5 * (data[IX(1, side-1)] + data[IX(0, side-2)]);
        data[IX(side-1, 0)]      = 0.5 * (data[IX(side-2, 0)] + data[IX(side-1, 1)]);
        data[IX(side-1, side-1)] = 0.5 * (data[IX(side-2, side-1)] + data[IX(side-1, side-2)]);
    }

    void diffuse(int b, float diffCoeff, float dt) {
        float a = dt * diffCoeff * (side - 2) * (side - 2);
        for (int k = 0; k < iterGS; k++) {
            for (int j = 1; j <= side - 2; j++) {
                for (int i = 1; i <= side - 2; i++) {
                    int id = IX(i, j);
                    data[id] = (tmp[id] + a * (
                        data[IX(i - 1, j)] + data[IX(i + 1, j)] +
                        data[IX(i, j - 1)] + data[IX(i, j + 1)]
                    )) / (1 + 4 * a);
                }
            }
            setBound(b);
        }
    }

    void advect(int b, const ScalarField& u, const ScalarField& v, float dt) {
        float dt0 = dt * (side - 2);
	const float* udat;
	const float* vdat;
	if (b == 0) {
	    udat = u.data.data();
	    vdat = v.data.data();
	} else {
	    udat = u.tmp.data();  // we are required to read tmp when
	    vdat = v.tmp.data();  // updating u or v themselves
	}
        for (int j = 1; j <= side - 2; j++) {
            for (int i = 1; i <= side - 2; i++) {
                int id = IX(i, j);
                float f = i - dt0 * udat[id];
                float g = j - dt0 * vdat[id];
                f = clamp(f, 0.5f, side - 1.5f);
                g = clamp(g, 0.5f, side - 1.5f);
                int i0 = (int) f, j0 = (int) g;
                int i1 = i0 + 1, j1 = j0 + 1;
                float s1 = f - i0, s0 = 1 - s1;
                float t1 = g - j0, t0 = 1 - t1;
                data[id] = s0 * (t0 * tmp[IX(i0, j0)] + t1 * tmp[IX(i0, j1)])
			 + s1 * (t0 * tmp[IX(i1, j0)] + t1 * tmp[IX(i1, j1)]);
            }
        }
        setBound(b);
    } // advect
}; // ScalarField


struct FluidSim {
    ScalarField u, v;

    FluidSim() : u(N + 2), v(N + 2) {}

    void project() {
        float h = 1.0f / N;
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                v.tmpAt(i, j) = -0.5 * h * (
                    (u(i + 1, j) - u(i - 1, j)) +
                    (v(i, j + 1) - v(i, j - 1))
                );
                u.tmpAt(i, j) = 0;
            }
        }
        v.swapTemp(); v.setBound(0); v.swapTemp(); // v.tmp set bound
        u.swapTemp(); u.setBound(0); u.swapTemp(); // u.tmp set bound

        for (int k = 0; k < iterGS; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    u.tmpAt(i, j) = (
                        v.tmpAt(i, j) +
                        u.tmpAt(i - 1, j) + u.tmpAt(i + 1, j) +
                        u.tmpAt(i, j - 1) + u.tmpAt(i, j + 1)
                    ) / 4;
                }
            }
	    u.swapTemp(); u.setBound(0); u.swapTemp();
        }

        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                u(i, j) -= 0.5 * (u.tmpAt(i + 1, j) - u.tmpAt(i - 1, j)) / h;
                v(i, j) -= 0.5 * (u.tmpAt(i, j + 1) - u.tmpAt(i, j - 1)) / h;
            }
        }
	u.setBound(1);
	v.setBound(2);
    }

    void velocityStep(float dt) {
	u.swapTemp(); u.diffuse(1, visc, dt);
	v.swapTemp(); v.diffuse(2, visc, dt);
        project();
        u.swapTemp(); v.swapTemp();
	u.advect(1, u, v, dt); // interpolates from
	v.advect(2, u, v, dt); // tmp buffer as b != 0.
        project();
    }
}; // FluidSim

struct DyeSim {
    ScalarField dr, dg, db;

    DyeSim() : dr(N + 2), dg(N + 2), db(N + 2) {}

    void densityStep(const ScalarField& u, const ScalarField& v, float dt) {
        dr.swapTemp(); dr.diffuse(0, diff, dt); dr.swapTemp(); dr.advect(0, u, v, dt);
        dg.swapTemp(); dg.diffuse(0, diff, dt); dg.swapTemp(); dg.advect(0, u, v, dt);
        db.swapTemp(); db.diffuse(0, diff, dt); db.swapTemp(); db.advect(0, u, v, dt);
    } // densityStep
}; // DyeSim

struct Geyser {
    float next = 0.0f;
    bool par1 = true, par2 = true;
    int dist = 1;
    float col[3] = {1.0f, 1.0f, 1.0f};

    const int flowTime = 3000;
    const float pow = 1.0f;

    void resetGeyser(float now) {
	next = now + 500 + (rand() / (float) RAND_MAX) * 1500;
	par1 = rand() % 2; // random numbers not guaranteed to be uniform here
	par2 = rand() % 2;
	dist = rand() % N + 1;
	int i = rand() % 6;
	float f = (rand() / (float) RAND_MAX);
	float s = 0.85f;
	float pal[3] = {1 - s, 1 - s * f, 1 - s * (1 - f)};
	col[0] = pal[i % 3];
	col[1] = pal[i > 2 ? (i + 1) % 3 : (i + 2) % 3];
	col[2] = pal[i > 2 ? (i + 2) % 3 : (i + 1) % 3];
    } // resetGeyser
};

void updateGeyserDye(FluidSim& sim, DyeSim& dye, Geyser& g, float now) {
    if (now < g.next) return;
    if (now >= g.next + g.flowTime) {
	g.resetGeyser(now);
	return;
    }

    float vol = 155.0f;
    int nf = g.par2 ? 2 : N - 1;
    int sgn = g.par2 ? 1 : -1;
    float jet = g.pow * sgn;
    int jetSize = 3;

    for (int k = -jetSize; k <= jetSize; k++) {
	int dd = clamp(g.dist + k, 1, N);
	if (g.par1) {
	    dye.dr(dd, nf) += vol * g.col[0];
	    dye.dg(dd, nf) += vol * g.col[1];
	    dye.db(dd, nf) += vol * g.col[2];
	    sim.v(dd, nf) += jet;
	} else {
	    dye.dr(nf, dd) += vol * g.col[0];
	    dye.dg(nf, dd) += vol * g.col[1];
	    dye.db(nf, dd) += vol * g.col[2];
	    sim.u(nf, dd) += jet;
	}
    }
} // updateGeyserDye

void updateGeyser(FluidSim& sim, Geyser& g, float now) {
    if (now < g.next) return;
    if (now >= g.next + g.flowTime) {
	g.resetGeyser(now);
	return;
    }

    int nf = g.par2 ? 2 : N - 1;
    int sgn = g.par2 ? 1 : -1;
    float jet = g.pow * sgn;
    int jetSize = 3;

    for (int k = -jetSize; k <= jetSize; k++) {
	int dd = clamp(g.dist + k, 1, N);
	if (g.par1) sim.v(dd, nf) += jet;
	else        sim.u(nf, dd) += jet;
    }
} // updateGeyser


struct MouseState {
    float pmx = -1, pmy = -1;
    float ppmx = -1, ppmy = -1;
};

void pointerMotion(SDL_Event& e, float winW, float winH, FluidSim& sim, MouseState& m) {
    float gx = 1 + (e.motion.x / winW) * (N - 1);
    float gy = 1 + (e.motion.y / winH) * (N - 1);

    if (m.pmx < 0) {
        m.pmx = m.ppmx = gx;
        m.pmy = m.ppmy = gy;
        return;
    }

    float dx = (gx - m.ppmx) / 2;
    float dy = (gy - m.ppmy) / 2;
    m.ppmx = m.pmx; m.pmx = gx;
    m.ppmy = m.pmy; m.pmy = gy;

    int ii = clamp((int) gx, 1, N);
    int jj = clamp((int) gy, 1, N);
    int r = 2;

    int i0 = clamp(ii - r, 1, N), i1 = clamp(ii + r, 1, N);
    int j0 = clamp(jj - r, 1, N), j1 = clamp(jj + r, 1, N);

    for (int i = i0; i <= i1; i++) {
	for (int j = j0; j <= j1; j++) {
	    float dist2 = (ii - i) * (ii - i) + (jj - j) * (jj - j);
	    float wgt = exp(-2 * dist2 / r);
	    sim.u(i, j) += dx * wgt;
	    sim.v(i, j) += dy * wgt;
	}
    }
}

void initialise(SDL_Window*& win, SDL_Renderer*& ren, SDL_Texture*& tex) {
    SDL_Init(SDL_INIT_VIDEO);

    win = SDL_CreateWindow("Smoke", 1280, 800, SDL_WINDOW_RESIZABLE);
//    win = SDL_CreateWindow("Smoke", 1280, 800, SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIGH_PIXEL_DENSITY);

    ren = SDL_CreateRenderer(win, NULL);
//    SDL_SetRenderVSync(ren, 1);
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);

    tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, N, N);

//    SDL_SetTextureScaleMode(tex, SDL_SCALEMODE_PIXELART);
    SDL_SetTextureScaleMode(tex, SDL_SCALEMODE_LINEAR);
}

struct Particle {
    float x, y;
    float px, py;

    Particle() : x(0), y(0), px(0), py(0) {reset();}

    void reset() {
	x = px = 0.5f + N * (rand() / (float) RAND_MAX);
	y = py = 0.5f + N * (rand() / (float) RAND_MAX);
    } // reset

    void update(const ScalarField& u, const ScalarField& v, float dt0) {
        px = x; py = y;
	float f = clamp(x, 0.5f, N + 0.5f);
	float g = clamp(y, 0.5f, N + 0.5f);
	int i0 = (int) f, i1 = i0 + 1;
	int j0 = (int) g, j1 = j0 + 1;
	float s1 = f - i0, s0 = 1.0f - s1;
	float t1 = g - j0, t0 = 1.0f - t1;
	float ui = s0 * (t0 * u(i0, j0) + t1 * u(i0, j1))
		 + s1 * (t0 * u(i1, j0) + t1 * u(i1, j1));
	float vi = s0 * (t0 * v(i0, j0) + t1 * v(i0, j1))
		 + s1 * (t0 * v(i1, j0) + t1 * v(i1, j1));
	x += dt0 * ui;
	y += dt0 * vi;
    } // update
}; // Particle

struct ParticleSystem {
    vector<Particle> particles;

    const float resetPerc = 0.01;

    ParticleSystem(int size) : particles(size) {}

    void step(ScalarField& u, ScalarField& v, float dt) {
	float dt0 = dt * N;
	for (auto& p : particles) {
	    bool reset = rand() / (float) RAND_MAX < resetPerc;
	    if (reset) {
		p.reset();
	    } else {
	        p.update(u, v, dt0);
	    }
	}
    } // step
}; // ParticleSystem

void renderDye(SDL_Renderer* ren, SDL_Texture* tex, float winW, float winH, const DyeSim& dye, vector<unsigned char>& pixels) {
    for (int j = 1; j <= N; j++) {
	for (int i = 1; i <= N; i++) {
	    int k = ((j - 1) * N + (i - 1)) * 4;
	    pixels[k]     = clamp((int) dye.dr(i, j), 0, 255);
	    pixels[k + 1] = clamp((int) dye.dg(i, j), 0, 255);
	    pixels[k + 2] = clamp((int) dye.db(i, j), 0, 255);
	    pixels[k + 3] = 255;
	}
    }

    SDL_UpdateTexture(tex, NULL, pixels.data(), N * 4);

    SDL_FRect dst = {0, 0, winW, winH};

    SDL_RenderClear(ren);
    SDL_RenderTexture(ren, tex, NULL, &dst);
    SDL_RenderPresent(ren);
}

void renderParticles(SDL_Renderer* ren, float winW, float winH, const ParticleSystem& ps) {
    float sx = winW / N;
    float sy = winH / N;

    SDL_SetRenderDrawColor(ren, 0, 0, 0, 60);
    SDL_RenderFillRect(ren, NULL);

    SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);

    for (const auto& p : ps.particles) {
	SDL_RenderLine(ren,
	    p.px * sx,
	    p.py * sy,
	    p.x  * sx,
	    p.y  * sy);
    }

    SDL_RenderPresent(ren);
}

enum class RenderMode {
    Dye,
    Particles
}; // RenderMode

RenderMode mode = RenderMode::Dye;

int main(int argc, char** argv) {
    if (argc > 1 && string(argv[1]) == "--particles") {
	mode = RenderMode::Particles;
    }

    SDL_Window* win;
    SDL_Renderer* ren;
    SDL_Texture* tex;
    initialise(win, ren, tex);
    int w, h;
    SDL_GetWindowSize(win, &w, &h);
    float winW = (float) w;
    float winH = (float) h;

    FluidSim sim;
    DyeSim dye;
    vector<unsigned char> pixels(N * N * 4);

    MouseState mouse;
    Geyser geyser;
    geyser.resetGeyser((float) SDL_GetTicks());

    ParticleSystem particles(particleCount);

    float last = SDL_GetTicks() / 1000.0;
    bool running = true;

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
	    if (e.type == SDL_EVENT_WINDOW_RESIZED) {
		winW = (float) e.window.data1;
		winH = (float) e.window.data2;
	    }
            if (e.type == SDL_EVENT_MOUSE_MOTION)
                pointerMotion(e, winW, winH, sim, mouse);
        }

        float now = SDL_GetTicks() / 1000.0;
        float dt = min(now - last, 1.0f / 30);
        last = now;

	if (mode == RenderMode::Dye) {
	    updateGeyserDye(sim, dye, geyser, now * 1000);
	    sim.velocityStep(dt);
	    dye.densityStep(sim.u, sim.v, dt);
	    renderDye(ren, tex, winW, winH, dye, pixels);
	} else if (mode == RenderMode::Particles) {
	    updateGeyser(sim, geyser, now * 1000);
	    sim.velocityStep(dt);
	    particles.step(sim.u, sim.v, dt);
	    renderParticles(ren, winW, winH, particles);
	}

    }

    SDL_Quit();
}
