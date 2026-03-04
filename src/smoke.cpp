#include <SDL3/SDL.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const int N = 164; // simulation size is N x N voxels
const int iterGS = 10; // number of iterations of the Gauss--Seidel linear solvers

inline int IX(int i, int j) { return i + (N + 2) * j; }

struct FluidSim {
    const int size = (N + 2) * (N + 2); // Size of vectors containing velocities/densities. (includes a buffer of one pixel around whole grid)
    const float diff = 0.000001f;       // diffusion for dyes
    const float visc = 0.00005f;        // diffusion for velocities of fluid

    vector<float> u, v, u0, v0;
    vector<float> dr, dg, db, dr0, dg0, db0;

    FluidSim() :
        u(size, 0), v(size, 0), u0(size, 0), v0(size, 0),
        dr(size, 0), dr0(size, 0), dg(size, 0), dg0(size, 0), db(size, 0), db0(size, 0) {}

    void setBound(int b, vector<float>& x) {
        for (int i = 1; i <= N; i++) {
            x[IX(0, i)]     = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)]     = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
        }
        x[IX(0, 0)]         = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)]     = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)]     = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }

    void diffuse(int b, vector<float>& x, vector<float>& x0, float diffCoeff, float dt) {
        float a = dt * diffCoeff * N * N;
        for (int k = 0; k < iterGS; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    int id = IX(i, j);
                    x[id] = (x0[id] + a * (
                        x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                        x[IX(i, j - 1)] + x[IX(i, j + 1)]
                    )) / (1 + 4 * a);
                }
            }
            setBound(b, x);
        }
    }

    void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& u, vector<float>& v, float dt) {
        float dt0 = dt * N;
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                int id = IX(i, j);
                float f = i - dt0 * u[id];
                float g = j - dt0 * v[id];
                f = clamp(f, 0.5f, N + 0.5f);
                g = clamp(g, 0.5f, N + 0.5f);
                int i0 = floor(f), j0 = floor(g);
                int i1 = i0 + 1, j1 = j0 + 1;
                float s1 = f - i0, s0 = 1 - s1;
                float t1 = g - j0, t0 = 1 - t1;
                d[id] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)])
                      + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        setBound(b, d);
    } 

    void project(vector<float>& u, vector<float>& v, vector<float>& p, vector<float>& div) {
        float h = 1.0f / N;
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                div[IX(i, j)] = -0.5 * h * (
                    (u[IX(i + 1, j)] - u[IX(i - 1, j)]) +
                    (v[IX(i, j + 1)] - v[IX(i, j - 1)])
                );
                p[IX(i, j)] = 0;
            }
        }
        setBound(0, div);
        setBound(0, p);

        for (int k = 0; k < iterGS; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    p[IX(i, j)] = (
                        div[IX(i, j)] +
                        p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                        p[IX(i, j - 1)] + p[IX(i, j + 1)]
                    ) / 4;
                }
            }
            setBound(0, p);
        }

        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
                v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
            }
        }
        setBound(1, u);
        setBound(2, v);
    }

    void velocityStep(float dt) {
        u.swap(u0); diffuse(1, u, u0, visc, dt);
        v.swap(v0); diffuse(2, v, v0, visc, dt);
        project(u, v, u0, v0);
        u.swap(u0); v.swap(v0);
        advect(1, u, u0, u0, v0, dt);
        advect(2, v, v0, u0, v0, dt);
        project(u, v, u0, v0);
    }

    void densityStep(float dt) {
        dr.swap(dr0); diffuse(0, dr, dr0, diff, dt); dr.swap(dr0); advect(0, dr, dr0, u, v, dt);
        dg.swap(dg0); diffuse(0, dg, dg0, diff, dt); dg.swap(dg0); advect(0, dg, dg0, u, v, dt);
        db.swap(db0); diffuse(0, db, db0, diff, dt); db.swap(db0); advect(0, db, db0, u, v, dt);
    }
};

struct GeyserState {
    float next = 0.0f;
    bool par1 = true, par2 = true;
    int dist = 1;
    float col[3] = {1.0f, 1.0f, 1.0f};
};

void resetGeyser(GeyserState& g) { 
    g.par1 = rand() % 2; // random numbers not guaranteed to be uniform here
    g.par2 = rand() % 2;
    g.dist = rand() % N + 1; 
    int i = rand() % 6;
    float f = (rand() / (float) RAND_MAX);
    float s = 0.85f;
    float pal[3] = {1 - s, 1 - s * f, 1 - s * (1 - f)};
    g.col[0] = pal[i % 3];
    g.col[1] = pal[i > 2 ? (i + 1) % 3 : (i + 2) % 3];
    g.col[2] = pal[i > 2 ? (i + 2) % 3 : (i + 1) % 3];
} // resetGeyser

void updateGeyser(FluidSim& sim, GeyserState& g, float now) {
    if (now < g.next) return;
    if (now >= g.next + 3000) {
	g.next = now + 500 + (rand() / (float) RAND_MAX) * 1500;
	resetGeyser(g);
	return;
    }

    float pow = 1.0f, vol = 155.0f;
    int nf = g.par2 ? 2 : N - 1;
    int sgn = g.par2 ? 1 : -1;
    float jet = pow * sgn;
    int jetSize = 3;

    for (int k = -jetSize; k <= jetSize; k++) {
	int dd = clamp(g.dist + k, 1, N);
	int id = g.par1 ? IX(dd, nf) : IX(nf, dd);

	sim.dr[id] += vol * g.col[0];
	sim.dg[id] += vol * g.col[1];
	sim.db[id] += vol * g.col[2];

	if (g.par1) sim.v[id] += jet;
	else	    sim.u[id] += jet;
    }
} // updateGeyser

struct MouseState {
    float pmx = -1, pmy = -1;
    float ppmx = -1, ppmy = -1;
};

void pointerMotion(SDL_Event& e, SDL_Window* win, FluidSim& sim, MouseState& m) {
    int w, h;
    SDL_GetWindowSize(win, &w, &h);

    float gx = 1 + (e.motion.x / (float)w) * (N - 1);
    float gy = 1 + (e.motion.y / (float)h) * (N - 1);

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
	    int id = IX(i, j);
	    float dist = hypot(ii - i, jj - j);
	    float wgt = exp(-2 * dist / r);
	    sim.u[id] += dx * wgt;
	    sim.v[id] += dy * wgt;
	}
    }
}

void initialise(SDL_Window*& win, SDL_Renderer*& ren, SDL_Texture*& tex) {
    SDL_Init(SDL_INIT_VIDEO);

    win = SDL_CreateWindow("Smoke", 1280, 800, SDL_WINDOW_RESIZABLE);
//    win = SDL_CreateWindow("Smoke", 1280, 800, SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIGH_PIXEL_DENSITY);

    ren = SDL_CreateRenderer(win, NULL);
//    SDL_SetRenderVSync(ren, 1);
//    SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);

    tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, N, N);

//    SDL_SetTextureScaleMode(tex, SDL_SCALEMODE_PIXELART); 
    SDL_SetTextureScaleMode(tex, SDL_SCALEMODE_LINEAR); 
}

void render(SDL_Renderer* ren, SDL_Texture* tex, SDL_Window* win, FluidSim& sim, vector<unsigned char>& pixels) {
    for (int j = 1; j <= N; j++) {
	for (int i = 1; i <= N; i++) {
	    int id = IX(i, j);
	    int k = ((j - 1) * N + (i - 1)) * 4;
	    pixels[k]     = clamp((int) sim.dr[id], 0, 255);
	    pixels[k + 1] = clamp((int) sim.dg[id], 0, 255);
	    pixels[k + 2] = clamp((int) sim.db[id], 0, 255);
	    pixels[k + 3] = 255;
	}
    }

    SDL_UpdateTexture(tex, NULL, pixels.data(), N * 4);

    int w, h;
    SDL_GetWindowSize(win, &w, &h);
    SDL_FRect dst = {0, 0, (float) w, (float) h};

    SDL_RenderClear(ren);
    SDL_RenderTexture(ren, tex, NULL, &dst);
    SDL_RenderPresent(ren);
}

int main() {
    SDL_Window* win;
    SDL_Renderer* ren;
    SDL_Texture* tex;
    initialise(win, ren, tex);

    FluidSim sim;
    vector<unsigned char> pixels(N * N * 4);

    MouseState mouse;
    GeyserState geyser;
    resetGeyser(geyser);
    geyser.next = SDL_GetTicks() + 500 + (rand() / (float) RAND_MAX) * 1500;

    float last = SDL_GetTicks() / 1000.0;
    bool running = true;

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_MOUSE_MOTION)
                pointerMotion(e, win, sim, mouse);
        }

        float now = SDL_GetTicks() / 1000.0;
        float dt = min(now - last, 1.0f / 30);
        last = now;

        updateGeyser(sim, geyser, now * 1000);
        sim.velocityStep(dt);
        sim.densityStep(dt);

        render(ren, tex, win, sim, pixels);
    }

    SDL_Quit();
}
