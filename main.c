#include <SDL3/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define kinematic_viscosity 0.01
#define GRID_SIZE 100
#define CELL_SIZE 8
#define WINDOW_SIZE (GRID_SIZE * CELL_SIZE)

typedef struct {
    double **u;  // x-velocity
    double **v;  // y-velocity
    double **p;  // pressure
} Grid;

Grid* create_grid(int rows, int cols) {
    Grid* grid = (Grid*)malloc(sizeof(Grid));
    grid->u = (double**)malloc(rows * sizeof(double*));
    grid->v = (double**)malloc(rows * sizeof(double*));
    grid->p = (double**)malloc(rows * sizeof(double*));

    for (int i = 0; i < rows; i++) {
        grid->u[i] = (double*)malloc(cols * sizeof(double));
        grid->v[i] = (double*)malloc(cols * sizeof(double));
        grid->p[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            grid->u[i][j] = 0.0;
            grid->v[i][j] = 0.0;
            grid->p[i][j] = 0.0;
        }
    }
    return grid;
}

void free_grid(Grid* grid, int rows) {
    for (int i = 0; i < rows; i++) {
        free(grid->u[i]);
        free(grid->v[i]);
        free(grid->p[i]);
    }
    free(grid->u);
    free(grid->v);
    free(grid->p);
    free(grid);
}

void diffuse(Grid* grid, double dt) {
    double a = dt * kinematic_viscosity * GRID_SIZE * GRID_SIZE;
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i < GRID_SIZE-1; i++) {
            for (int j = 1; j < GRID_SIZE-1; j++) {
                grid->u[i][j] = (grid->u[i][j] + a * (grid->u[i-1][j] + grid->u[i+1][j] + 
                    grid->u[i][j-1] + grid->u[i][j+1])) / (1 + 4 * a);
                grid->v[i][j] = (grid->v[i][j] + a * (grid->v[i-1][j] + grid->v[i+1][j] + 
                    grid->v[i][j-1] + grid->v[i][j+1])) / (1 + 4 * a);
            }
        }
    }
}

void project(Grid* grid) {
    double h = 1.0 / GRID_SIZE;
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i < GRID_SIZE-1; i++) {
            for (int j = 1; j < GRID_SIZE-1; j++) {
                grid->p[i][j] = (grid->p[i-1][j] + grid->p[i+1][j] + 
                    grid->p[i][j-1] + grid->p[i][j+1] - 
                    h * (grid->u[i+1][j] - grid->u[i-1][j] + 
                    grid->v[i][j+1] - grid->v[i][j-1])) / 4;
            }
        }
    }
}

void advect(Grid* grid, double dt) {
    double dt0 = dt * GRID_SIZE;
    for (int i = 1; i < GRID_SIZE-1; i++) {
        for (int j = 1; j < GRID_SIZE-1; j++) {
            double x = i - dt0 * grid->u[i][j];
            double y = j - dt0 * grid->v[i][j];
            
            if (x < 0.5) x = 0.5;
            if (x > GRID_SIZE-1.5) x = GRID_SIZE-1.5;
            if (y < 0.5) y = 0.5;
            if (y > GRID_SIZE-1.5) y = GRID_SIZE-1.5;
            
            int i0 = (int)x;
            int j0 = (int)y;
            int i1 = i0 + 1;
            int j1 = j0 + 1;
            
            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;
            
            grid->u[i][j] = s0 * (t0 * grid->u[i0][j0] + t1 * grid->u[i0][j1]) +
                          s1 * (t0 * grid->u[i1][j0] + t1 * grid->u[i1][j1]);
            grid->v[i][j] = s0 * (t0 * grid->v[i0][j0] + t1 * grid->v[i0][j1]) +
                          s1 * (t0 * grid->v[i1][j0] + t1 * grid->v[i1][j1]);
        }
    }
}

void step(Grid* grid, double dt) {
    diffuse(grid, dt);
    project(grid);
    advect(grid, dt);
    project(grid);
}

void add_force(Grid* grid, int x, int y, double fx, double fy) {
    if (x >= 0 && x < GRID_SIZE && y >= 0 && y < GRID_SIZE) {
        grid->u[x][y] += fx;
        grid->v[x][y] += fy;
    }
}

void render_fluid(SDL_Renderer* renderer, Grid* grid) {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            double speed = sqrt(grid->u[i][j] * grid->u[i][j] + grid->v[i][j] * grid->v[i][j]);
            double normalized_speed = fmin(speed * 40, 1.0);
            
            int r, g, b;
            if (normalized_speed < 0.5) {
                r = (int)(255 * normalized_speed * 2);
                g = 0;
                b = 0;
            } else {
                r = 255;
                g = (int)(255 * (normalized_speed - 0.5) * 2);
                b = 0;
            }
            
            SDL_SetRenderDrawColor(renderer, r, g, b, 255);
            SDL_RenderPoint(renderer, i * CELL_SIZE, j * CELL_SIZE);
        }
    }
}

void add_vortex(Grid* grid, int center_x, int center_y, double strength) {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            double dx = i - center_x;
            double dy = j - center_y;
            double dist = sqrt(dx*dx + dy*dy);
            if (dist < 30) {
                double angle = atan2(dy, dx);
                double force = strength * (1 - dist/30);
                add_force(grid, i, j, -sin(angle) * force, cos(angle) * force);
            }
        }
    }
}

void initialize_fluid(Grid* grid) {
    // Add multiple strong vortices
    add_vortex(grid, GRID_SIZE/3, GRID_SIZE/2, 5.0);
    add_vortex(grid, 2*GRID_SIZE/3, GRID_SIZE/2, -5.0);
    add_vortex(grid, GRID_SIZE/2, GRID_SIZE/3, 4.0);
    add_vortex(grid, GRID_SIZE/2, 2*GRID_SIZE/3, -4.0);
}

int main(int argc, char *argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow(
        "SDL Fluid Simulation",
        WINDOW_SIZE,
        WINDOW_SIZE,
        SDL_WINDOW_RESIZABLE
    );    
    if (!window) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    if (!renderer) {
        printf("Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    Grid* grid = create_grid(GRID_SIZE, GRID_SIZE);
    initialize_fluid(grid);
    
    int quit = 0;
    SDL_Event e;
    int mouse_x = 0, mouse_y = 0;
    int last_mouse_x = 0, last_mouse_y = 0;
    int is_dragging = 0;

    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_EVENT_QUIT) {
                quit = 1;
            } else if (e.type == SDL_EVENT_MOUSE_MOTION) {
                mouse_x = e.motion.x / CELL_SIZE;
                mouse_y = e.motion.y / CELL_SIZE;
                
                if (is_dragging) {
                    // Calculate vortex strength based on mouse movement
                    double dx = mouse_x - last_mouse_x;
                    double dy = mouse_y - last_mouse_y;
                    double speed = sqrt(dx*dx + dy*dy);
                    double strength = speed * 5.0; // Scale factor for vortex strength
                    
                    // Create vortex at current mouse position
                    add_vortex(grid, mouse_x, mouse_y, strength);
                }
                
                last_mouse_x = mouse_x;
                last_mouse_y = mouse_y;
            } else if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN) {
                if (e.button.button == SDL_BUTTON_LEFT) {
                    is_dragging = 1;
                }
            } else if (e.type == SDL_EVENT_MOUSE_BUTTON_UP) {
                if (e.button.button == SDL_BUTTON_LEFT) {
                    is_dragging = 0;
                }
            }
        }

        step(grid, 0.2);

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        render_fluid(renderer, grid);
        SDL_RenderPresent(renderer);
        SDL_Delay(16);
    }

    free_grid(grid, GRID_SIZE);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}


