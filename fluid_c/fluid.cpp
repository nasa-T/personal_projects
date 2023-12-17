#include <cmath>
#include <iostream>
#include "fluid.h"
#include <SDL2/SDL.h>
#include "sample.cpp"
#include <cstdlib>

class VelocityVector {
    public:
        VelocityVector(float vx_ = 0, float vy_ = 0, float ax_ = 0, float ay_ = 0): vx(vx_), vy(vy_), ax(ax_), ay(ay_) {}
        VelocityVector operator+(VelocityVector& other) {
            VelocityVector added(this->vx + other.vx, this->vy + other.vy, other.ax, other.ay);
            return added;
        }
        float getVx() {
            return vx;
        }
        float getVy() {
            return vy;
        }
        float getAx() {
            return ax;
        }
        float getAy() {
            return ay;
        }
        void setVx(float v) {
            vx = v;
        }
        void setVy(float v) {
            vy = v;
        }
        void accelerate() {
            vx += ax * consts::dt;
            vy += ay * consts::dt;
        }
        
    private:
        float vx, vy, ax, ay;
};

class Neighbors {
    public:
        Neighbors() {
            top = NULL;
            bottom = NULL;
            left  = NULL;
            right = NULL;
        }
    FluidCell *top, *bottom, *left, *right;
};


class FluidCell {
    public:
        FluidCell(float mass, float width, float height, float temp): mass(mass), width(width), height(height), temperature(temp) {
            // size is the physical area
            size = width * height;
            density = mass/size;
            pressure = density/mass * consts::kb * temperature;
            vel = VelocityVector();
            neighbors = Neighbors();
        }
        float getMass() {
            return mass;
        }
        float getSize() {
            return size;
        }
        float getDensity() {
            return density;
        }
        float getTemp() {
            return temperature;
        }
        float getPressure() {
            return pressure;
        }
        void setMass(float m) {
            if (m >= 0) mass = m;
        }
        void transferMass(FluidCell *other, float m) {
            float otherMass = other->getMass();
            other->setMass(otherMass + m);
        }

    private:
        float mass, size, density, temperature, pressure, width, height;
        VelocityVector vel;
        Neighbors neighbors;
};

class FluidGrid {
    public:
        FluidGrid(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {
            int i, j;
            grid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            cellWidth = width / c;
            cellHeight = height / r;
            // Test sim(height, width);
            cells = rows * cols;
            for (i = 0; i < rows; i++) {
                grid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
                for (j = 0; j < cols; j++) {
                    //random mass for now
                    float mass = std::rand() % 256;
                    grid[i][j] = FluidCell(mass, cellWidth, cellHeight, 100);
                }
            }
        }
        FluidCell *getCell(int i, int j) {
            return &grid[i][j];
        }

        void freeGrid() {
            int i, j;
            for (i = 0; i < rows; i++) {
                free(grid[i]);
            }
            free(grid);
        }

        void update() {
            int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    //random mass for now
                    FluidCell *cell = getCell(i, j);
                    cell->setMass(cell->getMass()-5);
                }
            }
        }

    private:
        float width, height;
        float cellWidth, cellHeight;
        int rows, cols, cells;
        FluidCell **grid;
};

class Simulator {
    public:
        Simulator(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            cells = rows * cols;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(consts::GRID_WIDTH, consts::GRID_HEIGHT, rows, cols);
            
            SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
            SDL_CreateWindowAndRenderer(consts::GRID_WIDTH, consts::GRID_HEIGHT, 0, &window, &renderer);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
            SDL_RenderClear(renderer);
            
        }
        void drawGridLines() {
            int i, j;
            SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);
            for (i = 1; i < rows; i++) {
                SDL_RenderDrawLine(renderer, 0, i*consts::GRID_HEIGHT/rows, consts::GRID_WIDTH, i*consts::GRID_HEIGHT/rows);
            }
            for (j = 1; j < cols; j++) {
                SDL_RenderDrawLine(renderer, j*consts::GRID_WIDTH/cols, 0, j*consts::GRID_WIDTH/cols, consts::GRID_HEIGHT);
            }
            // Clear the newly created window
            SDL_RenderPresent(renderer);
        }
        void drawCells() {
            int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    float mass = grid->getCell(i,j)->getMass();
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                    SDL_RenderFillRect(renderer, &rect);
                }
            }
            SDL_RenderPresent(renderer);
        }
        
        void step() {
            grid->update();
            drawCells();
            // int i, j;
            // for (i = 0; i < rows; i++) {
            //     for (j = 0; j < cols; j++) {
            //         FluidCell *cell = grid->getCell(i, j);
            //         cell->
            //     }
            // }
        }

        void freeSim() {
            grid->freeGrid();
            free(grid);
        }
    private:
        float width, height;
        int rows, cols, cells;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        float SCALE_H, SCALE_W;
        FluidGrid *grid;
};

std::ostream& operator<<(std::ostream &s, VelocityVector &vec) {
    return s << "vx: " << vec.getVx() << " vy: " << vec.getVy();
}

int main(int argv, char **argc) {
    if (argv > 5) std::srand((unsigned) atoi(argc[5]));
    else std::srand((unsigned) std::time(NULL));
    Simulator sim(atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), atoi(argc[4]));
    sim.drawCells();
    // sim.drawGridLines();

    SDL_Event event;

    while(!(event.type == SDL_QUIT)){
        SDL_Delay(10);  // setting some Delay
        sim.step();
        SDL_PollEvent(&event);  // Catching the poll event.
    }
    sim.freeSim();
    return 1;
}
