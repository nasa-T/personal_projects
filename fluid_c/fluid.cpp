#include <cmath>
#include <iostream>
#include "fluid.h"
#include <SDL2/SDL.h>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
GLFWwindow* window;

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
using namespace glm;

void print(char *string) {
    std::cout << string << std::endl;
}

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
        float getMag() {
            return pow(vx*vx+vy*vy,1/2);
        }
        // std::ostream& operator<<(std::ostream &s, VelocityVector &vec) {
        //     return s << "vx: " << vec.getVx() << " vy: " << vec.getVy();
        // }
    private:
        float vx, vy, ax, ay;
};

class VelocityBox {
    public:
        float top, bottom, left, right;
        VelocityBox(float top = 0, float bottom = 0, float left = 0, float right = 0): top(top), bottom(bottom), left(left), right(right) {}
        void setVelocity(VelocityVector v) {
            left = v.getVx();
            right = v.getVx();
            top = v.getVy();
            bottom = v.getVy();
        }
};

class VelocityGrid {
    public:
        VelocityGrid(int rows, int cols): rows(rows), cols(cols) {
            int i, j;
            vyArray = (float **) malloc(sizeof(float *) * (rows + 1));
            vxArray = (float **) malloc(sizeof(float *) * rows);
            for (i = 0; i < rows + 1; i++) {
                vyArray[i] = (float *) malloc(sizeof(float) * cols);
                for (j = 0; j < cols; j++) {
                    vyArray[i][j] = 0;
                }
            }
            for (i = 0; i < rows; i++) {
                vxArray[i] = (float *) malloc(sizeof(float) * (cols + 1));
                for (j = 0; j < cols + 1; j++) {
                    vxArray[i][j] = 0;
                }
            }
        }

        VelocityBox getVelocityBox(int i, int j) {
            return VelocityBox(vyArray[i][j], vyArray[i+1][j], vxArray[i][j], vxArray[i][j+1]);
        }

        VelocityVector getVelocityVector(int i, int j) {
            return VelocityVector(vxArray[i][j], vyArray[i][j]);
        }
    private:
        int rows, cols;
        float **vyArray, **vxArray;
};

class Neighbors {
    public:
        FluidCell *top, *bottom, *left, *right;
        Neighbors() {
            top = NULL;
            bottom = NULL;
            left  = NULL;
            right = NULL;
        }
};


class FluidCell {
    public:
        FluidCell(float mass, float width, float height, float temp): mass(mass), width(width), height(height), temperature(temp) {
            // size is the physical area
            size = width * height;
            density = mass/size;
            pressure = density/mass * consts::kb * temperature;
            velocity = VelocityVector();
            // velocity.setVx(1);
            // velocity.setVx((std::rand() % (int)(2*width)) - width/2);
            // velocity.setVy((std::rand() % (int)(2*height)) - height/2);
            vBounds = VelocityBox(0,0,0,0);
            vBounds.setVelocity(velocity);
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
        void addMass(float m) {
            setMass(mass + m);
        }

        FluidCell *getRight() {
            return neighbors.right;
        }
        FluidCell *getLeft() {
            return neighbors.left;
        }
        FluidCell *getTop() {
            return neighbors.top;
        }
        FluidCell *getBottom() {
            return neighbors.bottom;
        }

        void setRight(FluidCell *neighbor) {
            neighbors.right = neighbor;
        }
        void setLeft(FluidCell *neighbor) {
            neighbors.left = neighbor;
        }
        void setTop(FluidCell *neighbor) {
            neighbors.top = neighbor;
        }
        void setBottom(FluidCell *neighbor) {
            neighbors.bottom = neighbor;
        }

        void transferMass(FluidCell *other, float m) {
            if (other != NULL) { // halt at the boundaries; NOTHING GETS OUT
                if (m > 0) {
                    if (this->getMass() < m) {
                        m = this->getMass();
                    }
                    float otherMass = other->getMass();
                    other->setMass(otherMass + m);
                    this->setMass(mass - m);
                } else if (m < 0) {
                    m = -m;
                    if (other->getMass() < m) {
                        m = other->getMass();
                    }
                    float otherMass = other->getMass();
                    other->setMass(otherMass - m);
                    this->setMass(mass + m);
                }
            }
        }
        // void addVelocity(VelocityVector vel2) {
        //     vel = vel + vel2;
        // }
        void setVelocity(VelocityVector vel2) {
            velocity = vel2;
            vBounds.setVelocity(vel2);
        }
        VelocityVector getVelocity() {
            return velocity;
        }

        // void solveVelocity() {
        //     // pressure difference first
            
        // }

        char isActive() {
            char b, t, l, r = 0;
            if (neighbors.bottom != NULL) b = !!(neighbors.bottom->getMass());
            if (neighbors.top != NULL) t = !!(neighbors.top->getMass());
            if (neighbors.left != NULL) l = !!(neighbors.left->getMass());
            if (neighbors.right != NULL) r = !!(neighbors.right->getMass());
            // printf("%d ", !!mass || b || t || l || r);
            return !!mass || b || t || l || r;
        }
        
        void setLoc(int r, int c) {
            row = r;
            col = c;
        }

        int getRow() {
            return row;
        }

        int getCol() {
            return col;
        }

    private:
        float mass, size, density, temperature, pressure, width, height;
        int row, col;
        VelocityVector velocity;
        VelocityBox vBounds;
        Neighbors neighbors;
};

class FluidGrid {
    public:
        FluidGrid(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {
            int i, j;
            cells = rows * cols;
            grid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            newGrid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            // vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
            // vGrid[0] = VelocityGrid(r, c);
            cellWidth = width / c;
            cellHeight = height / r;
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            maxV = 0;
            for (i = 0; i < rows; i++) {
                grid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
                newGrid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
                for (j = 0; j < cols; j++) {
                    //random mass for now
                    // float mass = std::rand() % 256;
                    // grid[i][j] = FluidCell(mass, cellWidth, cellHeight, 100);
                    if (i == rows/2 && j == cols/2) {
                        grid[i][j] = FluidCell(100, cellWidth, cellHeight, 100);
                        grid[i][j].setVelocity(VelocityVector(20,0));
                    } 
                    else {
                        grid[i][j] = FluidCell(0, cellWidth, cellHeight, 100);
                    }
                    if (i != 0) {
                        grid[i][j].setTop(&grid[i - 1][j]);
                        grid[i - 1][j].setBottom(&grid[i][j]);
                    }
                    if (j != 0) {
                        grid[i][j].setLeft(&grid[i][j - 1]);
                        grid[i][j - 1].setRight(&grid[i][j]);
                    }
                    if (grid[i][j].isActive()) {
                        activeCells.push_back(&grid[i][j]);
                    }
                    if (grid[i][j].getVelocity().getMag() > maxV) {
                        maxV = grid[i][j].getVelocity().getMag();
                    }
                    grid[i][j].setLoc(i, j);
                    newGrid[i][j] = grid[i][j];
                }
            }
            dt = 0.02;
            // 0.7/maxV;
            nActive = activeCells.size();
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

        void advect() {
            int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    FluidCell cell = *getCell(i,j);
                    VelocityVector v = cell.getVelocity();
                    float physX = cellWidth*j+cellWidth/2;
                    float physY = cellHeight*i+cellHeight/2;
                    float prevX = physX - dt*v.getVx();
                    float prevY = physY - dt*v.getVy();
                    // printf("%f ", prevX);
                    // std::cout << physX << " ";
                    if ((prevX < 0) || (prevY < 0) || (prevX > width) || (prevY > height)) {
                        newGrid[i][j].setMass(0);
                    } else {
                        int prev_i = prevY / cellHeight;
                        int prev_j = prevX / cellWidth;
                        // printf("%f; %f, %f: %d, %d\n", cell.getMass(), prevX, prevY, prev_j, prev_i);
                        
                        newGrid[i][j].setMass(grid[prev_i][prev_j].getMass());
                        //  = grid[prev_i][prev_j];
                        // newGrid[i][j].setLoc(i,j);
                    }
                }
            }
        }

        void update(SDL_Event event) {
            int i, j;
            if (event.type == SDL_MOUSEBUTTONDOWN) {
                SDL_MouseButtonEvent buttonEvent = event.button;
                Sint32 x = buttonEvent.x;
                Sint32 y = buttonEvent.y;
                i = y * SCALE_H / cellHeight;
                j = x * SCALE_W / cellWidth;
                FluidCell *clickedCell = getCell(i, j);
                // printf("%f\n", clickedCell->getMass());
                // if (!clickedCell->isActive()) activeCells.push_back(clickedCell);
                clickedCell->addMass(100);
                printf("%d\n", nActive);
            }
            // float totalMass = 0;
            nActive = activeCells.size();
            activeCells.clear();
            // int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    if (getCell(i, j)->isActive()) {
                        activeCells.push_back(getCell(i, j));
                    }
                }
            }
            // for (i = 0; i < nActive; i++) {
            //     FluidCell *cell = activeCells[i];
            //     float initialMass = cell->getMass();
                
            //     cell->transferMass(cell->getBottom(), 5);
            //     totalMass += cell->getMass();
            //     if (!cell->isActive()) activeCells.erase(activeCells.begin() + i);
            //     if (!initialMass && !!cell->getMass()) {
            //         if ((cell->getBottom() != NULL) && (cell->getBottom()->getMass() == 0)) activeCells.push_back(cell->getBottom());
            //         if ((cell->getTop() != NULL) && (cell->getTop()->getMass() == 0)) activeCells.push_back(cell->getTop());
            //         if ((cell->getLeft() != NULL) && (cell->getLeft()->getMass() == 0)) activeCells.push_back(cell->getLeft());
            //         if ((cell->getRight() != NULL) && (cell->getRight()->getMass() == 0)) activeCells.push_back(cell->getRight());
            //         // if ((cell->getBottom() != NULL) && (cell->getBottom()->isActive())) activeCells.push_back(cell->getBottom());
            //         // if ((cell->getTop() != NULL) && (cell->getTop()->isActive())) activeCells.push_back(cell->getTop());
            //         // if ((cell->getLeft() != NULL) && (cell->getLeft()->isActive())) activeCells.push_back(cell->getLeft());
            //         // if ((cell->getRight() != NULL) && (cell->getRight()->isActive())) activeCells.push_back(cell->getRight());
            //     }
                
            // }
            nActive = activeCells.size();
            grid = newGrid;
            
            // advect();
            // print("advected");
            // printf("%f\n", maxV);
            // for (i = 0; i < rows; i++) {
            //     for (j = 0; j < cols; j++) {
            //         //random mass for now
            //         FluidCell *cell = getCell(i, j);
            //         // cell->transferMass(cell->getRight(), 5);
            //         cell->transferMass(cell->getBottom(), 5);
            //         totalMass += cell->getMass();
            //         // cell->setMass(cell->getMass()-5);
            //     }
            // }
        }

        std::vector<FluidCell*> getActive() {
            return activeCells;
        }
        unsigned int nActive;
    private:
        float width, height;
        float cellWidth, cellHeight;
        float dt, maxV;
        int rows, cols, cells;
        FluidCell **grid;
        FluidCell **newGrid;
        std::vector<FluidCell*> activeCells;
        float SCALE_H, SCALE_W;
        // VelocityGrid *vGrid;
        // VelocityGrid *new_vGrid;
};

class Simulator {
    public:
        Simulator(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            cells = rows * cols;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(width, height, rows, cols);
            
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
            // for (i = 0; i < grid->nActive; i++) {
                    float mass = grid->getCell(i,j)->getMass();
                    // float mass = grid->getActive()[i]->getMass();
                    // FluidCell *cell = grid->getActive()[i];
                    
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    // SDL_Rect rect{cell->getCol()*consts::GRID_WIDTH/cols,cell->getRow()*consts::GRID_HEIGHT/rows,(cell->getCol()+1)*consts::GRID_WIDTH/cols,(cell->getRow()+1)*consts::GRID_HEIGHT/rows};
                    if (mass > 255) {
                        mass = 255;
                    }
                    // if (!cell->isActive()) SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
                    SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                    SDL_RenderFillRect(renderer, &rect);
            // }
                }
            }
            SDL_RenderPresent(renderer);
        }
        
        void step(SDL_Event event) {
            grid->update(event);
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


int main(int argv, char **argc) {
    if (argv > 5) std::srand((unsigned) atoi(argc[5]));
    else std::srand((unsigned) std::time(NULL));
    Simulator sim(atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), atoi(argc[4]));
    sim.drawCells();
    SDL_Event event;
    // 
    // sim.step(event);
    // sim.drawGridLines();

    while(!(event.type == SDL_QUIT)){
        SDL_Delay(20);  // setting some Delay
        sim.step(event);
        SDL_PollEvent(&event);  // Catching the poll event.
    }
    sim.freeSim();
    return 1;
}

