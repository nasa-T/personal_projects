#include <cmath>
#include <iostream>
#include "fluid.h"

#include <map>
#include <string>


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
        VelocityGrid(float width, float height, int rows, int cols): width(width), height(height), rows(rows), cols(cols) {
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
            cellHeight = height/rows;
            cellWidth = width/cols;
        }

        VelocityBox getVelocityBox(int i, int j) {
            return VelocityBox(vyArray[i][j], vyArray[i+1][j], vxArray[i][j], vxArray[i][j+1]);
        }
        VelocityVector sampleVelocityAtPoint(float x, float y) {
            int vyi = ceil((height - y) / cellHeight);
            int vyj = floor((x-cellWidth/2) / cellWidth);
            // the physical positions relative to nearest box of vys
            float cvyX, cvyY; 
            float cvxX, cvxY; // same for vxs
            cvyX = (fmod(x+cellWidth/2, cellWidth))/cellWidth;
            cvyY = fmod(y, cellHeight)/cellHeight;

            // will throw an error for when x and y are off the grid entirely
            float vy10, vy00, vy01, vy11 = 0;
            if ((vyj > -1)) {
                vy10 = vyArray[vyi][vyj];
            }
            if ((vyj > -1) && (vyi > 0)) {
                vy00 = vyArray[vyi-1][vyj];
            }
            if ((vyj < cols-1) && (vyi > 0)) {
                vy01 = vyArray[vyi-1][vyj+1];
            }
            if (vyj < cols-1) {
                vy11 = vyArray[vyi][vyj+1];
            }

            float vy = (1-cvyX)*(1-cvyY)*vy10 + cvyX*(1-cvyY)*vy11 + (1-cvyX)*cvyY*vy00 + cvyX*cvyY*vy01;

            int vxi = ceil((height - (y+cellHeight/2)) / cellHeight);
            int vxj = floor(x / cellWidth);

            cvxX = fmod(x, cellWidth)/cellWidth;
            cvxY = (fmod(y+cellHeight/2, cellHeight))/cellHeight;

            float vx10, vx00, vx01, vx11 = 0;
            if ((vxi < rows)) {
                vx10 = vxArray[vxi][vxj];
            }
            if ((vxi > 0)) {
                vx00 = vxArray[vxi-1][vxj];
            }
            if ((vxi > 0)) {
                vx01 = vxArray[vxi-1][vxj+1];
            }
            if (vxi < rows) {
                vx11 = vxArray[vxi][vxj+1];
            }
            float vx = (1-cvxX)*(1-cvxY)*vx10 + cvxX*(1-cvxY)*vx11 + (1-cvxX)*cvxY*vx00 + cvxX*cvxY*vx01;
            return VelocityVector(vx, vy);
        }
        VelocityVector getVelocityVector(int i, int j) {
            return VelocityVector(vxArray[i][j], vyArray[i][j]);
        }
        float getVy(int i, int j) {
            return vyArray[i][j];
        }
        float getVx(int i, int j) {
            return vxArray[i][j];
        }
        void setVy(int i, int j, float vy) {
            vyArray[i][j] = vy;
        }
        void setVx(int i, int j, float vx) {
            vxArray[i][j] = vx;
        }
        float getWidth() {
            return width;
        }
        float getHeight() {
            return height;
        }
        float getCellWidth() {
            return cellWidth;
        }
        float getCellHeight() {
            return cellHeight;
        }
    private:
        float width, height;
        float cellWidth, cellHeight;
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

        void setVelocity(float vx, float vy) {
            velocity.setVx(vx);
            velocity.setVy(vy);
        }

        VelocityVector getVelocity() {
            return velocity;
        }

        char isActive() {
            char b, t, l, r = 0;
            if (neighbors.bottom != NULL) b = !!(neighbors.bottom->getMass());
            if (neighbors.top != NULL) t = !!(neighbors.top->getMass());
            if (neighbors.left != NULL) l = !!(neighbors.left->getMass());
            if (neighbors.right != NULL) r = !!(neighbors.right->getMass());
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
        // VelocityBox vBounds;
        Neighbors neighbors;
};

class Source {
    public:
        Source(float _x, float _y, FluidGrid *g, VelocityGrid *vgrid, uint sourceType):  grid(g), vGrid(vgrid), type(sourceType) {
            setX(_x);
            setY(_y);
            switch (sourceType) {
                case SMOKEGUN:
                    i = floor((vgrid->getHeight() - ys) / vgrid->getCellHeight());
                    j = floor(xs / vgrid->getCellWidth());
                    grid->getCell(i,j)->addMass(100);
                    break;
                case FAN:
                    break;
                case POINTSOURCE:
                    i = floor((vgrid->getHeight() - ys) / vgrid->getCellHeight());
                    j = floor(xs / vgrid->getCellWidth());
                    grid->getCell(i,j)->addMass(100);
                    setVx(vgrid->getWidth()/4);
                    setVy(vgrid->getHeight()/4);
                    break;
            }
        }
        void setX(float x) {
            xs = x;
        }
        void setY(float y) {
            ys = y;
        }
        float getX() {
            return xs;
        }
        float getY() {
            return ys;
        }
        float getVx() {
            return vx;
        }
        float getVy() {
            return vy;
        }
        void setVx(float v) {
            int i = floor((vGrid->getHeight() - ys) / vGrid->getCellHeight());
            int j = floor(xs / vGrid->getCellWidth());
            vx = v;
            if (type == SMOKEGUN) {
                vGrid->setVx(i,j,vx);
            } else if (type == POINTSOURCE) {
                vGrid->setVx(i,j,-vx);
            }
            vGrid->setVx(i,j+1,vx);
        }
        void setVy(float v) {
            int i = floor((vGrid->getHeight() - ys) / vGrid->getCellHeight());
            int j = floor(xs / vGrid->getCellWidth());
            vy = v;
            vGrid->setVy(i,j,vy);
            if (type == SMOKEGUN) {
                vGrid->setVy(i+1,j,vy);
            } else if (type == POINTSOURCE) {
                vGrid->setVy(i+1,j,-vy);
            }
        }
        void addMass() {
            int i = floor((vGrid->getHeight() - ys) / vGrid->getCellHeight());
            int j = floor(xs / vGrid->getCellWidth());
            if (getType() == POINTSOURCE) {
                grid->getCell(i,j)->addMass(5);
            } else if (getType() == SMOKEGUN) {
                grid->getCell(i,j)->addMass(200);
            }
        }
        uint getType() {
            return type;
        }
    private:
        uint type;
        FluidGrid *grid;
        VelocityGrid *vGrid;
        float vx, vy = 0;
        float xs, ys;
        int i, j;
};

        FluidGrid::FluidGrid(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {
            int i, j;
            cells = rows * cols;
            grid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            newGrid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
            vGrid[0] = VelocityGrid(width, height, r, c);
            new_vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
            new_vGrid[0] = VelocityGrid(width, height, r, c);
            sourceArray = (Source *) malloc(MAXSOURCES*sizeof(Source));
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
                    grid[i][j] = FluidCell(0, cellWidth, cellHeight, 100);

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
                    new_vGrid->setVx(i,j,vGrid->getVx(i,j));
                    new_vGrid->setVy(i,j,vGrid->getVy(i,j));
                }
            }
            dt = 0.05;
            nActive = activeCells.size();
        }

        FluidCell *FluidGrid::getCell(int i, int j) {
            return &grid[i][j];
        }

        void FluidGrid::freeGrid() {
            int i, j;
            for (i = 0; i < rows; i++) {
                free(grid[i]);
            }
            free(grid);
        }

        std::map<uint, float> FluidGrid::sampleCellAtPoint(float x, float y) {
            int i = ceil((height - y - cellHeight/2) / cellHeight);
            int j = floor((x - cellWidth/2) / cellWidth);

            // the physical positions relative to nearest box of cells
            float X, Y; 
            X = fmod(x+cellWidth/2, cellWidth)/cellWidth;
            Y = fmod(y+cellHeight/2, cellHeight)/cellHeight;

            // will still throw an error for when x and y are off the grid entirely
            std::map<uint, float> p10, p00, p01, p11;
            if ((j > -1) && (i < rows)) {
                FluidCell *cell = getCell(i,j);
                p10[MASS] = cell->getMass();
                p10[TEMPERATURE] = cell->getTemp();
                p10[PRESSURE] = cell->getPressure();
            } else {
                p10[MASS] = 0;
                p10[TEMPERATURE] = 0;
                p10[PRESSURE] = 0;
            }
            if ((j > -1) && (i > 0)) {
                FluidCell *cell = getCell(i-1,j);
                p00[MASS] = cell->getMass();
                p00[TEMPERATURE] = cell->getTemp();
                p00[PRESSURE] = cell->getPressure();
            } else {
                p00[MASS] = 0;
                p00[TEMPERATURE] = 0;
                p00[PRESSURE] = 0;
            }
            if ((j < cols-1) && (i > 0)) {
                FluidCell *cell = getCell(i-1,j+1);
                p01[MASS] = cell->getMass();
                p01[TEMPERATURE] = cell->getTemp();
                p01[PRESSURE] = cell->getPressure();
            } else {
                p01[MASS] = 0;
                p01[TEMPERATURE] = 0;
                p01[PRESSURE] = 0;
            }
            if ((j < cols-1) && (i < rows)) {
                FluidCell *cell = getCell(i,j+1);
                p11[MASS] = cell->getMass();
                p11[TEMPERATURE] = cell->getTemp();
                p11[PRESSURE] = cell->getPressure();
            } else {
                p11[MASS] = 0;
                p11[TEMPERATURE] = 0;
                p11[PRESSURE] = 0;
            }
            std::map<uint, float> props;
            props[MASS] = (1-X)*(1-Y)*p10[MASS] + X*(1-Y)*p11[MASS] + (1-X)*Y*p00[MASS] + X*Y*p01[MASS];
            props[TEMPERATURE] = (1-X)*(1-Y)*p10[TEMPERATURE] + X*(1-Y)*p11[TEMPERATURE] + (1-X)*Y*p00[TEMPERATURE] + X*Y*p01[TEMPERATURE];
            props[PRESSURE] = (1-X)*(1-Y)*p10[PRESSURE] + X*(1-Y)*p11[PRESSURE] + (1-X)*Y*p00[PRESSURE] + X*Y*p01[PRESSURE];

            return props;
        }

        void FluidGrid::projection(int iters) {
            int i, j;
            for (int n = 0; n < iters; n++) {
                for (i = 0; i < rows; i++) {
                    for (j = 0; j < cols; j++) {
                        int sLeft = !!j;
                        int sRight = !!(cols - 1 - j);
                        int sTop = !!i;
                        int sBottom = !!(rows - 1 - i);                        
                        int s = sLeft + sRight + sTop + sBottom;

                        float left = vGrid->getVx(i, j);
                        float right = vGrid->getVx(i, j+1);
                        float top = vGrid->getVy(i, j);
                        float bottom = vGrid->getVy(i+1, j);
                        float d = -left + right + top - bottom;

                        vGrid->setVx(i, j, left+d*sLeft/s);
                        vGrid->setVx(i, j+1, right-d*sRight/s);
                        vGrid->setVy(i, j, top-d*sTop/s);
                        vGrid->setVy(i+1, j, bottom+d*sBottom/s);
                    }
                }
            }
        }

        void FluidGrid::advect() {
            int i, j;
            for (i = 0; i < rows+1; i++) {
                for (j = 0; j < cols+1; j++) {
                    if ((i < rows) && (j < cols)) {
                        FluidCell cell = *getCell(i,j);
                        float physX = cellWidth*j+cellWidth/2;
                        float physY = height - (cellHeight*i+cellHeight/2); // I want y pointed up
                        VelocityVector vCell = vGrid->sampleVelocityAtPoint(physX, physY);
                        float prevX = physX - dt*vCell.getVx();
                        float prevY = physY - dt*vCell.getVy();
                        std::map<uint, float> prevProps;
                        if (round(vCell.getVx()+vCell.getVy()) == 0) {
                            prevProps[MASS] = cell.getMass();
                            prevProps[TEMPERATURE] = cell.getTemp();
                            prevProps[PRESSURE] = cell.getPressure();
                        } else {
                            prevProps = sampleCellAtPoint(prevX, prevY);
                        }

                        if ((prevX < 0) || (prevY < 0) || (prevX > width) || (prevY > height)) {
                            newGrid[i][j].setMass(0);

                        } else {
                            newGrid[i][j].setMass(prevProps[MASS]);
                        }
                    }
                    if (i < rows) {
                        // vxs of vgrid
                        float vxPhysX = cellWidth*j;
                        float vxPhysY = height - cellHeight*i - cellHeight/2;
                        VelocityVector vVx = vGrid->sampleVelocityAtPoint(vxPhysX, vxPhysY);
                        float prevVxX = vxPhysX - dt*vVx.getVx();
                        float prevVxY = vxPhysY - dt*vVx.getVy();
                        new_vGrid->setVx(i,j,vGrid->sampleVelocityAtPoint(prevVxX, prevVxY).getVx());
                    }
                    if (j < cols) {
                        // vys of vgrid
                        float vyPhysX = cellWidth*j + cellWidth/2;
                        float vyPhysY = height - cellHeight*i;
                        VelocityVector vVy = vGrid->sampleVelocityAtPoint(vyPhysX, vyPhysY);
                        float prevVyX = vyPhysX - dt*vVy.getVx();
                        float prevVyY = vyPhysY - dt*vVy.getVy();
                        new_vGrid->setVy(i,j,vGrid->sampleVelocityAtPoint(prevVyX, prevVyY).getVy());
                        
                    }
                }
            }
            grid = newGrid;
            vGrid = new_vGrid;
        }

        void FluidGrid::update(SDL_Event event) {
            int i, j;
            
            if ((event.type == SDL_MOUSEBUTTONDOWN) || (event.type != SDL_MOUSEBUTTONUP)) {
                
                SDL_MouseButtonEvent buttonEvent = event.button;
                if (buttonEvent.button == SDL_BUTTON_RIGHT) {
                    mouseVelFlag = !mouseVelFlag;
                } else if (buttonEvent.button == SDL_BUTTON_LEFT) {
                    buttonHeld += 1;
                    Sint32 x = buttonEvent.x;
                    Sint32 y = buttonEvent.y;
                    if ((prevMouseX == 0) || (prevMouseY == 0)) {
                        prevMouseX = x;
                        prevMouseY = y;
                    }
                    i = y * SCALE_H / cellHeight;
                    j = x * SCALE_W / cellWidth;
                    if (mouseVelFlag == 1) {
                        float vxMouse = (x - prevMouseX)*SCALE_W;
                        float vyMouse = -(y - prevMouseY)*SCALE_H;
                        vGrid->setVx(i, j+1, vxMouse);
                        vGrid->setVx(i, j, vxMouse);
                        vGrid->setVy(i, j, vyMouse);
                        vGrid->setVy(i+1, j, vyMouse);
                    } else if (mouseVelFlag == 0) {
                        FluidCell *clickedCell = getCell(i, j);
                        clickedCell->addMass(100);
                    } else if ((mouseVelFlag == 2) && (buttonHeld < 2) && (sourceList.size() < MAXSOURCES)) {
                        sourceArray[sourceList.size()] = Source(x*SCALE_W, height-y*SCALE_H, this, vGrid, SMOKEGUN);
                        sourceList.push_back(&sourceArray[sourceList.size()]);
                    } else if ((mouseVelFlag == 3) && (buttonHeld < 2) && (sourceList.size() < MAXSOURCES)) {

                        sourceArray[sourceList.size()] = Source(x*SCALE_W, height-y*SCALE_H, this, vGrid, POINTSOURCE);
                        sourceList.push_back(&sourceArray[sourceList.size()]);
                    }
                } else if (event.key.type == SDL_KEYDOWN) {
                    switch (event.key.keysym.sym) {
                        case SDLK_1:
                            mouseVelFlag = 2;
                            break;
                        case SDLK_2:
                            mouseVelFlag = 3;
                            break;
                        case SDLK_3:
                            mouseVelFlag = 4;
                    }
                }

            } else if (event.type == SDL_MOUSEBUTTONUP) {
                buttonHeld = 0;
                if (mouseVelFlag == 2) {
                    SDL_MouseButtonEvent buttonEvent = event.button;
                    Sint32 x = buttonEvent.x;
                    Sint32 y = buttonEvent.y;
                    float vxMouse = (x - prevMouseX)*SCALE_W;
                    float vyMouse = -(y - prevMouseY)*SCALE_H;
                    Source *s = sourceList.back();
                    s->setVx(vxMouse);
                    s->setVy(vyMouse);
                }
                prevMouseX = 0;
                prevMouseY = 0;
            }
            for (int i = 0; i < sourceList.size(); i++) {
                Source *s = sourceList[i];
                if (s->getVx()*s->getVy() != 0) {
                    s->setVx(s->getVx());
                    s->setVy(s->getVy());
                    if ((s->getType() == SMOKEGUN) || (s->getType() == POINTSOURCE)) {
                        s->addMass();
                    }
                }
            }
            nActive = activeCells.size();
            
            projection(4);
            advect();

        }

        std::vector<FluidCell*> FluidGrid::getActive() {
            return activeCells;
        }

class Simulator {
    public:
        Simulator(float width, float height, int c, int r): width(width), height(height), rows(r), cols(c) {
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
                    float mass = grid->getCell(i,j)->getMass();
                    
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    if (mass > 255) {
                        mass = 255;
                    }
                    SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                    SDL_RenderFillRect(renderer, &rect);

                }
            }
            SDL_RenderPresent(renderer);
        }
        
        void step(SDL_Event event) {
            grid->update(event);
            drawCells();
        }

        void freeSim() {
            grid->freeGrid();
            free(grid);
        }
        float SCALE_H, SCALE_W;
        float width, height;
        FluidGrid *grid;
    private:
        
        int rows, cols, cells;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        
        
};


int main(int argv, char **argc) {
    if (argv > 5) std::srand((unsigned) atoi(argc[5]));
    else std::srand((unsigned) std::time(NULL));
    Simulator sim(atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), atoi(argc[4]));
    sim.drawCells();
    SDL_Event event;

    while(!(event.type == SDL_QUIT)){
        SDL_Delay(20);  // setting some Delay
        sim.step(event);
        SDL_PollEvent(&event);  // Catching the poll event.
    }
    sim.freeSim();
    return 1;
}

