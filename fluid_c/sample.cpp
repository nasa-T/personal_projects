#include <iostream>
#include <SDL2/SDL.h>
#include <cmath>
// #include "fluid.h"

class Test {
    public:
        Test(int height_, int width_): height(height_), width(width_) {
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
            SDL_CreateWindowAndRenderer(width, height, 0, &window, &renderer);
            // SDL_RenderSetScale(renderer, SCALE_W, SCALE_H);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
            SDL_RenderClear(renderer);
            screenSurface = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);
            SDL_Rect rect{100,100,200,200};
            SDL_SetRenderDrawColor(renderer, 155, 0, 155, 255);
            // SDL_FillRect(screenSurface, NULL, SDL_MapRGB(screenSurface->format, 0xFF, 0xFF, 0xFF));                        
            SDL_RenderDrawRect(renderer, &rect);
            SDL_SetRenderDrawColor(renderer, 155, 0, 0, 255);

            SDL_RenderFillRect(renderer, &rect);
            // SDL_UpdateWindowSurface(window);
            SDL_Rect rect2{150,150,100,150};

            SDL_SetRenderDrawColor(renderer, 0, 100, 0, 100);
            SDL_RenderFillRect(renderer, &rect2);

            // Clear the newly created window
            SDL_RenderPresent(renderer);
        }
            // if (event.type == SDL_MOUSEBUTTONDOWN) {
                //     bool mouseDown = true;
                //     int buttonType = event.button.button;
                //     while (mouseDown) {
                //         if (buttonType == SDL_BUTTON_LEFT) std::cout << "ur mom\n";
                //         else if (buttonType == SDL_BUTTON_RIGHT) std::cout << "ur dad\n";
                //         SDL_PollEvent(&event);
                //         if (event.type == SDL_MOUSEBUTTONUP) mouseDown = false;
                //         SDL_Delay(100);
                //         SDL_Delay(1000*consts::dt * consts::SCALE_TIME);
                //         draw_cells();
                //         grid.solveVelocity();
                //         grid.moveMass();
                //         SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
                //         SDL_RenderClear(renderer);
                //     }
                // }
    private:
        int height, width;
        float SCALE_H, SCALE_W;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
};


// int main(int argc, char* argv[]) {
//     Test sim(800, 800);
//     // sim.draw_grid();
//     SDL_Event event;

//     while(!(event.type == SDL_QUIT)){
//         SDL_Delay(10);  // setting some Delay
//         SDL_PollEvent(&event);  // Catching the poll event.
//     }
    
//     return 0;
// }