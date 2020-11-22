#include <stdint.h>

#define TILE_SIZE 32

unsigned int GetWidth();
unsigned int GetHeight();
void RedrawTile(int tile_x, int tile_y, int ylines=TILE_SIZE);
void SmallPTLoop();
extern int running;
