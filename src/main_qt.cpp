#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPixmap>
#include "smallpt.h"

uint32_t *workBuffer;
QImage *img;
QLabel *l;
QApplication *app;

void RedrawTile(int tile_x, int tile_y, int ylines)
{
    for (int y = 0; y < ylines; y++)
    {
        uint32_t *line = (uint32_t *)(img->scanLine(tile_y * TILE_SIZE + y));
        uint32_t w = GetWidth();
        for (unsigned x=0; x < w; x++)
        {
            line[x] = workBuffer[x + (tile_y*TILE_SIZE+y)*w];
        }
    }

    l->setPixmap(QPixmap::fromImage(*img));
    app->processEvents();
    (void)tile_x;
    (void)tile_y;
    (void)ylines;
}

int main( int argc, char **argv )
{
    QApplication a( argc, argv );
    QLabel label;
    app = &a;
    QImage image(GetWidth(), GetHeight(), QImage::Format_ARGB32);
    l = &label;
    workBuffer = new uint32_t[GetWidth() * GetHeight()];
    img = &image;
    label.setPixmap(QPixmap::fromImage(image));
    label.show();

    SmallPTLoop();

    return a.exec();
}
