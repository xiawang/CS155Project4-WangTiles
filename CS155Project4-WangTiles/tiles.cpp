#include "tiles.h"


Tile::Tile()
{
}

void Tile::genTexture(){
    
}

Tile::~Tile()
{
}

Tiles::Tiles(Image* scr, int hc, int vc)
{
    
    // sample hc + vc squares of 50/sqrt(2) * 50/sqrt(2) pixel
    // put into himage_, vimage_
    
    // generate model, put into tiles_;
    genTiles();
    
    // generate texture for each tile model
    for (Tile t : tiles_){
        t.genTexture();
    }
        
    
}

void Tiles::genTiles(){
    
}






Image* Tiles::tilePlain()
{
    Image* dest = new Image(2,2);
    return dest;
}