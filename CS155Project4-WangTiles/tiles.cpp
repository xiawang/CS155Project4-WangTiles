#include "tiles.h"

/*
 * constructor */
Tile::Tile()
{
}

/*
 * gen a tile with four colors specified, 
 * generate the corresponding texture and store in tex_
 *
 * look at himage_ and vimage_, merge corresponding squares into diamond,
 * and crop the middle to get the texture
 */
void Tile::genTexture(){
    
}

// a bunch of sets and gets to be defined


/*
 * given a N color and a W color, return true if the tile can be placed in the next spot
 * a color code of -1 indicates no restriction
 * 
 * will be called in Tiles::tilePlain()
 */
bool validNeighbor(int up, int left);


/*
 * destructor
 */
Tile::~Tile()
{
}





Tiles::Tiles(Image* scr, int hc, int vc)
{
    
    // sample hc + vc squares of 50/sqrt(2) * 50/sqrt(2) pixel
    // put into himage_, vimage_
    
    // generate model, put into tiles_;
    genTiles(hc, vc);
    
    // generate texture for each tile model
    for (Tile t : tiles_){
        t.genTexture();
    }
        
    
}


/*
 * given colors available for horizontal and vertical edges
 * generate at least 2 * hc * vc color coded tile models and put them into tiles_
 * models will contain all possible combinations of NW edges
 * colors for SE edges will be randomly selected
 */
void Tiles::genTiles(int hc, int vc){
    
}





/*
 * given desired width w and desired height h,
 * generate a tiled plain of that size
 */
Image* Tiles::tilePlain(int w, int h)
{
    Image* dest = new Image(2,2);
    return dest;
}