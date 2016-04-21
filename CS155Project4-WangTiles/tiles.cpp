#include "tiles.h"
#include <time.h>
#include <iostream>
#include "image.h"


/*
 * constructor */
Tile::Tile()
{
}

Tile::Tile(int up, int left){
    up_ = up;
    left_ = left;
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

// a bunch of sets and gets
void Tile::setUp(int u){
    up_ = u;
}

void Tile::setDown(int d){
    down_ = d;
}

void Tile::setLeft(int l){
    left_ = l;
}

void Tile::setRight(int r){
    right_ = r;
}

void Tile::setSize(int s){
    tilesize_ = s;
}

void Tile::setTexture(Image* src){
    tex_ = src;
}

int Tile::getUp(){
    return up_;
}

int Tile::getDown(){
    return down_;
}

int Tile::getRight(){
    return right_;
}

int Tile::getLeft(){
    return left_;
}

int Tile::getSize(){
    return tilesize_;
}

Image* Tile::getTexture(){
    return tex_;
}

/*
 * given a N color and a W color (that this tile should have), return true if the tile can be placed in the next spot
 * a color code of -1 indicates no restriction
 * 
 * will be called in Tiles::tilePlain()
 */
bool Tile::validNeighbor(int up, int left){
    if (up > 0 && up != up_){
        return false;
    }
    
    if (left > 0 && left != left_){
        return false;
    }
    
    return true;
}


/*
 * destructor
 */
Tile::~Tile()
{
}

void Tiles::genColors(){
    
    vector<Pixel> colors;
    for (int i=0; i<(hc_ + vc_); i++) {
        int digit0 = (i>>0)&1;
        int digit1 = (i>>1)&1;
        int digit2 = (i>>2)&1;
        cout << digit0 << digit1 << digit2 << endl;
        Pixel* tempp = new Pixel(digit0,digit1,digit2);
        colors.push_back(*tempp);
    }
    random_shuffle ( colors.begin(), colors.end() );
    colors_ = colors;
}

Image* Tiles::genDummyTile(int n, int e, int s, int w){
    
    vector<Pixel> vtiles;
    vector<Pixel> htiles;
    
    for (int i = 0; i < hc_ + vc_; i++){
        if (i < hc_){
            htiles.push_back(colors_[i]);
        } else {
            vtiles.push_back(colors_[i]);
        }
    }
    
    Image* dumb = new Image(tw_,th_);
    
    for (int i = 0; i < tw_; i++){
        for (int j = 0; j < th_; j++){
            if (i > j){
                if ((tw_-i) > j){
                    // North
                    dumb->setPixel(i, j, htiles[n-1]);
                } else {
                    // East
                    dumb->setPixel(i, j, vtiles[e-1]);
                }
            } else {
                if ((th_-i) < j){
                    // South
                    dumb->setPixel(i, j, htiles[s-1]);
                } else {
                    // West
                    dumb->setPixel(i, j, vtiles[w-1]);
                }
            }
        }
    }
    
    return dumb;
}


// Dummy constructor for testing using models in tiles folder
Tiles::Tiles(int hc, int vc, int tw, int th){
    
    hc_ = hc;
    vc_ = vc;
    
    tw_ = tw;
    th_ = th;
    
    genColors();
    genTiles();
    
    for (int i = 0; i < tiles_.size(); i++){
        Tile* t = &tiles_[i];
        t->setTexture(genDummyTile(t->getUp(), t->getRight(), t->getDown(), t->getLeft()));
    }
    

}



Tiles::Tiles(Image* scr, int hc, int vc, int tw, int th)
{
    
    hc_ = hc;
    vc_ = vc;
    tw_ = tw;
    th_ = th;
    
    // sample hc + vc squares of 50/sqrt(2) * 50/sqrt(2) pixel
    // put into himage_, vimage_
    
    // generate model, put into tiles_;
    genTiles();
    
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
void Tiles::genTiles(){
    
    for (int i = 1; i <= hc_; i ++){
        for (int j = 1; j <= vc_; j++) {

            int h1 = rand() % hc_;
            int h2 = rand() % hc_;
            int v1 = rand() % vc_;
            int v2 = rand() % vc_;

            while ((h1 == h2) && (v1 == v2)){
                h1 = rand() % hc_;
                h2 = rand() % hc_;
                v1 = rand() % vc_;
                v2 = rand() % vc_;
            }

            Tile* t1 = new Tile(i, j);
            Tile* t2 = new Tile(i, j);

            t1->setRight(v1+1);
            t1->setDown(h1+1);

            t2->setRight(v2+1);
            t2->setDown(h2+1);
            
            tiles_.push_back(*t1);
            tiles_.push_back(*t2);
            
        }
    }
    
}


int Tiles::getRandomTile(int up, int left){
    
    vector<int> valid;
    for (int i = 0; i < tiles_.size(); i++){
        if (tiles_[i].validNeighbor(up, left)){
            valid.push_back(i);
        }
    }
    
    int index = rand() % valid.size();
    
    return valid[index];
    
}


/*
 * given desired width w and desired height h,
 * generate a tiled plain of that size
 */
Image* Tiles::tilePlain(int w, int h)
{
    
    Image* dest = new Image(w, h);
    int nh = w/tw_ + 1;                // number of tiles horizontally
    int nv = h/th_ + 1;               // number of tiles vertically
    int tileModel[nv][nh];                   // a 2D array that specify types of tiles used
    
    tileModel[0][0] = rand() % tiles_.size();
    
    // generate an abstract valid tiling
    for (int i = 0; i < nv; i++){
        for (int j = 0; j < nh; j++){
            
            int up = (i == 0) ? -1 : tiles_[tileModel[i-1][j]].getDown();
            int left = (j == 0) ? -1 : tiles_[tileModel[i][j-1]].getRight();
            
            tileModel[i][j] = getRandomTile(up, left);
            cout << "i: " << i << endl;
            cout << "j: " << j << endl;
            cout << "model: " << tileModel[i][j] << endl;
            
        }
    }
    
    // fill in the corresponding textured tile according to the tile models
    for (int i = 0; i < h; i ++){
        for (int j = 0; j < w; j++){
            
            int tileIndex = tileModel[i / th_][j / tw_];
            Tile t = tiles_[tileIndex];
            Pixel p = t.getTexture()->getPixel(j % tw_,i % th_);
            dest->setPixel(j, i, p);
        }
    }
    
    
    return dest;
}