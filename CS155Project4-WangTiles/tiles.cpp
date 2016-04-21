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
    for (int i=0; i<8; i++) {
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

Image* Tiles::genDummyTile(int n, int e, int s, int w, int tw, int th){
    
    vector<Pixel> vtiles;
    vector<Pixel> htiles;
    
    htiles.push_back(colors_[0]);
    htiles.push_back(colors_[1]);
    
    vtiles.push_back(colors_[2]);
    vtiles.push_back(colors_[3]);
    
    Image* dumb = new Image(tw,th);
    
//    Pixel* r = new Pixel(1,0,0);
//    Pixel* g = new Pixel(0,1,0);
//    Pixel* y = new Pixel(1,1,0);
//    Pixel* b = new Pixel(0,0,1);
//    
//    htiles.push_back(*r);
//    htiles.push_back(*g);
//    
//    vtiles.push_back(*y);
//    vtiles.push_back(*b);
    
    
    for (int i = 0; i < tw; i++){
        for (int j = 0; j < th; j++){
            if (i > j){
                if ((tw-i) > j){
                    // North
                    dumb->setPixel(i, j, htiles[n-1]);
                } else {
                    // East
                    dumb->setPixel(i, j, vtiles[e-1]);
                }
            } else {
                if ((th-i) < j){
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
Tiles::Tiles(int tw, int th){
    
    genColors();
    
    //1111
    Tile* t1 = new Tile(1,1);
    t1->setRight(1);
    t1->setDown(1);
    t1->setTexture(genDummyTile(1,1,1,1,tw,th));
    tiles_.push_back(*t1);
    
    //1122
    Tile* t2 = new Tile(1,2);
    t2->setRight(1);
    t2->setDown(2);
    t2->setTexture(genDummyTile(1,1,2,2,tw,th));
    tiles_.push_back(*t2);
    
    //1212
    Tile* t3 = new Tile(1,2);
    t3->setRight(2);
    t3->setDown(1);
    t3->setTexture(genDummyTile(1,2,1,2,tw,th));
    tiles_.push_back(*t3);
    
    //1221
    Tile* t4 = new Tile(1,1);
    t4->setRight(2);
    t4->setDown(2);
    t4->setTexture(genDummyTile(1,2,2,1,tw,th));
    tiles_.push_back(*t4);
    
    //2112
    Tile* t5 = new Tile(2,2);
    t5->setRight(1);
    t5->setDown(1);
    t5->setTexture(genDummyTile(2,1,1,2,tw,th));
    tiles_.push_back(*t5);
    
    //2121
    Tile* t6 = new Tile(2,1);
    t6->setRight(1);
    t6->setDown(2);
    t6->setTexture(genDummyTile(2,1,2,1,tw,th));
    tiles_.push_back(*t6);
    
    //2211
    Tile* t7 = new Tile(2,1);
    t7->setRight(2);
    t7->setDown(1);
    t7->setTexture(genDummyTile(2,2,1,1,tw,th));
    tiles_.push_back(*t7);
    
    //2222
    Tile* t8 = new Tile(2,2);
    t8->setRight(2);
    t8->setDown(2);
    t8->setTexture(genDummyTile(2,2,2,2,tw,th));
    tiles_.push_back(*t8);
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
    
    for (int i = 0; i < hc; i ++){
        for (int j = 0; j < vc; j++) {

            int h1 = rand() % hc;
            int h2 = rand() % hc;
            int v1 = rand() % vc;
            int v2 = rand() % vc;

            while ((h1 == h2) && (v1 == v2)){
                h1 = rand() % hc;
                h2 = rand() % hc;
                v1 = rand() % vc;
                v2 = rand() % vc;
            }

            Tile* t1 = new Tile(i, j);
            Tile* t2 = new Tile(i, j);

            t1->setRight(v1);
            t1->setDown(h1);

            t2->setRight(v2);
            t2->setDown(h2);
            
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
Image* Tiles::tilePlain(int w, int h, int tw, int th)
{
    
    Image* dest = new Image(w, h);
    int nh = w/tw + 1;                // number of tiles horizontally
    int nv = h/th + 1;               // number of tiles vertically
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
            
            int tileIndex = tileModel[i / th][j / tw];
            Tile t = tiles_[tileIndex];
            Pixel p = t.getTexture()->getPixel(j % tw,i % th);
            dest->setPixel(j, i, p);
        }
    }
    
    
    return dest;
}