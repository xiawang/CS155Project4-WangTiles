#include "tiles.h"
#include <time.h>
#include <iostream>
#include "image.h"
#include "ip.h"
#include <cmath>


/*
 * constructor */
Tile::Tile()
{
}

Tile::Tile(int up, int left){
    up_ = up;
    left_ = left;
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

void Tiles::initColors(){
    
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


void Tiles::initTextures(Image* src){
    
    int patch_size = tw_/sqrt(2)*(1+OVERLAP);
    int width = src->getWidth();
    int height = src->getHeight();
    
    for (int i = 0; i < hc_; i++){

        Image* patch = new Image(patch_size,patch_size);
        int off_x = rand() % (width - patch_size);
        int off_y = rand() % (height - patch_size);
        
        for (int m=0; m<patch_size; m++) {
            for (int n=0; n<patch_size; n++) {
                
                // Wrap around edges
                double x = m + off_x;
                double y = n + off_y;
                
                double r = src->getPixel(x, y, 0);
                double g = src->getPixel(x, y, 1);
                double b = src->getPixel(x, y, 2);
                patch->setPixel(m, n, 0, r);
                patch->setPixel(m, n, 1, g);
                patch->setPixel(m, n, 2, b);
            }
        }
        himage_.push_back(patch);
    }
    
    for (int j = 0; j < vc_; j++){
        
        Image* patch = new Image(patch_size,patch_size);
        int off_x = rand() % (width - patch_size);
        int off_y = rand() % (height - patch_size);
        
        for (int m=0; m<patch_size; m++) {
            for (int n=0; n<patch_size; n++) {
                
                // Wrap around edges
                double x = m + off_x;
                double y = n + off_y;
                
                double r = src->getPixel(x, y, 0);
                double g = src->getPixel(x, y, 1);
                double b = src->getPixel(x, y, 2);
                patch->setPixel(m, n, 0, r);
                patch->setPixel(m, n, 1, g);
                patch->setPixel(m, n, 2, b);
            }
        }
        vimage_.push_back(patch);
        
    }
    
    
}

/*
 * Return an array of indices indicating the lowest energy seam
 * between the two overlapped subimage (NW vs. SE)
 */
vector<int> Tiles::getVerticalSeam(int n, int e, int s, int w){
    
    Image* north = himage_[n];
    Image* east = vimage_[e];
    Image* south = himage_[s];
    Image* west = vimage_[w];
    
    int wt = tw_/sqrt(2)* OVERLAP * 2;
    int ht = tw_/sqrt(2)*2;
    
    Image* nw_se_1 = new Image(wt, ht);
    Image* nw_se_2 = new Image(wt, ht);
    
    // generate vertical overlap
    for (int i = 0; i < wt; i++){
        for (int j = 0; j < ht; j++){
            
            Pixel p_nw;
            Pixel p_se;
            if (j < tw_/sqrt(2)){
                p_nw = north->getPixel(north->getWidth() - 1 - wt + i, j);
                p_se = east->getPixel(i, j);
            } else {
                p_nw = west->getPixel(west->getWidth() - 1 - wt + i, j - (int)(tw_/sqrt(2)));
                p_se = south->getPixel(i, j - (int)(tw_/sqrt(2)));
            }
            
            // nw_se_1 is the overlap contributed by NW
            nw_se_1->setPixel(i, j, 0, p_nw.getColor(0));
            nw_se_1->setPixel(i, j, 1, p_nw.getColor(1));
            nw_se_1->setPixel(i, j, 2, p_nw.getColor(2));
            
            // nw_se_2 is the overlap contributed by SE
            nw_se_2->setPixel(i, j, 0, p_se.getColor(0));
            nw_se_2->setPixel(i, j, 1, p_se.getColor(1));
            nw_se_2->setPixel(i, j, 2, p_se.getColor(2));
            
        }
    }
    
    // get the dividing seam with least difference between nw_se_1 and nw_se_2
    vector<int> v_carve = seamCarving(nw_se_1, nw_se_2);
    return v_carve;
    
}

/*
 * Return an array of indices indicating the lowest energy seam
 * between the two overlapped subimage (NE vs. SW)
 */
vector<int> Tiles::getHorizontalSeam(int n, int e, int s, int w){
    
    Image* north = himage_[n];
    Image* east = vimage_[e];
    Image* south = himage_[s];
    Image* west = vimage_[w];
    
    int wt = tw_/sqrt(2)*OVERLAP*2;
    int ht = tw_/sqrt(2)*2;
    
    Image* ne_sw_1 = new Image(wt, ht);
    Image* ne_sw_2 = new Image(wt, ht);
    
    // generate horizonal overlap (represented as vertical images)
    for (int i = 0; i < wt; i++){
        for (int j = 0; j < ht; j++){
            
            Pixel p_ne;
            Pixel p_sw;
            if (j < tw_/sqrt(2)){
                p_ne = north->getPixel(j, north->getHeight() - wt + i);
                p_sw = west->getPixel(j, i);
            } else {
                p_ne = east->getPixel(j - (int)(tw_/sqrt(2)), east->getHeight() - wt + i);
                p_sw = south->getPixel(j - (int)(tw_/sqrt(2)), i);
            }
            
            // ne_sw_1 is the overlap contributed by NE
            ne_sw_1->setPixel(i, j, 0, p_ne.getColor(0));
            ne_sw_1->setPixel(i, j, 1, p_ne.getColor(1));
            ne_sw_1->setPixel(i, j, 2, p_ne.getColor(2));
            
            // ne_sw_2 is the overlap contributed by SW
            ne_sw_1->setPixel(i, j, 0, p_sw.getColor(0));
            ne_sw_1->setPixel(i, j, 1, p_sw.getColor(1));
            ne_sw_1->setPixel(i, j, 2, p_sw.getColor(2));
            
        }
    }
    
    // get the dividing seam with least difference between ne_sw_1 and ne_sw_2
    vector<int> h_carve = seamCarving(ne_sw_1, ne_sw_2);
    
    return h_carve;
    
}

/*
 * Given indices of samples, generate the 
 * corresponding collage of the four samples
 *
 *   N | E
 *   W | S
 *
 * Keeps the orientation of the samples (i.e. (0,0) at upper left)
 */
Image* Tiles::genCollage(int n, int e, int s, int w){
    
    Image* north = himage_[n];
    Image* east = vimage_[e];
    Image* south = himage_[s];
    Image* west = vimage_[w];
    
    double perc = 1 - OVERLAP;
    
    // get the dividing seams
    vector<int> v_carve = getVerticalSeam(n, e, s, w);
    vector<int> h_carve = getHorizontalSeam(n, e, s, w);
    
    
    // generate collated image
    int col_size = tw_*sqrt(2);
    Image* collage = new Image(col_size, col_size);
    
    for (int i = 0; i < col_size; i++){ // x
        for (int j = 0; j < col_size; j++){  // y
            
            bool left = i < v_carve[j] + tw_/sqrt(2) * perc;
            bool up = j < h_carve[i] + tw_/sqrt(2) * perc;
            
            Pixel cp;
            if (left){
                if (up){
                    // use north
                    cp = north->getPixel(i, j);
                } else {
                    // use west
                    cp = west->getPixel(i, j - tw_/sqrt(2) * perc);
                }
            } else {
                if (up){
                    // use east
                    cp = east->getPixel(i - tw_/sqrt(2) * perc, j);
                } else {
                    // use south
                    cp = south->getPixel(i - tw_/sqrt(2) * perc, j - tw_/sqrt(2) * perc);
                }
            }
            
            collage->setPixel(i, j, cp);
            
        }
    }
    
    return collage;
    
    
}


/*
 * for each tile with four colors specified,
 * generate the corresponding textures and store in tex_
 *
 * look at himage_ and vimage_, merge corresponding squares into diamond,
 * and crop the middle to get the texture
 */
Image* Tiles::genTextures(int n, int e, int s, int w){
    

    Image* collage = genCollage(n-1, e-1, s-1, w-1);
    //Image* tex = new Image(tw_, th_);
    Image* rot = getMiddleSquare(collage, -45, collage->getWidth()/2, collage->getHeight()/2);
    
    
    return rot;
    
}

/*
 * Given a collage of four samples, rotate 45 degrees clockwise, 
 * and get the square in the middle as the corresponding texture
 */
Image* Tiles::getMiddleSquare (Image* src, double theta, int x, int y)
{
    float sin_val = sin(theta * M_PI / 180);
    float cos_val = cos(theta * M_PI / 180);
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width, height);
    // initialize as all black
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            int new_x = w - x;
            int new_y = h - y;
            // rotate point
            float rot_x = new_x * cos_val - new_y * sin_val;
            float rot_y = new_x * sin_val + new_y * cos_val;
            rot_x += x;
            rot_y += y;
            // resemble
            Pixel pix;
            if (rot_x < 0 || rot_x >= width || rot_y < 0 || rot_y >= height) {
                pix = Pixel(0, 0, 0);
            } else {
                pix = resamplePixel(src, rot_x, rot_y);
            }
            
            dest->setPixel_(w, h, pix);
        }
    }
    
    int start = (sqrt(2)-1)*tw_/2;
    
    int off = 16;
    
    Image* res = new Image(tw_ - off, tw_ - off);
    for (int i = 0; i < tw_ - off; i++){
        for (int j = 0; j < tw_ - off; j++){
            Pixel p = dest->getPixel(start + i + off, start + j + off);
            res->setPixel(i, j, p);
        }
    }
    return res;
}


/*
 * Given the mapped coordinates before rotation, 
 * sample a pixel from the orignal image,
 * using bilinear sampling technique
 */
Pixel Tiles::resamplePixel(Image* src, double x, double y){
    int width = src->getWidth();
    int height = src->getHeight();
    int x1 = floor(x);
    if(x1<0){
        x1 = 0;
    }else if(x1>width-2){
        x1 = width-2;
    }
    int x2 = x1+1;
    int y1 = floor(y);
    if(y1<0){
        y1 = 0;
    }else if(y1>height-2){
        y1 = height-2;
    }
    int y2 = y1+1;
    Pixel upleft = src->getPixel(x1, y1);
    Pixel downleft = src->getPixel(x1, y2);
    Pixel upright = src->getPixel(x2, y1);
    Pixel downright = src->getPixel(x2, y2);
    Pixel average1 = takeWeightedAverage(upleft, upright, x2-x, x-x1);
    Pixel average2 = takeWeightedAverage(downleft, downright, x2-x, x-x1);
    Pixel average = takeWeightedAverage(average1, average2, y2-y, y-y1);
    return average;
}

/*
 * Helper function for resamplePixel
 */
Pixel Tiles::takeWeightedAverage(Pixel pixel1, Pixel pixel2, double weight1, double weight2){
    double r1 = pixel1.getColor(0);
    double g1 = pixel1.getColor(1);
    double b1 = pixel1.getColor(2);
    double r2 = pixel2.getColor(0);
    double g2 = pixel2.getColor(1);
    double b2 = pixel2.getColor(2);
    double r = r1*weight1 + r2*weight2;
    double g = g1*weight1 + g2*weight2;
    double b = b1*weight1 + b2*weight2;
    return Pixel(r,g,b);
}



Image* Tiles::genDummyTexture(int n, int e, int s, int w){
    
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
    
    initColors();
    initTiles();
    
    for (int i = 0; i < tiles_.size(); i++){
        Tile* t = &tiles_[i];
        t->setTexture(genDummyTexture(t->getUp(), t->getRight(), t->getDown(), t->getLeft()));
    }
    

}



Tiles::Tiles(Image* src, int hc, int vc, int tw, int th)
{
    
    hc_ = hc;
    vc_ = vc;
    tw_ = tw;
    th_ = th;
    
    // sample hc + vc squares of 50/sqrt(2) * 50/sqrt(2) pixel
    // put into himage_, vimage_
    
    // generate model, put into tiles_;
    initTiles();
    initTextures(src);
    
    // generate texture for each tile model
    for (int i = 0; i < tiles_.size(); i++){
        Tile* t = &tiles_[i];
        t->setTexture(genTextures(t->getUp(), t->getRight(), t->getDown(), t->getLeft()));
    }
        
    
}


/*
 * given colors available for horizontal and vertical edges
 * generate at least 2 * hc * vc color coded tile models and put them into tiles_
 * models will contain all possible combinations of NW edges
 * colors for SE edges will be randomly selected
 */
void Tiles::initTiles(){
    
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

double Tiles::colorDiff(double r, double g, double b){
    double res = 0.0;
    res = sqrt(r*r + g*g + b*b);
    return res;
}

vector<int> Tiles::seamCarving(Image* src1, Image* src2){
    vector<int> res;
    int width = src1->getWidth();
    int height = src1->getHeight();
    
    // calculate pixel color difference
    Image* temp = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = abs(src1->getPixel_(w,h,0) - src2->getPixel_(w,h,0));
            double g = abs(src1->getPixel_(w,h,1) - src2->getPixel_(w,h,1));
            double b = abs(src1->getPixel_(w,h,2) - src2->getPixel_(w,h,2));
            temp->setPixel_(w, h, 0, r);
            temp->setPixel_(w, h, 1, g);
            temp->setPixel_(w, h, 2, b);
        }
    }
    
    // create dp table
    vector<vector<double>> dp_table;
    for (int i=0; i<height; i++) {
        vector<double> temp_vec(width);
        dp_table.push_back(temp_vec);
    }
    vector<vector<double>> bt_table;
    for (int i=0; i<height; i++) {
        vector<double> temp_vec(width);
        bt_table.push_back(temp_vec);
    }
//    double **dp_table = new double*[height];
//    for (int x=0; x<height; x++){
//        dp_table[x] = new double[width];
//    }
//    double **bt_table = new double*[height];
//    for (int x=0; x<height; x++){
//        bt_table[x] = new double[width];
//    }
    
    // initialize the frst row
    for (int i = 0; i < width; i++) {
        double r = temp->getPixel_(i,0,0);
        double g = temp->getPixel_(i,0,1);
        double b = temp->getPixel_(i,0,2);
        double temp_c = sqrt(r*r + g*g + b*b);
        dp_table[0][i] = temp_c;
    }
    
    for(int i=1; i<height; i++){
        for(int j=0; j<width; j++){
            double r = temp->getPixel_(j,i,0);
            double g = temp->getPixel_(j,i,1);
            double b = temp->getPixel_(j,i,2);
            
            // on the left
            if (j == 0) {
                double r2 = temp->getPixel_(j,i-1,0);
                double g2 = temp->getPixel_(j,i-1,1);
                double b2 = temp->getPixel_(j,i-1,2);
                double r3 = temp->getPixel_(j+1,i-1,0);
                double g3 = temp->getPixel_(j+1,i-1,1);
                double b3 = temp->getPixel_(j+1,i-1,2);
                dp_table[i][j] = colorDiff(r,g,b) + min(colorDiff(r2,g2,b2),colorDiff(r3,g3,b3));
                if (colorDiff(r2,g2,b2)<colorDiff(r3,g3,b3)) {
                    bt_table[i][j] = j;
                } else {
                    bt_table[i][j] = j+1;
                }
            } else if (j > 0 && j < width-1) {  // middle
                double r1 = temp->getPixel_(j-1,i-1,0);
                double g1 = temp->getPixel_(j-1,i-1,1);
                double b1 = temp->getPixel_(j-1,i-1,2);
                double r2 = temp->getPixel_(j,i-1,0);
                double g2 = temp->getPixel_(j,i-1,1);
                double b2 = temp->getPixel_(j,i-1,2);
                double r3 = temp->getPixel_(j+1,i-1,0);
                double g3 = temp->getPixel_(j+1,i-1,1);
                double b3 = temp->getPixel_(j+1,i-1,2);
                dp_table[i][j] = colorDiff(r,g,b) + min(min(colorDiff(r1,g1,b1),colorDiff(r2,g2,b2)),colorDiff(r3,g3,b3));
                if (colorDiff(r1,g1,b1)<colorDiff(r2,g2,b2)) {
                    if (colorDiff(r1,g1,b1)<colorDiff(r3,g3,b3)) {
                        bt_table[i][j] = j-1;
                    } else {
                        bt_table[i][j] = j+1;
                    }
                } else {
                    if (colorDiff(r2,g2,b2)<colorDiff(r3,g3,b3)) {
                        bt_table[i][j] = j;
                    } else {
                        bt_table[i][j] = j+1;
                    }
                }
            } else {  // on the right
                double r1 = temp->getPixel_(j-1,i-1,0);
                double g1 = temp->getPixel_(j-1,i-1,1);
                double b1 = temp->getPixel_(j-1,i-1,2);
                double r2 = temp->getPixel_(i-1,j,0);
                double g2 = temp->getPixel_(i-1,j,1);
                double b2 = temp->getPixel_(i-1,j,2);
                dp_table[i][j] = colorDiff(r,g,b) + min(colorDiff(r1,g1,b1),colorDiff(r2,g2,b2));
                if (colorDiff(r1,g1,b1)<colorDiff(r2,g2,b2)) {
                    bt_table[i][j] = j-1;
                } else {
                    bt_table[i][j] = j;
                }
            }
        }
    }
    
    int min_bt = dp_table[height-1][0];
    for (int i=1; i<width; i++) {
        if (min_bt >= dp_table[height-1][i]) {
            min_bt = dp_table[height-1][i];
        }
    }
    
    res.push_back(min_bt);
    for (int i=height-1; i>0; i--) {
        res.push_back(bt_table[i][min_bt]);
        min_bt = bt_table[i][min_bt];
    }
    
    reverse(res.begin(),res.end());
    
    
//    
//    // TODELETE: dummy
//    for (int i = 0; i < src1->getHeight(); i++){
//        res.push_back(0);
//    }
    
    return res;
}

Image* Tiles::testSeamCarving(){
    
    vector<int> h_carve = seamCarving(himage_[0], himage_[1]);
    
    Image* pic  = new Image(himage_[0]->getWidth(), himage_[0]->getHeight());
    for (int i = 0; i < himage_[0]->getWidth(); i++){
        for (int j = 0; j < himage_[0]->getHeight(); j++){
            
            Pixel p;
            if (i < h_carve[j]){
                p = himage_[0] -> getPixel(i, j);
            } else {
                p = himage_[1] -> getPixel(i, j);
            }
            
            pic->setPixel(i, j, p);
        }
    }
    
    
    return pic;
}


/*
 * given desired width w and desired height h,
 * generate a tiled plain of that size
 */
Image* Tiles::tilePlain(int w, int h)
{
//    tw_ = tw_*sqrt(2);
//    th_ = th_*sqrt(2);
    
    int rth = tiles_[0].getTexture()->getHeight();
    int rtw = tiles_[0].getTexture()->getWidth();
    
    
    Image* dest = new Image(w, h);
    int nh = w/rtw + 1;                      // number of tiles horizontally
    int nv = h/rth + 1;                      // number of tiles vertically
    int tileModel[nv][nh];                   // a 2D array that specify types of tiles used
    
    tileModel[0][0] = rand() % tiles_.size();
    
    // generate an abstract valid tiling
    for (int i = 0; i < nv; i++){
        for (int j = 0; j < nh; j++){
            
            int up = (i == 0) ? -1 : tiles_[tileModel[i-1][j]].getDown();
            int left = (j == 0) ? -1 : tiles_[tileModel[i][j-1]].getRight();
            
            tileModel[i][j] = getRandomTile(up, left);
            
        }
    }
    
    // fill in the corresponding textured tile according to the tile models
    for (int i = 0; i < h; i ++){
        for (int j = 0; j < w; j++){
            
            int tileIndex = tileModel[i / rth][j / rtw];
            Tile t = tiles_[tileIndex];
            Pixel p = t.getTexture()->getPixel(j % rtw,i % rth);
            dest->setPixel(j, i, p);
        }
    }
    
    return dest;
    //return genTextures(1,1,2,2);
}