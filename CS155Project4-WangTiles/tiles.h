#include "image.h"
#include <vector>

using namespace std;

const int tileWidth = 50;
const int tileHeight = 50;

class Tile
{
    public:
        Tile ();
        Tile (int up, int left);
    
        // gets and sets
        void setUp(int u);
        void setDown(int d);
        void setLeft(int l);
        void setRight(int r);
        void setSize(int s);
        void setTexture(Image* src);
        int getUp();
        int getDown();
        int getRight();
        int getLeft();
        int getSize();
        Image* getTexture();
    
        bool validNeighbor(int up, int left);
        void genTexture();
        ~Tile ();
        
        
        
    private:
        int up_;
        int down_;
        int left_;
        int right_;
        int tilesize_;
        Image* tex_;
};


class Tiles
{
    public:
        Tiles(int tw, int th);   // default constructor for testing, using dummy models in tiles folder
                   // Ignores the source image, but uses specified width and height
        Tiles(Image* scr, int hc, int vc);
        Image* tilePlain(int w, int h, int tw, int th); // generete FINAL RESULT!

    private:
        friend class Tile;
        vector<Tile> tiles_;     // up_ = 1 ==> corresponds to himage_.get(1)
                                 // left = 1 ==> corresponds to vimage_.get(1)
        vector<Image*> himage_;               // each image in himage_ corresponds to a horizontal color
        vector<Image*> vimage_;               // each image in vimage_ corresponds to a vertical color
        void genTiles(int hc, int vc);        // generate tile MODELS, bet on probability
        int getRandomTile(int up, int left);  // randomly select a tile with matching NW edge
        Image* genDummyTile(int n, int e, int s, int w, int tw, int th);



};