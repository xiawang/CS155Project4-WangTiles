#include "image.h"
#include <vector>

using namespace std;

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
        void setTexture(Image* src);
        int getUp();
        int getDown();
        int getRight();
        int getLeft();
        Image* getTexture();
    
        bool validNeighbor(int type, Tile n);
        void genTexture();
        ~Tile ();
        
        
        
    private:
        int up_;
        int down_;
        int left_;
        int right_;
        Image* tex_;
};


class Tiles
{
    public:
        Tiles(Image* scr, int hc, int vc);
        Image* tilePlain(int w, int h); // generete FINAL RESULT!

    private:
        friend class Tile;
        vector<Tile> tiles_;
        vector<Image*> himage_;
        vector<Image*> vimage_;
        void genTiles(int hc, int vc);   // generate tile MODELS, bet on probability



};