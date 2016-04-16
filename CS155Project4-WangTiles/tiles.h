#include "image.h"
#include <vector>

using namespace std;

class Tile
{
    public:
        Tile ();
        
        void setUp(int u);
        void setDown(int d);
        void setLeft(int l);
        void setRight(int r);
        void setTexture(Image* src);
        int getUp();
        int getDown();
        int getRight();
        int getLeft();
        bool validNeighbor(int type, Tile n);
        Image* getTexture();
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
        Image* tilePlain();

    private:
        vector<Tile> tiles_;
        void genTiles();



};