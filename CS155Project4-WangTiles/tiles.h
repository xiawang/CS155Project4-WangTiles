#include "image.h"

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
        Image* getTexture();
        ~Tile ();
        
        
        
    private:
        int up_;
        int down_;
        int left_;
        int right_;
        Image* tex_;
};