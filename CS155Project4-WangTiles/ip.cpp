#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "tiles.h"


/*
 * tiling a plain using Wang tiles
 */
Image* ip_tile (Image* src, int hc, int vc, int w, int h, bool source)
{
    cout << "hc is: "  << hc << endl;
    cout << "vc is: "  << vc << endl;
    Tiles* t  = new Tiles();
    return t->tilePlain(w, h);

}



/*
 * convolve with a box filter
 */
Image* ip_blur_box (Image* src, int size)
{
    double **kernel = new double*[size];
    for (int x=0; x<size; x++){
        kernel[x] = new double[size];
    }
    for(int i=0;i<size;i++){
        for(int j=0; j<size; j++){
            kernel[i][j] = 1.0/(size*size);
        }
    }
    return ip_convolve(src, size, kernel);
}

/*
 * convolve with a gaussian filter
 */
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    double **kernel = new double*[size];
    for (int x=0; x<size; x++){
        kernel[x] = new double[size];
    }
    double sum = 0.0;
    for(int i=0;i<size;i++){
        for(int j=0; j<size; j++){
            int span = size / 2;
            int a = i - span;
            int b = j - span;
            double pi = atan(1.0)*4;
            double n = -(a*a + b*b) / (2.0 * sigma * sigma);
            double g_value = (1 / sqrt(2.0 * pi * sigma * sigma)) * exp(n);
            kernel[i][j] = g_value;
            sum += kernel[i][j];
        }
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            kernel[i][j] /= sum;
        }
    }
    return ip_convolve(src, size, kernel);
}


/*
 * convolve with a triangle filter
 */
Image* ip_blur_triangle (Image* src, int size)
{
	
    cerr << "This filter is no longer in the assignment.\n";
    return NULL;
}


/*
 * interpolate with a black image
 */
Image* ip_brighten (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* black = new Image(width,height);
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            black->setPixel_(w, h, 0, 0);
            black->setPixel_(w, h, 1, 0);
            black->setPixel_(w, h, 2, 0);
        }
    }
    dest = ip_interpolate (src, black, alpha);
    return dest;
}

Image* ip_color_shift(Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            dest->setPixel(w, h, 0, b);
            dest->setPixel(w, h, 1, r);
            dest->setPixel(w, h, 2, g);
        }
    }
    return dest;
}

/*
 * interpolate with the average intensity of the src image
 */
Image* ip_contrast (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* grey = new Image(width,height);
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            grey->setPixel_(w, h, 0, 0.5);
            grey->setPixel_(w, h, 1, 0.5);
            grey->setPixel_(w, h, 2, 0.5);
        }
    }
    dest = ip_interpolate (src, grey, alpha);
    return dest;
}


/*
 * convolve an image with another image
 */
Image* ip_convolve (Image* src, int size, double** kernel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image(width,height);
    int hspan = (size-1)/2;
    int wspan = (size-1)/2;
    
    for (int w=0; w<width;w++){
        for (int h=0; h<height; h++){
            double r = 0; double g = 0; double b = 0;
            for(int dw=-wspan; dw<=wspan; dw++){
                for(int dh=-hspan; dh<=hspan; dh++){
                    int ww = w + dw; int hh = h + dh;
                        double weight = kernel[dw+wspan][dh+hspan];
                        double rr = src->getPixel_(ww,hh,0);
                        double gg = src->getPixel_(ww,hh,1);
                        double bb = src->getPixel_(ww,hh,2);
                        r = r + weight * rr;
                        g = g + weight * gg;
                        b = b + weight * bb;
                }
            }
            output->setPixel_(w,h,0,r);
            output->setPixel_(w,h,1,g);
            output->setPixel_(w,h,2,b);
        }
    }
    return output;
}


/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
Image* ip_composite(Image* src1, Image* src2, Image* mask)
{
    int width = src1->getWidth();
    int height = src1->getHeight();
    int width1 = src2->getWidth();
    int height1 = src2->getHeight();
    int width2 = mask->getWidth();
    int height2 = mask->getHeight();
    cout << "src1_w: " << width << " src1_h: " << height << endl;
    cout << "src2_w: " << width1 << " src2_h: " << height1 << endl;
    cout << "mask_w: " << width2 << " mask_h: " << height2 << endl;
    Image* dest = new Image(width1,height1);
    // background layer
    for (int w=0; w<width1; w++){
        for (int h=0;h<height1; h++){
            double r2 = src2->getPixel(w, h, 0);
            double g2 = src2->getPixel(w, h, 1);
            double b2 = src2->getPixel(w, h, 2);
            dest->setPixel_(w, h, 0, r2);
            dest->setPixel_(w, h, 1, g2);
            dest->setPixel_(w, h, 2, b2);
        }
    }
    // face layer
    int x;
    int y;
    cout << "face_x: " << endl;
    cin >> x;
    cout << "face_y: " << endl;
    cin >> y;
    
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r1 = src1->getPixel(w, h, 0);
            double g1 = src1->getPixel(w, h, 1);
            double b1 = src1->getPixel(w, h, 2);
            double r2 = src2->getPixel(w+x, h+y, 0);
            double g2 = src2->getPixel(w+x, h+y, 1);
            double b2 = src2->getPixel(w+x, h+y, 2);
            double alpha_r = mask->getPixel(w, h, 0);
            double alpha_g = mask->getPixel(w, h, 1);
            double alpha_b = mask->getPixel(w, h, 2);
            double rd = alpha_r*r1 + (1-alpha_r)*r2;
            double gd = alpha_g*g1 + (1-alpha_g)*g2;
            double bd = alpha_b*b1 + (1-alpha_b)*b2;
            dest->setPixel_(w+x, h+y, 0, rd);
            dest->setPixel_(w+x, h+y, 1, gd);
            dest->setPixel_(w+x, h+y, 2, bd);
        }
    }

    return dest;
}


/*
 * cut away all but a subset of the image
 */
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    Image* dest = new Image((x1-x0),(y1-y0));
    for (int w=x0; w<x1; w++){
        for (int h=y0;h<y1; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            dest->setPixel_(w-x0, h-y0, 0, r);
            dest->setPixel_(w-x0, h-y0, 1, g);
            dest->setPixel_(w-x0, h-y0, 2, b);
        }
    }
    return dest;
}


/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
    double **kernel = new double*[3];
    for (int x=0; x<3; x++){
        kernel[x] = new double[3];
    }
    for(int i=0;i<3;i++){
        for(int j=0; j<3; j++){
            kernel[i][j] = (i == 1 && j == 1) ? 8 : -1;
        }
    }
    return ip_convolve(src, 3, kernel);
}


/*
 * create a new image with color values from one channel of src
 */
Image* ip_extract (Image* src, int channel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0; h<height; h++){
            for (int c=0; c<3; c++) {
                if (c != channel) {
                    dest->setPixel(w, h, c, 0);
                } else {
                    double cvalue = src->getPixel(w, h, channel);
                    dest->setPixel(w, h, c, cvalue);
                }
            }
        }
    }
    return dest;
}

/*
 * perform do what you like
 * you must query user for parameters here
 */
Image* ip_fun_warp (Image* src,int samplingMethod)
{
    int x = 0;
    int y = 0;
    int r = 0;
    int w = 0;
    float t = 0;
    float s = 0;
    int size = 0;
    double sigma = 0.0;
    cout << "Ripple effect x value: " << endl;
    cin >> x;
    cout << "Ripple effect y value: " << endl;
    cin >> y;
    cout << "Ripple effect radius value: " << endl;
    cin >> r;
    cout << "Ripple effect wavelength value: " << endl;
    cin >> w;
    cout << "Ripple effect trainwidth value: " << endl;
    cin >> t;
    cout << "Ripple effect superphase value: " << endl;
    cin >> s;
    cout << "Size: " << endl;
    cin >> size;
    cout << "sigma: " << endl;
    cin >> sigma;
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    // initialize as all black
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            int new_x = 0;
            int new_y = 0;
            double dx = w - x;
            double dy = h - y;
            double rr = (sqrt(dx*dx+dy*dy)-r)/w ;
            double k = rr - (1-s)*r/w ;
            double a = 1 / (1.0 + (rr/t)*(rr/t));
            double depth = a * sin(k*2*M_PI);
            new_x = w + depth*20;
            new_y = h + depth*20;
            // resemble
            Pixel pix;
            if (new_x < 0 || new_x >= width || new_y < 0 || new_y >= height) {
                pix = Pixel(0, 0, 0);
            } else {
                if (samplingMethod == I_NEAREST) {
                    pix = ip_resample_nearest(src, new_x, new_y);
                } else if (samplingMethod == I_BILINEAR) {
                    pix = ip_resample_bilinear(src, new_x, new_y);
                } else {
                    pix = ip_resample_gaussian(src, new_x, new_y, size, sigma);
                }
            }
            
            dest->setPixel_(w, h, pix);
        }
    }
    return dest;
}


/*
 * create a new image with values equal to the psychosomatic intensities
 * of the source image
 */
Image* ip_grey (Image* src)
{
    double cr=  0.2126;
    double cg=  0.7152;
    double cb=  0.0722;
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double gvalue = cr*r + cg*g + cb*b;
            dest->setPixel(w, h, 0, gvalue);
            dest->setPixel(w, h, 1, gvalue);
            dest->setPixel(w, h, 2, gvalue);
        }
    }
    return dest;
}

/*
 * shift image by dx, dy
 *
*/
Image* ip_image_shift(Image* src, int dx, int dy)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            int nw = (w+dx) % width;
            int nh = (h+dy) % height;
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            dest->setPixel(nw, nh, 0, r);
            dest->setPixel(nw, nh, 1, g);
            dest->setPixel(nw, nh, 2, b);
        }
    }
    return dest;
}

/*
 * interpolate an image with another image
 */
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    int width = src1->getWidth();
    int height = src1->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r1 = src1->getPixel(w, h, 0);
            double g1 = src1->getPixel(w, h, 1);
            double b1 = src1->getPixel(w, h, 2);
            double r2 = src2->getPixel(w, h, 0);
            double g2 = src2->getPixel(w, h, 1);
            double b2 = src2->getPixel(w, h, 2);
            double rd = alpha*r1 + (1-alpha)*r2;
            double gd = alpha*g1 + (1-alpha)*g2;
            double bd = alpha*b1 + (1-alpha)*b2;
            dest->setPixel_(w, h, 0, rd);
            dest->setPixel_(w, h, 1, gd);
            dest->setPixel_(w, h, 2, bd);
        }
    }
    return dest;
}

/*
 * subtract the image from a white image
 */
Image* ip_invert (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    Image* grey = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            grey->setPixel_(w, h, 0, 0.5);
            grey->setPixel_(w, h, 1, 0.5);
            grey->setPixel_(w, h, 2, 0.5);
        }
    }
    dest = ip_interpolate (src, grey, -1);
    return dest;
}

/*
 * median filter
 */
Image*  ip_median(Image* src, int n)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width-2; w++){
        for (int h=0;h<height-2; h++){
            // construct vector for sorting
            vector<double> vec_r;
            vector<double> vec_g;
            vector<double> vec_b;
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    // get r,g,b values
                    double r = src->getPixel(w+i, h+j, 0);
                    double g = src->getPixel(w+i, h+j, 1);
                    double b = src->getPixel(w+i, h+j, 2);
                    vec_r.push_back(r);
                    vec_g.push_back(g);
                    vec_b.push_back(b);
                }
            }
            // sorting r,g,b vectors
            sort(vec_r.begin(), vec_r.end());
            sort(vec_g.begin(), vec_g.end());
            sort(vec_b.begin(), vec_b.end());
            dest->setPixel(w, h, 0, vec_r[4]);
            dest->setPixel(w, h, 1, vec_g[4]);
            dest->setPixel(w, h, 2, vec_b[4]);
        }
    }
    return dest;
}

/* 
 * misc - We implemented as the "wind effect"
 *        This method requires the input for intensity
 */
Image* ip_misc(Image* src) 
{
    int intensity = 0;
    cout << "Please input the intensity of the wind: " << endl;
    cin >> intensity;
    int width = src->getWidth();
    int height = src->getHeight();
    Image* newImage = new Image(width,height);
    for (int j = 0; j<height; j++){
        int max = j + intensity;
        int min = j - intensity;
        int randj = rand()%(max-min + 1) + min;
        if (randj<0){
            randj = 0;
        }else if(randj>=height){
            randj = height - 1;
        }
        for (int i=0; i<width; i++){
            Pixel pixel = src->getPixel(i,randj);
            newImage->setPixel_(i, j, pixel);
        }
    }
    return newImage;
}


/*
 * round each pixel to the nearest value in the new number of bits
 */

//round a number to the closest number that can be represented by bit values.
double bitRound(double number, int bit){
    int levels = pow(2,bit);
    double grid = 1.0/double(levels-1);
    double val;
    for (double cur=0.0; cur<=1.0; cur+=grid){
        if (cur >= 1.0 || fabs(cur-number)<fabs(cur+grid-number)) {
            val = cur;
            break;
        }
    }
    return val;
}

Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double rp = bitRound(r, bitsPerChannel);
            double gp = bitRound(g, bitsPerChannel);
            double bp = bitRound(b, bitsPerChannel);
            dest->setPixel_(w, h, 0, rp);
            dest->setPixel_(w, h, 1, gp);
            dest->setPixel_(w, h, 2, bp);
        }
    }
    return dest;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using a static 2x2 matrix
 */

//round a number to the closest number that can be represented by a 2x2 block.
double blockRound(double number, int bit){
    int levels = 4*(-1+pow(2,bit))+1;
    double grid = 1.0/double(levels-1);
    double val;
    for (double cur=0.0; cur<=1.0; cur+=grid){
        if (cur >= 1.0 || fabs(cur-number)<fabs(cur+grid-number)) {
            val = cur;
            break;
        }
    }
    return val;
}

Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    int levels = pow(2,bitsPerChannel);
    double grid = 1.0/double(levels-1);
    double offset0 = grid*3.0/8.0;
    double offset1 = grid*1.0/8.0;
    double offset2 = -grid*1.0/8.0;
    double offset3 = -grid*3.0/8.0;
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            if (w%2==0 && h%2 ==1){
                r += offset0; g += offset0; b += offset0;
            }else if(w%2==1 && h%2 ==0){
                r += offset1; g += offset1; b += offset1;
            }else if(w%2==1 && h%2 ==1){
                r += offset2; g += offset2; b += offset2;
            }else{
                r += offset3; g += offset3; b += offset3;
            }
            double rp = bitRound(r, bitsPerChannel);
            double gp = bitRound(g, bitsPerChannel);
            double bp = bitRound(b, bitsPerChannel);
            dest->setPixel_(w, h, 0, rp);
            dest->setPixel_(w, h, 1, gp);
            dest->setPixel_(w, h, 2, bp);
        }
    }
    return dest;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using error diffusion
 */
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    // floyd-steinberg factors
    float alpha = 7.0 / 16.0;
    float beta = 3.0 / 16.0;
    float chi = 5.0 / 16.0;
    float delta = 1.0 / 16.0;
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    
    // construct and initialize error float matrix
    float e_matrix[height][width][3];
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            e_matrix[i][j][0] = 0.0;
            e_matrix[i][j][1] = 0.0;
            e_matrix[i][j][2] = 0.0;
        }
    }
    
    for (int h=0;h<height; h++){
        for (int w=0; w<width; w++){
            // read in errors
            float error_r = e_matrix[h][w][0];
            float error_g = e_matrix[h][w][1];
            float error_b = e_matrix[h][w][2];
            
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double nr = blockRound(r+error_r,bitsPerChannel);
            double ng = blockRound(g+error_g,bitsPerChannel);
            double nb = blockRound(b+error_b,bitsPerChannel);
            dest->setPixel_(w, h, 0, nr);
            dest->setPixel_(w, h, 1, ng);
            dest->setPixel_(w, h, 2, nb);
            
            // store errors
            double cur_error_r = nr - r;
            double cur_error_g = ng - g;
            double cur_error_b = nb - b;
            
            // alpha
            if (w == width-1) {  // right most corner
                // do nothing
            } else {
                e_matrix[h][w][0] += alpha*cur_error_r;
                e_matrix[h][w][1] += alpha*cur_error_g;
                e_matrix[h][w][2] += alpha*cur_error_b;
            }
            
            // beta
            if (w == 0 || h == height-1) {  // lower left most corner
                // do nothing
            } else {
                e_matrix[h][w][0] += beta*cur_error_r;
                e_matrix[h][w][1] += beta*cur_error_g;
                e_matrix[h][w][2] += beta*cur_error_b;
            }
            
            // chi
            if (h == height-1) {  // bottom line
                // do nothing
            } else {
                e_matrix[h][w][0] += chi*cur_error_r;
                e_matrix[h][w][1] += chi*cur_error_g;
                e_matrix[h][w][2] += chi*cur_error_b;
            }
            
            // delta
            if (w == width-1 || h == height-1) {  // lower right most corner
                // do nothing
            } else {
                e_matrix[h][w][0] += delta*cur_error_r;
                e_matrix[h][w][1] += delta*cur_error_g;
                e_matrix[h][w][2] += delta*cur_error_b;
            }
        }
    }
    return dest;
}

/* helper functions you may find useful for resampling */

/*
 * nearest neighbor sample
 */

Pixel ip_resample_nearest(Image* src, double x, double y)
{
    double x_floor = floor(x);
    double y_floor = floor(y);
    double dist_1 = sqrt(pow((x-x_floor), 2.0) + pow((y-y_floor), 2.0));
    double dist_2 = sqrt(pow((x-x_floor-1), 2.0) + pow((y-y_floor), 2.0));
    double dist_3 = sqrt(pow((x-x_floor), 2.0) + pow((y-y_floor-1), 2.0));
    double dist_4 = sqrt(pow((x-x_floor-1), 2.0) + pow((y-y_floor-1), 2.0));
    vector<double> vec;
    vec.push_back(dist_1);
    vec.push_back(dist_2);
    vec.push_back(dist_3);
    vec.push_back(dist_4);
    double index = min_element( vec.begin(), vec.end() ) - vec.begin();
    Pixel my_pix;
    if (index == 0) {
        my_pix = src->getPixel_(x_floor, y_floor);
    } else if (index == 1) {
        my_pix = src->getPixel_(x_floor+1, y_floor);
    } else if (index == 2) {
        my_pix = src->getPixel_(x_floor, y_floor+1);
    } else {
        my_pix = src->getPixel_(x_floor+1, y_floor+1);
    }
    return my_pix;
}

/*
 * bilinear resample
 */

Pixel takeWeightedAverage(Pixel pixel1, Pixel pixel2, double weight1, double weight2){
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

Pixel ip_resample_bilinear(Image* src, double x, double y){
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
 * gausian sample
 */
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma) {
    int width = src->getWidth();
    int height = src->getHeight();
    double **kernel = new double*[size];
    for (int x=0; x<size; x++){
        kernel[x] = new double[size];
    }
    int x0 = floor(x)-(size/2-1); int y0 = floor(y)-(size/2-1);
    double sum = 0;
    for(int i=x0; i<x0+size; i++){
        for(int j=y0; j<y0+size;j++){
            kernel[i-x0][j-y0] = exp(-(pow(x-i,2)+pow(y-j,2))/(2*pow(sigma,2)));
            sum = sum + kernel[i-x0][j-y0];
        }
    }
    double r = 0; double g = 0; double b = 0;
    for(int i=x0; i<x0+size; i++){
        for(int j=y0; j<y0+size;j++){
            if(i>0 && i<width && j>0 && j<height){
                double weight = kernel[i-x0][j-y0]/sum;
                r = r + weight*(src->getPixel(i, j, 0));
                g = g + weight*(src->getPixel(i, j, 1));
                b = b + weight*(src->getPixel(i, j, 2));
            }
        }
    }
    return Pixel(r,g,b);
}


/*
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta, int x, int y, int mode, 
                  int size, double sigma)
{
    float sin_val = sin(theta * M_PI / 180);
    float cos_val = cos(theta * M_PI / 180);
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
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
                if (mode == I_NEAREST) {
                    pix = ip_resample_nearest(src, rot_x, rot_y);
                } else if (mode == I_BILINEAR) {
                    pix = ip_resample_bilinear(src, rot_x, rot_y);
                } else {
                    pix = ip_resample_gaussian(src, rot_x, rot_y, size, sigma);
                }
            }
            
            dest->setPixel_(w, h, pix);
        }
    }
    return dest;
}

/*
 * interpolate with the greyscale version of the image
 */
Image* ip_saturate (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* greyscale = ip_grey(src);
    Image* dest = new Image(width,height);
    dest = ip_interpolate (src, greyscale, alpha);
    return dest;
}


/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac, int mode, 
                 int size, double sigma)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width*xFac,height*yFac);
    // initialize as all black
    for (int w=0; w<width*xFac; w++){
        for (int h=0;h<height*yFac; h++){
            // rescaling
            int new_x = w * 1.0 / xFac;
            int new_y = h * 1.0 / yFac;
            // resemble
            Pixel pix;
            if (new_x < 0 || new_x >= width || new_y < 0 || new_y >= height) {
                pix = Pixel(0, 0, 0);
            } else {
                if (mode == I_NEAREST) {
                    pix = ip_resample_nearest(src, new_x, new_y);
                } else if (mode == I_BILINEAR) {
                    pix = ip_resample_bilinear(src, new_x, new_y);
                } else {
                    pix = ip_resample_gaussian(src, new_x, new_y, size, sigma);
                }
            }
            
            dest->setPixel_(w, h, pix);
        }
    }
    return dest;
}

/*
 * create a new sepia tones image
 */
Image* ip_sepia (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double rp = 0.393*r + 0.769*g + 0.189*b;
            double gp = 0.349*r + 0.686*g + 0.168*b;
            double bp = 0.272*r + 0.534*g + 0.131*b;
            dest->setPixel_(w, h, 0, rp);
            dest->setPixel_(w, h, 1, gp);
            dest->setPixel_(w, h, 2, bp);
        }
    }
    return dest;
}


/*
 * create a new one bit/channel image with the intensity 0 or 1
 * depending on whether the input value is above or below the 
 * threshold
 */
Image* ip_threshold (Image* src, double cutoff)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double rp = (r < cutoff) ? 0 : 1;
            double gp = (g < cutoff) ? 0 : 1;
            double bp = (b < cutoff) ? 0 : 1;
            dest->setPixel(w, h, 0, rp);
            dest->setPixel(w, h, 1, gp);
            dest->setPixel(w, h, 2, bp);
        }
    }
    return dest;
}




