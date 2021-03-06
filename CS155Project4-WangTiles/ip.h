#ifndef IP_H
#define IP_H


#include "common.h"
#include "image.h"


/*
 * IMPORTANT - DO NOT CHANGE THE INTERFACES DEFINED HERE - IMPORTANT
 *
 * It's very important that the interface to your program does not
 * change from what is provided, so that automated testing scripts
 * will work.  If you add additional requests for input, these scripts
 * will no longer work and it will negatively affect your grade.  Each
 * method has sufficient information to carry out the required task.
 *
 * The only method you are may get more user input in is ip_warp().
 * This option will not be tested automatically, since everyone's
 * will be different, so you may ask for whatever input is necessary.
 * To support multiple warps, have ip_warp() print a menu of options.
 *
 * Of course, you may add whatever other functions you may need.
 */


/*
 * definitions for the required sampling techniques
 */
#define I_NEAREST	0
#define I_BILINEAR	1
#define I_GAUSSIAN	2

Image*  ip_tile (Image* src, int hc, int vc, int w, int h, int tw, int th, bool source);
Image* ip_quilt (Image* src, int patch_size, double patch_w, double patch_h, int num_h, int num_w);
Image*	ip_blur_box (Image* src, int size);
Image*	ip_blur_gaussian (Image* src, int size, double sigma);
Image*	ip_blur_triangle (Image* src, int size);

Image*	ip_brighten (Image* src, double alpha);
Image*	ip_color_shift (Image* src);
Image*	ip_contrast (Image* src, double alpha);
Image*	ip_composite (Image* src1, Image* src2, Image* mask);
Image*	ip_crop (Image* src, int x0, int y0, int x1, int y1);
Image*	ip_edge_detect (Image* src);
Image*	ip_extract (Image* src, int channel);
Image*	ip_grey (Image* src);
Image*  ip_image_shift (Image* src, int dx, int dy);
Image*	ip_invert (Image* src);
double	ip_kernel_calc(Image* src, int i, int j, int k, double* kernel, int size);
Image*  ip_median(Image* src, int n);
Image*	ip_misc(Image* src);
Image*	ip_noisify (Image* src, double alpha);

Image*	ip_quantize_simple (Image* src, int bitsPerChannel);
Image*	ip_quantize_ordered (Image* src, int bitsPerChannel);
Image*	ip_quantize_fs (Image* src, int bitsPerChannel);

Pixel	ip_resample_nearest(Image* src, double x, double y);
double	ip_resample_nearest_channel(Image* src, double x, double y, int channel);
Pixel	ip_resample_bilinear(Image* src, double x, double y);
double	ip_resample_bilinear_channel(Image* src, double x, double y, int channel);
Pixel   ip_resample_gaussian(Image* src, double x, double y, int filtersize, double sigma);
double   ip_resample_gaussian_channel(Image* src, double x, double y,int filtersize, double sigma,  int channel);


Image*	ip_rotate (Image* src, double theta, int x, int y, int mode,
                   int size, double sigma);
Image*	ip_saturate (Image* src, double alpha);
Image*	ip_sepia (Image* src);
Image*	ip_scale (Image* src, double x, double y, int mode, int size,
                  double sigma);
Image*	ip_threshold (Image* src, double cutoff);
Image*	ip_fun_warp (Image* src, int samplingMethod);
Image*  ip_interpolate(Image* src1, Image* src2, double alpha);
Image* ip_convolve(Image* src, int size, double** kernel);
Image* ip_misc(Image* src);



#endif // IP_H
