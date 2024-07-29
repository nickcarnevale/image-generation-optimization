#include "image.h"
#include <time.h>
using ImagePtr = std::shared_ptr<Image>;

ImagePtr applyGamma(ImagePtr image_ptr, double gamma);
ImagePtr applyTint(ImagePtr image_ptr, const double *tints);
ImagePtr applyBlur(ImagePtr imag_ptr);
void writeImage(ImagePtr image_ptr);

//MY OPTIMIZED FUNCTIONS
ImagePtr applyGammaReduceFunctionCalls(ImagePtr image_ptr, double gamma);
ImagePtr applyGammaUnrolled2x1(ImagePtr image_ptr, double gamma);
ImagePtr applyGammaUnrolled4x1(ImagePtr image_ptr, double gamma);
ImagePtr applyBlurUnrolled2x1(ImagePtr image_ptr);
ImagePtr applyBlurUnrolled2x1withInnerLoop(ImagePtr image_ptr);
ImagePtr applyTintUnrolled2x1(ImagePtr image_ptr, const double *tints);


void process_images(const std::vector<ImagePtr>& image_vector) {
  const double tint_array[] = {0.75, 0, 0};
  for (ImagePtr img : image_vector) {
    writeImage(img);
    img = applyGamma(img, 1.4); 
    img = applyTint(img, tint_array);
    img = applyBlur(img);
    writeImage(img);
  }
  
}

ImagePtr applyGamma(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", 
                                                  IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  //const int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      //const Pixel& p = in_rows[i][j]; 
      double v = 0.3*in_rows[i][j].bgra[2] + 0.59*in_rows[i][j].bgra[1] + 0.11*in_rows[i][j].bgra[0];
      double res = pow(v, gamma);
      if(res > MAX_BGR_VALUE) res = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res, res, res);
    }
  }
  return output_image_ptr;
}

ImagePtr applyTint(ImagePtr image_ptr, const double *tints) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_tinted", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();

  for (unsigned long i = 0; i < image_ptr->rows().size(); ++i ) {
    for (int j = 0; j < image_ptr->rows()[1] - image_ptr->rows()[0]; ++j ) {
      double b = (double)in_rows[i][j].bgra[0] + (MAX_BGR_VALUE-in_rows[i][j].bgra[0])*tints[0];
      double g = (double)in_rows[i][j].bgra[1] + (MAX_BGR_VALUE-in_rows[i][j].bgra[1])*tints[1];
      double r = (double)in_rows[i][j].bgra[2] + (MAX_BGR_VALUE-in_rows[i][j].bgra[0])*tints[2];
      out_rows[i][j].bgra[0] = b > MAX_BGR_VALUE ? MAX_BGR_VALUE:b;
      out_rows[i][j].bgra[1] = g > MAX_BGR_VALUE ? MAX_BGR_VALUE:g;
      out_rows[i][j].bgra[2] = r > MAX_BGR_VALUE ? MAX_BGR_VALUE:r;
    }
  }
  return output_image_ptr;
}

ImagePtr applyBlur(ImagePtr image_ptr) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_blurred", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  double b, g, r;
  //int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      // Average = ([i-1][j-1] + [i-1][j] + [i-1][j+1] + [i][j-1] + [i][j] + [i][j+1] + [i+1][j-1] + [i+1][j] + [i+1][j+1])/ 9
      if (i == 0) {                        /* first row */
        if (j == 0) {                     /* first row, first column */
          b = (0 + 0 + 0 + 0 + in_rows[i][j].bgra[0] + in_rows[i+1][j].bgra[0] + 0 + in_rows[i][j+1].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (0 + 0 + 0 + 0 + in_rows[i][j].bgra[1] + in_rows[i+1][j].bgra[1] + 0 + in_rows[i][j+1].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (0 + 0 + 0 + 0 + in_rows[i][j].bgra[2] + in_rows[i+1][j].bgra[2] + 0 + in_rows[i][j+1].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {          /* first row, last column */
          b = (0 + 0 + 0 + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + 0 + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + 0) / 9;
          g = (0 + 0 + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + 0 + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + 0) / 9;
          r = (0 + 0 + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + 0 + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + 0) / 9;
        } 
        else {                          /* first row, middle columns */
          b = (0 + 0 + 0 + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (0 + 0 + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (0 + 0 + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      } 
      else if (i == in_rows.size() - 1) {        /* last row */
        if (j == 0) {             /* last row, first column */
          b = (0 + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + 0 + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + 0 + 0 + 0) / 9;
          g = (0 + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + 0 + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + 0 + 0 + 0) / 9;
          r = (0 + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + 0 + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + 0 + 0 + 0) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* last row, last column */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j+1].bgra[0] + 0 + in_rows[i][i-1].bgra[0] + in_rows[i][j].bgra[0] + 0 + 0 + 0 + 0) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j+1].bgra[1] + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + 0 + 0 + 0 + 0) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j+1].bgra[2] + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + 0 + 0 + 0 + 0) / 9;
        } 
        else {                          /* last row, middle columns */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + 0 + 0 + 0) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + 0 + 0 + 0) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + 0 + 0 + 0) / 9;
        }
      } 
      else {                            /* middle rows */
        if (j == 0) {                 /* middle row, first column */
          b = ( 0 + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + 0 + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + 0 + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = ( 0 + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + 0 + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + 0 + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = ( 0 + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + 0 + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + 0 + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* middle row, last column */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + 0 + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + 0 + in_rows[i+1][j-1].bgra[0]+ in_rows[i+1][j].bgra[0] + 0) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + 0 + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + 0) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + 0 + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + 0) / 9;
        } 
        else {                          /* middle row, middle columns */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      }
      out_rows[i][j].bgra[0] = (b > MAX_BGR_VALUE)? MAX_BGR_VALUE : b;
      out_rows[i][j].bgra[1] = (g > MAX_BGR_VALUE)? MAX_BGR_VALUE : g;
      out_rows[i][j].bgra[2] = (r > MAX_BGR_VALUE)? MAX_BGR_VALUE : r;
    }
  }
  return output_image_ptr;
}

//MY OPTIMIZED FUNCTIONS
ImagePtr applyGammaReduceFunctionCalls(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", 
                                                  IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  const int height = in_rows.size();
  const int width = in_rows[1] - in_rows[0];

  const double redWeight = 0.11;
  const double greenWeight = 0.59;
  const double blueWeight = 0.3;
  
  for (unsigned long i = 0; i < height; ++i ) {
    for (int j = 0; j < width; ++j ) {
      const Pixel& p = in_rows[i][j]; 
      double v = blueWeight*p.bgra[2] + greenWeight*p.bgra[1] + redWeight*p.bgra[0];
      double res = pow(v, gamma);
      if(res > MAX_BGR_VALUE) res = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res, res, res);
    }
  }
  return output_image_ptr;
}

ImagePtr applyGammaUnrolled2x1(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  const int height = in_rows.size();
  const int width = in_rows[1] - in_rows[0];
  int limit = width -1;
  // Move constant luminance coefficients outside the loop
  // Code motion
  const double redWeight = 0.11;
  const double greenWeight = 0.59;
  const double blueWeight = 0.3;
  int i, j;
  for (i = 0; i < height; ++i) {
    for (j = 0; j < limit; j += 2) { // 2x1 loop unrolling
      Pixel p1 = in_rows[i][j];
      Pixel p2 = in_rows[i][j + 1];

      double v1 = blueWeight * p1.bgra[2] + greenWeight * p1.bgra[1] + redWeight * p1.bgra[0];
      double v2 = blueWeight * p2.bgra[2] + greenWeight * p2.bgra[1] + redWeight * p2.bgra[0];

      double res1 = pow(v1, gamma);
      double res2 = pow(v2, gamma);

      if (res1 > MAX_BGR_VALUE) res1 = MAX_BGR_VALUE;
      if (res2 > MAX_BGR_VALUE) res2 = MAX_BGR_VALUE;

      out_rows[i][j] = Pixel(res1, res1, res1);
      out_rows[i][j + 1] = Pixel(res2, res2, res2);
    }
    for(; j<width; j++){
      Pixel p1 = in_rows[i][j];
      double v1 = blueWeight * p1.bgra[2] + greenWeight * p1.bgra[1] + redWeight * p1.bgra[0];
      double res1 = pow(v1, gamma);
      if (res1 > MAX_BGR_VALUE) res1 = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res1, res1, res1);
    }
  }
  return output_image_ptr;
}

ImagePtr applyGammaUnrolled4x1(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  const int height = in_rows.size();
  const int width = in_rows[1] - in_rows[0];
  int limit = width -3;
  // Move constant luminance coefficients outside the loop
  // Code motion
  const double redWeight = 0.11;
  const double greenWeight = 0.59;
  const double blueWeight = 0.3;
  int i, j;
  for (i = 0; i < height; ++i) {
    for (j = 0; j < limit; j += 4) { // 2x1 loop unrolling
      Pixel p1 = in_rows[i][j];
      Pixel p2 = in_rows[i][j + 1];
      Pixel p3 = in_rows[i][j + 2];
      Pixel p4 = in_rows[i][j + 3];

      double v1 = blueWeight * p1.bgra[2] + greenWeight * p1.bgra[1] + redWeight * p1.bgra[0];
      double v2 = blueWeight * p2.bgra[2] + greenWeight * p2.bgra[1] + redWeight * p2.bgra[0];
      double v3 = blueWeight * p3.bgra[2] + greenWeight * p3.bgra[1] + redWeight * p3.bgra[0];
      double v4 = blueWeight * p4.bgra[2] + greenWeight * p4.bgra[1] + redWeight * p4.bgra[0];

      double res1 = pow(v1, gamma);
      double res2 = pow(v2, gamma);
      double res3 = pow(v3, gamma);
      double res4 = pow(v4, gamma);

      if (res1 > MAX_BGR_VALUE) res1 = MAX_BGR_VALUE;
      if (res2 > MAX_BGR_VALUE) res2 = MAX_BGR_VALUE;
      if (res3 > MAX_BGR_VALUE) res3 = MAX_BGR_VALUE;
      if (res4 > MAX_BGR_VALUE) res4 = MAX_BGR_VALUE;

      out_rows[i][j] = Pixel(res1, res1, res1);
      out_rows[i][j + 1] = Pixel(res2, res2, res2);
      out_rows[i][j + 2] = Pixel(res3, res3, res3);
      out_rows[i][j + 3] = Pixel(res4, res4, res4);
 
    }
    for(; j<width; j++){
      Pixel p1 = in_rows[i][j];
      double v1 = blueWeight * p1.bgra[2] + greenWeight * p1.bgra[1] + redWeight * p1.bgra[0];
      double res1 = pow(v1, gamma);
      if (res1 > MAX_BGR_VALUE) res1 = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res1, res1, res1);
    }
  }
  return output_image_ptr;
}

ImagePtr applyBlurUnrolled2x1(ImagePtr image_ptr) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_blurred", IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  const int height = in_rows.size() - 1; // Avoid edges
  const int width = in_rows[1] - in_rows[0] - 1; // Avoid edges
  double b[2], g[2], r[2];

  for (size_t i = 1; i < height; ++i) {
    for (int j = 1; j < width; j += 2) {
      b[0] = b[1] = 0;
      g[0] = g[1] = 0;
      r[0] = r[1] = 0;

      for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
          Pixel& p1 = in_rows[i + x][j + y];
          Pixel& p2 = in_rows[i + x][j + y + 1];
          b[0] += p1.bgra[0];
          g[0] += p1.bgra[1];
          r[0] += p1.bgra[2];
          b[1] += p2.bgra[0];
          g[1] += p2.bgra[1];
          r[1] += p2.bgra[2];
        }
      }

      Pixel& out_pixel1 = out_rows[i][j];
      Pixel& out_pixel2 = out_rows[i][j + 1];

      out_pixel1.bgra[0] = (b[0] / 9);
      out_pixel1.bgra[1] = (g[0] / 9);
      out_pixel1.bgra[2] = (r[0] / 9);
      out_pixel2.bgra[0] = (b[1] / 9);
      out_pixel2.bgra[1] = (g[1] / 9);
      out_pixel2.bgra[2] = (r[1] / 9);
    }
  }

  return output_image_ptr;
}

ImagePtr applyBlurUnrolled2x1withInnerLoop(ImagePtr image_ptr) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_blurred", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  const int height = in_rows.size()-1; //avoid edges
  const int width = in_rows[1] - in_rows[0]-1; //avoid edges
  double b,g,r;

  for (size_t i = 1; i < height; ++i) { 
      for (int j = 1; j < width; j += 2) {  
          for (int k = 0; k < 2; ++k) {
            b = 0; 
            g = 0; 
            r = 0;
            for (int x = -1; x <= 1; ++x) {
                for (int y = -1; y <= 1; ++y) {
                  Pixel& p = in_rows[i + x][j + y + k];
                  b += p.bgra[0];
                  g += p.bgra[1];
                  r += p.bgra[2];
                }
            }
            Pixel& out_pixel = out_rows[i][j + k];
            out_pixel.bgra[0] = (b / 9);
            out_pixel.bgra[1] = (g / 9);
            out_pixel.bgra[2] = (r / 9);
          }
      }
  }
    return output_image_ptr;
}

ImagePtr applyTintUnrolled2x1(ImagePtr image_ptr, const double *tints) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_tinted", IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  double b1,b2,g1,g2,r1,r2;

  for (unsigned long i = 0; i < in_rows.size(); ++i) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; j += 2) {
      b1 = (double)in_rows[i][j].bgra[0];
      g1 = (double)in_rows[i][j].bgra[1];
      r1 = (double)in_rows[i][j].bgra[2];
      b2 = (double)in_rows[i][j + 1].bgra[0];
      g2 = (double)in_rows[i][j + 1].bgra[1];
      r2 = (double)in_rows[i][j + 1].bgra[2];

      b1 = b1 + (MAX_BGR_VALUE - b1) * tints[0];
      g1 = g1 + (MAX_BGR_VALUE - g1) * tints[1];
      r1 = r1 + (MAX_BGR_VALUE - r1) * tints[2];
      b2 = b2 + (MAX_BGR_VALUE - b2) * tints[0];
      g2 = g2 + (MAX_BGR_VALUE - g2) * tints[1];
      r2 = r2 + (MAX_BGR_VALUE - r2) * tints[2];

      out_rows[i][j].bgra[0] = (b1 > MAX_BGR_VALUE) ? MAX_BGR_VALUE : b1;
      out_rows[i][j].bgra[1] = (g1 > MAX_BGR_VALUE) ? MAX_BGR_VALUE : g1;
      out_rows[i][j].bgra[2] = (r1 > MAX_BGR_VALUE) ? MAX_BGR_VALUE : r1;
      out_rows[i][j + 1].bgra[0] = (b2 > MAX_BGR_VALUE) ? MAX_BGR_VALUE : b2;
      out_rows[i][j + 1].bgra[1] = (g2 > MAX_BGR_VALUE) ? MAX_BGR_VALUE : g2;
      out_rows[i][j + 1].bgra[2] = (r2 > MAX_BGR_VALUE) ? MAX_BGR_VALUE : r2;
    }
  }

  return output_image_ptr;
}



void writeImage(ImagePtr image_ptr) {
  image_ptr->write( (image_ptr->name() + ".bmp").c_str());
}

void checkCorrectness(std::vector<ImagePtr> image_vector, std::vector<ImagePtr> image_vector_original){
  for(unsigned long k=0; k<image_vector.size(); k++){
    ImagePtr img1 = image_vector[k];
    ImagePtr img2 = image_vector_original[k];
    auto in_rows1 = img1->rows();
    auto in_rows2 = img2->rows();
    int height1 = img1->height();
    int width1 = img1->width();
    int height2 = img2->height();
    int width2 = img2->width();
    if(height1 != height2 || width1 != width2){
      printf("The two images do not have the same dimensions");
      return;
    }
    for(int i=0; i<height1; i++){
      for(int j=0; j<width1; j++){
        if(in_rows1[i][j].value != in_rows2[i][j].value){
            printf("Correctness check failed for pixels image_%ld[%d][%d] = %d and image_refrence_%ld[%d][%d] = %d\n", (k+1), i, j, in_rows1[i][j].value, (k+1), i, i, in_rows2[i][j].value);
            return;
        }
      }
    }
    printf("Correctness check passed for image %ld\n", (k+1));
  }
}

int main() {
  const double tint_array[] = {0.75, 0, 0};
  std::vector<ImagePtr> image_vector, image_vector_reference;
  for ( int i = 2000; i <= 2000000; i *= 10 ) {
    image_vector.push_back(makeFractalImage(i));
    image_vector_reference.push_back(makeFractalImage(i));
  }
  // Store the output of the original functions applyGamma, applyTint, applyBlur in image_vector_reference
  // The output images will be used for checking the correctness of the optimized functions

  //Adding clock variables to this for loop to determine the time each function takes
  clock_t start, end;
  int counter = 1;
  int total;
  int gamma;
  int tint;
  int blur;
  float g;
  float t;
  float b;

  for(ImagePtr img:image_vector_reference){
    total = 0;
    gamma = 0;
    blur = 0;
    tint = 0;
    start = clock();
    img = applyGammaUnrolled4x1(img, 1.4);
    end = clock();
    gamma = (end-start);
    printf("Image %d Gamma 4x1 Function Time: %d cycles.\n", counter, gamma);

    total += gamma;

    start = clock();
    img = applyTintUnrolled2x1(img, tint_array);
    end = clock();
    tint = (end-start);
    printf("Image %d Tint 2x1 Function Time: %d cycles.\n", counter, tint);
    total += tint;

    start = clock();
    img = applyBlurUnrolled2x1withInnerLoop(img);
    end = clock();
    blur = (end-start);
    printf("Image %d Blur 2x1 IL Function Time: %d cycles =", counter, blur);
    total += blur;

    printf(" Total Time: %d cycles.\n", total);
    g = (float) gamma / total;  
    g *= 100;
    t = (float) tint / total;
    t *= 100;
    b = (float) blur / total;
    b *= 100;
    printf("Gamma Percentage: %.2f%%\n", g);
    printf("Tint Percentage: %.2f%%\n", t);
    printf("Blur Percentage: %.2f%%\n\n", b);
    counter++;
  }

  process_images(image_vector);
  checkCorrectness(image_vector, image_vector_reference);
  return 0;
}

// INITIAL RUN WITHOUT OMTIMIZATIONS:
// Image 1 Gamma Function Time: 28741 cycles.
// Image 1 Tint Function Time: 15100 cycles.
// Image 1 blur Function Time: 29778 cycles = Total Time: 73619 cycles.
// Gamma Percentage: 39.04%
// Tint Percentage: 20.51%
// Blur Percentage: 40.45%
// Image 2 Gamma Function Time: 30424 cycles.
// Image 2 Tint Function Time: 15800 cycles.
// Image 2 blur Function Time: 28924 cycles = Total Time: 75148 cycles.
// Gamma Percentage: 40.49%
// Tint Percentage: 21.03%
// Blur Percentage: 38.49%
// Image 3 Gamma Function Time: 32480 cycles.
// Image 3 Tint Function Time: 15915 cycles.
// Image 3 blur Function Time: 28913 cycles = Total Time: 77308 cycles.
// Gamma Percentage: 42.01%
// Tint Percentage: 20.59%
// Blur Percentage: 37.40%
// Image 4 Gamma Function Time: 30672 cycles.
// Image 4 Tint Function Time: 15700 cycles.
// Image 4 blur Function Time: 28914 cycles = Total Time: 75286 cycles.
// Gamma Percentage: 40.74%
// Tint Percentage: 20.85%
// Blur Percentage: 38.41%
// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4


// REDUCED FUNCTION CALLS ON GAMMA FUNCTION
// Image 1 Gamma RFC Function Time: 26358 cycles.
// Image 1 Tint Function Time: 15872 cycles.
// Image 1 blur Function Time: 29336 cycles = Total Time: 71566 cycles.
// Gamma Percentage: 36.83%
// Tint Percentage: 22.18%
// Blur Percentage: 40.99%
// Image 2 Gamma RFC Function Time: 28610 cycles.
// Image 2 Tint Function Time: 15958 cycles.
// Image 2 blur Function Time: 29141 cycles = Total Time: 73709 cycles.
// Gamma Percentage: 38.81%
// Tint Percentage: 21.65%
// Blur Percentage: 39.54%
// Image 3 Gamma RFC Function Time: 30915 cycles.
// Image 3 Tint Function Time: 15345 cycles.
// Image 3 blur Function Time: 29161 cycles = Total Time: 75421 cycles.
// Gamma Percentage: 40.99%
// Tint Percentage: 20.35%
// Blur Percentage: 38.66%
// Image 4 Gamma RFC Function Time: 28610 cycles.
// Image 4 Tint Function Time: 15873 cycles.
// Image 4 blur Function Time: 28932 cycles = Total Time: 73415 cycles.
// Gamma Percentage: 38.97%
// Tint Percentage: 21.62%
// Blur Percentage: 39.41%
// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4


// 2x1 LOOP UNROLLING ON GAMMA FUNCTION
// Image 1 Gamma 2x1 Function Time: 23183 cycles.
// Image 1 Tint Function Time: 15644 cycles.
// Image 1 blur Function Time: 29622 cycles = Total Time: 68449 cycles.
// Gamma Percentage: 33.87%
// Tint Percentage: 22.85%
// Blur Percentage: 43.28%
// Image 2 Gamma 2x1 Function Time: 23864 cycles.
// Image 2 Tint Function Time: 16037 cycles.
// Image 2 blur Function Time: 29951 cycles = Total Time: 69852 cycles.
// Gamma Percentage: 34.16%
// Tint Percentage: 22.96%
// Blur Percentage: 42.88%
// Image 3 Gamma 2x1 Function Time: 24984 cycles.
// Image 3 Tint Function Time: 15632 cycles.
// Image 3 blur Function Time: 29261 cycles = Total Time: 69877 cycles.
// Gamma Percentage: 35.75%
// Tint Percentage: 22.37%
// Blur Percentage: 41.88%
// Image 4 Gamma 2x1 Function Time: 25091 cycles.
// Image 4 Tant Fun Fon tine: 1546 cycles cles.
// Image 4 blur Function Time: 29254 cycles = Total Time: 68959 cycles.
// Gamma Percentage: 35.15%
// Tint Percentage: 22.42%
// Blur Percentage: 42.42%
// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4


// 4x1 LOOP UNROLLING ON GAMMA FUNCTION
// Image 1 Gamma 4x1 Function Time: 23744 cycles.
// Image 1 Tint Function Time: 16337 cycles.
// Image 1 blur Function Time: 31144 cycles = Total Time: 71225 cycles.
// Gamma Percentage: 33.34%
// Tint Percentage: 22.94%
// Blur Percentage: 43.73%
// Image 2 Gamma 4x1 Function Time: 23662 cycles.
// Image 2 Tint Function Time: 16745 cycles.
// Image 2 blur Function Time: 29985 cycles = Total Time: 70392 cycles.
// Gamma Percentage: 33.61%
// Tint Percentage: 23.79%
// Blur Percentage: 42.60%
// Image 3 Gamma 4x1 Function Time: 24041 cycles.
// Image 3 Tint Function Time: 15525 cycles.
// Image 3 blur Function Time: 29973 cycles = Total Time: 69539 cycles.
// Gamma Percentage: 34.57%
// Tint Percentage: 22.33%
// Blur Percentage: 43.10%
// Image 4 Gamma 4x1 Function Time: 24141 cycles.
// Image 4 Tint Function Time: 16650 cycles.
// Image 4 blur Function Time: 29814 cycles = Total Time: 70605 cycles.
// Gamma Percentage: 34.19%
// Tint Percentage: 23.58%
// Blur Percentage: 42.23%
// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4


// 2x1 LOOP UNROLLING ON BLUR FUNCTION
// Image 1 Gamma 4x1 Function Time: 21881 cycles.
// Image 1 Tint Function Time: 15213 cycles.
// Image 1 Blur 2x1 Function Time: 28043 cycles = Total Time: 65137 cycles.
// Gamma Percentage: 33.59%
// Tint Percentage: 23.36%
// Blur Percentage: 43.05%

// Image 2 Gamma 4x1 Function Time: 23396 cycles.
// Image 2 Tint Function Time: 15993 cycles.
// Image 2 Blur 2x1 Function Time: 29993 cycles = Total Time: 69382 cycles.
// Gamma Percentage: 33.72%
// Tint Percentage: 23.05%
// Blur Percentage: 43.23%

// Image 3 Gamma 4x1 Function Time: 23218 cycles.
// Image 3 Tint Function Time: 15726 cycles.
// Image 3 Blur 2x1 Function Time: 30091 cycles = Total Time: 69035 cycles.
// Gamma Percentage: 33.63%
// Tint Percentage: 22.78%
// Blur Percentage: 43.59%

// Image 4 Gamma 4x1 Function Time: 23345 cycles.
// Image 4 Tint Function Time: 15198 cycles.
// Image 4 Blur 2x1 Function Time: 28069 cycles = Total Time: 66612 cycles.
// Gamma Percentage: 35.05%
// Tint Percentage: 22.82%
// Blur Percentage: 42.14%

// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4


// 2x1 LOOP UNROLLING WITH INNER LOOP ON BLUR FUNCTION 
// Image 1 Gamma 4x1 Function Time: 23122 cycles.
// Image 1 Tint Function Time: 14962 cycles.
// Image 1 Blue 2x1 IL Function Time: 27390 cycles = Total Time: 65474 cycles.
// Gamma Percentage: 35.31%
// Tint Percentage: 22.85%
// Blur Percentage: 41.83%
// Image 2 Gamma 4x1 Function Time: 21713 cycles.
// Image 2 Tint Function Time: 14605 cycles.
// Image 2 Blue 2x1 IL Function Time: 27222 cycles = Total Time: 63540 cycles.
// Gamma Percentage: 34.17%
// Tint Percentage: 22.99%
// Blur Percentage: 42.84%
// Image 3 Gamma 4x1 Function Time: 22223 cycles.
// Image 3 Tint Function Time: 14621 cycles.
// Image 3 Blue 2x1 IL Function Time: 27205 cycles = Total Time: 64049 cycles.
// Gamma Percentage: 34.70%
// Tint Percentage: 22.83%
// Blur Percentage: 42.48%
// Image 4 Gamma 4x1 Function Time: 22659 cycles.
// Image 4 Tint Function Time: 14784 cycles.
// Image 4 Blue 2x1 IL Function Time: 27224 cycles = Total Time: 64667 cycles.
// Gamma Percentage: 35.04%
// Tint Percentage: 22.86%
// Blur Percentage: 42.10%
// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4


// 2x1 LOOP UNROLLING FOR TINT FUNCTION
// Image 1 Gamma 4x1 Function Time: 22081 cycles.
// Image 1 Tint 2x1 Function Time: 9629 cycles.
// Image 1 Blur 2x1 IL Function Time: 26649 cycles = Total Time: 58359 cycles.
// Gamma Percentage: 37.84%
// Tint Percentage: 16.50%
// Blur Percentage: 45.66%
// Image 2 Gamma 4x1 Function Time: 21959 cycles.
// Image 2 Tint 2x1 Function Time: 9101 cycles.
// Image 2 Blur 2x1 IL Function Time: 26470 cycles = Total Time: 57530 cycles.
// Gamma Percentage: 38.17%
// Tint Percentage: 15.82%
// Blur Percentage: 46.01%
// Image 3 Gamma 4x1 Function Time: 22230 cycles.
// Image 3 Tint 2x1 Function Time: 9098 cycles.
// Image 3 Blur 2x1 IL Function Time: 26528 cycles = Total Time: 57856 cycles.
// Gamma Percentage: 38.42%
// Tint Percentage: 15.73%
// Blur Percentage: 45.85%
// Image 4 Gamma 4x1 Function Time: 22158 cycles.
// Image 4 Tint 2x1 Function Time: 9051 cycles.
// Image 4 Blur 2x1 IL Function Time: 26639 cycles = Total Time: 57848 cycles.
// Gamma Percentage: 38.30%
// Tint Percentage: 15.65%
// Blur Percentage: 46.05%
// Correctness check passed for image 1
// Correctness check passed for image 2
// Correctness check passed for image 3
// Correctness check passed for image 4