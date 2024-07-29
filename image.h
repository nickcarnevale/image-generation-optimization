/*
Copyright (C) 2019 Intel Corporation

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.

SPDX-License-Identifier: MIT
*/
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 800;
const int MAX_BGR_VALUE = 255;

union Pixel {
    std::uint8_t bgra[4];
    std::uint32_t value;
    Pixel() {}
    template <typename T> Pixel(T b, T g, T r) {
      bgra[0] = (std::uint8_t)b, 
      bgra[1] = (std::uint8_t)g, 
      bgra[2] = (std::uint8_t)r, 
      bgra[3] = 0;
    }
};

class Image {
  private:
    std::string myName;
    int myWidth;
    int myHeight;
    int myPadSize;
    std::vector<Pixel> myData; //raw raster data
    std::vector<Pixel*> myRows;

    void reset(int w, int h);
    //don't allow copying
    Image(const Image&);
    void operator=(const Image&);
    //data structures 'file' and 'info' are using to store an Image as BMP file
    //for more details see https://en.wikipedia.org/wiki/BMP_file_format
    using BITMAPFILEHEADER = struct {
      std::uint16_t sizeRest; // field is not from specification, 
              // was added for alignemt. store size of rest of the fields
      std::uint16_t type;
      std::uint32_t size;
      std::uint32_t reserved; 
      std::uint32_t offBits;
    };
    BITMAPFILEHEADER file;

    using BITMAPINFOHEADER = struct {
      std::uint32_t size;
      std::int32_t width;
      std::int32_t height;
      std::uint16_t planes;
      std::uint16_t bitCount;
      std::uint32_t compression;
      std::uint32_t sizeImage;
      std::int32_t xPelsPerMeter;
      std::int32_t yPelsPerMeter;
      std::uint32_t clrUsed;
      std::uint32_t clrImportant;
    };
    BITMAPINFOHEADER info;
 public:
  // constructor
  Image(const std::string &n, int w, int h);
  std::string name() const;
  std::string setName(const std::string &n);
  int width() const;
  int height() const;
  void write(const char* fname) const;
  void fill(std::uint8_t r, 
            std::uint8_t g, 
            std::uint8_t b, 
            int x, 
            int y);
  template <typename F> void fill(F f);
  std::vector<Pixel*>& rows();

};

//! Fractal class
class Fractal {
    private:
        //! Size of the Fractal area
        const int mySize[2];
        //! Fractal properties
        double cx = -0.7436;
        const double cy = 0.1319;
        const double myMagn;
        const int maxIter = 1000;
    public:
        //! Constructor
        Fractal(int x, int y, double m);
        //! One pixel calculation routine
        double calcOnePixel(int x0, int y0); 
        
};
std::shared_ptr<Image> makeFractalImage(double magn);



