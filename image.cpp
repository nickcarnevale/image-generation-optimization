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

#include "image.h"


void Image::reset(int w, int h) {
      if(w <= 0 || h <= 0) {
        std::cout << "Warning: Invalid Image size.\n";
        return;
      }
      myWidth = w, myHeight = h;
      //reset raw data
      myData.resize(myWidth*myHeight);
      myRows.resize(myHeight);
      //reset rows
      for(unsigned long i = 0; i < myRows.size(); ++i)
        myRows[i] = &myData[0]+i*myWidth;

      myPadSize = (4-(w*sizeof(myData[0]))%4)%4;
      int sizeData = w*h*sizeof(myData[0]) + h*myPadSize;
      int sizeAll = sizeData + sizeof(file) + sizeof(info);

      //BITMAPFILEHEADER
      file.sizeRest = 14;
      file.type = 0x4d42; //same as 'BM' in ASCII
      file.size = sizeAll;
      file.reserved = 0;
      file.offBits = 54;

      //BITMAPINFOHEADER
      info.size = 40;
      info.width = w;
      info.height = h;
      info.planes = 1;
      info.bitCount = 32;
      info.compression = 0;
      info.sizeImage = sizeData;
      info.yPelsPerMeter = 0;
      info.xPelsPerMeter = 0;
      info.clrUsed = 0;
      info.clrImportant = 0;
  }
  // constructor
  Image::Image(const std::string &n, int w = 1920, int h = 1080) : myName(n) {
    reset(w, h);
  }

  std::string Image::name() const { return myName; }
  std::string Image::setName(const std::string &n) { return myName = n; }

  int Image::width() const { return myWidth; }
  int Image::height() const { return myHeight; }

  void Image::write(const char* fname) const {
    if(myData.empty()) {
      std::cout << "Warning: Image is empty.\n";
      return;
    }
    std::ofstream stream{fname};
    stream.write((char*)&file.type, file.sizeRest);
    stream.write((char*)&info, info.size);
    stream.write((char*)myData[0].bgra, myData.size()*sizeof(myData[0]));
  }

  void Image::fill(std::uint8_t r, std::uint8_t g, 
                   std::uint8_t b, int x = -1, int y = -1) {
    if(myData.empty())
      return;

    if(x < 0 && y < 0) //fill whole Image
      std::fill(myData.begin(), myData.end(), Pixel(b, g, r));
    else {
      auto& bgra = myData[myWidth*x + y].bgra;
      bgra[3] = 0, bgra[2] = r, bgra[1] = g, bgra[0] = b;
    }
  }

  template <typename F> void Image::fill(F f) {
    if(myData.empty())
      reset(myWidth, myHeight);

    int i = -1;
    int w = this->myWidth;
    std::for_each(myData.begin(), myData.end(), [&i, w, f](Pixel& p) {
      ++i;
      int x = i / w, y = i % w;
      auto val = f(x, y);
      if(val > 255)
        val = 255;
      p = Pixel(val, val, val);
    });
  }

  std::vector<Pixel*>& Image::rows() { return myRows; }

  //! Constructor
Fractal::Fractal(int x, int y, double m = 2000000.0): mySize{x, y}, myMagn(m) {}

//! One pixel calculation routine
double Fractal::calcOnePixel(int x0, int y0) {
    double fx0 = double(x0) - double(mySize[0]) / 2;
    double fy0 = double(y0) - double(mySize[1]) / 2;
    fx0 = fx0 / myMagn + cx;
    fy0 = fy0 / myMagn + cy;

    double res = 0, x = 0, y = 0;
    for(int iter = 0; x*x + y*y <= 4 && iter < maxIter; ++iter) {
        const double val = x*x - y*y + fx0;
        y = 2*x*y + fy0, x = val;
        res += exp(-sqrt(x*x+y*y));
    }
    return res;
}
std::shared_ptr<Image> makeFractalImage(double magn = 2000000) {
  static int index = 1;
  const std::string name = std::string("fractal_") + std::to_string(index);
  auto image_ptr = std::make_shared<Image>(name, IMAGE_WIDTH, IMAGE_HEIGHT);
  Fractal fr(image_ptr->width(), image_ptr->height(), magn);
  image_ptr->fill([&fr](int x, int y) { return fr.calcOnePixel(x, y); });
  index++;
  return image_ptr;
}

