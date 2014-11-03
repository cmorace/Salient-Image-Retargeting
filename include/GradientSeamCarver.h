//
//  GradientSeamCarver.h
//  ImageRetargeting
//
//  Created by Charles Morace on 11/2/14.
//
//

#ifndef ImageRetargeting_GradientSeamCarver_h
#define ImageRetargeting_GradientSeamCarver_h

#include <stdio.h>
#include "segment-image.h"

#include "opencv2/opencv.hpp"
#include "CinderOpenCV.h"

using namespace cv;


struct ImageSeam
{
    unsigned int index = 0;
    unsigned int preValue = 0;
    unsigned int totalValue = 0;
    bool isVerticalSeam = true;
    bool isForwardSeam = true;
    int* seamPath;
};

bool compareImageSeam(ImageSeam a, ImageSeam b)
{
    return a.totalValue < b.totalValue;
}


class GradientSeamCarver {
public:
    GradientSeamCarver();
    cinder::Surface getGradientImage(cinder::Surface imgData, bool isSobel);
    //cinder::Surface computeVerticalSeam(cinder::Surface imgData);
    cinder::Surface getVerticalSeams(cinder::Surface imgData, int nSeams = 1, bool returnGradient=false);
    cinder::Surface getHorizontalSeams(cinder::Surface imgData, int nSeams = 1, bool returnGradient = false);
    cinder::Surface deleteVerticalSeam(cinder::Surface imgData);
    int scale;
    int delta;
    
    double gradTime;
    double carveTime;
    unsigned int nbOfSeams;
    
private:
    cinder::Timer* timer;
    ImageSeam* getTopToBottom(cinder::Surface imgData);
    ImageSeam* getBottomToTop(cinder::Surface imgData);
    ImageSeam* getRightToLeft(cinder::Surface imgData);
    ImageSeam* getLeftToRight(cinder::Surface imgData);
    void drawSeams(cinder::Surface imgData, int nSeams, std::vector<ImageSeam> seamVector);
    bool isSobel = true;
};

GradientSeamCarver::GradientSeamCarver()
{
    timer = new cinder::Timer(false);
   
    scale = 1;
    delta = 50;
    
    gradTime = 0;
    gradTime = 0;
    nbOfSeams = 10;
    
}



cinder::Surface GradientSeamCarver::getGradientImage(cinder::Surface imgData, bool isSobel)
{
    this->isSobel = isSobel;
    Mat src, src_gray;
    
    Mat grad;
    
    int ddepth = CV_16S;
    
    src = cinder::toOcv(imgData);
    timer->start();
    GaussianBlur( src, src, cv::Size(3,3), 0, 0, BORDER_DEFAULT );
    
    /// Convert it to gray
    cvtColor( src, src_gray, CV_RGB2GRAY );
    
    /// Generate grad_x and grad_y
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;
    
    /// Gradient X
    if(isSobel){
      Sobel( src_gray, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
    }
    else{
       Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
    }
    //
    convertScaleAbs( grad_x, abs_grad_x );
    
    /// Gradient Y
    if(isSobel){
       Sobel( src_gray, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
    }
    else{
       Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
    }
    convertScaleAbs( grad_y, abs_grad_y );
    
    /// Total Gradient (approximate)
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
    timer->stop();
    gradTime = timer->getSeconds();
    return cinder::Surface(cinder::fromOcv(grad));
}



void GradientSeamCarver::drawSeams(cinder::Surface imgData, int nSeams, std::vector<ImageSeam> seamVector)
{
    ci::ColorT<uchar>* r = new ci::ColorT<uchar>(255,0,0);
    for (std::vector<ImageSeam>::iterator it=seamVector.begin(); it!=seamVector.begin() + nSeams; ++it)
    {
        
        int h = imgData.getHeight();
        int w = imgData.getWidth();
        if(it->isVerticalSeam){
            int x = it->index;
            if(it->isForwardSeam){
                for(int y=h-1; y>=0; y--){
                    imgData.setPixel(ci::Vec2i(x,y), *r);
                    x += it->seamPath[y*w+x];
                }
            }
            else{
                for(int y=0; y<h; y++){
                    imgData.setPixel(ci::Vec2i(x,y), *r);
                    x += it->seamPath[y*w+x];
                    
                }
            }
        }
        else{
            int y = it->index;
            if(it->isForwardSeam){
                for(int x=w-1; x>=0; x--){
                    imgData.setPixel(ci::Vec2i(x,y), *r);
                    y += it->seamPath[y*w+x];
                }
            }
            else{
                for(int x=0; x<w; x++){
                    imgData.setPixel(ci::Vec2i(x,y), *r);
                    y += it->seamPath[y*w+x];
                    
                }
            }
        }
    }
}


cinder::Surface GradientSeamCarver::getVerticalSeams(cinder::Surface img, int nSeams,bool returnGrad)
{
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    //
    timer->start();
    int w = imgData.getWidth();
    ImageSeam* verticalSeam_bt = getBottomToTop(imgData);
    ImageSeam* verticalSeam_tb = getTopToBottom(imgData);
    ImageSeam* verticalSeams = new ImageSeam[2*w];
    for(int i=0; i<w; i++){
        verticalSeams[2*i] = verticalSeam_bt[i];
        verticalSeams[2*i+1] = verticalSeam_tb[i];
    }
    
    std::vector<ImageSeam> c(verticalSeams, verticalSeams + 2*w);
    std::sort (c.begin(), c.end(), compareImageSeam);
    timer->stop();
    carveTime  = timer->getSeconds();
    //
    if(returnGrad){
        drawSeams(imgData, nSeams, c);
        return imgData;
    }
    else{
        drawSeams(img, nSeams, c);
        return img;
    }

    
    
}

cinder::Surface GradientSeamCarver::getHorizontalSeams(cinder::Surface img, int nSeams, bool returnGrad)
{
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    //
    timer->start();
    int h = imgData.getHeight();
    ImageSeam* horizSeam_lr = getLeftToRight(imgData);
    ImageSeam* horizSeam_rl = getRightToLeft(imgData);
    ImageSeam* horizSeams = new ImageSeam[2*h];
    
    for(int i=0; i<h; i++){
        horizSeams[2*i] = horizSeam_lr[i];
        horizSeams[2*i+1] = horizSeam_rl[i];
    }
    
    std::vector<ImageSeam> c(horizSeams, horizSeams + 2*h);
    std::sort (c.begin(), c.end(), compareImageSeam);
    
    timer->stop();
    carveTime  = timer->getSeconds();
    
    if(returnGrad){
        drawSeams(imgData, nSeams, c);
        return imgData;
    }
    else{
        drawSeams(img, nSeams, c);
        return img;
    }
}


cinder::Surface GradientSeamCarver::deleteVerticalSeam(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    int seamPath[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    
    for(int x=0; x<w; x++){
        seamPath[x] = 0;
        verticalSeam[x] = ImageSeam();
        verticalSeam[x].index = x;
        verticalSeam[x].preValue = (imgData.getPixel(ci::Vec2i(x,0))).b;
        verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,0))).b;
    }
    
    for(int y=1; y<h; y++){
        for(int x=0; x<w; x++){
            
            if(x==0){//First Column
                if(verticalSeam[x].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue <= verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int x=0; x<w; x++){
            verticalSeam[x].preValue = verticalSeam[x].totalValue;
        }
    }
    
    std::vector<ImageSeam> myvector (verticalSeam, verticalSeam+w);
    std::sort (myvector.begin(), myvector.end(), compareImageSeam);
    
    ci::Surface s = cinder::Surface(w-1,h,false);
    
    for (std::vector<ImageSeam>::iterator it=myvector.begin(); it!=myvector.begin() + 1; ++it)
    {
        int seamX = it->index;
        for(int y=h-1; y>=0; y--){
            //s.setPixel(ci::Vec2i(x,y), *r);
            
            for (int x = 0; x<w-1; x++) {
                if(x>=seamX){
                    s.setPixel(ci::Vec2i(x,y), imgData.getPixel(ci::Vec2i(x+1,y)));
                }
                else{
                  s.setPixel(ci::Vec2i(x,y), imgData.getPixel(ci::Vec2i(x,y)));  
                }
            }
            seamX += seamPath[y*w+seamX];
        }
    }
    return s;
}


ImageSeam* GradientSeamCarver::getTopToBottom(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    
    for(int x=0; x<w; x++){
        seamPath[x] = 0;
        verticalSeam[x] = ImageSeam();
        verticalSeam[x].index = x;
        verticalSeam[x].preValue = (imgData.getPixel(ci::Vec2i(x,0))).b;
        verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,0))).b;
        verticalSeam[x].isVerticalSeam = true;
        verticalSeam[x].isForwardSeam = true;
        verticalSeam[x].seamPath = seamPath;
    }
    
    for(int y=1; y<h; y++){
        for(int x=0; x<w; x++){
            
            if(x==0){//First Column
                if(verticalSeam[x].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue <= verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int x=0; x<w; x++){
            verticalSeam[x].preValue = verticalSeam[x].totalValue;
        }
    }
    return verticalSeam;
}

ImageSeam* GradientSeamCarver::getBottomToTop(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    
    for(int x=0; x<w; x++){
        seamPath[(w-1)*h + x] = 0;
        verticalSeam[x] = ImageSeam();
        verticalSeam[x].index = x;
        verticalSeam[x].preValue = (imgData.getPixel(ci::Vec2i(x,h-1))).b;
        verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,h-1))).b;
        verticalSeam[x].isVerticalSeam = true;
        verticalSeam[x].isForwardSeam = false;
        verticalSeam[x].seamPath = seamPath;
    }
    
    for(int y=h-2; y>=0; y--){
        for(int x=0; x<w; x++){
            
            if(x==0){//First Column
                if(verticalSeam[x].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue <= verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                    
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int x=0; x<w; x++){
            verticalSeam[x].preValue = verticalSeam[x].totalValue;
        }
    }
    return verticalSeam;
}

ImageSeam* GradientSeamCarver::getLeftToRight(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[h];
    
    for(int y=0; y<h; y++){
        seamPath[y] = 0;
        verticalSeam[y] = ImageSeam();
        verticalSeam[y].index = y;
        verticalSeam[y].preValue = (imgData.getPixel(ci::Vec2i(0,y))).b;
        verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(0,y))).b;
        verticalSeam[y].isVerticalSeam = false;
        verticalSeam[y].isForwardSeam = true;
        verticalSeam[y].seamPath = seamPath;
    }
    
    for(int x=1; x<w; x++){
        for(int y=0; y<h; y++){
            
            if(y==0){//First Row
                if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(verticalSeam[y].preValue <= verticalSeam[y-1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[y-1].preValue <= verticalSeam[y].preValue && verticalSeam[y-1].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int y=0; y<h; y++){
            verticalSeam[y].preValue = verticalSeam[y].totalValue;
        }
    }
    return verticalSeam;
}

ImageSeam* GradientSeamCarver::getRightToLeft(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[h];
    
    //fill right row with 0
    for(int y=0; y<h; y++){
        seamPath[w*y+w-1] = 0;
        verticalSeam[y] = ImageSeam();
        verticalSeam[y].index = y;
        verticalSeam[y].preValue = (imgData.getPixel(ci::Vec2i(w-1,y))).b;
        verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(w-1,y))).b;
        verticalSeam[y].isVerticalSeam = false;
        verticalSeam[y].isForwardSeam = false;
        verticalSeam[y].seamPath = seamPath;
    }
    
    //start with row 1 from right
    for(int x=w-2; x>=0; x--){
        for(int y=0; y<h; y++){
            
            if(y==0){//First Row
                if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(verticalSeam[y].preValue <= verticalSeam[y-1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                    
                }
                else{
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Rows
                if(verticalSeam[y-1].preValue <= verticalSeam[y].preValue && verticalSeam[y-1].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int y=0; y<h; y++){
            verticalSeam[y].preValue = verticalSeam[y].totalValue;
        }
    }
    return verticalSeam;
}

#endif
