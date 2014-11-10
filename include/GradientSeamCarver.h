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
    Mat getGradientImageCV(cinder::Surface imgData, bool isSobel);
    std::vector<ImageSeam> getVerticalSeams(cinder::Surface img);
    std::vector<ImageSeam> getHorizontalSeams(cinder::Surface img);
    cinder::Surface drawVerticalSeams(cinder::Surface imgData, int nSeams = 1, bool returnGradient=false);
    cinder::Surface drawHorizontalSeams(cinder::Surface imgData, int nSeams = 1, bool returnGradient = false);
    cinder::Surface deleteVerticalSeam(cinder::Surface imgData);
    cinder::Surface deleteHorizontalSeam(cinder::Surface imgData);
    cinder::Surface deleteMinEnergySeam(cinder::Surface imgData);
    cinder::Surface addVerticalSeam(cinder::Surface imgData, int nSeams);
    cinder::Surface addHorizontalSeam(cinder::Surface imgData, int nSeams);
    
    unsigned int scale;
    unsigned int delta;
    
    double gradTime;
    double carveTime;
    unsigned int nbOfSeams;
    
    unsigned int newWidth;
    unsigned int newHeight;
    
private:
    cinder::Timer* timer;
    ImageSeam* getTopToBottom(cinder::Surface imgData);
    ImageSeam* getBottomToTop(cinder::Surface imgData);
    ImageSeam* getRightToLeft(cinder::Surface imgData);
    ImageSeam* getLeftToRight(cinder::Surface imgData);
    
    ImageSeam* getTopToBottom(Mat imgData);
    ImageSeam* getBottomToTop(Mat imgData);
    ImageSeam* getRightToLeft(Mat imgData);
    ImageSeam* getLeftToRight(Mat imgData);
    
    void drawSeams(cinder::Surface imgData, int nSeams, std::vector<ImageSeam> seamVector);
    cinder::Surface carveSeam(cinder::Surface imgData, std::vector<ImageSeam>::iterator seam);
    bool isSobel = true;
};

GradientSeamCarver::GradientSeamCarver()
{
    timer = new cinder::Timer(false);
   
    scale = 1;
    delta = 50;
    
    gradTime = 0.0;
    carveTime = 0.0;
    nbOfSeams = 1;
    
    newWidth =100;
    newHeight =100;
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

Mat GradientSeamCarver::getGradientImageCV(cinder::Surface imgData, bool isSobel)
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
    return grad;
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

cinder::Surface GradientSeamCarver::carveSeam(cinder::Surface imgData, std::vector<ImageSeam>::iterator seam)
{
    
    //ci::ColorT<uchar>* g = new ci::ColorT<uchar>(0,255,0);
    int h = imgData.getHeight();
    int w = imgData.getWidth();
    cinder::Surface s;
    
    if(seam->isVerticalSeam){
        s = cinder::Surface(w-1,h,false);
        int seamX = seam->index;
        if(seam->isForwardSeam){
            for(int y=h-1; y>=0; y--){
                for (int x = 0; x<w-1; x++) {
                    if(x<seamX){
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y)));
                    }
                    else{
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x+1,y)));
                    }
                }
                seamX += seam->seamPath[y*w+seamX];
            }
        }
        else{
            for(int y=0; y<h; y++){
                for (int x = 0; x<w-1; x++) {
                    if(x<seamX){
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y)));
                    }
                    else{
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x+1,y)));
                    }
                }
                seamX += seam->seamPath[y*w+seamX];
            }
        }
    }
    else{
        s = cinder::Surface(w,h-1,false);
        int seamY = seam->index;
        if(seam->isForwardSeam){
            for(int x=w-1; x>=0; x--){
                for (int y = 0; y<h-1; y++) {
                    if(y<seamY){
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y)));
                    }
                    else{
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y+1)));
                    }
                }
                seamY += seam->seamPath[seamY*w+x];
            }
        }
        else{
            for(int x=0; x<w; x++){
                for (int y = 0; y<h-1; y++) {
                    if(y<seamY){
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y)));
                    }
                    else{
                        s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y+1)));
                    }
                }
                seamY += seam->seamPath[seamY*w+x];
            }
        }
    }
    return s;
}

std::vector<ImageSeam> GradientSeamCarver::getVerticalSeams(cinder::Surface imgData)
{
   //test
    //Mat img = cinder::toOcv(imgData);
    int w = imgData.getWidth();
    ImageSeam* vertSeam_tb = getTopToBottom(imgData);
    ImageSeam* vertSeam_bt = getBottomToTop(imgData);
    ImageSeam* vertSeam = new ImageSeam[2*w];
    
    for(int i=0; i<w; i++){
        vertSeam[2*i] = vertSeam_tb[i];
        vertSeam[2*i+1] = vertSeam_bt[i];
    }
    
    std::vector<ImageSeam> c(vertSeam, vertSeam + 2*w);
    std::sort (c.begin(), c.end(), compareImageSeam);
    return c;
}

std::vector<ImageSeam> GradientSeamCarver::getHorizontalSeams(cinder::Surface imgData)
{
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
    
    return c;
}

cinder::Surface GradientSeamCarver::drawVerticalSeams(cinder::Surface img, int nSeams,bool returnGrad)
{
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    //
    timer->start();
    std::vector<ImageSeam> c = getVerticalSeams(imgData);
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

cinder::Surface GradientSeamCarver::drawHorizontalSeams(cinder::Surface img, int nSeams, bool returnGrad)
{
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    //
    timer->start();
    std::vector<ImageSeam> c = getHorizontalSeams(imgData);
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

cinder::Surface GradientSeamCarver::deleteVerticalSeam(cinder::Surface img)
{
   
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    timer->start();
    std::vector<ImageSeam> c = getVerticalSeams(imgData);
    carveTime  = timer->getSeconds();
    img = carveSeam(img, c.begin());
    return img;
}

cinder::Surface GradientSeamCarver::deleteHorizontalSeam(cinder::Surface img)
{
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    //
    timer->start();
    std::vector<ImageSeam> c = getHorizontalSeams(imgData);
    timer->stop();
    carveTime  = timer->getSeconds();
    img = carveSeam(img, c.begin());
    return img;
}

cinder::Surface GradientSeamCarver::deleteMinEnergySeam(cinder::Surface img){
    cinder::Surface imgData = getGradientImage(img.clone(),this->isSobel);
    //
    timer->start();
    std::vector<ImageSeam> a = getHorizontalSeams(imgData);
    std::vector<ImageSeam> b = getVerticalSeams(imgData);
    timer->stop();
    carveTime  = timer->getSeconds();
    std::vector<ImageSeam>::iterator aSeam =  a.begin();
    std::vector<ImageSeam>::iterator bSeam =  b.begin();
    if(1.f*aSeam->totalValue/img.getWidth() > 1.f*bSeam->totalValue/img.getHeight()){
      img = carveSeam(img, aSeam);
    }
    else{
      img = carveSeam(img, bSeam);
    }
    
    return img;
}

cinder::Surface GradientSeamCarver::addVerticalSeam(cinder::Surface imgData, int nSeams)
{
    return imgData;
}

cinder::Surface GradientSeamCarver::addHorizontalSeam(cinder::Surface imgData, int nSeams)
{
    return imgData;
}

// Get Seams in all 4 direction
//========================================================================
ImageSeam* GradientSeamCarver::getTopToBottom(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    
    for(int x=0; x<w; x++){
        //seamPath[x] = 0;
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
                if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue < verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
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
       // seamPath[(w-1)*h + x] = 0;
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
                if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue < verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                    
                }
                else{
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
        }for(int x=0; x<w; x++){
            verticalSeam[x].preValue = verticalSeam[x].totalValue;
           // printf("\nverticalSeam[%d].totalValue = %d",x,verticalSeam[x].totalValue);
        }
        //printf("\n\n");
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
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(verticalSeam[y].preValue <= verticalSeam[y-1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[y-1].preValue <= verticalSeam[y].preValue && verticalSeam[y-1].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
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
        verticalSeam[y].preValue = (imgData.getPixel(ci::Vec2i(w-1,y))).b/255.0f;
        verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(w-1,y))).b/255.0f;
        verticalSeam[y].isVerticalSeam = false;
        verticalSeam[y].isForwardSeam = false;
        verticalSeam[y].seamPath = seamPath;
    }
    
    //start with row 1 from right
    for(int x=w-2; x>=0; x--){
        for(int y=0; y<h; y++){
            
            if(y==0){//First Row
                if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(verticalSeam[y].preValue <= verticalSeam[y-1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                    
                }
                else{
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Rows
                if(verticalSeam[y-1].preValue <= verticalSeam[y].preValue && verticalSeam[y-1].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue = (imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
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


///

ImageSeam* GradientSeamCarver::getTopToBottom(Mat img)
{
    int w = img.rows;
    int h = img.cols;
    uchar* raster = img.data;
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    
    for(int x=0; x<w; x++){
        //seamPath[x] = 0;
        verticalSeam[x] = ImageSeam();
        verticalSeam[x].index = x;
        verticalSeam[x].preValue = raster[x]/255.0f;
        verticalSeam[x].totalValue = raster[x]/255.0f;
        verticalSeam[x].isVerticalSeam = true;
        verticalSeam[x].isForwardSeam = true;
        verticalSeam[x].seamPath = seamPath;
    }
    
    for(int y=1; y<h; y++){
        for(int x=0; x<w; x++){
            
            if(x==0){//First Column
                if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x].preValue;
                    seamPath[y * w + x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue < verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue = raster[y * w + x] + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int x=0; x<w; x++){
            verticalSeam[x].preValue = verticalSeam[x].totalValue;
            printf("\nverticalSeam[x].totalValue = %f",verticalSeam[x].totalValue);
        }
        printf("\n\n");
    }
    return verticalSeam;
}


ImageSeam* GradientSeamCarver::getBottomToTop(Mat img)
{
    int w = img.rows;
    int h = img.cols;
    uchar* raster = img.data;
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    /*
    for(int x=0; x<w; x++){
        // seamPath[(w-1)*h + x] = 0;
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
                if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (x == w-1){ //Last Column
                if(verticalSeam[x].preValue < verticalSeam[x-1].preValue){
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                    
                }
                else{
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[x-1].preValue <= verticalSeam[x].preValue && verticalSeam[x-1].preValue <= verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[x].preValue < verticalSeam[x+1].preValue){
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[x].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[x+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
        }
        for(int x=0; x<w; x++){
            verticalSeam[x].preValue = verticalSeam[x].totalValue;
        }
    }
     */
    return verticalSeam;
}


ImageSeam* GradientSeamCarver::getLeftToRight(Mat img)
{
    int w = img.rows;
    int h = img.cols;
    uchar* raster = img.data;
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[h];
    /*
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
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(verticalSeam[y].preValue <= verticalSeam[y-1].preValue){
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(verticalSeam[y-1].preValue <= verticalSeam[y].preValue && verticalSeam[y-1].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int y=0; y<h; y++){
            verticalSeam[y].preValue = verticalSeam[y].totalValue;
        }
    }
     */
    return verticalSeam;
}


ImageSeam* GradientSeamCarver::getRightToLeft(Mat img)
{
    
    int w = img.rows;
    int h = img.cols;
    uchar* raster = img.data;
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[h];
    /*
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
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(verticalSeam[y].preValue <= verticalSeam[y-1].preValue){
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                    
                }
                else{
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Rows
                if(verticalSeam[y-1].preValue <= verticalSeam[y].preValue && verticalSeam[y-1].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(verticalSeam[y].preValue <= verticalSeam[y+1].preValue){
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    verticalSeam[y].totalValue += imgData.getPixel(ci::Vec2i(x,y))).b + verticalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
            
        }
        for(int y=0; y<h; y++){
            verticalSeam[y].preValue = verticalSeam[y].totalValue;
        }
    }
     */
    return verticalSeam;
}

#endif
