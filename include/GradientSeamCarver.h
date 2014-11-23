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
#include "opencv2/opencv.hpp"
#include "CinderOpenCV.h"


class GradientSeamCarver {
public:
    enum EdgeDetection
    {
        Sobel,
        Scharr,
        Canny, //todo
        Undefined
    };
    EdgeDetection prevGradientMethod = EdgeDetection::Undefined;
    
    GradientSeamCarver();
    // demo 1
    cinder::Surface getGradientImage(cinder::Surface imgData, EdgeDetection edgeDetect = EdgeDetection::Sobel);
    cinder::Surface drawVerticalSeamsGradient();
    cinder::Surface drawHorizontalSeamsGradient();
    void drawCurrentSeam(cinder::Surface img);
    cinder::Surface deleteCurrentSeam(cinder::Surface img);
    cinder::Surface deleteVerticalSeam(cinder::Surface imgData);
    cinder::Surface deleteHorizontalSeam(cinder::Surface imgData);
    void startGradientTimer();
    void stopGradientTimer();
    
    void startCarveTimer();
    void stopCarveTimer();
    //cinder::Surface deleteMinEnergySeam(cinder::Surface imgData);
    //cinder::Surface addVerticalSeam(cinder::Surface imgData, int nSeams);
    //cinder::Surface addHorizontalSeam(cinder::Surface imgData, int nSeams);
    //
    unsigned int scale;
    unsigned int delta;
    unsigned int nbOfSeams;
    unsigned int newWidth;
    unsigned int newHeight;
    double gradTime;
    double seamTime;
    double carveTime;
    
private:
    
    struct ImageSeam
    {
        enum SeamDirection
        {
            LeftToRight,
            RightToLeft,
            TopToBottom,
            BottomToTop,
            Undefined
        };
        unsigned int index = 0;
        unsigned int preValue = 0;
        unsigned int totalValue = 0;
        SeamDirection direction = SeamDirection::Undefined;
        int* seamPath = nullptr;
        bool operator < (const ImageSeam& seam) const
        {
            return totalValue < seam.totalValue;
        }
    };
    
    cv::Mat getGradientImageCV(cinder::Surface imgData, EdgeDetection edgeDetect);
    ImageSeam* getVerticalSeams(cv::Mat imgData);
    ImageSeam* getHorizontalSeams(cv::Mat imgData);
    std::vector<ImageSeam> getSortedVerticalSeams(cv::Mat imgData);
    std::vector<ImageSeam> getSortedHorizontalSeams(cv::Mat imgData);
    void deleteCurrentSeams();
    void drawSeams(cinder::Surface imgData, int nSeams, std::vector<ImageSeam> seamVector);
    
    cinder::Surface carveSeam(cinder::Surface imgData, ImageSeam seam);
    cinder::Timer* timer;
    cv::Mat currentGradient;
    std::vector<ImageSeam> currentSeams;
};

GradientSeamCarver::GradientSeamCarver()
{
    timer = new cinder::Timer(false);
    scale = 1;
    delta = 0;
    gradTime = 0.0;
    seamTime = 0.0;
    carveTime = 0.0;
    nbOfSeams = 1;
    newWidth = 800;
    newHeight = 600;
}

void GradientSeamCarver::startCarveTimer()
{
    timer->start();
    carveTime = 0;

}

void GradientSeamCarver::stopCarveTimer()
{
    timer->stop();
    carveTime = timer->getSeconds();
}

void GradientSeamCarver::startGradientTimer()
{
    timer->start();
    gradTime = 0;
    
}

void GradientSeamCarver::stopGradientTimer()
{
    timer->stop();
    gradTime = timer->getSeconds();
}


cinder::Surface GradientSeamCarver::getGradientImage(cinder::Surface imgData, EdgeDetection edgeDetect)
{
    currentGradient = getGradientImageCV(imgData, edgeDetect);
    prevGradientMethod = edgeDetect;
    return cinder::Surface( cinder::fromOcv( currentGradient ) );
}

cv::Mat GradientSeamCarver::getGradientImageCV(cinder::Surface imgData, EdgeDetection edgeDetect)
{
    //this->isSobel = isSobel;
    cv::Mat src, src_gray;
    int ddepth = CV_16S;
    
    src = cinder::toOcv(imgData);
    GaussianBlur( src, src, cv::Size(3,3), 0, 0, cv::BORDER_DEFAULT );
    
    /// Convert it to gray
    cvtColor( src, src_gray, CV_RGB2GRAY );
    
    /// Generate grad_x and grad_y
    cv::Mat grad_x, grad_y;
    cv::Mat abs_grad_x, abs_grad_y;
    
    /// Gradient X
    if(edgeDetect == EdgeDetection::Sobel){
        cv::Sobel( src_gray, grad_x, ddepth, 1, 0, 3, scale, delta, cv::BORDER_DEFAULT );
    }
    else if(edgeDetect == EdgeDetection::Scharr){
        cv::Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, cv::BORDER_DEFAULT );
    }
    //
    convertScaleAbs( grad_x, abs_grad_x );
    
    /// Gradient Y
    if(edgeDetect == EdgeDetection::Sobel){
        cv::Sobel( src_gray, grad_y, ddepth, 0, 1, 3, scale, delta, cv::BORDER_DEFAULT );
    }
    else if(edgeDetect == EdgeDetection::Scharr){
        cv::Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, cv::BORDER_DEFAULT );
    }
    convertScaleAbs( grad_y, abs_grad_y );
    
    /// Total Gradient (approximate)
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, currentGradient);
    return currentGradient;
}


void GradientSeamCarver::drawCurrentSeam(cinder::Surface img)
{
    drawSeams(img,1,currentSeams);
}

void GradientSeamCarver::drawSeams(cinder::Surface imgData, int nSeams, std::vector<ImageSeam> seamVector)
{
    if(!seamVector.empty()){
        ci::ColorT<uchar>* r = new ci::ColorT<uchar>(255,0,0);
        int h = imgData.getHeight();
        int w = imgData.getWidth();
        for (std::vector<ImageSeam>::iterator it=seamVector.begin(); it!=seamVector.begin() + nSeams; ++it)
        {
            if(it->direction == ImageSeam::SeamDirection::TopToBottom){
                int x = it->index;
                for(int y=h-1; y>=0; y--){
                    imgData.setPixel(ci::Vec2i(x,y), *r);
                    x += it->seamPath[y * w + x];
                }
            }
            else if(it->direction == ImageSeam::SeamDirection::LeftToRight){
                int y = it->index;
                for(int x=w-1; x>=0; x--){
                    imgData.setPixel(ci::Vec2i(x,y), *r);
                    y += it->seamPath[y * w + x];
                }
            }
        }
    }
    
}

cinder::Surface GradientSeamCarver::deleteCurrentSeam(cinder::Surface img)
{
    return carveSeam(img,currentSeams.at(0));
}

void GradientSeamCarver::deleteCurrentSeams()
{
    if(!currentSeams.empty()){
        delete[] currentSeams.at(0).seamPath;
        currentSeams.clear();
    }
    
}

cinder::Surface GradientSeamCarver::carveSeam(cinder::Surface imgData, ImageSeam seam)
{
    int h = imgData.getHeight();
    int w = imgData.getWidth();
    cinder::Surface s;
    if(seam.direction == ImageSeam::SeamDirection::TopToBottom){
        int i = seam.index;
        s = cinder::Surface(w-1,h,false);
        for(int y=h-1; y>=0; y--){
            for (int x = 0; x<w-1; x++) {
                if(x<i){
                    s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y)));
                }
                else{
                    s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x+1,y)));
                }
            }
            i += seam.seamPath[y * w + i];
        }
    }
    else if(seam.direction == ImageSeam::SeamDirection::LeftToRight){
        int i = seam.index;
        s = cinder::Surface(w,h-1,false);
        for(int x=w-1; x>=0; x--){
            for (int y = 0; y<h-1; y++) {
                if(y<i){
                    s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y)));
                }
                else{
                    s.setPixel(cinder::Vec2i(x,y), imgData.getPixel(cinder::Vec2i(x,y+1)));
                }
            }
            i += seam.seamPath[i * w + x];
        }
    }
    return s;
}

std::vector<GradientSeamCarver::ImageSeam> GradientSeamCarver::getSortedVerticalSeams(cv::Mat img)
{
    int w = img.cols;
    ImageSeam* vertSeam = getVerticalSeams(img);
    std::vector<ImageSeam> c(vertSeam, vertSeam + w);
    std::sort(c.begin(), c.end());
    return c;
}

std::vector<GradientSeamCarver::ImageSeam> GradientSeamCarver::getSortedHorizontalSeams(cv::Mat img)
{
    int h = img.rows;
    ImageSeam* horizSeam = getHorizontalSeams(img);
    std::vector<ImageSeam> c(horizSeam, horizSeam + h);
    std::sort (c.begin(), c.end());
    return c;
}

cinder::Surface GradientSeamCarver::drawVerticalSeamsGradient()
{
    cinder::Surface result = cinder::fromOcv(currentGradient);
    deleteCurrentSeams();
    currentSeams = getSortedVerticalSeams(currentGradient);
    drawSeams(result, nbOfSeams, currentSeams);
    return result;
}

cinder::Surface GradientSeamCarver::drawHorizontalSeamsGradient()
{
    cinder::Surface result = cinder::fromOcv(currentGradient);
    deleteCurrentSeams();
    currentSeams = getSortedHorizontalSeams(currentGradient);
    drawSeams(result, nbOfSeams, currentSeams);
    return result;
}

cinder::Surface GradientSeamCarver::deleteVerticalSeam(cinder::Surface img)
{
    deleteCurrentSeams();
    currentGradient = getGradientImageCV(img, prevGradientMethod);
    currentSeams = getSortedVerticalSeams(currentGradient);
    img = carveSeam(img, currentSeams.at(0));
    return img;
}

cinder::Surface GradientSeamCarver::deleteHorizontalSeam(cinder::Surface img)
{
    deleteCurrentSeams();
    currentGradient = getGradientImageCV(img, prevGradientMethod);
    currentSeams = getSortedHorizontalSeams(currentGradient);
    img = carveSeam(img, currentSeams.at(0));
    return img;
}


/*
cinder::Surface GradientSeamCarver::deleteMinEnergySeam(cinder::Surface img){
    cv::Mat imgData = getGradientImageCV(img.clone(),this->isSobel);
    //
    timer->start();
    std::vector<ImageSeam> a = getSortedHorizontalSeams(imgData);
    std::vector<ImageSeam> b = getSortedVerticalSeams(imgData);
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

 */


GradientSeamCarver::ImageSeam* GradientSeamCarver::getVerticalSeams(cv::Mat img)
{
    int w = img.cols;
    int h = img.rows;
    uchar* raster = img.data;
    
    //remember to clean up after use
    int* seamPath = new int[w*h];
    ImageSeam* verticalSeam = new ImageSeam[w];
    
    for(int x=0; x<w; x++){
        verticalSeam[x] = ImageSeam();
        verticalSeam[x].direction = ImageSeam::SeamDirection::TopToBottom;
        verticalSeam[x].index = x;
        verticalSeam[x].preValue = raster[x];
        verticalSeam[x].totalValue = raster[x];
        verticalSeam[x].seamPath = seamPath;
        seamPath[x] = 0;
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
                if(verticalSeam[x-1].preValue < verticalSeam[x].preValue && verticalSeam[x-1].preValue < verticalSeam[x+1].preValue){
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
        }
    }
    return verticalSeam;
}

GradientSeamCarver::ImageSeam* GradientSeamCarver::getHorizontalSeams(cv::Mat img)
{
    int w = img.cols;
    int h = img.rows;
    uchar* raster = img.data;
    
    //remember to clean up after use
    int* seamPath = new int[w*h];
    ImageSeam* horizontalSeam = new ImageSeam[h];
    
    for(int y=0; y<h; y++){
        horizontalSeam[y] = ImageSeam();
        horizontalSeam[y].direction = ImageSeam::SeamDirection::LeftToRight;
        horizontalSeam[y].index = y;
        horizontalSeam[y].preValue = raster[y*w];
        horizontalSeam[y].totalValue = raster[y*w];
        horizontalSeam[y].seamPath = seamPath;
        seamPath[y*w] = 0;
    }
    
    for(int x=1; x<w; x++){
        for(int y=0; y<h; y++){
            
            if(y==0){//First Row
                if(horizontalSeam[y].preValue < horizontalSeam[y+1].preValue){
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
                
            }
            else if (y == h-1){ //Last Row
                if(horizontalSeam[y].preValue < horizontalSeam[y-1].preValue){
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
            }
            else{ //Interior Column
                if(horizontalSeam[y-1].preValue < horizontalSeam[y].preValue && horizontalSeam[y-1].preValue < horizontalSeam[y+1].preValue){
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y-1].preValue;
                    seamPath[y*w+x] = -1;
                }
                else if(horizontalSeam[y].preValue < horizontalSeam[y+1].preValue){
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y].preValue;
                    seamPath[y*w+x] = 0;
                }
                else{
                    horizontalSeam[y].totalValue = raster[y * w + x] + horizontalSeam[y+1].preValue;
                    seamPath[y*w+x] = 1;
                }
            }
        }
        for(int y=0; y<h; y++){
            horizontalSeam[y].preValue = horizontalSeam[y].totalValue;
        }
    }
    return horizontalSeam;
}


// Get Seams in all 4 direction
//========================================================================
/*
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
*/

///




/*
ImageSeam* GradientSeamCarver::getBottomToTop(cv::Mat img)
{
    int w = img.rows;
    int h = img.cols;
    uchar* raster = img.data;
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
 
    return verticalSeam;
}
*/



/*

ImageSeam* GradientSeamCarver::getRightToLeft(cv::Mat img)
{
    
    int w = img.rows;
    int h = img.cols;
    uchar* raster = img.data;
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
 
    return verticalSeam;
}

*/


#endif
