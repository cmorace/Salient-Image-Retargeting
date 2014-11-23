//
//  SaliencySegmentor.h
//  SaliencyRetarget
//
//  Created by Charles Morace on 10/28/14.
//
//

#ifndef __SaliencyRetarget__SaliencySegmentor__
#define __SaliencyRetarget__SaliencySegmentor__

#include <stdio.h>
#include <algorithm>
#include "segment-image.h"
#include "opencv2/opencv.hpp"
#include "CinderOpenCV.h"

class SaliencySegmentor {
public:
    enum SaliencyMethod
    {
        Sobel,
        Scharr,
        Canny, //todo
        Undefined
    };
    SaliencyMethod edgeDetect = SaliencyMethod::Sobel;
    unsigned int scale = 1;
    unsigned int delta = 0;
    
    SaliencySegmentor();
    void setSegmentationParameters(float sigma, float k, float minSize);
    cinder::Surface getSegmentedImage(cinder::Surface imgData);
    cinder::Surface getSaliencyMap(cinder::Surface imgData, SaliencyMethod edgeDetect);
    cinder::Surface getSalientSegmentedImage(cinder::Surface imgData);

    cinder::Surface mergeCurrentSegmentAndSaliency();
    float sigma_Seg;
    float k_Seg;
    int minSize_Seg;
    int nbOfSegments = 0;
    double segTime = 0.0;
    double saliencyTime = 0.0;
    
private:
    struct SaliencyScore
    {
        unsigned long totalScore;
        unsigned long totalPixels;
        float normalScore;
    };
    
    cinder::Timer* timer;
    cv::Mat currentGradientCV;
    cinder::Surface currentGradientSurface;
    universe* currentSegmentGraph;
    image<rgb> *segmentImage(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    image<rgb> *segmentColorImage(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    cv::Mat mergeSegmentedAndSaliency(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    cv::Mat getSaliencyMapCV(cinder::Surface imgData, SaliencyMethod edgeDetect);
};

SaliencySegmentor::SaliencySegmentor()
{
    sigma_Seg = 0.8f;
    k_Seg = 800;
    minSize_Seg = 100;
    timer = new cinder::Timer(false);
    currentSegmentGraph = new universe(0);
}

cinder::Surface SaliencySegmentor::getSegmentedImage(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    
    // assuming rgb format no alpha channel
    image<rgb>* input = new image<rgb>(w,h,false);
    input->data = (rgb*)(imgData.getData());
    
    timer->start();
    image<rgb>* output = segmentImage(input, sigma_Seg, k_Seg, minSize_Seg, &nbOfSegments);
    timer->stop();
    segTime = timer->getSeconds();
    
    // assuming rgb format no alpha channel
    cinder::Surface result = cinder::Surface((uchar*)(output->data), w,h, w*3, cinder::SurfaceChannelOrder::RGB);
    delete output;
    return result;
}


ci::Surface SaliencySegmentor::getSaliencyMap(cinder::Surface imgData, SaliencyMethod edgeDetect)
{
    currentGradientCV = getSaliencyMapCV(imgData, edgeDetect);
    return cinder::Surface( cinder::fromOcv( currentGradientCV ) );
}

cv::Mat SaliencySegmentor::getSaliencyMapCV(cinder::Surface imgData, SaliencyMethod edgeDetect)
{
    cv::Mat src, src_gray;
    int ddepth = CV_16S;
    
    src = cinder::toOcv(imgData);
    timer->start();
    GaussianBlur( src, src, cv::Size(3,3), 0, 0, cv::BORDER_DEFAULT );
    
    /// Convert it to gray
    cvtColor( src, src_gray, CV_RGB2GRAY );
    
    /// Generate grad_x and grad_y
    cv::Mat grad_x, grad_y;
    cv::Mat abs_grad_x, abs_grad_y;
    
    /// Gradient X
    if(edgeDetect == SaliencyMethod::Sobel){
        cv::Sobel( src_gray, grad_x, ddepth, 1, 0, 3, scale, delta, cv::BORDER_DEFAULT );
    }
    else if(edgeDetect == SaliencyMethod::Scharr){
        cv::Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, cv::BORDER_DEFAULT );
    }
    //
    convertScaleAbs( grad_x, abs_grad_x );
    
    /// Gradient Y
    if(edgeDetect == SaliencyMethod::Sobel){
        cv::Sobel( src_gray, grad_y, ddepth, 0, 1, 3, scale, delta, cv::BORDER_DEFAULT );
    }
    else if(edgeDetect == SaliencyMethod::Scharr){
        cv::Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, cv::BORDER_DEFAULT );
    }
    convertScaleAbs( grad_y, abs_grad_y );
    currentGradientCV.release();
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, currentGradientCV);
    return currentGradientCV;
}


cinder::Surface SaliencySegmentor::getSalientSegmentedImage(cinder::Surface imgData)
{
    printf("\ngetSalientSegmentedImage");
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    // assuming rgb format no alpha channel
    image<rgb>* input = new image<rgb>(w,h,false);
    
    input->data = (rgb*)(imgData.getData());
    cv::Mat output = mergeSegmentedAndSaliency(input, sigma_Seg, k_Seg, minSize_Seg, &nbOfSegments);
    //delete input;
    cinder::Surface result = cinder::Surface(cinder::fromOcv(output));
    return result;
}


cv::Mat SaliencySegmentor::mergeSegmentedAndSaliency(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs)
{
    int width = im->width();
    int height = im->height();
    
    image<float> *r = new image<float>(width, height);
    image<float> *g = new image<float>(width, height);
    image<float> *b = new image<float>(width, height);
    
    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            imRef(r, x, y) = imRef(im, x, y).r;
            imRef(g, x, y) = imRef(im, x, y).g;
            imRef(b, x, y) = imRef(im, x, y).b;
        }
    }
    image<float> *smooth_r = smooth(r, sigma);
    image<float> *smooth_g = smooth(g, sigma);
    image<float> *smooth_b = smooth(b, sigma);
    delete r;
    delete g;
    delete b;
    
    // build graph
    edge *edges = new edge[width*height*4];
    int num = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (x < width-1) {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
                num++;
            }
            
            if (y < height-1) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + x;
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
                num++;
            }
            
            if ((x < width-1) && (y < height-1)) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
                num++;
            }
            
            if ((x < width-1) && (y > 0)) {
                edges[num].a = y * width + x;
                edges[num].b = (y-1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
                num++;
            }
        }
    }
    delete smooth_r;
    delete smooth_g;
    delete smooth_b;
    
    // segment
    universe *u = segment_graph(width*height, num, edges, c);
    
    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    *num_ccs = u->num_sets();
    
    std::map<unsigned int, SaliencyScore> patchIndices;
    
    // get gradient from OpenCV
    uchar* raster = currentGradientCV.data;
    //
    int w = currentGradientCV.cols;
    int h = currentGradientCV.rows;
    for(int x=0; x<w; x++)
    {
        for(int y=0; y<h; y++)
        {
            uchar val = raster[y * w + x];
            int segmentIndex = u->find(y * w + x);
            patchIndices[segmentIndex].totalScore += val;
            patchIndices[segmentIndex].totalPixels += 1;
        }
    }
    
    SaliencyScore firstScore = patchIndices.begin()->second;
    float highScore = 0;
    float lowScore = 1.f*firstScore.totalScore/firstScore.totalPixels;
    
    for(std::map<unsigned int, SaliencyScore>::iterator patchIter = patchIndices.begin(); patchIter != patchIndices.end(); patchIter++)
    {
        SaliencyScore score = patchIter->second;
        score.normalScore = 1.f*score.totalScore/score.totalPixels;
        
        if(score.normalScore > highScore){
            highScore = score.normalScore;
        }
        else if(score.normalScore < lowScore){
            lowScore = score.normalScore;
        }
        patchIter->second = score;
    }
    
    // norm saliency values [0,1]
    float range = highScore - lowScore;
    for(std::map<unsigned int, SaliencyScore>::iterator patchIter = patchIndices.begin(); patchIter != patchIndices.end(); patchIter++)
    {
        SaliencyScore score = patchIter->second;
        score.normalScore = (score.normalScore - lowScore)/range;
        patchIter->second = score;
    }
    
    cv::Mat result = cv::Mat(h,w,CV_8U);
    for(int x=0; x<w; x++)
    {
        for(int y=0; y<h; y++)
        {
            result.at<uchar>(cv::Point(x,y)) = (uchar)(255*patchIndices[u->find(y * w + x)].normalScore);
        }
    }
    
    delete currentSegmentGraph;
    currentSegmentGraph = u;
    return result;
}

image<rgb> * SaliencySegmentor::segmentImage(image<rgb> *im, float sigma, float c, int min_size,
                                                    int *num_ccs)
{
    //return segment_image(im, sigma, c, min_size, num_ccs);
    
    int width = im->width();
    int height = im->height();
    
    image<float> *r = new image<float>(width, height);
    image<float> *g = new image<float>(width, height);
    image<float> *b = new image<float>(width, height);
    
    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            imRef(r, x, y) = imRef(im, x, y).r;
            imRef(g, x, y) = imRef(im, x, y).g;
            imRef(b, x, y) = imRef(im, x, y).b;
        }
    }
    image<float> *smooth_r = smooth(r, sigma);
    image<float> *smooth_g = smooth(g, sigma);
    image<float> *smooth_b = smooth(b, sigma);
    delete r;
    delete g;
    delete b;
    
    // build graph
    edge *edges = new edge[width*height*4];
    int num = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (x < width-1) {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
                num++;
            }
            
            if (y < height-1) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + x;
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
                num++;
            }
            
            if ((x < width-1) && (y < height-1)) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
                num++;
            }
            
            if ((x < width-1) && (y > 0)) {
                edges[num].a = y * width + x;
                edges[num].b = (y-1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
                num++;
            }
        }
    }
    delete smooth_r;
    delete smooth_g;
    delete smooth_b;
    
    // segment
    
    universe *u = segment_graph(width*height, num, edges, c);
    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    *num_ccs = u->num_sets();
    
    //
    /* testing
    int nbPatches = currentSegmentGraph->num_sets();
    std::map<unsigned int, unsigned long> patchIndices;
    printf("\nnbPatches = %d",nbPatches);
    uchar* raster = currentGradientCV.data;
    
    for(int x=0; x<width; x++)
    {
        for(int y=0; y<height; y++)
        {
            int segmentIndex = u->find(y*height + x);
            patchIndices[segmentIndex]  += 1;
        }
    }
    
    printf("\npatchIndices size = %lu", patchIndices.size() );
    */
    ///////////////////
    
    image<rgb> *output = new image<rgb>(width, height);
    
    // pick random colors for each component
    rgb *colors = new rgb[width*height];
    for (int i = 0; i < width*height; i++){
        colors[i] = random_rgb();
    }
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            imRef(output, x, y) = colors[comp];
            //printf("\n index = %d",comp);
        }
    }
    
    delete [] colors;  
    //delete u;
    delete currentSegmentGraph;
    currentSegmentGraph = u;
    
    return output;
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
image<rgb> * SaliencySegmentor::segmentColorImage(image<rgb> *im, float sigma, float c, int min_size,
                                int *num_ccs) {
    
    int width = im->width();
    int height = im->height();
    
    image<float> *r = new image<float>(width, height);
    image<float> *g = new image<float>(width, height);
    image<float> *b = new image<float>(width, height);
    
    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            imRef(r, x, y) = imRef(im, x, y).r;
            imRef(g, x, y) = imRef(im, x, y).g;
            imRef(b, x, y) = imRef(im, x, y).b;
        }
    }
    image<float> *smooth_r = smooth(r, sigma);
    image<float> *smooth_g = smooth(g, sigma);
    image<float> *smooth_b = smooth(b, sigma);
    delete r;
    delete g;
    delete b;
    
    // build graph
    edge *edges = new edge[width*height*4];
    int num = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (x < width-1) {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
                num++;
            }
            
            if (y < height-1) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + x;
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
                num++;
            }
            
            if ((x < width-1) && (y < height-1)) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
                num++;
            }
            
            if ((x < width-1) && (y > 0)) {
                edges[num].a = y * width + x;
                edges[num].b = (y-1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
                num++;
            }
        }
    }
    delete smooth_r;
    delete smooth_g;
    delete smooth_b;
    
    // segment
    universe *u = segment_graph(width*height, num, edges, c);
    
    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    *num_ccs = u->num_sets();
    
    image<rgb> *output = new image<rgb>(width, height);
    rgb *colors = new rgb[u->num_sets()];
    int* setIndices =  new int[u->num_sets()];
    int setIndex = 1;
    setIndices[0] = u->find(0);
    colors[0] = imRef(im, 0, 0);
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            bool newComp = true;
            int colorIndex = 0;
            for (int i=0; i<setIndex; i++) {
                if(comp == setIndices[i]){
                    newComp = false;
                    colorIndex = i;
                    break;
                }
            }
            if(newComp){
                setIndices[setIndex] = comp;
                colors[setIndex] = imRef(im, x, y);
                colorIndex = setIndex;
                setIndex++;
                //printf("comp = %d\n",comp);
            }
            imRef(output, x, y) = colors[colorIndex];
        }
    }
    
    delete [] colors;
    delete [] setIndices;
    delete currentSegmentGraph;
    currentSegmentGraph = u;
    return output;
    
}

#endif /* defined(__SaliencyRetarget__SaliencySegmentor__) */
