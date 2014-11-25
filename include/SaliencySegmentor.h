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
//#include "MeshWarpRetargetter.h"



class SaliencySegmentor {
public:
    
    
    enum SaliencyMethod
    {
        Sobel, //using gradient for saliency for now
        Scharr,
        Canny, //todo
        Undefined
    };
    SaliencyMethod edgeDetect = SaliencyMethod::Sobel;
    
    struct Patch
    {
        unsigned long totalScore;
        unsigned long totalPixels;
        float normalScore;
        std::vector<ci::Vec2i> pixels;
        std::vector<int> edges;
    };
    typedef std::map<unsigned int, Patch> PatchMap;
    typedef std::map<unsigned int, Patch>::iterator PatchMapIterator;
    
    unsigned int scale = 1;
    unsigned int delta = 0;
    
    SaliencySegmentor();
    cinder::Surface getSegmentedImage(cinder::Surface imgData);
    cinder::Surface getSegmentedColorImage(cinder::Surface imgData);
    cinder::Surface getSaliencyMap(cinder::Surface imgData, SaliencyMethod edgeDetect);
    cinder::Surface getSegmentedSalientImage(cinder::Surface imgData);
    //void setPatchCenters(void);
    const PatchMap getPatchMap();
    universe* getUniverse();

    float segBlurDeviation;
    float segNeighborThreshold;
    int segMinSize;
    int nbOfSegments = 0;
    double segTime = 0.0;
    double saliencyTime = 0.0;
    
private:
    PatchMap patchMap;
    universe* currentUniverse;
    
    struct ColorScore
    {
        unsigned long totalRedScore;
        unsigned long totalGreenScore;
        unsigned long totalBlueScore;
        unsigned long totalPixels;
        rgb averageColor;
    };
    
    struct EdgeGraph
    {
        edge* edgeSet;
        int size;
    };
    
    
    cinder::Timer* timer;
    cv::Mat currentGradientCV;
    cinder::Surface currentGradientSurface;
    image<rgb> *segmentImage(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    image<rgb> *segmentColorImage(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    cv::Mat mergeSegmentedAndSaliency(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    cv::Mat getSaliencyMapCV(cinder::Surface imgData, SaliencyMethod edgeDetect);
    EdgeGraph getEdges(image<rgb> *im, float sigma);
};

SaliencySegmentor::SaliencySegmentor()
{
    segBlurDeviation = 0.8f;
    segNeighborThreshold = 800;
    segMinSize = 100;
    timer = new cinder::Timer(false);
    currentUniverse = new universe(0);
}

cinder::Surface SaliencySegmentor::getSegmentedImage(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    
    // assuming rgb format no alpha channel
    image<rgb>* input = new image<rgb>(w,h,false);
    input->data = (rgb*)(imgData.getData());
    
    timer->start();
    image<rgb>* output = segmentImage(input, segBlurDeviation, segNeighborThreshold, segMinSize, &nbOfSegments);
    timer->stop();
    segTime = timer->getSeconds();
    cinder::Surface result = cinder::Surface((uchar*)(output->data), w,h, w*3, cinder::SurfaceChannelOrder::RGB);
    delete output;
    return result;
}

cinder::Surface SaliencySegmentor::getSegmentedColorImage(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    
    // assuming rgb format no alpha channel
    image<rgb>* input = new image<rgb>(w,h,false);
    input->data = (rgb*)(imgData.getData());
    
    timer->start();
    image<rgb>* output = segmentColorImage(input, segBlurDeviation, segNeighborThreshold, segMinSize, &nbOfSegments);
    timer->stop();
    segTime = timer->getSeconds();
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
    //Fake Saliency
    GaussianBlur( currentGradientCV, currentGradientCV, cv::Size(5,5), 10, 10, cv::BORDER_DEFAULT );
    GaussianBlur( currentGradientCV, currentGradientCV, cv::Size(5,5), 10, 10, cv::BORDER_DEFAULT );
    GaussianBlur( currentGradientCV, currentGradientCV, cv::Size(5,5), 10, 10, cv::BORDER_DEFAULT );
    GaussianBlur( currentGradientCV, currentGradientCV, cv::Size(5,5), 10, 10, cv::BORDER_DEFAULT );
    GaussianBlur( currentGradientCV, currentGradientCV, cv::Size(5,5), 10, 10, cv::BORDER_DEFAULT );
    GaussianBlur( currentGradientCV, currentGradientCV, cv::Size(5,5), 10, 10, cv::BORDER_DEFAULT );
    timer->stop();
    saliencyTime = timer->getSeconds();
    return currentGradientCV;
}


cinder::Surface SaliencySegmentor::getSegmentedSalientImage(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    // assuming rgb format no alpha channel
    image<rgb>* input = new image<rgb>(w,h,false);
    
    input->data = (rgb*)(imgData.getData());
    timer->start();
    cv::Mat output = mergeSegmentedAndSaliency(input, segBlurDeviation, segNeighborThreshold, segMinSize, &nbOfSegments);
    timer->stop();
    segTime = timer->getSeconds();
    cinder::Surface result = cinder::Surface(cinder::fromOcv(output));
    return result;
}

SaliencySegmentor::EdgeGraph SaliencySegmentor::getEdges(image<rgb> *im, float sigma)
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
    
    EdgeGraph edgeGraph;
    edgeGraph.edgeSet = edges;
    edgeGraph.size = num;
    return edgeGraph;
}


cv::Mat SaliencySegmentor::mergeSegmentedAndSaliency(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs)
{
    int w = im->width();
    int h = im->height();
    
    EdgeGraph e = getEdges(im, sigma);
    edge *edges = e.edgeSet;
    int num = e.size;
    
    // segment
    universe *u = segment_graph(w*h, num, edges, c);
    
    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    patchMap.clear();
    
    // get gradient from OpenCV
    uchar* raster = currentGradientCV.data;
    for(int x=0; x<w; x++)
    {
        for(int y=0; y<h; y++)
        {
            uchar val = raster[y * w + x];
            int segmentIndex = u->find(y * w + x);
            patchMap[segmentIndex].totalScore += val;
            patchMap[segmentIndex].totalPixels += 1;
            patchMap[segmentIndex].pixels.push_back(ci::Vec2i(x,y));
        }
    }
    Patch firstScore = patchMap.begin()->second;
    float highScore = 0;
    float lowScore = 1.f*firstScore.totalScore/firstScore.totalPixels;
    for(PatchMapIterator patchIter = patchMap.begin(); patchIter != patchMap.end(); patchIter++)
    {
        Patch patch = patchIter->second;
        patch.normalScore = 1.f*patch.totalScore/patch.totalPixels;
        
        if(patch.normalScore > highScore){
            highScore = patch.normalScore;
        }
        else if(patch.normalScore < lowScore){
            lowScore = patch.normalScore;
        }
        patchIter->second = patch;
    }
    
    // norm saliency values [0.1,1]
    float range = highScore - lowScore;
    for(PatchMapIterator patchIter = patchMap.begin(); patchIter != patchMap.end(); patchIter++)
    {
        Patch patch = patchIter->second;
        patch.normalScore = (0.9f*(patch.normalScore - lowScore)/range) + 0.1f;
        patchIter->second = patch;
    }
    
    // make saliency map
    cv::Mat result = cv::Mat(h,w,CV_8U);
    for(int x=0; x<w; x++)
    {
        for(int y=0; y<h; y++)
        {
            result.at<uchar>(cv::Point(x,y)) = (uchar)(255*patchMap[u->find(y * w + x)].normalScore);
        }
    }
    *num_ccs = u->num_sets();
    delete currentUniverse;
    currentUniverse = u;
    return result;
}

const SaliencySegmentor::PatchMap SaliencySegmentor::getPatchMap()
{
    return patchMap;
}

universe* SaliencySegmentor::getUniverse()
{
    return currentUniverse;
}


image<rgb> * SaliencySegmentor::segmentImage(image<rgb> *im, float sigma, float c, int min_size,
                                                    int *num_ccs)
{
    // build graph
    EdgeGraph e = getEdges(im, sigma);
    edge* edges = e.edgeSet;
    int num = e.size;
    
    // segment
    int w = im->width();
    int h = im->height();
    universe *u = segment_graph(w*h, num, edges, c);
    
    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    image<rgb> *output = new image<rgb>(w, h);
    
    // pick random colors for each component
    rgb *colors = new rgb[w*h];
    for (int i = 0; i < w*h; i++){
        colors[i] = random_rgb();
    }
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int comp = u->find(y * w + x);
            imRef(output, x, y) = colors[comp];
        }
    }
    delete [] colors;
    *num_ccs = u->num_sets();
    delete u;
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
                                int *num_ccs)
{
    // build graph
    EdgeGraph e = getEdges(im, sigma);
    edge* edges = e.edgeSet;
    int num = e.size;
    int w = im->width();
    int h = im->height();
    
    // segment
    universe *u = segment_graph(w*h, num, edges, c);
    
    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    std::map<unsigned int, ColorScore> patchIndices;
    for(int x=0; x<w; x++)
    {
        for(int y=0; y<h; y++)
        {
            rgb color = imRef(im, x, y);
            int segmentIndex = u->find(y * w + x);
            patchIndices[segmentIndex].totalRedScore += color.r;
            patchIndices[segmentIndex].totalGreenScore += color.g;
            patchIndices[segmentIndex].totalBlueScore += color.b;
            patchIndices[segmentIndex].totalPixels += 1;
        }
    }
    for(std::map<unsigned int, ColorScore>::iterator patchIter = patchIndices.begin(); patchIter != patchIndices.end(); patchIter++)
    {
        ColorScore score = patchIter->second;
        rgb color;
        color.r = score.totalRedScore/score.totalPixels;
        color.g = score.totalGreenScore/score.totalPixels;
        color.b = score.totalBlueScore/score.totalPixels;
        score.averageColor = color;
        patchIter->second = score;
    }
    image<rgb> *output = new image<rgb>(w, h);
    for(int x=0; x<w; x++)
    {
        for(int y=0; y<h; y++)
        {
            imRef(output, x, y) = patchIndices[u->find(y * w + x)].averageColor;
        }
    }
    *num_ccs = u->num_sets();
    delete u;
    return output;
}

#endif /* defined(__SaliencyRetarget__SaliencySegmentor__) */
