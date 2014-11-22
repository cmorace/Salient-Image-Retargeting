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
#include "segment-image.h"
#include "opencv2/opencv.hpp"
#include "CinderOpenCV.h"

class SaliencySegmentor {
public:
    enum EdgeDetection
    {
        Sobel,
        Scharr,
        Canny, //todo
        Undefined
    };
    EdgeDetection edgeDetect = EdgeDetection::Sobel;
    unsigned int scale = 1;
    unsigned int delta = 0;
    
    SaliencySegmentor();
    void setSegmentationParameters(float sigma, float k, float minSize);
    cinder::Surface segmentImage(cinder::Surface imgData);
    cinder::Surface getSaliencyMap(cinder::Surface imgData);
    float sigma_Seg;
    float k_Seg;
    int minSize_Seg;
    int nbOfSegments;
    double segTime;
    double saliencyTime;
    
private:
    cinder::Timer* timer;
    cv::Mat currentGradient;
    image<rgb> *segment_color_image(image<rgb> *im, float sigma, float c, int min_size, int *num_ccs);
    cv::Mat getSaliencyMapCV(cinder::Surface imgData, EdgeDetection edgeDetect);
};

SaliencySegmentor::SaliencySegmentor()
{
    sigma_Seg = 0.8f;
    k_Seg = 1000;
    minSize_Seg = 200;
    nbOfSegments = 0;
    segTime = 0.0;
    timer = new cinder::Timer(false);
}

cinder::Surface SaliencySegmentor::segmentImage(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    
    // assuming rgb format no alpha channel
    image<rgb>* input = new image<rgb>(w,h,false);
    //
    
    input->data = (rgb*)(imgData.getData());
    timer->start();
    image<rgb>* output = segment_color_image(input, sigma_Seg, k_Seg, minSize_Seg, &nbOfSegments);
    timer->stop();
    segTime = timer->getSeconds();
    
    // assuming rgb format no alpha channel
    return cinder::Surface((uchar*)(output->data), w,h, w*3, cinder::SurfaceChannelOrder::RGB);
    //


}

ci::Surface SaliencySegmentor::getSaliencyMap(cinder::Surface imgData)
{
    currentGradient = getSaliencyMapCV(imgData, edgeDetect);
    return cinder::Surface( cinder::fromOcv( currentGradient ) );
}

cv::Mat SaliencySegmentor::getSaliencyMapCV(cinder::Surface imgData, EdgeDetection edgeDetect)
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
    timer->stop();
    saliencyTime = timer->getSeconds();
    return currentGradient;
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
image<rgb> * SaliencySegmentor::segment_color_image(image<rgb> *im, float sigma, float c, int min_size,
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
    delete u;
    
    return output;
    
}

#endif /* defined(__SaliencyRetarget__SaliencySegmentor__) */
