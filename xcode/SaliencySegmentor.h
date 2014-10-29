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

class SaliencySegmentor {
public:
    SaliencySegmentor();
    void setSegmentationParameters(float sigma, float k, float minSize);
    cinder::Surface segmentImage(cinder::Surface imgData);
    float sigma_Seg;
    float k_Seg;
    int minSize_Seg;
    int nbOfSegments;
    double segTime;
    
private:
    cinder::Timer* timer;
};

SaliencySegmentor::SaliencySegmentor()
{
    sigma_Seg = 1.0f;
    k_Seg = 315;
    minSize_Seg = 530;
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
    input->data = (rgb*)(imgData.getData());
    timer->start();
    image<rgb>* output = segment_image(input, sigma_Seg, k_Seg, minSize_Seg, &nbOfSegments);
    timer->stop();
    segTime = timer->getSeconds();
    // assuming rgb format no alpha channel
    return cinder::Surface((uchar*)(output->data), w,h, w*3, cinder::SurfaceChannelOrder::RGB);
}

#endif /* defined(__SaliencyRetarget__SaliencySegmentor__) */
