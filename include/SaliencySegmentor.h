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
//#include "segment-image.h"

class SaliencySegmentor {
public:
    SaliencySegmentor();
    
    void setSegmentationParameters(float sigma, float k, float minSize);
    cinder::Surface segmentImage(cinder::Surface imgData);
    
private:
    float sigma_Seg;
    float k_Seg;
    int minSize_Seg;
};

#endif /* defined(__SaliencyRetarget__SaliencySegmentor__) */
