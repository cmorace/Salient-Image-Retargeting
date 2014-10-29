//
//  SaliencySegmentor.cpp
//  SaliencyRetarget
//
//  Created by Charles Morace on 10/28/14.
//
//

#include "SaliencySegmentor.h"
#include "segment-image.h"


SaliencySegmentor::SaliencySegmentor()
{
    printf("constructor\n");
    sigma_Seg = 0.5f;
    k_Seg = 50;
    minSize_Seg = 10;
}

cinder::Surface SaliencySegmentor::segmentImage(cinder::Surface imgData)
{
    int w = imgData.getWidth();
    int h = imgData.getHeight();
    image<rgb>* input = new image<rgb>(w,h,true);
    
    input->data = (rgb*)imgData.getData();
    
    printf("%u\n",input->data[50].r);
    printf("%u\n",input->data[50].g);
    printf("%u\n",input->data[50].b);
    
    
    cinder::ColorAT<uchar> c  = imgData.getPixel(cinder::Vec2i(50,0));
    
    printf("%u\n",c.r);
    printf("%u\n",c.g);
    printf("%u\n",c.b);
    
    int num_ccs;
    image<rgb>* output = segment_image(input, sigma_Seg, k_Seg, minSize_Seg, &num_ccs);
    printf("nb segments = %d\n", num_ccs);
    
    printf("%u\n",output->data[50].r);
    printf("%u\n",output->data[50].g);
    printf("%u\n",output->data[50].b);
    
    // ussuming no alpha channel
    return cinder::Surface((uchar*)(output->data), w,h, w*3, cinder::SurfaceChannelOrder::RGB);
    
    
    //return cinder::Surface((uchar*)(input->data), w,h, w*3, cinder::SurfaceChannelOrder::RGB);
}
