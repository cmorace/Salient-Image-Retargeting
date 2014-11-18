//
//  MeshWarpRetargetter.h
//  ImageRetargeting
//
//  Created by Charles Morace on 11/18/14.
//
//

#ifndef ImageRetargeting_MeshWarpRetargetter_h
#define ImageRetargeting_MeshWarpRetargetter_h

#include "MeshWarper.h"

using namespace ph::warping;

class MeshWarpRetargetter {
    
public:
    MeshWarpRetargetter();
    
    //void setTexture(ci::gl::Texture texture);
    void draw(ci::gl::Texture texture);
    void resize(void);
    void setMesh(void);
    void setSize(const Area a);
    ci::gl::Texture	myTexture;
    MeshWarperRef meshWarper;
    unsigned int quadSize;
    
private:
    
};

MeshWarpRetargetter::MeshWarpRetargetter()
{
    printf("MeshWarpRetargetter");
    meshWarper = MeshWarperRef(new MeshWarper());
    meshWarper->enableEditMode();
    meshWarper->setNumControlX(3);
    meshWarper->setNumControlY(3);
    quadSize = 20;
}


void MeshWarpRetargetter::draw(ci::gl::Texture texture)
{
    WarpRef warp( meshWarper );
    warp->draw(texture);
}



void MeshWarpRetargetter::resize()
{
    meshWarper->resize();
}

void MeshWarpRetargetter::setMesh()
{
    meshWarper->setMesh(quadSize);
}

#endif
