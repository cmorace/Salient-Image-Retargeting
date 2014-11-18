//
//  MeshWarper.h
//  ImageRetargeting
//
//  Created by Charles Morace on 11/10/14.
//
//

#ifndef ImageRetargeting_MeshWarper_h
#define ImageRetargeting_MeshWarper_h

#include "WarpBilinear.h"

using namespace ph::warping;
using namespace cinder;

typedef std::shared_ptr<class MeshWarper>	MeshWarperRef;

class MeshWarper : public WarpBilinear
{
    
public:
    void setMesh(int quadSize);

protected:
    void draw(bool controls = true);
    void createMesh(int resolutionX=20, int resolutionY=20);
};


void MeshWarper::draw(bool controls)
{
    createBuffers();
    
    if(!mVboMesh) return;
    
    // save current texture mode, drawing color, line width and depth buffer state
    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
    
    gl::disableDepthRead();
    gl::disableDepthWrite();
    
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
    // adjust brightness
    if( mBrightness < 1.f )
    {
        ColorA currentColor;
        glGetFloatv(GL_CURRENT_COLOR, currentColor.ptr());
        
        ColorA drawColor = mBrightness * currentColor;
        drawColor.a = currentColor.a;
        
        gl::color( drawColor );
    }
    
    // draw textured mesh
    gl::draw( mVboMesh );
    
    // draw edit interface
    if( isEditModeEnabled() ) {
        glDisable( GL_TEXTURE_2D );
        
        // draw wireframe
        gl::color( ColorA(1, 1, 1, 0.5f) );
        gl::enableAlphaBlending();
        gl::enableWireframe();
        gl::draw( mVboMesh );
        gl::disableAlphaBlending();
        gl::disableWireframe();
        
        if(controls) {
            // draw control points
           // for(unsigned i=0;i<mPoints.size();i++)
            //    drawControlPoint( getControlPoint(i) * mWindowSize, i == mSelected );
        }
    }
    
    // restore states
    glPopAttrib();
}
void MeshWarper::setMesh(int quadSize)
{
    printf("\nsetMesh");
    printf("\nsetMesh = %d", quadSize);
    mResolution = quadSize;
    setNumControlX(getWidth()/quadSize);
    setNumControlY(getHeight()/quadSize);
    createMesh(mResolution,mResolution);
}

void MeshWarper::createMesh(int resolutionX, int resolutionY)
{
    // convert from number of quads to number of vertices
    printf("\ncreateMesh");
    ++resolutionX;	++resolutionY;
    
    // find a value for resolutionX and resolutionY that can be
    // evenly divided by mControlsX and mControlsY
    /*
    if(mControlsX < resolutionX) {
        int dx = (resolutionX-1) % (mControlsX-1);
        if(dx >= (mControlsX/2)) dx -= (mControlsX-1);
        resolutionX -= dx;
    } else {
        resolutionX = mControlsX;
    }
    
    if(mControlsY < resolutionY) {
        int dy = (resolutionY-1) % (mControlsY-1);
        if(dy >= (mControlsY/2)) dy -= (mControlsY-1);
        resolutionY -= dy;
    } else {
        resolutionY = mControlsY;
    }
    */
    //test
    resolutionX = mControlsX;
    resolutionY = mControlsY;
    //
    
    //
    mResolutionX = resolutionX;
    mResolutionY = resolutionY;
    
    //
    int numVertices = (resolutionX * resolutionY);
    int numQuads = (resolutionX - 1) * (resolutionY - 1);
    int numIndices = numQuads * 4;
    
    //
    gl::VboMesh::Layout	layout;
    layout.setStaticIndices();
    layout.setStaticTexCoords2d();
    layout.setDynamicPositions();
    
    //
    mVboMesh = gl::VboMesh( numVertices, numIndices, layout, GL_QUADS );
    if(!mVboMesh) return;
    
    // buffer static data
    int i = 0;
    int j = 0;
    std::vector<uint32_t>	indices( numIndices );
    std::vector<Vec2f>		texCoords( numVertices );
    for(int x=0; x<resolutionX; ++x) {
        for(int y=0; y<resolutionY; ++y) {
            // index
            if( ((x+1)<resolutionX) && ((y+1)<resolutionY) ) {
                indices[i++] = (x+0) * resolutionY + (y+0);
                indices[i++] = (x+1) * resolutionY + (y+0);
                indices[i++] = (x+1) * resolutionY + (y+1);
                indices[i++] = (x+0) * resolutionY + (y+1);
            }
            // texCoords
            float tx = lerp<float, float>( mX1, mX2, x / (float)(resolutionX-1) );
            float ty = lerp<float, float>( mY1, mY2, y / (float)(resolutionY-1) );
            texCoords[j++] = Vec2f(tx, ty );
        }
    }
    mVboMesh.bufferIndices( indices );
    mVboMesh.bufferTexCoords2d( 0, texCoords );
    
    //
    mIsDirty = true;
}




#endif
