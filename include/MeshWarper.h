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
    void setTexture(cinder::gl::Texture t);
    cinder::gl::Texture currentTexture;
    bool isdrawWireFrame = true;
    bool isdrawControlPoints = true;

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
    if( isdrawWireFrame ) {
        glDisable( GL_TEXTURE_2D );
        
        // draw wireframe
        gl::color( ColorA(1, 1, 1, 0.3f) );
        gl::enableAlphaBlending();
        gl::enableWireframe();
        gl::draw( mVboMesh );
        gl::disableAlphaBlending();
        gl::disableWireframe();
        
        if(isdrawControlPoints) {
            // draw control points
            for(unsigned i=0;i<mPoints.size();i++)
                drawControlPoint( getControlPoint(i) * mWindowSize, i == mSelected );
        }
    }
    
    // restore states
    glPopAttrib();
}
void MeshWarper::setMesh(int quadSize)
{
    mResolution = quadSize;
    createMesh(quadSize,quadSize);
}

void MeshWarper::createMesh(int quadResolutionX, int quadResolutionY)
{
    /*
    int totalVertices = VERTICES_X * VERTICES_Z;
    int totalQuads = ( VERTICES_X - 1 ) * ( VERTICES_Z - 1 );
    gl::VboMesh::Layout layout;
    layout.setStaticIndices();
    layout.setDynamicPositions();
    layout.setStaticTexCoords2d();
    mVboMesh = gl::VboMesh::create( numVertices, numQuads * 4, layout, GL_QUADS );
    
    // buffer our static data - the texcoords and the indices
    vector<uint32_t> indices;
    vector<Vec2f> texCoords;
    for( int x = 0; x < VERTICES_X; ++x ) {
        for( int z = 0; z < VERTICES_Z; ++z ) {
            // create a quad for each vertex, except for along the bottom and right edges
            if( ( x + 1 < VERTICES_X ) && ( z + 1 < VERTICES_Z ) ) {
                indices.push_back( (x+0) * VERTICES_Z + (z+0) );
                indices.push_back( (x+1) * VERTICES_Z + (z+0) );
                indices.push_back( (x+1) * VERTICES_Z + (z+1) );
                indices.push_back( (x+0) * VERTICES_Z + (z+1) );
            }
            // the texture coordinates are mapped to [0,1.0)
            texCoords.push_back( Vec2f( x / (float)VERTICES_X, z / (float)VERTICES_Z ) );
        }
    }
    
    mVboMesh->bufferIndices( indices );
    mVboMesh->bufferTexCoords2d( 0, texCoords );
    
    // make a second Vbo that uses the statics from the first
    mVboMesh2 = gl::VboMesh::create( totalVertices, totalQuads * 4, mVboMesh->getLayout(), GL_QUADS, &mVboMesh->getIndexVbo(), &mVboMesh->getStaticVbo(), NULL );
    mVboMesh2->setTexCoordOffset( 0, mVboMesh->getTexCoordOffset( 0 ) );
    
    mTexture = gl::Texture::create( loadImage( loadResource( RES_IMAGE ) ) );
    // convert from number of quads to number of vertices
    
    ++resolutionX;	++resolutionY;
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
     */
}




#endif
