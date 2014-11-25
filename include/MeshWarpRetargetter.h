//
//  MeshWarpRetargetter.h
//  ImageRetargeting
//
//  Created by Charles Morace on 11/18/14.
//
//

#ifndef ImageRetargeting_MeshWarpRetargetter_h
#define ImageRetargeting_MeshWarpRetargetter_h

//#include "SaliencySegmentor.h"
#include "cinder/Rand.h"

class MeshWarpRetargetter {
    
public:
    
    MeshWarpRetargetter();
    void initMesh(unsigned int imgWidth, unsigned int imgHeight);
    void initMesh(unsigned int imgWidth, unsigned int imgHeight, SaliencySegmentor* segmentor);
    //void setEdges(SaliencySegmentor::PatchMap patchMap);
    void drawMesh(ci::gl::Texture texture);
    void drawEdges(ci::gl::Texture texture);
    void getRigidTransformationTerms(void);
    //void drawCenterEdge(ci::gl::Texture texture);
    //void intitRepresentativeEdge();
    //void resize(void);
    //void setMesh();
    
    bool isDrawingWireFrame = true;
    bool isDrawingVertices = false;
    
    bool isLinearInterpolated = true;
    int quadSize = 50; // this is desired quad width set by user
    
    
private:
    struct MeshEdge{
        ci::Vec2f a,b; // store endpoints vertices
    };
    std::vector<MeshEdge>	meshEdges;
    typedef std::vector<MeshEdge>::iterator MeshEdgeIterator;
    
    struct MeshQuad{
        ci::Vec2f t,b,l,r; // store quad vertices
    };
    std::vector<MeshQuad>	meshQuads;
    SaliencySegmentor::PatchMap currentPatchMap;
    MeshEdge getEdge(int x, int y);
    
    //Vec2f getPoint(int col, int row, int yRes) const;
    //Vec2f cubicInterpolate( const std::vector<Vec2f> &knots, float t ) const;
    //void setNumControlX(int xRes, int yRes);
    //void setNumControlY(int xRes, int yRes);
    
    ci::gl::VboMeshRef vboMesh;
    float quadWidth, quadHeight;  // this is actual quad width and height
    unsigned int numVertices;
    unsigned int numEdges;
    unsigned int numQuads;
    unsigned int numPatches;
    unsigned int m;
    unsigned int mPrime;
    unsigned int n;
    unsigned int nPrime;
};


MeshWarpRetargetter::MeshWarpRetargetter()
{
    printf("\nMeshWarpRetargetter");
}

void MeshWarpRetargetter::initMesh(unsigned int imgWidth, unsigned int imgHeight)
{
    quadWidth = quadHeight = quadSize;
    m = imgHeight;
    n = imgWidth;
    int const VERTICES_Y = m / quadHeight + 1;
    int const VERTICES_X = n / quadWidth + 1;
    
    numVertices = VERTICES_X * VERTICES_Y;
    numEdges = 2 * (VERTICES_Y - 1) * ( VERTICES_X - 1 ) + (VERTICES_X + VERTICES_Y - 2);
    numQuads = ( VERTICES_X - 1 ) * ( VERTICES_Y - 1 );
    
    printf("\nm = %d", m);
    printf("\nn = %d", n);
    printf("\nVERTICES_X = %d", VERTICES_X);
    printf("\nVERTICES_Y = %d", VERTICES_Y);
    printf("\nnumVertices = %d", numVertices);
    printf("\nnumEdges = %d", numEdges);
    printf("\nnumQuads = %d", numQuads);
    
    quadHeight = m / (VERTICES_Y - 1.f);
    quadWidth = n / (VERTICES_X - 1.f);
    
    printf("\n(quadWidth,quadWidth) = (%f,%f)",quadWidth,quadHeight);
    
    cinder::gl::VboMesh::Layout layout;
    layout.setStaticIndices();
    layout.setDynamicPositions();
    layout.setStaticTexCoords2d();
    vboMesh = ci::gl::VboMesh::create( numVertices, numQuads * 4, layout, GL_QUADS );
    
    // buffer our static data - the texcoords and the indices
    std::vector<uint32_t> indices;
    std::vector<ci::Vec2f> texCoords;
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    // buffer edge and quad indices for retargetting
    meshEdges.clear();
    meshQuads.clear();
    
    for( int x = 0; x < VERTICES_X; ++x ) {
        for( int y = 0; y < VERTICES_Y; ++y ) {
            // set Texture coordinates
            ci::Vec2f v = Vec2f( x / (VERTICES_X-1.f), y / (VERTICES_Y-1.f) );
            // the texture coordinates are mapped to [0,1]x[0,1]
            texCoords.push_back( v );
            // the vertex coordinates are mapped to [0,n-1]x[0,m-1]
            v.x *= n-1;
            v.y *= m-1;
            iter.setPosition(v.x, v.y, 0.0f );
            ++iter;
            
            // create a quad for each vertex, except for along the bottom and right edges
            if( ( x + 1 < VERTICES_X ) && ( y + 1 < VERTICES_Y ) ) {
                int topLeft = (x+0) * VERTICES_Y + (y+0);
                int topRight = (x+1) * VERTICES_Y + (y+0);
                int bottomRight = (x+1) * VERTICES_Y + (y+1);
                int bottomLeft = (x+0) * VERTICES_Y + (y+1);
                indices.push_back(topLeft);
                indices.push_back(topRight);
                indices.push_back(bottomRight);
                indices.push_back(bottomLeft);
                
                Vec2f vTopRight = Vec2f(v.x+quadWidth , v.y);
                Vec2f vBottomRight = Vec2f(v.x+quadWidth , v.y+quadHeight);
                Vec2f vBottomLeft = Vec2f(v.x , v.y+quadHeight);
                
                MeshEdge topEdge = {v,vTopRight};
                MeshEdge leftEdge = {v,vBottomLeft};
                meshEdges.push_back(topEdge);
                meshEdges.push_back(leftEdge);
                
                MeshQuad quad = {v,vTopRight,vBottomRight,vBottomLeft};
                meshQuads.push_back(quad);
            }
            else if(y+1 == VERTICES_Y && x+1 != VERTICES_X)
            {
                Vec2f vTopRight = Vec2f(v.x+quadWidth , v.y);
                MeshEdge topEdge = {v,vTopRight};
                meshEdges.push_back(topEdge);
            }
            else if(x+1 == VERTICES_X && y+1 != VERTICES_Y)
            {
                Vec2f vBottomLeft = Vec2f(v.x, v.y + quadHeight);
                MeshEdge leftEdge = {v,vBottomLeft};
                meshEdges.push_back(leftEdge);
            }
            
        }
    }
    printf("\nnumEdges = %lu", meshEdges.size());
    printf("\nnumQuads = %lu", meshQuads.size());
    
    vboMesh->bufferIndices( indices );
    vboMesh->bufferTexCoords2d( 0, texCoords );
}


void MeshWarpRetargetter::initMesh(unsigned int imgWidth, unsigned int imgHeight, SaliencySegmentor* segmentor)
{
    quadWidth = quadHeight = quadSize;
    m = imgHeight;
    n = imgWidth;
    int const VERTICES_Y = m / quadHeight + 1;
    int const VERTICES_X = n / quadWidth + 1;
    
    numVertices = VERTICES_X * VERTICES_Y;
    numEdges = 2 * (VERTICES_Y - 1) * ( VERTICES_X - 1 ) + (VERTICES_X + VERTICES_Y - 2);
    numQuads = ( VERTICES_X - 1 ) * ( VERTICES_Y - 1 );
    
    printf("\nm = %d", m);
    printf("\nn = %d", n);
    printf("\nVERTICES_X = %d", VERTICES_X);
    printf("\nVERTICES_Y = %d", VERTICES_Y);
    printf("\nnumVertices = %d", numVertices);
    printf("\nnumEdges Expected = %d", numEdges);
    printf("\nnumQuads = %d", numQuads);
    
    quadHeight = m / (VERTICES_Y - 1.f);
    quadWidth = n / (VERTICES_X - 1.f);
    
    printf("\n(quadWidth,quadHeight) = (%f,%f)",quadWidth,quadHeight);
    
    gl::VboMesh::Layout layout;
    layout.setStaticIndices();
    layout.setDynamicPositions();
    layout.setStaticTexCoords2d();
    vboMesh = ci::gl::VboMesh::create( numVertices, numQuads * 4, layout, GL_QUADS );
    
    // buffer our static data - the texcoords and the indices
    std::vector<uint32_t> indices;
    std::vector<Vec2f> texCoords;
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    // buffer edge and quad indices for retargetting
    meshEdges.clear();
    meshQuads.clear();
    
    for( int x = 0; x < VERTICES_X; ++x )
    {
        for( int y = 0; y < VERTICES_Y; ++y )
        {
            // set Texture coordinates
            Vec2f v = Vec2f( x / (VERTICES_X-1.f), y / (VERTICES_Y-1.f) );
            // the texture coordinates are mapped to [0,1] x [0,1]
            texCoords.push_back( v );
            v.x *= n-1;
            v.y *= m-1;
            // the vertex coordinates are mapped to [0,n-1] x [0,m-1]
            iter.setPosition( v.x, v.y, 0.0f );
            ++iter;
            
            // create a quad for each vertex, except for along the bottom and right edges
            if( ( x + 1 < VERTICES_X ) && ( y + 1 < VERTICES_Y ) )
            {
                int topLeft = (x+0) * VERTICES_Y + (y+0);
                int topRight = (x+1) * VERTICES_Y + (y+0);
                int bottomRight = (x+1) * VERTICES_Y + (y+1);
                int bottomLeft = (x+0) * VERTICES_Y + (y+1);
                indices.push_back(topLeft);
                indices.push_back(topRight);
                indices.push_back(bottomRight);
                indices.push_back(bottomLeft);
                
                Vec2f vTopRight = Vec2f(v.x+quadWidth , v.y);
                Vec2f vBottomRight = Vec2f(v.x+quadWidth , v.y+quadHeight);
                Vec2f vBottomLeft = Vec2f(v.x , v.y+quadHeight);
                
                MeshEdge topEdge = {v,vTopRight};
                MeshEdge leftEdge = {v,vBottomLeft};
                meshEdges.push_back(topEdge);
                meshEdges.push_back(leftEdge);
                
                MeshQuad quad = {v,vTopRight,vBottomRight,vBottomLeft};
                meshQuads.push_back(quad);
            }
            else if(y+1 == VERTICES_Y && x+1 != VERTICES_X)
            {
                Vec2f vTopRight = Vec2f(v.x+quadWidth , v.y);
                MeshEdge topEdge = {v,vTopRight};
                meshEdges.push_back(topEdge);
            }
            else if(x+1 == VERTICES_X && y+1 != VERTICES_Y)
            {
                Vec2f vBottomLeft = Vec2f(v.x, v.y + quadHeight);
                MeshEdge leftEdge = {v,vBottomLeft};
                meshEdges.push_back(leftEdge);
            }
            
        }
    }
    printf("\nnumEdges = %lu", meshEdges.size());
    printf("\nnumQuads = %lu", meshQuads.size());
    
    universe* u = segmentor->getUniverse();
    currentPatchMap.clear();
    currentPatchMap = segmentor->getPatchMap();
    
    int edgeIndex = 0;
    for(std::vector<MeshEdge>::iterator edgeIter = meshEdges.begin(); edgeIter != meshEdges.end(); edgeIter++, edgeIndex++)
    {
        Vec2f midPoint = 0.5f * (edgeIter->a + edgeIter->b);
        int edgeX = round(midPoint.x);
        int edgeY = round(midPoint.y);
        //printf("\n(x,y) = (%d,%d)",edgeX,edgeY);
        int patchID = u->find(n * edgeY + edgeX);
        currentPatchMap[patchID].edges.push_back(edgeIndex);
    }
    vboMesh->bufferIndices( indices );
    vboMesh->bufferTexCoords2d( 0, texCoords );
}





/*
 
 // from http://www.paulinternet.nl/?page=bicubic : fast catmull-rom calculation
 Vec2f MeshWarpRetargetter::cubicInterpolate( const std::vector<Vec2f> &knots, float t ) const
 {
 assert( knots.size() >= 4 );
 
 return knots[1] + 0.5f * t*(knots[2] - knots[0] +
 t*(2.0f*knots[0] - 5.0f*knots[1] +
 4.0f*knots[2] - knots[3] +
 t*(3.0f*(knots[1] - knots[2]) +
 knots[3] - knots[0])));
 }
 
 
 Vec2f MeshWarpRetargetter::getPoint(int col, int row, int yRes) const
 {
 return vertices[(col * yRes) + row];
 }
 
 
void MeshWarpRetargetter::setNumControlX(int xRes, int yRes)
{
    // there should be a minimum of 2 control points
    //n = math<int>::max(2, n);
    
    // create a list of new points
    std::vector<Vec2f> temp(xRes * yRes);
    
    // perform spline fitting
    for(int row=0;row<yRes;++row) {
        std::vector<Vec2f> points;
        if(isLinearInterpolated) {
            // construct piece-wise linear spline
            for(int col=0;col<xRes;++col) {
                points.push_back( getPoint(col,row,yRes) );
            }
            
            BSpline2f s( points, 1, false, true );
            
            // calculate position of new control points
            float length = s.getLength(0.0f, 1.0f);
            float step = 1.0f / (xRes-1);
            for(int col=0;col<xRes;++col) {
                temp[(col * yRes) + row] = s.getPosition( s.getTime( length * col * step ) );
            }
        }
        else {
            // construct piece-wise catmull-rom spline
            for(int col=0;col<xRes;++col) {
                Vec2f p0 = getPoint(col-1, row,yRes);
                Vec2f p1 = getPoint(col, row,yRes);
                Vec2f p2 = getPoint(col+1, row,yRes);
                Vec2f p3 = getPoint(col+2, row,yRes);
                
                // control points according to an optimized Catmull-Rom implementation
                Vec2f b1 = p1 + (p2 - p0) / 6.0f;
                Vec2f b2 = p2 - (p3 - p1) / 6.0f;
                
                points.push_back(p1);
                
                if(col < (xRes-1)) {
                    points.push_back(b1);
                    points.push_back(b2);
                }
            }
            
            BSpline2f s(points, 3, false, true );
            
            // calculate position of new control points
            float length = s.getLength(0.0f, 1.0f);
            float step = 1.0f / (n-1);
            for(int col=0;col<n;++col) {
                temp[(col * yRes) + row] = s.getPosition( s.getTime( length * col * step ) );
            }
        }
    }
    
    // copy new control points
    vertices = temp;
}

void MeshWarpRetargetter::setNumControlY(int xRes, int yRes)
{
    // there should be a minimum of 2 control points
    //n = math<int>::max(2, n);
    
    // create a list of new points
    std::vector<Vec2f> temp(xRes * yRes);
    
    // perform spline fitting
    for(int col=0;col<xRes;++col) {
        std::vector<Vec2f> points;
        if(isLinearInterpolated) {
            // construct piece-wise linear spline
            for(int row=0;row<yRes;++row)
                points.push_back( getPoint(col, row, yRes) );
            
            BSpline2f s( points, 1, false, true );
            
            // calculate position of new control points
            float length = s.getLength(0.0f, 1.0f);
            float step = 1.0f / (yRes-1);
            for(int row=0;row<yRes;++row) {
                temp[(col * n) + row] = s.getPosition( s.getTime( length * row * step ) );
            }
        }
        else {
            // construct piece-wise catmull-rom spline
            for(int row=0;row<yRes;++row) {
                Vec2f p0 = getPoint(col, row-1,yRes);
                Vec2f p1 = getPoint(col, row,yRes);
                Vec2f p2 = getPoint(col, row+1,yRes);
                Vec2f p3 = getPoint(col, row+2,yRes);
                
                // control points according to an optimized Catmull-Rom implementation
                Vec2f b1 = p1 + (p2 - p0) / 6.0f;
                Vec2f b2 = p2 - (p3 - p1) / 6.0f;
                
                points.push_back(p1);
                
                if(row < (yRes-1)) {
                    points.push_back(b1);
                    points.push_back(b2);
                }
            }
            
            BSpline2f s( points, 3, false, true );
            
            // calculate position of new control points
            float length = s.getLength(0.0f, 1.0f);
            float step = 1.0f / (n-1);
            for(int row=0;row<n;++row) {
                temp[(col * n) + row] = s.getPosition( s.getTime( length * row * step ) );
            }
        }
    }
    
    // copy new verices
    vertices = temp;
    
}
*/



void MeshWarpRetargetter::drawMesh(ci::gl::Texture texture)
{
    if(!vboMesh) return;
    texture.enableAndBind();
    // save current texture mode, drawing color, line width and depth buffer state
    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
    
    gl::disableDepthRead();
    gl::disableDepthWrite();
    
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
    gl::draw( vboMesh );
    
    if( isDrawingWireFrame ) {
        glDisable( GL_TEXTURE_2D );
        
        // draw wireframe
        gl::color( ColorA(1, 1, 1, 0.5f) );
        gl::enableAlphaBlending();
        gl::enableWireframe();
        gl::draw( vboMesh );
        gl::disableAlphaBlending();
        gl::disableWireframe();
        
        if(isDrawingVertices) {
            // draw control points
            //for(unsigned i=0;i<vertices.size();i++)
                //drawControlPoint( vertices[i] * mWindowSize, i == mSelected );
        }
    }
    glPopAttrib();
}

void MeshWarpRetargetter::drawEdges(ci::gl::Texture texture)
{
    gl::draw(texture);
    //for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter != currentPatchMap.end()
    //glPushMatrix();
    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
    gl::disableDepthRead();
    gl::disableDepthWrite();
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glDisable( GL_TEXTURE_2D );
    gl::enableAlphaBlending();
    
    
    glBegin( GL_LINES );
    
    for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter!= currentPatchMap.end(); patchIter++){
        std::vector<int> patchEdgeIndices = patchIter->second.edges;
        float patchSaliencyScore = patchIter->second.normalScore;
        for(std::vector<int>::iterator edgeIndexIter = patchEdgeIndices.begin(); edgeIndexIter!=patchEdgeIndices.end();edgeIndexIter++)
        {
            gl::color( ColorA(1, 0, 0, patchSaliencyScore) );
            MeshEdge edge = meshEdges[*edgeIndexIter];
            ci::Vec2f a = edge.a;
            ci::Vec2f b = edge.b;
            glVertex2f(a.x,a.y);
            glVertex2f(b.x,b.y);
        }
    }
    glEnd( );
    gl::disableAlphaBlending();
    glPopAttrib();
}

void MeshWarpRetargetter::getRigidTransformationTerms(void)
{
    std::vector<float> rigidTransformationTerms;
    
    for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter!= currentPatchMap.end(); patchIter++){
        std::vector<int> patchEdgeIndices = patchIter->second.edges;
        float patchSaliencyScore = patchIter->second.normalScore;
        float patchRigidTransformationTerm;
        
        for(std::vector<int>::iterator edgeIndexIter = patchEdgeIndices.begin(); edgeIndexIter!=patchEdgeIndices.end();edgeIndexIter++)
        {
            //patchRigidTransformationTerm +=
        }
    }
}



#endif
