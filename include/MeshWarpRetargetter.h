//
//  MeshWarpRetargetter.h
//  ImageRetargeting
//
//  Created by Charles Morace on 11/18/14.
//
//

#ifndef ImageRetargeting_MeshWarpRetargetter_h
#define ImageRetargeting_MeshWarpRetargetter_h


class MeshWarpRetargetter {
    
public:
    MeshWarpRetargetter();
    void initMesh(unsigned int imgWidth, unsigned int imgHeight);
    void draw(ci::gl::Texture texture);
    void resize(void);
    void setMesh();
    
    bool isDrawingWireFrame = true;
    bool isDrawingVertices = false;
    
    bool isLinearInterpolated = true;
    int quadSize = 50; // this is desired quad width set by user
    float quadWidth, quadHeight;  // this is actual quad width and height
    
private:
    Vec2f getPoint(int col, int row, int yRes) const;
    Vec2f cubicInterpolate( const std::vector<Vec2f> &knots, float t ) const;
    void setNumControlX(int xRes, int yRes);
    void setNumControlY(int xRes, int yRes);
    
    ci::gl::VboMeshRef vboMesh;
    unsigned int numVertices;
    unsigned int numEdges;
    unsigned int numQuads;
    unsigned int numPatches;
    unsigned int m;
    unsigned int mPrime;
    unsigned int n;
    unsigned int nPrime;
    
    std::vector<ci::Vec2f>	vertices;
    
    struct Edge{
        unsigned int a,b; // save indices in vbo
    };
    struct Quad{
        unsigned int t,b,l,r; // save indices in vbo
    };
    
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
    numEdges = VERTICES_Y * ( VERTICES_X - 1 ) + VERTICES_X * ( VERTICES_Y - 1 );
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
    
    gl::VboMesh::Layout layout;
    layout.setStaticIndices();
    layout.setDynamicPositions();
    layout.setStaticTexCoords2d();
    vboMesh = ci::gl::VboMesh::create( numVertices, numQuads * 4, layout, GL_QUADS );
    
    // buffer our static data - the texcoords and the indices
    std::vector<uint32_t> indices;
    std::vector<Vec2f> texCoords;
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    for( int x = 0; x < VERTICES_X; ++x ) {
        for( int y = 0; y < VERTICES_Y; ++y ) {
            // create a quad for each vertex, except for along the bottom and right edges
            if( ( x + 1 < VERTICES_X ) && ( y + 1 < VERTICES_Y ) ) {
                indices.push_back( (x+0) * VERTICES_Y + (y+0) );
                indices.push_back( (x+1) * VERTICES_Y + (y+0) );
                indices.push_back( (x+1) * VERTICES_Y + (y+1) );
                indices.push_back( (x+0) * VERTICES_Y + (y+1) );
            }
            // set Texture coordinates
            Vec2f v = Vec2f( x / (VERTICES_X-1.f), y / (VERTICES_Y-1.f) );
            // the texture coordinates are mapped to [0,1.0]
            texCoords.push_back( v );
            // set Vertex Position
            iter.setPosition( n * v.x, m * v.y, 0.0f );
            ++iter;
        }
    }
    vboMesh->bufferIndices( indices );
    vboMesh->bufferTexCoords2d( 0, texCoords );
    
    /*
    vertices.clear();
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    for(int x=0;x<VERTICES_X;x++) {
        for(int y=0;y<VERTICES_Y;y++) {
            Vec2f v = Vec2f(n * x / float(VERTICES_X-1), m * y / float(VERTICES_Y-1));
            vertices.push_back( v );
            iter.setPosition( v.x, v.y, 0.0f );
            ++iter;
        }
    }
     */
    
    //setNumControlX(VERTICES_X,VERTICES_Y);
    printf("\ngot control x");
    
   // setNumControlY(VERTICES_X,VERTICES_Y);
    printf("\ngot control y");
    
    
    
    
    /*
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    Vec2f			p;
    float			u, v;
    int				col, row;
    
    std::vector<Vec2f>	cols, rows;
    
    for(int x=0; x<VERTICES_X; ++x) {
        for(int y=0; y<VERTICES_Y; ++y) {
            // transform coordinates to [0..numControls]
            u = x * (VERTICES_X - 1) / (float)(n - 1);
            v = y * (VERTICES_Y - 1) / (float)(m - 1);
            
            // determine col and row
            col = (int)(u);
            row = (int)(v);
            
            // normalize coordinates to [0..1]
            u -= col;
            v -= row;
            
            
            
            if(isLinearInterpolated) {
                Vec2f p1 = (1.0f - u) * getPoint(col, row,VERTICES_Y) + u * getPoint(col+1, row,VERTICES_Y);
                Vec2f p2 = (1.0f - u) * getPoint(col, row+1,VERTICES_Y) + u * getPoint(col+1, row+1,VERTICES_Y);
                p = ((1.0f - v) * p1 + v * p2) * Vec2f(n,m);
                
            }
            else {
                // perform bicubic interpolation
                rows.clear();
                for(int i=-1;i<3;++i) {
                    cols.clear();
                    for(int j=-1;j<3;++j) {
                        cols.push_back( getPoint(col + i, row + j,VERTICES_Y) );
                    }
                    rows.push_back( cubicInterpolate( cols, v ) );
                }
                p = cubicInterpolate( rows, u ) * Vec2f(n,m);
            }
            
            iter.setPosition( p.x, p.y, 0.0f );
            ++iter;
        }
    }
     */
}

Vec2f MeshWarpRetargetter::getPoint(int col, int row, int yRes) const
{
    return vertices[(col * yRes) + row];
}

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


void MeshWarpRetargetter::draw(ci::gl::Texture texture)
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



void MeshWarpRetargetter::resize()
{
    //meshWarper->resize();
}

void MeshWarpRetargetter::setMesh()
{
    
}



#endif
