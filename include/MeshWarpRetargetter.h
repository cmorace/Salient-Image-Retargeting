//
//  MeshWarpRetargetter.h
//  ImageRetargeting
//
//  Created by Charles Morace on 11/18/14.
//
//

#ifndef ImageRetargeting_MeshWarpRetargetter_h
#define ImageRetargeting_MeshWarpRetargetter_h

#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseQR"
#include "cinder/Rand.h"

class MeshWarpRetargetter {
    
public:
    
    MeshWarpRetargetter();
    void initMesh(unsigned int imgWidth, unsigned int imgHeight);
    void initMesh(unsigned int imgWidth, unsigned int imgHeight, SaliencySegmentor* segmentor);
    void drawMesh(ci::gl::Texture texture);
    void drawEdges(ci::gl::Texture texture);
    void computeOptimizationMatrix(int newWidth, int newHeight);
    void resizeMesh(int newWidth, int newHeight);
    
    int quadSize = 30;
    double transformationAlpha = 0.8;
    bool isDrawingWireFrame = true;
    
private:
    struct MeshEdge{
        // store endpoint vertices
        ci::Vec2f a;
        ci::Vec2f b;
        
        // store indices in optimization vector
        int aX_Index;
        int bX_Index;
        int aY_Index;
        int bY_Index;
    };
    std::vector<MeshEdge>	meshEdges;
    
    struct MeshQuad{
        // store quad vertices
        ci::Vec2f tl;
        ci::Vec2f tr;
        ci::Vec2f br;
        ci::Vec2f bl;
        
        // store indices in optimization vector
        int tlX_Index;
        int blX_Index;
        int trX_Index;
        int brX_Index;
        int tlY_Index;
        int trY_Index;
        int blY_Index;
        int brY_Index;
        
    };
    std::vector<MeshQuad>	meshQuads;
    
    std::vector<int>	topBoundaryIndices;
    std::vector<int>	bottomBoundaryIndices;
    std::vector<int>	leftBoundaryIndices;
    std::vector<int>	rightBoundaryIndices;
    
    ci::gl::VboMeshRef vboMesh;
    
    struct MeshPatch{
        SaliencySegmentor::Patch p;
        MeshEdge c;
        std::vector<MeshEdge> patchEdges;
        std::vector<Eigen::Matrix2d> transformation;
        std::vector<Eigen::Matrix2d> linearTransformation;
    };
    std::vector<MeshPatch> meshPatches;
    
    SaliencySegmentor::PatchMap currentPatchMap;
    
    Eigen::VectorXd vertexVectorX;
    Eigen::VectorXd answerVectorB;
    
    
    
    float quadWidth, quadHeight;
    int numXVertices;
    int numYVertices;
    int numVertices;
    int numEdges;
    int numQuads;
    int mOriginal;
    int nOriginal;
    

    Eigen::Matrix2d computeTransformation(MeshEdge c, MeshEdge e);
    Eigen::Matrix2d computeLinearTransformation(Eigen::Matrix2d, int newWidth, int newHeight, int oldWidth, int oldHeight);
};


MeshWarpRetargetter::MeshWarpRetargetter()
{
    printf("\nMeshWarpRetargetter");
    vertexVectorX = Eigen::VectorXd();
}


void MeshWarpRetargetter::initMesh(unsigned int imgWidth, unsigned int imgHeight)
{
    mOriginal = imgHeight;
    nOriginal = imgWidth;
    numXVertices = nOriginal / quadSize + 1;
    numYVertices = mOriginal / quadSize + 1;
    numVertices = numXVertices * numYVertices;
    numEdges = 2 * (numYVertices - 1) * ( numXVertices - 1 ) + (numXVertices + numYVertices - 2);
    numQuads = ( numXVertices - 1 ) * ( numYVertices - 1 );
    
    quadHeight = mOriginal / (numYVertices - 1.f);
    quadWidth = nOriginal / (numXVertices - 1.f);
    
    //setup vbo and texture map
    std::vector<uint32_t> indices;
    std::vector<ci::Vec2f> texCoords;
    cinder::gl::VboMesh::Layout layout;
    layout.setStaticIndices();
    layout.setDynamicPositions();
    layout.setStaticTexCoords2d();
    vboMesh = ci::gl::VboMesh::create( numVertices, numQuads * 4, layout, GL_QUADS );
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    
    //setup mesh optimization
    vertexVectorX.resize(numVertices * 2);
    meshEdges.clear();
    meshQuads.clear();
    
    int vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            
            // texture coordinates mapped to [0,1]x[0,1]
            ci::Vec2f v = Vec2f( x / (numXVertices-1.f),
                                 y / (numYVertices-1.f) );
            
            texCoords.push_back( v );
            
            // the vertex coordinates mapped to [0,n-1]x[0,m-1]
            v.x *= nOriginal-1;
            v.y *= mOriginal-1;
            iter.setPosition(v.x, v.y, 0.0f );
            ++iter;
            
            // save in vector for initial guess for optimization solver
            vertexVectorX(vertexCounter) = v.x;
            vertexVectorX(numVertices + vertexCounter) = v.y;
            vertexCounter++;
            
            // create a quad and 2 edges for each vertex, except along the bottom and right edges
            if( ( x + 1 < numXVertices ) && ( y + 1 < numYVertices ) ) {
                int topLeft = (x+0) * numYVertices + (y+0);
                int topRight = (x+1) * numYVertices + (y+0);
                int bottomRight = (x+1) * numYVertices + (y+1);
                int bottomLeft = (x+0) * numYVertices + (y+1);
                indices.push_back(topLeft);
                indices.push_back(topRight);
                indices.push_back(bottomRight);
                indices.push_back(bottomLeft);
                
                //calculate quad vertices
                Vec2f vTopRight = Vec2f(v.x + quadWidth , v.y);
                Vec2f vBottomRight = Vec2f(v.x + quadWidth , v.y + quadHeight);
                Vec2f vBottomLeft = Vec2f(v.x , v.y + quadHeight);
                
                // save edges
                MeshEdge topEdge = {v, vTopRight,
                                    topLeft, topRight,
                                    topLeft+numVertices, topRight+numVertices};
                
                MeshEdge leftEdge = {v, vBottomLeft,
                                     topLeft, bottomLeft,
                                     topLeft+numVertices, bottomLeft+numVertices,};
                
                meshEdges.push_back(topEdge);
                meshEdges.push_back(leftEdge);
                
                //save quad
                MeshQuad quad = {v, vTopRight, vBottomRight, vBottomLeft,
                                 topLeft, topRight, bottomRight, bottomLeft,
                                 topLeft+numVertices, topRight+numVertices, bottomRight+numVertices, bottomLeft+numVertices};
                
                meshQuads.push_back(quad);
            }
            
            // bottom boundary but not bottom right corner
            else if(y+1 == numYVertices && x+1 != numXVertices)
            {
                int topLeft = (x+0) * numYVertices + (y+0);
                int topRight = (x+1) * numYVertices + (y+0);
                
                Vec2f vTopRight = Vec2f(v.x+quadWidth , v.y);
                
                MeshEdge topEdge = {v, vTopRight,
                                    topLeft, topRight,
                                    topLeft+numVertices, topRight+numVertices};
                
                meshEdges.push_back(topEdge);

            }
            // right boundary but not bottom right corner
            else if(x+1 == numXVertices && y+1 != numYVertices)
            {
                int topLeft = (x+0) * numYVertices + (y+0);
                int bottomLeft = (x+0) * numYVertices + (y+1);
                Vec2f vBottomLeft = Vec2f(v.x, v.y + quadHeight);
                MeshEdge leftEdge = {v, vBottomLeft,
                                     topLeft, bottomLeft,
                                     topLeft+numVertices, bottomLeft+numVertices,};
                meshEdges.push_back(leftEdge);
            }
        }
    }
    vboMesh->bufferIndices( indices );
    vboMesh->bufferTexCoords2d( 0, texCoords );
    
    //set boundary indices
    topBoundaryIndices.clear();
    bottomBoundaryIndices.clear();
    leftBoundaryIndices.clear();
    rightBoundaryIndices.clear();
    
    //top and bottom boundary indices in optimization answerVectorB
    //choose y coordinate in our optimized vector
    for( int x = 0; x < numXVertices; ++x ) {
        topBoundaryIndices.push_back(x * numYVertices + numVertices);
        bottomBoundaryIndices.push_back((x+1) * numYVertices -1 + numVertices);
    }
    //left and right boundary indices in optimization answerVectorB
    for( int y = 0; y < numYVertices; ++y ) {
        leftBoundaryIndices.push_back(y);
        rightBoundaryIndices.push_back(numVertices - numYVertices + y);
    }
}


void MeshWarpRetargetter::initMesh(unsigned int imgWidth, unsigned int imgHeight, SaliencySegmentor* segmentor)
{
    initMesh(imgWidth, imgHeight);
    universe* u = segmentor->getUniverse();
    
    currentPatchMap.clear(); // todo:: should also clean up currentPatchMap contents here
    
    currentPatchMap = segmentor->getPatchMap();
    
    // for every edge find its patch (we use middle pixel)
    int edgeIndex = 0;
    for(std::vector<MeshEdge>::iterator edgeIter = meshEdges.begin(); edgeIter != meshEdges.end(); edgeIter++, edgeIndex++)
    {
        Vec2f midPoint = 0.5f * (edgeIter->a + edgeIter->b);
        int edgeX = round(midPoint.x);
        int edgeY = round(midPoint.y);
        //printf("\n(x,y) = (%d,%d)",edgeX,edgeY);
        int patchID = u->find(nOriginal * edgeY + edgeX);
        currentPatchMap[patchID].edges.push_back(edgeIndex);
    }
    
    meshPatches.clear();  // todo:: should also clean up meshPatches contents here
    
    for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter!= currentPatchMap.end(); patchIter++)
    {
        SaliencySegmentor::Patch p = patchIter->second;
        MeshEdge c;
        std::vector<MeshEdge> patchEdges;
        std::vector<Eigen::Matrix2d> transformations;
        std::vector<Eigen::Matrix2d> linearTransformations;
        if(p.edges.size() > 0)
        {
            c = meshEdges[p.edges[0]];
            for (std::vector<int>::iterator edgeIndexIter = p.edges.begin(); edgeIndexIter != p.edges.end(); edgeIndexIter++) {
                MeshEdge patchEdge = meshEdges[*edgeIndexIter];
                patchEdges.push_back(patchEdge);
                Eigen::Matrix2d L = computeTransformation(c,patchEdge);
                transformations.push_back(L);
               linearTransformations.push_back(L); //just reserve space, will replace linearTransformation when resizing
            }
        }
        MeshPatch meshPatch = {p,c,patchEdges,transformations,linearTransformations};
        meshPatches.push_back(meshPatch);
    }
}

Eigen::Matrix2d MeshWarpRetargetter::computeTransformation(MeshEdge c, MeshEdge e)
{
    double cx = c.a.x - c.b.x;
    double cy = c.a.y - c.b.y;
    double ex = e.a.x - e.b.x;
    double ey = e.a.y - e.b.y;
    double denom = cx*cx + cy*cy;
    double s = (cx*ex + cy*ey) / denom;
    double r = (cy*ex - cx*ey) / denom;
    Eigen::Matrix2d T;
    
    T << s , r,
        -r , s;
    //std::cout << "T = " << std::endl << T << std::endl;
    return T;
}

Eigen::Matrix2d MeshWarpRetargetter::computeLinearTransformation(Eigen::Matrix2d T, int newWidth, int newHeight, int oldWidth, int oldHeight)
{
    Eigen::Matrix2d L;
    
    L << (1.0 * newHeight) / oldHeight , 0,
          0 , (1.0 * newWidth) / oldWidth;
    //std::cout << "L*T = " << std::endl << L*T << std::endl;
    return L*T;
}


void MeshWarpRetargetter::computeOptimizationMatrix(int newWidth, int newHeight)
{
    //compute the edge transformations for every patch
    std::vector<Eigen::Triplet<double>> matrixA_Entries;
    
    /*
    for (std::vector<MeshPatch>::iterator iter = meshPatches.begin(); iter != meshPatches.end(); iter++) {
        MeshPatch p = *iter;
        
        if(p.transformation.size() > 0)
        {
            p.linearTransformation.clear();
            printf("\np size before = %lu", p.linearTransformation.size());
            int matCounter = 0;
            for (std::vector<Eigen::Matrix2d>::iterator matIter = p.transformation.begin(); matIter != p.transformation.end(); matIter++,matCounter++) {
                Eigen::Matrix2d T = *matIter;
                Eigen::Matrix2d LT = computeLinearTransformation(T,newWidth,newHeight,nOriginal,mOriginal);
                double s = LT(0,0);
                double r = LT(0,1);
                //printf("\n(s,r) = (%f,%f)",s,r);
                //p.linearTransformation.push_back(LT);
            }
            //printf("\nlinearTransformations.size() = %lu",p.linearTransformation.size());
        }
    }
     */
    
    
    /////TEST
    
    Eigen::VectorXd x2(2*numVertices);
    
    int vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            x2(vertexCounter) = vertexVectorX(vertexCounter);
            x2(numVertices + vertexCounter) = vertexVectorX(numVertices + vertexCounter);
            vertexCounter++;
        }
    }
    
    printf("\n compute transformation terms");
    //compute transformation terms
    
    for (std::vector<MeshPatch>::iterator iter = meshPatches.begin(); iter != meshPatches.end(); iter++) {
        MeshPatch p = *iter;
        
        if(p.patchEdges.size() > 0)
        {
            MeshEdge c = p.c;
            int edgeCounter = 0;
            for (std::vector<Eigen::Matrix2d>::iterator iter = p.transformation.begin(); iter != p.transformation.end(); iter++,edgeCounter++) {
                
                Eigen::Matrix2d T = *iter;
                Eigen::Matrix2d LT = computeLinearTransformation(T,newWidth,newHeight,nOriginal,mOriginal);
                
                // transformation equations
                // (we could compute this first transformation during preprocessing)
                // how to to account for identical vertices among edges
                /*
                double s = T(0,0);
                double r = T(0,1);
                MeshEdge edgeI = p.patchEdges[edgeCounter];
                double cx = x2(c.aX_Index) - x2(c.bX_Index);
                double cy = x2(c.aY_Index) - x2(c.bY_Index);
                
                x2(edgeI.aX_Index) = x2(edgeI.bX_Index) + s*cx + r*cy;
                x2(edgeI.aY_Index) = x2(edgeI.bY_Index) - r*cx + s*cy;
                
                x2(edgeI.bX_Index) = x2(edgeI.aX_Index) - s*cx - r*cy;
                x2(edgeI.bY_Index) = x2(edgeI.aY_Index) + r*cx - s*cy;
                */
                
                // linear transformation equations
                
                double lt00 = T(0,0);
                double lt01 = T(0,1);
                double lt10 = T(1,0);
                double lt11 = T(1,1);
                
                double xScale = 1.0*newWidth/nOriginal;
                double yScale = 1.0*newHeight/mOriginal;
                
                MeshEdge edgeI = p.patchEdges[edgeCounter];
                double cx = vertexVectorX(c.aX_Index) - vertexVectorX(c.bX_Index);
                double cy = vertexVectorX(c.aY_Index) - vertexVectorX(c.bY_Index);
                
                //test linear scaling
                x2(edgeI.aX_Index) = xScale*(vertexVectorX(edgeI.bX_Index) + lt00*cx + lt01*cy);
                x2(edgeI.aY_Index) = yScale*(vertexVectorX(edgeI.bY_Index) + lt10*cx + lt11*cy);
                
                x2(edgeI.bX_Index) = xScale*(vertexVectorX(edgeI.aX_Index) - lt00*cx - lt01*cy);
                x2(edgeI.bY_Index) = yScale*(vertexVectorX(edgeI.aY_Index) - lt10*cx - lt11*cy);
                
            }
        }
    }
    
    /*
    printf("\n compute grid orientation terms");
    for(std::vector<MeshQuad>::iterator quadIter = meshQuads.begin(); quadIter != meshQuads.end(); quadIter++)
    {
        MeshQuad quad = *quadIter;
        
        x2(quad.tlX_Index) = x2(quad.blX_Index);
        x2(quad.trX_Index) = x2(quad.brX_Index);
        x2(quad.tlY_Index) = x2(quad.trY_Index);
        x2(quad.brY_Index) = x2(quad.blY_Index);
        
    }
    
    
    printf("\n compute boundary condition terms");
    std::vector<int>::iterator botIter = bottomBoundaryIndices.begin();
    for (std::vector<int>::iterator topIter = topBoundaryIndices.begin(); topIter!=topBoundaryIndices.end(); topIter++,botIter++)
    {
        x2(*topIter) = 0;
        x2(*botIter) = newHeight;
    }
    
    std::vector<int>::iterator leftIter = leftBoundaryIndices.begin();
    for (std::vector<int>::iterator rightIter = rightBoundaryIndices.begin(); rightIter!=rightBoundaryIndices.end(); rightIter++,leftIter++)
    {
        x2(*leftIter) = 0;
        x2(*rightIter) = newWidth;
    }
    */
    
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            float vX = x2(vertexCounter);
            float vY = x2((numVertices + vertexCounter));
            iter.setPosition(vX, vY, 0.0f );
            ++iter;
            vertexCounter++;
        }
    }
    
    /*
    Eigen::SparseMatrix<double> B(66790,2*numVertices);
    
    // D_st Matrix
    
    for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter!= currentPatchMap.end(); patchIter++)
    {
        SaliencySegmentor::Patch p = patchIter->second;
        std::vector<int> edgeIndices = p.edges;
        if(edgeIndices.size() > 0)
        {
            int patchEdgeIndex = 0;
            double patchSaliency = p.normalScore * transformationAlpha;
            std::vector<int>::iterator edgeIter = edgeIndices.begin();
            MeshEdge edgeC = meshEdges[*edgeIter];
            edgeIter++;
            for(; edgeIter!= edgeIndices.end(); edgeIter++)
            {
                MeshEdge edgeI = meshEdges[*edgeIter];
                
                Eigen::Matrix2d T = patchTransformations[patchIndex].transformation[patchEdgeIndex];
                
                double s = T(0,0)*patchSaliency;
                double r = T(0,1)*patchSaliency;
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.aX_Index,patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.bX_Index,-patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aX_Index,-s));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aX_Index,-r));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bX_Index,s));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bX_Index,r));
                printf("\n\ntest\n");
                B.insert(rowIndex,edgeI.aX_Index) = patchSaliency;
                B.insert(rowIndex,edgeI.bX_Index) = -patchSaliency;
                B.insert(rowIndex,edgeC.aX_Index) = -s;
                B.insert(rowIndex,edgeC.aX_Index) = -r;
                B.insert(rowIndex,edgeC.bX_Index) = s;
                B.insert(rowIndex,edgeC.bX_Index) = r;
                
                rowIndex++;
                
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.aY_Index,patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.bY_Index,-patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aX_Index,r));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bX_Index,-r));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aY_Index,-s));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bY_Index,s));
                
                B.insert(rowIndex,edgeI.aY_Index) = patchSaliency;
                B.insert(rowIndex,edgeI.bY_Index) = -patchSaliency;
                B.insert(rowIndex,edgeC.aX_Index) = r;
                B.insert(rowIndex,edgeC.bX_Index) = -r;
                B.insert(rowIndex,edgeC.aY_Index) = -s;
                B.insert(rowIndex,edgeC.bY_Index) = s;
                
                rowIndex++;
                patchEdgeIndex++;
            }
        }
        patchIndex++;
    }
    
    patchIndex = 0;
    
    // D_lt Matrix
    for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter!= currentPatchMap.end(); patchIter++)
    {
        SaliencySegmentor::Patch p = patchIter->second;
        std::vector<int> edgeIndices = p.edges;
        if(edgeIndices.size() > 0)
        {
            int patchEdgeIndex = 0;
            double patchSaliency = p.normalScore * (1-transformationAlpha);
            std::vector<int>::iterator edgeIter = edgeIndices.begin();
            MeshEdge edgeC = meshEdges[*edgeIter];
            edgeIter++;
            for(; edgeIter!= edgeIndices.end(); edgeIter++)
            {
                MeshEdge edgeI = meshEdges[*edgeIter];
                
                Eigen::Matrix2d T = patchTransformations[patchIndex].linearTransformation[patchEdgeIndex];
                
                double s = T(0,0)*patchSaliency;
                double r = T(0,1)*patchSaliency;
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.aX_Index,patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.bX_Index,-patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aX_Index,-s));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aX_Index,-r));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bX_Index,s));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bX_Index,r));
                
                B.insert(rowIndex,edgeI.aX_Index) = patchSaliency;
                B.insert(rowIndex,edgeI.bX_Index) = -patchSaliency;
                B.insert(rowIndex,edgeC.aX_Index) = -s;
                B.insert(rowIndex,edgeC.aX_Index) = -r;
                B.insert(rowIndex,edgeC.bX_Index) = s;
                B.insert(rowIndex,edgeC.bX_Index) = r;
                
                rowIndex++;
                
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.aY_Index,patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeI.bY_Index,-patchSaliency));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aX_Index,r));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bX_Index,-r));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.aY_Index,-s));
                matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,edgeC.bY_Index,s));
                
                B.insert(rowIndex,edgeI.aY_Index) = patchSaliency;
                B.insert(rowIndex,edgeI.bY_Index) = -patchSaliency;
                B.insert(rowIndex,edgeC.aX_Index) = r;
                B.insert(rowIndex,edgeC.bX_Index) = -r;
                B.insert(rowIndex,edgeC.aY_Index) = -s;
                B.insert(rowIndex,edgeC.bY_Index) = s;
                
                rowIndex++;
                patchEdgeIndex++;
            }
        }
        patchIndex++;
    }
    
    //D_or Matrix
    for (std::vector<MeshQuad>::iterator quadIter = meshQuads.begin(); quadIter!=meshQuads.end(); quadIter++) {
        MeshQuad currentQuad = *quadIter;
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.tlY_Index,1));
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.trY_Index,-1));
        
        B.insert(rowIndex,currentQuad.tlY_Index) = 1;
        B.insert(rowIndex,currentQuad.trY_Index) = -1;
        
        rowIndex++;
        
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.blY_Index,1));
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.brY_Index,-1));
        
        B.insert(rowIndex,currentQuad.blY_Index) = 1;
        B.insert(rowIndex,currentQuad.brY_Index) = -1;
        
        rowIndex++;
        
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.tlX_Index,1));
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.blX_Index,-1));
        
        B.insert(rowIndex,currentQuad.tlX_Index) = 1;
        B.insert(rowIndex,currentQuad.blX_Index) = -1;
        
        rowIndex++;
        
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.trX_Index,1));
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,currentQuad.brX_Index,-1));
        
        B.insert(rowIndex,currentQuad.trX_Index) = 1;
        B.insert(rowIndex,currentQuad.brX_Index) = -1;
        
        rowIndex++;
    }
    
    for (std::vector<int>::iterator boundaryIter = topBoundaryIndices.begin(); boundaryIter!=topBoundaryIndices.end(); boundaryIter++) {
        int topVertexY_Index = *boundaryIter;
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,topVertexY_Index,1));
        
        B.insert(rowIndex,topVertexY_Index) = 1;
        
        rowIndex++;
    }
    int bottomIndexStart = rowIndex;
    for (std::vector<int>::iterator boundaryIter = bottomBoundaryIndices.begin(); boundaryIter!=bottomBoundaryIndices.end(); boundaryIter++) {
        int bottomVertexY_Index = *boundaryIter;
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,bottomVertexY_Index,numVertices));
        
        B.insert(rowIndex,bottomVertexY_Index) = numVertices;
        
        rowIndex++;
    }
    int bottomIndexEnd = rowIndex;
    
    for (std::vector<int>::iterator boundaryIter = leftBoundaryIndices.begin(); boundaryIter!=leftBoundaryIndices.end(); boundaryIter++) {
        int leftVertexX_Index = *boundaryIter;
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,leftVertexX_Index,1));
        
        B.insert(rowIndex,leftVertexX_Index) = 1;

        rowIndex++;
    }
    int rightIndexStart = rowIndex;
    for (std::vector<int>::iterator boundaryIter = rightBoundaryIndices.begin(); boundaryIter!=rightBoundaryIndices.end(); boundaryIter++) {
        int rightVertexX_Index = *boundaryIter;
        matrixA_Entries.push_back(Eigen::Triplet<double>(rowIndex,rightVertexX_Index,numVertices));
        
        B.insert(rowIndex,rightVertexX_Index) = numVertices;
        
        rowIndex++;
        //
    }
    int rightIndexEnd = rowIndex;
    
    printf("\nnbRows = %d",rowIndex);
    Eigen::SparseMatrix<double> A(rowIndex,2*numVertices);
    Eigen::VectorXd x2(2*numVertices);
    Eigen::VectorXd b(rowIndex);
    
    //set A
    A.setFromTriplets(matrixA_Entries.begin(), matrixA_Entries.end());
    
    // set b
    for (int i=0; i<rowIndex; i++) {
        if(i >= bottomIndexStart && i <bottomIndexEnd)
        {
           b(i) = 1.0*numVertices*newHeight;
        }
        else if(i >= rightIndexStart && i <rightIndexEnd)
        {
            b(i) = 1.0*numVertices*newWidth;
        }
        else
        {
           b(i) = 0.0;
        }
    }
    
    // check A's values
    for (std::vector<Eigen::Triplet<double>>::iterator i = matrixA_Entries.begin();
         i!= matrixA_Entries.end(); i++) {
        Eigen::Triplet<double> t = *i;
        
        double mv = B.coeff(t.row(), t.col());
        double tv = t.value();
       
        printf("\n%f = %f",mv,tv);
    }
    
    */
    
    /*
    A.makeCompressed();
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    x2 = solver.solve(b);
    
    
    
    Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > cg;
    cg.compute(A.transpose()*A);
    x2 = cg.solveWithGuess(A.transpose()*b, vertexVectorX);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    
    
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    int vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            float vX = x2(vertexCounter);
            float vY = x2((numVertices + vertexCounter));
            iter.setPosition(vX, vY, 0.0f );
            ++iter;
            vertexCounter++;
        }
    }
     */
    
    
    
}


void MeshWarpRetargetter::resizeMesh(int newWidth, int newHeight)
{
    
    
    /*
    for(int i=0; i<2*numVertices; i++){
        
        vertexVectorX(i) += 0.5-Rand::randFloat(1.f);
    }
    */
    
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    int vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            float vX = vertexVectorX(vertexCounter);
            float vY = vertexVectorX((numVertices + vertexCounter));
            iter.setPosition(vX, vY, 0.0f );
            ++iter;
            vertexCounter++;
        }
    }
    // fill A and b
    
    /*
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
    cg.compute(A);
    //x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    // update b, and solve again
    x = cg.solve(b);
     */
}


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
        gl::color( ColorA(1, 1, 1, 0.5f) );
        gl::enableAlphaBlending();
        gl::enableWireframe();
        gl::draw( vboMesh );
        gl::disableAlphaBlending();
        gl::disableWireframe();
    }
    glPopAttrib();
}

void MeshWarpRetargetter::drawEdges(ci::gl::Texture texture)
{
    //gl::draw(texture);
    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
    gl::disableDepthRead();
    gl::disableDepthWrite();
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glDisable( GL_TEXTURE_2D );
    gl::enableAlphaBlending();
    
    glLineWidth(2.0f);
    glBegin( GL_LINES );
    
    
    for(SaliencySegmentor::PatchMapIterator patchIter = currentPatchMap.begin(); patchIter!= currentPatchMap.end(); patchIter++){
        std::vector<int> patchEdgeIndices = patchIter->second.edges;
        float patchSaliencyScore = patchIter->second.normalScore;
        for(std::vector<int>::iterator edgeIndexIter = patchEdgeIndices.begin(); edgeIndexIter!=patchEdgeIndices.end();edgeIndexIter++)
        {
            // gl::color( ColorA(1, 0, 0, patchSaliencyScore) );
            float red = 0.2f+patchSaliencyScore;
            float green = 0.2f+patchSaliencyScore;
            float blue = 0.2f+patchSaliencyScore;
            float alpha = 0.2f+patchSaliencyScore;
            gl::color( ColorA(red,green,blue,alpha) );
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



#endif
