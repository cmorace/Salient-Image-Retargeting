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
    void resizeMeshRect(int newWidth, int newHeight);
    void resizeMeshEllipse(int newWidth, int newHeight);
    
    
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
        int trX_Index;
        int brX_Index;
        int blX_Index;
        int tlY_Index;
        int trY_Index;
        int brY_Index;
        int blY_Index;
        
    };
    std::vector<MeshQuad>	meshQuads;
    
    std::vector<int>	topBoundaryIndices;
    std::vector<int>	bottomBoundaryIndices;
    std::vector<int>	leftBoundaryIndices;
    std::vector<int>	rightBoundaryIndices;
    
    std::vector<int> circleBoundaryX_Indices;
    std::vector<int> circleBoundaryY_Indices;
    
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
    int computeTransformationTerms(std::vector< Eigen::Triplet<double>> &terms, double saliencyWeight, int rowIndex);
    int computeLinearTransformationTerms(std::vector< Eigen::Triplet<double> > &terms, double saliencyWeight, int rowIndex, int newWidth, int newHeight);
    int computeGridOrientationTerms(std::vector< Eigen::Triplet<double>> &terms, int rowIndex);
    int computeBoundaryConditionTerms(std::vector< Eigen::Triplet<double>> &terms, Eigen::VectorXd &b, double weight, int rowIndex, int newWidth, int newHeight);
    int computeCircleBoundaryConditionTerms(std::vector< Eigen::Triplet<double>> &terms, Eigen::VectorXd &b, double weight, int rowIndex, int newWidth, int newHeight);
    


    void testEnergyTerms(int newWidth, int newHeight);
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
    numEdges = (numXVertices - 1) * numYVertices + (numYVertices - 1) * numXVertices;
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
    
    
    circleBoundaryX_Indices.clear();
    circleBoundaryY_Indices.clear();
    
    for( int y = 0; y < numYVertices; ++y ) { //2*y
        circleBoundaryX_Indices.push_back(y);
        circleBoundaryY_Indices.push_back(y + numVertices);
    }
    for( int x = 1; x < numXVertices; ++x ) { //2x-2)
        circleBoundaryX_Indices.push_back((x+1) * numYVertices -1);
        circleBoundaryY_Indices.push_back((x+1) * numYVertices -1 + numVertices);
    }
    for( int y = 1; y < numYVertices; ++y ) { //2y-2)
        circleBoundaryX_Indices.push_back(numVertices - y - 1);
        circleBoundaryY_Indices.push_back(numVertices -  y - 1 + numVertices);
    }
    for( int x = numXVertices-2; x > 0; --x ) { //2x - 4)
        circleBoundaryX_Indices.push_back(x * numYVertices);
        circleBoundaryY_Indices.push_back(x * numYVertices + numVertices);
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

void MeshWarpRetargetter::resizeMeshRect(int newWidth, int newHeight)
{
    int rows = 2*numEdges                                                           // transformation terms
    + 2*numEdges                                                           // linear transformation terms
    + 2*(numXVertices-2)*(numYVertices-2)+4*(numXVertices+numYVertices-3)  // grid orientation terms
    + 2*(numXVertices+numYVertices);                                       // boundary condition terms
    
    printf("\ncalculated rows = %d",rows);
    
    Eigen::SparseMatrix<double> A(rows,2*numVertices);
    A.setZero();
    
    Eigen::VectorXd b(rows);
    b.setZero();
    
    int rowIndex = 0;
    double w1 = 0.8; //saliency weight
    double w2 = numVertices; //boundary weight
    
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(rows);
    rowIndex = computeTransformationTerms(tripletList, w1, rowIndex);
    rowIndex = computeLinearTransformationTerms(tripletList, w1, rowIndex, newWidth, newHeight);
    rowIndex = computeGridOrientationTerms(tripletList, rowIndex);
    rowIndex = computeBoundaryConditionTerms(tripletList, b, w2, rowIndex, newWidth, newHeight);
    
    printf("\ntripletList.size = %lu",tripletList.size());
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    printf("\ntotal rows = %d\n",rowIndex);
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > cg;
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> ATA = AT*A;
    Eigen::VectorXd ATb = AT*b;
    cg.compute(ATA);
    
    Eigen::VectorXd residual;
    
    
    Eigen::VectorXd x2(2*numVertices);
    int vertexCounter = 0;
    
    
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            x2(vertexCounter) = vertexVectorX(vertexCounter);
            x2(numVertices + vertexCounter) = vertexVectorX(numVertices + vertexCounter);
            vertexCounter++;
        }
    }
    
    
    x2 = cg.solveWithGuess(ATb,x2);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    
    
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            float vX = x2(vertexCounter);
            float vY = x2(numVertices + vertexCounter);
            iter.setPosition(vX, vY, 0.0f );
            ++iter;
            vertexCounter++;
        }
    }
}


void MeshWarpRetargetter::resizeMeshEllipse(int newWidth, int newHeight)
{
    int rows = 2*numEdges                                                           // transformation terms
             + 2*numEdges                                                           // linear transformation terms
             + 2*(numXVertices-2)*(numYVertices-2)+4*(numXVertices+numYVertices-3)  // grid orientation terms
             + 2*(numXVertices+numYVertices);                                       // boundary condition terms
    
    //cicle boundary testing
    rows += 2*(numXVertices+numYVertices)-8;
    printf("\ncalculated rows = %d",rows);
    
    Eigen::SparseMatrix<double> A(rows,2*numVertices);
    A.setZero();
    
    Eigen::VectorXd b(rows);
    b.setZero();
    
    
    
    int rowIndex = 0;
    double w1 = 0.8; //saliency weight
    double w2 = numVertices; //boundary weight
    
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(rows);
    rowIndex = computeTransformationTerms(tripletList, w1, rowIndex);
    rowIndex = computeLinearTransformationTerms(tripletList, w1, rowIndex, newWidth, newHeight);
    rowIndex = computeGridOrientationTerms(tripletList, rowIndex);
    //rowIndex = computeBoundaryConditionTerms(tripletList, b, w2, rowIndex, newWidth, newHeight);
    rowIndex = computeCircleBoundaryConditionTerms(tripletList, b, w2, rowIndex, newWidth, newHeight);
    
    printf("\ntripletList.size = %lu",tripletList.size());
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    printf("\ntotal rows = %d\n",rowIndex);
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > cg;
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> ATA = AT*A;
    Eigen::VectorXd ATb = AT*b;
    cg.compute(ATA);
    
    Eigen::VectorXd residual;
    
    
    Eigen::VectorXd x2(2*numVertices);
    int vertexCounter = 0;
    
    
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            x2(vertexCounter) = vertexVectorX(vertexCounter);
            x2(numVertices + vertexCounter) = vertexVectorX(numVertices + vertexCounter);
            vertexCounter++;
        }
    }
    
    
    x2 = cg.solveWithGuess(ATb,x2);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    
    
    gl::VboMesh::VertexIter iter = vboMesh->mapVertexBuffer();
    vertexCounter = 0;
    for( int x = 0; x < numXVertices; ++x ) {
        for( int y = 0; y < numYVertices; ++y ) {
            float vX = x2(vertexCounter);
            float vY = x2(numVertices + vertexCounter);
            iter.setPosition(vX, vY, 0.0f );
            ++iter;
            vertexCounter++;
        }
    }
}

int MeshWarpRetargetter::computeBoundaryConditionTerms(std::vector< Eigen::Triplet<double> > &terms, Eigen::VectorXd &b, double weight, int rowIndex, int newWidth, int newHeight)
{
    std::vector<int>::iterator botIter = bottomBoundaryIndices.begin();
    for (std::vector<int>::iterator topIter = topBoundaryIndices.begin(); topIter!=topBoundaryIndices.end(); topIter++,botIter++)
    {
        //we could precompute these
        //x2(*topIter) = 0;
        terms.push_back(Eigen::Triplet<double>(rowIndex, *topIter, weight));
        rowIndex++;
        
        //x2(*botIter) = newHeight;
        terms.push_back(Eigen::Triplet<double>(rowIndex, *botIter, weight));
        b(rowIndex) = weight*newHeight;
        rowIndex++;
    }
    
    std::vector<int>::iterator leftIter = leftBoundaryIndices.begin();
    for (std::vector<int>::iterator rightIter = rightBoundaryIndices.begin(); rightIter!=rightBoundaryIndices.end(); rightIter++,leftIter++)
    {
        //we could precompute these
        //x2(*leftIter) = 0;
        terms.push_back(Eigen::Triplet<double>(rowIndex, *leftIter, weight));
        rowIndex++;
        
        //x2(*rightIter) = newWidth;
        terms.push_back(Eigen::Triplet<double>(rowIndex, *rightIter, weight));
        b(rowIndex) = weight*newWidth;
        rowIndex++;
    }
    return rowIndex;
}

int MeshWarpRetargetter::computeCircleBoundaryConditionTerms(std::vector< Eigen::Triplet<double> > &terms, Eigen::VectorXd &b, double weight, int rowIndex, int newWidth, int newHeight)
{
    
    cinder::Vec2i c = Vec2i(newWidth/2,newHeight/2);
    float PI = 3.14159f;
    float theta3 = atan((1.0*c.y)/c.x);
    float theta0 = PI-theta3;
    float dty = 2*theta3/(numYVertices-1);
    float dtx = (PI - 2*theta3)/(numXVertices-1);
    int vertexCounter = 0;
    std::vector<int>::iterator xIter = circleBoundaryX_Indices.begin();
    for (std::vector<int>::iterator yIter = circleBoundaryY_Indices.begin(); yIter!=circleBoundaryY_Indices.end(); xIter++,yIter++)
    {
        if(vertexCounter < numYVertices){
            float theta = theta0 + vertexCounter * dty;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *xIter, weight));
            b(rowIndex) = weight*(c.x + c.x*cos(theta));
            rowIndex++;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *yIter, weight));
            b(rowIndex) = weight*(c.y - c.y*sin(theta));
            rowIndex++;
            vertexCounter++;
        }
        else if( vertexCounter < numYVertices + numXVertices - 2)
        {
            
            float theta = theta0 + 2*theta3 + dtx + (vertexCounter - numYVertices)*dtx;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *xIter, weight));
            b(rowIndex) = weight*(c.x + c.x*cos(theta));
            rowIndex++;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *yIter, weight));
            b(rowIndex) = weight*(c.y - c.y*sin(theta));
            rowIndex++;
            vertexCounter++;
        }
        else if( vertexCounter < 2*numYVertices + numXVertices - 3)
        {
            float theta = theta0 + PI + dty + (vertexCounter - (numYVertices + numXVertices - 2))*dty;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *xIter, weight));
            b(rowIndex) = weight*(c.x + c.x*cos(theta));
            rowIndex++;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *yIter, weight));
            b(rowIndex) = weight*(c.y - c.y*sin(theta));
            rowIndex++;
            vertexCounter++;
        }
        else{
            float theta = theta0 + 2*theta3 + PI + dtx + (vertexCounter - (2* numYVertices + numXVertices - 3))*dtx;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *xIter, weight));
            b(rowIndex) = weight*(c.x + c.x*cos(theta));
            rowIndex++;
            terms.push_back(Eigen::Triplet<double>(rowIndex, *yIter, weight));
            b(rowIndex) = weight*(c.y - c.y*sin(theta));
            rowIndex++;
            vertexCounter++;
        }
    }
    return rowIndex;
}


int MeshWarpRetargetter::computeGridOrientationTerms(std::vector< Eigen::Triplet<double>> &terms, int rowIndex){
    std::vector<MeshQuad>::iterator quadIter = meshQuads.begin();
    
    for(int x=0; x<numXVertices-1; x++)
    {
        for(int y=0; y<numYVertices-1; y++,quadIter++)
        {
            MeshQuad quad = *quadIter;
            
            if(y < numYVertices-2 && x < numXVertices-2)
            {
                //x2(quad.tlY_Index) = x2(quad.trY_Index);
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.tlY_Index, 1));
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.trY_Index, -1));
                rowIndex++;
                
                //x2(quad.tlX_Index) = x2(quad.blX_Index);
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.tlX_Index, 1));
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.blX_Index, -1));
                rowIndex++;
            }
            else    //  bottom row and left column set all 4 vertices
            {
                // x2(quad.tlY_Index) = x2(quad.trY_Index);
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.tlY_Index, 1));
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.trY_Index, -1));
                rowIndex++;
                
                // x2(quad.tlX_Index) = x2(quad.blX_Index);
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.tlX_Index, 1));
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.blX_Index, -1));
                rowIndex++;
                
                // x2(quad.trX_Index) = x2(quad.brX_Index);
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.trX_Index, 1));
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.brX_Index, -1));
                rowIndex++;
                
                // x2(quad.blY_Index) = x2(quad.brY_Index);
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.blY_Index, 1));
                terms.push_back(Eigen::Triplet<double>(rowIndex, quad.blY_Index, -1));
                rowIndex++;
            }
        }
    }
    return rowIndex;
}


int MeshWarpRetargetter::computeLinearTransformationTerms(std::vector< Eigen::Triplet<double> > &terms, double saliencyWeight, int rowIndex, int newWidth, int newHeight)
{
    for (std::vector<MeshPatch>::iterator iter = meshPatches.begin(); iter != meshPatches.end(); iter++) {
        MeshPatch p = *iter;
        
        if(p.patchEdges.size() > 0)
        {
            MeshEdge c = p.c;
            int edgeCounter = 0;
            for (std::vector<Eigen::Matrix2d>::iterator iter = p.transformation.begin(); iter != p.transformation.end(); iter++,edgeCounter++) {
                
                Eigen::Matrix2d T = *iter;
                Eigen::Matrix2d LT = computeLinearTransformation(T,newWidth,newHeight,nOriginal,mOriginal);
                
                int w = (1-saliencyWeight) * p.p.normalScore;
                // Papers equations
                double lt00 = w*LT(0,0);
                double lt01 = w*LT(0,1);
                double lt10 = w*LT(1,0);
                double lt11 = w*LT(1,1);
                
                MeshEdge edgeI = p.patchEdges[edgeCounter];
                
                if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
                {
                    if (c.bX_Index == edgeI.bX_Index
                        && c.bY_Index == edgeI.bY_Index) //c = e
                    {
                        //printf("\n c = e (same edge)");
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt00-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, w-lt00));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt01));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt01));
                        rowIndex++;
                        
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt10));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt11-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, w-lt11));
                        rowIndex++;
                    }
                    else //c.a = e.a
                    {
                        terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bX_Index, w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt00-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt00));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt01));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt01));
                        rowIndex++;
                        
                        terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bY_Index, w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt10));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt11-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt11));
                        rowIndex++;
                    }
                }
                else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index) //  c.b = e.a
                {
                    //printf("\n c.b = e.a (bottom right)");
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bX_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt00));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10-w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt01));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt01));
                    rowIndex++;
                    
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bY_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt11));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt11-w));
                    rowIndex++;
                    
                }
                else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index) //  c.a = e.b
                {
                    // printf("\n c.a = e.b (top left)");
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aX_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, w+lt00));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt00));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt01));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt01));
                    rowIndex++;
                    
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aY_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, w+lt11));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt11));
                    rowIndex++;
                    
                }
                else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index) //  c.b = e.b
                {
                    //printf("\n c.b = e.b (bottom left)");
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aX_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt00));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, w-lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt01));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt01));
                    rowIndex++;
                    
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aY_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt11));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, w-lt11));
                    rowIndex++;
                }
                else{ // edges not connected
                    //printf("\n edges not connected");
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aX_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bX_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt00));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt00));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt01));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt01));
                    rowIndex++;
                    
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aY_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bY_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -lt10));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, lt11));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -lt11));
                    rowIndex++;
                }
            }
        }
    }
    return rowIndex;
}


int MeshWarpRetargetter::computeTransformationTerms(std::vector< Eigen::Triplet<double>> &terms, double saliencyWeight, int rowIndex)
{
    
    for (std::vector<MeshPatch>::iterator iter = meshPatches.begin(); iter != meshPatches.end(); iter++) {
        MeshPatch p = *iter;
        
        if(p.patchEdges.size() > 0)
        {
            MeshEdge c = p.c;
            double w = saliencyWeight*p.p.normalScore; //saliency weight
            int edgeCounter = 0;
            for (std::vector<Eigen::Matrix2d>::iterator iter = p.transformation.begin(); iter != p.transformation.end(); iter++,edgeCounter++) {
                
                Eigen::Matrix2d T = *iter;
                // transformation equations
                // (we could compute this first transformation during preprocessing)
                
                double s = w*T(0,0);
                double r = w*T(0,1);
                MeshEdge edgeI = p.patchEdges[edgeCounter];
                
                // we need to account for edges that share vertices with the patch representative edge
                if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
                {
                    if (c.bX_Index == edgeI.bX_Index
                        && c.bY_Index == edgeI.bY_Index)
                        //  c = e
                    {
                        //x2(c.aX_Index) = (1-s)*x2(c.bX_Index) + s*x2(c.aX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, s-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, w-s));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, r));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -r));
                        rowIndex++;
                        
                        //x2(c.aY_Index) = (1-s)*x2(c.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index);
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, -r));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, r));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, s-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, w-s));
                        rowIndex++;
                    }
                    else                                     // c.a = e.a
                    {
                        //x2(c.aX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                        terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bX_Index, w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, s-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -s));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, r));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -r));
                        rowIndex++;
                        
                        //x2(c.aY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
                        terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bY_Index, w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, -r));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, r));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, s-w));
                        terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -s));
                        rowIndex++;
                    }
                }
                else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index)
                    //  c.b = e.a
                {
                    //printf("\n c.b = e.a (bottom right)");
                    //x2(c.bX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bX_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -w-s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -r));
                    rowIndex++;
                    
                    //x2(c.bY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bY_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, -r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -w-s));
                    rowIndex++;
                    
                }
                else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index)
                    //  c.a = e.b
                {
                    //printf("\n c.a = e.b (top left)");
                    //x2(edgeI.aX_Index) = (1+s)*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aX_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, w+s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -r));
                    rowIndex++;
                    
                    //x2(edgeI.aY_Index) = (1+s)*x2(c.aY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) - s*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aY_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, -r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, w+s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -s));
                    rowIndex++;
                }
                else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index)
                    //  c.b = e.b
                {
                    //printf("\n c.b = e.b (bottom left)");
                    //x2(edgeI.aX_Index) = (1-s)*x2(c.bX_Index) + s*x2(c.aX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aX_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, w-s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -r));
                    rowIndex++;
                    
                    //x2(edgeI.aY_Index) = (1-s)*x2(c.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aY_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, -r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, w-s));
                    rowIndex++;
                }
                else{                                       // edges not connected
                    //printf("\n edges not connected");
                    //x2(edgeI.aX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aX_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bX_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, -s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -r));
                    rowIndex++;
                    
                    //x2(edgeI.aY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.aY_Index, -w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, edgeI.bY_Index, w));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aX_Index, -r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bX_Index, r));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.aY_Index, s));
                    terms.push_back(Eigen::Triplet<double>(rowIndex, c.bY_Index, -s));
                    rowIndex++;
                }
            }
        }
    }
    return rowIndex;
}


void MeshWarpRetargetter::testEnergyTerms(int newWidth, int newHeight)
{
    /////                               TESTING ONLY!!!
    
    /*
     printf("\n compute transformation terms");
     
     //compute transformation terms
     for (std::vector<MeshPatch>::iterator iter = meshPatches.begin(); iter != meshPatches.end(); iter++) {
     MeshPatch p = *iter;
     
     if(p.patchEdges.size() > 0)
     {
     MeshEdge c = p.c;
     double w = w1*p.p.normalScore; //saliency weight
     int edgeCounter = 0;
     for (std::vector<Eigen::Matrix2d>::iterator iter = p.transformation.begin(); iter != p.transformation.end(); iter++,edgeCounter++) {
     
     Eigen::Matrix2d T = *iter;
     // transformation equations
     // (we could compute this first transformation during preprocessing)
     
     double s = w*T(0,0);
     double r = w*T(0,1);
     MeshEdge edgeI = p.patchEdges[edgeCounter];
     
     // we need to account for edges that share vertices with the patch representative edge
     if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
     {
     if (c.bX_Index == edgeI.bX_Index
     && c.bY_Index == edgeI.bY_Index)
     //  c = e
     {
     //x2(c.aX_Index) = (1-s)*x2(c.bX_Index) + s*x2(c.aX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
     A.insert(rowIndex, c.aX_Index) = s-w;
     A.insert(rowIndex, c.bX_Index) = w-s;
     A.insert(rowIndex, c.aY_Index) = r;
     A.insert(rowIndex, c.bY_Index) = -r;
     rowIndex++;
     
     //x2(c.aY_Index) = (1-s)*x2(c.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index);
     A.insert(rowIndex, c.aX_Index) = -r;
     A.insert(rowIndex, c.bX_Index) = r;
     A.insert(rowIndex, c.aY_Index) = s-w;
     A.insert(rowIndex, c.bY_Index) = w-s;
     rowIndex++;
     }
     else                                     // c.a = e.a
     {
     //x2(c.aX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.bX_Index) = w;
     A.insert(rowIndex, c.aX_Index) = s-w;
     A.insert(rowIndex, c.bX_Index) = -s;
     A.insert(rowIndex, c.aY_Index) = r;
     A.insert(rowIndex, c.bY_Index) = -r;
     rowIndex++;
     
     //x2(c.aY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.bY_Index) = w;
     A.insert(rowIndex, c.aX_Index) = -r;
     A.insert(rowIndex, c.bX_Index) = r;
     A.insert(rowIndex, c.aY_Index) = s-w;
     A.insert(rowIndex, c.bY_Index) = -s;
     rowIndex++;
     }
     }
     else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index)
     //  c.b = e.a
     {
     //printf("\n c.b = e.a (bottom right)");
     //x2(c.bX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.bX_Index) = w;
     A.insert(rowIndex, c.aX_Index) = s;
     A.insert(rowIndex, c.bX_Index) = -w-s;
     A.insert(rowIndex, c.aY_Index) = r;
     A.insert(rowIndex, c.bY_Index) = -r;
     rowIndex++;
     
     //x2(c.bY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.bY_Index) = w;
     A.insert(rowIndex, c.aX_Index) = -r;
     A.insert(rowIndex, c.bX_Index) = r;
     A.insert(rowIndex, c.aY_Index) = s;
     A.insert(rowIndex, c.bY_Index) = -w-s;
     rowIndex++;
     
     }
     else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index)
     //  c.a = e.b
     {
     //printf("\n c.a = e.b (top left)");
     //x2(edgeI.aX_Index) = (1+s)*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.aX_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = w+s;
     A.insert(rowIndex, c.bX_Index) = -s;
     A.insert(rowIndex, c.aY_Index) = r;
     A.insert(rowIndex, c.bY_Index) = -r;
     rowIndex++;
     
     //x2(edgeI.aY_Index) = (1+s)*x2(c.aY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) - s*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.aY_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = -r;
     A.insert(rowIndex, c.bX_Index) = r;
     A.insert(rowIndex, c.aY_Index) = w+s;
     A.insert(rowIndex, c.bY_Index) = -s;
     rowIndex++;
     }
     else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index)
     //  c.b = e.b
     {
     //printf("\n c.b = e.b (bottom left)");
     //x2(edgeI.aX_Index) = (1-s)*x2(c.bX_Index) + s*x2(c.aX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.aX_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = s;
     A.insert(rowIndex, c.bX_Index) = w-s;
     A.insert(rowIndex, c.aY_Index) = r;
     A.insert(rowIndex, c.bY_Index) = -r;
     rowIndex++;
     
     //x2(edgeI.aY_Index) = (1-s)*x2(c.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index);
     A.insert(rowIndex, edgeI.aY_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = -r;
     A.insert(rowIndex, c.bX_Index) = r;
     A.insert(rowIndex, c.aY_Index) = s;
     A.insert(rowIndex, c.bY_Index) = w-s;
     rowIndex++;
     }
     else{                                       // edges not connected
     //printf("\n edges not connected");
     //x2(edgeI.aX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.aX_Index) = -w;
     A.insert(rowIndex, edgeI.bX_Index) = w;
     A.insert(rowIndex, c.aX_Index) = s;
     A.insert(rowIndex, c.bX_Index) = -s;
     A.insert(rowIndex, c.aY_Index) = r;
     A.insert(rowIndex, c.bY_Index) = -r;
     rowIndex++;
     
     //x2(edgeI.aY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
     A.insert(rowIndex, edgeI.aY_Index) = -w;
     A.insert(rowIndex, edgeI.bY_Index) = w;
     A.insert(rowIndex, c.aX_Index) = -r;
     A.insert(rowIndex, c.bX_Index) = r;
     A.insert(rowIndex, c.aY_Index) = s;
     A.insert(rowIndex, c.bY_Index) = -s;
     rowIndex++;
     }
     }
     }
     }
     
     printf("\n nb of transformation rows = %d",rowIndex);
     
     for (std::vector<MeshPatch>::iterator iter = meshPatches.begin(); iter != meshPatches.end(); iter++) {
     MeshPatch p = *iter;
     
     if(p.patchEdges.size() > 0)
     {
     MeshEdge c = p.c;
     int edgeCounter = 0;
     for (std::vector<Eigen::Matrix2d>::iterator iter = p.transformation.begin(); iter != p.transformation.end(); iter++,edgeCounter++) {
     
     Eigen::Matrix2d T = *iter;
     Eigen::Matrix2d LT = computeLinearTransformation(T,newWidth,newHeight,nOriginal,mOriginal);
     
     int w = (1-w1) * p.p.normalScore;
     // Papers equations
     double lt00 = w*LT(0,0);
     double lt01 = w*LT(0,1);
     double lt10 = w*LT(1,0);
     double lt11 = w*LT(1,1);
     
     MeshEdge edgeI = p.patchEdges[edgeCounter];
     
     if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
     {
     if (c.bX_Index == edgeI.bX_Index
     && c.bY_Index == edgeI.bY_Index) //c = e
     {
     //printf("\n c = e (same edge)");
     A.insert(rowIndex, c.aX_Index) = lt00-w;
     A.insert(rowIndex, c.bX_Index) = w-lt00;
     A.insert(rowIndex, c.aY_Index) = lt01;
     A.insert(rowIndex, c.bY_Index) = -lt01;
     rowIndex++;
     
     A.insert(rowIndex, c.aX_Index) = lt10;
     A.insert(rowIndex, c.bX_Index) = -lt10;
     A.insert(rowIndex, c.aY_Index) = lt11-w;
     A.insert(rowIndex, c.bY_Index) = w-lt11;
     rowIndex++;
     
     }
     else //c.a = e.a
     {
     
     A.insert(rowIndex, edgeI.bX_Index) = w;
     A.insert(rowIndex, c.aX_Index) = lt00-w;
     A.insert(rowIndex, c.bX_Index) = -lt00;
     A.insert(rowIndex, c.aY_Index) = lt01;
     A.insert(rowIndex, c.bY_Index) = -lt01;
     rowIndex++;
     
     A.insert(rowIndex, edgeI.bY_Index) = w;
     A.insert(rowIndex, c.aX_Index) = lt10;
     A.insert(rowIndex, c.bX_Index) = -lt10;
     A.insert(rowIndex, c.aY_Index) = lt11-w;
     A.insert(rowIndex, c.bY_Index) = -lt11;
     rowIndex++;
     }
     }
     else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index) //  c.b = e.a
     {
     //printf("\n c.b = e.a (bottom right)");
     
     A.insert(rowIndex, edgeI.bX_Index) = w;
     A.insert(rowIndex, c.aX_Index) = lt00;
     A.insert(rowIndex, c.bX_Index) = -lt00-w;
     A.insert(rowIndex, c.aY_Index) = lt01;
     A.insert(rowIndex, c.bY_Index) = -lt01;
     rowIndex++;
     
     A.insert(rowIndex, edgeI.bY_Index) = w;
     A.insert(rowIndex, c.aX_Index) = lt10;
     A.insert(rowIndex, c.bX_Index) = -lt10;
     A.insert(rowIndex, c.aY_Index) = lt11;
     A.insert(rowIndex, c.bY_Index) = -lt11-w;
     rowIndex++;
     }
     else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index) //  c.a = e.b
     {
     // printf("\n c.a = e.b (top left)");
     
     A.insert(rowIndex, edgeI.aX_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = w+lt00;
     A.insert(rowIndex, c.bX_Index) = -lt00;
     A.insert(rowIndex, c.aY_Index) = +lt01;
     A.insert(rowIndex, c.bY_Index) = -lt01;
     rowIndex++;
     
     A.insert(rowIndex, edgeI.aY_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = lt10;
     A.insert(rowIndex, c.bX_Index) = -lt10;
     A.insert(rowIndex, c.aY_Index) = w+lt11;
     A.insert(rowIndex, c.bY_Index) = -lt11;
     rowIndex++;
     
     }
     else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index) //  c.b = e.b
     {
     //printf("\n c.b = e.b (bottom left)");
     
     A.insert(rowIndex, edgeI.aX_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = lt00;
     A.insert(rowIndex, c.bX_Index) = w-lt00;
     A.insert(rowIndex, c.aY_Index) = lt01;
     A.insert(rowIndex, c.bY_Index) = -lt01;
     rowIndex++;
     
     A.insert(rowIndex, edgeI.aY_Index) = -w;
     A.insert(rowIndex, c.aX_Index) = lt10;
     A.insert(rowIndex, c.bX_Index) = -lt10;
     A.insert(rowIndex, c.aY_Index) = lt11;
     A.insert(rowIndex, c.bY_Index) = w-lt11;
     rowIndex++;
     }
     else{ // edges not connected
     //printf("\n edges not connected");
     
     A.insert(rowIndex, edgeI.aX_Index) = -w;
     A.insert(rowIndex, edgeI.bX_Index) = w;
     A.insert(rowIndex, c.aX_Index) = lt00;
     A.insert(rowIndex, c.bX_Index) = -lt00;
     A.insert(rowIndex, c.aY_Index) = lt01;
     A.insert(rowIndex, c.bY_Index) = -lt01;
     rowIndex++;
     
     A.insert(rowIndex, edgeI.aY_Index) = -w;
     A.insert(rowIndex, edgeI.bY_Index) = w;
     A.insert(rowIndex, c.aX_Index) = lt10;
     A.insert(rowIndex, c.bX_Index) = -lt10;
     A.insert(rowIndex, c.aY_Index) = lt11;
     A.insert(rowIndex, c.bY_Index) = -lt11;
     rowIndex++;
     }
     }
     }
     }
     
     
     printf("\n compute grid orientation terms");
     std::vector<MeshQuad>::iterator quadIter = meshQuads.begin();
     
     for(int x=0; x<numXVertices-1; x++)
     {
     for(int y=0; y<numYVertices-1; y++,quadIter++)
     {
     MeshQuad quad = *quadIter;
     
     if(y < numYVertices-2 && x < numXVertices-2)
     {
     //x2(quad.tlY_Index) = x2(quad.trY_Index);
     A.insert(rowIndex, quad.tlY_Index) = 1;
     A.insert(rowIndex, quad.trY_Index) = -1;
     rowIndex++;
     
     //x2(quad.tlX_Index) = x2(quad.blX_Index);
     A.insert(rowIndex, quad.tlX_Index) = 1;
     A.insert(rowIndex, quad.blX_Index) = -1;
     rowIndex++;
     }
     else{
     // x2(quad.tlY_Index) = x2(quad.trY_Index);
     A.insert(rowIndex, quad.tlY_Index) = 1;
     A.insert(rowIndex, quad.trY_Index) = -1;
     rowIndex++;
     
     // x2(quad.tlX_Index) = x2(quad.blX_Index);
     A.insert(rowIndex, quad.tlX_Index) = 1;
     A.insert(rowIndex, quad.blX_Index) = -1;
     rowIndex++;
     
     // x2(quad.trX_Index) = x2(quad.brX_Index);
     A.insert(rowIndex, quad.trX_Index) = 1;
     A.insert(rowIndex, quad.brX_Index) = -1;
     rowIndex++;
     
     // x2(quad.blY_Index) = x2(quad.brY_Index);
     A.insert(rowIndex, quad.blY_Index) = 1;
     A.insert(rowIndex, quad.brY_Index) = -1;
     rowIndex++;
     }
     }
     }
     
     
     printf("\n compute boundary condition terms");
     std::vector<int>::iterator botIter = bottomBoundaryIndices.begin();
     for (std::vector<int>::iterator topIter = topBoundaryIndices.begin(); topIter!=topBoundaryIndices.end(); topIter++,botIter++)
     {
     //we could precompute these
     //x2(*topIter) = 0;
     A.insert(rowIndex, *topIter) = w2;
     rowIndex++;
     
     //x2(*botIter) = newHeight;
     A.insert(rowIndex, *botIter) = w2;
     b(rowIndex) = w2*newHeight;
     rowIndex++;
     }
     
     std::vector<int>::iterator leftIter = leftBoundaryIndices.begin();
     for (std::vector<int>::iterator rightIter = rightBoundaryIndices.begin(); rightIter!=rightBoundaryIndices.end(); rightIter++,leftIter++)
     {
     //we could precompute these
     //x2(*leftIter) = 0;
     A.insert(rowIndex, *leftIter) = w2;
     rowIndex++;
     
     //x2(*rightIter) = newWidth;
     A.insert(rowIndex, *rightIter) = w2;
     b(rowIndex) = w2*newWidth;
     rowIndex++;
     }
     */
    
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
                
                
                double s = T(0,0);
                double r = T(0,1);
                MeshEdge edgeI = p.patchEdges[edgeCounter];
                
                
                if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
                {
                    if (c.bX_Index == edgeI.bX_Index
                        && c.bY_Index == edgeI.bY_Index) //c = e
                    {
                        printf("\n c = e (same edge)");
                        x2(c.aX_Index) = (1-s)*x2(c.bX_Index) + s*x2(c.aX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                        x2(c.aY_Index) = (1-s)*x2(c.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index);
                        
                        x2(c.bX_Index) = (1-s)*x2(c.aX_Index) + s*x2(c.bX_Index) - r*x2(c.aY_Index) + r*x2(c.bY_Index);
                        x2(c.bY_Index) = (1-s)*x2(c.aY_Index) + r*x2(c.aX_Index) - r*x2(c.bX_Index) + s*x2(c.bY_Index);
                    }
                    else //c.a = e.a
                    {
                        printf("\n c.a = e.a (top right)");
                        x2(c.aX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                        x2(c.aY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
                        
                        x2(edgeI.bX_Index) = (1-s)*x2(c.aX_Index) + s*x2(c.bX_Index) - r*x2(c.aY_Index) + r*x2(c.bY_Index);
                        x2(edgeI.bY_Index) = (1-s)*x2(c.aY_Index) + r*x2(c.aX_Index) - r*x2(c.bX_Index) + s*x2(c.bY_Index);
                    }
                }
                else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index) //  c.b = e.a
                {
                    printf("\n c.b = e.a (bottom right)");
                    x2(c.bX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    x2(c.bY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
                    
                    x2(edgeI.bX_Index) = (1+s)*x2(c.bX_Index) - s*x2(c.aX_Index) - r*x2(c.aY_Index) + r*x2(c.bY_Index);
                    x2(edgeI.bY_Index) = (1+s)*x2(c.bY_Index) + r*x2(c.aX_Index) - r*x2(c.bX_Index) - s*x2(c.aY_Index);
                    
                }
                else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index) //  c.a = e.b
                {
                    printf("\n c.a = e.b (top left)");
                    x2(edgeI.aX_Index) = (1+s)*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    x2(edgeI.aY_Index) = (1+s)*x2(c.aY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) - s*x2(c.bY_Index);
                    
                    x2(c.aX_Index) = x2(edgeI.aX_Index) - s*x2(c.aX_Index) + s*x2(c.bX_Index) - r*x2(c.aY_Index) + r*x2(c.bY_Index);
                    x2(c.aY_Index) = x2(edgeI.aY_Index) + r*x2(c.aX_Index) - r*x2(c.bX_Index) - s*x2(c.aY_Index) + s*x2(c.bY_Index);
                }
                else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index) //  c.b = e.b
                {
                    printf("\n c.b = e.b (bottom left)");
                    x2(edgeI.aX_Index) = (1-s)*x2(c.bX_Index) + s*x2(c.aX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    x2(edgeI.aY_Index) = (1-s)*x2(c.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index);
                    
                    x2(c.bX_Index) = x2(edgeI.aX_Index) - s*x2(c.aX_Index) + s*x2(c.bX_Index) - r*x2(c.aY_Index) + r*x2(c.bY_Index);
                    x2(c.bY_Index) = x2(edgeI.aY_Index) + r*x2(c.aX_Index) - r*x2(c.bX_Index) - s*x2(c.aY_Index) + s*x2(c.bY_Index);
                }
                else
                { // edges not connected
                    //printf("\n edges not connected");
                    x2(edgeI.aX_Index) = x2(edgeI.bX_Index) + s*x2(c.aX_Index) - s*x2(c.bX_Index) + r*x2(c.aY_Index) - r*x2(c.bY_Index);
                    x2(edgeI.aY_Index) = x2(edgeI.bY_Index) - r*x2(c.aX_Index) + r*x2(c.bX_Index) + s*x2(c.aY_Index) - s*x2(c.bY_Index);
                    
                    x2(edgeI.bX_Index) = x2(edgeI.aX_Index) - s*x2(c.aX_Index) + s*x2(c.bX_Index) - r*x2(c.aY_Index) + r*x2(c.bY_Index);
                    x2(edgeI.bY_Index) = x2(edgeI.aY_Index) + r*x2(c.aX_Index) - r*x2(c.bX_Index) - s*x2(c.aY_Index) + s*x2(c.bY_Index);
                }
                //
                
                
                
                // Papers equations
                double lt00 = LT(0,0);
                double lt01 = LT(0,1);
                double lt10 = LT(1,0);
                double lt11 = LT(1,1);
                
                //MeshEdge edgeI = p.patchEdges[edgeCounter];
                
                if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
                {
                    if (c.bX_Index == edgeI.bX_Index
                        && c.bY_Index == edgeI.bY_Index) //c = e
                    {
                        printf("\n c = e (same edge)");
                        
                        x2(c.aX_Index) = (1-lt00)*vertexVectorX(c.bX_Index)
                        + lt00*vertexVectorX(c.aX_Index)
                        + lt01*vertexVectorX(c.aY_Index)
                        - lt01*vertexVectorX(c.bY_Index);
                        
                        x2(c.aY_Index) = (1-lt11)*vertexVectorX(c.bY_Index)
                        + lt10*vertexVectorX(c.aX_Index)
                        - lt10*vertexVectorX(c.bX_Index)
                        + lt11*vertexVectorX(c.aY_Index);
                        
                        x2(c.bX_Index) = (1-lt00)*vertexVectorX(c.aX_Index)
                        + lt00*vertexVectorX(c.bX_Index)
                        - lt01*vertexVectorX(c.aY_Index)
                        + lt01*vertexVectorX(c.bY_Index);
                        
                        x2(c.bY_Index) = (1-lt11)*vertexVectorX(c.aY_Index)
                        - lt10*vertexVectorX(c.aX_Index)
                        + lt10*vertexVectorX(c.bX_Index)
                        + lt11*vertexVectorX(c.bY_Index);
                    }
                    else //c.a = e.a
                    {
                        printf("\n c.a = e.a (top right)");
                        x2(c.aX_Index) = vertexVectorX(edgeI.bX_Index)
                        + lt00*vertexVectorX(c.aX_Index)
                        - lt00*vertexVectorX(c.bX_Index)
                        + lt01*vertexVectorX(c.aY_Index)
                        - lt01*vertexVectorX(c.bY_Index);
                        
                        x2(c.aY_Index) = vertexVectorX(edgeI.bY_Index)
                        + lt10*vertexVectorX(c.aX_Index)
                        - lt10*vertexVectorX(c.bX_Index)
                        + lt11*vertexVectorX(c.aY_Index)
                        - lt11*vertexVectorX(c.bY_Index);
                        
                        x2(edgeI.bX_Index) = (1-lt00)*vertexVectorX(c.aX_Index)
                        + lt00*vertexVectorX(c.bX_Index)
                        - lt01*vertexVectorX(c.aY_Index)
                        + lt01*vertexVectorX(c.bY_Index);
                        
                        x2(edgeI.bY_Index) = (1-lt11)*vertexVectorX(c.aY_Index)
                        - lt10*vertexVectorX(c.aX_Index)
                        + lt10*vertexVectorX(c.bX_Index)
                        + lt11*vertexVectorX(c.bY_Index);
                    }
                }
                else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index) //  c.b = e.a
                {
                    printf("\n c.b = e.a (bottom right)");
                    x2(c.bX_Index) = vertexVectorX(edgeI.bX_Index)
                    + lt00*vertexVectorX(c.aX_Index)
                    - lt00*vertexVectorX(c.bX_Index)
                    + lt01*vertexVectorX(c.aY_Index)
                    - lt01*vertexVectorX(c.bY_Index);
                    
                    x2(c.bY_Index) = vertexVectorX(edgeI.bY_Index)
                    + lt10*vertexVectorX(c.aX_Index)
                    - lt10*vertexVectorX(c.bX_Index)
                    + lt11*vertexVectorX(c.aY_Index)
                    - lt11*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.bX_Index) = (1+lt00)*vertexVectorX(c.bX_Index)
                    - lt00*vertexVectorX(c.aX_Index)
                    - lt01*vertexVectorX(c.aY_Index)
                    + lt01*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.bY_Index) = (1+lt11)*vertexVectorX(c.bY_Index)
                    - lt10*vertexVectorX(c.aX_Index)
                    + lt10*vertexVectorX(c.bX_Index)
                    - lt11*vertexVectorX(c.aY_Index);
                }
                else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index) //  c.a = e.b
                {
                    printf("\n c.a = e.b (top left)");
                    x2(edgeI.aX_Index) = (1+lt00)*vertexVectorX(c.aX_Index)
                    - lt00*vertexVectorX(c.bX_Index)
                    + lt01*vertexVectorX(c.aY_Index)
                    - lt01*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.aY_Index) = (1+lt11)*vertexVectorX(c.aY_Index)
                    + lt10*vertexVectorX(c.aX_Index)
                    - lt10*vertexVectorX(c.bX_Index)
                    - lt11*vertexVectorX(c.bY_Index);
                    
                    x2(c.aX_Index) = vertexVectorX(edgeI.aX_Index)
                    - lt00*vertexVectorX(c.aX_Index)
                    + lt00*vertexVectorX(c.bX_Index)
                    - lt01*vertexVectorX(c.aY_Index)
                    + lt01*vertexVectorX(c.bY_Index);
                    
                    x2(c.aY_Index) = vertexVectorX(edgeI.aY_Index)
                    - lt10*vertexVectorX(c.aX_Index)
                    + lt10*vertexVectorX(c.bX_Index)
                    - lt11*vertexVectorX(c.aY_Index)
                    + lt11*vertexVectorX(c.bY_Index);
                }
                else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index) //  c.b = e.b
                {
                    printf("\n c.b = e.b (bottom left)");
                    x2(edgeI.aX_Index) = (1-lt00)*vertexVectorX(c.bX_Index)
                    + lt00*vertexVectorX(c.aX_Index)
                    + lt01*vertexVectorX(c.aY_Index)
                    - lt01*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.aY_Index) = (1-lt11)*vertexVectorX(c.bY_Index)
                    + lt10*vertexVectorX(c.aX_Index)
                    - lt10*vertexVectorX(c.bX_Index)
                    + lt11*vertexVectorX(c.aY_Index);
                    
                    x2(c.bX_Index) = vertexVectorX(edgeI.aX_Index)
                    - lt00*vertexVectorX(c.aX_Index)
                    + lt00*vertexVectorX(c.bX_Index)
                    - lt01*vertexVectorX(c.aY_Index)
                    + lt01*vertexVectorX(c.bY_Index);
                    
                    x2(c.bY_Index) = vertexVectorX(edgeI.aY_Index)
                    - lt10*vertexVectorX(c.aX_Index)
                    + lt10*vertexVectorX(c.bX_Index)
                    - lt11*vertexVectorX(c.aY_Index)
                    + lt11*vertexVectorX(c.bY_Index);
                }
                else{ // edges not connected
                    //printf("\n edges not connected");
                    x2(edgeI.aX_Index) = vertexVectorX(edgeI.bX_Index)
                    + lt00*vertexVectorX(c.aX_Index)
                    - lt00*vertexVectorX(c.bX_Index)
                    + lt01*vertexVectorX(c.aY_Index)
                    - lt01*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.aY_Index) = vertexVectorX(edgeI.bY_Index)
                    + lt10*vertexVectorX(c.aX_Index)
                    - lt10*vertexVectorX(c.bX_Index)
                    + lt11*vertexVectorX(c.aY_Index)
                    - lt11*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.bX_Index) = vertexVectorX(edgeI.aX_Index)
                    - lt00*vertexVectorX(c.aX_Index)
                    + lt00*vertexVectorX(c.bX_Index)
                    - lt01*vertexVectorX(c.aY_Index)
                    + lt01*vertexVectorX(c.bY_Index);
                    
                    x2(edgeI.bY_Index) = vertexVectorX(edgeI.aY_Index)
                    - lt10*vertexVectorX(c.aX_Index)
                    + lt10*vertexVectorX(c.bX_Index)
                    - lt11*vertexVectorX(c.aY_Index)
                    + lt11*vertexVectorX(c.bY_Index);
                }
                
                
                
                
                //test linear scaling
                // my own scaling (pure linear scaling using original matrix T)
                /*
                double lt00 = T(0,0);
                double lt01 = T(0,1);
                double lt10 = T(1,0);
                double lt11 = T(1,1);
                
                MeshEdge edgeI = p.patchEdges[edgeCounter];
                double xScale = 1.0*newWidth/nOriginal;
                double yScale = 1.0*newHeight/mOriginal;
                
                if(c.aX_Index == edgeI.aX_Index && c.aY_Index == edgeI.aY_Index) //2 cases (c=e or ca=ea (top))
                {
                    if (c.bX_Index == edgeI.bX_Index
                        && c.bY_Index == edgeI.bY_Index) //c = e
                    {
                        printf("\n c = e (same edge)");
                        
                        x2(c.aX_Index) = xScale*((1-lt00)*vertexVectorX(c.bX_Index)
                                                 + lt00*vertexVectorX(c.aX_Index)
                                                 + lt01*vertexVectorX(c.aY_Index)
                                                 - lt01*vertexVectorX(c.bY_Index));
                        
                        x2(c.aY_Index) = yScale*((1-lt11)*vertexVectorX(c.bY_Index)
                                                 + lt10*vertexVectorX(c.aX_Index)
                                                 - lt10*vertexVectorX(c.bX_Index)
                                                 + lt11*vertexVectorX(c.aY_Index));
                        
                        x2(c.bX_Index) = xScale*((1-lt00)*vertexVectorX(c.aX_Index)
                                                 + lt00*vertexVectorX(c.bX_Index)
                                                 - lt01*vertexVectorX(c.aY_Index)
                                                 + lt01*vertexVectorX(c.bY_Index));
                        
                        x2(c.bY_Index) = yScale*((1-lt11)*vertexVectorX(c.aY_Index)
                                                 - lt10*vertexVectorX(c.aX_Index)
                                                 + lt10*vertexVectorX(c.bX_Index)
                                                 + lt11*vertexVectorX(c.bY_Index));
                    }
                    else //c.a = e.a
                    {
                        printf("\n c.a = e.a (top right)");
                        x2(c.aX_Index) = xScale*(vertexVectorX(edgeI.bX_Index)
                                                 + lt00*vertexVectorX(c.aX_Index)
                                                 - lt00*vertexVectorX(c.bX_Index)
                                                 + lt01*vertexVectorX(c.aY_Index)
                                                 - lt01*vertexVectorX(c.bY_Index));
                        
                        x2(c.aY_Index) = yScale*(vertexVectorX(edgeI.bY_Index)
                                                 + lt10*vertexVectorX(c.aX_Index)
                                                 - lt10*vertexVectorX(c.bX_Index)
                                                 + lt11*vertexVectorX(c.aY_Index)
                                                 - lt11*vertexVectorX(c.bY_Index));
                        
                        x2(edgeI.bX_Index) = xScale*((1-lt00)*vertexVectorX(c.aX_Index)
                                                     + lt00*vertexVectorX(c.bX_Index)
                                                     - lt01*vertexVectorX(c.aY_Index)
                                                     + lt01*vertexVectorX(c.bY_Index));
                        
                        x2(edgeI.bY_Index) = yScale*((1-lt11)*vertexVectorX(c.aY_Index)
                                                     - lt10*vertexVectorX(c.aX_Index)
                                                     + lt10*vertexVectorX(c.bX_Index)
                                                     + lt11*vertexVectorX(c.bY_Index));
                    }
                }
                else if(c.bX_Index == edgeI.aX_Index && c.bY_Index == edgeI.aY_Index) //  c.b = e.a
                {
                    printf("\n c.b = e.a (bottom right)");
                    x2(c.bX_Index) = xScale*(vertexVectorX(edgeI.bX_Index)
                                             + lt00*vertexVectorX(c.aX_Index)
                                             - lt00*vertexVectorX(c.bX_Index)
                                             + lt01*vertexVectorX(c.aY_Index)
                                             - lt01*vertexVectorX(c.bY_Index));
                    
                    x2(c.bY_Index) = yScale*(vertexVectorX(edgeI.bY_Index)
                                             + lt10*vertexVectorX(c.aX_Index)
                                             - lt10*vertexVectorX(c.bX_Index)
                                             + lt11*vertexVectorX(c.aY_Index)
                                             - lt11*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.bX_Index) = xScale*((1+lt00)*vertexVectorX(c.bX_Index)
                                                 - lt00*vertexVectorX(c.aX_Index)
                                                 - lt01*vertexVectorX(c.aY_Index)
                                                 + lt01*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.bY_Index) = yScale*((1+lt11)*vertexVectorX(c.bY_Index)
                                                 - lt10*vertexVectorX(c.aX_Index)
                                                 + lt10*vertexVectorX(c.bX_Index)
                                                 - lt11*vertexVectorX(c.aY_Index));
                }
                else if(c.aX_Index == edgeI.bX_Index && c.aY_Index == edgeI.bY_Index) //  c.a = e.b
                {
                    printf("\n c.a = e.b (top left)");
                    x2(edgeI.aX_Index) = xScale*((1+lt00)*vertexVectorX(c.aX_Index)
                                                 - lt00*vertexVectorX(c.bX_Index)
                                                 + lt01*vertexVectorX(c.aY_Index)
                                                 - lt01*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.aY_Index) = yScale*((1+lt11)*vertexVectorX(c.aY_Index)
                                                 + lt10*vertexVectorX(c.aX_Index)
                                                 - lt10*vertexVectorX(c.bX_Index)
                                                 - lt11*vertexVectorX(c.bY_Index));
                    
                    x2(c.aX_Index) = xScale*(vertexVectorX(edgeI.aX_Index)
                                             - lt00*vertexVectorX(c.aX_Index)
                                             + lt00*vertexVectorX(c.bX_Index)
                                             - lt01*vertexVectorX(c.aY_Index)
                                             + lt01*vertexVectorX(c.bY_Index));
                    
                    x2(c.aY_Index) = yScale*(vertexVectorX(edgeI.aY_Index)
                                             - lt10*vertexVectorX(c.aX_Index)
                                             + lt10*vertexVectorX(c.bX_Index)
                                             - lt11*vertexVectorX(c.aY_Index)
                                             + lt11*vertexVectorX(c.bY_Index));
                }
                else if(c.bX_Index == edgeI.bX_Index && c.bY_Index == edgeI.bY_Index) //  c.b = e.b
                {
                    printf("\n c.b = e.b (bottom left)");
                    x2(edgeI.aX_Index) = xScale*((1-lt00)*vertexVectorX(c.bX_Index)
                                                 + lt00*vertexVectorX(c.aX_Index)
                                                 + lt01*vertexVectorX(c.aY_Index)
                                                 - lt01*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.aY_Index) = yScale*((1-lt11)*vertexVectorX(c.bY_Index)
                                                 + lt10*vertexVectorX(c.aX_Index)
                                                 - lt10*vertexVectorX(c.bX_Index)
                                                 + lt11*vertexVectorX(c.aY_Index));
                    
                    x2(c.bX_Index) = xScale*(vertexVectorX(edgeI.aX_Index)
                                             - lt00*vertexVectorX(c.aX_Index)
                                             + lt00*vertexVectorX(c.bX_Index)
                                             - lt01*vertexVectorX(c.aY_Index)
                                             + lt01*vertexVectorX(c.bY_Index));
                    
                    x2(c.bY_Index) = yScale*(vertexVectorX(edgeI.aY_Index)
                                             - lt10*vertexVectorX(c.aX_Index)
                                             + lt10*vertexVectorX(c.bX_Index)
                                             - lt11*vertexVectorX(c.aY_Index)
                                             + lt11*vertexVectorX(c.bY_Index));
                }
                else{ // edges not connected
                    //printf("\n edges not connected");
                    x2(edgeI.aX_Index) = xScale*(vertexVectorX(edgeI.bX_Index)
                                                 + lt00*vertexVectorX(c.aX_Index)
                                                 - lt00*vertexVectorX(c.bX_Index)
                                                 + lt01*vertexVectorX(c.aY_Index)
                                                 - lt01*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.aY_Index) = yScale*(vertexVectorX(edgeI.bY_Index)
                                                 + lt10*vertexVectorX(c.aX_Index)
                                                 - lt10*vertexVectorX(c.bX_Index)
                                                 + lt11*vertexVectorX(c.aY_Index)
                                                 - lt11*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.bX_Index) = xScale*(vertexVectorX(edgeI.aX_Index)
                                                 - lt00*vertexVectorX(c.aX_Index)
                                                 + lt00*vertexVectorX(c.bX_Index)
                                                 - lt01*vertexVectorX(c.aY_Index)
                                                 + lt01*vertexVectorX(c.bY_Index));
                    
                    x2(edgeI.bY_Index) = yScale*(vertexVectorX(edgeI.aY_Index)
                                                 - lt10*vertexVectorX(c.aX_Index)
                                                 + lt10*vertexVectorX(c.bX_Index)
                                                 - lt11*vertexVectorX(c.aY_Index)
                                                 + lt11*vertexVectorX(c.bY_Index));
                }
                 */
            }
        }
        
    }
    
    
    printf("\n compute grid orientation terms");
    std::vector<MeshQuad>::iterator quadIter = meshQuads.begin();
    
    for(int x=0; x<numXVertices-1; x++)
    {
        for(int y=0; y<numYVertices-1; y++,quadIter++)
        {
            MeshQuad quad = *quadIter;
            
            if(y < numYVertices-2 && x < numXVertices-2)
            {
                x2(quad.tlX_Index) = x2(quad.blX_Index);
                x2(quad.tlY_Index) = x2(quad.trY_Index);
            }
            else{
                x2(quad.tlX_Index) = x2(quad.blX_Index);
                x2(quad.tlY_Index) = x2(quad.trY_Index);
                x2(quad.trX_Index) = x2(quad.brX_Index);
                x2(quad.blY_Index) = x2(quad.brY_Index);
            }
            
        }
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
