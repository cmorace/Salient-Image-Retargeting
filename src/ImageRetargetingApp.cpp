#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Texture.h"
#include "cinder/params/Params.h"
#include "SaliencySegmentor.h"
#include "MeshWarper.h"
#include "GradientSeamCarver.h"
#include "MeshWarpRetargetter.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class ImageRetargetingApp : public AppNative {
public:
    void setup();
    void keyDown( KeyEvent event );
    void fileDrop( FileDropEvent event );
    void mouseDownRetarget( MouseEvent event );
    void mouseDragRetarget( MouseEvent event );
    
    
private:
    enum class SeamCarvingState
    {
        ShowImage,
        ShowGradient,
        SeamCarving,
        Undefined
    };
    SeamCarvingState seamCarvingState = SeamCarvingState::ShowImage;
    
    enum class MeshWarpingState
    {
        ShowImage,
        ShowSegmentedImage,
        ShowSaliencyMap,
        ShowMesh,
        ShowMeshWarping,
        ShowPatchCenterEdge,
        ShowImageRetargeting,
        Undefined
    };
    MeshWarpingState meshWarpingState = MeshWarpingState::ShowImage;
    
    int originalWidth = 0;
    int originalHeight = 0;
    int linearScaleWidth;
    int linearScaleHeight;
    SaliencySegmentor*   saliencySegmentor;
    MeshWarpRetargetter* meshWarpRetargetter;
    GradientSeamCarver*  seamCarver;
    
    // Surfaces are raster images living on cpu;
    Surface             originalImage;
    Surface             segmentedImage;
    Surface             saliencyImage;
    Surface             gradientImage;
    Surface             seamCarvedImage;
    Surface             seamCarvedImageCopy;
    
    // gl::Textures are raster images living on gpu
    gl::Texture         originalTexture;
    gl::Texture         segmentedTexture;
    gl::Texture         saliencyTexture;
    gl::Texture         gradientTexture;
    gl::Texture         seamCarvedTexture;
    
    params::InterfaceGlRef  linearScalingParams;
    params::InterfaceGlRef	meshWarpingParams;
    params::InterfaceGlRef	seamCarvingParams;
    
    Rectf*              linearScaleRec;
    WindowRef           linearScaleWindow;
    WindowRef           meshWarpingWindow;
    WindowRef           seamCarvingWindow;
    
    void initData();
    void initWindows();
    void initTextures(fs::path path);
    void initGUI();
    void resetWindowOriginalSize(WindowRef window);
    void resetAllWindowsOriginalSize();
    
    // Linear Resizing
    void linearResizeResetButtonClick();
    void linearResizeButtonClick();
    
    //Warping
    void meshWarperResetButtonClick();
    void segmentRandomButtonClick();
    void segmentColorButtonClick();
    void segmentSaliencyButtonClick();
    void sobelSaliencyButtonClick();
    void toggleWireFrameButtonClick();
    void getPatchEdgeClick();
    void resizeMeshRect();
    void resizeMeshEllipse();
    
    //SeamCarving
    void seamCarveResetButtonClick();
    void sobelGradientButtonClick();
    void verticalSeamGradientButtonClick();
    void horizontalSeamGradientButtonClick();
    void showCurrentSeamButtonClick();
    void deleteCurrentSeamButtonClick();
    void addCurrentSeamButtonClick();
    void resizeSeamButtonClick();
    
    void drawLinearScaleWindow();
    void drawMeshWarpingWindow();
    void drawSeamCarvingWindow();
    
    void linearResizeWindowResize();
    void seamCarvingWindowResize();
    void meshWarpingWindowResize();
};


void ImageRetargetingApp::setup()
{
    printf("setup");
    initData();
    initWindows();
    try {
        fs::path path = getOpenFilePath( "", ImageIo::getLoadExtensions() );
        if( ! path.empty() ) {
            initTextures(path);
            resetAllWindowsOriginalSize();
        }
    }
    catch( ... ) {
        console() << "unable to load the texture file!" << std::endl;
    }
    initGUI();
}


void ImageRetargetingApp::initData()
{
    linearScaleRec = new Rectf(0,0,0,0);
    saliencySegmentor = new SaliencySegmentor();
    meshWarpRetargetter = new MeshWarpRetargetter();
    seamCarver = new GradientSeamCarver();
}


void ImageRetargetingApp::initWindows()
{
    linearScaleWindow = this->getWindow();
    linearScaleWindow->setTitle("Linear Scaling");
    linearScaleWindow->connectDraw(&ImageRetargetingApp::drawLinearScaleWindow, this);
    linearScaleWindow->connectResize(&ImageRetargetingApp::linearResizeWindowResize, this);
    linearScalingParams = params::InterfaceGl::create( linearScaleWindow, "Original Image Data", toPixels( ci::Vec2i( 200, 400 ) ) );
    
    seamCarvingWindow = createWindow();
    seamCarvingWindow->setTitle("Seam Carving");
    seamCarvingWindow->connectDraw(&ImageRetargetingApp::drawSeamCarvingWindow, this);
    seamCarvingWindow->connectResize(&ImageRetargetingApp::seamCarvingWindowResize, this);
    seamCarvingParams = params::InterfaceGl::create( seamCarvingWindow, "Seam Carving", toPixels( ci::Vec2i( 200, 400 ) ) );
    
    meshWarpingWindow = createWindow();
    meshWarpingWindow->setTitle("Mesh Warping");
    meshWarpingWindow->connectDraw(&ImageRetargetingApp::drawMeshWarpingWindow, this);
    meshWarpingWindow->connectResize(&ImageRetargetingApp::meshWarpingWindowResize, this);
    meshWarpingParams = params::InterfaceGl::create( meshWarpingWindow, "Mesh Warping", toPixels( ci::Vec2i( 200, 400 ) ) );
}


void ImageRetargetingApp::initTextures(fs::path path)
{
    originalTexture = gl::Texture(loadImage( path ));
    originalImage =  Surface(originalTexture);
    
    gradientImage = seamCarver->getGradientImage(originalImage.clone());
    gradientTexture = gl::Texture(gradientImage);;
    
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    
    segmentedImage = originalImage.clone();
    segmentedTexture = gl::Texture(segmentedImage);
    
    saliencyImage = saliencySegmentor->getSaliencyMap(originalImage.clone(), SaliencySegmentor::SaliencyMethod::Sobel);
    saliencyImage = saliencySegmentor->getSegmentedSalientImage(saliencyImage);
    saliencyTexture = gl::Texture(saliencyImage);
    
    originalWidth = originalImage.getWidth();
    originalHeight = originalImage.getHeight();
    
    meshWarpRetargetter->initMesh(originalWidth, originalHeight, saliencySegmentor);
    linearScaleRec->set(0, 0, originalWidth, originalHeight);
}

void ImageRetargetingApp::initGUI()
{
    linearScalingParams->addButton( "Reset Image", std::bind( &ImageRetargetingApp::linearResizeResetButtonClick, this ) );
    linearScalingParams->addSeparator();
    linearScalingParams->addParam("Original Width: ", &originalWidth, true);
    linearScalingParams->addParam("Original Height: ", &originalHeight, true);
    linearScalingParams->addParam( "Resize Width", &linearScaleWidth ).min( 1 ).max( 1200 ).step( 1 );
    linearScalingParams->addParam( "Resize Height", &linearScaleHeight ).min( 1 ).max( 750 ).step( 1 );
    linearScalingParams->addButton( "Resize", std::bind( &ImageRetargetingApp::linearResizeButtonClick, this ) );
    linearScalingParams->addSeparator();
    
    meshWarpingParams->addButton( "Reset Mesh", std::bind( &ImageRetargetingApp::meshWarperResetButtonClick, this ) );
    meshWarpingParams->addParam( "Quad Size", &(meshWarpRetargetter->quadSize) ).min( 5 ).max( 100 ).step( 1 );
    meshWarpingParams->addButton( "Toggle Wire Frame", std::bind( &ImageRetargetingApp::toggleWireFrameButtonClick, this ) );
    meshWarpingParams->addSeparator();
    meshWarpingParams->addText("Segmentation");
    meshWarpingParams->addParam( "Blur", &saliencySegmentor->segBlurDeviation ).min( 0.01f ).max( 1.00f ).step( 0.01f );
    meshWarpingParams->addParam( "K",&saliencySegmentor->segNeighborThreshold ).min( 0.0f ).max( 1500.0f ).step( 10.f );
    meshWarpingParams->addParam( "Min Size", &saliencySegmentor->segMinSize).min( 0 ).max( 1500 ).step( 10 );
    meshWarpingParams->addButton( "Segment Random", std::bind( &ImageRetargetingApp::segmentRandomButtonClick, this ) );
    meshWarpingParams->addButton( "Segment Color", std::bind( &ImageRetargetingApp::segmentColorButtonClick, this ) );
    meshWarpingParams->addSeparator();
    meshWarpingParams->addText("Saliency");
    meshWarpingParams->addParam( "Scale", &(saliencySegmentor->scale) ).min( 1 ).max( 10 ).step( 1 );
    meshWarpingParams->addButton( "Sobel Saliency", std::bind( &ImageRetargetingApp::sobelSaliencyButtonClick, this ) );
    meshWarpingParams->addButton( "Segment Saliency", std::bind( &ImageRetargetingApp::segmentSaliencyButtonClick, this ) );
    meshWarpingParams->addButton( "Edge Saliency", std::bind( &ImageRetargetingApp::getPatchEdgeClick, this ) );
    meshWarpingParams->addSeparator();
    meshWarpingParams->addText("Quad Warping");
    meshWarpingParams->addParam( "Resize Width", &(meshWarpRetargetter->resizeWidth) ).min( 1 ).max( 1200 ).step( 1 );
    meshWarpingParams->addParam( "Resize Height", &(meshWarpRetargetter->resizeHeight) ).min( 1 ).max( 750 ).step( 1 );
    meshWarpingParams->addParam( "Linear Weight", &(meshWarpRetargetter->transformationWeight) ).min( 0.f ).max( 1.f ).step( 0.1f );
    meshWarpingParams->addButton( "Resize Rect", std::bind( &ImageRetargetingApp::resizeMeshRect, this ) );
    meshWarpingParams->addButton( "Resize Ellipse", std::bind( &ImageRetargetingApp::resizeMeshEllipse, this ) );
    meshWarpingParams->addParam("Resize Time: ", &(meshWarpRetargetter->resizeTime), true);
    meshWarpingParams->addSeparator();
    
    seamCarvingParams->addButton( "Reset Image", std::bind( &ImageRetargetingApp::seamCarveResetButtonClick, this ) );
    seamCarvingParams->addSeparator();
    seamCarvingParams->addParam( "Scale", &(seamCarver->scale) ).min( 1 ).max( 20 ).step( 1 );
    seamCarvingParams->addButton( "Sobel Gradient", std::bind( &ImageRetargetingApp::sobelGradientButtonClick, this ) );
    seamCarvingParams->addSeparator();
    seamCarvingParams->addButton( "Show Vertical Seam", std::bind( &ImageRetargetingApp::verticalSeamGradientButtonClick, this ) );
    seamCarvingParams->addButton( "Show Horizontal Seam", std::bind( &ImageRetargetingApp::horizontalSeamGradientButtonClick, this ) );
    seamCarvingParams->addSeparator();
    seamCarvingParams->addButton( "Get Gradient", std::bind( &ImageRetargetingApp::sobelGradientButtonClick, this ) );
    seamCarvingParams->addButton( "Get Seam", std::bind( &ImageRetargetingApp::verticalSeamGradientButtonClick, this ) );
    seamCarvingParams->addButton( "Delete Seam", std::bind( &ImageRetargetingApp::deleteCurrentSeamButtonClick, this ) );
    //seamCarvingParams->addButton( "Add Seam", std::bind( &ImageRetargetingApp::addCurrentSeamButtonClick, this ) );
    seamCarvingParams->addSeparator();
    seamCarvingParams->addParam( "Resize Width", &(seamCarver->newWidth) ).min( 1 ).max( 1200 ).step( 1 );
    seamCarvingParams->addParam( "Resize Height", &(seamCarver->newHeight) ).min( 1 ).max( 750 ).step( 1 );
    seamCarvingParams->addButton( "Resize", std::bind( &ImageRetargetingApp::resizeSeamButtonClick, this ) );
    seamCarvingParams->addParam("Resize Time: ", &(seamCarver->carveTime), true);
    seamCarvingParams->addSeparator();
}


void ImageRetargetingApp::resetWindowOriginalSize(WindowRef window)
{
    window->setSize(originalTexture.getWidth(), originalTexture.getHeight());
}

void ImageRetargetingApp::resetAllWindowsOriginalSize()
{
    resetWindowOriginalSize(linearScaleWindow);
    resetWindowOriginalSize(seamCarvingWindow);
    resetWindowOriginalSize(meshWarpingWindow);
}
    


//                              WINDOW DRAWING
//==============================================================================

void ImageRetargetingApp::drawLinearScaleWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    
    if( originalTexture ) {
        gl::draw(originalTexture, *linearScaleRec);
    }
    linearScalingParams->draw();
}


void ImageRetargetingApp::drawMeshWarpingWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    
    switch (meshWarpingState) {
        case MeshWarpingState::ShowImage :
            
            if (originalTexture) {
                meshWarpRetargetter->drawMesh(originalTexture);
            }
            break;
            
        case MeshWarpingState::ShowSegmentedImage :
            if( segmentedTexture ) {
                meshWarpRetargetter->drawMesh(segmentedTexture);
            }
            break;
            
        case MeshWarpingState::ShowSaliencyMap :
            if( saliencyTexture ) {
                meshWarpRetargetter->drawMesh(saliencyTexture);
            }
            break;
          
        case MeshWarpingState::ShowMesh :
            if( originalTexture ) {
                meshWarpRetargetter->drawMesh(originalTexture);
            }
            break;
            
        case MeshWarpingState::ShowPatchCenterEdge :
            if( saliencyTexture ) {
                meshWarpRetargetter->drawEdges(saliencyTexture);
            }
            break;
          
        case MeshWarpingState::ShowMeshWarping :
            if( originalTexture ) {
                meshWarpRetargetter->drawMesh(originalTexture);
            }
            break;
            
        default:
            break;
    }
    meshWarpingParams->draw();
}

void ImageRetargetingApp::drawSeamCarvingWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    switch(seamCarvingState)
    {
        case SeamCarvingState::ShowImage:
            if( seamCarvedTexture ) {
                gl::draw(seamCarvedTexture);
            }
            break;
            
        case SeamCarvingState::ShowGradient:
            if( gradientTexture ) {
                gl::draw(gradientTexture);
            }
            break;
            
        case SeamCarvingState::SeamCarving:
            if(seamCarvedTexture){
                int dw = seamCarver->newWidth - seamCarvedImage.getWidth();
                int dh = seamCarver->newHeight - seamCarvedImage.getHeight();
                
                if (dw<0 && dh<0){
                    seamCarvingWindow->setSize(seamCarvedTexture.getWidth()-1, seamCarvedTexture.getHeight());
                    seamCarvedImage = seamCarver->deleteVerticalSeam(seamCarvedImage);
                }
                else if (dw<0){
                    seamCarvingWindow->setSize(seamCarvedTexture.getWidth()-1, seamCarvedTexture.getHeight());
                    seamCarvedImage = seamCarver->deleteVerticalSeam(seamCarvedImage);
                }
                else if (dh<0){
                    seamCarvingWindow->setSize(seamCarvedTexture.getWidth(), seamCarvedTexture.getHeight()-1);
                    seamCarvedImage = seamCarver->deleteHorizontalSeam(seamCarvedImage);
                }
                //TODO:: Add Seam
                else if (dw>=0 && dh>=0) {
                    seamCarvingState = SeamCarvingState::ShowImage;
                    seamCarver->stopCarveTimer();
                }
                seamCarvedTexture = gl::Texture(seamCarvedImage);
                gl::draw(seamCarvedTexture);
            }
            
            break;
        case SeamCarvingState::Undefined:
            break;
    }
    seamCarvingParams->draw();
}


//                        RESET ORIGINAL IMAGE
//==============================================================================
void ImageRetargetingApp::seamCarveResetButtonClick()
{
    gradientImage =  seamCarver->getGradientImage(originalImage.clone()) ;
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    gradientTexture = gl::Texture(gradientImage);
    seamCarvingState = SeamCarvingState::ShowImage;
    resetWindowOriginalSize(seamCarvingWindow);
}

void ImageRetargetingApp::meshWarperResetButtonClick()
{
    meshWarpingState = MeshWarpingState::ShowImage;
    saliencyImage = saliencySegmentor->getSaliencyMap(originalImage.clone(), SaliencySegmentor::SaliencyMethod::Sobel);
    saliencyImage = saliencySegmentor->getSegmentedSalientImage(saliencyImage);
    saliencyTexture = gl::Texture(saliencyImage);
    
    meshWarpRetargetter->initMesh(originalImage.getWidth(), originalImage.getHeight(), saliencySegmentor);
    meshWarpingWindow->setSize(originalImage.getWidth(), originalImage.getHeight());
}

//                        EDGE DETECTION
//==============================================================================
void ImageRetargetingApp::sobelGradientButtonClick()
{
    seamCarvingState = SeamCarvingState::ShowGradient;
    gradientTexture = gl::Texture(seamCarver->getGradientImage(seamCarvedImage, GradientSeamCarver::EdgeDetection::Sobel));
}


//                        SALIENCY
//==============================================================================
void ImageRetargetingApp::sobelSaliencyButtonClick()
{
    meshWarpingState = MeshWarpingState::ShowSaliencyMap;
    saliencyImage = saliencySegmentor->getSaliencyMap(originalImage.clone(), SaliencySegmentor::SaliencyMethod::Sobel);
    saliencyTexture = gl::Texture(saliencyImage);
}

//                        SEGMENTATION
//==============================================================================
void ImageRetargetingApp::segmentRandomButtonClick()
{
    segmentedImage = Surface(originalTexture);
    segmentedImage = saliencySegmentor->getSegmentedImage(segmentedImage);
    segmentedTexture = gl::Texture(segmentedImage);
    meshWarpingState = MeshWarpingState::ShowSegmentedImage;
}

void ImageRetargetingApp::segmentColorButtonClick()
{
    segmentedImage = Surface(originalTexture);
    segmentedImage = saliencySegmentor->getSegmentedColorImage(segmentedImage);
    segmentedTexture = gl::Texture(segmentedImage);
    meshWarpingState = MeshWarpingState::ShowSegmentedImage;
}

void ImageRetargetingApp::segmentSaliencyButtonClick()
{
    saliencyImage = Surface(originalTexture);
    saliencyImage = saliencySegmentor->getSegmentedSalientImage(saliencyImage);
    saliencyTexture = gl::Texture(saliencyImage);
    meshWarpingState = MeshWarpingState::ShowSaliencyMap;
}

//                        MESH WARPING
//==============================================================================

void ImageRetargetingApp::toggleWireFrameButtonClick()
{
    meshWarpRetargetter->isDrawingWireFrame = !meshWarpRetargetter->isDrawingWireFrame;
}

void ImageRetargetingApp::getPatchEdgeClick()
{
    meshWarpRetargetter->initMesh(saliencyImage.getWidth(),saliencyImage.getHeight(), saliencySegmentor);
    meshWarpingState = MeshWarpingState::ShowPatchCenterEdge;
}

void ImageRetargetingApp::resizeMeshRect()
{
    meshWarpRetargetter->startTimer();
    meshWarpRetargetter->resizeMeshRect(meshWarpRetargetter->resizeWidth , meshWarpRetargetter->resizeHeight);
    meshWarpRetargetter->stopTimer();
    meshWarpingState = MeshWarpingState::ShowMeshWarping;
    meshWarpingWindow->setSize(meshWarpRetargetter->resizeWidth , meshWarpRetargetter->resizeHeight);
}

void ImageRetargetingApp::resizeMeshEllipse()
{
    meshWarpRetargetter->resizeMeshEllipse(meshWarpRetargetter->resizeWidth , meshWarpRetargetter->resizeHeight);
    meshWarpingState = MeshWarpingState::ShowMeshWarping;
    meshWarpingWindow->setSize(meshWarpRetargetter->resizeWidth , meshWarpRetargetter->resizeHeight);
}



//                        SEAM CARVING BUTTON EVENTS
//==============================================================================

void ImageRetargetingApp::verticalSeamGradientButtonClick()
{
    gradientImage = seamCarver->drawVerticalSeamsGradient();
    gradientTexture = gl::Texture(gradientImage);
    seamCarvingState = SeamCarvingState::ShowGradient;
}

void ImageRetargetingApp::horizontalSeamGradientButtonClick()
{
    gradientImage = seamCarver->drawHorizontalSeamsGradient();
    gradientTexture = gl::Texture(gradientImage);
    seamCarvingState = SeamCarvingState::ShowGradient;
}

void ImageRetargetingApp::showCurrentSeamButtonClick()
{
    seamCarver->drawCurrentSeam(seamCarvedImageCopy = seamCarvedImage.clone());
    seamCarvedTexture = gl::Texture(seamCarvedImageCopy);
    seamCarvingState = SeamCarvingState::ShowImage;
}

void ImageRetargetingApp::deleteCurrentSeamButtonClick()
{
    seamCarvedImage = seamCarver->deleteCurrentSeam(seamCarvedImage);
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    seamCarvingState = SeamCarvingState::ShowImage;
}

//-------TODO:
void ImageRetargetingApp::addCurrentSeamButtonClick()
{
    
}

void ImageRetargetingApp::resizeSeamButtonClick()
{
    seamCarver->startCarveTimer();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    seamCarvingState = SeamCarvingState::SeamCarving;
}

//                          LINEAR RESIZE BUTTON EVENTS
//==============================================================================

void ImageRetargetingApp::linearResizeResetButtonClick()
{
    linearScaleRec->set(0,0,originalImage.getWidth(),originalImage.getHeight());
    linearScaleWindow->setSize(originalImage.getWidth(),originalImage.getHeight());
}


void ImageRetargetingApp::linearResizeButtonClick()
{
    linearScaleRec->set(0,0,linearScaleWidth,linearScaleHeight);
    linearScaleWindow->setSize(linearScaleWidth,linearScaleHeight);
}



//                          KEY EVENTS
//==============================================================================

void ImageRetargetingApp::keyDown( KeyEvent event )
{
    
    if( event.getChar() == 'o' ) {
        fs::path path = getOpenFilePath( "", ImageIo::getLoadExtensions() );
        if( ! path.empty() ){
            initTextures(path);
            resetAllWindowsOriginalSize();
        }
    }
    /*
    todo: use fbo to save images
    else if( event.getChar() == 's' ) {
        fs::path path = getSaveFilePath();
        if( ! path.empty() ) {
            Surface s8( retargetedTexture );
            writeImage( writeFile( path ), s8 );
        }
    }
     */
}


//                          MOUSE EVENTS
//==============================================================================

void ImageRetargetingApp::fileDrop( FileDropEvent event )
{
    try {
        initTextures(event.getFile(0));
        resetAllWindowsOriginalSize();
    }
    catch( ... ) {
        console() << "unable to load the texture file!" << std::endl;
    };
}


//                          WINDOW EVENTS
//==============================================================================

void ImageRetargetingApp::linearResizeWindowResize()
{
    linearScaleWidth = linearScaleWindow->getWidth();
    linearScaleHeight = linearScaleWindow->getHeight();
    linearScaleRec->set(0,0,linearScaleWidth,linearScaleHeight);
}

void ImageRetargetingApp::seamCarvingWindowResize()
{
    if(seamCarvingState != SeamCarvingState::SeamCarving)
    {
        seamCarver->newWidth = seamCarvingWindow->getWidth();
        seamCarver->newHeight = seamCarvingWindow->getHeight();
    }
}

void ImageRetargetingApp::meshWarpingWindowResize()
{
    meshWarpRetargetter->resizeWidth = meshWarpingWindow->getWidth();
    meshWarpRetargetter->resizeHeight = meshWarpingWindow->getHeight();
    if(meshWarpRetargetter->getNumVertices() < 500)
    {
        meshWarpRetargetter->resizeMeshRect(meshWarpRetargetter->resizeWidth , meshWarpRetargetter->resizeHeight);
        meshWarpingState = MeshWarpingState::ShowMeshWarping;
    }
}



CINDER_APP_NATIVE( ImageRetargetingApp, RendererGl )
