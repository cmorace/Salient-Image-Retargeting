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
    void resizeWarp();
    
private:
    enum class SeamCarvingState
    {
        ShowImage,
        ShowGradient,
        ShowSeamCarving,
        Undefined
    };
    SeamCarvingState seamCarvingState = SeamCarvingState::ShowImage;
    
    enum class MeshWarpingState
    {
        ShowImage,
        ShowSegmentedImage,
        ShowSaliencyMap,
        ShowMeshWarping,
        Undefined
    };
    MeshWarpingState meshWarpingState = MeshWarpingState::ShowImage;
    
    void initData();
    void initWindows();
    void initTextures(fs::path path);
    
    void updateData();
    void updateWindows();
    void updateApplication();
    
    //Warping
    void segmentButtonClick();
    void getMeshButtonClick();
    //SeamCarving
    void seamCarveResetButtonClick();
    void sobelGradientButtonClick();
    void scharrGradientButtonClick();
    void sobelSaliencyButtonClick();
    void scharrSaliencyButtonClick();
    
    void verticalSeamGradientButtonClick();
    void horizontalSeamGradientButtonClick();
    void showCurrentSeamButtonClick();
    void deleteCurrentSeamButtonClick();
    void addCurrentSeamButtonClick();
    void resizeSeamButtonClick();
    ////////////////////////////////////////////////
    
    void drawOriginalImageWindow();
    void drawSegmentedImageWindow();
    void drawGradientImageWindow();
    void drawRetargettedImageWindow();
    
    SaliencySegmentor*  segmentor;
    MeshWarpRetargetter* meshWarpRetargetter;
    GradientSeamCarver*  seamCarver;
    
    Surface             originalImage;
    Surface             segmentedImage;
    Surface             saliencyImage;
    Surface             gradientImage;
    Surface             seamCarvedImage;
    Surface             seamCarvedImageCopy;
    
    gl::Texture         originalTexture;
    gl::Texture         segmentedTexture;
    gl::Texture         saliencyTexture;
    gl::Texture         gradientTexture;
    gl::Texture         seamCarvedTexture;
    gl::Texture         retargetedTexture;
    
    params::InterfaceGlRef  origParams;
    params::InterfaceGlRef	segParams;
    params::InterfaceGlRef	gradParams;
    params::InterfaceGlRef	retargParams;
    
    Rectf*              retargRec;
    Rectf*              thumbnailRect;
    WindowRef           originalImageWindow;
    WindowRef           segmentedImageWindow;
    WindowRef           gradientImageWindow;
    WindowRef           retargetedImageWindow;
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
            updateApplication();
        }
    }
    catch( ... ) {
        console() << "unable to load the texture file!" << std::endl;
    }
}


void ImageRetargetingApp::initData()
{
    retargRec = new Rectf (0,0,0,0);
    thumbnailRect = new Rectf (0,0,0,0);
    segmentor = new SaliencySegmentor();
    meshWarpRetargetter = new MeshWarpRetargetter();
    seamCarver = new GradientSeamCarver();
    
}


void ImageRetargetingApp::initWindows()
{
    retargetedImageWindow = this->getWindow();
    retargetedImageWindow->setTitle("Retargetted Image");
    retargetedImageWindow->connectDraw(&ImageRetargetingApp::drawRetargettedImageWindow, this);
    retargetedImageWindow->connectMouseDown(&ImageRetargetingApp::mouseDownRetarget, this);
    retargetedImageWindow->connectMouseDrag(&ImageRetargetingApp::mouseDragRetarget, this);
    
    originalImageWindow = createWindow();
    originalImageWindow->setTitle("Original Image");
    originalImageWindow->connectDraw(&ImageRetargetingApp::drawOriginalImageWindow, this);
    origParams = params::InterfaceGl::create( originalImageWindow, "Original Image Data", toPixels( ci::Vec2i( 200, 400 ) ) );
    
    gradientImageWindow = createWindow();
    gradientImageWindow->setTitle("Gradient / Seam Carving");
    gradientImageWindow->connectDraw(&ImageRetargetingApp::drawGradientImageWindow, this);
    gradParams = params::InterfaceGl::create( gradientImageWindow, "Gradient Parameters", toPixels( ci::Vec2i( 200, 400 ) ) );
    
    segmentedImageWindow = createWindow();
    segmentedImageWindow->setTitle("Segmented Image");
    segmentedImageWindow->connectDraw(&ImageRetargetingApp::drawSegmentedImageWindow, this);
    segmentedImageWindow->connectResize(&ImageRetargetingApp::resizeWarp, this);
    segParams = params::InterfaceGl::create( segmentedImageWindow, "Segmentation Parameters", toPixels( ci::Vec2i( 200, 400 ) ) );
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
    
    retargetedTexture = gl::Texture(originalImage.clone());
    retargRec->set(0, 0, originalTexture.getWidth(), originalTexture.getHeight());
}


void ImageRetargetingApp::updateData()
{
    origParams->clear();
    origParams->addText("Image Width: " + to_string(originalImage.getWidth()));
    origParams->addText("Image Height: " + to_string(originalImage.getHeight()));
    origParams->addText("Has Alpha: " + to_string(originalImage.hasAlpha()));
    origParams->addText("Row Bytes: " + to_string(originalImage.getRowBytes()));
    
    segParams->clear();
    segParams->addParam( "Smoothing", &(segmentor->sigma_Seg) ).min( 0.01f ).max( 1.00f ).step( 0.01f );
    segParams->addParam( "K",&(segmentor->k_Seg) ).min( 0.0f ).max( 1500.0f ).step( 2.0f );
    segParams->addParam( "Min Size", &(segmentor->minSize_Seg)).min( 0 ).max( 1500 ).step( 2 );
    segParams->addButton( "Segment", std::bind( &ImageRetargetingApp::segmentButtonClick, this ) );
    segParams->addText("NB of Segments: " + to_string(segmentor->nbOfSegments));
    segParams->addText("Time: " + to_string(segmentor->segTime));
    segParams->addSeparator();
    segParams->addParam( "Scale", &(segmentor->scale) ).min( 1 ).max( 10 ).step( 1 );
    segParams->addParam( "Delta", &(segmentor->delta) ).min( 0 ).max( 100 ).step( 1 );
    segParams->addButton( "Sobel Gradient", std::bind( &ImageRetargetingApp::sobelSaliencyButtonClick, this ) );
    segParams->addButton( "Scharr Gradient", std::bind( &ImageRetargetingApp::scharrSaliencyButtonClick, this ) );
    segParams->addText("Time: " + to_string(segmentor->saliencyTime));
    segParams->addSeparator();
    segParams->addParam( "Quad Size", &(meshWarpRetargetter->quadSize) ).min( 10 ).max( 100 ).step( 1 );
    segParams->addButton( "GetMesh", std::bind( &ImageRetargetingApp::getMeshButtonClick, this ) );
    segParams->addSeparator();
    segParams->addText("Image Width: " + to_string(segmentedTexture.getWidth()));
    segParams->addText("Image Height: " + to_string(segmentedTexture.getHeight()));
    
    gradParams->clear();
    gradParams->addButton( "Original Image", std::bind( &ImageRetargetingApp::seamCarveResetButtonClick, this ) );
    gradParams->addSeparator();
    gradParams->addParam( "Scale", &(seamCarver->scale) ).min( 1 ).max( 10 ).step( 1 );
    gradParams->addParam( "Delta", &(seamCarver->delta) ).min( 0 ).max( 100 ).step( 1 );
    gradParams->addButton( "Sobel Gradient", std::bind( &ImageRetargetingApp::sobelGradientButtonClick, this ) );
    gradParams->addButton( "Scharr Gradient", std::bind( &ImageRetargetingApp::scharrGradientButtonClick, this ) );
    gradParams->addText("Time: " + to_string(seamCarver->gradTime));
    gradParams->addSeparator();
    gradParams->addParam( "Nb of Seams", &(seamCarver->nbOfSeams)).min( 1 ).max( 500 ).step( 1 );
    gradParams->addButton( "Get Vertical Seam", std::bind( &ImageRetargetingApp::verticalSeamGradientButtonClick, this ) );
    gradParams->addButton( "Get Horizontal Seam", std::bind( &ImageRetargetingApp::horizontalSeamGradientButtonClick, this ) );
    gradParams->addText("Time: " + to_string(seamCarver->seamTime));
    gradParams->addSeparator();
    
    gradParams->addButton( "Show Seam", std::bind( &ImageRetargetingApp::showCurrentSeamButtonClick, this ) );
    gradParams->addButton( "Delete Seam", std::bind( &ImageRetargetingApp::deleteCurrentSeamButtonClick, this ) );
    gradParams->addButton( "Add Seam", std::bind( &ImageRetargetingApp::addCurrentSeamButtonClick, this ) );
    gradParams->addSeparator();
    
    gradParams->addText("Image Width: " + to_string(seamCarvedTexture.getWidth()));
    gradParams->addText("Image Height: " + to_string(seamCarvedTexture.getHeight()));
    gradParams->addParam( "Desired Width", &(seamCarver->newWidth) ).min( 1 ).max( 1000 ).step( 1 );
    gradParams->addParam( "Desired Height", &(seamCarver->newHeight) ).min( 1 ).max( 1000 ).step( 1 );
    gradParams->addButton( "Resize", std::bind( &ImageRetargetingApp::resizeSeamButtonClick, this ) );
    
}

void ImageRetargetingApp::updateWindows()
{
    if(originalImageWindow){
      originalImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
    }
    if(retargetedImageWindow){
      retargetedImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
    }
    if(segmentedImageWindow){
        segmentedImageWindow->setSize(segmentedTexture.getWidth(), segmentedTexture.getHeight());
    }
    if(gradientImageWindow){
        gradientImageWindow->setSize(seamCarvedTexture.getWidth(), seamCarvedTexture.getHeight());
    }
}

void ImageRetargetingApp::updateApplication()
{
    updateWindows();
    updateData();
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// WINDOW DRAWING
/////////////////////////////////////////////////////////////////////////////////////////////////

void ImageRetargetingApp::drawOriginalImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    if( originalTexture ) {
        gl::draw(originalTexture);
    }
    origParams->draw();
}


void ImageRetargetingApp::drawSegmentedImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    
    switch (meshWarpingState) {
        case MeshWarpingState::ShowImage :
            if( segmentedTexture ) {
                gl::draw(segmentedTexture);
            }
            break;
            
        case MeshWarpingState::ShowSegmentedImage :
            if( segmentedTexture ) {
                gl::draw(segmentedTexture);
            }
            break;
            
        case MeshWarpingState::ShowMeshWarping :
            if( segmentedTexture ) {
                meshWarpRetargetter->draw(segmentedTexture);
            }
            break;
            
        default:
            break;
    }
    
    /*
    if( segmentedTexture ) {
        gl::draw(segmentedTexture);
    }
    */
    segParams->draw();
}

void ImageRetargetingApp::drawRetargettedImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    if( retargetedTexture ) {
        gl::draw(retargetedTexture, *retargRec);
        thumbnailRect->set(10, (this->getWindowSize()).y-50, 50, (this->getWindowSize()).y-10);
        gl::draw(retargetedTexture, *thumbnailRect);
    }
}

void ImageRetargetingApp::drawGradientImageWindow()
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
            
        case SeamCarvingState::ShowSeamCarving:
            if(seamCarvedTexture){
                int dw = seamCarver->newWidth-seamCarvedImage.getWidth();
                int dh = seamCarver->newHeight-seamCarvedImage.getHeight();
                if (dw<0 && dh<0){
                    //seamCarvedImage = seamCarver->deleteMinEnergySeam(seamCarvedImage);
                    seamCarvedImage = seamCarver->deleteVerticalSeam(seamCarvedImage);
                }
                else if (dw<0){
                    seamCarvedImage = seamCarver->deleteVerticalSeam(seamCarvedImage);
                }
                else if (dh<0){
                    seamCarvedImage = seamCarver->deleteHorizontalSeam(seamCarvedImage);
                }
                //TODO:: Add Seam
                else if (dw>=0 && dh>=0) {
                    updateApplication();
                    seamCarvingState = SeamCarvingState::ShowImage;
                }
                seamCarvedTexture = gl::Texture(seamCarvedImage);
                updateData();
                gl::draw(seamCarvedTexture);
            }
            break;
        case SeamCarvingState::Undefined:
            break;
    }
    gradParams->draw();
}


//                        Segmentation GUI
//==============================================================================
void ImageRetargetingApp::segmentButtonClick()
{
    segmentedImage =  originalImage.clone();
    segmentedTexture = gl::Texture(segmentor->segmentImage(segmentedImage));
    meshWarpingState = MeshWarpingState::ShowSegmentedImage;
    updateApplication();
}

void ImageRetargetingApp::getMeshButtonClick()
{
    //meshWarpRetargetter->setMesh();
    //meshWarpRetargetter->setTexture(segmentedImage);
    meshWarpRetargetter->initMesh(segmentedImage.getWidth(),segmentedImage.getHeight());
    meshWarpingState = MeshWarpingState::ShowMeshWarping;
    updateApplication();
}


//                        RESET
//===============================================================
void ImageRetargetingApp::seamCarveResetButtonClick()
{
    seamCarvingState = SeamCarvingState::ShowImage;
    gradientImage =  seamCarver->getGradientImage(originalImage.clone()) ;
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    gradientTexture = gl::Texture(gradientImage);
    updateApplication();
}

//                        EDGE DETECTION
//===============================================================
void ImageRetargetingApp::sobelGradientButtonClick()
{
    seamCarvingState = SeamCarvingState::ShowGradient;
    gradientTexture = gl::Texture(seamCarver->getGradientImage(seamCarvedImage, GradientSeamCarver::EdgeDetection::Sobel));
    updateApplication();
}

void ImageRetargetingApp::scharrGradientButtonClick()
{
    seamCarvingState = SeamCarvingState::ShowGradient;
    gradientTexture = gl::Texture(seamCarver->getGradientImage(seamCarvedImage, GradientSeamCarver::EdgeDetection::Scharr));
    updateApplication();
}

//                        SALIENCY
//===============================================================
void ImageRetargetingApp::sobelSaliencyButtonClick()
{
   // seamCarvingState = SeamCarvingState::;
    gradientTexture = gl::Texture(seamCarver->getGradientImage(seamCarvedImage, GradientSeamCarver::EdgeDetection::Sobel));
    updateApplication();
}

void ImageRetargetingApp::scharrSaliencyButtonClick()
{
    seamCarvingState = SeamCarvingState::ShowGradient;
    gradientTexture = gl::Texture(seamCarver->getGradientImage(seamCarvedImage, GradientSeamCarver::EdgeDetection::Scharr));
    updateApplication();
}

//                        SEAM CARVING
//===============================================================

void ImageRetargetingApp::verticalSeamGradientButtonClick(){
    seamCarvingState = SeamCarvingState::ShowGradient;
    gradientImage = seamCarver->drawVerticalSeamsGradient();
    gradientTexture = gl::Texture(gradientImage);
    updateApplication();
}

void ImageRetargetingApp::horizontalSeamGradientButtonClick(){
    seamCarvingState = SeamCarvingState::ShowGradient;
    gradientImage = seamCarver->drawHorizontalSeamsGradient();
    gradientTexture = gl::Texture(gradientImage);
    updateApplication();
}

void ImageRetargetingApp::showCurrentSeamButtonClick(){
    seamCarvingState = SeamCarvingState::ShowImage;
    seamCarver->drawCurrentSeam(seamCarvedImageCopy = seamCarvedImage.clone());
    seamCarvedTexture = gl::Texture(seamCarvedImageCopy);
    updateApplication();
}

void ImageRetargetingApp::deleteCurrentSeamButtonClick(){
    seamCarvingState = SeamCarvingState::ShowImage;
    seamCarvedImage = seamCarver->deleteCurrentSeam(seamCarvedImage);
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    updateApplication();
}

//-------TODO:
void ImageRetargetingApp::addCurrentSeamButtonClick(){
    updateApplication();
}

void ImageRetargetingApp::resizeSeamButtonClick()
{
    seamCarvingState = SeamCarvingState::ShowSeamCarving;
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    updateApplication();
}

//                        KEY
//===============================================================

void ImageRetargetingApp::keyDown( KeyEvent event )
{
    if( event.getChar() == 'f' ) {
        setFullScreen( ! isFullScreen() );
    }
    else if( event.getCode() == app::KeyEvent::KEY_ESCAPE ) {
        setFullScreen( false );
    }
    else if( event.getChar() == 'o' ) {
        fs::path path = getOpenFilePath( "", ImageIo::getLoadExtensions() );
        if( ! path.empty() ){
            initTextures(path);
            updateApplication();
        }
    }
    else if( event.getChar() == 's' ) {
        fs::path path = getSaveFilePath();
        if( ! path.empty() ) {
            Surface s8( retargetedTexture );
            writeImage( writeFile( path ), s8 );
        }
    }
}

//                        MOUSE
//===============================================================

void ImageRetargetingApp::fileDrop( FileDropEvent event )
{
    try {
        initTextures(event.getFile(0));
        updateApplication();
    }
    catch( ... ) {
        console() << "unable to load the texture file!" << std::endl;
    };
}

void ImageRetargetingApp::mouseDownRetarget( MouseEvent event )
{
    int x = event.getX();
    int y = event.getY();
    retargRec->set(x,y,x,y);
}

void ImageRetargetingApp::mouseDragRetarget( MouseEvent event )
{
    int x = event.getX();
    int y = event.getY();
    retargRec->set(retargRec->x1, retargRec->y1, x, y);
}

void ImageRetargetingApp::resizeWarp()
{
    // tell the warps our window has been resized, so they properly scale up or down
    
    meshWarpRetargetter->resize();
}


CINDER_APP_NATIVE( ImageRetargetingApp, RendererGl )
