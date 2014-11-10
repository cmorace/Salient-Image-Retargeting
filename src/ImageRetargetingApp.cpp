#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Texture.h"
#include "cinder/params/Params.h"

#include "SaliencySegmentor.h"
#include "GradientSeamCarver.h"


using namespace ci;
using namespace ci::app;
using namespace std;

class ImageRetargettingApp : public AppNative {
public:
    void setup();
    void keyDown( KeyEvent event );
    void fileDrop( FileDropEvent event );
    void mouseDownRetarget( MouseEvent event );
    void mouseDragRetarget( MouseEvent event );
    
private:
    void initData();
    void initWindows();
    void initTextures(fs::path path);
    
    void updateData();
    void updateWindows();
    void updateApplication();
    
    void segmentButtonClick();
    void sobelGradientButtonClick();
    void scharrGradientButtonClick();
    void verticalSeamGradientButtonClick();
    void horizontalSeamGradientButtonClick();
    void verticalSeamImageButtonClick();
    void horizontalSeamImageButtonClick();
    void resizeSeamButtonClick();
    
    void getVerticalSeamButtonClick();
    void deleteVerticalSeamButtonClick();
    void getHorizontalSeamButtonClick();
    void deleteHorizontalSeamButtonClick();
    
    void drawOriginalImageWindow();
    void drawSegmentedImageWindow();
    void drawGradientImageWindow();
    void drawRetargettedImageWindow();
    
    SaliencySegmentor*  segmentor;
    GradientSeamCarver*  seamCarver;
    
    Surface             originalImage;
    Surface             segmentedImage;
    Surface             gradientImage;
    Surface             seamCarvedImage;
    Surface             seamCarvedImageCopy;
    
    
    gl::Texture         originalTexture;
    gl::Texture         segmentedTexture;
    gl::Texture         gradientTexture;
    gl::Texture         seamCarvedTexture;
    gl::Texture         seamCarvedTextureCopy;
    
    bool showGradient = false;
    bool isSeamCarving = false;
    
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

void ImageRetargettingApp::setup()
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


void ImageRetargettingApp::initData()
{
    retargRec = new Rectf (0,0,0,0);
    thumbnailRect = new Rectf (0,0,0,0);
    segmentor = new SaliencySegmentor();
    seamCarver = new GradientSeamCarver();
}


void ImageRetargettingApp::initWindows()
{
    retargetedImageWindow = this->getWindow();
    retargetedImageWindow->setTitle("Retargetted Image");
    retargetedImageWindow->connectDraw(&ImageRetargettingApp::drawRetargettedImageWindow, this);
    retargetedImageWindow->connectMouseDown(&ImageRetargettingApp::mouseDownRetarget, this);
    retargetedImageWindow->connectMouseDrag(&ImageRetargettingApp::mouseDragRetarget, this);
    
    originalImageWindow = createWindow();
    originalImageWindow->setTitle("Original Image");
    originalImageWindow->connectDraw(&ImageRetargettingApp::drawOriginalImageWindow, this);
    origParams = params::InterfaceGl::create( originalImageWindow, "Original Image Data", toPixels( ci::Vec2i( 200, 400 ) ) );
    
    segmentedImageWindow = createWindow();
    segmentedImageWindow->setTitle("Segmented Image");
    segmentedImageWindow->connectDraw(&ImageRetargettingApp::drawSegmentedImageWindow, this);
    segParams = params::InterfaceGl::create( segmentedImageWindow, "Segmentation Parameters", toPixels( ci::Vec2i( 200, 400 ) ) );
    
    gradientImageWindow = createWindow();
    gradientImageWindow->setTitle("Gradient / Seam Carving");
    gradientImageWindow->connectDraw(&ImageRetargettingApp::drawGradientImageWindow, this);
    gradParams = params::InterfaceGl::create( gradientImageWindow, "Gradient Parameters", toPixels( ci::Vec2i( 200, 400 ) ) );
}


void ImageRetargettingApp::initTextures(fs::path path)
{
    originalTexture = gl::Texture(loadImage( path ));
    originalImage =  Surface(originalTexture);
    
    gradientImage = originalImage.clone();
    gradientTexture = gl::Texture(gradientImage);;
    
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    
    segmentedImage = originalImage.clone();
    segmentedTexture = gl::Texture(segmentedImage);
    
    retargetedTexture = gl::Texture(originalImage.clone());
    retargRec->set(0, 0, originalTexture.getWidth(), originalTexture.getHeight());
    
    
}


void ImageRetargettingApp::updateData()
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
    segParams->addButton( "Segment", std::bind( &ImageRetargettingApp::segmentButtonClick, this ) );
    segParams->addText("NB of Segments: " + to_string(segmentor->nbOfSegments));
    segParams->addText("Time: " + to_string(segmentor->segTime));
    segParams->addSeparator();
    segParams->addText("Image Width: " + to_string(segmentedTexture.getWidth()));
    segParams->addText("Image Height: " + to_string(segmentedTexture.getHeight()));
    
    gradParams->clear();
    gradParams->addParam( "Scale", &(seamCarver->scale) ).min( 1 ).max( 50 ).step( 1 );
    gradParams->addParam( "Delta", &(seamCarver->delta) ).min( 0 ).max( 100 ).step( 1 );
    gradParams->addButton( "Sobel Gradient", std::bind( &ImageRetargettingApp::sobelGradientButtonClick, this ) );
    gradParams->addButton( "Scharr Gradient", std::bind( &ImageRetargettingApp::scharrGradientButtonClick, this ) );
    gradParams->addText("Time: " + to_string(seamCarver->gradTime));
    
    gradParams->addSeparator();
    
    gradParams->addParam( "Nb of Seams", &(seamCarver->nbOfSeams)).min( 1 ).max( 200 ).step( 1 );
    gradParams->addButton( "Vertical Seam Gradient", std::bind( &ImageRetargettingApp::verticalSeamGradientButtonClick, this ) );
    gradParams->addButton( "Horizontal Seam Gradient", std::bind( &ImageRetargettingApp::horizontalSeamGradientButtonClick, this ) );
    gradParams->addButton( "Vertical Seam Image", std::bind( &ImageRetargettingApp::verticalSeamImageButtonClick, this ) );
    gradParams->addButton( "Horizontal Seam Image", std::bind( &ImageRetargettingApp::horizontalSeamImageButtonClick, this ) );
    gradParams->addText("Time: " + to_string(seamCarver->carveTime));
    
    gradParams->addSeparator();
    
    gradParams->addButton( "Get Vertical Seam", std::bind( &ImageRetargettingApp::getVerticalSeamButtonClick, this ) );
    gradParams->addButton( "Delete Vertical Seam", std::bind( &ImageRetargettingApp::deleteVerticalSeamButtonClick, this ) );
    gradParams->addButton( "Get Horizontal Seam", std::bind( &ImageRetargettingApp::getHorizontalSeamButtonClick, this ) );
    gradParams->addButton( "Delete Horizontal Seam", std::bind( &ImageRetargettingApp::deleteHorizontalSeamButtonClick, this ) );
    
    gradParams->addSeparator();
    
    gradParams->addText("Image Width: " + to_string(seamCarvedTexture.getWidth()));
    gradParams->addText("Image Height: " + to_string(seamCarvedTexture.getHeight()));
    gradParams->addParam( "Desired Width", &(seamCarver->newWidth) ).min( 1 ).max( 1000 ).step( 1 );
    gradParams->addParam( "Desired Height", &(seamCarver->newHeight) ).min( 1 ).max( 1000 ).step( 1 );
    gradParams->addButton( "Resize", std::bind( &ImageRetargettingApp::resizeSeamButtonClick, this ) );
    
}

void ImageRetargettingApp::updateWindows()
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


void ImageRetargettingApp::updateApplication()
{
    updateWindows();
    updateData();
}


void ImageRetargettingApp::drawOriginalImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    if( originalTexture ) {
        gl::draw(originalTexture);
        origParams->draw();
    }
}


void ImageRetargettingApp::drawSegmentedImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    if( segmentedTexture ) {
        gl::draw(segmentedTexture);
        segParams->draw();
    }
}

void ImageRetargettingApp::drawRetargettedImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    if( retargetedTexture ) {
        gl::draw(retargetedTexture, *retargRec);
        thumbnailRect->set(10, (this->getWindowSize()).y-50, 50, (this->getWindowSize()).y-10);
        gl::draw(retargetedTexture, *thumbnailRect);
    }
}

void ImageRetargettingApp::drawGradientImageWindow()
{
    gl::clear( Color( 0.f, 0.f, 0.f ) );
    if(showGradient){
        if( gradientTexture ) {
            gl::draw(gradientTexture);
            gradParams->draw();
        }
    }
    else if( seamCarvedTexture ) {
        if(isSeamCarving){
            printf("carving");
            int dw = seamCarver->newWidth-seamCarvedImage.getWidth();
            int dh = seamCarver->newHeight-seamCarvedImage.getHeight();
            
            if (dw<0 && dh<0){
                seamCarvedImage = seamCarver->deleteMinEnergySeam(seamCarvedImage);
            }
            else if (dw<0){
               seamCarvedImage = seamCarver->deleteVerticalSeam(seamCarvedImage);
            }
            else if (dh<0){
               seamCarvedImage = seamCarver->deleteHorizontalSeam(seamCarvedImage);
            }
            /*
            else if(dw>0){
                seamCarvedImage = seamCarver->addVerticalSeam(seamCarvedImage,dw);
            }
             */
            
            /*
            else if (dh>0){
                seamCarvedImage = seamCarver->addHorizontalSeam(seamCarvedImage,dw);
            }
             */
            
            if (seamCarvedImage.getWidth() == seamCarver->newWidth//){
                && seamCarvedImage.getHeight() == seamCarver->newHeight) {
                
                updateApplication();
                isSeamCarving = false;
            }
            seamCarvedTexture = gl::Texture(seamCarvedImage);
            updateData();
        }
        gl::draw(seamCarvedTexture);
        gradParams->draw();
    }
}


//EVENT LISTENERS
void ImageRetargettingApp::segmentButtonClick()
{
    segmentedImage =  originalImage.clone();
    segmentedTexture = gl::Texture(segmentor->segmentImage(segmentedImage));
    updateApplication();
}

void ImageRetargettingApp::sobelGradientButtonClick()
{
    showGradient = true;
    gradientImage =  originalImage.clone();
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    gradientTexture = gl::Texture(seamCarver->getGradientImage(gradientImage, true));
    updateApplication();
}

void ImageRetargettingApp::scharrGradientButtonClick()
{
    showGradient = true;
    gradientImage =  originalImage.clone();
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    gradientTexture = gl::Texture(seamCarver->getGradientImage(gradientImage, false));
    updateApplication();
}

void ImageRetargettingApp::verticalSeamGradientButtonClick(){
    showGradient = true;
    gradientImage = Surface(originalImage.clone());
    //seamCarvedImage = originalImage.clone();
    gradientTexture = gl::Texture(seamCarver->drawVerticalSeams(seamCarvedImage,seamCarver->nbOfSeams,showGradient));
    updateApplication();
}

void ImageRetargettingApp::horizontalSeamGradientButtonClick(){
    showGradient = true;
    //gradientImage = Surface(gradientTexture);
    seamCarvedImage = originalImage.clone();
    gradientTexture = gl::Texture(seamCarver->drawHorizontalSeams(seamCarvedImage,seamCarver->nbOfSeams,showGradient));
    updateApplication();
}

void ImageRetargettingApp::verticalSeamImageButtonClick(){
    showGradient = false;
    //seamCarvedImage = Surface(seamCarvedTexture);
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarver->drawVerticalSeams(seamCarvedImage,seamCarver->nbOfSeams,showGradient));
    updateApplication();
}

void ImageRetargettingApp::horizontalSeamImageButtonClick(){
    showGradient = false;
    //seamCarvedImage = Surface(seamCarvedTexture);
    seamCarvedImage = originalImage.clone();
    seamCarvedTexture = gl::Texture(seamCarver->drawHorizontalSeams(seamCarvedImage,seamCarver->nbOfSeams,showGradient));
    updateApplication();
}

void ImageRetargettingApp::getVerticalSeamButtonClick(){
    showGradient = false;
    seamCarvedImageCopy = seamCarvedImage.clone();
    seamCarvedTexture = gl::Texture(seamCarver->drawVerticalSeams(seamCarvedImage,1,showGradient));
    seamCarvedImage = Surface(seamCarvedTexture);
    updateApplication();
}

void ImageRetargettingApp::deleteVerticalSeamButtonClick()
{
    showGradient = false;
    seamCarvedTexture = gl::Texture(seamCarver->deleteVerticalSeam(seamCarvedImageCopy));
    seamCarvedImage = Surface(seamCarvedTexture);
    updateApplication();
}

void ImageRetargettingApp::getHorizontalSeamButtonClick(){
    showGradient = false;
    seamCarvedImageCopy = seamCarvedImage.clone();
    seamCarvedTexture = gl::Texture(seamCarver->drawHorizontalSeams(seamCarvedImage,1,showGradient));
    seamCarvedImage = Surface(seamCarvedTexture);
    updateApplication();
}

void ImageRetargettingApp::deleteHorizontalSeamButtonClick()
{
    showGradient = false;
    seamCarvedTexture = gl::Texture(seamCarver->deleteHorizontalSeam(seamCarvedImageCopy));
    seamCarvedImage = Surface(seamCarvedTexture);
    updateApplication();
}

void ImageRetargettingApp::resizeSeamButtonClick()
{
    showGradient = false;
    //int maxSeams = seamCarvedImage.getWidth() - seamCarver->newWidth;
    isSeamCarving = true;
    seamCarvedTexture = gl::Texture(seamCarvedImage);
    printf("carving Call");
    updateApplication();
}

void ImageRetargettingApp::keyDown( KeyEvent event )
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

void ImageRetargettingApp::fileDrop( FileDropEvent event )
{
    try {
        initTextures(event.getFile(0));
        updateApplication();
    }
    catch( ... ) {
        console() << "unable to load the texture file!" << std::endl;
    };
}


void ImageRetargettingApp::mouseDownRetarget( MouseEvent event )
{
    int x = event.getX();
    int y = event.getY();
    retargRec->set(x,y,x,y);
}

void ImageRetargettingApp::mouseDragRetarget( MouseEvent event )
{
    int x = event.getX();
    int y = event.getY();
    retargRec->set(retargRec->x1, retargRec->y1, x, y);
}


CINDER_APP_NATIVE( ImageRetargettingApp, RendererGl )
