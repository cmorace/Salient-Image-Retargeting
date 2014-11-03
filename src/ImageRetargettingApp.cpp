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
    void gradientButtonClick();
    void verticalSeamButtonClick();
    void horizontalSeamButtonClick();
    
    void deleteVerticalSeamButtonClick();
    
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
    
    gl::Texture         originalTexture;
    gl::Texture         segmentedTexture;
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
    gradientTexture = gl::Texture(seamCarver->getGradientImage(gradientImage));
    
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
    gradParams->addButton( "Gradient", std::bind( &ImageRetargettingApp::gradientButtonClick, this ) );
    gradParams->addText("Time: " + to_string(seamCarver->gradTime));
    gradParams->addSeparator();
    
    gradParams->addParam( "Nb of Seams", &(seamCarver->nbOfSeams)).min( 1 ).max( 200 ).step( 1 );
    gradParams->addButton( "Vertical Seam", std::bind( &ImageRetargettingApp::verticalSeamButtonClick, this ) );
    gradParams->addButton( "Horizontal Seam", std::bind( &ImageRetargettingApp::horizontalSeamButtonClick, this ) );

    gradParams->addText("Time: " + to_string(seamCarver->carveTime));
    gradParams->addSeparator();
    gradParams->addButton( "Delete Vertical Seam", std::bind( &ImageRetargettingApp::deleteVerticalSeamButtonClick, this ) );
    gradParams->addSeparator();
    gradParams->addText("Image Width: " + to_string(gradientTexture.getWidth()));
    gradParams->addText("Image Height: " + to_string(gradientTexture.getHeight()));
    
}

void ImageRetargettingApp::updateWindows()
{
    originalImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
    retargetedImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
    segmentedImageWindow->setSize(segmentedTexture.getWidth(), segmentedTexture.getHeight());
    gradientImageWindow->setSize(gradientTexture.getWidth(), gradientTexture.getHeight());
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
    if( gradientTexture ) {
        gl::draw(gradientTexture);
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

void ImageRetargettingApp::gradientButtonClick()
{
    gradientImage =  seamCarvedImage.clone();
    gradientTexture = gl::Texture(seamCarver->getGradientImage(gradientImage));
    updateApplication();
}

void ImageRetargettingApp::verticalSeamButtonClick(){
    gradientImage = Surface(gradientTexture);
    gradientTexture = gl::Texture(seamCarver->getVerticalSeams(gradientImage,seamCarver->nbOfSeams));
    updateApplication();
}

void ImageRetargettingApp::horizontalSeamButtonClick(){
    gradientImage = Surface(gradientTexture);
    gradientTexture = gl::Texture(seamCarver->getHorizontalSeams(gradientImage,seamCarver->nbOfSeams));
    updateApplication();
}

void ImageRetargettingApp::deleteVerticalSeamButtonClick()
{
    gradientTexture = gl::Texture(seamCarver->deleteVerticalSeam(Surface(gradientTexture)));
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
