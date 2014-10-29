#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Texture.h"
#include "SaliencySegmentor.h"
#include "cinder/params/Params.h"



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
    
    void drawOriginalImageWindow();
    void drawSegmentedImageWindow();
    void drawRetargettedImageWindow();
    
    SaliencySegmentor*  segmentor;
    
    gl::Texture         originalTexture;
    gl::Texture         segmentedTexture;
    gl::Texture         retargettedTexture;
    
    params::InterfaceGlRef  origParams;
    params::InterfaceGlRef	segParams;
    params::InterfaceGlRef	retargParams;
    
    Surface             originalImage;
    
    Rectf*              retargRec;
    Rectf*              thumbnailRect;
    WindowRef           originalImageWindow;
    WindowRef           segmentedImageWindow;
    WindowRef           retargetedImageWindow;

};

void ImageRetargettingApp::setup()
{
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
    origParams = params::InterfaceGl::create( originalImageWindow, "Original Image Data", toPixels( Vec2i( 200, 400 ) ) );
    
    segmentedImageWindow = createWindow();
    segmentedImageWindow->setTitle("Segmented Image");
    segmentedImageWindow->connectDraw(&ImageRetargettingApp::drawSegmentedImageWindow, this);
    segParams = params::InterfaceGl::create( segmentedImageWindow, "Segmentation Parameters", toPixels( Vec2i( 200, 400 ) ) );
}

void ImageRetargettingApp::initTextures(fs::path path)
{
    originalTexture = gl::Texture(loadImage( path ));
    originalImage =  Surface(originalTexture);
    
    segmentedTexture = gl::Texture(segmentor->segmentImage(originalImage));
    
    retargettedTexture = gl::Texture(originalTexture);
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
    segParams->addParam( "K",&(segmentor->k_Seg) ).min( 0.0f ).max( 1200.0f ).step( 2.0f );
    segParams->addParam( "Min Size", &(segmentor->minSize_Seg));
    segParams->addButton( "Segment", std::bind( &ImageRetargettingApp::segmentButtonClick, this ) );
    segParams->addSeparator();
    segParams->addText("Image Width: " + to_string(segmentedTexture.getWidth()));
    segParams->addText("Image Height: " + to_string(segmentedTexture.getHeight()));
    segParams->addText("NB of Segments: " + to_string(segmentor->nbOfSegments));
    segParams->addText("Time: " + to_string(segmentor->segTime));
}

void ImageRetargettingApp::updateWindows()
{
    originalImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
    retargetedImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
    segmentedImageWindow->setSize(originalTexture.getWidth(), originalTexture.getHeight());
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
    if( retargettedTexture ) {
        gl::draw(retargettedTexture, *retargRec);
        thumbnailRect->set(10, (this->getWindowSize()).y-50, 50, (this->getWindowSize()).y-10);
        gl::draw(retargettedTexture, *thumbnailRect);
    }
}


//EVENT LISTENERS
void ImageRetargettingApp::segmentButtonClick()
{
    originalImage =  Surface(originalTexture);
    segmentedTexture = gl::Texture(segmentor->segmentImage(originalImage));
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
            Surface s8( retargettedTexture );
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
