////////////////////////////////////////////////////////////////////////
//
//   KAIST, Spring 2023
//   CS380: Introduction to Computer Graphics
//   Instructor: Minhyuk Sung (mhsung@kaist.ac.kr)
//   Last Update: Juil Koo (63days@kaist.ac.kr)
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <list>
#include <fstream>
#include <sstream>
#include <iostream>
// If your OS is LINUX, uncomment the line below.
//#include <tr1/memory>

#include <GL/glew.h>
#ifdef __MAC__
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

#include "cvec.h"
#include "matrix4.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"

#include "rigtform.h"
#include "arcball.h"

#include "asstcommon.h"
#include "drawer.h"
#include "picker.h"
#include "scenegraph.h"

#include "sgutils.h"
#include "interpolation.h"

using namespace std;      // for string, vector, iostream, and other standard C++ stuff
// If your OS is LINUX, uncomment the line below.
//using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.0 or
// OpenGL 3.x+ with GLSL 1.3.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.3. Make sure that your machine supports the version of GLSL you
// are using. In particular, on Mac OS X currently there is no way of using
// OpenGL 3.x with GLSL 1.3 when GLUT is used.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------
const bool g_Gl2Compatible = true;


static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = -2.0;      // y coordinate of the ground
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;

static const int PICKING_SHADER = 2; // index of the picking shader is g_shaerFiles
static const int g_numShaders = 3;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"}
};
static const char* const g_shaderFilesGl2[g_numShaders][2] = {
    {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
    {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"},
    {"./shaders/basic-gl2.vshader", "./shaders/pick-gl2.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz)
    : p(x,y,z), n(nx, ny, nz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
  }
};

typedef SgGeometryShapeNode<Geometry> MyShapeNode;

// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_arcball, g_point, g_sphere;

static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_planeNode, g_robot2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do picking

// --------- Scene

static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 7.0));
static RigTForm g_objectRbt[4] = {RigTForm(Cvec3(-1, 0, 0)), RigTForm(Cvec3(1, 0, 0)), RigTForm(Cvec3(0, 0, 0)), RigTForm(Cvec3(0, 0, 0))}; // currently only 4 objs are defined
static Cvec3f g_objectColors[2] = {Cvec3f(0, 0, 1), Cvec3f(1, 0, 0)};
static Cvec3f g_arcballColor = Cvec3f(0, 1, 0);

int g_vcount = 0;
int g_mcount = 0;
static RigTForm g_orjSkyRbt = RigTForm(Cvec3(0.0, 0.25, 8.0));
static RigTForm g_orjObjectRbt[2] = {RigTForm(Cvec3(-1, 0, 0)), RigTForm(Cvec3(1, 0, 0))}; // currently only 2 objs are defined

double g_arcballScreenRadius;
double g_arcballScale;

bool g_isPicking = false;
bool g_isArcball = false;

static const bool isDebug = false;

// --------- Animation

std::list<std::vector<RigTForm> >g_keyFrames = std::list<std::vector<RigTForm > >();
std::list<std::vector<RigTForm> >::iterator g_keyFramesIter = g_keyFrames.begin();
std::string g_filename = "keyframes.txt"; 

bool g_isAnimation = false;

int g_msBetweenKeyFrames = 2000; // 2000 if 2 seconds between keyframes
static int g_msOffset = 100;
static int g_animateFramesPerSecond = 60; // frames to render per second during animation playback

static bool g_isLinear = false;

///////////////// END OF G L O B A L S //////////////////////////////////////////////////

static shared_ptr<SgRbtNode> returnEyeNode (int cnt) {
  switch (cnt % 2){
  case 0:
    return g_skyNode;
    break;
  
  case 1:
    return g_planeNode;
    break;

  default:
    throw std::runtime_error("returnEyeNode; fail to handle all cases");
    break;
  }
}

// helper function to print object
void printObject(int n) {
  cout << "Current selected object: ";
  switch(n % 3) {
  case 0:
    cout << "Camera" << endl;
    break;

  case 1:
    cout << "Cube1" << endl;
    break;

  case 2:
    cout << "Cube2" << endl;
    break;

  default:
    throw std::runtime_error("printObject; fail to handle all cases");
    break;
  }
}

// helper function to print view
void printView(int n)
{
  cout << "Current view: ";
  switch (n % 2)
  {
  case 0:
    cout << "Camera" << endl;
    break;

  case 1:
    cout << "Plane" << endl;
    break;

  default:
    throw std::runtime_error("printObject; fail to handle all cases");
    break;
  }
}

// helper function to print Matrix4 objects
void printMatrix(const string str, const Matrix4 &matrix) {
  cout << "printMatrix - " << str << ": " << endl;
  // print out the elements of the matrix
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      cout << matrix[i * 4 + j] << " ";
    }
    cout << std::endl;
  }
}

// helper function to print Matrix4 objects
void printRigTForm(const string str, const RigTForm &tform) {
  cout << "printRigTForm - " << str << ": " << endl;
  cout << "t_ : " << tform.getTranslation()[0] << ", " << tform.getTranslation()[1] << ", " << tform.getTranslation()[2] << endl;
  cout << "r_ : " << tform.getRotation()[0] << ", " << tform.getRotation()[1] << ", " << tform.getRotation()[2] << ", " << tform.getRotation()[3] << endl;
}

// helper function to print Cvec3 objects
void printCvec3(const string str, const Cvec3 &cvec) {
  cout << "printCvec3 - " << str << ": " << cvec[0] << ", " << cvec[1] << ", " << cvec[2] << endl;
}

// helper function to print Cvec3 objects
void printCvec4(const string str, const Cvec3 &cvec) {
  cout << "printCvec4 - " << str << ": " << cvec[0] << ", " << cvec[1] << ", " << cvec[2] << ", " << cvec[3] << endl;
}

// helper function to print the length of current keyframes
void printKeyFramesNum() {
  cout << "Total of "<< distance(g_keyFrames.begin(), g_keyFrames.end()) << " keyframes" << endl;
}

// helper function to print current iterator position
void printCurrentKeyFrame() {
  if (g_keyFrames.empty()) {
    cout << "Keyframes empty" << endl;
  } else {
    cout << "currently " << distance(g_keyFrames.begin(), g_keyFramesIter) + 1 << " keyframe out of " << distance(g_keyFrames.begin(), g_keyFrames.end()) << " keyframes" << endl;
  }
}
// helper function to print current interpolation mode
void printInterpolationMode(){
  if (g_isLinear) {
    cout << "Linear Interpolation" << endl;
  } else {
    cout << "Catmull-Rom Interpolation" << endl;
  }
}

// save current scene graph to vector, return vector
std::vector<RigTForm > saveToKeyFrame() {
  std::vector<std::shared_ptr<SgRbtNode> > l_nodes;
  std::vector<RigTForm > l_keyframe = std::vector<RigTForm > ();
  dumpSgRbtNodes(g_world, l_nodes);
  for (int i = 0; i < l_nodes.size(); i++) {
    l_keyframe.push_back(l_nodes[i]->getRbt());
  }
  return l_keyframe;
}

// load from vector to current scene graph
void loadFromKeyFrame(std::vector<RigTForm> l_keyframe) {
  std::vector<std::shared_ptr<SgRbtNode> > l_nodes;
  dumpSgRbtNodes(g_world, l_nodes);
  assert(l_nodes.size() == l_keyframe.size());
  for (int i = 0; i < l_nodes.size(); i++) {
    l_nodes[i]->setRbt(l_keyframe[i]);
  }
  return;
}

void save_keyframes(const std::list<std::vector<RigTForm> >& l_keyFrames, const std::string& filename) {
  std::ofstream out(filename);
  Cvec3 temp_t;
  Quat temp_r;
  for(const auto& keyframe : l_keyFrames) {
    for(const auto& rigtform : keyframe) {
      temp_t = rigtform.getTranslation();
      temp_r = rigtform.getRotation();
      for(int i = 0; i < 3; i++) 
          out << temp_t[i] << ",";
      for(int i = 0; i < 4; i++) 
          out << temp_r[i] << (i < 3 ? "," : "\n");
    }
    out << "---\n";  // Delimiter between different keyframes
  }
}

std::list<std::vector<RigTForm> > load_keyframes(const std::string& filename) {
  std::list<std::vector<RigTForm> > l_keyFrames;
  std::ifstream in(filename);
  std::string line;
  Cvec3 temp_t;
  Quat temp_r;
  std::vector<RigTForm> currentKeyFrame;
  cout << filename << endl;
  while (std::getline(in, line)) {
    if (line == "---") {
      l_keyFrames.push_back(currentKeyFrame);
      currentKeyFrame.clear();
    } else {
      std::istringstream iss(line);
      RigTForm rigtform;
      for (int i = 0; i < 3; i++) {
          std::string value;
          std::getline(iss, value, ',');
          temp_t[i] = std::stod(value);
      }
      for (int i = 0; i < 4; i++) {
          std::string value;
          std::getline(iss, value, i < 3 ? ',' : '\n');
          temp_r[i] = std::stod(value);
      }
      rigtform.setRotation(temp_r);
      rigtform.setTranslation(temp_t);
      currentKeyFrame.push_back(rigtform);
    }
  }

  return l_keyFrames;
}

// Given t in the range [0, n], perform interpolation and draw the scene
// for the particular t. Returns true if we are at the end of the animation
// sequence, or false otherwise.
bool interpolateAndDisplay(float t) {
  if (t < g_keyFrames.size() - 3) {
    std::list<std::vector<RigTForm> >::iterator it = g_keyFrames.begin();
    for (int j = -1; j < std::floor(t); j++) {
      it++;
    } 
    it--;
    std::vector<RigTForm> base = *it;
    it++;
    std::vector<RigTForm> front = *it;
    it++;
    std::vector<RigTForm> back = *it;
    it++;
    std::vector<RigTForm> next = *it;

    std::vector<RigTForm> interp_keyFrame = std::vector<RigTForm>();
    double alpha = t - std::floor(t);
    assert(front.size() == back.size());
    RigTForm interpolated = RigTForm();
    for (int i = 0; i < front.size(); i++) {
      if (g_isLinear) {
        interpolated = Interpolation::Linear(front[i], back[i], alpha);
      } else {
        interpolated = Interpolation::CatmullRom(base[i], front[i], back[i], next[i], alpha);
      }
      interp_keyFrame.push_back(interpolated);
    }
    loadFromKeyFrame(interp_keyFrame);
    glutPostRedisplay();
    return false;
  } else {
    std::list<std::vector<RigTForm> >::iterator it = g_keyFrames.begin();
    std::advance(it, static_cast<size_t>(g_keyFrames.size()-2));
    std::vector<RigTForm> front = *it;
    loadFromKeyFrame(front);
    glutPostRedisplay();
    return true;
  }
}

// Interpret "ms" as milliseconds into the animation
static void animateTimerCallback(int ms) {
  if (!g_isAnimation) {
    return;
  }
  float t = (float)ms / (float)g_msBetweenKeyFrames;
  bool endReached = interpolateAndDisplay(t);
  // cout << endReached << endl;
  // cout << ms << endl;
  if (!endReached)
    glutTimerFunc(1000 / g_animateFramesPerSecond, animateTimerCallback,
                  ms + 1000 / g_animateFramesPerSecond);
  else {
    cout << "Animation playback finished" << endl;
    g_isAnimation = false;
  }
}

static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));

  makeCube(0.05, vtx.begin(), idx.begin());
  g_point.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initSphere() {
  int ibLen, vbLen;
  int radius = 1;
  int slices = 30;
  int stacks = 30;
  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  // Temporary storage for sphere geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeSphere(radius, slices, stacks, vtx.begin(), idx.begin());
  g_sphere.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initArcball() { 
  int ibLen, vbLen;
  int radius = 1;
  int slices = 10;
  int stacks = 10;
  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  // Temporary storage for arcball geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeSphere(radius, slices, stacks, vtx.begin(), idx.begin());
  g_arcball.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS380_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void drawStuff(const ShaderState& curSS, bool picking) {
  // short hand for current shader state
  // const ShaderState& curSS = *g_shaderStates[g_activeShader];

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // use the skyRbt as the eyeRbt
  const RigTForm eyeRbt = getPathAccumRbt(g_world, returnEyeNode(g_vcount));
  // printMatrix("eyeRbt", RigTFormToMatrix(eyeRbt));
  const RigTForm invEyeRbt = inv(eyeRbt);
  const Cvec3 eyeLight1 = Cvec3(RigTFormToMatrix(invEyeRbt) * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(RigTFormToMatrix(invEyeRbt) * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  const RigTForm groundRbt = RigTForm();  // identity
  Matrix4 MVM;
  Matrix4 NMVM;

  // ===========
  // draw ground
  //
  if (!picking) {
    Drawer drawer(invEyeRbt, curSS);
    g_world->accept(drawer);

    // draw arcball as part of asst3
    // draw sphere (arcball)
    g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);

    Matrix4 l_scale;
    double radius;

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // if manipulating a cube but is not the current eyeframe
    if (g_currentPickedRbtNode != shared_ptr<SgRbtNode>()) {
      MVM = RigTFormToMatrix(invEyeRbt * getPathAccumRbt(g_world, g_currentPickedRbtNode));
      radius = g_arcballScreenRadius * g_arcballScale;

      MVM *= l_scale.makeScale(radius);
      NMVM = normalMatrix(MVM);
      sendModelViewNormalMatrix(curSS, MVM, NMVM);
      safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[2]);
      g_arcball->draw(curSS);
      g_isArcball = true;
    } else {
      if (!(g_mouseMClickButton ||
            (g_mouseLClickButton && g_mouseRClickButton))) {
        MVM = RigTFormToMatrix(invEyeRbt * groundRbt);
        g_arcballScale =
            getScreenToEyeScale(-abs(norm(Cvec3(MVM[3], MVM[7], MVM[11]))),
                                g_frustFovY, g_windowHeight);
      }
      if (g_isArcball) {
        MVM = RigTFormToMatrix(invEyeRbt * groundRbt);
        radius = g_arcballScreenRadius * g_arcballScale;

        MVM *= l_scale.makeScale(radius);
        NMVM = normalMatrix(MVM);
        sendModelViewNormalMatrix(curSS, MVM, NMVM);
        safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1],
                         g_arcballColor[2]);
        g_arcball->draw(curSS);
      }
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else {
    cout << "Picked..." << endl;
    Picker picker(invEyeRbt, curSS);
    g_world->accept(picker);
    glFlush();
    g_currentPickedRbtNode =
        picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
    // cout << "(" << g_mouseClickX << ", " <<  g_mouseClickY << ")" << endl;
    if (g_currentPickedRbtNode == g_groundNode || g_currentPickedRbtNode == shared_ptr<SgRbtNode>()){
      g_currentPickedRbtNode = shared_ptr<SgRbtNode>();  // set to NULL
      g_isArcball = false;
    }
  }
}

static void pick() {
  // We need to set the clear color to black, for pick rendering.
  // so let's save the clear color
  GLdouble clearColor[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);

  glClearColor(0, 0, 0, 0);

  // using PICKING_SHADER as the shader
  glUseProgram(g_shaderStates[PICKING_SHADER]->program);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff(*g_shaderStates[PICKING_SHADER], true);

  // Uncomment below and comment out the glutPostRedisplay in mouse(...) call
  // back to see result of the pick rendering pass 
  // glutSwapBuffers();

  // Now set back the clear color
  glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);

  checkGlErrors();
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff(*g_shaderStates[g_activeShader], false);

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}

static void motion(const int x, const int y) {
  const double dx = x - g_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY;

  // coordinates of p1(x1, y1), p2(x2, y2)
  const double x1 = g_mouseClickX;
  const double y1 = g_mouseClickY;
  const double x2 = x;
  const double y2 = g_windowHeight - y - 1;

  // d1 and d2
  RigTForm eyeRbt = getPathAccumRbt(g_world, returnEyeNode(g_vcount));
  RigTForm invEyeRbt = inv(eyeRbt);

  // printMatrix("eyeRbt", RigTFormToMatrix(eyeRbt));

  Cvec3 d1 = getModelViewRay(Cvec2(x1, y1), g_frustFovY, g_windowWidth, g_windowHeight);
  Cvec3 d2 = getModelViewRay(Cvec2(x2, y2), g_frustFovY, g_windowWidth, g_windowHeight);

  // origin of sphere
  const double radius = g_arcballScale * g_arcballScreenRadius;
  const RigTForm groundRbt = RigTForm(); // identity

  RigTForm m;
  Cvec3 org = Cvec3();

  if (g_isArcball && g_currentPickedRbtNode != shared_ptr<SgRbtNode>()) {
    org = getPathAccumRbt(g_world, g_currentPickedRbtNode).getTranslation() - eyeRbt.getTranslation();
  } else {
    if (g_mcount % 2 == 1) {
      org = Cvec3() - eyeRbt.getTranslation();
    } else {
      org = Cvec3();
    }
  }

  org = invEyeRbt.getRotation() * org;

  if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
    
    // discriminant 
    const double dis1 = dot(org, d1) * dot(org, d1) - norm2(d1) * (norm2(org) - radius * radius);
    const double dis2 = dot(org, d2) * dot(org, d2) - norm2(d2) * (norm2(org) - radius * radius);

    Cvec3 v1;
    Cvec3 v2;
    double temp;
    double temp1;
    double temp2;

    if (dis1 < 0) {
      // no intersection
      v1 = normalize(d1 * dot(org - Cvec3(0, 0, 0), d1) - org) * radius;
    } else {
      // one or two intersection
      temp1 = (-dot(org, d1) + sqrt(dis1)) / norm2(d1);
      temp2 = (-dot(org, d1) - sqrt(dis1)) / norm2(d1);

      temp = fabs(temp1) < fabs(temp2) ? temp1 : temp2;
      v1 = - org - d1 * temp;
    }

    if (dis2 < 0) {
      // no intersection
      v2 = normalize(d2 * dot(org - Cvec3(0, 0, 0), d2) - org) * radius;
    } else {
      // one or two intersection
      temp1 = (-dot(org, d2) + sqrt(dis2)) / norm2(d2);
      temp2 = (-dot(org, d2) - sqrt(dis2)) / norm2(d2);

      temp = fabs(temp1) < fabs(temp2) ? temp1 : temp2;
      v2 = - org - d2 * temp;
    }

    g_objectRbt[2] = RigTForm(eyeRbt.getRotation() * (org + v1) + eyeRbt.getTranslation());
    g_objectRbt[3] = RigTForm(eyeRbt.getRotation() * (org + v2) + eyeRbt.getTranslation());

    m = RigTForm(Quat(0, normalize(v2)) * Quat(0, normalize(-v1)));                                                 
  } else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    if (g_vcount % 3 == 0 && g_mcount % 2 == 1) {
      // perform ego motion - arcball not in use
      m = RigTForm(Cvec3(dx, dy, 0) * 0.01);
    } else {
      m = RigTForm(Cvec3(dx, dy, 0) * g_arcballScale);
    }
  } else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    if (g_vcount % 3 == 0 && g_mcount % 2 == 1) {
      // perform ego motion - arcball not in use
      m = RigTForm(Cvec3(dx, dy, 0) * 0.01);
    } else {
      m = RigTForm(Cvec3(0, 0, -dy) * g_arcballScale);
    }
  }

  if (g_mouseClickDown) {
    RigTForm eyeRbt = getPathAccumRbt(g_world, returnEyeNode(g_vcount));
    RigTForm a;
    RigTForm temp;
  
    // cout << g_isArcball << endl;
    if (g_isArcball) {
      if (g_vcount % 3 == 0 && g_mcount % 2 == 1) {
        // perform ego motion
        RigTForm l_orgRbt = RigTForm(Cvec3(0.0, 0.0, 0.0));
        a = makeMixedFrame(l_orgRbt, eyeRbt); 
        temp = doMtoOwrtA(inv(m), eyeRbt, a);
        returnEyeNode(g_vcount)->setRbt(temp);
      } else {  
        // move the object
        a = makeMixedFrame(getPathAccumRbt(g_world, g_currentPickedRbtNode), eyeRbt);
        temp = doMtoOwrtA(m, g_currentPickedRbtNode->getRbt(), inv(getPathAccumRbt(g_world, g_currentPickedRbtNode, 1)) * a);
        g_currentPickedRbtNode->setRbt(temp);
      }
    } else {
      // move the eye
      temp = doMtoOwrtA(m, eyeRbt, eyeRbt);
      returnEyeNode(g_vcount)->setRbt(temp);
    }
    glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;

  if(g_isPicking) {
    pick();
    g_mcount = 0;
    g_isPicking = false;
  }
  glutPostRedisplay();
}

static void constructPlane(shared_ptr<SgTransformNode> base,
                           const Cvec3& color) {
  const double WING_LEN = 1.5, WING_THICK = 0.2, WING_WIDTH = 0.5;
  const double TAIL_LEN = 0.8, TAIL_THICK = 0.2, TAIL_WIDTH = 0.2;
  const double BODY_LEN = 1.5, BODY_THICK = 0.2, BODY_WIDTH = 0.3;
  const int NUM_JOINTS = 3, NUM_SHAPES = 3;

  struct JointDesc {
    int parent;
    float x, y, z;
  };

  JointDesc jointDesc[NUM_JOINTS] = {
      {-1},                                       // body
      {0, 0, 0.3, 0},   // wing
      {0, 0, -(BODY_LEN / 2), 0},                         // tail
  };

  struct ShapeDesc {
    int parentJointId;
    float x, y, z, sx, sy, sz;
    shared_ptr<Geometry> geometry;
  };

  ShapeDesc shapeDesc[NUM_SHAPES] = {
      {0, 0, 0, 0, BODY_WIDTH, BODY_LEN, BODY_THICK, g_cube},  // torso
      {1, 0, 0, 0.1, WING_LEN, WING_WIDTH, WING_THICK,
       g_cube},  // wing
      {2, 0, 0, 0, TAIL_LEN, TAIL_WIDTH, TAIL_THICK,
       g_cube},  // tail
  };

  shared_ptr<SgTransformNode> jointNodes[NUM_JOINTS];

  for (int i = 0; i < NUM_JOINTS; ++i) {
    if (jointDesc[i].parent == -1)
      jointNodes[i] = base;
    else {
      jointNodes[i].reset(new SgRbtNode(
          RigTForm(Cvec3(jointDesc[i].x, jointDesc[i].y, jointDesc[i].z))));
      jointNodes[jointDesc[i].parent]->addChild(jointNodes[i]);
    }
  }
  for (int i = 0; i < NUM_SHAPES; ++i) {
    shared_ptr<MyShapeNode> shape(new MyShapeNode(
        shapeDesc[i].geometry, color,
        Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z), Cvec3(0, 0, 0),
        Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
    jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
  }
}

static void initScene() {
  g_world.reset(new SgRootNode());

  g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 4.0))));

  g_groundNode.reset(new SgRbtNode());
  g_groundNode->addChild(shared_ptr<MyShapeNode>(
      new MyShapeNode(g_ground, Cvec3(0.1, 0.95, 0.1))));

  g_planeNode.reset(new SgRbtNode(RigTForm(Cvec3(-1.5, 1, 0))));

  constructPlane(g_planeNode, Cvec3(1, 0, 0));  // plane

  g_world->addChild(g_skyNode);
  g_world->addChild(g_groundNode);
  g_world->addChild(g_planeNode);
}

static void keyboard(const unsigned char key, const int x, const int y) {
  std::vector<RigTForm> temp = std::vector<RigTForm>();
  std::list<std::vector<RigTForm> >::iterator temp_iter = g_keyFramesIter;
  std::vector<RigTForm> curr = std::vector<RigTForm>();
  switch (key) {
  case 27:  // ESC
    exit(0);                                  
    cout << "Exit by ESC key" << endl;
  case 'h':
    cout << " ============== H E L P ==============\n\n"
         << "Esc\t\t Exit\n"
         << "h\t\tHelp menu\n"
         << "s\t\tSave screenshot\n"
         << "f\t\tToggle flat shading on/off.\n"
         << "v\t\tCycle view\n"
         << "m\t\tCycle the manipulated frame (only if both is sky camera)\n"
         << "p\t\tPart Picking\n"
         << "r\t\tReset as initial state\n"
         << "========== About Animation ==========\n\n"
         << "Space\t Display current keyframe\n"
         << "u\t\t Update current keyframe\n"
         << ">\t\t Advance to next keyframe\n"
         << "<\t\t Retreat to previous keyframe\n"
         << "d\t\t Delete current keyframe\n"
         << "n\t\t Create new keyframe\n"
         << "i\t\t Import(read) a keyframe list file\n"
         << "w\t\t Export(write) a keyframe list file\n"
         << "y\t\t Start/stop animation\n"
         << "+\t\t Faster animation\n"
         << "-\t\t Slower animation\n"
         << "l\t\t Switch between linear and Catmull-Rom\n"
         << "drag left mouse to rotate\n"
         << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    cout << "Screenshot saved" <<endl;
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  case 'v':
    g_vcount += 1;
    printView(g_vcount);
    g_currentPickedRbtNode = shared_ptr<SgRbtNode>();
    break;
  case 'm':
    if (g_vcount % 3 == 0) {
      if (g_mcount == 1) {
        cout << "Sky-sky frame" << endl;
        g_isArcball = false;
        g_mcount = 0;
      } else {
        cout << "World-sky frame" << endl;
        g_isArcball = true;
        g_mcount = 1;
      }
    } else {
      cout << "m is only for when view and object is all sky-camera" << endl;
    }
    g_currentPickedRbtNode = shared_ptr<SgRbtNode>(); 
    break;
  case 'p':
    g_isPicking = true;
    cout << "Enabling Picking..." << endl;
    break;
  case 'r':
    g_isArcball = false;
    g_isPicking = false;
    g_mcount = 0;
    g_vcount = 0;
    initScene();
    cout << "Reset" << endl;
    break;
  case 32:  //Spacebar
    if(g_keyFrames.size() == 0){
      cout << "Keyframes are empty, Add new keyframe to display" <<endl;
    } else {
      loadFromKeyFrame(*g_keyFramesIter);
      cout << "Display current keyframe" << endl;
      printCurrentKeyFrame();
    }
    break;
  case 'u':
    cout << "Update current keyframe" << endl;
    temp = saveToKeyFrame();
    if (g_keyFrames.empty()) {
      g_keyFrames.insert(g_keyFramesIter, temp);
      g_keyFramesIter--;
    } else {
      g_keyFramesIter = g_keyFrames.erase(g_keyFramesIter);
      g_keyFrames.insert(g_keyFramesIter, temp);
    }
    g_keyFramesIter--;
    printCurrentKeyFrame();
    break;
  case '>':
    if (g_keyFrames.empty()) {
      cout << "Keyframes are empty" << endl;
    } else {
      if (++temp_iter == g_keyFrames.end()) {
        cout << "Last keyframe" << endl;
      } else {
        cout << "Advance to next keyframe" << endl;
        g_keyFramesIter++;
      }
      loadFromKeyFrame(*g_keyFramesIter);
    printCurrentKeyFrame();
    }
    break;
  case '<':
    if (g_keyFrames.empty()) {
    cout << "Keyframes are empty" << endl;
    } else {
      if (g_keyFramesIter == g_keyFrames.begin()) {
        cout << "First keyframe" << endl;
      } else {
        cout << "Retreat to previous keyframe" << endl;
        g_keyFramesIter--;
      }
      loadFromKeyFrame(*g_keyFramesIter);
      printCurrentKeyFrame();
    }
    break;
  case 'd':
    cout << "Delete current keyframe" << endl;
    if (g_keyFrames.empty()) {
      cout << "Keyframe list is empty" << endl;
    } else {
      g_keyFramesIter = g_keyFrames.erase(g_keyFramesIter);
      if (distance(g_keyFrames.begin(), g_keyFramesIter) != 0) {
        g_keyFramesIter--;
      }
      if (!g_keyFrames.empty()) {
        loadFromKeyFrame(*g_keyFramesIter);
      }
    }
    printKeyFramesNum();
    printCurrentKeyFrame();
    break;
  case 'n':
    cout << "Create new keyframe" << endl;
    temp = saveToKeyFrame();
    g_keyFrames.insert(++temp_iter, temp);
    if (g_keyFrames.size() == 1) {
      g_keyFramesIter--;
    } else {
      g_keyFramesIter++;
    }
    printKeyFramesNum();
    printCurrentKeyFrame();
    break;
  case 'i':
    cout << "Import(read) a keyframe list file" << endl;
    g_keyFrames.clear();
    for (const auto& vector : load_keyframes(g_filename)) {
      g_keyFrames.push_back(vector);
    }
    g_keyFramesIter = g_keyFrames.begin();
    printKeyFramesNum();
    printCurrentKeyFrame();
    break;
  case 'w':
    cout << "Export(write) a keyframe list file" << endl;
    save_keyframes(g_keyFrames, g_filename);
    printKeyFramesNum();
    break;
  case 'y':
    if (!g_isAnimation) {
      if (g_keyFrames.size() < 4) {
        cout << "Need at least 4 keyframes to play animation" << endl;
      } else {
        cout << "Start animation" << endl;
        printInterpolationMode();
        cout << "g_msBetweenKeyFrames: " << g_msBetweenKeyFrames << endl;
        g_isAnimation = true;
        animateTimerCallback(0);
      }
    } else {
      cout << "Stop animation" << endl;
      g_isAnimation = false;
    }
    break;
  case '+':
    cout << "Increase animation speed" << endl;
    if (g_msBetweenKeyFrames == g_msOffset) {
      cout << "Maximum Speed Approached" << endl;
    } else {
      g_msBetweenKeyFrames -= g_msOffset;
    }
    cout << "g_msBetweenKeyFrames: " << g_msBetweenKeyFrames << endl;
    break;
  case '-':
    cout << "Decrease animation speed" << endl;
    g_msBetweenKeyFrames += g_msOffset;
    cout << "g_msBetweenKeyFrames: " << g_msBetweenKeyFrames << endl;
    break;
  case 'l':
    g_isLinear = !g_isLinear;
    printInterpolationMode();
  }
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("CS380: Assignment 4");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
      g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
    else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() {
  initGround();
  initCubes();
  initArcball();
  initSphere();
}

int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

    glewInit(); // load the OpenGL extensions

    cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3") << endl;
    if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");

    initGLState();
    initShaders();
    initGeometry();
    initScene();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}