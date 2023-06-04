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

#include "geometry.h"

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

// --------- Materials
// This should replace all the contents in the Shaders section, e.g.,
// g_numShaders, g_shaderFiles, and so on
static shared_ptr<Material> g_redDiffuseMat, g_blueDiffuseMat, g_bumpFloorMat,
    g_arcballMat, g_pickingMat, g_lightMat;
shared_ptr<Material> g_shadowMapMat, g_shadowMat;

shared_ptr<Material> g_overridingMaterial;


// --------- Geometry

typedef SgGeometryShapeNode MyShapeNode;

// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_arcball, g_point, g_sphere;

static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node;
static shared_ptr<SgRbtNode> g_light1Node, g_light2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do picking

// --------- Scene

static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 7.0));
static RigTForm g_objectRbt[4] = {RigTForm(Cvec3(-1, 0, 0)), RigTForm(Cvec3(1, 0, 0)), RigTForm(Cvec3(0, 0, 0)), RigTForm(Cvec3(0, 0, 0))}; // currently only 4 objs are defined

int g_vcount = 0;
int g_mcount = 0;
static RigTForm g_orjSkyRbt = RigTForm(Cvec3(0.0, 0.25, 8.0));
static RigTForm g_orjObjectRbt[2] = {RigTForm(Cvec3(-1, 0, 0)), RigTForm(Cvec3(1, 0, 0))}; // currently only 2 objs are defined

double g_arcballScreenRadius;
double g_arcballScale;

bool g_isPicking = false;
bool g_isArcball = false;

static const bool isDebug = false;

// --------- Shadow

GLuint loadShaders(const char* vertex_file_path,
                   const char* fragment_file_path) {
  // Create the shaders
  GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
  GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

  // Read the Vertex Shader code from the file
  std::string VertexShaderCode;
  FILE* vertexShaderFile = fopen(vertex_file_path, "r");
  if (vertexShaderFile != NULL) {
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), vertexShaderFile))
      VertexShaderCode += buffer;
    fclose(vertexShaderFile);
  }

  // Read the Fragment Shader code from the file
  std::string FragmentShaderCode;
  FILE* fragmentShaderFile = fopen(fragment_file_path, "r");
  if (fragmentShaderFile != NULL) {
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), fragmentShaderFile))
      FragmentShaderCode += buffer;
    fclose(fragmentShaderFile);
  }

  GLint Result = GL_FALSE;
  int InfoLogLength;

  // Compile Vertex Shader
  printf("Compiling shader : %s\n", vertex_file_path);
  char const* VertexSourcePointer = VertexShaderCode.c_str();
  glShaderSource(VertexShaderID, 1, &VertexSourcePointer, NULL);
  glCompileShader(VertexShaderID);

  // Check Vertex Shader
  glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
  glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
  if (InfoLogLength > 0) {
    std::vector<char> VertexShaderErrorMessage(InfoLogLength + 1);
    glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL,
                       &VertexShaderErrorMessage[0]);
    printf("%s\n", &VertexShaderErrorMessage[0]);
  }

  // Compile Fragment Shader
  printf("Compiling shader : %s\n", fragment_file_path);
  char const* FragmentSourcePointer = FragmentShaderCode.c_str();
  glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer, NULL);
  glCompileShader(FragmentShaderID);

  // Check Fragment Shader
  glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
  glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
  if (InfoLogLength > 0) {
    std::vector<char> FragmentShaderErrorMessage(InfoLogLength + 1);
    glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL,
                       &FragmentShaderErrorMessage[0]);
    printf("%s\n", &FragmentShaderErrorMessage[0]);
  }

  // Link the program
  printf("Linking program\n");
  GLuint ProgramID = glCreateProgram();
  glAttachShader(ProgramID, VertexShaderID);
  glAttachShader(ProgramID, FragmentShaderID);
  glLinkProgram(ProgramID);

  // Check the program
  glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
  glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
  if (InfoLogLength > 0) {
    std::vector<char> ProgramErrorMessage(InfoLogLength + 1);
    glGetProgramInfoLog(ProgramID, InfoLogLength, NULL,
                        &ProgramErrorMessage[0]);
    printf("%s\n", &ProgramErrorMessage[0]);
  }

  glDeleteShader(VertexShaderID);
  glDeleteShader(FragmentShaderID);

  return ProgramID;
}

int depthMapFBO;
int g_frameBufferWidth = 512;
int g_frameBufferHeight = 512;
GLuint g_framebufferName = 0;
GLuint g_depthProgramID;
GLuint g_depthMatrixID;
GLuint g_programID;
GLuint g_matrixID;
GLuint g_depthBiasID;
GLuint g_shadowMapID;
GLuint g_vertexPosition_modelspaceID;
GLuint g_vertexPosition_modelspace2ID;
GLuint g_vertexBuffer;

///////////////// END OF G L O B A L S //////////////////////////////////////////////////

static shared_ptr<SgRbtNode> returnEyeNode (int cnt) {
  switch (cnt % 3){
  case 0:
    return g_skyNode;
    break;
  
  case 1:
    return g_robot1Node;
    break;

  case 2:
    return g_robot2Node;
    break;

  default:
    throw std::runtime_error("returnEyeNode; fail to handle all cases");
    break;
  }
}

static RigTForm& returnManipulateRbt(int cnt){
  switch (cnt % 3){
  case 0:
    return g_skyRbt;
    break;

  case 1:
    return g_objectRbt[0];
    break;

  case 2:
    return g_objectRbt[1];
    break;

  default:
    throw std::runtime_error("returnManipulateRbt; fail to handle all cases");
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
  switch (n % 3)
  {
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

void print(int i) {
  cout << "check: " << to_string(i) << endl; 
  return;
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


static void initGround() {
  int ibLen, vbLen;
  getPlaneVbIbLen(vbLen, ibLen);

  // Temporary storage for cube Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makePlane(g_groundSize * 2, vtx.begin(), idx.begin());
  g_ground.reset(
      new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
}

GLfloat& toFloat(Matrix4 abc) {
  GLfloat hi[16];
  for (int i = 0; i < 16; ++i) {
    hi[i] = static_cast<GLfloat>(abc[i]);
  }
  return hi[0];
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initSphere() {
  int ibLen, vbLen;
  getSphereVbIbLen(20, 10, vbLen, ibLen);

  // Temporary storage for sphere Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  makeSphere(1, 20, 10, vtx.begin(), idx.begin());
  g_sphere.reset(
      new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vtx.size(), idx.size()));
}

// takes a projection matrix and send to the the shaders
inline void sendProjectionMatrix(Uniforms& uniforms, const Matrix4& projMatrix) {
  uniforms.put("uProjMatrix", projMatrix);
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

static void drawStuff(bool picking) {
  // short hand for current shader state
  // const ShaderState& curSS = *g_shaderStates[g_activeShader];

  Uniforms uniforms;

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(uniforms, projmat);

  // use the skyRbt as the eyeRbt
  const RigTForm eyeRbt = getPathAccumRbt(g_world, returnEyeNode(g_vcount));
  const RigTForm invEyeRbt = inv(eyeRbt);

  Cvec3 eyeLight1 = getPathAccumRbt(g_world, g_light1Node).getTranslation();
  Cvec3 eyeLight2 = getPathAccumRbt(g_world, g_light2Node).getTranslation();

  uniforms.put("uLight", Cvec3(invEyeRbt * Cvec4(eyeLight1, 1)));
  uniforms.put("uLight2", Cvec3(invEyeRbt * Cvec4(eyeLight2, 1)));

  const RigTForm groundRbt = RigTForm();  // identity
  Matrix4 MVM;
  Matrix4 NMVM;

  // ===========
  // draw ground
  //
  if (!picking) {
    Drawer drawer(invEyeRbt, uniforms);
    g_world->accept(drawer);

    // draw arcball as part of asst3
    // draw sphere (arcball)
    g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);

    Matrix4 l_scale;
    double radius;

    // if manipulating a cube but is not the current eyeframe
    if (g_currentPickedRbtNode != shared_ptr<SgRbtNode>()) {
      MVM = RigTFormToMatrix(invEyeRbt * getPathAccumRbt(g_world, g_currentPickedRbtNode));
      radius = g_arcballScreenRadius * g_arcballScale;

      MVM *= l_scale.makeScale(radius);
      NMVM = normalMatrix(MVM);
      sendModelViewNormalMatrix(uniforms, MVM, NMVM);
      // safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[2]);
      g_arcballMat->draw(*g_sphere, uniforms);
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
        sendModelViewNormalMatrix(uniforms, MVM, NMVM);
        // safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[2]);
        g_arcballMat->draw(*g_sphere, uniforms);
      }
    }
  } else {
    cout << "Picked..." << endl;
    Picker picker(invEyeRbt, uniforms);
    g_overridingMaterial = g_pickingMat;
    g_world->accept(picker);
    g_overridingMaterial.reset();
    glFlush();
    g_currentPickedRbtNode =
        picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
    // cout << "(" << g_mouseClickX << ", " <<  g_mouseClickY << ")" << endl;
    if (g_currentPickedRbtNode == g_groundNode || g_currentPickedRbtNode == shared_ptr<SgRbtNode>()){
      g_currentPickedRbtNode = shared_ptr<SgRbtNode>();  // set to NULL
      g_isArcball = false;
    }
  }


  // // draw shadow
 
  // // Assuming we have the projection, view and model matrices as:
  // Matrix4 projectionMatrix =
  //     makeProjectionMatrix();  // Projection matrix is already provided by
  //                              // `makeProjectionMatrix` method
  // Matrix4 viewMatrix =
  //     inv(RigTFormToMatrix(eyeRbt));  // View matrix from eyeRbt
  // Matrix4 modelMatrix =
  //     RigTFormToMatrix(groundRbt);  // Model matrix from groundRbt

  // // Construct the MVP matrix by multiplication
  // Matrix4 mvpMatrix = projectionMatrix * viewMatrix * modelMatrix;

  // // Construct the depthMVP matrix for shadows. This involves projection, view,
  // // and model matrices for the light source
  // Matrix4 depthProjectionMatrix = Matrix4::makeProjection(
  //     g_frustFovY, g_windowWidth / static_cast<double>(g_windowHeight),
  //     g_frustNear, g_frustFar);  // Adjust as per light source
  // Matrix4 depthViewMatrix = inv(RigTFormToMatrix(getPathAccumRbt(
  //     g_world, g_light1Node)));  // Assuming light1Node is the light source
  // Matrix4 depthModelMatrix = modelMatrix;  // Shadow of the same object

  // Matrix4 depthMVP = depthProjectionMatrix * depthViewMatrix * depthModelMatrix;
  
  // Uniforms extraUniforms;
  // // Set vertexPosition_modelspace attribute
  // glBindBuffer(GL_ARRAY_BUFFER,
  //              g_vertexBuffer);  // Assuming vertexBuffer is your vertex buffer
  // glEnableVertexAttribArray(0);
  // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  // extraUniforms.put("depthMVP", depthMVP);

  // g_shadowMapMat->draw(*g_sphere, extraUniforms);
}

static void pick() {
  // We need to set the clear color to black, for pick rendering.
  // so let's save the clear color
  GLdouble clearColor[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);

  glClearColor(0, 0, 0, 0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff(true);

  // Uncomment below and comment out the glutPostRedisplay in mouse(...) call
  // back to see result of the pick rendering pass 
  // glutSwapBuffers();

  // Now set back the clear color
  glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);

  checkGlErrors();
}

static void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff(false);

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

static void constructRobot(shared_ptr<SgTransformNode> base, shared_ptr<Material> material) {
  const double ARM_LEN = 0.7, ARM_THICK = 0.25, TORSO_LEN = 1.5,
               TORSO_THICK = 0.25, TORSO_WIDTH = 1, HEAD_RADIUS = 0.5;
  const int NUM_JOINTS = 10, NUM_SHAPES = 10;

  struct JointDesc {
    int parent;
    float x, y, z;
  };

  JointDesc jointDesc[NUM_JOINTS] = {
      {-1},                                       // torso
      {0, (TORSO_WIDTH) / 2, TORSO_LEN / 2, 0},   // upper right arm
      {1, ARM_LEN, 0, 0},                         // lower right arm
      {0, -(TORSO_WIDTH) / 2, TORSO_LEN / 2, 0},  // upper left arm
      {3, -ARM_LEN, 0, 0},                        // lower left arm
      {0, (TORSO_WIDTH) / 4, -(TORSO_LEN / 2),
       0},                  // upper right leg
      {5, 0, -ARM_LEN, 0},  // lower right leg
      {0, -(TORSO_WIDTH) / 4, -(TORSO_LEN / 2),
       0},                                      // upper left leg
      {7, 0, -ARM_LEN, 0},                      // lower left leg
      {0, 0, (TORSO_LEN / 2), 0}  // head
  };

  struct ShapeDesc {
    int parentJointId;
    float x, y, z, sx, sy, sz;
    shared_ptr<Geometry> geometry;
  };

  ShapeDesc shapeDesc[NUM_SHAPES] = {
      {0, 0, 0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube},  // torso
      {1, ARM_LEN / 2, 0, 0, ARM_LEN / 2, ARM_THICK / 2, ARM_THICK / 2,
       g_sphere},  // upper right arm
      {2, ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK,
       g_cube},  // lower right arm
      {3, -ARM_LEN / 2, 0, 0, ARM_LEN / 2, ARM_THICK / 2, ARM_THICK / 2,
       g_sphere},  // upper left arm
      {4, -ARM_LEN / 2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK,
       g_cube},  // lower left arm
      {5, 0, -ARM_LEN / 2, 0, ARM_THICK / 2, ARM_LEN / 2, ARM_THICK / 2,
       g_sphere},  // upper right leg
      {6, 0, -ARM_LEN / 2, 0, ARM_THICK, ARM_LEN, ARM_THICK,
       g_cube},  // lower right leg
      {7, 0, -ARM_LEN / 2, 0, ARM_THICK / 2, ARM_LEN / 2, ARM_THICK / 2,
       g_sphere},  // upper left leg
      {8, 0, -ARM_LEN / 2, 0, ARM_THICK, ARM_LEN, ARM_THICK,
       g_cube},  // lower left leg
      {9, 0, HEAD_RADIUS, 0, HEAD_RADIUS, HEAD_RADIUS, HEAD_RADIUS, g_sphere},
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
  // The new MyShapeNode takes in a material as opposed to color
  for (int i = 0; i < NUM_SHAPES; ++i) {
    shared_ptr<SgGeometryShapeNode> shape(new MyShapeNode(
        shapeDesc[i].geometry, material,  // USE MATERIAL as opposed to color
        Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z), Cvec3(0, 0, 0),
        Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
    jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
  }
}

static void initScene() {
  g_world.reset(new SgRootNode());

  g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 7.0))));

  g_light1Node.reset(new SgRbtNode(RigTForm(Cvec3(0, 1.0, 8.0))));
  g_light2Node.reset(new SgRbtNode(RigTForm(Cvec3(0, 1.0, -6.0))));
  g_light1Node->addChild(shared_ptr<MyShapeNode>(
      new MyShapeNode(g_sphere, g_lightMat, Cvec3(0, 0, 0))));
  g_light2Node->addChild(shared_ptr<MyShapeNode>(
      new MyShapeNode(g_sphere, g_lightMat, Cvec3(0, 0, 0))));

  g_groundNode.reset(new SgRbtNode());
  g_groundNode->addChild(shared_ptr<MyShapeNode>(
      new MyShapeNode(g_ground, g_bumpFloorMat, Cvec3(0, g_groundY, 0))));

  g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-3, 1, 0))));
  g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(3, 1, 0))));

  constructRobot(g_robot1Node, g_redDiffuseMat);  // a Red robot
  constructRobot(g_robot2Node, g_blueDiffuseMat);  // a Blue robot

  g_world->addChild(g_skyNode);
  g_world->addChild(g_groundNode);
  g_world->addChild(g_light1Node);
  g_world->addChild(g_light2Node);
  g_world->addChild(g_robot1Node);
  g_world->addChild(g_robot2Node);
}

static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "v\t\tCycle view\n"
    << "m\t\tCycle the manipulated frame (only if both is sky camera)\n"
    << "p\t\tPart Picking\n"
    << "r\t\tReset as initial state\n"
    << "drag left mouse to rotate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
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

static void initMaterials() {
  // Create some prototype materials
  Material diffuse("./shaders/basic-gl2.vshader",
                   "./shaders/diffuse-gl2.fshader");
  Material solid("./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader");

  // copy diffuse prototype and set red color
  g_redDiffuseMat.reset(new Material(diffuse));
  g_redDiffuseMat->getUniforms().put("uColor", Cvec3f(1, 0, 0));

  // copy diffuse prototype and set blue color
  g_blueDiffuseMat.reset(new Material(diffuse));
  g_blueDiffuseMat->getUniforms().put("uColor", Cvec3f(0, 0, 1));

  // normal mapping material
  g_bumpFloorMat.reset(new Material("./shaders/normal-gl2.vshader",
                                    "./shaders/normal-gl2.fshader"));
  g_bumpFloorMat->getUniforms().put(
      "uTexColor",
      shared_ptr<ImageTexture>(new ImageTexture("Fieldstone.ppm", true)));
  g_bumpFloorMat->getUniforms().put(
      "uTexNormal", shared_ptr<ImageTexture>(
                        new ImageTexture("FieldstoneNormal.ppm", false)));
  // if the compile errors occur, use below instead (intel mac)
  // g_bumpFloorMat.reset(new Material("./shaders/normal-gl2.vshader",
  // "./shaders/normal-gl2.fshader"));
  // g_bumpFloorMat->getUniforms().put_tex("uTexColor",
  // shared_ptr<ImageTexture>(new ImageTexture("Fieldstone.ppm", true)));
  // g_bumpFloorMat->getUniforms().put_tex("uTexNormal",
  // shared_ptr<ImageTexture>(new ImageTexture("FieldstoneNormal.ppm", false)));

  // copy solid prototype, and set to wireframed rendering
  g_arcballMat.reset(new Material(solid));
  g_arcballMat->getUniforms().put("uColor", Cvec3f(0.27f, 0.82f, 0.35f));
  g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // copy solid prototype, and set to color white
  g_lightMat.reset(new Material(solid));
  g_lightMat->getUniforms().put("uColor", Cvec3f(1, 1, 1));

  // pick shader
  g_pickingMat.reset(new Material("./shaders/basic-gl2.vshader",
                                  "./shaders/pick-gl2.fshader"));

};


// static bool initShadow() {
//   GLuint depthTexture;

//   g_depthProgramID = loadShaders("./shaders/shadow-gl2.vshader",
//                                  "./shaders/shadow-gl2.fshader");
//   g_programID = loadShaders("./shaders/simple-gl2.vshader",
//                             "./shaders/simple-gl2.fshader");

//   g_depthMatrixID = glGetUniformLocation(g_depthProgramID, "depthMVP");
//   g_vertexPosition_modelspaceID =
//       glGetAttribLocation(g_depthProgramID, "vertexPosition_modelspace");
//   g_vertexPosition_modelspace2ID =
//       glGetAttribLocation(g_programID, "vertexPosition_modelspace");
//   g_matrixID = glGetUniformLocation(g_programID, "MVP");
//   g_depthBiasID = glGetUniformLocation(g_programID, "DepthBiasMVP");
//   g_shadowMapID = glGetUniformLocation(g_programID, "shadowMap");

//   glGenFramebuffers(1, &g_framebufferName);
//   glBindFramebuffer(GL_FRAMEBUFFER, g_framebufferName);

//   // Create the depth texture.
//   glGenTextures(1, &depthTexture);
//   glBindTexture(GL_TEXTURE_2D, depthTexture);
//   glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT16, g_windowWidth,
//                g_windowHeight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE,
//                   GL_COMPARE_R_TO_TEXTURE);

//   // Attach the depth texture to the Framebuffer Object.
//   glBindFramebuffer(GL_FRAMEBUFFER, g_framebufferName);
//   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
//                          depthTexture, 0);
//   glDrawBuffer(GL_NONE);  // No color buffer is drawn to.
//   glReadBuffer(GL_NONE);

//   // Always check that our framebuffer is ok
//   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
//     cout << "false" << endl;
//     return false;

//   return true;
// }

static void initGeometry() {
  initGround();
  initCubes();
  // initArcball();
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
    initMaterials();
    initGeometry();
    initScene();
    // initShadow();
    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
