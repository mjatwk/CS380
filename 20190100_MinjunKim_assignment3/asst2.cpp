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
static const bool g_Gl2Compatible = true;


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

struct ShaderState {
  GlProgram program;

  // Handles to uniform variables
  GLint h_uLight, h_uLight2;
  GLint h_uProjMatrix;
  GLint h_uModelViewMatrix;
  GLint h_uNormalMatrix;
  GLint h_uColor;

  // Handles to vertex attributes
  GLint h_aPosition;
  GLint h_aNormal;

  ShaderState(const char* vsfn, const char* fsfn) {
    readAndCompileShader(program, vsfn, fsfn);

    const GLuint h = program; // short hand

    // Retrieve handles to uniform variables
    h_uLight = safe_glGetUniformLocation(h, "uLight");
    h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
    h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
    h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
    h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");
    h_uColor = safe_glGetUniformLocation(h, "uColor");

    // Retrieve handles to vertex attributes
    h_aPosition = safe_glGetAttribLocation(h, "aPosition");
    h_aNormal = safe_glGetAttribLocation(h, "aNormal");

    if (!g_Gl2Compatible)
      glBindFragDataLocation(h, 0, "fragColor");
    checkGlErrors();
  }

};

static const int g_numShaders = 2;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"}
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


// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_arcball, g_point;

// --------- Scene

static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 4.0));
static RigTForm g_objectRbt[4] = {RigTForm(Cvec3(-1, 0, 0)), RigTForm(Cvec3(1, 0, 0)), RigTForm(Cvec3(0, 0, 0)), RigTForm(Cvec3(0, 0, 0))}; // currently only 4 objs are defined
static Cvec3f g_objectColors[2] = {Cvec3f(0, 0, 1), Cvec3f(1, 0, 0)};
static Cvec3f g_arcballColor = Cvec3f(0, 1, 0);

int g_vcount = 0;
int g_ocount = 1;
int g_mcount = 0;
static RigTForm g_orjSkyRbt = RigTForm(Cvec3(0.0, 0.25, 4.0));
static RigTForm g_orjObjectRbt[2] = {RigTForm(Cvec3(-1, 0, 0)), RigTForm(Cvec3(1, 0, 0))}; // currently only 2 objs are defined

double g_arcballScreenRadius;
double g_arcballScale;

static const bool isDebug = false;

///////////////// END OF G L O B A L S //////////////////////////////////////////////////

static RigTForm returnEyeRbt (int cnt) {
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
    throw std::runtime_error("returnEyeRbt; fail to handle all cases");
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

// takes MVM and its normal matrix to the shaders
static void sendModelViewNormalMatrix(const ShaderState& curSS, const Matrix4& MVM, const Matrix4& NMVM) {
  GLfloat glmatrix[16];
  MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
  safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

  NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
  safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);
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

static void drawStuff() {
  // short hand for current shader state
  const ShaderState& curSS = *g_shaderStates[g_activeShader];

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // use the skyRbt as the eyeRbt
  const RigTForm eyeRbt = returnEyeRbt(g_vcount);
  const RigTForm invEyeRbt = inv(eyeRbt);
  const Cvec3 eyeLight1 = Cvec3(rigTFormToMatrix(invEyeRbt) * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(rigTFormToMatrix(invEyeRbt) * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  // draw ground
  // ===========
  //
  const RigTForm groundRbt = RigTForm();  // identity
  Matrix4 MVM = rigTFormToMatrix(invEyeRbt * groundRbt);
  Matrix4 NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
  g_ground->draw(curSS);

  // draw cubes
  // ==========

  // draw red cube
  MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[0]);
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, g_objectColors[0][0], g_objectColors[0][1], g_objectColors[0][2]);
  g_cube->draw(curSS);

  // draw green cube
  MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[1]);
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, g_objectColors[1][0], g_objectColors[1][1], g_objectColors[1][2]);
  g_cube->draw(curSS);

  // draw sphere (arcball)
  g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);

  Matrix4 l_scale;
  double radius;

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  if (g_ocount % 3 != 0) {
    if (g_ocount % 3 != g_vcount % 3) {
      // if manipulating a cube but is not the current eyeframe
      MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[g_ocount % 3 - 1]);

      if (!(g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton))) {
        g_arcballScale = getScreenToEyeScale(-abs(norm(Cvec3(MVM[3], MVM[7], MVM[11]))), g_frustFovY, g_windowHeight);
      }
      radius = g_arcballScreenRadius * g_arcballScale;

      MVM *= l_scale.makeScale(radius);
      NMVM = normalMatrix(MVM);
      sendModelViewNormalMatrix(curSS, MVM, NMVM);
      safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[2]);
      g_arcball->draw(curSS);
    }
  } else {
    if (g_mcount % 2 == 1) {
      // if manipulating sky camera w.r.t the world-sky coordinate
      MVM = rigTFormToMatrix(invEyeRbt * groundRbt);

      if (!(g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton))) {
        g_arcballScale = getScreenToEyeScale(-abs(norm(Cvec3(MVM[3], MVM[7], MVM[11]))), g_frustFovY, g_windowHeight);
      }
      radius = g_arcballScreenRadius * g_arcballScale;

      MVM *= l_scale.makeScale(radius);
      NMVM = normalMatrix(MVM);
      sendModelViewNormalMatrix(curSS, MVM, NMVM);
      safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[2]);
      g_arcball->draw(curSS);
    }
  }
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // draw points
  if (isDebug) {
    MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[2]);
    NMVM = normalMatrix(MVM);
    sendModelViewNormalMatrix(curSS, MVM, NMVM);
    safe_glUniform3f(curSS.h_uColor, g_objectColors[0][0], g_objectColors[0][1], g_objectColors[0][2]);
    g_point->draw(curSS);

    MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[3]);
    NMVM = normalMatrix(MVM);
    sendModelViewNormalMatrix(curSS, MVM, NMVM);
    safe_glUniform3f(curSS.h_uColor, g_objectColors[1][0], g_objectColors[1][1], g_objectColors[1][2]);
    g_point->draw(curSS);
  }
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff();

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
  RigTForm eyeRbt = returnEyeRbt(g_vcount);
  RigTForm invEyeRbt = inv(returnEyeRbt(g_vcount));

  Cvec3 d1 = getModelViewRay(Cvec2(x1, y1), g_frustFovY, g_windowWidth, g_windowHeight);
  Cvec3 d2 = getModelViewRay(Cvec2(x2, y2), g_frustFovY, g_windowWidth, g_windowHeight);

  // origin of sphere
  const double radius = g_arcballScale * g_arcballScreenRadius;
  const RigTForm groundRbt = RigTForm(); // identity

  RigTForm m;
  Cvec3 org = Cvec3();

  // if the sphere exists
  if (g_ocount % 3 != 0)
  {
    if (g_ocount % 3 != g_vcount % 3)
    {
      // if manipulating a cube but is not the current eyeframe
      org = returnManipulateRbt(g_ocount).getTranslation() - returnEyeRbt(g_vcount).getTranslation();
    }
  }
  else
  {
    if (g_mcount % 2 == 1)
    {
      // if manipulating sky camera w.r.t the world-sky coordinate
      org = groundRbt.getTranslation() - returnEyeRbt(g_vcount).getTranslation();
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

    g_objectRbt[2] = RigTForm(eyeRbt.getRotation() * (org + v1) + returnEyeRbt(g_vcount).getTranslation());
    g_objectRbt[3] = RigTForm(eyeRbt.getRotation() * (org + v2) + returnEyeRbt(g_vcount).getTranslation());

    m = RigTForm(Quat(0, normalize(v2)) * Quat(0, normalize(-v1)));                                                             
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    m = RigTForm(Cvec3(dx, dy, 0) * g_arcballScale);
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    m = RigTForm(Cvec3(0, 0, -dy) * g_arcballScale);
  }

  if (g_mouseClickDown) {
    RigTForm& l_objRbt = returnManipulateRbt(g_ocount);
    RigTForm& l_eyeRbt = returnManipulateRbt(g_vcount); 
    RigTForm a = makeMixedFrame(l_objRbt, l_eyeRbt);
    RigTForm temp;

    if (g_ocount % 3 == g_vcount % 3) {
      if (g_ocount % 3 == 0 && g_mcount % 2 == 1) {
        // perform ego motion
        RigTForm l_orgRbt = RigTForm(Cvec3(0.0, 0.0, 0.0));
        a = makeMixedFrame(l_orgRbt, l_eyeRbt); 
        temp = doMtoOwrtA(inv(m), l_eyeRbt, a);
        softBind_rigtform(temp, l_objRbt);
      } else {  
        // move the eye
        temp = doMtoOwrtA(m, l_objRbt, a);
        softBind_rigtform(temp, l_objRbt);
      }
    } else {
      // move the object
      temp = doMtoOwrtA(m, l_objRbt, a);
      softBind_rigtform(temp, l_objRbt);
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

  glutPostRedisplay();
}


static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "o\t\tCycle object to edit\n"
    << "v\t\tCycle view\n"
    << "m\t\tCycle the manipulated frame (only if both is sky camera)\n"
    << "r\t\tReset as initial state\n"
    << "drag left mouse to rotate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  case 'o':
    switch (g_vcount % 3) {
    case 0:
      g_ocount += 1;
      break;
    case 1:
    case 2:
      if (g_ocount % 3 == 0) {
        throw runtime_error("Error: pressed o button, cannot control sky while cube view");
      }  
      g_ocount += (g_ocount % 3 == 1) ? 1 : 2;
      break;

    default:
      throw runtime_error("Error: pressed o button, the switch-case statement doesnt handle all cases");
      break;
    }
    printObject(g_ocount);
    break;
  case 'v':
    g_vcount += 1;
    if (g_vcount % 3 == 1 && g_ocount % 3 == 0) {
      g_ocount = 1;
    } 
    printView(g_vcount);
    break;
  case 'm':
    if (g_ocount % 3 == 0 && g_vcount % 3 == 0) {
      g_mcount += 1;
      if (g_mcount % 2 == 0) {
        cout << "Sky-sky frame" << endl;
      } else {
        cout << "World-sky frame" << endl;
      }
    } else {
      cout << "m is only for when view and object is all sky-camera" << endl;
    }
    break;
  case 'r':
    g_skyRbt = g_orjSkyRbt;
    memcpy(&g_objectRbt,&g_orjObjectRbt, sizeof(RigTForm)*2);
    cout << "Reset to initial position" << endl;
    break;
  }
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("CS380: Assignment 2");                       // title the window

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

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
