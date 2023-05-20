#include <GL/glew.h>

#include "picker.h"

using namespace std;

Picker::Picker(const RigTForm& initialRbt, const ShaderState& curSS)
  : drawer_(initialRbt, curSS)
  , idCounter_(0)
  , srgbFrameBuffer_(!g_Gl2Compatible) {}

bool Picker::visit(SgTransformNode& node) {
  // TODO_T
  drawer_.visit(node);
  nodeStack_.push_back(node.shared_from_this());
  return true;
}

bool Picker::postVisit(SgTransformNode& node) {
  // TODO_T
  drawer_.postVisit(node);
  nodeStack_.pop_back();
  return true;
}   

bool Picker::visit(SgShapeNode& node) {
  // TODO_T
  addToMap(++idCounter_, dynamic_pointer_cast<SgRbtNode>(nodeStack_.back()));
  safe_glUniform3f(drawer_.getCurSS().h_uIdColor, idToColor(idCounter_)[0], idToColor(idCounter_)[1], idToColor(idCounter_)[2]);
  drawer_.visit(node);
  return true;
}

bool Picker::postVisit(SgShapeNode& node) {
  // TODO_T
  drawer_.postVisit(node);
  return true;
}

shared_ptr<SgRbtNode> Picker::getRbtNodeAtXY(int x, int y) {
  // TODO_T
  PackedPixel color;
  glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &color);
  return find(colorToId(color));
}

//------------------
// Helper functions
//------------------
//
void Picker::addToMap(int id, shared_ptr<SgRbtNode> node) {
  idToRbtNode_[id] = node;
}

shared_ptr<SgRbtNode> Picker::find(int id) {
  IdToRbtNodeMap::iterator it = idToRbtNode_.find(id);
  if (it != idToRbtNode_.end())
    return it->second;
  else
    return shared_ptr<SgRbtNode>(); // set to null
}

// encode 2^4 = 16 IDs in each of R, G, B channel, for a total of 16^3 number of objects
static const int NBITS = 4, N = 1 << NBITS, MASK = N-1;

Cvec3 Picker::idToColor(int id) {
  assert(id > 0 && id < N * N * N);
  Cvec3 framebufferColor = Cvec3(id & MASK, (id >> NBITS) & MASK, (id >> (NBITS+NBITS)) & MASK);
  framebufferColor = framebufferColor / N + Cvec3(0.5/N);

  if (!srgbFrameBuffer_)
    return framebufferColor;
  else {
    // if GL3 is used, the framebuffer will be in SRGB format, and the color we supply needs to be in linear space
    Cvec3 linearColor;
    for (int i = 0; i < 3; ++i) {
      linearColor[i] = framebufferColor[i] <= 0.04045 ? framebufferColor[i]/12.92 : pow((framebufferColor[i] + 0.055)/1.055, 2.4);
    }
    return linearColor;
  }
}

int Picker::colorToId(const PackedPixel& p) {
  const int UNUSED_BITS = 8 - NBITS;
  int id = p.r >> UNUSED_BITS;
  id |= ((p.g >> UNUSED_BITS) << NBITS);
  id |= ((p.b >> UNUSED_BITS) << (NBITS+NBITS));
  return id;
}