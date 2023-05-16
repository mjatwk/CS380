#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS380_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r) {
    //TODO_T
    softBind_cvec3(t, t_);
    softBind_quat(r, r_);
  }

  explicit RigTForm(const Cvec3& t) {
    // TODO_T
    softBind_cvec3(t, t_);
  }

  explicit RigTForm(const Quat& r) {
    // TODO_T
    softBind_quat(r, r_);
  }

  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  Cvec4 operator * (const Cvec4& a) const;

  RigTForm operator * (const RigTForm& a) const;
};

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline Matrix4 RigTFormToMatrix(const RigTForm& tform) {
  // TODO_T
  Matrix4 m = quatToMatrix(tform.getRotation());

  m[3] = tform.getTranslation()[0];
  m[7] = tform.getTranslation()[1];
  m[11] = tform.getTranslation()[2];
  m[15] = 1;

  return m;
}

inline RigTForm inv(const RigTForm &tform)
{
  // TODO_T
  Matrix4 r_inv = inv(RigTFormToMatrix(RigTForm(tform.getRotation())));
  Cvec3 rt(r_inv[0] * tform.getTranslation()[0] + r_inv[1] * tform.getTranslation()[1] + r_inv[2] * tform.getTranslation()[2], r_inv[4] * tform.getTranslation()[0] + r_inv[5] * tform.getTranslation()[1] + r_inv[6] * tform.getTranslation()[2], r_inv[8] * tform.getTranslation()[0] + r_inv[9] * tform.getTranslation()[1] + r_inv[10] * tform.getTranslation()[2]);
  return RigTForm(rt * (-1), inv(tform.getRotation()));
}

inline Cvec4 RigTForm::operator*(const Cvec4 &a) const
{
  // TODO_T
  return RigTFormToMatrix(*this) * a;
}

inline RigTForm RigTForm::operator*(const RigTForm &a) const
{
  // TODO_T
  Matrix4 mat = RigTFormToMatrix(*this);
  Cvec3 rt(mat[0] * a.getTranslation()[0] + mat[1] * a.getTranslation()[1] + mat[2] * a.getTranslation()[2], mat[4] * a.getTranslation()[0] + mat[5] * a.getTranslation()[1] + mat[6] * a.getTranslation()[2], mat[8] * a.getTranslation()[0] + mat[9] * a.getTranslation()[1] + mat[10] * a.getTranslation()[2]);
  return RigTForm(t_ + rt, r_ * a.getRotation());
}

// -- CUSTOM FUNCTION FROM HERE -- //

// from ppt
// do M to O with respect to A
static RigTForm doMtoOwrtA(const RigTForm &M, const RigTForm &O, const RigTForm &A)
{
  return A * M * inv(A) * O;
}

// from ppt
// make affine matrix A
static RigTForm makeMixedFrame(const RigTForm &O, const RigTForm &E)
{
  return RigTForm(O.getTranslation(), E.getRotation());
}

// softbind function for RigTForm type
inline void softBind_rigtform(const RigTForm &from, RigTForm &to)
{
  to.setRotation(from.getRotation());
  to.setTranslation(from.getTranslation());
  return;
}

#endif
