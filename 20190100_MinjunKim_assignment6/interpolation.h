#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <cmath>
#include <iostream>
#include "rigtform.h"

namespace Interpolation {

// Linear interpolation of two coordinate vectors
// TODO: compute the linear interpolation of two Cvec3 inputs and return the
// result
inline Cvec3 lerp(const Cvec3& c0, const Cvec3& c1, const double& alpha) {
                Cvec3 result = Cvec3();
                result = c0 * (1 - alpha) + c1 * alpha;
                return result;
}

// Spherical interpolation of two quaternions
// TODO: compute the spherical interpolation of two Quat inputs and return the
// result
inline Quat slerp(const Quat& q0, const Quat& q1, const double& alpha) {
                Quat result = Quat();
                Quat temp = q1 * inv(q0);
                if (temp[0] < 0) {
                  temp *= -1;
                }
                double p = atan2(
                    std::sqrt(std::pow(temp[1], 2) + std::pow(temp[2], 2) +
                              std::pow(temp[3], 2)),
                    temp[0]);
                double ap = alpha * p;

                const double threshold = 0.001;
                if (p > CS380_PI - threshold) {  // Large angle
                  // Use lerp instead
                  result = normalize(q0 * (1.0 - alpha) + q1 * alpha);
                } else if (p < threshold) {  // Small angle
                  // Use lerp instead
                  result = normalize(q0 * (1.0 - alpha) + q1 * alpha);
                } else {
                  result[0] = cos(ap);
                  result[1] = temp[1] / sin(p) * sin(ap);
                  result[2] = temp[2] / sin(p) * sin(ap);
                  result[3] = temp[3] / sin(p) * sin(ap);
                  result = result * q0;
                }
                return result;
}

// Linear interpolation of two RigTForms
// TODO: compute the linear interpolation of two RigTForm inputs and return the
// result Note: you should use the lerp and slerp functions you implemented
// above
inline RigTForm Linear(const RigTForm& rbt0,
                       const RigTForm& rbt1,
                       const double& alpha) {
                RigTForm interpolate_rbt = RigTForm();
                interpolate_rbt.setTranslation(
                    lerp(rbt0.getTranslation(), rbt1.getTranslation(), alpha));
                interpolate_rbt.setRotation(
                    slerp(rbt0.getRotation(), rbt1.getRotation(), alpha));
                return interpolate_rbt;
}

// Catmull-Rom interpolation of two RigTForms
// TODO: compute the Catmull-Rom interpolation of four RigTForm inputs and
// return the result Note: To Catmull-Rom interpolate two RigTFrom rbt0 and
// rbt1, we need 4 Keyframes keyframe rbt_1, rbt0, rbt1, rbt2. 		 You can use the
// lerp and slerp functions you implemented above.
inline RigTForm CatmullRom(const RigTForm& rbt_1,
                           const RigTForm& rbt0,
                           const RigTForm& rbt1,
                           const RigTForm& rbt2,
                           const double& alpha) {
 // For translation
  Cvec3 c_i = rbt0.getTranslation();
  Cvec3 c_i1 = rbt1.getTranslation();
  Cvec3 c_i_1 = rbt_1.getTranslation();
  Cvec3 c_i2 = rbt2.getTranslation();

  Cvec3 d_i = (c_i1 - c_i_1) * (1/6) + c_i;
  Cvec3 e_i = (c_i2 - c_i) * (-1/6) + c_i1;

  // For rotation
  Quat q_i = rbt0.getRotation();
  Quat q_i1 = rbt1.getRotation();
  Quat q_i_1 = rbt_1.getRotation();
  Quat q_i2 = rbt2.getRotation();

  Quat q_i_prime = slerp(Quat(), q_i1 * inv(q_i_1), 1.0 / 6.0) * q_i;
  Quat q_i1_prime = slerp(Quat(), q_i2 * inv(q_i), 1.0 / 6.0) * q_i1;

  // Interpolate the translation and rotation separately
  Cvec3 interpolated_translation = lerp(lerp(c_i, d_i, alpha), lerp(e_i, c_i1, alpha), alpha);
  Quat interpolated_rotation = slerp(slerp(q_i, q_i_prime, alpha), slerp(q_i1_prime, q_i1, alpha), alpha);

  // Return the interpolated RigTForm
  return RigTForm(interpolated_translation, interpolated_rotation);
}

}  // namespace Interpolation
#endif