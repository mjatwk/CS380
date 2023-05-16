#ifndef INTERPOLATION_H
#define INTERPOLATION_H


#include <iostream>
#include <cmath>
#include "rigtform.h"

namespace Interpolation {

	// Linear interpolation of two coordinate vectors
	// TODO: compute the linear interpolation of two Cvec3 inputs and return the result
	inline Cvec3 lerp(const Cvec3& c0, const Cvec3& c1, const double& alpha) {
		Cvec3 result = Cvec3();
		result = c0 * (1-alpha) + c1 * alpha;
		return result;
	}

	// Spherical interpolation of two quaternions
	// TODO: compute the spherical interpolation of two Quat inputs and return the result
	inline Quat slerp(const Quat& q0, const Quat& q1, const double& alpha) {
		Quat result = Quat();
		Quat temp = q1 * inv(q0);
		if (temp[0] < 0) {
			temp *= -1;
		} 
		double p = atan2(std::sqrt(std::pow(temp[1], 2) + std::pow(temp[2], 2) + std::pow(temp[3], 2)), temp[0]);
		double ap = alpha * p;
		
		const double threshold = 0.001;
		if (p > CS380_PI - threshold) { // Large angle
			// Use lerp instead
			result = normalize(q0 * (1.0 - alpha) + q1 * alpha );
		} else if (p < threshold) { // Small angle
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
	// TODO: compute the linear interpolation of two RigTForm inputs and return the result
	// Note: you should use the lerp and slerp functions you implemented above
	inline RigTForm Linear(const RigTForm& rbt0, const RigTForm& rbt1, const double& alpha) {
		RigTForm interpolate_rbt = RigTForm();
		interpolate_rbt.setTranslation(lerp(rbt0.getTranslation(), rbt1.getTranslation(), alpha));
		interpolate_rbt.setRotation(slerp(rbt0.getRotation(), rbt1.getRotation(), alpha));
		return interpolate_rbt;	
	}

}
#endif