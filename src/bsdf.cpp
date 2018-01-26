#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CMU462 {

	void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

		Vector3D z = Vector3D(n.x, n.y, n.z);
		Vector3D h = z;
		if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
		else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
		else h.z = 1.0;

		z.normalize();
		Vector3D y = cross(h, z);
		y.normalize();
		Vector3D x = cross(z, y);
		x.normalize();

		o2w[0] = x;
		o2w[1] = y;
		o2w[2] = z;
	}

	// Diffuse BSDF //

	Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
		return albedo * (1.0 / PI);
	}

	Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
		*wi = sampler.get_sample(pdf);
		return f(wo, *wi);
	}

	// Mirror BSDF //

	Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
		Vector3D N(0, 0, 1);
		Vector3D H = (wo + wi).unit();
		return reflectance * pow(std::max(dot(N, H), 0.0), roughness);
	}

	Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

		// TODO:
		// Implement MirrorBSDF
		reflect(wo, wi);
		*pdf = 1.0f;
		return reflectance * (1.0 / abs_cos_theta(wo));
	}

	// Glossy BSDF //

	/*
	   Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
	   return Spectrum();
	   }

	   Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
	 *pdf = 1.0f;
	 return reflect(wo, wi, reflectance);
	 }
	 */

	 // Refraction BSDF //

	Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
		return Spectrum();
	}

	Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

		// TODO:
		// Implement RefractionBSDF

		return Spectrum();
	}

	// Glass BSDF //

	Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
		return Spectrum();
		double no = 1.0;
		double ni = ior;
		double cos_theta_i = wi.z;
		if (cos_theta_i > 0.0)
		{
			cos_theta_i = -cos_theta_i;
			no = ior;
			ni = 1.0;
		}

		double Ro = ((ni - no) / (ni + no)) * ((ni - no) / (ni + no));
		double temp = 1 - cos_theta_i;
		double R = Ro + (1 - Ro) * temp * temp * temp * temp * temp;
		if (wo.z * wi.z >= 0.0)
		{
			// reflect
			Vector3D N(0, 0, 1);
			Vector3D H = (wo + wi).unit();
			return R * reflectance * pow(std::max(dot(N, H), 0.0), roughness);
		}
		else
		{
			// refract
			return (1 - R) * transmittance * abs_cos_theta(wo);
		}
	}

	Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
		// TODO:
		// Compute Fresnel coefficient and either reflect or refract based on it.
		if (refract(wo, wi, ior))
		{
			double no = 1.0;
			double ni = ior;
			double cos_theta_i = wi->z;
			if (cos_theta_i > 0.0)
			{
				no = ior;
				ni = 1.0;
			}
			cos_theta_i = fabs(cos_theta_i);
			double Ro = ((ni - no) / (ni + no)) * ((ni - no) / (ni + no));
			double temp = 1 - cos_theta_i;
			double R = Ro + (1 - Ro) * temp * temp * temp * temp * temp;
			double random = (double)(std::rand()) / RAND_MAX;
			double P = (R + 0.5) / 2.0;
			if (random <= P)
			{
				// reflect
				reflect(wo, wi);
				*pdf = R / P;
				return R*reflectance * (1.0 / abs_cos_theta(wo));
 			}
 			else
 			{
 				// refract
 				*pdf = 1.0 - P;
 				return (1.0 - R) * transmittance * ((ni / no) * (ni / no) / abs_cos_theta(wo));
 			}
		}
		else
		{
			// reflect
			reflect(wo, wi);
			*pdf = 1.0f;
			return reflectance * (1.0 / abs_cos_theta(wo));
		}
	}

	void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

		// TODO:
		// Implement reflection of wo about normal (0,0,1) and store result in wi.
		*wi = Vector3D(0, 0, 2 * wo[2]) - wo;
	}

	bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

		// TODO:
		// Use Snell's Law to refract wo surface and store result ray in wi.
		// Return false if refraction does not occur due to total internal reflection
		// and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
		// ray entering the surface through vacuum.

		double cos_theta_o = wo[2];
		double sin_theta_o = sin_theta(wo);
		double sin_theta_i = sin_theta_o / ior;
		if (cos_theta_o < 0.0)
		{
			double sin_theta_o = sin_theta(wo);
			double sin_theta_i = sin_theta_o * ior;
			if (sin_theta_i > 1.0)
				return false;
		}

		double _cos_phi = cos_phi(wo);
		double _sin_phi = sin_phi(wo);
		double xi = sin_theta_i * _cos_phi;
		double yi = sin_theta_i * _sin_phi;
		double zi = sqrtf(1.0 - sin_theta_i*sin_theta_i);
		*wi = Vector3D(-xi, -yi, -zi);

		return true;
	}

	// Emission BSDF //

	Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
		return Spectrum();
	}

	Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
		*wi = sampler.get_sample(pdf);
		return Spectrum();
	}

} // namespace CMU462
