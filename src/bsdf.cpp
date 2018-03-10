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
		return reflectance * (1.0 / abs_cos_theta(*wi));
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

	double GlassBSDF::smithG1(const Vector3D &v, const Vector3D &m, double roughness) const {
		const double tanTheta = std::abs(tan_theta(v));
		/* perpendicular incidence -- no shadowing/masking */
		if (tanTheta == 0.0f)
			return 1.0f;

		/* Can't see the back side from the front and vice versa */
		if (dot(v, m) * cos_theta(v) <= 0)
			return 0.0f;

		roughness = std::sqrt(0.5f * roughness + 1) / tanTheta;
			const double a = 1.0f / (roughness * tanTheta);
			const double aSqr = a * a;
			if (a >= 1.6f)
				return 1.0f;

			return (3.535f * a + 2.181f * aSqr)
				/ (1.0f + 2.276f * a + 2.577f * aSqr);
	}

	double GlassBSDF::fresnel(double cosThetaI) const {
		double etaI = ior, etaT = 1.0;

		/* Swap the indices of refraction if the interaction starts
		at the inside of the object */
		if (cosThetaI > 0.0)
			std::swap(etaI, etaT);

		/* Using Snell's law, calculate the sine of the angle
		between the transmitted ray and the surface normal */
		double sinThetaT = etaI / etaT *
			std::sqrt(std::max(0.0, 1.0 - cosThetaI*cosThetaI));

		if (sinThetaT > 1.0)
			return 1.0;  /* Total internal reflection! */

		double cosThetaT = std::sqrt(1.0 - sinThetaT*sinThetaT);

		/* Finally compute the reflection coefficient */
		cosThetaI = fabs(cosThetaI);
		double Rs = (etaI * cosThetaI - etaT * cosThetaT) / (etaI * cosThetaI + etaT * cosThetaT);
		double Rp = (etaT * cosThetaI - etaI * cosThetaT) / (etaT * cosThetaI + etaI * cosThetaT);
		return (Rs * Rs + Rp * Rp) / 2.0;
	}

	// Glass BSDF //

	Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
		bool reflect = (wi.z * wo.z) > 0.0;
		if (reflect)
		{
			// reflect
			Vector3D N(0, 0, 1);
			Vector3D H = ((wo + wi)* abs_cos_theta(wo)).unit();
			if (cos_theta(H) <= 0.0)
				return Spectrum();
			double D = (roughness + 2) / (2 * PI)* std::pow(cos_theta(H), roughness);
			double F = fresnel(dot(wi, H));
			double G = smithG1(wi, H, roughness) * smithG1(wo, H, roughness);
			return D * F * G / (4 * wo.z * wi.z) * reflectance;
		}
		else
		{
			// refract
			double etaI = ior, etaT = 1.0;
			if (wi.z > 0.0)
				std::swap(etaI, etaT);
			Vector3D H = (ior > 1.0 ? 1.0 : -1.0)* (wi*etaI + wo*etaT).unit();
			if (cos_theta(H) <= 0.0)
				return Spectrum();
			double D = (roughness + 2) / (2 * PI)* std::pow(std::max(0.0,cos_theta(H)), roughness);
			double F = fresnel(dot(wi, H));
			double G = smithG1(wi, H, roughness) * smithG1(wo, H, roughness);

			double sqrtDenom = etaI * dot(wi, H) + etaT * dot(wo, H);
			double value = ((1 - F) * D * G * etaT * etaT * dot(wi, H)*dot(wo, H)) /
				(cos_theta(wi) * cos_theta(wo) * sqrtDenom * sqrtDenom);


			return value * transmittance ;
		}
	}

	Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
		// TODO:
		// Compute Fresnel coefficient and either reflect or refract based on it.
		if (refract(wo, wi, ior))
		{
			
			double etaI = ior, etaT = 1.0;
			double cosThetaI = cos_theta(*wi);
			if (cosThetaI > 0.0)
				std::swap(etaI, etaT);
			double R = fresnel(cosThetaI);
			double random = (double)(std::rand()) / RAND_MAX;
			if (random <= R)
			{
				// reflect
				reflect(wo, wi);
				*pdf = R ;
				return R * reflectance * (1.0 / std::max(std::fabs((*wi)[2]), 1e-8));
 			}
 			else
 			{
 				// refract
 				*pdf = (1.0 - R);
 				return (1.0 - R) * transmittance * ((etaI / etaT) * (etaI / etaT) / std::max(std::fabs((*wi)[2]), 1e-8));
 			}
		}
		else
		{
			// reflect
			reflect(wo, wi);
			*pdf = 1.0;
			return reflectance * (1.0 / std::max(std::fabs((*wi)[2]), 1e-8));
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
		double etaI = ior, etaT = 1;
		if (cos_theta_o < 0.0)
			std::swap(etaI, etaT);

		double sin_theta_i = sin_theta_o * etaT / etaI;
		if (sin_theta_i > 1.0)
			return false;

		Vector3D N(0, 0, 1);
		Vector3D TX = (wo - N * cos_theta_o).unit();
		*wi = N * (cos_theta_o < 0.0 ? 1 : -1) * sqrt(1.0 - sin_theta_i* sin_theta_i) - TX * sin_theta_i;
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
