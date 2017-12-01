#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
	Vector3D temp = r.o - o;
	double twoa = 2.0f*dot(r.d, r.d);
	double b = 2.0f*dot(r.d, temp);
	double c = dot(temp, temp) - this->r*this->r;
	double discriminant = b*b - 2.0f*twoa*c;
	if (discriminant > 0.0f)
	{
		discriminant = sqrt(discriminant);
		double t = (-b - discriminant) / twoa;
		if (t < r.min_t) t = (-b + discriminant) / twoa;
		if (t < r.min_t || t > r.max_t) return false;

		t1 = (-b - discriminant) / twoa;
		t2 = (-b + discriminant) / twoa;
		return true;
	}

	return false;
}

bool Sphere::intersect(const Ray& r) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.

	double t1, t2;
	if (test(r,t1,t2))
	{
		if (t1 < r.min_t)
			r.max_t = t2;
		else
			r.max_t = t1;

		return true;
	}

    return false;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.

	double t1, t2;
	if (test(r, t1, t2))
	{
		if (t1 < r.min_t)
			r.max_t = t2;
		else
			r.max_t = t1;

		i->t = r.max_t;
		i->primitive = this;
		i->bsdf = get_bsdf();
		Vector3D hit = r.o + i->t * r.d;
		Vector3D vNormal = hit - o;
		vNormal.normalize();
		i->n = vNormal;
		return true;
	}

	return false;
}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CMU462
