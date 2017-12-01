#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, vector<size_t>& v) :
    mesh(mesh), v(v) { }
Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {
  
  // TODO: 
  // compute the bounding box of the triangle
  
  return BBox();
}

bool Triangle::intersect(const Ray& r) const {
  
  // TODO: implement ray-triangle intersection
	Vector3D p0 = mesh->positions[v1];
	Vector3D p1 = mesh->positions[v2];
	Vector3D p2 = mesh->positions[v3];
	double t;
	double A = p0.x - p1.x;
	double B = p0.y - p1.y;
	double C = p0.z - p1.z;

	double D = p0.x - p2.x;
	double E = p0.y - p2.y;
	double F = p0.z - p2.z;

	double G = r.d.x;
	double H = r.d.y;
	double I = r.d.z;

	double J = p0.x - r.o.x;
	double K = p0.y - r.o.y;
	double L = p0.z - r.o.z;

	double EIHF = E*I - H*F;
	double GFDI = G*F - D*I;
	double DHEG = D*H - E*G;

	double denom = (A*EIHF + B*GFDI + C*DHEG);

	double beta = (J*EIHF + K*GFDI + L*DHEG) / denom;

	if (beta <= 0.0 || beta >= 1)
	{
		return false;
	}

	double AKJB = A*K - J*B;
	double JCAL = J*C - A*L;
	double BLKC = B*L - K*C;

	double gamma = (I*AKJB + H*JCAL + G*BLKC) / denom;
	if (gamma <= 0.0 || beta + gamma >= 1.0) return false;

	t = -(F*AKJB + E*JCAL + D*BLKC) / denom;
	if (t >= r.min_t && t <= r.max_t)
	{
		r.max_t = t;
		return true;
	}
	return false;
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // TODO: 
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
	Vector3D p0 = mesh->positions[v1];
	Vector3D p1 = mesh->positions[v2];
	Vector3D p2 = mesh->positions[v3];
	double t;
	double A = p0.x - p1.x;
	double B = p0.y - p1.y;
	double C = p0.z - p1.z;

	double D = p0.x - p2.x;
	double E = p0.y - p2.y;
	double F = p0.z - p2.z;

	double G = r.d.x;
	double H = r.d.y;
	double I = r.d.z;

	double J = p0.x - r.o.x;
	double K = p0.y - r.o.y;
	double L = p0.z - r.o.z;

	double EIHF = E*I - H*F;
	double GFDI = G*F - D*I;
	double DHEG = D*H - E*G;

	double denom = (A*EIHF + B*GFDI + C*DHEG);

	double beta = (J*EIHF + K*GFDI + L*DHEG) / denom;

	if (beta <= 0.0 || beta >= 1)
	{
		return false;
	}

	double AKJB = A*K - J*B;
	double JCAL = J*C - A*L;
	double BLKC = B*L - K*C;

	double gamma = (I*AKJB + H*JCAL + G*BLKC) / denom;
	if (gamma <= 0.0 || beta + gamma >= 1.0) return false;

	t = -(F*AKJB + E*JCAL + D*BLKC) / denom;
	if (t >= r.min_t && t <= r.max_t)
	{
		r.max_t = t;
		isect->t = t;
		isect->primitive = this;
		isect->bsdf = get_bsdf();
		Vector3D vNormal = cross((p1 - p0), (p2 - p0));
		vNormal.normalize();
		if (dot(vNormal, r.d) > 0.0) 
			vNormal = -vNormal;
		isect->n = vNormal;
		return true;
	}
	return false;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CMU462
