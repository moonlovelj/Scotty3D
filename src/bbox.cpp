#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CMU462 {

  bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

    // TODO:
    // Implement ray - bounding box intersection test
    // If the ray intersected the bouding box within the range given by
    // t0, t1, update t0 and t1 with the new intersection times.

	  double interval_min = r.min_t;
	  double interval_max = r.max_t;
	  Vector3D pp[2];
	  pp[0] = min;
	  pp[1] = max;

	  double tt0 = (pp[r.sign[0]].x - r.o.x) * r.inv_d.x;
	  double tt1 = (pp[r.sign[0]^1].x - r.o.x) * r.inv_d.x;
	  if (tt0 > interval_min) interval_min = tt0;
	  if (tt1 < interval_max) interval_max = tt1;
	  if (interval_min > interval_max) return false;

	  tt0 = (pp[r.sign[1]].y - r.o.y) * r.inv_d.y;
	  tt1 = (pp[r.sign[1]^1].y - r.o.y) * r.inv_d.y;
	  if (tt0 > interval_min) interval_min = tt0;
	  if (tt1 < interval_max) interval_max = tt1;
	  if (interval_min > interval_max) return false;

	  tt0 = (pp[r.sign[2]].z - r.o.z) * r.inv_d.z;
	  tt1 = (pp[r.sign[2]^1].z - r.o.z) * r.inv_d.z;
	  if (tt0 > interval_min) interval_min = tt0;
	  if (tt1 < interval_max) interval_max = tt1;
	  if (interval_min > interval_max) return false;
	  
	  t0 = interval_min;
	  t1 = interval_max;
	  return true;
  }

  void BBox::draw(Color c) const {

    glColor4f(c.r, c.g, c.b, c.a);

    // top
    glBegin(GL_LINE_STRIP);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(max.x, max.y, max.z);
    glEnd();

    // bottom
    glBegin(GL_LINE_STRIP);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, min.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glEnd();

    // side
    glBegin(GL_LINES);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(min.x, min.y, max.z);
    glEnd();

  }

  std::ostream& operator<<(std::ostream& os, const BBox& b) {
    return os << "BBOX(" << b.min << ", " << b.max << ")";
  }

} // namespace CMU462
