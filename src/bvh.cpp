#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {

	size_t qsplit(std::vector<Primitive *> &_primitives, size_t start, size_t range, double pivot_val, size_t axis)
	{
		BBox bbox;
		double centroid;
		size_t ret_val = start;
		size_t size = start + range;

		for (int i = start; i < size; i++)
		{
			bbox = _primitives[i]->get_bbox();
			if (axis == 0)
				centroid = bbox.centroid().x;
			else if (axis == 1)
				centroid = bbox.centroid().y;
			else
				centroid = bbox.centroid().z;

			if (centroid < pivot_val)
			{
				std::swap(_primitives[i], _primitives[ret_val]);
				ret_val++;
			}
		}
		if (ret_val == start || ret_val == size) ret_val = (start + size) / 2;

		return ret_val;
	}

  BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
      size_t max_leaf_size) {

    this->primitives = _primitives;

    // TODO:
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration. The starter code build a BVH aggregate with a
    // single leaf node (which is also the root) that encloses all the
    // primitives.

    BBox bb;
    for (size_t i = 0; i < primitives.size(); ++i) {
      bb.expand(primitives[i]->get_bbox());
    }

    root = new BVHNode(bb, 0, primitives.size());

	if (primitives.size() > max_leaf_size)
	{
		size_t mid_point = qsplit(primitives, 0, primitives.size(), bb.centroid().x, 0);
		root->l = buildBranch(primitives, 0, mid_point, 1, max_leaf_size);
		root->r = buildBranch(primitives, mid_point, primitives.size() - mid_point, 1, max_leaf_size);
	}
	else
	{
		root->l = new BVHNode(bb, 0, primitives.size());
		root->r = new BVHNode(bb, 0, primitives.size());
	}
  }

  BVHAccel::~BVHAccel() {

    // TODO:
    // Implement a proper destructor for your BVH accelerator aggregate
	  if (root)
	  {
		  if (root->l)
			  deleteBranch(root->l);
		  if (root->r)
			  deleteBranch(root->r);

		  delete root;
	  }
  }

  BBox BVHAccel::get_bbox() const {
    return root->bb;
  }

  bool BVHAccel::intersect(const Ray &ray) const {

    // TODO:
    // Implement ray - bvh aggregate intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.
	  if (NULL == root)
		  return false;
	  return intersect(root, ray);
  }

  bool BVHAccel::intersect(const BVHNode* pNode, const Ray& ray) const
  {
	  double t0, t1;
	  if (!pNode->bb.intersect(ray, t0, t1))
		  return false;
	  if (pNode->l)
	  {
		  if (intersect(pNode->l, ray))
			  return true;
	  }

	  if (pNode->r)
	  {
		  if (intersect(pNode->r, ray))
			  return true;
	  }

	  if (pNode->isLeaf())
	  {
		  size_t size = pNode->start + pNode->range;
		  for (size_t p = pNode->start; p < size; ++p) 
		  {
			  if (primitives[p]->intersect(ray))
				  return true;
		  }
	  }

	  return false;
  }

  bool BVHAccel::intersect(const Ray &ray, Intersection *i) const {

    // TODO:
    // Implement ray - bvh aggregate intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate. When an intersection does happen.
    // You should store the non-aggregate primitive in the intersection data
    // and not the BVH aggregate itself.

	  if (NULL == root)
		  return false;
	  return intersect(root, ray, i);
  }

  bool BVHAccel::intersect(BVHNode* pNode, const Ray& ray, Intersection* i) const
  {
	  double t0, t1;
	  if (!pNode->bb.intersect(ray, t0, t1))
		  return false;

	  if (pNode->isLeaf())
	  {
		  bool hit = false;
		  size_t size = pNode->start + pNode->range;
		  for (size_t p = pNode->start; p < size; ++p)
		  {
			  if (primitives[p]->intersect(ray, i))
				  hit = true;
		  }

		  return hit;
	  }

	  bool hitL = false, hitR = false;
	  if (pNode->l)
		  hitL = intersect(pNode->l, ray, i);

	  if (pNode->r)
		  hitR = intersect(pNode->r, ray, i);

	  return (hitL || hitR);
  }

  BVHNode* BVHAccel::buildBranch(std::vector<Primitive *> &_primitives, size_t start, size_t range, size_t axis, size_t max_leaf_size)
  {
	  BBox bb;
	  size_t size = start + range;
	  for (size_t i = start; i < size; ++i)
		  bb.expand(primitives[i]->get_bbox());
		
	  BVHNode* pNode = new BVHNode(bb, start, range);
	  if (range <= max_leaf_size)
		  return pNode;

	  size_t mid_point = qsplit(primitives, start, range, bb.centroid().x, axis);
	  pNode->l = buildBranch(primitives, start, mid_point - start, (axis + 1) % 3, max_leaf_size);
	  pNode->r = buildBranch(primitives, mid_point, size - mid_point, (axis + 1) % 3, max_leaf_size);

	  return pNode;
  }

  void BVHAccel::deleteBranch(BVHNode* pNode)
  {
	  if (pNode->l)
		  deleteBranch(pNode->l);
	  if (pNode->r)
		  deleteBranch(pNode->r);
	  delete pNode;
  }

}  // namespace StaticScene
}  // namespace CMU462
