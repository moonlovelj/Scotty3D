#include <float.h>
#include <assert.h>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"

namespace CMU462 {

  VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {

    // TODO: (meshEdit)
    // This method should split the given edge and return an iterator to the
    // newly inserted vertex. The halfedge of this vertex should point along
    // the edge that was split, rather than the new edges.
	  if (e0->isBoundary())
		  return e0->halfedge()->vertex();

	  // collect elements
	  HalfedgeIter h0 = e0->halfedge();
	  HalfedgeIter h1 = h0->next();
	  HalfedgeIter h2 = h1->next();
	  HalfedgeIter h3 = h0->twin();
	  HalfedgeIter h4 = h3->next();
	  HalfedgeIter h5 = h4->next();
	  VertexIter v0 = h1->vertex();
	  VertexIter v1 = h2->vertex();
	  VertexIter v2 = h4->vertex();
	  VertexIter v3 = h5->vertex();
	  FaceIter f0 = h0->face();
	  FaceIter f1 = h3->face();

	  //allocate new elements
	  VertexIter v4 = newVertex();
	  EdgeIter e1 = newEdge();
	  EdgeIter e2 = newEdge();
	  EdgeIter e3 = newEdge();
	  FaceIter f2 = newFace();
	  FaceIter f3 = newFace();
	  HalfedgeIter h6 = newHalfedge();
	  HalfedgeIter h7 = newHalfedge();
	  HalfedgeIter h8 = newHalfedge();
	  HalfedgeIter h9 = newHalfedge();
	  HalfedgeIter h10 = newHalfedge();
	  HalfedgeIter h11 = newHalfedge();

	  // reassign elements
	  v0->halfedge() = h1;
	  v1->halfedge() = h2;
	  v2->halfedge() = h4;
	  v3->halfedge() = h5;
	  v4->halfedge() = h0;
	  e0->halfedge() = h0;
	  e1->halfedge() = h7;
	  e2->halfedge() = h9;
	  e3->halfedge() = h11;
	  f0->halfedge() = h0;
	  f1->halfedge() = h7;
	  f2->halfedge() = h9;
	  f3->halfedge() = h11;
	  h0->setNeighbors(h1, h3, v4, e0, f0);
	  h1->setNeighbors(h6, h1->twin(), v0, h1->edge(), f0);
	  h2->setNeighbors(h8, h2->twin(), v1, h2->edge(), f1);
	  h3->setNeighbors(h11, h0, v0, e0, f3);
	  h4->setNeighbors(h10, h4->twin(), v2, h4->edge(), f2);
	  h5->setNeighbors(h3, h5->twin(), v3, h5->edge(), f3);
	  h6->setNeighbors(h0, h7, v1, e1, f0);
	  h7->setNeighbors(h2, h6, v4, e1, f1);
	  h8->setNeighbors(h7, h9, v2, e2, f1);
	  h9->setNeighbors(h4, h8, v4, e2, f2);
	  h10->setNeighbors(h9, h11, v3, e3, f2);
	  h11->setNeighbors(h5, h10, v4, e3, f3);
	  v4->position = 0.5 * (v0->position + v2->position);
	  e0->isNew = e2->isNew = false;
	  e1->isNew = e3->isNew = true;

    return v4;
  }

  VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {

    // TODO: (meshEdit)
    // This method should collapse the given edge and return an iterator to
    // the new vertex created by the collapse.
	  if (e->isBoundary())
		  return e->halfedge()->vertex();

	  if (e->halfedge()->face()->isBoundary() || e->halfedge()->twin()->face()->isBoundary())
		  return e->halfedge()->vertex();

	  // collect elements
	  HalfedgeIter h0 = e->halfedge();
	  HalfedgeIter h1 = h0->next();
	  HalfedgeIter h2 = h1->next();
	  HalfedgeIter h3 = h0->twin();
	  HalfedgeIter h4 = h3->next();
	  HalfedgeIter h5 = h4->next();
	  HalfedgeIter h6 = h5->twin();
	  HalfedgeIter h7 = h4->twin();
	  HalfedgeIter h8 = h2->twin();
	  HalfedgeIter h9 = h1->twin();
	  HalfedgeIter h10 = h9->next();
	  HalfedgeIter h11 = h10->twin();
	  HalfedgeIter h12 = h11->next();
	  HalfedgeIter h13 = h12->twin();
	  HalfedgeIter h14 = h13->next();
	  HalfedgeIter h15 = h14->twin();
	  HalfedgeIter h16 = h7->next();
	  HalfedgeIter h17 = h16->twin();
	  HalfedgeIter h18 = h17->next();
	  HalfedgeIter h19 = h18->twin();
	  HalfedgeIter h20 = h19->next();
	  HalfedgeIter h21 = h20->twin();

	  VertexIter v0 = h5->vertex();
	  VertexIter v1 = h2->vertex();
	  VertexIter v2 = h3->vertex();
	  VertexIter v3 = h0->vertex();

	  EdgeIter e1 = h5->edge();
	  EdgeIter e2 = h1->edge();
	  EdgeIter e3 = h2->edge();
	  EdgeIter e4 = h4->edge();

	  FaceIter f0 = h0->face();
	  FaceIter f1 = h3->face();

	  // reassign elements
	  h6->setNeighbors(h6->next(), h7, v2, e1, h6->face());
	  h7->setNeighbors(h7->next(), h6, v0, e1, h7->face());
	  h8->setNeighbors(h8->next(), h9, v2, e2, h8->face());
	  h9->setNeighbors(h9->next(), h8, v1, e2, h9->face());
	  h10->setNeighbors(h10->next(), h11, v2, h10->edge(), h10->face());
	  h11->setNeighbors(h11->next(), h10, h11->vertex(), h11->edge(), h11->face());
	  h12->setNeighbors(h12->next(), h13, v2, h12->edge(), h12->face());
	  h13->setNeighbors(h13->next(), h12, h13->vertex(), h13->edge(), h13->face());
	  h14->setNeighbors(h14->next(), h15, v2, h14->edge(), h14->face());
	  h15->setNeighbors(h6, h14, h15->vertex(), h15->edge(), h15->face());
	  h16->setNeighbors(h16->next(), h17, v2, h16->edge(), h16->face());
	  h17->setNeighbors(h17->next(), h16, h17->vertex(), h17->edge(), h17->face());
	  h18->setNeighbors(h18->next(), h19, v2, h18->edge(), h18->face());
	  h19->setNeighbors(h19->next(), h18, h19->vertex(), h19->edge(), h19->face());
	  h20->setNeighbors(h20->next(), h21, v2, h20->edge(), h20->face());
	  h21->setNeighbors(h8, h20, h21->vertex(), h21->edge(), h21->face());

	  v0->halfedge() = h7;
	  v1->halfedge() = h9;
	  v2->halfedge() = h6;
	  v2->position = 0.5f * (v2->position + v3->position);

	  e1->halfedge() = h6;
	  e2->halfedge() = h8;


	  // delete unused elements
	  deleteHalfedge(h0);
	  deleteHalfedge(h1);
	  deleteHalfedge(h2);
	  deleteHalfedge(h3);
	  deleteHalfedge(h4);
	  deleteHalfedge(h5);

	  deleteVertex(v3);

	  deleteEdge(e);
	  deleteEdge(e3);
	  deleteEdge(e4);

	  deleteFace(f0);
	  deleteFace(f1);

    return v2;
  }

  VertexIter HalfedgeMesh::collapseFace(FaceIter f) {

    // TODO: (meshEdit)
    // This method should collapse the given face and return an iterator to
    // the new vertex created by the collapse.
	
	  return VertexIter();
  }

  FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {

    // TODO: (meshEdit)
    // This method should replace the given vertex and all its neighboring
    // edges and faces with a single face, returning the new face.
    return FaceIter();
  }

  FaceIter HalfedgeMesh::eraseEdge( EdgeIter e ) {
    // TODO: (meshEdit)
    // This method should erase the given edge and return an iterator to the
    // merged face.
    return FaceIter();
  }

  EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {

    // TODO: (meshEdit)
    // This method should flip the given edge and return an iterator to the
    // flipped edge.

	  if (e0->isBoundary())
		  return e0;

	  // collect elements
	  HalfedgeIter h1 = e0->halfedge();
	  HalfedgeIter h2 = h1->next();
	  HalfedgeIter h0 = h2->next();
	  HalfedgeIter h3 = h1->twin();
	  HalfedgeIter h4 = h3->next();
	  HalfedgeIter h5 = h4->next();
	  VertexIter v0 = h0->vertex();
	  VertexIter v1 = h4->vertex();
	  VertexIter v2 = h5->vertex();
	  VertexIter v3 = h2->vertex();
	  FaceIter f0 = h1->face();
	  FaceIter f1 = h3->face();

	  // reassign elements
	  v0->halfedge() = h0;
	  v1->halfedge() = h4;
	  v2->halfedge() = h5;
	  v3->halfedge() = h2;
	  e0->halfedge() = h1;	
	  f0->halfedge() = h1;
	  f1->halfedge() = h3;
	  h0->setNeighbors(h4, h0->twin(), v0, h0->edge(), f0);
	  h1->setNeighbors(h0, h3, v2, e0, f0);
	  h2->setNeighbors(h3, h2->twin(), v3, h2->edge(), f1);
	  h3->setNeighbors(h5, h1, v0, e0, f1);
	  h4->setNeighbors(h1, h4->twin(), v1, h4->edge(), f0);
	  h5->setNeighbors(h2, h5->twin(), v2, h5->edge(), f1);

	  return e0;
  }

  void HalfedgeMesh::subdivideQuad( bool useCatmullClark )
  {
    // Unlike the local mesh operations (like bevel or edge flip), we will perform
    // subdivision by splitting *all* faces into quads "simultaneously."  Rather
    // than operating directly on the halfedge data structure (which as you've seen
    // is quite difficult to maintain!) we are going to do something a bit nicer:
    //
    //    1. Create a raw list of vertex positions and faces (rather than a full-
    //       blown halfedge mesh).
    //
    //    2. Build a new halfedge mesh from these lists, replacing the old one.
    //
    // Sometimes rebuilding a data structure from scratch is simpler (and even more
    // efficient) than incrementally modifying the existing one.  These steps are
    // detailed below.

    // TODO Step I: Compute the vertex positions for the subdivided mesh.  Here we're
    // going to do something a little bit strange: since we will have one vertex in
    // the subdivided mesh for each vertex, edge, and face in the original mesh, we
    // can nicely store the new vertex *positions* as attributes on vertices, edges,
    // and faces of the original mesh.  These positions can then be conveniently
    // copied into the new, subdivided mesh.
    // [See subroutines for actual "TODO"s]
    if( useCatmullClark )
    {
      computeCatmullClarkPositions();
    }
    else
    {
      computeLinearSubdivisionPositions();
    }

    // TODO Step II: Assign a unique index (starting at 0) to each vertex, edge, and
    // face in the original mesh.  These indices will be the indices of the vertices
    // in the new (subdivided mesh).  They do not have to be assigned in any particular
    // order, so long as no index is shared by more than one mesh element, and the
    // total number of indices is equal to V+E+F, i.e., the total number of vertices
    // plus edges plus faces in the original mesh.  Basically we just need a one-to-one
    // mapping between original mesh elements and subdivided mesh vertices.
    // [See subroutine for actual "TODO"s]
    assignSubdivisionIndices();

    // TODO Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
    // the element indices defined above.  In other words, each new quad should be of
    // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
    // original mesh elements.  Note that it is essential to get the orientation right
    // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
    // circulate in the same direction as old faces (think about the right-hand rule).
    // [See subroutines for actual "TODO"s]
    vector< vector<Index> > subDFaces;
    vector< Vector3D > subDVertices;
    buildSubdivisionFaceList( subDFaces );
    buildSubdivisionVertexList( subDVertices );

    // TODO Step IV: Pass the list of vertices and quads to a routine that clears the
    // internal data for this halfedge mesh, and builds new halfedge data from scratch,
    // using the two lists.
    rebuild( subDFaces, subDVertices );
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * simple linear interpolation, e.g., the edge midpoints and face
   * centroids.
   */
  void HalfedgeMesh::computeLinearSubdivisionPositions()
  {
    // TODO For each vertex, assign Vertex::newPosition to
    // its original position, Vertex::position.
	  for (VertexIter it = vertices.begin(); it != vertices.end(); it++)
	  {
		  it->newPosition = it->position;
	  }

    // TODO For each edge, assign the midpoint of the two original
    // positions to Edge::newPosition.
	  for (EdgeIter it = edges.begin(); it != edges.end(); it++)
	  {
		  it->newPosition = 0.5f * (it->halfedge()->vertex()->position + it->halfedge()->twin()->vertex()->position);
	  }

    // TODO For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::newPosition.  Note
    // that in general, NOT all faces will be triangles!
	  for (FaceIter it = faces.begin(); it != faces.end(); it++)
	  {
		  int n = 0;
		  HalfedgeIter hIt = it->halfedge();
		  Vector3D allPosition;
		  do 
		  {
			  VertexIter vIt = hIt->vertex();
			  allPosition += vIt->position;
			  hIt = hIt->next();
			  ++n;
		  } while (hIt != it->halfedge());

		  it->newPosition = allPosition / n;
	  }
  }

  /**
   * Compute new vertex positions for a mesh that splits each polygon
   * into quads (by inserting a vertex at the face midpoint and each
   * of the edge midpoints).  The new vertex positions will be stored
   * in the members Vertex::newPosition, Edge::newPosition, and
   * Face::newPosition.  The values of the positions are based on
   * the Catmull-Clark rules for subdivision.
   */
  void HalfedgeMesh::computeCatmullClarkPositions()
  {
    // TODO The implementation for this routine should be
    // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules.  (These rules are outlined in the Developer Manual.)

    // TODO face
	  for (FaceIter it = faces.begin(); it != faces.end(); it++)
	  {
		  int n = 0;
		  HalfedgeIter hIt = it->halfedge();
		  Vector3D allPosition;
		  do
		  {
			  VertexIter vIt = hIt->vertex();
			  allPosition += vIt->position;
			  hIt = hIt->next();
			  ++n;
		  } while (hIt != it->halfedge());

		  it->newPosition = allPosition / n;
	  }
    // TODO edges
	  for (EdgeIter it = edges.begin(); it != edges.end(); it++)
	  {
		  it->newPosition = (it->halfedge()->vertex()->position + it->halfedge()->twin()->vertex()->position
								+ it->halfedge()->face()->newPosition + it->halfedge()->twin()->face()->newPosition) * 0.25;
	  }

    // TODO vertices
	  for (VertexIter it = vertices.begin(); it != vertices.end(); it++)
	  {
		  Size n = it->degree();
		  Vector3D Q, R, S = it->position;
		  HalfedgeIter h = it->halfedge();
		  do 
		  {
			  Q += h->face()->newPosition;
			  R += h->edge()->newPosition;
			  h = h->twin()->next();
		  } while (h != it->halfedge());
		  Q /= n;
		  R /= n;
		  it->newPosition = (Q + 2 * R + (n - 3)*S) / n;
	  }

  }

  /**
   * Assign a unique integer index to each vertex, edge, and face in
   * the mesh, starting at 0 and incrementing by 1 for each element.
   * These indices will be used as the vertex indices for a mesh
   * subdivided using Catmull-Clark (or linear) subdivision.
   */
  void HalfedgeMesh::assignSubdivisionIndices()
  {
    // TODO Start a counter at zero; if you like, you can use the
    // "Index" type (defined in halfedgeMesh.h)
	  Index index = 0;
    // TODO Iterate over vertices, assigning values to Vertex::index
	for (VertexIter it = vertices.begin(); it != vertices.end(); it++)
	{
		it->index = index;
		++index;
	}
    // TODO Iterate over edges, assigning values to Edge::index
	for (EdgeIter it = edges.begin(); it != edges.end(); it++)
	{
		it->index = index;
		++index;
	}
    // TODO Iterate over faces, assigning values to Face::index
	for (FaceIter it = faces.begin(); it != faces.end(); it++)
	{
		it->index = index;
		++index;
	}
  }

  /**
   * Build a flat list containing all the vertex positions for a
   * Catmull-Clark (or linear) subdivison of this mesh.  The order of
   * vertex positions in this list must be identical to the order
   * of indices assigned to Vertex::newPosition, Edge::newPosition,
   * and Face::newPosition.
   */
  void HalfedgeMesh::buildSubdivisionVertexList( vector<Vector3D>& subDVertices )
  {
    // TODO Resize the vertex list so that it can hold all the vertices.

    // TODO Iterate over vertices, assigning Vertex::newPosition to the appropriate
    // location in the new vertex list.
	  for (VertexIter it = vertices.begin(); it != vertices.end(); it++)
	  {
		  subDVertices.push_back(it->newPosition);
	  }
    // TODO Iterate over edges, assigning Edge::newPosition to the appropriate
    // location in the new vertex list.
	  for (EdgeIter it = edges.begin(); it != edges.end(); it++)
	  {
		  subDVertices.push_back(it->newPosition);
	  }
    // TODO Iterate over faces, assigning Face::newPosition to the appropriate
    // location in the new vertex list.
	  for (FaceIter it = faces.begin(); it != faces.end(); it++)
	  {
		  subDVertices.push_back(it->newPosition);
	  }
  }

  /**
   * Build a flat list containing all the quads in a Catmull-Clark
   * (or linear) subdivision of this mesh.  Each quad is specified
   * by a vector of four indices (i,j,k,l), which come from the
   * members Vertex::index, Edge::index, and Face::index.  Note that
   * the ordering of these indices is important because it determines
   * the orientation of the new quads; it is also important to avoid
   * "bowties."  For instance, (l,k,j,i) has the opposite orientation
   * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
   * will look like a bowtie.
   */
  void HalfedgeMesh::buildSubdivisionFaceList( vector< vector<Index> >& subDFaces )
  {
    // TODO This routine is perhaps the most tricky step in the construction of
    // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
    // vertex positions).  Basically what you want to do is iterate over faces,
    // then for each for each face, append N quads to the list (where N is the
    // degree of the face).  For this routine, it may be more convenient to simply
    // append quads to the end of the list (rather than allocating it ahead of
    // time), though YMMV.  You can of course iterate around a face by starting
    // with its first halfedge and following the "next" pointer until you get
    // back to the beginning.  The tricky part is making sure you grab the right
    // indices in the right order---remember that there are indices on vertices,
    // edges, AND faces of the original mesh.  All of these should get used.  Also
    // remember that you must have FOUR indices per face, since you are making a
    // QUAD mesh!

    // TODO iterate over faces
    // TODO loop around face
    // TODO build lists of four indices for each sub-quad
    // TODO append each list of four indices to face list
	  for (FaceIter fIt = faces.begin(); fIt != faces.end(); fIt++)
	  {
		  HalfedgeIter hIt = fIt->halfedge();
		  do 
		  {
			  vector<Index> temp;
			  temp.push_back(hIt->edge()->index);
			  temp.push_back(hIt->next()->vertex()->index);
			  temp.push_back(hIt->next()->edge()->index);
			  temp.push_back(fIt->index);
			  subDFaces.push_back(temp);
			  hIt = hIt->next();
		  } while (hIt != fIt->halfedge());
	  }
  }

  void HalfedgeMesh::_bevel_fc_reposition_with_dist( vector<Vector3D>& orig, // list of vertex positions of the original face (before bevel)
                                                     vector<HalfedgeIter>& hs, // list of halfedges pointing from the vertices of the new, beveled face to the old, original face
                                                     double shift, // user-requested amount to shift the face in the normal direction
                                                     double inset ) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled face.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //

	  int N = hs.size();
	  for (int i = 0; i < hs.size(); i++)
	  {
		  //Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
		  int a = (i + N - 1) % N;
		  int b = i;
		  int c = (i + 1) % N;
		  // Get the actual 3D vertex coordinates at these vertices
		  Vector3D pa = orig[a];
		  Vector3D pb = orig[b];
		  Vector3D pc = orig[c];
		  Vector3D v1 = pa - pb;
		  Vector3D v2 = pc - pb;
		  Vector3D v = v1 + v2;
		  v.normalize();
		  hs[i]->vertex()->position = pb + v * inset;
	  }
	  
	  if (N == 0)
		  return;

	  Vector3D fNormal = hs[0]->twin()->next()->twin()->face()->normal();
	  Vector3D fNormal1 = hs[1]->twin()->next()->twin()->face()->normal();
	  for (int i = 0; i < hs.size(); i++)
	  {
		  hs[i]->vertex()->position = hs[i]->vertex()->position +shift * fNormal;
	  }

	  for (int i = 0; i < orig.size(); i++)
	  {
		  orig[i] = hs[i]->vertex()->position;
	  }
	  
  }

  void HalfedgeMesh::_bevel_vtx_reposition_with_dist( Vector3D orig, // original vertex position, before the bevel
                                                      vector<HalfedgeIter>& hs, // list of halfedges pointing from the vertices of the new, beveled face to the neighbors of the original vertex
                                                      double inset ) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled vertex.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //
  }

  void HalfedgeMesh::_bevel_edge_reposition_with_dist( vector<Vector3D>& origs,  // list of vertex positions of the neighbors of the two endpoints of the edge, before the bevel
                                                       vector<HalfedgeIter>& hs,  // list of halfedges pointing from the vertices of the new, beveled face to the neighbors of the endpoints of the old, original edge
                                                       double inset) // user-requested amount by which to inset (i.e., grow/shrink) the beveled face
  {
    // TODO Compute new vertex positions for the vertices of the beveled edge.
    //
    // These vertices can be accessed via hs[i]->vertex()->position for i = 1, ..., hs.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the orig array) to compute an offset vertex position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in hs and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < hs.size(); hs++ )
    // {
    //    Vector3D pi = orig[i]; // get the original vertex position correponding to vertex i
    // }
    //
  }

  FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {

    // TODO This method should replace the vertex v with a face, corresponding to a bevel operation.
    // It should return the new face.  NOTE: This method is responsible for updating the *connectivity*
    // of the mesh only---it does not need to update the vertex positions.  These positions will be
    // updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to implement!)

    return facesBegin();
  }

  FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {

    // TODO This method should replace the edge e with a face, corresponding to a bevel operation.
    // It should return the new face.  NOTE: This method is responsible for updating the *connectivity*
    // of the mesh only---it does not need to update the vertex positions.  These positions will be
    // updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to implement!)

    return facesBegin();
  }

  FaceIter HalfedgeMesh::bevelFace(FaceIter f) {

    // TODO This method should replace the face f with an additional, inset face (and ring of faces around it),
    // corresponding to a bevel operation. It should return the new face.  NOTE: This method is responsible for
    // updating the *connectivity* of the mesh only---it does not need to update the vertex positions.  These
    // positions will be updated in HalfedgeMesh::_bevel_vtx_reposition_with_dist (which you also have to
    // implement!)
	  HalfedgeIter h0 = f->halfedge();
	  HalfedgeIter h1 = h0->next();
	  HalfedgeIter h2 = h1->next();
	  VertexIter v0 = h0->vertex();
	  VertexIter v1 = h1->vertex();
	  VertexIter v2 = h2->vertex();

	  VertexIter v3 = newVertex();
	  VertexIter v4 = newVertex();
	  VertexIter v5 = newVertex();
	  EdgeIter e0 = newEdge();
	  EdgeIter e1 = newEdge();
	  EdgeIter e2 = newEdge();
	  EdgeIter e3 = newEdge();
	  EdgeIter e4 = newEdge();
	  EdgeIter e5 = newEdge();
	  FaceIter f0 = newFace();
	  FaceIter f1 = newFace();
	  FaceIter f2 = newFace();
	  HalfedgeIter h3 = newHalfedge();
	  HalfedgeIter h4 = newHalfedge();
	  HalfedgeIter h5 = newHalfedge();
	  HalfedgeIter h6 = newHalfedge();
	  HalfedgeIter h7 = newHalfedge();
	  HalfedgeIter h8 = newHalfedge();
	  HalfedgeIter h9 = newHalfedge();
	  HalfedgeIter h10 = newHalfedge();
	  HalfedgeIter h11 = newHalfedge();
	  HalfedgeIter h12 = newHalfedge();
	  HalfedgeIter h13 = newHalfedge();
	  HalfedgeIter h14 = newHalfedge();
	  
	  v3->halfedge() = h5;
	  v4->halfedge() = h11;
	  v5->halfedge() = h8;
	  e0->halfedge() = h3;
	  e1->halfedge() = h9;
	  e2->halfedge() = h6;
	  e3->halfedge() = h14;
	  e4->halfedge() = h12;
	  e5->halfedge() = h13;
	  f->halfedge() = h12;
	  f0->halfedge() = h0;
	  f1->halfedge() = h1;
	  f2->halfedge() = h2;
	  h0->setNeighbors(h9, h0->twin(), v0, h0->edge(), f0);
	  h1->setNeighbors(h6, h1->twin(), v1, h1->edge(), f1);
	  h2->setNeighbors(h3, h2->twin(), v2, h2->edge(), f2);
	  h3->setNeighbors(h4, h11, v0, e0, f2);
	  h4->setNeighbors(h5, h14, v4, e3, f2);
	  h5->setNeighbors(h2, h6, v3, e2, f2);
	  h6->setNeighbors(h7, h5, v2, e2, f1);
	  h7->setNeighbors(h8, h13, v3, e5, f1);
	  h8->setNeighbors(h1, h9, v5, e1, f1);
	  h9->setNeighbors(h10, h8, v1, e1, f0);
	  h10->setNeighbors(h11, h12, v5, e4, f0);
	  h11->setNeighbors(h0, h3, v4, e0, f0);
	  h12->setNeighbors(h13, h10, v4, e4, f);
	  h13->setNeighbors(h14, h7, v5, e5, f);
	  h14->setNeighbors(h12, h4, v3, e3, f);

	  return f;
  }

  void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
    for (auto f : fcs) splitPolygon(f);
  }

  void HalfedgeMesh::splitPolygon(FaceIter f) {
    // TODO triangulation
	  if (f->degree() <= 3)
		  return;

	  HalfedgeIter hBegin = f->halfedge();
	  HalfedgeIter hLast = hBegin;
	  HalfedgeIter h = hBegin->next();
	  FaceIter curF = f;
	  FaceIter nextF = newFace();
	  bool bContinue = true;
	  do 
	  {
		  EdgeIter e = newEdge();
		  HalfedgeIter h1 = newHalfedge();
		  HalfedgeIter h2 = newHalfedge();
		  e->halfedge() = h1;
		  curF->halfedge() = hLast;
		  h1->setNeighbors(hLast, h2, h->next()->vertex(), e, curF);
		  h2->setNeighbors(h->next(), h1, hBegin->vertex(), e, nextF);
		  hLast->face() = curF;
		  h->face() = curF;
		  hLast = h2;
		  curF = nextF;
		  HalfedgeIter hTemp = h;
		  h = h->next();
		  hTemp->next() = h1;
		  if (h->next()->next() != hBegin)
			  nextF = newFace();
		  else
		  {
			  nextF->halfedge() = h2;
			  h->face() = nextF;
			  h->next()->face() = nextF;
			  h->next()->next() = h2;
			  bContinue = false;
		  }

	  } while (bContinue);

  }

  EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {

    // TODO: (meshEdit)
    // Compute the combined quadric from the edge endpoints.
    // -> Build the 3x3 linear system whose solution minimizes the quadric error
    //    associated with these two endpoints.
    // -> Use this system to solve for the optimal position, and store it in
    //    EdgeRecord::optimalPoint.
    // -> Also store the cost associated with collapsing this edg in
    //    EdgeRecord::Cost.

  }

  void MeshResampler::upsample(HalfedgeMesh& mesh)
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
  {

    // TODO: (meshEdit)
    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::newPosition.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::newPosition.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::isNew. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::position.

    // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order.  For future
    // reference, we're also going to store some information about which subdivided
    // edges come from splitting an edge in the original mesh, and which edges are new.
    // In this loop, we only want to iterate over edges of the original mesh---otherwise,
    // we'll end up splitting edges that we just split (and the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
	  for (VertexIter it = mesh.verticesBegin();it != mesh.verticesEnd(); ++it)
	  {
		  it->isNew = false;
		  HalfedgeIter h = it->halfedge();
		  Size n = it->degree();
		  double u = 0;
		  if (n == 3) 
			  u = 3.0 / 16;
		  else 
			  u = 3.0 / (8.0*n);

		  Vector3D newPosition;
		  do 
		  {
			  newPosition += h->next()->vertex()->position;
			  h = h->twin()->next();
		  } while (h != it->halfedge());

		  it->newPosition = (1.0 - n * u) * it->position + u * newPosition;
	  }

	  for (EdgeIter it = mesh.edgesBegin(); it != mesh.edgesEnd(); it++)
	  {
		  HalfedgeIter eH = it->halfedge();
		  it->newPosition= (eH->vertex()->position + eH->twin()->vertex()->position) * 3.0 / 8 +
			  (eH->next()->next()->vertex()->position + eH->twin()->next()->next()->vertex()->position) * 1.0 / 8;
	  }
	  // iterate over all edges in the mesh
	  int n = mesh.nEdges();
	  EdgeIter e = mesh.edgesBegin();
	  for (int i = 0; i < n; i++) {

		  // get the next edge NOW!
		  EdgeIter nextEdge = e;
		  nextEdge++;

		  e->isNew = false;
		  // Compute new positions associated with the vertices that will be inserted at edge midpoints,
		  // and store them in Edge::newPosition.
		  Vector3D newVertexPosition = e->newPosition;
		  VertexIter insertedVertex = mesh.splitEdge(e);
		  insertedVertex->isNew = true;
		  insertedVertex->newPosition = newVertexPosition;
		  e = nextEdge;
	  }

	  // Flip any new edge that connects an old and new vertex.
	  for (EdgeIter eIt = mesh.edgesBegin(); eIt != mesh.edgesEnd(); eIt++)
	  {
		  if (eIt->isNew)
		  {
			  if (( !eIt->halfedge()->vertex()->isNew && eIt->halfedge()->twin()->vertex()->isNew) ||
				  ( eIt->halfedge()->vertex()->isNew && !eIt->halfedge()->twin()->vertex()->isNew))
			  {
				  mesh.flipEdge(eIt);
			  }
		  }
	  }

	  for (VertexIter it = mesh.verticesBegin(); it != mesh.verticesEnd(); it++)
	  {
			it->position = it->newPosition;
	  }
  }

  void MeshResampler::downsample(HalfedgeMesh& mesh)
  {

    // TODO: (meshEdit)
    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in Face::quadric
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in Vertex::quadric
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an EdgeRecord for each edge and sticking it in the
    //    queue.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

  }

  void MeshResampler::resample(HalfedgeMesh& mesh) {

    // TODO: (meshEdit)
    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

  }

} // namespace CMU462
