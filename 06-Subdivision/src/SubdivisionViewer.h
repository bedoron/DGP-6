//=============================================================================
//
//  CLASS ParamViewer
//
//=============================================================================


#ifndef SUBDIVISION_VIEWER_HH
#define SUBDIVISION_VIEWER_HH


//== INCLUDES =================================================================

#include "MeshViewer.hh"
#include <gmm.h>

//== CLASS DEFINITION =========================================================

	      

class SubdivisionViewer : public MeshViewer
{
public:
  

   
  /// default constructor
  SubdivisionViewer(const char* _title, int _width, int _height);

  // destructor
  ~SubdivisionViewer();

  /// open mesh
  virtual bool open_mesh(const char* _meshfilename);


private:
  
  virtual void init();
  virtual void draw(const std::string& _draw_mode);

  virtual void keyboard(int key, int x, int y);

  void Perform_CatmullClark();
  void Perform_Loop();

private:
  
  //OpenMesh::VPropHandleT<OpenMesh::Vec2f>   vparam_u_, vparam_h_, texcoord_u_, texcoord_h_  ;
  OpenMesh::VPropHandleT<Mesh::Point>       vpos_;
  OpenMesh::FPropHandleT<Mesh::Point>       fcentroid_;
  OpenMesh::EPropHandleT<Mesh::Point>       enewpoint_;
  OpenMesh::EPropHandleT<VertexHandle>      enewvertex_;
  //OpenMesh::EPropHandleT<Mesh::Scalar>      eweight_;

  Mesh::Point _bbMin3D, _bbMax3D;//, _bbMin2D, _bbMax2D;

  bool Chaikin(Mesh::EdgeIter& eh,Mesh::Point& edgePoint);
  bool Chaikin(Mesh::VertexIter& vh, Mesh::Point& newPos);
  bool Chaikin(Mesh::EdgeHandle& eh,Mesh::Point& edgePoint);
  bool Chaikin(Mesh::VertexHandle& vh, Mesh::Point& newPos);
};


//=============================================================================
#endif // SUBDIVISION_VIEWER_HH defined
//=============================================================================

