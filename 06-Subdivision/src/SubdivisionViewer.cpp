//=============================================================================
//
//  CLASS SubdivisionViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================
#include "SubdivisionViewer.h"
#include <windows.h>

//== IMPLEMENTATION ========================================================== 

static bool _SubdivisioneterizationComputed_u = false, _SubdivisioneterizationComputed_h = false;
static bool _BoundingBox2DComputed = false;


SubdivisionViewer::
	SubdivisionViewer(const char* _title, int _width, int _height)
	: MeshViewer(_title, _width, _height)
{ 

	mesh_.request_face_status();
	mesh_.request_halfedge_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	init();
}

SubdivisionViewer::
	~SubdivisionViewer()
{ 
}

void
	SubdivisionViewer::
	init()
{
	// base class first
	MeshViewer::init();

}

//-----------------------------------------------------------------------------
typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;
typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;

static Mesh *qmesh = 0;
static TriMesh *tmesh = 0;

void SubdivisionViewer::Perform_CatmullClark()
{
	mesh_.add_property(fcentroid_);
	mesh_.add_property(enewpoint_);
	mesh_.add_property(enewvertex_);


	/***********Catmull-Clark Subdivision************
	1. Find face centroids and fill them in fcentroid_
	2. Find edge newpoints and fill them in enewpoint_
	3. Allocate new edge vertices in enewvertex_
	4. Delete the old faces, and enter the new faces

	NOTICE: make sure to call mesh_.garbage_collection after deleting new faces, or you
	might get mixed up when running over all faces\halfedges\vertices
	**********************************************************/
	//1. Find face centroids and fill them in fcentroid_
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		Mesh::Point p(0,0,0);
		
		int total_points = 0;
		for(Mesh::FVIter fviter = mesh_.fv_iter(fiter); fviter; ++fviter) {
			++total_points;
			p += mesh_.point(fviter);
		}
		p /= total_points;

		mesh_.property(fcentroid_, fiter) = p;
	}
	//2. Find edge newpoints and fill them in enewpoint_
	for(Mesh::EdgeIter eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		Mesh::HalfedgeHandle heh0 = mesh_.halfedge_handle(eiter, 0);
		Mesh::HalfedgeHandle heh1 = mesh_.halfedge_handle(eiter, 0);

		Mesh::FaceHandle fh0 = mesh_.face_handle(heh0);
		Mesh::FaceHandle fh1 = mesh_.face_handle(heh1);

		Mesh::Point fp0 = mesh_.property(fcentroid_, fh0);
		Mesh::Point fp1 = mesh_.property(fcentroid_, fh1);

		Mesh::Point p0 = mesh_.point(mesh_.from_vertex_handle(heh0));
		Mesh::Point p1 = mesh_.point(mesh_.to_vertex_handle(heh0));

		Mesh::Point ep = (p0+p1+fp0+fp1)/4.0;

		mesh_.property(enewpoint_, eiter) = ep;
	}
	//3. Allocate new edge vertices in enewvertex_
	for(Mesh::EdgeIter eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		Mesh::Point ep = mesh_.property(enewpoint_, eiter);
		Mesh::VertexHandle vp = mesh_.add_vertex(ep);
		mesh_.property(enewvertex_, eiter) = vp;
	}
	
	if(qmesh != 0)
		delete qmesh;
	qmesh = new Mesh(mesh_);
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		qmesh->delete_face(fiter, false);
	}
	// qmesh now contains only vertices, connect them by the rules:
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		Mesh::Point faceCentroid = mesh_.property(fcentroid_, fiter);
		Mesh::FaceEdgeIter feiter = mesh_.fe_iter(fiter);
		std::vector<Mesh::VertexHandle> edgeCenters;
		for(;feiter; ++feiter) {
			edgeCenters.push_back(mesh_.property(enewvertex_, feiter));
		}
	}	
	

	//4. Delete the old faces, and enter the new faces

	mesh_.remove_property(fcentroid_);
	mesh_.remove_property(enewpoint_);
	mesh_.remove_property(enewvertex_);

	mesh_.update_normals();

	update_face_indices();

	std::cerr << mesh_.n_vertices() << " vertices, "
		<< mesh_.n_faces()    << " faces\n";
}

void SubdivisionViewer::Perform_Loop()
{
	mesh_.add_property(enewpoint_);
	mesh_.add_property(enewvertex_);


	/************Loop subdivision***************
	1. Find New Edge midpoints, and enter them in enewpoint_
	2. Allocate new edge vertices in enewvertex_
	3. Delete the old faces, and enter the new faces
	**********************************************/

	mesh_.remove_property(enewpoint_);
	mesh_.remove_property(enewvertex_);

	mesh_.update_normals();

	update_face_indices();

	std::cerr << mesh_.n_vertices() << " vertices, "
		<< mesh_.n_faces()    << " faces\n";
}

//-----------------------------------------------------------------------------





bool
	SubdivisionViewer::
	open_mesh(const char* _meshfilename)
{
	// load mesh
	if (MeshViewer::open_mesh(_meshfilename))
	{
		// store vertex initial positions and 3D mesh bounding box
		Mesh::VertexIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end());
		_bbMin3D = _bbMax3D = mesh_.point(v_it);
		for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
		{
			//mesh_.property(vpos_,v_it) = mesh_.point(v_it);
			_bbMin3D.minimize(mesh_.point(v_it));
			_bbMax3D.maximize(mesh_.point(v_it));

		}
		return true;
	}
	return false;
}


void 
	SubdivisionViewer::
	draw(const std::string& _draw_mode)
{
	MeshViewer::draw(_draw_mode);
}


//-----------------------------------------------------------------------------

void
	SubdivisionViewer::
	keyboard(int key, int x, int y)
{
	switch (toupper(key))
	{ 
	case 'O':
		{
			OPENFILENAME ofn={0};
			char szFileName[MAX_PATH]={0};
			ofn.lStructSize=sizeof(OPENFILENAME);
			ofn.Flags=OFN_ALLOWMULTISELECT|OFN_EXPLORER;
			ofn.lpstrFilter="All Files (*.*)\0*.*\0";
			ofn.lpstrFile=szFileName;
			ofn.nMaxFile=MAX_PATH;
			if(GetOpenFileName(&ofn))
				open_mesh(szFileName);
		}
		break;
	case 'C':
		Perform_CatmullClark();
		glutPostRedisplay();
		break;
	case 'L':
		Perform_Loop();
		glutPostRedisplay();
		break;
	default:
		MeshViewer::keyboard(key, x, y);
		break;
	}
}
//-----------------------------------------------------------------------------

//=============================================================================
