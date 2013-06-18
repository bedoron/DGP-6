//=============================================================================
//
//  CLASS SubdivisionViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================
#include "SubdivisionViewer.h"
#include <windows.h>
#include <algorithm>
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
	cout << "Finding face centroids and filling the in fcentroid_\n";
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
	cout << "Finding edge newpoints and filling them in enewpoint_\n";
	for(Mesh::EdgeIter eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {

		Mesh::HalfedgeHandle heh0 = mesh_.halfedge_handle(eiter, 0);
		Mesh::HalfedgeHandle heh1 = mesh_.halfedge_handle(eiter, 1);

		Mesh::FaceHandle fh0 = mesh_.face_handle(heh0);
		Mesh::FaceHandle fh1 = mesh_.face_handle(heh1);

		Mesh::Point p0 = mesh_.point(mesh_.from_vertex_handle(heh0));
		Mesh::Point p1 = mesh_.point(mesh_.to_vertex_handle(heh0));		

		Mesh::Point fp0(0,0,0);
		Mesh::Point fp1(0,0,0);

		float divider = 2;

		if(!mesh_.is_boundary(heh0)) {
			fp0 = mesh_.property(fcentroid_, fh0);
			++divider;
		}

		if(!mesh_.is_boundary(heh1)) {
			fp1 = mesh_.property(fcentroid_, fh1);
			++divider;
		}

		Mesh::Point ep = (p0+p1+fp0+fp1)/divider;

		mesh_.property(enewpoint_, eiter) = ep;
	}

	// Calcualte to where old points need to move and store them in oldP
	cout << "Calcualte to where old points need to move and store them in oldP\n";
	std::map<Mesh::VertexHandle, Mesh::Point> oldP;
	for(Mesh::VertexIter vi = mesh_.vertices_begin(); vi != mesh_.vertices_end(); ++vi) {
		// Calculate the average face centers around VI
		Mesh::VertexHandle vh = vi;
		Mesh::Point P = mesh_.point(vi);
		Mesh::Point F(0,0,0);
		int incidents = 0;
		for(Mesh::VFIter vfIt = mesh_.vf_iter(vi); vfIt; ++vfIt) {
			++incidents;
			F += mesh_.property(fcentroid_, vfIt);
		}
		F /= ((float)incidents);

		// Calculate the average center point 
		Mesh::Point R(0,0,0);
		int n = 0;
		for(Mesh::VOHIter vohit = mesh_.voh_iter(vi); vohit; ++vohit) {
			Mesh::VertexHandle to =  mesh_.to_vertex_handle(vohit);
			Mesh::Point avg = (P + mesh_.point(to))/2.0;
			R += avg;
			++n;
		}
		R /= ((float)n);

		Mesh::Point originalMove(0,0,0);
		originalMove = (F + R*2 + P*(n-3))/((float)n);
		oldP[vh] = originalMove;
	}
	// Now move old points
	 cout << "Moving old points\n";
	for(std::map<Mesh::VertexHandle, Mesh::Point>::iterator oldPit = oldP.begin(); oldPit != oldP.end(); ++oldPit) {
		mesh_.set_point(oldPit->first, oldPit->second);
	}

	//3. Allocate new edge vertices in enewvertex_
	cout << "Allocating new vertices in enewvertex_\n";
	for(Mesh::EdgeIter eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		Mesh::Point ep = mesh_.property(enewpoint_, eiter);
		Mesh::VertexHandle vp = mesh_.add_vertex(ep);
		mesh_.property(enewvertex_, eiter) = vp;
	}
	
	// Build new faces
	cout << "Building new faces and storing them inside a map\n";
	std::vector< std::vector<Mesh::VertexHandle> > faces;
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		Mesh::Point faceCentroid = mesh_.property(fcentroid_, fiter);
		Mesh::VertexHandle centroidvh = mesh_.add_vertex(faceCentroid);	

		Mesh::FaceHandle fh = fiter.handle();

		std::map< Mesh::VertexHandle, std::vector<Mesh::VertexHandle> > new_faces;
		for(Mesh::FaceHalfedgeIter fheiter = mesh_.fh_iter(fiter); fheiter; ++fheiter) {
			Mesh::VertexHandle enewvh = mesh_.property(enewvertex_, mesh_.edge_handle(fheiter));

			Mesh::VertexHandle key0 = mesh_.from_vertex_handle(fheiter);
			Mesh::VertexHandle key1 = mesh_.to_vertex_handle(fheiter);

			bool push_centroid = false;
			if(std::find(new_faces[key0].begin(),new_faces[key0].end(),key0)==new_faces[key0].end()) {
				new_faces[key0].push_back(key0);
				push_centroid = true;
			}

			new_faces[key0].push_back(enewvh);
			if(push_centroid)
				new_faces[key0].push_back(centroidvh);

			/// -------------
			bool push_vertex = false;
			if(std::find(new_faces[key1].begin(),new_faces[key1].end(),centroidvh)==new_faces[key1].end()) {
				new_faces[key1].push_back(centroidvh);
				push_vertex = true;
			}
			new_faces[key1].push_back(enewvh);

			if(push_vertex)
				new_faces[key1].push_back(key1);

		}

		std::map< Mesh::VertexHandle, std::vector<Mesh::VertexHandle> >::iterator nfIter = new_faces.begin();
		for(;nfIter!=new_faces.end(); ++nfIter) {
			faces.push_back(nfIter->second);
		}
	}
	//4. Delete the old faces, and enter the new faces
	cout << "killing old faces *muhahaha* \n";
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		mesh_.delete_face(fiter,false);
	}

	mesh_.garbage_collection();


	for(int i = 0 ; i < faces.size(); ++i) {
		mesh_.add_face(faces[i]);
	}

	mesh_.remove_property(enewpoint_);
	mesh_.remove_property(enewvertex_);

	//mesh_ = (*qmesh);

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
	//1. Find New Edge midpoints, and enter them in enewpoint_
	for(Mesh::EdgeIter eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		Mesh::HalfedgeHandle heh0 = mesh_.halfedge_handle(eiter, 0);
		Mesh::HalfedgeHandle heh1 = mesh_.halfedge_handle(eiter, 1);

		Mesh::Point p0 = mesh_.point(mesh_.from_vertex_handle(heh0));
		Mesh::Point p1 = mesh_.point(mesh_.to_vertex_handle(heh0));

		Mesh::HalfedgeHandle topHeh = mesh_.next_halfedge_handle(heh1);
		Mesh::HalfedgeHandle bottomHeh = mesh_.next_halfedge_handle(heh0);

		Mesh::Point top = mesh_.point(mesh_.to_vertex_handle(topHeh));
		Mesh::Point bottom = mesh_.point(mesh_.to_vertex_handle(bottomHeh));

		Mesh::Point edgePoint(0,0,0);

		edgePoint = (p0 + p1)*0.375 + (top + bottom)*0.125;
		mesh_.property(enewpoint_, eiter) = edgePoint;
	}
	//2. Allocate new edge vertices in enewvertex_
	for(Mesh::EdgeIter eiter = mesh_.edges_begin(); eiter != mesh_.edges_end(); ++eiter) {
		Mesh::Point pt = mesh_.property(enewpoint_, eiter);
		Mesh::VertexHandle vh = mesh_.add_vertex(pt);
		mesh_.property(enewvertex_, eiter) = vh;
	}
	//3. Delete the old faces, and enter the new faces
	std::map<Mesh::VertexHandle, Mesh::Point> originals;
	for(Mesh::VertexIter vh = mesh_.vertices_begin(); vh != mesh_.vertices_end(); ++vh) {
		Mesh::Point P = mesh_.point(vh);
		int valence = 0;
		Mesh::Point neighbors(0,0,0);
		for(Mesh::VVIter vvit = mesh_.vv_iter(vh); vvit ; ++vvit) {
			++valence;
			neighbors += mesh_.point(vvit);
		}
		// Build Beta
		float square = ((3.0/8.0)+0.25*cos(2.8*M_PI/(float)valence));
		float beta = (1.0/(float)valence)*((5.0/8.0)-(square*square));

		originals[vh.handle()] = P*(1-((float)valence)*beta) + neighbors*beta;
//		originals[vh.handle()] = Mesh::Point(0,0,0);//P*(1-((float)valence)*beta) + neighbors*beta;
	}
	// Move Existing vertices
	for(std::map<Mesh::VertexHandle, Mesh::Point>::iterator oit = originals.begin(); oit != originals.end(); ++oit) {
		mesh_.set_point(oit->first, oit->second);
	}
	// Build faces
	std::vector< std::vector<Mesh::VertexHandle> > faces;
	for(Mesh::FaceIter fit = mesh_.faces_begin(); fit != mesh_.faces_end(); ++fit) {

		for(Mesh::FaceHalfedgeIter fheit = mesh_.fh_iter(fit); fheit; ++fheit ) {

			Mesh::HalfedgeHandle prev, tentative;
			prev = mesh_.next_halfedge_handle(mesh_.next_halfedge_handle(fheit));
			tentative = mesh_.next_halfedge_handle(prev);

			prev = (tentative == fheit.handle())?prev:tentative;

			Mesh::EdgeHandle currentEdge = mesh_.edge_handle(fheit);			
			Mesh::EdgeHandle prevEdge = mesh_.edge_handle(prev);

			Mesh::VertexHandle v0, currCenter, prevCenter;

			v0 = mesh_.from_vertex_handle(fheit);
			currCenter = mesh_.property(enewvertex_, currentEdge);
			prevCenter = mesh_.property(enewvertex_, prevEdge);

			std::vector< Mesh::VertexHandle > face;
			face.push_back(v0);
			face.push_back(currCenter);
			face.push_back(prevCenter);

			faces.push_back(face);
		}
	}

	//4. Delete the old faces, and enter the new faces
	cout << "killing old faces *muhahaha* \n";
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		mesh_.delete_face(fiter,false);
	}

	mesh_.garbage_collection();


	for(int i = 0 ; i < faces.size(); ++i) {
		mesh_.add_face(faces[i]);
	}

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
