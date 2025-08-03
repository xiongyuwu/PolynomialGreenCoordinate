#ifndef MESHPROCESSING_MAIN_VIEWWE_WIDGET_H
#define MESHPROCESSING_MAIN_VIEWWE_WIDGET_H

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <QtGui>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
//main widget
#include "InteractiveViewerWidget.h"
#include "MeshParamDialog.h"

class MainViewerWidget : public QDialog
{
	Q_OBJECT
public:
	MainViewerWidget(QWidget* _parent = 0);
	~MainViewerWidget();

	void setDrawMode(int dm)
	{
		MeshViewer->setDrawMode(dm);
	}
	void setMouseMode(int dm)
	{
		MeshViewer->setMouseMode(dm);
	}

	void set_show_BBox()
	{
		MeshViewer->set_draw_bbox_ok();
	}
	void set_show_mesh_boundary()
	{
		MeshViewer->set_draw_mesh_boundary_ok();
	}
	void change_cage_mode()
	{
		MeshViewer->PGC.set_change_cage_mode();
	}
	void set_show_control_point()
	{
		MeshViewer->set_draw_control_point_ok();
	}
	void set_show_move_vertices()
	{
		MeshViewer->set_draw_move_vertices_ok();
	}
	void set_show_cage()
	{
		MeshViewer->set_draw_cage_ok();
	}
	void generate_deformation()
	{
		MeshViewer->PGC.do_deformation();
		MeshViewer->setmesh(MeshViewer->PGC.getmesh());
	}
	void generate_variational_deformation()
	{
		MeshViewer->PGC.do_variational_deformation();
		MeshViewer->setmesh(MeshViewer->PGC.getmesh());
	}
	void generate_deformation_by_control_points()
	{
		Mesh control_points_mesh;
		OpenMesh::IO::read_mesh(control_points_mesh, "deformed_control_points.obj");
		std::vector<Mesh::Point> represent_control_points;
		for (auto& vh : control_points_mesh.vertices())
			represent_control_points.push_back(control_points_mesh.point(vh));

		MeshViewer->PGC.set_represent_control_points(represent_control_points);
		MeshViewer->PGC.control_point_is_initialed = true;
		MeshViewer->PGC.do_deformation();
		MeshViewer->setmesh(MeshViewer->PGC.getmesh());
	}
	void optimize_control_points()
	{
		MeshViewer->PGC.do_optimize_control_points();
		MeshViewer->setmesh(MeshViewer->PGC.getmesh());
	}
	void record_origin_mesh()
	{
		std::string prefix;
		std::array<bool, 3> compute_level;
		std::vector<Mesh::Point>sample_points;
		prefix = "";
		compute_level = { true,false,false };
		auto& mesh = MeshViewer->mesh_ref();
		for (const auto& vh : mesh.vertices())
			sample_points.push_back(mesh.point(vh));
		MeshViewer->PGC.compute_integral_para_for_sample_points(prefix, compute_level, sample_points);
	}
	void read_integral_para()
	{
		MeshViewer->PGC.do_read_integral_para();
	}
	void read_cage()
	{
		MeshViewer->PGC.do_read_cage();
	}
	void compute_para_for_control_points()
	{
		std::string prefix;
		std::array<bool, 3> compute_level;
		std::vector<Mesh::Point>sample_points;
		prefix = "";
		compute_level = { true,true,true };
		auto& mesh = MeshViewer->mesh_ref();
		for (const auto& vh : mesh.vertices())
			sample_points.push_back(mesh.point(vh));
		MeshViewer->PGC.do_compute_para_for_control_points(prefix, compute_level, sample_points);
	}
	void generate_para_complete_flow()
	{
		MeshViewer->PGC.do_generate_para_complete_flow();
		MeshViewer->setmesh(MeshViewer->PGC.getmesh());
	}
	void clear_select_points() {
		MeshViewer->clear_select_points();
	}
	void openMesh_fromMain(char* filename)
	{
		QString str(filename);
		open_mesh_gui(filename);
	}

	void edit_undo()
	{
		MeshViewer->edit_undo_viewer();
	}

	void edit_redo()
	{
		MeshViewer->edit_redo_viewer();
	}

public slots:
	void open_mesh_query()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open mesh file"),
			tr("../models/"),
			tr("OFF Files (*.off);;"
				"OBJ Files (*.obj);;"
				"PLY Files (*.ply);;"
				"STL Files (*.stl);;"
				"All Files (*)"));
		if (!fileName.isEmpty())
		{
			open_mesh_gui(fileName);
		}
	}
	void save_mesh_query()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save mesh file"),
			tr("../models/untitled.off"),
			tr("OFF Files (*.off);;"
				"OBJ Files (*.obj);;"
				"PLY Files (*.ply);;"
				"STL Files (*.stl);;"
				"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_mesh_gui(fileName);
		}
	}
	void saveOpenGLScreen()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			("Save screen as image file"),
			("../Results/untitled.png"),
			("PNG Files (*.png);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
				"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_screen_gui(fileName);
		}
	}
	void save_opengl_screen(QString str)
	{
		MeshViewer->saveScreen(str.toLocal8Bit());
	}
	virtual void update_mesh()
	{
		if (MeshViewer->mesh_ref().n_vertices() != 0)
		{
			MeshViewer->updateMesh();
		}
	}

	virtual void clear_all_mesh()
	{
		if (LoadMeshSuccess)
		{
			LoadMeshSuccess = false;
			MeshViewer->clearAllMesh();
		}
	}

	virtual void clear_all_selected()
	{
		if (LoadMeshSuccess)
		{
			MeshViewer->clearSelectedData();
			MeshViewer->updateGL();
		}
	}

	void LoadMeshFromInner(bool OK, QString fname)
	{
		LoadMeshSuccess = OK;
		if (LoadMeshSuccess)
		{
			SetMeshForALL();
		}
		emit(haveLoadMesh(fname));
	};

	//
	void print_info();

signals:
	void haveLoadMesh(QString filePath);
	void setMouseMode_signal_main(int);
	void setDrawMode_signal_main(int);

	void set_edit_undo_enable_signal(bool);
	void set_edit_redo_enable_signal(bool);

protected:
	virtual void initViewerWindow();
	virtual void createParamIDialog();
	virtual void createViewerDialog();
	virtual void save_mesh_gui(QString fname);
	virtual void open_mesh_gui(QString fname);
	virtual void save_screen_gui(QString fname);

protected:
	bool LoadMeshSuccess;

private:
	InteractiveViewerWidget* MeshViewer;
	MeshParamDialog* MeshParam;

	void SetMeshForALL()
	{
	}

#pragma region Auxiliary Function
public:
	void aux_inverse_mesh_connectivity()
	{
		MeshViewer->inverse_mesh_connectivity();
	}

	void aux_scale_mesh_using_BBox(int max_len)
	{
		MeshViewer->scale_mesh_using_BBox(max_len);
	}

	void aux_split_quad_mesh()
	{
		MeshViewer->split_quad_mesh();
	}

	void transform_mesh(const std::vector<double>& m)
	{
		MeshViewer->transform_mesh(m);
	}

	void aux_find_vertex_by_id(int id)
	{
		MeshViewer->find_vertex_by_id(id);
	}

	void aux_find_face_by_id(int id)
	{
		MeshViewer->find_face_by_id(id);
	}

	void aux_find_edge_by_id(int id)
	{
		MeshViewer->find_edge_by_id(id);
	}

	void aux_find_vertex_by_valance(int valance)
	{
		MeshViewer->find_vertex_by_valance(valance);
	}

	void aux_delete_vertex_valence_four()
	{
		MeshViewer->delete_vertex_valence_four();
	}

	void aux_delete_vertex_valence_three()
	{
		MeshViewer->delete_vertex_valence_three();
	}

	void aux_split_vertex_valence_eight()
	{
		MeshViewer->split_vertex_valence_eight();
	}

#pragma endregion

};


#endif