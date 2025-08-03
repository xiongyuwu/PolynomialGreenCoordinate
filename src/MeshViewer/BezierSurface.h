#pragma once
#include <QString>
#include <QMessageBox>
#include <QFileDialog>

#include <Eigen/Dense>

#include "QGLViewerWidget.h"
#include "MeshDefinition.h"

class BezierSurface {
public:
	std::vector<Mesh::Point> ControlPoints = { {-1.0, -1.0, 0.0},{1.0, -1.0, 0.0},{-1.0, 1.0, 0.0},{1.0, 1.0, 0.0} };
	bool cage_mode = false;
	int control_point_per_edge;

	double step_length = 0.01;
	Mesh::Point ComputeSurfacePoint(double u, double v);
	void DrawBezierSurface();
	void SetControlPoints(const std::vector<Mesh::Point>& control_points, const bool cage_mode_, const int control_point_per_edge_) {
		ControlPoints = control_points;
		cage_mode = cage_mode_;
		control_point_per_edge = control_point_per_edge_;
	};
	int factorial(int n);
};
