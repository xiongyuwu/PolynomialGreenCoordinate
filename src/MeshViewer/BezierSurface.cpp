#include <iostream>
#include <algorithm>
#include <cmath>

#include <qapplication.h>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "BezierSurface.h"
#include "../Common/CommonDefinitions.h"

Mesh::Point BezierSurface::ComputeSurfacePoint(double u, double v) {
	Mesh::Point result = { 0,0,0 };
	//quad cage
	if (cage_mode == false) {
		int idx = 0;
		for (int i = 0; i < control_point_per_edge; ++i) {
			for (int j = 0; j < control_point_per_edge; ++j) {
				int factor = factorial(control_point_per_edge - 1) / factorial(i) /
					factorial(control_point_per_edge - 1 - i)*factorial(control_point_per_edge - 1) /
					factorial(j) / factorial(control_point_per_edge - 1 - j);
				result += std::pow(u, j)*std::pow(1 - u, control_point_per_edge - 1 - j)*std::pow(v, i)*
					std::pow(1 - v, control_point_per_edge - 1 - i)*factor*ControlPoints[idx];
				idx++;
			}
		}
		/*if (ControlPoints.size() == 4) {
			result = (1 - u) * (1 - v) * ControlPoints[0] + (1 - u) * v * ControlPoints[1] + u * v * ControlPoints[2] + u * (1 - v) * ControlPoints[3];
		}*/
	}
	//tri cage
	else {
		double w = 1 - u - v;
		int idx = 0;
		for (int i = 0; i < control_point_per_edge; ++i) {
			for (int j = 0; j < control_point_per_edge - i; ++j) {
				int k = control_point_per_edge - 1 - i - j;
				int factor = factorial(control_point_per_edge - 1) / factorial(i) / factorial(j) / factorial(k);
				result += std::pow(u, j)*std::pow(v, i)*std::pow(w, k)*factor*ControlPoints[idx];
				idx++;
			}
		}
		/*double w = 1 - u - v;
		if (ControlPoints.size() == 10) {
			result += ControlPoints[0] * w*w*w;
			result += ControlPoints[1] * 3 * u*w*w;
			result += ControlPoints[2] * 3 * u*u*w;
			result += ControlPoints[3] * u*u*u;
			result += ControlPoints[4] * 3 * v*w*w;
			result += ControlPoints[5] * 6 * u*v*w;
			result += ControlPoints[6] * 3 * v*u*u;
			result += ControlPoints[7] * 3 * v*v*w;
			result += ControlPoints[8] * 3 * v*v*u;
			result += ControlPoints[9] * v*v*v;
		}*/
	}
	return result;
}

void BezierSurface::DrawBezierSurface() {
	//draw facet
	glDepthMask(GL_FALSE);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	GLfloat global_ambient[] = { 0.8f, 0.8f, 0.8f, 1.0f }; // 全局环境光
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

	// 关闭光源的漫反射和镜面反射，只保留环境光
	GLfloat light_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };  // 环境光关闭
	GLfloat light_diffuse[] = { 0.0f, 0.0f, 0.0f, 1.0f };  // 漫反射关闭
	GLfloat light_specular[] = { 0.0f, 0.0f, 0.0f, 1.0f }; // 镜面反射关闭
	GLfloat light_position[] = { 0.0f, 0.0f, 1.0f, 0.0f }; // 光源方向（无效）

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	double r = 150;
	double g = 150;
	double b = 250;
	GLfloat mat_a[] = { r / 255,g / 255, b / 255, 0.1f }; // 环境光
	GLfloat mat_d[] = { r / 255 * 0.5,g / 255 * 0.5, b / 255 * 0.5, 0.1f }; // 漫反射
	GLfloat mat_s[] = { 0, 0, 0, 0.1f };  // 镜面反射
	GLfloat shine[] = { 0.0f }; // 光泽度

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);

	if (cage_mode == false) {
		glBegin(GL_QUADS);
		for (double u = 0; u < 1 - step_length * 0.1; u += step_length) {
			for (double v = 0; v < 1 - step_length * 0.1; v += step_length) {
				Mesh::Point p1 = ComputeSurfacePoint(u, v);
				Mesh::Point p2 = ComputeSurfacePoint(u + step_length, v);
				Mesh::Point p3 = ComputeSurfacePoint(u, v + step_length);
				Mesh::Point p4 = ComputeSurfacePoint(u + step_length, v + step_length);

				glVertex3dv(p1.data());
				glVertex3dv(p2.data());
				glVertex3dv(p4.data());
				glVertex3dv(p3.data());
			}
		}
		glEnd();
	}
	else {
		glBegin(GL_TRIANGLES);
		for (double u = 0; u < 1 - step_length * 0.1; u += step_length) {
			for (double v = 0; v < 1 - u - step_length * 0.1; v += step_length) {
				Mesh::Point p1 = ComputeSurfacePoint(u, v);
				Mesh::Point p2 = ComputeSurfacePoint(u + step_length, v);
				Mesh::Point p3 = ComputeSurfacePoint(u, v + step_length);

				glVertex3dv(p1.data());
				glVertex3dv(p2.data());
				glVertex3dv(p3.data());

				if (u + v + step_length * 1.1 <= 1) {
					Mesh::Point p4 = ComputeSurfacePoint(u + step_length, v + step_length);
					glVertex3dv(p2.data());
					glVertex3dv(p3.data());
					glVertex3dv(p4.data());
				}
			}
		}
		glEnd();
	}

	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_LIGHTING);

	//draw warframe
	glLineWidth(3);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	r = 150 / 2;
	g = 150 / 2;
	b = 250 / 2;
	glColor4f(r / 255, g / 255, b / 255, 0.3f);
	if (cage_mode == false) {
		/*glBegin(GL_QUADS);
		for (double u = 0; u < 1 - step_length * 0.1; u += step_length) {
			for (double v = 0; v < 1 - step_length * 0.1; v += step_length) {
				//glBegin(GL_POLYGON);
				Mesh::Point p1 = ComputeSurfacePoint(u, v);
				Mesh::Point p2 = ComputeSurfacePoint(u + step_length, v);
				Mesh::Point p3 = ComputeSurfacePoint(u, v + step_length);
				Mesh::Point p4 = ComputeSurfacePoint(u + step_length, v + step_length);

				glVertex3dv(p1.data());
				glVertex3dv(p2.data());
				glVertex3dv(p4.data());
				glVertex3dv(p3.data());
				//glEnd();
			}
		}
		glEnd();*/
		glBegin(GL_LINE_STRIP);
		for (double u = 0; u < 1 - step_length * 0.1; u += step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(u, 0);
			glVertex3dv(p1.data());
		}
		for (double v = 0; v < 1; v += step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(1, v);
			glVertex3dv(p1.data());
		}
		for (double u = 1; u > 0; u -= step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(u, 1);
			glVertex3dv(p1.data());
		}
		for (double v = 1; v >= 0; v -= step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(0, v);
			glVertex3dv(p1.data());
		}
		glEnd();
	}
	else {
		/*glBegin(GL_TRIANGLES);
		for (double u = 0; u < 1 - step_length * 0.1; u += step_length) {
			for (double v = 0; v < 1 - u - step_length * 0.1; v += step_length) {
				Mesh::Point p1 = ComputeSurfacePoint(u, v);
				Mesh::Point p2 = ComputeSurfacePoint(u + step_length, v);
				Mesh::Point p3 = ComputeSurfacePoint(u, v + step_length);

				glVertex3dv(p1.data());
				glVertex3dv(p2.data());
				glVertex3dv(p3.data());

				if (u + v + step_length * 1.1 <= 1) {
					Mesh::Point p4 = ComputeSurfacePoint(u + step_length, v + step_length);
					glVertex3dv(p2.data());
					glVertex3dv(p3.data());
					glVertex3dv(p4.data());
				}
			}
		}*/
		glBegin(GL_LINE_STRIP);
		for (double u = 0; u < 1 - step_length * 0.1; u += step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(u, 0);
			glVertex3dv(p1.data());
		}
		for (double u = 1; u > 0; u -= step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(u, 1 - u);
			glVertex3dv(p1.data());
		}
		for (double v = 1; v >= 0; v -= step_length) {
			Mesh::Point p1 = ComputeSurfacePoint(0, v);
			glVertex3dv(p1.data());
		}
		glEnd();
	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDepthMask(GL_TRUE);
}

int BezierSurface::factorial(int n) {
	if (n <= 1) return 1;
	return n * factorial(n - 1);
}