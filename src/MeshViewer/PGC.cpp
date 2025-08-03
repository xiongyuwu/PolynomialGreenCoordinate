#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include <queue>
#include <unordered_set>
#include <iomanip>

#include <qapplication.h>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "PGC.h"

PGC::PGC(const Mesh& mesh_) :mesh(mesh_) {
	cage_mode = false;
	control_point_size_per_facet = control_point_per_edge * control_point_per_edge;
	updateControlPoints();
}

void PGC::draw_represent_control_points() {
	glColor3f(0, 0, 0);
	glPointSize(control_point_size_per_facet);
	glBegin(GL_POINTS);
	for (const auto& point : represent_control_points) {
		glVertex3d(point[0], point[1], point[2]);
	}
	glEnd();
	//std::cout << "represent_control_points.size()=" << represent_control_points.size() << "\n";
}

void PGC::draw_cage() {
	if (cage_mode == true) {
		BezierSurface facet;
		for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
			std::vector<Mesh::Point> control_points_on_facet(control_points.begin() + facet_idx * control_point_size_per_facet,
				control_points.begin() + (facet_idx + 1) * control_point_size_per_facet);
			facet.SetControlPoints(control_points_on_facet, cage_mode, control_point_per_edge);
			facet.DrawBezierSurface();
		}
	}
	else {
		BezierSurface facet;
		for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
			std::vector<Mesh::Point> control_points_on_facet(control_points.begin() + facet_idx * control_point_size_per_facet,
				control_points.begin() + (facet_idx + 1) * control_point_size_per_facet);
			facet.SetControlPoints(control_points_on_facet, cage_mode, control_point_per_edge);
			facet.DrawBezierSurface();
		}
	}
}

void PGC::draw_position_constraints() {
	read_position_constraints();

	glColor3f(1.0f, 0.5f, 0.0f);
	double radius = cage_mesh_diag_length*0.01;
	for (auto& position_constraint : position_constraints) {
		auto point = mesh.point(mesh.vertex_handle(position_constraint.first));
		glPushMatrix();
		glTranslatef(point[0], point[1], point[2]);
		GLUquadric *quadric = gluNewQuadric();
		gluSphere(quadric, radius, 30, 30);
		gluDeleteQuadric(quadric);
		glPopMatrix();
	}
	glColor3f(0.0f, 0.0f, 1.0f);
	for (auto& fixed_vertex : fixed_vertices) {
		auto point = mesh.point(mesh.vertex_handle(fixed_vertex));
		glPushMatrix();
		glTranslatef(point[0], point[1], point[2]);
		GLUquadric *quadric = gluNewQuadric();
		gluSphere(quadric, radius, 30, 30);
		gluDeleteQuadric(quadric);
		glPopMatrix();
	}
}

void PGC::do_generate_para_complete_flow() {
	generate_smooth_sample_points();

	read_position_constraints();
	read_orientation_constraints();
	read_rigidity_constraints();
	read_smoothness_constraints();

	std::string prefix;
	std::array<bool, 3> compute_level;
	std::vector<Mesh::Point>sample_points;
	prefix = "";
	compute_level = { true,true,true };
	for (const auto& vh : mesh.vertices())
		sample_points.push_back(mesh.point(vh));
	compute_integral_para_for_sample_points(prefix, compute_level, sample_points);

	prefix = "rigidity_";
	compute_level = { true ,true,false };
	compute_integral_para_for_sample_points(prefix, compute_level, rigidity_constraints);

	sample_points.clear();
	prefix = "smoothness_";
	compute_level = { true,false,true };
	compute_integral_para_for_sample_points(prefix, compute_level, smoothness_constraints);

	constraint_is_set = true;
	do_optimize_control_points();
}

void PGC::do_read_cage() {
	origin_cage_is_initialed = true;
	bool read_ok = OpenMesh::IO::read_mesh(cage_mesh, "cage.obj");
	if (!read_ok) {
		std::cout << "read cage error!\n";
	}
	origin_cage_vertices.clear();
	map_facet_to_point_idxs.clear();

	OpenMesh::Vec3d bbox_min(DBL_MAX, DBL_MAX, DBL_MAX), bbox_max(-DBL_MAX,-DBL_MAX,-DBL_MAX);
	for (const auto& vh : cage_mesh.vertices()) {
		auto& point = cage_mesh.point(vh);
		origin_cage_vertices.push_back(point);
		for (int i = 0; i < 3; ++i) {
			if (point[i] > bbox_max[i])
				bbox_max[i] = point[i];
			if (point[i] < bbox_min[i])
				bbox_min[i] = point[i];
		}
	}
	cage_mesh_diag_length = (bbox_max - bbox_min).length();

	for (const auto&fh : cage_mesh.faces()) {
		std::vector<size_t> vertices;
		for (auto& it = cage_mesh.cfv_iter(fh); it.is_valid(); ++it) {
			vertices.push_back((*it).idx());
		}
		map_facet_to_point_idxs.push_back(vertices);
	}

	std::cout << "read cage complete!\n";
	updateControlPoints();
	control_point_is_initialed = true;
	std::ofstream represent_control_points_out("control_points.obj");
	for (auto& point : represent_control_points) {
		represent_control_points_out << "v " << point << "\n";
	}
	represent_control_points_out.close();
}

void PGC::updateControlPoints() {
	if (!control_point_is_initialed) {
		control_points.clear();
	}
	facet_and_control_points_pairs.clear();
	if (cage_mode == false) {
		if (!origin_cage_is_initialed) {
			origin_cage_vertices = { {-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1} };
			map_facet_to_point_idxs = { {0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7} };
		}

		for (int facet_idx = 0; facet_idx < map_facet_to_point_idxs.size(); ++facet_idx) {
			auto& facet_point_idxs = map_facet_to_point_idxs[facet_idx];
			auto& P0 = origin_cage_vertices[facet_point_idxs[0]];
			auto& P1 = origin_cage_vertices[facet_point_idxs[1]];
			auto& P2 = origin_cage_vertices[facet_point_idxs[2]];
			auto& P3 = origin_cage_vertices[facet_point_idxs[3]];
			if (!control_point_is_initialed) {
				//默认平行四边形
				auto u_vec = P1 - P0;
				auto v_vec = P3 - P0;
				for (int i = 0; i < control_point_per_edge; ++i) {
					for (int j = 0; j < control_point_per_edge; ++j) {
						control_points.push_back(P0 + (j * u_vec + i * v_vec) / (control_point_per_edge - 1));
					}
				}
				//control_points.push_back(origin_cage_vertices[facet_point_idxs[0]]);
				//control_points.push_back(origin_cage_vertices[facet_point_idxs[1]]);
				//control_points.push_back(origin_cage_vertices[facet_point_idxs[2]]);
				//control_points.push_back(origin_cage_vertices[facet_point_idxs[3]]);
			}

			std::vector<Mesh::Point> facet_vertices = { P0,P1,P2,P3 };
			std::vector<Mesh::Point> control_points_on_facet(control_points.begin() + facet_idx * control_point_size_per_facet,
				control_points.begin() + (facet_idx + 1) * control_point_size_per_facet);
			facet_and_control_points_pairs.push_back({ facet_vertices,facet_vertices });
		}
	}
	if (cage_mode == true) {
		if (!origin_cage_is_initialed) {
			double sqrt_2 = std::sqrt(2);
			double sqrt_6 = std::sqrt(6);
			origin_cage_vertices = { {sqrt_2,-sqrt_6,-1},{sqrt_2,sqrt_6,-1},{-2 * sqrt_2,0,-1},{0,0,3} };
			map_facet_to_point_idxs = { {0,2,1},{0,1,3},{1,2,3},{2,0,3} };
		}

		for (int facet_idx = 0; facet_idx < map_facet_to_point_idxs.size(); ++facet_idx) {
			auto& facet_point_idxs = map_facet_to_point_idxs[facet_idx];
			auto& P003 = origin_cage_vertices[facet_point_idxs[0]];
			auto& P030 = origin_cage_vertices[facet_point_idxs[1]];
			auto& P300 = origin_cage_vertices[facet_point_idxs[2]];

			if (!control_point_is_initialed) {
				for (int i = 0; i < control_point_per_edge; ++i) {
					for (int j = 0; j < control_point_per_edge - i; ++j) {
						control_points.push_back((i*P300 + j * P030 + (control_point_per_edge - 1 - i - j)*P003)
							/ (control_point_per_edge - 1));
					}
				}
				//auto& P012 = 2.0 / 3 * P003 + 1.0 / 3 * P030;
				//auto& P021 = 1.0 / 3 * P003 + 2.0 / 3 * P030;
				//auto& P102 = 2.0 / 3 * P003 + 1.0 / 3 * P300;
				//auto& P111 = 1.0 / 3 * P300 + 1.0 / 3 * P030 + 1.0 / 3 * P003;
				//auto& P120 = 2.0 / 3 * P030 + 1.0 / 3 * P300;
				//auto& P201 = 2.0 / 3 * P300 + 1.0 / 3 * P003;
				//auto& P210 = 2.0 / 3 * P300 + 1.0 / 3 * P030;

				//control_points.push_back(P003);
				//control_points.push_back(P012);
				//control_points.push_back(P021);
				//control_points.push_back(P030);
				//control_points.push_back(P102);
				//control_points.push_back(P111);
				//control_points.push_back(P120);
				//control_points.push_back(P201);
				//control_points.push_back(P210);
				//control_points.push_back(P300);
			}

			std::vector<Mesh::Point> facet_vertices = { P003,P030,P300 };
			std::vector<Mesh::Point> control_points_on_facet(control_points.begin() + facet_idx * control_point_size_per_facet,
				control_points.begin() + (facet_idx + 1) * control_point_size_per_facet);
			facet_and_control_points_pairs.push_back({ facet_vertices,control_points_on_facet });
		}
	}
	if (!control_point_is_initialed) {
		processPoints();
	}
}

void PGC::generate_smooth_sample_points() {
	double retract_length = DBL_MAX;
	for (auto& eh : cage_mesh.edges()) {
		auto& heh = cage_mesh.halfedge_handle(eh);
		auto& p0 = cage_mesh.point(cage_mesh.to_vertex_handle(heh));
		auto& p1 = cage_mesh.point(cage_mesh.from_vertex_handle(heh));
		double length = (p0 - p1).length();
		if (retract_length > length)
			retract_length = length;
	}
	retract_length /= 10;
	std::cout << "retract_length:" << retract_length << "\n";
	std::vector<Mesh::Point> retract_vertices;
	std::vector<Mesh::Point> face_normals;
	for (auto& fh : cage_mesh.faces()) {
		face_normals.push_back(cage_mesh.calc_face_normal(fh));
	}

	std::ofstream smoothness_sample_points_out("smoothness_sample_points.txt");
	std::ofstream smoothness_sample_points("smoothness_sample_points.obj");

	for (const auto& vh : cage_mesh.vertices()) {
		Mesh::Point retract_direction = { 0,0,0 };
		std::vector<Mesh::Point> normals;
		for (auto& fh : cage_mesh.vf_range(vh)) {
			normals.push_back(face_normals[fh.idx()]);
		}
		retract_direction = find_positive_inner_product_vector(normals);

		retract_vertices.push_back(cage_mesh.point(vh) - retract_direction * retract_length);
		if (retract_direction.length() < 1e-10)
			continue;
		smoothness_sample_points_out << retract_vertices[retract_vertices.size() - 1] << "\n";
		smoothness_sample_points << "v " << retract_vertices[retract_vertices.size() - 1] << "\n";
	}

	int segment_number = 7;
	for (const auto& eh : cage_mesh.edges()) {
		auto& heh = cage_mesh.halfedge_handle(eh, 0);
		auto& vh1 = cage_mesh.to_vertex_handle(heh);
		auto& vh2 = cage_mesh.from_vertex_handle(heh);

		for (int i = 1; i < segment_number; ++i) {
			auto& point = retract_vertices[vh1.idx()] * i / segment_number + retract_vertices[vh2.idx()] * (segment_number - i) / segment_number;
			bool valid_flag = true;
			for (auto& fh : cage_mesh.vf_range(vh1)) {
				double distance = std::abs((point - cage_mesh.point(vh1)).dot(face_normals[fh.idx()]));
				if (distance < retract_length*0.25) {
					valid_flag = false;
					break;
				}
			}
			if (!valid_flag)continue;

			valid_flag = true;
			for (auto& fh : cage_mesh.vf_range(vh2)) {
				double distance = std::abs((point - cage_mesh.point(vh2)).dot(face_normals[fh.idx()]));
				if (distance < retract_length*0.25) {
					valid_flag = false;
					break;
				}
			}
			if (!valid_flag)continue;
			smoothness_sample_points_out << point << "\n";
			smoothness_sample_points << "v " << point << "\n";
		}
	}
	for (const auto&fh : cage_mesh.faces()) {
		std::vector<Mesh::VertexHandle> facet_vertices;
		for (auto& it = cage_mesh.cfv_iter(fh); it.is_valid(); ++it) {
			facet_vertices.push_back((*it));
		}
		for (int i = 1; i < segment_number - 1; ++i) {
			for (int j = 1; j < segment_number - i; ++j) {
				int k = segment_number - i - j;
				auto& point = retract_vertices[facet_vertices[0].idx()] * i / segment_number +
					retract_vertices[facet_vertices[1].idx()] * j / segment_number +
					retract_vertices[facet_vertices[2].idx()] * k / segment_number;

				bool valid_flag = true;
				for (auto& vh : facet_vertices) {
					if (!valid_flag)break;
					for (auto& fh : cage_mesh.vf_range(vh)) {
						double distance = std::abs((point - cage_mesh.point(vh)).dot(face_normals[fh.idx()]));
						if (distance < retract_length*0.25) {
							valid_flag = false;
							break;
						}
					}
				}
				if (!valid_flag)continue;

				smoothness_sample_points_out << point << "\n";
				smoothness_sample_points << "v " << point << "\n";
			}
		}
	}
	smoothness_sample_points_out.close();
	smoothness_sample_points.close();
	std::cout << "generate smoothness sample points complete!" << "\n";
}

void PGC::read_position_constraints() {
	if (position_constraints_is_set)return;
	position_constraints_is_set = true;

	std::string line;
	double number;
	std::ifstream position_constraints_in("position_constraints.txt");
	position_constraints.clear();
	while (std::getline(position_constraints_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		position_constraints.push_back({ (int)numbers[0],{ numbers[1],numbers[2] ,numbers[3] } });
	}
	position_constraints_in.close();

	std::ifstream fixed_vertices_in("fixed_vertices.txt");
	while (fixed_vertices_in >> number) {
		fixed_vertices.push_back(number);
	}
	fixed_vertices_in.close();
	std::cout << "read position constraints complete!\n";
}

void PGC::read_orientation_constraints() {
	if (orientation_constraints_is_set)return;
	orientation_constraints_is_set = true;

	std::string line;
	double number;
	std::ifstream orientation_constraints_in("orientation_constraints.txt");
	orientation_constraints.clear();
	while (std::getline(orientation_constraints_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		orientation_constraints.push_back({ (int)numbers[0], Eigen::Map<Eigen::Matrix3d>(numbers.data() + 1) });
	}
	orientation_constraints_in.close();
	std::cout << "read orientation constraints complete!\n";
}

void PGC::read_rigidity_constraints() {
	if (rigidity_constraints_is_set)return;
	rigidity_constraints_is_set = true;
	Mesh polyline;
	OpenMesh::IO::read_mesh(polyline, "polylines.off");
	for (const auto& vh : polyline.vertices())
		rigidity_constraints.push_back(polyline.point(vh));
	std::cout << "read rigidity constraints complete!\n";
}

void PGC::read_smoothness_constraints() {
	if (smoothness_constraints_is_set)return;
	smoothness_constraints_is_set = true;
	std::ifstream smoothness_sample_points("smoothness_sample_points.txt");

	std::string line;
	double number;
	while (std::getline(smoothness_sample_points, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		smoothness_constraints.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	smoothness_sample_points.close();
	std::cout << "read smoothness constraints complete!\n";
}

void  PGC::compute_integral_para_for_sample_points(const std::string& prefix,
	const std::array<bool, 3>& compute_level, const std::vector<Mesh::Point>&sample_points) {
	std::string filename;
	if (cage_mode == false) {
		filename = "quad";
	}
	else {
		filename = "tri";
	}
	/*std::ofstream dirichlet_out(prefix + "dirichlet_para_" + filename + ".txt");
	std::ofstream neumann_out(prefix + "neumann_para_" + filename + ".txt");
	std::ofstream dirichlet_gradient_out(prefix + "dirichlet_gradient_para_" + filename + ".txt");
	std::ofstream neumann_gradient_out(prefix + "neumann_gradient_para_" + filename + ".txt");
	std::ofstream dirichlet_hessian_out(prefix + "dirichlet_hessian_para_" + filename + ".txt");
	std::ofstream neumann_hessian_out(prefix + "neumann_hessian_para_" + filename + ".txt");*/

	DirichletParametersForMeshPoints.clear();
	NeumannParametersForMeshPoints.clear();
	DDirichletParametersForMeshPoints.clear();
	DDDirichletParametersForMeshPoints.clear();
	DirichletGradientParametersForMeshPoints.clear();
	NeumannGradientParametersForMeshPoints.clear();
	DirichletHessianParametersForMeshPoints.clear();
	NeumannHessianParametersForMeshPoints.clear();


	std::clock_t start, end;
	start = std::clock();
	for (vertex_idx = 0; vertex_idx < sample_points.size(); ++vertex_idx) {
		if (vertex_idx % 1000 == 0) {
			std::cout << "curr vertex_idx:" << vertex_idx << "\n";
		}
		for (facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
			auto& facet = facet_and_control_points_pairs[facet_idx].first;

			auto& v1 = facet[1] - facet[0];
			auto& v2 = facet[facet.size() - 1] - facet[0];
			auto& v3 = facet[2] - facet[0];
			auto& point_to_p0 = facet[0] - sample_points[vertex_idx];

			/*std::cout << "\n";
			std::cout << "vertex_idx:" << vertex_idx << " facet_idx:" << facet_idx << "\n";
			std::cout << "facet:" << facet[0] << "," << facet[1] << "," << facet[2] << "\n";
			std::cout << "v1:" << v1 << "\n";
			std::cout << "v2:" << v2 << "\n";
			std::cout << "v3:" << v3 << "\n";
			std::cout << "point_to_p0:" << point_to_p0 << "\n";*/

			v1_cross_v2 = (v1%v2);
			projection_height = (point_to_p0 | v1_cross_v2);
			normal = v1_cross_v2.normalized();
			projection_to_p0 = point_to_p0 - (point_to_p0 | normal)*normal;
			projection_to_p1 = projection_to_p0 + v1;
			projection_to_p2 = projection_to_p0 + v3;
			projection_to_p3 = projection_to_p0 + v2;

			Eigen::Matrix3d M00 = convertToEigen(point_to_p0)*convertToEigen(point_to_p0).transpose();
			Eigen::Matrix3d M10 = convertToEigen(v1)*convertToEigen(point_to_p0).transpose();
			Eigen::Matrix3d M20 = convertToEigen(v2)*convertToEigen(point_to_p0).transpose();
			Eigen::Matrix3d M11 = convertToEigen(v1)*convertToEigen(v1).transpose();
			Eigen::Matrix3d M12 = convertToEigen(v1)*convertToEigen(v2).transpose();
			Eigen::Matrix3d M22 = convertToEigen(v2)*convertToEigen(v2).transpose();
			Eigen::Matrix3d M30 = convertToEigen(v1_cross_v2)*convertToEigen(point_to_p0).transpose();
			Eigen::Matrix3d M31 = convertToEigen(v1_cross_v2)*convertToEigen(v1).transpose();
			Eigen::Matrix3d M32 = convertToEigen(v1_cross_v2)*convertToEigen(v2).transpose();
			M10 += M10.transpose().eval();
			M20 += M20.transpose().eval();
			M12 += M12.transpose().eval();
			M30 += M30.transpose().eval();
			M31 += M31.transpose().eval();
			M32 += M32.transpose().eval();

			//(r,theta)坐标下p0到点(r,theta)的向量为b1*r*cos(theta)+b2*r*sin(theta)+b3
			b1 = projection_to_p0.normalized();
			b2 = (normal%b1);
			b3 = -projection_to_p0;
			//std::cout << "b1:" << b1 << "\n";
			cord_1 = ComputeCordWithBasis(v1, v2, b1);
			cord_2 = ComputeCordWithBasis(v1, v2, b2);
			cord_3 = ComputeCordWithBasis(v1, v2, b3);
			b11 = cord_1.first;
			b21 = cord_1.second;
			b12 = cord_2.first;
			b22 = cord_2.second;
			b13 = cord_3.first;
			b23 = cord_3.second;

			beta1 = ComputeAngle(projection_to_p0, projection_to_p1, normal);
			beta2 = ComputeAngle(projection_to_p1, projection_to_p2, normal);
			beta3 = ComputeAngle(projection_to_p2, projection_to_p3, normal);
			beta4 = ComputeAngle(projection_to_p3, projection_to_p0, normal);
			alpha1 = ComputeAngle(v1, -projection_to_p0, normal);
			alpha2 = ComputeAngle(v3 - v1, -projection_to_p1, normal);
			alpha3 = ComputeAngle(v2 - v3, -projection_to_p2, normal);
			alpha4 = ComputeAngle(-v2, -projection_to_p3, normal);
			c0 = (point_to_p0 - projection_to_p0).length();

			s1 = -(v1 | v1);
			s2 = -(v1 | v2);
			s3 = -(v1 | point_to_p0);
			t1 = -(v2 | v1);
			t2 = -(v2 | v2);
			t3 = -(v2 | point_to_p0);
			s = (s3 - t3 / t2 * s2) / (s1 - t1 / t2 * s2);
			t = (t3 - s3 / s1 * t1) / (t2 - s2 / s1 * t1);

			cordinate = ComputeCordWithBasis(v1, v2, v3);
			a = cordinate.first;
			b = cordinate.second;

			c_u = { (a - 1) / b, a / (b - 1), 0 };
			d_u = { 1 ,-a / (b - 1),0 };
			e_u = { {0,b},{b,1},{1,0} };

			c_v = { 0, b / (a - 1),(b - 1) / a };
			d_v = { 0,-b / (a - 1), 1 };
			e_v = { {0,1},{1,a},{a,0} };



			if (cage_mode == false) {
				for (int i = 0; i <= 2 * control_point_per_edge - 4; ++i) {
					for (int j = 0; j <= 2 * control_point_per_edge - 4; ++j) {
						NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
							NeumannSubIntegral(i, j, v1, v2, v3, point_to_p0);
						/*neumann_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";*/
					}
				}

				for (int i = 0; i <= control_point_per_edge - 1; ++i) {
					for (int j = 0; j <= control_point_per_edge - 1; ++j) {
						DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
							DirichletSubIntegral(i, j, v1, v2, v3, point_to_p0)* projection_height;
						/*dirichlet_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";*/
					}
				}
			}
			else {
				if (compute_level[1] == false && compute_level[2] == false) {
					for (int i = 0; i <= 2 * control_point_per_edge - 4; ++i) {
						for (int j = 0; j <= 2 * control_point_per_edge - 4 - i; ++j) {
							NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
								NeumannSubIntegral(i, j, v1, v2, v3, point_to_p0);
							//neumann_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							//	NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";

							//std::cout << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							//	NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
						}
					}

					for (int i = 0; i <= control_point_per_edge - 1; ++i) {
						for (int j = 0; j <= control_point_per_edge - 1 - i; ++j) {
							DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
								DirichletSubIntegral(i, j, v1, v2, v3, point_to_p0)* projection_height;
							//dirichlet_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							//	DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";

							//std::cout << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							//	DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";】
						}
					}

					////for analytic
					//for (int i = 0; i <= 2 * control_point_per_edge - 4; ++i) {
					//	for (int j = 0; j <= 2 * control_point_per_edge - 4 - i; ++j) {
					//		double analytic_result = integrate_triangle([i, j, v1, v2, point_to_p0](double x, double y) {
					//			return Neumann_integrate_func(i, j, v1, v2, point_to_p0, x, y);
					//		});
					//		NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] = analytic_result;
					//		//std::cout << "analytic result:" << analytic_result << "\n";
					//		//double diff = std::abs(NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] - analytic_result);
					//		//if (diff > 1e-3)std::cout << "N " << vertex_idx << " " << facet_idx << " " << i << " " << j << " " << " diff:" << diff << "\n";
					//	}
					//}
					////for analytic
					//for (int i = 0; i <= control_point_per_edge - 1; ++i) {
					//	for (int j = 0; j <= control_point_per_edge - 1 - i; ++j) {
					//		double analytic_result = integrate_triangle([i, j, v1, v2, point_to_p0](double x, double y) {
					//			return Dirichlet_integrate_func(i, j, v1, v2, point_to_p0, x, y);
					//		});
					//		DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] = analytic_result * projection_height;
					//		//std::cout << "analytic result:" << analytic_result << "\n";
					//		//double diff = std::abs(DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] / projection_height - analytic_result);
					//		//if (diff > 1e-3)std::cout << "D " << vertex_idx << " " << facet_idx << " " << i << " " << j << " " << " diff:" << diff << "\n";
					//	}
					//}
				}
				else {
					//for variation, use these
					for (int i = 0; i <= 2 * control_point_per_edge - 4; ++i) {
						for (int j = 0; j <= 2 * control_point_per_edge - 4 - i; ++j) {
							NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
								NeumannSubIntegral(i, j, v1, v2, v3, point_to_p0);
							//neumann_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							//	NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
						}
					}

					//for variation, use these
					for (int i = 0; i <= 2 * control_point_per_edge - 3; ++i) {
						for (int j = 0; j <= 2 * control_point_per_edge - 3 - i; ++j) {
							DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
								DirichletSubIntegral(i, j, v1, v2, v3, point_to_p0)* projection_height;
							//dirichlet_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
							//	DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
						}
					}
				}
			}

			if (compute_level[1] == false && compute_level[2] == false)continue;

			for (int i = 0; i <= 2 * control_point_per_edge - 4; ++i) {
				for (int j = 0; j <= 2 * control_point_per_edge - 4 - i; ++j) {
					NeumannGradientParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
						point_to_p0 * DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] / projection_height +
						v1 * DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j}] / projection_height +
						v2 * DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 1}] / projection_height;
					//neumann_gradient_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
					//	NeumannGradientParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
				}
			}

			for (int i = 0; i <= 2 * control_point_per_edge - 2; ++i) {
				for (int j = 0; j <= 2 * control_point_per_edge - 2 - i; ++j) {
					DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
						DDirichletSubIntegral(i, j, v1, v2, v3, point_to_p0);
					//std::cout << "i=" << i << " j=" << j << " DD=" <<
					//	DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
				}
			}

			for (int i = 0; i <= control_point_per_edge - 1; ++i) {
				for (int j = 0; j <= control_point_per_edge - 1 - i; ++j) {
					DirichletGradientParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
						-(v1%v2)*DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] / projection_height +
						3 * projection_height*(point_to_p0*DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] +
							v1 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j}] +
							v2 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 1}]);
					//dirichlet_gradient_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
					//	DirichletGradientParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
				}
			}

			if (compute_level[2] == false)continue;

			for (int i = 0; i <= control_point_per_edge + 1; ++i) {
				for (int j = 0; j <= control_point_per_edge + 1 - i; ++j) {
					DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] =
						DDDirichletSubIntegral(i, j, v1, v2, v3, point_to_p0);
					//std::cout << "i=" << i << " j=" << j << " DDD=" <<
					//	DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] << "\n";
				}
			}

			for (int i = 0; i <= 2 * control_point_per_edge - 4; ++i) {
				for (int j = 0; j <= 2 * control_point_per_edge - 4 - i; ++j) {
					Eigen::Matrix3d mat =
						-Eigen::Matrix3d::Identity()*DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] / projection_height +
						3 * (M00*DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] +
							M10 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j}] +
							M20 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 1}] +
							M12 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j + 1}] +
							M11 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 2, j}] +
							M22 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 2}]);
					NeumannHessianParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] = mat;
					//neumann_hessian_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
					//	Eigen::Map<Eigen::VectorXd>(mat.data(), mat.size()).transpose() << "\n";
				}
			}

			for (int i = 0; i <= control_point_per_edge - 1; ++i) {
				for (int j = 0; j <= control_point_per_edge - 1 - i; ++j) {
					Eigen::Matrix3d mat = -3 * ((M30 + Eigen::Matrix3d::Identity()*projection_height)*
						DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] +
						M31 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j}] +
						M32 * DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 1}]) +
						15 * projection_height*
						(M00*DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] +
							M10 * DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j}] +
							M20 * DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 1}] +
							M12 * DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 1, j + 1}] +
							M11 * DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i + 2, j}] +
							M22 * DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j + 2}]);
					DirichletHessianParametersForMeshPoints[{vertex_idx, facet_idx}][{i, j}] = mat;
					//dirichlet_hessian_out << vertex_idx << " " << facet_idx << " " << i << " " << j << " " <<
					//	Eigen::Map<Eigen::VectorXd>(mat.data(), mat.size()).transpose() << "\n";
				}
			}
		}
	}
	//dirichlet_out.close();
	//neumann_out.close();
	//dirichlet_gradient_out.close();
	//neumann_gradient_out.close();
	//dirichlet_hessian_out.close();
	//neumann_hessian_out.close();
	end = std::clock();

	std::cout << "avg per point per face time cost:" << double(end - start) / CLOCKS_PER_SEC / sample_points.size() / facet_and_control_points_pairs.size() << "s" << "\n";
	std::cout << "Complete computing " << prefix << " integral parameters!" << "\n";
	do_compute_para_for_control_points(prefix, compute_level, sample_points);
}

void PGC::do_compute_para_for_control_points(const std::string& prefix,
	const std::array<bool, 3>& compute_level, const std::vector<Mesh::Point>&sample_points) {
	std::string filename;
	if (cage_mode == false) {
		filename = "quad";
	}
	else {
		filename = "tri";
	}
	std::ofstream control_point_deformation_para_out(prefix + "control_point_deformation_para_" + filename + ".txt");
	std::ofstream control_point_pair_deformation_para_out(prefix + "control_point_pair_deformation_para_" + filename + ".txt");
	std::ofstream control_point_gradient_deformation_para_out(prefix + "control_point_gradient_deformation_para_" + filename + ".txt");
	std::ofstream control_point_pair_gradient_deformation_para_out(prefix + "control_point_pair_gradient_deformation_para_" + filename + ".txt");
	std::ofstream control_point_hessian_deformation_para_out(prefix + "control_point_hessian_deformation_para_" + filename + ".txt");
	std::ofstream control_point_pair_hessian_deformation_para_out(prefix + "control_point_pair_hessian_deformation_para_" + filename + ".txt");

	std::vector<double> map_sample_point_idx_and_control_point_to_deformation_para;
	std::vector<double> map_sample_point_idx_and_control_point_pair_to_deformation_para;
	std::vector<Mesh::Point> map_sample_point_idx_and_control_point_to_gradient_deformation_para;
	std::vector<Mesh::Point> map_sample_point_idx_and_control_point_pair_to_gradient_deformation_para;
	std::vector<Eigen::Matrix3d> map_sample_point_idx_and_control_point_to_hessian_deformation_para;
	std::vector<Eigen::Matrix3d> map_sample_point_idx_and_control_point_pair_to_hessian_deformation_para;

	if (cage_mode == true) {
		//存储方式为u,v
		std::vector<TernaryPolynomial> P;
		std::vector<TernaryPolynomial> P_u;
		std::vector<TernaryPolynomial> P_v;
		//P.resize(control_point_size_per_facet);
		P_u.resize(control_point_size_per_facet);
		P_v.resize(control_point_size_per_facet);

		TernaryPolynomial u;
		TernaryPolynomial v;
		TernaryPolynomial w;
		u.addTerm(1, 0, 0, 1);
		v.addTerm(0, 1, 0, 1);
		w.addTerm(0, 0, 0, 1);
		w.addTerm(1, 0, 0, -1);
		w.addTerm(0, 1, 0, -1);
		for (int i = 0; i < control_point_per_edge; ++i) {
			for (int j = 0; j < control_point_per_edge - i; ++j) {
				int k = control_point_per_edge - 1 - i - j;
				int factor = factorial(control_point_per_edge - 1) / factorial(i) / factorial(j) / factorial(k);
				P.push_back(u.power(j)*v.power(i)*w.power(k)*factor);
			}
		}
		//P[0] = w.power(3);
		//P[1] = u * w*w * 3;
		//P[2] = u * u*w * 3;
		//P[3] = u * u*u;
		//P[4] = v * w*w * 3;
		//P[5] = u * v*w * 6;
		//P[6] = v * u*u * 3;
		//P[7] = v * v*w * 3;
		//P[8] = v * v*u * 3;
		//P[9] = v * v*v;

		for (int i = 0; i < control_point_size_per_facet; ++i) {
			for (const auto& term : P[i].terms) {
				if (term.first[0] > 0)
					P_u[i].addTerm(term.first[0] - 1, term.first[1], term.first[2], term.second*term.first[0]);
				if (term.first[1] > 0)
					P_v[i].addTerm(term.first[0], term.first[1] - 1, term.first[2], term.second*term.first[1]);
			}
		}

		for (int vertex_idx = 0; vertex_idx < sample_points.size(); ++vertex_idx) {
			for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
				if (compute_level[0]) {
					auto& DirichletParameter = DirichletParametersForMeshPoints.at({ vertex_idx,facet_idx });
					auto& NeumannParameter = NeumannParametersForMeshPoints.at({ vertex_idx,facet_idx });

					//Dirichlet term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						double para = 0;
						for (const auto& term : P[i].terms) {
							para += term.second*DirichletParameter[{term.first[0], term.first[1]}];
						}
						para /= 4 * M_PI;
						map_sample_point_idx_and_control_point_to_deformation_para.push_back(para);
						//control_point_deformation_para_out << para << "\n";
					}

					//Neumann term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						for (int j = 0; j < control_point_size_per_facet; ++j) {
							auto&temp = P_u[i] * P_v[j];
							double para = 0;
							for (const auto& term : temp.terms) {
								para += term.second*NeumannParameter[{term.first[0], term.first[1]}];
							}
							para /= 4 * M_PI;
							map_sample_point_idx_and_control_point_pair_to_deformation_para.push_back(para);
							//control_point_pair_deformation_para_out << para << "\n";
						}
					}
				}

				if (compute_level[1]) {
					auto& DirichletGradientParameter = DirichletGradientParametersForMeshPoints.at({ vertex_idx,facet_idx });
					auto& NeumannGradientParameter = NeumannGradientParametersForMeshPoints.at({ vertex_idx,facet_idx });

					//Dirichlet gradient term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						Mesh::Point para = { 0,0,0 };
						for (const auto& term : P[i].terms) {
							para += term.second*DirichletGradientParameter[{term.first[0], term.first[1]}];
						}
						para /= 4 * M_PI;
						map_sample_point_idx_and_control_point_to_gradient_deformation_para.push_back(para);
						//control_point_gradient_deformation_para_out << para << "\n";
					}

					//Neumann gradient term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						for (int j = 0; j < control_point_size_per_facet; ++j) {
							auto&temp = P_u[i] * P_v[j];
							Mesh::Point para = { 0,0,0 };
							for (const auto& term : temp.terms) {
								para += term.second*NeumannGradientParameter[{term.first[0], term.first[1]}];
							}
							para /= 4 * M_PI;
							map_sample_point_idx_and_control_point_pair_to_gradient_deformation_para.push_back(para);
							//control_point_pair_gradient_deformation_para_out << para << "\n";
						}
					}
				}

				if (compute_level[2]) {
					auto& DirichletHessianParameter = DirichletHessianParametersForMeshPoints.at({ vertex_idx,facet_idx });
					auto& NeumannHessianParameter = NeumannHessianParametersForMeshPoints.at({ vertex_idx,facet_idx });

					//Dirichlet hessian term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						Eigen::Matrix3d para = Eigen::Matrix3d::Zero();
						for (const auto& term : P[i].terms) {
							para += term.second*DirichletHessianParameter[{term.first[0], term.first[1]}];
						}
						para /= 4 * M_PI;
						map_sample_point_idx_and_control_point_to_hessian_deformation_para.push_back(para);
						/*control_point_hessian_deformation_para_out <<
							Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";*/
					}

					//Neumann hessian term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						for (int j = 0; j < control_point_size_per_facet; ++j) {
							auto&temp = P_u[i] * P_v[j];
							Eigen::Matrix3d para = Eigen::Matrix3d::Zero();
							for (const auto& term : temp.terms) {
								para += term.second*NeumannHessianParameter[{term.first[0], term.first[1]}];
							}
							para /= 4 * M_PI;
							map_sample_point_idx_and_control_point_pair_to_hessian_deformation_para.push_back(para);
							/*control_point_pair_hessian_deformation_para_out <<
								Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";*/
						}
					}
				}
			}
		}
		for (auto& para : map_sample_point_idx_and_control_point_to_deformation_para)
			control_point_deformation_para_out << std::setprecision(20) << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_pair_to_deformation_para)
			control_point_pair_deformation_para_out << std::setprecision(20) << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_to_gradient_deformation_para)
			control_point_gradient_deformation_para_out << std::setprecision(20) << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_pair_to_gradient_deformation_para)
			control_point_pair_gradient_deformation_para_out << std::setprecision(20) << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_to_hessian_deformation_para)
			control_point_hessian_deformation_para_out << std::setprecision(20) << Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_pair_to_hessian_deformation_para)
			control_point_pair_hessian_deformation_para_out << std::setprecision(20) << Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";

		control_point_deformation_para_out.close();
		control_point_pair_deformation_para_out.close();
		control_point_gradient_deformation_para_out.close();
		control_point_pair_gradient_deformation_para_out.close();
		control_point_hessian_deformation_para_out.close();
		control_point_pair_hessian_deformation_para_out.close();
	}
	else if (cage_mode == false) {
		//存储方式为u,v
		std::vector<TernaryPolynomial> P;
		std::vector<TernaryPolynomial> P_u;
		std::vector<TernaryPolynomial> P_v;
		//P.resize(control_point_size_per_facet);
		P_u.resize(control_point_size_per_facet);
		P_v.resize(control_point_size_per_facet);

		TernaryPolynomial u;
		TernaryPolynomial v;
		TernaryPolynomial one_minus_u;
		TernaryPolynomial one_minus_v;
		u.addTerm(1, 0, 0, 1);
		v.addTerm(0, 1, 0, 1);
		one_minus_u.addTerm(0, 0, 0, 1);
		one_minus_u.addTerm(1, 0, 0, -1);
		one_minus_v.addTerm(0, 0, 0, 1);
		one_minus_v.addTerm(0, 1, 0, -1);

		for (int i = 0; i < control_point_per_edge; ++i) {
			for (int j = 0; j < control_point_per_edge; ++j) {
				int factor = factorial(control_point_per_edge - 1) / factorial(i) /
					factorial(control_point_per_edge - 1 - i)*factorial(control_point_per_edge - 1) /
					factorial(j) / factorial(control_point_per_edge - 1 - j);
				P.push_back(u.power(j)*one_minus_u.power(control_point_per_edge - 1 - j)*v.power(i)*
					one_minus_v.power(control_point_per_edge - 1 - i)*factor);
			}
		}

		//P[0].addTerm(0, 0, 0, 1);
		//P[0].addTerm(1, 0, 0, -1);
		//P[0].addTerm(0, 1, 0, -1);
		//P[0].addTerm(1, 1, 0, 1);
		//P[1].addTerm(1, 0, 0, 1);
		//P[1].addTerm(1, 1, 0, -1);
		//P[2].addTerm(1, 1, 0, 1);
		//P[3].addTerm(0, 1, 0, 1);
		//P[3].addTerm(1, 1, 0, -1);

		for (int i = 0; i < control_point_size_per_facet; ++i) {
			for (const auto& term : P[i].terms) {
				if (term.first[0] > 0)
					P_u[i].addTerm(term.first[0] - 1, term.first[1], term.first[2], term.second*term.first[0]);
				if (term.first[1] > 0)
					P_v[i].addTerm(term.first[0], term.first[1] - 1, term.first[2], term.second*term.first[1]);
			}
		}

		for (int vertex_idx = 0; vertex_idx < sample_points.size(); ++vertex_idx) {
			for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
				if (compute_level[0]) {
					auto& DirichletParameter = DirichletParametersForMeshPoints.at({ vertex_idx,facet_idx });
					auto& NeumannParameter = NeumannParametersForMeshPoints.at({ vertex_idx,facet_idx });

					//Dirichlet term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						double para = 0;
						for (const auto& term : P[i].terms) {
							para += term.second*DirichletParameter[{term.first[0], term.first[1]}];
						}
						para /= 4 * M_PI;
						map_sample_point_idx_and_control_point_to_deformation_para.push_back(para);
						//control_point_deformation_para_out << para << "\n";
					}

					//Neumann term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						for (int j = 0; j < control_point_size_per_facet; ++j) {
							auto&temp = P_u[i] * P_v[j];
							double para = 0;
							for (const auto& term : temp.terms) {
								para += term.second*NeumannParameter[{term.first[0], term.first[1]}];
							}
							para /= 4 * M_PI;
							map_sample_point_idx_and_control_point_pair_to_deformation_para.push_back(para);
							//control_point_pair_deformation_para_out << para << "\n";
						}
					}
				}

				if (compute_level[1]) {
					auto& DirichletGradientParameter = DirichletGradientParametersForMeshPoints.at({ vertex_idx,facet_idx });
					auto& NeumannGradientParameter = NeumannGradientParametersForMeshPoints.at({ vertex_idx,facet_idx });

					//Dirichlet gradient term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						Mesh::Point para = { 0,0,0 };
						for (const auto& term : P[i].terms) {
							para += term.second*DirichletGradientParameter[{term.first[0], term.first[1]}];
						}
						para /= 4 * M_PI;
						map_sample_point_idx_and_control_point_to_gradient_deformation_para.push_back(para);
						//control_point_gradient_deformation_para_out << para << "\n";
					}

					//Neumann gradient term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						for (int j = 0; j < control_point_size_per_facet; ++j) {
							auto&temp = P_u[i] * P_v[j];
							Mesh::Point para = { 0,0,0 };
							for (const auto& term : temp.terms) {
								para += term.second*NeumannGradientParameter[{term.first[0], term.first[1]}];
							}
							para /= 4 * M_PI;
							map_sample_point_idx_and_control_point_pair_to_gradient_deformation_para.push_back(para);
							//control_point_pair_gradient_deformation_para_out << para << "\n";
						}
					}
				}

				if (compute_level[2]) {
					auto& DirichletHessianParameter = DirichletHessianParametersForMeshPoints.at({ vertex_idx,facet_idx });
					auto& NeumannHessianParameter = NeumannHessianParametersForMeshPoints.at({ vertex_idx,facet_idx });

					//Dirichlet hessian term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						Eigen::Matrix3d para = Eigen::Matrix3d::Zero();
						for (const auto& term : P[i].terms) {
							para += term.second*DirichletHessianParameter[{term.first[0], term.first[1]}];
						}
						para /= 4 * M_PI;
						map_sample_point_idx_and_control_point_to_hessian_deformation_para.push_back(para);
						//control_point_hessian_deformation_para_out <<
						//	Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";
					}

					//Neumann hessian term
					for (int i = 0; i < control_point_size_per_facet; ++i) {
						for (int j = 0; j < control_point_size_per_facet; ++j) {
							auto&temp = P_u[i] * P_v[j];
							Eigen::Matrix3d para = Eigen::Matrix3d::Zero();
							for (const auto& term : temp.terms) {
								para += term.second*NeumannHessianParameter[{term.first[0], term.first[1]}];
							}
							para /= 4 * M_PI;
							map_sample_point_idx_and_control_point_pair_to_hessian_deformation_para.push_back(para);
							//control_point_pair_hessian_deformation_para_out <<
							//	Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";
						}
					}
				}
			}
		}
		for (auto& para : map_sample_point_idx_and_control_point_to_deformation_para)
			control_point_deformation_para_out << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_pair_to_deformation_para)
			control_point_pair_deformation_para_out << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_to_gradient_deformation_para)
			control_point_gradient_deformation_para_out << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_pair_to_gradient_deformation_para)
			control_point_pair_gradient_deformation_para_out << para << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_to_hessian_deformation_para)
			control_point_hessian_deformation_para_out << Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";
		for (auto& para : map_sample_point_idx_and_control_point_pair_to_hessian_deformation_para)
			control_point_pair_hessian_deformation_para_out << Eigen::Map<Eigen::VectorXd>(para.data(), para.size()).transpose() << "\n";

		control_point_deformation_para_out.close();
		control_point_pair_deformation_para_out.close();
		control_point_gradient_deformation_para_out.close();
		control_point_pair_gradient_deformation_para_out.close();
		control_point_hessian_deformation_para_out.close();
		control_point_pair_hessian_deformation_para_out.close();
	}
	std::cout << "Compute " << prefix << " para for control points completed!" << "\n";
}

void PGC::do_optimize_control_points() {
	if (!constraint_is_set) {
		read_position_constraints();
		read_orientation_constraints();
		read_rigidity_constraints();

		std::string line;
		std::ifstream smoothness_sample_points("smoothness_sample_points.txt");

		while (std::getline(smoothness_sample_points, line)) {
			std::istringstream iss(line);
			std::vector<double> numbers;

			double number;
			while (iss >> number) {
				numbers.push_back(number);
			}
			smoothness_constraints.push_back({ numbers[0],numbers[1] ,numbers[2] });
		}
		smoothness_sample_points.close();
	}

	std::vector<double> map_vertex_idx_and_control_point_to_deformation_para;
	std::vector<double> map_vertex_idx_and_control_point_pair_to_deformation_para;

	std::vector<double> map_vertex_idx_and_control_point_to_rigidity_deformation_para;
	std::vector<double> map_vertex_idx_and_control_point_pair_to_rigidity_deformation_para;

	std::vector<double> map_vertex_idx_and_control_point_to_smoothness_deformation_para;
	std::vector<double> map_vertex_idx_and_control_point_pair_to_smoothness_deformation_para;

	std::vector<Mesh::Point> map_vertex_idx_and_control_point_to_rigidity_gradient_deformation_para;
	std::vector<Mesh::Point> map_vertex_idx_and_control_point_pair_to_rigidity_gradient_deformation_para;

	std::vector<Eigen::Matrix3d> map_vertex_idx_and_control_point_to_hessian_deformation_para;
	std::vector<Eigen::Matrix3d> map_vertex_idx_and_control_point_pair_to_hessian_deformation_para;

	std::vector<Eigen::Matrix3d> map_vertex_idx_and_control_point_to_smoothness_hessian_deformation_para;
	std::vector<Eigen::Matrix3d> map_vertex_idx_and_control_point_pair_to_smoothness_hessian_deformation_para;

	std::string filename;
	if (cage_mode == false) {
		filename = "quad";
	}
	else {
		filename = "tri";
	}

	std::ifstream control_point_deformation_para_in("control_point_deformation_para_" + filename + ".txt");
	std::ifstream control_point_pair_deformation_para_in("control_point_pair_deformation_para_" + filename + ".txt");

	std::ifstream rigidity_control_point_deformation_para_in("rigidity_control_point_deformation_para_" + filename + ".txt");
	std::ifstream rigidity_control_point_pair_deformation_para_in("rigidity_control_point_pair_deformation_para_" + filename + ".txt");

	std::ifstream smoothness_control_point_deformation_para_in("smoothness_control_point_deformation_para_" + filename + ".txt");
	std::ifstream smoothness_control_point_pair_deformation_para_in("smoothness_control_point_pair_deformation_para_" + filename + ".txt");

	std::ifstream control_point_gradient_deformation_para_in("control_point_gradient_deformation_para_" + filename + ".txt");
	std::ifstream control_point_pair_gradient_deformation_para_in("control_point_pair_gradient_deformation_para_" + filename + ".txt");

	std::ifstream rigidity_control_point_gradient_deformation_para_in("rigidity_control_point_gradient_deformation_para_" + filename + ".txt");
	std::ifstream rigidity_control_point_pair_gradient_deformation_para_in("rigidity_control_point_pair_gradient_deformation_para_" + filename + ".txt");

	std::ifstream control_point_hessian_deformation_para_in("control_point_hessian_deformation_para_" + filename + ".txt");
	std::ifstream control_point_pair_hessian_deformation_para_in("control_point_pair_hessian_deformation_para_" + filename + ".txt");

	std::ifstream smoothness_control_point_hessian_deformation_para_in("smoothness_control_point_hessian_deformation_para_" + filename + ".txt");
	std::ifstream smoothness_control_point_pair_hessian_deformation_para_in("smoothness_control_point_pair_hessian_deformation_para_" + filename + ".txt");

	double number;
	std::string line;
	while (control_point_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_to_deformation_para.push_back(number);
	}
	control_point_deformation_para_in.close();

	while (control_point_pair_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_pair_to_deformation_para.push_back(number);
	}
	control_point_pair_deformation_para_in.close();

	while (rigidity_control_point_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_to_rigidity_deformation_para.push_back(number);
	}
	rigidity_control_point_deformation_para_in.close();

	while (rigidity_control_point_pair_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_pair_to_rigidity_deformation_para.push_back(number);
	}
	rigidity_control_point_pair_deformation_para_in.close();

	while (smoothness_control_point_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_to_smoothness_deformation_para.push_back(number);
	}
	smoothness_control_point_deformation_para_in.close();

	while (smoothness_control_point_pair_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_pair_to_smoothness_deformation_para.push_back(number);
	}
	smoothness_control_point_pair_deformation_para_in.close();

	//gradient

	while (std::getline(control_point_gradient_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;
		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_to_gradient_deformation_para.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	control_point_gradient_deformation_para_in.close();

	while (std::getline(control_point_pair_gradient_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;
		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_pair_to_gradient_deformation_para.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	control_point_pair_gradient_deformation_para_in.close();

	while (std::getline(rigidity_control_point_gradient_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;
		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_to_rigidity_gradient_deformation_para.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	rigidity_control_point_gradient_deformation_para_in.close();

	while (std::getline(rigidity_control_point_pair_gradient_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;
		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_pair_to_rigidity_gradient_deformation_para.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	rigidity_control_point_pair_gradient_deformation_para_in.close();

	//hessian
	while (std::getline(control_point_hessian_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_to_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
	}
	control_point_hessian_deformation_para_in.close();

	while (std::getline(control_point_pair_hessian_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_pair_to_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
	}
	control_point_pair_hessian_deformation_para_in.close();

	while (std::getline(smoothness_control_point_hessian_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_to_smoothness_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
	}
	smoothness_control_point_hessian_deformation_para_in.close();

	while (std::getline(smoothness_control_point_pair_hessian_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_pair_to_smoothness_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
	}
	smoothness_control_point_pair_hessian_deformation_para_in.close();

	int control_point_pairs_per_facet = control_point_size_per_facet * control_point_size_per_facet;
	std::vector<int> map_var_to_represent_var(control_points.size() +
		facet_and_control_points_pairs.size() * control_point_pairs_per_facet, -1);
	std::vector<int> map_var_to_represent_var_sign(control_points.size() +
		facet_and_control_points_pairs.size() *control_point_pairs_per_facet, 1);
	int var_size = represent_control_points.size();
	for (int i = 0; i < control_points.size(); ++i) {
		map_var_to_represent_var[i] = map_control_point_to_represent_control_point_idx[i];
	}
	int not_appeared_flag = map_var_to_represent_var.size();

	std::vector<int> represent_control_pair_is_appeared(represent_control_points.size()*represent_control_points.size(), not_appeared_flag);
	for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
		for (int control_point_idx_0 = 0; control_point_idx_0 < control_point_size_per_facet; ++control_point_idx_0) {
			for (int control_point_idx_1 = 0; control_point_idx_1 < control_point_size_per_facet; ++control_point_idx_1) {
				int index = facet_idx * control_point_pairs_per_facet + control_point_idx_0 * control_point_size_per_facet + control_point_idx_1;
				int represent_control_point_idx_0 = map_control_point_to_represent_control_point_idx[control_point_idx_0 + facet_idx * control_point_size_per_facet];
				int represent_control_point_idx_1 = map_control_point_to_represent_control_point_idx[control_point_idx_1 + facet_idx * control_point_size_per_facet];

				if (represent_control_point_idx_0 > represent_control_point_idx_1) {
					int i = represent_control_point_idx_0 * represent_control_points.size() + represent_control_point_idx_1;
					if (represent_control_pair_is_appeared[i] == not_appeared_flag) {
						map_var_to_represent_var[control_points.size() + index] = var_size;
						represent_control_pair_is_appeared[i] = var_size;
						var_size++;
					}
					else {
						map_var_to_represent_var[control_points.size() + index] = represent_control_pair_is_appeared[i];
					}
				}
				else if (represent_control_point_idx_0 < represent_control_point_idx_1) {
					int i = represent_control_point_idx_1 * represent_control_points.size() + represent_control_point_idx_0;
					if (represent_control_pair_is_appeared[i] == not_appeared_flag) {
						map_var_to_represent_var[control_points.size() + index] = var_size;
						map_var_to_represent_var_sign[control_points.size() + index] = -1;
						represent_control_pair_is_appeared[i] = var_size;
						var_size++;
					}
					else {
						map_var_to_represent_var[control_points.size() + index] = represent_control_pair_is_appeared[i];
						map_var_to_represent_var_sign[control_points.size() + index] = -1;
					}
				}
				else {
					map_var_to_represent_var[control_points.size() + index] = 0;
					map_var_to_represent_var_sign[control_points.size() + index] = 0;
				}
			}
		}
	}
	var_size *= 3;

#if DEBUG_MODE
		double total_hessian_norm = 0;
		for (int i = 0; i < smoothness_constraints.size(); ++i)
		{
			Tensor hessian = {};

			for (int control_point_idx = 0; control_point_idx < control_points.size(); ++control_point_idx) {
				int index = i * control_points.size() + control_point_idx;

				sum_point_product(hessian, control_points[control_point_idx], map_vertex_idx_and_control_point_to_smoothness_hessian_deformation_para[index]);
			}

			for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
				for (int control_point_idx_0 = 0; control_point_idx_0 < control_point_size_per_facet; ++control_point_idx_0) {
					for (int control_point_idx_1 = 0; control_point_idx_1 < control_point_size_per_facet; ++control_point_idx_1) {
						int index = i * facet_and_control_points_pairs.size() * control_point_pairs_per_facet + facet_idx * control_point_pairs_per_facet +
							control_point_idx_0 * control_point_size_per_facet + control_point_idx_1;
						auto& P0_cross_P1 = (control_points[control_point_idx_0 +
							facet_idx * control_point_size_per_facet] % control_points[control_point_idx_1 + facet_idx * control_point_size_per_facet]);

						sum_point_product(hessian, P0_cross_P1, map_vertex_idx_and_control_point_pair_to_smoothness_hessian_deformation_para[index]);
					}
				}
			}
			//std::cout << "before:" << smoothness_constraints[i] << " after:" << deformed_point << "\n";
			double hessian_norm = 0;
			for (const auto& matrix : hessian) {
				for (const auto& row : matrix) {
					for (const auto& element : row) {
						hessian_norm += element * element;
					}
				}
			}
			total_hessian_norm += hessian_norm;
			//std::cout << i << " hessian norm:" << hessian_norm << "\n";
		}
		std::cout << "smoothness sample points total hessian norm:" << total_hessian_norm << " avg hessian norm:"
			<< total_hessian_norm / smoothness_constraints.size() << "\n";

		total_hessian_norm = 0;
		for (int i = 0; i < mesh.n_vertices(); ++i)
		{
			Tensor hessian = {};

			for (int control_point_idx = 0; control_point_idx < control_points.size(); ++control_point_idx) {
				int index = i * control_points.size() + control_point_idx;

				sum_point_product(hessian, control_points[control_point_idx], map_vertex_idx_and_control_point_to_hessian_deformation_para[index]);
			}

			for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
				for (int control_point_idx_0 = 0; control_point_idx_0 < control_point_size_per_facet; ++control_point_idx_0) {
					for (int control_point_idx_1 = 0; control_point_idx_1 < control_point_size_per_facet; ++control_point_idx_1) {
						int index = i * facet_and_control_points_pairs.size() * control_point_pairs_per_facet +
							facet_idx * control_point_pairs_per_facet + control_point_idx_0 * control_point_size_per_facet + control_point_idx_1;
						auto& P0_cross_P1 = (control_points[control_point_idx_0 +
							facet_idx * control_point_size_per_facet] % control_points[control_point_idx_1 +
							facet_idx * control_point_size_per_facet]);

						sum_point_product(hessian, P0_cross_P1, map_vertex_idx_and_control_point_pair_to_hessian_deformation_para[index]);
					}
				}
			}
			//std::cout << "before:" << smoothness_constraints[i] << " after:" << deformed_point << "\n";
			double hessian_norm = 0;
			for (const auto& matrix : hessian) {
				for (const auto& row : matrix) {
					for (const auto& element : row) {
						hessian_norm += element * element;
					}
				}
			}
			total_hessian_norm += hessian_norm;
			//std::cout << i << " hessian norm:" << hessian_norm << "\n";
		}
		std::cout << "mesh total hessian norm:" << total_hessian_norm << " avg hessian norm:"
			<< total_hessian_norm / mesh.n_vertices() << "\n";

		//Eigen::VectorXd x(var_size);
		//std::ifstream optimized_control_points_in("optimized_control_points.obj");
		//int i = 0;
		//while (optimized_control_points_in >> number) {
		//	x[i] = number;
		//	i++;
		//}
		//optimized_control_points_in.close();

		//std::ifstream optimized_control_point_pairs_in("optimized_control_point_pairs.obj");
		//while (optimized_control_point_pairs_in >> number) {
		//	x[i] = number;
		//	i++;
		//}
		//optimized_control_point_pairs_in.close();s
#endif

	std::cout << "rigidity_constraints.size():" << rigidity_constraints.size() <<
		" smoothness_constraints.size():" << smoothness_constraints.size() << "\n";
	int equation_size = 3 * position_constraints.size() + 9 * orientation_constraints.size() +
		9 * rigidity_constraints.size() + 15 * smoothness_constraints.size() + 3 * fixed_vertices.size();
	std::cout << "equation_num:" << equation_size << ",var_size:" << var_size << "\n";

	Eigen::MatrixXd global_solver_mat(equation_size, var_size);
	Eigen::VectorXd global_solver_vec(equation_size);
	global_solver_mat.setZero();
	global_solver_vec.setZero();

	int mat_line_begin = 0;
	for (int i = 0; i < position_constraints.size(); ++i) {
		int vertex_idx = position_constraints[i].first;
		for (int j = 0; j < control_points.size(); ++j) {
			for (int k = 0; k < 3; ++k) {
				global_solver_mat(3 * i + k, 3 * map_var_to_represent_var[j] + k) +=
					position_weight * map_vertex_idx_and_control_point_to_deformation_para[vertex_idx*control_points.size() + j];
			}
		}
		for (int j = 0; j < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++j) {
			for (int k = 0; k < 3; ++k) {
				global_solver_mat(3 * i + k, 3 * map_var_to_represent_var[control_points.size() + j] + k) += position_weight *
					map_var_to_represent_var_sign[control_points.size() + j] *
					map_vertex_idx_and_control_point_pair_to_deformation_para
					[facet_and_control_points_pairs.size() * control_point_pairs_per_facet * vertex_idx + j];
			}
		}
		global_solver_vec(3 * i) += position_weight * position_constraints[i].second[0];
		global_solver_vec(3 * i + 1) += position_weight * position_constraints[i].second[1];
		global_solver_vec(3 * i + 2) += position_weight * position_constraints[i].second[2];
	}

	mat_line_begin += 3 * position_constraints.size();

	for (int i = 0; i < orientation_constraints.size(); ++i) {
		int vertex_idx = orientation_constraints[i].first;
		for (int j = 0; j < control_points.size(); ++j) {
			auto& temp = orientation_weight * map_vertex_idx_and_control_point_to_gradient_deformation_para[vertex_idx*control_points.size() + j];
			for (int ii = 0; ii < 3; ++ii) {
				for (int jj = 0; jj < 3; ++jj) {
					global_solver_mat(mat_line_begin + 9 * i + 3 * ii + jj, 3 * map_var_to_represent_var[j] + ii) += temp[jj];
				}
			}
		}
		for (int j = 0; j < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++j) {
			auto& temp = orientation_weight * map_vertex_idx_and_control_point_pair_to_gradient_deformation_para
				[facet_and_control_points_pairs.size() * control_point_pairs_per_facet * vertex_idx + j];
			for (int ii = 0; ii < 3; ++ii) {
				for (int jj = 0; jj < 3; ++jj) {
					global_solver_mat(mat_line_begin + 9 * i + 3 * ii + jj, 3 * map_var_to_represent_var[control_points.size() + j] + ii) +=
						map_var_to_represent_var_sign[control_points.size() + j] * temp[jj];
				}
			}
		}
		auto& matrix = orientation_weight * orientation_constraints[i].second;
		global_solver_vec(mat_line_begin + 9 * i) += matrix(0,0);
		global_solver_vec(mat_line_begin + 9 * i + 1) += matrix(0, 1);
		global_solver_vec(mat_line_begin + 9 * i + 2) += matrix(0, 2);
		global_solver_vec(mat_line_begin + 9 * i + 3) += matrix(1, 0);
		global_solver_vec(mat_line_begin + 9 * i + 4) += matrix(1, 1);
		global_solver_vec(mat_line_begin + 9 * i + 5) += matrix(1, 2);
		global_solver_vec(mat_line_begin + 9 * i + 6) += matrix(2, 0);
		global_solver_vec(mat_line_begin + 9 * i + 7) += matrix(2, 1);
		global_solver_vec(mat_line_begin + 9 * i + 8) += matrix(2, 2);
	}

	mat_line_begin += 9 * orientation_constraints.size();
	for (int i = 0; i < rigidity_constraints.size(); ++i) {
		for (int j = 0; j < control_points.size(); ++j) {
			auto& temp = rigidity_weight * map_vertex_idx_and_control_point_to_rigidity_gradient_deformation_para[i*control_points.size() + j];
			for (int ii = 0; ii < 3; ++ii) {
				for (int jj = 0; jj < 3; ++jj) {
					global_solver_mat(mat_line_begin + 9 * i + 3 * ii + jj, 3 * map_var_to_represent_var[j] + ii) += temp[jj];
				}
			}
		}
		for (int j = 0; j < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++j) {
			auto& temp = rigidity_weight * map_vertex_idx_and_control_point_pair_to_rigidity_gradient_deformation_para
				[facet_and_control_points_pairs.size() * control_point_pairs_per_facet * i + j];
			for (int ii = 0; ii < 3; ++ii) {
				for (int jj = 0; jj < 3; ++jj) {
					global_solver_mat(mat_line_begin + 9 * i + 3 * ii + jj, 3 * map_var_to_represent_var[control_points.size() + j] + ii) +=
						map_var_to_represent_var_sign[control_points.size() + j] * temp[jj];
				}
			}
		}
		//global_solver_vec(mat_line_begin + 9 * i) += rigidity_weight;
		//global_solver_vec(mat_line_begin + 9 * i + 1) += 0;
		//global_solver_vec(mat_line_begin + 9 * i + 2) += 0;
		//global_solver_vec(mat_line_begin + 9 * i + 3) += 0;
		//global_solver_vec(mat_line_begin + 9 * i + 4) += rigidity_weight;
		//global_solver_vec(mat_line_begin + 9 * i + 5) += 0;
		//global_solver_vec(mat_line_begin + 9 * i + 6) += 0;
		//global_solver_vec(mat_line_begin + 9 * i + 7) += 0;
		//global_solver_vec(mat_line_begin + 9 * i + 8) += rigidity_weight;
	}

	mat_line_begin += 9 * rigidity_constraints.size();
	for (int i = 0; i < smoothness_constraints.size(); ++i) {
		for (int j = 0; j < control_points.size(); ++j) {
			auto& temp = smoothness_weight * map_vertex_idx_and_control_point_to_smoothness_hessian_deformation_para
				[i * control_points.size() + j];
			for (int k = 0; k < 3; ++k) {
				for (int ii = 0; ii < 3; ++ii) {
					for (int jj = 0; jj <= ii && jj != 2; ++jj) {
						global_solver_mat(mat_line_begin + 15 * i + 5 * k + (ii * ii + ii) / 2 + jj, 3 * map_var_to_represent_var[j] + k) += temp(ii, jj);
					}
				}
			}
		}
		for (int j = 0; j < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++j) {
			auto& temp = smoothness_weight * map_vertex_idx_and_control_point_pair_to_smoothness_hessian_deformation_para
				[i * facet_and_control_points_pairs.size() * control_point_pairs_per_facet + j];
			for (int k = 0; k < 3; ++k) {
				for (int ii = 0; ii < 3; ++ii) {
					for (int jj = 0; jj <= ii && jj != 2; ++jj) {
						global_solver_mat(mat_line_begin + 15 * i + 5 * k + (ii * ii + ii) / 2 + jj,
							3 * map_var_to_represent_var[control_points.size() + j] + k) +=
							map_var_to_represent_var_sign[control_points.size() + j] * temp(ii, jj);
					}
				}
			}
		}
		//for (int j = 0; j < 15; ++j) {
		//	global_solver_vec(mat_line_begin + 15 * i + j) += 0;
		//}
	}

	mat_line_begin += 15 * smoothness_constraints.size();
	for (int i = 0; i < fixed_vertices.size(); ++i) {
		for (int j = 0; j < control_points.size(); ++j) {
			for (int k = 0; k < 3; ++k) {
				global_solver_mat(mat_line_begin + 3 * i + k, 3 * map_var_to_represent_var[j] + k) += fixed_vertices_weight *
					map_vertex_idx_and_control_point_to_deformation_para[fixed_vertices[i] * control_points.size() + j];
			}
		}
		for (int j = 0; j < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++j) {
			for (int k = 0; k < 3; ++k) {
				global_solver_mat(mat_line_begin + 3 * i + k, 3 * map_var_to_represent_var[control_points.size() + j] + k) +=
					fixed_vertices_weight * map_var_to_represent_var_sign[control_points.size() + j] *
					map_vertex_idx_and_control_point_pair_to_deformation_para
					[facet_and_control_points_pairs.size() * control_point_pairs_per_facet * fixed_vertices[i] + j];
			}
		}
		auto& point = fixed_vertices_weight * mesh.point(mesh.vertex_handle(fixed_vertices[i]));
		global_solver_vec(mat_line_begin + 3 * i) += point[0];
		global_solver_vec(mat_line_begin + 3 * i + 1) += point[1];
		global_solver_vec(mat_line_begin + 3 * i + 2) += point[2];
	}

	std::ofstream ofs("global_solver_mat.dat", std::ios::binary);
	int rows = global_solver_mat.rows();
	int cols = global_solver_mat.cols();
	ofs.write(reinterpret_cast<const char*>(&rows), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&cols), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(global_solver_mat.data()), rows * cols * sizeof(double));
	ofs.close();

	std::ofstream global_solver_vec_out("global_solver_vec.txt");
	for (int i = 0; i < equation_size; ++i) {
		global_solver_vec_out << global_solver_vec(i) << "\n";
	}
	global_solver_vec_out.close();

	// SVD
	std::cout << "begin SVD" << "\n";
	std::clock_t start, end;
	start = std::clock();
	Eigen::BDCSVD<Eigen::MatrixXd> solver(global_solver_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
	solver.setThreshold(1e-10);
	saveBDCSVD(solver, "SVD.dat");
	end = std::clock();

	std::cout << "decomposition complete, spend " << double(end - start) / CLOCKS_PER_SEC << "s" << "\n";
	std::cout << "threshold:" << solver.threshold() << " tolerance:"
		<< solver.singularValues().coeff(0)*solver.threshold() << " rank:" << solver.rank() << "\n";

	loadBDCSVD();

	int rigidity_start = 3 * position_constraints.size() + 9 * orientation_constraints.size();
	double energy = DBL_MAX;
	double energy_diff = DBL_MAX;
	int iter = 0;
	Eigen::VectorXd global_solver_rotation_vec(equation_size);
	global_solver_rotation_vec.setZero();
	for (int i = 0; i < rigidity_constraints.size(); ++i) {
		Eigen::Vector3d axis = Eigen::Vector3d::Random().normalized();
		double angle = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
		Eigen::AngleAxisd angleAxis(angle, axis);
		auto& rotation_mat = angleAxis.toRotationMatrix();
		global_solver_rotation_vec[rigidity_start + 9 * i] = rotation_mat(0, 0);
		global_solver_rotation_vec[rigidity_start + 9 * i + 1] = rotation_mat(1, 0);
		global_solver_rotation_vec[rigidity_start + 9 * i + 2] = rotation_mat(2, 0);
		global_solver_rotation_vec[rigidity_start + 9 * i + 3] = rotation_mat(0, 1);
		global_solver_rotation_vec[rigidity_start + 9 * i + 4] = rotation_mat(1, 1);
		global_solver_rotation_vec[rigidity_start + 9 * i + 5] = rotation_mat(2, 1);
		global_solver_rotation_vec[rigidity_start + 9 * i + 6] = rotation_mat(0, 2);
		global_solver_rotation_vec[rigidity_start + 9 * i + 7] = rotation_mat(1, 2);
		global_solver_rotation_vec[rigidity_start + 9 * i + 8] = rotation_mat(2, 2);
	}
	Eigen::VectorXd x;

	int section_1 = 3 * position_constraints.size();
	int section_2 = section_1 + 9 * orientation_constraints.size();
	int section_3 = section_2 + 9 * rigidity_constraints.size();
	int section_4 = section_3 + 15 * smoothness_constraints.size();

	// Local/Global
	start = std::clock();
	while ((energy_diff > 1e-7 || iter < 5) && iter < 1000) {
		std::cout << "Iteration " << iter + 1 << "\n";
		//Global
		x = solver.solve(global_solver_vec + global_solver_rotation_vec);
#if DEBUG_MODE
			auto& x_bar = SVD_mat * (global_solver_vec + global_solver_rotation_vec);
			std::cout << "diff between x and x_bar:" << (x - x_bar).squaredNorm() << "\n";
#endif
		Eigen::VectorXd x_star = global_solver_mat * x;

		//Local
		for (int i = 0; i < rigidity_constraints.size(); ++i) {
			Eigen::Matrix3d gradient = Eigen::Map<Eigen::Matrix3d>(x_star.data() + rigidity_start + 9 * i) / rigidity_weight;
			auto& rotation_mat = optimalRotationMatrix(gradient)*rigidity_weight;
			global_solver_rotation_vec[rigidity_start + 9 * i] = rotation_mat(0, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 1] = rotation_mat(1, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 2] = rotation_mat(2, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 3] = rotation_mat(0, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 4] = rotation_mat(1, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 5] = rotation_mat(2, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 6] = rotation_mat(0, 2);
			global_solver_rotation_vec[rigidity_start + 9 * i + 7] = rotation_mat(1, 2);
			global_solver_rotation_vec[rigidity_start + 9 * i + 8] = rotation_mat(2, 2);
		}

		auto& delta = x_star - global_solver_vec - global_solver_rotation_vec;
		double E_total = delta.squaredNorm();
		double E_position = 0, E_orientation = 0, E_rigidity = 0, E_smoothness = 0, E_fixed_vertices = 0;
		for (int i = 0; i < section_1; ++i) {
			E_position += delta(i)*delta(i);
		}
		for (int i = section_1; i < section_2; ++i) {
			E_orientation += delta(i)*delta(i);
		}
		for (int i = section_2; i < section_3; ++i) {
			E_rigidity += delta(i)*delta(i);
		}
		for (int i = section_3; i < section_4; ++i) {
			E_smoothness += delta(i)*delta(i);
		}
		for (int i = section_4; i < delta.size(); ++i) {
			E_fixed_vertices += delta(i)*delta(i);
		}

		energy_diff = std::abs(energy - E_total);
		energy = E_total;
		std::cout << "total energy=" << E_total << " avg E_position=" << E_position / position_constraints.size()
			<< "avg E_orientation=" << E_orientation / orientation_constraints.size()
			<< " avg E_rigidity=" << E_rigidity / rigidity_constraints.size()
			<< " avg E_smoothness=" << E_smoothness / smoothness_constraints.size()
			<< " avg E_fixed_vertices=" << E_fixed_vertices / fixed_vertices.size()
			<< " energy_diff=" << energy_diff << "\n";

		iter++;
	}
	end = std::clock();
	std::cout << "Complete Optimizing, spend " << double(end - start) / CLOCKS_PER_SEC * 1000 << "ms" << "\n";
	std::cout << "average time cost per step:" << double(end - start) / CLOCKS_PER_SEC * 1000 / (iter + 1) << "ms" << "\n";

	/*Eigen::VectorXd origin_x(var_size);
	int index = 0;
	for (int i = 0; i < control_points.size(); ++i) {
		int j = map_var_to_represent_var[i];
		origin_x[3 * j] = control_points[i][0];
		origin_x[3 * j + 1] = control_points[i][1];
		origin_x[3 * j + 2] = control_points[i][2];
	}
	for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
		for (int control_point_idx_0 = 0; control_point_idx_0 < control_point_size_per_facet; ++control_point_idx_0) {
			for (int control_point_idx_1 = 0; control_point_idx_1 < control_point_size_per_facet; ++control_point_idx_1) {
				auto& P0_cross_P1 = (control_points[control_point_idx_0 +
				facet_idx * control_point_size_per_facet] % control_points[control_point_idx_1 +
				facet_idx * control_point_size_per_facet]);
				int index = control_points.size() + facet_idx * control_point_pairs_per_facet +
				control_point_idx_0 * control_point_size_per_facet + control_point_idx_1;
				int j = map_var_to_represent_var[index];
				origin_x[3 * j] = P0_cross_P1[0] * map_var_to_represent_var_sign[index];
				origin_x[3 * j + 1] = P0_cross_P1[1] * map_var_to_represent_var_sign[index];
				origin_x[3 * j + 2] = P0_cross_P1[2] * map_var_to_represent_var_sign[index];
			}
		}
	}

	//check origin_x
	for (const auto& vh : mesh.vertices())
	{
		Mesh::Point deformed_point(0, 0, 0);
		for (int i = 0; i < control_points.size(); ++i) {
			int j = map_var_to_represent_var[i];
			int index = vh.idx()*control_points.size() + i;
			Mesh::Point control_point = { origin_x[3 * j], origin_x[3 * j + 1], origin_x[3 * j + 2] };

			deformed_point += control_point *
				map_vertex_idx_and_control_point_to_deformation_para[index];
		}
		for (int i = 0; i < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++i) {
			int j = map_var_to_represent_var[control_points.size() + i];
			int index = vh.idx()*facet_and_control_points_pairs.size() * control_point_pairs_per_facet + i;
			Mesh::Point control_point_pair = { origin_x[3 * j], origin_x[3 * j + 1], origin_x[3 * j + 2] };

			deformed_point += control_point_pair * map_var_to_represent_var_sign[control_points.size() + i] *
				map_vertex_idx_and_control_point_pair_to_deformation_para[index];
		}
		std::cout << vh.idx() << " before:" << mesh.point(vh) << " after:" << deformed_point << "\n";
	}

	Eigen::VectorXd origin_x_star = global_solver_mat * origin_x;
		auto& origin_delta = origin_x_star - global_solver_vec;
		double origin_E_total = origin_delta.squaredNorm();
		double origin_E_position = 0, origin_E_orientation = 0, origin_E_rigidity = 0, origin_E_smoothness = 0, origin_E_fixed_vertices = 0;
		for (int i = 0; i < section_1; ++i) {
			origin_E_position += origin_delta(i)*origin_delta(i);
		}
		for (int i = section_1; i < section_2; ++i) {
			origin_E_orientation += origin_delta(i)*origin_delta(i);
		}
		for (int i = section_2; i < section_3; ++i) {
			origin_E_rigidity += origin_delta(i)*origin_delta(i);
		}
		for (int i = section_3; i < section_4; ++i) {
			origin_E_smoothness += origin_delta(i)*origin_delta(i);
		}
		for (int i = section_4; i < origin_delta.size(); ++i) {
			origin_E_fixed_vertices += origin_delta(i)*origin_delta(i);
		}
		std::cout << "origin_total energy=" << origin_E_total << " avg origin_E_position=" << origin_E_position / position_constraints.size()
			<< "avg origin_E_orientation=" << origin_E_orientation / orientation_constraints.size()
			<< " avg origin_E_rigidity=" << origin_E_rigidity / rigidity_constraints.size()
			<< " avg origin_E_smoothness=" << origin_E_smoothness / smoothness_constraints.size()
			<< " avg origin_E_fixed_vertices=" << origin_E_fixed_vertices / fixed_vertices.size() << "\n";*/

#if DEBUG_MODE
	std::ofstream optimized_var_out("optimized_var.txt");
	for (int i = 0; i < var_size / 3; ++i) {
		optimized_var_out << x[3 * i] << " " << x[3 * i + 1] << " " << x[3 * i + 2] << "\n";
	}
	optimized_var_out.close();

	std::ofstream optimized_control_points_out("optimized_control_points.obj");
	for (int i = 0; i < represent_control_points.size(); ++i) {
		represent_control_points[i] = { x[3 * i], x[3 * i + 1], x[3 * i + 2] };
		optimized_control_points_out << "v " << represent_control_points[i] << "\n";
	}
	optimized_control_points_out.close();

	std::ofstream optimized_control_point_pairs_out("optimized_control_point_pairs.obj");
	for (int i = represent_control_points.size(); i < var_size / 3; ++i) {
		optimized_control_point_pairs_out << "v " << x[3 * i] << " " << x[3 * i + 1] << " " << x[3 * i + 2] << "\n";
	}
	optimized_control_point_pairs_out.close();

	total_hessian_norm = 0;
	double max_hessian = 0;
	std::ofstream optimized_smoothness_sample_points_out("optimized_smoothness_sample_points.obj");
	for (int vertex_idx = 0; vertex_idx < smoothness_constraints.size(); ++vertex_idx) {
		Mesh::Point deformed_point(0, 0, 0);
		Tensor hessian = {};
		for (int i = 0; i < control_points.size(); ++i) {
			int j = map_var_to_represent_var[i];
			int index = vertex_idx * control_points.size() + i;
			Mesh::Point control_point = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point *
				map_vertex_idx_and_control_point_to_smoothness_deformation_para[index];
			sum_point_product(hessian, control_point, map_vertex_idx_and_control_point_to_smoothness_hessian_deformation_para[index]);
		}
		for (int i = 0; i < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++i) {
			int j = map_var_to_represent_var[control_points.size() + i];
			int index = vertex_idx * facet_and_control_points_pairs.size() * control_point_pairs_per_facet + i;
			Mesh::Point control_point_pair = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point_pair * map_var_to_represent_var_sign[control_points.size() + i] *
				map_vertex_idx_and_control_point_pair_to_smoothness_deformation_para[index];
			sum_point_product(hessian, control_point_pair, map_vertex_idx_and_control_point_pair_to_smoothness_hessian_deformation_para[index]
				* map_var_to_represent_var_sign[control_points.size() + i]);
		}
		optimized_smoothness_sample_points_out << "v " << deformed_point << "\n";
		/*std::cout << "hessian:";
		printTensor(hessian);
		std::cout << "\n";*/
		double hessian_norm = 0;
		for (const auto& matrix : hessian) {
			for (const auto& row : matrix) {
				for (const auto& element : row) {
					hessian_norm += element * element;
				}
			}
		}
		total_hessian_norm += hessian_norm;
		if (hessian_norm > max_hessian) {
			max_hessian = hessian_norm;
		}
	}
	optimized_smoothness_sample_points_out.close();
	std::cout << "smoothness avg hessian norm:" << total_hessian_norm / smoothness_constraints.size() << " max hessian norm:" << max_hessian << "\n";

	std::ofstream optimized_rigidity_sample_points_out("optimized_polyline.obj");
	for (int vertex_idx = 0; vertex_idx < rigidity_constraints.size(); ++vertex_idx) {
		Mesh::Point deformed_point(0, 0, 0);
		Eigen::Matrix3d gradient = Eigen::Matrix3d::Zero();
		for (int i = 0; i < control_points.size(); ++i) {
			int j = map_var_to_represent_var[i];
			int index = vertex_idx * control_points.size() + i;
			Mesh::Point control_point = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point *
				map_vertex_idx_and_control_point_to_rigidity_deformation_para[index];
			gradient += convertToEigen(control_point)*
				convertToEigen(map_vertex_idx_and_control_point_to_rigidity_gradient_deformation_para[index]).transpose();
		}
		for (int i = 0; i < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++i) {
			int j = map_var_to_represent_var[control_points.size() + i];
			int index = vertex_idx * facet_and_control_points_pairs.size() * control_point_pairs_per_facet + i;
			Mesh::Point control_point_pair = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point_pair * map_var_to_represent_var_sign[control_points.size() + i] *
				map_vertex_idx_and_control_point_pair_to_rigidity_deformation_para[index];
			gradient += convertToEigen(control_point_pair) *
				convertToEigen(map_vertex_idx_and_control_point_pair_to_rigidity_gradient_deformation_para[index]).transpose()
				* map_var_to_represent_var_sign[control_points.size() + i];
		}
		optimized_rigidity_sample_points_out << "v " << deformed_point << "\n";
		//std::cout << "rigidity gradient:" << gradient << "\n";
	}
	optimized_rigidity_sample_points_out.close();

	total_hessian_norm = 0;
	int max_hessian_idx = -1;
	max_hessian = 0;
#endif

	for (const auto& vh : mesh.vertices())
	{
		Mesh::Point deformed_point(0, 0, 0);
#if DEBUG_MODE
		Eigen::Matrix3d gradient = Eigen::Matrix3d::Zero();
		Tensor hessian = {};
#endif
		for (int i = 0; i < control_points.size(); ++i) {
			int j = map_var_to_represent_var[i];
			int index = vh.idx()*control_points.size() + i;
			Mesh::Point control_point = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point *
				map_vertex_idx_and_control_point_to_deformation_para[index];
#if DEBUG_MODE
			gradient += convertToEigen(control_point)*
				convertToEigen(map_vertex_idx_and_control_point_to_gradient_deformation_para[index]).transpose();
			sum_point_product(hessian, control_point, map_vertex_idx_and_control_point_to_hessian_deformation_para[index]);
#endif
		}
		for (int i = 0; i < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++i) {
			int j = map_var_to_represent_var[control_points.size() + i];
			int index = vh.idx()*facet_and_control_points_pairs.size() * control_point_pairs_per_facet + i;
			Mesh::Point control_point_pair = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point_pair * map_var_to_represent_var_sign[control_points.size() + i] *
				map_vertex_idx_and_control_point_pair_to_deformation_para[index];
#if DEBUG_MODE
			gradient += convertToEigen(control_point_pair) *
			convertToEigen(map_vertex_idx_and_control_point_pair_to_gradient_deformation_para[index]).transpose()
				* map_var_to_represent_var_sign[control_points.size() + i];
			sum_point_product(hessian, control_point_pair, map_vertex_idx_and_control_point_pair_to_hessian_deformation_para[index]
				* map_var_to_represent_var_sign[control_points.size() + i]);
#endif
		}
		/*std::cout << vh.idx() << " before:" << mesh.point(vh) << " after:" << deformed_point << "\n";
		std::cout << "gradient:" << gradient << "\n";
		std::cout << "hessian:";
		printTensor(hessian);
		std::cout << "\n";*/

		mesh.set_point(vh, deformed_point);
#if DEBUG_MODE
		double hessian_norm = 0;
		for (const auto& matrix : hessian) {
			for (const auto& row : matrix) {
				for (const auto& element : row) {
					hessian_norm += element * element;
				}
			}
		}
		total_hessian_norm += hessian_norm;
		if (hessian_norm > max_hessian) {
			max_hessian = hessian_norm;
			max_hessian_idx = vh.idx();
		}
		if (hessian_norm > 100)std::cout << "idx:" << vh.idx() << " hessian norm:" << hessian_norm << "\n";
#endif
	}
	updateGL();
	OpenMesh::IO::write_mesh(mesh, "deformed_mesh.obj");
#if DEBUG_MODE
	std::cout << "mesh avg hessian norm:" << total_hessian_norm / mesh.n_vertices() << " max hessian norm:" << max_hessian
		<< " idx:" << max_hessian_idx << "\n";
#endif
	std::cout << "Generate deformation completed!" << "\n";
}

void PGC::do_read_integral_para() {
	std::string filename;
	if (cage_mode == false) {
		filename = "quad";
	}
	else {
		filename = "tri";
	}

	/*std::ifstream dirichlet_in("dirichlet_para_" + filename + ".txt");
	std::ifstream neumann_in("neumann_para_" + filename + ".txt");
	std::ifstream dirichlet_gradient_in("dirichlet_gradient_para_" + filename + ".txt");
	std::ifstream neumann_gradient_in("neumann_gradient_para_" + filename + ".txt");
	std::ifstream dirichlet_hessian_in("dirichlet_hessian_para_" + filename + ".txt");
	std::ifstream neumann_hessian_in("neumann_hessian_para_" + filename + ".txt");*/

	std::ifstream control_point_deformation_para_in("control_point_deformation_para_" + filename + ".txt");
	std::ifstream control_point_pair_deformation_para_in("control_point_pair_deformation_para_" + filename + ".txt");
	/*std::ifstream control_point_gradient_deformation_para_in("control_point_gradient_deformation_para_" + filename + ".txt");
	std::ifstream control_point_pair_gradient_deformation_para_in("control_point_pair_gradient_deformation_para_" + filename + ".txt");
	std::ifstream control_point_hessian_deformation_para_in("control_point_hessian_deformation_para_" + filename + ".txt");
	std::ifstream control_point_pair_hessian_deformation_para_in("control_point_pair_hessian_deformation_para_" + filename + ".txt");*/
	std::string line;

	/*while (std::getline(dirichlet_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		double number;
		while (iss >> number) {
			numbers.push_back(number);
		}
		DirichletParametersForMeshPoints[{ numbers[0], numbers[1] }][{ numbers[2], numbers[3] }] = numbers[4];
	}
	dirichlet_in.close();

	while (std::getline(neumann_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		double number;
		while (iss >> number) {
			numbers.push_back(number);
		}
		NeumannParametersForMeshPoints[{ numbers[0], numbers[1] }][{ numbers[2], numbers[3] }] = numbers[4];
	}
	neumann_in.close();

	while (std::getline(dirichlet_gradient_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		double number;
		while (iss >> number) {
			numbers.push_back(number);
		}
		DirichletGradientParametersForMeshPoints[{ numbers[0], numbers[1] }][{ numbers[2], numbers[3] }] =
		{ numbers[4],numbers[5],numbers[6] };
	}
	dirichlet_gradient_in.close();

	while (std::getline(neumann_gradient_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		double number;
		while (iss >> number) {
			numbers.push_back(number);
		}
		NeumannGradientParametersForMeshPoints[{ numbers[0], numbers[1] }][{ numbers[2], numbers[3] }] =
		{ numbers[4],numbers[5],numbers[6] };
	}
	neumann_gradient_in.close();

	while (std::getline(dirichlet_hessian_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		double number;
		while (iss >> number) {
			numbers.push_back(number);
		}
		DirichletHessianParametersForMeshPoints[{ numbers[0], numbers[1] }][{ numbers[2], numbers[3] }] =
			Eigen::Map<Eigen::Matrix3d>(numbers.data() + 4);
	}
	dirichlet_hessian_in.close();

	while (std::getline(neumann_hessian_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		double number;
		while (iss >> number) {
			numbers.push_back(number);
		}
		NeumannHessianParametersForMeshPoints[{ numbers[0], numbers[1] }][{ numbers[2], numbers[3] }] =
			Eigen::Map<Eigen::Matrix3d>(numbers.data() + 4);
	}
	neumann_hessian_in.close();*/

	double number;
	while (control_point_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_to_deformation_para.push_back(number);
	}
	control_point_deformation_para_in.close();

	while (control_point_pair_deformation_para_in >> number) {
		map_vertex_idx_and_control_point_pair_to_deformation_para.push_back(number);
	}
	control_point_pair_deformation_para_in.close();

	//gradient
	/*while (std::getline(control_point_gradient_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;
		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_to_gradient_deformation_para.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	control_point_gradient_deformation_para_in.close();

	while (std::getline(control_point_pair_gradient_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;
		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_pair_to_gradient_deformation_para.push_back({ numbers[0],numbers[1] ,numbers[2] });
	}
	control_point_pair_gradient_deformation_para_in.close();

	//hessian
	while (std::getline(control_point_hessian_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_to_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
	}
	control_point_hessian_deformation_para_in.close();

	while (std::getline(control_point_pair_hessian_deformation_para_in, line)) {
		std::istringstream iss(line);
		std::vector<double> numbers;

		while (iss >> number) {
			numbers.push_back(number);
		}
		map_vertex_idx_and_control_point_pair_to_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
	}
	control_point_pair_hessian_deformation_para_in.close();*/
	std::cout << "read integral para complete!" << "\n";
}

void PGC::do_deformation() {
	//double total_hessian_norm = 0;
	//int max_hessian_idx = -1;
	//double max_hessian = 0;
#if DEBUG_MODE
	double max_diff = 0;
	double avg_diff = 0;
#endif

	std::clock_t start, end;
	start = std::clock();
	for (const auto& vh : mesh.vertices())
	{
		vertex_idx = vh.idx();
		Mesh::Point deformed_point(0, 0, 0);
		//Eigen::Matrix3d gradient = Eigen::Matrix3d::Zero();
		//Tensor hessian = {};

		for (int control_point_idx = 0; control_point_idx < control_points.size(); ++control_point_idx) {
			int index = vertex_idx * control_points.size() + control_point_idx;

			deformed_point += control_points[control_point_idx] * map_vertex_idx_and_control_point_to_deformation_para[index];
			//gradient += convertToEigen(control_points[control_point_idx])*
			//	convertToEigen(map_vertex_idx_and_control_point_to_gradient_deformation_para[index]).transpose();
			//sum_point_product(hessian, control_points[control_point_idx], map_vertex_idx_and_control_point_to_hessian_deformation_para[index]);
		}

		for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
			for (int control_point_idx_0 = 0; control_point_idx_0 < control_point_size_per_facet; ++control_point_idx_0) {
				for (int control_point_idx_1 = 0; control_point_idx_1 < control_point_size_per_facet; ++control_point_idx_1) {
					int index = vertex_idx * facet_and_control_points_pairs.size() * control_point_size_per_facet*control_point_size_per_facet
						+ facet_idx * control_point_size_per_facet*control_point_size_per_facet +
						control_point_idx_0 * control_point_size_per_facet + control_point_idx_1;
					auto& P0_cross_P1 = (control_points[control_point_idx_0 + facet_idx * control_point_size_per_facet] %
						control_points[control_point_idx_1 + facet_idx * control_point_size_per_facet]);

					deformed_point += P0_cross_P1 * map_vertex_idx_and_control_point_pair_to_deformation_para[index];
					//gradient += convertToEigen(P0_cross_P1) * convertToEigen(map_vertex_idx_and_control_point_pair_to_gradient_deformation_para[index]).transpose();
					//sum_point_product(hessian, P0_cross_P1, map_vertex_idx_and_control_point_pair_to_hessian_deformation_para[index]);
				}
			}
		}
		//std::cout << "before:" << mesh.point(vh) << " after:" << deformed_point << "\n";
#if DEBUG_MODE
		double diff = (mesh.point(vh) - deformed_point).length();
		if (max_diff < diff) max_diff = diff;
		avg_diff += diff;
#endif
		//
		//gradient -= Eigen::Matrix3d::Identity();
		//double hessian_norm = 0;
		//for (const auto& matrix : hessian) {
		//	for (const auto& row : matrix) {
		//		for (const auto& element : row) {
		//			hessian_norm += element * element;
		//		}
		//	}
		//}
		//total_hessian_norm += hessian_norm;
		//if (hessian_norm > max_hessian) {
		//	max_hessian = hessian_norm;
		//	max_hessian_idx = vh.idx();
		//}
		//double gradient_norm = gradient.norm();
		//if (gradient_norm > 1)std::cout << "idx:" << vh.idx() << " gradient norm:" << gradient_norm << "\n";
		//if (hessian_norm > 1)std::cout << "idx:" << vh.idx() << " hessian norm:" << hessian_norm << "\n";
		mesh.set_point(vh, deformed_point);
	}
	//std::cout << "mesh avg hessian norm:" << total_hessian_norm / mesh.n_vertices() << " max hessian norm:" << max_hessian
	//	<< " idx:" << max_hessian_idx << "\n";
#if DEBUG_MODE
	std::cout << "max diff*10e6:" << max_diff * 10e6 << " avg_diff*10e6:" << avg_diff / mesh.n_vertices()*10e6 << "\n";
#endif
	end = std::clock();
	std::cout << "deformation time:" << double(end - start) / CLOCKS_PER_SEC * 1000 << "ms" << "\n";

	updateGL();
	std::cout << "Generate deformation completed!" << "\n";
}

void PGC::do_variational_deformation() {
	int control_point_pairs_per_facet = control_point_size_per_facet * control_point_size_per_facet;
	if (!variational_is_initialed) {
		loadBDCSVD();
		loadglobalsolver();
		read_position_constraints();
		read_orientation_constraints();
		read_rigidity_constraints();
		read_smoothness_constraints();

		std::string filename;
		if (cage_mode == false) {
			filename = "quad";
		}
		else {
			filename = "tri";
		}

		std::string line;
		double number;
		map_vertex_idx_and_control_point_to_deformation_para.clear();
		map_vertex_idx_and_control_point_pair_to_deformation_para.clear();
		std::ifstream control_point_deformation_para_in("control_point_deformation_para_" + filename + ".txt");
		std::ifstream control_point_pair_deformation_para_in("control_point_pair_deformation_para_" + filename + ".txt");
		while (control_point_deformation_para_in >> number) {
			map_vertex_idx_and_control_point_to_deformation_para.push_back(number);
		}
		control_point_deformation_para_in.close();

		while (control_point_pair_deformation_para_in >> number) {
			map_vertex_idx_and_control_point_pair_to_deformation_para.push_back(number);
		}
		control_point_pair_deformation_para_in.close();

#if DEBUG_MODE
		//read hessian paras for debug
		std::ifstream control_point_hessian_deformation_para_in("control_point_hessian_deformation_para_" + filename + ".txt");
		std::ifstream control_point_pair_hessian_deformation_para_in("control_point_pair_hessian_deformation_para_" + filename + ".txt");
		while (std::getline(control_point_hessian_deformation_para_in, line)) {
			std::istringstream iss(line);
			std::vector<double> numbers;

			while (iss >> number) {
				numbers.push_back(number);
			}
			map_vertex_idx_and_control_point_to_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
		}
		control_point_hessian_deformation_para_in.close();

		while (std::getline(control_point_pair_hessian_deformation_para_in, line)) {
			std::istringstream iss(line);
			std::vector<double> numbers;

			while (iss >> number) {
				numbers.push_back(number);
			}
			map_vertex_idx_and_control_point_pair_to_hessian_deformation_para.push_back(Eigen::Map<Eigen::Matrix3d>(numbers.data()));
		}
		control_point_pair_hessian_deformation_para_in.close();
#endif

		map_var_to_represent_var.resize(control_points.size() + facet_and_control_points_pairs.size() * control_point_pairs_per_facet, -1);
		map_var_to_represent_var_sign.resize(control_points.size() + facet_and_control_points_pairs.size() * control_point_pairs_per_facet, 1);
		var_size = represent_control_points.size();
		for (int i = 0; i < control_points.size(); ++i) {
			map_var_to_represent_var[i] = map_control_point_to_represent_control_point_idx[i];
		}
		not_appeared_flag = map_var_to_represent_var.size();

		represent_control_pair_is_appeared.resize(represent_control_points.size()*represent_control_points.size(), not_appeared_flag);
		for (int facet_idx = 0; facet_idx < facet_and_control_points_pairs.size(); ++facet_idx) {
			for (int control_point_idx_0 = 0; control_point_idx_0 < control_point_size_per_facet; ++control_point_idx_0) {
				for (int control_point_idx_1 = 0; control_point_idx_1 < control_point_size_per_facet; ++control_point_idx_1) {
					int index = facet_idx * control_point_pairs_per_facet + control_point_idx_0 * control_point_size_per_facet + control_point_idx_1;
					int represent_control_point_idx_0 = map_control_point_to_represent_control_point_idx
						[control_point_idx_0 + facet_idx * control_point_size_per_facet];
					int represent_control_point_idx_1 = map_control_point_to_represent_control_point_idx
						[control_point_idx_1 + facet_idx * control_point_size_per_facet];

					if (represent_control_point_idx_0 > represent_control_point_idx_1) {
						int i = represent_control_point_idx_0 * represent_control_points.size() + represent_control_point_idx_1;
						if (represent_control_pair_is_appeared[i] == not_appeared_flag) {
							map_var_to_represent_var[control_points.size() + index] = var_size;
							represent_control_pair_is_appeared[i] = var_size;
							var_size++;
						}
						else {
							map_var_to_represent_var[control_points.size() + index] = represent_control_pair_is_appeared[i];
						}
					}
					else if (represent_control_point_idx_0 < represent_control_point_idx_1) {
						int i = represent_control_point_idx_1 * represent_control_points.size() + represent_control_point_idx_0;
						if (represent_control_pair_is_appeared[i] == not_appeared_flag) {
							map_var_to_represent_var[control_points.size() + index] = var_size;
							map_var_to_represent_var_sign[control_points.size() + index] = -1;
							represent_control_pair_is_appeared[i] = var_size;
							var_size++;
						}
						else {
							map_var_to_represent_var[control_points.size() + index] = represent_control_pair_is_appeared[i];
							map_var_to_represent_var_sign[control_points.size() + index] = -1;
						}
					}
					else {
						map_var_to_represent_var[control_points.size() + index] = 0;
						map_var_to_represent_var_sign[control_points.size() + index] = 0;
					}
				}
			}
		}

		equation_size = 3 * position_constraints.size() + 9 * orientation_constraints.size() + 9 * rigidity_constraints.size() +
			15 * smoothness_constraints.size() + 3 * fixed_vertices.size();
		rigidity_start = 3 * position_constraints.size() + 9 * orientation_constraints.size();

		global_solver_rotation_vec.resize(equation_size);
		global_solver_rotation_vec.setZero();

		for (int i = 0; i < rigidity_constraints.size(); ++i) {
			Eigen::Vector3d axis = Eigen::Vector3d::Random().normalized();
			double angle = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
			Eigen::AngleAxisd angleAxis(angle, axis);
			auto& rotation_mat = angleAxis.toRotationMatrix();
			global_solver_rotation_vec[rigidity_start + 9 * i] = rotation_mat(0, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 1] = rotation_mat(1, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 2] = rotation_mat(2, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 3] = rotation_mat(0, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 4] = rotation_mat(1, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 5] = rotation_mat(2, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 6] = rotation_mat(0, 2);
			global_solver_rotation_vec[rigidity_start + 9 * i + 7] = rotation_mat(1, 2);
			global_solver_rotation_vec[rigidity_start + 9 * i + 8] = rotation_mat(2, 2);
		}

		section_1 = 3 * position_constraints.size();
		section_2 = section_1 + 9 * orientation_constraints.size();
		section_3 = section_2 + 9 * rigidity_constraints.size();
		section_4 = section_3 + 15 * smoothness_constraints.size();

		variational_is_initialed = true;
		std::cout << "variation initial complete\n";
	}

	for (int i = 0; i < position_constraints.size(); ++i) {
		global_solver_vec(3 * i) = position_weight * position_constraints[i].second[0];
		global_solver_vec(3 * i + 1) = position_weight * position_constraints[i].second[1];
		global_solver_vec(3 * i + 2) = position_weight * position_constraints[i].second[2];
	}

#if DEBUG_MODE
	//compute the mesh generated by previous optimized var 
	std::ifstream vec_in("optimized_var.txt");
	double number;
	std::vector<double> numbers;
	while (vec_in >> number) {
		numbers.push_back(number);
	}
	x = Eigen::VectorXd::Map(numbers.data(), numbers.size());

	vec_in.close();
	for (const auto& vh : mesh.vertices())
	{
		Mesh::Point deformed_point(0, 0, 0);
		for (int i = 0; i < control_points.size(); ++i) {
			int j = map_var_to_represent_var[i];
			int index = vh.idx()*control_points.size() + i;
			Mesh::Point control_point = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point *
				map_vertex_idx_and_control_point_to_deformation_para[index];
		}
		for (int i = 0; i < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++i) {
			int j = map_var_to_represent_var[control_points.size() + i];
			int index = vh.idx()*facet_and_control_points_pairs.size() * control_point_pairs_per_facet + i;
			Mesh::Point control_point_pair = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point_pair * map_var_to_represent_var_sign[control_points.size() + i] *
				map_vertex_idx_and_control_point_pair_to_deformation_para[index];
		}
		mesh.set_point(vh, deformed_point);
	}
	OpenMesh::IO::write_mesh(mesh, "previous_variational_deformed_mesh.obj");
#endif

	energy = DBL_MAX;
	energy_diff = DBL_MAX;
	iter = 0;

	// Local/Global
	std::clock_t start, end;
	start = std::clock();
	while ((energy_diff > 1e-7 || iter < 5) && iter < 1000) {
		std::cout << "Iteration " << iter + 1 << "\n";
		//Global
		x = SVD_mat * (global_solver_vec + global_solver_rotation_vec);
		Eigen::VectorXd x_star = global_solver_mat * x;

		//Local
		for (int i = 0; i < rigidity_constraints.size(); ++i) {
			Eigen::Matrix3d gradient = Eigen::Map<Eigen::Matrix3d>(x_star.data() + rigidity_start + 9 * i) / rigidity_weight;
			auto& rotation_mat = optimalRotationMatrix(gradient)*rigidity_weight;

			global_solver_rotation_vec[rigidity_start + 9 * i] = rotation_mat(0, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 1] = rotation_mat(1, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 2] = rotation_mat(2, 0);
			global_solver_rotation_vec[rigidity_start + 9 * i + 3] = rotation_mat(0, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 4] = rotation_mat(1, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 5] = rotation_mat(2, 1);
			global_solver_rotation_vec[rigidity_start + 9 * i + 6] = rotation_mat(0, 2);
			global_solver_rotation_vec[rigidity_start + 9 * i + 7] = rotation_mat(1, 2);
			global_solver_rotation_vec[rigidity_start + 9 * i + 8] = rotation_mat(2, 2);
		}

		auto& delta = x_star - global_solver_vec - global_solver_rotation_vec;
		double E_total = delta.squaredNorm();
		double E_position = 0, E_orientation = 0, E_rigidity = 0, E_smoothness = 0, E_fixed_vertices = 0;
		for (int i = 0; i < section_1; ++i) {
			E_position += delta(i)*delta(i);
		}
		for (int i = section_1; i < section_2; ++i) {
			E_orientation += delta(i)*delta(i);
		}
		for (int i = section_2; i < section_3; ++i) {
			E_rigidity += delta(i)*delta(i);
		}
		for (int i = section_3; i < section_4; ++i) {
			E_smoothness += delta(i)*delta(i);
		}
		for (int i = section_4; i < delta.size(); ++i) {
			E_fixed_vertices += delta(i)*delta(i);
		}

		energy_diff = std::abs(energy - E_total);
		energy = E_total;
		std::cout << "total energy=" << E_total << " avg E_position=" << E_position / position_constraints.size()
			<< "avg E_orientation=" << E_orientation / orientation_constraints.size()
			<< " avg E_rigidity=" << E_rigidity / rigidity_constraints.size()
			<< " avg E_smoothness=" << E_smoothness / smoothness_constraints.size()
			<< " avg E_fixed_vertices=" << E_fixed_vertices / fixed_vertices.size()
			<< " energy_diff=" << energy_diff << "\n";

		iter++;
	}
	end = std::clock();
	std::cout << "Complete Optimizing, spend " << double(end - start) / CLOCKS_PER_SEC * 1000 << "ms" << "\n";
	std::cout << "average time cost per step:" << double(end - start) / CLOCKS_PER_SEC * 1000 / (iter + 1) << "ms" << "\n";

	for (const auto& vh : mesh.vertices())
	{
		Mesh::Point deformed_point(0, 0, 0);
#if DEBUG_MODE
		Tensor hessian = {};
#endif
		for (int i = 0; i < control_points.size(); ++i) {
			int j = map_var_to_represent_var[i];
			int index = vh.idx()*control_points.size() + i;
			Mesh::Point control_point = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point *
				map_vertex_idx_and_control_point_to_deformation_para[index];
#if DEBUG_MODE
			sum_point_product(hessian, control_point, map_vertex_idx_and_control_point_to_hessian_deformation_para[index]);
#endif
		}
		for (int i = 0; i < facet_and_control_points_pairs.size() * control_point_pairs_per_facet; ++i) {
			int j = map_var_to_represent_var[control_points.size() + i];
			int index = vh.idx()*facet_and_control_points_pairs.size() * control_point_pairs_per_facet + i;
			Mesh::Point control_point_pair = { x[3 * j], x[3 * j + 1], x[3 * j + 2] };

			deformed_point += control_point_pair * map_var_to_represent_var_sign[control_points.size() + i] *
				map_vertex_idx_and_control_point_pair_to_deformation_para[index];
#if DEBUG_MODE
			sum_point_product(hessian, control_point_pair, map_vertex_idx_and_control_point_pair_to_hessian_deformation_para[index]
				* map_var_to_represent_var_sign[control_points.size() + i]);
#endif
		}
		mesh.set_point(vh, deformed_point);

#if DEBUG_MODE
		double hessian_norm = 0;
		for (const auto& matrix : hessian) {
			for (const auto& row : matrix) {
				for (const auto& element : row) {
					hessian_norm += element * element;
				}
			}
		}
		if (hessian_norm > 10)std::cout << "idx:" << vh.idx() << " hessian norm:" << hessian_norm << "\n";
#endif
	}
	updateGL();
	std::cout << "Generate deformation completed!" << "\n";
}

Mesh::Point PGC::DeformedPositionByFacet(const int vertex_idx, const int facet_idx, const std::vector<Mesh::Point>& control_point) {
	Mesh::Point dirichlet_part = { 0,0,0 }, neumann_part = { 0,0,0 };

	auto& DirichletParameter = DirichletParametersForMeshPoints.at({ vertex_idx,facet_idx });
	auto& NeumannParameter = NeumannParametersForMeshPoints.at({ vertex_idx,facet_idx });
	if (cage_mode == false) {
		dirichlet_part += control_point[0] * DirichletParameter[{0, 0}];
		dirichlet_part += (-control_point[0] + control_point[1]) * DirichletParameter[{1, 0}];
		dirichlet_part += (-control_point[0] + control_point[3]) *  DirichletParameter[{0, 1}];
		dirichlet_part += (control_point[0] - control_point[1] + control_point[2] - control_point[3]) *  DirichletParameter[{1, 1}];

		neumann_part += ((-control_point[0] + control_point[1]) % (-control_point[0] + control_point[3]))*NeumannParameter[{0, 0}];
		neumann_part += ((-control_point[0] + control_point[1]) % (-control_point[3] + control_point[2])) *NeumannParameter[{1, 0}];
		neumann_part += ((-control_point[1] + control_point[2]) % (-control_point[0] + control_point[3]))*NeumannParameter[{0, 1}];
	}
	else {
		auto& P003 = control_point[0];
		auto& P012 = control_point[1];
		auto& P021 = control_point[2];
		auto& P030 = control_point[3];
		auto& P102 = control_point[4];
		auto& P111 = control_point[5];
		auto& P120 = control_point[6];
		auto& P201 = control_point[7];
		auto& P210 = control_point[8];
		auto& P300 = control_point[9];

		std::map<std::pair<int, int>, Mesh::Point> P;
		P[{0, 0}] = P003;
		P[{1, 0}] = -3 * P003 + 3 * P012;
		P[{0, 1}] = -3 * P003 + 3 * P102;
		P[{2, 0}] = 3 * P003 - 6 * P012 + 3 * P021;
		P[{1, 1}] = 6 * P003 - 6 * P012 - 6 * P102 + 6 * P111;
		P[{0, 2}] = 3 * P003 - 6 * P102 + 3 * P201;
		P[{3, 0}] = -P003 + 3 * P012 - 3 * P021 + P030;
		P[{2, 1}] = -3 * P003 + 6 * P012 - 3 * P021 + 3 * P102 - 6 * P111 + 3 * P120;
		P[{1, 2}] = -3 * P003 + 3 * P012 + 6 * P102 - 6 * P111 - 3 * P201 + 3 * P210;
		P[{0, 3}] = -P003 + 3 * P102 - 3 * P201 + P300;

		for (const auto& term : P) {
			dirichlet_part += term.second * DirichletParameter[{term.first.first, term.first.second}];
		}

		std::map<std::pair<int, int>, Mesh::Point> P_u_cross_P_v;
		for (int i = 0; i <= 4; ++i) {
			for (int j = 0; j <= 4 - i; ++j) {
				P_u_cross_P_v[{i, j}] = { 0,0,0 };
			}
		}
		for (const auto& term1 : P) {
			if (term1.first.first == 0)continue;
			for (const auto& term2 : P) {
				if (term2.first.second == 0)continue;
				P_u_cross_P_v[{term1.first.first + term2.first.first - 1, term1.first.second + term2.first.second - 1}] +=
					term1.first.first*term2.first.second*(term1.second%term2.second);
			}
		}

		for (const auto& term : P_u_cross_P_v) {
			neumann_part += term.second * NeumannParameter[{term.first.first, term.first.second}];
		}
	}

	return (dirichlet_part + neumann_part) / 4 / M_PI;
}

Eigen::Matrix3d PGC::ComputeGradientByFacet(const int point_idx, const int facet_idx, const std::vector<Mesh::Point>& control_point) {
	Eigen::Matrix3d dirichlet_part, neumann_part;
	dirichlet_part.setZero();
	neumann_part.setZero();

	auto& DirichletGradientParameter = DirichletGradientParametersForMeshPoints.at({ point_idx,facet_idx });
	auto& NeumannGradientParameter = NeumannGradientParametersForMeshPoints.at({ point_idx,facet_idx });
	if (cage_mode == false) {
		dirichlet_part += convertToEigen(control_point[0]) * convertToEigen(DirichletGradientParameter[{0, 0}]).transpose();
		dirichlet_part += convertToEigen(-control_point[0] + control_point[1]) * convertToEigen(DirichletGradientParameter[{1, 0}]).transpose();
		dirichlet_part += convertToEigen(-control_point[0] + control_point[3]) *  convertToEigen(DirichletGradientParameter[{0, 1}]).transpose();
		dirichlet_part += convertToEigen(control_point[0] - control_point[1] + control_point[2] - control_point[3]) *
			convertToEigen(DirichletGradientParameter[{1, 1}]).transpose();

		neumann_part += convertToEigen((-control_point[0] + control_point[1]) % (-control_point[0] + control_point[3]))*
			convertToEigen(NeumannGradientParameter[{0, 0}]).transpose();
		neumann_part += convertToEigen((-control_point[0] + control_point[1]) % (-control_point[3] + control_point[2])) *
			convertToEigen(NeumannGradientParameter[{1, 0}]).transpose();
		neumann_part += convertToEigen((-control_point[1] + control_point[2]) % (-control_point[0] + control_point[3]))*
			convertToEigen(NeumannGradientParameter[{0, 1}]).transpose();
	}
	else {
		auto& P003 = control_point[0];
		auto& P012 = control_point[1];
		auto& P021 = control_point[2];
		auto& P030 = control_point[3];
		auto& P102 = control_point[4];
		auto& P111 = control_point[5];
		auto& P120 = control_point[6];
		auto& P201 = control_point[7];
		auto& P210 = control_point[8];
		auto& P300 = control_point[9];

		std::map<std::pair<int, int>, Mesh::Point> P;
		P[{0, 0}] = P003;
		P[{1, 0}] = -3 * P003 + 3 * P012;
		P[{0, 1}] = -3 * P003 + 3 * P102;
		P[{2, 0}] = 3 * P003 - 6 * P012 + 3 * P021;
		P[{1, 1}] = 6 * P003 - 6 * P012 - 6 * P102 + 6 * P111;
		P[{0, 2}] = 3 * P003 - 6 * P102 + 3 * P201;
		P[{3, 0}] = -P003 + 3 * P012 - 3 * P021 + P030;
		P[{2, 1}] = -3 * P003 + 6 * P012 - 3 * P021 + 3 * P102 - 6 * P111 + 3 * P120;
		P[{1, 2}] = -3 * P003 + 3 * P012 + 6 * P102 - 6 * P111 - 3 * P201 + 3 * P210;
		P[{0, 3}] = -P003 + 3 * P102 - 3 * P201 + P300;

		for (const auto& term : P) {
			dirichlet_part += convertToEigen(term.second) * convertToEigen(DirichletGradientParameter[{term.first.first, term.first.second}]).transpose();
		}

		std::map<std::pair<int, int>, Mesh::Point> P_u_cross_P_v;
		for (int i = 0; i <= 4; ++i) {
			for (int j = 0; j <= 4 - i; ++j) {
				P_u_cross_P_v[{i, j}] = { 0,0,0 };
			}
		}
		for (const auto& term1 : P) {
			if (term1.first.first == 0)continue;
			for (const auto& term2 : P) {
				if (term2.first.second == 0)continue;
				P_u_cross_P_v[{term1.first.first + term2.first.first - 1, term1.first.second + term2.first.second - 1}] +=
					term1.first.first*term2.first.second*(term1.second%term2.second);
			}
		}

		for (const auto& term : P_u_cross_P_v) {
			neumann_part += convertToEigen(term.second)*convertToEigen(NeumannGradientParameter[{term.first.first, term.first.second}]).transpose();
		}
	}
	return (dirichlet_part + neumann_part) / 4 / M_PI;
}

PGC::Tensor PGC::ComputeHessianByFacet
(const int point_idx, const int facet_idx, const std::vector<Mesh::Point>& control_point) {
	Tensor dirichlet_part, neumann_part, result;
	dirichlet_part = {};
	neumann_part = {};

	auto& DirichletHessianParameter = DirichletHessianParametersForMeshPoints.at({ point_idx,facet_idx });
	auto& NeumannHessianParameter = NeumannHessianParametersForMeshPoints.at({ point_idx,facet_idx });
	if (cage_mode == false) {
		sum_point_product(dirichlet_part, control_point[0], DirichletHessianParameter[{0, 0}]);
		sum_point_product(dirichlet_part, -control_point[0] + control_point[1], DirichletHessianParameter[{1, 0}]);
		sum_point_product(dirichlet_part, -control_point[0] + control_point[3], DirichletHessianParameter[{0, 1}]);
		sum_point_product(dirichlet_part, control_point[0] - control_point[1] + control_point[2] - control_point[3],
			DirichletHessianParameter[{1, 1}]);

		sum_point_product(neumann_part, (-control_point[0] + control_point[1]) % (-control_point[0] + control_point[3]),
			NeumannHessianParameter[{0, 0}]);
		sum_point_product(neumann_part, (-control_point[0] + control_point[1]) % (-control_point[3] + control_point[2]),
			NeumannHessianParameter[{1, 0}]);
		sum_point_product(neumann_part, (-control_point[1] + control_point[2]) % (-control_point[0] + control_point[3]),
			NeumannHessianParameter[{0, 1}]);
	}
	else {
		auto& P003 = control_point[0];
		auto& P012 = control_point[1];
		auto& P021 = control_point[2];
		auto& P030 = control_point[3];
		auto& P102 = control_point[4];
		auto& P111 = control_point[5];
		auto& P120 = control_point[6];
		auto& P201 = control_point[7];
		auto& P210 = control_point[8];
		auto& P300 = control_point[9];

		std::map<std::pair<int, int>, Mesh::Point> P;
		P[{0, 0}] = P003;
		P[{1, 0}] = -3 * P003 + 3 * P012;
		P[{0, 1}] = -3 * P003 + 3 * P102;
		P[{2, 0}] = 3 * P003 - 6 * P012 + 3 * P021;
		P[{1, 1}] = 6 * P003 - 6 * P012 - 6 * P102 + 6 * P111;
		P[{0, 2}] = 3 * P003 - 6 * P102 + 3 * P201;
		P[{3, 0}] = -P003 + 3 * P012 - 3 * P021 + P030;
		P[{2, 1}] = -3 * P003 + 6 * P012 - 3 * P021 + 3 * P102 - 6 * P111 + 3 * P120;
		P[{1, 2}] = -3 * P003 + 3 * P012 + 6 * P102 - 6 * P111 - 3 * P201 + 3 * P210;
		P[{0, 3}] = -P003 + 3 * P102 - 3 * P201 + P300;

		for (const auto& term : P) {
			sum_point_product(dirichlet_part, term.second, DirichletHessianParameter[{term.first.first, term.first.second}]);
		}

		std::map<std::pair<int, int>, Mesh::Point> P_u_cross_P_v;
		for (int i = 0; i <= 4; ++i) {
			for (int j = 0; j <= 4 - i; ++j) {
				P_u_cross_P_v[{i, j}] = { 0,0,0 };
			}
		}
		for (const auto& term1 : P) {
			if (term1.first.first == 0)continue;
			for (const auto& term2 : P) {
				if (term2.first.second == 0)continue;
				P_u_cross_P_v[{term1.first.first + term2.first.first - 1, term1.first.second + term2.first.second - 1}] +=
					term1.first.first*term2.first.second*(term1.second%term2.second);
			}
		}

		for (const auto& term : P_u_cross_P_v) {
			sum_point_product(neumann_part, term.second, NeumannHessianParameter[{term.first.first, term.first.second}]);
		}
	}

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				result[i][j][k] = (dirichlet_part[i][j][k] + neumann_part[i][j][k]) / 4 / M_PI;
			}
		}
	}

	return result;
}

double PGC::DirichletSubIntegral
(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2, const Mesh::Point& v3, const Mesh::Point& point_to_p0) {
	double result = 0;
	if (i == 0 && j == 0) {
		result += DirichletRThetaSubIntegral(c0, alpha1, beta1 + alpha1, projection_to_p0.length()*std::sin(alpha1));
		result += DirichletRThetaSubIntegral(c0, alpha2, beta2 + alpha2, projection_to_p1.length()*std::sin(alpha2));
		result += DirichletRThetaSubIntegral(c0, alpha3, beta3 + alpha3, projection_to_p2.length()*std::sin(alpha3));
		result += DirichletRThetaSubIntegral(c0, alpha4, beta4 + alpha4, projection_to_p3.length()*std::sin(alpha4));
		return result * (b11*b22 - b12 * b21);
	}

	//std::cout << "\n" << "i:" << i << ",j:" << j << "\n";

	//存储方式为x^a*y^b*z^0对应u^a*v^b
	TernaryPolynomial f0, f1;
	for (int k = 0; k < i; ++k) {
		f0.addTerm(k, j, 0, std::pow(-s, i - k - 1));
	}
	for (int k = 0; k < j; ++k) {
		f1.addTerm(0, k, 0, std::pow(-t, j - k - 1));
	}
	//u^i*v^j=(s1*u+s2*v+s3)*P_fun+(t1*u+t2*v+t3)*Q_fun+(-s)^i*(-t)^j
	auto&P_fun = f0 * (1.0 / (s1 - t1 / t2 * s2)) + f1 * (std::pow(-s, i) / (s2 - t2 / t1 * s1));
	auto&Q_fun = f0 * (1.0 / (t1 - s1 / s2 * t2)) + f1 * (std::pow(-s, i) / (t2 - s2 / s1 * t1));
	//std::cout << "P_fun:";
	//P_fun.display();
	//std::cout << "Q_fun:";
	//Q_fun.display();

	//计算(s1*u+s2*v+s3)*P_fun部分的积分
	for (const auto& term : P_fun.terms) {
		if (term.first[0] > 0) {
			result += term.second*(-term.first[0])*
				NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{term.first[0] - 1, term.first[1]}];
		}
	}

	for (int edge_index = 0; edge_index < 3; ++edge_index) {
		//边为u=c0*v+d0
		if (std::abs(b - 1) < 1e-6 && edge_index == 1)continue;
		if (std::abs(b) < 1e-6 && edge_index == 0)continue;
		double c0 = c_u[edge_index];
		double d0 = d_u[edge_index];

		//积分分母为sqrt(v_1*v^2+v_2*v+v_3)
		double v_1 = ((c0*v1 + v2) | (c0*v1 + v2));
		double v_2 = 2 * ((c0*v1 + v2) | (d0*v1 + point_to_p0));
		double v_3 = ((d0*v1 + point_to_p0) | (d0*v1 + point_to_p0));

		//u=c0*v+d0代入P_fun后分子
		Polynomial v_edge_integral_0(0);
		for (const auto& term : P_fun.terms) {
			v_edge_integral_0 += term.second*polynomial_power(Polynomial{ d0, c0 }, term.first[0])*
				polynomial_power(Polynomial{ 0, 1 }, term.first[1]);
		}

		//换元为v=p*v_bar+q
		double p = std::sqrt(v_3 / v_1 - v_2 * v_2 / 4 / v_1 / v_1);
		double q = -v_2 / v_1 / 2;
		double lower_bound = (e_u[edge_index].first - q) / p;
		double upper_bound = (e_u[edge_index].second - q) / p;

		Polynomial v_bar_edge_integral(0);
		for (int k = 0; k < v_edge_integral_0.size(); ++k) {
			v_bar_edge_integral += v_edge_integral_0[k] * polynomial_power(Polynomial{ q, p }, k);
		}

		//std::cout << "v_bar_edge_integral:" << v_bar_edge_integral << "\n";

		for (int k = 0; k < v_bar_edge_integral.size(); ++k) {
			result += v_bar_edge_integral[k].real() / std::sqrt(v_1)*
				(SubsidiaryIntegral(k, upper_bound) - SubsidiaryIntegral(k, lower_bound));
		}
	}

	//计算(t1*u+t2*v+t3)*Q_fun部分的积分
	for (const auto& term : Q_fun.terms) {
		if (term.first[1] > 0) {
			result += term.second*(-term.first[1])*
				NeumannParametersForMeshPoints[{vertex_idx, facet_idx}][{term.first[0], term.first[1] - 1}];
		}
	}

	for (int edge_index = 0; edge_index < 3; ++edge_index) {
		//边为v=c0*u+d0
		if (std::abs(a) < 1e-6 && edge_index == 2)continue;
		if (std::abs(a - 1) < 1e-6 && edge_index == 1)continue;
		double c0 = c_v[edge_index];
		double d0 = d_v[edge_index];

		//积分分母为sqrt(u_1*u^2+u_2*u+u_3)
		double u_1 = ((v1 + c0 * v2) | (v1 + c0 * v2));
		double u_2 = 2 * ((v1 + c0 * v2) | (d0*v2 + point_to_p0));
		double u_3 = ((d0*v2 + point_to_p0) | (d0*v2 + point_to_p0));

		//v=c0*u+d0代入Q_fun后分子
		Polynomial u_edge_integral_0(0);
		for (const auto& term : Q_fun.terms) {
			u_edge_integral_0 += term.second*polynomial_power(Polynomial{ 0, 1 }, term.first[0])*
				polynomial_power(Polynomial{ d0, c0 }, term.first[1]);
		}

		//换元为u=p*u_bar+q
		double p = std::sqrt(u_3 / u_1 - u_2 * u_2 / 4 / u_1 / u_1);
		double q = -u_2 / u_1 / 2;
		double lower_bound = (e_v[edge_index].first - q) / p;
		double upper_bound = (e_v[edge_index].second - q) / p;

		Polynomial u_bar_edge_integral(0);
		for (int k = 0; k < u_edge_integral_0.size(); ++k) {
			u_bar_edge_integral += u_edge_integral_0[k] * polynomial_power(Polynomial{ q, p }, k);
		}

		//std::cout << "u_bar_edge_integral:" << u_bar_edge_integral << "\n";

		for (int k = 0; k < u_bar_edge_integral.size(); ++k) {
			result -= u_bar_edge_integral[k].real() / std::sqrt(u_1)*
				(SubsidiaryIntegral(k, upper_bound) - SubsidiaryIntegral(k, lower_bound));
		}
	}

	//计算(-s)^i*(-t)^j部分的积分
	result += std::pow(-s, i)*std::pow(-t, j)* DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{0, 0}] / projection_height;

	return result;
}

double PGC::DDirichletSubIntegral(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
	const Mesh::Point& v3, const Mesh::Point& point_to_p0) {
	double result = 0;
	if (i == 0 && j == 0) {
		result += DDirichletRThetaSubIntegral(c0, alpha1, beta1 + alpha1, projection_to_p0.length()*std::sin(alpha1));
		result += DDirichletRThetaSubIntegral(c0, alpha2, beta2 + alpha2, projection_to_p1.length()*std::sin(alpha2));
		result += DDirichletRThetaSubIntegral(c0, alpha3, beta3 + alpha3, projection_to_p2.length()*std::sin(alpha3));
		result += DDirichletRThetaSubIntegral(c0, alpha4, beta4 + alpha4, projection_to_p3.length()*std::sin(alpha4));
		return result * (b11*b22 - b12 * b21);
	}

	//std::cout << "\n" << "i:" << i << ",j:" << j << "\n";

	//d/du的分子为s1*u+s2*v+s3,d/dv的分子为t1*u+t2*v+t3,可分别线性表示为u+s,v+t
	double s1 = -3 * (v1 | v1);
	double s2 = -3 * (v1 | v2);
	double s3 = -3 * (v1 | point_to_p0);
	double t1 = -3 * (v2 | v1);
	double t2 = -3 * (v2 | v2);
	double t3 = -3 * (v2 | point_to_p0);

	//存储方式为x^a*y^b*z^0对应u^a*v^b
	TernaryPolynomial f0, f1;
	for (int k = 0; k < i; ++k) {
		f0.addTerm(k, j, 0, std::pow(-s, i - k - 1));
	}
	for (int k = 0; k < j; ++k) {
		f1.addTerm(0, k, 0, std::pow(-t, j - k - 1));
	}
	//u^i*v^j=(s1*u+s2*v+s3)*P_fun+(t1*u+t2*v+t3)*Q_fun+(-s)^i*(-t)^j
	auto&P_fun = f0 * (1.0 / (s1 - t1 / t2 * s2)) + f1 * (std::pow(-s, i) / (s2 - t2 / t1 * s1));
	auto&Q_fun = f0 * (1.0 / (t1 - s1 / s2 * t2)) + f1 * (std::pow(-s, i) / (t2 - s2 / s1 * t1));
	//std::cout << "P_fun:";
	//P_fun.display();
	//std::cout << "Q_fun:";
	//Q_fun.display();

	//计算(s1*u+s2*v+s3)*P_fun部分的积分
	for (const auto& term : P_fun.terms) {
		if (term.first[0] > 0) {
			result += term.second*(-term.first[0])*
				DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{term.first[0] - 1, term.first[1]}] / projection_height;
		}
	}

	for (int edge_index = 0; edge_index < 3; ++edge_index) {
		//边为u=c0*v+d0
		if (std::abs(b - 1) < 1e-6 && edge_index == 1)continue;
		if (std::abs(b) < 1e-6 && edge_index == 0)continue;
		double c0 = c_u[edge_index];
		double d0 = d_u[edge_index];

		//积分分母为sqrt(v_1*v^2+v_2*v+v_3)^3
		double v_1 = ((c0*v1 + v2) | (c0*v1 + v2));
		double v_2 = 2 * ((c0*v1 + v2) | (d0*v1 + point_to_p0));
		double v_3 = ((d0*v1 + point_to_p0) | (d0*v1 + point_to_p0));

		//u=c0*v+d0代入P_fun后分子
		Polynomial v_edge_integral_0(0);
		for (const auto& term : P_fun.terms) {
			v_edge_integral_0 += term.second*polynomial_power(Polynomial{ d0, c0 }, term.first[0])*
				polynomial_power(Polynomial{ 0, 1 }, term.first[1]);
		}

		//换元为v=p*v_bar+q
		double p = std::sqrt(v_3 / v_1 - v_2 * v_2 / 4 / v_1 / v_1);
		double q = -v_2 / v_1 / 2;
		double lower_bound = (e_u[edge_index].first - q) / p;
		double upper_bound = (e_u[edge_index].second - q) / p;

		Polynomial v_bar_edge_integral(0);
		for (int k = 0; k < v_edge_integral_0.size(); ++k) {
			v_bar_edge_integral += v_edge_integral_0[k] * polynomial_power(Polynomial{ q, p }, k);
		}

		//std::cout << "v_bar_edge_integral:" << v_bar_edge_integral << "\n";

		for (int k = 0; k < v_bar_edge_integral.size(); ++k) {
			result += v_bar_edge_integral[k].real() / p / p / v_1 / std::sqrt(v_1)*
				(DSubsidiaryIntegral(k, upper_bound) - DSubsidiaryIntegral(k, lower_bound));
		}
	}

	//计算(t1*u+t2*v+t3)*Q_fun部分的积分
	for (const auto& term : Q_fun.terms) {
		if (term.first[1] > 0) {
			result += term.second*(-term.first[1])*
				DirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{term.first[0], term.first[1] - 1}] / projection_height;
		}
	}

	for (int edge_index = 0; edge_index < 3; ++edge_index) {
		//边为v=c0*u+d0
		if (std::abs(a) < 1e-6 && edge_index == 2)continue;
		if (std::abs(a - 1) < 1e-6 && edge_index == 1)continue;
		double c0 = c_v[edge_index];
		double d0 = d_v[edge_index];

		//积分分母为sqrt(u_1*u^2+u_2*u+u_3)^3
		double u_1 = ((v1 + c0 * v2) | (v1 + c0 * v2));
		double u_2 = 2 * ((v1 + c0 * v2) | (d0*v2 + point_to_p0));
		double u_3 = ((d0*v2 + point_to_p0) | (d0*v2 + point_to_p0));

		//v=c0*u+d0代入Q_fun后分子
		Polynomial u_edge_integral_0(0);
		for (const auto& term : Q_fun.terms) {
			u_edge_integral_0 += term.second*polynomial_power(Polynomial{ 0, 1 }, term.first[0])*
				polynomial_power(Polynomial{ d0, c0 }, term.first[1]);
		}

		//换元为u=p*u_bar+q
		double p = std::sqrt(u_3 / u_1 - u_2 * u_2 / 4 / u_1 / u_1);
		double q = -u_2 / u_1 / 2;
		double lower_bound = (e_v[edge_index].first - q) / p;
		double upper_bound = (e_v[edge_index].second - q) / p;

		Polynomial u_bar_edge_integral(0);
		for (int k = 0; k < u_edge_integral_0.size(); ++k) {
			u_bar_edge_integral += u_edge_integral_0[k] * polynomial_power(Polynomial{ q, p }, k);
		}

		//std::cout << "u_bar_edge_integral:" << u_bar_edge_integral << "\n";

		for (int k = 0; k < u_bar_edge_integral.size(); ++k) {
			result -= u_bar_edge_integral[k].real() / p / p / u_1 / std::sqrt(u_1)*
				(DSubsidiaryIntegral(k, upper_bound) - DSubsidiaryIntegral(k, lower_bound));
		}
	}

	//计算(-s)^i*(-t)^j部分的积分
	result += std::pow(-s, i)*std::pow(-t, j)* DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{0, 0}];

	return result;
}

double PGC::DDDirichletSubIntegral(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
	const Mesh::Point& v3, const Mesh::Point& point_to_p0) {
	double result = 0;
	if (i == 0 && j == 0) {
		result += DDDirichletRThetaSubIntegral(c0, alpha1, beta1 + alpha1, projection_to_p0.length()*std::sin(alpha1));
		result += DDDirichletRThetaSubIntegral(c0, alpha2, beta2 + alpha2, projection_to_p1.length()*std::sin(alpha2));
		result += DDDirichletRThetaSubIntegral(c0, alpha3, beta3 + alpha3, projection_to_p2.length()*std::sin(alpha3));
		result += DDDirichletRThetaSubIntegral(c0, alpha4, beta4 + alpha4, projection_to_p3.length()*std::sin(alpha4));
		return result * (b11*b22 - b12 * b21);
	}

	//std::cout << "\n" << "i:" << i << ",j:" << j << "\n";

	//d/du的分子为s1*u+s2*v+s3,d/dv的分子为t1*u+t2*v+t3,可分别线性表示为u+s,v+t
	double s1 = -5 * (v1 | v1);
	double s2 = -5 * (v1 | v2);
	double s3 = -5 * (v1 | point_to_p0);
	double t1 = -5 * (v2 | v1);
	double t2 = -5 * (v2 | v2);
	double t3 = -5 * (v2 | point_to_p0);

	//存储方式为x^a*y^b*z^0对应u^a*v^b
	TernaryPolynomial f0, f1;
	for (int k = 0; k < i; ++k) {
		f0.addTerm(k, j, 0, std::pow(-s, i - k - 1));
	}
	for (int k = 0; k < j; ++k) {
		f1.addTerm(0, k, 0, std::pow(-t, j - k - 1));
	}
	//u^i*v^j=(s1*u+s2*v+s3)*P_fun+(t1*u+t2*v+t3)*Q_fun+(-s)^i*(-t)^j
	auto&P_fun = f0 * (1.0 / (s1 - t1 / t2 * s2)) + f1 * (std::pow(-s, i) / (s2 - t2 / t1 * s1));
	auto&Q_fun = f0 * (1.0 / (t1 - s1 / s2 * t2)) + f1 * (std::pow(-s, i) / (t2 - s2 / s1 * t1));
	//std::cout << "P_fun:";
	//P_fun.display();
	//std::cout << "Q_fun:";
	//Q_fun.display();

	//计算(s1*u+s2*v+s3)*P_fun部分的积分
	for (const auto& term : P_fun.terms) {
		if (term.first[0] > 0) {
			result += term.second*(-term.first[0])*
				DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{term.first[0] - 1, term.first[1]}];
		}
	}

	for (int edge_index = 0; edge_index < 3; ++edge_index) {
		//边为u=c0*v+d0
		if (std::abs(b - 1) < 1e-6 && edge_index == 1)continue;
		if (std::abs(b) < 1e-6 && edge_index == 0)continue;
		double c0 = c_u[edge_index];
		double d0 = d_u[edge_index];

		//积分分母为sqrt(v_1*v^2+v_2*v+v_3)^3
		double v_1 = ((c0*v1 + v2) | (c0*v1 + v2));
		double v_2 = 2 * ((c0*v1 + v2) | (d0*v1 + point_to_p0));
		double v_3 = ((d0*v1 + point_to_p0) | (d0*v1 + point_to_p0));

		//u=c0*v+d0代入P_fun后分子
		Polynomial v_edge_integral_0(0);
		for (const auto& term : P_fun.terms) {
			v_edge_integral_0 += term.second*polynomial_power(Polynomial{ d0, c0 }, term.first[0])*
				polynomial_power(Polynomial{ 0, 1 }, term.first[1]);
		}

		//换元为v=p*v_bar+q
		double p = std::sqrt(v_3 / v_1 - v_2 * v_2 / 4 / v_1 / v_1);
		double q = -v_2 / v_1 / 2;
		double lower_bound = (e_u[edge_index].first - q) / p;
		double upper_bound = (e_u[edge_index].second - q) / p;

		Polynomial v_bar_edge_integral(0);
		for (int k = 0; k < v_edge_integral_0.size(); ++k) {
			v_bar_edge_integral += v_edge_integral_0[k] * polynomial_power(Polynomial{ q, p }, k);
		}

		//std::cout << "v_bar_edge_integral:" << v_bar_edge_integral << "\n";

		for (int k = 0; k < v_bar_edge_integral.size(); ++k) {
			result += v_bar_edge_integral[k].real() / std::pow(p, 4) / v_1 / v_1 / std::sqrt(v_1)*
				(DDSubsidiaryIntegral(k, upper_bound) - DDSubsidiaryIntegral(k, lower_bound));
		}
	}

	//计算(t1*u+t2*v+t3)*Q_fun部分的积分
	for (const auto& term : Q_fun.terms) {
		if (term.first[1] > 0) {
			result += term.second*(-term.first[1])*
				DDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{term.first[0], term.first[1] - 1}];
		}
	}

	for (int edge_index = 0; edge_index < 3; ++edge_index) {
		//边为v=c0*u+d0
		if (std::abs(a) < 1e-6 && edge_index == 2)continue;
		if (std::abs(a - 1) < 1e-6 && edge_index == 1)continue;
		double c0 = c_v[edge_index];
		double d0 = d_v[edge_index];

		//积分分母为sqrt(u_1*u^2+u_2*u+u_3)^3
		double u_1 = ((v1 + c0 * v2) | (v1 + c0 * v2));
		double u_2 = 2 * ((v1 + c0 * v2) | (d0*v2 + point_to_p0));
		double u_3 = ((d0*v2 + point_to_p0) | (d0*v2 + point_to_p0));

		//v=c0*u+d0代入Q_fun后分子
		Polynomial u_edge_integral_0(0);
		for (const auto& term : Q_fun.terms) {
			u_edge_integral_0 += term.second*polynomial_power(Polynomial{ 0, 1 }, term.first[0])*
				polynomial_power(Polynomial{ d0, c0 }, term.first[1]);
		}

		//换元为u=p*u_bar+q
		double p = std::sqrt(u_3 / u_1 - u_2 * u_2 / 4 / u_1 / u_1);
		double q = -u_2 / u_1 / 2;
		double lower_bound = (e_v[edge_index].first - q) / p;
		double upper_bound = (e_v[edge_index].second - q) / p;

		Polynomial u_bar_edge_integral(0);
		for (int k = 0; k < u_edge_integral_0.size(); ++k) {
			u_bar_edge_integral += u_edge_integral_0[k] * polynomial_power(Polynomial{ q, p }, k);
		}

		//std::cout << "u_bar_edge_integral:" << u_bar_edge_integral << "\n";

		for (int k = 0; k < u_bar_edge_integral.size(); ++k) {
			result -= u_bar_edge_integral[k].real() / std::pow(p, 4) / u_1 / u_1 / std::sqrt(u_1)*
				(DDSubsidiaryIntegral(k, upper_bound) - DDSubsidiaryIntegral(k, lower_bound));
		}
	}

	//计算(-s)^i*(-t)^j部分的积分
	result += std::pow(-s, i)*std::pow(-t, j)* DDDirichletParametersForMeshPoints[{vertex_idx, facet_idx}][{0, 0}];

	return result;
}

double PGC::NeumannSubIntegral
(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2, const Mesh::Point& v3, const Mesh::Point& point_to_p0) {
	double result = 0;
	//u=b11*r*cos(theta)+b12*r*sin(theta)+b13
	//u=b21*r*cos(theta)+b22*r*sin(theta)+b23
	//存储方式为(r,sin(theta),cos(theta))
	TernaryPolynomial u, v, result_poly;
	u.addTerm(1, 0, 1, b11);
	u.addTerm(1, 1, 0, b12);
	u.addTerm(0, 0, 0, b13);
	v.addTerm(1, 0, 1, b21);
	v.addTerm(1, 1, 0, b22);
	v.addTerm(0, 0, 0, b23);
	result_poly.addTerm(1, 0, 0, b11*b22 - b12 * b21);
	result_poly = result_poly * u.power(i)*v.power(j);
	//std::cout << "u:";
	//u.display();
	//std::cout << "v:";
	//v.display();
	//std::cout << "u.power(i):";
	//u.power(i).display();
	//std::cout << "v.power(i):";
	//v.power(i).display();
	//std::cout << "result_poly:";
	//result_poly.display();

	for (const auto& term : result_poly.terms) {
		result += term.second*NeumannRThetaSubIntegral(term.first[1], term.first[2], term.first[0], c0,
			0, beta1, projection_to_p0.length()*std::sin(alpha1), alpha1);
		result += term.second* NeumannRThetaSubIntegral(term.first[1], term.first[2], term.first[0], c0,
			beta1, beta1 + beta2, projection_to_p1.length()*std::sin(alpha2), alpha2);
		result += term.second*NeumannRThetaSubIntegral(term.first[1], term.first[2], term.first[0], c0,
			beta1 + beta2, beta1 + beta2 + beta3, projection_to_p2.length()*std::sin(alpha3), alpha3);
		result += term.second*NeumannRThetaSubIntegral(term.first[1], term.first[2], term.first[0], c0,
			beta1 + beta2 + beta3, beta1 + beta2 + beta3 + beta4, projection_to_p3.length()*std::sin(alpha4), alpha4);
	}
	return result;
}

double PGC::ComputeAngle(const Mesh::Point& start, const Mesh::Point& end, const Mesh::Point& normal) {
	double start_length = start.length();
	double end_length = end.length();
	if (start_length < 1e-6 || end_length < 1e-6)return M_PI / 2;
	double cosTheta = (start | end) / start_length / end_length;
	cosTheta = std::fmax(-1.0, std::fmin(1.0, cosTheta));
	double angle = std::acos(cosTheta);
	auto&crossProd = (start%end);

	if ((crossProd | normal) < 0) {
		angle = -angle;
	}
	return angle;
}

//checked
double PGC::NeumannRThetaSubIntegral(const int m, const int n, const int l, const double c0,
	const double theta_lower_bound, const double theta_upper_bound, const double c1, const double alpha) {
	double result = 0;
	if (std::abs(c1) < 1e-6)return 0;
	if (std::abs(theta_upper_bound - theta_lower_bound) < 1e-6)return 0;
	if (l != m + n + 1) {
		std::cout << "l!=m+n+1,error!" << "\n";
	}
	auto&B_func = BIntegral(m, n);
	//std::cout << "B:";
	//B_func.display();
	double A0 = AIntegral(l, 0, c0);
	result += B_func.evalue(theta_upper_bound, std::sin(theta_upper_bound), std::cos(theta_upper_bound))*
		(AIntegral(l, c1 / std::sin(alpha + theta_upper_bound - theta_lower_bound), c0) - A0) -
		B_func.evalue(theta_lower_bound, std::sin(theta_lower_bound), std::cos(theta_lower_bound))*
		(AIntegral(l, c1 / std::sin(alpha), c0) - A0);

	TernaryPolynomial translated_B_func;
	for (const auto& term : B_func.terms) {
		TernaryPolynomial translated_theta;
		TernaryPolynomial translated_sin_theta;
		TernaryPolynomial translated_cos_theta;
		double translated_distance = alpha - theta_lower_bound;
		translated_theta.addTerm(0, 0, 0, -translated_distance);
		translated_theta.addTerm(1, 0, 0, 1);
		translated_sin_theta.addTerm(0, 1, 0, std::cos(translated_distance));
		translated_sin_theta.addTerm(0, 0, 1, -std::sin(translated_distance));
		translated_cos_theta.addTerm(0, 1, 0, std::sin(translated_distance));
		translated_cos_theta.addTerm(0, 0, 1, std::cos(translated_distance));
		translated_B_func = translated_B_func + translated_theta.power(term.first[0])*
			translated_sin_theta.power(term.first[1])*	translated_cos_theta.power(term.first[2])*term.second;
	}
	//std::cout << "trans_B:";
	//translated_B_func.display();

	double new_theta_lower_bound = alpha;
	double new_theta_upper_bound = theta_upper_bound - theta_lower_bound + alpha;
	double k = c0 / c1;

	double partial_result = 0;
	for (const auto& term : translated_B_func.terms) {
		double partial_result_of_term = 0;
		int a = term.first[0];
		int b = term.first[2] + 1;
		int c = l - term.first[1] + 1;

		//if (std::abs(term.second) < 1e-6)continue;

		//compute Integrate[theta^a*cos(theta)^b/sin(theta)^c/sqrt(1+k^2*sin(theta)^2)]
		if (term.first[0] > 1 || (term.first[0] == 1 && term.first[1] + term.first[2] != 0) ||
			(term.first[0] == 1 && l % 2 == 0) || ((l - term.first[1] - term.first[2]) % 2 == 0)) {
			std::cout << "NeumannRThetaSubIntegral error!!" << "\n";
			continue;
		}
		else if (term.first[0] == 1) {
			//use partial integral
			partial_result_of_term += new_theta_upper_bound * NeumannSinCosIntegral(b, c, new_theta_upper_bound, k) -
				new_theta_lower_bound * NeumannSinCosIntegral(b, c, new_theta_lower_bound, k);
			partial_result_of_term -= NeumannRThetaSubIntegralOfConstIntegral(c, new_theta_upper_bound, new_theta_lower_bound, k);
		}
		else {
			partial_result_of_term += NeumannSinCosIntegral(b, c, new_theta_upper_bound, k) -
				NeumannSinCosIntegral(b, c, new_theta_lower_bound, k);
		}
		partial_result += partial_result_of_term * term.second;
	}
	result += partial_result * std::pow(c1, l);
	return result;
}

//checked
//Integrate[r^l/sqrt(r*r+c0*c0)]
double PGC::AIntegral(const int l, const double r, const double c0) {
	return std::pow(c0, l)*SubsidiaryIntegral(l, r / c0);
}

//checked
//Integrate[Integrate[cos(theta)/sin(theta)^c/sqrt(1+k*k*sin(theta)^2)]],c为偶数
double PGC::NeumannRThetaSubIntegralOfConstIntegral(const int c, const double theta_upper_bound, const double theta_lower_bound, const double k) {
	if (c % 2 == 1) {
		std::cout << "NeumannRThetaSubIntegralOfConstIntegral error" << "\n";
		return 0;
	}
	double result = 0;
	auto&subsidiary_poly = std::pow(k, c - 1)*SubsidiaryIntegralAsPolynomial(-c);
	//std::cout << "subsidiary_poly:" << subsidiary_poly << "\n";
	for (int i = 0; i < subsidiary_poly.size(); ++i) {
		if (i % 2 == 0)continue;
		result += subsidiary_poly[i].real()*(1.0 / std::pow(k, i)*SubsidiaryIntegral2(i, k, theta_upper_bound) +
			1.0 / std::pow(k, i - 2)*SubsidiaryIntegral2(i - 2, k, theta_upper_bound));
		result -= subsidiary_poly[i].real()*(1.0 / std::pow(k, i)*SubsidiaryIntegral2(i, k, theta_lower_bound) +
			1.0 / std::pow(k, i - 2)*SubsidiaryIntegral2(i - 2, k, theta_lower_bound));
	}
	return result;
}

//checked
//Integrate[1/sqrt(1+k^2*sin(theta)^2)/sin(theta)^i]
double PGC::SubsidiaryIntegral2(const int i, const double k, const double theta) {
	if (i % 2 == 0 || i < -1) {
		std::cout << "SubsidiaryIntegral2 error!" << "\n";
		return 0;
	}
	if (i == -1) {
		double t = std::tan(theta / 2);
		return std::atan(k*(-1 + t * t) / std::sqrt(1 + (2 + 4 * k*k)*t*t + std::pow(t, 4))) / k;
	}
	double temp = std::sqrt(1 + k * k*std::sin(theta)*std::sin(theta));
	if (i == 1) {
		return -std::atanh(std::cos(theta) / temp);
	}
	if (i == 3) {
		return (k*k - 1) / 2 * std::atanh(std::cos(theta) / temp) -
			temp * std::cos(theta) / 2 / std::sin(theta) / std::sin(theta);
	}
	return 1.0 / (i - 1)*(-std::cos(theta) / std::pow(std::sin(theta), i - 1)*temp +
		(i - 2)*(1 - k * k)*SubsidiaryIntegral2(i - 2, k, theta) +
		k * k*(i - 3)*SubsidiaryIntegral2(i - 4, k, theta));
}

//checked
//Integrate[cos(theta)^b/sin(theta)^c/sqrt(1+k*k*sin(theta)^2)],c>=b+1且c-b为奇数，b>=0
double PGC::NeumannSinCosIntegral(const int b, const int c, const double theta, const double k) {
	if (c <= b || (c - b) % 2 == 0) {
		std::cout << "NeumannSinCosIntegral error!" << "\n";
		return 0;
	}
	if (b == 0) {
		return SubsidiaryIntegral2(c, k, theta);
	}
	if (b == 1) {
		return std::pow(k, c - 1)*SubsidiaryIntegral(-c, std::sin(theta) * k);
	}
	return NeumannSinCosIntegral(b - 2, c, theta, k) - NeumannSinCosIntegral(b - 2, c - 2, theta, k);
}

//存储方式为(theta,sin(theta),cos(theta))
PGC::TernaryPolynomial PGC::BIntegral(const int m, const int n) {
	TernaryPolynomial result;
	if (m == 0 && n == 0) {
		result.addTerm(1, 0, 0, 1);
		return result;
	}
	if (m == 1) {
		result.addTerm(0, 0, n + 1, -1.0 / (n + 1));
		return result;
	}
	if (n == 1) {
		result.addTerm(0, m + 1, 0, 1.0 / (m + 1));
		return result;
	}
	if (m == 0) {
		auto&temp = BIntegral(0, n - 2);
		for (const auto& term : temp.terms) {
			result.addTerm(term.first[0], term.first[1], term.first[2], term.second*(n - 1) / n);
		}
		result.addTerm(0, 1, n - 1, 1.0 / n);
		return result;
	}
	if (n == 0) {
		auto&temp = BIntegral(m - 2, 0);
		for (const auto& term : temp.terms) {
			result.addTerm(term.first[0], term.first[1], term.first[2], term.second*(m - 1) / m);
		}
		result.addTerm(0, m - 1, 1, -1.0 / m);
		return result;
	}
	auto&temp = BIntegral(m + 2, n - 2);
	for (const auto& term : temp.terms) {
		result.addTerm(term.first[0], term.first[1], term.first[2], term.second*(n - 1) / (m + 1));
	}
	result.addTerm(0, m + 1, n - 1, 1.0 / (m + 1));
	return result;
}

//checked
double PGC::DirichletRThetaSubIntegral(const double c0, const double theta_lower_bound,
	const double theta_upper_bound, const double c1) {
	double result = 0;
	if (std::abs(c1) < 1e-6)return 0;
	if (std::abs(theta_upper_bound - theta_lower_bound) < 1e-6)return 0;
	result += (theta_upper_bound - theta_lower_bound) / c0;

	/*double new_lower_bound = std::tan(theta_lower_bound / 2);
	double new_upper_bound = std::tan(theta_upper_bound / 2);

	result -= std::atan(c0*(-1 + new_upper_bound * new_upper_bound) / std::sqrt(4 * c0*c0*new_upper_bound*new_upper_bound +
		c1 * c1*(1 + 2 * new_upper_bound*new_upper_bound + std::pow(new_upper_bound, 4)))) / c0;
	result += std::atan(c0*(-1 + new_lower_bound * new_lower_bound) / std::sqrt(4 * c0*c0*new_lower_bound*new_lower_bound +
		c1 * c1*(1 + 2 * new_lower_bound*new_lower_bound + std::pow(new_lower_bound, 4)))) / c0;*/

	double temp_upper = std::sin(theta_upper_bound);
	double temp_lower = std::sin(theta_lower_bound);
	double sign_upper = 1;
	double sign_lower = 1;
	if (temp_upper < 0)sign_upper = -1;
	if (temp_lower < 0)sign_lower = -1;

	result -= sign_upper * std::atan2(std::sqrt(c1*c1 + c0 * c0*std::pow(temp_upper, 2)), c0 * std::cos(theta_upper_bound)) / c0;
	result += sign_lower * std::atan2(std::sqrt(c1*c1 + c0 * c0*std::pow(temp_lower, 2)), c0 * std::cos(theta_lower_bound)) / c0;

	return result;
}

//checked
double PGC::DDirichletRThetaSubIntegral(const double c0, const double theta_lower_bound,
	const double theta_upper_bound, const double c1) {
	double result = 0;
	if (std::abs(c1) < 1e-6)return 0;
	if (std::abs(theta_upper_bound - theta_lower_bound) < 1e-6)return 0;
	result += (theta_upper_bound - theta_lower_bound) / 3 / std::pow(c0, 3);

	double k = c0 / c1;
	double temp_upper = std::cos(theta_upper_bound) / std::sqrt(1 + k * k*std::pow(std::sin(theta_upper_bound), 2));
	double temp_lower = std::cos(theta_lower_bound) / std::sqrt(1 + k * k*std::pow(std::sin(theta_lower_bound), 2));
	result -= 1.0 / 3 / std::pow(c1, 3)*(temp_upper / k / k / (1 + k * k) - 1.0 / std::pow(k, 3)*std::atan(k*temp_upper));
	result += 1.0 / 3 / std::pow(c1, 3)*(temp_lower / k / k / (1 + k * k) - 1.0 / std::pow(k, 3)*std::atan(k*temp_lower));
	return result;
}

double PGC::DDDirichletRThetaSubIntegral(const double c0, const double theta_lower_bound,
	const double theta_upper_bound, const double c1) {
	double result = 0;
	if (std::abs(c1) < 1e-6)return 0;
	if (std::abs(theta_upper_bound - theta_lower_bound) < 1e-6)return 0;

	double c0_2 = c0 * c0;
	double c1_2 = c1 * c1;
	double temp_upper = std::sqrt(c1*c1 + c0 * c0*std::pow(std::sin(theta_upper_bound), 2));
	double temp_lower = std::sqrt(c1*c1 + c0 * c0*std::pow(std::sin(theta_lower_bound), 2));
	double sign_upper = 1;
	double sign_lower = 1;
	if (std::sin(theta_upper_bound) < 0)sign_upper = -1;
	if (std::sin(theta_lower_bound) < 0)sign_lower = -1;

	result += (3 * theta_upper_bound - 3 * std::atan2(temp_upper, c0*std::cos(theta_upper_bound))*sign_upper +
		c0 * c1_2*std::cos(theta_upper_bound)*sign_upper*
		(-3 * c0_2*c0_2 - 7 * c0_2* c1_2 - 3 * c1_2*c1_2 + c0_2 * (3 * c0_2 + 2 * c1_2) *std::cos(2 * theta_upper_bound)) /
		(c0_2 + c1_2) / (c0_2 + c1_2) / std::pow(temp_upper, 3)) / 15 / std::pow(c0, 5);
	result -= (3 * theta_lower_bound - 3 * std::atan2(temp_lower, c0*std::cos(theta_lower_bound))*sign_lower +
		c0 * c1_2*std::cos(theta_lower_bound)*sign_lower*
		(-3 * c0_2*c0_2 - 7 * c0_2* c1_2 - 3 * c1_2*c1_2 + c0_2 * (3 * c0_2 + 2 * c1_2) *std::cos(2 * theta_lower_bound)) /
		(c0_2 + c1_2) / (c0_2 + c1_2) / std::pow(temp_lower, 3)) / 15 / std::pow(c0, 5);

	return result;
}

//checked
//Integrate[u^m/sqrt(u^2+1)],m为非负数或负偶数
double PGC::SubsidiaryIntegral(const int m, const double u) {
	if (m == 0) {
		return std::asinh(u);
	}
	if (m == 1) {
		return std::sqrt(1 + u * u);
	}
	if (m >= 2) {
		return std::pow(u, m - 1)*std::sqrt(1 + u * u) / m - (m - 1)*SubsidiaryIntegral(m - 2, u) / m;
	}
	//if (m == -1) {
	//	return -std::atanh(std::sqrt(1 + u * u));
	//}
	if (m == -2) {
		return -std::sqrt(1 + u * u) / u;
	}
	return std::pow(u, m + 1)*std::sqrt(1 + u * u) / (m + 1) - (m + 2)*SubsidiaryIntegral(m + 2, u) / (m + 1);
}

//Integrate[u^m/sqrt(u^2+1)^3],m为非负数或负偶数
double PGC::DSubsidiaryIntegral(const int m, const double u) {
	if (m == 0) {
		return u / std::sqrt(1 + u * u);
	}
	if (m == 1) {
		return -1.0 / std::sqrt(1 + u * u);
	}
	return SubsidiaryIntegral(m - 2, u) - DSubsidiaryIntegral(m - 2, u);
}

//Integrate[u^m/sqrt(u^2+1)^5],m为非负数或负偶数
double PGC::DDSubsidiaryIntegral(const int m, const double u) {
	if (m == 0) {
		return u * (3 + 2 * u*u) / 3 / std::pow(std::sqrt(1 + u * u), 3);
	}
	if (m == 1) {
		return -1.0 / 3 / std::pow(std::sqrt(1 + u * u), 3);
	}
	return DSubsidiaryIntegral(m - 2, u) - DDSubsidiaryIntegral(m - 2, u);
}

//Integrate[u^m/sqrt(u^2+1)],其中m为负偶数，返回不定积分的多项式，存储方式为单项式x^i对应1/u^i*sqrt(u^2+1)
PGC::Polynomial PGC::SubsidiaryIntegralAsPolynomial(const int m) {
	if (m >= 0 || m % 2 == 1) {
		std::cout << "SubsidiaryIntegralAsPolynomial error!" << "\n";
		return Polynomial{ 0 };
	}
	if (m == -2) {
		return Polynomial{ 0,-1 };
	}
	auto&result_poly = SubsidiaryIntegralAsPolynomial(m + 2);
	return -((m + 2.0) / (m + 1.0))*result_poly + 1.0 / (m + 1)*polynomial_power(Polynomial{ 0.0,1.0 }, -m - 1);
}

PGC::Polynomial PGC::polynomial_power(const Polynomial& p, int n) {
	Polynomial result(1.0);
	Polynomial base = p;

	while (n > 0) {
		if (n % 2 == 1) {
			result = result * base;
		}
		base = base * base;
		n /= 2;
	}
	return result;
}

//v3=a*v1+b*v2,return (a,b)
std::pair<double, double> PGC::ComputeCordWithBasis(const Mesh::Point& v1, const Mesh::Point& v2, const Mesh::Point& v3) {
	auto&v1_cross_v2 = (v1%v2);
	auto&v3_cross_v2 = (v3%v2);
	double a = v3_cross_v2.length() / v1_cross_v2.length();
	if ((v3_cross_v2 | v1_cross_v2) < 0)a = -a;
	auto&bv2 = v3 - a * v1;
	double b = bv2.length() / v2.length();
	if ((bv2 | v2) < 0)b = -b;
	return { a,b };
}

Eigen::Vector3d PGC::convertToEigen(const OpenMesh::Vec3d& vec) {
	return Eigen::Vector3d(vec[0], vec[1], vec[2]);
}

void PGC::sum_point_product(Tensor& tensor, const OpenMesh::Vec3d& vec, const Eigen::Matrix3d& mat) {
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				tensor[i][j][k] += vec[i] * mat(j, k);
			}
		}
	}
}

void PGC::addTensors(Tensor& t1, const Tensor& t2) {
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			for (size_t k = 0; k < 3; ++k) {
				t1[i][j][k] += t2[i][j][k];
			}
		}
	}
}

void PGC::printTensor(const Tensor& tensor) {
	for (const auto& matrix : tensor) {
		for (const auto& row : matrix) {
			for (const auto& element : row) {
				std::cout << element << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
}

Eigen::Matrix3d PGC::optimalRotationMatrix(const Eigen::Matrix3d& A) {
	// 奇异值分解 (SVD)
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV();

	// 构造 D 矩阵
	Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
	if ((U * V.transpose()).determinant() < 0) {
		D(2, 2) = -1; // 如果行列式为负，调整最后一个奇异值
	}

	// 计算旋转矩阵 S
	Eigen::Matrix3d S = U * D * V.transpose();
	return S;
}

void PGC::processPoints() {
	// 1. 获取不重复的代表点
	struct PointHash {
		size_t operator()(const Mesh::Point& point) const {
			// 使用点的模值来生成哈希
			return std::hash<double>{}(point.norm());
		}
	};

	struct PointEqual {
		bool operator()(const Mesh::Point& a, const Mesh::Point& b) const {
			return (a - b).norm() < 1e-6; // 判断两点是否近似相等
		}
	};

	std::unordered_set<Mesh::Point, PointHash, PointEqual> unique_points;

	for (const auto& point : control_points) {
		unique_points.insert(point);
	}

	// 将去重后的点保存回原始容器
	represent_control_points.assign(unique_points.begin(), unique_points.end());

	// 2. 构建映射关系：point -> represent_point
	map_control_point_to_represent_control_point_idx.resize(control_points.size());
	for (size_t i = 0; i < control_points.size(); ++i) {
		// 找到当前点在represent_points中的位置
		auto it = std::find_if(represent_control_points.begin(), represent_control_points.end(), [&](const Mesh::Point& rp) {
			return (control_points[i] - rp).norm() < 1e-6;
		});
		map_control_point_to_represent_control_point_idx[i] = std::distance(represent_control_points.begin(), it);
	}

	// 3. 构建映射关系：represent_point -> 子点索引
	map_represent_control_point_idx_to_control_point_idxs.clear();
	for (int i = 0; i < control_points.size(); ++i) {
		int represent_idx = map_control_point_to_represent_control_point_idx[i];
		map_represent_control_point_idx_to_control_point_idxs[represent_idx].push_back(i);
	}
}

// 寻找满足条件的向量
Mesh::Point PGC::find_positive_inner_product_vector(const std::vector<Mesh::Point>& vectors) {
	Mesh::Point u(0, 0, 0);
	double max_c = 0;

	for (double ux = -1.0; ux <= 1.0; ux += 0.1) {
		for (double uy = -1.0; uy <= 1.0; uy += 0.1) {
			for (double uz = -1.0; uz <= 1.0; uz += 0.1) {
				Mesh::Point candidate_u = { ux, uy, uz };
				if (ux == 0 && uy == 0 && uz == 0)continue;
				candidate_u.normalize();
				double min_dot = 1;
				bool checked = true;
				for (const auto& v : vectors) {
					double dot = candidate_u.dot(v);
					if (dot < max_c) {
						checked = false;
						break;
					}
					if (min_dot > dot)
						min_dot = dot;
				}
				if (!checked)continue;
				u = candidate_u;
				max_c = min_dot;
			}
		}
	}
	/*for (const auto& v : vectors) {
		if (u.dot(v) < max_c) {
			std::cout << "error" << "\n";
		}
	}
	std::cout << "max_c:" << max_c << "\n";*/
	return u;
}

void PGC::saveBDCSVD(const Eigen::BDCSVD<Eigen::MatrixXd>& solver, const std::string& filename) {
	Eigen::VectorXd singularValues = solver.singularValues();
	Eigen::MatrixXd U = solver.matrixU();
	Eigen::MatrixXd V = solver.matrixV();

	// 设置阈值来处理小的奇异值
	double tolerance = solver.singularValues().coeff(0)*solver.threshold();

	// 计算有效秩 l_rank
	int l_rank = (singularValues.array() > tolerance).count();

	// 计算 V * Σ^+ * U^T
	Eigen::MatrixXd V_sigma_plus_UT;
	V_sigma_plus_UT.noalias() = V.leftCols(l_rank) * singularValues.head(l_rank).asDiagonal().inverse()
		* U.leftCols(l_rank).adjoint();

	// 保存矩阵 V * Σ^+ * U^T
	std::ofstream ofs(filename, std::ios::binary);
	int rows = V_sigma_plus_UT.rows();
	int cols = V_sigma_plus_UT.cols();
	ofs.write(reinterpret_cast<const char*>(&rows), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(&cols), sizeof(int));
	ofs.write(reinterpret_cast<const char*>(V_sigma_plus_UT.data()), rows * cols * sizeof(double));
	ofs.close();

	// 1. 保存奇异值
	std::ofstream ofs_0("singularValues.dat", std::ios::binary);
	int singularSize = singularValues.size();
	ofs_0.write(reinterpret_cast<const char*>(&singularSize), sizeof(int));
	ofs_0.write(reinterpret_cast<const char*>(singularValues.data()), singularSize * sizeof(double));

	// 2. 保存 U 矩阵
	std::ofstream ofs_1("U.dat", std::ios::binary);
	int rowsU = U.rows(), colsU = U.cols();
	ofs_1.write(reinterpret_cast<const char*>(&rowsU), sizeof(int));
	ofs_1.write(reinterpret_cast<const char*>(&colsU), sizeof(int));
	ofs_1.write(reinterpret_cast<const char*>(U.data()), rowsU * colsU * sizeof(double));

	// 3. 保存 V 矩阵
	std::ofstream ofs_2("V.dat", std::ios::binary);
	int rowsV = V.rows(), colsV = V.cols();
	ofs_2.write(reinterpret_cast<const char*>(&rowsV), sizeof(int));
	ofs_2.write(reinterpret_cast<const char*>(&colsV), sizeof(int));
	ofs_2.write(reinterpret_cast<const char*>(V.data()), rowsV * colsV * sizeof(double));

	ofs.close();
}

void PGC::loadBDCSVD() {
	// 读取 V * Σ^+ * U^T
	std::ifstream ifs("SVD.dat", std::ios::binary);

	int rows, cols;
	ifs.read(reinterpret_cast<char*>(&rows), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&cols), sizeof(int));

	SVD_mat.resize(rows, cols);
	ifs.read(reinterpret_cast<char*>(SVD_mat.data()), rows * cols * sizeof(double));
	ifs.close();
	std::cout << "rows=" << rows << " cols=" << cols << " read SVD.dat complete\n";
}

void PGC::loadglobalsolver() {
	std::ifstream mat_in("global_solver_mat.dat", std::ios::binary);
	std::ifstream vec_in("global_solver_vec.txt");

	int rows, cols;
	mat_in.read(reinterpret_cast<char*>(&rows), sizeof(int));
	mat_in.read(reinterpret_cast<char*>(&cols), sizeof(int));

	global_solver_mat.resize(rows, cols);
	mat_in.read(reinterpret_cast<char*>(global_solver_mat.data()), rows * cols * sizeof(double));
	mat_in.close();

	double number;
	std::vector<double> numbers;
	while (vec_in >> number) {
		numbers.push_back(number);
	}
	global_solver_vec = Eigen::VectorXd::Map(numbers.data(), numbers.size());

	vec_in.close();
	std::cout << "load global solver complete\n";
}

int PGC::factorial(int n) {
	if (n <= 1) return 1;
	return n * factorial(n - 1);
}

double PGC::integrate_triangle(const std::function<double(double, double)>& f) {
	// 使用15个节点的高斯-勒让德积分规则
	constexpr int nodes = 150;
	boost::math::quadrature::gauss<double, nodes> integrator;

	// 外层积分：对x积分，内层对y积分
	auto inner_integral = [&](double x) {
		return integrator.integrate(
			[x, &f](double y) { return f(x, y); },  // 内层被积函数（关于y）
			0.0, 1.0 - x                                // y的积分区间[0,1]
		);
	};

	// 计算外层积分（关于x）
	double result = integrator.integrate(
		inner_integral,  // 外层被积函数（关于x）
		0.0, 1.0         // x的积分区间[0,1]
	);

	return result;
}

double Neumann_integrate_func(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
	const Mesh::Point& point_to_p0, const double u, const double v) {
	auto distance = (point_to_p0 + v1 * u + v2 * v).length();
	return std::pow(u, i)*std::pow(v, j) / distance;
}

double Dirichlet_integrate_func(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
	const Mesh::Point& point_to_p0, const double u, const double v) {
	auto distance = (point_to_p0 + v1 * u + v2 * v).length();
	return std::pow(u, i)*std::pow(v, j) / std::pow(distance, 3);
}