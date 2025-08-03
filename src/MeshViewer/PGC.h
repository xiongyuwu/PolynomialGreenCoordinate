#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

#include <boost/math/tools/polynomial.hpp>

#include "QGLViewerWidget.h"
#include "BezierSurface.h"

class PGC : public QGLViewerWidget
{
public:
	PGC(const Mesh& mesh_);
	PGC() {
		cage_mode = false;
		control_point_size_per_facet = control_point_per_edge * control_point_per_edge;
		updateControlPoints();
	};

	void setmesh(const Mesh& mesh_) {
		mesh = mesh_;
	}
	Mesh& getmesh() {
		return mesh;
	}
	void set_change_cage_mode()
	{
		cage_mode = cage_mode ? false : true;
		std::cout << "curr cage mode:";
		if (cage_mode) {
			control_point_size_per_facet = control_point_per_edge * (control_point_per_edge + 1) / 2;
			std::cout << "tri control point num per edge:" << control_point_per_edge << "\n";
		}
		else {
			control_point_size_per_facet = control_point_per_edge * control_point_per_edge;
			std::cout << "quad control point num per edge:" << control_point_per_edge << "\n";
		}
		updateControlPoints();
		updateGL();
	}
	void draw_represent_control_points();
	void draw_cage();
	void draw_position_constraints();
	void updateControlPoints();
	void do_read_cage();
	void do_generate_para_complete_flow();
	Mesh::Point find_nearest(const std::vector<Mesh::Point>& source, const Mesh::Point& target) {
		double min_distance = 1e5;
		Mesh::Point pos;
		for (auto& point : source) {
			double distance = (point - target).sqrnorm();
			if (distance < min_distance) {
				min_distance = distance;
				pos = point;
			}
		}
		return pos;
	}

	std::vector<Mesh::Point>& get_represent_control_points() {
		return represent_control_points;
	}
	void set_represent_control_points(const int idx, const Mesh::Point& new_represent_control_point) {
		represent_control_points[idx] = new_represent_control_point;
		for (auto& i : map_represent_control_point_idx_to_control_point_idxs[idx]) {
			control_points[i] = new_represent_control_point;
		}
	}
	void set_represent_control_points(const std::vector<Mesh::Point>& new_represent_control_points) {
		represent_control_points = new_represent_control_points;
		for (int i = 0; i < control_points.size(); ++i) {
			control_points[i] = represent_control_points[map_control_point_to_represent_control_point_idx[i]];
		}
	}

	std::vector<Mesh::Point> get_position_constraints() {
		std::vector<Mesh::Point> points;
		for (auto& position_constraint : position_constraints) {
			points.push_back(mesh.point(mesh.vertex_handle(position_constraint.first)));
		}
		return points;
	}
	void set_position_constraints(const int idx, const Mesh::Point& new_pos) {
		position_constraints[idx].second = new_pos;
	}

	void do_deformation();
	void do_optimize_control_points();
	void do_variational_deformation();
	//compute level:0 for normal, 1 for gradient, 2 for hessian
	void compute_integral_para_for_sample_points(const std::string& prefix,
		const std::array<bool, 3>& compute_level, const std::vector<Mesh::Point>&sample_points);
	void do_read_integral_para();
	bool control_point_is_initialed = false;
	bool origin_cage_is_initialed = false;
	void generate_smooth_sample_points();
	void read_position_constraints();
	void read_orientation_constraints();
	void read_rigidity_constraints();
	void read_smoothness_constraints();

	//compute level:0 for normal, 1 for gradient, 2 for hessian
	void do_compute_para_for_control_points(const std::string& prefix,
		const std::array<bool, 3>& compute_level, const std::vector<Mesh::Point>&sample_points);

	const Mesh::Point ERROR_POINT = { DBL_MAX, DBL_MAX, DBL_MAX };
private:
	//0 for quad cage,1 for tri cage
	bool cage_mode;
	Mesh mesh;
	Mesh cage_mesh;
	double cage_mesh_diag_length = 0;

	class TernaryPolynomial {
	public:
		std::map<std::vector<int>, double> terms;

		void addTerm(int x, int y, int z, double coeff) {
			terms[{x, y, z}] += coeff;
		}

		TernaryPolynomial operator+(const TernaryPolynomial& other) const {
			TernaryPolynomial result;
			result.terms = terms;
			for (const auto& term : other.terms) {
				result.addTerm(term.first[0], term.first[1], term.first[2], term.second);
			}
			return result;
		}

		TernaryPolynomial operator-(const TernaryPolynomial& other) const {
			TernaryPolynomial result;
			result.terms = terms;
			for (const auto& term : other.terms) {
				result.addTerm(term.first[0], term.first[1], term.first[2], -term.second);
			}
			return result;
		}

		TernaryPolynomial operator*(const double k) const {
			TernaryPolynomial result;
			for (const auto& term : terms) {
				result.addTerm(term.first[0], term.first[1], term.first[2], term.second * k);
			}
			return result;
		}

		TernaryPolynomial operator*(const TernaryPolynomial& other) const {
			TernaryPolynomial result;
			for (const auto& term1 : terms) {
				for (const auto& term2 : other.terms) {
					result.addTerm(term1.first[0] + term2.first[0], term1.first[1] + term2.first[1]
						, term1.first[2] + term2.first[2], term1.second * term2.second);
				}
			}
			return result;
		}

		TernaryPolynomial power(int n) const {
			TernaryPolynomial result;
			result.addTerm(0, 0, 0, 1);
			TernaryPolynomial base = *this;

			while (n > 0) {
				if (n % 2 == 1) {
					result = result * base;
				}
				base = base * base;
				n /= 2;
			}
			return result;
		}

		double evalue(double x, double y, double z)const {
			double result = 0;
			for (const auto& term : terms) {
				result += term.second*std::pow(x, term.first[0])*std::pow(y, term.first[1])*std::pow(z, term.first[2]);
			}
			return result;
		}

		void display() const {
			for (const auto& term : terms) {
				if (term.second != 0) {
					std::cout << term.second << "x^" << term.first[0] << "y^" << term.first[1] << "z^" << term.first[2] << " + ";
				}
			}
			std::cout << "\b\b \n";
		}
	};

	typedef boost::math::tools::polynomial<std::complex<double>> Polynomial;
	typedef std::array<std::array<std::array<double, 3>, 3>, 3> Tensor;

	int control_point_per_edge = 4;
	int control_point_size_per_facet;
	std::vector<Mesh::Point> control_points;
	std::vector<std::pair<std::vector<Mesh::Point>, std::vector<Mesh::Point>>> facet_and_control_points_pairs;
	std::vector<Mesh::Point> origin_cage_vertices;
	std::vector<std::vector<size_t>> map_facet_to_point_idxs;
	std::vector<Mesh::Point> represent_control_points;
	std::vector<int> map_control_point_to_represent_control_point_idx;
	std::map<int, std::vector<int>> map_represent_control_point_idx_to_control_point_idxs;

	void processPoints();

	//constraints
	std::vector<std::pair<int, Mesh::Point>> position_constraints;
	std::vector<std::pair<int, Eigen::Matrix3d>> orientation_constraints;
	std::vector<Mesh::Point> rigidity_constraints;
	std::vector<Mesh::Point> smoothness_constraints;
	std::vector<int> fixed_vertices;
	double position_weight = 10000 / std::sqrt(3);
	double orientation_weight = 10000.0 / std::sqrt(9);
	double rigidity_weight = 1.0 / std::sqrt(9);
	double smoothness_weight = 1.0 / std::sqrt(27) * 2;
	double fixed_vertices_weight = 10000 / std::sqrt(3);

	bool position_constraints_is_set=false;
	bool orientation_constraints_is_set = false;
	bool rigidity_constraints_is_set = false;
	bool smoothness_constraints_is_set = false;
	bool constraint_is_set = false;

	Mesh::Point DeformedPositionByFacet(const int point_idx, const int facet_idx, const std::vector<Mesh::Point>& control_point);
	Eigen::Matrix3d ComputeGradientByFacet(const int point_idx, const int facet_idx, const std::vector<Mesh::Point>& control_point);
	Tensor ComputeHessianByFacet(const int point_idx, const int facet_idx, const std::vector<Mesh::Point>& control_point);

	Eigen::Matrix3d optimalRotationMatrix(const Eigen::Matrix3d& A);

	double DirichletSubIntegral(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
		const Mesh::Point& v3, const Mesh::Point& point_to_p0);
	double DDirichletSubIntegral(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
		const Mesh::Point& v3, const Mesh::Point& point_to_p0);
	double DDDirichletSubIntegral(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
		const Mesh::Point& v3, const Mesh::Point& point_to_p0);
	double DirichletRThetaSubIntegral(const double c0, const double theta_lower_bound, const double theta_upper_bound, const double c1);
	double DDirichletRThetaSubIntegral(const double c0, const double theta_lower_bound, const double theta_upper_bound, const double c1);
	double DDDirichletRThetaSubIntegral(const double c0, const double theta_lower_bound, const double theta_upper_bound, const double c1);

	double NeumannSubIntegral(const int i, const int j, const Mesh::Point& v1,
		const Mesh::Point& v2, const Mesh::Point& v3, const Mesh::Point& point_to_p0);
	double NeumannRThetaSubIntegral(const int m, const int n, const int l, const double c0, const double theta_lower_bound,
		const double theta_upper_bound, const double c1, const double alpha);
	double AIntegral(const int l, const double r, const double c0);
	double NeumannRThetaSubIntegralOfConstIntegral(const int c, const double theta_upper_bound, const double theta_lower_bound, const double k);
	double NeumannSinCosIntegral(const int b, const int c, const double theta, const double k);

	TernaryPolynomial BIntegral(const int m, const int n);

	double SubsidiaryIntegral(const int m, const double u);
	double DSubsidiaryIntegral(const int m, const double u);
	double DDSubsidiaryIntegral(const int m, const double u);
	Polynomial SubsidiaryIntegralAsPolynomial(const int m);
	double SubsidiaryIntegral2(const int i, const double k, const double theta);

	std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> DirichletParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> NeumannParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> DDirichletParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> DDDirichletParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, Mesh::Point>> DirichletGradientParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, Mesh::Point>> NeumannGradientParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, Eigen::Matrix3d>> DirichletHessianParametersForMeshPoints;
	std::map<std::pair<int, int>, std::map<std::pair<int, int>, Eigen::Matrix3d>> NeumannHessianParametersForMeshPoints;

	std::vector<double> map_vertex_idx_and_control_point_to_deformation_para;
	std::vector<double> map_vertex_idx_and_control_point_pair_to_deformation_para;
	std::vector<Mesh::Point> map_vertex_idx_and_control_point_to_gradient_deformation_para;
	std::vector<Mesh::Point> map_vertex_idx_and_control_point_pair_to_gradient_deformation_para;
	std::vector<Eigen::Matrix3d> map_vertex_idx_and_control_point_to_hessian_deformation_para;
	std::vector<Eigen::Matrix3d> map_vertex_idx_and_control_point_pair_to_hessian_deformation_para;

	Polynomial polynomial_power(const Polynomial& p, int n);
	double ComputeAngle(const Mesh::Point& start, const Mesh::Point& end, const Mesh::Point& normal);
	std::pair<double, double> ComputeCordWithBasis(const Mesh::Point& v1, const Mesh::Point& v2, const Mesh::Point& v3);
	Eigen::Vector3d convertToEigen(const OpenMesh::Vec3d& vec);
	void sum_point_product(Tensor& tensor, const OpenMesh::Vec3d& vec, const Eigen::Matrix3d& mat);
	void addTensors(Tensor& t1, const Tensor& t2);
	void printTensor(const Tensor& tensor);
	int factorial(int n);

	//计算同点、同面、不同i,j的对应积分时使用的公共参数
	int vertex_idx, facet_idx;
	double projection_height;
	Mesh::Point v1_cross_v2, normal, projection_to_p0, projection_to_p1, projection_to_p2, projection_to_p3, b1, b2, b3;
	std::pair<double, double> cord_1, cord_2, cord_3;

	//(r,theta)坐标下p0到点(r,theta)的向量为b1*r*cos(theta)+b2*r*sin(theta)+b3
	double b11, b21, b12, b22, b13, b23;
	double beta1, beta2, beta3, beta4, alpha1, alpha2, alpha3, alpha4, c0;
	double s1, s2, s3, t1, t2, t3, s, t;

	std::pair<double, double> cordinate;
	double a, b;

	std::vector<double> c_u, d_u, c_v, d_v;
	std::vector<std::pair<double, double>> e_u, e_v;

	Mesh::Point find_positive_inner_product_vector(const std::vector<Mesh::Point>& vectors);
	Eigen::MatrixXd SVD_mat, global_solver_mat;
	Eigen::VectorXd global_solver_vec;
	void saveBDCSVD(const Eigen::BDCSVD<Eigen::MatrixXd>& solver, const std::string& filename);
	void loadBDCSVD();
	void loadglobalsolver();

	bool variational_is_initialed = false;
	std::vector<int> map_var_to_represent_var;
	std::vector<int> map_var_to_represent_var_sign;
	int var_size, equation_size;
	int not_appeared_flag;

	std::vector<int> represent_control_pair_is_appeared;
	double energy;
	double energy_diff;
	int iter;
	Eigen::VectorXd global_solver_rotation_vec;
	Eigen::VectorXd x;

	int section_1, section_2, section_3, section_4, rigidity_start;

	//数值积分
	double integrate_triangle(const std::function<double(double, double)>& f);
};

double Neumann_integrate_func(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
	const Mesh::Point& point_to_p0, const double u, const double v);
double Dirichlet_integrate_func(const int i, const int j, const Mesh::Point& v1, const Mesh::Point& v2,
	const Mesh::Point& point_to_p0, const double u, const double v);