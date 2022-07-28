
#ifndef ELASTIC_PROBLEM_H_
#define ELASTIC_PROBLEM_H_

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/sparse_direct.h>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <filesystem>
#include <deal.II/lac/lapack_full_matrix.h>
#include "read_file.h"

#define DOMAIN_DIM 2
#define OUTPUT_DIM 3
#define pi 3.1415926535897932384626433832795028841971

namespace DEAL{
using namespace dealii;



class ElasticProblem
{
public:
	ElasticProblem();
	void run();

private:

	Triangulation<DOMAIN_DIM> triangulation;
	FESystem<DOMAIN_DIM>          fe;
	DoFHandler<DOMAIN_DIM>    dof_handler;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	SparseMatrix<double> constraint_matrix;

	AffineConstraints<double> constraints;

	read_file rf;
	some_data dat;

	int quadegadd = 10;
	Vector<double> solution;
	Vector<double> displacement_vector;


	int refinelevel = 2;
	double L = 1.0;


	void make_grid();
	void setup_system();
	void setup_constraints();
	void extract_params();
	void initialize_configurations();
	void output_data(Vector<double> &);
	void calculate_displacement_vector(Vector<double> &);

	std::vector<bool> boolean_vector(int,int);




};
}



#endif /* ELASTIC_PROBLEM_H_ */
