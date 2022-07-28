
#ifndef ELASTIC_PROBLEM_CC_
#define ELASTIC_PROBLEM_CC_

#include "../include/elastic_problem.h"

#include <deal.II/base/iterator_range.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>
#include <deal.II/dofs/dof_accessor.templates.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_q1.h>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>  // std::setprecision()
#include<bits/stdc++.h>




namespace DEAL
{
using namespace dealii;



ElasticProblem::ElasticProblem()
: fe(FESystem<DOMAIN_DIM>(FE_Q<DOMAIN_DIM>(2), OUTPUT_DIM),1,FESystem<DOMAIN_DIM>(FE_Q<DOMAIN_DIM>(2), DOMAIN_DIM),OUTPUT_DIM,FESystem<DOMAIN_DIM>(FE_DGQ<DOMAIN_DIM>(0), DOMAIN_DIM),OUTPUT_DIM)
, dof_handler(triangulation){}



void ElasticProblem::extract_params(){
	std::string str = "../paramfiles/params.dat";
	char *cstr = new char[str.length() + 1];
	strcpy(cstr, str.c_str());
	rf.readInputFile(cstr, dat);

	refinelevel = dat.refinelevel; std::cout << "	Refine Level: " << refinelevel << std::endl;
	L = dat.L; std::cout << "	L: " << L << std::endl;



}


void ElasticProblem::make_grid()
{
	GridGenerator::hyper_cube(triangulation, 0.0, L);
	triangulation.refine_global(refinelevel);

	std::cout << "   Number of active cells: " << triangulation.n_active_cells()
																																																																																																																																																																																																																																															<< std::endl << "   Total number of cells: "
																																																																																																																																																																																																																																															<< triangulation.n_cells() << std::endl;
}

void ElasticProblem::setup_system(){
	dof_handler.distribute_dofs(fe);
	solution.reinit(dof_handler.n_dofs());
	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()<< std::endl;


	QGauss<DOMAIN_DIM> quadrature_formula(fe.degree + quadegadd);
	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler,
			dsp,
			constraints,
			/*keep_constrained_dofs = */ false);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);
	//std::ofstream out("sparsity_pattern.svg");
	//sparsity_pattern.print_svg(out);
}

std::vector<bool> ElasticProblem::boolean_vector(int n,int N){
	std::vector<bool> bvec(N);
	for (int i = 0; i < N; i++){
		bvec[i] = false;
	}
	bvec[n] = true;
	return bvec;
}

void ElasticProblem::setup_constraints(){
	constraints.clear();

	DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	int number_dofs = dof_handler.n_dofs();


	std::vector<bool> u1_components = boolean_vector(0,15);
	ComponentMask u1_mask(u1_components);

	std::vector<bool> u2_components = boolean_vector(1,15);;
	ComponentMask u2_mask(u2_components);

	std::vector<bool> u3_components = boolean_vector(2,15);;
	ComponentMask u3_mask(u3_components);


	std::vector<bool> is_u1_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u1_mask, is_u1_comp);

	std::vector<bool> is_u2_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u2_mask, is_u2_comp);

	std::vector<bool> is_u3_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u3_mask, is_u3_comp);



	std::vector<Point<DOMAIN_DIM>> support_points(number_dofs);
	MappingQ1<DOMAIN_DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);
	/*
	unsigned int closestzeropointind = 0;
	double tempdist = sqrt(pow(support_points[0][0],2.0) + pow(support_points[0][1],2.0));

	unsigned int furthestrightmidpointind = 0;
	double rightdist = sqrt(pow(support_points[0][0] - L/2.0,2.0) + pow(support_points[0][1],2.0));

	for (unsigned int i = 0; i < support_points.size(); i++) {
		double idist = sqrt(pow(support_points[i][0],2.0) + pow(support_points[i][1],2.0));
		if (idist < tempdist){
			tempdist = idist;
			closestzeropointind = i;
		}

		double temprdist = sqrt(pow(support_points[i][0] - L/2.0,2.0) + pow(support_points[i][1],2.0));

		if (temprdist < rightdist) {
			rightdist = temprdist;
			furthestrightmidpointind = i;
		}
	}

	for (unsigned int i = 0; i < number_dofs; i++) {

		if (fabs(support_points[i][0] - support_points[closestzeropointind][0]) < 1.0e-6
				&& fabs(support_points[i][1]- support_points[closestzeropointind][1]) < 1.0e-6
				&& (is_u1_comp[i] || is_u2_comp[i] || is_u3_comp[i] || is_w1_comp[i] || is_w2_comp[i])) {
			constraints.add_line(i);
		}

		if (fabs(support_points[i][0] - support_points[furthestrightmidpointind][0]) < 1.0e-6
				&& fabs(support_points[i][1]- support_points[furthestrightmidpointind][1]) < 1.0e-6
				&& is_u2_comp[i]) {
			constraints.add_line(i);
		}



	}
	 */

	constraints.close();
}

void ElasticProblem::initialize_configurations(){

	// Initializing the plate to be flat

	DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	int number_dofs = dof_handler.n_dofs();


	std::vector<bool> u1_components = boolean_vector(0,15);
	ComponentMask u1_mask(u1_components);

	std::vector<bool> u2_components = boolean_vector(1,15);;
	ComponentMask u2_mask(u2_components);

	std::vector<bool> u3_components = boolean_vector(2,15);;
	ComponentMask u3_mask(u3_components);


	std::vector<bool> is_u1_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u1_mask, is_u1_comp);

	std::vector<bool> is_u2_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u2_mask, is_u2_comp);

	std::vector<bool> is_u3_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u3_mask, is_u3_comp);



	std::vector<Point<DOMAIN_DIM>> support_points(number_dofs);
	MappingQ1<DOMAIN_DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	for (int i = 0; i < number_dofs; i++){

		if (is_u1_comp[i]) {
			solution[i] = support_points[i][0];
		} else if (is_u2_comp[i]) {
			solution[i] = support_points[i][1];
		}
	}


}
void ElasticProblem::calculate_displacement_vector(Vector<double> & vin){
	displacement_vector.reinit(vin.size());

	DoFTools::make_hanging_node_constraints(dof_handler, constraints);
		int number_dofs = dof_handler.n_dofs();


		std::vector<bool> u1_components = boolean_vector(0,15);
		ComponentMask u1_mask(u1_components);

		std::vector<bool> u2_components = boolean_vector(1,15);;
		ComponentMask u2_mask(u2_components);

		std::vector<bool> u3_components = boolean_vector(2,15);;
		ComponentMask u3_mask(u3_components);


		std::vector<bool> is_u1_comp(number_dofs, false);
		DoFTools::extract_dofs(dof_handler, u1_mask, is_u1_comp);

		std::vector<bool> is_u2_comp(number_dofs, false);
		DoFTools::extract_dofs(dof_handler, u2_mask, is_u2_comp);

		std::vector<bool> is_u3_comp(number_dofs, false);
		DoFTools::extract_dofs(dof_handler, u3_mask, is_u3_comp);



		std::vector<Point<DOMAIN_DIM>> support_points(number_dofs);
		MappingQ1<DOMAIN_DIM> mapping;
		DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

		for (int i = 0; i < number_dofs; i++){

			if (is_u1_comp[i]) {
				displacement_vector[i] = solution[i] - support_points[i][0];
			} else if (is_u2_comp[i]) {
				displacement_vector[i] = solution[i] - support_points[i][1];
			} else {
				displacement_vector[i] = solution[i];
			}
		}

}

void ElasticProblem::output_data(Vector<double> & out){
	std::vector<std::string> solution_names;


	solution_names.emplace_back("u1");
	solution_names.emplace_back("u2");
	solution_names.emplace_back("u3");

	solution_names.emplace_back("xi1");
	solution_names.emplace_back("xi1");
	solution_names.emplace_back("xi1");

	solution_names.emplace_back("xi2");
	solution_names.emplace_back("xi2");
	solution_names.emplace_back("xi2");

	solution_names.emplace_back("xi3");
	solution_names.emplace_back("xi3");
	solution_names.emplace_back("xi3");

	solution_names.emplace_back("lambda1");
	solution_names.emplace_back("lambda1");
	solution_names.emplace_back("lambda1");

	solution_names.emplace_back("lambda2");
	solution_names.emplace_back("lambda2");
	solution_names.emplace_back("lambda2");

	solution_names.emplace_back("lambda3");
	solution_names.emplace_back("lambda3");
	solution_names.emplace_back("lambda3");


	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	data_component_interpretation(
			OUTPUT_DIM, DataComponentInterpretation::component_is_part_of_vector);
	for (int i = 0; i < 12; i++){
		data_component_interpretation.emplace_back(DataComponentInterpretation::component_is_part_of_vector);
	}

	/*
	DataOut<DOMAIN_DIM> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution,
			solution_names,
			DataOut<DOMAIN_DIM>::type_dof_data,
			data_component_interpretation);
	data_out.build_patches();

	std::ofstream output(
			"solution.vtk");
	std::cout << "gothere" << std::endl;
	data_out.write_vtk(output);
	 */

	DataOut<DOMAIN_DIM> data_out;
	data_out.attach_dof_handler(dof_handler);
	calculate_displacement_vector(solution);
	data_out.add_data_vector(displacement_vector, solution_names);
	data_out.build_patches();

	std::ofstream output("solutions.vtk");
	data_out.write_vtk(output);
}

void ElasticProblem::run()
{


	extract_params();
	make_grid();
	setup_system();
	setup_constraints();
	initialize_configurations();
	output_data(solution);

}


}

#endif // ELASTIC_PROBLEM_CC_
