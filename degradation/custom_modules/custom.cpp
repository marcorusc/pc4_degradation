/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

#include <iostream>
#include <vector>
#include <omp.h>
#include <functional>

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	for (auto* pCD: *getElongatedCellDefinitions()){
		pCD->functions.instantiate_cell = instantiate_elongated_cell;
		pCD->functions.update_phenotype = phenotype_function;
		pCD->functions.custom_cell_rule = custom_function;
		pCD->functions.contact_function = contact_function;
	}

	/*
       Cell rule definitions 
	*/

	//setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");

	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			//position[0] = Xmin + UniformRandom()*Xrange; 
			//position[1] = Ymin + UniformRandom()*Yrange; 
			//position[2] = Zmin + UniformRandom()*Zrange; 

			position[0] = 0.0; 
			position[1] = 0.0; 
			position[2] = 0.0; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
			int voxel_index = pC->get_current_voxel_index();
			microenvironment.density_vector(voxel_index)[ecm_index] = 0.0;
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

/////////////// Custom functions for the cell phenotype //////////////////

std::vector<double> calculate_migration_direction(Cell_elongated* cell)
{
    std::vector<double>& velocity = cell->velocity;

    double magnitude = norm(velocity);

    std::vector<double> migration_direction(3, 0.0);
    if (magnitude > 0)
    {
        migration_direction[0] = velocity[0] / magnitude;
        migration_direction[1] = velocity[1] / magnitude;
        migration_direction[2] = velocity[2] / magnitude;
    }
    else
    {
        // Generate a random unit vector if the cell is stationary
        migration_direction = PhysiCell::UniformOnUnitSphere();
    }
    return migration_direction;
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 	
	Cell_elongated* Elongated = dynamic_cast<Cell_elongated*>(pCell);

	std::cout << "Phenotype function" << std::endl;


	 // Calculate the migration direction
    	std::vector<double> migration_direction;

	migration_direction.assign(3, 0.0);

	migration_direction[0] = pCell->phenotype.motility.motility_vector[0];
	migration_direction[1] = pCell->phenotype.motility.motility_vector[1];
	migration_direction[2] = pCell->phenotype.motility.motility_vector[2];

	std::cout << "migration direction" << migration_direction << std::endl;


    	double voxel_size = microenvironment.mesh.dx; // Assuming cubic voxels with size dx

	// Calculate positions one voxel away in the migration direction
	std::vector<double> position_ahead = {
		Elongated->position[0] + migration_direction[0] * voxel_size,
		Elongated->position[1] + migration_direction[1] * voxel_size,
		Elongated->position[2] + migration_direction[2] * voxel_size
	};

	std::vector<double> position_behind = {
		Elongated->position[0] - migration_direction[0] * voxel_size,
		Elongated->position[1] - migration_direction[1] * voxel_size,
		Elongated->position[2] - migration_direction[2] * voxel_size
	};

	if (parameters.bools("cathepsine"))
	{
		int h_plus_index = BioFVM::microenvironment.find_density_index("h_plus");
		int cathepsine_index = BioFVM::microenvironment.find_density_index("cathepsine");

		// I need to deactivate the secretion of H+ and set the secretion of cathepsine in the voxel ahead

		Elongated->phenotype.secretion.secretion_rates[h_plus_index] = 0.0;
		Elongated->phenotype.secretion.secretion_rates[cathepsine_index] = parameters.doubles("cathepsine_secretion_rate");

		std::cout << "position ahead" << position_ahead << std::endl;

		if (microenvironment.mesh.is_position_valid(position_ahead[0], position_ahead[1], position_ahead[2]))
		{
		Elongated->simulate_secretion_and_uptake_in_voxel(&microenvironment, dt, microenvironment.nearest_voxel_index(position_ahead));
		}

		Elongated->phenotype.secretion.secretion_rates[cathepsine_index] = 0.0;
		Elongated->phenotype.secretion.secretion_rates[h_plus_index] = parameters.doubles("h_plus_secretion_rate");

		std::cout << "position behind" << position_behind << std::endl;

		if (microenvironment.mesh.is_position_valid(position_behind[0], position_behind[1], position_behind[2]))
		{
		Elongated->simulate_secretion_and_uptake_in_voxel(&microenvironment, dt, microenvironment.nearest_voxel_index(position_behind));
		}

		Elongated->phenotype.secretion.secretion_rates[cathepsine_index] = 0.0;
		Elongated->phenotype.secretion.secretion_rates[h_plus_index] = 0.0;
	}
	if (parameters.bools("MMP14"))
	{
		int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
		Elongated->phenotype.secretion.uptake_rate("ecm") = parameters.doubles("ecm_uptake_rate");

		if (microenvironment.mesh.is_position_valid(position_ahead[0], position_ahead[1], position_ahead[2]))
		{
		Elongated->simulate_secretion_and_uptake_in_voxel(&microenvironment, dt, microenvironment.nearest_voxel_index(position_ahead));
		}
		Elongated->phenotype.secretion.uptake_rate("ecm") = 0.0;
	}
	
	return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 	
	// ADDING ECM PHYSICAL INTERACTION AND ADHESION

	pCell->custom_data["ecm_contact"] = 0.0;
	pCell->custom_data["nucleus_deform"] = 0.0;

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 ){
		add_ecm_interaction( pCell, ecm_index, pCell->get_current_voxel_index() );
		//add_TGFbeta_interaction(pCell, pCell->get_current_mechanics_voxel_index());
		std::vector<int>::iterator neighbor_voxel_index;
		std::vector<int>::iterator neighbor_voxel_index_end = 
		microenvironment.mesh.moore_connected_voxel_indices[pCell->get_current_voxel_index()].end();

		for( neighbor_voxel_index = 
			microenvironment.mesh.moore_connected_voxel_indices[pCell->get_current_voxel_index()].begin();
			neighbor_voxel_index != neighbor_voxel_index_end; 
			++neighbor_voxel_index )
		{
			add_ecm_interaction( pCell, ecm_index, *neighbor_voxel_index );
			
		}

		/* pCell->update_motility_vector(dt); 
		pCell->velocity += phenotype.motility.motility_vector; */
	}
	
	return; 	
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

////////////////// Functions for the reaction-diffusion system //////////////////

// Define the reaction system
void reactionSystem(double t, std::vector<double> &dydt,
                    double k_reaction, bool consume_H, bool consume_cathepsine, int current_voxel) {

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	int h_plus = BioFVM::microenvironment.find_density_index("h_plus");
	int cathepsine = BioFVM::microenvironment.find_density_index("cathepsine");

	std::vector<double> y(3, 0.0);

	y[0] = microenvironment.density_vector(current_voxel)[ecm_index];
	y[1] = microenvironment.density_vector(current_voxel)[h_plus];
	y[2] = microenvironment.density_vector(current_voxel)[cathepsine];

    double A = y[0]; // ecm concentration
    double B = y[1]; // H+ concentration
    double C = y[2]; // cathepsine concentration

    // Reaction rate
    double reaction_rate = k_reaction * A * B * C;

    // Differential equations
	dydt[0] = -reaction_rate;                // d[C]/dt
    dydt[1] = consume_H ? -reaction_rate : 0; // d[B]/dt
	dydt[2] = consume_cathepsine ? -reaction_rate : 0; // d[A]/dt
}

// Parallelized function to update reactions in the grid
void updateReactionsInGrid(double k_reaction, double dt, bool consume_H, bool consume_cathepsine) {

	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		std::vector<double> dydt(3, 0.0);
		// Call the reaction system
		reactionSystem(0.0, dydt, k_reaction, consume_H, consume_cathepsine, n);

		 // Update concentrations using Euler's method for simplicity
		for (size_t density_index = 0; density_index < microenvironment.density_vector(n).size(); ++density_index) {
                    microenvironment.density_vector(n)[density_index] += dydt[density_index] * dt; // y = y + dydt * dt
                }
		
	}


}

////////////////// Functions for the ECM interaction //////////////////

/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void add_ecm_interaction(Cell* pC, int index_ecm, int index_voxel )
{
	Cell_elongated* pCell = dynamic_cast<Cell_elongated*>(pCell);
	// Check if there is ECM material in given voxel
	//double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = pC->get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	double ecmrad = sqrt(3.0) * pC->get_microenvironment()->mesh.dx * 0.5;
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		pC->displacement = pC->position - microenvironment.mesh.voxels[index_voxel].center;
		double distance = norm(pC->displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, EPSILON);
		
		double dd = pC->phenotype.geometry.radius + ecmrad;  
		double dnuc = pC->phenotype.geometry.nuclear_radius + ecmrad;  

		double tmp_r = 0;
		// Cell overlap with ECM node, add a repulsion term
		if ( distance < dd )
		{
			// repulsion stronger if nucleii overlap, see Macklin et al. 2012, 2.3.1
			if ( distance < dnuc )
			{
				double M = 1.0;
				double c = 1.0 - dnuc/dd;
				c *= c;
				c -= M;
				tmp_r = c*distance/dnuc + M;
				pC->custom_data["nucleus_deform"] += (dnuc-distance);
			}
			else
			{
				tmp_r = ( 1 - distance / dd );
				tmp_r *= tmp_r;
			}
			tmp_r *= dens * parameters.doubles("cell_ecm_repulsion_strength");
		}

		// Cell adherence to ECM through integrins
		double max_interactive_distance = (pC->phenotype.mechanics.relative_maximum_adhesion_distance) + ecmrad;
		if ( distance < max_interactive_distance ) 
		{	
			double temp_a = 1 - distance/max_interactive_distance; 
			temp_a *= temp_a; 
			/* \todo change dens with a maximal density ratio ? */

			pC->custom_data["ecm_contact"] += dens * (max_interactive_distance-distance);

			temp_a *= dens * parameters.doubles("cell_ecm_adhesion_strength");
			
			tmp_r -= temp_a;
		}
		
		/////////////////////////////////////////////////////////////////
		if(tmp_r==0)
			return;
		tmp_r/=distance;

		axpy( &pC->velocity , tmp_r , pC->displacement ); 
	}

}

Cell* instantiate_elongated_cell() { return new Cell_elongated; }
