
#include "Cell_elongated.h"


bool isElongated(PhysiCell::Cell_Definition * cellDef)
{
    const auto agentname = std::string(cellDef->name);
    const auto ecm = std::string("elongated");
    const auto matrix = std::string("elongated_cell");
    const auto fiber = std::string("cell_elongated");
    const auto fibre = std::string("long_cell");

    return (agentname.find(ecm) != std::string::npos ||
        agentname.find(matrix) != std::string::npos ||
        agentname.find(fiber) != std::string::npos ||
        agentname.find(fibre) != std::string::npos 
    );
}

std::vector<PhysiCell::Cell_Definition*>* getElongatedCellDefinitions() {
    std::vector<PhysiCell::Cell_Definition*>* result = new std::vector<PhysiCell::Cell_Definition*>();
    PhysiCell::Cell_Definition* pCD;
    
    
    for (auto& cd_name: PhysiCell::cell_definitions_by_name) {
        if (isElongated(cd_name.second)) {
            result->push_back(cd_name.second);
        }
    }
    
    return result;
}

Cell_elongated::Cell_elongated() : Cell()
{
    // Custom initialization for Cell_elongated
    std::cout << "Cell_elongated constructor" << std::endl;
}


void Cell_elongated::simulate_secretion_and_uptake_in_voxel(Microenvironment* pS, double dt, int voxel_index )
{
    if(!this->is_active)
    { return; }

    std::cout << "Cell_elongated::simulate_secretion_and_uptake_in_voxel: " << voxel_index << std::endl;

    this->set_internal_uptake_constants(dt);
    
    
    if( default_microenvironment_options.track_internalized_substrates_in_each_agent == true )
    {

        total_extracellular_substrate_change.assign( total_extracellular_substrate_change.size() , 1.0 ); // 1

        total_extracellular_substrate_change -= cell_source_sink_solver_temp2; // 1-c2
        total_extracellular_substrate_change *= (*pS)(voxel_index); // (1-c2)*rho 
        total_extracellular_substrate_change += cell_source_sink_solver_temp1; // (1-c2)*rho+c1 
        total_extracellular_substrate_change /= cell_source_sink_solver_temp2; // ((1-c2)*rho+c1)/c2
        total_extracellular_substrate_change *= pS->voxels(voxel_index).volume; // W*((1-c2)*rho+c1)/c2 
        
        *internalized_substrates -= total_extracellular_substrate_change; // opposite of net extracellular change 	
    }
    
    (*pS)(voxel_index) += cell_source_sink_solver_temp1; 
    (*pS)(voxel_index) /= cell_source_sink_solver_temp2; 

    std::cout << "Cell_elongated::simulate_secretion_and_uptake_in_voxel: " << cell_source_sink_solver_temp1 << std::endl;
    std::cout << "Cell_elongated::simulate_secretion_and_uptake_in_voxel: " << cell_source_sink_solver_temp2 << std::endl;

    
    // now do net export 
    (*pS)(voxel_index) += cell_source_sink_solver_temp_export2; 

    return; 
}











