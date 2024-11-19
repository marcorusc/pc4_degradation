#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;


bool isElongated(PhysiCell::Cell_Definition * cellDef);

std::vector<PhysiCell::Cell_Definition*>* getElongatedCellDefinitions();

class Cell_elongated : public Cell
{
public:
    // Constructor
    Cell_elongated();
    ~Cell_elongated() {};

    // Any other custom methods or properties

    void custom_cell_elongated_update_phenotype(PhysiCell::Cell *pCell, PhysiCell::Phenotype &phenotype, double dt);

    void simulate_secretion_and_uptake_in_voxel(Microenvironment* pS, double dt, int voxel_index );
};