// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Thivin Anandh
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <MainUtilities.h>

#include <stdlib.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>


// Include Helper Files
#include "../HelperFunctions/ScalarHelperFunctions.h"

// Include Your example File 
#include "../Examples/CD_2D/SineLaplace.h"

int main(int argc, char *argv[])
{
	// Initialise the FE Libraries
	TDatabase *Database;
	TFEDatabase2D *FEDatabase;
	InitialiseFELibraries(Database,FEDatabase);

	;
	TFESpace2D *Scalar_FeSpace, *fesp[1];
	TFEFunction2D *Scalar_FeFunction;
	TSystemCD2D *SystemMatrix;
	TOutput2D *Output;
	MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

	std::ostringstream os;
	os << " ";

	// Create a Domain with the values given by the Parmeter files 
	TDomain *Domain = new TDomain(argv[1]);

	// Read the Mesh File and refine them 
	GenerateMeshAndRefine(Domain);


	// Get the Collection of Geometric Cells from the mesh 
	TCollection *coll = Domain->GetCollection(It_Finest, 0);

	// Obtain all the finite element Parameters 
	int ORDER  =  TDatabase::ParamDB->ANSATZ_ORDER;
	int N_Cells = coll->GetN_Cells();

	// Declare a Finite Element Space with the Given Order and Boundary Condition
	// Scalar_FeSpace = 




}
