// Functions which are usefull for the scalar convection equations

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



//This Function Initialises the Finite Element Data Structures for a 2D Problem
void InitialiseFELibraries(TDatabase *Database,TFEDatabase2D *FEDatabase)
{
    Database = new TDatabase();
	FEDatabase = new TFEDatabase2D();
}




// This Function reads the mesh and converts them into the Finite Element DataStructures
void GenerateMeshAndRefine(TDomain *Domain)
{
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
        OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
        OutPut("GMSH used for meshing !!!" << endl);
    }                                            // gmsh mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 2) // triangle mesh
    {
        OutPut("Triangle.h used for meshing !!!" << endl);
        // TriaReMeshGen(Domain);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }


    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();  
    
    if(TDatabase::ParamDB->WRITE_PS)
        Domain->PS("Domain.ps", It_Finest, 0);
    
    if(TDatabase::ParamDB->WRITE_VTK)
        mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

}


void GenerateScalarFESpace(TFESpace2D* Scalar_FeSpace, TCollection* coll, BoundCondFunct2D *BoundaryCondition,int order)
{
    Scalar_FeSpace = new TFESpace2D(coll, (char*)"ScalarFESpace", (char*)"FESapace for Scalar Problem in 2D", BoundaryCondition, order, NULL);
}