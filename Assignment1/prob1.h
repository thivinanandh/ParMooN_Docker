// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
    OutPut("Example: DrivenCavity.h" << endl);
}

// ========================================================================
// exact solution   -- NO need to take into account for this Assignment
// ========================================================================
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

///////////// BOUNDARY COMPONENT ID  /////////////////////////////////
/* 
    FOr a Square Domain, The edge y= 0    ID  = 0
                         The Edge x = 1   ID  = 1
                         The edge y = 1   ID  = 2
                         The edge x = 0   ID  = 3   
                         Inner Circle     ID  = 4
*/
//////


/////////// HEMKER PROBLEM /////////////////////////
// Simple perturbed Eliptic Equation 
//                       Eps( grad)^u + Ux  = 0 
// Similar to Convection diffusion advection equation ( CD2D) with b1 = 1 , b2 = 0, c = 0 
// Boundary Conditions, at the unit circle, the Value of unknown wil be "1" ( Diriclet )
// Its zero on inlet and walls, No flux on outlet



// ========================================================================
// boundary conditions
// Provide The Boundary Condition as "DIRICHLET" or "NEUMANN"
// ========================================================================
// TO DO
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
   /*
    Function provides Boundary Condition on the edges of the domain. 
    for eg, if you want to provide DIRICHLET for Egde where y = 1 ( ID  = 2) and Neuman every where
    Then use either "if - else if  -else" ( as below ) or Switch case based on value "i"
        if(i == 0)
            cond  =  DIRICHLET
        else
            cond  =  NEUMANN
    Note : Make sure to give Bound cond for all Edges, else the execution will fail. 
    */
}

// Boundary  value at the Boundary Edges  // TO DO
void BoundValue(int BdComp, double Param, double &value)
{
    /*
    Give the Value at the Boundary based on the BdComp (Boundary COnponent ID)
    For eg : For the Edge y = 1, if I need to provide value =  1, then use, 
    either Switch condition or an Ifelse cond to apply bd Values. 
     if(BdComp == 0)
        value = 1.0;
    Note : Make sure to give Bound Value for all Edges, else the execution will fail. 
 */
}


// DO NOT  CHANGE THIS FUNCTION for Assignment
void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff, *param;

  v1 = 1.0;
  v2 = 0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;     // ( 1/Pe )
    coeff[1] = v1;      // ( b1     )
    coeff[2] = v2;      // b2
    coeff[3] = 0;       // c

    coeff[4] = 0;       // f
  }
}


// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
void HemkerAssembly(double quad_wt, double *coeff, double *param,
                    double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{

	// The below values N_x, N_y , etc are ARRAY of Values 
    // Which provides Value of a Particular Shape function or its derivative ( Based on a DOF )
    // at the given Quadrature POint. 
    //
    // Here the Size of the Array is equal to the NUmber of DOF in the cell 
    
    // For Eg : Orig0[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
    // at the given Quadrature Point.

    double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2];
    
	
	double **A11, *F1;
    double val = 0.;
    
	
	A11 = LocMatrices[0];  // Local Stiffenss matrix

    F1 = LocRhs[0];


    double c0 = coeff[0]; // nu
    double b1 = coeff[1]; // b1
    double b2 = coeff[2]; // b2
    double c  = coeff[3]; // c
    double f  = coeff[4]; //f

    int N_DOF_perCell = N_BaseFuncts[0];

    for (int i = 0; i < N_DOF_perCell; i++)
    {
        // Assemble RHS 
        F1[i] =  ; // TO DO
        for (int j = 0; j < N_DOF_perCell; j++)
        {

            A11[i][j]  = //   TO DO 
        } // endfor j
    }     // endfor i

}


// Ignore this Function for assignment
void NSType1Galerkin(double Mult, double *coeff,
                     double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
}