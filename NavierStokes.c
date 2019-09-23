/*========================================================*\
|							   |
|		 Navier-Stokes solver			   |
|		 finite volume method			   |
|							   |
\*========================================================*/

#include <iostream>
#include <fstream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;



// =======================================================
// 			case parameters
// =======================================================
// domain
#define	Lx		2.5
#define	Ly		0.5
#define	Lz		0.01

// grid
#define	xCells		50
#define	yCells		10
#define	zCells		1

// material properties
#define rho		1000.0
#define nu		1e-6

// velocity field
#define	Ux		0.0
#define	Uy		0.0

// boundary conditions for x-velocity
#define westBCtype	"Dirichlet"
#define westBCvalue	1.0
#define eastBCtype	"Neumann"
#define eastBCvalue	0.0
#define southBCtype	"Dirichlet"
#define southBCvalue	0.0
#define northBCtype	"Dirichlet"
#define northBCvalue	0.0

// numerical methods
#define linearSolver	BiCGSTAB
#define maxIter		1000
#define tolerance	1e-8


// =======================================================
// 		     global variables
// =======================================================
int	cells=xCells*yCells;
double	dx=Lx/xCells, dy=Ly/yCells, dz=Lz/zCells;
SparseMatrix<double> 	M(cells,cells);
VectorXd 		B(cells), X(cells);
double	Ap, Fw, Fe, Fs, Fn, Dw, De, Ds, Dn, Su;
int i=0, j, k;
double u[yCells][xCells], v[yCells][xCells], p[yCells][xCells];



// =======================================================
// 			 main code 
// =======================================================
int main()
{
	// -------------------------------------------
	//		intialize fields
	// -------------------------------------------
	for (int i=0; i<yCells; i++)
		for (int j=0; j<xCells; j++)
		{
			u[i][j] = Ux;
			v[i][j] = Uy;
			p[i][j] = double(xCells-j) / double(xCells);
		};

	// -------------------------------------------
	//	        print mesh info
	// -------------------------------------------
	cout << endl << "Domain: " << Lx << " X " << Ly << " X " << Lz << " m^3." << endl;
	cout << "Grid: "  << xCells << " X " << yCells << " X " << zCells << " = " <<
		 cells << " cells." << endl;
	cout << "Cell: "  << dx << " X " << dy << " X " << dz << " m^3." << endl;


	// -------------------------------------------
	//		assemble matrix
	// -------------------------------------------
	cout << "Assembling matrix..." << flush;

	// loop all cells
	for (i=0; i<cells; i++)
	{
		// coord's for current cell
		j = i / xCells; // y
		k = i % xCells; // x

		// diffusion
		Dw = nu / dx * dy * dz;
		De = nu / dx * dy * dz;
		Ds = nu / dy * dx * dz;
		Dn = nu / dy * dx * dz;

		// face fluxes
		if (k==0)		Fw = u[j][k] * dy * dz;
		else			Fw = ( u[j][k-1] + u[j][k] ) / 2.0 * dy * dz;
		if (k==xCells-1)	Fe = u[j][k] * dy * dz;	
		else			Fe = ( u[j][k] + u[j][k+1] ) / 2.0 * dy * dz;
		if (j==yCells-1)	Fs = v[j][k] * dx * dz;
		else			Fs = ( v[j+1][k] + v[j][k] ) / 2.0 * dx * dz;
		if (j==0)		Fn = v[j][k] * dx * dz;
		else			Fn = ( v[j][k] + v[j-1][k] ) / 2.0 * dx * dz;

		// reset source term
		Su = -1.0 / rho * 100.0;//( p[j][k+1] - p[j][k-1] );

		// reset contribution to current cell
		Ap = 0;

		// west boundary?
		if (k==0)
		{
			// Dirichet
			Su += +Fw * westBCvalue  +2.0*Dw * westBCvalue;
			Ap += +2.0*Dw;
		}
		else
		{
			Ap += -Fw/2.0 + Dw;
			M.insert(i,i-1) = -Fw/2.0 -Dw;
		}
	    
		// east boundary?
		if (k==xCells-1)
		{
			// Neumann
			Su += eastBCvalue;
		}
		else
		{
			Ap += +Fe/2.0 + De;
			M.insert(i,i+1) = +Fe/2.0 -De;
		}
    
		// south boundary?
		if (j==yCells-1)
		{
			// Dirichet
			Su += +Fs * southBCvalue +2.0*Ds * southBCvalue;
			Ap += -2.0*Ds;
		}
		else
		{
			Ap += -Fs/2.0 + Ds;
			M.insert(i,i+yCells) = -Fw/2.0 -Ds;
		}
	    
		// north boundary?
		if (j==0)
		{
			// Dirichet
			Su += -Fn * northBCvalue +2.0*Dn * northBCvalue;
			Ap += +2.0*Dn;
		}
		else
		{
			Ap += +Fn/2.0 + Dn;
			M.insert(i,i-yCells) = +Fn/2.0 -Dn;
		}
		
		// current cell
		M.insert(i,i) = Ap;

		B(i) = Su;
	};	
	cout << "done." << endl;
	

	// -------------------------------------------
	//		solve M*X=B
	// -------------------------------------------
	cout << "Solving matrix..." << flush;
	linearSolver< SparseMatrix<double> >  solver;
	solver.setMaxIterations(maxIter);
	solver.setTolerance(tolerance);
	solver.compute(M);
	X = solver.solve(B);
	cout << "done " << " (" << solver.error() << ", " << solver.iterations() << ")." << endl;


	// -------------------------------------------
	//		write fields
	// -------------------------------------------
	cout << "Writing results..." << flush;
	ofstream outFile;

	// -------------  p field  -------------------
	outFile.open("p");
	for (int j=0; j<yCells; j++)
	{
		for (int k=0; k<xCells; k++)
		{
			outFile << p[j][k];
			if (k<xCells-1) outFile << " ";
		}
		outFile << endl;
	}
	outFile << endl;
	outFile.close();

	// ------------  Ux field  -------------------
	outFile.open("Ux");
	for (int j=0; j<yCells; j++)
	{
		for (int k=0; k<xCells; k++)
		{
			outFile << X[k*yCells+j];
			if (k<xCells-1) outFile << " ";
		}
		outFile << endl;
	}
	outFile << endl;
	outFile.close();

	// ----------------- plot --------------------
	system("gnuplot plotFields.gp");
};

