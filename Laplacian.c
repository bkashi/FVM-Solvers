/*========================================================*\
|							   |
|		   Laplacian solver			   |
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
#define	Lx 0.3
#define	Ly 0.4
#define	Lz 0.01

// grid
#define	xCells 30
#define	yCells 40
#define	zCells 1
int	cells=xCells*yCells;
double	dx=Lx/xCells, dy=Ly/yCells, dz=Lz/zCells;

// material properties
double	k = 1000.0;

// numerical methods
#define linearSolver	BiCGSTAB
#define maxIter		1000
#define tolerance	1e-6

// boundary conditions	
#define westBCtype "Dirichlet"
#define westBCvalue 100.0
#define eastBCtype "Dirichlet"
#define eastBCvalue 100.0
#define southBCtype "Neumann"
#define southBCvalue 500000.0
#define northBCtype "Dirichlet"
#define northBCvalue 100.0


// =======================================================
// 		     global variables
// =======================================================
SparseMatrix<double> 	M(cells,cells);
VectorXd 		B(cells), X(cells);
double	W, E, S, N, P, C;
int i=0;



// =======================================================
// 			classes
// =======================================================
class BoundaryCondition
{
		string	location;
		string	type;
		double	value;

	public:
		void setup(string aLocation, string aType, double aValue)
		{
			type = aType;
			value = aValue;
			location = aLocation;
		};

		void apply()
		{
			double C=0;

			if (location=="west")	C = W;
			if (location=="east")	C = E;
			if (location=="south")	C = S;
			if (location=="north")	C = N;

			if (type=="Neumann")
			{
				B(i) += value * dy * dz;
				P = P - C;
			};

			if (type=="Dirichlet")
			{
				B(i) += 2.0 * C * value;
				P = P + C;
			};
		};
};



// =======================================================
// 			 main code 
// =======================================================
int main()
{
	// setup boundary conditions	
	BoundaryCondition westBC, eastBC, southBC, northBC;
	westBC.setup("west", westBCtype, westBCvalue);
	eastBC.setup("east", eastBCtype, eastBCvalue);
	southBC.setup("south", southBCtype, southBCvalue);
	northBC.setup("north", northBCtype, northBCvalue);

	// print info
	cout << endl << "Domain: " << Lx << " X " << Ly << " X " << Lz << " m^3." << endl;
	cout << "Grid: "  << xCells << " X " << yCells << " X " << zCells << " = " <<
		 cells << " cells." << endl;
	cout << "Cell: "  << dx << " X " << dy << " X " << dz << " m^3." << endl;


	// -------------------------------------------
	//		build matrix
	// -------------------------------------------
	cout << "Building matrix..." << flush;

	// diffusion coefficients
	W = k * dy * dz / dx; 
	E = k * dy * dz / dx;
	S = k * dx * dz / dy;
	N = k * dx * dz / dy;

	// loop all cells
	//for (int i=0; i<cells; i++)
	while (i<cells)
	{
	        P = W + E + S + N;

		// west domain boundary?
		if (i<yCells) westBC.apply();
		else M.insert(i,i-yCells) = -W;    
	    
		// east domain boundary?
		if (i>=(xCells-1)*yCells) eastBC.apply();
		else M.insert(i,i+yCells) = -E;
    
		// south domain boundary?
		if ((i+1) % yCells == 1) southBC.apply();
		else M.insert(i,i-1) = -S;
	    
		// north domain boundary?
		if ((i+1) % yCells == 0) northBC.apply();
		else M.insert(i,i+1) = -N;

		// current cell
		M.insert(i,i) = P;

		// next cell
		i++;
	};	

	cout << "done." << endl;
	//cout << M<< endl;
	//cout << B << endl;
	

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
	//cout << X << endl;


	// -------------------------------------------
	//		write solution
	// -------------------------------------------
	cout << "Writing to file..." << flush;
	ofstream outFile;
	outFile.open("output.txt");

	for (int i=0; i<cells; i++)
	{
		outFile << int(X(i));
		if ((i+1) % yCells == 0) outFile << endl;
			else outFile << " ";
	}
	outFile << endl;

	cout << "done." << endl << endl;
	outFile.close();

};

