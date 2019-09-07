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
#define	Lx 0.6
#define	Ly 0.4
#define	Lz 0.01

// grid
#define	xCells 60
#define	yCells 40
#define	zCells 1

// material properties
#define	K 1.5e+1
#define rho 1000.0

// velocity field
#define	Ux 1.0
#define	Uy 0.0

// numerical methods
#define linearSolver	BiCGSTAB
#define maxIter		1000
#define tolerance	1e-6

// boundary conditions	
#define westBCtype "Dirichlet"
#define westBCvalue 0.0
#define eastBCtype "Neumann"
#define eastBCvalue 0.0
#define southBCtype "Dirichlet"
#define southBCvalue 100.0
#define northBCtype "Dirichlet"
#define northBCvalue 100.0


// =======================================================
// 		     global variables
// =======================================================
int	cells=xCells*yCells;
double	dx=Lx/xCells, dy=Ly/yCells, dz=Lz/zCells;
SparseMatrix<double> 	M(cells,cells);
VectorXd 		B(cells), X(cells);
double	Ap, Fw, Fe, Fs, Fn, Dw, De, Ds, Dn, Aw, Ae, As, An, Sp, Su;
int i=0;



// =======================================================
// 			classes
// =======================================================

class BoundaryCondition
{
	private:
		string	location;
		string	type;
		double	value;

		void Dirichlet()
		{
			if (location=="west")
			{
				Aw  =	0;
				Sp +=	-2*Dw - Fw;
				Su +=	(2*Dw + Fw) * westBCvalue;
			}

			if (location=="east")
			{
				Ae  =	0;
				Sp +=	-2*De + Fe;
				Su +=	(2*De - Fe) * eastBCvalue;
			}

			if (location=="south")
			{
				As  =	0;
				Sp +=	-2*Ds - Fs;
				Su +=	(2*Ds + Fs) * southBCvalue;
			}

			if (location=="north")
			{
				An  =	0;
				Sp +=	-2*Dn + Fn;
				Su +=	(2*Dn - Fn) * northBCvalue;
			}
		};

		void Neumann()
		{
			if (location=="west")
			{
				Aw  =	0;
				Su +=	dy * dz * westBCvalue;
			}

			if (location=="east")
			{
				Ae  =	0;
				Su +=	dy * dz * eastBCvalue;
			}

			if (location=="south")
			{
				As  =	0;
				Su +=	dx * dz * southBCvalue;
			}

			if (location=="north")
			{
				An  =	0;
				Su +=	dx * dz * northBCvalue;
			}
		};

	public:
		void setup(string aLocation, string aType, double aValue)
		{
			type =		aType;
			value =		aValue;
			location =	aLocation;
		};

		void apply()
		{
			if (type=="Dirichlet")	Dirichlet();
			if (type=="Neumann")	Neumann();
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

	// setup velocity field
	// ...

	// print info
	cout << endl << "Domain: " << Lx << " X " << Ly << " X " << Lz << " m^3." << endl;
	cout << "Grid: "  << xCells << " X " << yCells << " X " << zCells << " = " <<
		 cells << " cells." << endl;
	cout << "Cell: "  << dx << " X " << dy << " X " << dz << " m^3." << endl;


	// -------------------------------------------
	//		build matrix
	// -------------------------------------------
	cout << "Building matrix..." << flush;


	Fw = rho * Ux * dy * dz;
	Fe = rho * Ux * dy * dz;
	Fs = rho * Uy * dx * dz;
	Fn = rho * Uy * dx * dz;

	Dw = K / dx * dy * dz;
	De = K / dx * dy * dz;
	Ds = K / dy * dx * dz;
	Dn = K / dy * dx * dz;


	// loop all cells
	while (i<cells)
	{
		Aw = Dw + Fw/2;
		Ae = De - Fe/2;
		As = Ds + Fs/2;
		An = Dn - Fn/2;

		Sp = 0;
		Su = 0;

		// west domain boundary?
		if (i<yCells)			westBC.apply();
		else				M.insert(i,i-yCells) = -Aw;    
	    
		// east domain boundary?
		if (i>=(xCells-1)*yCells)	eastBC.apply();
		else				M.insert(i,i+yCells) = -Ae;
    
		// south domain boundary?
		if ((i+1) % yCells == 1)	southBC.apply();
		else				M.insert(i,i-1) = -As;
	    
		// north domain boundary?
		if ((i+1) % yCells == 0)	northBC.apply();
		else				M.insert(i,i+1) = -An;

		// current cell
		Ap = Aw + Ae + As + An + Fe - Fw + Fn - Fs - Sp;
		M.insert(i,i) = Ap;

		B(i) = Su;

		// next cell
		i++;
	};	

	cout << "done." << endl;
	//cout << M << endl;
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

	system("octave plot.m");

};

