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
#define	Lx		1.5
#define	Ly		0.5
#define	Lz		0.01

// grid
#define	xCells		90
#define	yCells		30
#define	zCells		1

// material properties
#define	K		1.5e+1
#define rho		1000.0

// velocity field
#define	Ux		1.0
#define	UxProfile	"parabolic"
#define	Uy		0.0
#define	UyProfile	"uniform"

// boundary conditions for scalar
#define westBCtype	"Dirichlet"
#define westBCvalue	20.0
#define eastBCtype	"Neumann"
#define eastBCvalue	0.0
#define southBCtype	"Neumann"
#define southBCvalue	350.0
#define northBCtype	"Neumann"
#define northBCvalue	200.0

// numerical methods
#define linearSolver	BiCGSTAB
#define maxIter		1000
#define tolerance	1e-6


// =======================================================
// 		     global variables
// =======================================================
int	cells=xCells*yCells;
double	dx=Lx/xCells, dy=Ly/yCells, dz=Lz/zCells;
SparseMatrix<double> 	M(cells,cells);
VectorXd 		B(cells), X(cells);
double	Ap, Fw, Fe, Fs, Fn, Dw, De, Ds, Dn, Aw, Ae, As, An, Sp, Su;
int i=0, j, k;
double u[xCells][yCells], v[xCells][yCells];



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
	// -------------------------------------------
	//	     setup boundary conditions
	// -------------------------------------------
	BoundaryCondition westBC, eastBC, southBC, northBC;
	westBC.setup("west", westBCtype, westBCvalue);
	eastBC.setup("east", eastBCtype, eastBCvalue);
	southBC.setup("south", southBCtype, southBCvalue);
	northBC.setup("north", northBCtype, northBCvalue);


	// -------------------------------------------
	//	          print info
	// -------------------------------------------
	cout << endl << "Domain: " << Lx << " X " << Ly << " X " << Lz << " m^3." << endl;
	cout << "Grid: "  << xCells << " X " << yCells << " X " << zCells << " = " <<
		 cells << " cells." << endl;
	cout << "Cell: "  << dx << " X " << dy << " X " << dz << " m^3." << endl;


	// -------------------------------------------
	//	     setup velocity field
	// -------------------------------------------
	for (int k=0; k<xCells; k++)
		for (int j=0; j<yCells; j++)
		{
			if (UxProfile=="uniform")	u[k][j] = Ux;
			else	u[k][j] = Ux * (1 - pow( double(yCells+1)/2.0 - double(j+1), 2 )
					/ pow( double(yCells+1)/2.0, 2 ) );

			if (UyProfile=="uniform")	v[k][j] = Uy;
			else	v[k][j] = Uy  * (1 - pow( double(xCells+1)/2.0 - double(k+1), 2 )
					/ pow( double(xCells+1)/2.0, 2 ) );
		};
	cout << endl;


	// -------------------------------------------
	//		 build matrix
	// -------------------------------------------
	cout << "Building matrix..." << flush;

	// loop all cells
	for (i=0; i<cells; i++)
	{
		// coord's for current cell
		j = i % yCells;
		k = i / yCells;

		// west domain boundary?
		if (i<yCells)			Fw = rho * u[k][j] * dy * dz;
		else				Fw = rho * ( u[k-1][j] + u[k][j] ) / 2.0 * dy * dz;
	    
		// east domain boundary?
		if (i>=(xCells-1)*yCells)	Fe = rho * u[k][j] * dy * dz;
		else				Fe = rho * ( u[k][j] + u[k+1][j] ) / 2.0 * dy * dz;
    
		// south domain boundary?
		if ((i+1) % yCells == 1)	Fs = rho * v[k][j] * dx * dz;
		else				Fs = rho * ( v[k][j-1] + v[k][j] ) / 2.0 * dx * dz;
	    
		// north domain boundary?
		if ((i+1) % yCells == 0)	Fn = rho * v[k][j] * dx * dz;
		else				Fn = rho * ( v[k][j] + v[k][j+1] ) / 2.0 * dx * dz;

		Dw = K / dx * dy * dz;
		De = K / dx * dy * dz;
		Ds = K / dy * dx * dz;
		Dn = K / dy * dx * dz;

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
	//		write results
	// -------------------------------------------
	cout << "Writing to files..." << flush;
	ofstream outFile;

	// ------------  Ux field  -------------------
	outFile.open("Ux");

	for (int j=yCells-1; j>=0; j--)
	{
		for (int k=0; k<xCells; k++)
		{
			outFile << u[k][j];
			if (k<xCells-1) outFile << " ";
		}
		outFile << endl;
	}
	outFile << endl;
	outFile.close();

	// ------------  Uy field  -------------------
	outFile.open("Uy");

	for (int j=0; j<yCells; j++)
	{
		for (int k=0; k<xCells; k++)
		{
			outFile << v[k][j];;
			if (k<xCells-1) outFile << " ";
		}
		outFile << endl;
	}
	outFile << endl;
	outFile.close();


	// -------------  T field  -------------------
	outFile.open("T");

	for (int i=0; i<cells; i++)
	{
		outFile << X(i);
		if ((i+1) % yCells == 0) outFile << endl;
			else outFile << " ";
	}
	outFile << endl;

	cout << "done." << endl << endl;
	outFile.close();

	system("octave plotFields.m");

};

