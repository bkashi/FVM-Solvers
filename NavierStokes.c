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
#define westBCvalue	1e-4
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
double	Ap, Fw, Fe, Fs, Fn, Dw, De, Ds, Dn, Aw, Ae, As, An, Sp, Su;
int i=0, j, k;
double u[yCells][xCells], v[yCells][xCells], p[yCells][xCells];



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
	//	     setup boundary conditions
	// -------------------------------------------
	BoundaryCondition westBC, eastBC, southBC, northBC;
	westBC.setup("west", westBCtype, westBCvalue);
	eastBC.setup("east", eastBCtype, eastBCvalue);
	southBC.setup("south", southBCtype, southBCvalue);
	northBC.setup("north", northBCtype, northBCvalue);


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

		// west domain boundary?
		if (k==0)			Fw = u[j][k] * dy * dz;
		else				Fw = ( u[j][k-1] + u[j][k] ) / 4.0 * dy * dz;
	    
		// east domain boundary?
		if (k==xCells-1)		Fe = u[j][k] * dy * dz;
		else				Fe = ( u[j][k] + u[j][k+1] ) / 4.0 * dy * dz;
    
		// south domain boundary?
		if (j==yCells-1)		Fs = v[j][k] * dx * dz;
		else				Fs = ( v[j+1][k] + v[j][k] ) / 4.0 * dx * dz;
	    
		// north domain boundary?
		if (j==0)			Fn = v[j][k] * dx * dz;
		else				Fn = ( v[j][k] + v[j-1][k] ) / 4.0 * dx * dz;

		Dw = nu / dx * dy * dz;
		De = nu / dx * dy * dz;
		Ds = nu / dy * dx * dz;
		Dn = nu / dy * dx * dz;

		Aw = -Dw -Fw;
		Ae = -De +Fe;
		As = -Ds -Fs;
		An = -Dn +Fn;

		Sp = 0;
		Su = -1.0 * ( p[j][k+1] - p[j][k-1] ) / rho;

		// west domain boundary?
		if (k==0)			westBC.apply();
		else				M.insert(i,i-1) = Aw;    
	    
		// east domain boundary?
		if (k==xCells-1)		eastBC.apply();
		else				M.insert(i,i+1) = Ae;
    
		// south domain boundary?
		if (j==yCells-1)		southBC.apply();
		else				M.insert(i,i+xCells) = As;
	    
		// north domain boundary?
		if (j==0)			northBC.apply();
		else				M.insert(i,i-xCells) = An;

		// current cell
		Ap = Aw + Ae + As + An - Sp;
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

