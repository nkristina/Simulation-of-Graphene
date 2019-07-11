// THis code generates the configuration for MD simulation
// It writes atomic coordinates into the DATA.FILE, which is the initial configuration for the LAMMPS script
// It writes output .txt files which contain information about all input parameters
// It writes the PARM.FILE included in the LAMMPS script



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


// Global structure which represents the position of an atom (x, y, z)
struct atom{
	double x;
	double y;
	double z;
};

// Global structure which represents one group of atoms when dividing the plate into different groups with different atom_types
struct Group{
	double B;
	double T;
};


/* Prototypes of functions */

int main(); // main function

void MakeFCC111(int n1y, int nL, double D, double LL[], struct atom r1[]); // makes FCC (111) plate

void MakeGraphene(double Lx, double Ly, double a, struct atom r3[]); // makes flat layer of graphene

void MakeHydrogen(int N4,struct atom r4[], double Lx, double Ly, double Lz);

double maxAr(int len, double array[]); // finds a maximum of a given array

double minAr(int len, double array[]); // finds a minimum of a given array

double fmax(double a, double b); // finds a maximum out of two numbers

double fmin(double a, double b); // finds a minimum out of two numbers

double dist(double r1[], double r2[]); // calculates the distance between two points

double rand2(double a, double b); // finds a random number between two  numbers

double grapheneLength(double Ly); // calculates the length of graphene along the y-direction for a given length of system: L = Ly/2, D = L/3, H = D/3

void Output(int n1x, int n1y, int N1, int nL, int N3, double D, double sigmaP, double charge, double mP, double epsPP, double K, int PD, int nLPD[], double delz, double Lx, double Ly, struct atom r1[], struct atom r3[], int N4, struct atom r4[]); // writes the output files

/****************************************************************************************************************************************/
void MakeHydrogen(int N4, struct atom r4[],double Lx, double Ly, double Lz)
{
	int i;
	const double pi=3.141592654;

	double D = Ly/3;
	double H = Ly/9;

	for (i=0; i<N4; i++)
	{
		r4[i].x = rand2(0,Lx);
		r4[i].y = rand2(0,Ly);
		if ((r4[i].y >= (Ly/2 - D)) && (r4[i].y <= (Ly/2 + D))){
			r4[i].z = 3 + rand2(H*(1 + cos((r4[i].y - Ly/2)*pi/D)),25*D+Lz);
		}
		else r4[i].z = 3 + rand2(0,25*D+Lz);

	}
}

int main(){

	// global variables
	int nL; // total number of PLATE layers

	int n1x, n1y, n1z; // number of atoms along all three directions of the BOTTOM plate
	int N1; // total number of atoms in the BOTTOM plate

	int PD; // number of different groups of atoms in a plate => PlateDivision (PD)

	int tmp; // temporary variable

	double D; // center-to-center distance between consecutive atoms for PLATE atoms
	double Lx, Ly, Lz; // length of the BOTTOM plate in all three directions
	double LL[3]; // array which contains information about the plate dimensions
	double charge; // absolute value of the charge

	double sigmaP; // size of the LJ_sigma parameter: P - plates

	double mP; // mass of plate atoms

	double K; // bond_coeff for atom-ghost_atom pairs connected with elastic springs

	int i, j, k; // counter

	double epsPP; // values of the epsilon parameters of corresponding LJ interactions
	// epsPP = eps for atoms within the plates

	int N4; // number of hydrogen atoms

	FILE *fp; // file to which we output the values of input parameters
	char filename[64]; // name of that file



/****************************************************************************************************************************************/

	/* Read the input parameters */

	printf("INPUT PROCEDURE\n");

	printf("\n\n\n");
	printf("ENTER THE VALUES OF PARAMETERS\n");
	printf("\n\n\n");

	printf("Number of atoms along the y-axis of BOTTOM plate: \n");
	scanf("%d", &n1y);

	printf("Size of the LJ_sigma parameter in angstroms (PLATES): \n");
	scanf("%lf", &sigmaP);

	printf("Absolute value of the charge: \n");
	scanf("%lf", &charge);

	printf("Mass of PLATES' ATOMS (mP): \n");
	scanf("%lf", &mP);

	printf("Epsilon parameter for interactions within the plates, epsPP: \n");
	scanf("%lf", &epsPP);

	printf("Number of groups to which each plate is divided, PD: \n");
	scanf("%d", &PD);

	int nLPD[PD]; // define the array of PD members - they contain the number of layers in each group of a plate

	for (i=0; i<PD; i++){
		printf("Number of layers in %d group of a plate:\n", i+1);
		scanf("%d", &tmp);
		nLPD[i] = tmp;

	}
	printf("\n");

	printf("Bonf coeff, K: \n");
	scanf("%lf", &K);

	printf("\n");

	printf("Number of hydrogen atoms, N4: \n");
	scanf("%d",&N4);

/****************************************************************************************************************************************/

	/* Write the input parameters to the screen as a check */

	printf("\n\n\n\n");

	printf("Number of atoms along the y-axis of BOTTOM plate: %d\n", n1y);
	printf("\n");

	printf("Size of the LJ_sigma parameter in angstroms (PLATES' ATOMS): %.3lf\n", sigmaP);
	printf("\n");

	printf("Absolute value of the charge: %.3lf\n", charge);
	printf("\n");

	printf("Mass of PLATES' ATOMS (mP): %.3lf\n", mP);
	printf("\n");

	printf("Epsilon parameter epsPP: %.3lf\n", epsPP);
	printf("\n");

	printf("Bonf coeff, K: %.3lf\n", K);
	printf("\n");


/*************************************************************************************************************************************/

	/* Write the values of all input parameters into a .txt file */

    sprintf(filename, "./InputParameters.txt");
	fp = fopen(filename, "w");

	fprintf(fp, "Number of atoms along the y-axis of BOTTOM plate: %d\n", n1y);
	fprintf(fp, "\n");

	fprintf(fp, "Size of the LJ_sigma parameter in angstroms (PLATES' ATOMS): %.3lf\n", sigmaP);
	fprintf(fp, "\n");

	printf("Absolute value of the charge: %.3lf\n", charge);
	printf("\n");

	fprintf(fp, "Mass of PLATES' ATOMS (mP): %.3lf\n", mP);
	fprintf(fp, "\n");

	fprintf(fp, "Number of groups to which each plate is divided, PD: %d\n", PD);
	fprintf(fp, "\n");

	for (i=0; i<PD; i++){
		fprintf(fp, "Number of layers in %d group of a plate: %d\n", i+1, nLPD[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Bond coeff K: %.3lf\n", K);
	fprintf(fp, "\n");

	fclose(fp);


/****************************************************************************************************************************************/

	D = 1.5*sigmaP; // lattice constant of Cu plate

	/* FCC_111 LATTICE - calculation of the total number of layers in both plates */
	nL = 0;
	for (i=0; i<PD; i++){
		nL += nLPD[i];
	}
	printf("Total number of layers in the plates: %d\n", nL);


	// BOTTOM plate => struct atom r1[N1]
	n1x = n1y;
	n1z = floor(nL/3.0);
	N1 = 6*n1x*n1y*n1z; // total number of atoms in the BOTTOM plate
	struct atom r1[N1]; // BOTTOM plate atoms

	MakeFCC111(n1y, nL, D, LL, r1);

	Lx = LL[0];
	Ly = LL[1];
	Lz = LL[2];
	printf("Lengths of the BOTTOM plate (Lx, Ly, Lz): %3.3lf %3.3lf %3.3lf\n", Lx, Ly, Lz);


	// Length of graphene
	double Ll = Ly;
	double Dd = Ll/6;
	double Hh = Dd/3;
	double sgraph = grapheneLength(Ly)/Dd;

	printf("Length of the curved part of graphene s = %10.10lf\n", sgraph);
	double graphenL = sgraph + (Ll - 2*Dd);
	printf("Length of graphene graphenL = %10.10lf\n", graphenL);


	// Graphene => struct atom r3[N3]
	int nnx, nny, N3;
	double aa;
	double aax, aay;

	aa = D; // lattice constant of graphene
	aax = aa*sqrt(3);
	aay = 3*aa;

	nnx = round(Lx/aax);
	nny = round(Ly/aay);

	N3 = 4*nnx*nny; // total number of atoms in graphene
	struct atom r3[N3]; // graphene atoms

	MakeGraphene(Lx, Ly, aa, r3);

	struct atom r4[N4]; // hydrogen atoms
	for(i=0;i<N4;i++)
	{
		r4[i].x=0;
		r4[i].y=0;
		r4[i].z=0;
	}

	MakeHydrogen( N4,  r4, Lx, Ly, Lz);



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double delz = D*sqrt(3)/3;
	printf("Value of the delz parameter: %3.3lf\n", delz);


	Output(n1x, n1y, N1, nL, N3, D, sigmaP, charge, mP, epsPP, K, PD, nLPD, delz, Lx, Ly, r1, r3, N4, r4);

	return 0;
};

				              /* DEFINITIONS OF THE FUNCTIONS */
/****************************************************************************************************************************************/

void MakeFCC111(int n1y, int nL, double D, double LL[], struct atom r1[]){

	int i, j, k, layer, count;
	double ax, ay, az, x2, y2, y3, y4, y5, y6;
	double tx, ty, tz;
	double x1min, x1max, y1min, y1max, z1min, z1max; // borders of the BOTTOM plate
	int n1x, n1z, N1;


	// BOTTOM plate
	n1x = n1y;
	n1z = floor(nL/3.0);
	N1 = 6*n1x*n1y*n1z; // total number of atoms in the BOTTOM plate

	printf("Number of atoms in the BOTTOM plate: %d\n", N1);

	double arrayx1[N1];
	double arrayy1[N1];
	double arrayz1[N1];

	ax = D*sqrt(2)/2;
	ay = D*sqrt(6)/2;
	az = D*sqrt(3);

	x2 = D*sqrt(2)/4;
	y2 = D*sqrt(6)/4;
	y3 = D*sqrt(6)/6;
	y4 = D*sqrt(6)*5/12;
	y5 = D*sqrt(6)*2/6;
	y6 = D*sqrt(6)/12;


	/* Formation of the BOTTOM plate */
	count = -1;
	for(i=0; i<n1x; i++){
    		for(j=0; j<n1y; j++){
      			layer = 0;
                	for(k=0; k<n1z; k++){

				count++;
                		tx = i*ax;
				ty = j*ay;
				tz = layer/3.0 * az;

				r1[count].x = tx;
				r1[count].y = ty;
				r1[count].z = tz;
				/************************************************************************/
				count++;
                		tx = x2 + i*ax;
				ty = y2 + j*ay;
				tz = layer/3.0 * az;

				r1[count].x = tx;
				r1[count].y = ty;
				r1[count].z = tz;

                		layer = layer + 1;
				/***********************************************************************/
				count++;
                		tx = i*ax;
				ty = y3 + j*ay;
				tz = layer/3.0 * az;

				r1[count].x = tx;
				r1[count].y = ty;
				r1[count].z = tz;
				/***********************************************************************/
      				count++;
                		tx = x2 + i*ax;
				ty = y4 + j*ay;
				tz = layer/3.0 * az;

				r1[count].x = tx;
				r1[count].y = ty;
				r1[count].z = tz;

                		layer = layer + 1;
				/***********************************************************************/
				count++;
                		tx = i*ax;
				ty = y5 + j*ay;
				tz = layer/3.0 * az;

				r1[count].x = tx;
				r1[count].y = ty;
				r1[count].z = tz;
				/***********************************************************************/
				count++;
                		tx = x2 + i*ax;
				ty = y6 + j*ay;
				tz = layer/3.0 * az;

				r1[count].x = tx;
				r1[count].y = ty;
				r1[count].z = tz;

                		layer = layer + 1;
				/************************************************************************/
			}
		}
	}

	printf("Formation of the BOTTOM plate is done!\n");

	for (i=0; i<N1; i++){
		arrayx1[i] = r1[i].x;
		arrayy1[i] = r1[i].y;
		arrayz1[i] = r1[i].z;
	}

	x1min = minAr(N1, arrayx1);
	x1max = maxAr(N1, arrayx1);
	y1min = minAr(N1, arrayy1);
	y1max = maxAr(N1, arrayy1);
	z1min = minAr(N1, arrayz1);
	z1max = maxAr(N1, arrayz1);
	printf("(x1min x1max y1min y1max z1min z1max) = %3.3lf %3.3lf %3.3lf %3.3lf %3.3lf %3.3lf\n", x1min, x1max, y1min, y1max, z1min, z1max);

	LL[0] = x1max - x1min;
	LL[1] = y1max - y1min;
	LL[2] = z1max - z1min;
	printf("Lengths of the BOTTOM plate (Lx, Ly, Lz): %3.3lf %3.3lf %3.3lf\n", LL[0], LL[1], LL[2]);


	return;
};


void MakeGraphene(double Lx, double Ly, double a, struct atom r3[]){

	int i, j, count;
	double ax, ay;
	double x1, y1, x2, y2, x3, y3, x4, y4;
	double tx, ty;
	int nx, ny;
	double D, H;
	const pi = 3.141592654;

	ax = a*sqrt(3);
	ay = 3*a;

	nx = round(Lx/ax);
	ny = round(Ly/ay);

	x1 = 0;
	y1 = 0;
	x2 = 0;
	y2 = a;
	x3 = a*sqrt(3)/2;
	y3 = 3*a/2;
	x4 = a*sqrt(3)/2;
	y4 = 5*a/2;

	D = Ly/3; // length parameter of graphene along the y-direction
	H = Ly/9; // height parameter of the graphene deformation along the z-direction

	/* Graphene formation using 4 points of the unit cell */
	count = -1;
	for (i=0; i<nx; i++){
		for (j=0; j<ny; j++){
			count++;
			tx = x1 + i*ax;
			ty = y1 + j*ay;

			r3[count].x = tx;
			r3[count].y = ty;
			if ((r3[count].y >= (Ly/2 - D)) && (r3[count].y <= (Ly/2 + D))){
				r3[count].z = a + H*(1 + cos((r3[count].y - Ly/2)*pi/D));
			}
			else r3[count].z = a;
	/*******************************************************************************************************/
			count++;
			tx = x2 + i*ax;
			ty = y2 + j*ay;

			r3[count].x = tx;
			r3[count].y = ty;
			if ((r3[count].y >= (Ly/2 - D)) && (r3[count].y <= (Ly/2 + D))){
				r3[count].z = a + H*(1 + cos((r3[count].y - Ly/2)*pi/D));
			}
			else r3[count].z = a;
	/******************************************************************************************************/
			count++;
			tx = x3 + i*ax;
			ty = y3 + j*ay;

			r3[count].x = tx;
			r3[count].y = ty;
			if ((r3[count].y >= (Ly/2 - D)) && (r3[count].y <= (Ly/2 + D))){
				r3[count].z = a + H*(1 + cos((r3[count].y - Ly/2)*pi/D));
			}
			else r3[count].z = a;
	/******************************************************************************************************/
			count++;
			tx = x4 + i*ax;
			ty = y4 + j*ay;

			r3[count].x = tx;
			r3[count].y = ty;
			if ((r3[count].y >= (Ly/2 - D)) && (r3[count].y <= (Ly/2 + D))){
				r3[count].z = a + H*(1 + cos((r3[count].y - Ly/2)*pi/D));
			}
			else r3[count].z = a;
		}
	}


	return;
};

double rand2(double a, double b)
{
	return  a +(b-a)*((double)rand()/RAND_MAX);
}


double maxAr(int len, double array[]){

	int i;
	double out;

	out = -10000000;
	for (i=0; i<len; i++){
		if (array[i] > out){
			 out = array[i];
		}
	}

	return out;
};


double minAr(int len, double array[]){

	int i;
	double out;

	out = 10000000;
	for (i=0; i<len; i++){
		if (array[i] < out){
			out = array[i];
		}
	}

	return out;
};


double fmax(double a, double b){

	if (a >= b){
	   return a;
	}
	else {
	   return b;
	}
}


double fmin(double a, double b){

	if (a <= b){
	   return a;
	}
	else {
	   return b;
	}
}


double dist(double r1[], double r2[]){

  	double d;
  	int i;
  	double dr[3];

  	d = 0.0;
  	for (i = 0; i<3; i++){
    		dr[i] = r1[i] - r2[i];
    		d = d + dr[i] * dr[i];
  	}
  	d = sqrt(d);

  	return d;
};


double grapheneLength(double Ly){

	int i;
	int np = 1000;
	double L, D, H, dy, s;
	double y[np];
	const pi = 3.141592654;

	L = Ly;
	D = L/4;
	H = D/3;

	printf("L = %3.3lf, D = %3.3lf, H = %3.3lf\n", L, D, H);

	for (i=0; i<np; i++){
		y[i] = -D + (i*(2*D))/np;
		/* if (i % 100 == 0){
			printf("In step %d y[i] = %3.3lf\n", i, y[i]);
		} */
	}

	dy = y[1] - y[0];

	s = 0;
	for (i=0; i<np; i++){
		s += sqrt(1 + (H*H*pi*pi*sin(pi*y[i]/D)*sin(pi*y[i]/D))/(D*D))*dy;
		/* if (i % 100 == 0){
			printf("In step %d s = %3.3lf\n", i, s);
		} */
	}

	return s;
};


void Output(int n1x, int n1y, int N1, int nL, int N3, double D, double sigmaP, double charge, double mP, double epsPP, double K, int PD, int nLPD[], double delz, double Lx, double Ly, struct atom r1[], struct atom r3[], int N4, struct atom r4[]){

	int i, j, k, ind, count1; // counters

	int Ntypes = 4; // number of different types of atoms (BOTTOM plate)
	int mol_tag; // molecule_tag of the atoms
	int atom_type; // atom_type of the atoms

	int bond_type; // bond_type of the corresponding bond
	int BONDtypes = 1; // there is only one bond_type which is the same for every atom - ghost_atom pair
	int Nbonds = N1/nL; // total number of bonds: it is equal to the number of atom - ghost_atom pairs

	double KB[BONDtypes]; // array which contains the K coefficient of the bond types
	double r0[BONDtypes]; // array which contains the r0 coefficient of the bond types

	double q[Ntypes]; // array which contains the charge of the atom types

	double xlow, xhigh, ylow, yhigh, zlow, zhigh;
	double xc, yc, zc, xmin, xmax, ymin, ymax, zmin, zmax, S;
	double arrayx[N1], arrayy[N1], arrayz[N1], arrayz1[N1];
	double MM; // current value of the mass needed for the output to PARM.FILE


	struct atom r[N1];
	struct atom rtot[N1];
	struct Group groups1[PD];

	FILE *fp;
	FILE *ffB;
	char filename[64];


	// Setting the charge of each atom type and its index
	for (i=0; i<Ntypes; i++){
		q[i] = 0.0; // BOTTOM plate
	}

	// Seting the bond_coeff for every bond type
	for (i=0; i<BONDtypes; i++){
		KB[i] = K;
		r0[i] = sigmaP;
	}

	// Unification of the coordinates of two plates r1[N1] and r2t[N2t] into one array r[N1 + N2t]
	for (i=0; i<N1; i++){
		r[i].x = r1[i].x;
		r[i].y = r1[i].y;
		r[i].z = r1[i].z;
	}

	for (i=0; i<N1; i++){
		arrayx[i] = r[i].x;
		arrayy[i] = r[i].y;
		arrayz[i] = r[i].z;
	}

	// Set the simulation box dimensions according to the PBC (PERIODIC: x, y; NONPERIODIC: z)
	// In this folder periodic directions are x and y, while z is a nonperiodic direction

	xmin = minAr(N1, arrayx);
	xmax = maxAr(N1, arrayx);

	ymin = minAr(N1, arrayy);
	ymax = maxAr(N1, arrayy);

	zmin = minAr(N1, arrayz);
	zmax = maxAr(N1, arrayz);

	S = (xmax - xmin + D)*(ymax - ymin + D/2); // Calculate the cross-section area of the plates

	sprintf(filename, "./Area_of_Plates");
	fp = fopen(filename, "w");
	fprintf(fp, "DIMENSIONS OF THE PLATES\n\n\n");
	fprintf(fp, "(xmin, xmax, ymin, ymax, zmin, zmax) = (%3.3lf, %3.3lf, %3.3lf, %3.3lf, %3.3lf, %3.3lf)\n\n", xmin, xmax, ymin, ymax, zmin, zmax);
	fprintf(fp, "Cross-section area of the plates is: S = %3.3lf\n", S);
	fclose(fp);


	xlow = xmin; xhigh = n1x*D*sqrt(2)/2;
	ylow = ymin; yhigh = n1y*D*sqrt(6)/2;
	zlow = (zmin - 10*D); zhigh = (zmax + 25*D);


	/********************************************************************************************************************************/

	// Unification of (N1) atoms in the BOTTOM plate with other atoms
	// Total number of atoms in the system: Ntot = 2*N1 + N3

	for (i=0; i<N1; i++){
		rtot[i].x = r1[i].x;
		rtot[i].y = r1[i].y;
		rtot[i].z = r1[i].z;
	}


/****************************************************************************************************************************************/

	/* Output to the LAMMPS file */
	sprintf(filename, "./DATA.FILE");
	fp = fopen(filename, "w");

	fprintf(fp, "LAMMPS description\n");
	fprintf(fp, "\n");
	fprintf(fp, "%d		atoms\n", (2*N1/nL + N3 + N4));
	fprintf(fp, "\n");
	fprintf(fp, "%d		bonds\n", Nbonds);
	fprintf(fp, "\n");
	fprintf(fp, "%d		atom types\n", Ntypes);
	fprintf(fp, "%d		bond types\n", BONDtypes);
	fprintf(fp, "\n");
	fprintf(fp, "%3.3lf %3.3lf xlo xhi\n", xlow, xhigh);
	fprintf(fp, "%3.3lf %3.3lf ylo yhi\n", ylow, yhigh);
	fprintf(fp, "%3.3lf %3.3lf zlo zhi\n", zlow, zhigh);
	fprintf(fp, "\n");
	fprintf(fp, "Atoms\n");
	fprintf(fp, "\n");

	/* Additional arrays in order to process the selection of atom_types along the z-axis */
	for (i=0; i<N1; i++){
		arrayz1[i] = r1[i].z;
	}

	/* Selection of atom_types */

	/* Division of BOTTOM plate */
	groups1[0].B = minAr(N1, arrayz1);
	groups1[0].T = groups1[0].B + (nLPD[0] - 1)*delz;

	for (j=1; j<PD; j++){
		groups1[j].B = groups1[j - 1].T;
		groups1[j].T = groups1[j].B + (nLPD[j]*delz);
	}

	/* Test the division of BOTTOM plate */
	for (j=0; j<PD; j++){
		printf("groups1.B, groups1.T: %3.3lf, %3.3lf \n", groups1[j].B, groups1[j].T);
	}

	sprintf(filename, "./BOTTOM_atomtypes.txt");
	ffB = fopen(filename, "w");

	count1 = 0;

	/* Selection of atom groups along the z-axis */

	/* Multi layers in the plate */
	/*
	for (i=0; i<N1; i++){

		ind = i + 1;

		if (ind <= N1){
			for (j=0; j<PD; j++){
				if (j == 0){
					if (rtot[i].z >= groups1[j].B && rtot[i].z <= groups1[j].T){
						atom_type = j + 1;
						fprintf(ffB, "BOTTOM plate atom_type: %d\n", atom_type);
						count1++;
						mol_tag = atom_type;
					}
				}
				else if (j != 0){
					if (rtot[i].z > groups1[j].B && rtot[i].z <= groups1[j].T){
						atom_type = j + 1;
						fprintf(ffB, "BOTTOM plate atom_type: %d\n", atom_type);
						count1++;
						mol_tag = atom_type;
					}
				}
			}
		}

	fprintf(fp, "%d %d %d %3.3lf %3.3lf %3.3lf %3.3lf\n", ind, mol_tag, atom_type, q[atom_type - 1], rtot[i].x, rtot[i].y, rtot[i].z);
	}
	*/

	/* Keep just the lowest layer in the plate */

	ind = 0;
	for (i=0; i<N1; i++){
		if (rtot[i].z == zmin){
			ind++;
			atom_type = 1;
			mol_tag = atom_type;
			fprintf(fp, "%d %d %d %3.3lf %3.3lf %3.3lf %3.3lf\n", ind, mol_tag, atom_type, q[atom_type - 1], rtot[i].x, rtot[i].y, rtot[i].z);
			fprintf(fp, "%d %d %d %3.3lf %3.3lf %3.3lf %3.3lf\n", ind+1, mol_tag+1, atom_type+1, q[atom_type - 1], rtot[i].x, rtot[i].y, rtot[i].z);
			ind++;
		}
	}

	for (i=0; i<N3; i++){
		ind++;
		atom_type = 3;
		mol_tag = atom_type;
		fprintf(fp, "%d %d %d %3.3lf %3.3lf %3.3lf %3.3lf\n", ind, mol_tag, atom_type, q[atom_type - 1], r3[i].x, r3[i].y, r3[i].z);
	}

	for (i=0; i<N4; i++){
		ind++;
		atom_type = 4;
		mol_tag = atom_type;
		fprintf(fp, "%d %d %d %3.3lf %3.3lf %3.3lf %3.3lf\n", ind, mol_tag, atom_type, q[atom_type - 1], r4[i].x, r4[i].y, r4[i].z);
	}


	fprintf(ffB, "Number of atoms in BOTTOM plate is: %d\n", ind);
	fclose(ffB);


	fprintf(fp, "\n");
	fprintf(fp, "Bonds\n");
	fprintf(fp, "\n");

	for (i=0; i<Nbonds; i++){
		ind = i + 1;
		bond_type = 1;
	fprintf(fp, "%d %d %d %d\n", ind, bond_type, 2*i + 1, 2*i + 2);
	}

	fclose(fp);

/*****************************************************************************************************************************/

	/* Output to the PARM.FILE */

	/*
	sprintf(filename, "./PARM.FILE");
	fp = fopen(filename, "w");

	fprintf(fp, "# THis is the parameter file included in LAMMPS script\n");
	fprintf(fp, "# It contains information about the pair_style and pair_coeff\n");
	fprintf(fp, "# of all the interactions in the system\n\n\n");

	fprintf(fp, "pair_style		lj/cut/coul/msm 	%3.3lf\n\n\n", 3*fmax(fmax(sigmaC, sigmaA), sigmaP));

	for (i=0; i<Ntypes; i++){

		ind = i + 1;

		if (ind == (2*PD + 1)){
			MM = mC;
		}
		else if (ind == (2*PD + 2)){
			MM = mA;
		}
		else {
			MM = mP;
		}

		fprintf(fp, "mass	%d		%3.3lf\n", ind, MM);
	}

	fprintf(fp, "\n\n\n");

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", 1, PD, 1, PD, epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", (PD + 1), (2*PD), (PD + 1), (2*PD), epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", 1, PD, (2*PD + 3), (2*PD + 4), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", (PD + 1), (2*PD), (2*PD + 5), (2*PD + 6), epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 3), (2*PD + 3), epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 4), (2*PD + 4), epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 5), (2*PD + 5), epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 6), (2*PD + 6), epsPP, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", 1, PD, (PD + 1), (2*PD), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", 1, PD, (2*PD + 5), (2*PD + 6), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", (PD + 1), (2*PD), (2*PD + 3), (2*PD + 4), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d*%d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 3), (2*PD + 4), (2*PD + 5), (2*PD + 6), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 3), (2*PD + 4), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 5), (2*PD + 6), epsBT, sigmaP, 3*sigmaP);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 1), (2*PD + 1), epsII, sigmaC, 3*sigmaC);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 2), (2*PD + 2), epsII, sigmaA, 3*sigmaA);

	fprintf(fp, "pair_coeff %d	%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 1), (2*PD + 2), epsII, (sigmaC + sigmaA)/2, 3*(sigmaC + sigmaA)/2);

	fprintf(fp, "pair_coeff %d*%d	%d	%3.3lf	%3.3lf	%3.3lf\n", 1, (2*PD), (2*PD + 1), epsIP, (sigmaC + sigmaP)/2, 3*(sigmaC + sigmaP)/2);

	fprintf(fp, "pair_coeff %d*%d	%d	%3.3lf	%3.3lf	%3.3lf\n", 1, (2*PD), (2*PD + 2), epsIP, (sigmaA + sigmaP)/2, 3*(sigmaA + sigmaP)/2);

	fprintf(fp, "pair_coeff %d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 1), (2*PD + 3), (2*PD + 6), epsIP, (sigmaC + sigmaP)/2, 3*(sigmaC + sigmaP)/2);

	fprintf(fp, "pair_coeff %d	%d*%d	%3.3lf	%3.3lf	%3.3lf\n", (2*PD + 2), (2*PD + 3), (2*PD + 6), epsIP, (sigmaA + sigmaP)/2, 3*(sigmaA + sigmaP)/2);


	fclose(fp);
	*/

/*********************************************************************************************/



	return;
};


