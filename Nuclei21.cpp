#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <list>
#include<vector>
#include "subboxes.cpp"
#include "geompack.hpp"
#include "geompack.cpp"
#include "lodepng.h"
using namespace std;
bool boxes_avail=false;
#define FILENAME_LENGTH 100
#define LINELENGTH 10000
#define DELIMS " ,\t"
#define exponent_of_r 2.0
#define MAXDIM 3
#define floats double
#define twoPi 6.2831853071795864769
#define SIX 6.0

#define dx0_to_b(dx0, b, TYPE) \
  switch (TYPE) {             \
    case 1:                   \
      b = dx0; break;         \
    case 2:                   \
      b = dx0; break;         \
    case 3:                   \
      b = dx0; break;     \
    case 4:                   \
      b = dx0; break; \
    default:                  \
      b = dx0; break;         \
}


#define CALC_dX0(XX, iD, jD, dX)   \
      for(nu=0; nu<D; nu++) { \
        dX[nu] = XX[iD + nu] - XX[jD + nu]; \
      } 


#define CALC_dV0(VV, iD, jD, dV)   \
      for(nu=0; nu<D; nu++) { \
        dV[nu] = VV[iD + nu] - VV[jD + nu]; \
      } 

#define ADD_COMPONENTS(ADD_WHAT)          \
      for(j=0; j<Wall_AtomCount0; j++) {  \
        k  = Wall_AtomIndex0[j];          \
        kD = k*D;                         \
        for(nu=0; nu<pseudoD; nu++) {     \
          kDpnu = kD+nu;                  \
          ADD_WHAT                        \
        }                                 \
      }                                   \





#define ADD_BOUNDARYFORCE(FORCE_TYPE, bi, nu, Fnu) \
    if(BoundaryConnQ[bi] != 0) {  	   \
      if(BoundaryConnQ[bi] > 0)   	   \
	FORCE_TYPE[bi][nu] -= Fnu; 	   \
      else				   \
	FORCE_TYPE[bi][nu] += Fnu; 	   \
    }





struct BoundaryOsciStruct {
  floats time;
  floats amplitude[MAXDIM][MAXDIM];
  floats frequency[MAXDIM][MAXDIM];
  floats timeShift[MAXDIM][MAXDIM];
};

floats dt,dt_h;
floats SQRT3DT;
floats GAMMA,GAMMAi;
floats SIGMA,SIGMAoGAMMA;
long N,D,ND;
long pseudoD = -1;
floats L[MAXDIM];
floats Lh[MAXDIM];
long NSubBoxes[MAXDIM];

floats vBoundary[MAXDIM][MAXDIM]; // current direction of the boundaries
BoundaryOsciStruct BoundaryOsci;  // current oscillation velocity of the boundaries

floats L_LE[MAXDIM][MAXDIM]; 
floats L_LEinv[MAXDIM][MAXDIM];
int LeesEdwardsTrue;
int BoundaryOsciTrue;
int OverdampedTrue = 0;
int PeriodicBoundaryTrue[MAXDIM];
floats fmod_OffsetHalf(floats x, floats y);

struct interaction {
  int i1;      // first particle of the interaction
  int i2;      // second particle of the interaction
  short type;  // interaction type
  short damp;  // damping type ( 0 for no damping )
  floats dx0;   // the actual rest length
  floats dxbb1; // bond breaks for dx < dxbb1
  floats dxbb2; // bond breaks for dx > dxbb2
  floats dx2bb1,dx2bb2; // square of the above
  floats a;     // the energy scale
  floats b;     // typically sets the rest length
  floats c;     // the damping coefficient. Only considered for interactions with damp>0.
  floats c1,c2; // damping factors: 
                // [damping force] = force_damp(c,v) * (c1*[total velocity] + c2*[normal velocity])
  int life_time;//added by Lutz, to survey middle life time of Interactions
  double work;//added by Lutz, keep track of work done by connection
};

struct connectionCutStruct {
  long i1;      // first particle of the interaction
  long i2;      // second particle of the interaction
  short cutQ;   // cut connection ?
};

struct WallMoveStruct {
      floats time;
      int WallNumber;
      floats direction[MAXDIM];
};

struct BoundaryMoveStruct {
  floats time;
  floats vel[MAXDIM][MAXDIM];
};



void CALC_Xinv(floats *XX, floats *XXinv) {
  static int nu,bi;
  static floats dXinv;

  if (LeesEdwardsTrue) {
    for (bi=0; bi<D; bi++) {
      dXinv = 0.;
      for (nu=0; nu<D; nu++)
        dXinv += L_LEinv[bi][nu] * XX[nu];
      
      XXinv[bi] = dXinv - floor(dXinv);
    }
  }
  else {
    for(nu=0; nu<D; nu++) {
      if (L[nu]>0.) {
	dXinv = XX[nu]/L[nu]; 
	XXinv[nu] = dXinv - floor(dXinv);
      }
      else {
	XXinv[nu] = XX[nu];
      }
    }
  }
}

void CALC_dX(floats *dX, long *PeriodicCount) {
  static int nu,bi;
  static floats dXinv;
  
  if (LeesEdwardsTrue) {
    for (bi=0; bi<D; bi++) {
      dXinv = 0.;
      for (nu=0; nu<D; nu++) {
        dXinv += L_LEinv[bi][nu] * dX[nu];
      }
      PeriodicCount[bi] = (long)floor(dXinv + .5);
    }
    
    for (bi=0; bi<D; bi++) {
      for (nu=0; nu<D; nu++) {
        dX[bi] -= L_LE[bi][nu] * PeriodicCount[nu]; 
      }
    }
  }
  else {
    for(nu=0; nu<D; nu++) {
      if (L[nu]>0.) {
	dXinv = floor(dX[nu]/L[nu] + .5);
	dX[nu] = dX[nu] - L[nu]*dXinv;  // = fmod_OffsetHalf(dX[nu], L[nu]); 
	PeriodicCount[nu] = (long)dXinv;
      }
      else
	PeriodicCount[nu] = 0;
    }
  }
}


void CALC_dV(floats *dV, long *PeriodicCount) {
  static int nu,bi;
  
  if (LeesEdwardsTrue) {
    for (bi=0; bi<D; bi++) {
      for (nu=0; nu<D; nu++) {
        dV[nu] -= vBoundary[bi][nu] * PeriodicCount[nu]; 
      }
    }
  }
  else {
    for(nu=0; nu<D; nu++) 
      if (L[nu]>0.) 
        dV[nu] -= vBoundary[nu][nu] * PeriodicCount[nu]; 
  }
}

floats SQR(floats value) {
  return value*value;
}

floats CALC_SQR_VEC(floats *VEC) {
  static floats VEC2;
  static int nu;
  
  VEC2 = VEC[0]*VEC[0];
  for(nu=1;nu<D;nu++)
    VEC2 += VEC[nu]*VEC[nu];
  
  return VEC2;
}



int Inverse (double aIn[][MAXDIM], double ainv[][MAXDIM], int n) {
  int   i, j;                    // Zeile, Spalte
  int   s;                       // Elimininationsschritt
  int   pzeile;                  // Pivotzeile
  int   fehler = 0;              // Fehlerflag
  double f;                      // Multiplikationsfaktor
  const double Epsilon = 0.001;  // Genauigkeit
  double Maximum;                // Zeilenpivotisierung
  //extern FILE *fout;
  int pivot = 1;
  double a[MAXDIM][2*MAXDIM];
  

  // ergänze die Matrix a um eine Einheitsmatrix (rechts anhängen)
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      a[i][j] = aIn[i][j];

    for (j = 0; j < n; j++) {
      if (i==j) a[i][n+j] = 1.0;
      else      a[i][n+j] = 0.0;
    }
  }
#if DEBUG
  MatOut (stdout, a, n, 2*n);
#endif
  
  // die einzelnen Eliminationsschritte
  s = 0;
  do {
	  // Pivotisierung vermeidet unnötigen Abbruch bei einer Null in der Diagnonalen und
		// erhöht die Rechengenauigkeit
    Maximum = fabs(a[s][s]);
    if (pivot)
    {
      pzeile = s ; 
      for (i = s+1; i < n; i++)
        if (fabs(a[i][s]) > Maximum) {
        Maximum = fabs(a[i][s]) ;
        pzeile = i;
        }
    }
    fehler = (Maximum < Epsilon);

    if (fehler) break;           // nicht lösbar 

    if (pivot)
    {
      if (pzeile != s)  // falls erforderlich, Zeilen tauschen
      { double h;
      for (j = s ; j < 2*n; j++) {
        h = a[s][j];
        a[s][j] = a[pzeile][j];
        a[pzeile][j]= h;
      }
      }
    }

    // Eliminationszeile durch Pivot-Koeffizienten f = a[s][s] dividieren
    f = a[s][s];
    for (j = s; j < 2*n; j++)
      a[s][j] = a[s][j] / f;

    // Elimination --> Nullen in Spalte s oberhalb und unterhalb der Diagonalen
    // durch Addition der mit f multiplizierten Zeile s zur jeweiligen Zeile i
    for (i = 0; i < n; i++ ) {
      if (i != s) 
      {
        f = -a[i][s];                 // Multiplikationsfaktor
        for (j = s; j < 2*n ; j++)    // die einzelnen Spalten
          a[i][j] += f*a[s][j];       // Addition der Zeilen i, s
      }
    }
#if DEBUG
    fprintf(stdout, "Nach %1i-tem Eliminationschritt:\n", s+1);
    MatOut (stdout, a, n, 2*n);
#endif
    s++;
  } while ( s < n ) ;

  if (fehler) 
  {
    fprintf (stdout, "Inverse: Matrix ist singulär\n");
    return 0; 
  }
  // Die angehängte Einheitsmatrix Matrix hat sich jetzt in die inverse Matrix umgewandelt
  // Umkopieren auf die Zielmatrix
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      ainv[i][j] = a[i][n+j];
    }
  }
  return 1;  
}

void PrintMatrix(floats MAT[][MAXDIM]) {
      int bi,nu; 
      for (bi=0; bi<D; bi++) {
        cout << "| ";
        for (nu=0; nu<D; nu++)
          cout << MAT[bi][nu] << "\t";
        
        cout << "|" << endl;
      }
}

void PrintOsciMatrix(BoundaryOsciStruct OSCI) {
  int bi,nu; 
  for (bi=0; bi<D; bi++) {
    cout << "|";
    for (nu=0; nu<D; nu++) {
      cout << " " << OSCI.amplitude[bi][nu] << "\tcos( 2 pi " << OSCI.frequency[bi][nu] << "\t";
      if (OSCI.timeShift[bi][nu] == 0.)
        cout << " t";
      else
        cout << "( t - " << OSCI.timeShift[bi][nu] << " )";
      cout << " )\t";
    }
    cout << "|" << endl;
  }
}

floats fmodt(floats x, floats y) {
  return x - y*floor(x/y);
}

floats fmod_OffsetHalf(floats x, floats y) {
  return x - y*floor(x/y + .5);
}

/*int fmod_LeesEdwards(floats *DeltaX) {
  static floats dX_PeriodicCount[MAXDIM];
  static int nu,bi;
  
  for (nu=0; nu<D; nu++)
    dX_PeriodicCount[nu] = floor(DeltaX[nu]/L[nu] + .5);
  
  for (nu=0; nu<D; nu++) {
    DeltaX[nu] -= L[nu] * dX_PeriodicCount[nu];
    for (bi=0; bi<D; bi++) {
      if (nu != bi)
        DeltaX[nu] -= LEshift[bi][nu] * dX_PeriodicCount[bi];
    }
  }
  
  return 0;
}*/


int comp0(const void *a,const void *b) {  //compare the first matrix elements
  return round(( **(floats**)a - **(floats**)b ));
}

int compWM(const void *a,const void *b) {  
  return ( (*(WallMoveStruct*)a).time > (*(WallMoveStruct*)b).time );
}

int compBM(const void *a,const void *b) {  
  return ( (*(BoundaryMoveStruct*)a).time > (*(BoundaryMoveStruct*)b).time );
}

int compBO(const void *a,const void *b) {  
  return ( (*(BoundaryOsciStruct*)a).time > (*(BoundaryOsciStruct*)b).time );
}

void CheckParticleNumber(const char *POSFILE, long *N, long *D){
  int d,k;
  floats dummy;
  ifstream POS0(POSFILE);
    
  if(!POS0.is_open()) {
    cout << "Position file \"" << POSFILE << "\" essential, but not found." << endl;
    exit(42);
  }
  
    *D=0;
    while(POS0.peek()!=13 && POS0.peek()!=10 && !POS0.eof()) {   // if next char is newline "\n" (ASCII 13) or carriage return "\r" (ASCII 10)
      POS0 >> dummy;
      (*D)++;
      //cout << dummy << "  " << *D << endl;
    }

    k=1;
    while(!POS0.eof()){
      POS0 >> dummy;
      d=1;
      while(POS0.peek()!=13 && POS0.peek()!=10 && !POS0.eof()) {
	POS0 >> dummy;
	d++;
      }
      if (d==1)
	break;
      k++;
      if (d!=(*D)) {printf("Dimension not consistent! Line %d in %s\n",k,POSFILE);exit(1);}
    }
    POS0.close();
    *N=k;
}

long CheckNumberOfLines(const char *FILE) {
  long k=0;
  char dummy[1000];
  ifstream stream(FILE);
  if (stream.is_open()) {
    while(!stream.eof()) {
      stream.getline(dummy,1000);
      //if (dummy!="")
      if ( strcmp(dummy, "") != 0)
        k++;
    }
    stream.close();
    return k;
  }
  else
    return -1;
}


string addNDot(string st, int k) {
  int kk = st.rfind(".",st.length());
  if(kk==(int)string::npos) 
    kk=st.length();
  
  char StHelp[8];
  sprintf(StHelp,"%d",k);
  return st.insert(kk,StHelp);
}

floats randomKick() {
  floats rand = drand48();
  if(rand <= .16666666667) return -SQRT3DT;
  if(rand >  .83333333333) return SQRT3DT;
  return 0.;
  
  /* *** this would be a gaussian: ***
  double z1, z2, w, y1, y2;
  do {
    z1 = 2.0 * drand48() - 1.0;
    z2 = 2.0 * drand48() - 1.0;
    w = z1 * z1 + z2 * z2;
  } while ( w >= 1.0 );
  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = z1 * w;
  //y2 = z2 * w;
  
  return(y1);
  */
}


  
floats force(short type, floats x2, floats a, floats b) {
  /** 
  * This function "force" can be used for central forces. 
  * The input parameter x2 must be x^2
  * For a potential V(x) it must return V'(x)/x
  * Hence, it is made such that the actual force vector is:   
  * [Fx Fy] = force(...) [x, y]       or   [Fx Fy Fz] = force(...) [x, y, z]
  * hence, case 1, a "constant" force(...) would correspond to a Hookian spring 
  *
  * special cases: 
  * V(x)       = b/n * x^n
  * force(...) = b   * x^(n-2) = b * x2^(n/2 - 1)
  * 
  * V(x)       = b/n * (x-x0)^n
  * force(...) = b   * (x-x0)^(n-1) / x
  **/
  switch (type)
  {
    case 1://stepwise
      if(sqrt(x2) > b/2.0)
        return a/2.0;
      else
        return a;
      break;
    case 2://fast devaying
      return a*pow(b/(2.0*sqrt(x2)),3.0);
      break;
    case 3://slow decaying
      //theta force is older now we do this semi soft hard core. The name
      //sounds
      //extremely stupid
      if(sqrt(x2) <= b/2.0)
      {
          return a*pow(b,3)/(8.0*x2);
          break;
      }
      if(sqrt(x2) <= b)
      {
          return a*(b-sqrt(x2));
          break;
      }
      return 0;
      break;
      //this last case must not be reached if everything works well
    case 4://constant
      return a;//now we have a theta like force
      break;
      //return a*exp(-sqrt(x2)/b);
    case 5:
      return 0;
    default: 
      cout << "invalid force type: " << type << endl;
      //return 0;
      exit(1);
      break;
  }
}

floats Fcrit_to_dxbb1(short type, floats a, floats dx0, floats Fcrit) {
  //static floats x2i,x6i;
  if (Fcrit < 0.)
    return -1.;
  else {
    switch (type) {
      case 1:
        return 0.; // Hookian spring
        break;
      case 2: 
        return 0.; // V(x) = a * ln(x)
        break;
      case 3: 
        return 0.;
        break;
      case 4:
        //x2i = 1. / x2;
        //x6i = x2i * x2i * x2i;        // Lenard Jones 
        return 0.;   // hence x = 1/b^(1/6) is the equilibrium distance
        break;
      case 5:
        return dx0 - Fcrit/(2.*a); // Hookian Spring with rest length dx0
        break;
  
      default: 
        cout << "invalid force type: " << type << endl;
        //return 0;
        exit(1);
        break;
    }
  }
}

floats Fcrit_to_dxbb2(short type, floats a, floats dx0, floats Fcrit) {
  //static floats x2i,x6i;
  if (Fcrit < 0.)
    return -1.;
  else {
    switch (type) {
      case 1:
        return 0.; // Hookian spring
        break;
      case 2: 
        return 0.; // V(x) = a * ln(x)
        break;
      case 3: 
        return 0.;
        break;
      case 4:
        //x2i = 1. / x2;
        //x6i = x2i * x2i * x2i;        // Lenard Jones 
        return 0.;   // hence x = 1/b^(1/6) is the equilibrium distance
        break;
      case 5:
        return dx0 + Fcrit/(2.*a); // Hookian Spring with rest length dx0
        break;
  
      default: 
        cout << "invalid force type: " << type << endl;
        //return 0;
        exit(1);
        break;
    }
  }
}


floats force_damp(short type, floats v2, floats c) {
  /** 
  * This function "force" can be used for central forces. 
  * The input parameter x2 must be x^2
  * For a potential V(x) it must return V'(x)/x
  * Hence, it is made such that the actual force vector is:   
  * [Fx Fy] = force(...) [x, y]       or   [Fx Fy Fz] = force(...) [x, y, z]
  * hence, case 1, a "constant" force(...) would correspond to a Hookian spring 
  *
  * special cases: 
  * V(x)       = b/n * x^n
  * force(...) = b   * x^(n-2) = b * x2^(n/2 - 1)
  * 
  * V(x)       = b/n * (x-x0)^n
  * force(...) = b   * (x-x0)^(n-1) / x
  **/
  switch (type)
  {
    case 1:
      return -c; // corresponds to F_damp(v) = c*v;
      break;

    default: 
      cout << "invalid damping force type: " << type << endl;
      exit(1);
      break;
  }
}
list<int> get_neigbours(int cell,box* map,double* X,double lc,int D,double* L);
double angle_to_uv(double* vec,int D);
int get_nmnpeh(int n)
{
    return (n+1)*(n+2)/2;
}
double get_shielded_force(double f,int i,int j,double* X,box* map,double lc,int D,double* L,double cell_radius,int* shield_matrix,double kappa)
{
    //just for konstant force to get a shielded output. Connections should be
    //reduced to neighbourhood
    int i1,i2;
    if(i>j)
    {
        i1=j;
        i2=i;
    }
    else
    {
        i1=i;
        i2=j;
    }
    if(boxes_avail)
    {
    if (shield_matrix[N*i1+i2-get_nmnpeh(i1)]==0)
    {
    double dX[D],pdX[D],pdX2,h,dX2;
    long dX_PeriodicCount[D];
    int nu;
    list<int> neigbours;
    neigbours=get_neigbours(i1,(box*)map,X,lc,D,L);
    list<int>::iterator neig;
	CALC_dX0(X,i1*D,i2*D,dX);
    CALC_dX(dX, dX_PeriodicCount);
    dX2=sqrt(CALC_SQR_VEC(dX));
    for (neig=neigbours.begin();neig!=neigbours.end();++neig)
    {
        if((*neig)!=i2)
        {
	    CALC_dX0(X,i1*D,(*neig)*D,pdX);
        CALC_dX(pdX, dX_PeriodicCount);
        pdX2=sqrt(CALC_SQR_VEC(pdX));
        if(pdX2<dX2)
        {
        h=abs(pdX2*sin(angle_to_uv(dX,D)-angle_to_uv(pdX,D)))/2.0;
        if (h<cell_radius)
        {
            neigbours.clear();
            shield_matrix[N*i1+i2-get_nmnpeh(i1)]=1;
            return kappa*f;
        }
        }
        }
    }
    neigbours.clear();
    shield_matrix[N*i1+i2-get_nmnpeh(i1)]=-1;
    return f;
    }
    if(shield_matrix[N*i1+i2-get_nmnpeh(i1)]==1)
        return kappa*f;
    if(shield_matrix[N*i1+i2-get_nmnpeh(i1)]==-1)
    {
        return f;
    }
    cout << "undefined value " << shield_matrix[N*i1+i2-get_nmnpeh(i1)] << " at " << i1 << " " << i2 << " " << N*i1+i2-get_nmnpeh(i1) << " in shield matrix\n"; 
    return 0;
    }
    return f;
}
double cut_off_power(double in)
{
    int out=round(pow(in,2));
    out=round(pow(out,2));
    return (double)out;
}
void calc_F(floats *X, floats *V, long Wall_AtomCount0, long *Wall_AtomIndex0, list<interaction> connections, floats *F,box* map,double lc,double cell_radius,int* shield_matrix,double shielding,double box_size,double cellr2) {
  long k,nu;
  list<int>::iterator int_it;
  int i;
  double excl_pow;
  long conn_i1D, conn_i2D;
  //long kD;
  long  dX_PeriodicCount[D];
  short conn_damp;
  list<int> tempN;
  floats dX[D];
  floats dV[D];
  floats dX2,dV2;
  floats Fh,Fnu,c1Fh,c2Fh,dXdV;
  
  interaction conn;
  list<interaction>::iterator it;
  
  // initialize all forces
  for(k=0; k<ND; k++)
    F[k]=0.;
  
  
  //long j=0;
  for(it = connections.begin(); it != connections.end(); ++it) {
    conn = *it;
    conn_i1D = conn.i1*D;
    conn_i2D = conn.i2*D;
  
    CALC_dX0(X, conn_i1D, conn_i2D, dX);
    CALC_dX(dX, dX_PeriodicCount);
    dX2 = CALC_SQR_VEC(dX);
  
    /*if(dX2>2500.) {
      cout << j << " ";
      cout << conn.i1 << ": [" << X[conn_i1D] << ", " << X[conn_i1D+1] << "] ";
      cout << conn.i2 << ": [" << X[conn_i2D] << ", " << X[conn_i2D+1] << "] ";
      cout << " [" << dX2 << "]" << endl;
    }
      j++;*/
  
    // calculate contribution of central force
    Fh = get_shielded_force(force(conn.type, dX2, conn.a, conn.b),conn.i1,conn.i2,X,map,lc,D,L,cell_radius,shield_matrix,shielding);
    //Fh = force(conn.type, dX2, conn.a, conn.b);
    for(nu=0; nu<D; nu++) {
      Fnu = dX[nu]*Fh;
      F[conn_i1D + nu] += Fnu;
      F[conn_i2D + nu] -= Fnu;
    }
  
    // calculate contribution of damping force
    conn_damp = conn.damp;
    if(conn_damp > 0) {
      CALC_dV0(V, conn_i1D, conn_i2D, dV);
      CALC_dV(dV, dX_PeriodicCount);
  
      if(conn_damp == 1) {
        Fh = - conn.c;
      }
      else {
        dV2 = CALC_SQR_VEC(dV);
        Fh = force_damp(conn_damp, dV2, conn.c);
      }
  
      if (conn.c1 != 0.) {
        c1Fh = conn.c1*Fh;
        for(nu=0; nu<D; nu++) {
          Fnu = dV[nu]*c1Fh;
          F[conn_i1D + nu] += Fnu;
          F[conn_i2D + nu] -= Fnu;
        }
      }
  
      if (conn.c2 != 0.) {
        dXdV = dX[0]*dV[0];
        for(nu=1; nu<D; nu++)
          dXdV += dX[nu]*dV[nu];
  
        c2Fh = conn.c2 * dXdV / dX2;
        for(nu=0; nu<D; nu++) {
          Fnu = dX[nu]*c2Fh;
          F[conn_i1D + nu] += Fnu;
          F[conn_i2D + nu] -= Fnu;
        }
      }
  
    }
  }
  //calculate Excluded Volume Forces
  if(boxes_avail)
  {
        for (int jota=0;jota<N;jota++)
        {
            tempN=get_neigbours(jota,map,X,box_size,D,L);
            for (int_it=tempN.begin();int_it!=tempN.end();++int_it)
            {
                i=(*int_it);
                if (jota!=i)
                {
                CALC_dX0(X,i*D,jota*D,dX);
                CALC_dX(dX, dX_PeriodicCount);
                dX2 = CALC_SQR_VEC(dX);
                //dX2 = sqrt(dX2)/2.0;
                //if(dX2 < cellr2)
                //{
                excl_pow=(pow(cellr2/dX2,3)*pow(cellr2/dX2,3))/sqrt(dX2);
                for (int erta=0;erta<D;erta++)
                {
                    F[i*D+erta]=F[i*D+erta]+dX[erta]*excl_pow;
                    F[jota*D+erta]=F[jota*D+erta]-dX[erta]*excl_pow;
                }
                //}
                }
            }
        }
  }
}

//Measurements, added by Lutz Kuenneke
list<int> mymerge(list<int> list1,list<int> list2,int cell)
{
    list<int>::iterator it,at;
    bool isin;
    for (at=list2.begin();at!=list2.end();++at)
    {
        isin=false;
        for (it=list1.begin();it!=list1.end();++it)
        {
            if (*it==*at)
                isin=true;
        }
        if(isin==false && *at != cell)
            list1.push_back(*at);
    }
    return list1;
}
string type_to_name(int type)
{
    //cosmetic function
    switch (type) {
        case 1:
            return "Stepwise decaying force\n";
            break;
        case 2:
            return "Fast decaying force\n";
            break;
        case 3:
            return "Slow decaying force\n";
            break;
        case 4:
            return "constant force\n";
            break;
        default :
            return "unknown force\n";
    }
}
double periodic_middle(double x1,double x2,double l)
{
    //this function shall find the middle of two points respecting periodic
    //boundary conditions
    while (x1<0)
        x1+=l;
    while (x1>=l)
        x1-=l;
    while (x2<0)
        x2+=l;
    while (x2>=l)
        x2-=l;
    double dist[3];
    double temp;
    dist[0]=abs(x1-x2);//both in one period
    dist[1]=abs(x1+l-x2);//x1 is one period ahead
    dist[2]=abs(x1-l-x2);//x2 is one period ahead
    if (dist[0]<=dist[1] && dist[0] <= dist[2])
    {
        return (x1+x2)/2.0;
    }
    else
    {
        if (dist[1]<dist[2])
        {
            temp=(x1+l+x2)/2.0;
            while (temp>l)
                temp-=l;
            while(temp<0)
                temp+=l;
            return temp;
        }
        else
        {  
            temp=(x1-l+x2)/2.0;
            while (temp>l)
                temp-=l;
            while(temp<0)
                temp+=l;
            return temp;
        }
    }
}
int plot_work_hist(int* histogramm,double dw,double max,double t)
{
    stringstream ss;
	ss<< "./run/OUT_work_hist" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
    double sum=0;
    for (int kappa=-ceil(max/dw);kappa<=ceil(max/dw);kappa++)
    {
        sum+=histogramm[kappa+(int)ceil(max/dw)]*dw;
    }
    if (sum>0)
    {
    for (int kappa=-ceil(max/dw);kappa<=ceil(max/dw);kappa++)
    {
        out << dw*(double)kappa << " " << histogramm[kappa+(int)ceil(max/dw)]/sum << endl;
    }
    }
    //reset histogramm now
    //for (int erta=-ceil(max/dw);erta<=ceil(max/dw);erta++)
    //    histogramm[erta+(int)ceil(max/dw)]=0;
    return 0;
}
int distribution_of_connections(double mean,boost::mt19937 rng)
{
    return mean;
    //return mean;
    //the purpose of this function is to get a distribution of connections
    //added at distance r around a given mean value.
    //Currently in use: Gaussian Distribution around mean n with deviation
    //sqrt(n)
    if (mean==0)
        return 0;
    boost::normal_distribution<> nd(mean, sqrt(mean));
    boost::variate_generator<boost::mt19937&, 
    boost::normal_distribution<> > var_nor(rng, nd);
    double d = var_nor();
    if(d<0)
        return 0;
    return d;//this should do the job
}
double mean_of_conn(double r,double lc,double eta,double kappa)
{
    return eta*(1.0/pow(r,kappa)-1.0/pow(lc,kappa));
}
void init_boxes(double* L,box* map,double lc,int D,double* X,double N)
{
    list<int> members;
    int boxes[D];
    for (int kappa=0;kappa<D;kappa++)
        boxes[kappa]=(int)ceil(L[kappa]/lc);
    cout << "-- creating " << boxes[0] << "x" << boxes[1] << " boxes on boxsize " << lc << " " << L[0] << " " << L[1] << "--\n";
    //box map[boxes[0]][boxes[1]];
    //support for higher dimensions not yet implemented
    for(int jota=0;jota<boxes[0];jota++)
    {
        for(int eta=0;eta<boxes[1];eta++)
        {
            //find members in corresponding area
            members.clear();
            for (int i=0;i<N;i++)
            {
                while(X[i*D] >= L[0])
                    X[i*D]-=L[0];
                if(X[i*D+1] >= L[1])
                    X[i*D+1]-=L[1];
                while(X[i*D]<0)
                    X[i*D]+=L[0];
                while(X[i*D+1]<0)
                    X[i*D+1]+=L[1];
                if ((int)floor(X[i*D]/lc) == jota && (int)floor(X[i*D+1]/lc) == eta)
                {
      //              cout << "adding member " << i << endl;
                    members.push_back(i);
                }
            }
            map[eta+boxes[1]*jota]=map[eta+boxes[1]*jota].new_box(members);
        }
    }
    boxes_avail=true;
}
list<int> get_neigbours(int cell,box* map,double* X,double lc,int D,double* L)
{
    list<int> neigbours;
    double x,y;
    x=X[cell*D];
    y=X[cell*D+1];
    while (x>=L[0])
        x-=L[0];
    while (x<0)
        x+=L[0];
    while (y>=L[1])
        y-=L[1];
    while (y<0)
        y+=L[1];
    neigbours.sort();
    list<int> spacer;
    int xbc=(int)floor(x/lc);
    int ybm;
    if(L[1]/lc == floor(L[1]/lc))
        ybm=(int)floor(L[1]/lc);
    else
        ybm=(int)floor(L[1]/lc)+1;
    int xbm;
    if(L[0]/lc == floor(L[0]/lc))
        xbm=(int)floor(L[0]/lc);
    else
        xbm=(int)floor(L[1]/lc)+1;
    int ybc=(int)floor(y/lc);
    int xb,yb;
    for(int jota=-1;jota<=1;jota++)
    {
        for(int eta=-1;eta<=1;eta++)
        {
            xb=xbc+jota;
            yb=ybc+eta;
            if(xb<0)
                xb=(int)floor(L[0]/lc);
            if(yb<0)
                yb=(int)floor(L[1]/lc);
            if(xb>=xbm)
                xb=0;
            if(yb>=ybm)
                yb=0;
           // cout << "requesting " << yb << "+" << ybm << "*" << xb << endl;
            spacer=map[yb+ybm*xb].get_members();
            //cout << "succes\n";
            spacer.sort();
            neigbours=mymerge(neigbours,spacer,cell);
            neigbours.sort();
            spacer.clear();
        }
    }
    return neigbours;
}
double get_total_energy(double* X,double *V,double N,int D,double rm,double lc,double f,double a,double r,double rho)//this is now the potential energy
{
    double ekin=0,epot=0,m=1,dX2,emin=1.0;
    double dX[D];
    long dX_PeriodicCount[D];
    int nu;
    //assume mass=1 for simplicicity
    for(int kappa=0;kappa<N*D;kappa++)
        ekin+=0.5*V[kappa]*V[kappa]*m;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(i!=j)
            {
	        CALC_dX0(X,i*D,j*D,dX);
            CALC_dX(dX, dX_PeriodicCount);
            dX2 = CALC_SQR_VEC(dX);
            dX2=sqrt(dX2);
            if(dX2<=rm)
                epot+=a/r*f*(dX2/(lc*lc)+1.0/lc);
            }
        }
    }
    return (ekin+epot)/(N*emin);
}
int get_life_time(double* life_time_distr,double t,int size,double step,double norm)
{
    stringstream ss;
	ss<< "./run/OUT_life_time" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
    for(int kappa=0;kappa<size;kappa++)
    {
        if(norm>0)
            out << step*(double)kappa+0.5*step << " " << life_time_distr[kappa]/norm << endl;
        else
            out << step*(double)kappa+0.5*step << " 0" << endl;
    }
    return 0;
}
int get_energy_distr(double* energy_distr,double* energy_normator,double t,int size,double step,double norm)
{
    stringstream ss;
    double enorm=0;
    for (int kappa=0;kappa<size;kappa++)
        enorm+=energy_distr[kappa]*step;
	ss<< "./run/OUT_energy_distr" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
    if (norm>0)
    {
    for(int kappa=0;kappa<size;kappa++)
    {
            out << step*(double)kappa+0.5*step << " " << energy_distr[kappa]<< endl;
    }
    }
    else
    {
        out << "0 0" << endl;
    }
    return 0;
}
double mod_by_per(double x,double l)
{
    while (x<0)
        x+=l;
    while(x>=l)
        x-=l;
    return x;
}
void to_period(double* X,double* L,int D,int N,int* map)
{
    int walker=N*D;
    for(int kappa=0;kappa<N;kappa++)
    {
        //setting up map and forcing inner points to boundary C
        map[kappa]=kappa;
	    for (int n=0;n<D;n++)
	    {
	    	while (X[kappa*D+n]>=L[n])
	    		X[kappa*D+n]=X[kappa*D+n]-L[n];
	    	while (X[kappa*D+n]<0)
	    		X[kappa*D+n]=X[kappa*D+n]+L[n];
	    }
    }
    //now add 3N more particels
    for (int kappa=0;kappa<N;kappa++)
    {
        for(int eta=0;eta<3;eta++)
            map[N+3*kappa+eta]=kappa;
        //every point gets three corresponding points outside the boundary
        if(X[kappa*D] > L[0]/2.0)
        {
            if(X[kappa*D+1] >= L[1]/2.0)
            {
                X[walker]=X[kappa*D];
                X[walker+1]=X[kappa*D+1]-L[1];
                walker+=2;
                X[walker]=X[kappa*D]-L[0];
                X[walker+1]=X[kappa*D+1];
                walker+=2;
                X[walker]=X[kappa*D]-L[0];
                X[walker+1]=X[kappa*D+1]-L[1];
                walker+=2;
            }
            if(X[kappa*D+1] < L[1]/2.0)
            {
                X[walker]=X[kappa*D];
                X[walker+1]=X[kappa*D+1]+L[1];
                walker+=2;
                X[walker]=X[kappa*D]-L[0];
                X[walker+1]=X[kappa*D+1];
                walker+=2;
                X[walker]=X[kappa*D]-L[0];
                X[walker+1]=X[kappa*D+1]+L[1];
                walker+=2;
            }
        }
        if(X[kappa*D] < L[0]/2.0)
        {
            if(X[kappa*D+1] >= L[1]/2.0)
            {
                X[walker]=X[kappa*D];
                X[walker+1]=X[kappa*D+1]-L[1];
                walker+=2;
                X[walker]=X[kappa*D]+L[0];
                X[walker+1]=X[kappa*D+1];
                walker+=2;
                X[walker]=X[kappa*D]+L[0];
                X[walker+1]=X[kappa*D+1]-L[1];
                walker+=2;
            }
            if(X[kappa*D+1] < L[1]/2.0)
            {
                X[walker]=X[kappa*D];
                X[walker+1]=X[kappa*D+1]+L[1];
                walker+=2;
                X[walker]=X[kappa*D]+L[0];
                X[walker+1]=X[kappa*D+1];
                walker+=2;
                X[walker]=X[kappa*D]+L[0];
                X[walker+1]=X[kappa*D+1]+L[1];
                walker+=2;
            }
        }
    }
    //at last perform dict ordering, means (x,y) <= (x',y') < = > x<x' or x=x'
    //and y<=y'
    //make next neigbours ordering, 4*N times comp in O(N^2)
    double x,y;
    int tmp;
    for(int tot=0;tot<4*N;tot++)
    {
    for (int ccount=0;ccount<4*N-1;ccount++)
    {
        if(X[ccount*D+D] < X[ccount*D])
        {
            x=X[ccount*D];
            X[ccount*D]=X[ccount*D+D];
            X[ccount*D+D]=x;
            y=X[ccount*D+1];
            X[ccount*D+1]=X[ccount*D+D+1];
            X[ccount*D+D+1]=y;
            tmp=map[ccount];
            map[ccount]=map[ccount+1];
            map[ccount+1]=tmp;
        }
        else { if (X[ccount*D+D] == X[ccount*D] && X[ccount*D+D+1] < X[ccount*D+1])
        {
            x=X[ccount*D];
            X[ccount*D]=X[ccount*D+D];
            X[ccount*D+D]=x;
            y=X[ccount*D+1];
            X[ccount*D+1]=X[ccount*D+D+1];
            X[ccount*D+D+1]=y;
            tmp=map[ccount];
            map[ccount]=map[ccount+1];
            map[ccount+1]=tmp;
        }
        }
    }
    }
}

void get_msd(double *msd_val,double *X,double *X0,int D,double N,double t,double* mean,double dt,double small_dt)
{
	//calculate mean
	int nu,kappa,jota;
	for (nu=0;nu<2;nu++)
		mean[nu]=0;
	for (kappa=0;kappa<N;kappa++)
	{
		for (jota=0;jota<D;jota++)
		{
			msd_val[kappa]+=pow(X[jota+kappa*D]-X0[jota+kappa*D],2)*dt;
		}
		mean[0]+=msd_val[kappa]/(t+small_dt);//in fact the timestep small_dt was made before measurement is performed
	}
	mean[0]=mean[0]/N;
	//calculate variance
	for (int kappa=0;kappa<N;kappa++)
	{
		mean[1]+=pow(msd_val[kappa]-mean[0],2);
	}
	mean[1]=sqrt(mean[0]/(N*(N-1)));
	//return mean;
}
int get_msd_distr(double *msd_val,double t,double N,double step)
{
    stringstream ss;
    double mean=0;
    //for (int kappa=0;kappa<N;kappa++)
    //{
    //    if(msd_val[kappa]>max)
    //        max=msd_val[kappa];
    //}
    //if (max==0)
    //{
    //    cout << "there is a weird max=0 in msd distr please check whats up\n";
    //    return 1;
    //}
	ss<< "./run/OUT_msd_dist" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
    int size=11,place;
    double arr[size];
    for (int kappa=0;kappa<size;kappa++)
        arr[kappa]=0;
    for (int kappa=0;kappa<N;kappa++)
        mean+=msd_val[kappa];
    mean=mean/N;
    for (int kappa=0;kappa<N;kappa++)
    {
        //center the msd distr on mean
        place=round(((double)size-1.0)*msd_val[kappa]/(2.0*mean));
        if (place<size)
            arr[place]+=1;
    }
    //write out
    for (int kappa=0;kappa<size;kappa++)
        out << 2.0/((double)size-1.0)*mean*(double)kappa << " " << arr[kappa]/N << endl;
    return 0;
}
//help function, absolute value
double myabs(double in)
{
	return sqrt(in*in);
}
//help function, get vector absolute value
double vabs(double *vec,int D)
{
	double vabs=0;
	for (int i=0;i<D;i++)
	{
		vabs+=pow(vec[i],2);
	}
	return sqrt(vabs);
}
//help function, get angle versus a certain unity vector

double angle_to_uv(double *vec,int D)
{
	//get angle to Vector v=(1,0,...,0)
//    double temp=acos(vec[0]/vabs(vec,D));
//    return temp;
    double p1[D],p2[D],p3[D];
    for(int kappa=0;kappa<D;kappa++)
    {
        p1[kappa]=0;//set vec, normally dX to origin
        p2[kappa]=vec[kappa];
        p3[kappa]=0;
    }
    p3[0]=vec[0];
    return angle_rad_2d(p2,p1,p3);//basically counterclockwise, see doc(geompack)
}
//help function, gaussian
double gaussian(double x0,double s,double p)//my stoneage gaussian, to be replaced some day, x0 mean, s deviation, higher p means higher computation time, and better gaussian distribution. p=100 is reasonable
{
       x0=x0-0.5;
       double gauss=0;
       int k=(int)p;
       if (k==0)
             k=1;
      for (int n=0;n<k;n++)
          gauss+=(drand48()-0.5)*s*sqrt(p)+0.5;
    return gauss/(p)+x0;
}
double get_op_simple(double *X,double *dX,int D,double N,double dX2Add,list<interaction> connections)
{
	list<interaction>::iterator it;
	double opam=0;
	double dX2,norm=0;//get norm by counting the number of summands
	int j,k,jD,kD,nu;
	//We assume particles nearer than dXAdd as neighboured
      for(it = connections.begin(); it != connections.end();++it ) 	{
        j = (*it).i1; jD=j*D;
        k = (*it).i2; kD=k*D;
	CALC_dX0(X,jD,kD,dX);
	dX2=dX[0]*dX[0];
	for(nu=1;nu<D;nu++)
		dX2+=dX[nu]*dX[nu];
		  opam+=cos(SIX*angle_to_uv(dX,D));
		  norm++;
          //cout << opam << " " << norm <<  " with " << j << " and " << k << endl;
	 }
      if (norm > 0)
	return myabs((double)opam/norm);
      else
	      return -1;
}
//more dedicated
double get_op(double *X,double *dX,int D,double N,double dX2Add,list<interaction> connections)
{
	list<interaction>::iterator it,at;
	complex<double> opam=0,im_unit(0,SIX);//im_unit is 6i
	double dX2,odX2,pdX[D],norm=0;
	int j,k,jD,kD,nu,lD,mD,l,m;
    long dX_PeriodicCount[D];
	//We assume particles nearer than dXAdd as neighboured
      for(it = connections.begin(); it != connections.end(); ++it) 	{
        j = (*it).i1; jD=j*D;
        k = (*it).i2; kD=k*D;
		 CALC_dX0(X,kD,jD,dX);
         CALC_dX(dX, dX_PeriodicCount);
		 dX2=dX[0]*dX[0];
		 pdX[0]=dX[0];
          for(nu=1;nu<D;nu++)
	  {
		  dX2+=dX[nu]*dX[nu];
		  pdX[nu]=dX[nu];
	  }
	
		  //in this case we are in the first pairs neighbourhood ant want to look for any other neighbourhoods
      for(at = connections.begin(); at != connections.end(); ++at) 	{
        l = (*at).i1; lD=l*D;
        m = (*at).i2; mD=m*D;
				CALC_dX0(X,mD,lD,dX);
                CALC_dX(dX, dX_PeriodicCount);
				odX2=dX[0]*dX[0];
				for(nu=1;nu<D;nu++)
					odX2+=dX[nu]*dX[nu];
				if (j==l)//look for local OP
				{
						// in this case we are in both neighbourhoods
						opam+=exp(im_unit*(angle_to_uv(pdX,D)-angle_to_uv(dX,D)));
						norm+=1.0;
				}
			  }
		  }
      if (norm>0)
	return abs(opam)/norm;
      else
	      return -1;
}
double get_sigma_over_mu(list<interaction> connections,double* X,int D,int N)
{
    //return deviation of mean distance over mean distance. Scalar measure for
    //spatial oredring
    //get mean
    list<interaction>::iterator it;
    int i,j,nu;
    double mean=0,dev=0,dX2,norm=0;
    long dX_PeriodicCount[D];
    double dX[D];
    for (it=connections.begin();it!=connections.end();++it)
    {
        i=it->i1;
        j=it->i2;
		 CALC_dX0(X,i*D,j*D,dX);
         CALC_dX(dX,dX_PeriodicCount);
		 dX2=dX[0]*dX[0];
         for (int erta=1;erta<D;erta++)
             dX2+=dX[erta]*dX[erta];
         dX2=sqrt(dX2);
         mean+=dX2;
         norm++;
    }
    mean=mean/norm;
    for (it=connections.begin();it!=connections.end();++it)
    {
        i=it->i1;
        j=it->i2;
		 CALC_dX0(X,i*D,j*D,dX);
         CALC_dX(dX,dX_PeriodicCount);
		 dX2=dX[0]*dX[0];
         for (int erta=1;erta<D;erta++)
             dX2+=dX[erta]*dX[erta];
         dX2=sqrt(dX2);
         dev+=pow(mean-dX2,2);
    }
    dev=sqrt(dev/(norm-1.0));
    return dev/norm;
}
//get the local op. Here we will have to make a time dependent output file
int get_op_local(double *X,double *dX,int D,double N,double dX2Add,double t,list<interaction> connections,double rho,double ending)
{
	list<interaction>::iterator it,at;
	stringstream ss;
	complex<double> opam(0,0),im_unit(0,6.0);
	double dX2,odX2,pdX[D],norm;
	int j,k,jD,kD,nu,l,lD,m,mD;
	double x=0.0001;
    long dX_PeriodicCount[D];
	ss<< "./run/OUT_loc_OP" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
	while(x<ending)
	{
		norm=0;
        opam=0;
      for(at = connections.begin(); at != connections.end(); ++at) 	{
        j = (*at).i1; jD=j*D;
        k = (*at).i2; kD=k*D;
		CALC_dX0(X,jD,kD,pdX);
        CALC_dX(pdX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(pdX);
		//in this case we are in the first pairs neighbourhood and want to look for any other neighbourhoods
        for(it = connections.begin(); it != connections.end(); ++it) {
          l = (*it).i1; lD=l*D;
          m = (*it).i2; mD=m*D;
          //	rewrite when the boxes are there
          //	|j-l|
	      CALC_dX0(X,jD,lD,dX);
          CALC_dX(dX, dX_PeriodicCount);
          odX2 = CALC_SQR_VEC(dX);
	      CALC_dX0(X,lD,mD,dX);
          CALC_dX(dX, dX_PeriodicCount);
						if (sqrt(odX2) <= x)
						{
						    opam+=exp(im_unit*(angle_to_uv(pdX,D)-angle_to_uv(dX,D)));
                           // if(j==m && k==l)
                           // cout << "added a value of " << exp(im_unit*abs((angle_to_uv(pdX,D)-angle_to_uv(dX,D)))) << " on an angle of " << abs(abs(angle_to_uv(pdX,D))-abs(angle_to_uv(dX,D))/(2.0*pi)) << " just as i looked at " << j << "->" << k << " versus " << l << "->" << m << endl;
						    norm++;
						}
					}
				}
      if (norm > 0 && x <= ending)
      {
	out << x << " " << (abs(opam)/norm) << endl;
    //cout << "just found opam " << abs(opam) << " on norm " << norm << endl;
      }
      else
          return 1;
      if(x==0.0001)
          x=0;
		x+=0.01*ending;
	}
    return 0;
}
int get_pkf(double *X,double *dX,int D,double N,double t,double dX2Add,double rho,double ending)
{
	stringstream ss;
    long dX_PeriodicCount[D];
	//double dXAdd=sqrt(dX2Add);
	double dXAdd=0.1*sqrt(dX2Add),pi=2.0*asin(1);
    int size=(int)round(ending/dXAdd)+1;
	int pkf[size],place;
    for (int kappa=0;kappa<size;kappa++)
        pkf[kappa]=0;
	//list<int>::iterator it;
	double dX2,x;
	int i,j,iD,jD,nu;
	ss<< "./run/OUT_pkf" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
		nu=0;
		//y=0;
	for (i=0;i<N;i++)
	{
		for (j=i+1;j<N;j++)
		{
			jD=j*D;
			iD=i*D;
		    CALC_dX0(X,iD,jD,dX);
            CALC_dX(dX, dX_PeriodicCount);
            dX2 = CALC_SQR_VEC(dX);
			dX2=sqrt(dX2);
            place=(int)round(dX2/dXAdd);
            if (place<size)
                pkf[place]++;
			//cout << "placing at " << place << " with dX2: " << dX2 << " and dX2Add: " << dX2Add << endl;
		}
	}
	//write out
	x=0.5*dXAdd;
	for(int kappa=0;kappa<size;kappa++)
	{
		out << x << " " << ((double)pkf[kappa])/(rho*N/2.0*pi*(pow(x+0.5*dXAdd,2)-pow(x-0.5*dXAdd,2))) << endl;
		x+=dXAdd;
	}
    return 0;
}
//get the mean length of a connection
void get_mean_length(double *X,double *dX,int D,list<interaction> connections,double* mean)
{
	list<interaction>::iterator it;
	int iD,jD,nu,n=0,i,j;
	double dX2;
    long dX_PeriodicCount[D];
	for(nu=0;nu<2;nu++)
		mean[nu]=0;
	//calculate mean value
	for (it=connections.begin();it!=connections.end();++it)
	{
		n++;
		i=(*it).i1;
		j=(*it).i2;
        iD=i*D;
        jD=j*D;
		CALC_dX0(X,iD,jD,dX);
        CALC_dX(dX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(dX);
		mean[0]+=sqrt(dX2);
	}
	if (n >= 1)
	{
	mean[0]=mean[0]/((float)n);
	//calculate mean deviation
	for (it=connections.begin();it!=connections.end();++it)
	{
		i=(*it).i1;
		j=(*it).i2;
        iD=i*D;
        jD=j*D;
		CALC_dX0(X,iD,jD,dX);
        CALC_dX(dX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(dX);
		mean[1]+=pow(sqrt(dX2)-mean[0],2);
	}
    if(n>=2)
	mean[1]=sqrt(mean[1]/((int)(n*(n-1))));
    else
	mean[1]=0;
	}
	else
	{
		mean[0]=0;
		mean[1]=0;
	//return mean;
	}
}
void get_mean_length_distr(double* X,list<interaction> connections,double t,double* dX,double max,int D)
{
    stringstream ss;
	ss<< "./run/OUT_length_distr" << t << ".txt";
	ofstream out(((string)ss.str()).c_str());
    double dX2;
    long dX_PeriodicCount[D];
    int i,iD,j,jD,nu;
    list<interaction>::iterator it;
    //for (it=connections.begin();it!=connections.end();++it)
    //{
    //        i=(*it).a;
    //        j=(*it).b;
	//		jD=j*D;
	//		iD=i*D;
	//	    CALC_dX0(X,iD,jD,dX);
    //        CALC_dX(dX, dX_PeriodicCount);
    //        dX2 = CALC_SQR_VEC(dX);
	//		dX2=sqrt(dX2);
    //        if(dX2 > max)
    //            max=dX2;
    //}
    int size=20,place;
    double arr[size];
    for (int kappa=0;kappa<size;kappa++)
        arr[kappa]=0;
    for (it=connections.begin();it!=connections.end();++it)
    {
            i=(*it).i1;
            j=(*it).i2;
			jD=j*D;
			iD=i*D;
		    CALC_dX0(X,iD,jD,dX);
            CALC_dX(dX, dX_PeriodicCount);
            dX2 = CALC_SQR_VEC(dX);
			dX2=sqrt(dX2);
            if(dX2>max)
            {
                cout << "detected breakout with length " << dX2 << endl;
            }
            else
            {
            place=(int)(dX2/max*(size-1));
            arr[place]++;
            }
    }
    for (int kappa=0;kappa<size;kappa++)
        out << max/(2.0*(double)(size-1))+max/((double)(size-1))*(double)kappa << " " << arr[kappa]/connections.size() << endl;
}
//get number of connections
double get_num_conn(list<interaction> connections,double N)
{
//	int counter=0;
//	list<interaction>::iterator it;
//		for (it=connections.begin();it!=connections.end();++it)
//		{
//				counter++;
//		}
//		return (((double)counter)/((double)N));
	return (connections.size()/N);
}
//get mean velocity with deviation
double get_mean_vel(double* V,double N,int D)//this is now used to measure middle temperature
{
    double temp=0;
    for(int kappa=0;kappa<D*N;kappa++)
        temp+=V[kappa]*V[kappa];
    return 0.5*temp/N;
}
int force_on_shell(double* V,double kT,double D,double N)
{
    double ekin=get_mean_vel(V,N,(int)D);
    double kappa=1.0/(sqrt(ekin/kT));
    for (int erta=0;erta<N*D;erta++)
        V[erta]=V[erta]*kappa;
    return 0;
}
int plot_delaunay(list<interaction> connections,double t,int D,double* X,double* L,string ID,double radius,list<int> marked,int force_type)
{
    //the array marked is a list giving numbers to mark with color on the plane
    //open up an image file
    int color[3*N];
    for (int erta=0;erta<N;erta++)
    {
        color[3*erta]=255;
        color[3*erta+1]=0;
        color[3*erta+2]=0;
    }
    list<int>::iterator intit;
    for (intit=marked.begin();intit!=marked.end();++intit)
    {
        if (intit==marked.begin())
        {
            color[3*(*intit)]=0;
            color[3*(*intit)+1]=255;
        }
        else
        {
            color[3*(*intit)]=0;
            color[3*(*intit)+2]=255;
        }
    }
    if (radius==0)
        radius=0.1;
    const int w=1024;
    const int h=1024;
    vector<unsigned char> image;
    image.resize(w*h*4);
    stringstream ss;
    int i,iD,j,jD,nu,px,py;
    double dX[D];
    long dX_PeriodicCount[D];
    double dX2,x,y,ymin,ymax,xmin,xmax;//stepwidth
	ss<< "./run/OUT_" << ID << t << ".png";
//	ofstream out(((string)ss.str()).c_str());
    list<interaction>::iterator it;
            for (x=0;x<w;x++)
            {
                for(y=0;y<h;y++)
                {
            image[4*w*y+4*x+0]=(unsigned char)255;
            image[4*w*y+4*x+1]=(unsigned char)255;
            image[4*w*y+4*x+2]=(unsigned char)255;
            image[4*w*y+4*x+3]=(unsigned char)255;
                }
            }
     //       int current=0;
    if (force_type!=5)//5 is no force. Other implementation causes problems blablabla, it works
    {
    for (it=connections.begin();it!=connections.end();++it)
    {
    //    current++;
    //    if(current>12)
    //        break;
        i=it->i1; iD=i*D;
        j=it->i2; jD=j*D;
        if (i>N || j>N)
            cout << "detected index out of bounds i:" << i << " j:" << j << endl;
		CALC_dX0(X,jD,iD,dX);
        CALC_dX(dX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(dX);
        dX2=sqrt(dX2);
        xmin=X[iD];
        xmax=X[iD]+dX[0];
        ymin=X[iD+1];
        ymax=X[iD+1]+dX[1];
        if(xmax != xmin || ymax != ymin)
        {
        for(int kappa=0;kappa<D;kappa++)
            dX[kappa]=dX[kappa]/(100.0*dX2);
        x=xmin;
        y=ymin;
        while (abs(x-xmin)<abs(xmax-xmin) || abs(y-ymin)<abs(ymax-ymin))
        {
            x+=dX[0];
            y+=dX[1];
            px=w*(mod_by_per(x,L[0])/L[0]);
            py=h*(mod_by_per(y,L[1])/L[1]);
            if(px>=w)
                px=w-1;
            if(py>=h)
                py=h-1;
            //cout << "plotting at " << px << " " << py << endl;
            image[4*w*py+4*px+0]=(unsigned char)0;
            image[4*w*py+4*px+1]=(unsigned char)0;
            image[4*w*py+4*px+2]=(unsigned char)255;
            image[4*w*py+4*px+3]=(unsigned char)255;
            //cout << x << " " << y << " " << dX[0] << " " << dX[1] << endl;
        }
        }
        else
        {
            cout << "Error at " << xmax << "=" << xmin << " and " << ymax << "=" << ymin << endl;
        }
    }
    }
        //add particle positions
        for(int nu=0;nu<N;nu++)
        {
            for(int ddx=-ceil(w*radius/L[0]);ddx<=ceil(w*radius/L[0]);ddx++)
            {
                for(int ddy=-ceil(h*radius/L[1]);ddy<=ceil(h*radius/L[1]);ddy++)
                {
                    if ((pow(ddx,2)+pow(ddy,2))<(pow(w*radius/L[0],2)))
                    {
            px=w*(mod_by_per(X[nu*D],L[0])/L[0])+ddx;
            py=h*(mod_by_per(X[nu*D+1],L[1])/L[1])+ddy;
            if(px>=w)
                px=w-1;
            if(py>=h)
                py=h-1;
            if(py<0)
                py=0;
            if(px<0)
                px=0;
            image[4*w*py+4*px+0]=(unsigned char)color[3*nu];
            image[4*w*py+4*px+1]=(unsigned char)color[3*nu+1];
            image[4*w*py+4*px+2]=(unsigned char)color[3*nu+2];
            image[4*w*py+4*px+3]=(unsigned char)255;
                    }
                }
            }
        }
        //include some nicers now, for example the bar indicating length scale
//        double bar_x=round(L[0]/10.0),bar_y=round(L[1]/10.0),bar_length=ceil(L[0]/10.0),ddx,ddy;
//        ddy=bar_y;
//        for(ddx=bar_x;ddx<bar_x+bar_length;ddx+=L[0]/w)
//        {
//            if(abs(ddx-bar_x-floor(ddx-bar_x)) < L[0]/w)
//            {
//                for(py=bar_y-1;py<bar_y+2;py++)
//                {
//            px=w*(mod_by_per(X[nu*D],L[0])/L[0])+ddx;
//            py=h*(mod_by_per(X[nu*D+1],L[1])/L[1])+ddy;
//            if(px>=w)
//                px=w-1;
//            if(py>=h)
//                py=h-1;
//            if(py<0)
//                py=0;
//            if(px<0)
//                px=0;
//            image[4*w*py+4*px+0]=(unsigned char)0;
//            image[4*w*py+4*px+1]=(unsigned char)0;
//            image[4*w*py+4*px+2]=(unsigned char)255;
//            image[4*w*py+4*px+3]=(unsigned char)255;
//                }
//            }
//            else
//            {
//            px=w*(mod_by_per(X[nu*D],L[0])/L[0])+ddx;
//            py=h*(mod_by_per(X[nu*D+1],L[1])/L[1])+ddy;
//            if(px>=w)
//                px=w-1;
//            if(py>=h)
//                py=h-1;
//            if(py<0)
//                py=0;
//            if(px<0)
//                px=0;
//            image[4*w*py+4*px+0]=(unsigned char)0;
//            image[4*w*py+4*px+1]=(unsigned char)0;
//            image[4*w*py+4*px+2]=(unsigned char)255;
//            image[4*w*py+4*px+3]=(unsigned char)255;
//            }
//        }

    LodePNG::Encoder encoder;
    encoder.addText("comment","Delaunay image, created with lodepng");
    encoder.getSettings().zlibsettings.windowSize=2048;
    vector<unsigned char> buffer;
    encoder.encode(buffer,image.empty() ? 0 : &image[0],w,h);
    if(image.empty())
        cout << "Empty image at " << t << endl;
    LodePNG::saveFile(buffer,((string)ss.str()).c_str());
    return 0;
}
int main(int argc, char** argv){
  //if(argc!=2 && argc!=4){
  list<interaction> neigbours;//usage of this list is good for completeness and bad for speed,sorry for the placement
  if( argc % 2 != 0 ) {
    printf("usage:\n");
    printf("./%s <Settings-File>\n\n",argv[0]);
    printf("<Settings-File> is a file that consequtively contains:\n");
    printf("1. File-Name for Positions of Atoms\n");
    printf("   must be lines containing continuouus particle numbers and particle coordinates and velocities\n");
    printf("2. \n");
    printf(".\n");
    printf(".\n");
    //WALLaFILE="lll.l";
    //cout addNDot(WALLaFILE,2) << "   " << WALLaFILE << endl;
    exit(1);
  }
  
  
  // ******************* prepare variables for Settings-File ******************* 
  char PosFILE[FILENAME_LENGTH];     strcpy(PosFILE, "POS.txt");
  char VelFILE[FILENAME_LENGTH];     strcpy(VelFILE, "VEL.txt");
  char ExtWallFILE[FILENAME_LENGTH]; strcpy(ExtWallFILE,  "EXT.txt");
  char ConnFILE[FILENAME_LENGTH];    strcpy(ConnFILE,     "CONN.txt");
  
  char PosOutFILE[FILENAME_LENGTH];     strcpy(PosOutFILE,   "run/OUTx.txt");
  char VelOutFILE[FILENAME_LENGTH];     strcpy(VelOutFILE,   "run/OUTv.txt");
  char VelFluctOutFILE[FILENAME_LENGTH];strcpy(VelFluctOutFILE,"none");
  char ConnOutFILE[FILENAME_LENGTH];    strcpy(ConnOutFILE,  "run/OUTc.txt");
  char SigmaOutFILE[FILENAME_LENGTH];   strcpy(SigmaOutFILE, "none");
  
  char BoundaryMoveFILE[FILENAME_LENGTH];    strcpy(BoundaryMoveFILE   , "BOUNDARYm.txt");
  char BoundaryOsciFILE[FILENAME_LENGTH];    strcpy(BoundaryOsciFILE   , "BOUNDARYo.txt");
  char WallMoveFILE[FILENAME_LENGTH];        strcpy(WallMoveFILE       , "WALLm.txt");
  char BoundaryPosOutFILE[FILENAME_LENGTH];   strcpy(BoundaryPosOutFILE  , "run/BOUNDARYpos.txt");
  char BoundaryForceOutFILE[FILENAME_LENGTH]; strcpy(BoundaryForceOutFILE, "run/BOUNDARYf.txt");
  char WallForceOutFILE[FILENAME_LENGTH];    strcpy(WallForceOutFILE   , "run/WALLf.txt");
  char WallPosOutFILE[FILENAME_LENGTH];      strcpy(WallPosOutFILE     , "run/WALLpos.txt");
  
  char IntegrationMethod[100];               strcpy(IntegrationMethod  , "default");
  char cuttingMethod[100];                   strcpy(cuttingMethod      , "moreThan");
  long cuttingMoreThanCoordNum = 3;
  floats AvCoordNum = -1.;
  
  int    weakSpringsQ = 0;
  floats weakSpringsa = 0.;
  floats weakSpringsc = 0.;
  
  long continueAt = -1;
  char lineIN[LINELENGTH];
  char *split;
  

  floats TMAX;
  floats dt_pos, dt_force;
  floats dX2Add;
  floats dX2Remove;
  floats AddRate; 
  floats AddNumber;
  floats RemoveRate,RemoveProb;
  floats Fmax,Fmax2;
  floats defaultDamping = 0.1;
  long OUTPRECISION = 6;
  short preventOverwriteQ = 0;
  //declarations by lutz Kuenneke
  ofstream msd("./MSD.txt"),op_simple("./OP_simple.txt"),op("./OP.txt"),length("./MEAN_LENGTH.txt"),conns("./MEAN_CONN.txt"),mean_vel("./MEAN_VEL.txt"),iface,total_energy("TOTAL_ENERGY.txt"),sigma_out("SIGMA_OVER_MU.txt");
  double dt_msd=0,dt_op_simple=0,dt_op=0,dt_op_local=0,dt_mean_length=0,dt_mean_num=0,dt_pkf=0,dt_mean_vel=0,dt_msd_distr=0,dt_mean_length_distr=0,dt_plot_delaunay=0,dt_voronoy=0,force_scaling=1,dt_energy_distr=0,dt_life_time=0,dt_total_energy=0,dt_plot_work_hist=0,dt_sigma_over_mu=0;
  int do_force_on_shell=0;
  int force_type=4;//default is constant force
  double mean[2],dt,shielding=1;
  double cell_radius=0;
  //end of declarations
  long i, j, k, nu;
  long conn_i1, conn_i2;
  long conn_i1D,conn_i2D, kD, kDpnu;
  
  for(nu=0; nu<MAXDIM; nu++)
    NSubBoxes[nu] = 1;
  
  
  
  
  // ******************* read Settings-File ******************* 
  ifstream SET(argv[1]);
  if (SET.is_open()) {
    cout << ":::::::::::::::::::::::::: reading input file " << argv[1] << " ::::::::::::::::::::::::::" << endl;
    cout << "   contents: " << endl;
    k = 0;
    while(!SET.eof()) {
      
      SET.getline(lineIN, LINELENGTH);
      //cout << lineIN << endl;
      k++;
      cout << ":" << k << ":\t";
      if ((split = strtok (lineIN, DELIMS)) != NULL ) {
        
        //cout << split << endl;
        
        if (strcmp(split, "InFilePos") == 0) {
          strcpy(PosFILE, strtok (NULL, DELIMS));
          cout << "File with initial particle positions: \t" << PosFILE;
        }
        
        else if (strcmp(split, "InFileVel") == 0) {
          strcpy(VelFILE, strtok (NULL, DELIMS)); 
          cout << "File with initial particle velocities: \t" << VelFILE;
        }
	//input from Lutz Kuenneke
	else if (strcmp(split, "force_type") == 0) {
          force_type = atof(strtok(NULL, DELIMS));
          cout << "Force type is = " << type_to_name(force_type);
        }
	else if (strcmp(split, "dt_mean_vel") == 0) {
          dt_mean_vel = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of mean velocity data: dt_mean_vel = " << dt_mean_vel;
        }
	else if (strcmp(split, "dt_life_time") == 0) {
          dt_life_time = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of lifetime data: dt_life_time = " << dt_life_time;
        }
	else if (strcmp(split, "dt_energy_distr") == 0) {
          dt_energy_distr = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of work - lifetime: dt_energy_distr = " << dt_energy_distr;
        }
	else if (strcmp(split, "force_scaling") == 0) {
          force_scaling = atof(strtok(NULL, DELIMS));
          cout << "Absolute strength of force is now set to = " << force_scaling;
        }
	else if (strcmp(split, "shielding") == 0) {
          shielding = atof(strtok(NULL, DELIMS));
          cout << "shielding constant is = " << shielding;
        }
	else if (strcmp(split, "force_on_shell") == 0) {
          do_force_on_shell = atof(strtok(NULL, DELIMS));
          if(do_force_on_shell)
                cout << "forcing particels on shell\n";
          else
              cout << "not forcing particles on shell\n";
        }
	else if (strcmp(split, "dt_sigma_over_mu") == 0) {
          dt_sigma_over_mu= atof(strtok(NULL, DELIMS));
          cout << "Time step for output of sigma over mu data: dt_sigma_over_mu= " << dt_sigma_over_mu;
        }
	else if (strcmp(split, "dt_msd") == 0) {
          dt_msd = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of msd data: dt_msd = " << dt_msd;
        }
	else if (strcmp(split, "cell_radius") == 0) {
          cell_radius = atof(strtok(NULL, DELIMS));
          cout << "Hard Core Cell Radius = " << cell_radius;
        }
	else if (strcmp(split, "dt_total_energy") == 0) {
          dt_total_energy = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of Total ENergy: dt_total_energy= " << dt_total_energy;
        }
	else if (strcmp(split, "dt_voronoy") == 0) {
          dt_voronoy= atof(strtok(NULL, DELIMS));
          cout << "Time step for output of msd data: dt_voronoy= " << dt_voronoy;
        }
	else if (strcmp(split, "dt_plot_delaunay") == 0) {
          dt_plot_delaunay = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of Delaunay Pictures: dt_delaunay = " << dt_plot_delaunay;
        }
	else if (strcmp(split, "dt_plot_work_hist") == 0) {
          dt_plot_work_hist = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of work hist: dt_plot_work_list= " << dt_plot_work_hist;
        }
	else if (strcmp(split, "dt_mean_length_distr") == 0) {
          dt_mean_length_distr = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of Mean length distribution data: dt_mean_length_distr = " << dt_mean_length_distr;
        }
	else if (strcmp(split, "dt_msd_distr") == 0) {
          dt_msd_distr = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of msd distribution data: dt_msd_distr = " << dt_msd_distr;
        }

        else if (strcmp(split, "dt_pkf") == 0) {
          dt_pkf = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of pkf data: dt_pkf = " << dt_pkf;
        }

        else if (strcmp(split, "dt_op_simple") == 0) {
          dt_op_simple = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of simple op data: dt_op_simple = " << dt_op_simple;
        }

        else if (strcmp(split, "dt_mean_length") == 0) {
          dt_mean_length = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of mean connection length: dt_mean_length = " << dt_mean_length;
        }

        else if (strcmp(split, "dt_mean_num") == 0) {
          dt_mean_num = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of the mean number of bounds: dt_mean_num = " << dt_mean_num;
        }

        else if (strcmp(split, "dt_op") == 0) {
          dt_op = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of op data: dt_op = " << dt_op;
        }

        else if (strcmp(split, "dt_op_local") == 0) {
          dt_op_local = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of local op data: dt_op_local = " << dt_op_local;
        }
	//input end
        
        else if (strcmp(split, "InFileConn") == 0) {
          strcpy(ConnFILE, strtok (NULL, DELIMS)); 
          cout << "File with initial particle bonds: \t" << ConnFILE;
        }
        
        else if (strcmp(split, "InFileExt") == 0) {
          strcpy(ExtWallFILE, strtok (NULL, DELIMS)); 
          cout << "File with extra particle properties: \t" << ExtWallFILE;
        }
        
        else if (strcmp(split, "OutFilePos") == 0) {
          strcpy(PosOutFILE, strtok (NULL, DELIMS)); 
          cout << "Output file for particle positions: \t" << PosOutFILE;
        }
        
        else if (strcmp(split, "OutFileVel") == 0) {
          strcpy(VelOutFILE, strtok (NULL, DELIMS)); 
          cout << "Output file for particle velocities: \t" << VelOutFILE;
        }

        else if (strcmp(split, "OutFileVFluct") == 0) {
          strcpy(VelFluctOutFILE, strtok (NULL, DELIMS)); 
          cout << "Output file for particle velocities relative to mean field: \t" << VelFluctOutFILE;
        }

        else if (strcmp(split, "OutFileConn") == 0) {
          strcpy(ConnOutFILE, strtok (NULL, DELIMS)); 
          cout << "Output file for particle bonds: \t" << ConnOutFILE;
        }
        
        else if (strcmp(split, "OutFileSigma") == 0) {
          strcpy(SigmaOutFILE, strtok (NULL, DELIMS));
	  if (strcmp(SigmaOutFILE, "none") == 0) cout << "OutFileSigma = none: no output of local stress tensor";
          else 				  	 cout << "Output file for local stress tensor: \t" << SigmaOutFILE;
        }
        
        else if (strcmp(split, "InFileBoundaryMove") == 0) {
          strcpy(BoundaryMoveFILE, strtok (NULL, DELIMS)); 
          cout << "Input  file for boundary moves: \t" << BoundaryMoveFILE;
        }
        
        else if (strcmp(split, "InFileBoundaryOsci") == 0) {
          strcpy(BoundaryOsciFILE, strtok (NULL, DELIMS)); 
          cout << "Input  file for boundary oscillations: \t" << BoundaryOsciFILE;
        }
        
        else if (strcmp(split, "OutFileBoundaryPos") == 0) {
          strcpy(BoundaryPosOutFILE, strtok(NULL, DELIMS));
          cout << "Output file for pos of boundaries: \t" << BoundaryPosOutFILE;
        }
        
        else if (strcmp(split, "OutFileBoundaryForce") == 0) {
          strcpy(BoundaryForceOutFILE, strtok(NULL, DELIMS));
          cout << "Output file for forces on boundaries: \t" << BoundaryForceOutFILE;
        }
        
        else if (strcmp(split, "InFileWallMove") == 0) {
          strcpy(WallMoveFILE, strtok(NULL, DELIMS));
          cout << "Input  file for wall moves: \t\t" << WallMoveFILE;
        }
        
        else if (strcmp(split, "OutFileWallForce") == 0) {
          strcpy(WallForceOutFILE, strtok(NULL, DELIMS));
          cout << "Output file for forces on walls: \t" << WallForceOutFILE;
        }
        
        else if (strcmp(split, "OutFileWallPos") == 0) {
          strcpy(WallPosOutFILE, strtok(NULL, DELIMS));
          cout << "Output file for positions of walls: \t" << WallPosOutFILE;
        }

        else if (strcmp(split, "dt") == 0) {
          dt = (double)atof(strtok(NULL, DELIMS));
          dt_h = dt*0.5;
          cout << "dt =\t" << dt;
        }
        
        else if (strcmp(split, "tmax") == 0) {
          TMAX = atof(strtok(NULL, DELIMS));
          cout << "tmax =\t" << TMAX;
        }
        
        else if (strcmp(split, "gamma") == 0) {
          GAMMA  = atof(strtok(NULL, DELIMS));
          cout << "gamma =\t" << GAMMA;
        }
        
        else if (strcmp(split, "sigma") == 0) {
          SIGMA = atof(strtok(NULL, DELIMS));
          cout << "sigma =\t" << SIGMA;
        }
        
        else if (strcmp(split, "dXAdd") == 0) {
          dX2Add = atof(strtok(NULL, DELIMS));
          cout << "maximum particle distance for adding bonds (dXAdd) =\t" << dX2Add;
          dX2Add *= dX2Add;
          cout << "\t\tdXAdd^2 = " << dX2Add;
        }
        
        else if (strcmp(split, "dXRemove") == 0) {
          dX2Remove = atof(strtok(NULL, DELIMS));
          cout << "particle distance beyond which bonds break (dXRemove) =\t" << dX2Remove;
          dX2Remove *= dX2Remove;
          cout << "\t\tdXRemove^2 = " << dX2Remove;
        }
        
        else if (strcmp(split, "FRemove") == 0) {
          Fmax = atof(strtok(NULL, DELIMS));
          if (Fmax >= 0.) {
            cout << "Force beyond which bonds break (FRemove) =\t" << Fmax;
            Fmax2 = Fmax*Fmax;
          }
          else {
            cout << "Bonds do not break due to maximum force.";
            Fmax2 = -1.;
          }
        }
        
        else if (strcmp(split, "AddRate") == 0) {
          AddRate   = (double)atof(strtok(NULL, DELIMS));
          AddNumber = (double)(AddRate*dt);
          cout << "rate, at which particles closer than dXAdd are bonded =\t" << AddRate;
          cout << "   --> AddNumber per time step = " << AddNumber;
        }
        
        else if (strcmp(split, "RemoveRate") == 0) {
          RemoveRate = (double)atof(strtok(NULL, DELIMS));
          RemoveProb = RemoveRate*dt;//formerly multiplied with  dt
          cout << "rate, at which particles bonds break spontaneously =\t" << RemoveRate;
          cout << "   --> RemoveProb per time step = " << RemoveProb;
        }
        
        else if (strcmp(split, "DampingDef") == 0) {
          defaultDamping = atof(strtok(NULL, DELIMS));
          cout << "default damping coefficient (if entry 7 in \"" << ConnFILE << "\" is empty) c = " << defaultDamping;
        }
	
        else if (strcmp(split, "AvCoordNum") == 0) {
          AvCoordNum = atof(strtok(NULL, DELIMS));
	  if(AvCoordNum >= 0.) 
	    cout << "bonds are randomly cut to reach average coordination number: z = " << AvCoordNum;
	  else
	    cout << "AvCoordNum < 0: no bonds being cut";
        }
        
        else if (strcmp(split, "cutMethod") == 0) {
          strcpy(cuttingMethod, strtok(NULL, DELIMS));
          cout << "Cutting Method: \t";
          if(strcmp(cuttingMethod, "moreThan") == 0) {
	    cuttingMoreThanCoordNum = atol( strtok(NULL, DELIMS) );
            cout << "cut bonds coordination with number more than " << cuttingMoreThanCoordNum;
	  }

          else if(strcmp(cuttingMethod, "highestFirst") == 0) {
            cout << "cut bonds with highest coordination first";
	  }

          else {
            strcpy(cuttingMethod, "moreThan");
	    cuttingMoreThanCoordNum = 3;
            cout << split << " not recognized. Using default method: \"moreThan 3\"";
          }
        }
	
        else if (strcmp(split, "weakSprings") == 0) {
	  if( (split = strtok(NULL, DELIMS)) != NULL) weakSpringsa = atof(split);
	  if( (split = strtok(NULL, DELIMS)) != NULL) weakSpringsc = atof(split);
	  if( weakSpringsa != 0. || weakSpringsc != 0. ) {
	    weakSpringsQ = 1;
	    cout << "Instead of cutting bonds: use \"weak springs\" with spring constant a = " << weakSpringsa << ", and damping coefficient c = " << weakSpringsc << ".";
	  }
	  else {
	    cout << "Cutting bonds. Not using \"weak springs\".";
	  }
        }
        
	
        else if (strcmp(split, "dt_pos") == 0) {
          dt_pos = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of particle data: dt_pos = " << dt_pos;
        }
        
        else if (strcmp(split, "dt_force") == 0) {
          dt_force = atof(strtok(NULL, DELIMS));
          cout << "Time step for output of force data:  dt_force = " << dt_force;
        }
        
        else if (strcmp(split, "BoxSize") == 0) {
          nu = 0;
          cout << "Box Size:   ";
          while( (split = strtok(NULL, DELIMS)) != NULL) {
            L[nu]  = atof( split );
            Lh[nu] = L[nu] * 0.5;
            cout << "L[" << nu << "] = " << L[nu] << "; ";
            nu++;
          }
        }
	
	else if (strcmp(split, "SubBoxes") == 0) {
          nu = 0;
          cout << "Number of subboxes:   ";
          while( (split = strtok(NULL, DELIMS)) != NULL) {
            NSubBoxes[nu]  = atol( split );
            cout << "N[" << nu << "] = " << NSubBoxes[nu] << "; ";
            nu++;
          }
        }
	
	else if (strcmp(split, "pseudoD") == 0) {
          pseudoD = atol(strtok(NULL, DELIMS));
          cout << "pseudo dimension: " << pseudoD;
        }
	
        else if (strcmp(split, "method") == 0) {
          strcpy(IntegrationMethod, strtok(NULL, DELIMS));
          cout << "Integration Method: \t";
          if(strcmp(IntegrationMethod, "LeapfrogNoDamp") == 0) 
            cout << "Leapfrog WITHOUT efficient damping";

          else if(strcmp(IntegrationMethod, "Leapfrog") == 0) 
            cout << "Leapfrog with efficient damping";

          else if(strcmp(IntegrationMethod, "VelocityVerlet") == 0) 
            cout << "Velocity Verlet with efficient damping";

          else if(strcmp(IntegrationMethod, "Euler") == 0) 
            cout << "Euler integration";

          else if(strcmp(IntegrationMethod, "Overdamped") == 0) {
            cout << "Overdamped dynamics, using Heun method. Particle masses are ignored.";
	    OverdampedTrue = 1;
	  }

          else {
            strcpy(IntegrationMethod, "Heun");
            cout << "default algorithm: Heun method";
          }
        }
	
	else if (strcmp(split, "outPrecisn") == 0) {
          OUTPRECISION = atol( strtok(NULL, DELIMS) );
          cout << "Precision (number of digits) for pos, vel and conn output files: " << OUTPRECISION;
        }
	
	else if (strcmp(split, "continueAt") == 0) {
          continueAt = atol( strtok(NULL, DELIMS) );
	  if(continueAt >= 0) {
	    cout << "Trying to continue the current run at output file # " << continueAt;
	    preventOverwriteQ = 1;
	  }
	  else
	    cout << "starting a new run at time 0";
        }


        else {
          cout << "keyword \"" << split << "\" ignored.";
        }
        
        
      }
      cout << endl;
    
    }
    SET.close();
    cout << "::::::::::::::::::::::::::        reading done       ::::::::::::::::::::::::::" << endl << endl;
  }
  
  if (argc >= 4){
    
    for (long argcounter = 2; argcounter < argc; argcounter += 2) {
      cout << "Overwriting the following settings from Settings-file \"" << argv[1] << "\":" << endl;
      if(strcmp(argv[argcounter], "continueAt") == 0) {
	continueAt = atol( argv[argcounter+1] );
	cout << "\"continueAt\": ";
	if(continueAt >= 0) {
	  cout << "Trying to continue the current run at output file # " << continueAt << endl;
	  preventOverwriteQ = 1;
	}
	else
	  cout << "starting a new run at time 0." << endl;
	//cout << "\"continueAt\"-enrty in " << argv[1] << "is being overwritten." << endl;
      }
      
      else if(strcmp(argv[argcounter], "tmax") == 0) {
	  TMAX = atof( argv[argcounter+1] );
          cout << "\"tmax\": tmax =\t" << TMAX << endl;
      }
    }
  }
  

  if(continueAt >= 0) {
    strcpy(PosFILE,  addNDot(PosOutFILE,continueAt).c_str());  cout << "Changing InFilePos to:  " << PosFILE << endl;
    strcpy(VelFILE,  addNDot(VelOutFILE,continueAt).c_str());  cout << "Changing InFileVel to:  " << VelFILE << endl;
    strcpy(ConnFILE, addNDot(ConnOutFILE,continueAt).c_str()); cout << "Changing InFileConn to: " << ConnFILE << endl;
    cout << endl;
  }

  CheckParticleNumber(PosFILE, &N, &D);
  ND=N*D;
  cout << "Summary:" << endl;
  cout << "========" << endl;
  cout << "Number of particles: " << N << endl;
  cout << "Spatial dimension:   " << D << endl;
  
  
  cout << "System Size:         L    = (" << L[0];
  for(nu=1; nu<D; nu++) cout << ", " << L[nu];
  cout << ")" << endl;
  
  if(0 <= pseudoD && pseudoD <= D)
    cout << "pseudo-dimension = " << pseudoD << ". POS- and VEL-components [nu > " << pseudoD << "] are fixed." << endl;
  else
    pseudoD = D;
  
  /*cout << "Sub Boxes:           N_sb = (" << NSubBoxes[0];
  for(nu=1; nu<D; nu++) cout << ", " << NSubBoxes[nu];
  cout << ")  (currently ignored)" << endl;*/
  
  if(AvCoordNum >= 0.) {
    cout << "bonds are randomly cut to reach average coordination number: z = " << AvCoordNum << endl;
    cout << "Cutting Method: \t";
    if(strcmp(cuttingMethod, "moreThan") == 0) 
      cout << "cut bonds with participants' coordination number more than " << cuttingMoreThanCoordNum << "." << endl;
    else if(strcmp(cuttingMethod, "highestFirst") == 0) 
      cout << "cut bonds with highest coordination first"  << endl;
  }
  else
    cout << "no bonds being cut" << endl;

  cout << "Trying to add " << AddNumber << " connections per time step" << endl;
  cout << "Removing connections with probability " << RemoveProb << " each time step" << endl;
  
  
  // ******************* initialize variables ******************* 
  AddNumber=AddNumber*(N*(N-1)/2);
  floats t;
  floats nextPosOut_t, nextForceOut_t;
  long count;
  long PosOutCount, ForceOutCount;
  //added by Lutz
  //bool its_okay;
  //bool not_moved[N];
  //double Xhard[D*N];
  //for (int erta=0;erta<D*N;erta++)
  //    Xhard[erta]=0;
  time_t seconds;
  double stats_time[4];
  for (int erta=0;erta<4;erta++)
      stats_time[erta]=0;
  list<int> marked;
  int is_moved[N];
  double work,work_max=1.0/(RemoveRate)*sqrt(dX2Add);
  if(force_scaling)
    work_max=force_scaling*work_max;
  if(AddRate)
    work_max=AddRate*work_max;
  //double work_max=10; //enable for no force debugging
  //double hard_scaling=force_scaling*AddRate/RemoveRate;
  double dw=work_max*0.001;
  int some_counter=(int)round(2*ceil(work_max/dw)+1);
  int work_histogramm[some_counter];
  for (int kappa=0;kappa<some_counter;kappa++)
      work_histogramm[kappa]=0;
  double diff_X[N*D];
  int conn_stats[N];
  for (int eta=0;eta<N;eta++)
      conn_stats[eta]=0;
  double xpos,ypos;
  list<int>::iterator int_it;
  int boxes[D];
  double cellr2=pow(2.0*cell_radius,2);
  //double excl_pow;
  double lc=sqrt(max(dX2Remove,cellr2));
  double box_size=lc;
  for (int eta=0;eta<D;eta++)
  {
      if(L[eta]/box_size == floor(L[eta]/box_size))
        boxes[eta]=(int)floor(L[eta]/(box_size));
      else
        boxes[eta]=(int)floor(L[eta]/(box_size))+1;
  }
  box boxmap[boxes[0]][boxes[1]];
  boxes_avail=true;
  cout << "using " << boxes[0] << " " << boxes[1] << " on boxsize " << box_size << endl;
  int nmnmeh=N*(N-1)/2;
  int shield_matrix[nmnmeh];
  for (int erta=0;erta<nmnmeh;erta++)
      shield_matrix[erta]=0;
  list<int> tempN;
  //double proj_1,proj_2,center;//fragment of my hard core times
  double norm=0;
  list<interaction>::iterator at;
  list<interaction> voronoy;
  int drop_isin,errorcode;
  int map[4*N],place;
  double div_op=0,div_op_local=0,div_neigbours=0;
  bool contain;
  int verts[3],tmp;
  double Xtemp[8*N];//this Vector is meant to be extended in every direction by the half box length
  int Tri_num;
  //double XEdges[24*N];
  int Tri_vert[24*N];
  int Tri_nabe[24*N];
  double msd_val[N],X0[N*D],rho,vol=1.0;
  for  (nu=0;nu<D;nu++)
      vol=vol*L[nu];
  rho=N/vol;
for (count=0;count<N;count++)
	msd_val[count]=0;
    double life_time_distr[101];
    for(int kappa=0;kappa<101;kappa++)
        life_time_distr[kappa]=0;
    double life_mean=1.0/RemoveRate;
    double life_max=2.0*life_mean;
    double life_step=0.01*life_max;
    double energy_distr[101];
    for(int kappa=0;kappa<101;kappa++)
        energy_distr[kappa]=0;
//end of Lutz
  SQRT3DT = sqrt(3.*dt);
  floats ONEminusGAMMAdt   = 1. - GAMMA*dt;
  floats ONEminusGAMMAdt_h = 1. - GAMMA*dt_h;
  //floats dtdtHalf = dt*dt*0.5;
  //floats SIGMAodt = SIGMA / dt;
  
  

  floats *X 	     = new floats[ND];	 // positions
  floats *V 	     = new floats[ND];	 // velocities
  floats *Vrand      = new floats[ND];	 // random velocity kicks
  floats *previous_X = new floats[ND];	 // positions from the previous time step
  floats *mi 	     = new floats[N];	 // inverse masses
  floats dX[D];
  floats dV[D];
  
  
  long dX_PeriodicCount[D];  // count how many periodic boundaries are crossed
  int  BoundaryConnQ[D];     // stores, whether the connection crosses the system boundary
  floats VplusAccDtH;
  floats dX2,dX0;
  floats dV2;
  floats dxbb;
  floats Fh,Fnu,c1Fh,c2Fh,dXdV;
  short conn_damp;
  
  floats *R1_X = new floats[ND];
  floats *R1_V = new floats[ND];
  floats *R1_F = new floats[ND];
  floats *R2_X = new floats[ND];
  floats *R2_V = new floats[ND];
  floats *R2_F = new floats[ND];
  
  

  floats X1_inv[D];
  floats X2_inv[D];
  floats dX_inv;
  
  int wi,wiD;
  int bi;//,biD;
  int Wallnk;
  int *Walln = new int[N];
  int WallN; //number of walls
  floats vWall[D]; //direction of a wall
  floats BoundaryForceElas[D][D];   //force on boundary ("stress tensor"), elastic contribution
  floats BoundaryForceDiss[D][D];   //force on boundary ("stress tensor"), dissipativ contribution
  floats SystemForceElas[D][D];     //summed forces ("stress tensor") in the system, elastic contribution
  floats SystemForceDiss[D][D];     //summed forces ("stress tensor") in the system, dissipativ contribution
  //floats SystemForceKine[D][D];     //summed forces ("stress tensor") in the system, kinetic contribution
  
  long Nsigma;
  if(strcmp(SigmaOutFILE, "none") == 0) Nsigma = 1;
  else 				 	Nsigma = N;

  floats ***sigmaLocalElas = new floats**[Nsigma];
  floats ***sigmaLocalDiss = new floats**[Nsigma];
  for(i=0; i<Nsigma; i++) {
    sigmaLocalElas[i] = new floats*[D];
    sigmaLocalDiss[i] = new floats*[D];
    for(j=0; j<D; j++) {
      sigmaLocalElas[i][j] = new floats[D];
      sigmaLocalDiss[i][j] = new floats[D];
    }
  }
  
  
  long NvNeighbors;
  if(strcmp(VelFluctOutFILE, "none") == 0) NvNeighbors = 1;
  else 				 	   NvNeighbors = N;
  long   *NNeighors     = new   long[NvNeighbors];
  floats *vNeighborsSum = new floats[NvNeighbors*D];
  floats  vField[D];


  long NConns     = 0;
  long NConns0    = 0;
  long NConnsGoal = -1;
  
  
  for (bi=0; bi<D; bi++) {
    for (nu=0; nu<D; nu++) {
      vBoundary[bi][nu] = 0.;
      if(bi==nu) L_LE[bi][nu] = L[nu];
      else       L_LE[bi][nu] = 0.;
      
      BoundaryOsci.amplitude[bi][nu] = 0.;
      BoundaryOsci.frequency[bi][nu] = 0.;
      BoundaryOsci.timeShift[bi][nu] = 0.;
    }
  }
  
  Inverse( L_LE, L_LEinv, D );
  
  
  interaction conn;
  list<interaction> connections;
  list<interaction>::iterator it;
  
  
  // ******************* read initial Position and Velocities ******************* 
  ifstream POS(PosFILE);
  for(k=0;k<ND;k++) 
	POS >> X[k];
  POS.close();
  
for (count=0;count<ND;count++)
	X0[count]=X[count];
  
  
  if( OverdampedTrue ) {
    cout << "velocity file \"" << VelFILE << "\" ignored due to use of overdamped dynamics." << endl;
    for(k=0;k<ND;k++)
      V[k] = 0.;
    
    if(GAMMA == 0.) {
      cout << "Warning: gamma = 0 invalid for overdamped dynamics. Setting gamma = 1." << endl;
      GAMMA = 1.;
    }
    GAMMAi      = 1./GAMMA;
    SIGMAoGAMMA = SIGMA/GAMMA;
    
  }
  else {
    ifstream VEL(VelFILE);
    if (VEL.is_open()) {
      for(k=0;k<ND;k++)
	    VEL >> V[k];
      VEL.close();
    }
    else {
      cout << "velocity file \"" << VelFILE << "\" not found. Setting initital velocities to 0." << endl;
      for(k=0;k<ND;k++)
	V[k] = 0.;
    }
  }
  //set initial kT for onshell forcing
    double kT=get_mean_vel(V,N,D);
  
  if(SIGMA==0.) {
    cout << "Sigma = 0. No temperature kicks." << endl;
    for(k=0;k<ND;k++)
      Vrand[k] = 0.;
  }
  
  
  
  
  

  // ******************* read Masses and Wall-Atoms ******************* 
  WallN = 0;
  ifstream EXTRAS(ExtWallFILE);
  if(EXTRAS.is_open()) {
    for(k=0; k<N; k++) {
      EXTRAS >> Wallnk; Walln[k] = Wallnk;    	// wall number of Atom k
      EXTRAS >> mi[k];   				// inverse mass
      if (Wallnk>0) {
	mi[k]=0.;
	if (Wallnk > WallN)
	  WallN = Wallnk;
      }
    }
    EXTRAS.close();
  }
  else {
    cout << "Ext-File \"" << ExtWallFILE << "\" not found. Not using walls and setting all masses to 1." << endl;
    for(k=0; k<N; k++) {
      Walln[k] = 0;    	// wall number of Atom k
      mi[k]    = 1.;    // inverse mass
    }
  }

  
  
  printf("Number of WALLS WallN = %d\n",WallN);
  
  long WallMoveSteps;
  long WallMoveStep;
  floats Wall_Force[(WallN+1)*D];
  long Wall_AtomMax[WallN+1];
  long Wall_AtomCount[WallN+1];
  long **Wall_AtomIndex;
  Wall_AtomIndex = new long*[WallN+1];
  
  long Wall_AtomCount0;     //number of "movable" atoms
  long *Wall_AtomIndex0;
      //pointer to Wall_AtomIndex[0], to save computation time
  
  long BoundaryMoveSteps;
  long BoundaryMoveStep;
  
  long BoundaryOsciSteps;
  long BoundaryOsciStep;
  
  
  for(wi=0; wi<=WallN; wi++) {
    Wall_AtomMax[wi] = 0;
    Wall_AtomCount[wi] = 0;
  }
  
  for(k=0; k<N; k++) 
    Wall_AtomMax[Walln[k]]++;
  
  Wall_AtomIndex = new long*[WallN+1];
  for(wi=0; wi<=WallN; wi++)
    Wall_AtomIndex[wi] = new long[Wall_AtomMax[wi]];
  
  for(k=0; k<N; k++) {
    wi=Walln[k];
    Wall_AtomIndex[wi][Wall_AtomCount[wi]++]=k;
  }
  
  Wall_AtomCount0 = Wall_AtomCount[0];  //number of "movable" atoms
  Wall_AtomIndex0 = Wall_AtomIndex[0];  
	  //pointer to Wall_AtomIndex[0], to save computation time
  
    
  // ******************* read connections between atoms ******************* 

  NConns = 0;
  ifstream CONN(ConnFILE);
  if(CONN.is_open()) {
    while(!CONN.eof()){
      CONN >> j;
      CONN >> k;
      if(j<k) { conn.i1 = j-1; conn.i2 = k-1; }
      else    { conn.i1 = k-1; conn.i2 = j-1; }
      
      if( !CONN.eof() ) {
	CONN >> conn.type;
	CONN >> conn.a;
	CONN >>    dX0;
	conn.dx0 = dX0;
	dx0_to_b(dX0, conn.b, conn.type);
	
	if(CONN.peek()==13 || CONN.peek()==10) {
	  if(defaultDamping == 0.) {
	    conn.damp  = 0;
	    conn.c     = 0.;
	  }
	  else {
	    conn.damp  = 1;
	    conn.c     = defaultDamping;
	  }
	    
	  conn.c1    = 1.;
	  conn.c2    = 0.;
	  dxbb = Fcrit_to_dxbb1(conn.type, conn.a, conn.dx0, Fmax); conn.dxbb1 = dxbb; conn.dx2bb1 = (dxbb>=0. ? SQR(dxbb) : -1.);
	  dxbb = Fcrit_to_dxbb2(conn.type, conn.a, conn.dx0, Fmax); conn.dxbb2 = dxbb; conn.dx2bb2 = (dxbb>=0. ? SQR(dxbb) : -1.);
	}
	else {
	  CONN >> conn.damp;
	  CONN >> conn.c;
	  CONN >> conn.c1;
	  CONN >> conn.c2;
	  if(CONN.peek()==13 || CONN.peek()==10 || CONN.eof()) {
	    dxbb = Fcrit_to_dxbb1(conn.type, conn.a, conn.dx0, Fmax); conn.dxbb1 = dxbb; conn.dx2bb1 = (dxbb>=0. ? SQR(dxbb) : -1.);
	    dxbb = Fcrit_to_dxbb2(conn.type, conn.a, conn.dx0, Fmax); conn.dxbb2 = dxbb; conn.dx2bb2 = (dxbb>=0. ? SQR(dxbb) : -1.);
	  }
	  else {
	    CONN >> dxbb; conn.dxbb1 = dxbb; conn.dx2bb1 = (dxbb>=0. ? SQR(dxbb) : -1.);
	    CONN >> dxbb; conn.dxbb2 = dxbb; conn.dx2bb2 = (dxbb>=0. ? SQR(dxbb) : -1.);
	  }
	}
    conn.life_time=0;
	connections.push_back(conn);
	NConns++;
      }
    }
    CONN.close();
    cout << "Number of connections read: " << NConns << endl;
  }
  else {
    cout << "Connection file \"" << ConnFILE << "\" not found. Initital number of connections: " << NConns << endl;
    //if no initial connections just set everything to equibrilibrium. Means we
    //add a/w*N(N-1)/2 connections
    //
	}
 //   cout << round(AddNumber) << " = " << AddNumber << "/" << RemoveProb << "*" << N*(N-1)/2 << endl;
 //   for(i=0; i<round(AddNumber/RemoveProb*N*(N-1)/2); i++) {
 //     j = (long)floor( drand48()*(double)N ); conn_i1D=j*D;
 //     k = (long)floor( drand48()*(double)N ); conn_i2D=k*D;
 //     if (j!=k) {
 //       
 //       CALC_dX0(X, conn_i1D, conn_i2D, dX);
 //       CALC_dX(dX, dX_PeriodicCount);
 //       dX2 = CALC_SQR_VEC(dX);
 //   
 //   if(dX2 < dX2Add) {
 //     if(j<k) {
 //       conn.i1 = j;
 //       conn.i2 = k;
 //     }
 //     else {
 //   	    conn.i1 = k;
 //       conn.i2 = j;
 //     }
 //     conn.type = force_type;
 //         conn.damp = 0;
 //     conn.a    = force_scaling;
 //         conn.dx0  = sqrt(dX2Add);
 //         dx0_to_b(conn.dx0, conn.b, conn.type);
 //         conn.c    = 0.;
 //         conn.c1   = 0.;
 //         conn.c2   = 0.;
 //         conn.dxbb1 = -1.;
 //         conn.dxbb2 = -1.;
 //         conn.life_time=0;
 //         conn.work=0;
 //         conn_stats[conn.i1]++;
 //         conn_stats[conn.i2]++;
 //     connections.push_back(conn);
 //     }
 //   }
 // }

  
  NConnsGoal = (long)floor((AvCoordNum*(floats)N)*.5 + .5);
  if(AvCoordNum >= 0. && NConns > NConnsGoal) {
    cout << "starting to cut connections until " << NConnsGoal << " are left... " << flush;
    
    
    long ParticleCoordNum[N];
    for(k=0; k<N; k++)
      ParticleCoordNum[k]=0;
    
    for(it = connections.begin(); it != connections.end(); ++it) {
      conn = *it;
      ParticleCoordNum[conn.i1]++;
      ParticleCoordNum[conn.i2]++;
    }
    
    connectionCutStruct connectionCut[NConns];
    k=0;
    for(it = connections.begin(); it != connections.end(); ++it) {
      conn = *it;
      connectionCut[k].i1   = conn.i1;
      connectionCut[k].i2   = conn.i2;
      connectionCut[k].cutQ = 0;
      k++;
    }
    
    if(k != NConns)
      cout << "Uhoh!" << endl;
    
    NConns0 = NConns;
    while(NConns > NConnsGoal) {
      k = (long)floor( drand48()*(double)NConns0 );
      if( ParticleCoordNum[connectionCut[k].i1]>cuttingMoreThanCoordNum && ParticleCoordNum[connectionCut[k].i2]>cuttingMoreThanCoordNum && connectionCut[k].cutQ == 0) {
	connectionCut[k].cutQ = 1;
	ParticleCoordNum[connectionCut[k].i1]--;
	ParticleCoordNum[connectionCut[k].i2]--;
	NConns--;
      }
    }
    
    
    if( weakSpringsQ == 1 ) {
      k=0; i=0;
      for(it = connections.begin(); it != connections.end(); ++it) {
	if( connectionCut[k].cutQ == 1 ) {
	  (*it).a = weakSpringsa;
	  (*it).c = weakSpringsc;
	}
	else 
	  i++;


	k++;
      }
      cout << "done. Number of leftover (non-weak) connections: " << i << endl;

    }
    else {
      k=0;
      for(it = connections.begin(); it != connections.end(); ) {
	if( connectionCut[k].cutQ == 1 )
    {
	  it = connections.erase(it); 
    }
	else
	  ++it;

	k++;
      }
      
      k=0;
      for(it = connections.begin(); it != connections.end(); ++it)
	k++;
      
      cout << "done. Number of connections: " << k << endl;
    }

  }
  
   // ******************* read wall movements ******************* 
  WallMoveSteps = WallN>0 ? CheckNumberOfLines(WallMoveFILE) : 0;
  WallMoveStruct WallMoves[WallMoveSteps];
  if (WallN>0) {
    
    ifstream WALLMOVE(WallMoveFILE);
    for(WallMoveStep=0; WallMoveStep<WallMoveSteps; WallMoveStep++){
      WALLMOVE >> WallMoves[WallMoveStep].time;
      WALLMOVE >> WallMoves[WallMoveStep].WallNumber;
      for(nu=0;nu<D;nu++)
        WALLMOVE >> WallMoves[WallMoveStep].direction[nu];
      
      
      if(WALLMOVE.peek()!=13 && WALLMOVE.peek()!=10 && !WALLMOVE.eof()) {
        cout << "Wallmoves have the wrong dimension";
        exit(1);
      }
      
      //cout << WallMoveStep << endl;
    }
    WALLMOVE.close();


    qsort (WallMoves, WallMoveSteps, sizeof(WallMoveStruct), compWM);

    cout << "Sorted Wall Moves:" << endl;
    for(WallMoveStep=0; WallMoveStep<WallMoveSteps; WallMoveStep++) {
      cout << WallMoves[WallMoveStep].time << " ";
      cout << WallMoves[WallMoveStep].WallNumber << " ";
      for(nu=0; nu<D; nu++)
	cout << WallMoves[WallMoveStep].direction[nu] << " ";
      cout << endl;
    }
  }
   
     // ******************* read boundary movements ******************* 
    
  LeesEdwardsTrue   = 0;  // do not consider Lees-Edwards boundary conditions, unless set to 1 later
  BoundaryMoveSteps = CheckNumberOfLines(BoundaryMoveFILE);
  
  if (BoundaryMoveSteps==0) 
    cout << "Warning: BoundaryMoveFile " << BoundaryMoveFILE << " empty. Hence, no boundary movments." << endl;

  if (BoundaryMoveSteps<0) {
    cout << "Warning: BoundaryMoveFile " << BoundaryMoveFILE << " does not exist. Hence, no boundary movments." << endl;
    BoundaryMoveSteps = 0;
  }
  
  BoundaryMoveStruct BoundaryMoves[BoundaryMoveSteps];
  if (BoundaryMoveSteps>0) {
    
    ifstream BOUNDARYMOVE(BoundaryMoveFILE);
    
    for(BoundaryMoveStep=0; BoundaryMoveStep<BoundaryMoveSteps; BoundaryMoveStep++){
      //BOUNDARYMOVE >> bi; BoundaryMoves[BoundaryMoveStep].Number = bi;
      for(bi=0; bi<D; bi++) {
        for(nu=0; nu<D; nu++) {
          BOUNDARYMOVE >>  BoundaryMoves[BoundaryMoveStep].vel[bi][nu];
          if( (bi!=nu) && (BoundaryMoves[BoundaryMoveStep].vel[bi][nu] != 0.)) 
            LeesEdwardsTrue = 1;
        }
      }
      BOUNDARYMOVE >> BoundaryMoves[BoundaryMoveStep].time;

      if(BOUNDARYMOVE.peek()!=13 && BOUNDARYMOVE.peek()!=10 && !BOUNDARYMOVE.eof()) {
        cout << "Boundary moves have the wrong format";
        exit(1);
      }
    }
    BOUNDARYMOVE.close();
    
    /*if (LeesEdwardsTrue)
      cout << "using Lees-Edwards boundary conditions" << endl;
    else
      cout << "NOT using Lees-Edwards boundary conditions" << endl;*/

    
    qsort (BoundaryMoves, BoundaryMoveSteps, sizeof(BoundaryMoveStruct), compBM);
    
    
    cout << "Sorted Boundary Moves:" << endl;
    for(BoundaryMoveStep=0; BoundaryMoveStep<BoundaryMoveSteps; BoundaryMoveStep++) {
      cout << "time: " << BoundaryMoves[BoundaryMoveStep].time << ":" << endl;
      PrintMatrix(BoundaryMoves[BoundaryMoveStep].vel);
      /*cout << BoundaryMoves[BoundaryMoveStep].Number << " ";
      for(nu=0;nu<D;nu++)
        cout << BoundaryMoves[BoundaryMoveStep].direction[nu] << " ";*/
      cout << endl;
    }
  }
  
       // ******************* read boundary oscillations ******************* 
    
  BoundaryOsciTrue  = 0;  
  BoundaryOsciSteps = CheckNumberOfLines(BoundaryOsciFILE);
  
  if (BoundaryOsciSteps==0) 
    cout << "Warning: BoundaryOsciFile " << BoundaryOsciFILE << " empty. Hence, no boundary oscillations." << endl;

  if (BoundaryOsciSteps<0) {
    cout << "Warning: BoundaryOsciFile " << BoundaryOsciFILE << " does not exist. Hence, no boundary oscillations." << endl;
    BoundaryOsciSteps = 0;
  }
  
  BoundaryOsciStruct BoundaryOscis[BoundaryOsciSteps];
  if (BoundaryOsciSteps>0) {
    BoundaryOsciTrue = 1;
    ifstream BOUNDARYOSCI(BoundaryOsciFILE);
    
    for(BoundaryOsciStep=0; BoundaryOsciStep<BoundaryOsciSteps; BoundaryOsciStep++){
      
      for(bi=0; bi<D; bi++) {
        for(nu=0; nu<D; nu++) {
          BOUNDARYOSCI >>  BoundaryOscis[BoundaryOsciStep].amplitude[bi][nu];
          if( (bi!=nu) && (BoundaryOscis[BoundaryOsciStep].amplitude[bi][nu] != 0.)) 
            LeesEdwardsTrue = 1;
        }
      }
      for(bi=0; bi<D; bi++) 
        for(nu=0; nu<D; nu++) 
          BOUNDARYOSCI >>  BoundaryOscis[BoundaryOsciStep].frequency[bi][nu];
      
      for(bi=0; bi<D; bi++) 
        for(nu=0; nu<D; nu++) 
          BOUNDARYOSCI >>  BoundaryOscis[BoundaryOsciStep].timeShift[bi][nu];

      BOUNDARYOSCI >> BoundaryOscis[BoundaryOsciStep].time;

      if(BOUNDARYOSCI.peek()!=13 && BOUNDARYOSCI.peek()!=10 && !BOUNDARYOSCI.eof()) {
        cout << "Boundary oscillations have the wrong format";
        exit(1);
      }
    }
    BOUNDARYOSCI.close();
    
    
    qsort (BoundaryOscis, BoundaryOsciSteps, sizeof(BoundaryOsciStruct), compBO);
    
    
    cout << "Sorted Boundary Oscillations:" << endl;
    for(BoundaryOsciStep=0; BoundaryOsciStep<BoundaryOsciSteps; BoundaryOsciStep++) {
      cout << "time: " << BoundaryOscis[BoundaryOsciStep].time << ":" << endl;
      PrintOsciMatrix(BoundaryOscis[BoundaryOsciStep]);
      cout << endl;
    }
  }

  // ******************* finished reading ******************* 
  
  if (LeesEdwardsTrue)
    cout << "using Lees-Edwards boundary conditions" << endl;
  else
    cout << "NOT using Lees-Edwards boundary conditions" << endl;

    //Create box mapping
  cout << "--- Doing boxes init ---\n";
    init_boxes(L,(box*)boxmap,box_size,D,X,N);
  
  // ******************* start the simulation ******************* 
  count = 0;
    boost::mt19937 rng(0); // I don't seed it on purpouse (it's not relevant)
  t = (floats)count*dt;
  
  PosOutCount    = 0;
  ForceOutCount  = 0;
  nextPosOut_t   = (floats)PosOutCount   * dt_pos;
  nextForceOut_t = (floats)ForceOutCount * dt_force;

  WallMoveStep = 0;
  BoundaryMoveStep = 0;
  BoundaryOsciStep = 0;
  
  // ******************* open streams for writing ******************* 
  ios_base::openmode mode = ios::out;
  if(continueAt >= 0) 
    mode = ios::in|ios::out;
  
  fstream BoundaryPosOutStream  (BoundaryPosOutFILE,  mode);
  fstream BoundaryForceOutStream(BoundaryForceOutFILE,mode);
  fstream WallPosOutStream      (WallPosOutFILE,      mode);
  fstream WallForceOutStream    (WallForceOutFILE,    mode);
  
  if(continueAt >= 0) {
    cout << "skipping dynamics until PosOutCount = " << continueAt << endl;
    while(PosOutCount <= continueAt) {

      if( t >= nextPosOut_t ) {
	if( PosOutCount >= continueAt ) 
	  break;
	PosOutCount++;  nextPosOut_t = (floats)PosOutCount * dt_pos;
	cout << " -- time: " << t << "\tcount: " << count << "\tnew PosOutCount: " << PosOutCount << "\tnextPosOut_t: " << nextPosOut_t << endl;
      }

      if( t >= nextForceOut_t ) {
	BoundaryPosOutStream.getline  (lineIN,LINELENGTH);
	//cout << "--- " << lineIN << endl;
	//cout << BoundaryPosOutStream.tellg() << endl;
	BoundaryForceOutStream.getline(lineIN,LINELENGTH);
	//cout << "--- " << lineIN << endl;
	//cout << BoundaryForceOutStream.tellg() << endl;
	WallPosOutStream.getline      (lineIN,LINELENGTH);
	WallForceOutStream.getline    (lineIN,LINELENGTH);

	ForceOutCount++;  nextForceOut_t = (floats)ForceOutCount * dt_force;
      }

      /******************************************\
      | update velocities of WALLS, if necessary |
      \******************************************/
      if (WallN>0) {
	while(WallMoveStep<WallMoveSteps && t>=WallMoves[WallMoveStep].time) {
	  wi = WallMoves[WallMoveStep].WallNumber;
	  
	  for(nu=0;nu<D;nu++)
	    vWall[nu] = WallMoves[WallMoveStep].direction[nu];
	  
	  for(j=0;j<Wall_AtomCount[wi];j++) {
	    kD = Wall_AtomIndex[wi][j]*D;
	    V[kD+nu] = vWall[nu];
	  }
	  
	  cout << "changing vel of wall " << wi << " to v = (";
	  cout << vWall[0];
	  for(nu=1; nu<D; nu++)
	    cout << "," << vWall[nu];
	  cout << ")" <<endl;
	  WallMoveStep++;
	}
      }
      
      /***********************************************\
      | update velocities of boundaries, if necessary |
      \***********************************************/
      if (BoundaryMoveSteps>0) {
	while(BoundaryMoveStep<BoundaryMoveSteps && t>=BoundaryMoves[BoundaryMoveStep].time) {
	  //bi = BoundaryMoves[BoundaryMoveStep].Number;
	  
	  for(bi=0; bi<D; bi++)
	    for(nu=0; nu<D; nu++)
	      vBoundary[bi][nu] = BoundaryMoves[BoundaryMoveStep].vel[bi][nu];
	  

	  cout << "changing boundary velocity to: " << endl;
	  PrintMatrix(vBoundary); cout << endl;
	  BoundaryMoveStep++;
	}
      }
      
      /*************************************************\
      | update oscillations of boundaries, if necessary |
      \*************************************************/
      if (BoundaryOsciSteps>0) {
	while(BoundaryOsciStep<BoundaryOsciSteps && t>=BoundaryOscis[BoundaryOsciStep].time) {
	  
	  for(bi=0; bi<D; bi++)
	    for(nu=0; nu<D; nu++)
	      BoundaryOsci = BoundaryOscis[BoundaryMoveStep];
	  

	  cout << "changing boundary oscillation velocity to: " << endl;
	  PrintOsciMatrix(BoundaryOsci); cout << endl;
	  BoundaryOsciStep++;
	}
      }
      

      
      // update position of boundaries
      if(BoundaryMoveSteps > 0) 
	for(bi=0; bi<D; bi++) 
	  for(nu=0; nu<D; nu++) 
	    L_LE[bi][nu] += dt * vBoundary[bi][nu];
      
      if(BoundaryOsciSteps > 0) 
	for(bi=0; bi<D; bi++) 
	  for(nu=0; nu<D; nu++) 
	    L_LE[bi][nu] += dt * BoundaryOsci.amplitude[bi][nu] * cos( twoPi * BoundaryOsci.frequency[bi][nu] * (t - BoundaryOsci.timeShift[bi][nu]) );
      
      for(bi=0; bi < D; bi++) {
	L[bi]  = L_LE[bi][bi];
	Lh[bi] = L[bi]*0.5;
      }

      // update position of Lees Edwards boundaries
      if (LeesEdwardsTrue) {
	Inverse( L_LE, L_LEinv, D ); 
      }


      count++;
      t = (floats)count*dt;
    }
    
    cout << endl;
    cout << "skipping done... continuing the dynamics at: ";
    cout << "time t = " << t << "\t count = " << count << endl;
    //cout << BoundaryPosOutStream.tellg() << " " << BoundaryForceOutStream.tellg() << " " << WallPosOutStream.tellg() << " " << WallForceOutStream.tellg() << " " <<endl;
    cout << endl;
  }





  short firstStep = 1;
  while(t<=TMAX) {
    /****************************************************************************************\
    | write particle positions, velocities, connections and local stress tensor if specified |
    \****************************************************************************************/
    if( t >= nextPosOut_t ) {
      if (preventOverwriteQ == 1) {
	cout << "Not overwriting output file " << addNDot(PosOutFILE,PosOutCount).c_str() << "." << endl;
	preventOverwriteQ = 0;
      }
      else {
	ofstream PosOutStream(addNDot(PosOutFILE,PosOutCount).c_str(),ios::out);
	ofstream VelOutStream(addNDot(VelOutFILE,PosOutCount).c_str(),ios::out);
	PosOutStream << setprecision(OUTPRECISION);
	VelOutStream << setprecision(OUTPRECISION);
	
	for(k=0;k<N;k++) {
	  kD=k*D;
	  PosOutStream << X[kD];
	  VelOutStream << V[kD];
	  for(nu=1;nu<D;nu++) {
	    PosOutStream << "\t" << X[kD+nu];
	    VelOutStream << "\t" << V[kD+nu];
	  }
	  PosOutStream << endl;
	  VelOutStream << endl;
	}


	PosOutStream.close();
	VelOutStream.close();
	
	
	
	// write connection file
	ofstream ConnOutStream(addNDot(ConnOutFILE,PosOutCount).c_str(),ios::out);
	ConnOutStream << setprecision(OUTPRECISION);
	
	for(it = connections.begin(); it != connections.end(); ++it) {
	  conn = *it;
	  conn.i1++;
	  conn.i2++;
	  ConnOutStream << conn.i1;
	  ConnOutStream << "\t" << conn.i2;
	  ConnOutStream << "\t" << conn.type;
	  ConnOutStream << "\t" << conn.a;
	  ConnOutStream << "\t" << conn.dx0;
	  ConnOutStream << "\t" << conn.damp;
	  ConnOutStream << "\t" << conn.c;
	  ConnOutStream << "\t" << conn.c1;
	  ConnOutStream << "\t" << conn.c2;
	  if(conn.dxbb1 >= 0. || conn.dxbb2 >= 0.) {
	    ConnOutStream << "\t" << conn.dxbb1;
	    ConnOutStream << "\t" << conn.dxbb2;
	  }
	  ConnOutStream << endl;
	}
	ConnOutStream.close();
	
	
	// write particle velocities relative to their neighbors
	if( strcmp(VelFluctOutFILE, "none") != 0 ) {

	  for(k=0; k<ND; k++)
	    vNeighborsSum[k] = 0.;

	  ofstream VelFluctOutStream(addNDot(VelFluctOutFILE,PosOutCount).c_str(),ios::out);
	  for(it = connections.begin(); it != connections.end(); ++it) {
	    conn = *it;
	    conn_i1 = conn.i1; conn_i1D = conn_i1*D;
	    conn_i2 = conn.i2; conn_i2D = conn_i2*D;
	    NNeighors[conn_i1]++;
	    NNeighors[conn_i2]++;
	    for(nu=0; nu<D; nu++) {
	      vNeighborsSum[conn_i1D+nu] += V[conn_i2D+nu];
	      vNeighborsSum[conn_i2D+nu] += V[conn_i1D+nu];
	    }
	  }
	  
	  for(k=0; k<N; k++) {
	    kD=k*D;
	    for(nu=0; nu<D; nu++) 
	      vField[nu] = V[kD+nu] - (vNeighborsSum[kD+nu]+V[kD+nu])/(floats)(NNeighors[k]+1);


	    VelFluctOutStream << vField[0];
	    for(nu=1; nu<D; nu++) 
	      VelFluctOutStream << "\t" << vField[nu];

	    VelFluctOutStream << endl;
	  }
	  VelFluctOutStream.close();
	}
	
	// write local stress tensor, if filename specified
	if(strcmp(SigmaOutFILE, "none") != 0) {

	  for(k=0; k<N; k++) {
	    for(nu=0; nu<D; nu++) {
	      for(bi=0; bi<D; bi++) {
		    sigmaLocalElas[k][bi][nu] = 0.;
		    sigmaLocalDiss[k][bi][nu] = 0.;
	      }
	    }
	  }

	  for(it = connections.begin(); it != connections.end(); ++it) {
	    conn = *it;
	    conn_i1  = conn.i1;
	    conn_i2  = conn.i2;
	    conn_i1D = conn.i1*D;
	    conn_i2D = conn.i2*D;
	  
	    CALC_dX0(X, conn_i1D, conn_i2D, dX);
	    CALC_dX(dX, dX_PeriodicCount);
	    dX2 = CALC_SQR_VEC(dX);
    
	    // calculate contribution of central force
	    Fh = get_shielded_force(force(conn.type, dX2, conn.a, conn.b),conn.i1,conn.i2,X,(box*)boxmap,lc,D,L,cell_radius,shield_matrix,shielding);
	    for(nu=0; nu<D; nu++) {
	      Fnu = dX[nu]*Fh;
	      for(bi=0; bi<D; bi++) {
		sigmaLocalElas[conn_i1][bi][nu] += Fnu * dX[bi];
		sigmaLocalElas[conn_i2][bi][nu] += Fnu * dX[bi];
	      }
	    }
	  
	    // calculate contribution of damping force
	    conn_damp = conn.damp;
	    if(conn_damp > 0) {
	      CALC_dV0(V, conn_i1D, conn_i2D, dV);
	      CALC_dV(dV, dX_PeriodicCount);
	  
	      if(conn_damp == 1) {
		Fh = - conn.c;
	      }
	      else {
		dV2 = CALC_SQR_VEC(dV);
		Fh = force_damp(conn_damp, dV2, conn.c);
	      }
	  
	      if (conn.c1 != 0.) {
		c1Fh = conn.c1*Fh;
		for(nu=0; nu<D; nu++) {
		  Fnu = dV[nu]*c1Fh;
		  for(bi=0; bi<D; bi++) {
		    sigmaLocalDiss[conn_i1][bi][nu] += Fnu * dX[bi];
		    sigmaLocalDiss[conn_i2][bi][nu] += Fnu * dX[bi];
		  }
		}
	      }
	  
	      if (conn.c2 != 0.) {
		dXdV = dX[0]*dV[0];
		for(nu=1; nu<D; nu++)
		  dXdV += dX[nu]*dV[nu];
	  
		c2Fh = conn.c2 * dXdV / dX2;
		for(nu=0; nu<D; nu++) {
		  Fnu = dX[nu]*c2Fh;
		  for(bi=0; bi<D; bi++) {
		    sigmaLocalDiss[conn_i1][bi][nu] += Fnu * dX[bi];
		    sigmaLocalDiss[conn_i2][bi][nu] += Fnu * dX[bi];
		  }
		}
	      }
	    }
	  }
	  
	  // start writing
	  ofstream SigmaOutStream(addNDot(SigmaOutFILE,PosOutCount).c_str(),ios::out);
	  
	  for(k=0; k<N; k++) {

	    for (bi=0; bi<D; bi++) 
	      for (nu=0; nu<D; nu++) {
		if(bi==0 && nu==0) SigmaOutStream << 0.5*sigmaLocalElas[k][bi][nu];
		else 		 SigmaOutStream << "\t" << 0.5*sigmaLocalElas[k][bi][nu];
	      }
	
	    for (bi=0; bi<D; bi++) 
	      for (nu=0; nu<D; nu++) 
		SigmaOutStream << "\t" << 0.5*sigmaLocalDiss[k][bi][nu];

	    SigmaOutStream << endl;
	  }
	  SigmaOutStream.close();
	}
      }
      
      PosOutCount++;  nextPosOut_t = (floats)PosOutCount * dt_pos;
      cout << " -- time: " << t << "\tcount: " << count << "\tnew PosOutCount: " << PosOutCount << "\tnextPosOut_t: " << nextPosOut_t << endl;

      //cout << "L:" << endl;         PrintMatrix(L_LE);    ;
      //cout << "inverse L:" << endl; PrintMatrix(L_LEinv); cout << endl;

    }

    
    /*****************************************************\
    | calculate and write forces on Walls and boundaries  |
    \*****************************************************/
    
    /** 
    | note on indices:
    | BOUNDARIES MOVEMENT:
    | L_LE[nu][bi] is the velocity of the "bi-boundary" (boundary with normal in bi-direction) in nu-direction
    | Keep in mind: the first index is the direction, the second index the boundary.
    | e.g. L_LE[0][1] = 1. is a shear of the top-boundary one length-unit to the right.
    | By Matrix multiplication, L_LE transforms the unit square to the actual simulation box, such that
    | " (X,Y) = L_LE . (Xinv, Yinv) " yields a position (X,Y) in the simulation box, for 0 <= Xinv,Yiny < 1.
    | 
    | BOUNDARY FORCES: 
    | BoundaryForceElas[bi][nu] is the elastic contribution to the stress tensor
    | BoundaryForceDiss[bi][nu] is the dissipative contribution to the stress tensor
    | The sum of both are the actual forces
    | BoundaryForce____[bi][nu] contains the "raw" forces. 
    | To obtain the stress tensor, the forces must be normalized with the system size by hand.
    | 
    | BoundaryForce____[bi][nu] is the force on the bi-boundary in nu-direction.
    |
    | note the switched indices compared to L_LE.
    \ **/
    
    if( t >= nextForceOut_t ) {
      // forces on Walls
      //cout << " " << L[0] << "  " << L[1] << endl;
      for(wi=1;wi<=WallN;wi++) {
        wiD = wi*D;
        for(nu=0;nu<D;nu++)
          Wall_Force[wiD+nu] = 0.;
	
        for(j=0;j<Wall_AtomCount[wi];j++) {
          kD = Wall_AtomIndex[wi][j]*D;
          for(nu=0;nu<D;nu++) {
            kDpnu = kD+nu;
            Wall_Force[wiD+nu] += 0.5*(R1_F[kDpnu] + R2_F[kDpnu]);
          }
        }
	
        for(nu=0;nu<D;nu++)
          WallForceOutStream << Wall_Force[wiD+nu] << "\t";
	
      }
      WallForceOutStream << t << endl;
      
      
      // write position of first wall particle
      for(wi=1; wi<=WallN; wi++) {
        kD = Wall_AtomIndex[wi][0]*D;
        for(nu=0; nu<D; nu++)
          WallPosOutStream << X[kD+nu] << "\t";
      }
      WallPosOutStream << t << endl;
      
      
      // calculate forces on boundaries and in the whole system
      for(bi=0; bi<D; bi++) {
        for(nu=0; nu<D; nu++) {
          BoundaryForceElas[bi][nu] = 0.;
          BoundaryForceDiss[bi][nu] = 0.;
	  SystemForceElas[bi][nu]   = 0.;
          SystemForceDiss[bi][nu]   = 0.;
        }
      }
      
      
      
      
      for(it = connections.begin(); it != connections.end(); ++it) {
        conn = *it;
        conn_i1D = conn.i1*D;
        conn_i2D = conn.i2*D;
	
	CALC_Xinv(X+conn_i1D, X1_inv);
        CALC_Xinv(X+conn_i2D, X2_inv);
	
	for(bi=0; bi<D; bi++) {
	  dX_inv = X1_inv[bi] - X2_inv[bi];
	  if(dX_inv < -0.5)
	    BoundaryConnQ[bi] = -1;
	  else if(dX_inv > 0.5)
	    BoundaryConnQ[bi] = 1;
	  else
	    BoundaryConnQ[bi] = 0;
	}
	
	
	CALC_dX0(X, conn_i1D, conn_i2D, dX);
	CALC_dX(dX, dX_PeriodicCount);
	dX2 = CALC_SQR_VEC(dX);
        
	// calculate contribution of central force
	Fh = get_shielded_force(force(conn.type, dX2, conn.a, conn.b),conn.i1,conn.i2,X,(box*)boxmap,lc,D,L,cell_radius,shield_matrix,shielding);

        for(nu=0; nu<D; nu++) {
          Fnu = dX[nu]*Fh;

	  for(bi=0; bi<D; bi++) {
	    // calculate tensor product contribution [  F (x) (R_i - R_j)  ],
	    // where "(x)" denotes the tensor product
	    SystemForceElas[bi][nu] += Fnu * dX[bi];

	    // add force contribution, if connection crosses the boundary in bi-direction
	    ADD_BOUNDARYFORCE(BoundaryForceElas, bi, nu, Fnu);
	  }
        }


        // calculate contribution of damping force
	conn_damp = conn.damp;
	if(conn_damp > 0) {
	  CALC_dV0(V, conn_i1D, conn_i2D, dV);
	  CALC_dV(dV, dX_PeriodicCount);
    
	  if(conn_damp == 1) {
	    Fh = - conn.c;
	  }
	  else {
	    dV2 = CALC_SQR_VEC(dV);
	    Fh = force_damp(conn_damp, dV2, conn.c);
	  }
    
	  if (conn.c1 != 0.) {
	    c1Fh = conn.c1*Fh;
	    for(nu=0; nu<D; nu++) {
	      Fnu = dV[nu]*c1Fh;
	      for(bi=0; bi<D; bi++) {
		SystemForceDiss[bi][nu] += Fnu * dX[bi];
		ADD_BOUNDARYFORCE(BoundaryForceDiss, bi, nu, Fnu);
	      }
	    }
	  }
    
	  if (conn.c2 != 0.) {
	    dXdV = dX[0]*dV[0];
	    for(nu=1; nu<D; nu++)
	      dXdV += dX[nu]*dV[nu];
      
	    c2Fh = conn.c2 * dXdV / dX2;
	    for(nu=0; nu<D; nu++) {
	      Fnu = dX[nu]*c2Fh;
	      for(bi=0; bi<D; bi++) {
		SystemForceDiss[bi][nu] += Fnu * dX[bi];
		ADD_BOUNDARYFORCE(BoundaryForceDiss, bi, nu, Fnu);
	      }
	    }
	  }
	}

      }
      // write forces on boundaries
      
      for (bi=0; bi<D; bi++) 
        for (nu=0; nu<D; nu++) 
          BoundaryForceOutStream << BoundaryForceElas[bi][nu] << "\t";
      
      for (bi=0; bi<D; bi++) 
        for (nu=0; nu<D; nu++) 
          BoundaryForceOutStream << BoundaryForceDiss[bi][nu] << "\t";

      for (bi=0; bi<D; bi++) 
        for (nu=0; nu<D; nu++) 
          BoundaryForceOutStream << SystemForceElas[bi][nu] << "\t";
      
      for (bi=0; bi<D; bi++) 
        for (nu=0; nu<D; nu++) 
          BoundaryForceOutStream << SystemForceDiss[bi][nu] << "\t";

	
      BoundaryForceOutStream << t << endl;

      // write boundary position
      for (bi=0; bi<D; bi++) {
        for (nu=0; nu<D; nu++) {
          if (bi == nu) 
            BoundaryPosOutStream << L[nu] << "\t";
          else
            BoundaryPosOutStream << L_LE[bi][nu] << "\t";
        }
      }
      BoundaryPosOutStream << t << endl;
      
      ForceOutCount++;  nextForceOut_t = (floats)ForceOutCount * dt_force;
    }
    
    
    /******************************************\
    | update velocities of WALLS, if necessary |
    \******************************************/
    if (WallN>0) {
      while(WallMoveStep<WallMoveSteps && t>=WallMoves[WallMoveStep].time) {
	wi = WallMoves[WallMoveStep].WallNumber;
	
	for(nu=0;nu<D;nu++)
	  vWall[nu] = WallMoves[WallMoveStep].direction[nu];
	
	for(j=0;j<Wall_AtomCount[wi];j++) {
	  kD = Wall_AtomIndex[wi][j]*D;
	  V[kD+nu] = vWall[nu];
	}
	
	cout << "changing vel of wall " << wi << " to v = (";
	cout << vWall[0];
	for(nu=1; nu<D; nu++)
	  cout << "," << vWall[nu];
	cout << ")" <<endl;
	WallMoveStep++;
      }
    }
    
    /***********************************************\
    | update velocities of boundaries, if necessary |
    \***********************************************/
    if (BoundaryMoveSteps>0) {
      while(BoundaryMoveStep<BoundaryMoveSteps && t>=BoundaryMoves[BoundaryMoveStep].time) {
        //bi = BoundaryMoves[BoundaryMoveStep].Number;
	
        for(bi=0; bi<D; bi++)
          for(nu=0; nu<D; nu++)
            vBoundary[bi][nu] = BoundaryMoves[BoundaryMoveStep].vel[bi][nu];
	

        cout << "changing boundary velocity to: " << endl;
        PrintMatrix(vBoundary); cout << endl;
        BoundaryMoveStep++;
      }
    }
    
    /*************************************************\
    | update oscillations of boundaries, if necessary |
    \*************************************************/
    if (BoundaryOsciSteps>0) {
      while(BoundaryOsciStep<BoundaryOsciSteps && t>=BoundaryOscis[BoundaryOsciStep].time) {
	
        for(bi=0; bi<D; bi++)
          for(nu=0; nu<D; nu++)
            BoundaryOsci = BoundaryOscis[BoundaryMoveStep];
	

        cout << "changing boundary oscillation velocity to: " << endl;
        PrintOsciMatrix(BoundaryOsci); cout << endl;
        BoundaryOsciStep++;
      }
    }
    
    /******************************************\
    | update positions of walls and boundaries |
    \******************************************/
    
    // update positions of walls
    for(wi=1; wi<=WallN; wi++) {
      for(j=0; j<Wall_AtomCount[wi]; j++) {
        kD = Wall_AtomIndex[wi][j]*D;
        for(nu=0; nu<D; nu++) {
          kDpnu = kD+nu;
          X[kDpnu] += dt*V[kDpnu];
        }
      }
    }
    
    // update position of boundaries
    if(BoundaryMoveSteps > 0) 
      for(bi=0; bi<D; bi++) 
        for(nu=0; nu<D; nu++) 
          L_LE[bi][nu] += dt * vBoundary[bi][nu];
    
    if(BoundaryOsciSteps > 0) 
      for(bi=0; bi<D; bi++) 
        for(nu=0; nu<D; nu++) 
          L_LE[bi][nu] += dt * BoundaryOsci.amplitude[bi][nu] * cos( twoPi * BoundaryOsci.frequency[bi][nu] * (t - BoundaryOsci.timeShift[bi][nu]) );
    
    for(bi=0; bi < D; bi++) {
      L[bi]  = L_LE[bi][bi];
      Lh[bi] = L[bi]*0.5;
    }

     // update position of Lees Edwards boundaries
    if (LeesEdwardsTrue) {
      Inverse( L_LE, L_LEinv, D ); 
      //cout << "L:" << endl;         PrintMatrix(L_LE);    cout << endl;
      //cout << "inverse L:" << endl; PrintMatrix(L_LEinv); cout << endl;
    }


    /*****************\
    | do the dynamics |
    \*****************/
    // calculate random kicks
    if(SIGMA>0.) {
      for(j=0; j<Wall_AtomCount0; j++) {
        kD = Wall_AtomIndex0[j]*D;
        for(nu=0; nu<D; nu++) 
          Vrand[kD+nu] = randomKick();
      }
    }
    seconds=time(NULL);
    
    
    
    // choose the method
    if( strcmp(IntegrationMethod, "LeapfrogNoDamp") == 0 )  { /* Leapfrog */
      // calculate R1
      ADD_COMPONENTS(R1_X[kDpnu] = X[kDpnu] + dt_h * V[kDpnu];);
      
      calc_F(R1_X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R1_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      
      // update positions of movable atoms
      ADD_COMPONENTS( V[kDpnu] = ONEminusGAMMAdt * V[kDpnu]  +  mi[k] * (dt*R1_F[kDpnu] + SIGMA*Vrand[kDpnu]);
                      X[kDpnu] = R1_X[kDpnu] + dt_h * V[kDpnu]; );

    }
    else if( strcmp(IntegrationMethod, "Leapfrog") == 0 )  { /* Leapfrog with damping*/
      if(firstStep == 1) {
        ADD_COMPONENTS( previous_X[kDpnu] = X[kDpnu]; );
	firstStep = 0;
      }
      
      // calculate R1
      calc_F(X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R1_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      
      // calculate R2
      ADD_COMPONENTS(R1_V[kDpnu] = V[kDpnu] + mi[k]*dt*R1_F[kDpnu];
                     R2_V[kDpnu] = 0.5*( (X[kDpnu] - previous_X[kDpnu])/dt   +   R1_V[kDpnu] ); );
      
      calc_F(X, R2_V, Wall_AtomCount0, Wall_AtomIndex0, connections, R2_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      
      // remember current positions for next step and 
      // update positions of movable atoms
      ADD_COMPONENTS( previous_X[kDpnu] = X[kDpnu];
                      V[kDpnu] = ONEminusGAMMAdt * V[kDpnu]  +  mi[k] * (dt*R2_F[kDpnu] + SIGMA*Vrand[kDpnu]);
                      X[kDpnu] = X[kDpnu] + dt*V[kDpnu]; );
      

    }
    else if( strcmp(IntegrationMethod, "VelocityVerlet") == 0 )  { /* Velocity Verlet with efficient damping */
      /* PSEUDO CODE: see http://wiki.vdrift.net/Old_Numerical_Integration
      
      if (not oldaccel)
        oldaccel = acceleration(state, t+dt)

      x += v*dt + 0.5*oldaccel*dt*dt
      v += 0.5*oldaccel*dt
      a = acceleration(state, t+dt)
      v += 0.5*a*dt

      oldaccel = a
      */
      
      if (count==0) {
        calc_F(X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R1_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      }
      
      // for count>0, R1_F is remembered from the previous tmie step
      ADD_COMPONENTS( //VplusAccDtH = V[kDpnu]                 + mi[k] *  dt_h*R1_F[kDpnu];
                      VplusAccDtH = ONEminusGAMMAdt*V[kDpnu]  +  mi[k] * (dt_h*R1_F[kDpnu] + SIGMA*Vrand[kDpnu]);
                      X[kDpnu] = X[kDpnu]  +  VplusAccDtH * dt;
                      V[kDpnu] =              VplusAccDtH;       );
      
      calc_F(X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R2_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      
      
      // update velocities of movable atoms and
      // remember current forces for next step 
      ADD_COMPONENTS( V[kDpnu] = ONEminusGAMMAdt*V[kDpnu]  +  mi[k] * (dt_h*R2_F[kDpnu] + SIGMA*Vrand[kDpnu]); 
                      R1_F[kDpnu] = R2_F[kDpnu]; );


    }
    else if( strcmp(IntegrationMethod, "Euler") == 0 )  {  /* Euler integration */
      // calculate R1 
      calc_F(X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R1_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      
      // update positions of movable atoms
      ADD_COMPONENTS( X[kDpnu] = X[kDpnu] + dt*V[kDpnu];
                      V[kDpnu] = ONEminusGAMMAdt * V[kDpnu]  +  mi[k] * (dt*R1_F[kDpnu] + SIGMA*Vrand[kDpnu]); );
      
    }
    else if( strcmp(IntegrationMethod, "Overdamped") == 0 ){  /* Overdamped Heun method */
      
      // calculate R1 using X and V
      calc_F(X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R1_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      ADD_COMPONENTS( R1_V[kDpnu] = GAMMAi * R1_F[kDpnu];
		      R1_X[kDpnu] = X[kDpnu] + dt*R1_V[kDpnu] + SIGMAoGAMMA*Vrand[kDpnu]; );
      
      // calculate R2 using R1
      calc_F(R1_X, R1_V, Wall_AtomCount0, Wall_AtomIndex0, connections, R2_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      ADD_COMPONENTS( R2_V[kDpnu] = GAMMAi * R2_F[kDpnu];
		      /*R2_X[kDpnu] = X[kDpnu] + dth*( R1_V[kDpnu] + R2_V[kDpnu]) + SIGMA*Vrand[kDpnu]; */);

      
      // update positions of movable atoms
      ADD_COMPONENTS( V[kDpnu] = 0.5*(R1_V[kDpnu] + R2_V[kDpnu]);
                      X[kDpnu] = X[kDpnu] + dt*V[kDpnu] + SIGMAoGAMMA*Vrand[kDpnu]; );

    }
    else {  /* Default algorithm: Heun method */
      
      // calculate R1 using X and V
      calc_F(X, V, Wall_AtomCount0, Wall_AtomIndex0, connections, R1_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      ADD_COMPONENTS( R1_X[kDpnu] = X[kDpnu] + dt_h*V[kDpnu];
                      R1_V[kDpnu] = ONEminusGAMMAdt_h * V[kDpnu]  +  mi[k] * (dt_h*R1_F[kDpnu] + SIGMA*Vrand[kDpnu]); );
      
      // calculate R2 using R1
      calc_F(R1_X, R1_V, Wall_AtomCount0, Wall_AtomIndex0, connections, R2_F,(box*)boxmap,lc,cell_radius,shield_matrix,shielding,box_size,cellr2);
      ADD_COMPONENTS( R2_X[kDpnu] = R1_X[kDpnu] + dt_h*R1_V[kDpnu];
                      R2_V[kDpnu] = ONEminusGAMMAdt_h * R1_V[kDpnu]  +  mi[k] * (dt_h*R2_F[kDpnu] + SIGMA*Vrand[kDpnu]); );
      
      //code fragments added by Lutz to keep track of work delivered by single
      //connections. Remark this only works on default integration method
      for (int erta=0;erta<ND;erta++)
          diff_X[erta]=X[erta];
      // update positions of movable atoms
      ADD_COMPONENTS( X[kDpnu] = 0.5*(R1_X[kDpnu] + R2_X[kDpnu]);
                      V[kDpnu] = 0.5*(R1_V[kDpnu] + R2_V[kDpnu]); );
      for (int erta=0;erta<ND;erta++)
          diff_X[erta]=X[erta]-diff_X[erta];
      //keep track of the energy
      if (dt_energy_distr > 0 || dt_plot_work_hist > 0)
      {
      for (it=connections.begin();it!=connections.end();++it)
      {
		CALC_dX0(X,it->i1*D,it->i2*D,dX);
        CALC_dX(dX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(dX);
        work=get_shielded_force(force(it->type,dX2,it->a,it->b),it->i1,it->i2,X,(box*)boxmap,lc,D,L,cell_radius,shield_matrix,shielding);
        work=work*(dX[0]*diff_X[it->i1*D]+dX[1]*diff_X[it->i1*D+1])-work*(dX[0]*diff_X[it->i2*D]+dX[1]*diff_X[it->i2*D+1]);
        it->work+=work;
        //if (abs(work) <=work_max)
        //{
        //    work=work+work_max;
        //    place=round(work/dw);
        //    work_histogramm[place]++;
        //    if(work_histogramm[place]<0)
        //        cout << "we are in serious trouble, hist is " << work_histogramm[place] << " at " << place << endl;
        //}
      }
      }
    }
seconds=time(NULL)-seconds;
stats_time[0]+=(double)seconds;
seconds=time(NULL);
    /*****************************\
    | remove and add interactions |
    \*****************************/
    // remove interactions
    if (RemoveProb > 0.) {
      for(it = connections.begin(); it != connections.end(); ) {
        //if( (*it).type==1 ) {
          conn = *it;
          conn_i1D = conn.i1*D;
          conn_i2D = conn.i2*D;
          
          CALC_dX0(X, conn_i1D, conn_i2D, dX);
          CALC_dX(dX, dX_PeriodicCount);
          dX2 = CALC_SQR_VEC(dX);
          if( dX2 > dX2Remove || dX2 < conn.dx2bb1 || (conn.dx2bb2 > 0. && dX2 > conn.dx2bb2) || round(it->a) <= 0 || drand48() < RemoveProb)
	  {
          if(it->life_time>0)
          {
          place=dt*((double)it->life_time)/life_step;
          if(dt*it->life_time<=life_max)
          {
              life_time_distr[place]++;
              energy_distr[place]+=it->work;
          }
          norm++;
          }
        if (dt_plot_work_hist > 0 && abs(it->work) <=work_max)
        {
            work=it->work+work_max;
            place=round(work/dw);
            work_histogramm[place]++;
            if(work_histogramm[place]<0)
                cout << "we are in serious trouble, hist is " << work_histogramm[place] << " at " << place << endl;
        }
            it = connections.erase(it); 
	  }
          else
            ++it;
      }
    }
    else {
      // RemoveProb=0  --> no need to spend time on random number generation:
      for(it = connections.begin(); it != connections.end(); ) 	{
        conn = *it;
        conn_i1D = conn.i1*D;
        conn_i2D = conn.i2*D;
          
        CALC_dX0(X, conn_i1D, conn_i2D, dX);
        CALC_dX(dX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(dX);
        // Fh = force((*it).type, dX2, (*it).a, (*it).b);
        
        if(shield_matrix[conn.i1*N+conn.i2-get_nmnpeh(conn.i1)]==1 || dX2 > dX2Remove || dX2 < conn.dx2bb1 || (conn.dx2bb2 > 0. && dX2 > conn.dx2bb2) )
	{
        if(it->life_time>0)
        {
          place=dt*((double)it->life_time)/life_step;
          if(dt*it->life_time<=life_max)
          {
              life_time_distr[place]++;
              energy_distr[place]+=it->work;
          }
          norm++;
        }
        if (dt_plot_work_hist > 0 && abs(it->work) <=work_max)
        {
            work=it->work+work_max;
            place=round(work/dw);
            work_histogramm[place]++;
            if(work_histogramm[place]<0)
                cout << "we are in serious trouble, hist is " << work_histogramm[place] << " at " << place << endl;
        }
            it = connections.erase(it); 
	}
        else
          ++it;
  
      }
    }
    //increment life times
    if (dt_life_time > 0 || dt_energy_distr > 0)
    {
    for (it=connections.begin();it!=connections.end();++it)
    {
        it->life_time++;
    }
    }
seconds=time(NULL)-seconds;
stats_time[1]+=(double)seconds;
seconds=time(NULL);
    //update boxes
    for(int erta=0;erta<N;erta++)
        is_moved[erta]=0;
    for (int eta=0;eta<ceil(L[0]/box_size);eta++)
    {
        for (int jota=0;jota<ceil(L[1]/box_size);jota++)
        {
            for(int kappa=0;kappa<boxmap[eta][jota].get_size();kappa++)
            {
                i=boxmap[eta][jota].get_member();
                if(is_moved[i]==0)
                {
                if(i!=-1)
                {
                xpos=X[i*D]-floor(X[i*D]/L[0])*L[0];
                ypos=X[i*D+1]-floor(X[i*D+1]/L[1])*L[1];
                while (xpos >= L[0])
		{
                        xpos-=L[0];
			cout << xpos << endl;
		}
                while (xpos <0)
		{
                        xpos+=L[0];
			cout << xpos << endl;
		}
                while (ypos >= L[1])
		{
                        ypos-=L[1];
			cout << ypos << endl;
		}
                while (ypos<0)
		{
                        ypos+=L[1];
			cout << ypos << endl;
		}
                xpos=floor(xpos/box_size);
                ypos=floor(ypos/box_size);
                boxmap[(int)xpos][(int)ypos].add_member(i);
                is_moved[i]=1;
                }
                }
                else
                {
                    boxmap[eta][jota].add_member(i);
                }
            }
        }
    }
seconds=time(NULL)-seconds;
stats_time[2]+=(double)seconds;
seconds=time(NULL);
   // its_okay=false;
   // while(!its_okay)
   // {
   //     its_okay=true;
   // for (int erta=0;erta<D*N;erta++)
   //     Xhard[erta]=X[erta];
   // for (int erta=0;erta<N;erta++)
   //     not_moved[erta]=true;
   //         //hard core
   //     for (int jota=0;jota<N;jota++)
   //     {
   //         tempN=get_neigbours(jota,(box*)boxmap,X,box_size,D,L);
   //         for (int_it=tempN.begin();int_it!=tempN.end();++int_it)
   //         {
   //             i=(*int_it);
   //             if (jota>i)
   //             {
   //             CALC_dX0(X,i*D,jota*D,dX);
   //             CALC_dX(dX, dX_PeriodicCount);
   //             dX2 = CALC_SQR_VEC(dX);
   //             //dX2 = sqrt(dX2)/2.0;
   //             //if(dX2 < cellr2)
   //             //{
   //             excl_pow=pow(cellr2/dX2,2)*hard_scaling/dX2;
   //             for (int erta=0;erta<D;erta++)
   //             {
   //                 V[i*D+erta]=V[i*D+erta]+dX[erta]*excl_pow;
   //                 V[jota*D+erta]=V[jota*D+erta]-dX[erta]*excl_pow;
   //             }
                //}
                //if (dX2 < 0.99*cell_radius)
                //{
                //    its_okay=false;
                //    if (not_moved[i] && not_moved[jota])
                //    {
                //    not_moved[i]=false;
                //    not_moved[jota]=false;
                //    proj_1=0;
                //    proj_2=0;
                //    for (int kappa=0;kappa<D;kappa++)
                //    {
                //        proj_1+=V[jota*D+kappa]*dX[kappa]/dX2;
                //    }
                //    for (int kappa=0;kappa<D;kappa++)
                //    {
                //        proj_2+=V[i*D+kappa]*dX[kappa]/dX2;
                //    }
                //    for (int kappa=0;kappa<D;kappa++)
                //    {
                //        V[jota*D+kappa]=V[jota*D+kappa]-proj_1*dX[kappa]/dX2+proj_2*dX[kappa]/dX2;//in the case m1=m2 tangent vels are switched
                //        V[i*D+kappa]=V[i*D+kappa]-proj_2*dX[kappa]/dX2+proj_1*dX[kappa]/dX2;
                //        if (dX2==0)
                //        {
                //            cout << "fatal error, dX=0\n";
                //            return 1;
                //        }
                //        center=X[i*D+kappa];
                //        while (center < 0)
                //            center+=L[kappa];
                //        while (center>=L[kappa])
                //            center-=L[kappa];
                //        //center=(center+X[jota*D+kappa])/2.0;
                //        center=periodic_middle(center,X[jota*D+kappa],L[kappa]);
                //        center=center+floor(X[jota*D+kappa]/L[kappa])*L[kappa];
                //        Xhard[jota*D+kappa]=center-dX[kappa]/dX2*(2.0*cell_radius-dX2);
                //        center=X[jota*D+kappa];
                //        while (center < 0)
                //            center+=L[kappa];
                //        while (center>=L[kappa])
                //            center-=L[kappa];
                //        center=periodic_middle(center,X[i*D+kappa],L[kappa]);
                //        center=center+floor(X[i*D+kappa]/L[kappa])*L[kappa];
                //        //center=(center+X[i*D+kappa])/2.0;
                //        Xhard[i*D+kappa]=center+dX[kappa]/dX2*(2.0*cell_radius-dX2);
                //        //cout << center << " " << dX[kappa]/dX2*cell_radius << endl;
                //    }
                //    }
                //}
   //             }
   //         }
   //     }
   // for (int erta=0;erta<D*N;erta++)
   //     X[erta]=Xhard[erta];
   // }
    // add interactions
    for(i=0; i<AddNumber; i++) {
      j = (long)floor( drand48()*(double)N ); conn_i1D=j*D;
      k = (long)floor( drand48()*(double)N ); conn_i2D=k*D;
      if (j!=k) {
        
        CALC_dX0(X, conn_i1D, conn_i2D, dX);
        CALC_dX(dX, dX_PeriodicCount);
        dX2 = CALC_SQR_VEC(dX);
	
	if(dX2 < dX2Add) {
	  if(j<k) {
	    conn.i1 = j;
	    conn.i2 = k;
	  }
	  else {
    	    conn.i1 = k;
	    conn.i2 = j;
	  }
      if(shield_matrix[conn.i1*N+conn.i2-get_nmnpeh(conn.i1)]==0)
      {
	  conn.type = force_type;
          conn.damp = 0;
	  conn.a    = force_scaling;
          conn.dx0  = sqrt(dX2Add);
          dx0_to_b(conn.dx0, conn.b, conn.type);
          conn.c    = 0.;
          conn.c1   = 0.;
          conn.c2   = 0.;
          conn.dxbb1 = -1.;
          conn.dxbb2 = -1.;
          conn.life_time=0;
          conn.work=0;
          conn_stats[conn.i1]++;
          conn_stats[conn.i2]++;
	  connections.push_back(conn);
      }
      }
        }
	}
    //}
    //}
    //reset shield matrix
    if(do_force_on_shell)
        force_on_shell(V,kT,(double)D,N);
    for(int erta=0;erta<nmnmeh;erta++)
        shield_matrix[erta]=0;
seconds=time(NULL)-seconds;
stats_time[3]+=(double)seconds;
	//initialize neigbours
if (round(t) == t)
{
cout << "count: " << count << " Mean dyn time " << stats_time[0]/count << " Mean rem " << stats_time[1]/count << " Mean box " << stats_time[2]/count << " Mean add " << stats_time[3]/count << " with " << connections.size() << " connections" <<endl;
    cout << t << endl;
    if ((dt_op_simple != 0 && ((int)t)%((int)dt_op_simple) == 0)  || (dt_op != 0 && ((int)t)%((int)dt_op) == 0) || (dt_op_local!= 0 && ((int)t)%((int)dt_op_local) == 0 ) ||  (dt_plot_delaunay != 0 && ((int)t)%((int)dt_plot_delaunay) == 0))
    {
	neigbours.clear();
    //order of X is changed in dtris2
    for (int kappa=0;kappa<2*N;kappa++)
    {
        Xtemp[kappa]=X[kappa];
    }
    to_period(Xtemp,L,D,N,map);//this function does not really clean to period
    errorcode=dtris2(4*N,Xtemp,&Tri_num,Tri_vert,Tri_nabe);
    if(errorcode != 0)
    {
        cout << "errorcode " << errorcode << "  in delaunay triangulation, exit\n";
        return 1;
    }
    //now change the delaunay Tri_vert due to map. If we use some day Tri_nabe
    //change that too
    for(int eta=0;eta<3*Tri_num;eta++)
    {
        //every cell that should be ignored is set to -1
        if(Xtemp[D*(Tri_vert[eta]-1)]> (L[0]*1.25) || Xtemp[D*(Tri_vert[eta]-1)] < -0.25*L[0] || Xtemp[D*(Tri_vert[eta]-1)+1] > 1.25*L[1] || Xtemp[D*(Tri_vert[eta]-1)+1] < -0.25*L[1])
            Tri_vert[eta]=-1;
        else
            Tri_vert[eta]=map[Tri_vert[eta]-1];//particle indizes are one above array indices, here Tri_vert struct is changing
    }
    ////now find neighbours by going through Tri_vert
    drop_isin=0;
    for (int kappa=0;kappa<Tri_num;kappa++)
    {
        for (int u=0;u<3;u++)
            verts[u]=Tri_vert[3*kappa+u];
        for(int u=0;u<3;u++)
        {
            if(u<2)
                k=u+1;
            else
                k=0;
            if(verts[u]<verts[k])
            {
                conn.i1=verts[u];
                conn.i2=verts[k];
            }
            else
            {
                conn.i2=verts[u];
                conn.i1=verts[k];
            }
            contain=false;
            for(it=neigbours.begin();it!=neigbours.end();++it)
            {
                if (it->i1==conn.i1 && it->i2==conn.i2)//for we do index sanitization there is no need to exclude points violating the boundarys
                {
                    contain=true;//this is OC
                    drop_isin++;
                }
            }
            if(contain==false && conn.i1 != conn.i2 && conn.i1!=-1 && conn.i2 != -1)
            {
                neigbours.push_back(conn);
                //push also reverse
                tmp=conn.i1;
                conn.i1=conn.i2;
                conn.i2=tmp;
                neigbours.push_back(conn);
            }
        }
    }
}
    //perform some consistency checks
    //first every i->j, j->i should be contained twice
    div_neigbours++;
    //in first timestep determine neigbours of cell 0
    if(t==0)
    {
        tmp=round(sqrt(N)/2.0+N/2.0);
        marked.push_back(tmp);
        for (it=neigbours.begin();it!=neigbours.end();++it)
        {
            if (it->i1 == tmp)
            {
                marked.push_back(it->i2);
            }
        }
    }
    //cout << Tri_num << endl;
	if (dt_msd != 0 && ((int)t)%((int)dt_msd) == 0 )
	{
		get_msd(msd_val,X,X0,D,N,t,mean,dt_msd,dt);
		//cout << "starting msd\n";
		msd << t << " " << mean[0] << " " << mean[1] << endl;//Currently just round reporting times are supported. This should not cause a problem
	}
	if (dt_total_energy != 0 && ((int)t)%((int)dt_total_energy) == 0 )
    {
        total_energy << t << " " << get_total_energy(X,V,N,D,sqrt(dX2Remove),sqrt(dX2Remove),force_scaling,AddNumber,RemoveProb,rho) << endl;
    }
	if (dt_op_simple != 0 && ((int)t)%((int)dt_op_simple) == 0 )
	{
		//cout << "Starting OP simple\n";
	op_simple << t << " " << get_op_simple(X,dX,D,N,dX2Add,neigbours) << endl;
	}
	if (dt_op != 0 && ((int)t)%((int)dt_op) == 0 )
	{
		//cout << "starting OP \n";
	op << t << " " << get_op(X,dX,D,N,dX2Add,neigbours) << endl;
    div_op++;
	}
	if (dt_op_local!= 0 && ((int)t)%((int)dt_op_local) == 0 )
	{
    	get_op_local(X,dX,D,N,dX2Add,t,neigbours,rho,min(L[0],L[1])/3.0);
        div_op_local++;
	}
	if (dt_mean_length != 0 && ((int)t)%((int)dt_mean_length) == 0)
	{
		//cout << "starting mean length\n";
		get_mean_length(X,dX,D,connections,mean);
		length << t << " " << mean[0] << " " << mean[1] << endl;
	}
	if (dt_mean_num != 0 && ((int)t)%((int)dt_mean_num) == 0)
	{
//		cout << "starting mean conn\n";
	conns << t << " " << get_num_conn(connections,N) << endl;
	}
	if (dt_pkf != 0 && ((int)t)%((int)dt_pkf) == 0 )
	{
	//	cout << "starting pkf\n";
		get_pkf(X,dX,D,N,t,dX2Add,rho,min(L[0],L[1])/2.0);
	}
	if (dt_mean_vel != 0 && ((int)t)%((int)dt_mean_vel) == 0 )
	{
		mean_vel << t << " " <<  get_mean_vel(V,N,D) << endl;
    }
	if (dt_msd_distr != 0 && ((int)t)%((int)dt_msd_distr) == 0 )
	{
        get_msd_distr(msd_val,t,N,sqrt(dX2Remove));
    }
	if (dt_mean_length_distr != 0 && ((int)t)%((int)dt_mean_length_distr) == 0 )
	{
        get_mean_length_distr(X,connections,t,dX,sqrt(dX2Remove),D);
    }
	if (dt_life_time!= 0 && ((int)t)%((int)dt_life_time) == 0 )
	{
        get_life_time(life_time_distr,t,101,life_step,norm);
    }
	if (dt_energy_distr!= 0 && ((int)t)%((int)dt_energy_distr) == 0 )
	{
        get_energy_distr(energy_distr,life_time_distr,t,101,life_step,norm);
    }
	if (dt_plot_delaunay != 0 && ((int)t)%((int)dt_plot_delaunay) == 0 )
	{
        plot_delaunay(connections,t,D,X,L,"connections",cell_radius,marked,force_type);
		//cout << "starting msd\n";
	}
	if (dt_plot_work_hist!= 0 && ((int)t)%((int)dt_plot_work_hist) == 0 )
	{
        plot_work_hist(work_histogramm,dw,work_max,t);
    }
	if (dt_sigma_over_mu!= 0 && ((int)t)%((int)dt_sigma_over_mu) == 0 )
    {
        sigma_out << get_sigma_over_mu(neigbours,X,D,N) << endl;
    }
}
    count++;
    t = (floats)count*dt;
    //cout << "time: " << t << " pos: " << X[0] << " " << X[1] << endl;
} 


  
  // free memory
  for (wi=0; wi <= WallN; wi++)
    delete[] Wall_AtomIndex[wi];
  delete[] Wall_AtomIndex;
  
  // close writing streams
  BoundaryPosOutStream.close();
  BoundaryForceOutStream.close();
  WallPosOutStream.close();
  WallForceOutStream.close();
}
