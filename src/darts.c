#include <R.h>
#include <Rmath.h>

/******************************************************************************************/

void EM(int *x, int *np, double *sinit, double *s, double *ll, int *niterp, 
	double *R);
double EMStep(int *x, int n, double s, double *R);
double ComputeExp(int x, double *a, double *b);
void ComputeExpConstants(double s, double *R, double *a);
double Integ1(double s, double r1, double r2);
void ComputeProbConstants(double s, double *R, double *b);
double Integ2(double s, double r1, double r2);
double Loglik(int *x, int n, double s, double *R);
double ComputeProb(int x, double *b);

/******************************************************************************************/

void EMCov(int *x, int *np, double *Sinit, double *S1, double *S2, double *S3, 
	   double *ll, int *niterp, int *computellp, double *R, double *ar, int *ii);
void EMCovStep(int *x, int n, double* S, double *R, double *ar, int *ii, double *A);
void SimulateExp(int x, double *S, double *R, double *ar, int *ii, double *B);
void RandomPt(int x, double *R, double *ar, int *ii, double *z);
void RandomSlicePt(int x, double r1, double r2, int *ii, double *z);
void RandomCirclePt(double r1, double r2, double *z);
double RandomR(double r1, double r2);
double LoglikCov(int *x, int n, double *S, double *R, double *ar);

/******************************************************************************************/

void BuildScoreMatrix(int *A, double *R, int *S);
int Score(double x, double y, double *R, int *S);
double MyMod(double x, double y);

/******************************************************************************************/

void EM(int *x, int *np, double *sinit, double *s, double *ll, int *niterp, 
	double *R) {
  int n = np[0];
  int niter = niterp[0];
  double scur = sinit[0];
  int i;

  for (i=0; i<niter; i++) {
    s[i] = EMStep(x,n,scur,R);
    ll[i] = Loglik(x,n,s[i],R);
    scur = s[i];
  }
}  

double EMStep(int *x, int n, double s, double *R) {
  double a[6];
  double b[6];
  ComputeExpConstants(s,R,a);
  ComputeProbConstants(s,R,b);
  double e = 0;
  int i;

  for (i=0; i<n; i++) {
    e += ComputeExp(x[i],a,b);
  }

  return e/(2*n);
}

double ComputeExp(int x, double *a, double *b) {
  if (x==1 || x==5 || x==7 || x==11 || x==13 || x==17 || x==19) {
    return a[2]/b[2];
  }
  else if (x==2 || x==4 || x==8 || x==10 || x==14 || x==16 || x==20) {
    return (a[2]+a[3])/(b[2]+b[3]);
  }
  else if (x==3 || x==9 || x==15) {
    return (a[2]+a[4])/(b[2]+b[4]);
  }
  else if (x==6 || x==12 || x==18) {
    return (a[2]+a[3]+a[4])/(b[2]+b[3]+b[4]);
  }
  else if (x==24 || x==30 || x==36) {
    return (a[3]+a[4])/(b[3]+b[4]);
  }
  else if (x==22 || x==26 || x==28 || x==32 || x==34 || x==38 || x==40) {
    return a[3]/b[3];
  }
  else if (x==21 || x==27 || x==33 || x==39 || x==42 || x==45 || x==48 ||
           x==51 || x==54 || x==57 || x==60) {
    return a[4]/b[4];
  }
  else if (x==25) {
    return a[1]/b[1];
  }
  else if (x==50) {
    return a[0]/b[0];
  }
  else {
    return a[5]/b[5];
  }
}

void ComputeExpConstants(double s, double *R, double *a) {
  a[0] = Integ1(s,0,R[0]); 
  a[1] = Integ1(s,R[0],R[1]); 
  a[2] = Integ1(s,R[1],R[2])/20 + Integ1(s,R[3],R[4])/20;
  a[3] = Integ1(s,R[4],R[5])/20;
  a[4] = Integ1(s,R[2],R[3])/20;
  a[5] = Integ1(s,R[5],-1);
}

double Integ1(double s, double r1, double r2) {
  if (r2==-1) {
    // r2 is assumed to be infinity
    return (r1*r1+2*s)*exp(-r1*r1/(2*s));
  }
  else {
    return (r1*r1+2*s)*exp(-r1*r1/(2*s)) - (r2*r2+2*s)*exp(-r2*r2/(2*s));
  }
}

void ComputeProbConstants(double s, double *R, double *b) {
  b[0] = Integ2(s,0,R[0]); 
  b[1] = Integ2(s,R[0],R[1]); 
  b[2] = Integ2(s,R[1],R[2])/20 + Integ2(s,R[3],R[4])/20;
  b[3] = Integ2(s,R[4],R[5])/20;
  b[4] = Integ2(s,R[2],R[3])/20;
  b[5] = Integ2(s,R[5],-1);
}

double Integ2(double s, double r1, double r2) {
  if (r2==-1) {
    // r2 is assumed to be infinity
    return exp(-r1*r1/(2*s));
  }
  else {
    return exp(-r1*r1/(2*s)) - exp(-r2*r2/(2*s));
  }
}

double Loglik(int *x, int n, double s, double *R) {
  double b[6];
  ComputeProbConstants(s,R,b);
  double p = 0;
  int i;

  for (i=0; i<n; i++) {
    p += log(ComputeProb(x[i],b));
  }
  
  return p;
}

double ComputeProb(int x, double *b) {
  if (x==1 || x==5 || x==7 || x==11 || x==13 || x==17 || x==19) {
    return b[2];
  }
  else if (x==2 || x==4 || x==8 || x==10 || x==14 || x==16 || x==20) {
    return b[2]+b[3];
  }
  else if (x==3 || x==9 || x==15) {
    return b[2]+b[4];
  }
  else if (x==6 || x==12 || x==18) {
    return b[2]+b[3]+b[4];
  }
  else if (x==24 || x==30 || x==36) {
    return b[3]+b[4];
  }
  else if (x==22 || x==26 || x==28 || x==32 || x==34 || x==38 || x==40) {
    return b[3];
  }
  else if (x==21 || x==27 || x==33 || x==39 || x==42 || x==45 || x==48 ||
           x==51 || x==54 || x==57 || x==60) {
    return b[4];
  }
  else if (x==25) {
    return b[1];
  }
  else if (x==50) {
    return b[0];
  }
  else {
    return b[5];
  }
}

/******************************************************************************************/

void EMCov(int *x, int *np, double *Sinit, double *S1, double *S2, double *S3, 
	   double *ll, int *niterp, int *computellp, double *R, double *ar, int *ii) {
  int n = np[0];
  int niter = niterp[0];
  int computell = computellp[0];
  double Scur[3];
  Scur[0] = Sinit[0]; 
  Scur[1] = Sinit[1]; 
  Scur[2] = Sinit[2];
  double A[3];
  int i;
  
  for (i=0; i<niter; i++) {
    EMCovStep(x,n,Scur,R,ar,ii,A);
    S1[i] = A[0]; 
    S2[i] = A[1]; 
    S3[i] = A[2];
    if (computell==0) ll[i] = 0;
    else ll[i] = LoglikCov(x,n,A,R,ar);
    Scur[0] = S1[i];
    Scur[1] = S2[i];
    Scur[2] = S3[i];
  }
}

void EMCovStep(int *x, int n, double *S, double *R, double *ar, int *ii, double *A) {
  double B[3];
  double C[3];
  C[0] = 0; 
  C[1] = 0;
  C[2] = 0;  
  int i;

  for (i=0; i<n; i++) {
    SimulateExp(x[i],S,R,ar,ii,B);
    C[0] += B[0];
    C[1] += B[1];
    C[2] += B[2];
  }

  A[0] = C[0]/n;
  A[1] = C[1]/n;
  A[2] = C[2]/n;
}

void SimulateExp(int x, double *S, double *R, double *ar, int *ii, double *B) {
  double det = S[0]*S[1]-S[2]*S[2];
  double z[2];
  double w, W = 0;
  B[0] = 0;
  B[1] = 0;
  B[2] = 0;
  int i;

  for (i=0; i<1000; i++) {
    RandomPt(x,R,ar,ii,z);
    w = exp(-(S[1]*z[0]*z[0] - 2*S[2]*z[0]*z[1] + S[0]*z[1]*z[1])/(2*det));
    B[0] += z[0]*z[0]*w;
    B[1] += z[1]*z[1]*w;
    B[2] += z[0]*z[1]*w;
    W += w;
  }

  B[0] /= W;
  B[1] /= W;
  B[2] /= W;
}

void RandomPt(int x, double *R, double *ar, int *ii, double *z) {
  double u = unif_rand();

  if (x==1 || x==5 || x==7 || x==11 || x==13 || x==17 || x==19) {
    if (u <= ar[2]/(ar[2]+ar[3])) {
      RandomSlicePt(x,R[1],R[2],ii,z);
    }
    else {
      RandomSlicePt(x,R[3],R[4],ii,z);
    }
  }

  else if (x==2 || x==4 || x==8 || x==10 || x==14 || x==16 || x==20) {
    if (u <= ar[4]/(ar[4]+ar[2]+ar[3])) {
      RandomSlicePt(x/2,R[4],R[5],ii,z);
    }
    else if (u <= (ar[4]+ar[2])/(ar[4]+ar[2]+ar[3])) {
      RandomSlicePt(x,R[1],R[2],ii,z);
    }
    else {
      RandomSlicePt(x,R[3],R[4],ii,z);
    }
  }

  else if (x==3 || x==9 || x==15) {
    if (u <= ar[5]/(ar[5]+ar[2]+ar[3])) {
      RandomSlicePt(x/3,R[2],R[3],ii,z);
    }
    else if (u <= (ar[5]+ar[2])/(ar[5]+ar[2]+ar[3])) {
      RandomSlicePt(x,R[1],R[2],ii,z);
    }
    else {
      RandomSlicePt(x,R[3],R[4],ii,z);
    }
  }

  else if (x==6 || x==12 || x==18) {
    if (u <= ar[5]/(ar[5]+ar[4]+ar[2]+ar[3])) {
      RandomSlicePt(x/3,R[2],R[3],ii,z);
    }
    else if (u <= (ar[5]+ar[4])/
             (ar[5]+ar[4]+ar[2]+ar[3])) {
      RandomSlicePt(x/2,R[4],R[5],ii,z); 
    }
    else if (u <= (ar[5]+ar[4]+ar[2])/
             (ar[5]+ar[4]+ar[2]+ar[3])) {
      RandomSlicePt(x,R[1],R[2],ii,z);
    }
    else {
      RandomSlicePt(x,R[3],R[4],ii,z); 
     } 
  }

  else if (x==24 || x==30 || x==36) {
    if (u <= ar[5]/(ar[5]+ar[4])) {
      RandomSlicePt(x/3,R[2],R[3],ii,z);
    }
    else {
      RandomSlicePt(x/2,R[4],R[5],ii,z);
    }
  }

  else if (x==22 || x==26 || x==28 || x==32 || x==34 || x==38 || x==40) {
    RandomSlicePt(x/2,R[4],R[5],ii,z);
   }

  else if (x==21 || x==27 || x==33 || x==39 || x==42 || x==45 || x==48 ||
           x==51 || x==54 || x==57 || x==60) {
    RandomSlicePt(x/3,R[2],R[3],ii,z);
  }

  else if (x==25) {
    RandomCirclePt(R[0],R[1],z);
  }

  else if (x==50) {
    RandomCirclePt(0,R[0],z);
  }
}

void RandomSlicePt(int x, double r1, double r2, int *ii, double *z) {
  int k = ii[x-1];
  double th = -2*PI/40 + (k-1)*2*PI/20 + 2*PI/20*unif_rand();
  th = PI/2 - th;
  double r = RandomR(r1,r2);
  z[0] = r*cos(th);
  z[1] = r*sin(th);
}

void RandomCirclePt(double r1, double r2, double *z) {
  double th = 2*PI*unif_rand();
  double r = RandomR(r1, r2);
  z[0] = r*cos(th);
  z[1] = r*sin(th);
}

double RandomR(double r1, double r2) {
  return sqrt(r1*r1 + (r2*r2-r1*r1)*unif_rand());
}

double LoglikCov(int *x, int n, double *S, double *R, double *ar) {
  return 0;
}

/******************************************************************************************/

void BuildScoreMatrix(int *A, double *R, int *S) {
  int i,j;

  for (i=0; i<681; i++) {
    for (j=0; j<681; j++) {
      A[i+681*j] = Score((double)(i-340),(double)(j-340),R,S);
    }
  }
}

int Score(double x, double y, double *R, int *S) {
  // Compute the radius
  double r = sqrt(x*x+y*y);

  // Check if it's off the board (do this for speed)
  if (r > R[5]) return 0;
  
  // Check for a center bullseye
  if (r <= R[0]) return 50;

  // Check for a donut bullseye
  if (r <= R[1]) return 25;

  // Now get the angle
  double theta = atan2(y, x);
  double phi = MyMod(PI/2 - theta + 2*PI/40, 2*PI);

  // Now get the number
  int i = (int)floor(phi/(2*PI)*20) + 1;
  if (i > 20) i = 20;
  int n = S[i-1];  

  // Check for a single
  if (r <= R[2]) return n;

  // Check for a triple
  if (r <= R[3]) return 3*n;

  // Check for a single
  if (r <= R[4]) return n;

  // If we got here, it must be a double
  return 2*n;
}

double MyMod(double x, double y) {
  return x-y*floor(x/y);
}
