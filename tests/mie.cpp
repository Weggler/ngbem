#include <comp.hpp>
#include <python_comp.hpp>

using namespace ngcomp;

double KAPPA = 5.;
const Complex I = Complex(0, 1);

// auxiliary routine evaluate the mie series
// r         [m]   = radius of the sphere
// w         [m/s] = radian frequency
// kappa           = wave number: kappa^2 = w^2 mu epsilon => kappa= w/c 
// theta,phi [rad] = Scattered field angles
// lmax            = Maximum mode number for computing Bessel functions
Complex*
mie_serie_current (double ctheta, double stheta, 
		   double radius, double kappa, int lmax)
{
  double pi = 3.1415926535897932;  // pi

  int     l;
  double  ll;
  double  s, ds, x;
  double  rbessel, ibessel;
  double  plk[lmax + 1] = {0.}, dplk[lmax + 1] = {0.};
  double  j, dj, jlk[4];
  double  p, dp, pl1[3];
  Complex h, dh, hlk[4];
  Complex S1p, S2p;
  Complex S1n, S2n;
  double  sgn1;
  Complex sgnI;
  Complex A, D, next;
  Complex *field;
  field = (Complex*) malloc(2*sizeof(Complex));

  int ierr, nz, c1, c2, nmax;
  double nu;
  Complex zx;
  double besselr[lmax+1];
  double besseli[lmax+1];
  double   dummy[lmax+1];
  Complex zhankel[lmax+1];
  Complex zbesselj[lmax+1];
  double cwrkr[lmax+1],  cwrki[lmax+1];
  Complex zwrk[lmax+1];
  
  // quasi evaluation point (in spherical coordinates)
  x  = kappa*radius;

  // evalutation of Associated Legendre Polynomials
  if(fabs(ctheta) < 0.999999)
    {
      for (l = 0; l < lmax; l++) {
	double pml1 = -std::assoc_legendre(l + 1, 1, ctheta);
	double pml = -std::assoc_legendre(l, 1, ctheta);
	plk[l] = pml1;
	dplk[l] = ((2. + l) * pml - (1. + l) * ctheta * pml1) / ((1. + ctheta) * (1. - ctheta));
      }
      p  =  plk[0];
      dp = dplk[0];
    }

  // evaluation of Hankel functions of second type
  // transform cylindrical Bessel functions into spherical Bessels!
  s  = sqrt(pi*x/2.);
  ds = 0.5*sqrt(pi/(2.*x));

  // gsl library
  // n == 1
  rbessel = std::cyl_bessel_j(1.5, x);
  ibessel = std::cyl_neumann(1.5, x);
  jlk[0]  = rbessel;
  hlk[0]  = rbessel - I*ibessel;

  // n == 2
  rbessel = std::cyl_bessel_j(2.5, x);
  ibessel = std::cyl_neumann(2.5, x);
  jlk[1]  = rbessel;
  hlk[1]  = rbessel - I*ibessel;

  // n=1 spherical Bessel&Hankel Balanis (11-240) p.655 
  //     and derivatives Balanis (IV-18) p.936 
  j  = s *jlk[0];
  dj = ds*jlk[0]+s*(-jlk[1]+1.5/x*jlk[0]);
  h  = s *hlk[0];
  dh = ds*hlk[0]+s*(-hlk[1]+1.5/x*hlk[0]);
  
  S1p = 0.0;
  S2p = 0.0;
  S1n = 0.0;
  S2n = 0.0;
  D   = 0.0;
  sgnI = -I; // => (-I)^l
  sgn1 = -1.;   // => (-1)^l

  for(l=1;l<=lmax;l++)
    {
      if(fabs(ctheta) < 0.999999)
	{
	  //     D != (0.1) due to Balanis, p.936 (IV-20)
	  D  = dj*h-j*dh;
	  A  = sgnI* (2.*l+1.)/(l*(l+1.));
	  next = A*((stheta*dp)/dh + I*p/(stheta*h));
	  if(std::imag(next) > 0.0) 
	    S1p += I*std::imag(next);
	  else
	    S1n += I*std::imag(next);
	  if(std::real(next) > 0.0) 
	    S1p += std::real(next);
	  else
	    S1n += std::real(next);

	  next = A*(p/(stheta*dh) - (stheta*dp)/(I*h));
	  if(std::imag(next) > 0.0) 
	    S2p += I*std::imag(next);
	  else
	    S2n += I*std::imag(next);
	  if(std::real(next) > 0.0) 
	    S2p += std::real(next);
	  else
	    S2n += std::real(next);
	}
      else 
	{
	  A  = sgnI* (2.*l+1.)/(l*(l+1.));
	  next = A*sgn1*(l*(l+1.)/2.)*(1./dh + I/h);
	  if(std::imag(next) > 0.0) 
	    S1p += I*std::imag(next);
	  else
	    S1n += I*std::imag(next);
	  if(std::real(next) > 0.0) 
	    S1p += std::real(next);
	  else
	    S1n += std::real(next);

	  next = A*sgn1*(l*(l+1.)/2.)*(1./dh - 1./(I*h));
	  if(std::imag(next) > 0.0) 
	    S2p += I*std::imag(next);
	  else
	    S2n += I*std::imag(next);
	  if(std::real(next) > 0.0) 
	    S2p += std::real(next);
	  else
	    S2n += std::real(next);
	}
     
      // recurrence relation for next Associated Legendre Polynomial
      // and next Hankel function
      p  =  plk[l];
      dp = dplk[l];

      // gsl library
      hlk[0]  = hlk[1];
      jlk[0]  = jlk[1];
      // n = l+1
      rbessel = std::cyl_bessel_j(l + 2.5, x);
      ibessel = std::cyl_neumann(l + 2.5, x);
      jlk[1]  = rbessel;
      hlk[1]  = rbessel - I*ibessel;

      // n = l+1
      j       = s *jlk[0];
      dj      = ds*jlk[0]+s*(-jlk[1]+(l+1.5)/x*jlk[0]);
      h       = s *hlk[0];
      dh      = ds*hlk[0]+s*(-hlk[1]+(l+1.5)/x*hlk[0]);

      sgnI *= -I;
      sgn1 *= -1.;
    }

  field[0] = S1p + S1n;
  field[1] = S2p + S2n;

  return field;
}

void mie_current(double *akappa, double *X, Complex *field)
{
  Complex *f;
  double x,y,z;
  double r, theta, phi;
  double ctheta, stheta;
  double kappa;
  kappa = *akappa;

  f = (Complex*) malloc(2*sizeof(Complex));
  
  x     = *X;
  y     = *(X+1);
  z     = *(X+2);
  
  r      = sqrt(x*x+y*y+z*z);
  phi    = atan2(y,x);
  theta  = acos(z/r);
  ctheta = cos(theta);
  stheta = sqrt(1.-ctheta*ctheta);

  f = mie_serie_current(ctheta,stheta,r,kappa,20);

  // Explanation of the factor, in order to not further 
  // change the hp-solution, I shift all on the mie solution, i.e.
  //  I/(c*mu*r*kappa) * -I*(mu*w)            = 1/r
  //  ----------------   -------------------   ------------------- 
  //  == I/(eta*x),see | the BIE solves for  | factor below
  //  p.294 Harrington | j_hp ~ -I*w*mu j_mie |

  f[0] *= -1./r*cos(phi);
  f[1] *= -1./r*sin(phi);

  // Neumann trace of total electric field component \gamma_N E^t 
  // in cartesian coordinates
  *(field  )=( f[0]*ctheta*cos(phi) - f[1]*sin(phi));
  *(field+1)=( f[0]*ctheta*sin(phi) + f[1]*cos(phi));
  *(field+2)=(-f[0]*stheta);

  free(f);
}

/*
  CoefficinetFunction returning mie current
  Base class is in ngsolve/fem/coefficient.hpp 
*/

class MieCurrent : public CoefficientFunction
{
public:
  MieCurrent() : CoefficientFunction(/*dimension = */ 3) { ; }

  // evaluation in one mapped point:
  virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
  {
    // cout << "MyXY - single point evaluation called" << endl;
    FlatVector<double> pnt = mip.GetPoint();
    double X[3] = {pnt(0), pnt(1), pnt(2)};
    Complex f[3] = {0.};
    mie_current(&KAPPA, X, f);
    //return std::norm(f[0]) * std::norm(f[0]) + std::norm(f[1]) * std::norm(f[1]) + std::norm(f[2]) * std::norm(f[2]);
    return std::real(f[2]);
  }

  // evaluate in all points of the integration-rule
  // values are stores as N x dim matrix
  virtual void Evaluate (const BaseMappedIntegrationRule & mir, BareSliceMatrix<Complex> values) const override
  {
    // cout << "MyXY - integration rule evaluation called" << endl;

    // all point coordinates of integration points (N x 2 matrix)
    SliceMatrix<double> pnts = mir.GetPoints();
    // cout << "pnts = " << pnts << endl;

    for (int i = 0; i < mir.Size(); i++) {
      double X[3] = {pnts(i, 0), pnts(i, 1), pnts(i, 2)};
      Complex f[3] = {0.};
      mie_current(&KAPPA, X, f);
      for (int j = 0; j < 3; j++)
	values(i,j) = f[j];
    }

    // more compact:
    // values.Col(0).Range(mir.Size()) = pw_mult(pnts.Col(0), pnts.Col(1));
  }  
};

extern "C" void Mie(py::object & res) {
  cout << "called mie series current" << endl;

  auto ngs = py::module::import("ngsolve");    
  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    

  
  py::class_<MieCurrent, shared_ptr<MieCurrent>, CoefficientFunction>
    (m, "MieCurrent", "CoefficientFunction that computes Mie series current.")
    .def(py::init<>())
    ;
  
  res = m;    
}    

