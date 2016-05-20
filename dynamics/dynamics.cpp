#include "../emacss.h"

/*---------------------------------------------------------------------------*/
// Purely dynamical differential equations (see AG2012, s2 & GALB2013, s2)

double dynamics::dNdt(){
  double dNdt = 0;
//FilippoC
  if (mynode->galaxy.f == 0)
    dNdt -= falpha()*xi()*mynode->N/mynode->t_rhp;		//Equation (6) AG2012
  
//FContenta
  return dNdt;
}

double dynamics::dmmdt(){			//Equation (??) AGLB2014
  double dmmdt = 0;
  if (mynode->galaxy.f == 0)
    dmmdt = gamma()*mynode->mm/mynode->t_rhp;
  return dmmdt;
}

double dynamics::drdt(){			//Equation (7) AG2012
  double drdt = 0;
  if (mynode->galaxy.f == 0)
    drdt += mu()*mynode->rh/mynode->t_rhp;
  return drdt;
}

double dynamics::dkdt(){			//Equation (5) GALB2013
  double dkdt = 0;
  dkdt += lambda()*mynode->kappa/mynode->t_rhp;
  return dkdt;
}

double dynamics::dtrhdt(){			//Counts t_rh
  double dtrhdt = 0;
  dtrhdt += 1.0/mynode->t_rh;
  return dtrhdt;			
}

double dynamics::dtrhpdt(){			//Counts t_rh'
  double dtrhdt = 0;
  dtrhdt += 1.0/mynode->t_rhp;
  return dtrhdt;			
}

double dynamics::drcdt(){			//Equation (10) GALB2013
  double drcdt = 0;

  drcdt += delta()*mynode->rc/mynode->t_rhp;
  return drcdt;                          
}

double dynamics::falpha(){
  double falpha = 0;

  falpha =  2.0*alpha()/3.0 - alpha()*alpha()/12.0;

  return falpha;
}

/*---------------------------------------------------------------------------*/
//Dynamical Dimensionless Parameter Equations (see AG2012, s2 & GALB2013, s2,3)

double dynamics::xi(){                          //Equation (26) AG2012
  double xi = 0;				//Equation (22) GALB2013
						//Equation (??) AGLB2013
  xi = xi_i()+xi_e();	
  
  return xi;
}
    
double dynamics::xi_e(){			//Equation (??) AGLB2013
  double xi = 0;	  	

  xi += F()*xi1*(1.0-P())+(f+(1-f)*F())*(3.0/5.0)*mynode->E.zeta*P();
  
  return xi;
}

double dynamics::xi_i(){
  double xi = 0;
  
  if (mynode->galaxy.type > 0 && mynode->s == 1)
      xi += f_ind()*myse->gamma_se();
  
  return xi;
}

double dynamics::gamma(){
  double gamma = 0;
  
  gamma = gamma_dyn()-myse->gamma_se();
  
  return gamma;  
}

double dynamics::gamma_dyn(){                   //Equation (??) AGL2013
  double gamma = 0, Fmin = 0, X = 0;
  
   if (mynode->s == 1 && mynode->galaxy.type > 0){
 
    Fmin = m_ej()/mynode->mm;
  
    X = (1.0-Fmin)*pow((mynode->MS-3.0)/(myse->MS_1-3.0),q);
    X *= fabs(mynode->m_max-mynode->mm)/mynode->m_max;

    gamma += X*xi(); 
  }
 return gamma;
}

double dynamics::lambda(){			//Equation (12) GALB2013
  double lambda = 0;
  
  if (mynode->s == 0){
    lambda += K()*(delta() - mu()); 
    if (mynode->E.source == 1) lambda += (k() - mynode->kappa)/k();
  }
  else if (mynode->trhpelapsed > 0.5*nc){
    lambda += (mynode->k1-mynode->kappa)*((2.0*mynode->trhpelapsed/nc)-1.0);
  }
   return lambda;
}

double dynamics::mu(){				//Equation (19) GALB2013e
  double mu = 0.0, rhrj = 0;
  
  if (mynode->s ==0 ){
    if (mynode->galaxy.type > 0) rhrj = mynode->Rhj;
  
    if (mynode->E.source == 0){
      mu = K()*delta();                       //Increase due to shrinking core
      mu += (rhrj/mynode->kappa - 2.0)*xi();//Shrink due to removal of stars
      mu /= (1.0 + K());// Correction because rh appears in expression for kappa
     }

     if (mynode->E.source == 1){
       mu = mynode->E.zeta - 2.0*xi();         // Make the balance (AG2012)
       mu += 2.0/3.0*K()*xi()/(1.0+mynode->N/N3);// Extra terms
     }
  }
  else{
      mu = myse->epsilon() - 2.0*xi() + lambda() + 2.0*gamma();  
  }
  return mu;
 }

 double dynamics::delta(){			//Equation (17) GALB2013
   double delta = 0;
   
   if (mynode->E.source == 0)
     delta += delta_1 + delta_2*mynode->t_rhp/mynode->t_rc;
   if (mynode->E.source == 1)
     delta += mu() + 2.0/3.0*xi()/(1.0 + mynode->N/N3);

   return delta;
 }

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
//Extra required parameter equations (see AG2012, s2 & AGLB2013, s2,3)

double dynamics::P(){				 //Equation (25) AG2012
  double f = 0, g = 0;
  
  f += pow((mynode->N*log(mynode->gamma*N1))/  \
		 (N1*log(mynode->gamma*mynode->N)),(1.0-x));
  if (mynode->galaxy.type > 0) g += pow((mynode->rv/mynode->rj)/R1,z);
  
  return f*g;
}

double dynamics::F(){				 //Equation (25) AG2012
  double F = 1;
  
  if (mynode->s == 0) F = Rchmin()/mynode->Rch; 
  else if (mynode->E.source < 2){
    F = 0;
//    F = pow(mynode->trhpelapsed/nc,1.5);
    if (mynode->trhpelapsed > 0.5*nc) 
      F = (2.0*mynode->trhpelapsed)/nc-1.0;
  }
  if (mynode->tcc == 0 && F == 1) mynode->tcc = mynode->time;
  
  return F;
}

double dynamics::f_ind(){
  double f= 0;
  
  if (mynode->Rhj > R1) 
     f =Y*pow(mynode->Rhj-R1,b);

  return f;
}

double dynamics::m_ej(){
  double m = 0;

  m = Fej*(mynode->mm-mynode->m_min)+mynode->m_min;
  
  return m;
}

double dynamics::K(){				 //Equation (29) GALB2013          
  double K = 0;
  
  K = 2.0*mynode->Rch*(mynode->k0-mynode->k1)*\
	exp(-pow(mynode->Rch/mynode->Rch0,2.0))/(sqrt(M_PI)*mynode->Rch0*k());
  
  return K;
}

double dynamics::k(){				 //Breakdown, GALB2013           
  double k = 0;
  
  k = (mynode->k1+(mynode->k0-mynode->k1)*erf(mynode->Rch/mynode->Rch0));
  
  return k;
}

void dynamics::tdf(double M_star, double R_star, double pcMyr, double T_star){                                  //Equation (4) AG2012

t_df = 225*pow((mynode->galaxy.R*R_star)/1000,2)*
	(mynode->galaxy.v/pcMyr)/(T_star/R_star)/
	((mynode->mm*mynode->N*M_star)/1e5);
t_df = t_df/T_star;
}

double dynamics::Rchmin(){          
  // Definition of the minimum ratio r_c/r_h that defines core collapse
  double Rchmin = 0;

  if (mynode->E.source == 0) Rchmin += pow(N2/mynode->N + N2/N3,2.0/3.0);
  else Rchmin += mynode->Rch;
  
  return Rchmin;
}

void dynamics::reset_K_constants(){          
  //Changes parameters to post-collapse values.
    
  if (mynode->s == 0){ 
    mynode->Rch0 = 0.2;
    mynode->k0 = 0.19;
    mynode->k1 = 0.27;
    mynode->E.source = 1;                          //Changes energy source
  }
}

double dynamics::fv(double vs) {
  int gamma = mynode->galaxy.gamma;
  double mg = mynode->galaxy.M;
  double a = mynode->galaxy.scale;
  double r = mynode->galaxy.R;

  double pot = -(mg/a) * ( 1.0/ (2.0-gamma));
  double fac = (r/(r+a));
  pot *= 1.0 - (gamma == 1 ? fac : fac*fac); 

  
  //Calculate f(vs) at r.
  double E = -(0.5*vs*vs + pot) * (a/mg);
  double E2 = E*E;
				   
    
  //Calculate f(E) (Dehnen 1993)
  double fE = 0.0;

  if (mynode->galaxy.gamma == 0) {
    fE = sqrt(2.0*E)*((3.0-4.0*E)/(1.0-2.0*E));
    fE -= 3.0*asinh(sqrt(2.0*E/(1.0 - 2.0*E)));
    fE *= 3.0*mg/(2.0*pow(M_PI,3) * pow(mg*a, 1.5));
  }
  else if (gamma == 1.0) {
    double E3 = E2*E;
    double E05 = sqrt(E);
    
    fE = 3.0*asin(E05) - sqrt(E*(1.0-E)) * (3.0 + 2.0*E - 24.0*E2 + 16.0*E3);

    const double EXPVAL = 11.313708498984761; // 2^3.5
    fE /= EXPVAL * pow(M_PI,3) * pow(1.0-E,2.5);
  }
  else{
    //NEED PRESCRIPTION
    cerr << "No prescription for f(E) for Gamma = " << gamma << endl;
    exit(1);
  }
  
  fE = fE*4.0*M_PI*vs*vs;
  return fE;
}

double dynamics::integrate(double a, double b) {
  const double h=1.0e-1;

  double inf=a;
  double sup=a+h;
  double tot=0.0;

  // Simpson's Rule
  while (inf < b) {
    double loc_h = sup-inf;

    tot += loc_h / 6.0 * (fv(inf) + 4.0 * fv((inf+sup) * 0.5) + fv(sup));
    

    inf += h;
    sup += h;
    sup = min(sup, b);
  }

  return tot;
}

double dynamics::dV_dR(){                        //computed by Mark and James: (dV/dR + V/R)
  if (mynode->galaxy.gamma == 0) 
    return sqrt(mynode->galaxy.M)*(mynode->galaxy.scale + mynode->galaxy.R*0.5)*pow(mynode->galaxy.scale + mynode->galaxy.R,-2.5);
  else
    return sqrt(mynode->galaxy.M)*(3*mynode->galaxy.scale + mynode->galaxy.R)/(2.0*pow(mynode->galaxy.R,0.5)*(mynode->galaxy.scale + mynode->galaxy.R));
}


double dynamics::dr2dt() { // Dynamical Friction (James)
  if (mynode->galaxy.f == 1 && mynode->galaxy.type == 3) {
    int gamma = mynode->galaxy.gamma;
    double mg = mynode->galaxy.M;
    double a = mynode->galaxy.scale;
    double r = mynode->galaxy.R;
    double Ms = mynode->N*mynode->mm;

    double ratio = (r/(r+a));
    if (gamma == 0.0)
      ratio *= ratio;
    
    double vesc = 2.0*(mg/a) * (1.0/(2.0-gamma)) * (1.0 - ratio);
    vesc = sqrt(vesc);

    double rhor = ((3.0-gamma)*mg*a) / (4.0*M_PI * (gamma == 0 ? 1.0 : r) * pow(r+a, 4.0-gamma));

    //assume circular orbit and calculate bmax and v
    double bmax = min(r*(r+a)/(a*gamma + 4.0*r),r);
    double v = vcirc();

    //calculate f(v*<v)
    double fnorm = integrate(1e-12,vesc-1.e-12);
    double fvlow = integrate(1e-12,v-1.e-12);
    double rt = mynode->rj;
    
    double Mrt = (4.0*M_PI*rhor/3.0)*pow(rt, 3);
    
    double r_in = r-rt;
    double m_in = mg*pow(r_in/(r_in+a),(3.0-gamma)); //menc(r_in)    

    double Lambda;
    if (Mrt >= m_in)
      Lambda = bmax/max(rt,(max((Ms/(v*v)), mynode->rh)));
    else
      Lambda = bmax/(max((Ms/(v*v)), mynode->rh));

    double Cfric = 2.0*M_PI*Ms*rhor*log(Lambda*Lambda+1.0);
    double fr = -Cfric*fvlow/fnorm/(v*v);

    //cout << "T_star = " << mynode->T_star << endl;
    //cout << "R_star = " << mynode->R_star << endl;
    //cout << "fr = " << fr / mynode->T_star / mynode->T_star * mynode->R_star << endl;
    cerr << fvlow << " " << fnorm << endl;
    
    return fr / dV_dR();
  }
  
  return 0.0;
}

double dynamics::vcirc(){                         //circular velocity Dehnen's models
    double vcirc = 0.0;
    vcirc = pow(mynode->galaxy.M*pow(mynode->galaxy.R , 2-mynode->galaxy.gamma) /
        pow(mynode->galaxy.R + mynode->galaxy.scale, 3-mynode->galaxy.gamma),0.5);

  return vcirc;
}


double dynamics::alpha(){                      //-dln(rho)/dlnr
  double alpha = 0;
  if (mynode->galaxy.type == 3) {
    alpha = mynode->galaxy.gamma - mynode->galaxy.R*(mynode->galaxy.gamma - 4.0)/(mynode->galaxy.R + mynode->galaxy.scale);
    //alpha = min(alpha,3.0);
  }
  else
    alpha = 2;


  //  cout << "alpha(" << mynode->galaxy.R << ") = " << alpha << endl;

  return alpha;
}
