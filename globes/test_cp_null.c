/* GLoBES -- General LOng Baseline Experiment Simulator */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "TFile.h"
#include "TH1.h"

#include <globes/globes.h>   /* GLoBES library */
#define RULE_T2K_NUMU         0
#define RULE_T2K_NUE          1
#define RULE_T2K_NUMU_NOSC    1
//#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="t2kUDChi_sensi_cp_null_all.root";

//char MYFILE1[]="output/test_cp_null_1.dat";
//char MYFILE2[]="output/test_cp_null_2.dat";
//char MYFILE3[]="output/test_cp_null_3.dat";
char AEDLFILE[]="t2k/t2k-fd.glb";
static inline double square(double x)
{
  return x*x;
}

static inline double gauss_likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}
/* Poisson likelihood */
static inline double poisson_likelihood(double true_rate, double fit_rate)
{
  double res;
  res = fit_rate - true_rate;
  if (true_rate > 0)
  {
    if (fit_rate <= 0.0)
      res = 1e100;
    else
      res += true_rate * log(true_rate/fit_rate);
  }
  else
    res = fabs(res);

  return 2.0 * res;
}

// ***************************************************************************
// * this use near detector
// * Calculate chi^2 for T2K, including the following systematical errors:   *
// *   x[ 0]: Correlated beam flux normalization error                       *
// *   x[ 1]: Energy scale error for mu-like events                          *
// *   x[ 2]: Energy scale error for e-like events                           *
// *   x[ 3]: Spectral tilt for mu-like events                               *
// *   x[ 4]: Spectral tilt for e-like events                                *
// *   x[ 5]: Normalization mu-like events                                   *
// *   x[ 6]: Normalization e-like events                                    *
// * The energy scale errors are included even in SYS_OFF simulations to     *
// * help resolving degeneracies between systematics params and osc. params  *
// ***************************************************************************/
double chiT2K(int exp, int rule, int n_params, double *x, double *errors,
              void *user_data)
{
  double exp_F = exp;
  double exp_N = exp - 1;
  if (exp_N < 0)
    return NAN;

  // T2K far detector data from Run 1-4
  // from Lorena Escudero's thesis, http://www.t2k.org/docs/thesis/070, fig 5.15 (left)
  static const double data_mu_F[] = {
    0, 0, 0, 0, 0,  3, 3, 8, 6, 4,
    4, 2, 3, 1, 4,  1, 1, 4, 4, 2,
    2, 1, 4, 5, 2,  2, 3, 1, 0, 5,
    3, 0, 3, 0, 1,  0, 1, 1, 0, 0,
    1, 1, 0, 2, 1,  0, 1, 2, 0, 1,
    0, 0, 0, 1, 0,  3, 0, 1, 1, 2,
    3, 1, 3, 2,
    3, 4,  1, 1,
    1
  };
  static const double data_e_F[] = {
    0, 0, 1, 0, 0,  0, 0, 2, 2, 2,
    3, 3, 3, 0, 1,  4, 2, 2, 1, 1,
    0, 0, 1, 0, 0
  };

  int n_bins_F = glbGetNumberOfBins(exp_F);
  int n_bins_N = glbGetNumberOfBins(exp_N);
  double emin_F, emax_F, emin_N, emax_N;
  int ew_low, ew_high;
  double fit_rate;
  double chi2 = 0.0;

  glbGetEminEmax(exp_F, &emin_F, &emax_F);
  glbGetEminEmax(exp_N, &emin_N, &emax_N);

  // Systematics ON
  // --------------
  if (n_params > 0)
  {
    double *bin_centers_F = glbGetBinCentersListPtr(exp_F);
    double E_center_F = 0.5 * (emax_F + emin_F);
    double signal_mu_F[n_bins_F], signal_e_F[n_bins_F];
    double bg_mu_F[n_bins_F], bg_e_F[n_bins_F];
    double signal_mu_N[n_bins_N], bg_mu_N[n_bins_N];
    double signal_mu_nosc_N[n_bins_N], bg_mu_nosc_N[n_bins_N];
    glbShiftEnergyScale(x[1], glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUMU),
                        signal_mu_F, n_bins_F, emin_F, emax_F);
    glbShiftEnergyScale(x[1], glbGetBGFitRatePtr(exp_F, RULE_T2K_NUMU),
                        bg_mu_F, n_bins_F, emin_F, emax_F);
    glbShiftEnergyScale(x[2], glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUE),
                        signal_e_F, n_bins_F, emin_F, emax_F);
    glbShiftEnergyScale(x[2], glbGetBGFitRatePtr(exp_F, RULE_T2K_NUE),
                        bg_e_F, n_bins_F, emin_F, emax_F);

    glbShiftEnergyScale(0.0, glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU),
                        signal_mu_N, n_bins_N, emin_N, emax_N);
    glbShiftEnergyScale(0.0, glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU),
                        bg_mu_N, n_bins_N, emin_N, emax_N);
    glbShiftEnergyScale(0.0, glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC),
                        signal_mu_nosc_N, n_bins_N, emin_N, emax_N);
    glbShiftEnergyScale(0.0, glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC),
                        bg_mu_nosc_N, n_bins_N, emin_N, emax_N);

    // Alter far detector prediction according to near detector measurement
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      signal_mu_F[i] *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
      signal_e_F[i]  *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
    }

    // Compare FD data to ND-corrected prediction
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[5]) * ( signal_mu_F[i] + bg_mu_F[i]
                 + x[3] * (bin_centers_F[i] - E_center_F) * (signal_mu_F[i] + bg_mu_F[i]) );
      chi2    += poisson_likelihood(data_mu_F[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[6]) * (signal_e_F[i] + bg_e_F[i]
                 + x[4] * (bin_centers_F[i] - E_center_F) * (signal_e_F[i] + bg_e_F[i]) );
      chi2    += poisson_likelihood(data_e_F[i], fit_rate);
    }

    // Systematics priors
    for (int i=0; i < n_params; i++)
      chi2 += square(x[i] / errors[i]);
  }

  // Systematics OFF
  // ---------------
  else
  {
    double *signal_mu_F      = glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUMU);
    double *bg_mu_F          = glbGetBGFitRatePtr(exp_F, RULE_T2K_NUMU);
    double *signal_e_F       = glbGetSignalFitRatePtr(exp_F, RULE_T2K_NUE);
    double *bg_e_F           = glbGetBGFitRatePtr(exp_F, RULE_T2K_NUE);
    double *signal_mu_N      = glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU);
    double *bg_mu_N          = glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU);
    double *signal_mu_nosc_N = glbGetSignalFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC);
    double *bg_mu_nosc_N     = glbGetBGFitRatePtr(exp_N, RULE_T2K_NUMU_NOSC);

    // Alter far detector prediction according to near detector measurement
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      signal_mu_F[i] *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
      signal_e_F[i]  *= (signal_mu_nosc_N[i]+bg_mu_nosc_N[i]) / (signal_mu_N[i]+bg_mu_N[i]);
    }

    // Compare FD data to ND-corrected prediction
    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = signal_mu_F[i] + bg_mu_F[i];
      chi2    += poisson_likelihood(data_mu_F[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp_F, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = signal_e_F[i] + bg_e_F[i];
      chi2    += poisson_likelihood(data_e_F[i], fit_rate);
    }
  }

  return chi2;
}

// ***************************************************************************
// * Calculate chi^2 for T2K, including the following systematical errors:   *
// *   x[ 0]: Correlated beam flux normalization error                       *
// *   x[ 1]: Energy scale error for mu-like events                          *
// *   x[ 2]: Energy scale error for e-like events                           *
// *   x[ 3]: Spectral tilt for mu-like events                               *
// *   x[ 4]: Spectral tilt for e-like events                                *
// *   x[ 5]: Normalization mu-like events                                   *
// *   x[ 6]: Normalization e-like events                                    *
// * The energy scale errors are included even in SYS_OFF simulations to     *
// * help resolving degeneracies between systematics params and osc. params  *
// ***************************************************************************/
double chiT2K_FDonly(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
{
  // T2K far detector data from Run 1-4
  // from Lorena Escudero's thesis, http://www.t2k.org/docs/thesis/070, fig 5.15 (left)
  static const double data_mu[] = {
    0, 0, 0, 0, 0,  3, 3, 8, 6, 4,
    4, 2, 3, 1, 4,  1, 1, 4, 4, 2,
    2, 1, 4, 5, 2,  2, 3, 1, 0, 5,
    3, 0, 3, 0, 1,  0, 1, 1, 0, 0,
    1, 1, 0, 2, 1,  0, 1, 2, 0, 1,
    0, 0, 0, 1, 0,  3, 0, 1, 1, 2,
    3, 1, 3, 2,
    3, 4,  1, 1,
    1
  };
  static const double data_e[] = {
    0, 0, 1, 0, 0,  0, 0, 2, 2, 2,
    3, 3, 3, 0, 1,  4, 2, 2, 1, 1,
    0, 0, 1, 0, 0
  };

  int n_bins = glbGetNumberOfBins(exp);
  double emin, emax;
  int ew_low, ew_high;
  double fit_rate;
  double chi2 = 0.0;

  glbGetEminEmax(exp, &emin, &emax);

  // Systematics ON
  // --------------
  if (n_params > 0)
  {
    double *bin_centers = glbGetBinCentersListPtr(exp);
    double E_center = 0.5 * (emax + emin);
    double signal_mu[n_bins], signal_e[n_bins];
    double bg_mu[n_bins], bg_e[n_bins];
    glbShiftEnergyScale(x[1], glbGetSignalFitRatePtr(exp, RULE_T2K_NUMU),
                        signal_mu, n_bins, emin, emax);
    glbShiftEnergyScale(x[1], glbGetBGFitRatePtr(exp, RULE_T2K_NUMU),
                        bg_mu, n_bins, emin, emax);
    glbShiftEnergyScale(x[2], glbGetSignalFitRatePtr(exp, RULE_T2K_NUE),
                        signal_e, n_bins, emin, emax);
    glbShiftEnergyScale(x[2], glbGetBGFitRatePtr(exp, RULE_T2K_NUE),
                        bg_e, n_bins, emin, emax);

    // Compare FD data to prediction
    glbGetEnergyWindowBins(exp, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[5]) * ( signal_mu[i] + bg_mu[i]
                 + x[3] * (bin_centers[i] - E_center) * (signal_mu[i] + bg_mu[i]) );
      chi2    += poisson_likelihood(data_mu[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[6]) * (signal_e[i] + bg_e[i]
                 + x[4] * (bin_centers[i] - E_center) * (signal_e[i] + bg_e[i]) );
      chi2    += poisson_likelihood(data_e[i], fit_rate);
    }

    // Systematics priors
    for (int i=0; i < n_params; i++)
      chi2 += square(x[i] / errors[i]);
  }

  // Systematics OFF
  // ---------------
  else
  {
    double *signal_mu = glbGetSignalFitRatePtr(exp, RULE_T2K_NUMU);
    double *bg_mu     = glbGetBGFitRatePtr(exp, RULE_T2K_NUMU);
    double *signal_e  = glbGetSignalFitRatePtr(exp, RULE_T2K_NUE);
    double *bg_e      = glbGetBGFitRatePtr(exp, RULE_T2K_NUE);

    // Compare FD data to ND-corrected prediction
    glbGetEnergyWindowBins(exp, RULE_T2K_NUMU, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = signal_mu[i] + bg_mu[i];
      chi2    += poisson_likelihood(data_mu[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = signal_e[i] + bg_e[i];
      chi2    += poisson_likelihood(data_e[i], fit_rate);
    }
  }

  return chi2;
}


/*void init_t2k(double *errors){
    glbDefineChiFunction(&chiT2K_FDonly,    7, "chiT2K",       NULL);
    glbDefineChiFunction(&chiT2K_FDonly,    0, "chiT2K-nosys", NULL);
    glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps);
    glbSetChiFunction(glb_num_of_exps-1, GLB_ALL, GLB_ON, "chiZero", NULL);
    //int i;
    //for (i=0; i < 29; i++)
     // { re_sys_errors[i] = errors[i];}
}*/


int main(int argc, char *argv[])
{ 
    /* Initialize libglobes */
    glbInit(argv[0]);
    
    /* Initialize experiment T2K */
    
    glbDefineChiFunction(&chiT2K_FDonly,    7, "chiT2K",       NULL);
    glbDefineChiFunction(&chiT2K_FDonly,    0, "chiT2K-nosys", NULL);
    glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps);
    //glbSetChiFunction(glb_num_of_exps-1, GLB_ALL, GLB_ON, "chiZero", NULL);
   
    
    
    /* Intitialize output */
    TFile *poutfile = new TFile(MYFILE,"RECREATE");
    
    
    
    /* Define standard oscillation parameters */
    //http://www.nu-fit.org/?q=node/211 fit wo atmospheric
    double theta12 = asin(sqrt(0.310));
    double theta13 = asin(sqrt(0.02241));
    double theta23 = M_PI/4;
    //double theta23 = M_PI*49.5/(4*45.0);
    //double theta23 = M_PI*43.5/(4*45.0);
    double deltacp = 0;
    double sdm = 7.39e-5;
    double ldm = 2.523e-3;
    
    double delta_true ;
    double delta_test ;
    
    //to define the mass hierarchy
    double signdm13 = -1. ;
    
    int count, i ;
    
    //locate 50 points of delta
    const int NDELTApoint = 50;//between [-M_PI, M_PI]
    double deltastep = 2./NDELTApoint;// in unit of M_PI
    float	delta[NDELTApoint] ;
    
    //chisquare array
    double chi0n_sys[NDELTApoint],chipn_sys[NDELTApoint],chi0i_sys[NDELTApoint],chipi_sys[NDELTApoint], chimin_sys ;
    double chi0n_cor[NDELTApoint],chipn_cor[NDELTApoint],chi0i_cor[NDELTApoint],chipi_cor[NDELTApoint], chimin_cor ;
    
    
    double thetheta13,x,y, res1,res2;
    
    double chimh_sys, chimh_cor, chimh_sys_min, chimh_cor_min;
    double chi_sys_min, chi_cor_min, chi_sys_min_mh, chi_cor_min_mh;
    
    chi_sys_min = 999. ;
    chi_cor_min = 999. ;
    chi_sys_min_mh = 999. ;
    chi_cor_min_mh = 999. ;
    
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params fit_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params minimum = 	glbAllocParams();
    glb_params central_values = glbAllocParams();
    
    //glb_projection theta13_projection = glbAllocProjection();
    
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetCentralValues(true_values);
    
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    
    /* Set starting values and input errors for all projections */
    //this to be understand
    glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0.0,sdm*0.03,ldm*0.05);
    //	glbDefineParams(input_errors,theta12*0.0003,theta13*0.0003,theta23*0.0011,0.0,sdm*0.0003,ldm*0.0003);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetInputErrors(input_errors);
    
    /* Set true value delta = 0. , NH */
    delta_true = 0.*M_PI ;
    glbDefineParams(true_values,theta12,theta13,theta23,delta_true,sdm,ldm);
    glbSetOscParams(true_values,+1.*ldm,GLB_DM_31);
    glbSetOscParams(true_values,delta_true,GLB_DELTA_CP);
    glbSetCentralValues(true_values);
    
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    TH1D *hchisq_truedcp0_nh_sys = new TH1D("hchisq_truedcp0_nh_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcp0_nh_proj= new TH1D("hchisq_truedcp0_nh_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcp0_ih_sys = new TH1D("hchisq_truedcp0_ih_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcp0_ih_proj= new TH1D("hchisq_truedcp0_ih_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    // Loop over test values of delta for both normal and inverted hierarchy
    count = 0;
    for (y = -1. ; y < 1.0 ; y=y+deltastep){
        delta_test = (y+0.5*deltastep)*M_PI ;
        // for normal hierarchy and delta
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,1.*ldm);
        glbSetOscParams(test_values,1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        //chisqr with systematic only
        chi0n_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chi0n_sys[count]<0.){chi0n_sys[count]=0. ;}
        hchisq_truedcp0_nh_sys->SetBinContent(count+1,chi0n_sys[count]);
        
        //to make project on dCP, need to put input error from disappearance?
        glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0,sdm*0.03,ldm*0.05);
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        
        //projection on ChiDelta
        chi0n_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chi0n_cor[count]<0.){chi0n_cor[count]=0. ;}
        hchisq_truedcp0_nh_proj->SetBinContent(count+1, chi0n_cor[count]);
        
        //for invered hierarchy; why need to invert delta as well?
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,-1.*ldm);
        glbSetOscParams(test_values,-1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,-1.*delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        chi0i_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chi0i_sys[count]<0.){chi0i_sys[count]=0. ;}
        hchisq_truedcp0_ih_sys->SetBinContent(count+1,chi0i_sys[count]);
        
        //projection on ChiDelta
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        chi0i_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chi0i_cor[count]<0.){chi0i_cor[count]=0. ;}
        hchisq_truedcp0_ih_proj->SetBinContent(count+1, chi0i_cor[count]);
        
        delta[count]=y;
        
        printf("Computing w/ dcp=0 true: %d/%d  points \n",count, NDELTApoint) ;
        count++;
        	
    }
    
    
    /* Set true value delta = pi , NH */
    delta_true = 1.*M_PI ;
    glbDefineParams(true_values,theta12,theta13,theta23,delta_true,sdm,1.*ldm);
    glbSetOscParams(true_values,+1.*ldm,GLB_DM_31);
    
    glbSetOscParams(true_values,delta_true,GLB_DELTA_CP);
    glbSetCentralValues(true_values);
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    // Loop over test values of delta for both normal and inverted hierarchy
    count = 0;
    TH1D *hchisq_truedcppi_nh_sys = new TH1D("hchisq_truedcppi_nh_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcppi_nh_proj= new TH1D("hchisq_truedcppi_nh_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcppi_ih_sys = new TH1D("hchisq_truedcppi_ih_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_truedcppi_ih_proj= new TH1D("hchisq_truedcppi_ih_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    
    for (y = -1. ; y < 1.0 ; y=y+deltastep){
        delta_test = (y+0.5*deltastep)*M_PI ;
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,ldm);
        glbSetOscParams(test_values,1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        chipn_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chipn_sys[count]<0.){chipn_sys[count]=0. ;}
        hchisq_truedcppi_nh_sys->SetBinContent(count+1,chipn_sys[count]);
        
        //projection on ChiDelta
        glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0,sdm*0.03,ldm*0.05);
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        
        chipn_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chipn_cor[count]<0.){chipn_cor[count]=0. ;}
        hchisq_truedcppi_nh_proj->SetBinContent(count+1,chipn_cor[count]);
        
        //chisqr with systematic only
        glbDefineParams(test_values,theta12,theta13,theta23,delta_test,sdm,-1.*ldm);
        glbSetOscParams(test_values,-1.*ldm,GLB_DM_31);
        glbSetOscParams(test_values,-1.*delta_test,GLB_DELTA_CP);
        glbSetOscParams(test_values,theta13,GLB_THETA_13);
        
        chipi_sys[count]=glbChiSys(test_values,GLB_ALL,GLB_ALL);
        if(chipi_sys[count]<0.){chipi_sys[count]=0. ;}
        hchisq_truedcppi_ih_sys->SetBinContent(count+1,chipi_sys[count]);
        
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        glbSetCentralValues(test_values);
        
        chipi_cor[count]=glbChiDelta(test_values,minimum,GLB_ALL);
        if(chipi_cor[count]<0.){chipi_cor[count]=0. ;}
        hchisq_truedcppi_ih_proj->SetBinContent(count+1,chipi_cor[count]);
        
        printf("Computing w/ dcp=pi true: %d/%d  points \n",count, NDELTApoint) ;
        count++;
    }
    
    //Test for minimum chi^2 for each value of delta
    TH1D *hchisq_min_glob_sys = new TH1D("hchisq_min_glob_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_min_nh_nh_sys= new TH1D("hchisq_min_nh_nh_sys","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_min_glob_proj = new TH1D("hchisq_min_glob_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    TH1D *hchisq_min_nh_proj= new TH1D("hchisq_min_nh_proj","",NDELTApoint,-1.0*M_PI,1.0*M_PI);
    
    for (i = 0; i < NDELTApoint ; i++){
        chi_sys_min = 999. ;
        chi_cor_min = 999. ;
        chi_sys_min_mh = 999. ;
        chi_cor_min_mh = 999. ;
        //global minimum
        if(chi0n_sys[i]<chi_sys_min){chi_sys_min=chi0n_sys[i] ;}
        if(chipn_sys[i]<chi_sys_min){chi_sys_min=chipn_sys[i] ;}
        if(chi0i_sys[i]<chi_sys_min){chi_sys_min=chi0i_sys[i] ;}
        if(chipi_sys[i]<chi_sys_min){chi_sys_min=chipi_sys[i] ;}
        hchisq_min_glob_sys->SetBinContent(i+1,sqrt(chi_sys_min));
        
        //if mh is known
        if(chi0n_sys[i]<chi_sys_min_mh){chi_sys_min_mh=chi0n_sys[i] ;}
        if(chipn_sys[i]<chi_sys_min_mh){chi_sys_min_mh=chipn_sys[i] ;}
        hchisq_min_nh_nh_sys->SetBinContent(i+1,sqrt(chi_sys_min_mh));
        
        //global mimum of projection
        if(chi0n_cor[i]<chi_cor_min){chi_cor_min=chi0n_cor[i] ;}
        if(chipn_cor[i]<chi_cor_min){chi_cor_min=chipn_cor[i] ;}
        if(chi0i_cor[i]<chi_cor_min){chi_cor_min=chi0i_cor[i] ;}
        if(chipi_cor[i]<chi_cor_min){chi_cor_min=chipi_cor[i] ;}
        hchisq_min_glob_proj->SetBinContent(i+1,sqrt(chi_cor_min));
        
        //if MH is known
        if(chi0n_cor[i]<chi_cor_min_mh){chi_cor_min_mh=chi0n_cor[i] ;}
        if(chipn_cor[i]<chi_cor_min_mh){chi_cor_min_mh=chipn_cor[i] ;}
        hchisq_min_nh_proj->SetBinContent(i+1,sqrt(chi_cor_min_mh));
        
        
    }
//output
    poutfile->cd();
    hchisq_truedcp0_nh_sys ->Write();
    hchisq_truedcp0_nh_proj ->Write();
    hchisq_truedcp0_ih_sys ->Write();
    hchisq_truedcp0_ih_proj ->Write();
    
    hchisq_truedcppi_nh_sys ->Write();
    hchisq_truedcppi_nh_proj ->Write();
    hchisq_truedcppi_ih_sys ->Write();
    hchisq_truedcppi_ih_proj ->Write();
    
    hchisq_min_glob_sys ->Write();
    hchisq_min_nh_nh_sys ->Write();
    hchisq_min_glob_proj ->Write();
    hchisq_min_nh_proj ->Write();
    
    
    poutfile->Close();
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors);
    
    exit(0);
}
