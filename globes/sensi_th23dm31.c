/* GLoBES -- General LOng Baseline Experiment Simulator */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"

#include <globes/globes.h>   /* GLoBES library */
#define RULE_T2K_NUMU         0
#define RULE_T2K_NUE          1
#define RULE_T2K_NUMU_NOSC    1
Bool_t isFitdata = false;
//#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
//char MYFILE[]="t2kUDChi_data_th23dm31";
//char MYFILE[]="t2kUDChi_sensi_th23dm31";
char MYFILE[]="t2kUDChi_sensi_th23dm31_trueth23e0d6";

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
    //for sensitivity
    double *true_rates_NUMU = glbGetRuleRatePtr(exp,RULE_T2K_NUMU);
    double *true_rates_NUE = glbGetRuleRatePtr(exp,RULE_T2K_NUE);
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
      if(isFitdata) chi2    += poisson_likelihood(data_mu[i], fit_rate);
      else chi2    += poisson_likelihood(true_rates_NUMU[i], fit_rate);
    }

    glbGetEnergyWindowBins(exp, RULE_T2K_NUE, &ew_low, &ew_high);
    for (int i=ew_low; i <= ew_high; i++)
    {
      fit_rate = (1 + x[0] + x[6]) * (signal_e[i] + bg_e[i]
                 + x[4] * (bin_centers[i] - E_center) * (signal_e[i] + bg_e[i]) );
        if(isFitdata) chi2    += poisson_likelihood(data_e[i], fit_rate);
        else chi2    += poisson_likelihood(true_rates_NUE[i], fit_rate);
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

int main(int argc, char *argv[])
{ 
    int nbins2RUN=-1;
    for(int i=1; i<argc; i++){
        if((strcmp(argv[i],"-n")==0) && i<argc){
            i++;
            nbins2RUN=atoi(argv[i]);
            continue;
        }
    }
    
    /* Initialize libglobes */
    glbInit(argv[0]);
    
    /* Initialize experiment NFstandard.glb */
    /* Initialize experiment T2K */
    
    glbDefineChiFunction(&chiT2K_FDonly,    7, "chiT2K",       NULL);
    glbDefineChiFunction(&chiT2K_FDonly,    0, "chiT2K-nosys", NULL);
    glbInitExperiment(AEDLFILE,&glb_experiment_list[0],&glb_num_of_exps);
    
    /* Intitialize output */
    TFile *poutfile = new TFile(Form("output/%s_bin%d.root",MYFILE,nbins2RUN),"RECREATE");
    TTree *SensiTree = new TTree("SensiTree","tree");
    double Chisq;
    double ChisqSys;
    double dm31test;
    double th23test;
    int binIndexx;
    int binIndexy;
    SensiTree->Branch("Chisq",&Chisq ,"Chisq/D");
    SensiTree->Branch("ChisqSys",&ChisqSys ,"ChisqSys/D");
    SensiTree->Branch("dm31test",&dm31test ,"dm31test/D");
    SensiTree->Branch("th23test",&th23test ,"th23test/D");
    SensiTree->Branch("binIndexx",&binIndexx ,"binIndexx/I");
    SensiTree->Branch("binIndexy",&binIndexy ,"binIndexy/I");
    
    
    
    /* Define standard oscillation parameters */
    //http://www.nu-fit.org/?q=node/211 fit wo atmospheric
    double theta12 = asin(sqrt(0.310));
    double theta13 = asin(sqrt(0.02241));
    //double theta23 = M_PI/4;
    //TMath::ASin(sqrt(0.6))*180/TMath::Pi()
    double theta23 = M_PI*50.768480/(4*45.0);
    //double theta23 = M_PI*43.5/(4*45.0);
    double deltacp = 0;
    double sdm = 7.39e-5;
    double ldm = 2.523e-3;
    
    double dm31_test ;
    double th23_test ;
    
    double signdm13 = -1. ;
    
    double chi0n_sys,chipn_sys,chi0i_sys,chipi_sys, chimin_sys ;
    double chi0n_cor,chipn_cor,chi0i_cor,chipi_cor, chimin_cor ;
    double chimh_sys, chimh_cor/*, chimh_sys_min, chimh_cor_min*/ ;
    
    double thetheta13,x,y, res1,res2;
    
    double chi_sys_nuapp ;
    double chi_sys_nudis ;
    double chi_sys_anuapp ;
    double chi_sys_anudis ;
    //double chi_sys_nuapp_t2k ;
    double chi_sys_nova ;
    
    
    /* Initialize parameter and projection vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params minimum =     glbAllocParams();
    
    
    //no need
    glb_projection myprojection = glbAllocProjection();
    //free theta_13 and dcp
    glbDefineProjection(myprojection,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FREE,GLB_FIXED,GLB_FIXED) ;
       glbSetDensityProjectionFlag(myprojection,GLB_FIXED,GLB_ALL) ;
       glbSetProjection(myprojection) ;
    
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbSetCentralValues(true_values);
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    
    /* Set starting values and input errors for all projections */
    glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0.0,sdm*0.03,ldm*0.05);
    glbSetDensityParams(input_errors,0.05,GLB_ALL);
    glbSetInputErrors(input_errors);
    
    double sinsqth23min = 0.3;
    double sinsqth23max = 0.7;
    double th23min = asin(sqrt(sinsqth23min));
    double th23max = asin(sqrt(sinsqth23max));
    int NPOINTTH23 = 80;
    double th13step = (th23max-th23min)/(1.*NPOINTTH23);
    double sinsqth13step = (sinsqth23max-sinsqth23min)/(1.*NPOINTTH23);
    
    //double theta13 = asin(sqrt(0.02241));
    
    /*double dm31min = 2.4e-3;
    int NPOINT = 40;
    double dm31step = 0.005e-3;*/
    double dm31min = 2.0e-3;
    int NPOINT = 40;
    double dm31step = 0.025e-3;
    double dm31max = dm31min+(NPOINT)*dm31step;
    TH2D *hchi_sys = new TH2D("hchi_sys","",NPOINTTH23,sinsqth23min,sinsqth23max,NPOINT,dm31min,dm31max);
    TH2D *hchi_cor = new TH2D("hchi_cor","",NPOINTTH23,sinsqth23min,sinsqth23max,NPOINT,dm31min,dm31max);
    
    // Loop over true values of delta
    // for (y = -1. ; y < 1.05 ; y=y+0.05){
    int county = nbins2RUN;
    int countx = 0;
    
    /*for (y = dm31min ; y < dm31max ; y=y+dm31step){*/
    dm31_test = dm31min+(nbins2RUN+1/2.)*dm31step ;
    county =nbins2RUN+1;
    
    countx = 0;//reset here
    for (x = sinsqth23min ; x < sinsqth23max ; x=x+sinsqth13step) {
        
        //glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,-ldm);
        //glbSetOscParams(test_values,-1.*ldm,GLB_DM_31);
        th23_test = asin(sqrt(x+sinsqth13step/2.)) ;
        countx += 1;
        glbSetOscParams(test_values,dm31_test,GLB_DM_31);
        glbSetOscParams(test_values,th23_test,GLB_THETA_23);
        glbSetCentralValues(test_values);
        glbDefineParams(input_errors,theta12*0.03,theta13*0.03,theta23*0.11,0.0,sdm*0.03,ldm*0.05);
        glbSetDensityParams(input_errors,0.05,GLB_ALL);
        glbSetInputErrors(input_errors);
        /* Compute Chi^2 for two-parameter correlation: minimize over deltacp only */
        glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON) ;
        
        //chi_sys_nuapp=glbChiSys(test_values,1,0);//for rule #1
        //chi_sys_nudis=glbChiSys(test_values,1,1);// for rule #2
        //chi_sys_anuapp=glbChiSys(test_values,2,0);// for rule #3
        //chi_sys_anudis=glbChiSys(test_values,2,1);//for rule #4
        // chi_sys_nuapp_t2k=glbChiSys(test_values,1,0);
        //chi_sys_nova = chi_sys_nuapp + chi_sys_nudis+ chi_sys_anuapp + chi_sys_anudis ;
        chimh_sys =   glbChiSys(test_values,0,GLB_ALL);//+glbChiSys(test_values,3,GLB_ALL);
        hchi_sys->SetBinContent(countx, county, chimh_sys);//fill histogram
        //printf("%g %g %g %g %g %g %g \n", y,x, chi_sys_nuapp, chi_sys_anuapp, chi_sys_nova, chi_sys_nuapp_t2k,chimh_sys) ;
        
        
        /* Compute Chi^2 for full correlation: minimize over all but delta */
        chimh_cor=glbChiNP(test_values,minimum,GLB_ALL);
        //         glbPrintParams(stdout,minimum) ;
        //         printf("%g %g %g %g %g %g %g \n", y,x, th23_test, chimh_sys, sqrt(chimh_sys), chimh_cor, sqrt(chimh_cor)) ;
        
        //         printf("%g %g %g %g  \n",y,x,chimh_sys, chimh_cor);
        hchi_cor->SetBinContent(countx, county, chimh_cor);//fill histogram
        //if(chimh_sys<chimh_sys_min){chimh_sys_min = chimh_sys ; }
        //if(chimh_cor<chimh_cor_min){chimh_cor_min = chimh_cor ; }
        printf("Processing truepoint %d/%d testpoint %d/%d \n",county,NPOINT,countx,NPOINTTH23);
        poutfile->cd();
        Chisq = chimh_cor;
        ChisqSys = chimh_sys;
        dm31test = dm31_test;
        th23test = x+sinsqth13step/2.;//sin^2
        binIndexy = county;
        binIndexx = countx;
        SensiTree->Fill();
        
    }
    
    
    /* }*/
    //printf("chi IH syst minimum %g \n",chimh_sys_min);
    //printf("chi IH cor minimum %g \n",chimh_cor_min);
    
    poutfile->cd();
    hchi_sys->Write();
    hchi_cor->Write();
    SensiTree->Write();
    
    
    poutfile->Close();
    /* Destroy parameter and projection vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors);
    glbFreeProjection(myprojection);
    
    exit(0);
    
 
}
