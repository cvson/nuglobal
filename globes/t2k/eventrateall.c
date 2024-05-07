/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * Example: Correlation between s22th13 and deltacp
 * Compile with ``make example1''
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"

#include <globes/globes.h>   /* GLoBES library */
//#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
//char MYFILE[]="eventrate_all_45d0_hyperktdr_wsmear_newNC_res.root";
char MYFILE[]="t2k2016_eventrate.root";
//char MYFILE[]="eventrate_app_49d5.root";
//char MYFILE[]="eventrate_app_43d5.root";

int main(int argc, char *argv[])
{ 
    /* Initialize libglobes */
    glbInit(argv[0]);
    
    /* Initialize experiment NFstandard.glb */
    glbInitExperiment("t2k-fd_chiOFF.glb",&glb_experiment_list[0],&glb_num_of_exps);
    
    /* Intitialize output */
    //InitOutput(MYFILE,"Format: Log(10,s22th13)   deltacp   chi^2 \n");
    
    TFile *poutfile = new TFile(MYFILE,"RECREATE");
    /*Retrieve the flux info.*/
    const int NFLAVOR = 3;
    const int NPOLARITY = 2;//neutrino and anti-neutrino
    const double expDistance = 295;//km
    const int NBINFLUX = 100;
    const double EMAXFLUX = 5.0;//GeV
    
    int polVal[NPOLARITY]={1,-1};//1 for neutrino and -1 for antineutrino
    
    int NFLUX = glbGetNumberOfFluxes(0);
    printf("Number of Flux %d\n",NFLUX);
    TH1D *hflux[NFLUX][NPOLARITY][NFLAVOR];
    //this flux is applied oscillation probability?
    for (int iflux=0; iflux<NFLUX; ++iflux) {
        for (int ipol=0; ipol<NPOLARITY; ++ipol) {
            for (int iflav=0; iflav<NFLAVOR; ++iflav) {
                hflux[iflux][ipol][iflav] = new TH1D(Form("flux%d_pol%d_flav%d",iflux,ipol,iflav+1),"",NBINFLUX,0,EMAXFLUX);//10 GeV
                for (Int_t ibin=1; ibin<=NBINFLUX; ++ibin) {
                    double energyval = hflux[iflux][ipol][iflav]->GetBinCenter(ibin);
                    double fluxval = glbFlux(0,iflux,energyval,expDistance,iflav+1,polVal[ipol]);
                    hflux[iflux][ipol][iflav]->SetBinContent(ibin,fluxval);
                }
            }
        }
    }
    /*Retrieve the cross section info*/
    int NXSEC = 3;//this is manual input, check your globes file
    TH1D *hxsec[NXSEC][NPOLARITY][NFLAVOR];
    //this flux is applied oscillation probability?
    for (int ixsec=0; ixsec<NXSEC; ++ixsec) {
        for (int ipol=0; ipol<NPOLARITY; ++ipol) {
            for (int iflav=0; iflav<NFLAVOR; ++iflav) {
                hxsec[ixsec][ipol][iflav] = new TH1D(Form("XSEC%d_pol%d_flav%d",ixsec,ipol,iflav+1),"",NBINFLUX,0,EMAXFLUX);//10 GeV
                for (Int_t ibin=1; ibin<=NBINFLUX; ++ibin) {
                    double energyval = hxsec[ixsec][ipol][iflav]->GetBinCenter(ibin);
                    double xsecval = glbXSection(0,ixsec,energyval,iflav+1,polVal[ipol]);
                    hxsec[ixsec][ipol][iflav]->SetBinContent(ibin,xsecval);
                }
            }
        }
    }
    
    
    /* Define standard oscillation parameters */
    //T2HK TDR
    double theta12 = asin(sqrt(0.8704))/2.0;
    double theta13 = asin(sqrt(0.1))/2.0;
    double theta23 = M_PI/4;
    //double theta23 = M_PI*49.5/(4*45.0);
    //double theta23 = M_PI*43.5/(4*45.0);
    //double deltacp = 0;
    //double deltacp = M_PI;
    double deltacp = M_PI/2;
    double sdm = 7.6e-5;
    double ldm = 2.4e-3;
    
    /* Initialize parameter vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(true_values,1.0,GLB_ALL);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbSetDensityParams(test_values,1.0,GLB_ALL);
    
    /* The simulated data are computed */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    
    //get event rate
    int i;
    int n_bins = glbGetNumberOfBins(0);
    double emin, emax;
    int ew_low, ew_high;
    glbGetEminEmax(0,&emin, &emax);
    glbGetEnergyWindowBins(0, 0, &ew_low, &ew_high);
    
    ///////////////////////////////////////////////////////////////////////////////////
    //Get information in channel-based
    ///////////////////////////////////////////////////////////////////////////////////
    int NCHANNELALL = glbGetNumberOfChannels(0);//
    printf("Number of Channel %d\n",NCHANNELALL);
    //channel rate pre-smearing
    TH1D **hrate_channel_pre = new TH1D*[NCHANNELALL];
    //channel rate post-smearing
    TH1D **hrate_channel_post = new TH1D*[NCHANNELALL];
    
    for (Int_t ichannel=0; ichannel<NCHANNELALL; ++ichannel) {
        hrate_channel_pre[ichannel] = new TH1D(Form("hrate_pre_channel%d",ichannel),"",n_bins,emin,emax);
        hrate_channel_post[ichannel] = new TH1D(Form("hrate_post_channel%d",ichannel),"",n_bins,emin,emax);
        double *channelrate_pre = glbGetChannelRatePtr(0,ichannel,GLB_PRE);
        double *channelrate_post =glbGetChannelRatePtr(0,ichannel,GLB_POST);
        for (i=0; i<n_bins; i++) {
            hrate_channel_pre[ichannel]->SetBinContent(i+1,channelrate_pre[i]);
            hrate_channel_post[ichannel]->SetBinContent(i+1,channelrate_post[i]);
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    //get information in ruled-based
    ///////////////////////////////////////////////////////////////////////////////////
    int NRULES = glbGetNumberOfRules(0);
    printf("Number of Rules %d\n",NRULES);
    
    TH1D **hrate_all_rule = new TH1D*[NRULES];//
    TH1D **hrate_sig_rule = new TH1D*[NRULES];//
    TH1D **hrate_bkg_rule = new TH1D*[NRULES];//
    
    
    for (Int_t irule=0; irule<NRULES; ++irule) {
        hrate_all_rule[irule] = new TH1D(Form("hrate_all_%d",irule),"",n_bins,emin,emax);
        hrate_sig_rule[irule] = new TH1D(Form("hrate_sig_%d",irule),"",n_bins,emin,emax);
        hrate_bkg_rule[irule] = new TH1D(Form("hrate_bkg_%d",irule),"",n_bins,emin,emax);
    }
    //fill the histogram
    
    
    double *aEbin = new double[n_bins];
    double *aRateTotal = new double[n_bins];
    printf("Emin % g - Emax % g",emin, emax);
    printf("\n");
    for (Int_t irule=0; irule<NRULES; ++irule) {
        printf("Simulated rates, experiment 0, rule %d: \n",irule);
        
        double *true_rates_rule0 = glbGetRuleRatePtr(0,irule);
        double *true_sigrates_rule0 = glbGetSignalRatePtr(0,irule);
        double *true_bkgrates_rule0 = glbGetBGRatePtr(0,irule);
        for (i=0; i<n_bins; i++) {
            printf("% g",true_rates_rule0[i]);
            printf("\n");
            aEbin[i] = emin + (emax-emin)*(i+0.5)/(n_bins*1.0);
            aRateTotal[i] = true_rates_rule0[i];
            hrate_all_rule[irule]->SetBinContent(i+1,true_rates_rule0[i]);
            hrate_sig_rule[irule]->SetBinContent(i+1,true_sigrates_rule0[i]);
            hrate_bkg_rule[irule]->SetBinContent(i+1,true_bkgrates_rule0[i]);
        }//end ibin
        
        //to get channel rate
        //int nchannelSig = glbGetChannelInRule(0,irule,-1,GLB_SIG);
        //int nchannelBkg = glbGetChannelInRule(0,irule,-1,GLB_BG);
        
        
        //delete true_rates_rule0;
        //delete true_sigrates_rule0;
        //delete true_bkgrates_rule0;
        
    }
    
    
    TGraph *pgrRateTotal = new TGraph(n_bins,aEbin,aRateTotal);
    poutfile->cd();
    pgrRateTotal->Write("pgrRateTotal");
    for (Int_t ichannel=0; ichannel<NCHANNELALL; ++ichannel) {
        hrate_channel_pre[ichannel]->Write(Form("hrate_pre_channel%d",ichannel));
        hrate_channel_post[ichannel]->Write(Form("hrate_post_channel%d",ichannel));
    }
    for (Int_t irule=0; irule<NRULES; ++irule) {
        hrate_all_rule[irule]->Write(Form("hrate_all_rule%d",irule));
        hrate_sig_rule[irule]->Write(Form("hrate_sig_rule%d",irule));
        hrate_bkg_rule[irule]->Write(Form("hrate_bkg_rule%d",irule));
    }
    //fill channels in rules
    for (Int_t irule=0; irule<NRULES; ++irule) {
        //for the signal
        int nchannelSig = glbGetLengthOfRule(0,irule,GLB_SIG);
        for (int ichansig=0; ichansig<nchannelSig; ++ichansig) {
            int channelinexinrule =glbGetChannelInRule(0,irule,ichansig,GLB_SIG);
            double channelcoeff = glbGetCoefficientInRule(0,irule,ichansig,GLB_SIG);
            printf("rule %d, signal chan%d, coff %g: \n",irule,channelinexinrule,channelcoeff);
            
            TH1D* hruthchannel= (TH1D*)hrate_channel_post[channelinexinrule]->Clone(Form("rule%d_sig_post_channel%d",irule,ichansig));
            hruthchannel->Scale(channelcoeff);
            hruthchannel->Write(Form("hrate_sig_rule%d_chan%d_post",irule,channelinexinrule));
            delete hruthchannel;
            
        }
        //for the background
        int nchannelBkg = glbGetLengthOfRule(0,irule,GLB_BG);
        for (int ichanbkg=0; ichanbkg<nchannelBkg; ++ichanbkg) {
            int channelinexinrule =glbGetChannelInRule(0,irule,ichanbkg,GLB_BG);
            double channelcoeff = glbGetCoefficientInRule(0,irule,ichanbkg,GLB_BG);
            printf("rule %d, b background chan%d, coff %g: \n",irule,channelinexinrule,channelcoeff);
            TH1D* hruthchannel= (TH1D*)hrate_channel_post[channelinexinrule]->Clone(Form("rule%d_bkg_post_channel%d",irule,ichanbkg));
            hruthchannel->Scale(channelcoeff);
            hruthchannel->Write(Form("hrate_bkg_rule%d_chan%d_post",irule,channelinexinrule));
            delete hruthchannel;
        }
        
        
        printf("For Rule %d, Number of channel in Signal %d, in Background %d \n",irule,nchannelSig, nchannelBkg);
    }
    
    
    //Save flux info
    for (int iflux=0; iflux<NFLUX; ++iflux) {
        for (int ipol=0; ipol<NPOLARITY; ++ipol) {
            for (int iflav=0; iflav<NFLAVOR; ++iflav) {
                hflux[iflux][ipol][iflav]->Write(Form("flux%d_pol%d_flav%d",iflux,ipol,iflav+1));
                
            }
        }
    }
    
    //save Xsec info.
    for (int ixsec=0; ixsec<NXSEC; ++ixsec) {
        for (int ipol=0; ipol<NPOLARITY; ++ipol) {
            for (int iflav=0; iflav<NFLAVOR; ++iflav) {
                hxsec[ixsec][ipol][iflav]->Write(Form("xsec%d_pol%d_flav%d",ixsec,ipol,iflav+1));
                
            }
        }
    }
    
    poutfile->Close();
    delete pgrRateTotal;
    //delete hrate_all_0;
    delete poutfile;
    /* Iteration over all values to be computed */
    /*double thetheta13,x,y,res;
     
     for(x=-4.0;x<-2.0+0.01;x=x+2.0/50)
     for(y=0.0;y<200.0+0.01;y=y+200.0/50)
     {
     // Set vector of test values
     thetheta13=asin(sqrt(pow(10,x)))/2;
     glbSetOscParams(test_values,thetheta13,GLB_THETA_13);
     glbSetOscParams(test_values,y*M_PI/180.0,GLB_DELTA_CP);
     
     // Compute Chi^2 for all loaded experiments and all rules
     res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
     AddToOutput(x,y,res);
     }*/
    
    /* Destroy parameter vector(s) */
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    
    exit(0);
}
