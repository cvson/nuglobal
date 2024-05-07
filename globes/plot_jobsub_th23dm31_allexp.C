#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TChain.h>
#include <TMath.h>
#include <TROOT.h>
#include <TGaxis.h>
#include <TList.h>
#include <TObjArray.h>
#include <TLegend.h>

#include <stdlib.h>
#include <iostream>
#include <vector>

//Subfunctions
void PrintHisto(TCanvas *canvas, TString name);
std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int twotheta, int graphNum);
void PreparePlot(TH2D* plot);
void PrepareBestPoint(TGraph* gr);
void SetStyleVariables(TStyle *t2kStyle);
void plot_jobsub_th23dm31_allexp(){
    //**** Set Style for Plots ****
    TStyle *t2kstyle = new TStyle("NU","NU approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("NU");
    
    TChain *Sensi = new TChain("SensiTree");
    TString subname="t2kUDChi_sensi_th23dm31_trueth23e0d6";
    //TString subname="t2kUDChi_sensi_th23dm31";
    //Sensi->Add("output/t2kUDChi_data_th23dm31_bin*.root");
    Sensi->Add(Form("output/%s_bin*.root",subname.Data()));
    std::cout << "Sensi Tree: " << Sensi->GetEntries() << std::endl;
    
    double dm31test = 0;
    double th23test = 0;
    int binIndexx = 0;
    int binIndexy = 0;
    double Chisq = 0.;
    Sensi->SetBranchAddress("dm31test", &dm31test);
    Sensi->SetBranchAddress("th23test", &th23test);
    Sensi->SetBranchAddress("binIndexx", &binIndexx);
    Sensi->SetBranchAddress("binIndexy", &binIndexy);
    Sensi->SetBranchAddress("Chisq", &Chisq);
    //the bining here is different from app/FitterSetting
    //TH2D *SA_th23dm31 = new TH2D("SA_th23dm31", "AsimovA Sensitivity", 51, -0.01, 1.01, 101, -0.01e-2, 2.01e-2);
    //TH2D *SA_th23dm31 = new TH2D("SA_th23dm31", "AsimovA Sensitivity", 51,0.796 , 1.204, 21, 1.975e-3, 3.025e-3);
     double th23min = asin(sqrt(0.3));
    double th23max = asin(sqrt(0.7));
    int NPOINTTH23 = 80;//40
    double th23step = (th23max-th23min)/(1.*NPOINTTH23);

    /*double dm31min = 2.4e-3;
    int NPOINT = 40;
    double dm31step = 0.005e-3;*/
    double dm31min = 2.e-3;
    int NPOINT = 40;
    double dm31step = 0.025e-3;
    
    double dm31max = dm31min+(NPOINT)*dm31step;

    //double M_PI = TMath::Pi();
    
    TH2D *SA_th23dm31 = new TH2D("SA_th23dm31", "", NPOINTTH23,0.3,0.7,NPOINT,dm31min,dm31max);
    //find minimize
    double globalmin = 999999.;
    double th23Testmin = 0;
    double dm31Testmin = 0;
    for(int i=0; i<Sensi->GetEntries(); i++){
        Sensi->GetEntry(i);
	if(Chisq<0) Chisq=0;
        if(Chisq<globalmin) {globalmin = Chisq;
        th23Testmin = th23test;
        dm31Testmin = dm31test;
        }
    }
    cout <<"min chisquare "<<globalmin<<" at dm31 "<<dm31Testmin<<" th23 "<<th23Testmin<<endl;
    
    for(int i=0; i<Sensi->GetEntries(); i++){
        Sensi->GetEntry(i);
        //debug
        //cout<<"at dm23" <<OscParams[2]<<" sin23 "<<OscParams[3]<<" likelihood "<<Chisq<<endl;
        //SA_th23dm31->SetBinContent(SA_th23dm31->FindBin(OscParams[3], OscParams[2]), Chisq);
        //after subtracting the minimal values
        SA_th23dm31->SetBinContent(SA_th23dm31->FindBin(th23test, dm31test), (Chisq-globalmin));
        
    }
    TCanvas *c1 = new TCanvas();
    c1->cd();
    SA_th23dm31->Draw("colz");
    PrintHisto(c1, Form("%s_allexp",subname.Data()));
    
    TFile *outFile = new TFile(Form("plots/%s.root",subname.Data()), "RECREATE");
    outFile->cd();
    //to make 68, 90, 99% C.L. contour
    int nlevels = 2;
    //double levels[2] = {1.1395, 2.3025};
    //double levels[2] = {4.61,11.83};//90% and 3sigma
    double levels[2] = {2.30,4.61};//68%, 90%
    TString gr_numbers[10] = {"0","1","2","3","4","5","6","7","8","9"};
    TCanvas *c2 = new TCanvas();
    c2->cd();
    TH2D* SA_th23dm31_cont = (TH2D*)SA_th23dm31->Clone();
    SA_th23dm31_cont->SetContour(nlevels, levels);
    SA_th23dm31_cont->Draw("cont list");
    c2->Update();
    TObjArray *contours_th23Test_th23True = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    if(!contours_th23Test_th23True){ std::cerr << "quitting!" << std::endl; return; }
    std::vector<TGraph*> gVect_th23Test_th23True_90;
    std::vector<TGraph*> gVect_th23Test_th23True_68;
    std::vector<TGraph*> gVect_th23Test_th23True_conv_90;
    std::vector<TGraph*> gVect_th23Test_th23True_conv_68;
    gVect_th23Test_th23True_90 = GetContourGraphs(contours_th23Test_th23True, "white", 0, 1);
    gVect_th23Test_th23True_68 = GetContourGraphs(contours_th23Test_th23True, "white", 0, 0);
    gVect_th23Test_th23True_conv_90 = GetContourGraphs(contours_th23Test_th23True, "black", 0, 1);
    gVect_th23Test_th23True_conv_68 = GetContourGraphs(contours_th23Test_th23True, "black", 0, 0);
    PreparePlot(SA_th23dm31);
    SA_th23dm31->Draw("colz");
    SA_th23dm31->GetZaxis()->SetRangeUser(0,30); 
    SA_th23dm31->Write();
    for(int j=0; j<gVect_th23Test_th23True_90.size(); j++){
        gVect_th23Test_th23True_90[j]->Draw("Lsame");
        gVect_th23Test_th23True_90[j]->Write("SA_th23dm31_cont_90_"+gr_numbers[j]);
        gVect_th23Test_th23True_conv_90[j]->Write("SA_th23dm31_cont_conv_90_"+gr_numbers[j]);
    }
    for(int j=0; j<gVect_th23Test_th23True_68.size(); j++){
        gVect_th23Test_th23True_68[j]->Draw("Lsame");
        gVect_th23Test_th23True_68[j]->Write("SA_th23dm31_cont_68_"+gr_numbers[j]);
        gVect_th23Test_th23True_conv_68[j]->Write("SA_th23dm31_cont_conv_68_"+gr_numbers[j]);
    }
    //  best_all->Draw("Psame");
    PrintHisto(c2, Form("%s_contour",subname.Data()));
    
}

void PrintHisto(TCanvas *canvas, TString name){
    //canvas->Print("plots/"+name+".eps");
    canvas->Print("plots/"+name+".pdf");
    //canvas->Print("plots/"+name+".png");
    //canvas->Print("plots/"+name+".C");
}

//#############################################################################
std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int twotheta, int graphNum){
    std::vector<TGraph*> gVect;
    double x,y;
    int ncontours = contours->GetEntries();
    std::cout << ncontours << std::endl;
    TList *list = (TList*)contours->At(graphNum);
    std::cout << list->GetEntries() << std::endl;
    
    for(int j=0; j<list->GetEntries(); j++){
        TGraph *gc = new TGraph();
        
        TGraph *gc_temp = (TGraph*)list->At(j);
        int count=0;
        std::cout << "points" << gc_temp->GetN() << std::endl;
        for(int i=0; i<gc_temp->GetN(); i++){
            gc_temp->GetPoint(i, x, y);
            if(twotheta==0){
                gc->SetPoint(i, x, y);
            }
            else{ //switch to sin^2(2 t23)
                //      if(x > 0.5){
                gc->SetPoint(count, 1-TMath::Power(1-2*x, 2), y);
                count++;
                //    }
            }
        }
        gc->SetLineWidth(2);
        gc->SetFillColor(0);
        gc->SetMarkerSize(0);
        
        if(graphNum==0){ gc->SetLineStyle(7); }
        
        if(twotheta==0){
            gc->GetXaxis()->SetLimits(0.4, 0.6);
            gc->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
        }
        else{
            gc->GetXaxis()->SetLimits(0.8, 1.0);
            gc->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
        }
        gc->GetYaxis()->SetRangeUser(0.0015, 0.004);
        
        gc->GetYaxis()->SetTitle("#Delta m^{2}_{32}");
        
        if(color.Contains("lue")){
            gc->SetLineColor(kBlue);
            gc->SetMarkerColor(kBlue);
        }
        else if(color.Contains("ed")){
            gc->SetLineColor(kRed);
            gc->SetMarkerColor(kRed);
        }
        else if(color.Contains("ack")){
            gc->SetLineColor(kBlack);
            gc->SetMarkerColor(kBlack);
        }
        else if(color.Contains("ite")){
            gc->SetLineColor(kWhite);
            gc->SetMarkerColor(kWhite);
        }
        
        gVect.push_back(gc);
    }
    return gVect;
}

//#############################################################################
void PreparePlot(TH2D* plot){
    plot->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    plot->GetYaxis()->SetTitle("#Delta m^{2}_{31} [eV^{2}/c^{4}]");
    TGaxis* gx = (TGaxis*)plot->GetXaxis();
    gx->SetMaxDigits(3);
    TGaxis* gy = (TGaxis*)plot->GetYaxis();
    gy->SetMaxDigits(3);
}
//#############################################################################
void PrepareBestPoint(TGraph* gr){
    gr->SetMarkerStyle(34);
    gr->SetMarkerColor(kRed);
    gr->SetMarkerSize(0.5);
}

void SetStyleVariables(TStyle *t2kStyle){
    
    t2kStyle->SetFrameBorderMode(0);
    t2kStyle->SetCanvasBorderMode(0);
    t2kStyle->SetPadBorderMode(0);
    t2kStyle->SetPadColor(0);
    t2kStyle->SetCanvasColor(0);
    t2kStyle->SetStatColor(0);
    t2kStyle->SetFillColor(0);
    t2kStyle->SetLegendBorderSize(1);
    
    t2kStyle->SetPaperSize(20,26);
    t2kStyle->SetPadTopMargin(0.05);
    t2kStyle->SetPadRightMargin(0.15); //0.05
    t2kStyle->SetPadBottomMargin(0.16);
    t2kStyle->SetPadLeftMargin(0.13);
    
    t2kStyle->SetTextFont(132);
    t2kStyle->SetTextSize(0.08);
    t2kStyle->SetLabelFont(132,"x");
    t2kStyle->SetLabelFont(132,"y");
    t2kStyle->SetLabelFont(132,"z");
    t2kStyle->SetLabelSize(0.05,"x");
    t2kStyle->SetTitleSize(0.06,"x");
    t2kStyle->SetLabelSize(0.05,"y");
    t2kStyle->SetTitleSize(0.06,"y");
    t2kStyle->SetLabelSize(0.05,"z");
    t2kStyle->SetTitleSize(0.06,"z");
    t2kStyle->SetLabelFont(132,"t");
    t2kStyle->SetTitleFont(132,"x");
    t2kStyle->SetTitleFont(132,"y");
    t2kStyle->SetTitleFont(132,"z");
    t2kStyle->SetTitleFont(132,"t");
    t2kStyle->SetTitleFillColor(0);
    t2kStyle->SetTitleX(0.25);
    t2kStyle->SetTitleFontSize(0.08);
    t2kStyle->SetTitleFont(132,"pad");
    
    t2kStyle->SetPadGridX(true);
    t2kStyle->SetPadGridY(true);
    
    
    t2kStyle->SetHistLineWidth(1.85);
    t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    
    
    t2kStyle->SetOptTitle(0);
    t2kStyle->SetOptStat(0);
    t2kStyle->SetOptFit(0);
    
    t2kStyle->SetPadTickX(1);
    t2kStyle->SetPadTickY(1);
    
    t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                     NCont);
    t2kStyle->SetNumberContours(NCont);
    
}


