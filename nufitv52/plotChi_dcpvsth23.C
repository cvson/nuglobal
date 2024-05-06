void titleStyleGraph(TGraph* h1);
#include <vector>
std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int graphNum);
void plotChi_dcpvsth23(){
    TString ntuplename= "chi_t23vsdcp";
   
    //for NuFit 5.2
    TFile *fnufit = new TFile("nufitv52/nufit52_SKoff_NO.root");
    //chi_t23;1
    TNtuple *ptuple = (TNtuple*)fnufit->Get(Form("%s",ntuplename.Data()));
    ptuple->Draw("dcp:th23:chisq","","goff");
    
    Int_t nentries =ptuple->GetSelectedRows();
    Double_t* a_dcp = ptuple->GetV1();
    Double_t* a_th23 = ptuple->GetV2();
    Double_t* a_chisq = ptuple->GetV3();
    //73x101
    Double_t dcp_max = TMath::MaxElement(nentries,a_dcp);
    Double_t dcp_min = TMath::MinElement(nentries,a_dcp);
    
    Int_t NbinY = 73;//(dcp_max-dcp_min)/5. + 1;
    Double_t StepY = (dcp_max-dcp_min)/(NbinY-1);
    cout<<"dcp max "<<dcp_max<<" min "<<dcp_min<<" step "<<StepY<<endl;
    
    Double_t sinsqth23_max = TMath::MaxElement(nentries,a_th23);
    Double_t sinsqth23_min = TMath::MinElement(nentries,a_th23);
    Int_t NbinX = 101;
    Double_t StepX = (sinsqth23_max-sinsqth23_min)/(NbinX-1);
    cout<<"sinsqth23 max "<<sinsqth23_max<<" min "<<sinsqth23_min<<" step "<<StepX<<endl;
    Double_t xrange_min = sinsqth23_min-StepX/2.;
    Double_t xrange_max = sinsqth23_max+StepX/2.;
    Double_t yrange_min = dcp_min-StepY/2.;
    Double_t yrange_max = dcp_max+StepY/2.;
    
    Double_t chisq_min = TMath::MinElement(nentries,a_chisq);
    cout<<"chisq_min "<<chisq_min<<endl;
    
    
    TH2D *hchisq = new TH2D("hchisq","",NbinX,xrange_min,xrange_max,NbinY,yrange_min,yrange_max);
    Int_t NBinCosdcp = 41;
    TH2D *hchisq_cosdcp = new TH2D("hchisq_cosdcp","",NbinX,xrange_min,xrange_max,41,-1.,1.);
    Double_t StepY_cosdcp = 2./(NBinCosdcp-1);
    for (Int_t ibinx=1; ibinx<=NbinX; ++ibinx) {
        for (Int_t ibiny=1; ibiny<=NBinCosdcp; ++ibiny) {
            hchisq_cosdcp->SetBinContent(ibinx,ibiny,99999);
        }
    }
    
    for (Int_t ibin=0; ibin<nentries; ++ibin) {
        Int_t thbinX =lround((a_th23[ibin]-sinsqth23_min)/StepX)+1;
        Int_t thbinY=lround((a_dcp[ibin]-dcp_min)/StepY)+1;
        hchisq->SetBinContent(thbinX,thbinY,a_chisq[ibin]);
        
        Int_t thbinY_dcp=lround((TMath::Cos(a_dcp[ibin]*TMath::Pi()/180.)+1)/StepY_cosdcp)+1;
        Double_t tmpBin = hchisq_cosdcp->GetBinContent(thbinX,thbinY_dcp);
        if (tmpBin>a_chisq[ibin]) {
            hchisq_cosdcp->SetBinContent(thbinX,thbinY_dcp,a_chisq[ibin]);
        }
    }
    
    //linear interpolation for the bin without filling
    for (Int_t ibinx=1; ibinx<=NbinX; ++ibinx) {
        for (Int_t ibiny=1; ibiny<=NBinCosdcp; ++ibiny) {
            Double_t tmpBin = hchisq_cosdcp->GetBinContent(ibinx,ibiny);
            Double_t tmpBin_dw = hchisq_cosdcp->GetBinContent(ibinx,ibiny-1);
            Double_t tmpBin_up = hchisq_cosdcp->GetBinContent(ibinx,ibiny+1);
            if(lround(tmpBin)==99999)hchisq_cosdcp->SetBinContent(ibinx,ibiny,(tmpBin_dw+tmpBin_up)/2.);
        }
    }
    
    
    new TCanvas;
    hchisq->Draw("colz");
    gPad->Print(Form("raw_%s.pdf",ntuplename.Data()));
    
    new TCanvas;
    hchisq_cosdcp->Draw("colz");
    hchisq_cosdcp->GetZaxis()->SetRangeUser(0,50);
    gPad->Print(Form("raw_%s_cosdcp.pdf",ntuplename.Data()));
    
    //3sigma contour
    const int nlevels = 2;
    double levels[2] = {6.18, 11.83};
    
    TH2D* hchisqCosdcp_clone=(TH2D*)hchisq_cosdcp->Clone("hchisqCosdcp_clone");
    hchisqCosdcp_clone->SetContour(nlevels,levels);
     std::vector<TGraph*> gVect_histCosdcp_0;
     std::vector<TGraph*> gVect_histCosdcp_1;
    TCanvas *c1 = new TCanvas();
     c1->cd();
    hchisqCosdcp_clone->Draw("cont list");
     c1->Update();
     TObjArray *contours_histCosdcp = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
     gVect_histCosdcp_0 = GetContourGraphs(contours_histCosdcp, "white",  0);
     gVect_histCosdcp_1 = GetContourGraphs(contours_histCosdcp, "black",  1);
     
    hchisqCosdcp_clone->Draw("colz");
     for(int j=0; j<gVect_histCosdcp_0.size(); j++){
         gVect_histCosdcp_0[j]->Draw("Lsame");
     }
     
     for(int j=0; j<gVect_histCosdcp_1.size(); j++){
         gVect_histCosdcp_1[j]->Draw("Lsame");
     }
    gPad->Print(Form("cont_%s_cosdcp.pdf",ntuplename.Data()));
    
    
    
    
    
    TH2D* hchisq_clone=(TH2D*)hchisq->Clone("hchisq_clone");
    hchisq_clone->SetContour(nlevels,levels);
     std::vector<TGraph*> gVect_hist_0;
     std::vector<TGraph*> gVect_hist_1;
    TCanvas *c2 = new TCanvas();
     c2->cd();
    hchisq_clone->Draw("cont list");
     c2->Update();
     TObjArray *contours_hist = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
     gVect_hist_0 = GetContourGraphs(contours_hist, "white",  0);
     gVect_hist_1 = GetContourGraphs(contours_hist, "black",  1);
     
    hchisq_clone->Draw("colz");
     for(int j=0; j<gVect_hist_0.size(); j++){
         gVect_hist_0[j]->Draw("Lsame");
     }
     
     for(int j=0; j<gVect_hist_1.size(); j++){
         gVect_hist_1[j]->Draw("Lsame");
     }
    gPad->Print(Form("cont_%s.pdf",ntuplename.Data()));
    
    //TH2D *hdcp_sinsqth23 = new TH2D("hdcp_sinsqth23","");
    //convert to cosdcp vs sinsqth23
    
    TH2D *hchisq_cosdcp_v2 = new TH2D("hchisq_cosdcp_v2","",NbinX,xrange_min,xrange_max,101,-1,1.);
    new TCanvas;
    hchisq_cosdcp_v2->Draw("AXIS");
    for(int j=0; j<gVect_hist_1.size(); j++){
        //gVect_hist_1[j]->Draw("Lsame");
        Double_t* exxcp2 =gVect_hist_1[j]->GetX();
        Double_t* eyycp2 =gVect_hist_1[j]->GetY();
        Double_t* eyycp2codcp =gVect_hist_1[j]->GetY();
        for (Int_t i=0; i<gVect_hist_1[j]->GetN(); ++i) {
            eyycp2codcp[i] = TMath::Cos(eyycp2[i]*TMath::Pi()/180.);
        }
        TGraph *pcosdcp = new TGraph(gVect_hist_1[j]->GetN(),exxcp2,eyycp2codcp);
        pcosdcp->Draw("Lsame");
    }
    gPad->Print(Form("cont_%s_cosdcp_test.pdf",ntuplename.Data()));
    
}

void titleStyleGraph(TGraph* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleOffset(0.9);
    h1->GetXaxis()->SetTitleOffset(0.9);
}

std::vector<TGraph*> GetContourGraphs(TObjArray *contours, TString color, int graphNum){
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
            gc->SetPoint(i, x, y);
        }
        gc->SetLineWidth(2);
        gc->SetFillColor(0);
        gc->SetMarkerSize(0);
        
        //if(graphNum==0){ gc->SetLineStyle(7); }
        
        
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
