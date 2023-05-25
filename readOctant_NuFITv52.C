void titleStyleGraph(TGraph* h1);
void readOctant_NuFITv52(){
    TString savename="nufitv52_octantdata";
    TFile *foutput = new TFile(Form("op_%s.root",savename.Data()),"RECREATE");
    
    const char *Optmassordering[2] = {"no","io"};
    const char *Optusereactor[2] = {"wreactor","noreactor"};
    const char *OptExp[4]={"minos","nova","t2k","combine"};
    const char *LegOptExp[4]={"MINOS","NOvA","T2K","Combined"};
    
    const char *colorcode[] = {
            "#000000",
            "#0072B2",
            "#D55E00",
            "#CC79A7",
            "#E69F00",
            "#009E73",
            "#56B4E9",
            "#F0E442"
        };
        Int_t ci;

    
    Int_t NMO =sizeof(Optmassordering)/sizeof(Optmassordering[0]);
    Int_t NRT =sizeof(Optusereactor)/sizeof(Optusereactor[0]);
    Int_t NEXP = sizeof(OptExp)/sizeof(OptExp[0]);
    
    Int_t NCONFIG = NMO*NRT*NEXP;

    TTree ****ptree = new TTree***[NMO];
    TGraph ****pvoltgraph_chi = new TGraph***[NMO];
   
    char pfilename[256];
    Long64_t nline_tmp;
    
    for (Int_t ioptMO=0; ioptMO<NMO; ++ioptMO) {
        pvoltgraph_chi[ioptMO]= new TGraph**[NRT];
        ptree[ioptMO] = new TTree**[NRT];
        
        for (Int_t ioptRT=0; ioptRT<NRT; ++ioptRT) {
            pvoltgraph_chi[ioptMO][ioptRT] = new TGraph*[NEXP];
            ptree[ioptMO][ioptRT] = new TTree*[NEXP];
            
            for (Int_t ioptEXP=0; ioptEXP<NEXP; ++ioptEXP) {
                ptree[ioptMO][ioptRT][ioptEXP] = new TTree(Form("tree%s_%s_%s",Optmassordering[ioptMO],Optusereactor[ioptRT],OptExp[ioptEXP]),"nufit52");
                
                sprintf(pfilename, "nufitv52/th23octant/%s_%s_%s.csv", Optmassordering[ioptMO],Optusereactor[ioptRT],OptExp[ioptEXP]);
                //FILE *fp = fopen(pfilename,"r");
                
                nline_tmp = ptree[ioptMO][ioptRT][ioptEXP]->ReadFile(pfilename,"th23/F:chi/F",',');
                foutput->cd();
                ptree[ioptMO][ioptRT][ioptEXP]->Write();
                ptree[ioptMO][ioptRT][ioptEXP]->Draw("chi:th23","","goff");
                pvoltgraph_chi[ioptMO][ioptRT][ioptEXP] = new TGraph(ptree[ioptMO][ioptRT][ioptEXP]->GetSelectedRows(),ptree[ioptMO][ioptRT][ioptEXP]->GetV2(), ptree[ioptMO][ioptRT][ioptEXP]->GetV1());
                ci = TColor::GetColor(colorcode[ioptEXP]);
                pvoltgraph_chi[ioptMO][ioptRT][ioptEXP] ->SetLineColor(ci);
                pvoltgraph_chi[ioptMO][ioptRT][ioptEXP] ->SetLineWidth(2);
                pvoltgraph_chi[ioptMO][ioptRT][ioptEXP]->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
                pvoltgraph_chi[ioptMO][ioptRT][ioptEXP]->GetYaxis()->SetTitle("#Delta #chi^{2}");
                //pvoltgraph_chi[ioptMO][ioptRT][ioptEXP] = ;
            }
        }
    }
    
    TGraph *pgr_no_wreactor_t2kclone = (TGraph*)pvoltgraph_chi[0][0][2]->Clone("pgr_no_wreactor_t2kclone");
    
    Int_t NGRPOINT =pgr_no_wreactor_t2kclone->GetN();
    Double_t *pX = pgr_no_wreactor_t2kclone->GetX();
    Double_t *pY = pgr_no_wreactor_t2kclone->GetY();
    for (Int_t ipoint=0; ipoint<NGRPOINT; ++ipoint) {
        for (Int_t ioptEXP=0; ioptEXP<2; ++ioptEXP) {
            pY[ipoint] +=pvoltgraph_chi[0][0][ioptEXP]->Eval(pX[ipoint]);
        }
        
    }
    double chisqmin = TMath::MinElement(NGRPOINT, pY);
    for (Int_t ipoint=0; ipoint<NGRPOINT; ++ipoint) {
        pY[ipoint] -=chisqmin;
    }
    
    TGraph *pgr_no_wreactor_combine_simple = new TGraph(NGRPOINT,pX,pY);
    ci = TColor::GetColor(colorcode[NEXP]);
    pgr_no_wreactor_combine_simple ->SetLineColor(ci);
    pgr_no_wreactor_combine_simple ->SetLineWidth(2);
    
    
    TLegend* leg0 = new TLegend(.38, .56, 0.65, .88);
           leg0->SetFillStyle(0);
           leg0->SetBorderSize(0);
           leg0->SetTextSize(24);
           leg0->SetTextFont(43);
           leg0->SetMargin(0.2);
    new TCanvas;
    pvoltgraph_chi[0][0][0]->Draw("APL");
    pvoltgraph_chi[0][0][0]->GetXaxis()->SetLimits(0.3,0.7);
    pvoltgraph_chi[0][0][0]->GetYaxis()->SetRangeUser(0.0,15.0);
    leg0->AddEntry(pvoltgraph_chi[0][0][0], "T2K, NuFIT 5.2");
    titleStyleGraph(pvoltgraph_chi[0][0][0]);
    for (Int_t ioptEXP=1; ioptEXP<NEXP; ++ioptEXP) {
        pvoltgraph_chi[0][0][ioptEXP]->Draw("L same");
        leg0->AddEntry( pvoltgraph_chi[0][0][ioptEXP], LegOptExp[ioptEXP]);
    }
    pgr_no_wreactor_combine_simple->SetLineStyle(8);
    pgr_no_wreactor_combine_simple->Draw("L same");
    leg0->AddEntry( pgr_no_wreactor_combine_simple, "#chi^{2} sum");
    leg0->Draw();
    gPad->Print(Form("plots/%s_no_wreactor.pdf",savename.Data()));
    
    // no without reactor
    TGraph *pgr_no_woreactor_t2kclone = (TGraph*)pvoltgraph_chi[0][1][2]->Clone("pgr_no_wreactor_t2kclone");
    
    Int_t NGRPOINT_woR =pgr_no_woreactor_t2kclone->GetN();
    Double_t *pX_woR = pgr_no_woreactor_t2kclone->GetX();
    Double_t *pY_woR = pgr_no_woreactor_t2kclone->GetY();
    for (Int_t ipoint=0; ipoint<NGRPOINT_woR; ++ipoint) {
        for (Int_t ioptEXP=0; ioptEXP<2; ++ioptEXP) {
            pY_woR[ipoint] +=pvoltgraph_chi[0][1][ioptEXP]->Eval(pX_woR[ipoint]);
        }
        
    }
    double chisqmin_woR = TMath::MinElement(NGRPOINT_woR, pY_woR);
    for (Int_t ipoint=0; ipoint<NGRPOINT_woR; ++ipoint) {
        pY_woR[ipoint] -=chisqmin_woR;
    }
    
    TGraph *pgr_no_woreactor_combine_simple = new TGraph(NGRPOINT_woR,pX_woR,pY_woR);
    ci = TColor::GetColor(colorcode[NEXP]);
    pgr_no_woreactor_combine_simple ->SetLineColor(ci);
    pgr_no_woreactor_combine_simple ->SetLineWidth(2);
    
    new TCanvas;
    pvoltgraph_chi[0][1][0]->Draw("APL");
    pvoltgraph_chi[0][1][0]->GetXaxis()->SetLimits(0.3,0.7);
    pvoltgraph_chi[0][1][0]->GetYaxis()->SetRangeUser(0.0,15.0);
    titleStyleGraph(pvoltgraph_chi[0][1][0]);
    for (Int_t ioptEXP=1; ioptEXP<NEXP; ++ioptEXP) {
        pvoltgraph_chi[0][1][ioptEXP]->Draw("L same");
    }
    pgr_no_woreactor_combine_simple->SetLineStyle(8);
    pgr_no_woreactor_combine_simple->Draw("L same");
    leg0->Draw();
    gPad->Print(Form("plots/%s_no_woreactor.pdf",savename.Data()));
    
    //io without reactor
    TGraph *pgr_io_woreactor_t2kclone = (TGraph*)pvoltgraph_chi[1][1][2]->Clone("pgr_io_wreactor_t2kclone");
    
    Int_t NGRPOINT_io_woR =pgr_io_woreactor_t2kclone->GetN();
    Double_t *pX_io_woR = pgr_io_woreactor_t2kclone->GetX();
    Double_t *pY_io_woR = pgr_io_woreactor_t2kclone->GetY();
    for (Int_t ipoint=0; ipoint<NGRPOINT_io_woR; ++ipoint) {
        for (Int_t ioptEXP=0; ioptEXP<2; ++ioptEXP) {
            pY_io_woR[ipoint] +=pvoltgraph_chi[1][1][ioptEXP]->Eval(pX_io_woR[ipoint]);
        }
        
    }
    double chisqmin_io_woR = TMath::MinElement(NGRPOINT_io_woR, pY_io_woR);
    for (Int_t ipoint=0; ipoint<NGRPOINT_io_woR; ++ipoint) {
        pY_io_woR[ipoint] -=chisqmin_io_woR;
    }
    
    TGraph *pgr_io_woreactor_combine_simple = new TGraph(NGRPOINT_io_woR,pX_io_woR,pY_io_woR);
    ci = TColor::GetColor(colorcode[NEXP]);
    pgr_io_woreactor_combine_simple ->SetLineColor(ci);
    pgr_io_woreactor_combine_simple ->SetLineWidth(2);
    
    new TCanvas;
    pvoltgraph_chi[1][1][0]->Draw("APL");
    pvoltgraph_chi[1][1][0]->GetXaxis()->SetLimits(0.3,0.7);
    pvoltgraph_chi[1][1][0]->GetYaxis()->SetRangeUser(0.,15.);
    titleStyleGraph(pvoltgraph_chi[1][1][0]);
    for (Int_t ioptEXP=1; ioptEXP<NEXP; ++ioptEXP) {
        pvoltgraph_chi[1][1][ioptEXP]->Draw("L same");
    }
    pgr_io_woreactor_combine_simple->SetLineStyle(8);
    pgr_io_woreactor_combine_simple->Draw("L same");
    leg0->Draw();
    gPad->Print(Form("plots/%s_io_woreactor.pdf",savename.Data()));
    
    //copy of T2K
    //TGraph *pgr_no_wreactor_combine_simple = (TGraph*)pvoltgraph_chi[0][0][2]->Clone("pgr_no_wreactor_combine_simple");
    
    //io with reactor
    TGraph *pgr_io_wreactor_t2kclone = (TGraph*)pvoltgraph_chi[1][0][2]->Clone("pgr_io_wreactor_t2kclone");
    
    Int_t NGRPOINT_io_wR =pgr_io_wreactor_t2kclone->GetN();
    Double_t *pX_io_wR = pgr_io_wreactor_t2kclone->GetX();
    Double_t *pY_io_wR = pgr_io_wreactor_t2kclone->GetY();
    for (Int_t ipoint=0; ipoint<NGRPOINT_io_wR; ++ipoint) {
        for (Int_t ioptEXP=0; ioptEXP<2; ++ioptEXP) {
            pY_io_wR[ipoint] +=pvoltgraph_chi[1][0][ioptEXP]->Eval(pX_io_wR[ipoint]);
        }
        
    }
    double chisqmin_io_wR = TMath::MinElement(NGRPOINT_io_wR, pY_io_wR);
    for (Int_t ipoint=0; ipoint<NGRPOINT_io_wR; ++ipoint) {
        pY_io_wR[ipoint] -=chisqmin_io_wR;
    }
    
    TGraph *pgr_io_wreactor_combine_simple = new TGraph(NGRPOINT_io_wR,pX_io_wR,pY_io_wR);
    ci = TColor::GetColor(colorcode[NEXP]);
    pgr_io_wreactor_combine_simple ->SetLineColor(ci);
    pgr_io_wreactor_combine_simple ->SetLineWidth(2);
    
    new TCanvas;
    pvoltgraph_chi[1][0][0]->Draw("APL");
    pvoltgraph_chi[1][0][0]->GetXaxis()->SetLimits(0.3,0.7);
    pvoltgraph_chi[1][0][0]->GetYaxis()->SetRangeUser(0.0,15.0);
    titleStyleGraph(pvoltgraph_chi[1][0][0]);
    for (Int_t ioptEXP=1; ioptEXP<NEXP; ++ioptEXP) {
        pvoltgraph_chi[1][0][ioptEXP]->Draw("L same");
    }
    pgr_io_wreactor_combine_simple->SetLineStyle(8);
    pgr_io_wreactor_combine_simple->Draw("L same");
    leg0->Draw();
    gPad->Print(Form("plots/%s_io_wreactor.pdf",savename.Data()));
    
    
    
    foutput->Close();
    
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
