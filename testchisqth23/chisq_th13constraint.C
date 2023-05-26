void titleStyleGraph(TGraph* h1);
void chisq_th13constraint(){
    TString subname= "test_cchisq_th13constraint";
    TFile *finput = new TFile("hyperk_wandwithout_external_constraint_th13.root");
    const char *ahistname[] = {"th23_nocontraint","th23_contraint_0p026sin2theta13","th23_contraint_0p01sin2theta13"};
    Int_t NGRAPH =sizeof(ahistname)/sizeof(ahistname[0]);
    TGraph **pgraph = new TGraph*[NGRAPH];
    for (Int_t igr=0; igr<NGRAPH; ++igr) {
        pgraph[igr] = (TGraph*)finput->Get(Form("%s",ahistname[igr]));
    }
    //
    double sinsqth23true = 0.572;
    double sinsqth13true = 0.02203;
    double th13true = TMath::ASin(TMath::Sqrt(sinsqth13true));
    double sigmath13_0p026 = th13true*0.013;// a half of sinsqth13 uncertainty
    double sigmath13_0p01 = th13true*0.005;
    
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
    //iso-probability sin2theta13*sinth23 = constant
    double isoProb = TMath::Sqrt(sinsqth23true)*2*TMath::Sqrt(sinsqth13true)*TMath::Sqrt(1-sinsqth13true);
    cout<<"iso prob value "<<isoProb<<endl;
    
    TGraph *grsub0p026 = (TGraph*)pgraph[1]->Clone("grsub0p026");
    Int_t NGRPOINT_0p026 =grsub0p026->GetN();
    Double_t *pX_0p026 = grsub0p026->GetX();
    Double_t *pY_0p026 = grsub0p026->GetY();
    
    Double_t *pY_0p026_chisq = new Double_t[NGRPOINT_0p026];
    Double_t tmpth13sol;
    Double_t tmpchisq;
    Double_t tmpchisqNom;
    for (Int_t ipoint=0; ipoint<NGRPOINT_0p026; ++ipoint) {
        tmpchisqNom = pgraph[0]->Eval(pX_0p026[ipoint]);
        pY_0p026[ipoint] -= tmpchisqNom;
        tmpth13sol = TMath::ASin(isoProb/TMath::Sqrt(pX_0p026[ipoint]))/2.;
        tmpchisq = (tmpth13sol-th13true)/sigmath13_0p026;
        //pY_0p026_chisq[ipoint] = pow(tmpchisq,2.);
        pY_0p026_chisq[ipoint] = 2*tmpchisq;
        //pY_0p026_chisq[ipoint] = 2*tmpchisq;//pow(tmpchisq,2);//*TMath::Prob(tmpchisqNom,1);
        //cout<<"sinsq "<< pX_0p026[ipoint]<<" chiNom "<<tmpchisqNom<<" ext "<<pY_0p026_chisq[ipoint]<<endl;
    }
    TGraph *grsub0p026_cp = new TGraph(NGRPOINT_0p026,pX_0p026,pY_0p026);
    ci = TColor::GetColor(colorcode[1]);
    grsub0p026_cp->SetLineColor(ci);
    grsub0p026_cp->SetLineWidth(2);
    
    TGraph *grchiExt0p026_cp = new TGraph(NGRPOINT_0p026,pX_0p026,pY_0p026_chisq);
    ci = TColor::GetColor(colorcode[2]);
    grchiExt0p026_cp->SetLineColor(ci);
    grchiExt0p026_cp->SetLineWidth(2);
    
    TGraph *grsub0p01 = (TGraph*)pgraph[2]->Clone("grsub0p01");
    Int_t NGRPOINT_0p01 =grsub0p01->GetN();
    Double_t *pX_0p01 = grsub0p01->GetX();
    Double_t *pY_0p01 = grsub0p01->GetY();
    
    Double_t *pY_0p01_chisq = new Double_t[NGRPOINT_0p026];
    
    for (Int_t ipoint=0; ipoint<NGRPOINT_0p01; ++ipoint) {
        tmpchisqNom = pgraph[0]->Eval(pX_0p01[ipoint]);
        pY_0p01[ipoint] -= tmpchisqNom;
        
        tmpth13sol = TMath::ASin(isoProb/TMath::Sqrt(pX_0p01[ipoint]))/2.;
        tmpchisq = (tmpth13sol-th13true)/sigmath13_0p01;
        //pY_0p01_chisq[ipoint] = pow(tmpchisq,2.);
        pY_0p01_chisq[ipoint] = 2*tmpchisq;
    }
    TGraph *grsub0p01_cp = new TGraph(NGRPOINT_0p01,pX_0p01,pY_0p01);
    ci = TColor::GetColor(colorcode[3]);
    grsub0p01_cp->SetLineColor(ci);
    grsub0p01_cp->SetLineWidth(2);
    
    TGraph *grchiExt0p01_cp = new TGraph(NGRPOINT_0p01,pX_0p01,pY_0p01_chisq);
    ci = TColor::GetColor(colorcode[4]);
    grchiExt0p01_cp->SetLineColor(ci);
    grchiExt0p01_cp->SetLineWidth(2);
    
    TLegend* leg0 = new TLegend(.35, .6, 0.85, .88);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(22);
    leg0->SetTextFont(43);
    leg0->SetMargin(0.15);
    
    
    new TCanvas;
    pgraph[0]->Draw("AL");
    leg0->AddEntry(pgraph[0],"#chi^{2}(HK only)","L");
    
    grsub0p01_cp->Draw("L same");
    leg0->AddEntry(grsub0p01_cp,"#chi^{2}(#sigma(sin^{2}#theta_{13})=1%)-#chi^{2}(HK only)");
    
    grsub0p026_cp->Draw("L same");
    leg0->AddEntry(grsub0p026_cp,"#chi^{2}(#sigma(sin^{2}#theta_{13})=2.6%)-#chi^{2}(HK only)");
    
    grchiExt0p01_cp->Draw("L same");
    leg0->AddEntry(grchiExt0p01_cp,"#chi^{2}(pred., #sigma(sin^{2}#theta_{13})=2.6%)");
    
    grchiExt0p026_cp->Draw("L same");
    leg0->AddEntry(grchiExt0p026_cp,"#chi^{2}(pred., #sigma(sin^{2}#theta_{13})=1%)");
    leg0->Draw();
    
    gPad->Print(Form("plots/%s.pdf",subname.Data()));
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
