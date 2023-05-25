void titleStyle(TGraph* h1);
void readOctant_frExps_v2(){
    Bool_t isIncludingSum = true;
    Bool_t isIncludingIcecube = true;
    Bool_t isIncludingSK = true;
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
    //for NuFit 5.2
    TFile *fnufit = new TFile("nufit52_SKoff_NO.root");
    //chi_t23;1
    TNtuple *ptuple1d_t23 = (TNtuple*)fnufit->Get("chi_t23");
    ptuple1d_t23->Draw("chisq:th23","","goff");
    TGraph *grnufit = new TGraph( ptuple1d_t23->GetSelectedRows(), ptuple1d_t23->GetV2(),  ptuple1d_t23->GetV1());
    
    
    //,"nufit52","th23:chisq"
    
    //read t2k 2023, using data until t2k 2020
    //https://zenodo.org/record/7929975#.ZGwsbS0RpbU
    //https://zenodo.org/record/7741399#.ZGwx8i0RpbW
    TFile *ft2k = new TFile("t2k2022summer/T2K_arxiv2303/Frequentist_DataRelease.root");
    const char *t2khistname[4] = {"h1D_th23chi2_woRC_NH","h1D_th23chi2_woRC_IH","h1D_th23chi2_wRC_NH","h1D_th23chi2_wRC_IH"};
    
    Int_t NT2KHIST =sizeof(t2khistname)/sizeof(t2khistname[0]);
    TH1D **ht2khist = new TH1D*[NT2KHIST];
    TGraph **ht2kgraph = new TGraph*[NT2KHIST];
    
    for (Int_t it2khist=0; it2khist<NT2KHIST; ++it2khist) {
        ht2khist[it2khist] = (TH1D*)ft2k->Get(Form("%s",t2khistname[it2khist]));
        ht2kgraph[it2khist] = new TGraph(ht2khist[it2khist]);
        ht2kgraph[it2khist]->SetLineWidth(2);
    }
    
    //read Nova, release with RC only
    //nova is significance, must be squared
    TFile *fnova = new TFile("NOvA_2020_data_release/NOvA_2020_official_slice_ssth23.root");
    const char *novagrname[2] = {"NO","IO"};
    
    Int_t NNOVAGRAPH=sizeof(novagrname)/sizeof(novagrname[0]);
    TGraph **grnova = new TGraph*[NNOVAGRAPH];
    TGraph **grnovaChisq = new TGraph*[NNOVAGRAPH];
    for (Int_t inovagr=0; inovagr<NNOVAGRAPH; ++inovagr) {
        grnova[inovagr] = (TGraph*)fnova->Get(Form("%s",novagrname[inovagr]));
        Int_t NGRPOINT =grnova[inovagr]->GetN();
        Double_t *pX = grnova[inovagr]->GetX();
        Double_t *pY = grnova[inovagr]->GetY();
        for (Int_t ipoint=0; ipoint<NGRPOINT; ++ipoint) {
            pY[ipoint] *=pY[ipoint];//chisquare
        }
        grnovaChisq[inovagr] = new TGraph(NGRPOINT,pX,pY);
        grnovaChisq[inovagr]->SetLineWidth(2);
    }
    
    
    //read MINOS
    const char *fminos[2] = {"minos2021_no_wreactor","minos2021_io_wreactor"};
    Int_t NMINOS_OPT =sizeof(fminos)/sizeof(fminos[0]);
    TTree **ptree_minos = new TTree*[NMINOS_OPT];
    TGraph **pchigraph_minos = new TGraph*[NMINOS_OPT];
    char pfilename[256];
    Long64_t nline_tmp;
    
    for (Int_t iminosopt=0; iminosopt<NMINOS_OPT; ++iminosopt) {
        ptree_minos[iminosopt] = new TTree(Form("tree%s",fminos[iminosopt]),"minos2021");
        sprintf(pfilename, "minos2021/%s.csv", fminos[iminosopt]);
        nline_tmp = ptree_minos[iminosopt]->ReadFile(pfilename,"th23/F:chi/F",',');

        ptree_minos[iminosopt]->Draw("chi:th23","","goff");
        pchigraph_minos[iminosopt] = new TGraph( ptree_minos[iminosopt]->GetSelectedRows(), ptree_minos[iminosopt]->GetV2(),  ptree_minos[iminosopt]->GetV1());
        
        pchigraph_minos[iminosopt]->SetLineWidth(2);
    }
    //Read Super-K
    const char *fsuperk[4] = {"sk2019_no_woreactor","sk2019_io_woreactor","sk2019_no_wreactor","sk2019_io_wreactor"};
    Int_t NSUPERK_OPT =sizeof(fsuperk)/sizeof(fsuperk[0]);
    TTree **ptree_superk = new TTree*[NSUPERK_OPT];
    TGraph **pchigraph_superk = new TGraph*[NSUPERK_OPT];
    
    for (Int_t isuperkopt=0; isuperkopt<NSUPERK_OPT; ++isuperkopt) {
        ptree_superk[isuperkopt] = new TTree(Form("tree%s",fsuperk[isuperkopt]),"superk2019");
        sprintf(pfilename, "superk2019/%s.csv", fsuperk[isuperkopt]);
        nline_tmp = ptree_superk[isuperkopt]->ReadFile(pfilename,"th23/F:chi/F",',');

        ptree_superk[isuperkopt]->Draw("chi:th23","","goff");
        pchigraph_superk[isuperkopt] = new TGraph( ptree_superk[isuperkopt]->GetSelectedRows(), ptree_superk[isuperkopt]->GetV2(),  ptree_superk[isuperkopt]->GetV1());
        
        pchigraph_superk[isuperkopt]->SetLineWidth(2);
        
    }
    
    //Read IceCube
    const char *ficecube[1] = {"icecube2023"};
    Int_t NICECUBE_OPT =sizeof(ficecube)/sizeof(ficecube[0]);
    TTree **ptree_icecube = new TTree*[NICECUBE_OPT];
    TGraph **pchigraph_icecube = new TGraph*[NICECUBE_OPT];
  
    for (Int_t iicecubeopt=0; iicecubeopt<NICECUBE_OPT; ++iicecubeopt) {
        ptree_icecube[iicecubeopt] = new TTree(Form("tree%s",ficecube[iicecubeopt]),"icecube2021");
        sprintf(pfilename, "icecube2023/%s.csv", ficecube[iicecubeopt]);
        nline_tmp = ptree_icecube[iicecubeopt]->ReadFile(pfilename,"th23/F:chi/F",',');

        ptree_icecube[iicecubeopt]->Draw("chi:th23","","goff");
        pchigraph_icecube[iicecubeopt] = new TGraph( ptree_icecube[iicecubeopt]->GetSelectedRows(), ptree_icecube[iicecubeopt]->GetV2(),  ptree_icecube[iicecubeopt]->GetV1());
        
        pchigraph_icecube[iicecubeopt]->SetLineWidth(2);
        
    }
    
    TLegend* leg0 = new TLegend(.38, .46, 0.65, .88);
           leg0->SetFillStyle(0);
           leg0->SetBorderSize(0);
           leg0->SetTextSize(22);
           leg0->SetTextFont(43);
           leg0->SetMargin(0.15);
    leg0->AddEntry((TObject*)0, "Normal m_{#nu} ordering", "");
    
    new TCanvas;
    gPad->SetBottomMargin(gPad->GetBottomMargin()*1.2);
    ht2kgraph[2]->Draw("APL");
    ht2kgraph[2]->GetYaxis()->SetTitle("#Delta #chi^{2}");
    ht2kgraph[2]->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    titleStyle(ht2kgraph[2]);
    
    ci = TColor::GetColor(colorcode[1]);
    ht2kgraph[2] ->SetLineColor(ci);
    //ht2kgraph[0]->Draw("L same");
    ht2kgraph[2]->GetYaxis()->SetRangeUser(0,15);
    ht2kgraph[2]->GetXaxis()->SetLimits(0.3,0.7);
    leg0->AddEntry(ht2kgraph[2], "T2K (2023)");
    
    //nova
    grnovaChisq[0]->Draw("L same");
    ci = TColor::GetColor(colorcode[2]);
    grnovaChisq[0] ->SetLineColor(ci);
    leg0->AddEntry(grnovaChisq[0], "NOvA (2020)");
    
    //minos
    pchigraph_minos[0]->Draw("L same");
    ci = TColor::GetColor(colorcode[3]);
    pchigraph_minos[0] ->SetLineColor(ci);
    leg0->AddEntry(pchigraph_minos[0], "MINOS (2020)");
    
    //superk
    pchigraph_superk[2]->Draw("L same");
    ci = TColor::GetColor(colorcode[4]);
    pchigraph_superk[2] ->SetLineColor(ci);
    leg0->AddEntry(pchigraph_superk[2], "Super-K (2019)");
    
    //icecube
    pchigraph_icecube[0]->Draw("L same");
    ci = TColor::GetColor(colorcode[5]);
    pchigraph_icecube[0] ->SetLineColor(ci);
    leg0->AddEntry(pchigraph_icecube[0], "IceCube (2023)");
    
    //make the sum
    TGraph *pgr_no_wreactor_t2kclone = (TGraph*)ht2kgraph[2]->Clone("pgr_no_wreactor_t2kclone");
    
    Int_t NGRPOINT =pgr_no_wreactor_t2kclone->GetN();
    Double_t *pX = pgr_no_wreactor_t2kclone->GetX();
    Double_t *pY = pgr_no_wreactor_t2kclone->GetY();
    for (Int_t ipoint=0; ipoint<NGRPOINT; ++ipoint) {
        //add nova
        pY[ipoint] +=grnovaChisq[0]->Eval(pX[ipoint]);
        //add minos
        pY[ipoint] +=pchigraph_minos[0]->Eval(pX[ipoint]);
        //add SK
        if (isIncludingSK) {
            pY[ipoint] +=pchigraph_superk[2]->Eval(pX[ipoint]);
        }
        //add Icecube
        if (isIncludingIcecube) {
            pY[ipoint] +=pchigraph_icecube[0]->Eval(pX[ipoint]);
        }
        
    }
    double chisqmin = TMath::MinElement(NGRPOINT, pY);
    for (Int_t ipoint=0; ipoint<NGRPOINT; ++ipoint) {
        pY[ipoint] -=chisqmin;
    }
    
    //TGraph *pgr_no_wreactor_combine_simple = new TGraph(NGRPOINT,pX,pY);
    ci = TColor::GetColor(colorcode[0]);
    grnufit ->SetLineColor(ci);
    grnufit ->SetLineWidth(2);
    grnufit->SetLineStyle(8);
    
    
    if(!isIncludingSum){
        leg0->Draw();
        gPad->Print("plots/octant_frExps.pdf");
    }
    else {
        grnufit->Draw("L same");
        leg0->AddEntry(grnufit, "NuFIT 5.2");
        leg0->Draw();
        if(isIncludingIcecube)gPad->Print("plots/octant_frExps_wsum_wIC.pdf");
        else {gPad->Print("plots/octant_frExps_wsum_woIC.pdf");
            if (!isIncludingSK) {
                gPad->Print("plots/octant_frExps_wsum_woIC_woSK.pdf");
            }
        }
    }
    
    //make a sum
}

void titleStyle(TGraph* h1){
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
