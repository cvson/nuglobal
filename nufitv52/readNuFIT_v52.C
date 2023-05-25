void titleStyle(TGraph* h1);
void titleStyle(TH1* h1);
void readNuFIT_v52(){
    TString savename="nufit52_2022nov";
    const char *oscfile[1] = {"v52.release-SKoff-NO"};
    Int_t NCOMP =sizeof(oscfile)/sizeof(oscfile[0]);
    //const Int_T NLINES3D_t23dmadcp=
    Double_t th23min = 0.250;
    Double_t th23max= 0.750;
    Double_t th23NBIN = 100;//should add1
    Double_t th23STEP = (th23max-th23min)/th23NBIN;
    Double_t th23min4hist = th23min - th23STEP/2.;
    Double_t th23max4hist = th23max + th23STEP/2.;
    
    //various bin
    Double_t DMAmin = 0.2;//e-3
    Double_t DMAmax= 7.0;//e-3
    Double_t DMANBIN = 100;//should add1
    Double_t DMASTEP = (DMAmax-DMAmin)/DMANBIN;
    Double_t DMAmin4hist = DMAmin - DMASTEP/2.;
    Double_t DMAmax4hist = DMAmax + DMASTEP/2.;
    
    const char *grpattern[22]={"T23/DMA/DCP","T13/T12","T13/DMS","T12/DMS","T13/T23","T13/DMA","T23/DMA","T13/DCP","T23/DCP","DMA/DCP","T12/T23","T12/DCP","T12/DMA","DMS/T23","DMS/DCP","DMS/DMA","T13","T12","T23","DCP","DMS","DMA"};
    Int_t NGRAPH =sizeof(grpattern)/sizeof(grpattern[0]);
    Int_t aLineIndex[NGRAPH];
    Int_t totalLine;
    //Long_t index;
    Float_t th23, th13, th12;
    Float_t dms, dma;
    Float_t dcp;
    Float_t chisq;
    TFile *foutput = new TFile("nufit52_SKoff_NO.root","RECREATE");
    
    TNtuple *ptuple3d = new TNtuple("chi","nufit52","th23:dma:dcp:chisq");
    
    TNtuple *ptuple2d_t13vst12 = new TNtuple("chi_t13vst12","nufit52","th13:th12:chisq");
    TNtuple *ptuple2d_t13vsdms = new TNtuple("chi_t13vsdms","nufit52","th13:dms:chisq");
    TNtuple *ptuple2d_t12vsdms = new TNtuple("chi_t12vsdms","nufit52","th12:dms:chisq");
    TNtuple *ptuple2d_t13vst23 = new TNtuple("chi_t13vsth23","nufit52","th13:th23:chisq");
    TNtuple *ptuple2d_t13vsdma = new TNtuple("chi_t13vsdma","nufit52","th13:dma:chisq");
    TNtuple *ptuple2d_t23vsdma = new TNtuple("chi_t23vsdma","nufit52","th23:dma:chisq");
    TNtuple *ptuple2d_t13vsdcp = new TNtuple("chi_t13vsdcp","nufit52","th13:dcp:chisq");
    TNtuple *ptuple2d_t23vsdcp = new TNtuple("chi_t23vsdcp","nufit52","th23:dcp:chisq");
    TNtuple *ptuple2d_dmavsdcp = new TNtuple("chi_dmavsdcp","nufit52","dma:dcp:chisq");
    TNtuple *ptuple2d_t12vst23 = new TNtuple("chi_t12vsth23","nufit52","th12:th23:chisq");
    TNtuple *ptuple2d_t12vsdcp = new TNtuple("chi_t12vsdcp","nufit52","th12:dcp:chisq");
    TNtuple *ptuple2d_t12vsdma = new TNtuple("chi_t12vsdma","nufit52","th12:dma:chisq");
    TNtuple *ptuple2d_dmsvst23 = new TNtuple("chi_dmsvsth23","nufit52","dms:th23:chisq");
    TNtuple *ptuple2d_dmsvsdcp = new TNtuple("chi_dmsvsdcp","nufit52","dms:dcp:chisq");
    TNtuple *ptuple2d_dmsvsdma = new TNtuple("chi_dmsvsdma","nufit52","dms:dma:chisq");
    
    TNtuple *ptuple1d_t13 = new TNtuple("chi_t13","nufit52","th13:chisq");
    TNtuple *ptuple1d_t12 = new TNtuple("chi_t12","nufit52","th12:chisq");
    TNtuple *ptuple1d_t23 = new TNtuple("chi_t23","nufit52","th23:chisq");
    TNtuple *ptuple1d_dcp = new TNtuple("chi_dcp","nufit52","dcp:chisq");
    TNtuple *ptuple1d_dms = new TNtuple("chi_dms","nufit52","dms:chisq");
    TNtuple *ptuple1d_dma = new TNtuple("chi_dma","nufit52","dma:chisq");
    
    char pfilename[256];
    char grpatternname[256];
    for (Int_t icomp=0; icomp<NCOMP; ++icomp) {
        sprintf(pfilename, "%s.txt", oscfile[icomp]);
        ifstream file("v52.release-SKoff-NO.txt");
        int line_counter = 0;
        string line;
        while (getline (file,line)) {
            line_counter++;
            
            for (int igr=0; igr<NGRAPH; ++igr) {
                sprintf(grpatternname, "# %s", grpattern[igr]);
                if (line.find(grpatternname) != string::npos) {
                    cout << "Match at line " << line_counter << endl;
                    aLineIndex[igr] = line_counter;
                }
            }
        }
        totalLine = line_counter;
        file.close();
    }

    for (Int_t igr2=0; igr2<NGRAPH; ++igr2) {
        cout<<"index "<<igr2<<" val "<<aLineIndex[igr2]<<endl;
    }
    
    cout<<"line end "<<totalLine<<endl;
    
    FILE *fp = fopen(pfilename,"r");
    char buffer[100];
    int ncols;
    //skip dumb line
    for (Int_t iline=0; iline<aLineIndex[0]; ++iline){
     fgets(buffer, 100, fp);
     }
     
     //T23/DMA/DCP projection
     for (Int_t iline=aLineIndex[0]; iline<aLineIndex[1]-2; ++iline) {
         ncols = fscanf(fp,"%f %f %f %f",&th23,&dma,&dcp,&chisq);
         ptuple3d->Fill(th23,dma,dcp,chisq);
     }
    fclose(fp);
    //------------end-----------
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[1]; ++iline){
     fgets(buffer, 100, fp);
     }
    
    for (Int_t iline=aLineIndex[1]; iline<aLineIndex[2]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th13,&th12,&chisq);
        ptuple2d_t13vst12->Fill(th13,th12,chisq);
    }
    fclose(fp);
    //------------end-----------
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[2]; ++iline){
     fgets(buffer, 100, fp);
     }
    
    for (Int_t iline=aLineIndex[2]; iline<aLineIndex[3]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th13,&dms,&chisq);
        ptuple2d_t13vsdms->Fill(th13,dms,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[3]; ++iline){
     fgets(buffer, 100, fp);
     }
    
    for (Int_t iline=aLineIndex[3]; iline<aLineIndex[4]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th12,&dms,&chisq);
        ptuple2d_t12vsdms->Fill(th12,dms,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[4]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[4]; iline<aLineIndex[5]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th13,&th23,&chisq);
        ptuple2d_t13vst23 ->Fill(th13,th23,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[5]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[5]; iline<aLineIndex[6]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th13,&dma,&chisq);
        ptuple2d_t13vsdma ->Fill(th13,dma,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[6]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[6]; iline<aLineIndex[7]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th23,&dma,&chisq);
        ptuple2d_t23vsdma ->Fill(th23,dma,chisq);
    }
    fclose(fp);
    //------------end-----------
     
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[7]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[7]; iline<aLineIndex[8]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th13,&dcp,&chisq);
        ptuple2d_t13vsdcp ->Fill(th13,dcp,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[8]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[8]; iline<aLineIndex[9]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th23,&dcp,&chisq);
        ptuple2d_t23vsdcp ->Fill(th23,dcp,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[9]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[9]; iline<aLineIndex[10]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&dma,&dcp,&chisq);
        ptuple2d_dmavsdcp ->Fill(dma,dcp,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[10]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[10]; iline<aLineIndex[11]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th12,&th23,&chisq);
        ptuple2d_t12vst23 ->Fill(th12,th23,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[11]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[11]; iline<aLineIndex[12]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th12,&dcp,&chisq);
        ptuple2d_t12vsdcp ->Fill(th12,dcp,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[12]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[12]; iline<aLineIndex[13]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&th12,&dma,&chisq);
        ptuple2d_t12vsdma ->Fill(th12,dma,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[13]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[13]; iline<aLineIndex[14]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&dms,&th23,&chisq);
        ptuple2d_dmsvst23 ->Fill(dms,th23,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[14]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[14]; iline<aLineIndex[15]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&dms,&dcp,&chisq);
        ptuple2d_dmsvsdcp ->Fill(dms,dcp,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Int_t iline=0; iline<aLineIndex[15]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Int_t iline=aLineIndex[15]; iline<aLineIndex[16]-2; ++iline) {
        ncols = fscanf(fp,"%f %f %f",&dms,&dma,&chisq);
        ptuple2d_dmsvsdma ->Fill(dms,dma,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Long_t iline=0; iline<aLineIndex[16]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Long_t iline=aLineIndex[16]; iline<aLineIndex[17]-2; ++iline) {
        ncols = fscanf(fp,"%f %f",&th13,&chisq);
        ptuple1d_t13 ->Fill(th13,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Long_t iline=0; iline<aLineIndex[17]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Long_t iline=aLineIndex[17]; iline<aLineIndex[18]-2; ++iline) {
        ncols = fscanf(fp,"%f %f",&th12,&chisq);
        ptuple1d_t12 ->Fill(th12,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Long_t iline=0; iline<aLineIndex[18]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Long_t iline=aLineIndex[18]; iline<aLineIndex[19]-2; ++iline) {
        ncols = fscanf(fp,"%f %f",&th23,&chisq);
        ptuple1d_t23 ->Fill(th23,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Long_t iline=0; iline<aLineIndex[19]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Long_t iline=aLineIndex[19]; iline<aLineIndex[20]-2; ++iline) {
        ncols = fscanf(fp,"%f %f",&dcp,&chisq);
        ptuple1d_dcp ->Fill(dcp,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Long_t iline=0; iline<aLineIndex[20]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Long_t iline=aLineIndex[20]; iline<aLineIndex[21]-2; ++iline) {
        ncols = fscanf(fp,"%f %f",&dms,&chisq);
        ptuple1d_dms ->Fill(dms,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    fp = fopen(pfilename,"r");
    //skip dummy line
    for (Long_t iline=0; iline<aLineIndex[21]; ++iline){
     fgets(buffer, 100, fp);
     }
    for (Long_t iline=aLineIndex[21]; iline<totalLine-1; ++iline) {
        ncols = fscanf(fp,"%f %f",&dma,&chisq);
        ptuple1d_dma ->Fill(dma,chisq);
    }
    fclose(fp);
    //------------end-----------
    
    
     
    ptuple3d->Print();
    foutput->cd();
    ptuple3d->Write();
    ptuple2d_t13vst12->Write();
    ptuple2d_t13vsdms->Write();
    ptuple2d_t12vsdms->Write();
    ptuple2d_t13vst23->Write();
    ptuple2d_t13vsdma->Write();
    ptuple2d_t23vsdma->Write();
    ptuple2d_t13vsdcp->Write();
    ptuple2d_t23vsdcp->Write();
    ptuple2d_dmavsdcp->Write();
    ptuple2d_t12vst23->Write();
    ptuple2d_t12vsdcp->Write();
    ptuple2d_t12vsdma->Write();
    ptuple2d_dmsvst23->Write();
    ptuple2d_dmsvsdcp->Write();
    ptuple2d_dmsvsdma->Write();
    ptuple1d_t13->Write();
    ptuple1d_t12->Write();
    ptuple1d_t23->Write();
    ptuple1d_dcp->Write();
    ptuple1d_dms->Write();
    ptuple1d_dma->Write();
    
    foutput->Close();
    
    
    
}

void titleStyle(TGraph* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.1);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.1);
    
    h1->GetYaxis()->SetTitleOffset(1.0);
    h1->GetXaxis()->SetTitleOffset(0.9);
}

void titleStyle(TH1* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.1);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.1);
    
    h1->GetYaxis()->SetTitleOffset(1.0);
    h1->GetXaxis()->SetTitleOffset(0.9);
}


