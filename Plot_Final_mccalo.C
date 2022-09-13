


#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TDirectory.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>                                                                                                                                                                                                         
#include <cstdio>   // needed for io       
using namespace std;
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TMath.h>
#include <vector>
#include <TGraphAsymmErrors.h>
#endif


Float_t getDPHI(Float_t phi1, Float_t phi2){
  Float_t dphi = phi1 - phi2;

  if(dphi > TMath::Pi())
    dphi = dphi - 2. * TMath::Pi();
  if(dphi <= -TMath::Pi())
    dphi = dphi + 2. * TMath::Pi();

  if(TMath::Abs(dphi) > TMath::Pi()) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
  }

  return dphi;
}


Float_t getDR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2){
  Float_t theDphi = getDPHI(phi1, phi2);
  Float_t theDeta = eta1 - eta2;
  return TMath::Sqrt(theDphi*theDphi + theDeta*theDeta);
}

int Plot(int indexdir=0,int indexfile=0,int centmin=0, int centmax=180){


  //int bitmin=274;
  // int bitmax=281;
  // int binnum=7;

  int bitmin=259;
  int bitmax=273;
  int binnum=14;

  
  const Double_t limits[14] = {8,16,24,28,32,36,40,44,48,56,60,64,72,80};
  TH1F* hjet_all_calo;
  TH1F* hjet_fired_calo[15];
  TH1F* hjet_fired_calo_ratio[15];
  TH1F* hcentrality;
 
  hcentrality=new TH1F("centrality","centrality",200,0,200);
  hjet_all_calo=new TH1F(Form("jetpt_all_calo"),Form("jetpt_all_calo"),150,0,300);

  hjet_all_calo->Sumw2();  
  for(int j=0;j<binnum;j++){

   
    hjet_fired_calo[j]=new TH1F(Form("jetpt_trig_calo%i",j),Form("jetpt_trig_calo%i",j),150,0,300);
    hjet_fired_calo[j]->Sumw2();}

 
   
    
  cout<<indexdir<<" "<<indexfile<<endl;  
  TFile *input = TFile::Open(Form("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/mitaylor/L1MenuStudies/MinBias_Hydjet_Drum5F_5p02TeV/MinBias_Hydjet_Drum5F_5p02TeV_Run3Winter22_122X_MiniAOD_1240_HF_14_16_MuShowers_v1/220726_160306/000%d/L1Ntuple_%d.root",indexdir,indexfile));
   
  if(!input || input->IsZombie()) { cout <<"The file could not be opened!"; return 1;}
   
    TDirectory* emuDir = input->GetDirectory("l1UpgradeEmuTree");
    TDirectoryFile *hiEvtAnalyzer=(TDirectoryFile*)input->Get("hiEvtAnalyzer");
    TDirectoryFile *jetbranch=(TDirectoryFile*)input->Get("akCs4PFJetAnalyzer");



    TTree *jets=(TTree*)jetbranch->Get("t");
    TTree *ev=(TTree*)hiEvtAnalyzer->Get("HiTree");
    TTree *emuTree=(TTree*)emuDir->Get("L1UpgradeTree");


    emuTree->AddFriend(jets);
    emuTree->AddFriend(ev);

    TTreeReader emuReader(emuTree);
    //TTreeReaderValue<std::vector<bool>> emuDecision(emuReader, "m_algoDecisionInitial");
    TTreeReaderValue<Int_t> mycent(emuReader,"hiBin");
    TTreeReaderValue<Float_t> weight(emuReader,"weight");   
    TTreeReaderValue<Int_t> nref(emuReader,"nref");
    TTreeReaderValue<Int_t> ncalo(emuReader,"ncalo");
    TTreeReaderArray<Float_t> jtpt(emuReader,"jtpt");
    TTreeReaderArray<Float_t> jteta(emuReader,"jteta");
    TTreeReaderArray<Float_t> jtphi(emuReader,"jtphi");
    TTreeReaderArray<Float_t> calopt(emuReader,"calopt");

    TTreeReaderArray<Float_t> caloeta(emuReader,"caloeta");
    TTreeReaderArray<Float_t> calophi(emuReader,"calophi");
    

    TTreeReaderValue<vector<float>> emuJetPt(emuReader, "jetEt");
     TTreeReaderValue<vector<float>> emuJetEta(emuReader, "jetEta");
    TTreeReaderValue<vector<float>> emuJetPhi(emuReader, "jetPhi");

       
     if(!emuTree)return 1;
    if(!jets) return 1;


    // vector<float>* emuJetPtVec;
    

    while (emuReader.Next()) {
      //emuJetPtVec = emuJetPt.Get();
      //cout<<emuJetPtVec<<endl;  
    hcentrality->Fill(*mycent,*weight);
      if(*mycent<centmin || *mycent>centmax) continue;
      double calomax=-999;
      double etacalomax=0;
      double phicalomax=0;
      for(int l=0;l<*ncalo;l++){
	
        if(calopt[l]>calomax){
	  calomax=calopt[l];
	  etacalomax=caloeta[l];
	  phicalomax=calophi[l];}}
     
      if(calomax==-999) continue;
      if(TMath::Abs(etacalomax)>2) continue; 
      hjet_all_calo->Fill(calomax,*weight);
 
      double l1max=-999;
      
      double dr=10;
      cout<<(*emuJetPt).size()<<"size"<<" "<<calomax<<endl;

      for(long unsigned n=0;n<(*emuJetPt).size();n++){
	dr=getDR((*emuJetEta)[n],(*emuJetPhi)[n],etacalomax,phicalomax); 
       if(dr>0.4) continue;
       if((*emuJetPt)[n]>l1max){
	 l1max=(*emuJetPt)[n];}
      }
       if(l1max==-999) continue;
   
      cout<<l1max<<" "<<calomax<<endl;
      for(int l=0;l<binnum;l++){

      if(l1max>limits[l]) hjet_fired_calo[l]->Fill(calomax,*weight);
      }
   

                                                               
   



             
                

	  
    }



    



    //    delete emuDir;
    //delete hiEvtAnalyzer;
    // delete jetbranch;
    delete jets;
    delete ev;
    delete emuTree;
    input->Close();
    delete input;


    TFile* fout = new TFile(Form("/grid_mnt/data__data.polcms/cms/lcunquei/L1StudiesFramework/ReadL1Ntuples/output_mc/l1output_selection_mc_indexfile%d_indexdir%d_centmin%d_centmax%d_calo_matched.root",indexfile,indexdir,centmin,centmax),"recreate");

    //  hjet_all->Write();
  hjet_all_calo->Write();
  hcentrality->Write();
  for(int k=bitmin;k<bitmax;k++){
    
    //  hjet_fired[k-bitmin]->Write();
    //hjet_fired_ratio[k-bitmin]=(TH1F*)hjet_fired[k-bitmin]->Clone(Form("hjet_fired_ratio[%i]",k-bitmin));

    //hjet_fired_ratio[k-bitmin]->Write();
   
    hjet_fired_calo[k-bitmin]->Write();
    hjet_fired_calo_ratio[k-bitmin]=(TH1F*)hjet_fired_calo[k-bitmin]->Clone(Form("hjet_fired_calo_ratio[%i]",k-bitmin));
    hjet_fired_calo_ratio[k-bitmin]->Write();
  }

  fout->Close();


  return 1;

}

int main(int argc,char* argv[]) {
  int indexdir=0;
  int indexfile=0;
  int centmin=0;
  int centmax=180;
  if(argc!=4) cout<<"hello more arguments than needed"<<endl;
  indexdir=atoi(argv[1]);
  indexfile=atoi(argv[2]); 
  centmin=atoi(argv[3]);
  centmax=atoi(argv[4]);
   
    return Plot(indexdir,indexfile,centmin,centmax);

 
}


