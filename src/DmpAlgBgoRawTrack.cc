#include "DmpAlgBgoRawTrack.h"
#include "DmpDataBuffer.h"
#include "DmpBgoBase.h"
#include "TH2F.h"
#include "TF1.h"
#include "stdio.h"
#include "TFile.h"
#include "TMath.h"
//-------------------------------------------------------------------
DmpAlgBgoRawTrack::DmpAlgBgoRawTrack()
 :DmpVAlg("BgoRawTrack"),
  fBgoHits(0)
{
}

//-------------------------------------------------------------------
DmpAlgBgoRawTrack::~DmpAlgBgoRawTrack(){
}

void DmpAlgBgoRawTrack::Reset(){
memset(x_lc,0,sizeof(x_lc));
memset(Ex_lc,0,sizeof(Ex_lc));
memset(LC_E,0,sizeof(LC_E));
Energy=0.;

}
//-------------------------------------------------------------------
bool DmpAlgBgoRawTrack::Initialize(){
  // read input data
  fBgoHits= new DmpEvtBgoHits();
  gDataBuffer->LinkRootFile("Event/Cal/Hits",fBgoHits);
  fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/Cal/Hits"));
  if(!fBgoHits){
    gDataBuffer->LinkRootFile("Event/MCTruth/BgoFDigit",fBgoHits);
    fBgoHits = dynamic_cast<DmpEvtBgoHits*>(gDataBuffer->ReadObject("Event/MCTruth/BgoFDigit"));
  } 
  Angle_xz=new TH1D("Angle_xz","Angle_xz;Angles;Counts",200,-60,60);
  Angle_yz=new TH1D("Angle_yz","Angle_yz;Angles;Counts",200,-60,60);
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoRawTrack::ProcessThisEvent(){
  Reset();
  short nHits=fBgoHits->fGlobalBarID.size();
  short gid=0,l=0,b=0;
  for(short n=0;n<nHits;++n){
  gid=fBgoHits->fGlobalBarID[n];
  l=DmpBgoBase::GetLayerID(gid);
  b=DmpBgoBase::GetBarID(gid);
  if(l%2==0)
  Ex_lc[l]+= fBgoHits->fEnergy[n]*fBgoHits->fPosition[n].y();  
  else
  Ex_lc[l]+= fBgoHits->fEnergy[n]*fBgoHits->fPosition[n].x();

  LC_E[l]+= fBgoHits->fEnergy[n];
  Energy+= fBgoHits->fEnergy[n];
  }
  for(short nl=0;nl<14;nl++){ 
    x_lc[nl]=Ex_lc[nl]/LC_E[nl];
  }
 // Angles
    TH2D *zx=new TH2D("trackzx","trackzx",700,-100,600,900,-450,450);
    TH2D *zy=new TH2D("trackzy","trackzy",700,-100,600,900,-450,450);
   double pi=3.1415926;
   int nPoints[2]={0,0};
   for(int i=0;i<14;i++){
 if(Energy>100&&LC_E[i]>7){//3Sigma=6.9MeV

 if(i%2==0){ 
 zy->Fill(i*29.0+35,x_lc[i],LC_E[i]/Energy); 
 nPoints[0]++;
  }
 else{
 zx->Fill(i*29.0+35,x_lc[i],LC_E[i]/Energy); 
 nPoints[1]++; 
 
  }
 }
}
  TF1 *linear=new TF1("linear","[0]+[1]*x",0,450);
  double par[2]={0.,0.};
 if(nPoints[0]>=2&&nPoints[1]>=2){
 
 zx->Fit(linear,"RQ0");
 linear->GetParameters(par);
 Angle_xz->Fill(-1*TMath::ATan(par[1])/2/pi*360);
 
 zy->Fit(linear,"RQ0");
 linear->GetParameters(par);
 Angle_yz->Fill(-1*TMath::ATan(par[1])/2/pi*360);
 }
 delete zx;
 delete zy;
 delete linear;
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoRawTrack::Finalize(){
  TF1 *mygaus=new TF1("mygaus","gaus",-45,45);
  std::string histFileName = "./Angle/"+gRootIOSvc->GetInputStem()+"_Angle.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  double mean=Angle_xz->GetMean();
  mygaus->SetRange(mean-15,mean+15);
  Angle_xz->Fit(mygaus,"RQ0");
  double par[3];
  mygaus->GetParameters(par);
  mygaus->SetRange(par[1]-3*par[2],par[1]+3*par[2]);
  Angle_xz->Fit(mygaus,"RQ");
  Angle_xz->Write();
  delete Angle_xz;
  
  mean=Angle_yz->GetMean();
  mygaus->SetRange(mean-15,mean+15);
  Angle_yz->Fit(mygaus,"RQ0");

  mygaus->GetParameters(par);
  mygaus->SetRange(par[1]-3*par[2],par[1]+3*par[2]);
  Angle_yz->Fit(mygaus,"RQ");
  Angle_yz->Write();
  delete Angle_yz;

  return true;
}

