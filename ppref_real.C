#include "uti.h"

const int nBins=10;
double ptBins[nBins+1] = {5., 10., 13., 17., 24., 30., 40., 50., 60., 70., 100.};
TString inputname_ref = "ppref_real.txt";
TString outputname_ref = "ppref_real_5p02TeV.root";
TString inputmcname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_addSamplePthat_pthatweight.root";

using namespace std;

void ppref_real()
{
  ifstream getdata(inputname_ref.Data());
  
  if(!getdata.is_open())
    {
      cout<<"Failure while opening the input file"<<endl;
    }
  
  // dsigmadpt units in [microbarn/GeV] 
  float ptmin[nBins], ptmax[nBins], yrange[nBins], dsigmadpt_real[nBins], e1[nBins], e2[nBins], e3[nBins], dsigmadpt_FONLL[nBins], e4[nBins], e5[nBins], dsigmadpt_FONLL_5p02[nBins], e6[nBins], e7[nBins], dsigmadpt_pythia[nBins], e_real[nBins], e_FONLL[nBins], e_FONLL_5p02[nBins], e_pythia[nBins], dsigmadpt_real_5p02[nBins];
  int i;
  
  for(i=0;i<nBins;i++)
    {
      getdata>>ptmin[i];
      getdata>>ptmax[i];
      getdata>>yrange[i];
      getdata>>dsigmadpt_real[i];
      getdata>>e1[i];
      getdata>>e2[i];
      getdata>>e3[i];
      getdata>>dsigmadpt_FONLL[i];
      getdata>>e4[i];
      getdata>>e5[i];
      getdata>>dsigmadpt_FONLL_5p02[i];
      getdata>>e6[i];
      getdata>>e7[i];
      getdata>>dsigmadpt_pythia[i];
      
      e_real[i] = sqrt(e1[i]*e1[i]+e2[i]*e2[i]+e3[i]*e3[i]);
      e_FONLL[i] = (e5[i]-e4[i])*0.5;
      e_FONLL_5p02[i] = (e7[i]-e6[i])*0.5;
      e_pythia[i] = 0.0;
      
      //dsigmadpt_real_5p02[i]=dsigmadpt_real[i]*dsigmadpt_FONLL_5p02[i]/dsigmadpt_FONLL[i];
    }
  
  TH1F* d_r = new TH1F("d_r","",nBins,ptBins[0],ptBins[nBins]);
  TH1F* d_f = new TH1F("d_f","",nBins,ptBins[0],ptBins[nBins]);
  TH1F* d_f5 = new TH1F("d_f5","",nBins,ptBins[0],ptBins[nBins]);
  TH1F* pythia_ref = new TH1F("pythia_ref","",nBins,ptBins[0],ptBins[nBins]);

  for(i=0;i<nBins;i++)
    {
      d_r->SetBinContent(i+1,dsigmadpt_real[i]);
      d_f->SetBinContent(i+1,dsigmadpt_FONLL[i]);
      d_f5->SetBinContent(i+1,dsigmadpt_FONLL_5p02[i]);
      pythia_ref->SetBinContent(i+1,dsigmadpt_pythia[i]);
      
      d_r->SetBinError(i+1,e_real[i]);
      d_f->SetBinError(i+1,e_FONLL[i]);
      d_f5->SetBinError(i+1,e_FONLL_5p02[i]);
      pythia_ref->SetBinError(i+1,e_pythia[i]);
    }
  
  TH1F* real_ref = (TH1F*) d_r->Clone("real_ref");
  real_ref->Sumw2();
  real_ref->Multiply(d_f5);
  real_ref->Divide(d_f);

  real_ref->SetTitle("Bpt distribution comparison");
  real_ref->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  real_ref->GetYaxis()->SetTitle("Normalized Distribution");
  real_ref->GetYaxis()->SetTitleOffset(1.0);
  
  float norm;
  norm = real_ref->Integral();
  real_ref->Scale(1.0/norm);
  norm = pythia_ref->Integral();
  pythia_ref->Scale(1.0/norm);

  TFile* inputMC = new TFile(inputmcname.Data());
  TTree* Gen = (TTree*) inputMC->Get("Bfinder/ntGen");
  TH1F* Genpt = new TH1F("Genpt","",nBins,ptBins[0],ptBins[nBins]);
  Gen->Project("Genpt","Gpt","(GisSignal==1)&&(GcollisionId==0)&&((Gpt>5&&Gpt<10&&TMath::Abs(Gy)<2.4)||(Gpt>10&&Gpt<17&&TMath::Abs(Gy)<1.45)||(Gpt>17&&Gpt<100&&TMath::Abs(Gy)<2.1))");
  
  norm = Genpt->Integral();
  Genpt->Scale(1.0/norm);
  
  TH1F* Ratiopt = (TH1F*) real_ref->Clone("Ratiopt");
  Ratiopt->Sumw2();
  Ratiopt->Divide(Genpt);
  
  TF1 *f = new TF1("f","[0]+[1]*x+[2]*x*x+[3]*x*x*x",5,100);
  f->SetParameters(1,0.1,0.1);
  f->SetParLimits(0,0,10);
  f->SetParLimits(1,-10,10);
  f->SetParLimits(2,-10,10);
  f->SetParLimits(3,-10,10);

  Ratiopt->Fit("f");

  printf("Bpt weight function: %f+%f*Bpt+%f*Bpt*Bpt+%f*Bpt*Bpt*Bpt\n",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3));
  
  Ratiopt->SetTitle("Ratio between Real data and Genpt");
  Ratiopt->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  Ratiopt->GetYaxis()->SetTitle("Real/Gen");
  Ratiopt->GetYaxis()->SetTitleOffset(1.0);
  
  TCanvas* c1 = new TCanvas("","",600,600);
  c1->cd();
  real_ref->SetLineColor(kRed);
  pythia_ref->SetLineColor(kGreen);
  Genpt->SetLineColor(kBlue);
  real_ref->Draw();
  pythia_ref->Draw("same");
  Genpt->Draw("same");

  TLegend *leg = new TLegend(0.45,0.75,0.75,0.90,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(real_ref,"real_reference","l");
  leg->AddEntry(pythia_ref,"pythia_reference","l");
  leg->AddEntry(Genpt,"Genpt_MC","l");
  leg->Draw("same");

  c1->SetLogy();
  c1->SaveAs("Bpt_comparison.png");

  TCanvas* c2 = new TCanvas("","",600,600);
  c2->cd();
  Ratiopt->Draw("");
  f->Draw("same");
  
  TLatex* fitftn = new TLatex(0.15,0.8,Form("Bpt weight: %f+%f*Bpt+%f*Bpt^2+%f*Bpt^3",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3)));
  fitftn->SetNDC();
  fitftn->SetTextFont(42);
  fitftn->SetTextSize(0.020);
  fitftn->SetLineWidth(1);
  fitftn->Draw("same");

  c2->SaveAs("Bpt_weight.png");

  TFile* output_ref = new TFile(outputname_ref.Data(),"recreate");
  output_ref->cd();
  real_ref->Write();
  pythia_ref->Write();
  Genpt->Write();
  Ratiopt->Write();
  output_ref->Close();
  
  return;
}

int main()
{
  ppref_real();
  return 0;
}
