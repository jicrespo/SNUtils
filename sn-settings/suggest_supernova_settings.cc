// Code by Nathaniel Tagg to determine thresholds for zero suppression in the SN stream from OM data
// Run it at near1

#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <map>


void suggest_supernova_settings( const char* argv )
{
  
  // int run = 17671; // original thresholds set for late July trials 2018
  // int run = 18468;  // Rerun aug 22 to try a reconfigured asic
  // int run = 21774;  // Rerun March 20, 2019

  int run = atoi( argv );

  double zero_fraction = 0.985; // 0.995
  //  jcrespo: 0.985 has a theoretical compression factor of 67
  bool symmetric = true;
  double pass_low_frac = (1.0-zero_fraction)/2;
  double pass_high_frac = (1.0-zero_fraction)/2;
 
  // read plex
  std::map<std::string,int> plex;
  ifstream in("plex.txt");
  while(in.good()) {
    std::string ccc;
    int w,plane,wirenum;
    in >> ccc >> w >> plane >> wirenum;
    // std::cout << ccc << " " << wirenum << std::endl;
    plex[ccc] = wirenum;
  }
 

  TH1* hccc = new TH1F("hccc",
                       "hccc",
                      10*20*64,
                      0,
                      10*20*64);
  hccc->GetXaxis()->SetTitle("Card address: (card + 20*crate)*64+channel");
  TH1* hwire = new TH1F("hwire",
                        "hwire",
                        8258,0,8258);
  hwire->GetXaxis()->SetTitle("Wire Number");
  TDirectory* dir;
  
  TFile* fout = new TFile(Form("suggested_supernova_thresholds_run%i.root", run),"RECREATE");
  TNtuple* nt = new TNtuple("nt","ntuple","crate:card:channel:wirenum:ped:low:high");
  fout->cd();
  
  dir = fout->mkdir("tpc/");
  dir->cd();
  TH1* hselect = new TH1D("hselect","ADCs selected with this cut",
                    0x2000,-0x1000,0x1000);
  hselect->GetXaxis()->SetTitle("Pedestal-subtracted ADC");
  hselect->GetYaxis()->SetTitle("Samples");
  hselect->SetDirectory(dir);
  TH1* hall = new TH1D("hall","All ADCs",
                    0x2000,-0x1000,0x1000);
  hall->GetXaxis()->SetTitle("Pedestal-subtracted ADC");
  hall->GetYaxis()->SetTitle("Samples");
  hall->SetDirectory(dir);

  
  fout->mkdir("tpc/mapccc/"); 
  dir = fout->GetDirectory("tpc/mapccc/");
  dir->cd();
  TH1* hccc_ped = (TH1*) (hccc->Clone("h_ped") ); hccc_ped ->SetTitle("Pedestal");
  TH1* hccc_high= (TH1*) (hccc->Clone("h_high")); hccc_high->SetTitle("Threshold-High");
  TH1* hccc_low = (TH1*) (hccc->Clone("h_low") ); hccc_low ->SetTitle("Threshold-Low");
  hccc_ped ->SetDirectory(dir);
  hccc_high->SetDirectory(dir);
  hccc_low ->SetDirectory(dir);

  fout->cd();
  fout->mkdir("tpc/mapwire/"); 
  dir = fout->GetDirectory("tpc/mapwire/");
  dir->cd();

  TH1* hwir_ped = (TH1*) (hwire->Clone("h_ped") );  hwir_ped ->SetTitle("Pedestal");
  TH1* hwir_high= (TH1*) (hwire->Clone("h_high"));  hwir_high->SetTitle("Threshold-High");
  TH1* hwir_low = (TH1*) (hwire->Clone("h_low") );  hwir_low ->SetTitle("Threshold-Low");
  hwir_ped ->SetDirectory(dir);
  hwir_high->SetDirectory(dir);
  hwir_low ->SetDirectory(dir);


  ofstream outtxt(Form("suggested_supernova_thresholds_run%i.txt", run));
  ofstream outfcl(Form("suggested_supernova_thresholds_run%i.fcl", run));
  outtxt << "Crate\tCard\tChannel\tpedestal\tthreshold_low\tthreshold_high" << std::endl;
  // TH2* h = new TH2F("h","",100);
  TFile fin(Form("/datalocal/om/run_%08d.om.root",run));
  for(int crate=1;crate<=9;crate++) {
    outfcl << Form("sebAppseb%02d: [ {",crate) << std::endl;
    outfcl << "  zero_suppression:{" << std::endl;
    outfcl << "    params:{ # Non-channel-wise parameters\n\
      channel_threshold : true\n\
      load_threshold_mean : 2\n\
      load_threshold_variance : 3\n\
      postsample : 7\n\
      presample : 7\n\
    }\n";
    for(int card=4;card<=18;card++) {
      
      for(int channel=0;channel<64;channel++) {
        TH1F* hraw = (TH1F*) fin.Get(Form("/tpc/crate%d/card%02d/chan%02d/h_01_raw",crate,card,channel));
        if(hraw) {
          // std::cout << crate << "|" << card << "|" << channel << " " << (hraw==0?"no":"yes") <<  " " <<std::endl;
          int pedestal_bin = hraw->GetMaximumBin();
          double pedestal = hraw->GetBinLowEdge(hraw->GetMaximumBin());
          // Now attempt to find rigorous inclusion rates.
          double tot = hraw->GetEntries();
          // Move in from LHS, adding bins until you get to the pass_low_frac
          double sum;
          int bin;
          double cut_low, cut_high;
          
          if(symmetric) {
            int cut;
            sum = hraw->GetBinContent(pedestal_bin);
            for(cut=1;cut<1000;cut++) {
              sum+=hraw->GetBinContent(pedestal_bin+cut);
              sum+=hraw->GetBinContent(pedestal_bin-cut);
              if(sum > tot*zero_fraction) break;
            }
            cut++;
            cut_high = cut;
            cut_low  = -cut;
            std::cout << sum/tot << "\t";            
          } else {
            // Asymmetric 
            sum = 0;
            for(bin=0;bin<pedestal_bin;bin++) {
              sum+=hraw->GetBinContent(bin);
              if(sum > tot*pass_low_frac) break;
            }
            cut_low = bin - pedestal_bin;

            sum = 0;
            for(bin=hraw->GetNbinsX()+1;bin>pedestal_bin;bin--) {
              sum+=hraw->GetBinContent(bin);
              if(sum > tot*pass_high_frac) break;
            }
            cut_high = bin - pedestal_bin;
          }
          
          for(bin=1;bin<hraw->GetNbinsX();bin++) 
          {
            double adc = hraw->GetBinLowEdge(bin);
            double x = adc - pedestal;
            hall->Fill(x,hraw->GetBinContent(bin));
            if(x<=cut_low || x>=cut_high) hselect->Fill(x,hraw->GetBinContent(bin));
          }
	  
          std::string ccc = Form("%1d|%02d|%02d",crate,card,channel);
          int wirenum = plex[ccc];
          
          int bin_ccc = channel + 64*(card + 20*crate)+1;
          int bin_wire= wirenum+1;
          hccc_ped ->SetBinContent(bin_ccc,pedestal);                     hccc_ped ->SetBinError(bin_ccc,1);
          hccc_high->SetBinContent(bin_ccc,cut_high);            hccc_high->SetBinError(bin_ccc,1);
          hccc_low ->SetBinContent(bin_ccc,-cut_low);             hccc_low ->SetBinError(bin_ccc,1);
          hwir_ped ->SetBinContent(bin_wire,pedestal);                     hwir_ped ->SetBinError(bin_wire,1);
          hwir_high->SetBinContent(bin_wire,cut_high);            hwir_high->SetBinError(bin_wire,1);
          hwir_low ->SetBinContent(bin_wire,-cut_low);             hwir_low ->SetBinError(bin_wire,1);
                                                                                    
          std::cout << Form("%1d|%02d|%02d %4.0lf %3.0lf %3.0lf",crate,card,channel,pedestal,cut_low,cut_high) << std::endl;
          outtxt    << Form("%1d\t%02d\t%02d\t%4.0lf\t%3.0lf\t%3.0lf",crate,card,channel,pedestal,cut_low,cut_high) << std::endl;
          nt->Fill(crate,card,channel,wirenum,pedestal,cut_low,cut_high);
          if(channel==0)
            outfcl <<  Form("    slot%d:{",card) << std::endl;
            
          outfcl << Form("      ch%d   : %d\n",channel,(int)cut_high);
          outfcl << Form("      ped%d  : %d\n",channel,(int)pedestal);
          outfcl << Form("      pol%d  : %d\n",channel,3); // 3= bipolar   
          if(channel==63)
            outfcl <<  Form("    }",card) << std::endl;

	  if( (wirenum%1000) == 0 ){           
	    TCanvas c(Form("c_hzs_wire%i", wirenum),"");
	    hraw->Draw();
	    TH1F* hzs = (TH1F*)hraw->Clone("hzs");
	    hzs->Reset();
	    for(int ibin = pedestal_bin+cut_low; ibin <= pedestal_bin+cut_high; ibin++){
	      hzs->SetBinContent(ibin, hraw->GetBinContent(ibin));
	    }
	    hzs->SetFillStyle(1);
	    hzs->SetFillColor(kRed);
	    hzs->Draw("same");
	    c.Print(".root");
	  }
        }
      }

    } // card
    
    outfcl << "  }" << std::endl;
    outfcl << "}  ]" << std::endl;
    outfcl << "\n\n";
  } //crate
  std::cout << "Selected data target: "<< 1.0-zero_fraction << " achieved: " << hselect->Integral()/hall->Integral() << std::endl;
  fout->cd();
  fout->Write(0,TObject::kOverwrite);
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){
  suggest_supernova_settings( argv[1] );
  return 0;
}
# endif
