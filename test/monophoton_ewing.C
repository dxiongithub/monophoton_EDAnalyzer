#define mvamonopho17_cxx
#include "mvamonopho17.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>


using namespace std;

//if there are only oot || if they don't match
Float_t newuMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newuMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newuMET;
}

//if they are both and they do match
Float_t newmMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t phoEt;
	Float_t phoPhi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t ITEtX = phoEt * cos(phoPhi);
	Float_t ITEtY = phoEt * sin(phoPhi);

	pfMETX += ITEtX;
	pfMETY += ITEtY;

	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newmMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newmMET;
}


//for checking matching

Float_t DeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2)
{

	Float_t ophoPhi = phi1;
	Float_t phoPhi = phi2;
	Float_t dphi = fabs(phoPhi - ophoPhi);
	Float_t tdphi = dphi;
    //phi wrap
	if (dphi > TMath::Pi()) tdphi = TMath::Pi()*2.0 - dphi;
	dphi = tdphi;

	Float_t ophoEta = eta1;
	Float_t phoEta = eta2;
	Float_t deta = fabs(phoEta - ophoEta);

	Float_t deltaR = sqrt(deta*deta + dphi * dphi);
	return deltaR;
}


//calulate invariant mass Z

Float_t InvariMass(Float_t Et1, Float_t Et2, Float_t Phi1, Float_t Phi2, Float_t Eta1, Float_t Eta2)
{
	Float_t Theta1 = 2 * atan(exp(-1.*Eta1));
	Float_t Theta2 = 2 * atan(exp(-1.*Eta2));
	Float_t phoPhi1 = Phi1;
	Float_t phoPhi2 = Phi2;
	Float_t Etot1 = Et1 / sin(Theta1); //E_tot1
	Float_t Etot2 = Et2 / sin(Theta2); //E_tot2
	
	//reconstruct the vectors for x, y and z

	
	Float_t phoX1 = Etot1 * cos(Phi1) * sin(Theta1); 
	Float_t phoY1 = Etot1 * sin(Phi1) * sin(Theta1);
	Float_t phoZ1 = Etot1 * cos(Theta1);

	Float_t phoX2 = Etot2 * cos(Phi2) * sin(Theta2);
	Float_t phoY2 = Etot2 * sin(Phi2) * sin(Theta2);
	Float_t phoZ2 = Etot2 * cos(Theta2);

	Float_t EX1 = Et1*sin(Phi1);
	Float_t EY1 = Et1*cos(Phi1);
	Float_t EZ1 = Etot1*cos(Theta1);

	Float_t EX2 = Et2*sin(Phi2);
	Float_t EY2 = Et2*cos(Phi2);
	Float_t EZ2 = Etot2*cos(Theta2);

	Float_t E1 = sqrt(Etot1*Etot1);
	Float_t E2 = sqrt(Etot2*Etot2);

	//Float_t InvariMassSq = 2 * E1*E2 - phoMag1 * phoMag1 - phoMag2 * phoMag2 - 2 * DotProd12;

	//1	
	Float_t InvariMassSq = (Etot1 + Etot2)*(Etot1 + Etot2) - (EX1 + EX2)*(EX1 + EX2) - (EY1 + EY2)*(EY1 + EY2) - (EZ1 + EZ2)*(EZ1 + EZ2);

	Float_t InvMass = sqrt(InvariMassSq);

	
	
	return InvMass;
	
}


void mvamonopho17::Loop()
{
//   In a ROOT session, you can do:
//      root> .L mvamonopho17.C
//      root> mvamonopho17 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    TFile *monophoseedtime = new TFile("monopho17_ewing.root", "RECREATE");

   
    
    //==================== ewing study =======================
    //candidate ewing > 0.01, in 3ns
    //spike ewing > 0.01, time < -12.5ns
    //prompt ewing > 0.01, in 3ns

    
    TH1F *ewing_candidate_3ns = new TH1F("ewing_candidate_3ns", "ewing_candidate_3ns", 2000, 0, 1.0);
	ewing_candidate_3ns->GetXaxis()->SetTitle("ewing");
	ewing_candidate_3ns->GetYaxis()->SetTitle("Entries");

    TH1F *timing_ewing_candidate_3ns = new TH1F("timing_ewing_candidate_3ns", "timing_ewing_candidate_3ns", 200, -25, 25);
	timing_ewing_candidate_3ns->GetXaxis()->SetTitle("Timing");
	timing_ewing_candidate_3ns->GetYaxis()->SetTitle("Entries"); 

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if ((jentry % 10000) == 0)
            // to print the number of processed entries
            std::cout << "Processed: " << jentry << std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        


        
        //============Selections==============

        /*
        ================================================================
        ================================================================
        ================================================================
        ================================================================
        ================================================================
        */

        //declare variables
        
		Bool_t misOOT = kFALSE; //if it's OOT photon, select the OOT, can be true or false but false makes more sense.
		Int_t mIDXSEL = -1;
        Int_t mZIDXSEL = -1;

        vector<Int_t> cellsEB_right;
        vector<Int_t> cellsEB_left;

        Int_t rightcellsIDX_only = -1;
        Int_t leftcellsIDX_only = -1;

        Int_t rightcellsIDX = -1;
        Int_t leftcellsIDX = -1;

        Float_t maxCellsE = 0;

        Float_t it_initimediff = 0;
        Float_t lictdit = 0;



        //=========================preselection============================

        //adding these tempaltes info for doing timing fit with LICTD info

        vector<Int_t> itpho;


        //in-time collection
		for (int ipho = 0; ipho < nPho; ipho++)
        {
            if ((*phoEt)[ipho] > 225 && fabs((*phoEta)[ipho]) < 1.442 && (*phoR9)[ipho] > 0.8 
            && (*phoHoverE)[ipho] < 0.042 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
            
				itpho.push_back(ipho);
                
			}
		}


        //===========================================================================


        //===========================analysis=====================================


        //---------------------------case 1, only it, no oot---------------------------
		if (itpho.size() != 0)
		{
			mIDXSEL = itpho[0];
        
            //candidate event
			if (pfMET > 210 && (*phohasPixelSeed)[mIDXSEL] == 0 && metFilters == 0 &&
            (*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
            //(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
            //(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
            (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015) //this cut will be replaced by the new Photon IDs (BDT)
			{      
                //searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itpho[0])
                            {
                                //cal aboslute timing diff
                                float it_abstimediff = abs((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (it_abstimediff > it_initimediff)
                                {
                                    it_initimediff = it_abstimediff;
                                    lictdit = ((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    //to see the distribution
                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        ewing_candidate_3ns->Fill(Ewing);
                    }
                    

                    //candidate for ewing & timing fit
                    if (Ewing > 0.01)
                    {
                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            timing_ewing_candidate_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itpho[0])
                            {
                                //cal aboslute timing diff
                                float it_abstimediff = abs((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (it_abstimediff > it_initimediff)
                                {
                                    it_initimediff = it_abstimediff;
                                    lictdit = ((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }

                    
                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        ewing_candidate_3ns->Fill(Ewing);
                    }
                    
                    //candidate for ewing & timing fit
                    if (Ewing > 0.01)
                    {
                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            timing_ewing_candidate_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    } 
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {
                        for (int clusidx = 0; clusidx < 30; clusidx++)
                        {
                            if (AllCellsE_EB[icell] > 1.0 && AllClusteredEB[icell][clusidx] == itpho[0])
                            {
                                //cal aboslute timing diff
                                float it_abstimediff = abs((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                //compare the timing difference (absolute value of the timing diff) with absolute_timingdiff
                                if (it_abstimediff > it_initimediff)
                                {
                                    it_initimediff = it_abstimediff;
                                    lictdit = ((*phoSeedTime)[mIDXSEL] - AllTimeEB[icell]);
                                }              
                            }                   
                        }                
                    }


                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        ewing_candidate_3ns->Fill(Ewing);
                    }
                    

                    //candidate for ewing & timing fit
                    if (Ewing > 0.01)
                    {
                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            timing_ewing_candidate_3ns->Fill((*phoSeedTime)[mIDXSEL]);                            
                        }
                    }
                }
            }
        }
    }
        

    //================= ewing ==============

    ewing_candidate_3ns->Write();
	

    timing_ewing_candidate_3ns->Write();
	
    //==============================test area========================================





    cout << "monopho17_ewing.root" << endl;
}
