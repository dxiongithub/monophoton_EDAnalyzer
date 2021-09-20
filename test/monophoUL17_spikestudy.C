#define mvamonopho17_mini2_cxx
#include "mvamonopho17_mini2.h"
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



//#include "extra_tools.cc"

//https://github.com/wadud001/UL17PhoEBID


#include "helpers.cc"
#include "/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/xgboost/1.3.3/include/xgboost/c_api.h"
R__ADD_INCLUDE_PATH(/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/xgboost/1.3.3/include/xgboost)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/xgboost/1.3.3/lib64/libxgboost.so)


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


void mvamonopho17_mini2::Loop()
{
//   In a ROOT session, you can do:
//      root> .L mvamonopho17_mini2.C
//      root> mvamonopho17_mini2 t
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
    TFile *monophoseedtime = new TFile("monophoUL17_spikesB.root", "RECREATE");


    //============ BDT setup ==================
    //variables
    /// BDT
    DMatrixHandle 		dTest;
    BoosterHandle 		phoBDT_h;
    XGBoosterCreate(NULL, 0, &phoBDT_h); 
    XGBoosterSetParam(phoBDT_h, "seed", "0"); 
    Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, "aNTGC_photon_BDT_EB_2021_08_26_09_39_52.model");

    //new IDs - Isolation corrections

    isoCorrMap ecalIsoRhoCorrMap("phoPFClusEcalIso_PtCorrections.txt", 2);
    isoCorrMap tkrIsoRhoCorrMap("phoTrkSumPtHollowConeDR03_RhoCorrections.txt", 2);
    isoCorrMap ecalIsoPtCorrMap("phoPFClusEcalIso_PtCorrections.txt", 2);

   

    TH1F *timing_beamhalo_iso_BDT_ewing1 = new TH1F("timing_beamhalo_iso_BDT_ewing1", "timing_beamhalo_iso_BDT_ewing1", 200, -25, 25);
	timing_beamhalo_iso_BDT_ewing1->GetXaxis()->SetTitle("Timing");
	timing_beamhalo_iso_BDT_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_candidate_iso_BDT_ewing1 = new TH1F("timing_candidate_iso_BDT_ewing1", "timing_candidate_iso_BDT_ewing1", 200, -25, 25);
	timing_candidate_iso_BDT_ewing1->GetXaxis()->SetTitle("timing");
	timing_candidate_iso_BDT_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_spike_iso_BDT_ewing1 = new TH1F("timing_spike_iso_BDT_ewing1", "timing_spike_iso_BDT_ewing1", 200, -25, 25);
	timing_spike_iso_BDT_ewing1->GetXaxis()->SetTitle("timing");
	timing_spike_iso_BDT_ewing1->GetYaxis()->SetTitle("Entries");

    TH1F *timing_promptZ_iso_BDT_ewing1 = new TH1F("timing_promptZ_iso_BDT_ewing1", "timing_promptZ_iso_BDT_ewing1", 200, -25, 25);
	timing_promptZ_iso_BDT_ewing1->GetXaxis()->SetTitle("timing");
	timing_promptZ_iso_BDT_ewing1->GetYaxis()->SetTitle("Entries");
    

    TH1F *candidate_BDTscore = new TH1F("candidate_BDTscore", "candidate_BDTscore", 200, 0, 1);
	candidate_BDTscore->GetXaxis()->SetTitle("BDT score");
	candidate_BDTscore->GetYaxis()->SetTitle("Entries");

    TH1F *spike_BDTscore = new TH1F("spike_BDTscore", "spike_BDTscore", 200, 0, 1);
	spike_BDTscore->GetXaxis()->SetTitle("BDT score");
	spike_BDTscore->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_BDTscore = new TH1F("promptZ_BDTscore", "promptZ_BDTscore", 200, 0, 1);
	promptZ_BDTscore->GetXaxis()->SetTitle("BDT score");
	promptZ_BDTscore->GetYaxis()->SetTitle("Entries");

    TH1F *beamhalo_BDTscore = new TH1F("beamhalo_BDTscore", "beamhalo_BDTscore", 200, 0, 1);
	beamhalo_BDTscore->GetXaxis()->SetTitle("BDT score");
	beamhalo_BDTscore->GetYaxis()->SetTitle("Entries");


    TH1F *candidate_ecaliso = new TH1F("candidate_ecaliso", "candidate_ecaliso", 2000, 0, 100);
	candidate_ecaliso->GetXaxis()->SetTitle("ecaliso");
	candidate_ecaliso->GetYaxis()->SetTitle("Entries");

    TH1F *spike_ecaliso = new TH1F("spike_ecaliso", "spike_ecaliso", 2000, 0, 100);
	spike_ecaliso->GetXaxis()->SetTitle("ecaliso");
	spike_ecaliso->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_ecaliso = new TH1F("promptZ_ecaliso", "promptZ_ecaliso", 2000, 0, 100);
	promptZ_ecaliso->GetXaxis()->SetTitle("ecaliso");
	promptZ_ecaliso->GetYaxis()->SetTitle("Entries");

    TH1F *beamhalo_ecaliso = new TH1F("beamhalo_ecaliso", "beamhalo_ecaliso", 2000, 0, 100);
	beamhalo_ecaliso->GetXaxis()->SetTitle("ecaliso");
	beamhalo_ecaliso->GetYaxis()->SetTitle("Entries");


    TH1F *candidate_trackiso = new TH1F("candidate_trackiso", "candidate_trackiso", 2000, 0, 100);
	candidate_trackiso->GetXaxis()->SetTitle("trackiso");
	candidate_trackiso->GetYaxis()->SetTitle("Entries");

    TH1F *spike_trackiso = new TH1F("spike_trackiso", "spike_trackiso", 2000, 0, 100);
	spike_trackiso->GetXaxis()->SetTitle("trackiso");
	spike_trackiso->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_trackiso = new TH1F("promptZ_trackiso", "promptZ_trackiso", 2000, 0, 100);
	promptZ_trackiso->GetXaxis()->SetTitle("trackiso");
	promptZ_trackiso->GetYaxis()->SetTitle("Entries");

    TH1F *beamhalo_trackiso = new TH1F("beamhalo_trackiso", "beamhalo_trackiso", 2000, 0, 100);
	beamhalo_trackiso->GetXaxis()->SetTitle("trackiso");
	beamhalo_trackiso->GetYaxis()->SetTitle("Entries");


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

        

        //declare variables
        Int_t nphotons = 0;
        Int_t nophotons = 0;
        
        Int_t mitIDXSEL = -1; //set index for photon, initial value could be something non-physical, negative something.
		Int_t mootIDXSEL = -1;
        Int_t mootZIDXSEL = -1;
        Int_t mitZIDXSEL = -1;
        
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



        //=========================preselection============================

        vector<Int_t> itpho;
		vector<Int_t> itootmatched;
		vector<Int_t> ootpho;
		vector<Int_t> ootitmatched;

        vector<Int_t> iiZpho;
		vector<Int_t> itootZmatched;
		vector<Int_t> ooZpho;
		vector<Int_t> ootitZmatched;

        //in-time collection
		for (int ipho = 0; ipho < nPho; ipho++)
        {
            if ((*phoEt)[ipho] > 225 && fabs((*phoEta)[ipho]) < 1.442 && (*phoR9)[ipho] > 0.8 
            && (*phoHoverE)[ipho] < 0.0346 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
				itpho.push_back(ipho);
				itootmatched.push_back(-1);

                //prompt Z in-time collection
                float PhoSep = 999;
				for (int iipho = ipho + 1; iipho < nPho; iipho++)
				{
					Float_t deltaR = DeltaR((*phoEta)[ipho], (*phoPhi)[ipho], (*phoEta)[iipho], (*phoPhi)[iipho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}
                
                    //collect the second photon
                    if ((*phoEt)[iipho] > 10.0 && fabs((*phoEta)[iipho]) < 1.442 && (*phoR9)[iipho] > 0.8 && (*phoHoverE)[iipho] < 0.0346 && PhoSep > 0.2)
                    {
                        iiZpho.push_back(iipho);
                        itootZmatched.push_back(-1);
                    }
                }
                
                nphotons++;
			}
		}

        //out-of-time collection
		for (int opho = 0; opho < onPho; opho++)
		{
            if ((*ophoEt)[opho] > 225 && fabs((*ophoEta)[opho]) < 1.442 
            && (*ophoR9)[opho] > 0.8 && (*ophoHoverE)[opho] < 0.0346 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
            
				ootpho.push_back(opho);
				ootitmatched.push_back(-1);


                //prompt Z out-of-time collection
                float PhoSep = 999;
				for (int oopho = opho + 1; oopho < nPho; oopho++)
				{
					Float_t deltaR = DeltaR((*ophoEta)[opho], (*ophoPhi)[opho], (*ophoEta)[oopho], (*ophoPhi)[oopho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}

                    //collect the second photon
                    if ((*ophoEt)[oopho] > 10.0 && fabs((*ophoEta)[oopho]) < 1.442 && (*ophoR9)[oopho] > 0.8 && (*ophoHoverE)[oopho] < 0.0346)
                    {
                        ooZpho.push_back(oopho);
                        ootitZmatched.push_back(-1);
                    }
                }

                nophotons++;

			}
		}

        //for getting the minimum deltaR
		//for (unsigned i = 0; i < itpho.size(); i++)
		for (unsigned i = 0; i < itpho.size(); ++i)
		{
			float PhoSep = 999;
			Int_t itootmatchedidx = -1; //defualt that they are not matched
			//Int_t itootseparateidx = 1;

			//for (unsigned j = 0; j < ootpho.size(); j++)
			for (unsigned j = 0; j < ootpho.size(); ++j)
			{
				Float_t deltaR = DeltaR((*phoEta)[i], (*phoPhi)[i], (*ophoEta)[j], (*ophoPhi)[j]);
				if (deltaR < PhoSep)
				{
					PhoSep = deltaR;
				}
				
				if (PhoSep < 0.2)
				{
					itootmatchedidx = j;
					itootmatched[i] = itootmatchedidx;
					ootitmatched[itootmatchedidx] = i;
				}
			}
		}

        //second photon matched for Z prompt
        for (unsigned i = 0; i < iiZpho.size(); ++i)
		{
			float PhoSep = 999;
			Int_t itootZmatchedidx = -1; //defualt that they are not matched
			//Int_t itootseparateidx = 1;

			//for (unsigned j = 0; j < ootpho.size(); j++)
			for (unsigned j = 0; j < ooZpho.size(); ++j)
			{
				Float_t deltaR = DeltaR((*phoEta)[i], (*phoPhi)[i], (*ophoEta)[j], (*ophoPhi)[j]);
				if (deltaR < PhoSep)
				{
					PhoSep = deltaR;
				}
				
				if (PhoSep < 0.2)
				{
					itootZmatchedidx = j;
					itootZmatched[i] = itootZmatchedidx;
					ootitZmatched[itootZmatchedidx] = i;
				}
			}
		}

        //===========================================================================


        //===========================analysis=====================================


        //---------------------------case 1, only it, no oot---------------------------
		if (itpho.size() != 0 && ootpho.size() == 0 && itootmatched.size() != 0)
		{
			mIDXSEL = itpho[0];
        
            //candidate event
			if (pfMET > 210 && (*phohasPixelSeed)[mIDXSEL] == 0 && metFilters == 0 &&
            (*phoMIPTotEnergy)[mIDXSEL] < 4.9)
			{      
                vector<Float_t> feats
                {
                    ((*phoE2x2Full5x5)[mIDXSEL]/((*phoR9Full5x5)[mIDXSEL]*(*phoSCRawE)[mIDXSEL])), 
                    ((*phoE1x3Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoE2ndFull5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoE2x5Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]),
                    ((*phoMaxEnergyXtal)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoSCEtaWidth)[mIDXSEL]/(*phoSCPhiWidth)[mIDXSEL]), 
                    (*phoSCEtaWidth)[mIDXSEL], 
                    (*phoSCPhiWidth)[mIDXSEL], 
                    (*phoR9Full5x5)[mIDXSEL], 
                    (*phoCalibEt)[mIDXSEL],
                    ((*phoE2x2Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]/(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL]), 
                    (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], 
                    (*phoSigmaIEtaIPhiFull5x5)[mIDXSEL], 
                    (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL]
                };

                XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                bst_ulong out_len;
                const float *prediction;
                XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                assert(out_len == 1);
                XGDMatrixFree(dTest);
                Float_t iPhoBDTpred = prediction[0];
                
                Float_t phoPFECALClusIsoCorr      = (*phoPFClusEcalIso)[mIDXSEL] - ecalIsoRhoCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), rho, 0) - ecalIsoPtCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), (*phoEt)[mIDXSEL], 1);
 		        Float_t phoTkrIsoCorr             = (*phoTrkSumPtHollowConeDR03)[mIDXSEL] - tkrIsoRhoCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), rho, 1);

                //cand1 searching adjencent cells
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

                    if (nphotons == 1) //on candidate only
                    {
                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                            {
                                candidate_BDTscore->Fill(iPhoBDTpred);
                            }

                            if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                            {
                                candidate_ecaliso->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                            {
                                candidate_trackiso->Fill(phoTkrIsoCorr);
                            }

                            if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                            {
                               timing_candidate_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //candidate timing fit
                            }
                        }
                    }

                    //spike
                    if (Ewing < 0.01)
                    {
                        if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            spike_BDTscore->Fill(iPhoBDTpred);
                        }

                        if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                        {
                            spike_ecaliso->Fill(phoPFECALClusIsoCorr);
                        }

                        if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                        {
                            spike_trackiso->Fill(phoTkrIsoCorr);
                        }

                        if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            timing_spike_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //spike timing fit
                        }
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    if (nphotons == 1) //on candidate only
                    {
                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                            {
                                candidate_BDTscore->Fill(iPhoBDTpred);
                            }

                            if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                            {
                                candidate_ecaliso->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                            {
                                candidate_trackiso->Fill(phoTkrIsoCorr);
                            }

                            if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                            {
                               timing_candidate_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //candidate timing fit
                            }
                        }
                    }

                    //spike
                    if (Ewing < 0.01)
                    {
                        if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            spike_BDTscore->Fill(iPhoBDTpred);
                        }

                        if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                        {
                            spike_ecaliso->Fill(phoPFECALClusIsoCorr);
                        }

                        if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                        {
                            spike_trackiso->Fill(phoTkrIsoCorr);
                        }

                        if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            timing_spike_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //spike timing fit
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

            
                    if (nphotons == 1) //on candidate only
                    {
                        if (Ewing > 0.01)
                        {
                            if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                            {
                                candidate_BDTscore->Fill(iPhoBDTpred);
                            }

                            if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                            {
                                candidate_ecaliso->Fill(phoPFECALClusIsoCorr);
                            }

                            if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                            {
                                candidate_trackiso->Fill(phoTkrIsoCorr);
                            }

                            if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                            {
                               timing_candidate_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //candidate timing fit
                            }
                        }
                    }

                    //spike
                    if (Ewing < 0.01)
                    {
                        if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            spike_BDTscore->Fill(iPhoBDTpred);
                        }

                        if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                        {
                            spike_ecaliso->Fill(phoPFECALClusIsoCorr);
                        }

                        if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                        {
                            spike_trackiso->Fill(phoTkrIsoCorr);
                        }

                        if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            timing_spike_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //spike timing fit
                        }
                    }
                }
            }
        

    
            //promptZ selection
            if (iiZpho.size() != 0 && ooZpho.size() == 0 && itootZmatched.size() != 0)
            {
                mZIDXSEL = iiZpho[0];
                //Prompt Z template:
                if ((*phohasPixelSeed)[mIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*phohasPixelSeed)[mZIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mZIDXSEL] < 4.9)
                {
                    vector<Float_t> feats
                    {
                        ((*phoE2x2Full5x5)[mIDXSEL]/((*phoR9Full5x5)[mIDXSEL]*(*phoSCRawE)[mIDXSEL])), 
                        ((*phoE1x3Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                        ((*phoE2ndFull5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                        ((*phoE2x5Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]),
                        ((*phoMaxEnergyXtal)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                        ((*phoSCEtaWidth)[mIDXSEL]/(*phoSCPhiWidth)[mIDXSEL]), 
                        (*phoSCEtaWidth)[mIDXSEL], 
                        (*phoSCPhiWidth)[mIDXSEL], 
                        (*phoR9Full5x5)[mIDXSEL], 
                        (*phoCalibEt)[mIDXSEL],
                        ((*phoE2x2Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                        ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]/(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL]), 
                        (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], 
                        (*phoSigmaIEtaIPhiFull5x5)[mIDXSEL], 
                        (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL]
                    };

                    XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                    bst_ulong out_len;
                    const float *prediction;
                    XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                    assert(out_len == 1);
                    XGDMatrixFree(dTest);
                    Float_t iPhoBDTpred = prediction[0];
                    
                    Float_t phoPFECALClusIsoCorr      = (*phoPFClusEcalIso)[mIDXSEL] - ecalIsoRhoCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), rho, 0) - ecalIsoPtCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), (*phoEt)[mIDXSEL], 1);
                    Float_t phoTkrIsoCorr             = (*phoTrkSumPtHollowConeDR03)[mIDXSEL] - tkrIsoRhoCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), rho, 1);

                    Float_t InvM = InvariMass((*phoEt)[mIDXSEL], (*phoEt)[mZIDXSEL], (*phoPhi)[mIDXSEL], (*phoPhi)[mZIDXSEL], (*phoEta)[mIDXSEL], (*phoEta)[mZIDXSEL]);
                    

                    if (InvM > 85 && InvM < 100)
                    {
                        //promptZ1 searching adjencent cells
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
                
                            if (Ewing > 0.01 && fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                                {
                                    promptZ_BDTscore->Fill(iPhoBDTpred);
                                }

                                if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                                {
                                    promptZ_ecaliso->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                                {
                                    promptZ_trackiso->Fill(phoTkrIsoCorr);
                                }

                                if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                                {
                                    timing_promptZ_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //promptZ timing
                                }
                            }
                        }
                        

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                            if (Ewing > 0.01 && fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                                {
                                    promptZ_BDTscore->Fill(iPhoBDTpred);
                                }

                                if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                                {
                                    promptZ_ecaliso->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                                {
                                    promptZ_trackiso->Fill(phoTkrIsoCorr);
                                }

                                if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                                {
                                    timing_promptZ_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //promptZ timing
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


                           if (Ewing > 0.01 && fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                                {
                                    promptZ_BDTscore->Fill(iPhoBDTpred);
                                }

                                if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                                {
                                    promptZ_ecaliso->Fill(phoPFECALClusIsoCorr);
                                }

                                if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                                {
                                    promptZ_trackiso->Fill(phoTkrIsoCorr);
                                }

                                if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                                {
                                    timing_promptZ_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]); //promptZ timing
                                }
                            }
                        }
                    }
                }
            }

            //BeamHalo template
            if (pfMET > 210 &&
                (*phohasPixelSeed)[mIDXSEL] == 0 &&
                (*phoMIPTotEnergy)[mIDXSEL] > 4.9 &&
                metFilters == 8)//106X beamhalo flag
            {
                vector<Float_t> feats
                {
                    ((*phoE2x2Full5x5)[mIDXSEL]/((*phoR9Full5x5)[mIDXSEL]*(*phoSCRawE)[mIDXSEL])), 
                    ((*phoE1x3Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoE2ndFull5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoE2x5Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]),
                    ((*phoMaxEnergyXtal)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoSCEtaWidth)[mIDXSEL]/(*phoSCPhiWidth)[mIDXSEL]), 
                    (*phoSCEtaWidth)[mIDXSEL], 
                    (*phoSCPhiWidth)[mIDXSEL], 
                    (*phoR9Full5x5)[mIDXSEL], 
                    (*phoCalibEt)[mIDXSEL],
                    ((*phoE2x2Full5x5)[mIDXSEL]/(*phoSCRawE)[mIDXSEL]), 
                    ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]/(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL]), 
                    (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL], 
                    (*phoSigmaIEtaIPhiFull5x5)[mIDXSEL], 
                    (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL]
                };

                XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                bst_ulong out_len;
                const float *prediction;
                XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                assert(out_len == 1);
                XGDMatrixFree(dTest);
                Float_t iPhoBDTpred = prediction[0];
                
                Float_t phoPFECALClusIsoCorr      = (*phoPFClusEcalIso)[mIDXSEL] - ecalIsoRhoCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), rho, 0) - ecalIsoPtCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), (*phoEt)[mIDXSEL], 1);
 		        Float_t phoTkrIsoCorr             = (*phoTrkSumPtHollowConeDR03)[mIDXSEL] - tkrIsoRhoCorrMap.getIsoCorr(fabs((*phoEta)[mIDXSEL]), rho, 1);

    
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

                    if (Ewing > 0.01)
                    {
                        if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            beamhalo_BDTscore->Fill(iPhoBDTpred);
                        }

                        if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                        {
                            beamhalo_ecaliso->Fill(phoPFECALClusIsoCorr);
                        }

                        if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                        {
                            beamhalo_trackiso->Fill(phoTkrIsoCorr);
                        }

                        if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            timing_beamhalo_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    if (Ewing > 0.01)
                    {
                        if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            beamhalo_BDTscore->Fill(iPhoBDTpred);
                        }

                        if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                        {
                            beamhalo_ecaliso->Fill(phoPFECALClusIsoCorr);
                        }

                        if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                        {
                            beamhalo_trackiso->Fill(phoTkrIsoCorr);
                        }

                        if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            timing_beamhalo_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
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

                    if (Ewing > 0.01)
                    {
                        if (phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            beamhalo_BDTscore->Fill(iPhoBDTpred);
                        }

                        if (iPhoBDTpred > 0.92 && phoTkrIsoCorr < 4.43)
                        {
                            beamhalo_ecaliso->Fill(phoPFECALClusIsoCorr);
                        }

                        if (phoPFECALClusIsoCorr < 6.65 && iPhoBDTpred > 0.92)
                        {
                            beamhalo_trackiso->Fill(phoTkrIsoCorr);
                        }

                        if (iPhoBDTpred > 0.92 && phoPFECALClusIsoCorr < 6.65 && phoTkrIsoCorr < 4.43)
                        {
                            timing_beamhalo_iso_BDT_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }
                }
            }
        }
            

        //---------------------------case 2, only oot, no it---------------------------
		if (itpho.size() == 0 && ootpho.size() != 0 && ootitmatched.size() != 0)
		{
			mIDXSEL = ootpho[0];
			misOOT = kTRUE;
			Float_t NewMET = newuMET(pfMET, pfMETPhi, (*ophoPhi)[mIDXSEL], (*ophoEt)[mIDXSEL]);

			//Candidate Events
			if (NewMET > 210 && (*ophohasPixelSeed)[mIDXSEL] == 0 && metFilters == 0 &&
            (*ophoMIPTotEnergy)[mIDXSEL] < 4.9)
			{
                vector<Float_t> feats
                {
                    ((*ophoE2x2Full5x5)[mIDXSEL]/((*ophoR9Full5x5)[mIDXSEL]*(*ophoSCRawE)[mIDXSEL])), 
                    ((*ophoE1x3Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoE2ndFull5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoE2x5Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]),
                    ((*ophoMaxEnergyXtal)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoSCEtaWidth)[mIDXSEL]/(*ophoSCPhiWidth)[mIDXSEL]), 
                    (*ophoSCEtaWidth)[mIDXSEL], 
                    (*ophoSCPhiWidth)[mIDXSEL], 
                    (*ophoR9Full5x5)[mIDXSEL], 
                    (*ophoEt)[mIDXSEL],
                    ((*ophoE2x2Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]/(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL]), 
                    (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], 
                    (*ophoSigmaIEtaIPhiFull5x5)[mIDXSEL], 
                    (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL]
                };

                XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                bst_ulong out_len;
                const float *prediction;
                XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                assert(out_len == 1);
                XGDMatrixFree(dTest);
                Float_t oPhoBDTpred = prediction[0];
                
                Float_t ophoPFECALClusIsoCorr      = (*ophoPFClusEcalIso)[mIDXSEL] - ecalIsoRhoCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), rho, 0) - ecalIsoPtCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), (*ophoEt)[mIDXSEL], 1);
 		        Float_t ophoTkrIsoCorr             = (*ophoTrkSumPtHollowConeDR03)[mIDXSEL] - tkrIsoRhoCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), rho, 1);

                        
                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
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

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (Ewing > 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            candidate_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            candidate_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            candidate_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                           timing_candidate_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //candidate timing fit
                        }
                    }

                    //spike
                    if (Ewing < 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            spike_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            spike_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            spike_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            timing_spike_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //spike timing fit
                        }
                    }
                }

                            
                            
                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (Ewing > 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            candidate_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            candidate_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            candidate_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                           timing_candidate_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //candidate timing fit
                        }
                    }

                    //spike
                    if (Ewing < 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            spike_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            spike_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            spike_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            timing_spike_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //spike timing fit
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

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

            
                    if (Ewing > 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            candidate_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            candidate_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            candidate_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                           timing_candidate_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //candidate timing fit
                        }
                    }

                    //spike
                    if (Ewing < 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            spike_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            spike_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            spike_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            timing_spike_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //spike timing fit
                        }
                    }
                }
            }
        
            //promptZ
            if (iiZpho.size() == 0 && ooZpho.size() != 0 && ootitZmatched.size() != 0)
            {
                mZIDXSEL = ooZpho[0];
                //Prompt Z template:
                if ((*ophohasPixelSeed)[mIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                    metFilters == 0 && //106X 256 for passing everything&&
                    (*ophohasPixelSeed)[mZIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mZIDXSEL] < 4.9)
                {
                    vector<Float_t> feats
                    {
                        ((*ophoE2x2Full5x5)[mIDXSEL]/((*ophoR9Full5x5)[mIDXSEL]*(*ophoSCRawE)[mIDXSEL])), 
                        ((*ophoE1x3Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                        ((*ophoE2ndFull5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                        ((*ophoE2x5Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]),
                        ((*ophoMaxEnergyXtal)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                        ((*ophoSCEtaWidth)[mIDXSEL]/(*ophoSCPhiWidth)[mIDXSEL]), 
                        (*ophoSCEtaWidth)[mIDXSEL], 
                        (*ophoSCPhiWidth)[mIDXSEL], 
                        (*ophoR9Full5x5)[mIDXSEL], 
                        (*ophoEt)[mIDXSEL],
                        ((*ophoE2x2Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                        ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]/(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL]), 
                        (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], 
                        (*ophoSigmaIEtaIPhiFull5x5)[mIDXSEL], 
                        (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL]
                    };

                    XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                    bst_ulong out_len;
                    const float *prediction;
                    XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                    assert(out_len == 1);
                    XGDMatrixFree(dTest);
                    Float_t oPhoBDTpred = prediction[0];
                    
                    Float_t ophoPFECALClusIsoCorr      = (*ophoPFClusEcalIso)[mIDXSEL] - ecalIsoRhoCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), rho, 0) - ecalIsoPtCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), (*ophoEt)[mIDXSEL], 1);
                    Float_t ophoTkrIsoCorr             = (*ophoTrkSumPtHollowConeDR03)[mIDXSEL] - tkrIsoRhoCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), rho, 1);

                    Float_t InvM = InvariMass((*ophoEt)[mIDXSEL], (*ophoEt)[mZIDXSEL], (*ophoPhi)[mIDXSEL], (*ophoPhi)[mZIDXSEL], (*ophoEta)[mIDXSEL], (*ophoEta)[mZIDXSEL]);

                    //promptZ for timing fit
                    if (InvM > 85 && InvM < 100)
                    {
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {       
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 
                            }

                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                            {               
                                cellsEB_left.push_back(icell);
                            }       

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
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

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                            if (Ewing > 0.01 && fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                                {
                                    promptZ_BDTscore->Fill(oPhoBDTpred);
                                }

                                if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                                {
                                    promptZ_ecaliso->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                                {
                                    promptZ_trackiso->Fill(ophoTkrIsoCorr);
                                }

                                if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                                {
                                    timing_promptZ_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //promptZ timing
                                }
                            }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        
                            if (Ewing > 0.01 && fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                                {
                                    promptZ_BDTscore->Fill(oPhoBDTpred);
                                }

                                if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                                {
                                    promptZ_ecaliso->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                                {
                                    promptZ_trackiso->Fill(ophoTkrIsoCorr);
                                }

                                if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                                {
                                    timing_promptZ_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //promptZ timing
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

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);


                            
                            if (Ewing > 0.01 && fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                                {
                                    promptZ_BDTscore->Fill(oPhoBDTpred);
                                }

                                if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                                {
                                    promptZ_ecaliso->Fill(ophoPFECALClusIsoCorr);
                                }

                                if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                                {
                                    promptZ_trackiso->Fill(ophoTkrIsoCorr);
                                }

                                if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                                {
                                    timing_promptZ_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]); //promptZ timing
                                }
                            }
                        }
                    }
                }
            }

            //BeamHalo template
            if (NewMET > 210 &&
                (*ophohasPixelSeed)[mIDXSEL] == 0 &&
                (*ophoMIPTotEnergy)[mIDXSEL] > 4.9 &&
                metFilters == 8)//106X beamhalo flag
            {

                vector<Float_t> feats
                {
                    ((*ophoE2x2Full5x5)[mIDXSEL]/((*ophoR9Full5x5)[mIDXSEL]*(*ophoSCRawE)[mIDXSEL])), 
                    ((*ophoE1x3Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoE2ndFull5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoE2x5Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]),
                    ((*ophoMaxEnergyXtal)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoSCEtaWidth)[mIDXSEL]/(*ophoSCPhiWidth)[mIDXSEL]), 
                    (*ophoSCEtaWidth)[mIDXSEL], 
                    (*ophoSCPhiWidth)[mIDXSEL], 
                    (*ophoR9Full5x5)[mIDXSEL], 
                    (*ophoEt)[mIDXSEL],
                    ((*ophoE2x2Full5x5)[mIDXSEL]/(*ophoSCRawE)[mIDXSEL]), 
                    ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]/(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL]), 
                    (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL], 
                    (*ophoSigmaIEtaIPhiFull5x5)[mIDXSEL], 
                    (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL]
                };

                XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                bst_ulong out_len;
                const float *prediction;
                XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                assert(out_len == 1);
                XGDMatrixFree(dTest);
                Float_t oPhoBDTpred = prediction[0];
                
                Float_t ophoPFECALClusIsoCorr      = (*ophoPFClusEcalIso)[mIDXSEL] - ecalIsoRhoCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), rho, 0) - ecalIsoPtCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), (*ophoEt)[mIDXSEL], 1);
 		        Float_t ophoTkrIsoCorr             = (*ophoTrkSumPtHollowConeDR03)[mIDXSEL] - tkrIsoRhoCorrMap.getIsoCorr(fabs((*ophoEta)[mIDXSEL]), rho, 1);


                for (int icell = 0; icell < nAllCellsEB; icell++)
                {       
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                    {               
                        cellsEB_left.push_back(icell);
                    }       

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
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

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);
                
                    if (Ewing > 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            beamhalo_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            beamhalo_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            beamhalo_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            timing_beamhalo_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (Ewing > 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            beamhalo_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            beamhalo_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            beamhalo_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            timing_beamhalo_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
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

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    if (Ewing > 0.01)
                    {
                        if (ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            beamhalo_BDTscore->Fill(oPhoBDTpred);
                        }

                        if (oPhoBDTpred > 0.92 && ophoTkrIsoCorr < 4.43)
                        {
                            beamhalo_ecaliso->Fill(ophoPFECALClusIsoCorr);
                        }

                        if (ophoPFECALClusIsoCorr < 6.65 && oPhoBDTpred > 0.92)
                        {
                            beamhalo_trackiso->Fill(ophoTkrIsoCorr);
                        }

                        if (oPhoBDTpred > 0.92 && ophoPFECALClusIsoCorr < 6.65 && ophoTkrIsoCorr < 4.43)
                        {
                            timing_beamhalo_iso_BDT_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                }
            }
        }
    }
        


    /*
    =================================================================================
    =================================================================================
    =================================================================================
    =================================================================================
    =================================================================================
    */


    //=================== timing BDT ========================
   


    timing_beamhalo_iso_BDT_ewing1->Write();

    timing_candidate_iso_BDT_ewing1->Write();

    timing_spike_iso_BDT_ewing1->Write();

    timing_promptZ_iso_BDT_ewing1->Write();
 

    //==============================test area========================================


    candidate_BDTscore->Write();
    candidate_ecaliso->Write();
    candidate_trackiso->Write();

    spike_BDTscore->Write();
    spike_ecaliso->Write();
    spike_trackiso->Write();

    promptZ_BDTscore->Write();
    promptZ_ecaliso->Write();
    promptZ_trackiso->Write();

    beamhalo_BDTscore->Write();
    beamhalo_ecaliso->Write();
    beamhalo_trackiso->Write();






    /*
    =================================================================================
    =================================================================================
    =================================================================================
    =================================================================================
    =================================================================================
    */






    
    cout << "monophoUL17_spikesB.root" << endl;



}
