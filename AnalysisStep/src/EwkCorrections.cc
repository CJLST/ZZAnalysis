#include "ZZAnalysis/AnalysisStep/interface/EwkCorrections.h"
#include "TLorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

namespace EwkCorrections
{
	//Read correction table
	std::vector<std::vector<float>> readFile_and_loadEwkTable(TString dtag){
		std::ifstream myReadFile;
		std::vector<float> Table_line;
		std::vector<std::vector<float>> Table_EWK;
//		TString name;
//		TString cmssw_path;
//		cmssw_path = getenv("CMSSW_BASE");
//		TString path = cmssw_path+"/src/ZZAnalysis/AnalysisStep/src/";
//		
//		if(dtag.Contains("ZZ")) name = path+"ZZ_EwkCorrections.dat";
		myReadFile.open(dtag);
		if(!myReadFile.is_open()) std::cout << "WARNING: "+dtag+" NOT FOUND" << std::endl;
		int Start=0;
		while (!myReadFile.eof()){
			Start++;
			std::string output;
			myReadFile >> output;
			if(Start%5!=0) Table_line.push_back(atof(output.c_str()));
			if(Start%5==0){
				Table_line.push_back(atof(output.c_str()));
				Table_EWK.push_back(Table_line);
				Table_line.clear();
			}
		}
		myReadFile.close();
		return Table_EWK;
	}


	//Find the right correction in the file	
	std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat){
		//find the range of sqrt s hat (each 200 lines it changes)
		unsigned int j = 0;
		float best = 0.8E+04; //highest value of sqrt s hat in the table
		if( sqrt_s_hat > best) j = 39800; //in the very rare case where we have bigger s than our table (table is for 8TeV and we run at 13TeV)
		else{
			for(unsigned int i = 0 ; i < 40000 ; i = i+200){
				if(fabs(sqrt_s_hat - Table_EWK[i][0]) < best){
					best = fabs(sqrt_s_hat - Table_EWK[i][0]);
					j = i;
				}
				else break ;
			}
		}
		best = Table_EWK[j+199][1];
		if(t_hat > best) j = j+199; //in the very rare case where we have bigger t than our table
		else{
			best = 0.1E+09;
			for(unsigned int k = j ; k < j + 200 ; k++){
				if(fabs(t_hat - Table_EWK[k][1]) < best){
					best = fabs(t_hat - Table_EWK[k][1]);
					j = k;
				}
				else break ;
			}
		}
		std::vector<float> EWK_w2_vec;
		EWK_w2_vec.push_back(Table_EWK[j][2]); //ewk corrections for quark u/c
		EWK_w2_vec.push_back(Table_EWK[j][3]); //ewk corrections for quark d/s
		EWK_w2_vec.push_back(Table_EWK[j][4]); //ewk corrections for quark b
		return EWK_w2_vec ;
	}


	//The main function, will return the kfactor
	double getEwkCorrections(const edm::Handle<edm::View<reco::Candidate> > & particles, 
	                         const std::vector<std::vector<float>> & Table, 
	                         const GenEventInfoProduct & eventInfo,
	                         TLorentzVector Z1, TLorentzVector Z2) {
	// , double & ewkCorrections_error){
		double kFactor = 1.;

                std::vector<const reco::Candidate*> genIncomingQuarks;
		std::vector<const reco::Candidate*> genIncomingGluons;
		std::vector<const reco::Candidate*> genLeptons;
		std::vector<const reco::Candidate*> genNeutrinos;

                for( edm::View<reco::Candidate>::const_iterator genParticle = particles->begin(); genParticle != particles->end(); ++ genParticle ) {
			if(fabs(genParticle->pdgId()) >= 1 && fabs(genParticle->pdgId()) <= 5 && genParticle->status() == 21) genIncomingQuarks.push_back(&*genParticle); //status 21 : incoming particles of hardest subprocess
			if(fabs(genParticle->pdgId()) == 21 && genParticle->status() == 21) genIncomingGluons.push_back(&*genParticle);
			if(genParticle->status()==1 && (fabs(genParticle->pdgId())==11 || fabs(genParticle->pdgId())==13 || fabs(genParticle->pdgId())==15)) genLeptons.push_back(&*genParticle); //status 1 : final state
			if(genParticle->status()==1 && (fabs(genParticle->pdgId())==12 || fabs(genParticle->pdgId())==14 || fabs(genParticle->pdgId())==16)) genNeutrinos.push_back(&*genParticle);
		}

		std::sort(genIncomingQuarks.begin(), genIncomingQuarks.end(), sort_CandidatesByPt);
		std::sort(genIncomingGluons.begin(), genIncomingGluons.end(), sort_CandidatesByPt);
		std::sort(genLeptons.begin(), genLeptons.end(), sort_CandidatesByPt);
		std::sort(genNeutrinos.begin(), genNeutrinos.end(), sort_CandidatesByPt);

		if(!eventInfo.pdf()) return 1; //no corrections can be applied because we need x1 and x2 
//		if( genLeptons.size() < 2 || genNeutrinos.size() < 2) return 1; //no corrections can be applied if we don't find our two Z's
		double x1 = eventInfo.pdf()->x.first; 
		double x2 = eventInfo.pdf()->x.second; 

//		LorentzVector Z1 = genLeptons[0]->p4() + genLeptons[1]->p4(); //First Z : charged leptons
//		LorentzVector Z2 = genNeutrinos[0]->p4() + genNeutrinos[1]->p4(); //Second Z : neutrinos
		TLorentzVector ZZ = Z1+Z2;
		TLorentzVector ZZ_t(ZZ.X(),ZZ.Y(),ZZ.Z(),ZZ.T()); //Need TLorentzVectors for several methods (boosts)
		TLorentzVector Z1_t(Z1.X(),Z1.Y(),Z1.Z(),Z1.T());
		TLorentzVector Z2_t(Z2.X(),Z2.Y(),Z2.Z(),Z2.T());

		double s_hat = pow(ZZ.M(),2); // s_hat = center-of-mass energy of 2 Z system

		//Boost quarks and Z1
		TLorentzVector Z1_b = Z1_t;
		TLorentzVector p1_b, p2_b;
		double energy = 6500. ; //13 TeV in total
		p1_b.SetXYZT(0.,0.,x1*energy,x1*energy); //x1 = fraction of momentum taken by the particle initiating the hard process
		p2_b.SetXYZT(0.,0.,-x2*energy,x2*energy);
		Z1_b.Boost( -ZZ_t.BoostVector()); //Inverse Lorentz transformation, to get to the center-of-mass frame
		p1_b.Boost( -ZZ_t.BoostVector());
		p2_b.Boost( -ZZ_t.BoostVector());

		//Unitary vectors
		TLorentzVector Z1_b_u = Z1_b*(1/Z1_b.P()); //Normalized to 1
		TLorentzVector p1_b_u = p1_b*(1/p1_b.P());
		TLorentzVector p2_b_u = p2_b*(1/p2_b.P());

		//Effective beam axis
		TLorentzVector diff_p = p1_b_u - p2_b_u;
		TLorentzVector eff_beam_axis = diff_p*(1./diff_p.P());
		double cos_theta = eff_beam_axis.X()*Z1_b_u.X() + eff_beam_axis.Y()*Z1_b_u.Y() + eff_beam_axis.Z()*Z1_b_u.Z();

		double m_z = 91.1876; //Z bosons assumed to be on-shell
		double t_hat = m_z*m_z - 0.5*s_hat + cos_theta * sqrt( 0.25*s_hat*s_hat - m_z*m_z*s_hat );

		int quark_type = 0; //Flavour of incident quark
		if(genIncomingQuarks.size() > 0) quark_type = fabs(genIncomingQuarks[0]->pdgId()); //Works unless if gg->ZZ process : it shouldn't be the case as we're using POWHEG

		std::vector<float> Correction_vec = findCorrection( Table, sqrt(s_hat), t_hat ); //Extract the corrections for the values of s and t computed
		//std::cout << Correction_vec[0] << " " << Correction_vec[1] << sqrt(s_hat) << 2*m_z << std::endl;
		
		if(quark_type==1) kFactor = 1. + Correction_vec[1]; //d
		if(quark_type==2) kFactor = 1. + Correction_vec[0]; //u
		if(quark_type==3) kFactor = 1. + Correction_vec[1]; //s as d
		if(quark_type==4) kFactor = 1. + Correction_vec[0]; //c as u
		if(quark_type==5) kFactor = 1. + Correction_vec[2]; //b

		if(sqrt(s_hat)< 2*m_z) kFactor = 1.; //Off-shell cases, not corrected to avoid non-defined values for t.

		//Computing the associated error:
		//Warning, several methods could be used.
		//In Run 1, CMS used (kFactor-1)*(kFactor_QCD -1) for all rho
		//And ATLAS used : 0 for rho < 0.3 and 1 for rho >0.3
		//
		//Here is an implementation that is using a mix of the two. It may change in the future (but the change won't be critical)
//		double kFactor_QCD = 15.99/9.89; //From arXiv1105.0020
		

		//At this point, we have the relative error on the delta_ewk ( = k_ewk -1 )
		//Let's - instead - return the absolute error on k: we do delta_ewk* the_relative_errir_on_it. This gives absolute error on delta, and so on k
//		ewkCorrections_error = fabs(ewkCorrections_error*(kFactor-1));
		
//		return 1;
		return kFactor;
	}

}

