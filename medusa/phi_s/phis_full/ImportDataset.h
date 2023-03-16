/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu,
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto,
 *                      Alessandro Maria Ricci
 *
 *   This file is part of Medusa Library.
 *
 *   Medusa is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Medusa is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Medusa.  If not, see <http://www.gnu.org/licenses/>.
 *
 *---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------
 *  Created: 17/02/2023
 *
 *  Author: Andrea Contu and Alessandro Maria Ricci
 *
 *  This library contains the function which reads the data from a tuple
 *  and separates them into different datasets.
 *---------------------------------------------------------------------------*/

#ifndef IMPORT_DATASET_H_
#define IMPORT_DATASET_H_

// std
#include <string>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

// Hydra
#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/multivector.h>
#include <hydra/Placeholders.h>
#include <hydra/Lambda.h>

// Medusa
#include <medusa/phi_s/Parameters.h>

// default namespace
using namespace ROOT;


namespace medusa {

    /* Function that imports the dataset for the full analytic model (time, helicity angles,
     * tagging, time resolution and acceptances) of the decay
     *
     *                          B0s -> J/psi (Phi -> K+ K-)
     *                                  |-> mu+ mu-
     * Parameters:
     * path             = path to the root file containing the dataset
     * dts_unbiased_SN  = containers which will contain the unbiased datasets
     * dts_biased_SN    = containers which will contain the biased datasets
     * wgt_unbiased_SN  = containers which will contain the unbiased weights
     * wgt_biased_SN    = containers which will contain the biased weights
     * mkk_bins         = bins of the masses of the K+K-
     * description      = description of the dataset
     * print            = true or false to have an additional printout
     */

    template <typename Container1,
              typename Container2,
              typename Container3>
    void ImportDataset(TString path,
                         Container1& dts_unbiased_S1, Container1& dts_unbiased_S2, Container1& dts_unbiased_S3,
                         Container1& dts_unbiased_S4, Container1& dts_unbiased_S5, Container1& dts_unbiased_S6,
                         Container1& dts_biased_S1,   Container1& dts_biased_S2,   Container1& dts_biased_S3,
                         Container1& dts_biased_S4,   Container1& dts_biased_S5,   Container1& dts_biased_S6,
                         Container2& wgt_unbiased_S1, Container2& wgt_unbiased_S2, Container2& wgt_unbiased_S3,
                         Container2& wgt_unbiased_S4, Container2& wgt_unbiased_S5, Container2& wgt_unbiased_S6,
                         Container2& wgt_biased_S1,   Container2& wgt_biased_S2,   Container2& wgt_biased_S3,
                         Container2& wgt_biased_S4,   Container2& wgt_biased_S5,   Container2& wgt_biased_S6,
                         Container3& mKK_bins, std::string description, bool print = false)
    {
        // open the root file
        TFile *infile=new TFile(path);
//        infile->ls();

        // import the tree
        TTree *tree = (TTree*)infile->Get("DecayTree");
//        tree->Print();

        // variables to be read
        dtime_t decay_time;
        costheta_h_t costheta_h;
        costheta_l_t costheta_l;
        phi_t phi_angle;
        qOS_t qOS;
        qSS_t qSS;
        etaOS_t etaOS;
        etaSS_t etaSS;
        delta_t delta_time;
        mKK_t mKK;
        weight_t weight;
        Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS Decision1;
        Jpsi_Hlt1DiMuonHighMassDecision_TOS Decision2;
        B_Hlt1TrackMuonDecision_TOS Decision3;
        B_Hlt1TwoTrackMVADecision_TOS Decision4;

        // set the addresses of the variables in the tree
        tree->SetBranchAddress("time", &decay_time);
        tree->SetBranchAddress("helcosthetaK", &costheta_h);
        tree->SetBranchAddress("helcosthetaL", &costheta_l);
        tree->SetBranchAddress("helphi", &phi_angle);
        tree->SetBranchAddress("OS_Combination_DEC", &qOS);
        tree->SetBranchAddress("B_SSKaonLatest_TAGDEC", &qSS);
        tree->SetBranchAddress("OS_Combination_ETA", &etaOS);
        tree->SetBranchAddress("B_SSKaonLatest_TAGETA", &etaSS);
        tree->SetBranchAddress("sigmat", &delta_time);
        tree->SetBranchAddress("X_M", &mKK);
        tree->SetBranchAddress("sw", &weight);
        tree->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS", &Decision1);
        tree->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_TOS", &Decision2);
        tree->SetBranchAddress("B_Hlt1TrackMuonDecision_TOS", &Decision3);
        tree->SetBranchAddress("B_Hlt1TwoTrackMVADecision_TOS", &Decision4);

        // lambda which computes the S-wave bin
        auto SwaveDecision = hydra::wrap_lambda(
            [mKK_bins]__hydra_dual__(mKK_t mKK)
            {
                if(mKK > mKK_bins[0][0] && mKK < mKK_bins[0][1]) return 1;
                else if(mKK > mKK_bins[1][0] && mKK < mKK_bins[1][1]) return 2;
                else if(mKK > mKK_bins[2][0] && mKK < mKK_bins[2][1]) return 3;
                else if(mKK > mKK_bins[3][0] && mKK < mKK_bins[3][1]) return 4;
                else if(mKK > mKK_bins[4][0] && mKK < mKK_bins[4][1]) return 5;
                else if(mKK > mKK_bins[5][0] && mKK < mKK_bins[5][1]) return 6;
                return 0;
            });

        // container to save the weights and variables to compute alpha (See formula on page 108
        // of Linn's thesis: https://www.physi.uni-heidelberg.de/Publications/linn_thesis.pdf)
        hydra::host::vector<weight_t> weights_unbiased_S1;
        hydra::host::vector<weight_t> weights_unbiased_S2;
        hydra::host::vector<weight_t> weights_unbiased_S3;
        hydra::host::vector<weight_t> weights_unbiased_S4;
        hydra::host::vector<weight_t> weights_unbiased_S5;
        hydra::host::vector<weight_t> weights_unbiased_S6;
        hydra::host::vector<weight_t> weights_biased_S1;
        hydra::host::vector<weight_t> weights_biased_S2;
        hydra::host::vector<weight_t> weights_biased_S3;
        hydra::host::vector<weight_t> weights_biased_S4;
        hydra::host::vector<weight_t> weights_biased_S5;
        hydra::host::vector<weight_t> weights_biased_S6;

        double alphaNum_unbiased_S1 = 0.;
        double alphaNum_unbiased_S2 = 0.;
        double alphaNum_unbiased_S3 = 0.;
        double alphaNum_unbiased_S4 = 0.;
        double alphaNum_unbiased_S5 = 0.;
        double alphaNum_unbiased_S6 = 0.;
        double alphaNum_biased_S1 = 0.;
        double alphaNum_biased_S2 = 0.;
        double alphaNum_biased_S3 = 0.;
        double alphaNum_biased_S4 = 0.;
        double alphaNum_biased_S5 = 0.;
        double alphaNum_biased_S6 = 0.;

        double alphaDenom_unbiased_S1 = 0.;
        double alphaDenom_unbiased_S2 = 0.;
        double alphaDenom_unbiased_S3 = 0.;
        double alphaDenom_unbiased_S4 = 0.;
        double alphaDenom_unbiased_S5 = 0.;
        double alphaDenom_unbiased_S6 = 0.;
        double alphaDenom_biased_S1 = 0.;
        double alphaDenom_biased_S2 = 0.;
        double alphaDenom_biased_S3 = 0.;
        double alphaDenom_biased_S4 = 0.;
        double alphaDenom_biased_S5 = 0.;
        double alphaDenom_biased_S6 = 0.;
        
        double alpha_unbiased_S1 = 0.;
        double alpha_unbiased_S2 = 0.;
        double alpha_unbiased_S3 = 0.;
        double alpha_unbiased_S4 = 0.;
        double alpha_unbiased_S5 = 0.;
        double alpha_unbiased_S6 = 0.;
        double alpha_biased_S1 = 0.;
        double alpha_biased_S2 = 0.;
        double alpha_biased_S3 = 0.;
        double alpha_biased_S4 = 0.;
        double alpha_biased_S5 = 0.;
        double alpha_biased_S6 = 0.;


        // read the variables and fill the appropriate containers
        for(int i = 0; i < tree->GetEntries(); i++)
        {
            // read a row of the dataset
            tree->GetEntry(i);

            // decide if unbiased, biased or nothing
            if(Decision1 && Decision2)  // unbiased
            {
                // create the tupla
                auto tupla = hydra::make_tuple(decay_time, costheta_h, costheta_l, phi_angle, qOS, qSS, etaOS, etaSS, delta_time);

                // compute the S-wave bin and fill the appropriate containers
                auto SwaveBin = SwaveDecision(mKK);
                switch (SwaveBin)
                {
                    case 1:
                    {
                        dts_unbiased_S1.push_back(tupla);
                        weights_unbiased_S1.push_back(weight);
                        alphaNum_unbiased_S1 += weight;
                        alphaDenom_unbiased_S1 += weight*weight;
                        break;
                    }
                    case 2:
                    {
                        dts_unbiased_S2.push_back(tupla);
                        weights_unbiased_S2.push_back(weight);
                        alphaNum_unbiased_S2 += weight;
                        alphaDenom_unbiased_S2 += weight*weight;
                        break;
                    }
                    case 3:
                    {
                        dts_unbiased_S3.push_back(tupla);
                        weights_unbiased_S3.push_back(weight);
                        alphaNum_unbiased_S3 += weight;
                        alphaDenom_unbiased_S3 += weight*weight;
                        break;
                    }
                    case 4:
                    {
                        dts_unbiased_S4.push_back(tupla);
                        weights_unbiased_S4.push_back(weight);
                        alphaNum_unbiased_S4 += weight;
                        alphaDenom_unbiased_S4 += weight*weight;
                        break;
                    }
                    case 5:
                    {
                        dts_unbiased_S5.push_back(tupla);
                        weights_unbiased_S5.push_back(weight);
                        alphaNum_unbiased_S5 += weight;
                        alphaDenom_unbiased_S5 += weight*weight;
                        break;
                    }
                    case 6:
                    {
                        dts_unbiased_S6.push_back(tupla);
                        weights_unbiased_S6.push_back(weight);
                        alphaNum_unbiased_S6 += weight;
                        alphaDenom_unbiased_S6 += weight*weight;
                        break;
                    }
                    default:
                        break;
                } // switch
            } // if
            else if(Decision1 && !Decision2 && (Decision3 || Decision4))  // biased
            {
                // create the tupla
                auto tupla = hydra::make_tuple(decay_time, costheta_h, costheta_l, phi_angle, qOS, qSS, etaOS, etaSS, delta_time);
                
                // compute the S-wave bin and fill the appropriate containers
                auto SwaveBin = SwaveDecision(mKK);
                switch (SwaveBin)
                {
                    case 1:
                    {
                        dts_biased_S1.push_back(tupla);
                        weights_biased_S1.push_back(weight);
                        alphaNum_biased_S1 += weight;
                        alphaDenom_biased_S1 += weight*weight;
                        break;
                    }
                    case 2:
                    {
                        dts_biased_S2.push_back(tupla);
                        weights_biased_S2.push_back(weight);
                        alphaNum_biased_S2 += weight;
                        alphaDenom_biased_S2 += weight*weight;
                        break;
                    }
                    case 3:
                    {
                        dts_biased_S3.push_back(tupla);
                        weights_biased_S3.push_back(weight);
                        alphaNum_biased_S3 += weight;
                        alphaDenom_biased_S3 += weight*weight;
                        break;
                    }
                    case 4:
                    {
                        dts_biased_S4.push_back(tupla);
                        weights_biased_S4.push_back(weight);
                        alphaNum_biased_S4 += weight;
                        alphaDenom_biased_S4 += weight*weight;
                        break;
                    }
                    case 5:
                    {
                        dts_biased_S5.push_back(tupla);
                        weights_biased_S5.push_back(weight);
                        alphaNum_biased_S5 += weight;
                        alphaDenom_biased_S5 += weight*weight;
                        break;
                    }
                    case 6:
                    {
                        dts_biased_S6.push_back(tupla);
                        weights_biased_S6.push_back(weight);
                        alphaNum_biased_S6 += weight;
                        alphaDenom_biased_S6 += weight*weight;
                        break;
                    }
                    default:
                        break;
                } // switch
            } // if
        } // for

        // compute alpha (See formula on page 108 of Linn's thesis:
        // https://www.physi.uni-heidelberg.de/Publications/linn_thesis.pdf)
        alpha_unbiased_S1 = alphaNum_unbiased_S1 / alphaDenom_unbiased_S1;
        alpha_unbiased_S2 = alphaNum_unbiased_S2 / alphaDenom_unbiased_S2;
        alpha_unbiased_S3 = alphaNum_unbiased_S3 / alphaDenom_unbiased_S3;
        alpha_unbiased_S4 = alphaNum_unbiased_S4 / alphaDenom_unbiased_S4;
        alpha_unbiased_S5 = alphaNum_unbiased_S5 / alphaDenom_unbiased_S5;
        alpha_unbiased_S6 = alphaNum_unbiased_S6 / alphaDenom_unbiased_S6;

        alpha_biased_S1 = alphaNum_biased_S1 / alphaDenom_biased_S1;
        alpha_biased_S2 = alphaNum_biased_S2 / alphaDenom_biased_S2;
        alpha_biased_S3 = alphaNum_biased_S3 / alphaDenom_biased_S3;
        alpha_biased_S4 = alphaNum_biased_S4 / alphaDenom_biased_S4;
        alpha_biased_S5 = alphaNum_biased_S5 / alphaDenom_biased_S5;
        alpha_biased_S6 = alphaNum_biased_S6 / alphaDenom_biased_S6;

        // insert the weights in the appropriate containers
        for(size_t i = 0; i < weights_unbiased_S1.size(); i++)
            wgt_unbiased_S1.push_back(alpha_unbiased_S1*weights_unbiased_S1[i]);
        for(size_t i = 0; i < weights_unbiased_S2.size(); i++)
            wgt_unbiased_S2.push_back(alpha_unbiased_S2*weights_unbiased_S2[i]);
        for(size_t i = 0; i < weights_unbiased_S3.size(); i++)
            wgt_unbiased_S3.push_back(alpha_unbiased_S3*weights_unbiased_S3[i]);
        for(size_t i = 0; i < weights_unbiased_S4.size(); i++)
            wgt_unbiased_S4.push_back(alpha_unbiased_S4*weights_unbiased_S4[i]);
        for(size_t i = 0; i < weights_unbiased_S5.size(); i++)
            wgt_unbiased_S5.push_back(alpha_unbiased_S5*weights_unbiased_S5[i]);
        for(size_t i = 0; i < weights_unbiased_S6.size(); i++)
            wgt_unbiased_S6.push_back(alpha_unbiased_S6*weights_unbiased_S6[i]);
        
        for(size_t i = 0; i < weights_biased_S1.size(); i++)
            wgt_biased_S1.push_back(alpha_biased_S1*weights_biased_S1[i]);
        for(size_t i = 0; i < weights_biased_S2.size(); i++)
            wgt_biased_S2.push_back(alpha_biased_S2*weights_biased_S2[i]);
        for(size_t i = 0; i < weights_biased_S3.size(); i++)
            wgt_biased_S3.push_back(alpha_biased_S3*weights_biased_S3[i]);
        for(size_t i = 0; i < weights_biased_S4.size(); i++)
            wgt_biased_S4.push_back(alpha_biased_S4*weights_biased_S4[i]);
        for(size_t i = 0; i < weights_biased_S5.size(); i++)
            wgt_biased_S5.push_back(alpha_biased_S5*weights_biased_S5[i]);
        for(size_t i = 0; i < weights_biased_S6.size(); i++)
            wgt_biased_S6.push_back(alpha_biased_S6*weights_biased_S6[i]);
        


        // close the root file
        infile->Close();

        if(print)
        {
            // output
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "-----------Import on Device -----------------------"  << std::endl;
            std::cout << "| B0s -> J/psi Phi -> mu+ mu- K+ K-"                  << std::endl;
            std::cout << "| Dataset: " << description                           << std::endl;
            std::cout << "| Unbiased:"                                          << std::endl;
            std::cout << "| Number of events for S1: "<< dts_unbiased_S1.size() << std::endl;
            std::cout << "| Number of events for S2: "<< dts_unbiased_S2.size() << std::endl;
            std::cout << "| Number of events for S3: "<< dts_unbiased_S3.size() << std::endl;
            std::cout << "| Number of events for S4: "<< dts_unbiased_S4.size() << std::endl;
            std::cout << "| Number of events for S5: "<< dts_unbiased_S5.size() << std::endl;
            std::cout << "| Number of events for S6: "<< dts_unbiased_S6.size() << std::endl;
            std::cout << "| biased:"                                            << std::endl;
            std::cout << "| Number of events for S1: "<< dts_biased_S1.size()   << std::endl;
            std::cout << "| Number of events for S2: "<< dts_biased_S2.size()   << std::endl;
            std::cout << "| Number of events for S3: "<< dts_biased_S3.size()   << std::endl;
            std::cout << "| Number of events for S4: "<< dts_biased_S4.size()   << std::endl;
            std::cout << "| Number of events for S5: "<< dts_biased_S5.size()   << std::endl;
            std::cout << "| Number of events for S6: "<< dts_biased_S6.size()   << std::endl;
            std::cout << "---------------------------------------------------"  << std::endl;
        }
    } // ImportDataset()
} // namespace medusa

#endif // IMPORT_DATASET_H_