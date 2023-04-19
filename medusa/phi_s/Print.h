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
/*--------------------------------------------------------
 *  Created on: 18/05/2022
 *
 *  Author: Alessandro Maria Ricci
 *
 *  This library contains the print functions.
 *--------------------------------------------------------*/

#ifndef MEDUSA_PRINT_H_
#define MEDUSA_PRINT_H_

// std
#include <iostream>
#include <string>
#include <vector>

// ROOT
#include <TROOT.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TString.h>
#include <TStyle.h>

//Hydra
#include <hydra/multivector.h>
#include <hydra/Tuple.h>
#include <hydra/Lambda.h>


namespace medusa {
    namespace print {

        /* Function that prints the Monte Carlo dataset (time and helicity angles)
         * for the decay
         *
         *                    B0s -> J/psi (Phi -> K+ K-)
         *                            |-> mu+ mu-
         * Parameters:
         * dataset      = container which contains the Monte Carlo dataset
         * description  = description of the dataset
         */

        template<typename Container>
        size_t PrintDataset(Container& dataset, std::string description)
        {
            // Print the first 10 lines of the dataset
            if(dataset.size() > 10)
            {
                std::cout << " " << std::endl;
                for( size_t i=0; i<10; i++ )
                    std::cout << "Dataset_" << description << ": {" << dataset[i]  << "}" << std::endl;
            }
            else std::cout << "The container " << description << " is empty!" << std::endl;

            return dataset.size();
        }


        /* Function that plots the Monte Carlo dataset (time and helicity angles)
         * for the decay
         *
         *                    B0s -> J/psi (Phi -> K+ K-)
         *                            |-> mu+ mu-
         * Parameters:
         * dataset      = container which contains the Monte Carlo dataset
         * description  = description of the dataset
         */

        template<typename Container>
        void PlotDataset(Container& dataset, std::string description)
        {
            // Create the histograms
            TH1D timedist("timedist","Decay Time; time (ps); Candidates / bin", 100, 0, 15);
            TH1D thetahdist("thetahdist","CosTheta_h; CosTheta_h; Candidates / bin", 100, -1, 1);
            TH1D thetaldist("thetaldist","CosTheta_l; CosTheta_l; Candidates / bin", 100, -1, 1);
            TH1D phidist("phidist","Phi Angle; phi (rad); Candidates / bin", 100, -PI, PI);

            for(auto x : dataset)
            {
                timedist.Fill( (double)hydra::get<0>(x) );
                thetahdist.Fill( (double)hydra::get<1>(x) );
                thetaldist.Fill( (double)hydra::get<2>(x) );
                phidist.Fill( (double)hydra::get<3>(x) );
            }

            TCanvas canvas1("canvas1","canvas1",3200,800);
            canvas1.Divide(4,1);

            canvas1.cd(1);
            gPad->SetLogy(1);
            timedist.Draw();

            canvas1.cd(2);
            thetahdist.Draw();

            canvas1.cd(3);
            thetaldist.Draw();

            canvas1.cd(4);
            phidist.Draw();

            std::string FileName = "Dataset_B0s_" + description + ".pdf";

            canvas1.SaveAs(FileName.c_str());
        }
    } // namespace print
} // namespace medusa

#endif // MEDUSA_PRINT_H