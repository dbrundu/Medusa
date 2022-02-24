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
 * benchmarks2_B0s_Jpsi_phi.inl
 *
 *  Created: 18/08/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  Simple benchmarks for fit_B0s_Jpsi_phi_full.inl
 *---------------------------------------------------------------------------*/

#ifndef BENCHMARKS_SIMPLE_B0S_JPSI_PHI_FULL_INL_
#define BENCHMARKS_SIMPLE_B0S_JPSI_PHI_FULL_INL_


//---------------------------------
//    Libraries and namespaces
//---------------------------------

// std
#include <iostream>
#include <chrono>

// command line arguments
#include <tclap/CmdLine.h>

// Hydra
#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/multivector.h>
#include <hydra/Parameter.h>
#include <hydra/Function.h>
#include <hydra/Lambda.h>
#include <hydra/Pdf.h>
#include <hydra/LogLikelihoodFCN.h>

// ROOT
#ifdef _ROOT_AVAILABLE_

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

#endif //_ROOT_AVAILABLE_

// Medusa
#include <medusa/phi_s/phis_full/FullAnalyticPhis.h>
#include <medusa/phi_s/phis_full/GenerateDataset.h>



//---------------------------------
//        Main program
//---------------------------------

int main(int argv, char** argc)
{


    //------------------------------------------------------
    //   TCLAP: read the arguments from the command line
    //------------------------------------------------------

    size_t nentries = 0;

	try {

		TCLAP::CmdLine cmd("Command line arguments for number of events and Vegas integrator", '=');

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events", false, 5e5, "size_t");
        cmd.add(EArg);

		// Parse the argv array.
		cmd.parse(argv, argc);

		// Get the value parsed by each arg.
        nentries = EArg.getValue();

	}
	catch (TCLAP::ArgException &e)  {
		std::cerr << " error: "  << e.error()
				  << " for arg " << e.argId()
				  << std::endl;
	}


    //---------------------------------
    //      Hydra print level
    //---------------------------------

    hydra::Print::SetLevel(hydra::WARNING);


    //---------------------------------
    //      Model generation
    //---------------------------------

    auto Model = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                            qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams, Spline_Knots,
                                                                                                            LowerLimit, UpperLimit);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_h;

    medusa::GenerateDataset_Full(Model, dataset_h, nentries, nentries, LowerLimit, UpperLimit);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h, dataset_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------
    
    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset_h: {"<< dataset_h[i]  << "}"<< std::endl;


    #ifdef _ROOT_AVAILABLE_

        // Plot of the dataset
        TH1D timedist("timedist","Decay Time; time (ps); Candidates / bin", 100, 0, 15);
        TH1D thetahdist("thetahdist","CosTheta_h; CosTheta_h; Candidates / bin", 100, -1, 1);
        TH1D thetaldist("thetaldist","CosTheta_l; CosTheta_l; Candidates / bin", 100, -1, 1);
        TH1D phidist("phidist","Phi angle; angle (rad); Candidates / bin", 100, -PI, PI);

        for(auto x : dataset_h)
        {
            timedist.Fill( (double)hydra::get<0>(x) );
            thetahdist.Fill( (double)hydra::get<1>(x) );
            thetaldist.Fill( (double)hydra::get<2>(x) );
            phidist.Fill( (double)hydra::get<3>(x) );
        }

        TCanvas canvas1("canvas","canvas",3200,800);
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

        canvas1.SaveAs("Dataset_B0s_JpsiPhi.pdf");


        // Plot of the cubic spline
        TCanvas canvas2("canvas","canvas",3200,800);
        canvas2.cd();
        Model.CreateHistogramPlot("Cubic Spline", "Cubic Spline", 100, 0, 20) -> Draw();
        canvas2.SaveAs("Cubic_Spline.pdf");

    #endif //_ROOT_AVAILABLE_


    //---------------------------------
    //          PDF generation
    //---------------------------------

    // integrator (it always returns the value 1.0, because the normalization is computed
    // in FullAnalyticPhis.h. This choice is justified by the fact that Hydra does not support
    // a normalization factor which depends from the experimental variables)
    auto integrator = hydra::AnalyticalIntegral< medusa::FullAnalyticPhis<CubicSpline,
                                                            dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto model_PDF = hydra::make_pdf(Model, integrator);


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn = hydra::make_loglikehood_fcn(model_PDF, dataset_d);

    fcn.SetFcnMaxValue(2.22507e+12);


    //---------------------------------
    //          Benchmarks
    //---------------------------------

    // print functor evaluation on 1 event
    auto x = dataset_d[0];

    auto start_functor = std::chrono::high_resolution_clock::now();
    Model(x);
    auto stop_functor = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_functor = stop_functor - start_functor;

    std::cout << "Functor = " << Model(x) << std::endl;
    std::cout << "Time (ms) = " << elapsed_functor.count() << std::endl;


    // print non-cached fcn evaluation
    auto start_fcn = std::chrono::high_resolution_clock::now();
    fcn(parameters);
    auto stop_fcn = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_fcn = stop_fcn - start_fcn;

    std::cout << "fcn = " << fcn(parameters) << std::endl;
    std::cout << "Time (ms) = " << elapsed_fcn.count() << std::endl;

    return 0;

} // main

#endif // BENCHMARKS_SIMPLE_B0S_JPSI_PHI_FULL_INL_