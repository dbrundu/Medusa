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
 *  Created: 29/10/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  B0s -> J/psi  (Phi -> K+ K-)
 *          |-> mu+ mu-
 * 
 *  Fit with full analytic model, i.e. signal plus experimental artifacts
 *  (tagging, time resolution and acceptances) with analytical convolution
 *  and integration.
 *---------------------------------------------------------------------------*/

#ifndef FIT_B0S_JPSI_PHI_FULL_INL_
#define FIT_B0S_JPSI_PHI_FULL_INL_


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

// Minuit2
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

// Medusa
#include <medusa/phi_s/phis_full/FullAnalyticPhis.h>
#include <medusa/phi_s/phis_full/GenerateDataset.h>

//default namespaces
using namespace ROOT::Minuit2;



//---------------------------------
//        Main program
//---------------------------------

int main(int argv, char** argc)
{

    //------------------------------------------------------
    //   TCLAP: read the arguments from the command line
    //------------------------------------------------------

    size_t nentries = 0;
    double edm = 0.;

	try {

		TCLAP::CmdLine cmd("Command line arguments for number of events and Vegas integrator", '=');

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events", false, 1e5, "size_t");
        cmd.add(EArg);

        TCLAP::ValueArg<double> EdmArg("e", "EDM", "Estimated vertical distance to minimum", false, 0.1, "double");
		cmd.add(EdmArg);

		// Parse the argv array.
		cmd.parse(argv, argc);

		// Get the value parsed by each arg.
        nentries = EArg.getValue();
        edm = EdmArg.getValue();

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

    auto Model_2015_unbiased = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams_2015_unbiased, Spline_Knots, LowerLimit, UpperLimit);

    auto Model_2015_biased = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams_2015_biased, Spline_Knots, LowerLimit, UpperLimit);

    auto Model_2016_unbiased = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams_2016_unbiased, Spline_Knots, LowerLimit, UpperLimit);

    auto Model_2016_biased = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams_2016_biased, Spline_Knots, LowerLimit, UpperLimit);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_2015_unbiased_h;

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_2015_biased_h;

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_2016_unbiased_h;

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_2016_biased_h;


    medusa::GenerateDataset_Full(Model_2015_unbiased, dataset_2015_unbiased_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased");
    
    medusa::GenerateDataset_Full(Model_2015_biased, dataset_2015_biased_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased");

    medusa::GenerateDataset_Full(Model_2016_unbiased, dataset_2016_unbiased_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased");
    
    medusa::GenerateDataset_Full(Model_2016_biased, dataset_2016_biased_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased");


    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_2015_unbiased_d(dataset_2015_unbiased_h.size());
    hydra::copy(dataset_2015_unbiased_h, dataset_2015_unbiased_d);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_2015_biased_d(dataset_2015_biased_h.size());
    hydra::copy(dataset_2015_biased_h, dataset_2015_biased_d);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_2016_unbiased_d(dataset_2016_unbiased_h.size());
    hydra::copy(dataset_2016_unbiased_h, dataset_2016_unbiased_d);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_2016_biased_d(dataset_2016_biased_h.size());
    hydra::copy(dataset_2016_biased_h, dataset_2016_biased_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------

    std::cout << " " << std::endl;
    for( size_t i=0; i<10; i++ )
        std::cout << "Dataset_2015_unbiased_h: {" << dataset_2015_unbiased_h[i]  << "}" << std::endl;

    std::cout << " " << std::endl;
    for( size_t i=0; i<10; i++ )
        std::cout << "Dataset_2015_biased_h: {" << dataset_2015_biased_h[i]  << "}" << std::endl;

    std::cout << " " << std::endl;
    for( size_t i=0; i<10; i++ )
        std::cout << "Dataset_2016_unbiased_h: {" << dataset_2016_unbiased_h[i]  << "}" << std::endl;

    std::cout << " " << std::endl;
    for( size_t i=0; i<10; i++ )
        std::cout << "Dataset_2016_biased_h: {" << dataset_2016_biased_h[i]  << "}" << std::endl;
    std::cout << " " << std::endl;


    #ifdef _ROOT_AVAILABLE_

        // Plot of the 2015 unbiased dataset
        TH1D timedist_2015_unbiased("timedist_2015_unbiased","Decay Time; time (ps); Candidates / bin", 100, 0, 15);
        TH1D thetahdist_2015_unbiased("thetahdist_2015_unbiased","CosTheta_h; CosTheta_h; Candidates / bin", 100, -1, 1);
        TH1D thetaldist_2015_unbiased("thetaldist_2015_unbiased","CosTheta_l; CosTheta_l; Candidates / bin", 100, -1, 1);
        TH1D phidist_2015_unbiased("phidist_2015_unbiased","Phi angle; angle (rad); Candidates / bin", 100, -PI, PI);

        for(auto x : dataset_2015_unbiased_h)
        {
            timedist_2015_unbiased.Fill( (double)hydra::get<0>(x) );
            thetahdist_2015_unbiased.Fill( (double)hydra::get<1>(x) );
            thetaldist_2015_unbiased.Fill( (double)hydra::get<2>(x) );
            phidist_2015_unbiased.Fill( (double)hydra::get<3>(x) );
        }

        TCanvas canvas1_2015_unbiased("canvas1_2015_unbiased","canvas1_2015_unbiased",3200,800);
        canvas1_2015_unbiased.Divide(4,1);

        canvas1_2015_unbiased.cd(1);
        gPad->SetLogy(1);
        timedist_2015_unbiased.Draw();

        canvas1_2015_unbiased.cd(2);
        thetahdist_2015_unbiased.Draw();

        canvas1_2015_unbiased.cd(3);
        thetaldist_2015_unbiased.Draw();

        canvas1_2015_unbiased.cd(4);
        phidist_2015_unbiased.Draw();

        canvas1_2015_unbiased.SaveAs("Dataset_B0s_JpsiPhi_2015_unbiased.pdf");


        // Plot of the 2015 biased dataset
        TH1D timedist_2015_biased("timedist_2015_biased","Decay Time; time (ps); Candidates / bin", 100, 0, 15);
        TH1D thetahdist_2015_biased("thetahdist_2015_biased","CosTheta_h; CosTheta_h; Candidates / bin", 100, -1, 1);
        TH1D thetaldist_2015_biased("thetaldist_2015_biased","CosTheta_l; CosTheta_l; Candidates / bin", 100, -1, 1);
        TH1D phidist_2015_biased("phidist_2015_biased","Phi angle; angle (rad); Candidates / bin", 100, -PI, PI);

        for(auto x : dataset_2015_biased_h)
        {
            timedist_2015_biased.Fill( (double)hydra::get<0>(x) );
            thetahdist_2015_biased.Fill( (double)hydra::get<1>(x) );
            thetaldist_2015_biased.Fill( (double)hydra::get<2>(x) );
            phidist_2015_biased.Fill( (double)hydra::get<3>(x) );
        }

        TCanvas canvas1_2015_biased("canvas1_2015_biased","canvas1_2015_biased",3200,800);
        canvas1_2015_biased.Divide(4,1);

        canvas1_2015_biased.cd(1);
        gPad->SetLogy(1);
        timedist_2015_biased.Draw();

        canvas1_2015_biased.cd(2);
        thetahdist_2015_biased.Draw();

        canvas1_2015_biased.cd(3);
        thetaldist_2015_biased.Draw();

        canvas1_2015_biased.cd(4);
        phidist_2015_biased.Draw();

        canvas1_2015_biased.SaveAs("Dataset_B0s_JpsiPhi_2015_biased.pdf");


        // Plot of the 2016 unbiased dataset
        TH1D timedist_2016_unbiased("timedist_2016_unbiased","Decay Time; time (ps); Candidates / bin", 100, 0, 15);
        TH1D thetahdist_2016_unbiased("thetahdist_2016_unbiased","CosTheta_h; CosTheta_h; Candidates / bin", 100, -1, 1);
        TH1D thetaldist_2016_unbiased("thetaldist_2016_unbiased","CosTheta_l; CosTheta_l; Candidates / bin", 100, -1, 1);
        TH1D phidist_2016_unbiased("phidist_2016_unbiased","Phi angle; angle (rad); Candidates / bin", 100, -PI, PI);

        for(auto x : dataset_2016_unbiased_h)
        {
            timedist_2016_unbiased.Fill( (double)hydra::get<0>(x) );
            thetahdist_2016_unbiased.Fill( (double)hydra::get<1>(x) );
            thetaldist_2016_unbiased.Fill( (double)hydra::get<2>(x) );
            phidist_2016_unbiased.Fill( (double)hydra::get<3>(x) );
        }

        TCanvas canvas1_2016_unbiased("canvas1_2016_unbiased","canvas1_2016_unbiased",3200,800);
        canvas1_2016_unbiased.Divide(4,1);

        canvas1_2016_unbiased.cd(1);
        gPad->SetLogy(1);
        timedist_2016_unbiased.Draw();

        canvas1_2016_unbiased.cd(2);
        thetahdist_2016_unbiased.Draw();

        canvas1_2016_unbiased.cd(3);
        thetaldist_2016_unbiased.Draw();

        canvas1_2016_unbiased.cd(4);
        phidist_2016_unbiased.Draw();

        canvas1_2016_unbiased.SaveAs("Dataset_B0s_JpsiPhi_2016_unbiased.pdf");


        // Plot of the 2016 biased dataset
        TH1D timedist_2016_biased("timedist_2016_biased","Decay Time; time (ps); Candidates / bin", 100, 0, 15);
        TH1D thetahdist_2016_biased("thetahdist_2016_biased","CosTheta_h; CosTheta_h; Candidates / bin", 100, -1, 1);
        TH1D thetaldist_2016_biased("thetaldist_2016_biased","CosTheta_l; CosTheta_l; Candidates / bin", 100, -1, 1);
        TH1D phidist_2016_biased("phidist_2016_biased","Phi angle; angle (rad); Candidates / bin", 100, -PI, PI);

        for(auto x : dataset_2016_biased_h)
        {
            timedist_2016_biased.Fill( (double)hydra::get<0>(x) );
            thetahdist_2016_biased.Fill( (double)hydra::get<1>(x) );
            thetaldist_2016_biased.Fill( (double)hydra::get<2>(x) );
            phidist_2016_biased.Fill( (double)hydra::get<3>(x) );
        }

        TCanvas canvas1_2016_biased("canvas1_2016_biased","canvas1_2016_biased",3200,800);
        canvas1_2016_biased.Divide(4,1);

        canvas1_2016_biased.cd(1);
        gPad->SetLogy(1);
        timedist_2016_biased.Draw();

        canvas1_2016_biased.cd(2);
        thetahdist_2016_biased.Draw();

        canvas1_2016_biased.cd(3);
        thetaldist_2016_biased.Draw();

        canvas1_2016_biased.cd(4);
        phidist_2016_biased.Draw();

        canvas1_2016_biased.SaveAs("Dataset_B0s_JpsiPhi_2016_biased.pdf");


        // Plot of the 2015 unbiased cubic spline
        TCanvas canvas2_2015_unbiased("canvas2_2015_unbiased","canvas2_2015_unbiased",3200,800);
        canvas2_2015_unbiased.cd();
        Model_2015_unbiased.CreateHistogramPlot("Cubic Spline_2015_unbiased", "Cubic Spline_2015_unbiased", 100, 0, 20) -> Draw();
        canvas2_2015_unbiased.SaveAs("Cubic_Spline_2015_unbiased.pdf");


        // Plot of the 2015 biased cubic spline
        TCanvas canvas2_2015_biased("canvas2_2015_biased","canvas2_2015_biased",3200,800);
        canvas2_2015_biased.cd();
        Model_2015_biased.CreateHistogramPlot("Cubic Spline_2015_biased", "Cubic Spline_2015_biased", 100, 0, 20) -> Draw();
        canvas2_2015_biased.SaveAs("Cubic_Spline_2015_biased.pdf");


        // Plot of the 2016 unbiased cubic spline
        TCanvas canvas2_2016_unbiased("canvas2_2016_unbiased","canvas2_2016_unbiased",3200,800);
        canvas2_2016_unbiased.cd();
        Model_2016_unbiased.CreateHistogramPlot("Cubic Spline_2016_unbiased", "Cubic Spline_2016_unbiased", 100, 0, 20) -> Draw();
        canvas2_2016_unbiased.SaveAs("Cubic_Spline_2016_unbiased.pdf");


        // Plot of the 2016 biased cubic spline
        TCanvas canvas2_2016_biased("canvas2_2016_biased","canvas2_2016_biased",3200,800);
        canvas2_2016_biased.cd();
        Model_2016_biased.CreateHistogramPlot("Cubic Spline_2016_biased", "Cubic Spline_2016_biased", 100, 0, 20) -> Draw();
        canvas2_2016_biased.SaveAs("Cubic_Spline_2016_biased.pdf");

    #endif //_ROOT_AVAILABLE_


    //---------------------------------
    //          PDF generation
    //---------------------------------

    // Integrator (it always returns the value 1.0, because the normalization is computed
    // in FullAnalyticPhis.h. This choice is justified by the fact that Hydra does not support
    // a normalization factor which depends from the experimental variables)
    auto Integrator = hydra::AnalyticalIntegral< medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                                                                        qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto PDF_2015_unbiased = hydra::make_pdf(Model_2015_unbiased, Integrator);
    auto PDF_2015_biased = hydra::make_pdf(Model_2015_biased, Integrator);
    auto PDF_2016_unbiased = hydra::make_pdf(Model_2016_unbiased, Integrator);
    auto PDF_2016_biased = hydra::make_pdf(Model_2016_biased, Integrator);


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn_2015_unbiased = hydra::make_loglikehood_fcn(PDF_2015_unbiased, dataset_2015_unbiased_d);
    fcn_2015_unbiased.SetFcnMaxValue(2.22507e+12);
    
    auto fcn_2015_biased = hydra::make_loglikehood_fcn(PDF_2015_biased, dataset_2015_biased_d);
    fcn_2015_biased.SetFcnMaxValue(2.22507e+12);
    
    auto fcn_2016_unbiased = hydra::make_loglikehood_fcn(PDF_2016_unbiased, dataset_2016_unbiased_d);
    fcn_2016_unbiased.SetFcnMaxValue(2.22507e+12);
    
    auto fcn_2016_biased = hydra::make_loglikehood_fcn(PDF_2016_biased, dataset_2016_biased_d);
    fcn_2016_biased.SetFcnMaxValue(2.22507e+12);

    auto sim_fcn = hydra::make_simultaneous_fcn(fcn_2015_unbiased, fcn_2015_biased, fcn_2016_unbiased, fcn_2016_biased);


    //---------------------------------
    //          fit by Minuit2
    //---------------------------------

    // print level
    MnPrint::SetLevel(3);

    // minimization strategy
	MnStrategy strategy(2);

    // create Migrad minimizer
	MnMigrad minimize(sim_fcn, sim_fcn.GetParameters().GetMnState(), strategy);

	// print parameters before fitting
	std::cout << sim_fcn.GetParameters().GetMnState() << std::endl;

	// minimize and profile the time
	auto start = std::chrono::high_resolution_clock::now();

	FunctionMinimum minimum = FunctionMinimum( minimize(std::numeric_limits<unsigned int>::max(), edm));

	auto stop = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed = stop - start;

	// print minuit result
	std::cout << " minimum: " << minimum << std::endl;

	//time
	std::cout << "-----------------------------------------"<< std::endl;
	std::cout << "| Time (ms) ="<< elapsed.count()          << std::endl;
	std::cout << "-----------------------------------------"<< std::endl;

    return 0;

} // main

#endif // FIT_B0S_JPSI_PHI_FULL_INL_