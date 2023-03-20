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
 *  Created: 18/08/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  Simple benchmarks for fit_sim_phis_full_model.inl
 *---------------------------------------------------------------------------*/

#ifndef BENCHMARKS_SIMPLE_SIM_PHIS_FULL_MODEL_INL_
#define BENCHMARKS_SIMPLE_SIM_PHIS_FULL_MODEL_INL_


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
#include <medusa/phi_s/Print.h>
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

    auto Model_2015_unbiased_S1 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S1, ExpParams_2015_unbiased_S1, Spline_Knots, LowerLimit, UpperLimit);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S1_h;

    medusa::GenerateDataset_Full(Model_2015_unbiased_S1, dts_2015_unbiased_S1_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S1");
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S1_d(dts_2015_unbiased_S1_h.size());
    hydra::copy(dts_2015_unbiased_S1_h, dts_2015_unbiased_S1_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------

    #ifdef _ROOT_AVAILABLE_

        // Check the datasets
        medusa::print::PrintDataset(dts_2015_unbiased_S1_h, "2015_unbiased_S1");

        // Plot the 2015 dataset with the S-wave in the first mass bin
        medusa::print::PlotDataset(dts_2015_unbiased_S1_h, "2015_unbiased_S1");

        // Plot of the 2015 unbiased cubic spline
        TCanvas canvas2_2015_unbiased_S1("canvas2_2015_unbiased_S1","canvas2_2015_unbiased_S1",3200,800);
        canvas2_2015_unbiased_S1.cd();
        Model_2015_unbiased_S1.CreateHistogramPlot("Cubic Spline_2015_unbiased_S1", "Cubic Spline_2015_unbiased_S1", 100, 0, 20) -> Draw();
        canvas2_2015_unbiased_S1.SaveAs("Cubic_Spline_2015_unbiased_S1.pdf");

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
    auto model_PDF = hydra::make_pdf(Model_2015_unbiased_S1, integrator);


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn = hydra::make_loglikehood_fcn(model_PDF, dts_2015_unbiased_S1_d);

    fcn.SetFcnMaxValue(2.22507e+12);


    //---------------------------------
    //          Benchmarks
    //---------------------------------

    // print functor evaluation on 1 event
    auto x = dts_2015_unbiased_S1_d[0];

    auto start_functor = std::chrono::high_resolution_clock::now();
    Model_2015_unbiased_S1(x);
    auto stop_functor = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_functor = stop_functor - start_functor;

    std::cout << "Functor = " << Model_2015_unbiased_S1(x) << std::endl;
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

#endif // BENCHMARKS_SIMPLE_SIM_PHIS_FULL_MODEL_INL_