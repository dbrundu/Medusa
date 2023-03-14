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
 *  Simple benchmarks for fit_sim_phis_only_signal.inl
 *---------------------------------------------------------------------------*/

#ifndef BENCHMARKS_SIMPLE_SIM_PHIS_ONLY_SIGNAL_INL_
#define BENCHMARKS_SIMPLE_SIM_PHIS_ONLY_SIGNAL_INL_


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
#include <medusa/generic/Print.h>
#endif //_ROOT_AVAILABLE_

// Medusa
#include <medusa/phi_s/phis_signal/PhisSignal.h>
#include <medusa/phi_s/phis_signal/GenerateDataset.h>



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

    auto Model = medusa::PhisSignal<B0sbar, dtime_t, costheta_h_t, costheta_l_t, phi_t>(ModelParams_S1);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t> , hydra::host::sys_t> dataset_h;

    medusa::GenerateDataset_SignalOnly(Model, dataset_h, nentries, nentries, LowerLimit, UpperLimit);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h, dataset_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------

    #ifdef _ROOT_AVAILABLE_

        // Plot the dataset with the S-wave in the first mass bin
        medusa::print::PrintDataset_B0s(dataset_h, "S1");

    #endif //_ROOT_AVAILABLE_


    //---------------------------------
    //          PDF generation
    //---------------------------------

    auto Integrator = hydra::AnalyticalIntegral<
                                medusa::PhisSignal< B0sbar, dtime_t, costheta_h_t, costheta_l_t, phi_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto model_PDF = hydra::make_pdf(Model, Integrator);

    // print the normalization factor
    std::cout << " "                                                            << std::endl;
    std::cout << "---------------------------------------------"                << std::endl;
    std::cout << "The norm factor is: " << model_PDF.GetNorm()                  << std::endl;
    std::cout << "The error in norm factor is: " << model_PDF.GetNormError()    << std::endl;
    std::cout << "---------------------------------------------"                << std::endl;


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn = hydra::make_loglikehood_fcn(model_PDF, dataset_d);


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


    // print non-cached integration with Vegas
    auto start_integrator = std::chrono::high_resolution_clock::now();
    auto NormFactor = Integrator.Integrate(Model);
    auto stop_integrator = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_integrator = stop_integrator - start_integrator;

    std::cout << "Integral = " << NormFactor.first << std::endl;
    std::cout << "Time (ms) = " << elapsed_integrator.count() << std::endl;


    // distortion of the parameters to avoid the cached integration
    for(size_t i=0; i<parameters.size(); i++)
    {
       parameters[i] *= 1.001;
    }


    // print non-cached fcn evaluation + non-cauched integration
    auto start_fcn = std::chrono::high_resolution_clock::now();
    fcn(parameters);
    auto stop_fcn = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_fcn = stop_fcn - start_fcn;

    std::cout << "fcn = " << fcn(parameters) << std::endl;
    std::cout << "Time (ms) = " << elapsed_fcn.count() << std::endl;

    return 0;

} // main

#endif // BENCHMARKS_SIMPLE_SIM_PHIS_ONLY_SIGNAL_INL_