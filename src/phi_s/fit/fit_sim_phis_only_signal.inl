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
/*----------------------------------------
 *  Created: 26/05/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  B0s -> J/psi  (Phi -> K+ K-)
 *          |-> mu+ mu-
 * 
 *  Fit with signal-only model
 *----------------------------------------*/

#ifndef FIT_SIM_PHIS_ONLY_SIGNAL_INL_
#define FIT_SIM_PHIS_ONLY_SIGNAL_INL_


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

// Minuit2
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

// Medusa
#include <medusa/phi_s/phis_signal/PhisSignal.h>
#include <medusa/phi_s/phis_signal/GenerateDataset.h>


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

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events", false, 5e5, "size_t");
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

    auto Model = medusa::PhisSignal<B0sbar, dtime_t, costheta_h_t, costheta_l_t, phi_t>(ModelParams_S1);


    //---------------------------------
    //      dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t> , hydra::host::sys_t> dataset_h;

    medusa::GenerateDataset_SignalOnly(Model, dataset_h, nentries, nentries, LowerLimit, UpperLimit);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h, dataset_d);


    //-----------------------------------------
    //      Print and plot the dataset
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
    auto Model_PDF = hydra::make_pdf(Model, Integrator);

    // print the normalization factor
    std::cout << " "                                                            << std::endl;
    std::cout << "---------------------------------------------"                << std::endl;
    std::cout << "The norm factor is: " << Model_PDF.GetNorm()                  << std::endl;
    std::cout << "The error in norm factor is: " << Model_PDF.GetNormError()    << std::endl;
    std::cout << "---------------------------------------------"                << std::endl;


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);

    fcn.SetFcnMaxValue(2.22507e+12);


    //---------------------------------
    //          fit by Minuit2
    //---------------------------------

    // print level (command for ROOT 6.22.08)
    MnPrint::SetLevel(3);
    // print level (command for ROOT 6.26)
//    MnPrint::SetGlobalLevel(2);

    // minimization strategy
	MnStrategy strategy(2);

    // create Migrad minimizer
	MnMigrad minimize(fcn, fcn.GetParameters().GetMnState(), strategy);

	// print parameters before fitting
	std::cout << fcn.GetParameters().GetMnState() << std::endl;

	// minimize and profile the time
	auto start = std::chrono::high_resolution_clock::now();

	FunctionMinimum minimum = FunctionMinimum( minimize(std::numeric_limits<unsigned int>::max(), edm));

	auto stop = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed = stop - start;

	// print minuit result
	std::cout << " minimum: " << minimum << std::endl;

	//time
	std::cout << "-----------------------------------------"<< std::endl;
	std::cout << "| Time (ms) = "<< elapsed.count()         << std::endl;
	std::cout << "-----------------------------------------"<< std::endl;

    return 0;

} // main

#endif // FIT_SIM_PHIS_ONLY_SIGNAL_INL_