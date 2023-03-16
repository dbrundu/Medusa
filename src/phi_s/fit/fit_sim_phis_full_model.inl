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
 *  Fit with Monte Carlo dataset and full analytic model, i.e. signal plus
 *  experimental artifacts (tagging, time resolution and acceptances) with
 *  analytical convolution and integration.
 *---------------------------------------------------------------------------*/

#ifndef FIT_SIM_PHIS_FULL_MODEL_INL_
#define FIT_SIM_PHIS_FULL_MODEL_INL_


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
#include <medusa/phi_s/phis_full/FullAnalyticPhis.h>
#include <medusa/phi_s/phis_full/GenerateDataset.h>

// default namespaces
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

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events per dataset", false, 2e4, "size_t");
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

    auto Model_2015_unbiased_S1 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S1, ExpParams_2015_unbiased_S1, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_unbiased_S2 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S2, ExpParams_2015_unbiased_S2, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_unbiased_S3 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S3, ExpParams_2015_unbiased_S3, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_unbiased_S4 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S4, ExpParams_2015_unbiased_S4, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_unbiased_S5 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S5, ExpParams_2015_unbiased_S5, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_unbiased_S6 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S6, ExpParams_2015_unbiased_S6, Spline_Knots, LowerLimit, UpperLimit);


    auto Model_2015_biased_S1 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S1, ExpParams_2015_biased_S1, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_biased_S2 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S2, ExpParams_2015_biased_S2, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_biased_S3 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S3, ExpParams_2015_biased_S3, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_biased_S4 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S4, ExpParams_2015_biased_S4, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_biased_S5 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S5, ExpParams_2015_biased_S5, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2015_biased_S6 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S6, ExpParams_2015_biased_S6, Spline_Knots, LowerLimit, UpperLimit);


    auto Model_2016_unbiased_S1 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S1, ExpParams_2016_unbiased_S1, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_unbiased_S2 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S2, ExpParams_2016_unbiased_S2, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_unbiased_S3 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S3, ExpParams_2016_unbiased_S3, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_unbiased_S4 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S4, ExpParams_2016_unbiased_S4, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_unbiased_S5 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S5, ExpParams_2016_unbiased_S5, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_unbiased_S6 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S6, ExpParams_2016_unbiased_S6, Spline_Knots, LowerLimit, UpperLimit);


    auto Model_2016_biased_S1 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S1, ExpParams_2016_biased_S1, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_biased_S2 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S2, ExpParams_2016_biased_S2, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_biased_S3 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S3, ExpParams_2016_biased_S3, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_biased_S4 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S4, ExpParams_2016_biased_S4, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_biased_S5 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S5, ExpParams_2016_biased_S5, Spline_Knots, LowerLimit, UpperLimit);
    auto Model_2016_biased_S6 = medusa::FullAnalyticPhis<CubicSpline, dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_S6, ExpParams_2016_biased_S6, Spline_Knots, LowerLimit, UpperLimit);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S1_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S2_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S3_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S4_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S5_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_unbiased_S6_h;


    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_biased_S1_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_biased_S2_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_biased_S3_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_biased_S4_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_biased_S5_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2015_biased_S6_h;


    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_unbiased_S1_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_unbiased_S2_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_unbiased_S3_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_unbiased_S4_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_unbiased_S5_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_unbiased_S6_h;


    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_biased_S1_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_biased_S2_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_biased_S3_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_biased_S4_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_biased_S5_h;
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dts_2016_biased_S6_h;


    medusa::GenerateDataset_Full(Model_2015_unbiased_S1, dts_2015_unbiased_S1_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S1");
    medusa::GenerateDataset_Full(Model_2015_unbiased_S2, dts_2015_unbiased_S2_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S2");
    medusa::GenerateDataset_Full(Model_2015_unbiased_S3, dts_2015_unbiased_S3_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S3");
    medusa::GenerateDataset_Full(Model_2015_unbiased_S4, dts_2015_unbiased_S4_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S4");
    medusa::GenerateDataset_Full(Model_2015_unbiased_S5, dts_2015_unbiased_S5_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S5");
    medusa::GenerateDataset_Full(Model_2015_unbiased_S6, dts_2015_unbiased_S6_h, nentries, nentries, LowerLimit, UpperLimit, "2015 unbiased S6");

    medusa::GenerateDataset_Full(Model_2015_biased_S1, dts_2015_biased_S1_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased S1");
    medusa::GenerateDataset_Full(Model_2015_biased_S2, dts_2015_biased_S2_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased S2");
    medusa::GenerateDataset_Full(Model_2015_biased_S3, dts_2015_biased_S3_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased S3");
    medusa::GenerateDataset_Full(Model_2015_biased_S4, dts_2015_biased_S4_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased S4");
    medusa::GenerateDataset_Full(Model_2015_biased_S5, dts_2015_biased_S5_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased S5");
    medusa::GenerateDataset_Full(Model_2015_biased_S6, dts_2015_biased_S6_h, nentries, nentries, LowerLimit, UpperLimit, "2015 biased S6");

    medusa::GenerateDataset_Full(Model_2016_unbiased_S1, dts_2016_unbiased_S1_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased S1");
    medusa::GenerateDataset_Full(Model_2016_unbiased_S2, dts_2016_unbiased_S2_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased S2");
    medusa::GenerateDataset_Full(Model_2016_unbiased_S3, dts_2016_unbiased_S3_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased S3");
    medusa::GenerateDataset_Full(Model_2016_unbiased_S4, dts_2016_unbiased_S4_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased S4");
    medusa::GenerateDataset_Full(Model_2016_unbiased_S5, dts_2016_unbiased_S5_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased S5");
    medusa::GenerateDataset_Full(Model_2016_unbiased_S6, dts_2016_unbiased_S6_h, nentries, nentries, LowerLimit, UpperLimit, "2016 unbiased S6");
    
    medusa::GenerateDataset_Full(Model_2016_biased_S1, dts_2016_biased_S1_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased S1");
    medusa::GenerateDataset_Full(Model_2016_biased_S2, dts_2016_biased_S2_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased S2");
    medusa::GenerateDataset_Full(Model_2016_biased_S3, dts_2016_biased_S3_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased S3");
    medusa::GenerateDataset_Full(Model_2016_biased_S4, dts_2016_biased_S4_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased S4");
    medusa::GenerateDataset_Full(Model_2016_biased_S5, dts_2016_biased_S5_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased S5");
    medusa::GenerateDataset_Full(Model_2016_biased_S6, dts_2016_biased_S6_h, nentries, nentries, LowerLimit, UpperLimit, "2016 biased S6");


    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S1_d(dts_2015_unbiased_S1_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S2_d(dts_2015_unbiased_S2_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S3_d(dts_2015_unbiased_S3_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S4_d(dts_2015_unbiased_S4_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S5_d(dts_2015_unbiased_S5_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_unbiased_S6_d(dts_2015_unbiased_S6_h.size());
    hydra::copy(dts_2015_unbiased_S1_h, dts_2015_unbiased_S1_d);
    hydra::copy(dts_2015_unbiased_S2_h, dts_2015_unbiased_S2_d);
    hydra::copy(dts_2015_unbiased_S3_h, dts_2015_unbiased_S3_d);
    hydra::copy(dts_2015_unbiased_S4_h, dts_2015_unbiased_S4_d);
    hydra::copy(dts_2015_unbiased_S5_h, dts_2015_unbiased_S5_d);
    hydra::copy(dts_2015_unbiased_S6_h, dts_2015_unbiased_S6_d);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_biased_S1_d(dts_2015_biased_S1_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_biased_S2_d(dts_2015_biased_S2_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_biased_S3_d(dts_2015_biased_S3_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_biased_S4_d(dts_2015_biased_S4_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_biased_S5_d(dts_2015_biased_S5_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2015_biased_S6_d(dts_2015_biased_S6_h.size());
    hydra::copy(dts_2015_biased_S1_h, dts_2015_biased_S1_d);
    hydra::copy(dts_2015_biased_S2_h, dts_2015_biased_S2_d);
    hydra::copy(dts_2015_biased_S3_h, dts_2015_biased_S3_d);
    hydra::copy(dts_2015_biased_S4_h, dts_2015_biased_S4_d);
    hydra::copy(dts_2015_biased_S5_h, dts_2015_biased_S5_d);
    hydra::copy(dts_2015_biased_S6_h, dts_2015_biased_S6_d);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_unbiased_S1_d(dts_2016_unbiased_S1_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_unbiased_S2_d(dts_2016_unbiased_S2_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_unbiased_S3_d(dts_2016_unbiased_S3_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_unbiased_S4_d(dts_2016_unbiased_S4_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_unbiased_S5_d(dts_2016_unbiased_S5_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_unbiased_S6_d(dts_2016_unbiased_S6_h.size());
    hydra::copy(dts_2016_unbiased_S1_h, dts_2016_unbiased_S1_d);
    hydra::copy(dts_2016_unbiased_S2_h, dts_2016_unbiased_S2_d);
    hydra::copy(dts_2016_unbiased_S3_h, dts_2016_unbiased_S3_d);
    hydra::copy(dts_2016_unbiased_S4_h, dts_2016_unbiased_S4_d);
    hydra::copy(dts_2016_unbiased_S5_h, dts_2016_unbiased_S5_d);
    hydra::copy(dts_2016_unbiased_S6_h, dts_2016_unbiased_S6_d);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_biased_S1_d(dts_2016_biased_S1_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_biased_S2_d(dts_2016_biased_S2_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_biased_S3_d(dts_2016_biased_S3_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_biased_S4_d(dts_2016_biased_S4_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_biased_S5_d(dts_2016_biased_S5_h.size());
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t,
                                etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dts_2016_biased_S6_d(dts_2016_biased_S6_h.size());
    hydra::copy(dts_2016_biased_S1_h, dts_2016_biased_S1_d);
    hydra::copy(dts_2016_biased_S2_h, dts_2016_biased_S2_d);
    hydra::copy(dts_2016_biased_S3_h, dts_2016_biased_S3_d);
    hydra::copy(dts_2016_biased_S4_h, dts_2016_biased_S4_d);
    hydra::copy(dts_2016_biased_S5_h, dts_2016_biased_S5_d);
    hydra::copy(dts_2016_biased_S6_h, dts_2016_biased_S6_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------

    #ifdef _ROOT_AVAILABLE_

        // Plot the 2015-2016 datasets with the S-wave in the first mass bin
        medusa::print::PrintDataset_B0s(dts_2015_unbiased_S1_h, "2015_unbiased_S1");

        // Plot of the 2015 unbiased cubic spline
        TCanvas canvas2_2015_unbiased_S1("canvas2_2015_unbiased_S1","canvas2_2015_unbiased_S1",3200,800);
        canvas2_2015_unbiased_S1.cd();
        Model_2015_unbiased_S1.CreateHistogramPlot("Cubic Spline_2015_unbiased_S1", "Cubic Spline_2015_unbiased_S1", 100, 0, 20) -> Draw();
        canvas2_2015_unbiased_S1.SaveAs("Cubic_Spline_2015_unbiased_S1.pdf");

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
    auto PDF_2015_unbiased_S1 = hydra::make_pdf(Model_2015_unbiased_S1, Integrator);
    auto PDF_2015_unbiased_S2 = hydra::make_pdf(Model_2015_unbiased_S2, Integrator);
    auto PDF_2015_unbiased_S3 = hydra::make_pdf(Model_2015_unbiased_S3, Integrator);
    auto PDF_2015_unbiased_S4 = hydra::make_pdf(Model_2015_unbiased_S4, Integrator);
    auto PDF_2015_unbiased_S5 = hydra::make_pdf(Model_2015_unbiased_S5, Integrator);
    auto PDF_2015_unbiased_S6 = hydra::make_pdf(Model_2015_unbiased_S6, Integrator);

    auto PDF_2015_biased_S1 = hydra::make_pdf(Model_2015_biased_S1, Integrator);
    auto PDF_2015_biased_S2 = hydra::make_pdf(Model_2015_biased_S2, Integrator);
    auto PDF_2015_biased_S3 = hydra::make_pdf(Model_2015_biased_S3, Integrator);
    auto PDF_2015_biased_S4 = hydra::make_pdf(Model_2015_biased_S4, Integrator);
    auto PDF_2015_biased_S5 = hydra::make_pdf(Model_2015_biased_S5, Integrator);
    auto PDF_2015_biased_S6 = hydra::make_pdf(Model_2015_biased_S6, Integrator);

    auto PDF_2016_unbiased_S1 = hydra::make_pdf(Model_2016_unbiased_S1, Integrator);
    auto PDF_2016_unbiased_S2 = hydra::make_pdf(Model_2016_unbiased_S2, Integrator);
    auto PDF_2016_unbiased_S3 = hydra::make_pdf(Model_2016_unbiased_S3, Integrator);
    auto PDF_2016_unbiased_S4 = hydra::make_pdf(Model_2016_unbiased_S4, Integrator);
    auto PDF_2016_unbiased_S5 = hydra::make_pdf(Model_2016_unbiased_S5, Integrator);
    auto PDF_2016_unbiased_S6 = hydra::make_pdf(Model_2016_unbiased_S6, Integrator);

    auto PDF_2016_biased_S1 = hydra::make_pdf(Model_2016_biased_S1, Integrator);
    auto PDF_2016_biased_S2 = hydra::make_pdf(Model_2016_biased_S2, Integrator);
    auto PDF_2016_biased_S3 = hydra::make_pdf(Model_2016_biased_S3, Integrator);
    auto PDF_2016_biased_S4 = hydra::make_pdf(Model_2016_biased_S4, Integrator);
    auto PDF_2016_biased_S5 = hydra::make_pdf(Model_2016_biased_S5, Integrator);
    auto PDF_2016_biased_S6 = hydra::make_pdf(Model_2016_biased_S6, Integrator);


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn_2015_unbiased_S1 = hydra::make_loglikehood_fcn(PDF_2015_unbiased_S1, dts_2015_unbiased_S1_d);
    auto fcn_2015_unbiased_S2 = hydra::make_loglikehood_fcn(PDF_2015_unbiased_S2, dts_2015_unbiased_S2_d);
    auto fcn_2015_unbiased_S3 = hydra::make_loglikehood_fcn(PDF_2015_unbiased_S3, dts_2015_unbiased_S3_d);
    auto fcn_2015_unbiased_S4 = hydra::make_loglikehood_fcn(PDF_2015_unbiased_S4, dts_2015_unbiased_S4_d);
    auto fcn_2015_unbiased_S5 = hydra::make_loglikehood_fcn(PDF_2015_unbiased_S5, dts_2015_unbiased_S5_d);
    auto fcn_2015_unbiased_S6 = hydra::make_loglikehood_fcn(PDF_2015_unbiased_S6, dts_2015_unbiased_S6_d);
    fcn_2015_unbiased_S1.SetFcnMaxValue(2.22507e+12);
    fcn_2015_unbiased_S2.SetFcnMaxValue(2.22507e+12);
    fcn_2015_unbiased_S3.SetFcnMaxValue(2.22507e+12);
    fcn_2015_unbiased_S4.SetFcnMaxValue(2.22507e+12);
    fcn_2015_unbiased_S5.SetFcnMaxValue(2.22507e+12);
    fcn_2015_unbiased_S6.SetFcnMaxValue(2.22507e+12);
    
    auto fcn_2015_biased_S1 = hydra::make_loglikehood_fcn(PDF_2015_biased_S1, dts_2015_biased_S1_d);
    auto fcn_2015_biased_S2 = hydra::make_loglikehood_fcn(PDF_2015_biased_S2, dts_2015_biased_S2_d);
    auto fcn_2015_biased_S3 = hydra::make_loglikehood_fcn(PDF_2015_biased_S3, dts_2015_biased_S3_d);
    auto fcn_2015_biased_S4 = hydra::make_loglikehood_fcn(PDF_2015_biased_S4, dts_2015_biased_S4_d);
    auto fcn_2015_biased_S5 = hydra::make_loglikehood_fcn(PDF_2015_biased_S5, dts_2015_biased_S5_d);
    auto fcn_2015_biased_S6 = hydra::make_loglikehood_fcn(PDF_2015_biased_S6, dts_2015_biased_S6_d);
    fcn_2015_biased_S1.SetFcnMaxValue(2.22507e+12);
    fcn_2015_biased_S2.SetFcnMaxValue(2.22507e+12);
    fcn_2015_biased_S3.SetFcnMaxValue(2.22507e+12);
    fcn_2015_biased_S4.SetFcnMaxValue(2.22507e+12);
    fcn_2015_biased_S5.SetFcnMaxValue(2.22507e+12);
    fcn_2015_biased_S6.SetFcnMaxValue(2.22507e+12);
    
    auto fcn_2016_unbiased_S1 = hydra::make_loglikehood_fcn(PDF_2016_unbiased_S1, dts_2016_unbiased_S1_d);
    auto fcn_2016_unbiased_S2 = hydra::make_loglikehood_fcn(PDF_2016_unbiased_S2, dts_2016_unbiased_S2_d);
    auto fcn_2016_unbiased_S3 = hydra::make_loglikehood_fcn(PDF_2016_unbiased_S3, dts_2016_unbiased_S3_d);
    auto fcn_2016_unbiased_S4 = hydra::make_loglikehood_fcn(PDF_2016_unbiased_S4, dts_2016_unbiased_S4_d);
    auto fcn_2016_unbiased_S5 = hydra::make_loglikehood_fcn(PDF_2016_unbiased_S5, dts_2016_unbiased_S5_d);
    auto fcn_2016_unbiased_S6 = hydra::make_loglikehood_fcn(PDF_2016_unbiased_S6, dts_2016_unbiased_S6_d);
    fcn_2016_unbiased_S1.SetFcnMaxValue(2.22507e+12);
    fcn_2016_unbiased_S2.SetFcnMaxValue(2.22507e+12);
    fcn_2016_unbiased_S3.SetFcnMaxValue(2.22507e+12);
    fcn_2016_unbiased_S4.SetFcnMaxValue(2.22507e+12);
    fcn_2016_unbiased_S5.SetFcnMaxValue(2.22507e+12);
    fcn_2016_unbiased_S6.SetFcnMaxValue(2.22507e+12);
    
    auto fcn_2016_biased_S1 = hydra::make_loglikehood_fcn(PDF_2016_biased_S1, dts_2016_biased_S1_d);
    auto fcn_2016_biased_S2 = hydra::make_loglikehood_fcn(PDF_2016_biased_S2, dts_2016_biased_S2_d);
    auto fcn_2016_biased_S3 = hydra::make_loglikehood_fcn(PDF_2016_biased_S3, dts_2016_biased_S3_d);
    auto fcn_2016_biased_S4 = hydra::make_loglikehood_fcn(PDF_2016_biased_S4, dts_2016_biased_S4_d);
    auto fcn_2016_biased_S5 = hydra::make_loglikehood_fcn(PDF_2016_biased_S5, dts_2016_biased_S5_d);
    auto fcn_2016_biased_S6 = hydra::make_loglikehood_fcn(PDF_2016_biased_S6, dts_2016_biased_S6_d);
    fcn_2016_biased_S1.SetFcnMaxValue(2.22507e+12);
    fcn_2016_biased_S2.SetFcnMaxValue(2.22507e+12);
    fcn_2016_biased_S3.SetFcnMaxValue(2.22507e+12);
    fcn_2016_biased_S4.SetFcnMaxValue(2.22507e+12);
    fcn_2016_biased_S5.SetFcnMaxValue(2.22507e+12);
    fcn_2016_biased_S6.SetFcnMaxValue(2.22507e+12);

    auto sim_fcn = hydra::make_simultaneous_fcn(fcn_2015_unbiased_S1, fcn_2015_biased_S1, fcn_2016_unbiased_S1, fcn_2016_biased_S1,
                                                fcn_2015_unbiased_S2, fcn_2015_biased_S2, fcn_2016_unbiased_S2, fcn_2016_biased_S2,
                                                fcn_2015_unbiased_S3, fcn_2015_biased_S3, fcn_2016_unbiased_S3, fcn_2016_biased_S3,
                                                fcn_2015_unbiased_S4, fcn_2015_biased_S4, fcn_2016_unbiased_S4, fcn_2016_biased_S4,
                                                fcn_2015_unbiased_S5, fcn_2015_biased_S5, fcn_2016_unbiased_S5, fcn_2016_biased_S5,
                                                fcn_2015_unbiased_S6, fcn_2015_biased_S6, fcn_2016_unbiased_S6, fcn_2016_biased_S6 );


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

#endif // FIT_SIM_PHIS_FULL_MODEL_INL_