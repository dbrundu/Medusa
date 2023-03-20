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
 *  Created on: 08/06/2020
 *      Author: Davide Brundu
 *      Updated by Alessandro Maria Ricci in 04/07/2021
 * 
 *  Complete benchmarks for fit_sim_phis_full_model.inl
 *---------------------------------------------------------------------------*/



#ifndef BENCHMARKS_COMPLETE_SIM_PHIS_FULL_MODEL_INL_
#define BENCHMARKS_COMPLETE_SIM_PHIS_FULL_MODEL_INL_

#define CATCH_CONFIG_ENABLE_BENCHMARKING


#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>


#define DEBUG(message)\
    std::cout<< "\033[1;33mDEBUG: \033[0m" << message << "\n";\



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


//-------------------------------------------
//              Benchmarks
//-------------------------------------------

TEST_CASE( "Benchmarks for B0s -> J/psi Phi -> mu+ mu- K+ K-")
{

    hydra::Print::SetLevel(hydra::WARNING);

    size_t  nentries   = 1000000;


    //---------------------------------
    //      model generation
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
    //      PDF and FCN generation
    //---------------------------------
    
    // Integrator (it always returns the value 1.0, because the normalization is computed
    // in FullAnalyticPhis.h. This choice is justified by the fact that Hydra does not support
    // a normalization factor which depends from the experimental variables)
    auto integrator = hydra::AnalyticalIntegral< medusa::FullAnalyticPhis<CubicSpline,
                                                            dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto Model_PDF = hydra::make_pdf(Model_2015_unbiased_S1, integrator);
    
    auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dts_2015_unbiased_S1_d);



    /*----------------------------------------------------------/
     * Benchmark for fcn evaluation with cached integration
     *---------------------------------------------------------*/
    BENCHMARK_ADVANCED( "Evaluation + cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        meter.measure([=] { return fcn(parameters); });
    };



    /*------------------------------------------------------------/
     * Benchmark for fcn evaluation with non-cached integration
     *-----------------------------------------------------------*/
    BENCHMARK_ADVANCED( "Evaluation + non-cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        // distortion of the parameters to avoid the cached integration
        for(size_t i=0; i<parameters.size(); i++) { parameters[i] *= 1.001; }

        meter.measure([=] { return fcn(parameters); });
    };



    /*--------------------------------------------------------/
     * Benchmark for fcn evaluation with all values cached
     *   This can be done because the FCN object is already
     *   created and evaluated previously, thus
     *   all the values are cached
     *--------------------------------------------------------*/
    BENCHMARK( "Cached Evaluation + cached Integration" )
    {    
        return fcn(parameters);
    };



    /*------------------------------------------------------/
     * Benchmark for direct functor call on 1 event
     *-----------------------------------------------------*/
    hydra::SeedRNG S{};
    auto rng = hydra::detail::RndUniform<size_t , hydra::default_random_engine >(S(), 0, dts_2015_unbiased_S1_d.size()-1);
    size_t index=0;
    
    BENCHMARK_ADVANCED( "Simple Functor call on 1 event" )(Catch::Benchmark::Chronometer meter)
    {
        const size_t i = rng(index++);
        auto x = dts_2015_unbiased_S1_d[i];

        meter.measure( [=] { return Model_2015_unbiased_S1(x); });
    };


} // TEST_CASE

#endif // BENCHMARKS_COMPLETE_SIM_PHIS_FULL_MODEL_INL_
