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
/*----------------------------------------------------------------
 *  benchmarks_B0s_Jpsi_phi.inl
 *
 *  Created on: 08/06/2020
 *      Author: Davide Brundu
 *      Updated by Alessandro Maria Ricci in 04/07/2021
 * 
 *  Complete benchmarks for fit_B0s_Jpsi_phi_signal.inl
 *----------------------------------------------------------------*/



#ifndef BENCHMARKS_COMPLETE_B0S_JPSI_PHI_SIGNAL_INL_
#define BENCHMARKS_COMPLETE_B0S_JPSI_PHI_SIGNAL_INL_

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
#include <medusa/phi_s/phis_signal/PhisSignal.h>
#include <medusa/phi_s/phis_signal/GenerateDataset.h>



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

    auto MODEL = medusa::PhisSignal<B0sbar, dtime_t, costheta_h_t, costheta_l_t, phi_t>(ModelParams);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t> , hydra::host::sys_t> dataset_h;

    medusa::GenerateDataset_SignalOnly(MODEL, dataset_h, nentries, nentries, LowerLimit, UpperLimit);
    
    hydra::multivector<hydra::tuple<dtime_t, costheta_h_t, costheta_l_t, phi_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h, dataset_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------
    
    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset_h: {"<< dataset_h[i]  << "}"<< std::endl;


    #ifdef _ROOT_AVAILABLE_

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

        TCanvas canvas("canvas","canvas",3200,800);
        canvas.Divide(4,1);

        canvas.cd(1);
        gPad->SetLogy(1);
        timedist.Draw();

        canvas.cd(2);
        thetahdist.Draw();

        canvas.cd(3);
        thetaldist.Draw();

        canvas.cd(4);
        phidist.Draw();

        canvas.SaveAs("Dataset_B0s_JpsiPhi.pdf");

    #endif //_ROOT_AVAILABLE_


    //---------------------------------
    //      PDF and FCN generation
    //---------------------------------

    auto Integrator = hydra::AnalyticalIntegral<
                                medusa::PhisSignal< B0sbar, dtime_t, costheta_h_t, costheta_l_t, phi_t> >(LowerLimit, UpperLimit);
    
    auto Model_PDF = hydra::make_pdf(MODEL, Integrator);
    
    auto fcn0 = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
    fcn0(parameters);



    /*------------------------------------------------------/
     * Benchmark for functor normalization
     *-----------------------------------------------------*/
    BENCHMARK( "Integration" )
    {
        return Integrator(MODEL) ; 
    };



    /*----------------------------------------------------------/
     * Benchmark for fcn evaluation with cached integration
     *   This can be done because the PDF object is already
     *   constructed, so it has its cached normalization,
     *   while the FCN object is recreated each time, 
     *   thus with non-cached fcn value.
     *---------------------------------------------------------*/
    BENCHMARK_ADVANCED( "Evaluation + cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
    
        meter.measure([=] { return fcn(parameters); });
    };



    /*------------------------------------------------------------/
     * Benchmark for fcn evaluation with non-cached integration
     *   This can be done because the FCN object is recreated 
     *   each time and the parameters are modified,
     *   in order to trigger the normalization in the PDF object.
     *-----------------------------------------------------------*/
    BENCHMARK_ADVANCED( "Evaluation + non-cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
        
        // distortion of the parameters to avoid the cached integration
        for(size_t i=0; i<parameters.size(); i++) { parameters[i] *= 1.001; }

        meter.measure([=] { return fcn(parameters); });
    };



    /*------------------------------------------------------/
     * Benchmark for fcn evaluation with all values cached
     *   This can be done because the 
     *   FCN object is already created and evaluated outside,
     *   thus all the values are cached
     *-----------------------------------------------------*/
    BENCHMARK( "Cached Evaluation + cached Integration" )
    {    
        return fcn0(parameters);
    };



    /*------------------------------------------------------/
     * Benchmark for direct functor call on 1 event
     *-----------------------------------------------------*/
    hydra::SeedRNG S{};
    auto rng = hydra::detail::RndUniform<size_t , hydra::default_random_engine >(S(), 0, dataset_d.size()-1);
    size_t index=0;
    
    BENCHMARK_ADVANCED( "Simple Functor call on 1 event" )(Catch::Benchmark::Chronometer meter)
    {
        const size_t i = rng(index++);
        auto x = dataset_d[i];

        meter.measure( [=] { return MODEL( x ); });
    };


} // TEST_CASE

#endif // BENCHMARKS_COMPLETE_B0S_JPSI_PHI_SIGNAL_INL_
