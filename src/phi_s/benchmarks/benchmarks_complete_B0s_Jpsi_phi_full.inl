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
 *  Complete benchmarks for fit_B0s_Jpsi_phi_full.inl
 *----------------------------------------------------------------*/



#ifndef BENCHMARKS_COMPLETE_B0S_JPSI_PHI_INL_
#define BENCHMARKS_COMPLETE_B0S_JPSI_PHI_INL_

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
    //      PDF and FCN generation
    //---------------------------------
    
    // Integrator (it always returns the value 1.0, because the normalization is computed
    // in FullAnalyticPhis.h. This choice is justified by the fact that Hydra does not support
    // a normalization factor which depends from the experimental variables)
    auto integrator = hydra::AnalyticalIntegral< medusa::FullAnalyticPhis<CubicSpline,
                                                            dtime_t, costheta_h_t, costheta_l_t, phi_t, qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto Model_PDF = hydra::make_pdf(Model, integrator);
    
    auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);



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
    auto rng = hydra::detail::RndUniform<size_t , hydra::default_random_engine >(S(), 0, dataset_d.size()-1);
    size_t index=0;
    
    BENCHMARK_ADVANCED( "Simple Functor call on 1 event" )(Catch::Benchmark::Chronometer meter)
    {
        const size_t i = rng(index++);
        auto x = dataset_d[i];

        meter.measure( [=] { return Model(x); });
    };


} // TEST_CASE

#endif // BENCHMARKS_COMPLETE_B0S_JPSI_PHI_INL_
