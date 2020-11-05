/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto
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
/*
 *  test_D2hhmumu.inl
 *
 *  Created on: 18/10/2020
 *      Author: Davide Brundu
 */



#ifndef TEST_D2HHMUMU_INL_
#define TEST_D2HHMUMU_INL_

#define CATCH_CONFIG_ENABLE_BENCHMARKING

#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>


#define DEBUG(message)\
    std::cout<< "\033[1;33mDEBUG: \033[0m" << message << "\n";\

/*---------------------------------
 * std
 * ---------------------------------
 */
#include <iostream>
#include <chrono>


/*---------------------------------
 * command line arguments
 *---------------------------------
 */
#include <tclap/CmdLine.h>

/*---------------------------------
 * Include hydra classes and
 * algorithms for
 *--------------------------------
 */
#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/PhaseSpace.h>
#include <hydra/Decays.h>
#include <hydra/Function.h>
#include <hydra/Lambda.h>
#include <hydra/multivector.h>
#include <hydra/Random.h>
#include <hydra/SeedRNG.h>
#include <hydra/Plain.h>
#include <hydra/VegasState.h>
#include <hydra/Vegas.h>
#include <hydra/Sobol.h>
#include <hydra/LogLikelihoodFCN.h>

/*-------------------------------------
 * Include classes from ROOT to fill
 * and draw histograms and plots.
 *-------------------------------------
 */
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

#include <medusa/models/Utils.h>
#include <medusa/models/D2hhmumu/D2hhmumuAngularDist.h>



using namespace hydra::placeholders;
using namespace hydra::arguments;



// D0 -> pi+ pi- mu+ mu-
// 


declarg(PionP, hydra::Vector4R)
declarg(PionM, hydra::Vector4R)
declarg(MuonP, hydra::Vector4R)
declarg(MuonM, hydra::Vector4R)

declarg(theta_l_t,  double)
declarg(chi_t,      double)




template<typename Backend, typename Model, typename Container>
size_t generate_dataset(Backend const& system, Model const& model, Container& final_dataset, size_t nevents, size_t bunch_size);



TEST_CASE( "D0->hhmumu Benchmarks")
{

    hydra::Print::SetLevel(hydra::WARNING);

    size_t  nentries   = 1000000; 
    bool do_benchmark  = true;
    
    /*----------------------------/
    *  Constants and Histograms
    *--------------------------- */

    const double D0_mass    = 1864.83;
    const double pi_mass    = 139.57061;
    const double mu_mass    = 105.6583745 ;
    
    #ifdef _ROOT_AVAILABLE_
    TH1D thetaldist("thetaldist","Theta_l Angle; angle (rad); Candidates / bin", 150, -1.0, 1.0);
    TH1D chidist("chidist","Chi angle;angle (rad);Candidates / bin", 150, 0.0, 2*PI);
    
    TH1D thetaldist_w("thetaldist_w","Theta_l Angle; angle (rad); Candidates / bin", 150, -1.0, 1.0);
    TH1D chidist_w("chidist_w","Chi angle;angle (rad);Candidates / bin", 150, 0.0, 2*PI);
    
    TH1D thetaldist_s("thetaldist_s","Theta_l Angle; angle (rad); Candidates / bin", 150, -1.0, 1.0);
    TH1D chidist_s("chidist_s","Chi angle;angle (rad);Candidates / bin", 150, 0.0, 2*PI);
    #endif


    /*----------------------------/
    * Angular Model
    *--------------------------- */
    const double H2 = 0.03; 
    const double H3 = -0.05; 
    const double H4 = 0.0 ;
    const double H5 = 0.0 ; 
    const double H6 = 0.0 ; 
    const double H7 = 0.0 ; 
    const double H8 = 0.0 ;
    const double H9 = 0.0 ;
    
    std::vector<double> parameters = {H2, H3, H4, H5, H6, H7, H8, H9};
    
    auto H2_p = hydra::Parameter::Create("H2").Value(H2).Error(0.0001);
    auto H3_p = hydra::Parameter::Create("H3").Value(H3).Error(0.0001);
    auto H4_p = hydra::Parameter::Create("H4").Value(H4).Error(0.0001);
    auto H5_p = hydra::Parameter::Create("H5").Value(H5).Error(0.0001);
    auto H6_p = hydra::Parameter::Create("H6").Value(H6).Error(0.0001);
    auto H7_p = hydra::Parameter::Create("H7").Value(H7).Error(0.0001);
    auto H8_p = hydra::Parameter::Create("H8").Value(H8).Error(0.0001);
    auto H9_p = hydra::Parameter::Create("H9").Value(H9).Error(0.0001);

    hydra::Parameter hydraparams[8] = {H2_p, H3_p, H4_p, H5_p, H6_p, H7_p, H8_p, H9_p};
    
    auto MODEL  = medusa::D2hhmumuAngularDist<theta_l_t, chi_t>(hydraparams);



    /*------------------------------------------------------/
     *  Unweighted dataset generation
     *-----------------------------------------------------*/
    hydra::multivector< hydra::tuple< theta_l_t, chi_t> , hydra::host::sys_t > dataset_h;

    generate_dataset(hydra::device::sys, MODEL, dataset_h, nentries, nentries);

    hydra::multivector<hydra::tuple<theta_l_t, chi_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
   
    hydra::copy(dataset_h , dataset_d);
    
    
    
    /*------------------------------------------------------/
     *  Print and plot
     *-----------------------------------------------------*/
    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset: {"<< dataset_h[i]  << "}"<< std::endl;

    #ifdef _ROOT_AVAILABLE_

    for(auto x : dataset_h){
         double ch = hydra::get<1>(x);
         double th = hydra::get<0>(x);
         thetaldist.Fill( ::cos(th) );
         chidist.Fill( ch );
    }
    
    TCanvas canvas("canvas","canvas",1600,800);
    canvas.Divide(2,1);
    canvas.cd(1);
    thetaldist.Draw();
    canvas.cd(2);
    chidist.Draw();
    canvas.SaveAs("test_D2hhmumu.pdf");

    #endif //_ROOT_AVAILABLE_



    /*------------------------------------------------------/
     *  Preparing the benchmark
     *-----------------------------------------------------*/
    const double min_theta = 0.0;
    const double max_theta = PI;
    const double min_chi   = 0.0;
    const double max_chi   = 2*PI;
    const size_t N         = 50;
        
    hydra::Plain<2,  hydra::device::sys_t> Integrator( {min_theta, min_chi}, {max_theta, max_chi}, 1000000);
    
    auto Model_PDF = hydra::make_pdf( MODEL, Integrator);
    
    auto fcn0 = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
    fcn0(parameters);




    /*------------------------------------------------------/
     * Benchmark for functor normalization
     *-----------------------------------------------------*/
    BENCHMARK( "Integration" )
    {
        return Integrator(MODEL) ; 
    };
    
    
     

    /*------------------------------------------------------/
     * Benchmark for fcn evaluation with cached integration
     *   This can be done because the PDF object is already
     *   constructed, so it has its cached normalization, while the
     *   FCN object is recreated each time, 
     *   thus with non-cached fcn value.
     *-----------------------------------------------------*/
    BENCHMARK_ADVANCED( "Evaluation + cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
    
        meter.measure([=] { return fcn(parameters); });
    };
    
    
    
    /*------------------------------------------------------/
     * Benchmark for fcn evaluation with non-cached integration
     *   This can be done because the 
     *   FCN object is recreated each time and the parameters
     *   are modified, in order to trigger the normalization
     *   in the PDF object.
     *-----------------------------------------------------*/
    BENCHMARK_ADVANCED( "Evaluation + non-cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
        for(auto& p : parameters ) p *= 1.01;
    
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
    auto rng = hydra::detail::RndUniform<size_t , hydra::default_random_engine >(S(), 0, dataset_h.size()-1);
    size_t index=0;
    
    BENCHMARK_ADVANCED( "Simple Functor call on 1 event" )(Catch::Benchmark::Chronometer meter)
    {
        const size_t i = rng(index++);
        auto x = dataset_h[i];
        meter.measure( [=] { return MODEL( x ); });
    };


}







template<typename Backend, typename Model, typename Container>
size_t generate_dataset(Backend const& system, Model const& model, Container& final, size_t nevents, size_t bunch_size)
{

    const double D0_mass    = 1864.83;
    const double pi_mass    = 139.57061;
    const double mu_mass    = 105.6583745 ;
    
    hydra::SeedRNG S{};
    
    auto CastToVariables  = hydra::wrap_lambda( [] __hydra_dual__ (PionP pip, PionM pim, MuonP mup, MuonM mum)
    {
        theta_l_t theta_l    = ::acos( medusa::cos_decay_angle(pip+pim+mup+mum, mup+mum,  mup) );
        chi_t     chiangle   = medusa::chi_plane_angle(pim, pip, mup, mum);
        return hydra::make_tuple( theta_l, chiangle) ;
    });
    
    
    auto model_for_rew = hydra::wrap_lambda(	[ model ] __hydra_dual__ ( PionP pip, PionM pim, MuonP mup, MuonM mum){
        theta_l_t theta_l    = ::acos( medusa::cos_decay_angle(pip+pim+mup+mum, mup+mum,  mup) );
        chi_t     chiangle   = medusa::chi_plane_angle(pim, pip, mup, mum);
        return model( theta_l, chiangle);
    } );
    
    
    hydra::Vector4R D0p(D0_mass, 0.0, 0.0, 0.0);
    const double ph_masses[4]{pi_mass, pi_mass, mu_mass, mu_mass };

    // Create phase space object D0 -> pi pi mu mu
    hydra::PhaseSpace<4> phsp(D0_mass, ph_masses);

    hydra::multivector< hydra::tuple< theta_l_t, chi_t> , Backend > dataset(bunch_size);

    auto start = std::chrono::high_resolution_clock::now();


    auto _data = hydra::Decays< hydra::tuple<PionP,PionM, MuonP, MuonM>, hydra::device::sys_t >( D0_mass, {pi_mass, pi_mass, mu_mass, mu_mass }, bunch_size);

    do {
        phsp.SetSeed( S() );
        phsp.Generate(D0p, _data );
        
        auto vars_unwgt = hydra::unweight(hydra::device::sys, _data, _data.GetEventWeightFunctor(model_for_rew), -1.0 , S() ) | CastToVariables ;
        
        hydra::multivector< hydra::tuple< theta_l_t, chi_t> , Backend > dataset_unwgt( vars_unwgt.size() );
        hydra::copy(vars_unwgt, dataset_unwgt);
        
        final.insert(final.size()==0? final.begin():final.end(), dataset_unwgt.begin(), dataset_unwgt.end() );
        
     } while(final.size()<=nevents );

    final.erase(final.begin()+nevents, final.end());
     
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    //output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-------------Generation on Device ---------"  << std::endl;
    std::cout << "| D0 -> pi pi mu mu                        "  << std::endl;
    std::cout << "| Number of events :"<< nevents             << std::endl;
    std::cout << "| Time (ms)        :"<< elapsed.count()     << std::endl;
    std::cout << "-------------------------------------------"  << std::endl;

    return final.size();


}

#endif
