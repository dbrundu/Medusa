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
 *  test_phis_JpsiKK.inl
 *
 *  Created on: 08/06/2020
 *      Author: Davide Brundu
 */



#ifndef TEST_PHIS_JPSIKK_INL_
#define TEST_PHIS_JPSIKK_INL_

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


// This library
#include <medusa/models/Utils.h>
#include <medusa/models/phi_s/PhisSignalOnly.h>



using namespace hydra::placeholders;
using namespace hydra::arguments;



// B0s -> J/psi  (Phi -> K+ K-)
//         |-> mu+ mu-

declarg(Jpsi,  hydra::Vector4R)
declarg(Phi,   hydra::Vector4R)
declarg(KaonP, hydra::Vector4R)
declarg(KaonM, hydra::Vector4R)
declarg(MuonP, hydra::Vector4R)
declarg(MuonM, hydra::Vector4R)

declarg(theta_h_t,  double)
declarg(theta_l_t,  double)
declarg(chi_t,      double)
declarg(dtime_t,    double)



template<typename Backend, typename Model, typename Container>
size_t generate_dataset(Backend const& system, Model const& model, Container& final_dataset, size_t nevents, size_t bunch_size);



TEST_CASE( "Phis Benchmarks")
{

    hydra::Print::SetLevel(hydra::WARNING);

    size_t  nentries   = 500000; 

    /*----------------------------/
    *  Constants
    *--------------------------- */
    const double B0_mass   = 5.336688;
    const double Jpsi_mass = 3.0969;
    const double Phi_mass  = 1.019462;
    const double K_mass    = 0.493677;
    const double pi_mass   = 0.13957061;
    const double mu_mass   = 0.1056583745;


    /*----------------------------/
    *  Model parameters
    *--------------------------- */
    const bool   B0bar      = false;

    const double A0         = ::sqrt(0.542);
    const double Aperp      = ::sqrt(0.206);
    const double AS         = ::sqrt(0.0037);

    const double phi0       = -0.082;
    const double phipar     = -0.043 + phi0;
    const double phiperp    = -0.074 + phi0;
    const double phiS       = 0.021 + phi0;

    const double lambda0    = 0.955;
    const double lambdapar  = 0.978*lambda0;
    const double lambdaperp = 1.23*lambda0;
    const double lambdaS    = 1.28*lambda0;

    const double delta0     = 0.0;
    const double deltapar   = 3.030 + delta0;
    const double deltaperp  = 2.60  + delta0;
    const double deltaS     = -0.30 + delta0;

    const double deltagammasd = -0.0044;
    const double deltagammas  = 0.0782;
    const double deltams      = 17.713;

    std::vector<double> parameters = {A0,     Aperp,     AS,         deltagammasd,
                                      deltagammas, deltams,
                                      phi0,    phipar,    phiperp,    phiS,
                                      lambda0, lambdapar, lambdaperp, lambdaS,
                                      delta0,  deltapar,  deltaperp,  deltaS};

    auto A0_p             = hydra::Parameter::Create("A0" ).Value(A0).Error(0.0001);
    auto Aperp_p          = hydra::Parameter::Create("Aperp").Value(Aperp).Error(0.0001);
    auto AS_p             = hydra::Parameter::Create("AS" ).Value(AS).Error(0.0001);

    auto DeltaGamma_sd_p  = hydra::Parameter::Create("DeltaGamma_sd" ).Value(deltagammasd).Error(0.0001);
    auto DeltaGamma_p     = hydra::Parameter::Create("DeltaGamma").Value(deltagammas).Error(0.0001);
    auto DeltaM_p         = hydra::Parameter::Create("DeltaM" ).Value(deltams).Error(0.0001);

    auto phi_0_p          = hydra::Parameter::Create("phi_0").Value(phi0).Error(0.0001);
    auto phi_par_p        = hydra::Parameter::Create("phi_par" ).Value(phipar).Error(0.0001);
    auto phi_perp_p       = hydra::Parameter::Create("phi_perp").Value(phiperp).Error(0.0001);
    auto phi_S_p          = hydra::Parameter::Create("phi_S" ).Value(phiS).Error(0.0001);

    auto lambda_0_p       = hydra::Parameter::Create("lambda_0").Value(lambda0).Error(0.0001);
    auto lambda_par_p     = hydra::Parameter::Create("lambda_par" ).Value(lambdapar).Error(0.0001);
    auto lambda_perp_p    = hydra::Parameter::Create("lambda_perp").Value(lambdaperp).Error(0.0001);
    auto lambda_S_p       = hydra::Parameter::Create("lambda_S" ).Value(lambdaS).Error(0.0001);

    auto delta_0_p        = hydra::Parameter::Create("delta_0").Value(delta0).Error(0.0001);
    auto delta_par_p      = hydra::Parameter::Create("delta_par").Value(deltapar).Error(0.0001);
    auto delta_perp_p     = hydra::Parameter::Create("delta_perp" ).Value(deltaperp).Error(0.0001);
    auto delta_S_p        = hydra::Parameter::Create("delta_S").Value(deltaS).Error(0.0001);

    hydra::Parameter hydraparams[18] = {A0_p, Aperp_p, AS_p,
                                        DeltaGamma_sd_p, DeltaGamma_p, DeltaM_p,
                                        phi_0_p,    phi_par_p,    phi_perp_p,    phi_S_p,
                                        lambda_0_p, lambda_par_p, lambda_perp_p, lambda_S_p,
                                        delta_0_p,  delta_par_p,  delta_perp_p,   delta_S_p};


    auto MODEL  = medusa::PhisSignalOnly< B0bar, dtime_t, theta_h_t, theta_l_t, chi_t>(hydraparams);



    /*------------------------------------------------------/
     *  Unweighted dataset generation
     *-----------------------------------------------------*/
    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, chi_t> , hydra::host::sys_t> dataset_h;

    generate_dataset(hydra::device::sys, MODEL, dataset_h, nentries, nentries);

    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, chi_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h , dataset_d);
    
    
    
    const double min_t     = 0.0;
    const double max_t     = 20.0;
    const double min_theta = 0.0;
    const double max_theta = PI;
    const double min_chi   = 0.0;
    const double max_chi   = 2*PI;
    const size_t N         = 50;

    hydra::Plain<4,  hydra::device::sys_t> Integrator( {min_t , min_theta, min_theta, min_chi},
                                                           {max_t , max_theta, max_theta, max_chi}, 1000000);
    
    auto Model_PDF = hydra::make_pdf( MODEL, Integrator);
    
    auto fcn0 = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
    fcn0(parameters);
    
    
    
    
    /*------------------------------------------------------/
     *  Print and plot
     *-----------------------------------------------------*/

    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset: {"<< dataset_h[i]  << "}"<< std::endl;


    #ifdef _ROOT_AVAILABLE_

    TH1D timedist("timedist","Decay Time; time (ps); Candidates / bin",100, 0, 20);
    TH1D thetahdist("thetahdist","Theta_h Angle; angle (rad); Candidates / bin",50, -1, 1);
    TH1D thetaldist("thetaldist","Theta_l Angle; angle (rad); Candidates / bin",100, -1, 1);
    TH1D chidist("chidist","Chi angle;angle (rad);Candidates / bin",50, 0, 2*PI);

    for(auto x : dataset_h){
         timedist.Fill( (double)hydra::get<0>(x) );
         thetahdist.Fill( ::cos((double)hydra::get<1>(x)) );
         thetaldist.Fill( ::cos((double)hydra::get<2>(x)) );
         chidist.Fill( (double)hydra::get<3>(x) );
    }

    TCanvas canvas("canvas","canvas",3200,800);
    canvas.Divide(4,1);

    canvas.cd(1);
    timedist.Draw();

    canvas.cd(2);
    thetahdist.Draw();

    canvas.cd(3);
    thetaldist.Draw();

    canvas.cd(4);
    chidist.Draw();

    canvas.SaveAs("test_phis_JpsiKK.pdf");

    #endif //_ROOT_AVAILABLE_




    /*------------------------------------------------------/
     *  BENCHMARKS
     *-----------------------------------------------------*/
     
    BENCHMARK( "Integration" )
    {
        return Integrator(MODEL) ; 
    };
    
    
    
    
    BENCHMARK_ADVANCED( "Evaluation + cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
    
        meter.measure([=] { return fcn(parameters); });
    };
    
    
    
    
    BENCHMARK_ADVANCED( "Evaluation + non-cached Integration" )(Catch::Benchmark::Chronometer meter)
    {
        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
        parameters[0] *= 1.01;
    
        meter.measure([=] { return fcn(parameters); });
    };
    
    
    
    
    BENCHMARK( "Cached Evaluation + cached Integration" )
    {    
        return fcn0(parameters); 
    };
    
    
    
    size_t index = 0;
    
    BENCHMARK( "Simple Functor call on 1 event" )
    {
        return MODEL( dataset_h[index++] );
    };


}





template<typename Backend, typename Model, typename Container>
size_t generate_dataset(Backend const& system, Model const& model, Container& final_dataset, size_t nevents, size_t bunch_size)
{


    const double B0_mass   = 5.336688;
    const double Jpsi_mass = 3.0969;
    const double Phi_mass  = 1.019462;
    const double K_mass    = 0.493677;
    const double pi_mass   = 0.13957061;
    const double mu_mass   = 0.1056583745;


    auto CastToVariables  = hydra::wrap_lambda(
            [] __hydra_dual__ (Jpsi jpsi, Phi phi, MuonP mup, MuonM mum, KaonP kaonp, KaonM kaonm, size_t n )
    {
        hydra_thrust::default_random_engine engine;
        hydra_thrust::uniform_real_distribution<double> uniDist(0.0, 20.0);
        engine.discard(n);
        dtime_t decay_time = uniDist(engine);

        theta_l_t theta_l    = ::acos( medusa::cos_decay_angle(jpsi + phi, phi,  kaonp) );
        theta_h_t theta_h    = ::acos( medusa::cos_decay_angle(jpsi + phi, jpsi, mup) );
        chi_t     chiangle   = medusa::chi_plane_angle(kaonm, kaonp, mup, mum);

        return hydra::make_tuple(decay_time, theta_h, theta_l, chiangle) ;

    });



    hydra::Vector4R B0(B0_mass, 0.0, 0.0, 0.0);

    hydra::PhaseSpace<2> B0_phsp(B0_mass, {Jpsi_mass, Phi_mass});

    hydra::PhaseSpace<2> Jpsi_phsp(Jpsi_mass, {mu_mass , mu_mass});

    hydra::PhaseSpace<2> Phi_phsp(Phi_mass, {K_mass , K_mass});

    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, chi_t> , Backend > dataset(bunch_size);

    hydra::SeedRNG S{};


    auto B0Decay   = hydra::Decays< hydra::tuple<Jpsi,Phi>,
                              hydra::device::sys_t >( B0_mass, {Jpsi_mass, Phi_mass }, bunch_size);

    auto JpsiDecay = hydra::Decays< hydra::tuple<MuonP,MuonM>,
                              hydra::device::sys_t >(Jpsi_mass, { mu_mass , mu_mass}, bunch_size);

    auto PhiDecay  = hydra::Decays< hydra::tuple<KaonP,KaonM>,
                              hydra::device::sys_t >(Phi_mass, { K_mass , K_mass}, bunch_size);

    size_t k = 0;

    auto start = std::chrono::high_resolution_clock::now();

    do {
        B0_phsp.SetSeed( S() );
        B0_phsp.Generate(B0, B0Decay);

        Jpsi_phsp.SetSeed( S() );
        Jpsi_phsp.Generate(B0Decay.GetDaugtherRange(_0), JpsiDecay);

        Phi_phsp.SetSeed( S() );
        Phi_phsp.Generate(B0Decay.GetDaugtherRange(_1), PhiDecay);

        auto range  = B0Decay.Meld( JpsiDecay, PhiDecay, hydra::range(0, B0Decay.size()) ) | CastToVariables;

        hydra::copy(range, dataset);

        auto dataset_unwgt = hydra::unweight( dataset, model);

        final_dataset.insert(final_dataset.end(), dataset_unwgt.begin(), dataset_unwgt.end() );

        ++k;

    } while(final_dataset.size()<nevents );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    final_dataset.erase(final_dataset.begin()+nevents, final_dataset.end());

    //output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------Generation on Device ---------"  << std::endl;
    std::cout << "| B0 -> J/psi Phi                        "  << std::endl;
    std::cout << "| Number of events :"<< nevents             << std::endl;
    std::cout << "| Number of events gen:"<< k*bunch_size     << std::endl;
    std::cout << "| Time (ms)        :"<< elapsed.count()     << std::endl;
    std::cout << "-----------------------------------------"  << std::endl;

    return final_dataset.size();


}

#endif