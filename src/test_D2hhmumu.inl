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



int main(int argv, char** argc)
{

    hydra::Print::SetLevel(hydra::WARNING);

    size_t  nentries   = 0; // number of events to generate, to be get from command line
    bool do_benchmark  = true;


    try {

        TCLAP::CmdLine cmd("Command line arguments for PHSP B0s -> J/psi Phi", '=');

        TCLAP::ValueArg<size_t> NArg("n", "nevents",
                "Number of events to generate. Default is [ 10e6 ].",
                true, 10e6, "unsigned long");
        cmd.add(NArg);

        TCLAP::ValueArg<bool> DoBenchmark("b", "benchmark",
                "Choose wether to do the benchmark or not. Default is true.",
                true, true, "bool");
        cmd.add(DoBenchmark);

        // Parse the argv array.
        cmd.parse(argv, argc);

        // Get the value parsed by each arg.
        nentries       = NArg.getValue();
        do_benchmark   = DoBenchmark.getValue();

    }
    catch (TCLAP::ArgException &e)  {
        std::cerr << "error: " << e.error() << " for arg " << e.argId()
                                                                << std::endl;
    }

    const double D0_mass    = 1864.83;
    const double pi_mass    = 139.57061;
    const double mu_mass    = 105.6583745 ;
    
#ifdef _ROOT_AVAILABLE_

    TH1D thetaldist("thetaldist","Theta_l Angle; angle (rad); Candidates / bin", 150, -1.0, 1.0);
    TH1D chidist("chidist","Chi angle;angle (rad);Candidates / bin", 150, 0.0, 2*PI);
    
    TH1D thetaldist_w("thetaldist","Theta_l Angle; angle (rad); Candidates / bin", 150, -1.0, 1.0);
    TH1D chidist_w("chidist","Chi angle;angle (rad);Candidates / bin", 150, 0.0, 2*PI);
    
#endif

    /*----------------------------/
    * Angular Model
    *--------------------------- */
    const double H2 = 0.03; // 0.03;
    const double H3 = -0.05; //-0.05 ;
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
    
    
    // Doing the reweight
    {
        hydra::Vector4R Parent(D0_mass, 0.0, 0.0, 0.0);

        double masses[4]{pi_mass, pi_mass, mu_mass, mu_mass };
        
        hydra::PhaseSpace<4> phsp{D0_mass, masses};
        
        auto model_for_rew = hydra::wrap_lambda(	[ MODEL ] __hydra_dual__ ( PionP pip, PionM pim, MuonP mup, MuonM mum){
                theta_l_t theta_l    = ::acos( medusa::cos_decay_angle(pip+pim+mup+mum, mup+mum,  mup) );
                chi_t     chiangle   = medusa::chi_plane_angle(pim, pip, mup, mum);
                return MODEL( theta_l, chiangle);
        } );
        
        auto CastToTheta_l  = hydra::wrap_lambda( [] __hydra_dual__ (PionP pip, PionM pim, MuonP mup, MuonM mum){
                theta_l_t costheta_l    = medusa::cos_decay_angle(pip+pim+mup+mum, mup+mum,  mup);
                return costheta_l;
            });
            
        auto CastToPhi  = hydra::wrap_lambda( [] __hydra_dual__ (PionP pip, PionM pim, MuonP mup, MuonM mum){
                chi_t     chiangle   = medusa::chi_plane_angle(pim, pip, mup, mum);
                return chiangle;
            });
            
        hydra::Decays< hydra::tuple<PionP,PionM,MuonP,MuonM>, hydra::device::sys_t > Events(D0_mass, masses, nentries);
        
        phsp.Generate(Parent, Events);
        
		auto theta_variables = Events | CastToTheta_l ;

		auto chi_variables   = Events | CastToPhi ;

		auto weights   = Events | Events.GetEventWeightFunctor(model_for_rew);
		
	
		auto Hist_theta = make_dense_histogram( hydra::device::sys, 150, -1.0, +1.0,	theta_variables,	weights);
		auto Hist_chi  = make_dense_histogram( hydra::device::sys, 150, 0.0, 2*PI,	chi_variables,	weights);
		
		
        for(size_t i=0; i< 150; i++) { 
            thetaldist_w.SetBinContent(i+1, Hist_theta.GetBinContent(i) );
            chidist_w.SetBinContent(i+1, Hist_chi.GetBinContent(i) );
        }
        
    
    }






    if(do_benchmark){

        const double min_theta = 0.0;
        const double max_theta = PI;
        const double min_chi   = 0.0;
        const double max_chi   = 2*PI;
        const size_t N         = 50;

        hydra::Plain<2,  hydra::device::sys_t> Integrator( {min_theta, min_chi},
                                                           {max_theta, max_chi}, 1000000);


        /*------------------------------------------------------/
         *  Benchmark of Integration
         *-----------------------------------------------------*/
        double sumtime = 0.0;
        for(size_t k=0 ; k<N; ++k)
        {
            auto start = std::chrono::high_resolution_clock::now();

            auto Model_PDF = hydra::make_pdf( MODEL, Integrator);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end - start;
            sumtime += elapsed.count();

        }

        std::cout << "| Integration - Time (ms)                        :"<< sumtime/N   << std::endl;



        /*------------------------------------------------------/
         *  Benchmark of Evaluation + cached integration
         *-----------------------------------------------------*/
       auto Model_PDF = hydra::make_pdf( MODEL, Integrator);

        sumtime = 0.0;

        for(size_t k=0 ; k<N; ++k)
        {
             auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);

             auto start = std::chrono::high_resolution_clock::now();

             fcn(parameters);

             auto end = std::chrono::high_resolution_clock::now();
             std::chrono::duration<double, std::milli> elapsed = end - start;
            sumtime += elapsed.count();

        }

        std::cout << "| Eval + Integration (cached) - Time (ms)        :"<< sumtime/N   << std::endl;



        /*------------------------------------------------------/
         *  Benchmark of Evaluation + Integration (not cached)
         *-----------------------------------------------------*/
        sumtime = 0.0;

        for(size_t k=0 ; k<N; ++k)
        {
            auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
            parameters[0] *= 1.01;

            auto start = std::chrono::high_resolution_clock::now();

            fcn(parameters);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end - start;
            sumtime += elapsed.count();

        }

        std::cout << "| Eval + Integration - Time (ms)                 :"<< sumtime/N   << std::endl;



        /*------------------------------------------------------/
         *  Benchmark of cached evaluation + cached integration
         *-----------------------------------------------------*/

        auto fcn = hydra::make_loglikehood_fcn(Model_PDF, dataset_d);
        fcn(parameters);

        auto start = std::chrono::high_resolution_clock::now();
        for(size_t k=0 ; k<N; ++k)  fcn(parameters);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;

        std::cout << "| Eval (Cached) - Time (ms)                      :"<< elapsed.count()/N   << std::endl << std::endl;

    } //end benchmark


    /*------------------------------------------------------/
     *  Print and plot
     *-----------------------------------------------------*/

    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset: {"<< dataset_h[i]  << "}"<< std::endl;


#ifdef _ROOT_AVAILABLE_

    TApplication *m_app = new TApplication("myapp",0,0);


    for(auto x : dataset_h){
         thetaldist.Fill( ::cos((double)hydra::get<0>(x)) );
         chidist.Fill( (double)hydra::get<1>(x) );
    }

    TCanvas canvas("canvas","canvas",1600,800);
    canvas.Divide(2,1);
    canvas.cd(1);
    thetaldist.Draw();
    canvas.cd(2);
    chidist.Draw();
    canvas.SaveAs("test_D2hhmumu.pdf");
    
    TCanvas canvas_w("canvas_w","canvas_w",1600,800);
    canvas_w.Divide(2,1);
    canvas_w.cd(1);
    thetaldist_w.Draw();
    canvas_w.cd(2);
    chidist_w.Draw();
    canvas_w.SaveAs("test_D2hhmumu_w.pdf");

    m_app->Run();

#endif //_ROOT_AVAILABLE_

    return 0;
}





template<typename Backend, typename Model, typename Container>
size_t generate_dataset(Backend const& system, Model const& model, Container& final, size_t nevents, size_t bunch_size)
{

    const double D0_mass    = 1864.83;
    const double pi_mass    = 139.57061;
    const double mu_mass    = 105.6583745 ;
    
    hydra::SeedRNG S{};
    
    auto CastToVariables  = hydra::wrap_lambda(
            [] __hydra_dual__ (PionP pip, PionM pim, MuonP mup, MuonM mum)
    {
        theta_l_t theta_l    = ::acos( medusa::cos_decay_angle(pip+pim+mup+mum, mup+mum,  mup) );
        chi_t     chiangle   = medusa::chi_plane_angle(pim, pip, mup, mum);
        return hydra::make_tuple( theta_l, chiangle) ;
    });
    
    
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
        
        auto range  = _data | CastToVariables;

        hydra::copy(range, dataset);

        auto dataset_unwgt = hydra::unweight( dataset, model );
        
        //auto dataset_2 = _data.Unweight() | CastToVariables ;
        //hydra::multivector< hydra::tuple< theta_l_t, chi_t> , Backend > dataset_unwgt(bunch_size);
        //hydra::copy(dataset_2, dataset_unwgt);
            
        final.insert(final.size()==0? final.begin():final.end(), dataset_unwgt.begin(), dataset_unwgt.end() );
        
     } while(final.size()<=nevents );

    final.erase(final.begin()+nevents, final.end());
     
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    //output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----------------- Device ----------------"  << std::endl;
    std::cout << "| D0 -> pi pi mu mu                        "  << std::endl;
    std::cout << "| Number of events :"<< nevents             << std::endl;
    std::cout << "| Time (ms)        :"<< elapsed.count()     << std::endl;
    std::cout << "-----------------------------------------"  << std::endl;

    return final.size();


}

#endif
