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
/*----------------------------------------------
 * benchmarks2_B0s_Jpsi_phi.inl
 *
 *  Created: 18/08/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  Other benchmarks for fit_B0s_Jpsi_phi.inl
 *----------------------------------------------*/

#ifndef BENCHMARKS2_B0S_JPSI_PHI_INL_
#define BENCHMARKS2_B0S_JPSI_PHI_INL_


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
#include <hydra/Vegas.h>

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
#include <medusa/phi_s/PhisSignal.h>
#include <medusa/phi_s/GenerateDataset.h>



//---------------------------------
//        Main program
//---------------------------------

int main(int argv, char** argc)
{


    //------------------------------------------------------
    //   TCLAP: read the arguments from the command line
    //------------------------------------------------------

    size_t nentries = 0;
    size_t calls = 0;
	size_t iterations = 0;
	double max_error = 0;

	try {

		TCLAP::CmdLine cmd("Command line arguments for number of events and Vegas integrator", '=');

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events", false, 5e5, "size_t");
        cmd.add(EArg);

		TCLAP::ValueArg<size_t> NCallsArg("c", "number-of-calls", "Number of call", false, 1e4, "size_t");
		cmd.add(NCallsArg);

		TCLAP::ValueArg<double> MaxErrorArg("e", "max-error", "Maximum error", false, 1.0e-5, "double");
		cmd.add(MaxErrorArg);

		TCLAP::ValueArg<size_t> IterationsArg("i", "max-iterations", "Maximum number of iterations",false, 150, "size_t");
		cmd.add(IterationsArg);

		// Parse the argv array.
		cmd.parse(argv, argc);

		// Get the value parsed by each arg.
        nentries = EArg.getValue();
        calls = NCallsArg.getValue();
		iterations = IterationsArg.getValue();
		max_error  = MaxErrorArg.getValue();

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

    // model parameters
    const bool   B0sbar     = false;    // boolean to specify wether B0s is B0sbar or not

    const double A0         = ::sqrt(0.542);
    const double Aperp      = ::sqrt(0.206);
    const double AS         = ::sqrt(0.0037);

    const double phi0       = -0.082;
    const double phipar     = -0.125;            // -0.043 + phi0
    const double phiperp    = -0.156;            // -0.074 + phi0
    const double phiS       = -0.061;            //  0.021 + phi0

    const double lambda0    = 0.955; 
    const double lambdapar  = 0.93399;           // 0.978*lambda0
    const double lambdaperp = 1.17465;           // 1.23*lambda0
    const double lambdaS    = 1.2224;            // 1.28*lambda0

    const double delta0     = 0.0;
    const double deltapar   = 3.030;            // + delta0
    const double deltaperp  = 2.60;             // + delta0
    const double deltaS     = -0.30;            // + delta0

    const double deltagammasd = -0.0044;
    const double deltagammas  = 0.0782;
    const double deltams      = 17.713;

    std::vector<double> parameters = {A0, Aperp, AS, 
                                      deltagammasd, deltagammas, deltams,
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

    auto Model = medusa::PhisSignal<B0sbar, dtime_t, theta_h_t, theta_l_t, phi_t>(hydraparams);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t> , hydra::host::sys_t> dataset_h;

    generate_dataset(Model, dataset_h, nentries, nentries);
    
    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h, dataset_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------
    
    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset_h: {"<< dataset_h[i]  << "}"<< std::endl;


    #ifdef _ROOT_AVAILABLE_

        TH1D timedist("timedist","Decay Time; time (ps); Candidates / bin",100, 0, 20);
        TH1D thetahdist("thetahdist","Theta_h Angle; angle (rad); Candidates / bin",50, -1, 1);
        TH1D thetaldist("thetaldist","Theta_l Angle; angle (rad); Candidates / bin",100, -1, 1);
        TH1D phidist("phidist","Phi angle; angle (rad); Candidates / bin",50, 0, 2*PI);

        for(auto x : dataset_h)
        {
            timedist.Fill( (double)hydra::get<0>(x) );
            thetahdist.Fill( ::cos((double)hydra::get<1>(x)) );
            thetaldist.Fill( ::cos((double)hydra::get<2>(x)) );
            phidist.Fill( (double)hydra::get<3>(x) );
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
        phidist.Draw();

        canvas.SaveAs("Dataset_B0s_JpsiPhi.pdf");

    #endif //_ROOT_AVAILABLE_


    //---------------------------------
    //          PDF generation
    //---------------------------------

    const dtime_t min_t         = 0.0;
    const dtime_t max_t         = 20.0;
    const theta_h_t min_theta_h = 0.0;
    const theta_h_t max_theta_h = PI;
    const theta_l_t min_theta_l = 0.0;
    const theta_l_t max_theta_l = PI;
    const phi_t min_phi         = 0.0;
    const phi_t max_phi         = 2*PI;
    const size_t N              = 4;

    // Vegas State_d holds the resources for performing the integration
    hydra::VegasState<N, hydra::device::sys_t> State_d({min_t, min_theta_h, min_theta_l, min_phi},
                                                            {max_t, max_theta_h, max_theta_l, max_phi});

    State_d.SetVerbose(-2);
    State_d.SetAlpha(1.5);
    State_d.SetIterations( iterations );
    State_d.SetUseRelativeError(1);
    State_d.SetMaxError( max_error );
    State_d.SetCalls( calls );
    State_d.SetTrainingCalls( calls/10 );
    State_d.SetTrainingIterations(2);

    // Vegas integrator
    hydra::Vegas<N, hydra::device::sys_t> Vegas_d(State_d);

    // make PDF
    auto model_PDF = hydra::make_pdf(Model, Vegas_d);

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

    // print functor evaluation
    auto x = dataset_d[0];

    auto start_functor = std::chrono::high_resolution_clock::now();
    Model(x);
    auto stop_functor = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_functor = stop_functor - start_functor;

    std::cout << "Functor = " << Model(x) << std::endl;
    std::cout << "Time (ms) = " << elapsed_functor.count() << std::endl;


    // print Vegas integration
    auto start_vegas = std::chrono::high_resolution_clock::now();
    Vegas_d.Integrate(Model);
    auto stop_vegas = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_vegas = stop_vegas - start_vegas;

    std::cout << "Integral = " << Vegas_d.GetState().GetResult() << " +/- " << Vegas_d.GetState().GetSigma() << std::endl; 
    std::cout << "Time (ms) = " << elapsed_vegas.count() << std::endl;


    // print FCN evaluation
    auto start_fcn = std::chrono::high_resolution_clock::now();
    fcn(parameters);
    auto stop_fcn = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_fcn = stop_fcn - start_fcn;

    std::cout << "fcn = " << fcn(parameters) << std::endl;
    std::cout << "Time (ms) = " << elapsed_fcn.count() << std::endl;

    return 0;

} // main

#endif // BENCHMARKS2_B0S_JPSI_PHI_INL_