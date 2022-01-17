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
 * benchmarks2_B0s_Jpsi_phi.inl
 *
 *  Created: 18/08/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  Simple benchmarks for fit_B0s_Jpsi_phi_full.inl
 *---------------------------------------------------------------------------*/

#ifndef BENCHMARKS_SIMPLE_B0S_JPSI_PHI_FULL_INL_
#define BENCHMARKS_SIMPLE_B0S_JPSI_PHI_FULL_INL_


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

    // temporal integration limits (ps)
    const dtime_t LowerLimit = 0.3;
    const dtime_t UpperLimit = 20.0;

    // enable the cubic spline
    const bool CubicSpline = true;

    // model parameters
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
                                      delta0,  deltapar,  deltaperp,  deltaS,
                                      b0, b1,
                         		      p0_OS, p1_OS, DeltaP0_OS, DeltaP1_OS, AvgEta_OS,
                         		      p0_SS, p1_SS, DeltaP0_SS, DeltaP1_SS, AvgEta_SS,
                        		      Omega[0], Omega[1], Omega[2], Omega[3], Omega[4],
                         		      Omega[5], Omega[6], Omega[7], Omega[8], Omega[9],
								      Spline_coeffs[0], Spline_coeffs[1], Spline_coeffs[2],
                                      Spline_coeffs[3], Spline_coeffs[4], Spline_coeffs[5],
                                      Spline_coeffs[6], Spline_coeffs[7], Spline_coeffs[8]};

    auto A0_pd             = hydra::Parameter::Create("A0" ).Value(A0).Error(0.0001);
    auto Aperp_pd          = hydra::Parameter::Create("Aperp").Value(Aperp).Error(0.0001);
    auto AS_pd             = hydra::Parameter::Create("AS" ).Value(AS).Error(0.0001);

    auto DeltaGamma_sd_pd  = hydra::Parameter::Create("DeltaGamma_sd" ).Value(deltagammasd).Error(0.0001);
    auto DeltaGamma_pd     = hydra::Parameter::Create("DeltaGamma").Value(deltagammas).Error(0.0001);
    auto DeltaM_pd         = hydra::Parameter::Create("DeltaM" ).Value(deltams).Error(0.0001);

    auto phi_0_pd          = hydra::Parameter::Create("phi_0").Value(phi0).Error(0.0001);
    auto phi_par_pd        = hydra::Parameter::Create("phi_par" ).Value(phipar).Error(0.0001);
    auto phi_perp_pd       = hydra::Parameter::Create("phi_perp").Value(phiperp).Error(0.0001);
    auto phi_S_pd          = hydra::Parameter::Create("phi_S" ).Value(phiS).Error(0.0001);

    auto lambda_0_pd       = hydra::Parameter::Create("lambda_0").Value(lambda0).Error(0.0001);
    auto lambda_par_pd     = hydra::Parameter::Create("lambda_par" ).Value(lambdapar).Error(0.0001);
    auto lambda_perp_pd    = hydra::Parameter::Create("lambda_perp").Value(lambdaperp).Error(0.0001);
    auto lambda_S_pd       = hydra::Parameter::Create("lambda_S" ).Value(lambdaS).Error(0.0001);

    auto delta_0_pd        = hydra::Parameter::Create("delta_0").Value(delta0).Error(0.0001);
    auto delta_par_pd      = hydra::Parameter::Create("delta_par").Value(deltapar).Error(0.0001);
    auto delta_perp_pd     = hydra::Parameter::Create("delta_perp" ).Value(deltaperp).Error(0.0001);
    auto delta_S_pd        = hydra::Parameter::Create("delta_S").Value(deltaS).Error(0.0001);

    hydra::Parameter ModelParams[18] = {A0_pd,             Aperp_pd,       AS_pd,
                                        DeltaGamma_sd_pd,  DeltaGamma_pd,  DeltaM_pd,
                                        phi_0_pd,          phi_par_pd,     phi_perp_pd,     phi_S_pd,
                                        lambda_0_pd,       lambda_par_pd,  lambda_perp_pd,  lambda_S_pd,
                                        delta_0_pd,        delta_par_pd,   delta_perp_pd,   delta_S_pd};

    auto Model = medusa::FullAnalyticPhis<CubicSpline, dtime_t, theta_h_t, theta_l_t, phi_t,
                                            qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams, Spline_Knots,
                                                                                                            LowerLimit, UpperLimit);


    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------

    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_h;

    GenerateDataset_Full(Model, dataset_h, nentries, nentries, LowerLimit, UpperLimit);
    
    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
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

    // integrator (it always returns the value 1.0, because the normalization is computed
    // in FullAnalyticPhis.h. This choice is justified by the fact that Hydra does not support
    // a normalization factor which depends from the experimental variables)
    auto integrator = hydra::AnalyticalIntegral< medusa::FullAnalyticPhis<CubicSpline,
                                                            dtime_t, theta_h_t, theta_l_t, phi_t, qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto model_PDF = hydra::make_pdf(Model, integrator);


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn = hydra::make_loglikehood_fcn(model_PDF, dataset_d);

    fcn.SetFcnMaxValue(2.22507e+12);


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


    // print non-cached fcn evaluation
    auto start_fcn = std::chrono::high_resolution_clock::now();
    fcn(parameters);
    auto stop_fcn = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed_fcn = stop_fcn - start_fcn;

    std::cout << "fcn = " << fcn(parameters) << std::endl;
    std::cout << "Time (ms) = " << elapsed_fcn.count() << std::endl;

    return 0;

} // main

#endif // BENCHMARKS_SIMPLE_B0S_JPSI_PHI_FULL_INL_