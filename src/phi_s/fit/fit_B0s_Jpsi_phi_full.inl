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
 *  Fit with full analytic model, i.e. signal plus experimental artifacts
 *  (tagging, time resolution and acceptances) with analytical convolution
 *  and integration.
 *---------------------------------------------------------------------------*/

#ifndef FIT_B0S_JPSI_PHI_FULL_INL_
#define FIT_B0S_JPSI_PHI_FULL_INL_


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

// Minuit2
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

// Medusa
#include <medusa/phi_s/phis_full/FullAnalyticPhis.h>
#include <medusa/phi_s/phis_full/GenerateDataset.h>


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
    double distortion = 0;
    double edm = 0;

	try {

		TCLAP::CmdLine cmd("Command line arguments for number of events and Vegas integrator", '=');

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events", false, 5e5, "size_t");
        cmd.add(EArg);

        TCLAP::ValueArg<double> DistortionArg("d", "distortion-factor", "Distortion factor with respect to the dataset", false, 1.0, "double");
        cmd.add(DistortionArg);

        TCLAP::ValueArg<double> EdmArg("e", "EDM", "Estimated vertical distance to minimum", false, 0.1, "double");
		cmd.add(EdmArg);

		// Parse the argv array.
		cmd.parse(argv, argc);

		// Get the value parsed by each arg.
        nentries = EArg.getValue();
        distortion = DistortionArg.getValue();
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

    // temporal integration limits (ps)
    const dtime_t LowerLimit = 0;
    const dtime_t UpperLimit = 20.0;

    // weight to improve the numerical fit
    const double Weight = 0.25;

    // model parameters
    const double A0_dataset         = ::sqrt(0.542);
    const double Aperp_dataset      = ::sqrt(0.206);
    const double AS_dataset         = ::sqrt(0.0037);

    const double phi0_dataset       = -0.082;
    const double phipar_dataset     = -0.125;            // -0.043 + phi0
    const double phiperp_dataset    = -0.156;            // -0.074 + phi0
    const double phiS_dataset       = -0.061;            //  0.021 + phi0

    const double lambda0_dataset    = 0.955; 
    const double lambdapar_dataset  = 0.93399;           // 0.978*lambda0
    const double lambdaperp_dataset = 1.17465;           // 1.23*lambda0
    const double lambdaS_dataset    = 1.2224;            // 1.28*lambda0

    const double delta0_dataset     = 0.0;
    const double deltapar_dataset   = 3.030;            // + delta0
    const double deltaperp_dataset  = 2.60;             // + delta0
    const double deltaS_dataset     = -0.30;            // + delta0

    const double deltagammasd_dataset = -0.0044;
    const double deltagammas_dataset  = 0.0782;
    const double deltams_dataset      = 17.713;

    std::vector<double> parameters_dataset = {A0_dataset, Aperp_dataset, AS_dataset, 
                                              deltagammasd_dataset, deltagammas_dataset, deltams_dataset,
                                              phi0_dataset,    phipar_dataset,    phiperp_dataset,    phiS_dataset,
                                              lambda0_dataset, lambdapar_dataset, lambdaperp_dataset, lambdaS_dataset,
                                              delta0_dataset,  deltapar_dataset,  deltaperp_dataset,  deltaS_dataset};

    auto A0_pd             = hydra::Parameter::Create("A0" ).Value(A0_dataset).Error(0.0001);
    auto Aperp_pd          = hydra::Parameter::Create("Aperp").Value(Aperp_dataset).Error(0.0001);
    auto AS_pd             = hydra::Parameter::Create("AS" ).Value(AS_dataset).Error(0.0001);

    auto DeltaGamma_sd_pd  = hydra::Parameter::Create("DeltaGamma_sd" ).Value(deltagammasd_dataset).Error(0.0001);
    auto DeltaGamma_pd     = hydra::Parameter::Create("DeltaGamma").Value(deltagammas_dataset).Error(0.0001);
    auto DeltaM_pd         = hydra::Parameter::Create("DeltaM" ).Value(deltams_dataset).Error(0.0001);

    auto phi_0_pd          = hydra::Parameter::Create("phi_0").Value(phi0_dataset).Error(0.0001);
    auto phi_par_pd        = hydra::Parameter::Create("phi_par" ).Value(phipar_dataset).Error(0.0001);
    auto phi_perp_pd       = hydra::Parameter::Create("phi_perp").Value(phiperp_dataset).Error(0.0001);
    auto phi_S_pd          = hydra::Parameter::Create("phi_S" ).Value(phiS_dataset).Error(0.0001);

    auto lambda_0_pd       = hydra::Parameter::Create("lambda_0").Value(lambda0_dataset).Error(0.0001);
    auto lambda_par_pd     = hydra::Parameter::Create("lambda_par" ).Value(lambdapar_dataset).Error(0.0001);
    auto lambda_perp_pd    = hydra::Parameter::Create("lambda_perp").Value(lambdaperp_dataset).Error(0.0001);
    auto lambda_S_pd       = hydra::Parameter::Create("lambda_S" ).Value(lambdaS_dataset).Error(0.0001);

    auto delta_0_pd        = hydra::Parameter::Create("delta_0").Value(delta0_dataset).Error(0.0001);
    auto delta_par_pd      = hydra::Parameter::Create("delta_par").Value(deltapar_dataset).Error(0.0001);
    auto delta_perp_pd     = hydra::Parameter::Create("delta_perp" ).Value(deltaperp_dataset).Error(0.0001);
    auto delta_S_pd        = hydra::Parameter::Create("delta_S").Value(deltaS_dataset).Error(0.0001);

    hydra::Parameter ModelParams_dataset[18] = {A0_pd,             Aperp_pd,       AS_pd,
                                                DeltaGamma_sd_pd,  DeltaGamma_pd,  DeltaM_pd,
                                                phi_0_pd,          phi_par_pd,     phi_perp_pd,     phi_S_pd,
                                                lambda_0_pd,       lambda_par_pd,  lambda_perp_pd,  lambda_S_pd,
                                                delta_0_pd,        delta_par_pd,   delta_perp_pd,   delta_S_pd};

    auto Model = medusa::FullAnalyticPhis<dtime_t, theta_h_t, theta_l_t, phi_t,
                                            qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams_dataset, ExpParams, LowerLimit, UpperLimit, Weight);

//    Model(2.54933, 1, 1, 1, 1, 1, 1, 1, 0.01);

    auto Spline = medusa::CubicSpline<7>(Spline_Knots, Spline_coeffs);

    std::cout << Spline.eval(0.4) << std::endl;
    std::cout << Spline.eval(1.0) << std::endl;
    std::cout << Spline.eval(6.99) << std::endl;
    std::cout << Spline.eval(8.0) << std::endl;

    double Gamma = deltagammasd_dataset + 0.65789;
    double HalfDeltaGamma = 0.5*deltagammas_dataset;

    double Spline_int_conv_exp_cosh[4];
    double Spline_int_conv_exp_sinh[4];
    double Spline_int_conv_exp_cos[4];
    double Spline_int_conv_exp_sin[4];

    for(size_t i=0; i<4; i++)
    {
        Spline_int_conv_exp_cosh[i] = Spline.Integrate_t_to_k_times_convolved_exp_sinhcosh(9.22784, i, Gamma, HalfDeltaGamma, 0, 0.0250024, LowerLimit, UpperLimit, 1);
        Spline_int_conv_exp_sinh[i] = Spline.Integrate_t_to_k_times_convolved_exp_sinhcosh(9.22784, i, Gamma, HalfDeltaGamma, 0, 0.0250024, LowerLimit, UpperLimit, -1);
        Spline_int_conv_exp_cos[i] = Spline.Integrate_t_to_k_times_convolved_exp_sincos(9.22784, i, Gamma, deltams_dataset, 0, 0.0250024, LowerLimit, UpperLimit, 1);
        Spline_int_conv_exp_sin[i] = Spline.Integrate_t_to_k_times_convolved_exp_sincos(9.22784, i, Gamma, deltams_dataset, 0, 0.0250024, LowerLimit, UpperLimit, -1);

        std::cout << "int_conv_cosh [" << i << "] = " << Spline_int_conv_exp_cosh[i] << std::endl;
        std::cout << "int_conv_sinh [" << i << "] = " << Spline_int_conv_exp_sinh[i] << std::endl;
        std::cout << "int_conv_cos [" << i << "] = " << Spline_int_conv_exp_cos[i] << std::endl;
        std::cout << "int_conv_sin [" << i << "] = " << Spline_int_conv_exp_sin[i] << std::endl;
    }
/*
    double int_conv_exp_cosh = medusa::functions::Integrate_convolved_exp_sinhcosh(9.22784, Gamma, HalfDeltaGamma, 0, 0.0250024, LowerLimit, UpperLimit, 1);
    double int_conv_exp_sinh = medusa::functions::Integrate_convolved_exp_sinhcosh(9.22784, Gamma, HalfDeltaGamma, 0, 0.0250024, LowerLimit, UpperLimit, -1);
    double int_conv_exp_cos = medusa::functions::Integrate_convolved_exp_sincos(9.22784, Gamma, deltams_dataset, 0, 0.0250024, LowerLimit, UpperLimit, 1);
    double int_conv_exp_sin = medusa::functions::Integrate_convolved_exp_sincos(9.22784, Gamma, deltams_dataset, 0, 0.0250024, LowerLimit, UpperLimit, -1);

    const double f = 0.2387324146378430; // 3./(4.*PI)

    double NormFactor = f * ( int_conv_exp_cosh + int_conv_exp_sinh + int_conv_exp_cos + int_conv_exp_sin );

    std::cout << int_conv_exp_cosh << std::endl;
    std::cout << int_conv_exp_sinh << std::endl;
    std::cout << int_conv_exp_cos << std::endl;
    std::cout << int_conv_exp_sin << std::endl;
    std::cout << NormFactor << std::endl;
*/
    //---------------------------------
    //  Unweighted dataset generation
    //---------------------------------
/*
    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::host::sys_t> dataset_h;

    GenerateDataset_FullAnalyticPhis(Model, dataset_h, nentries, nentries, LowerLimit, UpperLimit);
    
    hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t,
                                    qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> , hydra::device::sys_t> dataset_d(dataset_h.size());
    hydra::copy(dataset_h, dataset_d);


    //-----------------------------------------
    //  Print and plot the unweighted dataset
    //-----------------------------------------
    
    for( size_t i=0; i<10; i++ )
        std::cout <<"Dataset_h: {"<< dataset_h[i]  << "}"<< std::endl;


    #ifdef _ROOT_AVAILABLE_

        TH1D timedist("timedist","Decay Time; time (ps); Candidates / bin", 100, 0, 20);
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


    //---------------------------------------------
    //   Set the starting values for the fit
    //---------------------------------------------

    // model parameters
    const double A0         = ::sqrt(0.542)*distortion;
    const double Aperp      = ::sqrt(0.206)*distortion;
    const double AS         = ::sqrt(0.0037)*distortion;

    const double phi0       = -0.082*distortion;
    const double phipar     = -0.125*distortion;            // -0.043 + phi0
    const double phiperp    = -0.156*distortion;            // -0.074 + phi0
    const double phiS       = -0.061*distortion;            //  0.021 + phi0

    const double lambda0    = 0.955*distortion; 
    const double lambdapar  = 0.93399*distortion;           // 0.978*lambda0
    const double lambdaperp = 1.17465*distortion;           // 1.23*lambda0
    const double lambdaS    = 1.2224*distortion;            // 1.28*lambda0

    const double delta0     = 0.0;
    const double deltapar   = 3.030*distortion;            // + delta0
    const double deltaperp  = 2.60*distortion;             // + delta0
    const double deltaS     = -0.30*distortion;            // + delta0

    const double deltagammasd = -0.0044*distortion;
    const double deltagammas  = 0.0782*distortion;
    const double deltams      = 17.713*distortion;

    std::vector<double> parameters = {A0, Aperp, AS,
                                      deltagammasd, deltagammas, deltams,
                                      phi0,    phipar,    phiperp,    phiS,
                                      lambda0, lambdapar, lambdaperp, lambdaS,
                                      delta0,  deltapar,  deltaperp,  deltaS};

    auto A0_p             = hydra::Parameter::Create("A0" ).Value(A0).Error(0.0001).Limits(0.1, 0.9);
    auto Aperp_p          = hydra::Parameter::Create("Aperp").Value(Aperp).Error(0.0001).Limits(0.1, 0.9);
    auto AS_p             = hydra::Parameter::Create("AS" ).Value(AS).Error(0.0001).Limits(-0.1, 0.8);

    auto DeltaGamma_sd_p  = hydra::Parameter::Create("DeltaGamma_sd" ).Value(deltagammasd).Error(0.0001).Limits(-0.2, 0.2);
    auto DeltaGamma_p     = hydra::Parameter::Create("DeltaGamma").Value(deltagammas).Error(0.0001).Limits(0.03, 0.15);
    auto DeltaM_p         = hydra::Parameter::Create("DeltaM" ).Value(deltams).Error(0.0001).Limits(16.0, 20.0);

    auto phi_0_p          = hydra::Parameter::Create("phi_0").Value(phi0).Error(0.0001).Limits(-1.0, 1.0);
    auto phi_par_p        = hydra::Parameter::Create("phi_par" ).Value(phipar).Error(0.0001).Limits(-1.0, 1.0);
    auto phi_perp_p       = hydra::Parameter::Create("phi_perp").Value(phiperp).Error(0.0001).Limits(-1.0, 1.0);
    auto phi_S_p          = hydra::Parameter::Create("phi_S" ).Value(phiS).Error(0.0001).Limits(-1.0, 1.0);

    auto lambda_0_p       = hydra::Parameter::Create("lambda_0").Value(lambda0).Error(0.0001).Limits(0.7, 1.6);
    auto lambda_par_p     = hydra::Parameter::Create("lambda_par" ).Value(lambdapar).Error(0.0001).Limits(0.7, 1.6);
    auto lambda_perp_p    = hydra::Parameter::Create("lambda_perp").Value(lambdaperp).Error(0.0001).Limits(0.7, 1.6);
    auto lambda_S_p       = hydra::Parameter::Create("lambda_S" ).Value(lambdaS).Error(0.0001).Limits(0.7, 1.6);

    auto delta_0_p        = hydra::Parameter::Create("delta_0").Value(delta0).Error(0.0001).Limits(-6.0, 6.0);
    auto delta_par_p      = hydra::Parameter::Create("delta_par").Value(deltapar).Error(0.0001).Limits(-6.28, 6.28);
    auto delta_perp_p     = hydra::Parameter::Create("delta_perp" ).Value(deltaperp).Error(0.0001).Limits(-6.28, 6.28);
    auto delta_S_p        = hydra::Parameter::Create("delta_S").Value(deltaS).Error(0.0001).Limits(-6.0, 6.0);

    hydra::Parameter ModelParams[18] = {A0_p,            Aperp_p,      AS_p,
                                        DeltaGamma_sd_p, DeltaGamma_p, DeltaM_p,
                                        phi_0_p,         phi_par_p,    phi_perp_p,    phi_S_p,
                                        lambda_0_p,      lambda_par_p, lambda_perp_p, lambda_S_p,
                                        delta_0_p,       delta_par_p,  delta_perp_p,  delta_S_p};

    auto model = medusa::FullAnalyticPhis<dtime_t, theta_h_t, theta_l_t, phi_t,
                                            qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t>(ModelParams, ExpParams, LowerLimit, UpperLimit, Weight);

    //---------------------------------
    //          PDF generation
    //---------------------------------

    // integrator (it always returns the value 1.0, because the normalization is computed
    // in FullAnalyticPhis.h. This choice is justified by the fact that Hydra does not support
    // a normalization factor which depends from the experimental variables)
    auto integrator = hydra::AnalyticalIntegral<
                        medusa::FullAnalyticPhis<dtime_t, theta_h_t, theta_l_t, phi_t, qOS_t, qSS_t, etaOS_t, etaSS_t, delta_t> >(LowerLimit, UpperLimit);

    // make PDF
    auto model_PDF = hydra::make_pdf(model, integrator);


    //---------------------------------
    //          FCN generation
    //---------------------------------

    auto fcn = hydra::make_loglikehood_fcn(model_PDF, dataset_d);

    fcn.SetFcnMaxValue(2.22507e+12);


    //---------------------------------
    //          fit by Minuit2
    //---------------------------------

    // print level
    MnPrint::SetLevel(3);

    // minimization strategy
	MnStrategy strategy(2);

    // create Minimize minimizer
	MnMinimize minimize(fcn, fcn.GetParameters().GetMnState(), strategy);

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
	std::cout << "| Time (ms) ="<< elapsed.count()          << std::endl;
	std::cout << "-----------------------------------------"<< std::endl;
*/
    return 0;

} // main

#endif // FIT_B0S_JPSI_PHI_FULL_INL_