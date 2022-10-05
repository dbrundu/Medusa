/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu,
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto,
 * 						Alessandro Maria Ricci
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
 * 	Parameters.h
 *
 *  Created on: 04/06/2021
 *      Author: Alessandro Maria Ricci
 * 
 * 	This library contains all parameters used by PhisSignal.h
 *  and FullAnalyticPhis.h.
 *---------------------------------------------------------------------------*/

#ifndef PHIS_PARAMETERS_H_
#define PHIS_PARAMETERS_H_

// Hydra
#include <hydra/detail/Config.h>
#include <hydra/Vector4R.h>
#include <hydra/Parameter.h>

// Medusa
#include <medusa/generic/Constants.h>

// default namespaces
using namespace hydra::arguments;


//-------------------------------------
// New types for:
//  B0s -> J/psi  (Phi -> K+ K-)
//          |-> mu+ mu-
//-------------------------------------

// Particle types
declarg(Jpsi,  hydra::Vector4R)
declarg(Phi,   hydra::Vector4R)
declarg(KaonP, hydra::Vector4R)
declarg(KaonM, hydra::Vector4R)
declarg(MuonP, hydra::Vector4R)
declarg(MuonM, hydra::Vector4R)

// New types for time and helicity angles
declarg(costheta_h_t,  double)
declarg(costheta_l_t,  double)
declarg(phi_t,         double)
declarg(dtime_t,       double)

// New types for tagging variables
declarg(qOS_t, 	 int)
declarg(qSS_t, 	 int)
declarg(etaOS_t, double)
declarg(etaSS_t, double)
declarg(delta_t, double)


//-----------------------------------------------
// Constant parameters for:
//  B0s -> J/psi  (Phi -> K+ K-)
//          |-> mu+ mu-
//
// For FullAnalyticPhis.h or PhisSignal.h
//-----------------------------------------------

// Temporal integration limits (in ps) for FullAnalyticPhis.h
const dtime_t LowerLimit = 0.3;
const dtime_t UpperLimit = 15.0;

// Enable the cubic spline for FullAnalyticPhis.h
const bool CubicSpline = true;

// Specify wether B0s is B0sbar or not for PhisSignal.h
const bool B0sbar = false;


//-----------------------------------------------
// Model parameters for:
//  B0s -> J/psi  (Phi -> K+ K-)
//          |-> mu+ mu-
//
// For FullAnalyticPhis.h and PhisSignal.h
//-----------------------------------------------

const double A02         = 0.542;
const double Aperp2      = 0.206;
const double FS1         = 0.491;			 // FS = AS^2
const double FS2         = 0.0406;
const double FS3         = 0.0044;
const double FS4         = 0.0069;
const double FS5         = 0.073;
const double FS6         = 0.151;

const double phi0        = -0.082;
const double phipar0     = -0.043;            // phipar - phi0
const double phiperp0    = -0.074;            // phiperp - phi0
const double phiS0       =  0.021;            // phiS - phi0

const double lambda0     = 0.955; 
const double lambdapar0  = 0.978;            // lambdapar / lambda0
const double lambdaperp0 = 1.23;         	 // lambdaperp / lambda0
const double lambdaS0    = 1.28;             // lambdaS / lambda0

const double deltapar0   =  3.030;           // deltapar - delta0
const double deltaperp0  =  2.60;            // deltaperp - delta0
const double deltaS1perp =  2.21;            // deltaS1 - deltaperp
const double deltaS2perp =  1.55;
const double deltaS3perp =  1.07;
const double deltaS4perp = -0.28;
const double deltaS5perp = -0.536;
const double deltaS6perp = -1.10;

const double deltagammasd = -0.0044;
const double deltagammas  =  0.0782;
const double deltams      =  17.713;

auto A02_p            = hydra::Parameter::Create("A02" ).Value(A02).Error(0.0001).Limits(0., 0.9);
auto Aperp2_p         = hydra::Parameter::Create("Aperp2").Value(Aperp2).Error(0.0001).Limits(0., 0.9);
auto FS1_p            = hydra::Parameter::Create("FS1" ).Value(FS1).Error(0.0001).Limits(0., 0.8);
auto FS2_p            = hydra::Parameter::Create("FS2" ).Value(FS2).Error(0.0001).Limits(0., 0.8);
auto FS3_p            = hydra::Parameter::Create("FS3" ).Value(FS3).Error(0.0001).Limits(0., 0.8);
auto FS4_p            = hydra::Parameter::Create("FS4" ).Value(FS4).Error(0.0001).Limits(0., 0.8);
auto FS5_p            = hydra::Parameter::Create("FS5" ).Value(FS5).Error(0.0001).Limits(0., 0.8);
auto FS6_p            = hydra::Parameter::Create("FS6" ).Value(FS6).Error(0.0001).Limits(0., 0.8);

auto DeltaGamma_sd_p  = hydra::Parameter::Create("DeltaGamma_sd" ).Value(deltagammasd).Error(0.0001).Limits(-0.2, 0.2);
auto DeltaGamma_p     = hydra::Parameter::Create("DeltaGamma").Value(deltagammas).Error(0.0001).Limits(0., 0.15);
auto DeltaM_p         = hydra::Parameter::Create("DeltaM" ).Value(deltams).Error(0.0001).Limits(16.0, 20.0);

auto phi_0_p          = hydra::Parameter::Create("phi_0").Value(phi0).Error(0.0001).Limits(-1.0, 1.0);
auto phi_par0_p       = hydra::Parameter::Create("phi_par0" ).Value(phipar0).Error(0.0001).Limits(-1.0, 1.0);
auto phi_perp0_p      = hydra::Parameter::Create("phi_perp0").Value(phiperp0).Error(0.0001).Limits(-1.0, 1.0);
auto phi_S0_p         = hydra::Parameter::Create("phi_S0" ).Value(phiS0).Error(0.0001).Limits(-1.0, 1.0);

auto lambda_0_p       = hydra::Parameter::Create("lambda_0").Value(lambda0).Error(0.0001).Limits(0.7, 1.6);
auto lambda_par0_p    = hydra::Parameter::Create("lambda_par0" ).Value(lambdapar0).Error(0.0001).Limits(0.7, 1.6);
auto lambda_perp0_p   = hydra::Parameter::Create("lambda_perp0").Value(lambdaperp0).Error(0.0001).Limits(0.7, 1.6);
auto lambda_S0_p      = hydra::Parameter::Create("lambda_S0" ).Value(lambdaS0).Error(0.0001).Limits(0.7, 1.6);

auto delta_par0_p     = hydra::Parameter::Create("delta_par0").Value(deltapar0).Error(0.0001).Limits(-6.28, 6.28);
auto delta_perp0_p    = hydra::Parameter::Create("delta_perp0" ).Value(deltaperp0).Error(0.0001).Limits(-6.28, 6.28);
auto delta_S1perp_p   = hydra::Parameter::Create("delta_S1perp").Value(deltaS1perp).Error(0.0001).Limits(-6.0, 6.0);
auto delta_S2perp_p   = hydra::Parameter::Create("delta_S2perp").Value(deltaS2perp).Error(0.0001).Limits(-6.0, 6.0);
auto delta_S3perp_p   = hydra::Parameter::Create("delta_S3perp").Value(deltaS3perp).Error(0.0001).Limits(-6.0, 6.0);
auto delta_S4perp_p   = hydra::Parameter::Create("delta_S4perp").Value(deltaS4perp).Error(0.0001).Limits(-6.0, 6.0);
auto delta_S5perp_p   = hydra::Parameter::Create("delta_S5perp").Value(deltaS5perp).Error(0.0001).Limits(-6.0, 6.0);
auto delta_S6perp_p   = hydra::Parameter::Create("delta_S6perp").Value(deltaS6perp).Error(0.0001).Limits(-6.0, 6.0);

hydra::Parameter ModelParams_S1[17] = {A02_p,           Aperp2_p,      FS1_p,
                                       DeltaGamma_sd_p, DeltaGamma_p,  DeltaM_p,
                                       phi_0_p,         phi_par0_p,    phi_perp0_p,    phi_S0_p,
                                       lambda_0_p,      lambda_par0_p, lambda_perp0_p, lambda_S0_p,
                                       delta_par0_p,    delta_perp0_p, delta_S1perp_p};

hydra::Parameter ModelParams_S2[17] = {A02_p,           Aperp2_p,      FS2_p,
                                       DeltaGamma_sd_p, DeltaGamma_p,  DeltaM_p,
                                       phi_0_p,         phi_par0_p,    phi_perp0_p,    phi_S0_p,
                                       lambda_0_p,      lambda_par0_p, lambda_perp0_p, lambda_S0_p,
                                       delta_par0_p,    delta_perp0_p, delta_S2perp_p};

hydra::Parameter ModelParams_S3[17] = {A02_p,           Aperp2_p,      FS3_p,
                                       DeltaGamma_sd_p, DeltaGamma_p,  DeltaM_p,
                                       phi_0_p,         phi_par0_p,    phi_perp0_p,    phi_S0_p,
                                       lambda_0_p,      lambda_par0_p, lambda_perp0_p, lambda_S0_p,
                                       delta_par0_p,    delta_perp0_p, delta_S3perp_p};

hydra::Parameter ModelParams_S4[17] = {A02_p,           Aperp2_p,      FS4_p,
                                       DeltaGamma_sd_p, DeltaGamma_p,  DeltaM_p,
                                       phi_0_p,         phi_par0_p,    phi_perp0_p,    phi_S0_p,
                                       lambda_0_p,      lambda_par0_p, lambda_perp0_p, lambda_S0_p,
                                       delta_par0_p,    delta_perp0_p, delta_S4perp_p};

hydra::Parameter ModelParams_S5[17] = {A02_p,           Aperp2_p,      FS5_p,
                                       DeltaGamma_sd_p, DeltaGamma_p,  DeltaM_p,
                                       phi_0_p,         phi_par0_p,    phi_perp0_p,    phi_S0_p,
                                       lambda_0_p,      lambda_par0_p, lambda_perp0_p, lambda_S0_p,
                                       delta_par0_p,    delta_perp0_p, delta_S5perp_p};

hydra::Parameter ModelParams_S6[17] = {A02_p,           Aperp2_p,      FS6_p,
                                       DeltaGamma_sd_p, DeltaGamma_p,  DeltaM_p,
                                       phi_0_p,         phi_par0_p,    phi_perp0_p,    phi_S0_p,
                                       lambda_0_p,      lambda_par0_p, lambda_perp0_p, lambda_S0_p,
                                       delta_par0_p,    delta_perp0_p, delta_S6perp_p};


//-------------------------------------
// Experimental parameters for:
//  B0s -> J/psi  (Phi -> K+ K-)
//          |-> mu+ mu-
//
//  Only for FullAnalyticPhis.h
//-------------------------------------

const double b0 = 0.01297;
const double b1 = 0.8446;

const double p0_OS = 0.3890;
const double p1_OS = 0.8486;
const double DeltaP0_OS = 0.009;
const double DeltaP1_OS = 0.0143;
const double AvgEta_OS = 0.3602;

const double p0_SS = 0.4325;
const double p1_SS = 0.9241;
const double DeltaP0_SS = 0.;
const double DeltaP1_SS = 0.;
const double AvgEta_SS = 0.4167;


const double Spline_Knots[7] = {0.3, 0.58, 0.91, 1.35, 1.96, 3.01, 7.00};
double Spline_Coeffs_2015_unbiased[9] = {1.0, 1.05, 1.097, 0.969, 1.051, 1.05, 1.028, 1.094, 1.051};
double Spline_Coeffs_2015_biased[9] = {1.0, 1.69, 1.73, 1.85, 1.99, 1.92, 2.0, 2.19, 1.95};
double Spline_Coeffs_2016_unbiased[9] = {1.0, 1.008, 1.031, 1.001, 0.984, 1.0, 1.009, 0.989, 0.987};
double Spline_Coeffs_2016_biased[9] = {1.0, 1.49, 2.06, 2.12, 2.28, 2.29, 2.46, 2.25, 2.34};
double Spline_CoeffErr_2015_unbiased[9] = {0., 0.07, 0.049, 0.049, 0.051, 0.047, 0.064, 0.06, 0.05};
double Spline_CoeffErr_2015_biased[9] = {0., 0.31, 0.21, 0.27, 0.26, 0.26, 0.3, 0.31, 0.27};
double Spline_CoeffErr_2016_unbiased[9] = {0., 0.028, 0.018, 0.021, 0.018, 0.019, 0.024, 0.024, 0.02};
double Spline_CoeffErr_2016_biased[9] = {0., 0.16, 0.13, 0.17, 0.16, 0.17, 0.19, 0.18, 0.17};


const double Omega_2015_unbiased[10] = {1.0, 1.0434, 1.0442, -0.0026, -0.00142, 0.00139, 1.0156, -0.0014, 0.0006, -0.0171};
const double Omega_2015_biased[10] = {1.0, 1.0463, 1.0445, -0.0105, 0.0037, 0.0023, 1.0262, -0.0045, -0.0007, -0.0348};
const double Omega_2016_unbiased[10] = {1.0, 1.03788, 1.03765, -0.00079, 0.00026, 0.00023, 1.01022, 0.00004, 0.00010, -0.00362};
const double Omega_2016_biased[10] = {1.0, 1.0336, 1.0336, 0.0028, 0.00298, -0.0002, 1.0196, 0.00019, 0.00019, 0.0057};
const double OmegaErr_2015_unbiased[10] = {0.0, 0.0020, 0.0020, 0.0016, 0.00094, 0.00093, 0.0014, 0.0012, 0.0012, 0.0026};
const double OmegaErr_2015_biased[10] = {0.0, 0.0039, 0.0038, 0.0032, 0.0018, 0.0018, 0.0027, 0.0024, 0.0024, 0.0050};
const double OmegaErr_2016_unbiased[10] = {0.0, 0.00070, 0.00069, 0.00054, 0.00033, 0.00033, 0.00047, 0.00042, 0.00043, 0.00089};
const double OmegaErr_2016_biased[10] = {0.0, 0.0015, 0.0015, 0.0013, 0.00074, 0.00072, 0.0011, 0.00094, 0.00094, 0.0019};


const double Csp[6] = {0.8463, 0.8756, 0.8478, 0.8833, 0.9415, 0.9756};
const int KKmass[7] = {990, 1008, 1016, 1020, 1024, 1032, 1050};


auto b0_p = hydra::Parameter::Create("b0").Value(b0).Error(0.00022).Fixed();
auto b1_p = hydra::Parameter::Create("b1").Value(b1).Error(0.0057).Fixed();

auto p0_OS_p = hydra::Parameter::Create("p0_OS").Value(p0_OS).Error(0.0007 + 0.0028).Fixed();
auto p1_OS_p = hydra::Parameter::Create("p1_OS").Value(p1_OS).Error(0.0062 + 0.0265).Fixed();
auto DeltaP0_OS_p = hydra::Parameter::Create("DeltaP0_OS").Value(DeltaP0_OS).Error(0.0014).Fixed();
auto DeltaP1_OS_p = hydra::Parameter::Create("DeltaP1_OS").Value(DeltaP1_OS).Error(0.0124).Fixed();
auto AvgEta_OS_p = hydra::Parameter::Create("AvgEta_OS").Value(AvgEta_OS).Error(0.0).Fixed();

auto p0_SS_p = hydra::Parameter::Create("p0_SS").Value(p0_SS).Error(0.0108 + 0.003).Fixed();
auto p1_SS_p = hydra::Parameter::Create("p1_SS").Value(p1_SS).Error(0.1314 + 0.0196).Fixed();
auto DeltaP0_SS_p = hydra::Parameter::Create("DeltaP0_SS").Value(DeltaP0_SS).Error(0.03).Fixed();
auto DeltaP1_SS_p = hydra::Parameter::Create("DeltaP1_SS").Value(DeltaP1_SS).Error(0.03).Fixed();
auto AvgEta_SS_p = hydra::Parameter::Create("AvgEta_SS").Value(AvgEta_SS).Error(0.0).Fixed();

auto Csp1_p = hydra::Parameter::Create("Csp1").Value(Csp[0]).Error(0.0).Fixed();
auto Csp2_p = hydra::Parameter::Create("Csp2").Value(Csp[1]).Error(0.0).Fixed();
auto Csp3_p = hydra::Parameter::Create("Csp3").Value(Csp[2]).Error(0.0).Fixed();
auto Csp4_p = hydra::Parameter::Create("Csp4").Value(Csp[3]).Error(0.0).Fixed();
auto Csp5_p = hydra::Parameter::Create("Csp5").Value(Csp[4]).Error(0.0).Fixed();
auto Csp6_p = hydra::Parameter::Create("Csp6").Value(Csp[5]).Error(0.0).Fixed();

auto c0_2015_unbiased_p = hydra::Parameter::Create("c0_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[0]).Error(Spline_CoeffErr_2015_unbiased[0]).Fixed();
auto c1_2015_unbiased_p = hydra::Parameter::Create("c1_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[1]).Error(Spline_CoeffErr_2015_unbiased[1]).Fixed();
auto c2_2015_unbiased_p = hydra::Parameter::Create("c2_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[2]).Error(Spline_CoeffErr_2015_unbiased[2]).Fixed();
auto c3_2015_unbiased_p = hydra::Parameter::Create("c3_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[3]).Error(Spline_CoeffErr_2015_unbiased[3]).Fixed();
auto c4_2015_unbiased_p = hydra::Parameter::Create("c4_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[4]).Error(Spline_CoeffErr_2015_unbiased[4]).Fixed();
auto c5_2015_unbiased_p = hydra::Parameter::Create("c5_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[5]).Error(Spline_CoeffErr_2015_unbiased[5]).Fixed();
auto c6_2015_unbiased_p = hydra::Parameter::Create("c6_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[6]).Error(Spline_CoeffErr_2015_unbiased[6]).Fixed();
auto c7_2015_unbiased_p = hydra::Parameter::Create("c7_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[7]).Error(Spline_CoeffErr_2015_unbiased[7]).Fixed();
auto c8_2015_unbiased_p = hydra::Parameter::Create("c8_2015_unbiased").Value(Spline_Coeffs_2015_unbiased[8]).Error(Spline_CoeffErr_2015_unbiased[8]).Fixed();

auto c0_2015_biased_p = hydra::Parameter::Create("c0_2015_biased").Value(Spline_Coeffs_2015_biased[0]).Error(Spline_CoeffErr_2015_biased[0]).Fixed();
auto c1_2015_biased_p = hydra::Parameter::Create("c1_2015_biased").Value(Spline_Coeffs_2015_biased[1]).Error(Spline_CoeffErr_2015_biased[1]).Fixed();
auto c2_2015_biased_p = hydra::Parameter::Create("c2_2015_biased").Value(Spline_Coeffs_2015_biased[2]).Error(Spline_CoeffErr_2015_biased[2]).Fixed();
auto c3_2015_biased_p = hydra::Parameter::Create("c3_2015_biased").Value(Spline_Coeffs_2015_biased[3]).Error(Spline_CoeffErr_2015_biased[3]).Fixed();
auto c4_2015_biased_p = hydra::Parameter::Create("c4_2015_biased").Value(Spline_Coeffs_2015_biased[4]).Error(Spline_CoeffErr_2015_biased[4]).Fixed();
auto c5_2015_biased_p = hydra::Parameter::Create("c5_2015_biased").Value(Spline_Coeffs_2015_biased[5]).Error(Spline_CoeffErr_2015_biased[5]).Fixed();
auto c6_2015_biased_p = hydra::Parameter::Create("c6_2015_biased").Value(Spline_Coeffs_2015_biased[6]).Error(Spline_CoeffErr_2015_biased[6]).Fixed();
auto c7_2015_biased_p = hydra::Parameter::Create("c7_2015_biased").Value(Spline_Coeffs_2015_biased[7]).Error(Spline_CoeffErr_2015_biased[7]).Fixed();
auto c8_2015_biased_p = hydra::Parameter::Create("c8_2015_biased").Value(Spline_Coeffs_2015_biased[8]).Error(Spline_CoeffErr_2015_biased[8]).Fixed();

auto c0_2016_unbiased_p = hydra::Parameter::Create("c0_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[0]).Error(Spline_CoeffErr_2016_unbiased[0]).Fixed();
auto c1_2016_unbiased_p = hydra::Parameter::Create("c1_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[1]).Error(Spline_CoeffErr_2016_unbiased[1]).Fixed();
auto c2_2016_unbiased_p = hydra::Parameter::Create("c2_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[2]).Error(Spline_CoeffErr_2016_unbiased[2]).Fixed();
auto c3_2016_unbiased_p = hydra::Parameter::Create("c3_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[3]).Error(Spline_CoeffErr_2016_unbiased[3]).Fixed();
auto c4_2016_unbiased_p = hydra::Parameter::Create("c4_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[4]).Error(Spline_CoeffErr_2016_unbiased[4]).Fixed();
auto c5_2016_unbiased_p = hydra::Parameter::Create("c5_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[5]).Error(Spline_CoeffErr_2016_unbiased[5]).Fixed();
auto c6_2016_unbiased_p = hydra::Parameter::Create("c6_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[6]).Error(Spline_CoeffErr_2016_unbiased[6]).Fixed();
auto c7_2016_unbiased_p = hydra::Parameter::Create("c7_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[7]).Error(Spline_CoeffErr_2016_unbiased[7]).Fixed();
auto c8_2016_unbiased_p = hydra::Parameter::Create("c8_2016_unbiased").Value(Spline_Coeffs_2016_unbiased[8]).Error(Spline_CoeffErr_2016_unbiased[8]).Fixed();

auto c0_2016_biased_p = hydra::Parameter::Create("c0_2016_biased").Value(Spline_Coeffs_2016_biased[0]).Error(Spline_CoeffErr_2016_biased[0]).Fixed();
auto c1_2016_biased_p = hydra::Parameter::Create("c1_2016_biased").Value(Spline_Coeffs_2016_biased[1]).Error(Spline_CoeffErr_2016_biased[1]).Fixed();
auto c2_2016_biased_p = hydra::Parameter::Create("c2_2016_biased").Value(Spline_Coeffs_2016_biased[2]).Error(Spline_CoeffErr_2016_biased[2]).Fixed();
auto c3_2016_biased_p = hydra::Parameter::Create("c3_2016_biased").Value(Spline_Coeffs_2016_biased[3]).Error(Spline_CoeffErr_2016_biased[3]).Fixed();
auto c4_2016_biased_p = hydra::Parameter::Create("c4_2016_biased").Value(Spline_Coeffs_2016_biased[4]).Error(Spline_CoeffErr_2016_biased[4]).Fixed();
auto c5_2016_biased_p = hydra::Parameter::Create("c5_2016_biased").Value(Spline_Coeffs_2016_biased[5]).Error(Spline_CoeffErr_2016_biased[5]).Fixed();
auto c6_2016_biased_p = hydra::Parameter::Create("c6_2016_biased").Value(Spline_Coeffs_2016_biased[6]).Error(Spline_CoeffErr_2016_biased[6]).Fixed();
auto c7_2016_biased_p = hydra::Parameter::Create("c7_2016_biased").Value(Spline_Coeffs_2016_biased[7]).Error(Spline_CoeffErr_2016_biased[7]).Fixed();
auto c8_2016_biased_p = hydra::Parameter::Create("c8_2016_biased").Value(Spline_Coeffs_2016_biased[8]).Error(Spline_CoeffErr_2016_biased[8]).Fixed();

auto Omega1_2015_unbiased_p = hydra::Parameter::Create("Omega1_2015_unbiased").Value(Omega_2015_unbiased[0]).Error(OmegaErr_2015_unbiased[0]).Fixed();
auto Omega2_2015_unbiased_p = hydra::Parameter::Create("Omega2_2015_unbiased").Value(Omega_2015_unbiased[1]).Error(OmegaErr_2015_unbiased[1]).Fixed();
auto Omega3_2015_unbiased_p = hydra::Parameter::Create("Omega3_2015_unbiased").Value(Omega_2015_unbiased[2]).Error(OmegaErr_2015_unbiased[2]).Fixed();
auto Omega4_2015_unbiased_p = hydra::Parameter::Create("Omega4_2015_unbiased").Value(Omega_2015_unbiased[3]).Error(OmegaErr_2015_unbiased[3]).Fixed();
auto Omega5_2015_unbiased_p = hydra::Parameter::Create("Omega5_2015_unbiased").Value(Omega_2015_unbiased[4]).Error(OmegaErr_2015_unbiased[4]).Fixed();
auto Omega6_2015_unbiased_p = hydra::Parameter::Create("Omega6_2015_unbiased").Value(Omega_2015_unbiased[5]).Error(OmegaErr_2015_unbiased[5]).Fixed();
auto Omega7_2015_unbiased_p = hydra::Parameter::Create("Omega7_2015_unbiased").Value(Omega_2015_unbiased[6]).Error(OmegaErr_2015_unbiased[6]).Fixed();
auto Omega8_2015_unbiased_p = hydra::Parameter::Create("Omega8_2015_unbiased").Value(Omega_2015_unbiased[7]).Error(OmegaErr_2015_unbiased[7]).Fixed();
auto Omega9_2015_unbiased_p = hydra::Parameter::Create("Omega9_2015_unbiased").Value(Omega_2015_unbiased[8]).Error(OmegaErr_2015_unbiased[8]).Fixed();
auto Omega10_2015_unbiased_p = hydra::Parameter::Create("Omega10_2015_unbiased").Value(Omega_2015_unbiased[9]).Error(OmegaErr_2015_unbiased[9]).Fixed();

auto Omega1_2015_biased_p = hydra::Parameter::Create("Omega1_2015_biased").Value(Omega_2015_biased[0]).Error(OmegaErr_2015_biased[0]).Fixed();
auto Omega2_2015_biased_p = hydra::Parameter::Create("Omega2_2015_biased").Value(Omega_2015_biased[1]).Error(OmegaErr_2015_biased[1]).Fixed();
auto Omega3_2015_biased_p = hydra::Parameter::Create("Omega3_2015_biased").Value(Omega_2015_biased[2]).Error(OmegaErr_2015_biased[2]).Fixed();
auto Omega4_2015_biased_p = hydra::Parameter::Create("Omega4_2015_biased").Value(Omega_2015_biased[3]).Error(OmegaErr_2015_biased[3]).Fixed();
auto Omega5_2015_biased_p = hydra::Parameter::Create("Omega5_2015_biased").Value(Omega_2015_biased[4]).Error(OmegaErr_2015_biased[4]).Fixed();
auto Omega6_2015_biased_p = hydra::Parameter::Create("Omega6_2015_biased").Value(Omega_2015_biased[5]).Error(OmegaErr_2015_biased[5]).Fixed();
auto Omega7_2015_biased_p = hydra::Parameter::Create("Omega7_2015_biased").Value(Omega_2015_biased[6]).Error(OmegaErr_2015_biased[6]).Fixed();
auto Omega8_2015_biased_p = hydra::Parameter::Create("Omega8_2015_biased").Value(Omega_2015_biased[7]).Error(OmegaErr_2015_biased[7]).Fixed();
auto Omega9_2015_biased_p = hydra::Parameter::Create("Omega9_2015_biased").Value(Omega_2015_biased[8]).Error(OmegaErr_2015_biased[8]).Fixed();
auto Omega10_2015_biased_p = hydra::Parameter::Create("Omega10_2015_biased").Value(Omega_2015_biased[9]).Error(OmegaErr_2015_biased[9]).Fixed();

auto Omega1_2016_unbiased_p = hydra::Parameter::Create("Omega1_2016_unbiased").Value(Omega_2016_unbiased[0]).Error(OmegaErr_2016_unbiased[0]).Fixed();
auto Omega2_2016_unbiased_p = hydra::Parameter::Create("Omega2_2016_unbiased").Value(Omega_2016_unbiased[1]).Error(OmegaErr_2016_unbiased[1]).Fixed();
auto Omega3_2016_unbiased_p = hydra::Parameter::Create("Omega3_2016_unbiased").Value(Omega_2016_unbiased[2]).Error(OmegaErr_2016_unbiased[2]).Fixed();
auto Omega4_2016_unbiased_p = hydra::Parameter::Create("Omega4_2016_unbiased").Value(Omega_2016_unbiased[3]).Error(OmegaErr_2016_unbiased[3]).Fixed();
auto Omega5_2016_unbiased_p = hydra::Parameter::Create("Omega5_2016_unbiased").Value(Omega_2016_unbiased[4]).Error(OmegaErr_2016_unbiased[4]).Fixed();
auto Omega6_2016_unbiased_p = hydra::Parameter::Create("Omega6_2016_unbiased").Value(Omega_2016_unbiased[5]).Error(OmegaErr_2016_unbiased[5]).Fixed();
auto Omega7_2016_unbiased_p = hydra::Parameter::Create("Omega7_2016_unbiased").Value(Omega_2016_unbiased[6]).Error(OmegaErr_2016_unbiased[6]).Fixed();
auto Omega8_2016_unbiased_p = hydra::Parameter::Create("Omega8_2016_unbiased").Value(Omega_2016_unbiased[7]).Error(OmegaErr_2016_unbiased[7]).Fixed();
auto Omega9_2016_unbiased_p = hydra::Parameter::Create("Omega9_2016_unbiased").Value(Omega_2016_unbiased[8]).Error(OmegaErr_2016_unbiased[8]).Fixed();
auto Omega10_2016_unbiased_p = hydra::Parameter::Create("Omega10_2016_unbiased").Value(Omega_2016_unbiased[9]).Error(OmegaErr_2016_unbiased[9]).Fixed();

auto Omega1_2016_biased_p = hydra::Parameter::Create("Omega1_2016_biased").Value(Omega_2016_biased[0]).Error(OmegaErr_2016_biased[0]).Fixed();
auto Omega2_2016_biased_p = hydra::Parameter::Create("Omega2_2016_biased").Value(Omega_2016_biased[1]).Error(OmegaErr_2016_biased[1]).Fixed();
auto Omega3_2016_biased_p = hydra::Parameter::Create("Omega3_2016_biased").Value(Omega_2016_biased[2]).Error(OmegaErr_2016_biased[2]).Fixed();
auto Omega4_2016_biased_p = hydra::Parameter::Create("Omega4_2016_biased").Value(Omega_2016_biased[3]).Error(OmegaErr_2016_biased[3]).Fixed();
auto Omega5_2016_biased_p = hydra::Parameter::Create("Omega5_2016_biased").Value(Omega_2016_biased[4]).Error(OmegaErr_2016_biased[4]).Fixed();
auto Omega6_2016_biased_p = hydra::Parameter::Create("Omega6_2016_biased").Value(Omega_2016_biased[5]).Error(OmegaErr_2016_biased[5]).Fixed();
auto Omega7_2016_biased_p = hydra::Parameter::Create("Omega7_2016_biased").Value(Omega_2016_biased[6]).Error(OmegaErr_2016_biased[6]).Fixed();
auto Omega8_2016_biased_p = hydra::Parameter::Create("Omega8_2016_biased").Value(Omega_2016_biased[7]).Error(OmegaErr_2016_biased[7]).Fixed();
auto Omega9_2016_biased_p = hydra::Parameter::Create("Omega9_2016_biased").Value(Omega_2016_biased[8]).Error(OmegaErr_2016_biased[8]).Fixed();
auto Omega10_2016_biased_p = hydra::Parameter::Create("Omega10_2016_biased").Value(Omega_2016_biased[9]).Error(OmegaErr_2016_biased[9]).Fixed();

hydra::Parameter ExpParams_2015_unbiased_S1[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2015_unbiased_p, Omega2_2015_unbiased_p,
												   Omega3_2015_unbiased_p, Omega4_2015_unbiased_p,
												   Omega5_2015_unbiased_p, Omega6_2015_unbiased_p,
												   Omega7_2015_unbiased_p, Omega8_2015_unbiased_p,
												   Omega9_2015_unbiased_p, Omega10_2015_unbiased_p,
								  			       c0_2015_unbiased_p, c1_2015_unbiased_p, c2_2015_unbiased_p,
								  			       c3_2015_unbiased_p, c4_2015_unbiased_p, c5_2015_unbiased_p,
								  			       c6_2015_unbiased_p, c7_2015_unbiased_p, c8_2015_unbiased_p,
												   Csp1_p};

hydra::Parameter ExpParams_2015_unbiased_S2[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2015_unbiased_p, Omega2_2015_unbiased_p,
												   Omega3_2015_unbiased_p, Omega4_2015_unbiased_p,
												   Omega5_2015_unbiased_p, Omega6_2015_unbiased_p,
												   Omega7_2015_unbiased_p, Omega8_2015_unbiased_p,
												   Omega9_2015_unbiased_p, Omega10_2015_unbiased_p,
								  			       c0_2015_unbiased_p, c1_2015_unbiased_p, c2_2015_unbiased_p,
								  			       c3_2015_unbiased_p, c4_2015_unbiased_p, c5_2015_unbiased_p,
								  			       c6_2015_unbiased_p, c7_2015_unbiased_p, c8_2015_unbiased_p,
												   Csp2_p};

hydra::Parameter ExpParams_2015_unbiased_S3[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2015_unbiased_p, Omega2_2015_unbiased_p,
												   Omega3_2015_unbiased_p, Omega4_2015_unbiased_p,
												   Omega5_2015_unbiased_p, Omega6_2015_unbiased_p,
												   Omega7_2015_unbiased_p, Omega8_2015_unbiased_p,
												   Omega9_2015_unbiased_p, Omega10_2015_unbiased_p,
								  			       c0_2015_unbiased_p, c1_2015_unbiased_p, c2_2015_unbiased_p,
								  			       c3_2015_unbiased_p, c4_2015_unbiased_p, c5_2015_unbiased_p,
								  			       c6_2015_unbiased_p, c7_2015_unbiased_p, c8_2015_unbiased_p,
												   Csp3_p};

hydra::Parameter ExpParams_2015_unbiased_S4[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2015_unbiased_p, Omega2_2015_unbiased_p,
												   Omega3_2015_unbiased_p, Omega4_2015_unbiased_p,
												   Omega5_2015_unbiased_p, Omega6_2015_unbiased_p,
												   Omega7_2015_unbiased_p, Omega8_2015_unbiased_p,
												   Omega9_2015_unbiased_p, Omega10_2015_unbiased_p,
								  			       c0_2015_unbiased_p, c1_2015_unbiased_p, c2_2015_unbiased_p,
								  			       c3_2015_unbiased_p, c4_2015_unbiased_p, c5_2015_unbiased_p,
								  			       c6_2015_unbiased_p, c7_2015_unbiased_p, c8_2015_unbiased_p,
												   Csp4_p};

hydra::Parameter ExpParams_2015_unbiased_S5[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2015_unbiased_p, Omega2_2015_unbiased_p,
												   Omega3_2015_unbiased_p, Omega4_2015_unbiased_p,
												   Omega5_2015_unbiased_p, Omega6_2015_unbiased_p,
												   Omega7_2015_unbiased_p, Omega8_2015_unbiased_p,
												   Omega9_2015_unbiased_p, Omega10_2015_unbiased_p,
								  			       c0_2015_unbiased_p, c1_2015_unbiased_p, c2_2015_unbiased_p,
								  			       c3_2015_unbiased_p, c4_2015_unbiased_p, c5_2015_unbiased_p,
								  			       c6_2015_unbiased_p, c7_2015_unbiased_p, c8_2015_unbiased_p,
												   Csp5_p};

hydra::Parameter ExpParams_2015_unbiased_S6[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2015_unbiased_p, Omega2_2015_unbiased_p,
												   Omega3_2015_unbiased_p, Omega4_2015_unbiased_p,
												   Omega5_2015_unbiased_p, Omega6_2015_unbiased_p,
												   Omega7_2015_unbiased_p, Omega8_2015_unbiased_p,
												   Omega9_2015_unbiased_p, Omega10_2015_unbiased_p,
								  			       c0_2015_unbiased_p, c1_2015_unbiased_p, c2_2015_unbiased_p,
								  			       c3_2015_unbiased_p, c4_2015_unbiased_p, c5_2015_unbiased_p,
								  			       c6_2015_unbiased_p, c7_2015_unbiased_p, c8_2015_unbiased_p,
												   Csp6_p};

hydra::Parameter ExpParams_2015_biased_S1[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2015_biased_p, Omega2_2015_biased_p,
											     Omega3_2015_biased_p, Omega4_2015_biased_p,
											     Omega5_2015_biased_p, Omega6_2015_biased_p,
											     Omega7_2015_biased_p, Omega8_2015_biased_p,
											     Omega9_2015_biased_p, Omega10_2015_biased_p,
								  			  	 c0_2015_biased_p, c1_2015_biased_p, c2_2015_biased_p,
								  			  	 c3_2015_biased_p, c4_2015_biased_p, c5_2015_biased_p,
								  			  	 c6_2015_biased_p, c7_2015_biased_p, c8_2015_biased_p,
											  	 Csp1_p};

hydra::Parameter ExpParams_2015_biased_S2[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2015_biased_p, Omega2_2015_biased_p,
											     Omega3_2015_biased_p, Omega4_2015_biased_p,
											     Omega5_2015_biased_p, Omega6_2015_biased_p,
											     Omega7_2015_biased_p, Omega8_2015_biased_p,
											     Omega9_2015_biased_p, Omega10_2015_biased_p,
								  			     c0_2015_biased_p, c1_2015_biased_p, c2_2015_biased_p,
								  			     c3_2015_biased_p, c4_2015_biased_p, c5_2015_biased_p,
								  			     c6_2015_biased_p, c7_2015_biased_p, c8_2015_biased_p,
											     Csp2_p};

hydra::Parameter ExpParams_2015_biased_S3[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2015_biased_p, Omega2_2015_biased_p,
											     Omega3_2015_biased_p, Omega4_2015_biased_p,
											     Omega5_2015_biased_p, Omega6_2015_biased_p,
											     Omega7_2015_biased_p, Omega8_2015_biased_p,
											     Omega9_2015_biased_p, Omega10_2015_biased_p,
								  			     c0_2015_biased_p, c1_2015_biased_p, c2_2015_biased_p,
								  			     c3_2015_biased_p, c4_2015_biased_p, c5_2015_biased_p,
								  			     c6_2015_biased_p, c7_2015_biased_p, c8_2015_biased_p,
											     Csp3_p};

hydra::Parameter ExpParams_2015_biased_S4[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2015_biased_p, Omega2_2015_biased_p,
											     Omega3_2015_biased_p, Omega4_2015_biased_p,
											     Omega5_2015_biased_p, Omega6_2015_biased_p,
											     Omega7_2015_biased_p, Omega8_2015_biased_p,
											     Omega9_2015_biased_p, Omega10_2015_biased_p,
								  			     c0_2015_biased_p, c1_2015_biased_p, c2_2015_biased_p,
								  			     c3_2015_biased_p, c4_2015_biased_p, c5_2015_biased_p,
								  			     c6_2015_biased_p, c7_2015_biased_p, c8_2015_biased_p,
											     Csp4_p};

hydra::Parameter ExpParams_2015_biased_S5[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2015_biased_p, Omega2_2015_biased_p,
											     Omega3_2015_biased_p, Omega4_2015_biased_p,
											     Omega5_2015_biased_p, Omega6_2015_biased_p,
											     Omega7_2015_biased_p, Omega8_2015_biased_p,
											     Omega9_2015_biased_p, Omega10_2015_biased_p,
								  			     c0_2015_biased_p, c1_2015_biased_p, c2_2015_biased_p,
								  			     c3_2015_biased_p, c4_2015_biased_p, c5_2015_biased_p,
								  			     c6_2015_biased_p, c7_2015_biased_p, c8_2015_biased_p,
											     Csp5_p};

hydra::Parameter ExpParams_2015_biased_S6[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2015_biased_p, Omega2_2015_biased_p,
											     Omega3_2015_biased_p, Omega4_2015_biased_p,
											     Omega5_2015_biased_p, Omega6_2015_biased_p,
											     Omega7_2015_biased_p, Omega8_2015_biased_p,
											     Omega9_2015_biased_p, Omega10_2015_biased_p,
								  			     c0_2015_biased_p, c1_2015_biased_p, c2_2015_biased_p,
								  			     c3_2015_biased_p, c4_2015_biased_p, c5_2015_biased_p,
								  			     c6_2015_biased_p, c7_2015_biased_p, c8_2015_biased_p,
											     Csp6_p};

hydra::Parameter ExpParams_2016_unbiased_S1[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2016_unbiased_p, Omega2_2016_unbiased_p,
												   Omega3_2016_unbiased_p, Omega4_2016_unbiased_p,
												   Omega5_2016_unbiased_p, Omega6_2016_unbiased_p,
												   Omega7_2016_unbiased_p, Omega8_2016_unbiased_p,
												   Omega9_2016_unbiased_p, Omega10_2016_unbiased_p,
								  			       c0_2016_unbiased_p, c1_2016_unbiased_p, c2_2016_unbiased_p,
								  			       c3_2016_unbiased_p, c4_2016_unbiased_p, c5_2016_unbiased_p,
								  			       c6_2016_unbiased_p, c7_2016_unbiased_p, c8_2016_unbiased_p,
												   Csp1_p};

hydra::Parameter ExpParams_2016_unbiased_S2[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2016_unbiased_p, Omega2_2016_unbiased_p,
												   Omega3_2016_unbiased_p, Omega4_2016_unbiased_p,
												   Omega5_2016_unbiased_p, Omega6_2016_unbiased_p,
												   Omega7_2016_unbiased_p, Omega8_2016_unbiased_p,
												   Omega9_2016_unbiased_p, Omega10_2016_unbiased_p,
								  			       c0_2016_unbiased_p, c1_2016_unbiased_p, c2_2016_unbiased_p,
								  			       c3_2016_unbiased_p, c4_2016_unbiased_p, c5_2016_unbiased_p,
								  			       c6_2016_unbiased_p, c7_2016_unbiased_p, c8_2016_unbiased_p,
												   Csp2_p};

hydra::Parameter ExpParams_2016_unbiased_S3[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2016_unbiased_p, Omega2_2016_unbiased_p,
												   Omega3_2016_unbiased_p, Omega4_2016_unbiased_p,
												   Omega5_2016_unbiased_p, Omega6_2016_unbiased_p,
												   Omega7_2016_unbiased_p, Omega8_2016_unbiased_p,
												   Omega9_2016_unbiased_p, Omega10_2016_unbiased_p,
								  			       c0_2016_unbiased_p, c1_2016_unbiased_p, c2_2016_unbiased_p,
								  			       c3_2016_unbiased_p, c4_2016_unbiased_p, c5_2016_unbiased_p,
								  			       c6_2016_unbiased_p, c7_2016_unbiased_p, c8_2016_unbiased_p,
												   Csp3_p};

hydra::Parameter ExpParams_2016_unbiased_S4[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2016_unbiased_p, Omega2_2016_unbiased_p,
												   Omega3_2016_unbiased_p, Omega4_2016_unbiased_p,
												   Omega5_2016_unbiased_p, Omega6_2016_unbiased_p,
												   Omega7_2016_unbiased_p, Omega8_2016_unbiased_p,
												   Omega9_2016_unbiased_p, Omega10_2016_unbiased_p,
								  			       c0_2016_unbiased_p, c1_2016_unbiased_p, c2_2016_unbiased_p,
								  			       c3_2016_unbiased_p, c4_2016_unbiased_p, c5_2016_unbiased_p,
								  			       c6_2016_unbiased_p, c7_2016_unbiased_p, c8_2016_unbiased_p,
												   Csp4_p};

hydra::Parameter ExpParams_2016_unbiased_S5[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2016_unbiased_p, Omega2_2016_unbiased_p,
												   Omega3_2016_unbiased_p, Omega4_2016_unbiased_p,
												   Omega5_2016_unbiased_p, Omega6_2016_unbiased_p,
												   Omega7_2016_unbiased_p, Omega8_2016_unbiased_p,
												   Omega9_2016_unbiased_p, Omega10_2016_unbiased_p,
								  			       c0_2016_unbiased_p, c1_2016_unbiased_p, c2_2016_unbiased_p,
								  			       c3_2016_unbiased_p, c4_2016_unbiased_p, c5_2016_unbiased_p,
								  			       c6_2016_unbiased_p, c7_2016_unbiased_p, c8_2016_unbiased_p,
												   Csp6_p};

hydra::Parameter ExpParams_2016_unbiased_S6[32] = {b0_p, b1_p,
                         		  			       p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			       p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			       Omega1_2016_unbiased_p, Omega2_2016_unbiased_p,
												   Omega3_2016_unbiased_p, Omega4_2016_unbiased_p,
												   Omega5_2016_unbiased_p, Omega6_2016_unbiased_p,
												   Omega7_2016_unbiased_p, Omega8_2016_unbiased_p,
												   Omega9_2016_unbiased_p, Omega10_2016_unbiased_p,
								  			       c0_2016_unbiased_p, c1_2016_unbiased_p, c2_2016_unbiased_p,
								  			       c3_2016_unbiased_p, c4_2016_unbiased_p, c5_2016_unbiased_p,
								  			       c6_2016_unbiased_p, c7_2016_unbiased_p, c8_2016_unbiased_p,
												   Csp6_p};

hydra::Parameter ExpParams_2016_biased_S1[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2016_biased_p, Omega2_2016_biased_p,
											     Omega3_2016_biased_p, Omega4_2016_biased_p,
											     Omega5_2016_biased_p, Omega6_2016_biased_p,
											     Omega7_2016_biased_p, Omega8_2016_biased_p,
											     Omega9_2016_biased_p, Omega10_2016_biased_p,
								  			     c0_2016_biased_p, c1_2016_biased_p, c2_2016_biased_p,
								  			     c3_2016_biased_p, c4_2016_biased_p, c5_2016_biased_p,
								  			     c6_2016_biased_p, c7_2016_biased_p, c8_2016_biased_p,
											     Csp1_p};

hydra::Parameter ExpParams_2016_biased_S2[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2016_biased_p, Omega2_2016_biased_p,
											     Omega3_2016_biased_p, Omega4_2016_biased_p,
											     Omega5_2016_biased_p, Omega6_2016_biased_p,
											     Omega7_2016_biased_p, Omega8_2016_biased_p,
											     Omega9_2016_biased_p, Omega10_2016_biased_p,
								  			     c0_2016_biased_p, c1_2016_biased_p, c2_2016_biased_p,
								  			     c3_2016_biased_p, c4_2016_biased_p, c5_2016_biased_p,
								  			     c6_2016_biased_p, c7_2016_biased_p, c8_2016_biased_p,
											     Csp2_p};

hydra::Parameter ExpParams_2016_biased_S3[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2016_biased_p, Omega2_2016_biased_p,
											     Omega3_2016_biased_p, Omega4_2016_biased_p,
											     Omega5_2016_biased_p, Omega6_2016_biased_p,
											     Omega7_2016_biased_p, Omega8_2016_biased_p,
											     Omega9_2016_biased_p, Omega10_2016_biased_p,
								  			     c0_2016_biased_p, c1_2016_biased_p, c2_2016_biased_p,
								  			     c3_2016_biased_p, c4_2016_biased_p, c5_2016_biased_p,
								  			     c6_2016_biased_p, c7_2016_biased_p, c8_2016_biased_p,
											     Csp3_p};

hydra::Parameter ExpParams_2016_biased_S4[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2016_biased_p, Omega2_2016_biased_p,
											     Omega3_2016_biased_p, Omega4_2016_biased_p,
											     Omega5_2016_biased_p, Omega6_2016_biased_p,
											     Omega7_2016_biased_p, Omega8_2016_biased_p,
											     Omega9_2016_biased_p, Omega10_2016_biased_p,
								  			     c0_2016_biased_p, c1_2016_biased_p, c2_2016_biased_p,
								  			     c3_2016_biased_p, c4_2016_biased_p, c5_2016_biased_p,
								  			     c6_2016_biased_p, c7_2016_biased_p, c8_2016_biased_p,
											     Csp4_p};

hydra::Parameter ExpParams_2016_biased_S5[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2016_biased_p, Omega2_2016_biased_p,
											     Omega3_2016_biased_p, Omega4_2016_biased_p,
											     Omega5_2016_biased_p, Omega6_2016_biased_p,
											     Omega7_2016_biased_p, Omega8_2016_biased_p,
											     Omega9_2016_biased_p, Omega10_2016_biased_p,
								  			     c0_2016_biased_p, c1_2016_biased_p, c2_2016_biased_p,
								  			     c3_2016_biased_p, c4_2016_biased_p, c5_2016_biased_p,
								  			     c6_2016_biased_p, c7_2016_biased_p, c8_2016_biased_p,
											     Csp5_p};

hydra::Parameter ExpParams_2016_biased_S6[32] = {b0_p, b1_p,
                         		  			     p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  			     p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  			     Omega1_2016_biased_p, Omega2_2016_biased_p,
											     Omega3_2016_biased_p, Omega4_2016_biased_p,
											     Omega5_2016_biased_p, Omega6_2016_biased_p,
											     Omega7_2016_biased_p, Omega8_2016_biased_p,
											     Omega9_2016_biased_p, Omega10_2016_biased_p,
								  			     c0_2016_biased_p, c1_2016_biased_p, c2_2016_biased_p,
								  			     c3_2016_biased_p, c4_2016_biased_p, c5_2016_biased_p,
								  			     c6_2016_biased_p, c7_2016_biased_p, c8_2016_biased_p,
											     Csp6_p};

std::vector<double> parameters = {A02, Aperp2, FS1,
                                  deltagammasd, deltagammas, deltams,
                                  phi0,    phipar0,    phiperp0,    phiS0,
                                  lambda0, lambdapar0, lambdaperp0, lambdaS0,
                                  deltapar0,  deltaperp0,  deltaS1perp,
                                  b0, b1,
                         		  p0_OS, p1_OS, DeltaP0_OS, DeltaP1_OS, AvgEta_OS,
                         		  p0_SS, p1_SS, DeltaP0_SS, DeltaP1_SS, AvgEta_SS,
                        		  Omega_2016_biased[0], Omega_2016_biased[1], Omega_2016_biased[2], Omega_2016_biased[3], Omega_2016_biased[4],
                         		  Omega_2016_biased[5], Omega_2016_biased[6], Omega_2016_biased[7], Omega_2016_biased[8], Omega_2016_biased[9],
								  Spline_Coeffs_2016_biased[0], Spline_Coeffs_2016_biased[1], Spline_Coeffs_2016_biased[2],
                                  Spline_Coeffs_2016_biased[3], Spline_Coeffs_2016_biased[4], Spline_Coeffs_2016_biased[5],
                                  Spline_Coeffs_2016_biased[6], Spline_Coeffs_2016_biased[7], Spline_Coeffs_2016_biased[8],
								  Csp[0]};


//---------------------------------------------------------
// 		Structs for A_k, B_k, C_k, D_k, N_k
//  	coefficients and angular functions f_k
//---------------------------------------------------------

namespace medusa {
	namespace parameters {


		struct AngularTimeCoefficients
		{
			double k[10];
		};


		struct NFactors
		{
			double k[10];
		};


		struct AngularFunctions
		{
			// ctor of the angular functions f_k(costheta_h, costheta_l, phi)
			__hydra_dual__
			AngularFunctions(double costheta_h, double costheta_l, double phi)
			{
				const double ck = costheta_h;
				const double sk2 = 1 - costheta_h*costheta_h;
				const double sk = ::sqrt(sk2);
				const double cp = ::cos(phi);
				const double sp = ::sin(phi);
				const double cl = costheta_l;
				const double sl2 = 1 - costheta_l*costheta_l;
				const double sl = ::sqrt(sl2);

				//::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2)
				fk[0] = ck*ck*sl2;

				//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::cos(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
				fk[1] = 0.5*sk2*(1.0 - cp*cp * sl2);

				//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::sin(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
				fk[2] = 0.5*sk2*(1.0 - sp*sp * sl2);

				//::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(phi) * ::cos(phi)
				fk[3] = sk2*sl2*sp*cp;

				//sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
				fk[4] = M_Sqrt2* sk * ck * sl * cl * cp ;

				//-sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
				fk[5] = -M_Sqrt2* sk * ck * sl * cl * sp;

				//1./3. * ::pow( ::sin(theta_l) , 2 )
				fk[6] = M_1_3*sl2;

				//2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
				fk[7] = M_2_Sqrt6 * sk * sl * cl * cp;

				//-2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
				fk[8] = -M_2_Sqrt6* sk * sl * cl * sp;

				//2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 )
				fk[9] = M_2_Sqrt3 * ck * sl2;
			}

			double fk[10];

		};


	} // namespace parameters
}  // namespace medusa

#endif /* PHIS_PARAMETERS_H_ */
