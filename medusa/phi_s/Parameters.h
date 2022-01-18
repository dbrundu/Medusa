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

// particle types
declarg(Jpsi,  hydra::Vector4R)
declarg(Phi,   hydra::Vector4R)
declarg(KaonP, hydra::Vector4R)
declarg(KaonM, hydra::Vector4R)
declarg(MuonP, hydra::Vector4R)
declarg(MuonM, hydra::Vector4R)

// new types for time and helicity angles
declarg(theta_h_t,  double)
declarg(theta_l_t,  double)
declarg(phi_t,      double)
declarg(dtime_t,    double)

// new types for tagging variables
declarg(qOS_t, 	 int)
declarg(qSS_t, 	 int)
declarg(etaOS_t, double)
declarg(etaSS_t, double)
declarg(delta_t, double)


//-------------------------------------
// Model parameters for:
//  B0s -> J/psi  (Phi -> K+ K-)
//          |-> mu+ mu-
//-------------------------------------

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
const double deltapar   = 3.030;             // + delta0
const double deltaperp  = 2.60;              // + delta0
const double deltaS     = -0.30;             // + delta0

const double deltagammasd = -0.0044;
const double deltagammas  = 0.0782;
const double deltams      = 17.713;

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


//-------------------------------------
// Experimental parameters for:
//  B0s -> J/psi  (Phi -> K+ K-)
//          |-> mu+ mu-
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

const double Spline_Knots[7] = {0.3, 0.58, 0.91,
								1.35, 1.96, 3.01,
												7.00};

double Spline_coeffs[9] = {1.0, 1.49, 2.06,
							2.12, 2.28, 2.29,
							2.46, 2.25, 2.34};

const double Omega[10] = {1.0, 1.0336, 1.0336,
						  0.0028, 0.00298, -0.0002,
							1.0196, 0.00019, 0.00019,
													0.0057};

const double OmegaErr[10] = {0.0, 0.0015, 0.0015,
							 0.0013, 0.00074, 0.00072,
							 0.0011, 0.00094, 0.00094,
							 						0.0019};

const double Csp[6] = {0.8463, 0.8756, 0.8478, 0.8833, 0.9415, 0.9756};

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

auto Spline_c0_p = hydra::Parameter::Create("c0").Value(Spline_coeffs[0]).Error(0.0).Fixed();
auto Spline_c1_p = hydra::Parameter::Create("c1").Value(Spline_coeffs[1]).Error(0.16).Fixed();
auto Spline_c2_p = hydra::Parameter::Create("c2").Value(Spline_coeffs[2]).Error(0.13).Fixed();
auto Spline_c3_p = hydra::Parameter::Create("c3").Value(Spline_coeffs[3]).Error(0.17).Fixed();
auto Spline_c4_p = hydra::Parameter::Create("c4").Value(Spline_coeffs[4]).Error(0.16).Fixed();
auto Spline_c5_p = hydra::Parameter::Create("c5").Value(Spline_coeffs[5]).Error(0.17).Fixed();
auto Spline_c6_p = hydra::Parameter::Create("c6").Value(Spline_coeffs[6]).Error(0.19).Fixed();
auto Spline_c7_p = hydra::Parameter::Create("c7").Value(Spline_coeffs[7]).Error(0.18).Fixed();
auto Spline_c8_p = hydra::Parameter::Create("c8").Value(Spline_coeffs[8]).Error(0.17).Fixed();

auto Omega_1_p = hydra::Parameter::Create("Omega_1").Value(Omega[0]).Error(OmegaErr[0]).Fixed();
auto Omega_2_p = hydra::Parameter::Create("Omega_2").Value(Omega[1]).Error(OmegaErr[1]).Fixed();
auto Omega_3_p = hydra::Parameter::Create("Omega_3").Value(Omega[2]).Error(OmegaErr[2]).Fixed();
auto Omega_4_p = hydra::Parameter::Create("Omega_4").Value(Omega[3]).Error(OmegaErr[3]).Fixed();
auto Omega_5_p = hydra::Parameter::Create("Omega_5").Value(Omega[4]).Error(OmegaErr[4]).Fixed();
auto Omega_6_p = hydra::Parameter::Create("Omega_6").Value(Omega[5]).Error(OmegaErr[5]).Fixed();
auto Omega_7_p = hydra::Parameter::Create("Omega_7").Value(Omega[6]).Error(OmegaErr[6]).Fixed();
auto Omega_8_p = hydra::Parameter::Create("Omega_8").Value(Omega[7]).Error(OmegaErr[7]).Fixed();
auto Omega_9_p = hydra::Parameter::Create("Omega_9").Value(Omega[8]).Error(OmegaErr[8]).Fixed();
auto Omega_10_p = hydra::Parameter::Create("Omega_10").Value(Omega[9]).Error(OmegaErr[9]).Fixed();

hydra::Parameter ExpParams[31] = {b0_p, b1_p,
                         		  p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AvgEta_OS_p,
                         		  p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AvgEta_SS_p,
                        		  Omega_1_p, Omega_2_p, Omega_3_p, Omega_4_p, Omega_5_p,
                         		  Omega_6_p, Omega_7_p, Omega_8_p, Omega_9_p, Omega_10_p,
								  Spline_c0_p, Spline_c1_p, Spline_c2_p,
								  Spline_c3_p, Spline_c4_p, Spline_c5_p,
								  Spline_c6_p, Spline_c7_p, Spline_c8_p};

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
			// ctor of the angular functions f_k(theta_h, theta_l, phi)
			__hydra_dual__
			AngularFunctions(double theta_h, double theta_l, double phi)
			{
				const double cx = ::cos(theta_h);
				const double sx = ::sin(theta_h);
				const double cz = ::cos(phi);
				const double sz = ::sin(phi);
				const double cy = ::cos(theta_l);
				const double sy = ::sin(theta_l);

				//::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2)
				fk[0] = cx*sy; fk[0]*= fk[0];

				//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::cos(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
				fk[1] = 1.0 - cz*cz * sy*sy; fk[1] *= 0.5*sx*sx;

				//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::sin(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
				fk[2] = 1.0 - sz*sz * sy*sy; fk[2] *= 0.5*sx*sx;

				//::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(phi) * ::cos(phi)
				fk[3] = sx * sy; fk[3]*=cz*sz*fk[3];

				//sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
				fk[4] = M_Sqrt2* sx * cx * sy * cy * cz ;

				//-sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
				fk[5] = -M_Sqrt2* sx *  cx * sy * cy * sz;

				//1./3. * ::pow( ::sin(theta_l) , 2 )
				fk[6] = M_1_3*sy*sy;

				//2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
				fk[7] = M_2_Sqrt6 * sx * sy * cy * cz;

				//-2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
				fk[8] = -M_2_Sqrt6* sx * sy * cy * sz;

				//2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 )
				fk[9] = M_2_Sqrt3 * cx * sy * sy;
			}

			double fk[10];

		};


	} // namespace parameters
}  // namespace medusa

#endif /* PHIS_PARAMETERS_H_ */
