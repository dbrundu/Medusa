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
/*
 * 	Parameters.h
 *
 *  Created on: 04/06/2021
 *      Author: Alessandro Maria Ricci
 * 
 * 	This library contains all parameters used by PhisSignal.h amd PhisComplete.h.
 */


#ifndef PHIS_PARAMETERS_H_
#define PHIS_PARAMETERS_H_

// Hydra
#include <hydra/detail/Config.h>
#include <hydra/Vector4R.h>


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

const double Omega[10] = {1.0, 1.0208, 1.0208,
						  0.0024, 0.00321, -0.00018,
							1.0113, -0.00003, -0.00003,
													-0.0022};

const double OmegaErr[10] = {0, 0.0014, 0.0014,
							 0.0012, 0.00067, 0.00067,
							 0.001, 0.00087, 0.00088,
							 						0.0018};

auto b0_p = hydra::Parameter::Create("b0").Value(b0).Error(0.00022).Fixed();
auto b1_p = hydra::Parameter::Create("b1").Value(b1).Error(0.0057).Fixed();

auto p0_OS_p = hydra::Parameter::Create("p0_OS").Value(p0_OS).Error(0.0007 + 0.0028).Fixed();
auto p1_OS_p = hydra::Parameter::Create("p1_OS").Value(p1_OS).Error(0.0062 + 0.0265).Fixed();
auto DeltaP0_OS_p = hydra::Parameter::Create("DeltaP0_OS").Value(DeltaP0_OS).Error(0.0014).Fixed();
auto DeltaP1_OS_p = hydra::Parameter::Create("DeltaP1_OS").Value(DeltaP1_OS).Error(0.0124).Fixed();
auto AgvEta_OS_p = hydra::Parameter::Create("AvgEta_OS").Value(AvgEta_OS).Error(0).Fixed();

auto p0_SS_p = hydra::Parameter::Create("p0_SS").Value(p0_SS).Error(0.0108 + 0.003).Fixed();
auto p1_SS_p = hydra::Parameter::Create("p1_SS").Value(p1_SS).Error(0.1314 + 0.0196).Fixed();
auto DeltaP0_SS_p = hydra::Parameter::Create("DeltaP0_SS").Value(DeltaP0_SS).Error(0.03).Fixed();
auto DeltaP1_SS_p = hydra::Parameter::Create("DeltaP1_SS").Value(DeltaP1_SS).Error(0.03).Fixed();
auto AgvEta_SS_p = hydra::Parameter::Create("AvgEta_SS").Value(AvgEta_SS).Error(0).Fixed();

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

hydra::Parameter ExpParams[22] = {b0_p, b1_p,
                         		  p0_OS_p, p1_OS_p, DeltaP0_OS_p, DeltaP1_OS_p, AgvEta_OS_p,
                         		  p0_SS_p, p1_SS_p, DeltaP0_SS_p, DeltaP1_SS_p, AgvEta_SS_p,
                        		  Omega_1_p, Omega_2_p, Omega_3_p, Omega_4_p, Omega_5_p,
                         		  Omega_6_p, Omega_7_p, Omega_8_p, Omega_9_p, Omega_10_p};



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
				const static  double Sqrt2     = 1.414213562373095; //	\sqrt{2}
				const static  double OneThird  = 0.333333333333333; //	1./3.
				const static  double N2DSqrt6  = 0.816496580927726; //	2.0/sqrt6
				const static  double N2DSqrt3  = 1.154700538379252; //	2.0/sqrt3

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
				fk[4] = Sqrt2* sx * cx * sy * cy * cz ;

				//-sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
				fk[5] = -Sqrt2* sx *  cx * sy * cy * sz;

				//1./3. * ::pow( ::sin(theta_l) , 2 )
				fk[6] = OneThird*sy*sy;

				//2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
				fk[7] = N2DSqrt6 * sx * sy * cy * cz;

				//-2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
				fk[8] = -N2DSqrt6* sx * sy * cy * sz;

				//2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 )
				fk[9] = N2DSqrt3 * cx * sy * sy;
			}

			double fk[10];

		};


	} // namespace parameters

}  // namespace medusa

#endif /* PHIS_PARAMETERS_H_ */
