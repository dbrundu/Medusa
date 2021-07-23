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
 * 	This library contains all parameters used by PhisSignal.h.
 */


#ifndef MEDUSA_PARAMETERS_H_
#define MEDUSA_PARAMETERS_H_

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
			// ctor of the angular functions fk(omega)
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

#endif /* MEDUSA_PARAMETERS_H_ */
