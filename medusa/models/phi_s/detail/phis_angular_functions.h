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
 *
 *
 *  Created on: 11/05/2020
 *      Author: Davide Brundu
 *
 *  Modified on: 14/07/2020
 *       Author: Antonio Augusto Alves Jr.
 *          Log: Getting rid of the dozens of repeated
 *          and resource wasteful calls to pow, sin, cos.
 *
 */

#ifndef PHIS_ANGULAR_FUNCTIONS_H_
#define PHIS_ANGULAR_FUNCTIONS_H_

#include <hydra/Placeholders.h>



namespace medusa {

namespace detail {

struct __hydra_align__(16) AngularFactors {

	__hydra_dual__
	AngularFactors(double theta_h, double theta_l, double phi ) :
		ch(::cos(theta_h)),   sh(::sin(theta_h)),
		cp(::cos(phi)),       sp(::sin(phi)),
		cl(::cos(theta_l)),	  sl(::sin(theta_l))
	{
		const static  double Sqrt2     = 1.414213562373095; //sqrt{2}
		const static  double OneThird  = 0.333333333333333; //1./3.
		const static  double N2DSqrt6  = 0.816496580927726; //2.0/sqrt6
		const static  double N2DSqrt3  = 1.1547005383792515290; //2.0/sqrt3

		//::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2)
		fA[0] = ch*sl; fA[0] *= fA[0];

		//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::cos(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
		fA[1] = 1.0 - cp*cp * sl*sl; fA[1] *= 0.5*sh*sh;

		//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::sin(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
		fA[2] = 1.0 - sp*sp * sl*sl; fA[2] *= 0.5*sh*sh;

		//::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(phi) * ::cos(phi)
		fA[3] = sh * sl; fA[3] *= cp*sp*fA[3];

		//sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
		fA[4] =  Sqrt2* sh * ch * sl * cl * cp ;

		//-sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
		fA[5] = -Sqrt2 * sh * ch * sl * cl * sp;

		//1./3. * ::pow( ::sin(theta_l) , 2 )
		fA[6] = OneThird * sl * sl;

		//2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
		fA[7] = N2DSqrt6 * sh * sl * cl * cp;

		//-2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
		fA[8]= -N2DSqrt6* sh * sl * cl * sp;

		//2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 )
		fA[9]= N2DSqrt3 * ch * sl * sl;
	}

	double fA[10];

private:
	double ch ;
	double sh ;
	double cp ;
	double sp ;
	double cl ;
	double sl ;
};


} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_ANGULAR_FUNCTIONS_H_ */
