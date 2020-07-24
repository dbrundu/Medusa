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
 * AFactors.h
 *
 *  Created on: 24/07/2020
 *      Author: augalves
 */

#ifndef AFACTORS_H_
#define AFACTORS_H_

#include <hydra/detail/Config.h>

namespace medusa {

namespace detail {

struct AFactors
{

	__hydra_host__ __hydra_device__
	AFactors(double theta_h, double theta_l, double phi )
	{
		const static  double Sqrt2     = 1.414213562373095; //\sqrt{2}
		const static  double OneThird  = 0.333333333333333; //1./3.
		const static  double N2DSqrt6  = 0.816496580927726; //2.0/sqrt6
		const static  double N2DSqrt3  = 3.464101615137755; //2.0/sqrt3

		const double cx = ::cos(theta_h);
		const double sx = ::sin(theta_h);
		const double cz = ::cos(phi);
		const double sz = ::sin(phi);
		const double cy = ::cos(theta_l);
		const double sy = ::sin(theta_l);

		//::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2)
		fA[0] = cx*sy; fA[0]*= fA[0];

		//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::cos(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
		fA[1] = sx * sy * sz; fA[1]*=0.5*fA[1];

		//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::sin(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )
		fA[2] = sx * sy * cz; fA[2]*=0.5*fA[2];

		//::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(phi) * ::cos(phi)
		fA[3] = sx * sy; fA[3]*=cz*sz*fA[3];

		//sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
		fA[4] = Sqrt2* sx * cx * sy * cy * cz ;

		//-sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
		fA[5] = -Sqrt2* sx * sy * cy * sz;

		//1./3. * ::pow( ::sin(theta_l) , 2 )
		fA[6] = OneThird*sy*sy;

		//2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)
		fA[7] = N2DSqrt6 * sx * sy * cy * cz;

		//-2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)
		fA[8] = N2DSqrt6* sx * sy * cy * sz;

		//2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 )
		fA[9] = N2DSqrt3 * cx * sy * sy;
	}

	double fA[10];

};



} // namespace medusa::detail


}  // namespace medusa




#endif /* AFACTORS_H_ */
