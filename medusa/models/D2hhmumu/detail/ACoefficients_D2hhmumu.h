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
 * ACoefficients_D2hhmumu.h
 *
 *  Created on: 24/07/2020
 *      Author: augalves
 */

#ifndef ACOEFFICIENTS_H_
#define ACOEFFICIENTS_H_

#include <hydra/detail/Config.h>

namespace medusa {

namespace detail {

struct ACoefficients_D2hhmumu
{

	__hydra_host__ __hydra_device__
	ACoefficients_D2hhmumu(double const& theta_l, double const& phi )
	{
		const double cp = ::cos(phi);
		const double sp = ::sin(phi);
		const double cl = ::cos(theta_l);
		const double sl = ::sin(theta_l);
		

		//::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2)
		fC[0] = 1.0;

		//::cos(2*theta_l);
		fC[1] = cl*cl - sl*sl;

		//::pow(::sin(theta_l),2) * ::cos(2*chi);
		fC[2] = sl*sl*(cp*cp-sp*sp);

		//::sin(2*theta_l)*::cos(chi);
		fC[3] = 2.0*sl*cl*cp;

		//::sin(theta_l)*::cos(chi);
		fC[4] = sl*cp ;

		//::cos(theta_l);
		fC[5] = cl;

		//::sin(theta_l)*::sin(chi);
		fC[6] = sl*sp;

		//::sin(2*theta_l)*::sin(chi);
		fC[7] = 2.0*sl*cl*sp;

		//::pow(::sin(theta_l),2)*::sin(2*chi);
		fC[8] = 2.0*sp*cp*sl*sl;

	}

	double fC[9];

};



} // namespace medusa::detail


}  // namespace medusa




#endif /* ACOEFFICIENTS_H_ */
