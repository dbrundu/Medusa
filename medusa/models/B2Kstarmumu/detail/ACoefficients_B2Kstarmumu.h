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
 * ACoefficients_B2Kstarmumu.h
 *
 *  Created on: 18/10/2020
 *      Author: dbrundu
 */

#ifndef ACOEFFICIENTS_B2KSTARMUMU_H_
#define ACOEFFICIENTS_B2KSTARMUMU_H_

#include <hydra/detail/Config.h>

namespace medusa {

namespace detail {

struct ACoefficients_B2Kstarmumu
{

	__hydra_host__ __hydra_device__
	ACoefficients_B2Kstarmumu(double theta_h, double theta_l, double phi )
	{
		const double ch = ::cos(theta_h);
		const double sh = ::sin(theta_h);
		const double cp = ::cos(phi);
		const double sp = ::sin(phi);
		
		const double cl = ::cos(theta_l);
		sl = ::sin(theta_l);

		s2l = 2.0*sl*cl;
		c2l = cl*cl - sl*sl;
		
		const double s2h = 2.0*sh*ch;
		const double c2h = ch*ch - sh*sh;
		const double s2p = 2.0*sp*cp;
		const double c2p = cp*cp - sp*sp;


		//::pow( ::sin(theta_h) , 2)
		fC[0] = sh*sh;

		//::pow( ::cos(theta_h) , 2)
		fC[1] = ch*ch;

		//::pow( ::sin(theta_h) , 2) * ::cos(2*theta_l)
		fC[2] = sl*sl*c2l;

		//::pow( ::cos(theta_h) , 2) * ::cos(2*theta_l)
		fC[3] = ch*ch*c2l;

		//::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::cos(2*phi)
		fC[4] = sh*sh*sl*sl*c2p ;

		//::sin(2*theta_h) * ::sin(2*theta_l) * ::cos(phi)
		fC[5] = s2h*s2l*cp;

		// ::sin(2*theta_h) * ::sin(theta_l) * ::cos(phi)
		fC[6] = s2h*sl*cp;

		//::pow( ::sin(theta_h) , 2) * ::cos(theta_l)
		fC[7] = sh*sh*cl;

		//::sin(2*theta_h) * ::sin(theta_l) * ::sin(phi)
		fC[8] = s2h*sl*sp;

        //::sin(2*theta_h) * ::sin(2*theta_l) * ::sin(phi)
        fC[9] = s2h*s2l*sp;
        
        //::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(2*phi)
        fC[10] = sh*sh*sl*sl*s2p;
        
        //::pow( ::sin(theta_l) , 2)
        fC[11] = sl*sl;
        
        //::cos(theta_h)
        fC[12] = ch;
        
        //::sin(theta_h) * ::cos(phi)
        fC[13] = sh*cp;
        
        //::sin(theta_h) * ::sin(phi)
        fC[14] = sh*sp;
        
        
	}

	double fC[15];
	double sl, s2l, c2l; //--> used also externally

};



} // namespace medusa::detail


}  // namespace medusa




#endif /* ACOEFFICIENTS_B2KSTARMUMU_H_ */
