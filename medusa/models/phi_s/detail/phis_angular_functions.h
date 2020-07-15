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

struct __hydra_align__(16) AngularArgs {

	AngularArgs(double theta_h, double theta_l, double phi ) :
		cx(::cos(theta_h)),	sx(::sin(theta_h)),
		cz(::cos(phi)),	sz(::sin(phi)),
		cy(::cos(theta_l)),	sy(::sin(theta_l))
	{}

	double cx ;
	double sx ;
	double cz ;
	double sz ;
	double cy ;
	double sy ;
};


  __hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<1>, AngularArgs const& aa){

	//::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2)

   double cx = aa.cx*aa.cx;
   double sy = aa.sy*aa.sy;

	return cx*sy;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<2>, AngularArgs const& aa){

	//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::cos(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )

	double sx = aa.sx*aa.sx;
	double sz = aa.sz*aa.sz;
	double sy = aa.sy*aa.sy;

	return 0.5 * sx * sy * sz;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<3>, AngularArgs const& aa){

	//0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::sin(phi) , 2) *  ::pow( ::sin(theta_l) , 2) )

	double sx = aa.sx*aa.sx;
	double sy = aa.sy*aa.sy;
	double cz = aa.cz*aa.cz;

	return 0.5 * sx * sy * cz;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<4>, AngularArgs const& aa){

	//::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(phi) * ::cos(phi)
	double sx = aa.sx*aa.sx;
	double sy = aa.sy*aa.sy;

   return  sx * sy * aa.cz * aa.sz ;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<5>, AngularArgs const& aa){

	//sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)

	const static  double f = 1.414213562373095; //sqrt2

	return f * aa.sx * aa.cx * aa.sy * aa.cy * aa.cz ;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<6>, AngularArgs const& aa){

	//-sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)

	const static  double f = 1.414213562373095; //sqrt2

	return -f * aa.sx * aa.sy * aa.cy * aa.sz ;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<7>, AngularArgs const& aa){

	//1./3. * ::pow( ::sin(theta_l) , 2 )

	const static  double f = 0.333333333333333; //1./3.

	return  f * aa.sy*aa.sy;
}


__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<8>, AngularArgs const& aa){

	//2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi)

	const static  double f = 0.816496580927726; //2.0/sqrt6

	return f * aa.sx * aa.sy * aa.cy * aa.cz;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<9>, AngularArgs const& aa){

	//-2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi)

	const static  double f = 0.816496580927726; //2.0/sqrt6

	return -f * aa.sx * aa.sy * aa.cy * aa.sz;
}

__hydra_host__ __hydra_device__
inline double phis_angular_functions(hydra::placeholders::placeholder<10>, AngularArgs const& aa){

	//2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 )

	const static  double f = 3.464101615137755; //2./sqrt3

	return f * aa.cx * aa.sy * aa.sy;
}
    
} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_ANGULAR_FUNCTIONS_H_ */
