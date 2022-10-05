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
/*----------------------------------------
 *  Created on: 20/07/2022
 *
 *  Author: Davide Brundu and
 * 			Alessandro Maria Ricci
 * 
 *  Note: this is a fixed version of the
 *        corresponding hydra function.
 *----------------------------------------*/

#ifndef PLANES_DELTA_ANGLE_H_
#define PLANES_DELTA_ANGLE_H_

#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Tuple.h>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <cmath>

namespace medusa {

/**
 * \ingroup common_functions
 * \class PlanesDeltaAngle
 *
 *  This functor calculates the delta angle between decay plane of the particle with four-vector d2 and d3 (same plane)
 *  and h1 (other plane)
 */
class PlanesDeltaAngle: public hydra::BaseFunctor<PlanesDeltaAngle, double(hydra::Vector4R, hydra::Vector4R, hydra::Vector4R), 0>
{
	public:

	__hydra_dual__
	PlanesDeltaAngle(){};

	// ctor with other PlanesDeltaAngle instance (copy ctor)
	__hydra_dual__
	PlanesDeltaAngle( PlanesDeltaAngle const& other):
		BaseFunctor<PlanesDeltaAngle, double(hydra::Vector4R, hydra::Vector4R, hydra::Vector4R), 0>(other)
	{}

	__hydra_dual__
	inline PlanesDeltaAngle& operator=( PlanesDeltaAngle const& other)
	{
			if(this==&other) return  *this;
			BaseFunctor<PlanesDeltaAngle, double(hydra::Vector4R, hydra::Vector4R, hydra::Vector4R), 0>::operator=(other);
			return  *this;
	}

	__hydra_dual__
	inline double Evaluate(hydra::Vector4R const& d2, hydra::Vector4R const& d3, hydra::Vector4R const& h1) const
	{
		return chi_angle(d2, d3, h1);
	}


	private:

	__hydra_dual__
	inline hydra::GReal_t chi_angle(hydra::Vector4R const& d2, hydra::Vector4R const& d3, hydra::Vector4R const& h1) const
	{
		hydra::Vector4R D = d2 + d3;
     
		hydra::Vector4R d1_perp = d2 - (D.dot(d2) / D.dot(D)) * D;
		hydra::Vector4R h1_perp = h1 - (D.dot(h1) / D.dot(D)) * D;

		// orthogonal to both D and d1_perp
		hydra::Vector4R d1_prime = D.cross(d1_perp);

		d1_perp = d1_perp / d1_perp.d3mag();
		d1_prime = d1_prime / d1_prime.d3mag();

		hydra::GReal_t x, y;

		x = d1_perp.dot(h1_perp);
		y = d1_prime.dot(h1_perp);

		hydra::GReal_t chi = ::atan2(y, x);
      
		// Since any integer multiple of 2*PI can be added to the angle theta without changing either x or y,
		// implying an ambiguous value for the returned value, the principal value of the angle,
		// in the (left open, right closed) interval (−PI, PI] is returned.
		// Chi is signed with counterclockwise angles being positive and clockwise being negative.
		// In other words, ::atan2(y,x) is in the closed interval [0, PI] when y ≥ 0,
		// and in the open interval (−PI, 0) when y < 0.
		return chi;
	}

};

}  // namespace medusa

#endif // PLANES_DELTA_ANGLE_H_