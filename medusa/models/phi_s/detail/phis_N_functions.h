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
 *  Created on: 18/05/2020
 *      Author: Davide Brundu
 */

#ifndef PHIS_N_FUNCTIONS_H_
#define PHIS_N_FUNCTIONS_H_



namespace medusa {

namespace detail {

struct NFactors
{
	__hydra_dual__
	NFactors(double A_0, double A_perp,  double A_S, double A_par):
		fC1(A_0*A_0),
		fC2(A_par*A_par),
		fC3(A_perp*A_perp),
		fC4(A_perp*A_par),
		fC5(A_0*A_par),
		fC6(A_0*A_perp),
		fC7(A_S*A_S),
		fC8(A_S*A_par),
		fC9(A_S*A_perp),
		fC10(A_S*A_0)
		{}

	double fC1;
	double fC2;
	double fC3;
	double fC4;
	double fC5;
	double fC6;
	double fC7;
	double fC8;
	double fC9;
	double fC10;


};

    
} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_N_FUNCTIONS_H_ */
