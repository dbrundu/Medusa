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


#define PHIS_N_FUNCTION(N, formula)\
template<>\
__hydra_dual__ \
inline double phis_N_functions<N>(double const& A_0, double const& A_perp,  double const& A_S, double const& A_par){\
\
	return formula;\
}\



namespace medusa {

namespace detail {


    template<size_t N>
    inline double phis_N_functions(double const& A_0, double const& A_perp,  double const& A_S, double const& A_par);
    
    
    PHIS_N_FUNCTION(0, ::pow(A_0,2) )
    
    PHIS_N_FUNCTION(1, ::pow(A_par,2))
    
    PHIS_N_FUNCTION(2, ::pow(A_perp,2))
    
    PHIS_N_FUNCTION(3, A_perp*A_par)
    
    PHIS_N_FUNCTION(4, A_0*A_par)
    
    PHIS_N_FUNCTION(5, A_0*A_perp)
    
    PHIS_N_FUNCTION(6, ::pow(A_S,2))
    
    PHIS_N_FUNCTION(7, A_S*A_par)
    
    PHIS_N_FUNCTION(8, A_S*A_perp)
    
    PHIS_N_FUNCTION(9, A_S*A_0)

    
    
} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_N_FUNCTIONS_H_ */
