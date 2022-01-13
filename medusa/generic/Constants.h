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
 *  Created on: 17/12/2021
 *
 *  Author: Alessandro Maria Ricci
 *
 *  Constants.h
 *----------------------------------------*/

#ifndef MEDUSA_CONSTANTS_H
#define MEDUSA_CONSTANTS_H

// Hydra
#include <hydra/Types.h>
#include <hydra/Complex.h>

// use std::numeric_limits, since 1./0. and 0./0. fail with some compilers (MS)
#define Inf std::numeric_limits<double>::infinity()
#define NaN std::numeric_limits<double>::quiet_NaN()


namespace medusa {
    namespace math_constants {

		const double M_1_3      = 0.33333333333333333333;  // 1/3
        const double M_1_Sqrt8  = 0.35355339059327376220;  // 1/Sqrt(8)
        const double M_1_Sqrt32 = 0.17677669529663688110;  // 1/sqrt(32)
        const double M_2_Sqrt3  = 1.15470053837925152902;  // 2/sqrt(3)
		const double M_2_Sqrt6  = 0.81649658092772603273;  // 2/sqrt(6)
        const double M_1_SqrtPi = 0.56418958354775628695;  // 1/sqrt(Pi)

    } // namespace math_constants
} // namespace medusa


// macro
#define M_Sqrt2 hydra::math_constants::sqrt2
#define M_1_Sqrt2 hydra::math_constants::inverse_sqrt2
#define M_1_3 medusa::math_constants::M_1_3
#define M_1_Sqrt8 medusa::math_constants::M_1_Sqrt8
#define M_1_Sqrt32 medusa::math_constants::M_1_Sqrt32
#define M_2_Sqrt3 medusa::math_constants::M_2_Sqrt3
#define M_2_Sqrt6 medusa::math_constants::M_2_Sqrt6
#define M_1_SqrtPi medusa::math_constants::M_1_SqrtPi

#endif // MEDUSA_CONSTANTS_H
