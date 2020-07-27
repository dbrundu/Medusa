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
 *
 *  WARNING: This is a temporary
 *  and not optimized version, based on the augalves new design.
 *  Checkout augalves branch for the final refactory.
 */

#ifndef PHIS_TIME_FUNCTIONS_H_
#define PHIS_TIME_FUNCTIONS_H_





namespace medusa {

namespace detail {




struct TimeCoeficients
{
    double fC[10];
};




/*
struct __hydra_align__(16) TimeFactors {

    __hydra_host__ __hydra_device__
    TimeFactors(double p0, double p1, double p2, double p3,
               double p4, double p5, double p6, double p7,
               double p8, double p9, double p10, double p11,
               double H3,  double H4, double H5,
               double time, double CP) :
        p_0(p0), p_par(p1), p_perp(p2), p_S(p3),
        l_0(p4), l_par(p5), l_perp(p6), l_S(p7),
        d_0(p8), d_par(p9), d_perp(p10), d_S(p11),
        fH3(H3), fH4(H4), fH5(H5), ftime(time), fCP(CP)
    {
        const static double f  = 0.238732414637843003653;
        const static double dG = 0.65789;
        
        fA[0] = 0.5*(1.0 + l_0 * l_0);
        fA[1] = 0.5*(1.0 + l_par*l_par);
        fA[2] = 0.5*(1.0 + l_perp*l_perp);
        
        fA[3] = ::sin(d_perp - d_par - p_perp + p_par);
        fA[3] *= l_par * l_perp;
        fA[3] = ::sin(d_perp - d_par) - fA[3];
        fA[3] *= 0.5;
        
        fA[4] = ::cos(d_0 - d_par - p_0 + p_par);
        fA[4] *= l_0 * l_par;
        fA[4] = ::cos(d_0 - d_par) + fA[4];
        fA[4] *= 0.5;
        
        fA[5] = ::sin(d_0 - d_perp - p_0 + p_perp);
        fA[5] *= l_0 * l_perp;
        fA[5] = ::sin(d_0 - d_perp) - fA[5];
        fA[5] *= -0.5;
        
        fA[6] = 0.5*(1.0 + l_S*l_S);
        
        fA[7] = ::cos(d_S - d_par - p_S + p_par);
        fA[7] *= l_S * l_par;
        fA[7] = ::cos(d_S - d_par) - fA[7];
        fA[7] *= 0.5;
        
        fA[8] = ::sin(d_S - d_perp - p_S + p_perp);
        fA[8] *= l_S * l_perp;
        fA[8] = ::sin(d_S - d_perp) + fA[8] ;
        fA[8] *= -0.5;
        
        fA[9] = ::cos(d_S - d_0 - p_S + p_0);
        fA[9] *= l_S * l_0;
        fA[9] = ::cos(d_S - d_0) - fA[9];
        fA[9] *= 0.5;
        
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        
        
        f1B = -l_0 * ::cos(p_0);
        f2B = -l_par * ::cos(p_par);
        f3B = l_perp * ::cos(p_perp);
        
        f4B = l_perp * ::sin( d_perp - d_par - p_perp);
        f4B = l_par  * ::sin( d_par - d_perp - p_par) + f4B;
        f4B *= 0.5;
        
        f5B = l_0 * ::cos(d_0 - d_par - p_0);
        f5B = l_par * ::cos(d_par - d_0 - p_par) + f5B;
        f5B *= -0.5;
        
        f6B =  l_0 * ::sin(d_0 - d_perp - p_0);
        f6B = l_perp * ::sin(d_perp - d_0 - p_perp) + f6B;
        f6B *= 0.5;
        
        f7B = l_S * ::cos(p_S);
        
        f8B = l_S * ::cos(d_S - d_par - p_S);
        f8B = f8B - l_par * ::cos(d_par - d_S - p_par);
        f8B *= 0.5;
        
        f9B = l_S * ::sin( d_S - d_perp - p_S);
        f9B = f9B - l_perp * ::sin( d_perp - d_S - p_perp);
        f9B *= -0.5;
        
        f10B = l_S * ::cos(d_S - d_0 - p_S);
        f10B = f10B - l_0 * ::cos( d_0 - d_S - p_0);
        f10B *= 0.5;
        

        
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
       
        f1C = 0.5*(1.0 - l_0 * l_0);
        f2C = 0.5*(1.0 - l_par*l_par);
        f3C = 0.5*(1.0 - l_perp*l_perp);
        
        f4C = ::sin(d_perp - d_par - p_perp + p_par);
        f4C *= l_par * l_perp;
        f4C = ::sin(d_perp - d_par) + f4C;
        f4C *= 0.5;
        
        f5C = ::cos(d_0 - d_par - p_0 + p_par);
        f5C *= l_0 * l_par;
        f5C = ::cos(d_0 - d_par) - f5C;
        f5C *= 0.5;
        
        f6C = ::sin(d_0 - d_perp - p_0 + p_perp);
        f6C *= l_0 * l_perp;
        f6C = ::sin(d_0 - d_perp) + f6C;
        f6C *= -0.5;
        
        f7C = 0.5*(1.0 - l_S*l_S);
        
        f8C = ::cos(d_S - d_par - p_S + p_par);
        f8C *= l_S * l_par;
        f8C = ::cos(d_S - d_par) + f8C;
        f8C *= 0.5;
        
        f9C = ::sin(d_S - d_perp - p_S + p_perp);
        f9C *= l_S * l_perp;
        f9C = ::sin(d_S - d_perp) - f9C ;
        f9C *= -0.5;
        
        f10C = ::cos(d_S - d_0 - p_S + p_0);
        f10C *= l_S * l_0;
        f10C = ::cos(d_S - d_0) + f10C;
        f10C *= 0.5;
        
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        
        f1D = l_0 * ::sin(p_0);
        f2D = l_par * ::sin(p_par);
        f3D = -l_perp * ::sin(p_perp);
        
        f4D = l_perp * ::cos( d_perp - d_par - p_perp);
        f4D = l_par  * ::cos( d_par - d_perp - p_par) + f4D;
        f4D *= -0.5;
        
        f5D = l_0 * ::sin(d_0 - d_par - p_0);
        f5D = l_par * ::sin(d_par - d_0 - p_par) + f5D;
        f5D *= -0.5;
        
        f6D =  l_0 * ::cos(d_0 - d_perp - p_0);
        f6D = l_perp * ::cos(d_perp - d_0 - p_perp) + f6D;
        f6D *= -0.5;
        
        f7D = -l_S * ::sin(p_S);
        
        f8D = l_S * ::sin(d_S - d_par - p_S);
        f8D = f8D - l_par * ::sin(d_par - d_S - p_par);
        f8D *= 0.5;
        
        f9D = -l_S * ::cos( d_S - d_perp - p_S);
        f9D = f9D + l_perp * ::cos( d_perp - d_S - p_perp);
        f9D *= -0.5;
        
        f10D = l_S * ::sin(d_S - d_0 - p_S);
        f10D = f10D - l_0 * ::sin( d_0 - d_S - p_0);
        f10D *= 0.5;
        
        ////////////////////////////////////////////////////
        ////////////////////////////////////////////////////
        
        _a = ::cosh(0.5 * ftime * fH4);
        _b = ::sinh(0.5 * ftime * fH4);
        _c = ::cos(ftime * fH5)*fCP;
        _d = ::sin(ftime * fH5)*fCP;
        
        cons = f * ::exp( -(fH3 + dG) * ftime);
        
        fB1 = f1A*_a + f1B*_b + f1C*_c + f1D*_d;
        fB1 *= cons;
        
        fB2 = f2A*_a + f2B*_b + f2C*_c + f2D*_d;
        fB2 *= cons;
        
        fB3 = f3A*_a + f3B*_b + f3C*_c + f3D*_d;
        fB3 *= cons;
        
        fB4 = f4A*_a + f4B*_b + f4C*_c + f4D*_d;
        fB4 *= cons;
        
        fB5 = f5A*_a + f5B*_b + f5C*_c + f5D*_d;
        fB5 *= cons;
        
        fB6 = f6A*_a + f6B*_b + f6C*_c + f6D*_d;
        fB6 *= cons;
        
        fB7 = f7A*_a + f7B*_b + f7C*_c + f7D*_d;
        fB7 *= cons;
        
        fB8 = f8A*_a + f8B*_b + f8C*_c + f8D*_d;
        fB8 *= cons;
        
        fB9 = f9A*_a + f9B*_b + f9C*_c + f9D*_d;
        fB9 *= cons;
        
        fB10 = f10A*_a + f10B*_b + f10C*_c + f10D*_d;
        fB10 *= cons;
        
        
        
        }
        

        double fB1;
        double fB2;
        double fB3;
        double fB4;
        double fB5;
        double fB6;
        double fB7;
        double fB8;
        double fB9;
        double fB10;
        
private:

     double p_0;
     double p_par;
     double p_perp;
     double p_S;
     double l_0;
     double l_par;
     double l_perp;
     double l_S;
     double d_0;
     double d_par;
     double d_perp;
     double d_S; 
     
     double fH3, fH4, fH5;
     double ftime;
     double fCP;
     
     double fA[10];
     double fB[10];
     double fC[10];
     double fD[10];
     
     double f1A, f1B, f1C, f1D;
     double f2A, f2B, f2C, f2D;
     double f3A, f3B, f3C, f3D;
     double f4A, f4B, f4C, f4D;
     double f5A, f5B, f5C, f5D;
     double f6A, f6B, f6C, f6D;
     double f7A, f7B, f7C, f7D;
     double f8A, f8B, f8C, f8D;
     double f9A, f9B, f9C, f9D;
     double f10A, f10B, f10C, f10D;
     
     double _a, _b, _c, _d, cons;
    
};
*/
    
    
    


    
} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_TIME_FUNCTIONS_H_ */
        

