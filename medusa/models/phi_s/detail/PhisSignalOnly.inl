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
 * PhisSignalOnly.inl
 *
 *  Created on: 24/07/2020
 *      Author: augalves
 *
 *  WARNING: This is a temporary
 *  and not optimized version, based on the augalves new design.
 *  Checkout augalves branch for the final refactory.
 */

#ifndef PHISSIGNALONLY_INL_
#define PHISSIGNALONLY_INL_


namespace medusa {


template<bool B0bar,
typename ArgTypeTime,
typename ArgTypeThetah,
typename ArgTypeThetal,
typename ArgTypePhi,
typename Signature >
void PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi, Signature>::Update_NFactors()
{
    	/*
    	0: A_0,
        1: A_perp,
    	2: A_S,,
    	*/
    	double A_par = ::sqrt(1 - _par[0]*_par[0] -_par[1]*_par[1]);

    	fN.fC[0] = _par[0]*_par[0];//A_0*A_0 ;
    	fN.fC[1] =   A_par*A_par;
    	fN.fC[2] = _par[1]*_par[1]; //A_perp*A_perp;
    	fN.fC[3] = _par[1]*A_par;   //A_perp*A_par;
    	fN.fC[4] = _par[0]*A_par;   //A_0*A_par;
    	fN.fC[5] = _par[0]*_par[1]; //A_0*A_perp;
    	fN.fC[6] = _par[2]* _par[2];//A_S*A_S;
    	fN.fC[7] = _par[2]*A_par;   //A_S*A_par;
    	fN.fC[8] = _par[2]*_par[1]; //A_S*A_perp;
    	fN.fC[9] = _par[2]*_par[0]; //A_S*A_0;

}



template<bool B0bar,
typename ArgTypeTime,
typename ArgTypeThetah,
typename ArgTypeThetal,
typename ArgTypePhi,
typename Signature >
void PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi, Signature>::Update_TCoefficients()
{

	/*
	{
	 0: A_0,
	 1: A_perp,
	 2: A_S,
	 3: DeltaGamma_sd,
	 4: DeltaGamma,
	 5: DeltaM ,
	 6: phi_0,
	 7: phi_par,
	 8: phi_perp,
	 9: phi_S,
	 10: lambda_0,
	 11: lambda_par,
	 12: lambda_perp,
	 13: lambda_S,
	 14: delta_0,
	 15: delta_par,
	 16: delta_perp,
	 17: delta_S
	}
	*/

        double p_0    = _par[6];
        double p_par  = _par[7];
        double p_perp = _par[8];
        double p_S    = _par[9];
        double l_0    = _par[10];
        double l_par  = _par[11];
        double l_perp = _par[12];
        double l_S    = _par[13];
        double d_0    = _par[14];
        double d_par  = _par[15];
        double d_perp = _par[16];
        double d_S    = _par[17];

        fA.fC[0] = 0.5*(1.0 + l_0 * l_0);
        fA.fC[1] = 0.5*(1.0 + l_par*l_par);
        fA.fC[2] = 0.5*(1.0 + l_perp*l_perp);
        
        fA.fC[3] = ::sin(d_perp - d_par - p_perp + p_par);
        fA.fC[3] *= l_par * l_perp;
        fA.fC[3] = ::sin(d_perp - d_par) - fA.fC[3];
        fA.fC[3] *= 0.5;
        
        fA.fC[4] = ::cos(d_0 - d_par - p_0 + p_par);
        fA.fC[4] *= l_0 * l_par;
        fA.fC[4] = ::cos(d_0 - d_par) + fA.fC[4];
        fA.fC[4] *= 0.5;
        
        fA.fC[5] = ::sin(d_0 - d_perp - p_0 + p_perp);
        fA.fC[5] *= l_0 * l_perp;
        fA.fC[5] = ::sin(d_0 - d_perp) - fA.fC[5];
        fA.fC[5] *= -0.5;
        
        fA.fC[6] = 0.5*(1.0 + l_S*l_S);
        
        fA.fC[7] = ::cos(d_S - d_par - p_S + p_par);
        fA.fC[7] *= l_S * l_par;
        fA.fC[7] = ::cos(d_S - d_par) - fA.fC[7];
        fA.fC[7] *= 0.5;
        
        fA.fC[8] = ::sin(d_S - d_perp - p_S + p_perp);
        fA.fC[8] *= l_S * l_perp;
        fA.fC[8] = ::sin(d_S - d_perp) + fA.fC[8] ;
        fA.fC[8] *= -0.5;
        
        fA.fC[9] = ::cos(d_S - d_0 - p_S + p_0);
        fA.fC[9] *= l_S * l_0;
        fA.fC[9] = ::cos(d_S - d_0) - fA.fC[9];
        fA.fC[9] *= 0.5;
        
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        
        
        fB.fC[0] = -l_0 * ::cos(p_0);
        fB.fC[1] = -l_par * ::cos(p_par);
        fB.fC[2] = l_perp * ::cos(p_perp);
        
        fB.fC[3] = l_perp * ::sin( d_perp - d_par - p_perp);
        fB.fC[3] = l_par  * ::sin( d_par - d_perp - p_par) + fB.fC[3];
        fB.fC[3] *= 0.5;
        
        fB.fC[4] = l_0 * ::cos(d_0 - d_par - p_0);
        fB.fC[4] = l_par * ::cos(d_par - d_0 - p_par) + fB.fC[4];
        fB.fC[4] *= -0.5;
        
        fB.fC[5] =  l_0 * ::sin(d_0 - d_perp - p_0);
        fB.fC[5] = l_perp * ::sin(d_perp - d_0 - p_perp) + fB.fC[5];
        fB.fC[5] *= 0.5;
        
        fB.fC[6] = l_S * ::cos(p_S);
        
        fB.fC[7] = l_S * ::cos(d_S - d_par - p_S);
        fB.fC[7] = fB.fC[7] - l_par * ::cos(d_par - d_S - p_par);
        fB.fC[7] *= 0.5;
        
        fB.fC[8] = l_S * ::sin( d_S - d_perp - p_S);
        fB.fC[8] = fB.fC[8] - l_perp * ::sin( d_perp - d_S - p_perp);
        fB.fC[8] *= -0.5;
        
        fB.fC[9] = l_S * ::cos(d_S - d_0 - p_S);
        fB.fC[9] = fB.fC[9] - l_0 * ::cos( d_0 - d_S - p_0);
        fB.fC[9] *= 0.5;
        

        
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
       
        fC.fC[0] = 0.5*(1.0 - l_0 * l_0);
        fC.fC[1] = 0.5*(1.0 - l_par*l_par);
        fC.fC[2] = 0.5*(1.0 - l_perp*l_perp);
        
        fC.fC[3] = ::sin(d_perp - d_par - p_perp + p_par);
        fC.fC[3] *= l_par * l_perp;
        fC.fC[3] = ::sin(d_perp - d_par) + fC.fC[3];
        fC.fC[3] *= 0.5;
        
        fC.fC[4] = ::cos(d_0 - d_par - p_0 + p_par);
        fC.fC[4] *= l_0 * l_par;
        fC.fC[4] = ::cos(d_0 - d_par) - fC.fC[4];
        fC.fC[4] *= 0.5;
        
        fC.fC[5] = ::sin(d_0 - d_perp - p_0 + p_perp);
        fC.fC[5] *= l_0 * l_perp;
        fC.fC[5] = ::sin(d_0 - d_perp) + fC.fC[5];
        fC.fC[5] *= -0.5;
        
        fC.fC[6] = 0.5*(1.0 - l_S*l_S);
        
        fC.fC[7] = ::cos(d_S - d_par - p_S + p_par);
        fC.fC[7] *= l_S * l_par;
        fC.fC[7] = ::cos(d_S - d_par) + fC.fC[7];
        fC.fC[7] *= 0.5;
        
        fC.fC[8] = ::sin(d_S - d_perp - p_S + p_perp);
        fC.fC[8] *= l_S * l_perp;
        fC.fC[8] = ::sin(d_S - d_perp) - fC.fC[8] ;
        fC.fC[8] *= -0.5;
        
        fC.fC[9] = ::cos(d_S - d_0 - p_S + p_0);
        fC.fC[9] *= l_S * l_0;
        fC.fC[9] = ::cos(d_S - d_0) + fC.fC[9];
        fC.fC[9] *= 0.5;
        
        ///////////////////////////////////////////////////
        /////////////////////////////////////////////////// 
        
        fD.fC[0] = l_0 * ::sin(p_0);
        fD.fC[1] = l_par * ::sin(p_par);
        fD.fC[2] = -l_perp * ::sin(p_perp);
        
        fD.fC[3] = l_perp * ::cos( d_perp - d_par - p_perp);
        fD.fC[3] = l_par  * ::cos( d_par - d_perp - p_par) + fD.fC[3];
        fD.fC[3] *= -0.5;
        
        fD.fC[4] = l_0 * ::sin(d_0 - d_par - p_0);
        fD.fC[4] = l_par * ::sin(d_par - d_0 - p_par) + fD.fC[4];
        fD.fC[4] *= -0.5;
        
        fD.fC[5] =  l_0 * ::cos(d_0 - d_perp - p_0);
        fD.fC[5] = l_perp * ::cos(d_perp - d_0 - p_perp) + fD.fC[5];
        fD.fC[5] *= -0.5;
        
        fD.fC[6] = -l_S * ::sin(p_S);
        
        fD.fC[7] = l_S * ::sin(d_S - d_par - p_S);
        fD.fC[7] = fD.fC[7] - l_par * ::sin(d_par - d_S - p_par);
        fD.fC[7] *= 0.5;
        
        fD.fC[8] = -l_S * ::cos( d_S - d_perp - p_S);
        fD.fC[8] = fD.fC[8] + l_perp * ::cos( d_perp - d_S - p_perp);
        fD.fC[8] *= -0.5;
        
        fD.fC[9] = l_S * ::sin(d_S - d_0 - p_S);
        fD.fC[9] = fD.fC[9] - l_0 * ::sin( d_0 - d_S - p_0);
        fD.fC[9] *= 0.5;
        


}


}  // namespace medusa


#endif /* PHISSIGNALONLY_INL_ */
