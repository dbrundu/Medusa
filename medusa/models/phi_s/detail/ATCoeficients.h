
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
 * ATCoeficients.h
 *
 *  Created on: 18/07/2020
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef ATCoeficients_H_
#define ATCoeficients_H_

#include<map>
#include<array>

namespace medusa {

namespace detail {

struct AngularTimeCoeficients
{

	double fC[10];
};

/*
void fill_ATCoeficients( const double (&data)[20], AngularTimeCoeficients& A, AngularTimeCoeficients& B,
		AngularTimeCoeficients& C, AngularTimeCoeficients& D)
{

	//
	// data map:
	//
	// data[0]: sin(delta_perp) , data[1]: cos(delta_perp)
	// data[2]: sin(delta_par)  , data[3]: cos(delta_par)
	// data[4]: sin(delta_0)    , data[5]: cos(delta_0)
	// data[6]: sin(delta_S)    , data[7]: cos(delta_S)
	//
	// data[8]: sin(phi_perp)   ,  data[9]: cos(phi_perp)
	// data[10]: sin(phi_par)   , data[11]: cos(phi_par)
	// data[12]: sin(phi_0)     , data[13]: cos(phi_0)
	// data[14]: sin(phi_S)     , data[15]: cos(phi_S)
	//
	// data[16]: lambda_0       , data[17]: lambda_par
	// data[18]: lambda_perp    , data[19]: lambda_S

	double sdeltax = data[0]; double cdeltax = data[1];
	double sdeltay = data[2]; double cdeltay = data[3];
	double sdeltaz = data[4]; double cdeltaz = data[5];
	double sdeltaw = data[6]; double cdeltaw = data[7];

	double sphix = data[8] ;  double cphix = data[9];
	double sphiy = data[10];  double cphiy = data[11];
	double sphiz = data[12];  double cphiz = data[13];
	double sphiw = data[14];  double cphiw = data[14];

	lambda_0    = data[16];
	lambda_par  = data[17];
	lambda_perp = data[18];
	lambda_S    = data[19];

	//using the infamous identities
	//
	//sin(α – β) = sin(α) cos(β) – cos(α) sin(β)
	//sin(α + β) = sin(α) cos(β) + cos(α) sin(β)
	//
	//cos(α + β) = cos(α) cos(β) – sin(α) sin(β)
	//cos(α – β) = cos(α) cos(β) + sin(α) sin(β)
	//
	//rationale: we have over 50 calls to trigonometric
	//functions taking as arguments additive or substractive
	//of the 8 angles. The math is quite simple indeed:
	//with 16 calls to trigonometric, we can derive all other calls
	//as multiplications using the usual trigonometric
	//identities (see above).
	//Keep in mind that, even considering modern processors
	//have dedicated instructions for trigonometric functions
	//and this makes such calls ILP friendly into the pipeline,
	//a trigonometric call takes 15-20 cycles, while
	//a multiplication takes 1. So guys, make the math!

	//::sin(delta_perp - delta_par)
	double sdelta_xmy = sdeltax*cdeltay - cdeltax*sdeltay;
	//::cos(delta_perp - delta_par)
	double cdelta_xmy = cdeltax*cdeltay + sdeltax*sdeltay;

	//::cos(delta_0 - delta_par)
	double cdelta_zmy = cdeltaz*cdeltay + sdeltaz*sdeltay;
	//::sin(delta_0 - delta_par)
	double sdelta_zmy = cdeltaz*sdeltay - sdeltaz*cdeltay;

	//::sin(delta_0 - delta_perp)
	double sdelta_zmx = sdeltaz*cdeltax - cdeltaz*sdeltax;
	//::cos(delta_0 - delta_perp)
	double cdelta_zmx = cdeltaz*cdeltax + sdeltaz*sdeltax;

	//::sin(delta_S - delta_par)
	double sdelta_wmy = cdeltaw*sdeltay - sdeltaw*sdeltay;
	//::cos(delta_S - delta_par)
	double cdelta_wmy = cdeltaz*cdeltay + sdeltaz*sdeltay;

	//::sin(delta_S - delta_0)
	double sdelta_wmz = cdeltaw*sdeltaz - sdeltaw*cdeltaz;
	//::cos(delta_S - delta_0)
	double cdelta_wmz = cdeltaw*cdeltaz + sdeltaw*sdeltaz;

	//::sin(delta_perp + delta_par)
	double sdelta_xpy = sdeltax*cdeltay + cdeltax*sdeltay;
	//::cos(delta_0 + delta_par)
	double cdelta_zpy = cdeltaz*cdeltay - sdeltaz*sdeltay;
	//::sin(delta_0 + delta_perp)
	double sdelta_zpx = sdeltaz*cdeltax + cdeltaz*sdeltax;
	//::cos(delta_S + delta_par)
	double cdelta_wpy = cdeltaz*cdeltay - sdeltaz*sdeltay;
	//::sin(delta_S + delta_perp)
	double sdelta_wpx = sdeltaw*cdeltax + cdeltaw*sdeltax;
	//::cos(delta_S + delta_0)
	double cdelta_wpz = cdeltaw*cdeltaz - sdeltaw*sdeltaz;

	//::sin(phi_perp + phi_par)
	double sphi_xpy = sphix*cphiy + cphix*sphiy;
	//::cos(phi_perp + phi_par)
	double cphi_xpy = cphix*cphiy - sphix*sphiy;
	//::sin(phi_perp - phi_par)
	double sphi_xmy = sphix*cphiy - cphix*sphiy;
	//::cos(phi_perp - phi_par)
	double cphi_xmy = cphix*cphiy + sphix*sphiy;

	//::sin(phi_0 + phi_par)
	double sphi_zpy = sphiz*cphiy + cphiz*sphiy;
	//::cos(phi_0 + phi_par)
	double cphi_zpy = cphiz*cphiy - sphiz*sphiy;
	//::sin(phi_0 - phi_par)
	double sphi_zmy = sphiz*cphiy - cphiz*sphiy;
	//::cos(phi_0 - phi_par)
	double cphi_zmy = cphiz*cphiy + sphiz*sphiy;

	//::sin(phi_S + phi_par)
	double sphi_wpy = sphiw*cphiy + cphiw*sphiy;
	//::cos(phi_S + phi_par)
	double cphi_wpy = cphiw*cphiy - sphiw*sphiy;
	//::sin(phi_S - phi_par)
	double sphi_wmy = sphiw*cphiy - cphiw*sphiy;
	//::cos(phi_S - phi_par)
	double cphi_wmy = cphiw*cphiy + sphiw*sphiy;

	//::sin(phi_S + phi_0)
	double sphi_wpz = sphiw*cphiz + cphiw*sphiz;
	//::cos(phi_S + phi_0)
	double cphi_wpz = cphiw*cphiz - sphiw*sphiz;
	//::sin(phi_S - phi_0)
	double sphi_wmz = sphiw*cphiz - cphiw*sphiz;
	//::cos(phi_S - phi_0)
	double cphi_wmz = cphiw*cphiz + sphiw*sphiz;

	//::sin[(delta_perp - delta_par) - (phi_perp + phi_par) ]
	double sdelta_xmy_m_phi_xpy = sdelta_xmy*cphi_xpy - cdelta_xmy*sphi_xpy;

	//::cos(delta_0 - delta_par - phi_0 + phi_par)
	double cdelta_zmy_m_phi_zpy = cdelta_zmy*cphi_zpy + sdelta_xmy*sphi_zpy;

	//::sin(delta_0 - delta_perp - phi_0 + phi_perp)
	double sdelta_zmy_m_phi_zpy = cdelta_zmy*sphi_zpy - sdelta_xmy*cphi_zpy;

	//::cos(delta_S - delta_par - phi_S + phi_par)
	double cdelta_wmy_m_phi_wpy = cdelta_wmy*cphi_wpy + sdelta_wmy*sphi_wpy;

	//::sin(delta_S - delta_perp - phi_S + phi_perp)
	double sdelta_wmx_m_phi_wpx = cdelta_wmx*cphi_wpy - sdelta_wmy*sphi_wpy;

	//::cos(delta_S - delta_0 - phi_S + phi_0)
	double cdelta_wmz_m_phi_wpz = cdelta_wmz*cphi_wpz + sdelta_wmz*sphi_wpy;

	//::sin(delta_perp - delta_par - phi_perp)
	double sdelda_xmy_m_phix = sdelta_xmy*cphix - cdelta_xmy*sphix;
	//::cos(delta_perp - delta_par - phi_perp)
	double cdelda_xmy_m_phix = cdelta_xmy*cphix + sdelta_xmy*sphix;

	//::sin(delta_par - delta_perp - phi_par)= -::sin( phi_par + (delta_perp - delta_par) )
	double sdelta_ymx_m_phiy = -( sphiy*cdelta_xmy + cphiy*sdelta_xmy);
	//::cos(delta_par - delta_perp - phi_par)= ::cos( phi_par + (delta_perp - delta_par) )
	double cdelta_ymx_m_phiy =  cphiy*cdelta_xmy - sphiy*sdelta_xmy;

	//::cos(delta_0 - delta_par - phi_0)
	double cdelta_zmy_m_phiz = cdelta_zmy*cphiz + sdelta_zmy*sphiz;
	//::sin(delta_0 - delta_par - phi_0)
	double sdelta_zmy_m_phiz = sdelta_zmy*cphiz - cdelta_zmy*sphiz;

	//::cos(delta_par - delta_0 - phi_par)=::cos(delta_0 - delta_par + phi_par)
	double cdelta_ymz_m_phiy = cdelta_zmy*cphiy - sdelta_zmy*sphiy;
	//::sin(delta_par - delta_0 - phi_par)= -::sin(delta_0 - delta_par + phi_par)
	double sdelta_ymz_m_phiy = -(sdelta_zmy*cphiy + cdelta_zmy*sphiy);

	//::sin(delta_0 - delta_perp - phi_0)
	double sdelta_zmx_m_phiz = sdelta_zmx*cphiz - cdelta_zmx*sphiz;
	//::cos(delta_0 - delta_perp - phi_0)
	double cdelta_zmx_m_phiz = cdelta_zmx*cphiz + sdelta_zmx*sphiz;

	//::sin(delta_perp - delta_0 - phi_perp)=-::sin(delta_0 - delta_perp + phi_perp)
	double sdelta_xmz_m_phix = -(sdelta_zmx*cphix + cdelta_zmx*sphix);
	//::cos(delta_perp - delta_0 - phi_perp)=::cos(delta_0 - delta_perp + phi_perp)
	double cdelta_xmz_m_phix = cdelta_zmx*cphix - sdelta_zmx*sphix;

	//::cos(delta_S - delta_par - phi_S)
	double cdelta_wmy_m_phiw = cdelta_wmy*cphiw - sdelta_wmy*sphiw;
	//::sin(delta_S - delta_par - phi_S)
	double sdelta_wmy_m_phiw = sdelta_wmy*cphiw + cdelta_wmy*sphiw;

	//::cos(delta_par - delta_S - phi_par) = ::cos(delta_S - delta_par + phi_par )
	double cdelta_ymw_m_phiy = cdelta_wmy*cphiy -  sdelta_wmy*sphiy;
	//::sin(delta_par - delta_S - phi_par) = -::sin(delta_S - delta_par + phi_par )
	double sdelta_ymw_m_phiy = -(cdelta_wmy*sphiy +  sdelta_wmy*cphiy);

	//::sin(delta_S - delta_perp - phi_S)
	double sdelta_wmx_m_phiw = 	sdelta_wmx*cphiw - cdelta_wmx*sphiw;
	//::cos(delta_S - delta_perp - phi_S)
	double cdelta_wmx_m_phiw = 	cdelta_wmx*cphiw + sdelta_wmx*sphiw;

	//::sin(delta_perp - delta_S - phi_perp) = -::sin(delta_S - delta_perp + phi_perp)
	double sdelta_xmw_m_phix = -(sdelta_wmx*cphix + cdelta_wmx*sphix);
	//::cos(delta_perp - delta_S - phi_perp) = ::cos(delta_S - delta_perp + phi_perp)
	double cdelta_xmw_m_phix = -(cdelta_wmx*cphix + sdelta_wmx*sphix);

	//::cos(delta_S - delta_0 - phi_S)
	double cdelta_wmz_m_phiw = cdelta_wmz*cphiw + sdelta_wmz*sphiw;
	//::sin(delta_S - delta_0 - phi_S)
	double sdelta_wmz_m_phiw = cdelta_wmz*sphiw - sdelta_wmz*cphiw;

	//::cos(delta_0 - delta_S - phi_0) = ::cos(delta_S - delta_0 + phi_0)
	double cdelta_zmw_m_phiz =  cdelta_wmz*cphiz - sdelta_wmz*sphiw;
	//::sin(delta_0 - delta_S - phi_0) = -::sin(delta_S - delta_0 + phi_0)
	double sdelta_zmw_m_phiz =  -(cdelta_wmz*cphiz + cdelta_wmz*sphiw);

	//--------------------------------------------------------------------
	// A
	//--------------------------------------------------------------------

	A.fC[0] = 0.5*(1.0 + lambda_0*lambda_0);

	A.fC[1] = 0.5*(1.0 + lambda_par*lambda_par);

	A.fC[2] = 0.5*(1.0 + lambda_perp*lambda_perp);

	//0.5*(::sin(delta_perp - delta_par) - lambda_par*lambda_perp*::sin(delta_perp - delta_par - phi_perp + phi_par) )
	A.fC[3] = 0.5*(sdelta_xmy - lambda_par*lambda_perp*sdelta_xmy_m_phi_xpy);

	//0.5*(::cos(delta_0 - delta_par) + lambda_0*lambda_par*::cos(delta_0 - delta_par - phi_0 + phi_par) )
	A.fC[4] = 0.5*(cdelta_zmy + lambda_0*lambda_par*cdelta_zmy_m_phi_zpy );

	//-0.5*(::sin(delta_0 - delta_perp) - lambda_0*lambda_perp*::sin(delta_0 - delta_perp - phi_0 + phi_perp))
	A.fC[5] = -0.5*(sdelta_zmx - lambda_0*lambda_perp*sdelta_zmy_m_phi_zpy );

	//0.5*(1 + ::pow(lambda_S,2) )
	A.fC[6] = 0.5*(1 + lambda_S*lambda_S );

	//0.5*(::cos(delta_S - delta_par) - lambda_S*lambda_par*::cos(delta_S - delta_par - phi_S + phi_par))
	A.fC[7] = 0.5*(cdelta_wmy - lambda_S*lambda_par*cdelta_wmy_m_phi_wpy );

	//-0.5*(::sin(delta_S - delta_perp) + lambda_S*lambda_perp*::sin(delta_S - delta_perp - phi_S + phi_perp))
	A.fC[8] = -0.5*(sdelta_wmx + lambda_S*lambda_perp*sdelta_wmx_m_phi_wpx);

	//0.5*(::cos(delta_S - delta_0) - lambda_S*lambda_0*::cos(delta_S - delta_0 - phi_S + phi_0))
	A.fC[9] = 0.5*(cdelta_wmz - lambda_S*lambda_0*cdelta_wmz_m_phi_wpz);

	//--------------------------------------------------------------------
	// B
	//--------------------------------------------------------------------

	// -lambda_0*::cos(phi_0)
	B.fC[0] = -lambda_0*cphiz;

	// -lambda_par*::cos(phi_par)
	B.fC[1] = -lambda_par*cphiy;

	//lambda_perp*::cos(phi_perp)
	B.fC[2] = lambda_perp*cphix ;

	//0.5*(lambda_perp*::sin(delta_perp - delta_par - phi_perp) + lambda_par*::sin(delta_par - delta_perp - phi_par))
	B.fC[3] = 0.5*(lambda_perp*sdelda_xmy_m_phix + lambda_par*sdelta_ymx_m_phiy);

	//-0.5*(lambda_0*::cos(delta_0 - delta_par - phi_0) + lambda_par*::cos(delta_par - delta_0 - phi_par))
	B.fC[4] = -0.5*(lambda_0*cdelta_zmy_m_phiz + lambda_par*cdelta_ymz_m_phiy);

	//0.5*(lambda_0*::sin(delta_0 - delta_perp - phi_0) + lambda_perp*::sin(delta_perp - delta_0 - phi_perp))
	B.fC[5] = 0.5*(lambda_0*sdelta_zmx_m_phiz + lambda_perp*sdelta_xmz_m_phix);

	// lambda_S*::cos(phi_S)
	B.fC[6] =  lambda_S*cphiw;

	//0.5*(lambda_S*::cos(delta_S - delta_par - phi_S) - lambda_par*::cos(delta_par - delta_S - phi_par))
	B.fC[7] = 0.5*(lambda_S*cdelta_wmy_m_phiw - lambda_par*cdelta_ymw_m_phiy);

	//-0.5*(lambda_S*::sin(delta_S - delta_perp - phi_S) - lambda_perp*::sin(delta_perp - delta_S - phi_perp))
	B.fC[8] = -0.5*(lambda_S*sdelta_wmx_m_phiw - lambda_perp*sdelta_xmw_m_phix );

	// 0.5*(lambda_S*::cos(delta_S - delta_0 - phi_S) - lambda_0*::cos(delta_0 - delta_S - phi_0))
	B.fC[9] =  0.5*(lambda_S*cdelta_wmz_m_phiw - lambda_0*cdelta_zmw_m_phiz);

	//--------------------------------------------------------------------
	// C
	//--------------------------------------------------------------------

	//0.5*(1 - pow(lambda_0,2))
	C.fC[0] = 0.5*(1 - lambda_0*lambda_0);

	//0.5*(1 - pow(lambda_par,2))
	C.fC[1] = 0.5*(1 - lambda_par*lambda_0);

	//0.5*(1 - pow(lambda_perp,2))
	C.fC[2] = 0.5*(1 - lambda_perp*lambda_perp);

	//0.5*(::sin(delta_perp - delta_par) + lambda_par*lambda_perp*::sin(delta_perp - delta_par - phi_perp + phi_par) )
	C.fC[3] = 0.5*( sdelta_xmy + lambda_par*lambda_perp*sdelta_xmy_m_phi_xpy );

	//0.5*(::cos(delta_0 - delta_par) - lambda_0*lambda_par*::cos(delta_0 - delta_par - phi_0 + phi_par) )
	C.fC[4] = 0.5*(cdelta_zmy - lambda_0*lambda_par*cdelta_zmy_m_phi_zpy  );

	//-0.5*(::sin(delta_0 - delta_perp) + lambda_0*lambda_perp*::sin(delta_0 - delta_perp - phi_0 + phi_perp))
	C.fC[5] = -0.5*(sdelta_zmx + lambda_0*lambda_perp*sdelta_zmy_m_phi_zpy);

	//0.5*(1 - ::pow(lambda_S,2) )
	C.fC[6] = 0.5*(1 - lambda_S*lambda_S );

	//0.5*(::cos(delta_S - delta_par) + lambda_S*lambda_par*::cos(delta_S - delta_par - phi_S + phi_par))
	C.fC[7] = 0.5*(cdelta_wmy + lambda_S*lambda_par*cdelta_wmy_m_phi_wpy);

	//-0.5*(::sin(delta_S - delta_perp) - lambda_S*lambda_perp*::sin(delta_S - delta_perp - phi_S + phi_perp))
	C.fC[8] = -0.5*(sdelta_wmx - lambda_S*lambda_perp*sdelta_wmx_m_phi_wpx);

	//0.5*(::cos(delta_S - delta_0) + lambda_S*lambda_0*::cos(delta_S - delta_0 - phi_S + phi_0))
	C.fC[9] = 0.5*(cdelta_wmz + lambda_S*lambda_0*cdelta_wmz);

	//--------------------------------------------------------------------
	// D
	//--------------------------------------------------------------------

	//lambda_0*::sin(phi_0)
	D.fC[0] = lambda_0*sphiz;

	//lambda_par*::sin(phi_par)
	D.fC[1] = lambda_y*sphiy;

	//-lambda_perp*::sin(phi_perp)
	D.fC[2] = -lambda_x*sphiy;

	// -0.5*(lambda_perp*::cos(delta_perp - delta_par - phi_perp) + lambda_par*::cos(delta_par - delta_perp - phi_par))
	D.fC[3] = -0.5*(lambda_perp*cdelda_xmy_m_phix + lambda_par*cdelta_ymx_m_phiy);

	// -0.5*(lambda_0*::sin(delta_0 - delta_par - phi_0) + lambda_par*::sin(delta_par - delta_0 - phi_par))
	D.fC[4] = -0.5*(lambda_0*sdelta_zmy_m_phiz + lambda_par*::sin(delta_par - delta_0 - phi_par));

	//-0.5*(lambda_0*::cos(delta_0 - delta_perp - phi_0) + lambda_perp*::cos(delta_perp - delta_0 - phi_perp))
	D.fC[5] = -0.5*(lambda_0*cdelta_zmx_m_phiz + lambda_perp*cdelta_xmz_m_phix);

	//-lambda_S*::sin(phi_S)
	D.fC[6] = -lambda_S*sphiw;

	//0.5*(lambda_S*::sin(delta_S - delta_par - phi_S) - lambda_par*::sin(delta_par - delta_S - phi_par))
	D.fC[7] =  0.5*(lambda_S*sdelta_wmy_m_phiw - lambda_par*sdelta_ymw_m_phiy);

	//-0.5*(-lambda_S*::cos(delta_S - delta_perp - phi_S) + lambda_perp*::cos(delta_perp - delta_S - phi_perp))
	D.fC[8] = -0.5*(-lambda_S*cdelta_wmx_m_phiw + lambda_perp*cdelta_xmw_m_phix );

	//0.5*(lambda_S*::sin(delta_S - delta_0 - phi_S) - lambda_0*::sin(delta_0 - delta_S - phi_0))
	D.fC[9]  = 0.5*(lambda_S*sdelta_wmz_m_phiw - lambda_0*sdelta_zmw_m_phiz);

}
*/

} // namespace detail


}  // namespace medusa

#endif /* ATCoeficients_H_ */
