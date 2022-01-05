/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu,
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto,
 * 						Alessandro Maria Ricci
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
/*---------------------------------------------------------------------------
 *  FullAnalyticPhis.inl
 *
 *  Created on: 29/10/2020
 *      Author: Antonio Augusto Alves Junior
 * 		Updated by Alessandro Maria Ricci in 04/11/2021
 * 
 * 	This function updates the values of the angular coefficients a_k, b_k,
 *  c_k, d_k by using the formulas in Table 3 in arXiv:1906.08356v4.
 *---------------------------------------------------------------------------*/

#ifndef FULL_ANALYTIC_PHIS_INL_
#define FULL_ANALYTIC_PHIS_INL_

namespace medusa {

template<bool Spline,
		 typename ArgTypeTime,
		 typename ArgTypeThetah,
		 typename ArgTypeThetal,
		 typename ArgTypePhi,
		 typename ArgTypeQOS,
         typename ArgTypeQSS,
         typename ArgTypeEtaOS,
         typename ArgTypeEtaSS,
		 typename ArgTypeDelta,
		 typename Signature >
void FullAnalyticPhis<Spline, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
						ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta, Signature>::Update_ATCoefficients()
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

	//
	// data map:
	//
	// data[0]: sin(_par[16]), data[1]: cos(_par[16])
	// data[2]: sin(_par[15]), data[3]: cos(_par[15])
	// data[4]: sin(_par[14]), data[5]: cos(_par[14])
	// data[6]: sin(_par[17]), data[7]: cos(_par[17])
	//
	// data[8]:  sin(_par[8]), data[9]: cos(_par[8])
	// data[10]: sin(_par[7]), data[11]: cos(_par[7])
	// data[12]: sin(_par[6]), data[13]: cos(_par[6])
	// data[14]: sin(_par[9]), data[15]: cos(_par[9])
	//
	// data[16]: _par[10]    , data[17]: _par[11]
	// data[18]: _par[12]    , data[19]: _par[13]

	double sphiz   = ::sin(_par[6]);//data[12];
	double cphiz   = ::cos(_par[6]);//data[13];
	double sphiy   = ::sin(_par[7]);//data[10];
	double cphiy   = ::cos(_par[7]);//data[11];
	double sphix   = ::sin(_par[8]);//data[8] ;
	double cphix   = ::cos(_par[8]);//data[9];
	double sphiw   = ::sin(_par[9]);//data[14];
	double cphiw   = ::cos(_par[9]);//data[14];
	double lambda_0    = _par[10];//data[16];
	double lambda_par  = _par[11];//data[17];
	double lambda_perp = _par[12];//data[18];
	double lambda_S    = _par[13];//data[19];
	double sdeltaz = ::sin(_par[14]);//data[4];
	double cdeltaz = ::cos(_par[14]);//data[5];
	double sdeltay = ::sin(_par[15]);//data[2];
	double cdeltay = ::cos(_par[15]);//data[3];
	double sdeltax = ::sin(_par[16]);//data[0];
	double cdeltax = ::cos(_par[16]) ;//data[1];
	double sdeltaw = ::sin(_par[17]);//data[6];
	double cdeltaw = ::cos(_par[17]);//data[7];

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
	//======================================================
    //                !!! WARNNING !!!
	//------------------------------------------------------
	//DO NOT DELETE ANY OF THE FOLLOWING COMMENTED LINES
	//MAYBE THEY WILL BE NECESSARY...
	//======================================================
	//::sin(delta_perp - delta_par)
	double sdelta_xmy = sdeltax*cdeltay - cdeltax*sdeltay;//ok
	//::cos(delta_perp - delta_par)
	double cdelta_xmy = cdeltax*cdeltay + sdeltax*sdeltay;

	//::cos(delta_0 - delta_par)
	double cdelta_zmy = cdeltaz*cdeltay + sdeltaz*sdeltay;
	//::sin(delta_0 - delta_par)
	double sdelta_zmy = sdeltaz*cdeltay - cdeltaz*sdeltay;//ok

	//::sin(delta_0 - delta_perp)
	double sdelta_zmx = sdeltaz*cdeltax - cdeltaz*sdeltax;
	//::cos(delta_0 - delta_perp)
	double cdelta_zmx = cdeltaz*cdeltax + sdeltaz*sdeltax;

	//::sin(delta_S - delta_perp)
	double sdelta_wmx = sdeltaw*cdeltax - cdeltaw*sdeltax;//ok
	//::cos(delta_S - delta_perp)
	double cdelta_wmx = cdeltaw*cdeltax + sdeltaw*sdeltax;//ok


	//::sin(delta_S - delta_par)
	double sdelta_wmy = sdeltaw*cdeltay - cdeltaw*sdeltay;//ok
	//::cos(delta_S - delta_par)
	double cdelta_wmy = cdeltaw*cdeltay + sdeltaw*sdeltay;//ok

	//::sin(delta_S - delta_0)
	double sdelta_wmz = sdeltaw*cdeltaz - cdeltaw*sdeltaz;//ok
	//::cos(delta_S - delta_0)
	double cdelta_wmz = cdeltaw*cdeltaz + sdeltaw*sdeltaz;


	//::sin(delta_perp + delta_par)
	//double sdelta_xpy = sdeltax*cdeltay + cdeltax*sdeltay;
	//::cos(delta_0 + delta_par)
	//double cdelta_zpy = cdeltaz*cdeltay - sdeltaz*sdeltay;
	//::sin(delta_0 + delta_perp)
	//double sdelta_zpx = sdeltaz*cdeltax + cdeltaz*sdeltax;
	//::cos(delta_S + delta_par)
	//double cdelta_wpy = cdeltaz*cdeltay - sdeltaz*sdeltay;
	//::sin(delta_S + delta_perp)
	//double sdelta_wpx = sdeltaw*cdeltax + cdeltaw*sdeltax;
	//::cos(delta_S + delta_0)
	//double cdelta_wpz = cdeltaw*cdeltaz - sdeltaw*sdeltaz;

	//::sin(phi_perp + phi_par)
	double sphi_xpy = sphix*cphiy + cphix*sphiy;
	//::cos(phi_perp + phi_par)
	double cphi_xpy = cphix*cphiy - sphix*sphiy;
	//::sin(phi_perp - phi_par)
	double sphi_xmy = sphix*cphiy - cphix*sphiy;
	//::cos(phi_perp - phi_par)
	double cphi_xmy = cphix*cphiy + sphix*sphiy;
	//::sin(phi_0 - phi_perp)
	double sphi_zmx = sphiz*cphix - cphiz*sphix;
	//::cos(phi_0 - phi_perp)
	double cphi_zmx = cphiz*cphix + sphiz*sphix;




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
	//::sin(phi_S - phi_perp)
	double sphi_wmx = sphiw*cphix - cphiw*sphix;
	//::cos(phi_S - phi_perp)
	double cphi_wmx = cphiw*cphix + sphiw*sphix;

	//::sin(phi_S + phi_0)
	//double sphi_wpz = sphiw*cphiz + cphiw*sphiz;
	//::cos(phi_S + phi_0)
	double cphi_wpz = cphiw*cphiz - sphiw*sphiz;
	//::sin(phi_S - phi_0)
	double sphi_wmz = sphiw*cphiz - cphiw*sphiz;
	//::cos(phi_S - phi_0)
	double cphi_wmz = cphiw*cphiz + sphiw*sphiz;

	//::sin[(delta_perp - delta_par) - (phi_perp + phi_par) ]
	double sdelta_xmy_m_phi_xmy = sdelta_xmy*cphi_xmy - cdelta_xmy*sphi_xmy;

	//::cos(delta_0 - delta_par - phi_0 + phi_par) = ::cos( [delta_0 - delta_par] - [phi_0 - phi_par] )
	//
	double cdelta_zmy_m_phi_zmy = cdelta_zmy*cphi_zmy + sdelta_zmy*sphi_zmy;

	//::sin(delta_0 - delta_perp - phi_0 + phi_perp) = ::sin([delta_0 - delta_perp] - [phi_0 - phi_perp])
	double sdelta_zmx_m_phi_zmx = sdelta_zmx*cphi_zmx - cdelta_zmx*sphi_zmx;//ok

	//::cos(delta_S - delta_par - phi_S + phi_par) = ::cos([delta_S - delta_par] - [phi_S - phi_par])
	double cdelta_wmy_m_phi_wmy = cdelta_wmy*cphi_wmy + sdelta_wmy*sphi_wmy;//ok

	//::sin(delta_S - delta_perp - phi_S + phi_perp) = ::sin(delta_S - delta_perp - [phi_S - phi_perp])
	double sdelta_wmx_m_phi_wmx = sdelta_wmx*cphi_wmx - cdelta_wmx*sphi_wmx;//ok

	//::cos(delta_S - delta_0 - phi_S + phi_0) = ::cos(delta_S - delta_0 - [phi_S - phi_0])
	double cdelta_wmz_m_phi_wmz = cdelta_wmz*cphi_wmz + sdelta_wmz*sphi_wmz;//ok

	//::sin(delta_perp - delta_par - phi_perp)
	double sdelta_xmy_m_phix = sdelta_xmy*cphix - cdelta_xmy*sphix;
	//::cos(delta_perp - delta_par - phi_perp)
	double cdelta_xmy_m_phix = cdelta_xmy*cphix + sdelta_xmy*sphix;

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
	double cdelta_wmy_m_phiw = cdelta_wmy*cphiw + sdelta_wmy*sphiw;//ok
	//::sin(delta_S - delta_par - phi_S)
	double sdelta_wmy_m_phiw = sdelta_wmy*cphiw - cdelta_wmy*sphiw;//ok

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
	// A_k
	//--------------------------------------------------------------------

	A.k[0] = 0.5*(1.0 + lambda_0*lambda_0);

	A.k[1] = 0.5*(1.0 + lambda_par*lambda_par);

	A.k[2] = 0.5*(1.0 + lambda_perp*lambda_perp);

	//0.5*(::sin(delta_perp - delta_par) - lambda_par*lambda_perp*::sin(delta_perp - delta_par - phi_perp + phi_par) )
	A.k[3] = 0.5*(sdelta_xmy - lambda_par*lambda_perp*sdelta_xmy_m_phi_xmy);

	//0.5*(::cos(delta_0 - delta_par) + lambda_0*lambda_par*::cos(delta_0 - delta_par - phi_0 + phi_par) )
	A.k[4] = 0.5*(cdelta_zmy + lambda_0*lambda_par*cdelta_zmy_m_phi_zmy );

	//-0.5*(::sin(delta_0 - delta_perp) - lambda_0*lambda_perp*::sin(delta_0 - delta_perp - phi_0 + phi_perp))
	A.k[5] = -0.5*(sdelta_zmx - lambda_0*lambda_perp*sdelta_zmx_m_phi_zmx );

	//0.5*(1 + ::pow(lambda_S,2) )
	A.k[6] = 0.5*(1 + lambda_S*lambda_S );

	//0.5*(::cos(delta_S - delta_par) - lambda_S*lambda_par*::cos(delta_S - delta_par - phi_S + phi_par))
	A.k[7] = 0.5*(cdelta_wmy - lambda_S*lambda_par*cdelta_wmy_m_phi_wmy );

	//-0.5*(::sin(delta_S - delta_perp) + lambda_S*lambda_perp*::sin(delta_S - delta_perp - phi_S + phi_perp))
	A.k[8] = -0.5*(sdelta_wmx + lambda_S*lambda_perp*sdelta_wmx_m_phi_wmx);

	//0.5*(::cos(delta_S - delta_0) - lambda_S*lambda_0*::cos(delta_S - delta_0 - phi_S + phi_0))
	A.k[9] = 0.5*(cdelta_wmz - lambda_S*lambda_0*cdelta_wmz_m_phi_wmz);

	//--------------------------------------------------------------------
	// B_k
	//--------------------------------------------------------------------

	// -lambda_0*::cos(phi_0)
	B.k[0] = -lambda_0*cphiz;

	// -lambda_par*::cos(phi_par)
	B.k[1] = -lambda_par*cphiy;

	//lambda_perp*::cos(phi_perp)
	B.k[2] = lambda_perp*cphix ;

	//0.5*(lambda_perp*::sin(delta_perp - delta_par - phi_perp) + lambda_par*::sin(delta_par - delta_perp - phi_par))
	B.k[3] = 0.5*(lambda_perp*sdelta_xmy_m_phix + lambda_par*sdelta_ymx_m_phiy);

	//-0.5*(lambda_0*::cos(delta_0 - delta_par - phi_0) + lambda_par*::cos(delta_par - delta_0 - phi_par))
	B.k[4] = -0.5*(lambda_0*cdelta_zmy_m_phiz + lambda_par*cdelta_ymz_m_phiy);

	//0.5*(lambda_0*::sin(delta_0 - delta_perp - phi_0) + lambda_perp*::sin(delta_perp - delta_0 - phi_perp))
	B.k[5] = 0.5*(lambda_0*sdelta_zmx_m_phiz + lambda_perp*sdelta_xmz_m_phix);

	// lambda_S*::cos(phi_S)
	B.k[6] =  lambda_S*cphiw;

	//0.5*(lambda_S*::cos(delta_S - delta_par - phi_S) - lambda_par*::cos(delta_par - delta_S - phi_par))
	B.k[7] = 0.5*(lambda_S*cdelta_wmy_m_phiw - lambda_par*cdelta_ymw_m_phiy);

	//-0.5*(lambda_S*::sin(delta_S - delta_perp - phi_S) - lambda_perp*::sin(delta_perp - delta_S - phi_perp))
	B.k[8] = -0.5*(lambda_S*sdelta_wmx_m_phiw - lambda_perp*sdelta_xmw_m_phix );

	// 0.5*(lambda_S*::cos(delta_S - delta_0 - phi_S) - lambda_0*::cos(delta_0 - delta_S - phi_0))
	B.k[9] =  0.5*(lambda_S*cdelta_wmz_m_phiw - lambda_0*cdelta_zmw_m_phiz);

	//--------------------------------------------------------------------
	// C_k
	//--------------------------------------------------------------------

	//0.5*(1 - pow(lambda_0,2))
	C.k[0] = 0.5*(1 - lambda_0*lambda_0);

	//0.5*(1 - pow(lambda_par,2))
	C.k[1] = 0.5*(1 - lambda_par*lambda_par);

	//0.5*(1 - pow(lambda_perp,2))
	C.k[2] = 0.5*(1 - lambda_perp*lambda_perp);

	//0.5*(::sin(delta_perp - delta_par) + lambda_par*lambda_perp*::sin(delta_perp - delta_par - phi_perp + phi_par) )
	C.k[3] = 0.5*( sdelta_xmy + lambda_par*lambda_perp*sdelta_xmy_m_phi_xmy );

	//0.5*(::cos(delta_0 - delta_par) - lambda_0*lambda_par*::cos(delta_0 - delta_par - phi_0 + phi_par) )
	C.k[4] = 0.5*(cdelta_zmy - lambda_0*lambda_par*cdelta_zmy_m_phi_zmy  );

	//-0.5*(::sin(delta_0 - delta_perp) + lambda_0*lambda_perp*::sin(delta_0 - delta_perp - phi_0 + phi_perp))
	C.k[5] = -0.5*(sdelta_zmx + lambda_0*lambda_perp*sdelta_zmx_m_phi_zmx);

	//0.5*(1 - ::pow(lambda_S,2) )
	C.k[6] = 0.5*(1 - lambda_S*lambda_S );

	//0.5*(::cos(delta_S - delta_par) + lambda_S*lambda_par*::cos(delta_S - delta_par - phi_S + phi_par))
	C.k[7] = 0.5*(cdelta_wmy + lambda_S*lambda_par*cdelta_wmy_m_phi_wmy);

	//-0.5*(::sin(delta_S - delta_perp) - lambda_S*lambda_perp*::sin(delta_S - delta_perp - phi_S + phi_perp))
	C.k[8] = -0.5*(sdelta_wmx - lambda_S*lambda_perp*sdelta_wmx_m_phi_wmx);

	//0.5*(::cos(delta_S - delta_0) + lambda_S*lambda_0*::cos(delta_S - delta_0 - phi_S + phi_0))
	C.k[9] = 0.5*(cdelta_wmz + lambda_S*lambda_0*cdelta_wmz_m_phi_wmz);

	//--------------------------------------------------------------------
	// D_k
	//--------------------------------------------------------------------

	//lambda_0*::sin(phi_0)
	D.k[0] = lambda_0*sphiz;

	//lambda_par*::sin(phi_par)
	D.k[1] = lambda_par*sphiy;

	//-lambda_perp*::sin(phi_perp)
	D.k[2] = -lambda_perp*sphix;

	// -0.5*(lambda_perp*::cos(delta_perp - delta_par - phi_perp) + lambda_par*::cos(delta_par - delta_perp - phi_par))
	D.k[3] = -0.5*(lambda_perp*cdelta_xmy_m_phix + lambda_par*cdelta_ymx_m_phiy);

	// -0.5*(lambda_0*::sin(delta_0 - delta_par - phi_0) + lambda_par*::sin(delta_par - delta_0 - phi_par))
	D.k[4] = -0.5*(lambda_0*sdelta_zmy_m_phiz + lambda_par*sdelta_ymz_m_phiy);

	//-0.5*(lambda_0*::cos(delta_0 - delta_perp - phi_0) + lambda_perp*::cos(delta_perp - delta_0 - phi_perp))
	D.k[5] = -0.5*(lambda_0*cdelta_zmx_m_phiz + lambda_perp*cdelta_xmz_m_phix);

	//-lambda_S*::sin(phi_S)
	D.k[6] = -lambda_S*sphiw;

	//0.5*(lambda_S*::sin(delta_S - delta_par - phi_S) - lambda_par*::sin(delta_par - delta_S - phi_par))
	D.k[7] =  0.5*(lambda_S*sdelta_wmy_m_phiw - lambda_par*sdelta_ymw_m_phiy);

	//-0.5*(-lambda_S*::cos(delta_S - delta_perp - phi_S) + lambda_perp*::cos(delta_perp - delta_S - phi_perp))
	D.k[8] = -0.5*(-lambda_S*cdelta_wmx_m_phiw + lambda_perp*cdelta_xmw_m_phix );

	//0.5*(lambda_S*::sin(delta_S - delta_0 - phi_S) - lambda_0*::sin(delta_0 - delta_S - phi_0))
	D.k[9]  = 0.5*(lambda_S*sdelta_wmz_m_phiw - lambda_0*sdelta_zmw_m_phiz);

}


}  // namespace medusa


#endif /* FULL_ANALYTIC_PHIS_INL_ */
