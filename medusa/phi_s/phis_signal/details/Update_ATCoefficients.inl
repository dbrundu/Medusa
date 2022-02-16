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
 *  PhisSignal.inl
 *
 *  Created on: 24/07/2020
 *      Author: Antonio Augusto Alves Junior
 * 		Updated by Alessandro Maria Ricci in 04/06/2021
 * 
 * 	This function updates the values of the angular coefficients a_k, b_k,
 *  c_k, d_k by using the formulas in Table 3 in arXiv:1906.08356v4.
 *---------------------------------------------------------------------------*/

#ifndef PHIS_SIGNAL_INL_
#define PHIS_SIGNAL_INL_

namespace medusa {

template<bool B0sbar,
		 typename ArgTypeTime,
		 typename ArgTypeCosThetah,
		 typename ArgTypeCosThetal,
		 typename ArgTypePhi,
		 typename Signature >
void PhisSignal<B0sbar, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi, Signature>::Update_ATCoefficients()
{
	/*
	{
	 0: A_0^2,
	 1: A_perp^2,
	 2: A_S^2,
	 3: DeltaGamma_sd,
	 4: DeltaGamma,
	 5: DeltaM,
	 6: phi_0,
	 7: phi_par0: phi_par - phi_0,
	 8: phi_perp0: phi_perp - phi_0,
	 9: phi_S0: phi_S - phi_0,
	 10: lambda_0,
	 11: lambda_par0: lambda_par / lambda_0,
	 12: lambda_perp0: lambda_perp / lambda_0,
	 13: lambda_S0: lambda_S / lambda_0,
	 14: delta_par0: delta_par - delta_0,
	 15: delta_perp0: delta_perp - delta_0,
	 16: delta_Sperp: delta_S - delta_perp
	}
	*/

	//
	// data map:
	//
	// data[0]: sin(_par[15]), data[1]: cos(_par[15])
	// data[2]: sin(_par[14]), data[3]: cos(_par[14])
	// data[4]: sin(_par[16]), data[5]: cos(_par[16])
	//
	// data[6]:  sin(_par[8]), data[7]: cos(_par[8])
	// data[8]:  sin(_par[7]), data[9]: cos(_par[7])
	// data[10]: sin(_par[6]), data[11]: cos(_par[6])
	// data[12]: sin(_par[9]), data[13]: cos(_par[9])
	//
	// data[14]: _par[10]    , data[15]: _par[11]
	// data[16]: _par[12]    , data[17]: _par[13]

	double sphiz   = ::sin(_par[6]);//data[10];
	double cphiz   = ::cos(_par[6]);//data[11];
	double sphiy   = ::sin(_par[7]);//data[8];
	double cphiy   = ::cos(_par[7]);//data[9];
	double sphix   = ::sin(_par[8]);//data[6];
	double cphix   = ::cos(_par[8]);//data[7];
	double sphiw   = ::sin(_par[9]);//data[12];
	double cphiw   = ::cos(_par[9]);//data[13];
	double lambda_0     = _par[10];//data[14];
	double lambda_par0  = _par[11];//data[15];
	double lambda_perp0 = _par[12];//data[16];
	double lambda_S0    = _par[13];//data[17];
	double sdeltay = ::sin(_par[14]);//data[2];
	double cdeltay = ::cos(_par[14]);//data[3];
	double sdeltax = ::sin(_par[15]);//data[0];
	double cdeltax = ::cos(_par[15]) ;//data[1];
	double sdeltaw = ::sin(_par[16]);//data[4];
	double cdeltaw = ::cos(_par[16]);//data[5];

	// using the infamous identities
	//
	//sin(α – β) = sin(α) cos(β) – cos(α) sin(β)
	//sin(α + β) = sin(α) cos(β) + cos(α) sin(β)
	//
	//cos(α + β) = cos(α) cos(β) – sin(α) sin(β)
	//cos(α – β) = cos(α) cos(β) + sin(α) sin(β)
	//
	//cos(-α) =  cos(α)
	//sin(-α) = -sin(α)
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
	//::sin(delta_perp0 - delta_par0)
	double sdelta_xmy = sdeltax*cdeltay - cdeltax*sdeltay;//ok
	//::cos(delta_perp0 - delta_par0)
	double cdelta_xmy = cdeltax*cdeltay + sdeltax*sdeltay;

	//::cos(delta_Sperp + delta_perp0)
	double cdelta_wmx = cdeltaw*cdeltax - sdeltaw*sdeltax;
	//::sin(delta_Sperp + delta_perp0)
	double sdelta_wmx = sdeltaw*cdeltax + cdeltaw*sdeltax;


	//::cos(phi_par0 + phi_0)
	double cphi_ymz = cphiy*cphiz - sphiy*sphiz;//ok
	//::sin(phi_par0 + phi_0)
	double sphi_ymz = sphiy*cphiz + cphiy*sphiz;//ok

	//::cos(phi_perp0 + phi_0)
	double cphi_xmz = cphix*cphiz - sphix*sphiz;//ok
	//::sin(phi_perp0 + phi_0)
	double sphi_xmz = sphix*cphiz + cphix*sphiz;//ok

	//::sin(phi_par0 - phi_perp0)
	double sphi_ymx = sphiy*cphix - cphiy*sphix;
	//::cos(phi_par0 - phi_perp0)
	double cphi_ymx = cphiy*cphix + sphiy*sphix;

	//::sin(phi_par0 - phi_S0)
	double sphi_ymw = sphiy*cphiw - cphiy*sphiw;
	//::cos(phi_par0 - phi_S0)
	double cphi_ymw = cphiy*cphiw + sphiy*sphiw;

	//::sin(phi_perp0 - phi_S0)
	double sphi_xmw = sphix*cphiw - cphix*sphiw;
	//::cos(phi_perp0 - phi_S0)
	double cphi_xmw = cphix*cphiw + sphix*sphiw;

	//::cos(phi_S0 + phi_0)
	double cphi_wmz = cphiw*cphiz - sphiw*sphiz;
	//::sin(phi_S0 + phi_0)
	double sphi_wmz = sphiw*cphiz + cphiw*sphiz;


	//::sin(delta_Sperp + delta_perp0 - delta_par0) = ::sin( delta_Sperp + (delta_perp0 - delta_par0) )
	double sdelta_xmy_m_w = sdeltaw*cdelta_xmy + cdeltaw*sdelta_xmy;//ok
	//::cos(delta_Sperp + delta_perp0 - delta_par0) = ::cos( delta_Sperp + (delta_perp0 - delta_par0) )
	double cdelta_xmy_m_w = cdeltaw*cdelta_xmy - sdeltaw*sdelta_xmy;//ok

	//::sin(- delta_Sperp + delta_par0 - delta_perp0) = ::sin( (delta_par0 - delta_perp0) - deltaSperp )
	double sdelta_ymx_m_w = -sdelta_xmy*cdeltaw - cdelta_xmy*sdeltaw;//ok
	//::cos(- delta_Sperp + delta_par0 - delta_perp0) = ::cos( (delta_par0 - delta_perp0) - deltaSperp )
	double cdelta_ymx_m_w = cdelta_xmy*cdeltaw - sdelta_xmy*sdeltaw;//ok

	//::sin(delta_Sperp + phi_perp0 - phi_S0) = ::sin( delta_Sperp + ( phi_perp0 - phi_S0) )
	double sdeltaw_m_phi_xmw = sdeltaw*cphi_xmw + cdeltaw*sphi_xmw;//ok

	//::sin( (delta_perp0 - delta_par0) + (phi_par0 - phi_perp0) )
	double sdelta_xmy_m_phi_ymx = sdelta_xmy*cphi_ymx + cdelta_xmy*sphi_ymx;

	//::cos(phi_par0 - delta_par0)
	double cdeltay_m_phiy = cphiy*cdeltay + sdeltay*sphiy;

	//::sin(phi_perp0 - delta_perp0)
	double sdeltax_m_phix = sphix*cdeltax - cphix*sdeltax;//ok

	//::cos( delta_Sperp + ( delta_perp0 - delta_par0 ) + ( phi_par0 - phi_S0 ) )
	double cdelta_xmy_m_phi_xmy_m_w = cdelta_xmy_m_w*cphi_ymw - sdelta_xmy_m_w*sphi_ymw;

	//::cos( (delta_Sperp + delta_perp0) - phi_S0 )
	double cdelta_wmx_m_phiw = cdelta_wmx*cphiw + sdelta_wmx*sphiw;

	//::sin( (delta_perp0 - delta_par0) - (phi_perp0 + phi0) )
	double sdelta_xmy_m_phi_xmz = sdelta_xmy*cphi_xmz - cdelta_xmy*sphi_xmz;
	//::cos( (delta_perp0 - delta_par0) - (phi_perp0 + phi0) )
	double cdelta_xmy_m_phi_xmz = cdelta_xmy*cphi_xmz + sdelta_xmy*sphi_xmz;

	//::sin( (delta_par0 - delta_perp0) - (phi_par0 + phi0) )
	double sdelta_ymx_m_phi_ymz = -sdelta_xmy*cphi_ymz - cdelta_xmy*sphi_ymz;
	//::cos( (delta_par0 - delta_perp0) - (phi_par0 + phi0) )
	double cdelta_ymx_m_phi_ymz = cdelta_xmy*cphi_ymz - sdelta_xmy*sphi_ymz;

	//::sin(- delta_par0 - phi_0)
	double sdeltay_m_phiz = -sdeltay*cphiz - cdeltay*sphiz;
	//::cos(- delta_par0 - phi_0)
	double cdeltay_m_phiz = cdeltay*cphiz - sdeltay*sphiz;

	//::sin(- delta_perp0 - phi_0)
	double sdeltax_m_phiz = -sdeltax*cphiz - cdeltax*sphiz;
	//::cos(- delta_perp0 - phi_0)
	double cdeltax_m_phiz = cdeltax*cphiz - sdeltax*sphiz;

	//::sin( delta_par0 - (phi_par0 + phi_0) )
	double sdeltay_m_phi_ymz = sdeltay*cphi_ymz - cdeltay*sphi_ymz;
	//::cos( delta_par0 - (phi_par0 + phi_0) )
	double cdeltay_m_phi_ymz = cdeltay*cphi_ymz + sdeltay*sphi_ymz;

	//::sin(delta_perp0 - (phi_perp0 + phi_0) )
	double sdeltax_m_phi_xmz = sdeltax*cphi_xmz - cdeltax*sphi_xmz;
	//::cos(delta_perp0 - (phi_perp0 + phi_0) )
	double cdeltax_m_phi_xmz = cdeltax*cphi_xmz + sdeltax*sphi_xmz;

	//::sin( (delta_Sperp + (delta_perp0 - delta_par0)) - (phi_S0 + phi_0) )
	double sdeltaw_m_xmy_m_phi_wmz = sdelta_xmy_m_w*cphi_wmz - cdelta_xmy_m_w*sphi_wmz;
	//::cos( (delta_Sperp + (delta_perp0 - delta_par0)) - (phi_S0 + phi_0) )
	double cdeltaw_m_xmy_m_phi_wmz = cdelta_xmy_m_w*cphi_wmz + sdelta_xmy_m_w*sphi_wmz;

	//::sin( ((delta_par0 - delta_perp0) - delta_Sperp) - (phi_par0 + phi_0) )
	double sdeltaw_m_ymx_m_phi_ymz = sdelta_ymx_m_w*cphi_ymz - cdelta_ymx_m_w*sphi_ymz;
	//::cos( ((delta_par0 - delta_perp0) - delta_Sperp) - (phi_par0 + phi_0) )
	double cdeltaw_m_ymx_m_phi_ymz = cdelta_ymx_m_w*cphi_ymz + sdelta_ymx_m_w*sphi_ymz;

	//::sin( delta_Sperp - (phi_S0 + phi_0) )
	double sdeltaw_m_phi_wmz = sdeltaw*cphi_wmz - cdeltaw*sphi_wmz;
	//::cos( delta_Sperp - (phi_S0 + phi_0) )
	double cdeltaw_m_phi_wmz = cdeltaw*cphi_wmz + sdeltaw*sphi_wmz;

	//::sin(- delta_Sperp - (phi_perp0 + phi_0) )
	double sdeltaw_m_phi_xmz = -sdeltaw*cphi_xmz - cdeltaw*sphi_xmz;
	//::cos(- delta_Sperp - (phi_perp0 + phi_0) )
	double cdeltaw_m_phi_xmz = cdeltaw*cphi_xmz - sdeltaw*sphi_xmz;

	//::sin( (delta_Sperp + delta_perp0) - (phi_S0 + phi_0) )
	double sdelta_wmx_m_phi_wmz = sdelta_wmx*cphi_wmz - cdelta_wmx*sphi_wmz;
	//::cos( (delta_Sperp + delta_perp0) - (phi_S0 + phi_0) )
	double cdelta_wmx_m_phi_wmz = cdelta_wmx*cphi_wmz + sdelta_wmx*sphi_wmz;

	//::sin( -(delta_perp0 + delta_Sperp) - phi_0 )
	double sdelta_wmx_m_phiz = -sdelta_wmx*cphiz - cdelta_wmx*sphiz;
	//::cos( -(delta_perp0 + delta_Sperp) - phi_0 )
	double cdelta_wmx_m_phiz = cdelta_wmx*cphiz - sdelta_wmx*sphiz;


	//--------------------------------------------------------------------
	// A_k
	//--------------------------------------------------------------------

	A.k[0] = 0.5*(1.0 + lambda_0*lambda_0);

	A.k[1] = 0.5*(1.0 + lambda_par0*lambda_par0*lambda_0*lambda_0);

	A.k[2] = 0.5*(1.0 + lambda_perp0*lambda_perp0*lambda_0*lambda_0);

	//0.5*(::sin(delta_perp0 - delta_par0) - lambda_par0*lambda_perp0*lambda_0*lambda_0*
	//																	::sin( (delta_perp0 - delta_par0) + (phi_par0 - phi_perp0) ) )
	A.k[3] = 0.5*( sdelta_xmy - lambda_par0*lambda_perp0*lambda_0*lambda_0*sdelta_xmy_m_phi_ymx );

	//0.5*( ::cos(- delta_par0) + lambda_0*lambda_0*lambda_par0*::cos(phi_par0 - delta_par0) ) 
	A.k[4] = 0.5*( cdeltay + lambda_0*lambda_0*lambda_par0*cdeltay_m_phiy );

	//-0.5*(::sin(- delta_perp0) - lambda_0*lambda_0*lambda_perp0*::sin(phi_perp0 - delta_perp0))
	A.k[5] = -0.5*(-sdeltax - lambda_0*lambda_0*lambda_perp0*sdeltax_m_phix );

	//0.5*(1 + ::pow(lambda_S0*lambda_0,2) )
	A.k[6] = 0.5*(1 + lambda_S0*lambda_S0*lambda_0*lambda_0 );

	//0.5*( ::cos(delta_Sperp + delta_perp0 - delta_par0) - lambda_S0*lambda_par0*lambda_0*lambda_0*
	// 														::cos( delta_Sperp + (delta_perp0 - delta_par0) + (phi_par0 - phi_S0) ) )
	A.k[7] = 0.5*(cdelta_xmy_m_w - lambda_S0*lambda_par0*lambda_0*lambda_0*cdelta_xmy_m_phi_xmy_m_w );

	//-0.5*( ::sin(delta_Sperp) + lambda_S0*lambda_perp0*lambda_0*lambda_0*::sin( delta_Sperp + (phi_perp0  - phi_S0) ) )
	A.k[8] = -0.5*(sdeltaw + lambda_S0*lambda_perp0*lambda_0*lambda_0*sdeltaw_m_phi_xmw);

	//0.5*( ::cos(delta_Sperp + delta_perp0) - lambda_S0*lambda_0*lambda_0*::cos( (delta_Sperp + delta_perp0) - phi_S0 ) )
	A.k[9] = 0.5*(cdelta_wmx - lambda_S0*lambda_0*lambda_0*cdelta_wmx_m_phiw);


	//--------------------------------------------------------------------
	// B_k
	//--------------------------------------------------------------------

	// -lambda_0*::cos(phi_0)
	B.k[0] = -lambda_0*cphiz;

	// -lambda_par0*lambda_0*::cos(phi_par0 + phi_0)
	B.k[1] = -lambda_par0*lambda_0*cphi_ymz;

	//lambda_perp0*lambda_0*::cos(phi_perp0 + phi_0)
	B.k[2] = lambda_perp0*lambda_0*cphi_xmz;

	//0.5*( lambda_perp0*lambda_0*::sin( (delta_perp0 - delta_par0) - (phi_perp0 + phi0) ) + lambda_par0*lambda_0*
	//																			::sin( (delta_par0 - delta_perp0) - (phi_par0 + phi0) ) )
	B.k[3] = 0.5*(lambda_perp0*lambda_0*sdelta_xmy_m_phi_xmz + lambda_par0*lambda_0*sdelta_ymx_m_phi_ymz);

	//-0.5*(lambda_0*::cos(- delta_par0 - phi_0) + lambda_par0*lambda_0*::cos( delta_par0 - (phi_par0 + phi_0) ) )
	B.k[4] = -0.5*(lambda_0*cdeltay_m_phiz + lambda_par0*lambda_0*cdeltay_m_phi_ymz);

	//0.5*(lambda_0*::sin(- delta_perp0 - phi_0) + lambda_perp0*lambda_0*::sin(delta_perp0 - (phi_perp0 + phi_0) ) )
	B.k[5] = 0.5*(lambda_0*sdeltax_m_phiz + lambda_perp0*lambda_0*sdeltax_m_phi_xmz);

	// lambda_S0*lambda_0*::cos(phi_S0 + phi_0)
	B.k[6] = lambda_S0*lambda_0*cphi_wmz;

	//0.5*( lambda_S0*lambda_0*::cos( (delta_Sperp + (delta_perp0 - delta_par0)) - (phi_S0 + phi_0) ) -
	//											lambda_par0*lambda_0*::cos( ((delta_par0 - delta_perp0) - delta_Sperp) - (phi_par0 + phi_0) ) )
	B.k[7] = 0.5*(lambda_S0*lambda_0*cdeltaw_m_xmy_m_phi_wmz - lambda_par0*lambda_0*cdeltaw_m_ymx_m_phi_ymz);

	//-0.5*( lambda_S0*lambda_0*::sin( delta_Sperp - (phi_S0 + phi_0) ) - lambda_perp0*lambda_0*::sin(- delta_Sperp - (phi_perp0 + phi_0) ) )
	B.k[8] = -0.5*(lambda_S0*lambda_0*sdeltaw_m_phi_wmz - lambda_perp0*lambda_0*sdeltaw_m_phi_xmz );

	// 0.5*( lambda_S0*lambda_0*::cos( (delta_Sperp + delta_perp0) - (phi_S0 + phi_0) ) - lambda_0*::cos(- (delta_perp0 + delta_Sperp) - phi_0 ) )
	B.k[9] =  0.5*(lambda_S0*lambda_0*cdelta_wmx_m_phi_wmz - lambda_0*cdelta_wmx_m_phiz);


	//--------------------------------------------------------------------
	// C_k
	//--------------------------------------------------------------------

	//0.5*(1 - pow(lambda_0,2))
	C.k[0] = 0.5*(1 - lambda_0*lambda_0);

	//0.5*( 1 - pow(lambda_par0, 2) * pow(lambda_0, 2) )
	C.k[1] = 0.5*(1 - lambda_par0*lambda_par0*lambda_0*lambda_0);

	//0.5*( 1 - pow(lambda_perp0, 2) * pow(lambda_0, 2) )
	C.k[2] = 0.5*(1 - lambda_perp0*lambda_perp0*lambda_0*lambda_0);

	//0.5*( ::sin(delta_perp0 - delta_par0) + lambda_par0*lambda_perp0*lambda_0*lambda_0*
	//																		::sin( (delta_perp0 - delta_par0) + (phi_par0 - phi_perp0) ) )
	C.k[3] = 0.5*( sdelta_xmy + lambda_par0*lambda_perp0*lambda_0*lambda_0*sdelta_xmy_m_phi_ymx );

	//0.5*(::cos(- delta_par0) - lambda_0*lambda_0*lambda_par0*::cos(phi_par0 - delta_par0) )
	C.k[4] = 0.5*(cdeltay - lambda_0*lambda_0*lambda_par0*cdeltay_m_phiy  );

	//-0.5*(::sin(- delta_perp0) + lambda_0*lambda_0*lambda_perp0*::sin(phi_perp0 - delta_perp0) )
	C.k[5] = -0.5*(-sdeltax + lambda_0*lambda_0*lambda_perp0*sdeltax_m_phix);

	//0.5*(1 - ::pow(lambda_S0,2) * pow(lambda_0, 2) )
	C.k[6] = 0.5*(1 - lambda_S0*lambda_S0*lambda_0*lambda_0 );

	//0.5*( ::cos(delta_Sperp + delta_perp0 - delta_par0) + lambda_S0*lambda_par0*lambda_0*lambda_0*
	//															::cos( delta_Sperp + (delta_perp0 - delta_par0) + (phi_par0 - phi_S0) ) )
	C.k[7] = 0.5*(cdelta_xmy_m_w + lambda_S0*lambda_par0*lambda_0*lambda_0*cdelta_xmy_m_phi_xmy_m_w);

	//-0.5*( ::sin(delta_Sperp) - lambda_S0*lambda_perp0*lambda_0*lambda_0*::sin( delta_Sperp + (phi_perp0 - phi_S0) ) )
	C.k[8] = -0.5*(sdeltaw - lambda_S0*lambda_perp0*lambda_0*lambda_0*sdeltaw_m_phi_xmw);

	//0.5*( ::cos(delta_Sperp + delta_perp0) + lambda_S0*lambda_0*lambda_0*::cos( (delta_Sperp + delta_perp0) - phi_S0 ) )
	C.k[9] = 0.5*(cdelta_wmx + lambda_S0*lambda_0*lambda_0*cdelta_wmx_m_phiw);


	//--------------------------------------------------------------------
	// D_k
	//--------------------------------------------------------------------

	//lambda_0*::sin(phi_0)
	D.k[0] = lambda_0*sphiz;

	//lambda_par0*lambda_0*::sin(phi_par0 + phi_0)
	D.k[1] = lambda_par0*lambda_0*sphi_ymz;

	//-lambda_perp0*lambda_0*::sin(phi_perp0 + phi_0)
	D.k[2] = -lambda_perp0*lambda_0*sphi_xmz;

	// -0.5*( lambda_perp0*lambda_0*::cos( (delta_perp0 - delta_par0) - (phi_perp0 + phi0) ) + lambda_par*
	//																			::cos( (delta_par0 - delta_perp0) - (phi_par0 + phi0) ) )
	D.k[3] = -0.5*(lambda_perp0*lambda_0*cdelta_xmy_m_phi_xmz + lambda_par0*lambda_0*cdelta_ymx_m_phi_ymz);

	// -0.5*( lambda_0*::sin(- delta_par0 - phi_0) + lambda_par0*lambda_0*::sin( delta_par0 - (phi_par0 + phi_0) ) )
	D.k[4] = -0.5*(lambda_0*sdeltay_m_phiz + lambda_par0*lambda_0*sdeltay_m_phi_ymz);

	//-0.5*( lambda_0*::cos(- delta_perp0 - phi_0) + lambda_perp0*lambda_0*::cos( delta_perp0 - (phi_perp0 + phi_0) ) )
	D.k[5] = -0.5*(lambda_0*cdeltax_m_phiz + lambda_perp0*lambda_0*cdeltax_m_phi_xmz);

	//-lambda_S0*lambda_0*::sin(phi_S0 + phi_0)
	D.k[6] = -lambda_S0*lambda_0*sphi_wmz;

	//0.5*( lambda_S0*lambda_0*::sin( (delta_Sperp + (delta_perp0 - delta_par0)) - (phi_S0 + phi_0) ) - lambda_par0*lambda_0*
	//																::sin( ((delta_par0 - delta_perp0) - delta_Sperp) - (phi_par0 + phi_0) )
	D.k[7] =  0.5*(lambda_S0*lambda_0*sdeltaw_m_xmy_m_phi_wmz - lambda_par0*lambda_0*sdeltaw_m_ymx_m_phi_ymz);

	//-0.5*( -lambda_S0*lambda_0*::cos( delta_Sperp - (phi_S0 + phi_0) ) + lambda_perp0*lambda_0*::cos(- delta_Sperp - (phi_perp0 + phi_0) ) )
	D.k[8] = -0.5*(-lambda_S0*lambda_0*cdeltaw_m_phi_wmz + lambda_perp0*lambda_0*cdeltaw_m_phi_xmz );

	//0.5*( lambda_S0*lambda_0*::sin( (delta_Sperp + delta_perp0) - (phi_S0 + phi_0) ) - lambda_0*::sin( - (delta_perp0 + delta_Sperp) - phi_0 ) )
	D.k[9]  = 0.5*(lambda_S0*lambda_0*sdelta_wmx_m_phi_wmz - lambda_0*sdelta_wmx_m_phiz);

}


}  // namespace medusa


#endif /* PHIS_SIGNAL_INL_ */
