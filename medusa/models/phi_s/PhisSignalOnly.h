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
 *  Created on: 07/05/2020
 *      Author: Davide Brundu
 */


#ifndef PHISSIGNALONLY_H_
#define PHISSIGNALONLY_H_

#include <cmath>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <ratio>


#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/Tuple.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/Math.h>


#include<medusa/models/phi_s/detail/NFactors.h>
#include<medusa/models/phi_s/detail/AFactors.h>
#include<medusa/models/phi_s/detail/ATCoeficients.h>



namespace medusa {




/*
 *  @class PhisSignalOnly
 *  Funtor that provides the time dependent formula used in phi_s analysis
 *  The actual implementation is inside the detail/ folder
 *
 *  B0bar    = boolean to specify wether B is B0bar or not
 *  ArgTypes = argument types of the functor
 *
 */
template<bool B0bar, typename ArgTypeTime, typename ArgTypeThetah,
         typename ArgTypeThetal, typename ArgTypePhi,
         typename Signature=double(ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi) >
class PhisSignalOnly: public hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18>
{

    constexpr static int CP  =  (B0bar ? -1 : +1);
    
    using ThisBaseFunctor = hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18 >;

    using hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18 >::_par;
    
public:

    PhisSignalOnly() = delete;


    PhisSignalOnly(hydra::Parameter const& A_0,           hydra::Parameter const& A_perp,     hydra::Parameter const& A_S, 
                   hydra::Parameter const& DeltaGamma_sd, hydra::Parameter const& DeltaGamma, hydra::Parameter const& DeltaM,
                   hydra::Parameter const& phi_0,         hydra::Parameter const& phi_par,
                   hydra::Parameter const& phi_perp,      hydra::Parameter const& phi_S,
                   hydra::Parameter const& lambda_0,      hydra::Parameter const& lambda_par,
                   hydra::Parameter const& lambda_perp,   hydra::Parameter const& lambda_S,
                   hydra::Parameter const& delta_0,       hydra::Parameter const& delta_par,
                   hydra::Parameter const& delta_perp,  hydra::Parameter const& delta_S):
     ThisBaseFunctor({A_0, A_perp, A_S, DeltaGamma_sd, DeltaGamma, DeltaM ,
                      phi_0,       phi_par,    phi_perp,    phi_S,          lambda_0,   lambda_par, 
                      lambda_perp, lambda_S,   delta_0,     delta_par,      delta_perp, delta_S })
    {
       Update();
    }



    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit PhisSignalOnly( const hydra::Parameter (&Hs)[18] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17] }
    {
       Update();
    }
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit PhisSignalOnly( const double (&Hs)[18] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17] }
    {
       Update();
    }
    
    

    __hydra_dual__
    PhisSignalOnly( PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi> const& other):
    ThisBaseFunctor(other)
    {

    	#pragma unroll 10
    	for(size_t i=0; i<10;i++)
    	{
    		fA.fC[i] = 	other.GetA().fC[i];
    		fB.fC[i] = 	other.GetB().fC[i];
    		fC.fC[i] = 	other.GetC().fC[i];
    		fD.fC[i] = 	other.GetD().fC[i];
    		fN.fC[i] =  other.GetN().fC[i];

    	}

    }



    __hydra_dual__
    PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>& 
    operator=( PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);

		#pragma unroll 10
        for(size_t i=0; i<10;i++)
        {
        	fA.fC[i] = 	other.GetA().fC[i];
        	fB.fC[i] = 	other.GetB().fC[i];
        	fC.fC[i] = 	other.GetC().fC[i];
        	fD.fC[i] = 	other.GetD().fC[i];
        	fN.fC[i] =  other.GetN().fC[i];

        }

        return *this;
    }

    void Update(void) override {

    	Update_ATCoeficients();
    	Update_NFactors();

    }


    __hydra_dual__ 
    inline double Evaluate( ArgTypeTime  t, ArgTypeThetah theta_h, ArgTypeThetal theta_l, ArgTypePhi phi)  const
    {
                                       
        auto aF = detail::AFactors(theta_h, theta_l, phi);

        double result = 0;

        #pragma unroll 10
        for(size_t i=0; i<10; ++i)
        	result += aF.fA[i]*fN.fC[i]*Time_Factor(i, t);


        return  result;

        
    }

    __hydra_dual__
	const detail::AngularTimeCoeficients& GetA() const {
		return fA;
	}

    __hydra_dual__
	const detail::AngularTimeCoeficients& GetB() const {
		return fB;
	}

    __hydra_dual__
	const detail::AngularTimeCoeficients& GetC() const {
		return fC;
	}

    __hydra_dual__
	const detail::AngularTimeCoeficients& GetD() const {
		return fD;
	}

    __hydra_dual__
	const detail::NFactors& GetN() const {
		return fN;
	}

private:


    void Update_ATCoeficients();

    void Update_NFactors()
    {
    	/*
    	0: A_0,
        1: A_perp,
    	2: A_S,
    	3: DeltaGamma_sd,
    	4: DeltaGamma,
    	5: DeltaM ,
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

    __hydra_host__ __hydra_device__
    inline double Time_Factor(int i, double x) const
    {

		/*
		0: A_0,
		1: A_perp,
		2: A_S,
		3: DeltaGamma_sd,
		4: DeltaGamma,
		5: DeltaM ,
		*/

    	double X1=0.5*x*_par[4];
    	double X2=x*_par[5];
    	static const double f = 0.238732414638; //3./(4*PI);

        return f * ::exp( -(_par[3] + 0.65789) * x) *
                                ( fA.fC[i]*::cosh(X1)+ fB.fC[i]*::sinh(X1)+
                                  ( fC.fC[i]*::cos(X2) + fD.fC[i]*::sin(X2) )*CP );
    }


    detail::NFactors fN;
    detail::AngularTimeCoeficients fA;
    detail::AngularTimeCoeficients fB;
    detail::AngularTimeCoeficients fC;
    detail::AngularTimeCoeficients fD;

};

}  // namespace medusa

#include <medusa/models/phi_s/detail/PhisSignalOnly.inl>

#endif /* PHISSIGNALONLY_H_ */
