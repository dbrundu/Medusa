/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Davide Brundu, Antonio Augusto Alves Junior,
 *                      Piera Muzzetto et al.
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
 *  D2hhmumuAngularDist.h
 *
 *  Created on: 20/05/2020
 *      Author: Davide Brundu, Andrea Contu
 */

#ifndef D2HHMUMUANGULARDIST_H_
#define D2HHMUMUANGULARDIST_H_

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




namespace medusa {


/*
 *  @class D2hhmumuAngularDist
 *  Declaration of the B->K*mumu angular distribution class,
 *  the two parametrizations are for P-wave contribution only
 *  or both the S- and P-wave
 */
template<typename ArgType1, typename ArgType2, typename Signature=double(ArgType1, ArgType2)>
class D2hhmumuAngularDist : public hydra::BaseFunctor<D2hhmumuAngularDist<ArgType1, ArgType2>, Signature, 8>
{

    using ThisBaseFunctor = hydra::BaseFunctor<D2hhmumuAngularDist<ArgType1, ArgType2>, Signature, 8>;
    using ThisBaseFunctor::_par;


public:

    D2hhmumuAngularDist() = delete;

    D2hhmumuAngularDist(hydra::Parameter const& H2, hydra::Parameter const& H3,
                        hydra::Parameter const& H4, hydra::Parameter const& H5,
                        hydra::Parameter const& H6, hydra::Parameter const& H7,
                        hydra::Parameter const& H8, hydra::Parameter const& H9):
    ThisBaseFunctor{H2, H3, H4, H5, H6, H7, H8, H9}
    {}
    
    
    // ctor with array of hydra::Parameter
    // is needs to  respect the parameter order as the main ctor
    explicit D2hhmumuAngularDist( const hydra::Parameter (&Hs)[8] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3], Hs[4], Hs[5], Hs[6], Hs[7] }
    {}
    
    
    // ctor with array of double
    // is needs to  respect the parameter order as the main ctor
    explicit D2hhmumuAngularDist( const double (&Hs)[8] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3], Hs[4], Hs[5], Hs[6], Hs[7] }
    {}


    __hydra_dual__
    D2hhmumuAngularDist(const D2hhmumuAngularDist& other):
    ThisBaseFunctor(other)
    {}

    __hydra_dual__
    D2hhmumuAngularDist& operator=( D2hhmumuAngularDist const& other)
    {
        if(this==&other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ 
    inline double Evaluate(ArgType1 const& theta_l, ArgType2 const& phi)  const
    {
        const double  integral = 1.0;
        const double  I1fixed  = 0.5*integral + _par[0]/3.0;


        return 1.0/(2*PI) * ( C_1(theta_l, phi)*I1fixed + C_2(theta_l, phi)*_par[0] + C_3(theta_l, phi)*_par[1] + \
               C_4(theta_l, phi)*_par[2] + C_5(theta_l, phi)*_par[3] + C_6(theta_l, phi)*_par[4] + \
               C_7(theta_l, phi)*_par[5] + C_8(theta_l, phi)*_par[6] + C_9(theta_l, phi)*_par[7] );
    }


private:

    __hydra_dual__ inline
    double C_1(const double& theta_l, const double& chi) const {
                return 1.;
    }

    __hydra_dual__ inline
    double C_2(const double& theta_l, const double& chi) const {
                return ::cos(2*theta_l);
    }

    __hydra_dual__ inline
    double C_3(const double& theta_l, const double& chi) const {
                return ::pow(::sin(theta_l),2) * ::cos(2*chi);
    }

    __hydra_dual__ inline
    double C_4(const double& theta_l, const double& chi) const {
                return ::sin(2*theta_l)*::cos(chi);
    }

    __hydra_dual__ inline
    double C_5(const double& theta_l, const double& chi) const {
                return ::sin(theta_l)*::cos(chi);
    }

    __hydra_dual__ inline
    double C_6(const double& theta_l, const double& chi) const {
                return ::cos(theta_l);
    }

    __hydra_dual__ inline
    double C_7(const double& theta_l, const double& chi) const {
                return ::sin(theta_l)*::sin(chi);
    }

    __hydra_dual__ inline
    double C_8(const double& theta_l, const double& chi) const {
                return ::sin(2*theta_l)*::sin(chi);
    }

    __hydra_dual__ inline
    double C_9(const double& theta_l, const double& chi) const {
                return ::pow(::sin(theta_l),2)*::sin(2*chi);
    }
};




} // namespace medusa


// namespace hydra {
// 
//     template<>
//     struct IntegrationFormula< medusa::D2hhmumuAngularDist, 2>
//     {
// 
//     	inline std::pair<hydra::GReal_t, hydra::GReal_t>
//     	EvalFormula( medusa::D2hhmumuAngularDist const& , const double* , const double*  ) const
//     	{
//     		const double r =  1.0;
//     		return std::make_pair( r , 0.0);
//     	}
// 
//     };

// } // namespace hydra

#endif /* D2HHMUMUANGULARDIST_H_ */
