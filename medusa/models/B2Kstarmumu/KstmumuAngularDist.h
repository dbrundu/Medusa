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
 *  KstmumuAngularDist.h
 *
 *  Created on: 20/05/2020
 *      Author: Davide Brundu
 */


#ifndef KSTMUMUANGULARDIST_H_
#define KSTMUMUANGULARDIST_H_

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

#include <medusa/models/B2Kstarmumu/detail/kstarmumu_angular_functions.h>


namespace medusa {


enum KstmumuWave { PWave=0, PSWave };


/*
 *  @class KstmumuAngularDist
 *  Declaration of the B->K*mumu angular distribution class,
 *  the two parametrizations are for P-wave contribution only
 *  or both the S- and P-wave
 */
template< KstmumuWave L, typename ArgType1, typename ArgType2, typename ArgType3>
class KstmumuAngularDist{};




/*
 *  @class KstmumuAngularDist
 *  Specialization for P-wave contribution only
 *
 */
template< typename ArgType1, typename ArgType2, typename ArgType3>
class KstmumuAngularDist<PWave, ArgType1, ArgType2, ArgType3> : 
public hydra::BaseFunctor< KstmumuAngularDist<PWave, ArgType1, ArgType2, ArgType3>, double(ArgType1, ArgType2, ArgType3), 8>
{

    using ThisBaseFunctor = hydra::BaseFunctor< KstmumuAngularDist<PWave, ArgType1, ArgType2, ArgType3>, double(ArgType1, ArgType2, ArgType3), 8>;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;

public:

    KstmumuAngularDist()=delete;
    
    
    KstmumuAngularDist(Param const& FL,  Param const& S3,
                       Param const& S4,  Param const& S5,
                       Param const& AFB, Param const& S7,
                       Param const& S8,  Param const& S9):
    ThisBaseFunctor{FL, S3, S4, S5, AFB, S7, S8, S9}
    {}


    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit KstmumuAngularDist( const hydra::Parameter (&Hs)[8] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3], Hs[4], Hs[5], Hs[6], Hs[7] }
    {}
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit KstmumuAngularDist( const double (&Hs)[8] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3], Hs[4], Hs[5], Hs[6], Hs[7] }
    {}


    __hydra_dual__
    KstmumuAngularDist( KstmumuAngularDist<PWave, ArgType1, ArgType2, ArgType3> const& other):
    ThisBaseFunctor(other)
    {}


    __hydra_dual__
    KstmumuAngularDist& operator=( KstmumuAngularDist<PWave, ArgType1, ArgType2, ArgType3> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }


    __hydra_dual__ inline
    double Evaluate(ArgType1 const& theta_h, ArgType2 const& theta_l, ArgType3 const& phi)  const  {

        auto r = (3./4.)*(1-_par[0]) * detail::kstmumu_angular_fun<0>(theta_h, theta_l, phi)   +\
                             _par[0] * detail::kstmumu_angular_fun<1>(theta_h, theta_l, phi)   +\
                 (1./4.)*(1-_par[0]) * detail::kstmumu_angular_fun<2>(theta_h, theta_l, phi)   +\
                            -_par[0] * detail::kstmumu_angular_fun<3>(theta_h, theta_l, phi)   +\
                             _par[1] * detail::kstmumu_angular_fun<4>(theta_h, theta_l, phi)   +\
                             _par[2] * detail::kstmumu_angular_fun<5>(theta_h, theta_l, phi)   +\
                             _par[3] * detail::kstmumu_angular_fun<6>(theta_h, theta_l, phi)   +\
                       4./3.*_par[4] * detail::kstmumu_angular_fun<7>(theta_h, theta_l, phi)   +\
                             _par[5] * detail::kstmumu_angular_fun<8>(theta_h, theta_l, phi)   +\
                             _par[6] * detail::kstmumu_angular_fun<9>(theta_h, theta_l, phi)   +\
                             _par[7] * detail::kstmumu_angular_fun<10>(theta_h, theta_l, phi);
                             
        return r * 9./(32*PI) ;

    }
    
};






/*
 *  @class KstmumuAngularDist
 *  Specialization for S- plus P-wave contributions
 *
 */
template< typename ArgType1, typename ArgType2, typename ArgType3>
class KstmumuAngularDist<PSWave, ArgType1, ArgType2, ArgType3>: 
public hydra::BaseFunctor< KstmumuAngularDist<PSWave, ArgType1, ArgType2, ArgType3>, double(ArgType1, ArgType2, ArgType3), 15>
{

    using ThisBaseFunctor = hydra::BaseFunctor< KstmumuAngularDist<PSWave, ArgType1, ArgType2, ArgType3>, double(ArgType1, ArgType2, ArgType3), 15>;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;

public:

    KstmumuAngularDist()=delete;
    
    
    KstmumuAngularDist(Param const& FL,  Param const& S3,
                       Param const& S4,  Param const& S5,
                       Param const& AFB, Param const& S7,
                       Param const& S8,  Param const& S9,
                       Param const& FS,  Param const& S11,
                       Param const& S13, Param const& S14,
                       Param const& S15, Param const& S16,
                       Param const& S17):
    ThisBaseFunctor{FL, S3, S4, S5, AFB, S7, S8, S9, FS, S11, S13, S14, S15, S16, S17}
    {}
    
    
    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit KstmumuAngularDist( const hydra::Parameter (&Hs)[15] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14] }
    {}
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit KstmumuAngularDist( const double (&Hs)[15] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14] }
    {}



    __hydra_dual__
    KstmumuAngularDist( KstmumuAngularDist<PSWave, ArgType1, ArgType2, ArgType3> const& other):
    ThisBaseFunctor(other)
    {}


    __hydra_dual__
    KstmumuAngularDist& operator=( KstmumuAngularDist<PSWave, ArgType1, ArgType2, ArgType3> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }


    __hydra_dual__ inline
    double Evaluate(ArgType1 const& theta_h, ArgType2 const& theta_l, ArgType3 const& phi)  const  {

        auto r = (3./4.)*(1-_par[0]) * detail::kstmumu_angular_fun<0>(theta_h, theta_l, phi)   +\
                             _par[0] * detail::kstmumu_angular_fun<1>(theta_h, theta_l, phi)   +\
                 (1./4.)*(1-_par[0]) * detail::kstmumu_angular_fun<2>(theta_h, theta_l, phi)   +\
                            -_par[0] * detail::kstmumu_angular_fun<3>(theta_h, theta_l, phi)   +\
                             _par[1] * detail::kstmumu_angular_fun<4>(theta_h, theta_l, phi)   +\
                             _par[2] * detail::kstmumu_angular_fun<5>(theta_h, theta_l, phi)   +\
                             _par[3] * detail::kstmumu_angular_fun<6>(theta_h, theta_l, phi)   +\
                       4./3.*_par[4] * detail::kstmumu_angular_fun<7>(theta_h, theta_l, phi)   +\
                             _par[5] * detail::kstmumu_angular_fun<8>(theta_h, theta_l, phi)   +\
                             _par[6] * detail::kstmumu_angular_fun<9>(theta_h, theta_l, phi)   +\
                             _par[7] * detail::kstmumu_angular_fun<10>(theta_h, theta_l, phi);
        r *=  9./(32*PI);
                             
        auto res = (1 - _par[8])*r + (3./(16.0*PI)) * _par[8] * detail::kstmumu_angular_fun<11>(theta_h, theta_l, phi) +\
                   (9./(32.0*PI)) * (_par[9] + _par[10]*::cos(2*theta_l)) * detail::kstmumu_angular_fun<12>(theta_h, theta_l, phi) +\
                   (9./(32.0*PI)) * (_par[11]*::sin(2*theta_l) + _par[12]*::sin(theta_l)) * detail::kstmumu_angular_fun<13>(theta_h, theta_l, phi) +\
                   (9./(32.0*PI)) * (_par[13]*::sin(theta_l) + _par[14]*::sin(2*theta_l)) * detail::kstmumu_angular_fun<14>(theta_h, theta_l, phi);
                   
        return res;


    }
    
};






}  // namespace medusa


#endif /* KSTMUMUANGULARDIST_H_ */
