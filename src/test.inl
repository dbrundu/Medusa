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
 *  test.inl
 *
 *  Created on: 07/05/2020
 *      Author: Davide Brundu
 */



#ifndef TEST_INL_
#define TEST_INL_

#include <iostream>
#include <medusa/models/phi_s/PhisTimeRes.h>
#include <medusa/models/phi_s/PhisAngularDist.h>
#include <medusa/models/phi_s/PhisTimeDist.h>
#include <medusa/models/phi_s/PhisN.h>
#include <medusa/models/D2hhmumu/D2hhmumuAngularDist.h>
#include <medusa/models/B2Kstarmumu/KstmumuAngularDist.h>


using namespace hydra::arguments;

declarg(theta_h_type,  double)
declarg(theta_l_type,  double)
declarg(phi_type,      double)
declarg(time_type,     double)
declarg(Errtime_type,  double)


int main(int argv, char** argc)
{

    double pars15[15] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double pars3[3]   = {0.0, 0.0, 0.0};
    double pars4[4]   = {0.0, 0.0, 0.0, 0.0};
    double pars8[8]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // test compilation
    medusa::PhisTimeRes<time_type, Errtime_type>   test_timeres(pars4);
    medusa::PhisAngularDist<1, theta_h_type, theta_l_type, phi_type>   test_angular;
    medusa::PhisTimeDist<1, false, time_type> test_time(pars15);
    medusa::PhisN<1> test_N(pars3);
    medusa::KstmumuAngularDist<medusa::PWave , theta_h_type, theta_l_type, phi_type> test_angularB1(pars8);
    medusa::KstmumuAngularDist<medusa::PSWave, theta_h_type, theta_l_type, phi_type> test_angularB2(pars15);
    medusa::D2hhmumuAngularDist<theta_l_type, phi_type> test_angularD(pars8);
    
    std::cout << test_angular(0.3, 0.2, 0.4) << std::endl;
    
    return 0;
}

#endif
