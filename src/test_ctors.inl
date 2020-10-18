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
 *  test_ctors.inl
 *
 *  Created on: 07/05/2020
 *      Author: Davide Brundu
 */



#ifndef TEST_CTORS_INL_
#define TEST_CTORS_INL_

#include <iostream>
#include <medusa/models/phi_s/PhisSignalOnly.h>
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

    double pars18[18] = {0};
    double pars15[15] = {0};
    double pars3[3]   = {0};
    double pars4[4]   = {0};
    double pars8[8]   = {0};

    // test compilation of ctors

    medusa::PhisSignalOnly< 1, time_type, theta_h_type, theta_l_type, phi_type> test_phis(pars18);
    
    medusa::KstmumuAngularDist<medusa::PWave , theta_h_type, theta_l_type, phi_type> test_angularB1(pars8);

    medusa::KstmumuAngularDist<medusa::PSWave, theta_h_type, theta_l_type, phi_type> test_angularB2(pars15);

    medusa::D2hhmumuAngularDist<theta_l_type, phi_type> test_angularD(pars8);
    

    return 0;
}

#endif
