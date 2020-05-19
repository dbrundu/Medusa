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
 * 
 *
 *  Created on: 07/05/2020
 *      Author: Davide Brundu
 */



#ifndef TEST_INL_
#define TEST_INL_

#include <iostream>
#include <medusa/models/phi_s/PhisAngularDist.h>
#include <medusa/models/phi_s/PhisTimeDist.h>
#include <medusa/models/phi_s/PhisN.h>

int main(int argv, char** argc)
{

    // test compilation
    medusa::PhisAngularDist<1>   test_angular;
    medusa::PhisTimeDist<1, false, double> test_time(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    medusa::PhisN<1> test_N(0.0, 0.0, 0.0, 0.0);
    
    std::cout << test_angular(0.3, 0.2, 0.4) << std::endl;
    
    return 0;
}

#endif
