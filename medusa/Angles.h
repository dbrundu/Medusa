/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu,
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto,
 *                      Alessandro Maria Ricci
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

/*---------------------------------------------------------------
 *  Created: 26/06/2020
 *
 *  Authors: Davide Brundu
 *  
 *  Note: updated by Alessandro Maria Ricci in 28/05/2021
 *
 *  Functions to compute the decay angles theta and phi
 *  in the helicity basis.
 * 
 * The formulas used here can be found also in EvtGen 
 * [EvtGen/src/EvtKine.cpp] (https://evtgen.hepforge.org/).
 *---------------------------------------------------------------*/


#ifndef MEDUSA_ANGLES_H_
#define MEDUSA_ANGLES_H_


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
 *  This function returns the cosine of the decay angle theta.
 *  The decay angle calculated is that between
 *  the flight direction of the daughter meson, "D",
 *  in the rest frame of "Q" (the parent of "D"),
 *  with respect to "Q"'s flight direction in "P"'s
 *  (the parent of "Q") rest frame
 *  P == B0, Q = dimuon, D = muon
 *
 */

__hydra_dual__
inline double cos_decay_angle(hydra::Vector4R const& p, hydra::Vector4R const& q, hydra::Vector4R const& d){
  // P == B0, Q = Jpsi, D = muon
  double pd = p*d;
  double pq = p*q;
  double qd = q*d;
  double mp2 = p.mass2();
  double mq2 = q.mass2();
  double md2 = d.mass2();

  return (pd * mq2 - pq * qd) / ::sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));
}




/** Evaluate the angle phi beween two decay planes,
 *  formed by particles d2&d3 and h1&h2 correspondingly.
 *  The angle is evaluated in the rest frame
 *  of "mother" particles (defined as d2+d3+h1+h2)
 *  It is calculated as the angle formed by the h1 3vector projection
 *  on an x-y plane defined by d2(=x), h1+h2 (=z)
 *  For LHCb convention with B0->h+h-mu+mu- ==> d2 = h-, d3=h+, h1 = mu+, h2=mu-
 */

__hydra_dual__
inline double phi_plane_angle(hydra::Vector4R d2, hydra::Vector4R d3, hydra::Vector4R h1, hydra::Vector4R h2) {

  hydra::Vector4R Mother = d2 + d3 + h1 + h2;
  d2.applyBoostTo(Mother, /*inverse boost? == */ true);
  d3.applyBoostTo(Mother, /*inverse boost? == */ true);
  h1.applyBoostTo(Mother, /*inverse boost? == */ true);
  h2.applyBoostTo(Mother, /*inverse boost? == */ true);

  hydra::Vector4R D = d2 + d3;

  hydra::Vector4R d1_perp = d2 - (D.dot(d2) / D.dot(D)) * D; // d2 will be mu^+
  hydra::Vector4R h1_perp = h1 - (D.dot(h1) / D.dot(D)) * D;

  // orthogonal to both D and d1_perp
  // hydra::Vector4R d1_prime = D.cross(d1_perp);
  hydra::Vector4R d1_prime = d1_perp.cross(D);


  d1_perp  = d1_perp / d1_perp.d3mag();
  d1_prime = d1_prime / d1_prime.d3mag();

  double cos_phi, sin_phi;

  cos_phi = d1_perp.dot(h1_perp);   //cos_chi
  sin_phi = d1_prime.dot(h1_perp);  //sin_chi

  double phi = ::atan2(sin_phi, cos_phi);
  
  return (phi>=0)? phi : phi + 2*PI;
}




}  // namespace medusa


#endif /* MEDUSA_ANGLES_H_ */
