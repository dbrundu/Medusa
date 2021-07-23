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

/*----------------------------------------
 *  Generate_dataset.h
 *
 *  Created: 07/07/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *----------------------------------------*/

#ifndef GENERATE_DATASET_H_
#define GENERATE_DATASET_H_


// std
#include <iostream>
#include <chrono>

// Hydra
#include <hydra/device/System.h>
#include <hydra/detail/Config.h>
#include <hydra/Placeholders.h>
#include <hydra/multivector.h>
#include <hydra/PhaseSpace.h>
#include <hydra/Vector4R.h>
#include <hydra/SeedRNG.h>
#include <hydra/Lambda.h>
#include <hydra/Random.h>
#include <hydra/Tuple.h>

// Medusa
#include <medusa/phi_s/details/Parameters.h>
#include <medusa/Angles.h>



namespace medusa {


    /* Function that provides the Monte Carlo dataset (time and helicity angles) for the decay
     *                          B0s -> J/psi (Phi -> K+ K-)
     *                                  |-> mu+ mu-
     * Parameters:
     * model            = PDF according to which to generate the Monte Carlo dataset
     * final_dataset    = container where to save the Monte Carlo dataset
     * nevents          = number of events, namely dimension of final_dataset
     * bunch_size       = number of events generated from Hydra's Monte Carlo per filling cycle
     */

    template<typename Model, typename Container>
    size_t generate_dataset(Model const& model, Container& final_dataset, size_t nevents, size_t bunch_size)
    {

        // default namespaces
        using namespace hydra::placeholders;


        // constants
        const double B0s_mass  = 5.336688;
        const double Jpsi_mass = 3.0969;
        const double Phi_mass  = 1.019462;
        const double K_mass    = 0.493677;
        const double mu_mass   = 0.1056583745;


        // function that takes the four-vectors and returns the decay times and angles in the helicity basis
        auto CastToVariables  = hydra::wrap_lambda(
                [] __hydra_dual__ (Jpsi jpsi, Phi phi, MuonP mup, MuonM mum, KaonP kaonp, KaonM kaonm, size_t n )
        {
            hydra_thrust::default_random_engine engine;
            hydra_thrust::uniform_real_distribution<double> uniDist(0.0, 20.0);
            engine.discard(n);
            dtime_t decay_time = uniDist(engine);

            theta_l_t theta_l    = ::acos( medusa::cos_decay_angle(jpsi + phi, phi,  kaonp) );
            theta_h_t theta_h    = ::acos( medusa::cos_decay_angle(jpsi + phi, jpsi, mup) );
            phi_t     phiangle   = medusa::phi_plane_angle(kaonm, kaonp, mup, mum);

            return hydra::make_tuple(decay_time, theta_h, theta_l, phiangle) ;

        });


        //------------------------------------------
        //     Space-phase Monte Carlo
        //------------------------------------------

        // mother particle
        hydra::Vector4R B0s(B0s_mass, 0.0, 0.0, 0.0);


        // phase-space generator for B0s -> Jpsi Phi
        hydra::PhaseSpace<2> B0s_phsp(B0s_mass, {Jpsi_mass, Phi_mass});

        // phase-space generator for Jpsi -> mu+ mu-
        hydra::PhaseSpace<2> Jpsi_phsp(Jpsi_mass, {mu_mass , mu_mass});

        // phase-space generator for Phi -> K+ K-
        hydra::PhaseSpace<2> Phi_phsp(Phi_mass, {K_mass , K_mass});


        // container to hold the final state particles for B0s -> Jpsi Phi
        auto B0sDecay = hydra::Decays< hydra::tuple<Jpsi,Phi>,
                              hydra::device::sys_t >( B0s_mass, {Jpsi_mass, Phi_mass }, bunch_size);

        // container to hold the final state particles for Jpsi -> mu+ mu-
        auto JpsiDecay = hydra::Decays< hydra::tuple<MuonP,MuonM>,
                              hydra::device::sys_t >(Jpsi_mass, { mu_mass , mu_mass}, bunch_size);

        // container to hold the final state particles for Phi -> K+ K-
        auto PhiDecay = hydra::Decays< hydra::tuple<KaonP,KaonM>,
                              hydra::device::sys_t >(Phi_mass, { K_mass , K_mass}, bunch_size);


        // container to hold the times and helicity angles of the decay
        hydra::multivector<hydra::tuple<dtime_t, theta_h_t, theta_l_t, phi_t> , hydra::device::sys_t> dataset(bunch_size);


        // random seed generator inizialized to the default
        hydra::SeedRNG seeder{};

        // counter
        size_t k = 0;


        // start to profile the execution time
        auto start = std::chrono::high_resolution_clock::now();
        

        // fill final_dataset until nevents following the steps:
        // 1. generate an unmodeled dataset;
        // 2. erase events from the unmodeled dataset to model it in according to the given model;
        // 3. insert the modeled dataset in final_dataset.
        do {
            // generate the the final state particles for B0s -> Jpsi Phi
            B0s_phsp.SetSeed( seeder() );
            B0s_phsp.Generate(B0s, B0sDecay);

            // pass the list of J/psi to generate the final state particles for J/psi -> mu+ mu-
            Jpsi_phsp.SetSeed( seeder() );
            Jpsi_phsp.Generate(B0sDecay.GetDaugtherRange(_0), JpsiDecay);

            //pass the list of Phi to generate the final state particles for Phi -> K+ K-
            Phi_phsp.SetSeed( seeder() );
            Phi_phsp.Generate(B0sDecay.GetDaugtherRange(_1), PhiDecay);

            // merge the containers B0sDecay, JpsiDecay, PhiDecay and pass the global container
            // to CastToVariables to obtain time and angles of the decays
            auto range = B0sDecay.Meld( JpsiDecay, PhiDecay, hydra::range(0, B0sDecay.size()) ) | CastToVariables;

            hydra::copy(range, dataset);

            // erase events in the dataset to model it according to the given model
            auto dataset_unwgt = hydra::unweight( dataset, model);

            // insert dataset_unwgt in final_dataset
            final_dataset.insert(final_dataset.end(), dataset_unwgt.begin(), dataset_unwgt.end() );

            ++k;

        } while(final_dataset.size()<nevents );


        // end to profile the execution time
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;


        // erase final_dataset after nevents. This command is necessary, because we do not know
        // the number of events that have been erased by unweight() function
        final_dataset.erase(final_dataset.begin()+nevents, final_dataset.end());

        // output
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "-----------Generation on Device ---------"  << std::endl;
        std::cout << "| B0s -> J/psi Phi -> mu+ mu- K+ K-      "  << std::endl;
        std::cout << "| Number of events: "<< nevents             << std::endl;
        std::cout << "| Number of events gen: "<< k*bunch_size    << std::endl;
        std::cout << "| Time (ms):        "<< elapsed.count()     << std::endl;
        std::cout << "-----------------------------------------"  << std::endl;

        return final_dataset.size();


    } // generate_dataset()

} // namespace medusa

#endif // GENERATE_DATASET_H_