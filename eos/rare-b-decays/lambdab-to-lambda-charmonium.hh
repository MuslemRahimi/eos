/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Muslem Rahimi
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_CHARMONIUM_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_CHARMONIUM_HH 1

#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /**
     * Decay: Lambda_b -> Lambda psi
     * with psi a narrow charmonium: psi = { J/psi, psi(2S) }
     */
    class LambdabToLambdaCharmonium :
        public ParameterUser,
        public PrivateImplementationPattern<LambdabToLambdaCharmonium>
    {
        public:
            ///@name Basic operations.
            ///@{

            LambdabToLambdaCharmonium(const Parameters & parameters, const Options & options);
            ~LambdabToLambdaCharmonium();

            ///@}

            ///@name Observables
            ///@{

            /// Branching ratio
            double branching_ratio() const;

            /// Angular observables as detected in the decay Lambdab -> Lambda psi (-> l^+ l^-)
            double K1ss() const;
            double K1cc() const;
            double K2ss() const;
            double K3sc() const;
            double K4sc() const;

            ///@}
    };
}

#endif
