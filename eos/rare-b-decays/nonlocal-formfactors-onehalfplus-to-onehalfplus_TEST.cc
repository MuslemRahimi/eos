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

#include <test/test.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class NonlocalFormFactorBRvD2021Test :
    public TestCase
{
    public:
        NonlocalFormFactorBRvD2021Test() :
            TestCase("nonlocal_formfactor_BRvD2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 9.0;
                p["b->sccbar::t_s"]                          = -17.4724;

                Options o = { { "model", "WilsonScan" } };

                auto nc = NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::make("Lambda_b->Lambda::BRvD2021", p, o);


                auto diagnostics = nc->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {

                    std::make_pair( 0.0, eps),            // z(q2=10.0)
                    std::make_pair( 0.0, eps),            // z(q2=10.0)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

            }
        }
} nonlocal_formfactor_BRvD2021_test;