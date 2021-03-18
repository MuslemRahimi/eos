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
            static const double eps = 1e-4;

            {
                Parameters p = Parameters::Defaults();
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::Lambda"]                            = 1.115683;
                p["mass::Lambda_b"]                          = 5.61960;
                p["mass::D^0"]                               = 1.86483;
                p["b->sccbar::t_0"]                          = 9.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;

                p["Lambda_b->Lambdaccbar::Re{alpha_0^V_long}@BRvD2021"]  = 1.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_0^V_long}@BRvD2021"]  = 1.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_1^V_long}@BRvD2021"]  = 2.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_1^V_long}@BRvD2021"]  = 2.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_2^V_long}@BRvD2021"]  = 3.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_2^V_long}@BRvD2021"]  = 3.0;

                p["Lambda_b->Lambdaccbar::Re{alpha_0^V_perp}@BRvD2021"]  = 4.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_0^V_perp}@BRvD2021"]  = 4.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_1^V_perp}@BRvD2021"]  = 5.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_1^V_perp}@BRvD2021"]  = 5.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_2^V_perp}@BRvD2021"]  = 6.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_2^V_perp}@BRvD2021"]  = 6.0;

                p["Lambda_b->Lambdaccbar::Re{alpha_0^A_long}@BRvD2021"]  = 7.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_0^A_long}@BRvD2021"]  = 7.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_1^A_long}@BRvD2021"]  = 8.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_1^A_long}@BRvD2021"]  = 8.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_2^A_long}@BRvD2021"]  = 9.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_2^A_long}@BRvD2021"]  = 9.0;

                p["Lambda_b->Lambdaccbar::Re{alpha_0^A_perp}@BRvD2021"]  = 10.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_0^A_perp}@BRvD2021"]  = 10.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_1^A_perp}@BRvD2021"]  = 11.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_1^A_perp}@BRvD2021"]  = 11.0;
                p["Lambda_b->Lambdaccbar::Re{alpha_2^A_perp}@BRvD2021"]  = 12.0;
                p["Lambda_b->Lambdaccbar::Im{alpha_2^A_perp}@BRvD2021"]  = 12.0;

                Options o = { { "model", "SM" } };

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

                    std::make_pair( 0.752572, eps),         // Re{alpha_LbL}
                    std::make_pair( 0.0,  eps),             // Im{alpha_Lbl}

                    std::make_pair( -0.893623, eps),         // Re{blaschke_factor(q2 = 16)}
                    std::make_pair( 0.448818,  eps),         // Im{blaschke_factor(q2 = 16)}

                    std::make_pair( -0.402961, eps),        // real(z(q2 = 10.0))
                    std::make_pair( 0.915217, eps),         // img(z(q2 = 10.0))

                    //===================Testcase for outer function======================//
                    std::make_pair( 7.05783, eps),         // Re{phi_V_long(q2 = 16.0)}
                    std::make_pair( 33.5449,  eps),        // Im{phi_V_long(q2 = 16.0)}

                    std::make_pair( 4.19156, eps),        // Re{phi_V_perp(q2 = 16.0)}
                    std::make_pair( 19.9219,  eps),        // Im{phi_V_perp(q2 = 16.0)}

                    std::make_pair( 12.3545, eps),        // Re{phi_A_long(q2 = 16.0)}
                    std::make_pair( 58.7191,  eps),        // Im{phi_A_long(q2 = 16.0)}

                    std::make_pair( 10.9722, eps),        // Re{phi_A_perp(q2 = 16.0)}
                    std::make_pair( 52.1494,  eps),        // Im{phi_A_perp(q2 = 16.0)}

                    //===================Testcase for nonlocal FF======================//
                    std::make_pair( -0.275274,  eps),        // Re{H_V_long(q2 = 16.0)}
                    std::make_pair( 1.37506,  eps),        // Im{H_V_long(q2 = 16.0)}

                    std::make_pair( -0.877161,  eps),        // Re{H_V_long(q2 = 16.0)}
                    std::make_pair(  4.50042,  eps),        // Im{H_V_long(q2 = 16.0)}

                    std::make_pair( -0.437938,  eps),        // Re{H_A_long(q2 = 16.0)}
                    std::make_pair(  2.26822,  eps),        // Im{H_A_long(q2 = 16.0)}

                    std::make_pair( -0.651129,  eps),        // Re{H_A_long(q2 = 16.0)}
                    std::make_pair( 3.3887,  eps),          // Im{H_A_long(q2 = 16.0)}

                    //============Testcase for residue of nonlocal FF================//
                    std::make_pair( -0.696651,  eps),           // Re{H_V_long_residue_jpsi()}
                    std::make_pair( -0.696651,  eps),          // Im{H_V_long_residue_jpsi()}

                    std::make_pair( -3.00529,  eps),           // Re{H_V_perp_residue_jpsi()}
                    std::make_pair( -3.00529,  eps),          // Im{H_V_perp_residue_jpsi()}

                    std::make_pair( -1.69011,  eps),           // Re{H_A_long_residue_jpsi()}
                    std::make_pair( -1.69011,  eps),          // Im{H_A_long_residue_jpsi()}

                    std::make_pair( -3.27275,  eps),           // Re{H_A_perp_residue_jpsi()}
                    std::make_pair( -3.27275,  eps),          // Im{H_A_perp_residue_jpsi()}


                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

            }
        }
} nonlocal_formfactor_BRvD2021_test;