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
#include <eos/rare-b-decays/lambdab-to-lambda-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class LambdabToLambdaCharmoniumBRvD2021 :
    public TestCase
{
    public:
    LambdabToLambdaCharmoniumBRvD2021() :
            TestCase("lambdab_to_lambda_charmonium_BRvD2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::Lambda"]                            = 1.115683;
            p["mass::Lambda_b"]                          = 5.61960;
            p["mass::D^0"]                               = 1.86483;
            p["b->sccbar::t_0"]                          = 10.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
            p["life_time::Lambda_b"]                     = 1.471e-12;
            p["decay-constant::J/psi"]                   = 0.2773;
            p["hbar"]                                    = 6.5821e-25;
            p["Lambda::alpha"]                           = 0.750;

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

            Options oo;
            oo.set("model",          "SM");
            oo.set("formfactor",     "BRvD2021");
            oo.set("psi",            "J/psi");

            LambdabToLambdaCharmonium c(p, oo);

            //===============Angular-Observable===================//

            TEST_CHECK_NEARLY_EQUAL(c.branching_ratio(), 16858.471949, eps);
            TEST_CHECK_NEARLY_EQUAL(c.K1ss(),  0.286214, eps);
            TEST_CHECK_NEARLY_EQUAL(c.K1cc(),  0.427572, eps);
            TEST_CHECK_NEARLY_EQUAL(c.K2ss(), -0.198159, eps);
            TEST_CHECK_NEARLY_EQUAL(c.K2cc(), -0.319556, eps);
            TEST_CHECK_NEARLY_EQUAL(c.K3sc(),       0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(c.K4sc(), -0.0642176, eps);

            //===============Parameters===================//

            TEST_CHECK_NEARLY_EQUAL(c.abs_aplus(),  1.16454, eps);
            TEST_CHECK_NEARLY_EQUAL(c.abs_aminus(), 2.80841, eps);
            TEST_CHECK_NEARLY_EQUAL(c.abs_bplus(),  7.38048, eps);
            TEST_CHECK_NEARLY_EQUAL(c.abs_bminus(), 0.309115, eps);

            TEST_CHECK_NEARLY_EQUAL(c.arg_aplus(),  0.785398, eps);
            TEST_CHECK_NEARLY_EQUAL(c.arg_aminus(), -2.35619, eps);
            TEST_CHECK_NEARLY_EQUAL(c.arg_bplus(),  -2.35619, eps);
            TEST_CHECK_NEARLY_EQUAL(c.arg_bminus(), 0.785398, eps);

            TEST_CHECK_NEARLY_EQUAL(c.alpha_b(), 0.749799, eps);

        }
} lambdab_to_lambda_charmonium_BRvD2021_test;
