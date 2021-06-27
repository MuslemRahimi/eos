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
#include <eos/form-factors/parametric-brvd2021.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BRvD2021FormFactorsTest :
    public TestCase
{
    public:
        BRvD2021FormFactorsTest() :
            TestCase("BRvD2021_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-4;

            Parameters p = Parameters::Defaults();

            /* Outer function */
            {
                BRvD2021FormFactors ff(p, Options{ });
                //const double _t_p(LambdaBToLambda::tau_p);
                //const double _t_m(LambdaBToLambda::tau_m);
                const double _chi_1m( 5.131e-04); // fix values from Bhrarucha/Feldmann/Wick '10
                const double _chi_0p( 6.204e-03);
                const double _chi_1p( 3.894e-04);
                const double _chi_0m(19.421e-03);

                //===================Testcase for outer function======================//
                const unsigned phiParam_V_t[7] = {0, 0, 0, 1, 2, 0, 1};  // V_t polarization
                TEST_CHECK_NEARLY_EQUAL( 0.629133,  ff._phi(16, _chi_0p, phiParam_V_t) , eps);

                const unsigned phiParam_V_long[7] = {1, 0, 1, 0, 2, 1, 0};  // V_long polarization
                TEST_CHECK_NEARLY_EQUAL( 0.721545,  ff._phi(16, _chi_1m, phiParam_V_long) , eps);

                const unsigned phiParam_V_perp[7] = {0, 1, 0, 0, 1, 1, 0};  // V_perp polarization
                TEST_CHECK_NEARLY_EQUAL( 0.606014,  ff._phi(16, _chi_1m, phiParam_V_perp) , eps);

            }
            /* Local FF */
            {
                BRvD2021FormFactors ff(p, Options{ });
                const double _t_p(LambdaBToLambda::tau_p);
                const double _t_m(LambdaBToLambda::tau_m);
                const double _chi_1m( 5.131e-04); // fix values from Bhrarucha/Feldmann/Wick '10
                const double _chi_0p( 6.204e-03);
                const double _chi_1p( 3.894e-04);
                const double _chi_0m(19.421e-03);

                p["Lambda_b->Lambda::a^V,t_0@BRvD2021"]  = 1.0;
                p["Lambda_b->Lambda::a^V,t_1@BRvD2021"]  = 2.0;
                p["Lambda_b->Lambda::a^V,t_2@BRvD2021"]  = 3.0;

                p["Lambda_b->Lambda::a^V,long_0@BRvD2021"]  = 4.0;
                p["Lambda_b->Lambda::a^V,long_1@BRvD2021"]  = 5.0;
                p["Lambda_b->Lambda::a^V,long_2@BRvD2021"]  = 6.0;

                p["Lambda_b->Lambda::a^V,perp_0@BRvD2021"]  = 7.0;
                p["Lambda_b->Lambda::a^V,perp_1@BRvD2021"]  = 8.0;
                p["Lambda_b->Lambda::a^V,perp_2@BRvD2021"]  = 9.0;
                //===================Testcase for local ff======================//
                TEST_CHECK_NEARLY_EQUAL( 0.351287, ff.f_time_v(2.0) , eps);
                TEST_CHECK_NEARLY_EQUAL( 0.468319, ff.f_long_v(2.0) , eps);
                TEST_CHECK_NEARLY_EQUAL( 2.54744, ff.f_perp_v(2.0) , eps);
            }


        }
} BRvD2021_form_factor_test;
