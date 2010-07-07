/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/charm-loops.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace wf;

class OneLoopTest :
    public TestCase
{
    public:
        OneLoopTest() :
            TestCase("one_loop_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christoph Bobeth's result from May 2010 */

            /* One-Loop */
            {
                static const double mu = 4.2, s = 1.0, m_c = 1.4, m_b = 4.8, eps = 0.00001;
                TEST_CHECK_NEARLY_EQUAL(+1.57192, real(CharmLoops::h(mu, s)), eps);
                TEST_CHECK_NEARLY_EQUAL(+1.39626, imag(CharmLoops::h(mu, s)), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.58013, CharmLoops::h(mu, s, m_c), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.55926, CharmLoops::h(mu, s, m_b), eps);
            }
        }
} one_loop_test;

class HelperTest :
    public TestCase
{
    public:
        HelperTest() :
            TestCase("helper_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Mathematica results from July 2010 */

            /* C0 */
            {
                static const double m_b = 4.45, eps = 0.000001;
                // real parts
                TEST_CHECK_NEARLY_EQUAL(-1.64493406685, real(CharmLoops::C0( 0.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.648607,      real(CharmLoops::C0( 0.5,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.652304,      real(CharmLoops::C0( 1.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.659779,      real(CharmLoops::C0( 2.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.667360,      real(CharmLoops::C0( 3.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.690774,      real(CharmLoops::C0( 6.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.715257,      real(CharmLoops::C0( 9.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.740899,      real(CharmLoops::C0(12.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.767803,      real(CharmLoops::C0(15.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.796088,      real(CharmLoops::C0(18.0,  m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.807916,      real(CharmLoops::C0(19.21, m_b)), eps);

                // imag parts
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0( 1.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0( 6.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0(11.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0(16.0, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.000000000,  imag(CharmLoops::C0(19.0, m_b)), eps);
            }
        }
} helper_test;

class FormFactorsTest :
    public TestCase
{
    public:
        FormFactorsTest() :
            TestCase("form_factors_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christoph Bobeth's result from May 2010 */

            /* Formfactors, massless loops */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, eps = 0.0000001;
                TEST_CHECK_NEARLY_EQUAL(- 0.8832611, real(CharmLoops::F17(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.6937322, imag(CharmLoops::F17(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 5.2995666, real(CharmLoops::F27(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 4.1623936, imag(CharmLoops::F27(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 3.3632062, real(CharmLoops::F19(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 6.9078480, imag(CharmLoops::F19(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 3.4455298, real(CharmLoops::F29(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+24.6919276, imag(CharmLoops::F29(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(- 1.2486221, real(CharmLoops::F87(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 2.7925269, imag(CharmLoops::F87(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(- 3.2730189, real(CharmLoops::F89(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(  0.0000000, imag(CharmLoops::F89(mu, s, m_b)), eps);
            }

            /* Formfactors for O_8 are problematic near the zero recoil point */
            {
                static const double mu = 4.2, s = 19.2, m_b = 4.6, eps = 0.0000001;

                TEST_CHECK_NEARLY_EQUAL(- 0.9708796,  real(CharmLoops::F87(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 2.7925268,  imag(CharmLoops::F87(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(- 2.0208146,  real(CharmLoops::F89(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(  0.0000000,  imag(CharmLoops::F89(mu, s, m_b)), eps);
            }

            /* Formfactors, massive loops */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, m_c = 1.2, eps = 0.0000001;
                TEST_CHECK_NEARLY_EQUAL(- 0.73093991, real(CharmLoops::F17(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.17771334, imag(CharmLoops::F17(mu, s, m_b, m_c)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 4.38563254, real(CharmLoops::F27(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 1.06627403, imag(CharmLoops::F27(mu, s, m_b, m_c)), eps);

                TEST_CHECK_NEARLY_EQUAL(-34.40870331, real(CharmLoops::F19(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.25864665, imag(CharmLoops::F19(mu, s, m_b, m_c)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 6.27364439, real(CharmLoops::F29(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 1.55195807, imag(CharmLoops::F29(mu, s, m_b, m_c)), eps);
            }
        }
} two_loop_test;
