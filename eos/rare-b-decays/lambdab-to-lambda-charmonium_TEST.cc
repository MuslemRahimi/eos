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
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 9.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
            p["life_time::Lambda_b"]                     = 1.471e-12;
            p["decay-constant::J/psi"]                   = 0.2773;

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
            oo.set("model",          "WilsonScan");
            oo.set("q",              "d");
            oo.set("formfactor",     "GvDV2020");
            oo.set("psi",            "J/psi");

            LambdabToLambdaCharmonium c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.branching_ratio(),  0.1687825877, eps);
            /*
            TEST_CHECK_RELATIVE_ERROR(c.S_1c_LHCb(),  0.1687825877, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_1s_LHCb(),  0.6234130592, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_3_LHCb(),  -0.2413178232, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_4_LHCb(),   0.2354340047, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_8_LHCb(),  -0.0062319175, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_9_LHCb(),   0.0091313636, eps);
            */

        }
} lambdab_to_lambda_charmonium_BRvD2021_test;
