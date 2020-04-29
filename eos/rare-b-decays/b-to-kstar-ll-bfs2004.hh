/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_BFS2004_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_BFS2004_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::BFS2004> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_c;
            UsedParameter m_s_MSbar;

            UsedParameter f_B;
            UsedParameter f_Kstar_par;
            UsedParameter f_Kstar_perp;
            UsedParameter lambda_B_p;
            UsedParameter a_1_par;
            UsedParameter a_2_par;
            UsedParameter a_1_perp;
            UsedParameter a_2_perp;

            UsedParameter uncertainty_para;
            UsedParameter uncertainty_perp;
            UsedParameter uncertainty_long;

            UsedParameter uncertainty_xi_perp;
            UsedParameter uncertainty_xi_par;

            double e_q;

            char q;

            bool ccbar_resonance;

            bool use_nlo;

            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &)> qcdf_dilepton_massless_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_charm_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_bottom_case;

            std::string ff_relation;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes();

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            double m_b_PS() const;
            double mu_f() const;
            BToKstarDilepton::DipoleFormFactors dipole_form_factors(const double & q2, const WilsonCoefficients<BToS> & wc) const;
            double norm(const double & q2) const;
            double xi_perp(const double & q2) const;
            double xi_par(const double & q2) const;
    };
}

#endif
