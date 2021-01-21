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

#include <eos/utils/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options-impl.hh>
#include <cmath>
#include <eos/utils/options.hh>
#include <eos/utils/stringify.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

#include <map>

namespace eos
{

    namespace nc
    {
        struct LambdabToLambda
        {
            constexpr static const char * label = "Lambda_b->Lambda";
        };
        constexpr const char * LambdabToLambda::label;
    }

     namespace nc_p_to_p
    {
       template <typename Process_>
        class BRvD2021 :
            public NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>
        {
            private:

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_V_long;
                UsedParameter im_alpha_0_V_long;
                UsedParameter re_alpha_1_V_long;
                UsedParameter im_alpha_1_V_long;
                UsedParameter re_alpha_2_V_long;
                UsedParameter im_alpha_2_V_long;

                UsedParameter re_alpha_0_V_perp;
                UsedParameter im_alpha_0_V_perp;
                UsedParameter re_alpha_1_V_perp;
                UsedParameter im_alpha_1_V_perp;
                UsedParameter re_alpha_2_V_perp;
                UsedParameter im_alpha_2_V_perp;

                UsedParameter re_alpha_0_A_long;
                UsedParameter im_alpha_0_A_long;
                UsedParameter re_alpha_1_A_long;
                UsedParameter im_alpha_1_A_long;
                UsedParameter re_alpha_2_A_long;
                UsedParameter im_alpha_2_A_long;

                UsedParameter re_alpha_0_A_perp;
                UsedParameter im_alpha_0_A_perp;
                UsedParameter re_alpha_1_A_perp;
                UsedParameter im_alpha_1_A_perp;
                UsedParameter re_alpha_2_A_perp;
                UsedParameter im_alpha_2_A_perp;


                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // Lambda_b parameter
                UsedParameter m_LamB;

                // final state baryon parameters
                UsedParameter m_Lam;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

            public:
                BRvD2021(const Parameters & p, const Options & o) :

                re_alpha_0_V_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^V_long}@BRvD2021"], *this),
                im_alpha_0_V_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^V_long}@BRvD2021"], *this),
                re_alpha_1_V_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^V_long}@BRvD2021"], *this),
                im_alpha_1_V_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^V_long}@BRvD2021"], *this),
                re_alpha_2_V_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^V_long}@BRvD2021"], *this),
                im_alpha_2_V_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^V_long}@BRvD2021"], *this),

                re_alpha_0_V_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^V_perp}@BRvD2021"], *this),
                im_alpha_0_V_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^V_perp}@BRvD2021"], *this),
                re_alpha_1_V_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^V_perp}@BRvD2021"], *this),
                im_alpha_1_V_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^V_perp}@BRvD2021"], *this),
                re_alpha_2_V_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^V_perp}@BRvD2021"], *this),
                im_alpha_2_V_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^V_perp}@BRvD2021"], *this),

                re_alpha_0_A_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^A_long}@BRvD2021"], *this),
                im_alpha_0_A_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^A_long}@BRvD2021"], *this),
                re_alpha_1_A_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^A_long}@BRvD2021"], *this),
                im_alpha_1_A_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^A_long}@BRvD2021"], *this),
                re_alpha_2_A_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^A_long}@BRvD2021"], *this),
                im_alpha_2_A_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^A_long}@BRvD2021"], *this),

                re_alpha_0_A_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^A_perp}@BRvD2021"], *this),
                im_alpha_0_A_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^A_perp}@BRvD2021"], *this),
                re_alpha_1_A_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^A_perp}@BRvD2021"], *this),
                im_alpha_1_A_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^A_perp}@BRvD2021"], *this),
                re_alpha_2_A_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^A_perp}@BRvD2021"], *this),
                im_alpha_2_A_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^A_perp}@BRvD2021"], *this),

                m_Jpsi(p["mass::J/psi"], *this),
                m_psi2S(p["mass::psi(2S)"], *this),

                m_LamB(p["mass::Lambda_b"], *this),
                m_Lam(p["mass::Lambda"], *this),

                m_D0(p["mass::D^0"], *this),
                t_0(p["b->sccbar::t_0"], *this),

                t_s(p["b->sccbar::t_s"], *this),
                chiOPE(p["b->sccbar::chiOPE@BRvD2021"], *this)

                {
                }

                ~BRvD2021() = default;

                inline double beta(const double param) const
                {
                    const double s_plus = 4.0 * pow(m_D0, 2.0);

                    return pow(1.0 + param/s_plus, 0.5);
                }

                inline complex<double> phi(const double & q2, const unsigned phiParam[5]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF              a    b    c    d    e
                    // phi_(V,long)    1    0    4    1    0
                    // phi_(V,perp)    0    0    3    1    0
                    // phi_(A,long)    0    1    4    0    1
                    // phi_(A,perp)    0    0    3    0    1

                    const double m_Lam2  = pow(m_Lam, 2);
                    const double m_LamB2  = pow(m_LamB, 2);
                    const double m_D02 = pow(m_D0, 2);

                    const double s_0    = this->t_0();
                    const double s_plus = 4.0 * m_D02;
                    const double Q2     = this->t_s();
                    const double chi    = this->chiOPE();
                    const auto   z      = eos::nc_utils::z(q2, s_plus, s_0);

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3], e= phiParam[4];

                    const double Nlambda = 8*M_PI * pow((s_plus-s_0)/(3*chi),0.5) * m_LamB2 * pow(m_LamB + m_Lam, a) * pow(m_LamB - m_Lam, b); 
                    const complex<double> phi1 = pow((s_0*pow(1.0 + z,2.0)*(s_0*pow(1.0 + z,2.0) - 2.0*pow(-1.0 + z,2.0)*(m_Lam2 + m_LamB2)) + (-1.0 + pow(beta(pow(m_LamB2 - m_Lam2,2.0)),2.0)*pow(-1.0 + z,4.0) - pow(z,4.0) + z*(4.0 - 8.0*s_0 + 8.0*m_Lam2 + 8.0*m_LamB2) + pow(z,3.0)*(4.0 - 8.0*s_0 + 8.0*m_Lam2 + 8.0*m_LamB2) - 2.0*pow(z,2.0)*(3.0 + 8.0*s_0 + 8.0*m_Lam2 + 8.0*m_LamB2))*(s_plus) + 16.0*pow(z,2.0)*pow(s_plus,2.0))/pow(-1.0 + z,4.0),0.25);
                    const complex<double> phi2 = pow((s_0*pow(1. + z,2.0) - 4.0*z*(s_plus))/pow(-1. + z,2.0), 0.5);
                    const complex<double> phi3 = pow((-(s_0*pow(1.0 + z, 2.0)) + (-1.0 + pow(beta(pow(m_LamB - m_Lam, 2.0)),2.0) * pow(-1.0 + z,2.0) + 6.0*z -pow(z,2.0))*(s_plus))/pow(-1.0 + z, 2.0),0.5);
                    const complex<double> phi4 = pow((-(s_0*pow(1.0 + z, 2.0)) + (-1.0 + pow(beta(pow(m_LamB + m_Lam, 2.0)), 2.0) * pow(-1.0 + z,2.0) + 6.0*z -pow(z,2.0))*(s_plus))/pow(-1.0 + z, 2.0),0.5);
                    const complex<double> phi5 = pow((s_0*pow(1.0 + z, 2.0) + (1.0 - pow(beta(Q2),2.0)*pow(-1.0 + z, 2.0) - 6.0*z + pow(z,2.0))*(s_plus))/pow(-1.0 + z,2.0),-1.5);

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z,-1.5) * phi1 * pow(phi2, -c) * pow(phi3, d) * pow(phi4, e) * phi5;
                }


                inline complex<double> H_residue_jpsi(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                  const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zLbL    = eos::nc_utils::z(pow(m_LamB + m_Lam, 2.0), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,                  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,                 s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z_Jpsi, zLbL, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi*std::conj(z_psi2S)) / (z_Jpsi - z_psi2S);
                };

                inline complex<double> H_residue_psi2s(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                 const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zLbL    = eos::nc_utils::z(pow(m_LamB + m_Lam, 2.0), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,                  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,                 s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z_psi2S, zLbL, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S*std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi);
                };


                virtual complex<double> H_V_long(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_V_long, im_alpha_0_V_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_V_long, im_alpha_1_V_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_V_long, im_alpha_2_V_long);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,                       s_p, s_0);
                    const auto zLbL    = eos::nc_utils::z(pow(m_LamB + m_Lam, 2.0), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),           s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2),          s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {1,0,4,1,0};

                    return eos::nc_utils::PGvDV2020(z, zLbL, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> H_V_perp(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_V_perp, im_alpha_0_V_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_V_perp, im_alpha_1_V_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_V_perp, im_alpha_2_V_perp);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,                       s_p, s_0);
                    const auto zLbL    = eos::nc_utils::z(pow(m_LamB + m_Lam, 2.0), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),           s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2),          s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {0,0,3,1,0};

                    return eos::nc_utils::PGvDV2020(z, zLbL, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> H_A_long(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_A_long, im_alpha_0_A_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_A_long, im_alpha_1_A_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_A_long, im_alpha_2_A_long);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,                       s_p, s_0);
                    const auto zLbL    = eos::nc_utils::z(pow(m_LamB + m_Lam, 2.0), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),           s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2),          s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {0,1,4,0,1};

                    return eos::nc_utils::PGvDV2020(z, zLbL, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                    ;
                }

                virtual complex<double> H_A_perp(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_A_perp, im_alpha_0_A_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_A_perp, im_alpha_1_A_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_A_perp, im_alpha_2_A_perp);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,                       s_p, s_0);
                    const auto zLbL    = eos::nc_utils::z(pow(m_LamB + m_Lam, 2.0), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),           s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2),          s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {0,0,3,0,1};

                    return eos::nc_utils::PGvDV2020(z, zLbL, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }


                virtual complex<double> H_V_long_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_V_long, im_alpha_0_V_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_V_long, im_alpha_1_V_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_V_long, im_alpha_2_V_long);

                    const unsigned phiParam[5] = {1,0,4,1,0};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_V_long_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_V_long, im_alpha_0_V_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_V_long, im_alpha_1_V_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_V_long, im_alpha_2_V_long);

                    const unsigned phiParam[5] = {1,0,4,1,0};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_V_perp_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_V_perp, im_alpha_0_V_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_V_perp, im_alpha_1_V_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_V_perp, im_alpha_2_V_perp);

                    const unsigned phiParam[5] = {0,0,3,1,0};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_V_perp_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_V_perp, im_alpha_0_V_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_V_perp, im_alpha_1_V_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_V_perp, im_alpha_2_V_perp);

                    const unsigned phiParam[5] = {0,0,3,1,0};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_A_long_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_A_long, im_alpha_0_A_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_A_long, im_alpha_1_A_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_A_long, im_alpha_2_A_long);

                    const unsigned phiParam[5] = {0,1,4,0,1};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_A_long_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_A_long, im_alpha_0_A_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_A_long, im_alpha_1_A_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_A_long, im_alpha_2_A_long);

                    const unsigned phiParam[5] = {0,1,4,0,1};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_A_perp_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_A_long, im_alpha_0_A_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_A_long, im_alpha_1_A_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_A_long, im_alpha_2_A_long);

                    const unsigned phiParam[5] = {0,1,4,0,1};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_A_perp_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_A_perp, im_alpha_0_A_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_A_perp, im_alpha_1_A_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_A_perp, im_alpha_2_A_perp);

                    const unsigned phiParam[5] = {0,0,3,0,1};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };


                static NonlocalFormFactorPtr<nc::OneHalfPlusToOneHalfPlus> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::OneHalfPlusToOneHalfPlus>(new BRvD2021(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const double s_0   = this->t_0();
                    results.add({ real(eos::nc_utils::z(10.0, 4.0 * pow(m_D0, 2), s_0)), "real(z(q2 = 10.0))" });
                    results.add({ imag(eos::nc_utils::z(10.0, 4.0 * pow(m_D0, 2), s_0)), "imag(z(q2 = 10.0))" });

                    return results;
                }
        };
    }

    NonlocalFormFactorPtr<nc::OneHalfPlusToOneHalfPlus>
    NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<nc::OneHalfPlusToOneHalfPlus> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // parametrizations
            std::make_pair("Lambda_b->Lambda::BRvD2021",     &nc_p_to_p::BRvD2021<nc::LambdabToLambda>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<nc::OneHalfPlusToOneHalfPlus>(nullptr);
        }

        return i->second(p, o);
    }
}