/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BRVD2021_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BRVD2021_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/utils/model.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/complex.hh>


namespace eos
{
    class BRvD2021FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            const double _t_p, _t_m, _t_0;
            const double _chi_1m, _chi_0p;
            const double _chi_1p, _chi_0m;

            std::array<UsedParameter, 3> _a_V_time,  _a_V_long,  _a_V_perp;
            //std::array<UsedParameter, 4> _p_A_time,  _p_A_long,  _p_A_perp;
            //std::array<UsedParameter, 4> _p_T_long,  _p_T_perp;
            //std::array<UsedParameter, 4> _p_T5_long, _p_T5_perp;

            //UsedParameter _t_0;

            //UsedParameter m_Jpsi;
            //UsedParameter m_psi2S;

            // Lambda_b parameter
            UsedParameter m_LamB;

            // final state baryon parameters
            UsedParameter m_Lam;

            // Subtraction point for the dispersion relation...
            UsedParameter t_s;

            double beta(const double & param) const;
            double _z(const double & t, const double & t_0) const;
            //double _phi(const double & s, const double chi, const unsigned phiParam[7]) const;


            //Expansion in polynomials orthogonal on the arc of the unit circle (zXY, zXY*)
            double _PGvDV2020(const double z, const double & alpha_0, const double & alpha_1, const double & alpha_2) const;

            std::string _par_name(const std::string &) const;

        public:
            BRvD2021FormFactors(const Parameters &, const Options &);
            ~BRvD2021FormFactors();

            FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters &, const Options &);

            virtual double f_time_v(const double & s) const;
            virtual double f_long_v(const double & s) const;
            virtual double f_perp_v(const double & s) const;

            virtual double f_time_a(const double & s) const;
            virtual double f_long_a(const double & s) const;
            virtual double f_perp_a(const double & s) const;

            virtual double f_long_t(const double & s) const;
            virtual double f_perp_t(const double & s) const;

            virtual double f_long_t5(const double & s) const;
            virtual double f_perp_t5(const double & s) const;

            double _phi(const double & s, const double chi, const unsigned phiParam[7]) const;

    };

}

#endif