/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <eos/form-factors/parametric-brvd2021.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/complex.hh>



namespace eos
{
    BRvD2021FormFactors::BRvD2021FormFactors(const Parameters & p, const Options & o) :
        _t_p(LambdaBToLambda::tau_p),
        _t_m(LambdaBToLambda::tau_m),
        _t_0(LambdaBToLambda::tau_m),
        _chi_1m( 5.131e-04), // fix values from Bhrarucha/Feldmann/Wick '10
        _chi_0p( 6.204e-03),
        _chi_1p( 3.894e-04),
        _chi_0m(19.421e-03),
        _a_V_time{{
            UsedParameter(p[_par_name("V,t_0")],  *this),
            UsedParameter(p[_par_name("V,t_1")],  *this),
            UsedParameter(p[_par_name("V,t_2")],  *this)
        }},

        _a_V_long{{
            UsedParameter(p[_par_name("V,long_0")],  *this),
            UsedParameter(p[_par_name("V,long_1")],  *this),
            UsedParameter(p[_par_name("V,long_2")],  *this)        }},

        _a_V_perp{{
            UsedParameter(p[_par_name("V,perp_0")],  *this),
            UsedParameter(p[_par_name("V,perp_1")],  *this),
            UsedParameter(p[_par_name("V,perp_2")],  *this)        }},



        // here choice of t_0. Do we use t_0 = t_-?
        //_t_0(p["b->sccbar::t_0"], *this), //placeholder

        //m_Jpsi(p["mass::J/psi"], *this),
        //m_psi2S(p["mass::psi(2S)"], *this),

        m_LamB(p["mass::Lambda_b"], *this),
        m_Lam(p["mass::Lambda"], *this),
        t_s(p["b->sccbar::t_s"], *this) //placeholder

    {
    }

    BRvD2021FormFactors::~BRvD2021FormFactors() = default;

    double
    BRvD2021FormFactors::_z(const double & q2, const double & t_0) const
    {
        return (std::sqrt(_t_p - q2) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - q2) + std::sqrt(_t_p - t_0));
    }



    double
    BRvD2021FormFactors::beta(const double & param) const
    {
        //const double s_plus = pow(m_LamB + m_Lam, 2.0);
        //Check it again

        return pow(1.0 + param/_t_p, 0.5);
    }

    double
    BRvD2021FormFactors::_phi(const double & q2, const double chi, const unsigned phiParam[7]) const
    {
        // Values of a, b, c, d, e, f and g depends on the local form factor parametrization, see notes:
        // FF              a    b    c    d    e    f   g
        // phi_(V,t)       0    0    0    1    2    0   1
        // phi_(V,long)    1    0    1    0    2    1   0
        // phi_(V,perp)    0    1    0    0    1    1   0

        const double m_Lam2  = pow(m_Lam, 2.0 );
        const double m_LamB2  = pow(m_LamB, 2.0 );

        const double s_0    = _t_0;
        const double s_p    = _t_p;
        const double Q2     = this->t_s();
        const auto   z      = _z(q2, s_0);

        const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3], e= phiParam[4], f= phiParam[5], g= phiParam[6];

        const double Nlambda = pow( 4.0 * (s_p-s_0)/(2.0 * pow(M_PI,2.0) * pow(3.0, a) *chi) *pow( (2.0/3.0), b) ,0.5) *pow(m_LamB + m_Lam, c) * pow(m_LamB - m_Lam, d);
        const double phi1 = pow((s_0*pow(1.0 + z,2.0)*(s_0*pow(1.0 + z, 2.0) - 2.0*pow(-1.0 + z,2.0)*(m_Lam2 + m_LamB2)) + (-1.0 + pow(beta(pow(m_LamB2 - m_Lam2,2.0)),2.0)*pow(-1.0 + z,4.0) - pow(z,4.0) + z*(4.0 - 8.0*s_0 + 8.0*m_Lam2 + 8.0*m_LamB2) + pow(z,3.0)*(4.0 - 8.0*s_0 + 8.0*m_Lam2 + 8.0*m_LamB2) - 2.0*pow(z,2.0)*(3.0 + 8.0*s_0 + 8.0*m_Lam2 + 8.0*m_LamB2))*(s_p) + 16.0*pow(z,2.0)*pow(s_p,2.0))/pow(-1.0 + z,4.0),0.25);
        const double phi2 = pow((s_0*pow(1.0 + z,2.0) - 4.0*z*(s_p))/pow(-1. + z,2.0), 0.5);
        const double phi3 = pow((-(s_0*pow(1.0 + z, 2.0)) + (-1.0 + pow(beta(pow(m_LamB - m_Lam, 2.0)),2.0) * pow(-1.0 + z,2.0) + 6.0*z -pow(z,2.0))*(s_p))/pow(-1.0 + z, 2.0),0.5);
        const double phi4 = pow((-(s_0*pow(1.0 + z, 2.0)) + (-1.0 + pow(beta(pow(m_LamB + m_Lam, 2.0)), 2.0) * pow(-1.0 + z,2.0) + 6.0*z -pow(z,2.0))*(s_p))/pow(-1.0 + z, 2.0),0.5);
        const double phi5 = pow((s_0*pow(1.0 + z, 2.0) + (1.0 - pow(beta(Q2),2.0)*pow(-1.0 + z, 2.0) - 6.0*z + pow(z,2.0))*(s_p))/pow(-1.0 + z,2.0),-1.5);

        return Nlambda * pow(1.0 + z, 0.5) * pow(1.-z,-1.5) * phi1 * pow(phi2, -e) * pow(phi3, f) * pow(phi4, g) * phi5;
    }

    double
    BRvD2021FormFactors::_PGvDV2020(const double z, const double & alpha_0, const double & alpha_1, const double & alpha_2) const
    {

        const complex<double> m_B = 5.279;
        const complex<double> m_K = 0.493677;
        const complex<double> m_LamB = 5.61960;
        const complex<double> m_Lam = 1.115683;
        const complex<double> u = power_of<2>(m_LamB + m_Lam);
        const complex<double> _t_pp = power_of<2>(m_B+m_K);
        const complex<double> _t_00 = power_of<2>(m_LamB-m_Lam);

        const complex<double> zLbL = (std::sqrt(_t_pp - u ) - std::sqrt(_t_pp - _t_00)) / (std::sqrt(_t_pp - u) + std::sqrt(_t_pp - _t_00));
        const double alphaXY = std::abs(std::arg(zLbL));

        const double denom = 2.0*pow(alphaXY, 2) + cos(2.0*alphaXY) - 1;

        const double P0z = 1.0/sqrt(2*alphaXY);
        const double P1z = (z - sin(alphaXY)/alphaXY)*sqrt(alphaXY/denom);
        const double P2z = ( z*z + z*sin(alphaXY)*(sin(2*alphaXY)-2*alphaXY)/denom +
                                    2*sin(alphaXY)*(sin(alphaXY)-alphaXY*cos(alphaXY))/denom) *
                                    sqrt( 2*denom/(-9*alphaXY + 8*pow(alphaXY,3) + 8*alphaXY*cos(2*alphaXY) +
                                    alphaXY*cos(4*alphaXY) + 4*sin(2*alphaXY) - 2*sin(4*alphaXY)) );

        return alpha_0*P0z + alpha_1*P1z + alpha_2*P2z;
    }




    std::string
    BRvD2021FormFactors::_par_name(const std::string & ff_name) const
    {
        return std::string("Lambda_b->Lambda::a^") + ff_name + std::string("@BRvD2021");
    }

    FormFactors<OneHalfPlusToOneHalfPlus> *
    BRvD2021FormFactors::make(const Parameters & parameters, const Options & options)
    {
        return new BRvD2021FormFactors(parameters, options);
    }


    double
    BRvD2021FormFactors::f_time_v(const double & q2) const
    {
        const double z        = _z(q2, _t_0);

        // resonances for 0^+
        //const double blaschke = _z(q2, 6.329 * 6.329) * _z(q2, 6.910 * 6.910) * _z(q2, 7.020 * 7.020);
        const double blaschke = _z(q2, power_of<2>(5.367))*_z(q2, power_of<2>(5.416))*_z(q2, power_of<2>(5.711))*_z(q2, power_of<2>(5.750));
        const unsigned phiParam[7] = {0, 0, 0, 1, 2, 0, 1};

        return _PGvDV2020(z, _a_V_time[0], _a_V_time[1], _a_V_time[2]) / _phi(q2,_chi_0p, phiParam) / blaschke;
    }


    double
    BRvD2021FormFactors::f_long_v(const double & q2) const
    {

        const double z        = _z(q2, _t_0);

        // resonances for 0^+
        //const double blaschke = _z(q2, 6.329 * 6.329) * _z(q2, 6.910 * 6.910) * _z(q2, 7.020 * 7.020);
        const double blaschke = _z(q2, power_of<2>(5.367))*_z(q2, power_of<2>(5.416))*_z(q2, power_of<2>(5.711))*_z(q2, power_of<2>(5.750));
        const unsigned phiParam[7] = {1, 0, 1, 0, 2, 1, 0};

        return _PGvDV2020(z, _a_V_long[0], _a_V_long[1], _a_V_long[2]) / _phi(q2,_chi_1m, phiParam) / blaschke;
    }


    double
    BRvD2021FormFactors::f_perp_v(const double & q2) const
    {

        const double z        = _z(q2, _t_0);

        // resonances for 0^+
        //const double blaschke = _z(q2, 6.329 * 6.329) * _z(q2, 6.910 * 6.910) * _z(q2, 7.020 * 7.020);
        const double blaschke = _z(q2, power_of<2>(5.367))*_z(q2, power_of<2>(5.416))*_z(q2, power_of<2>(5.711))*_z(q2, power_of<2>(5.750));
        const unsigned phiParam[7] = {0, 1, 0, 0, 1, 1, 0};

        return _PGvDV2020(z, _a_V_perp[0], _a_V_perp[1], _a_V_perp[2]) / _phi(q2,_chi_1m, phiParam) / blaschke;
    }



    double
    BRvD2021FormFactors::f_time_a(const double & q2) const
    {

        return 0;
    }

    double
    BRvD2021FormFactors::f_long_a(const double & q2) const
    {

        return 0;
    }


    double
    BRvD2021FormFactors::f_perp_a(const double & q2) const
    {

        return 0;
    }




    double
    BRvD2021FormFactors::f_long_t(const double & q2) const
    {

        return 0;
    }


    double
    BRvD2021FormFactors::f_perp_t(const double & q2) const
    {

        return 0;
    }





    double
    BRvD2021FormFactors::f_long_t5(const double & q2) const
    {

        return 0;
    }


    double
    BRvD2021FormFactors::f_perp_t5(const double & q2) const
    {

        return 0;
    }


    /*
    Diagnostics
    BRvD2021FormFactors::diagnostics() const
    {
        Diagnostics results;

        {

            results.add({ f_time_v(10), "f_time_v(10)" });
        }


        return results;
    }
    */


}

