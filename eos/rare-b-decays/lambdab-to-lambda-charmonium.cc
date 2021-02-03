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

#include <eos/form-factors/mesonic.hh>
#include <eos/rare-b-decays/lambdab-to-lambda-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <complex>

namespace eos
{
    using std::abs;
    using std::arg;
    using std::conj;
    using std::norm;
    using std::real;
    using std::sqrt;

    
    template <>
    struct Implementation<LambdabToLambdaCharmonium>
    {
        UsedParameter g_fermi;

        UsedParameter hbar;

        std::shared_ptr<Model> model;

        SwitchOption q;

        UsedParameter m_LamB;

        UsedParameter tau_LamB;

        UsedParameter m_Lam;

        SwitchOption opt_formfactor;

        NonlocalFormFactorPtr<nc::OneHalfPlusToOneHalfPlus> nonlocal_formfactor;

        SwitchOption psi;

        UsedParameter m_psi;

        UsedParameter f_psi;

        UsedParameter alpha;

        std::function<complex<double> ()> residue_H_V_perp;
        std::function<complex<double> ()> residue_H_V_long;
        std::function<complex<double> ()> residue_H_A_perp;
        std::function<complex<double> ()> residue_H_A_long;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
        	// q(o, "q", { "d", "u" }, "d"), spectator quark
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u),
            model(Model::make(o.get("model", "SM"), p, o)),
            q(o, "q", { "d", "u" }, "d"),
            tau_LamB(p["life_time::Lambda_b"], u),
            m_Lam(p["mass::Lambda"], u),
            m_LamB(p["mass::Lambda_b"], u),
            opt_formfactor(o, "formfactor", { "GvDV2020", "GRvDV2021" }, "GvDV2020"),
            nonlocal_formfactor(NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::make("Lambda_b->Lambda::" + opt_formfactor.value(), p, o)),
            psi(o, "psi", { "J/psi", "psi(2S)" }, "J/psi"),
            m_psi(p["mass::" + psi.value()], u),
            f_psi(p["decay-constant::" + psi.value()], u),
            alpha(p["Lambda::alpha"], u)

        {
            if (! nonlocal_formfactor.get())
                throw InternalError("Cannot construct the nonlocal formfactor");

            if ("J/psi" == psi.value())
            {	
            	// Adapt to onehalfplus FF
                residue_H_V_perp = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_V_perp_residue_jpsi, nonlocal_formfactor);
                residue_H_V_long = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_V_long_residue_jpsi, nonlocal_formfactor);

                residue_H_A_perp = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_A_perp_residue_jpsi, nonlocal_formfactor);
                residue_H_A_long = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_A_long_residue_jpsi, nonlocal_formfactor);
            }
            else
            {
                residue_H_V_perp = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_V_perp_residue_psi2s, nonlocal_formfactor);
                residue_H_V_long = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_V_long_residue_psi2s, nonlocal_formfactor);

                residue_H_A_perp = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_A_perp_residue_psi2s, nonlocal_formfactor);
                residue_H_A_long = std::bind(&NonlocalFormFactor<nc::OneHalfPlusToOneHalfPlus>::H_A_long_residue_psi2s, nonlocal_formfactor);
            }

            u.uses(*model);
            u.uses(*nonlocal_formfactor);
        }

        ~Implementation() = default;

        struct Amplitudes
        {
            complex<double> A_V_perp, A_V_long, A_A_perp, A_A_long;
        };


        Amplitudes amplitudes() const
        {
           
            const double Q_c = 2.0/3.0;

            complex<double> A_V_perp = residue_H_V_perp() / (Q_c* f_psi * m_psi);
            complex<double> A_V_long = residue_H_V_long() / (Q_c* f_psi * m_psi);
            complex<double> A_A_perp = residue_H_A_perp() / (Q_c* f_psi * m_psi);
            complex<double> A_A_long = residue_H_A_long() / (Q_c* f_psi * m_psi);

            return { A_V_perp, A_V_long, A_A_perp, A_A_long };
        }
        
        double s_minus(const double & q) const
        {
            const double m_LamB = this->m_LamB();
            const double m_Lam = this->m_Lam();
            return pow(m_LamB - m_Lam, 2.0) - pow(q, 2.0); 
        }

        double s_plus(const double & q) const
        {
            const double m_LamB = this->m_LamB();
            const double m_Lam = this->m_Lam();
            return pow(m_LamB + m_Lam, 2.0) - pow(q ,2.0); 
        }
        double branching_ratio() const
        {
            const auto amps = amplitudes();
            const double m_LamB = this->m_LamB(); 
            const double m_Lam = this->m_Lam(); 
            const double tau_LamB = this->tau_LamB();
            const double m_psi = this->m_psi();
            const auto lambda = eos::lambda(pow(m_LamB, 2), pow(m_Lam, 2), pow(m_psi, 2));
            /*
            const auto prefactor = pow(g_fermi * abs(model->ckm_cb() * (model->ckm_cs())), 2)
                    * 6 * tau_LamB / (32 * M_PI) * m_LamB * sqrt(lambda);
            */
            const auto prefactor = 6 * tau_LamB / (32 * M_PI) * m_LamB * sqrt(lambda);

            const auto amps_res = s_minus(m_psi)* (pow(m_LamB + m_Lam, 2.0)/(2.0 * pow(m_psi, 2.0)) * norm(amps.A_V_long) + norm(amps.A_V_perp))
                                + s_plus(m_psi) *  (pow(m_LamB - m_Lam, 2.0)/(2.0 * pow(m_psi, 2.0)) * norm(amps.A_A_long) + norm(amps.A_A_perp));
            return prefactor * amps_res;
        }

        double residue_norm() const
        {
        
            return norm(residue_H_V_perp()) + norm(residue_H_V_long()) + norm(residue_H_A_perp()) +norm(residue_H_A_long());
        }
        double K1ss() const
        {

            return 1.0/residue_norm() * (norm(residue_H_V_perp()) + norm(residue_H_A_perp()) + 2.0 * norm(residue_H_V_long()) + 2.0 * norm(residue_H_A_long()));       
        }
        
        double K1cc() const
        {

            return 1.0/residue_norm() * (norm(residue_H_V_perp()) + norm(residue_H_A_perp()));
        }

        double K2ss() const
        {
            const double alpha = 0.20; 

            return alpha/residue_norm() * real(residue_H_V_perp() * conj(residue_H_A_perp()) + 2.0 * residue_H_V_long() * conj(residue_H_A_long()) );
        }

        double K2cc() const
        {
            
            const double alpha = 0.20;

            return 2.0 * alpha/residue_norm() * real(residue_H_V_perp() * conj(residue_H_A_perp()));
        }

        double K3sc() const
        {
            const double alpha = 0.20;

            return 2.0 *alpha/(pow(2.0, 0.5) * residue_norm()) * imag(residue_H_V_perp() * conj(residue_H_V_long()));
        }

        double K4sc() const
        {
            const double alpha = 0.20;

            return 2.0 *alpha/(pow(2.0, 0.5) * residue_norm()) * real(residue_H_V_perp() * conj(residue_H_A_long()) - residue_H_A_perp() * conj(residue_H_V_long()));
        }
    };

    LambdabToLambdaCharmonium::LambdabToLambdaCharmonium(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdabToLambdaCharmonium>(new Implementation<LambdabToLambdaCharmonium>(p, o, *this))
    {
    }

    LambdabToLambdaCharmonium::~LambdabToLambdaCharmonium() = default;

    double
    LambdabToLambdaCharmonium::branching_ratio() const
    {
        return _imp->branching_ratio();
    }
    

    double
    LambdabToLambdaCharmonium::K1ss() const
    {
        return _imp->K1ss();
    }

}
