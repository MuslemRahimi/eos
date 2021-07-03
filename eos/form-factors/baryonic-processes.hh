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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_PROCESSES_HH
#define EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_PROCESSES_HH 1

namespace eos
{
    /* 1/2^+ -> 1/2^+ Processes */

    struct LambdaBToLambda {
        using Transition = OneHalfPlusToOneHalfPlus;
        static constexpr const char * label = "Lambda_b->Lambda";
        static constexpr const double m_B = 5.279;
        static constexpr const double m_K = 0.493677;
        static constexpr const double m_LamB = 5.61960;
        static constexpr const double m_Lam = 1.115683;
        static constexpr const double m2_Br1m = 5.415 * 5.415; // B_s^*
        static constexpr const double m2_Br0p = 5.630 * 5.630; // B_s scalar
        static constexpr const double tau_p = (m_B + m_K) * (m_B + m_K);
        static constexpr const double tau_m = (m_LamB - m_Lam) * (m_LamB - m_Lam);
        static constexpr const bool uses_tensor_form_factors = true;
    };
}

#endif