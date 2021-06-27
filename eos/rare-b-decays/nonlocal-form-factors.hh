<<<<<<< HEAD
namespace nonlocal
{
    //struct PToP;
    //struct PToV;
    struct OneHalfPlusToOneHalfPlus;
}
    
template <typename Transition_>
class NonLocalFormFactors;

template
class NonLocalFormFactors<OneHalfPlusToOneHalfPlus> :
    public virtual ParameterUser
{
    public:
        virtual complex<double> H_V_perp(const double & q2) const = 0;
        virtual complex<double> H_V_0(const double & q2) const = 0;
        
        virtual complex<double> H_V_perp_residue_Jpsi() const = 0;
        virtual complex<double> H_V_0_residue_Jpsi() const = 0;
        
        
        virtual complex<double> H_A_perp(const double & q2) const = 0;
        virtual complex<double> H_A_0(const double & q2) const = 0;
        
        virtual complex<double> H_A_perp_residue_Jpsi() const = 0;
        virtual complex<double> H_A_0_residue_Jpsi() const = 0;

};
=======
#include <eos/utils/complex.hh>
#include <eos/utils/parameters.hh>

namespace eos
{

	namespace nonlocal
	{
	    //struct PToP;
	    //struct PToV;
	    struct OneHalfPlusToOneHalfPlus;
	}

	template <typename Transition_>
	class NonLocalFormFactors;

	template<>
	class NonLocalFormFactors<nonlocal::OneHalfPlusToOneHalfPlus> :
	    public virtual ParameterUser
	{
	    public:
	        virtual complex<double> H_V_perp(const double & q2) const = 0;
	        virtual complex<double> H_V_0(const double & q2) const = 0;

	        virtual complex<double> H_V_perp_residue_Jpsi() const = 0;
	        virtual complex<double> H_V_0_residue_Jpsi() const = 0;


	        virtual complex<double> H_A_perp(const double & q2) const = 0;
	        virtual complex<double> H_A_0(const double & q2) const = 0;

	        virtual complex<double> H_A_perp_residue_Jpsi() const = 0;
	        virtual complex<double> H_A_0_residue_Jpsi() const = 0;

	}; 
}
>>>>>>> [rare-b-decays] Add nonlocal formfactor base class
