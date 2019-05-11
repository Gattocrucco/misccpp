#ifndef UNCERTAINTIES_IO_HPP_6DDCDE20
#define UNCERTAINTIES_IO_HPP_6DDCDE20

/*! \file
\brief Defines stream operations on `UReal`s.
*/

#include <ostream>

namespace uncertainties {
    template<typename Real>
    class UReal;
    
    template<typename Real, typename CharT>
    std::basic_ostream<CharT> &operator<<(std::basic_ostream<CharT> &stream,
                                          const UReal<Real> &x) {
        return stream << x.format();
    }
}

#endif /* end of include guard: UNCERTAINTIES_IO_HPP_6DDCDE20 */
