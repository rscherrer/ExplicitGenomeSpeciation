//
// Created by p278834 on 9-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_ACCESSORY_H
#define EXPLICITGENOMESPECIATION_ACCESSORY_H

#include <string>
#include <vector>

// Forward declaration
typedef std::pair<double, double> TradeOffPt;

size_t divide_rounding_up( std::size_t dividend, std::size_t divisor )
{ return ( dividend + divisor - 1 ) / divisor; }

std::string to_string( std::vector< bool > const & bitvector ) {
std::string ret( divide_rounding_up( bitvector.size(), 8 ), 0 );
auto out = ret.begin();
int shift = 0;

for ( bool bit : bitvector ) {
* out |= bit << shift;

if ( ++ shift == 8 ) {
++ out;
shift = 0;
}
}
return ret;
}



#endif //EXPLICITGENOMESPECIATION_ACCESSORY_H
