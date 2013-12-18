#include "../interface/bitops.h"

bool test_bit( int mask, unsigned int iBit ) { return (mask >> iBit) & 1; }
void set_bit( int& mask, unsigned int iBit ) { mask |= (1<<iBit); }

bool test_bit_16( short mask, unsigned int iBit ) { return (mask >> iBit) & 1; }
void set_bit_16( short& mask, unsigned int iBit ) { mask |= (1<<iBit); }

bool test_bit_8( char mask, unsigned int iBit ) { return (mask >> iBit) & 1; }
void set_bit_8( char& mask, unsigned int iBit ) { mask |= (1<<iBit); }
