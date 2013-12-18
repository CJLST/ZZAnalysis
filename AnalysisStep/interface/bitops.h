#ifndef BITOPS_H
#define BITOPS_H

bool test_bit( int mask, unsigned int iBit );
void set_bit( int& mask, unsigned int iBit );

bool test_bit_16( short mask, unsigned int iBit );
void set_bit_16( short& mask, unsigned int iBit );

bool test_bit_8( char mask, unsigned int iBit );
void set_bit_8( char& mask, unsigned int iBit );


#endif
