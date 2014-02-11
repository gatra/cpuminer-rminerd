/*
 * Copyright 2011 ArtForz
 * Copyright 2011-2013 pooler
 * Copyright 2013 gatra
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.  See COPYING for more details.
 */

#include "cpuminer-config.h"
#include "miner.h"

#include <string.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>


extern void sha256d_80_swap(uint32_t *hash, const uint32_t *data);

static unsigned int largestPrimeInTable;
static unsigned int *primeTable;

static unsigned int primeTableAllocatedSize;
static unsigned int primeTableSize;


    int isPrime( int candidate )
    {
        for( unsigned int i = 0; i < primeTableSize; i++ )
        {
            const int prime = primeTable[i];
            if( prime * prime > candidate )
                return 1;
            if( (candidate % prime) == 0 )
                return 0;
        }
        return 1;
    }

int initPrimeTable( void )
{
	opt_sieve_size &= ~(unsigned int)32;
	opt_max_prime &= ~(unsigned int)32;
	largestPrimeInTable = opt_max_prime; //200*sqrt(opt_sieve_size);
	double d = largestPrimeInTable * 1.25506; // upper bound for PrimePi(n) - http://oeis.org/A209883
	d /= log(largestPrimeInTable);
	primeTableAllocatedSize = (unsigned int)(ceil(d) + 1); // upper bound on number of primes less than largestPrimeInTable
	primeTable = (unsigned int *) malloc( sizeof(unsigned int) * primeTableAllocatedSize );

	applog(LOG_INFO, "allocated space for %u primes in table", primeTableAllocatedSize);
	
	primeTableSize = 1;
	primeTable[0] = 2;
	for( int i = 3; i <= largestPrimeInTable; i += 2 )
	{
		if( isPrime(i) )
		{
			if( primeTableSize >= primeTableAllocatedSize )
			{
				applog(LOG_ERR, "primes don't fit allocated space"); // should never happen
				return 1;
			}
			primeTable[primeTableSize++] = i;
		}
	}
	applog(LOG_INFO, "using %u primes, largest prime in table is %u", primeTableSize, primeTable[primeTableSize-1]);
	return 0;
}
// end of init


void sieveReset( uint32_t *pSieve, int index )
{
	pSieve[index>>5] &= ~(  1U << (index & 0x1f)  );
}
uint32_t sieveGet( uint32_t *pSieve, int index )
{
	return pSieve[index>>5] & (  1U << (index & 0x1f)  );
}
void init( int *pSieve, const mpz_t base )
{
        memset( pSieve, 0xff, opt_sieve_size/8 );

        for( unsigned int primeIndex = 0; primeIndex < primeTableSize; primeIndex++ )
        {
            int prime = primeTable[primeIndex];
            
            for( int sieveIndex = prime - mpz_fdiv_ui(base, prime);
                 sieveIndex < opt_sieve_size;
                 sieveIndex += prime )
            {
        			sieveReset(pSieve, sieveIndex);

            }
        }
}

int getNext( int *pSieve, int _index )
{
        if( _index >= opt_sieve_size - 17 )
        {
            return -1;
        }
        while( 1 )
        {
            _index++;
            if( _index >= opt_sieve_size - 17 )
            {
                return -1;
            }
            if( sieveGet(pSieve, _index) && 
                sieveGet(pSieve, _index+4) && 
                sieveGet(pSieve, _index+6) && 
                sieveGet(pSieve, _index+10) &&
                sieveGet(pSieve, _index+12) && 
                sieveGet(pSieve, _index+16) )
                return _index;
        }
}


const int zeroesBeforeHashInPrime = 8;

void SetCompact(mpz_t p, uint32_t nCompact)
{
        unsigned int nSize = nCompact >> 24;
        //bool fNegative     =(nCompact & 0x00800000) != 0;
        unsigned int nWord = nCompact & 0x007fffff;
        if (nSize <= 3)
        {
            nWord >>= 8*(3-nSize);
            mpz_init_set_ui (p, nWord);
        }
        else
        {
            mpz_init_set_ui (p, nWord);
            mpz_mul_2exp(p, p, 8*(nSize-3));
        }
        //BN_set_negative(p, fNegative);
}

#define HASH_LEN_IN_BITS 256

unsigned int generatePrimeBase( mpz_t bnTarget, uint32_t *hash, uint32_t compactBits )
{
	int i;

    mpz_set_ui(bnTarget, 1);
    mpz_mul_2exp(bnTarget,bnTarget,zeroesBeforeHashInPrime);

    for ( i = 0; i < HASH_LEN_IN_BITS; i++ )
    {
	mpz_mul_2exp(bnTarget,bnTarget,1);
	mpz_add_ui(bnTarget,bnTarget,(hash[i/32] & 1));
        hash[i/32] >>= 1;
    }
    mpz_t nBits;
    SetCompact( nBits, compactBits );
    unsigned int trailingZeros = mpz_get_ui (nBits) - 1 - zeroesBeforeHashInPrime - HASH_LEN_IN_BITS;
    mpz_mul_2exp(bnTarget,bnTarget,trailingZeros);
    mpz_clear(nBits);
    return trailingZeros;
}

double riecoin_time_to_block( double hashrate, uint32_t compactBits, int primes )
{
/*

Name of prime k-tuple	Width	Pattern	Hardy-Littlewood constant Hk	Reciprocal Ck = 1/Hk
2-tuple, twin primes, twins	 2	0 2	1.32032	0.757392
3-tuple, triplet (type A)	 6	0 4 6	2.85825	0.349864
3-tuple, triplet (type B)	 6	0 2 6	2.85825	0.349864
4-tuple, quadruplet	 8	0 2 6 8	4.15118	0.240895
5-tuple, quintuplet (type A)	 12	0 4 6 10 12	10.13179	0.09869924
5-tuple, quintuplet (type B)	 12	0 2 6 8 12	10.13179	0.09869924
6-tuple, sextuplet	 16	0 4 6 10 12 16	17.29861	0.05780811
7-tuple, septuplet (type A)	 20	0 2 6 8 12 18 20	53.97195	0.01852814
7-tuple, septuplet (type B)	 20	0 2 8 12 14 18 20	53.97195	0.01852814
8-tuple, octuplet (type A)	 26	0 2 6 8 12 18 20 26	178.26195	0.005609722
8-tuple, octuplet (type B)	 26	  0 6 8 14 18 20 24 26  	178.26195	0.005609722
8-tuple, octuplet (type C)	 26	  0 2 6 12 14 20 24 26  	475.36521	0.002103646
*/
	mpz_t nBits;
	double f;
	static const double l2 = 0.69314718056; // ln(2)

	SetCompact( nBits, compactBits );
	f = mpz_get_ui (nBits);
	mpz_clear(nBits);
	
	f = pow( f * l2, primes ) / hashrate;
	if( primes == 6 )
		return f * 0.05780811; // reciprocal of Hardy-Littlewood constant H6
	if( primes == 5 )
		return f * 0.09869924; // reciprocal of Hardy-Littlewood constant H5 (type b)
	if( primes == 4 )
		return f * 0.240895; // reciprocal of Hardy-Littlewood constant H4 (delta=4)
	if( primes == 3 )
		return f * 0.349864; // reciprocal of Hardy-Littlewood constant H3 (type a)
	if( primes == 2 )
		return f * 0.757392; // reciprocal of Hardy-Littlewood constant H2 ???
	return 0;

}

#define MR_TESTS 12


int scanhash_riecoin(int thr_id, uint32_t *pdata, const int primes,
	uint32_t max_nonce, unsigned long *hashes_done, uint32_t *pSieve)
{
	uint32_t hash[8] __attribute__((aligned(32)));
	uint32_t n = pdata[RIECOIN_DATA_NONCE];
	const uint32_t first_nonce = n;
	mpz_t b;
	int sieveIndex = 0;

	mpz_t bnTarget;
	mpz_init( bnTarget );

	if( max_nonce <= opt_sieve_size )
	{
		max_nonce = 0;
	}
	else
	{
		max_nonce -= opt_sieve_size;
	}

	mpz_init(b);
	
	uint32_t aux = pdata[RIECOIN_DATA_NONCE];
	pdata[RIECOIN_DATA_NONCE] = 0x80000000UL;
	sha256d_80_swap(hash, pdata);
	pdata[RIECOIN_DATA_NONCE] = aux;

	generatePrimeBase( b, hash, swab32(pdata[RIECOIN_DATA_DIFF]) );

	mpz_add_ui( b, b, n );

	do {
		init( pSieve, b );
		
		while( (sieveIndex = getNext(pSieve, sieveIndex) ) > 0 )
                {
			mpz_set( bnTarget, b );
			mpz_add_ui( bnTarget, bnTarget, sieveIndex );
                	if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
                	{
				mpz_add_ui( bnTarget, bnTarget, 4 );
	                	if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
	                	{
				mpz_add_ui( bnTarget, bnTarget, 2 );
	                	if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
        	        	{
					mpz_add_ui( bnTarget, bnTarget, 4 );
                    			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
                    			{
						mpz_add_ui( bnTarget, bnTarget, 2 );
                        			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
                        			{
						mpz_add_ui( bnTarget, bnTarget, 4 );
                        			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) || primes < 6 )
                        			{
							pdata[RIECOIN_DATA_NONCE] = n + sieveIndex;
							pdata[RIECOIN_DATA_NONCE] = swab32(pdata[RIECOIN_DATA_NONCE]);
							*hashes_done = n + sieveIndex - first_nonce + 1;
							mpz_clear(bnTarget);
							mpz_clear(b);
							return 1;
						} }
					}
				} }
			}
		}
		
		n += opt_sieve_size;
		mpz_add_ui(b, b, opt_sieve_size );
	} while (n < max_nonce && !work_restart[thr_id].restart);
	
	*hashes_done = n - first_nonce + 1;
	pdata[RIECOIN_DATA_NONCE] = n;
	mpz_clear(bnTarget);
	mpz_clear(b);
	return 0;
}
