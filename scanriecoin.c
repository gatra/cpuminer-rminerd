/*
 * Copyright 2011 ArtForz
 * Copyright 2011-2013 pooler
 * Copyright 2013 gatra
 *
 * function single_modinv copied from jhPrimeminer-hg5fm
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


#undef REPORT_TESTS
//#define REPORT_TESTS

extern void sha256d_80_swap(uint32_t *hash, const uint32_t *data);


static unsigned int largestPrimeInTable;
static unsigned int *primeTable;
static unsigned int *primeTableInverses;
static unsigned int primeTableAllocatedSize;
static unsigned int primeTableSize;

#if 0
const unsigned int fixedMod = 210;
const unsigned int fixedModDelta = 97;

const unsigned int startingPrimeIndex = 4;

//const unsigned int fixedMod = 2310;
//const unsigned int fixedModDelta = 937;
#else

const unsigned int fixedMod = 2310;
const unsigned int fixedModDelta = 97;

const unsigned int startingPrimeIndex = 5;

const unsigned int efficiencyDivisor = 11-6;

#endif


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


int single_modinv (int a, int modulus)
{ /* start of single_modinv */
	int ps1, ps2, parity, dividend, divisor, rem, q, t;
	
	a = a % modulus;
	
	q = 1;
	rem = a;
	dividend = modulus;
	divisor = a;
	ps1 = 1;
	ps2 = 0;
	parity = 0;
	while (divisor > 1)
	{
		rem = dividend - divisor;
		t = rem - divisor;
		if (t >= 0) {
			q += ps1;
			rem = t;
			t -= divisor;
			if (t >= 0) {
				q += ps1;
				rem = t;
				t -= divisor;
				if (t >= 0) {
					q += ps1;
					rem = t;
					t -= divisor;
					if (t >= 0) {
						q += ps1;
						rem = t;
						t -= divisor;
						if (t >= 0) {
							q += ps1;
							rem = t;
							t -= divisor;
							if (t >= 0) {
								q += ps1;
								rem = t;
								t -= divisor;
								if (t >= 0) {
									q += ps1;
									rem = t;
									t -= divisor;
									if (t >= 0) {
										q += ps1;
										rem = t;
										if (rem >= divisor) {
											q = dividend/divisor;
											rem = dividend - q * divisor;
											q *= ps1;
										}}}}}}}}}
		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}

	if (parity == 0)
		return (ps1);
	else
		return (modulus - ps1);
} /* end of single_modinv */


int initPrimeTable( void )
{
	opt_sieve_size &= ~(unsigned int)32;
	opt_max_prime &= ~(unsigned int)32;
	largestPrimeInTable = opt_max_prime; //200*sqrt(opt_sieve_size);
	double d = largestPrimeInTable * 1.25506; // upper bound for PrimePi(n) - http://oeis.org/A209883
	d /= log(largestPrimeInTable);
	primeTableAllocatedSize = (unsigned int)(ceil(d) + 1); // upper bound on number of primes less than largestPrimeInTable
	primeTable = (unsigned int *) malloc( sizeof(unsigned int) * primeTableAllocatedSize );
	primeTableInverses = (unsigned int *) malloc( sizeof(unsigned int) * primeTableAllocatedSize );
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
			if( i > 7 )
			{
				primeTableInverses[primeTableSize] = i - single_modinv( fixedMod, i );
			}
			primeTable[primeTableSize++] = i;
		}
	}
	applog(LOG_INFO, "using %u primes, largest prime in table is %u", primeTableSize, primeTable[primeTableSize-1]);
	return 0;
}
// end of init


struct sieveData
{
	uint32_t *pSieve;
	int x;
};

void sieveReset( uint32_t *pSieve, int index )
{
	pSieve[index>>5] &= ~(  1U << (index & 0x1f)  );
}
uint32_t sieveGet( uint32_t *pSieve, int index )
{
	return pSieve[index>>5] & (  1U << (index & 0x1f)  );
}

uint32_t i2d( struct sieveData *pSieve, uint32_t index ) // sieve index to number offset (delta)
{
	return index * fixedMod + pSieve->x;
}

void init( struct sieveData *pSieve, const mpz_t base )
{
        memset( pSieve->pSieve, 0xff, opt_sieve_size/8 );
        
/*
        base = b30 (30)
        base = bp  (p)

        base + x = 11 (30)          base + x + q * 30 = 0 (p)   q= -(base+x) * 30^-1 (p)
        x = 11 - b30 (30) 
        let bx = base + x
        
        i2d = i*30 + x
        d2i = ( d - x ) / 30
        
        base + xx = 0 (p)
        xx = -bp (p)
        xx = x (30)

        base + x = 11 (30)
        x = 11 - b30 (30) 
        let bx = base + x
        
        bx + 30 * i = 0 (p)
        i = -bx * 30^-1 (p)
        
        i = q * p - bx * 30^-1 
        
        
        base + x + 30*i = 0 (p)

        d = 30^-1 * (-base-x) (p)
*/
	int x = (fixedMod + fixedModDelta - mpz_fdiv_ui(base, fixedMod)) % fixedMod;
	pSieve->x = x;

        for( unsigned int primeIndex = startingPrimeIndex; primeIndex < primeTableSize; primeIndex++ )
        {
		unsigned int prime = primeTable[primeIndex];
		unsigned int bp = mpz_fdiv_ui(base, prime);
		unsigned int k = primeTableInverses[primeIndex];
		unsigned int k2 = ((4 * k)%prime);
		unsigned int k3 = ((6 * k)%prime);
		unsigned int k4 = ((10 * k)%prime);
		unsigned int k5 = ((12 * k)%prime);
		unsigned int k6 = ((16 * k)%prime);
		unsigned int max = opt_sieve_size - prime;

		for( unsigned int sieveIndex = ((k * (x + bp)) % prime);
			sieveIndex < max;
			sieveIndex += prime )
		{
			sieveReset(pSieve->pSieve, sieveIndex);
			sieveReset(pSieve->pSieve, sieveIndex + k2);
			sieveReset(pSieve->pSieve, sieveIndex + k3);
			sieveReset(pSieve->pSieve, sieveIndex + k4);
			sieveReset(pSieve->pSieve, sieveIndex + k5);
			sieveReset(pSieve->pSieve, sieveIndex + k6);
		}
        }
}

int getNext( uint32_t *pSieve, int _index )
{
        if( _index >= opt_sieve_size )
        {
            return -1;
        }
        while( 1 )
        {
            _index++;
            if( _index >= opt_sieve_size )
            {
		return -1;
            }
            if( sieveGet(pSieve, _index) )
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
	mpz_t nBits;
	double f;
	static const double l2 = 0.69314718056; // ln(2)

	SetCompact( nBits, compactBits );
	f = mpz_get_ui (nBits);
	mpz_clear(nBits);
	
	f = pow( f * l2, primes ) / hashrate;

	if( primes == 6 )
		return f * 0.05780811; // reciprocal of Hardy-Littlewood constant H6
/*	if( primes == 5 )
		return f * 0.09869924; // reciprocal of Hardy-Littlewood constant H5 (type b)
	if( primes == 4 )
		return f * 0.240895; // reciprocal of Hardy-Littlewood constant H4 (delta=4)
	if( primes == 3 )
		return f * 0.349864; // reciprocal of Hardy-Littlewood constant H3 (type a)
	if( primes == 2 )
		return f * 0.757392; // reciprocal of Hardy-Littlewood constant H2 ???
		*/
	return 0;
}

#define MR_TESTS 1


int scanhash_riecoin(int thr_id, uint32_t *pdata, const int primes,
	uint64_t max_nonce, uint64_t *hashes_done, uint32_t *pSieve)
{
	int primes_found;
	uint32_t hash[8] __attribute__((aligned(32)));
	uint64_t n = *(uint64_t *)(pdata + RIECOIN_DATA_NONCE);
	const uint64_t first_nonce = n;
	#ifdef REPORT_TESTS
		int tests;
	#endif
	mpz_t b;
	int sieveIndex;
	struct sieveData mySieve;
	mySieve.pSieve = pSieve;

	mpz_t bnTarget;
	mpz_init( bnTarget );

	if( max_nonce <= (opt_sieve_size + 1) * fixedMod )
	{
		max_nonce = 0;
	}
	else
	{
		max_nonce -= (opt_sieve_size + 1) * fixedMod;
	}

	mpz_init(b);
	
	uint64_t aux = *(uint64_t *)&pdata[RIECOIN_DATA_NONCE];
	*(uint32_t *)&pdata[RIECOIN_DATA_NONCE] = 0x80000000UL;
	*(((uint32_t *)&pdata[RIECOIN_DATA_NONCE]) + 1) = 0;
	sha256d_80_swap(hash, pdata);
	*(uint64_t *)&pdata[RIECOIN_DATA_NONCE] = aux;

	generatePrimeBase( b, hash, swab32(pdata[RIECOIN_DATA_DIFF]) );
	mpz_tdiv_q_2exp( b, b, 32 );
	mpz_add_ui( b, b, n >> 32 );
	mpz_mul_2exp( b, b, 32 );
	mpz_add_ui( b, b, n );

	do {
               	#ifdef REPORT_TESTS
			tests = 0;
		#endif
		init( &mySieve, b );
		sieveIndex = 0;

		while( (sieveIndex = getNext(mySieve.pSieve, sieveIndex) ) >= 0 )
                {
                	#ifdef REPORT_TESTS
                		tests++;
                	#endif
			mpz_set( bnTarget, b );
			mpz_add_ui( bnTarget, bnTarget, i2d(&mySieve, sieveIndex) );
                	if( mpz_probab_prime_p ( bnTarget, MR_TESTS) == 0 )
                	{
                		continue;
                	}
                	primes_found = 1;
			mpz_add_ui( bnTarget, bnTarget, 4 );
			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
			{
				primes_found++;
			}
			else
			{
				if( primes_found + 4 < primes )
				{
					continue;
				}
			}

			mpz_add_ui( bnTarget, bnTarget, 2 );
			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
			{
				primes_found++;
			}
			else
			{
				if( primes_found + 3 < primes )
				{
					continue;
				}
			}
			mpz_add_ui( bnTarget, bnTarget, 4 );
			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
			{
				primes_found++;
			}
			else
			{
				if( primes_found + 2 < primes )
				{
					continue;
				}
			}
			mpz_add_ui( bnTarget, bnTarget, 2 );
			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
			{
				primes_found++;
			}
			else
			{
				if( primes_found + 1 < primes )
				{
					continue;
				}
			}
			mpz_add_ui( bnTarget, bnTarget, 4 );
			if( mpz_probab_prime_p ( bnTarget, MR_TESTS) )
			{
				primes_found++;
			}
			if( primes_found >= primes )
			{
				*(uint64_t *)(pdata + RIECOIN_DATA_NONCE) = n + i2d(&mySieve, sieveIndex);
				pdata[RIECOIN_DATA_NONCE] = swab32(pdata[RIECOIN_DATA_NONCE]);
				pdata[RIECOIN_DATA_NONCE+1] = swab32(pdata[RIECOIN_DATA_NONCE+1]);
				*hashes_done = (n + i2d(&mySieve, sieveIndex) - first_nonce + 1) / efficiencyDivisor;
				mpz_clear(bnTarget);
				mpz_clear(b);
				return 1;
			}
		}
		#ifdef REPORT_TESTS
			applog(LOG_INFO, "thread %d tests: %d   n = %"PRIu64",  i2d = %u   max = %"PRIu64", idx = %d, x = %d ", thr_id, tests, n, i2d(&mySieve, opt_sieve_size), max_nonce, sieveIndex, mySieve.x );
		#endif
		
		n += i2d(&mySieve, opt_sieve_size);
		mpz_add_ui( b, b, i2d(&mySieve, opt_sieve_size) );
	} while (n < max_nonce && !work_restart[thr_id].restart);
	
	*hashes_done = (n - first_nonce + 1) / efficiencyDivisor;
	*(uint64_t *)(pdata + RIECOIN_DATA_NONCE) = max_nonce;
	mpz_clear (bnTarget);
	mpz_clear (b);
	return 0;
}
