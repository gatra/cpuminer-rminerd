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

int initPrimeTable( void )
{
}

double riecoin_time_to_block( double hashrate, uint32_t compactBits, int primes )
{
	return 1;
}

int scanhash_riecoin(int thr_id, uint32_t *pdata, const int primes,
	uint32_t max_nonce, unsigned long *hashes_done, unsigned int *pSieve)
{
	return 0;
}
