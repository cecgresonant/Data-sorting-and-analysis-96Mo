#include "XIA_CFD.h"

#include <stdlib.h>
#include <stdint.h>

double XIA_CFD_Fraction_250MHz(uint16_t CFDvalue, char* fail)
{
	double correction;
	uint32_t cfdfailbit, cfdtrigsource, timecfd;

	cfdfailbit = ((CFDvalue & BIT15TO15) >> 15);
	cfdtrigsource = ((CFDvalue & BIT14TO14) >> 14);
	timecfd = ((CFDvalue & BIT13TO0) >> 0);
    
    if (cfdfailbit < 1) {
        correction = (((double)timecfd)/16384.0 - cfdtrigsource)*4.0;
        *fail = 0;
    } else {
        correction = 4.0*((double)rand()/RAND_MAX - 1);
        *fail = 1;
    }
    
	return correction;
}

double XIA_CFD_Fraction_500MHz(uint16_t CFDvalue, char* fail)
{
	double correction;
	uint32_t cfdtrigsource, timecfd;

	cfdtrigsource = ((CFDvalue & BIT15TO13) >> 13);
	timecfd = ((CFDvalue & BIT12TO0) >> 0);
    
    if (cfdtrigsource < 7) {
        correction = (((double)timecfd)/8192.0 + cfdtrigsource - 1.0)*2.0;
        *fail = 0;
    } else {
        correction = 2.0*(5.*(double)rand()/RAND_MAX - 1);
        *fail = 1;
    }
    
	return correction;
}
