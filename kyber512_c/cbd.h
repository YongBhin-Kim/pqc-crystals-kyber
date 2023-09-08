#ifndef CBD_H
#define CBD_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

void cbd_eta1(poly *res, const uint8_t buf[ETA1 * N/4]);
void cbd_eta2(poly *res, const uint8_t buf[ETA2 * N/4]);

#endif
