#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

int16_t csubq(int16_t x);
int16_t simple_reduce(int16_t x);
int16_t simple_reduce32(int32_t x);
int16_t barrett_reduce(int16_t a);

#endif