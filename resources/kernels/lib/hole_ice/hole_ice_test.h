#ifndef HOLE_ICE_TEST_H
#define HOLE_ICE_TEST_H

typedef double floating_t;

struct floating4_t {
  floating_t x;
  floating_t y;
  floating_t z;
  floating_t w;
};

extern inline floating_t my_sqrt(floating_t);
extern inline floating_t sqr(floating_t);
extern inline floating_t my_nan();
extern inline bool my_is_nan(floating_t);
extern inline floating_t min(floating_t, floating_t);

#endif