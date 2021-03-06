#ifndef RDTSC_H
#define RDTSC_H 

static __inline__ unsigned long long rdtsc (void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc":"=a" (lo), "=d" (hi));
  return ((unsigned long long) lo) | (((unsigned long long) hi) << 32);
}

#endif