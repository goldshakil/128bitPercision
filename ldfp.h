#include <stdio.h>

#define MAXLEN 300

typedef struct {
	unsigned char data[16];
	int IsZero;
} long_double;

long_double double_to_long_double(double op);
long_double FP_mul(long_double op1, long_double op2);
long_double FP_add(long_double op1, long_double op2);
char *long_double_print_bitseq(long_double op);
char *long_double_print_normalized(long_double op);