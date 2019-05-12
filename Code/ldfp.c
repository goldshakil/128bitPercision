#include <stdio.h>
#include <stdlib.h>
#include "ldfp.h"

typedef unsigned char* cast;
typedef unsigned char byte1;

struct Flags {
	union {
		struct {
			unsigned short a : 4;
			unsigned short b : 4;
			unsigned short c : 3;
			unsigned short d : 4;
			unsigned short : 1; //sign bit, but we don't want to access, so there is no name.
		};
		unsigned short e;
		unsigned char exp[2];
	};
};

union FRAC {
	unsigned long long Frac;
	unsigned char data[8];
};

int TwoOf(int n)
{
	if (n == 0)
		return 1;

	return TwoOf(n - 1) * 2;
}

long_double double_to_long_double(double op) //conversion
{
	long_double ret = { { 0x0 }, };
	ret.IsZero = 0; 

	cast Casted = (cast)&op;

	if (op == 0)
	{
		ret.IsZero = 1;
		return ret;
	}
	//Exp bit
	struct Flags Exp = { 0, }; // ALL -> zero - initialized 
	Exp.c = Casted[7] >> 4;
	Exp.b = Casted[7];
	Exp.a = Casted[6] >> 4;

	Exp.e = Exp.e - 1024 + 16384; // 64bit -> 128bit

	ret.data[15] = Exp.exp[1];
	ret.data[14] = Exp.exp[0];

	//sign bit
	byte1 finder_sign = 0x80;
	byte1 ucZero = 0x00;
	if ((byte1)(Casted[7] & finder_sign) != ucZero) // if the sign bit is 1 
		ret.data[15] = (byte1)(ret.data[15] | finder_sign); // ret.data[15] -> 10000000

															//frac bit
	for (int i = 13; i > 7; i--)
		ret.data[i] = (byte1)(Casted[i - 7] << 4) + (byte1)(Casted[i - 8] >> 4);

	ret.data[7] = (byte1)(Casted[0] << 4);

	return ret;
}

long_double FP_mul(long_double op1, long_double op2) {

	if (op1.IsZero)
		return op1;

	else if (op2.IsZero)
		return op2;

	/*Union Representing the Sign and Exponent bits*/
	typedef union exp_rep_128
	{
		unsigned char exp_byte[2];
		struct
		{
			unsigned short int exponent : 15;
			unsigned short int sign : 1;
		} field;
	}exp_rep_128;

	exp_rep_128 OP1, OP2, MUL;   /*MUL -> Multiplcation Answer*/

								 /*Extracting the E value from Both op1 and op2*/
	OP1.exp_byte[1] = op1.data[15];
	OP1.exp_byte[0] = op1.data[14];
	int E_op1 = OP1.field.exponent - 16383;

	OP2.exp_byte[1] = op2.data[15];
	OP2.exp_byte[0] = op2.data[14];
	int E_op2 = OP2.field.exponent - 16383;


	/*Calculating E and EXP for Multiplication Result*/
	int E_mul;
	if (E_op1 + E_op2 > 16384)   //Largest possible exp=32,767 / Largest possible E=16,384
		E_mul = 16384;

	else if (E_op1 + E_op2 < -16383)
	{
		op1.IsZero = 1;
		for (int i = 0; i < 16; i++)
			op1.data[i] = 0;

		return op1;
	}
	else
		E_mul = E_op2 + E_op1;  //in multiplcation we add the powers
								/*setting the Sign bit for the Multiplication Result */
	if (OP1.field.sign ^ OP2.field.sign)  // if any operand is negative
		MUL.field.sign = 1;

	else  //Both are positive
		MUL.field.sign = 0;

	MUL.field.exponent = E_mul + 16383;//setting MUL exp value

	long_double op = { { 0x0 }, };

	int mul[14] = { 0, };
	for (int i = 0; i<14; i++)
		mul[i] = 0;

	for (int i = 0; i<16; i++)
		op.data[i] = 0;

	for (int i = 7; i < 14; i++)
	{
		for (int j = 7; j < 14; j++)
		{
			int num = ((int)op1.data[j])*((int)op2.data[i]);
			mul[i + j - 14] += num;
		}
		mul[i] += (int)op2.data[i];
	}
	for (int i = 13; i > 6; i--) {
		mul[i] += (int)op1.data[i];
	}

	int q = 0;

	for (int i = 0; i < 14; i++)
	{
		q = mul[i] / 256;
		mul[i] %= 256;

		if (i < 13)
			mul[i + 1] += q;
	}

	q++;

	if (q / 2 == 1)
	{
		MUL.field.exponent++;
		for (int i = 0; i < 14; i++)
		{
			if (i != 13)
				op.data[i] = (byte1)((mul[i] >> 1) + (mul[i + 1] << 7));

			else
				op.data[i] = (byte1)((mul[i] >> 1) + (q % 2) * 128);
		}
		q /= 2;
	}

	else
	{
		for (int i = 0; i < 14; i++)
			op.data[i] = mul[i];
	}

	op.data[15] = MUL.exp_byte[1];
	op.data[14] = MUL.exp_byte[0];
	return op;
}

long_double FP_add(long_double op1, long_double op2)
{
	if (op1.IsZero == 1)
		return op2;

	if (op2.IsZero == 1)
		return op1;

	struct Flags Exp = { 0, };
	Exp.exp[1] = (byte1)(op1.data[15] & (byte1)(0x7f)); //We must make sign bit 0!
	Exp.exp[0] = op1.data[14];
	// Now, Exp1.e is the EXP of op1!
	signed short e1 = Exp.e - 16383; // E = EXP - bias(2^14 - 1)

	Exp.exp[1] = 0; Exp.exp[0] = 0;
	Exp.exp[1] = (byte1)(op2.data[15] & (byte1)(0x7f)); //We must make sign bit 0!
	Exp.exp[0] = op2.data[14];
	// Now, Exp2.e is the EXP of op2!
	signed short e2 = Exp.e - 16383; // E = EXP - bias(2^14 - 1)
									 //unsigned char mantissa[14] = { 0, };

	double m1 = 1, m2 = 1; //mantissa is 1.xxxxx! So, they're going to be added!
	double divider = 256;

	for (int i = 13; i > 6; i--) // find m1!
	{
		m1 += (double)((double)op1.data[i] / divider);
		divider *= 256;
	}

	divider = 256;

	for (int i = 13; i > 6; i--) // find m2!
	{
		m2 += (double)((double)op2.data[i] / divider);
		divider *= 256;
	}

	byte1 finder_sign = 0x80;
	byte1 ucZero = 0x00;
	double sign1 = 1, sign2 = 1;

	if ((byte1)(op1.data[15] & finder_sign) != ucZero) // if the sign bit is 1 
		sign1 = -1;

	if ((byte1)(op2.data[15] & finder_sign) != ucZero) // if the sign bit is 1 
		sign2 = -1;


	if (e1 == e2)
		m1 = sign1*m1 + sign2*m2;

	else if (e1 > e2) // e2-> e1! Instead, mantissa2 is divide by 2^e1-e2
	{
		if (e1 - e2 < 67)
			m1 = sign1*m1 + sign2*m2 / (double)TwoOf(e1 - e2);

		else
			m1 = sign1*m1;
	}
	else // e1 < e2 
	{
		if (e2 - e1 < 67)
			m1 = sign1*m1 / (double)TwoOf(e2 - e1) + sign2*m2;

		else
			m1 = sign2*m2;

		e1 = e2; // the final return value is op1 not op2! So, we change the value about op1!
	}

	if (m1 == (double)0)  //IF mantissa1 is zero -> zero
	{
		for (int i = 0; i < 16; i++)
			op1.data[i] = 0;
		op1.IsZero = 1;
		return op1;
	}

	//Sign bit
	if (m1 < 0)
	{
		sign1 = 1; //Now, we should consider sign1 as a sign bit! Ignore access to the value!
		m1 = -m1;
	}
	else
		sign1 = 0;

	if (m1 >= 2)
	{
		m1 /= 2;
		e1 += 1;
	}

	else if (m1 < 1)
	{
		while (m1 < 1)
		{
			m1 *= 2;
			e1 -= 1;
		}
	}

	//Exp sign bit. 
	Exp.e = (unsigned short)(e1 + 16383);

	if (sign1 == 1)
		op1.data[15] = (byte1)(Exp.exp[1] ^ finder_sign);
	else
		op1.data[15] = (byte1)(Exp.exp[1] & (byte1)(0x7f));

	op1.data[14] = Exp.exp[0];

	//Frac bit
	double frac = m1 - 1;
	union FRAC saver = { 0, }; // this is the thing that saves frac bits!

	int num = 0;
	while (frac != (double)0)
	{
		frac *= 2;

		if (frac >= 1) //this is 1
		{
			saver.Frac += 1;
			frac -= 1;
		}
		saver.Frac = saver.Frac << 1;

		num++;
	}

	while (num < 63) // we do shift -right until num is 64.
	{
		saver.Frac = saver.Frac << 1;
		num++;
	}

	for (int i = 13; i > 5; i--)
		op1.data[i] = saver.data[i - 6];

	return op1;
}

char *long_double_print_bitseq(long_double op)
{
	char *str = (char*)calloc(MAXLEN, sizeof(char)); // this memory needs to be freed from caller 

	byte1 byte;
	int k = 0;

	for (int i = 15; i >= 0; i--)
	{
		for (int j = 7; j >= 0; j--)
		{
			byte = (byte1)((byte1)(op.data[i] >> j) & (byte1)1);
			*(str + k) = byte + '0';
			k++;
		}
		str[k++] = ' ';
	}
	str[k++] = '\0';
	return str;
}

char *long_double_print_normalized(long_double op)
{
	char *str = (char*)calloc(MAXLEN, sizeof(char)); // this memory needs to be freed from caller 

	int index = 0;
	byte1 finder_sign = 0x80;
	byte1 ucZero = 0x00;

	// 1(sign)+ 114(1. +frac) + 1(blank) + 1(x) + 1(blank) + 2(2^) + 1(sign) + 5(E) + 1(NULL); 

	//mantissa part
	if (op.IsZero == 1) //   IF op == 0
	{
		for (int i = 0; i < 114; i++)
			str[i] = '0';
		str[1] = '.';
		index = 114;
	}

	else
	{
		//sign part
		if ((byte1)(op.data[15] & finder_sign) != ucZero) // if the sign bit is 1 -> negative
			str[index++] = '-';

		str[index++] = '1'; str[index++] = '.';

		byte1 byte;

		for (int i = 13; i >= 0; i--)
		{
			for (int j = 7; j >= 0; j--)
			{
				byte = (byte1)((byte1)(op.data[i] >> j) & (byte1)1);
				*(str + index) = byte + '0';
				index++;
			}
		}
	}
	//    3 blocks for 1(blank) + 1(x) + 1(blank)  + 2(2^)
	str[index++] = ' '; str[index++] = 'x'; str[index++] = ' '; str[index++] = '2'; str[index++] = '^';

	// 1(sign) + 5(E) + 1(NULL)
	if (op.IsZero == 1)
	{
		str[index++] = '-'; str[index++] = '1'; str[index++] = '6';
		str[index++] = '3'; str[index++] = '8'; str[index++] = '2';
		str[index++] = '\0';
	}
	else
	{
		struct Flags Exp = { 0, };
		Exp.exp[1] = (byte1)(op.data[15] & (byte1)(0x7f));
		Exp.exp[0] = op.data[14];
		int E = Exp.e - 16383; // E = EXP - bias(2^14 - 1)

		if (E == 0.000)
		{
			str[index++] = '0'; str[index++] = '\0';
			return str;
		}

		else if (E < 0)
		{
			str[index++] = '-';
			E -= 1;
		}

		int number = 0;
		int a = E;

		while (a != 0)
		{
			a /= 10;
			number++;
		}

		for (int i = number; i >0; i--)
		{
			str[index + i - 1] = '0' + (E % 10);
			E /= 10;
		}

		str[index + number] = '\0';
	}
	return str;
}
