#include<stdio.h>
#include <math.h>
// #include <string.h> // for memcpy

#include "fft.h"

void make_twiddle(my_complex *W)
{
   double delta = 2 * PI / N;
   int i;
   // team: divide and conquer 로 쓸려고하면 범위를 넘어가는 것 같다.
   for (i = 0; i < N; i++)
   {
      W[i].real = cos(i * delta);
      W[i].img = -sin(i * delta);
   }
}

int bit_revese(int x) // reordering function
{
   int a, b;
   int result = 0;
   for (int m = 0; m < M; m++)
   {
      a = (int)pow((double)2, (double)(m + 1));
      b = (int)pow((double)2, (double)m);
      result += (N / a) * ((x % a) / b);
   }
   return result;
}

my_complex mult(my_complex x, my_complex y) // multiply function of real and imaginary number
{
   my_complex result;
   result.real = x.real*y.real - x.img*y.img;
   result.img = x.real*y.img + x.img*y.real;
   return result;
}
my_complex add(my_complex x, my_complex y) // add function of real and imaginary number
{
   my_complex result;
   result.real = x.real + y.real;
   result.img = x.img + y.img;
   return result;
}
my_complex sub(my_complex x, my_complex y) // subtract function of real and imaginary number
{
   my_complex result;
   result.real = x.real - y.real;
   result.img = x.img - y.img;
   return result;
}

void dft_ref(my_complex* out, my_complex* input, my_complex* W_in)
{
   my_complex temp[N];

   int m, n;
   int k;

   /* comment out for fast testing */
   for(k = 0; k < N; k++)
   {
      temp[k].real = 0;
      temp[k].img = 0;
      for(n = 0; n < N; n++)
      {
         temp[k].real += input[n].real * ( cos((2 * PI / N) * k * n)) + input[n].img * ( sin((2 * PI / N) * k * n));
         temp[k].img += input[n].real * ( -sin((2 * PI / N) * k * n)) + input[n].img * ( cos((2 * PI / N) * k * n));
      }
   }


   for (n = 0; n < N; n++)
   {
      out[n] = temp[n];
   }


   return;
}



unsigned char reverse_table_2bits[64] ={
		0x0,0x10,0x20,0x30,0x4,0x14,0x24,0x34,
		0x8,0x18,0x28,0x38,0xc,0x1c,0x2c,0x3c,
		0x1,0x11,0x21,0x31,0x5,0x15,0x25,0x35,
		0x9,0x19,0x29,0x39,0xd,0x1d,0x2d,0x3d,
		0x2,0x12,0x22,0x32,0x6,0x16,0x26,0x36,
		0xa,0x1a,0x2a,0x3a,0xe,0x1e,0x2e,0x3e,
		0x3,0x13,0x23,0x33,0x7,0x17,0x27,0x37,
		0xb,0x1b,0x2b,0x3b,0xf,0x1f,0x2f,0x3f,
};


void fft_opt(my_complex* out, my_complex* input, my_complex* W_in)
{

    int
           iCnt1, iCnt2, iCnt3,
           iL,    iM,    iQ,
           iA,    iB,    iC,     iD;
       float
           fRealA,   fRealB,    fRealC,    fRealD,
           fImagA,   fImagB,    fImagC,    fImagD,
           fReal_Wq, fReal_W2q, fReal_W3q,
           fImag_Wq, fImag_W2q, fImag_W3q;

       int M_4 = (int)(log(N)/log(4));

       iL = 1;
       iM = N / 4;

       iL = 1;
       iM = N / 4;

       for (iCnt1 = 0; iCnt1 < M_4; ++iCnt1)
       {
           iQ = 0;
           for (iCnt2 = 0; iCnt2 < iM; ++iCnt2)
           {
               iA = iCnt2;
               fReal_Wq  = W_in[    iQ].real;
               fImag_Wq  = W_in[    iQ].img;
               fReal_W2q = W_in[2 * iQ].real;
               fImag_W2q = W_in[2 * iQ].img;
               fReal_W3q = W_in[3 * iQ].real;
               fImag_W3q = W_in[3 * iQ].img;

               for (iCnt3 = 0; iCnt3 < iL; ++iCnt3)
               {
                   iB = iA +     iM;
                   iC = iA + 2 * iM;
                   iD = iA + 3 * iM;

                   /* Butterfly: 42 FOP, 12 FMUL, 30 FADD */

                   fRealA = input[iA].real + input[iB].real
                          + input[iC].real + input[iD].real;
                   fImagA = input[iA].img + input[iB].img
                          + input[iC].img + input[iD].img;
                   fRealB = input[iA].real + input[iB].img
                          - input[iC].real - input[iD].img;
                   fImagB = input[iA].img - input[iB].real
                          - input[iC].img + input[iD].real;
                   fRealC = input[iA].real - input[iB].real
                          + input[iC].real - input[iD].real;
                   fImagC = input[iA].img - input[iB].img
                          + input[iC].img - input[iD].img;
                   fRealD = input[iA].real - input[iB].img
                          - input[iC].real + input[iD].img;
                   fImagD = input[iA].img + input[iB].real
                          - input[iC].img - input[iD].real;

                   input[iA].real = fRealA;
                   input[iA].img = fImagA;
                   input[iB].real = fRealB * fReal_Wq  - fImagB * fImag_Wq;
                   input[iB].img = fRealB * fImag_Wq  + fImagB * fReal_Wq;
                   input[iC].real = fRealC * fReal_W2q - fImagC * fImag_W2q;
                   input[iC].img = fRealC * fImag_W2q + fImagC * fReal_W2q;
                   input[iD].real = fRealD * fReal_W3q - fImagD * fImag_W3q;
                   input[iD].img = fRealD * fImag_W3q + fImagD * fReal_W3q;

                   iA = iA + 4 * iM;
               }
               iQ += iL;
           }
           iL *= 4;
           iM /= 4;
       }


   int n;
   /*
   for(n=0; n < N;n++)
      out[n] = input[reverse_table_radix2[n >> 6] | reverse_table_radix2[n & 0b111111] << 6];
      */
   for(n=0; n < N;n++){
      out[n] = input[reverse_table_2bits[n >> 6] | reverse_table_2bits[n & 0b111111] << 6];
   }




      return ;
}
