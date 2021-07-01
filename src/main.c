#include <stdlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fft.h"
#include <xtime_l.h>

#include "parameter.h"
#include "benchmarking.h"

#define TEST_CASE 2
#define TEST_ROUNDS 5

float W_I[N], W_Q[N];
float input_I[N], input_Q[N];
float output_I[N], output_Q[N];

XTime xstart, xstop;
float func_time;

unsigned int initializor_dummy(my_complex* uiParam0, my_complex* uiParam1, my_complex* uiParam2){return 1;}
unsigned int validator_dummy(my_complex* uiParam0, my_complex* uiParam1, my_complex* uiParam2){return 1;}

void make_input(my_complex *input)
{
	int i, j, k;
	for (i = 0; i < NOFDM; i++)
	{
		for (j = 0; j < NRX; j++)
		{
			for(k = 0; k < N; k++)
			{
				input[k].real = rand()%9;
				input[k].img = rand()%9;
			}
		}
	}
}

int main()
{

	xil_printf("\r\n");
	xil_printf("<%d-point Fourier Transform>\r\n", N);

	int i, j, mode_sel, ii = 0;

	BENCHMARK_CASE *pBenchmarkCase;
	BENCHMARK_STATISTICS *pStat;

	my_complex *in_FFT = 0;
	my_complex *out_FFT_ref = 0;
	my_complex *out_FFT_opt = 0;
	my_complex *W = 0;

	in_FFT = (my_complex*)malloc(sizeof(my_complex) * NOFDM * NRX * NSC); if (in_FFT == NULL) { xil_printf("Memory allocation error\r\n"); };
	out_FFT_ref = (my_complex*)malloc(sizeof(my_complex) * NOFDM * NRX * NSC); if (out_FFT_ref == NULL) { xil_printf("Memory allocation error\r\n"); };
	out_FFT_opt = (my_complex*)malloc(sizeof(my_complex) * NOFDM * NRX * NSC); if (out_FFT_opt == NULL) { xil_printf("Memory allocation error\r\n"); };
	W = (my_complex*)malloc(sizeof(my_complex) * N); if (W == NULL) { xil_printf("Memory allocation error\r\n"); };

	double error = 0;
	double signal = 0;
	double NSRdB = 0;
	double opt_time = 0;
	double ref_time = 0;

    make_twiddle(W);
    make_input(in_FFT);

	for (mode_sel = 0; mode_sel < _mode; mode_sel++)
	{
		if (mode_sel == 0)
		{
			for (i = 0; i < NOFDM; i++)
			{
				for (j = 0; j < NRX; j++)
				{
					dft_ref(&out_FFT_ref[i*NRX*NSC + j * NSC], &in_FFT[i*NRX*NSC + j * NSC], W);
				}
			}
		}
		else
		{
			for (i = 0; i < NOFDM; i++)
			{
				for (j = 0; j < NRX; j++)
				{
					fft_opt(&out_FFT_opt[i*NRX*NSC + j * NSC], &in_FFT[i*NRX*NSC + j * NSC], W);
				}
			}
		}

	}

	for (ii = 0; ii < NOFDM*NTX*NSC; ii++) {
		error += pow((out_FFT_ref[ii].real - out_FFT_opt[ii].real), 2) + pow((out_FFT_ref[ii].img - out_FFT_opt[ii].img), 2);;
		signal += pow((out_FFT_ref[ii].real), 2) + pow((out_FFT_ref[ii].img), 2);
	}

	NSRdB = 10 * log10(error / signal);
	printf("\nMeasured Accuracy: NSR(dB) = %0.3f \n", NSRdB);


	BENCHMARK_CASE BenchmarkCases[TEST_CASE] = {
		{"DFT Reference", TEST_ROUNDS, initializor_dummy, dft_ref,
			{out_FFT_ref, in_FFT, W},
			0, validator_dummy
		},
		{"FFT Optimization", TEST_ROUNDS, initializor_dummy, fft_opt,
				{out_FFT_opt, in_FFT, W},
			0, validator_dummy
		}
	};


	Xil_L1DCacheEnable();
	Xil_L2CacheDisable();
	Xil_L1ICacheEnable();

	xil_printf("\r\n");
	xil_printf("-----Benchmarking Start-----\r\n");
	for (i = 0; i < TEST_CASE; i++)
	{
		pBenchmarkCase = &BenchmarkCases[i];
		pStat = &(pBenchmarkCase->stat);
		printf("Case %d: %s\r\n", i, pBenchmarkCase->pName);
		run_benchmark_single(pBenchmarkCase);
		statistics_print(pStat);
		if (i == 0)
			ref_time = pStat->ullMax;
		else
			opt_time = pStat->ullMax;
	}

	xil_printf("----Benchmarking Complete----\r\n");
	xil_printf("\r\n");
	printf("Optimized FFT is x%.2f faster than Reference\r\n", (double)(ref_time/opt_time));
	xil_printf("\r\n");

	free(in_FFT);
	free(out_FFT_ref);
	free(out_FFT_opt);
	free(W);


    return 0;
}

