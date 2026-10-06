/*
 * osmo-fl2k, turns FL2000-based USB 3.0 to VGA adapters into
 * low cost DACs
 *
 * I/Q interpolating upconverter for fl2k_iq
 *
 * Copyright (C) 2026 by Krisztián Szöllősi <szollosi.krisztian.robert@gmail.com>
 *
 * SPDX-License-Identifier: GPL-2.0+
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "iq_upconv.h"

#ifndef M_PI
# define M_PI		3.14159265358979323846
#endif

/*
 * FIR stage: Kaiser-windowed sinc, cut-off at in_rate / 2, transition
 * band 0.4 ... 0.6 * in_rate, about 80 dB stop band attenuation.
 * The number of taps per polyphase branch does not depend on the
 * interpolation factor (the transition band scales with it).
 */
#define MAX_INTERPOLATION	64
#define TAPS_PER_PHASE		26
#define KAISER_BETA		7.857

/*
 * Linear interpolation between the f1 samples: its images around
 * multiples of f1 are attenuated by sinc^2, i.e. by about
 * 40 * log10(f1 / offset) dB. With f1 = 64 * 48 kHz and a signal edge at
 * 19.2 kHz this is 88 dB, for a 10 kHz wide DRM signal (5 kHz) 111 dB.
 */

#define FIFO_SAMPLES		(1 << 16)

/* The NCO phasor is recomputed from the exact phase every NCO_BLOCK samples. */
#define NCO_BLOCK		1024

#define TWO_POW_32		4294967296.0
#define TWO_POW_64		18446744073709551616.0

struct iq_upconv {
	struct iq_upconv_cfg cfg;

	/* FIR interpolator, polyphase layout: coef[p * TAPS_PER_PHASE + k] = h[k * interp + p] */
	int interp;
	float *coef;

	/* input FIFO (complex, full scale 1.0), keeps the filter history */
	float *fifo_i;
	float *fifo_q;
	uint64_t fifo_base;	/* absolute input index of fifo[0] */
	size_t fifo_len;
	uint64_t in_total;	/* samples in the FIFO so far, incl. flush padding */
	uint64_t in_real;	/* samples pushed by the caller */
	int flushed;

	/* fractional interpolator: position = m + frac / 2^32 in f1 samples */
	uint64_t m;
	uint64_t frac;
	uint64_t step;		/* f1 / out_rate, 32.32 fixed point */
	int pair_valid;
	float y0i, y0q, y1i, y1q;

	/* NCO: phase in cycles, 0.64 fixed point */
	uint64_t nco_step;
	double zr, zi;		/* exp(j * phase) of the current sample */
	double wr, wi;		/* exp(j * phase step) */

	/* output */
	double scale;
	uint32_t rng;
	uint64_t out_count;
	uint64_t clipped;

	/* input statistics */
	double power_sum;
	double peak;
};

/* Modified Bessel function of the first kind, order 0 (for the Kaiser window). */
static double bessel_i0(double x)
{
	double sum = 1.0, term = 1.0, half = x / 2.0;
	int k;

	for (k = 1; k < 200; k++) {
		term *= (half / k) * (half / k);
		sum += term;
		if (term < sum * 1e-17)
			break;
	}

	return sum;
}

static int choose_interpolation(const struct iq_upconv_cfg *cfg)
{
	int r = MAX_INTERPOLATION;

	/* f1 must stay below out_rate / 2, otherwise the linear interpolator aliases */
	while (r > 1 && r * cfg->in_rate > cfg->out_rate / 2)
		r /= 2;

	return r;
}

const char *iq_upconv_check(const struct iq_upconv_cfg *cfg)
{
	if (!cfg)
		return "no configuration";
	if (!isfinite(cfg->in_rate) || cfg->in_rate <= 0)
		return "input sample rate must be positive";
	if (!isfinite(cfg->out_rate) || cfg->out_rate <= 0)
		return "output sample rate must be positive";
	if (cfg->in_rate > cfg->out_rate / 2)
		return "input sample rate must be at most half the output sample rate";
	if (!isfinite(cfg->carrier) || cfg->carrier < 0 || cfg->carrier >= cfg->out_rate / 2)
		return "carrier must be between 0 and half the output sample rate";
	if (!isfinite(cfg->gain_db) || cfg->gain_db > 60 || cfg->gain_db < -60)
		return "gain must be between -60 and +60 dB";

	return NULL;
}

static int design_filter(iq_upconv_t *u)
{
	int n = u->interp * TAPS_PER_PHASE;
	double center = (n - 1) / 2.0;
	double norm = bessel_i0(KAISER_BETA);
	double *h, sum = 0.0;
	int i, p, k;

	h = malloc(n * sizeof(double));
	u->coef = malloc(n * sizeof(float));
	if (!h || !u->coef) {
		free(h);
		return -1;
	}

	for (i = 0; i < n; i++) {
		double x = (i - center) / u->interp;
		double sinc = (x == 0.0) ? 1.0 : sin(M_PI * x) / (M_PI * x);
		double r = (n > 1) ? 2.0 * i / (n - 1) - 1.0 : 0.0;
		double window = bessel_i0(KAISER_BETA * sqrt(1.0 - r * r)) / norm;

		h[i] = sinc * window;
		sum += h[i];
	}

	/* unity DC gain after interpolation: the taps sum to interp */
	for (p = 0; p < u->interp; p++)
		for (k = 0; k < TAPS_PER_PHASE; k++)
			u->coef[p * TAPS_PER_PHASE + k] = (float)(h[k * u->interp + p] * u->interp / sum);

	free(h);
	return 0;
}

iq_upconv_t *iq_upconv_create(const struct iq_upconv_cfg *cfg)
{
	iq_upconv_t *u;
	double f1;

	if (iq_upconv_check(cfg))
		return NULL;

	u = calloc(1, sizeof(*u));
	if (!u)
		return NULL;

	u->cfg = *cfg;
	u->interp = choose_interpolation(cfg);
	u->fifo_i = malloc(FIFO_SAMPLES * sizeof(float));
	u->fifo_q = malloc(FIFO_SAMPLES * sizeof(float));

	if (!u->fifo_i || !u->fifo_q || design_filter(u) < 0) {
		iq_upconv_destroy(u);
		return NULL;
	}

	f1 = cfg->in_rate * u->interp;
	u->step = (uint64_t)(f1 / cfg->out_rate * TWO_POW_32 + 0.5);
	u->nco_step = (uint64_t)(cfg->carrier / cfg->out_rate * TWO_POW_64);
	u->wr = cos(2.0 * M_PI * cfg->carrier / cfg->out_rate);
	u->wi = sin(2.0 * M_PI * cfg->carrier / cfg->out_rate);
	u->scale = 127.0 * pow(10.0, cfg->gain_db / 20.0);
	u->rng = 0x12345678u;

	return u;
}

void iq_upconv_destroy(iq_upconv_t *u)
{
	if (!u)
		return;

	free(u->coef);
	free(u->fifo_i);
	free(u->fifo_q);
	free(u);
}

int iq_upconv_interpolation(const iq_upconv_t *u)
{
	return u ? u->interp : 0;
}

int iq_format_from_name(const char *name, enum iq_format *fmt)
{
	if (!name || !fmt)
		return -1;

	if (!strcmp(name, "cs8"))
		*fmt = IQ_FMT_CS8;
	else if (!strcmp(name, "cs16"))
		*fmt = IQ_FMT_CS16;
	else if (!strcmp(name, "cf32"))
		*fmt = IQ_FMT_CF32;
	else
		return -1;

	return 0;
}

size_t iq_format_sample_size(enum iq_format fmt)
{
	switch (fmt) {
	case IQ_FMT_CS8:
		return 2;
	case IQ_FMT_CS16:
		return 4;
	case IQ_FMT_CF32:
		return 8;
	}

	return 0;
}

/* Drops the samples the filter will never need again. */
static void fifo_compact(iq_upconv_t *u)
{
	uint64_t oldest = u->m / u->interp;
	uint64_t keep_from, drop;

	keep_from = oldest >= TAPS_PER_PHASE - 1 ? oldest - (TAPS_PER_PHASE - 1) : 0;
	if (keep_from <= u->fifo_base)
		return;

	drop = keep_from - u->fifo_base;
	if (drop > u->fifo_len)
		drop = u->fifo_len;

	memmove(u->fifo_i, u->fifo_i + drop, (u->fifo_len - drop) * sizeof(float));
	memmove(u->fifo_q, u->fifo_q + drop, (u->fifo_len - drop) * sizeof(float));
	u->fifo_base += drop;
	u->fifo_len -= drop;
}

static float read_f32le(const uint8_t *b)
{
	uint32_t v = (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
		     ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
	float f;

	memcpy(&f, &v, sizeof(f));
	return f;
}

static void decode_sample(enum iq_format fmt, const uint8_t *b, float *i, float *q)
{
	switch (fmt) {
	case IQ_FMT_CS8:
		*i = (int8_t)b[0] / 127.0f;
		*q = (int8_t)b[1] / 127.0f;
		break;
	case IQ_FMT_CS16:
		*i = (int16_t)(b[0] | (b[1] << 8)) / 32767.0f;
		*q = (int16_t)(b[2] | (b[3] << 8)) / 32767.0f;
		break;
	case IQ_FMT_CF32:
		*i = read_f32le(b);
		*q = read_f32le(b + 4);
		break;
	}
}

size_t iq_upconv_push(iq_upconv_t *u, enum iq_format fmt, const void *buf, size_t len)
{
	const uint8_t *b = buf;
	size_t size = iq_format_sample_size(fmt);
	size_t count, n;

	if (!u || !buf || !size || u->flushed)
		return 0;

	fifo_compact(u);
	count = len / size;
	if (count > FIFO_SAMPLES - u->fifo_len)
		count = FIFO_SAMPLES - u->fifo_len;

	for (n = 0; n < count; n++) {
		float i, q, a;

		decode_sample(fmt, b + n * size, &i, &q);
		if (u->cfg.swap_iq) {
			float t = i;
			i = q;
			q = t;
		}

		u->power_sum += (double)i * i + (double)q * q;
		a = fabsf(i) > fabsf(q) ? fabsf(i) : fabsf(q);
		if (a > u->peak)
			u->peak = a;

		u->fifo_i[u->fifo_len] = i;
		u->fifo_q[u->fifo_len] = q;
		u->fifo_len++;
	}

	u->in_total += count;
	u->in_real += count;
	return count * size;
}

void iq_upconv_flush(iq_upconv_t *u)
{
	size_t n;

	if (!u || u->flushed)
		return;

	/* The last input sample reaches the output TAPS_PER_PHASE input samples later. */
	fifo_compact(u);
	for (n = 0; n < TAPS_PER_PHASE + 1 && u->fifo_len < FIFO_SAMPLES; n++) {
		u->fifo_i[u->fifo_len] = 0.0f;
		u->fifo_q[u->fifo_len] = 0.0f;
		u->fifo_len++;
		u->in_total++;
	}
	u->flushed = 1;
}

/* Is the f1 sample with index mm computable from the input so far? */
static inline int available(const iq_upconv_t *u, uint64_t mm)
{
	return mm / u->interp < u->in_total;
}

/* The f1 sample with index mm: one polyphase branch of the FIR filter. */
static inline void fir_value(const iq_upconv_t *u, uint64_t mm, float *oi, float *oq)
{
	uint64_t n = mm / u->interp;
	const float *c = u->coef + (mm % u->interp) * TAPS_PER_PHASE;
	float si = 0.0f, sq = 0.0f;
	int k;

	for (k = 0; k < TAPS_PER_PHASE; k++) {
		uint64_t idx;
		size_t j;

		if (n < (uint64_t)k)
			break;		/* before the first input sample: zero */
		idx = n - k;
		if (idx < u->fifo_base)
			break;		/* never happens: the history is kept */
		j = idx - u->fifo_base;
		si += c[k] * u->fifo_i[j];
		sq += c[k] * u->fifo_q[j];
	}

	*oi = si;
	*oq = sq;
}

/* Recomputes the NCO phasor from the exact phase of the current output sample. */
static void nco_sync(iq_upconv_t *u)
{
	uint64_t phase = u->nco_step * u->out_count;	/* wraps modulo one cycle */
	double angle = (double)phase * (2.0 * M_PI / TWO_POW_64);

	u->zr = cos(angle);
	u->zi = sin(angle);
}

static inline uint32_t xorshift32(uint32_t *state)
{
	uint32_t x = *state;

	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*state = x;
	return x;
}

/* Uniform in [-0.5, 0.5) */
static inline double uniform(uint32_t *state)
{
	return (xorshift32(state) >> 8) * (1.0 / 16777216.0) - 0.5;
}

static size_t produce(iq_upconv_t *u, int8_t *out8, float *outf, size_t n)
{
	size_t i;

	for (i = 0; i < n; i++) {
		float t, vi, vq;
		double y, nzr;
		uint64_t adv;

		if (!u->pair_valid) {
			if (!available(u, u->m + 1))
				break;
			fir_value(u, u->m, &u->y0i, &u->y0q);
			fir_value(u, u->m + 1, &u->y1i, &u->y1q);
			u->pair_valid = 1;
		}

		if ((u->out_count & (NCO_BLOCK - 1)) == 0)
			nco_sync(u);

		/* linear interpolation between the two f1 samples */
		t = (float)(u->frac * (1.0 / TWO_POW_32));
		vi = u->y0i + (u->y1i - u->y0i) * t;
		vq = u->y0q + (u->y1q - u->y0q) * t;

		/* Re{(I + jQ) * exp(j * phase)} */
		y = (vi * u->zr - vq * u->zi) * u->scale;
		nzr = u->zr * u->wr - u->zi * u->wi;
		u->zi = u->zr * u->wi + u->zi * u->wr;
		u->zr = nzr;

		if (outf) {
			outf[i] = (float)y;
		} else {
			int q;

			if (u->cfg.dither)
				y += uniform(&u->rng) + uniform(&u->rng);

			/* round half away from zero, clip symmetrically to +-127 */
			q = y >= 0 ? (int)(y + 0.5) : -(int)(-y + 0.5);
			if (q > 127) {
				q = 127;
				u->clipped++;
			} else if (q < -127) {
				q = -127;
				u->clipped++;
			}
			out8[i] = (int8_t)q;
		}
		u->out_count++;

		/* advance the fractional position */
		u->frac += u->step;
		adv = u->frac >> 32;
		u->frac &= 0xffffffffu;
		while (adv) {
			adv--;
			u->m++;
			u->y0i = u->y1i;
			u->y0q = u->y1q;
			if (available(u, u->m + 1)) {
				fir_value(u, u->m + 1, &u->y1i, &u->y1q);
			} else {
				/* recomputed when the input arrives */
				u->m += adv;
				adv = 0;
				u->pair_valid = 0;
			}
		}
	}

	return i;
}

size_t iq_upconv_pull(iq_upconv_t *u, int8_t *out, size_t n)
{
	if (!u || !out)
		return 0;

	return produce(u, out, NULL, n);
}

size_t iq_upconv_pull_float(iq_upconv_t *u, float *out, size_t n)
{
	if (!u || !out)
		return 0;

	return produce(u, NULL, out, n);
}

void iq_upconv_get_stats(const iq_upconv_t *u, struct iq_upconv_stats *stats)
{
	if (!u || !stats)
		return;

	stats->in_samples = u->in_real;
	stats->out_samples = u->out_count;
	stats->clipped = u->clipped;
	stats->in_power = u->in_real ? u->power_sum / u->in_real : 0.0;
	stats->in_peak = u->peak;
	stats->in_position = (u->m + u->frac / TWO_POW_32) / u->interp;
}
