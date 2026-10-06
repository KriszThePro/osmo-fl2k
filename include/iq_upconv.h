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

#ifndef IQ_UPCONV_H
#define IQ_UPCONV_H

#include <stddef.h>
#include <stdint.h>

/*
 * Converts a complex baseband signal (I/Q) at a low sample rate into a
 * real signal at the DAC sample rate, centred on a carrier frequency:
 *
 *   y(t) = Re{ x(t) * exp(j * 2 * pi * carrier * t) }
 *
 * The signal chain is:
 *
 *   input (cs8/cs16/cf32) -> complex float, full scale +-1.0
 *     -> polyphase FIR interpolator x R1        (in_rate -> f1 = R1 * in_rate)
 *     -> linear fractional interpolator         (f1 -> out_rate)
 *     -> mixer (NCO at the carrier), real part
 *     -> gain, optional TPDF dither, rounding and clipping to int8
 *
 * The module has no I/O, threads or device access, so it can be tested
 * on its own. The output does not depend on how the input is pushed or
 * how the output is pulled (chunk sizes): the same input always gives
 * the same output, bit for bit.
 */

enum iq_format {
	IQ_FMT_CS8,	/* signed 8 bit I, Q; full scale +-127 */
	IQ_FMT_CS16,	/* signed 16 bit little-endian I, Q; full scale +-32767 */
	IQ_FMT_CF32,	/* 32 bit float little-endian I, Q; full scale +-1.0 */
};

struct iq_upconv_cfg {
	double in_rate;		/* input sample rate [Hz] */
	double out_rate;	/* output (DAC) sample rate [Hz], the actual one */
	double carrier;		/* carrier frequency [Hz], 0 <= carrier < out_rate / 2 */
	double gain_db;		/* 0 dB: input full scale -> +-127 at the output */
	int dither;		/* add TPDF dither before rounding */
	int swap_iq;		/* swap I and Q (mirrors the spectrum) */
};

struct iq_upconv_stats {
	uint64_t in_samples;	/* complex input samples accepted */
	uint64_t out_samples;	/* output samples produced */
	uint64_t clipped;	/* output samples clipped to +-127 */
	double in_power;	/* mean |x|^2 of the input (full scale = 1.0) */
	double in_peak;		/* largest |I| or |Q| of the input */
	double in_position;	/* input samples consumed so far (fractional) */
};

typedef struct iq_upconv iq_upconv_t;

/* Returns NULL if the configuration is invalid (see iq_upconv_check). */
iq_upconv_t *iq_upconv_create(const struct iq_upconv_cfg *cfg);
void iq_upconv_destroy(iq_upconv_t *u);

/* NULL if the configuration is usable, otherwise a short explanation. */
const char *iq_upconv_check(const struct iq_upconv_cfg *cfg);

/* The interpolation factor of the FIR stage chosen for the configuration. */
int iq_upconv_interpolation(const iq_upconv_t *u);

/* "cs8", "cs16", "cf32" -> format; returns 0 on success, -1 otherwise. */
int iq_format_from_name(const char *name, enum iq_format *fmt);

/* Size of one complex sample (I and Q) in bytes. */
size_t iq_format_sample_size(enum iq_format fmt);

/*
 * Pushes raw input bytes. Only whole samples are accepted; returns the
 * number of bytes consumed, which is less than len if the internal FIFO
 * is full (pull some output first) or len is not a multiple of the
 * sample size.
 */
size_t iq_upconv_push(iq_upconv_t *u, enum iq_format fmt, const void *buf, size_t len);

/*
 * Produces up to n output samples, limited by the input available.
 * Returns the number of samples written.
 */
size_t iq_upconv_pull(iq_upconv_t *u, int8_t *out, size_t n);

/*
 * The same as iq_upconv_pull, but writes the value before rounding and
 * clipping (in output LSB units, i.e. +-127 is full scale). For tests
 * and measurements; do not mix with iq_upconv_pull on one instance.
 */
size_t iq_upconv_pull_float(iq_upconv_t *u, float *out, size_t n);

/*
 * End of input: appends enough zero samples for the interpolation
 * filter to deliver the last input samples. Pull the rest afterwards.
 */
void iq_upconv_flush(iq_upconv_t *u);

void iq_upconv_get_stats(const iq_upconv_t *u, struct iq_upconv_stats *stats);

#endif /* IQ_UPCONV_H */
