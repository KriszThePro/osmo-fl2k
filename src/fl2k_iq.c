/*
 * osmo-fl2k, turns FL2000-based USB 3.0 to VGA adapters into
 * low cost DACs
 *
 * fl2k_iq: transmits a complex baseband (I/Q) signal on a carrier
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

/*
 * The input is a raw I/Q stream (e.g. from an SDR application, GNU Radio
 * or a software modulator) at a low sample rate. It is interpolated to
 * the DAC sample rate, mixed up to the carrier frequency and sent to the
 * red channel of the FL2000 as a real signal:
 *
 *   modulator | fl2k_iq -c 7325000 -
 *
 * The DAC is the only clock: the input is read as fast as the device
 * consumes the samples, so a pipe throttles the producer by itself.
 *
 * With -w the samples go to a file instead of the device. No hardware is
 * needed then, which is useful for testing: the file can be analysed or
 * downconverted and fed to a receiver.
 */

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef _WIN32
#include <getopt.h>
#include <unistd.h>
#else
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#ifdef __MINGW32__
#include <getopt.h>
#else
#include "getopt/getopt.h"
#endif
#endif

#include <pthread.h>

#include "osmo-fl2k.h"
#include "iq_upconv.h"

#define DEFAULT_IN_RATE		48000
#define MIN_DEFAULT_RATE	10000000
#define MAX_DEFAULT_RATE	100000000
#define DEFAULT_BUFFERS		4
#define READ_CHUNK		16384

/* zero buffers sent after the end of the input before stopping (xfers in flight) */
#define DRAIN_BUFFERS		8

enum buf_state { BUF_FREE, BUF_FILLED };

static fl2k_dev_t *dev = NULL;
static volatile int do_exit = 0;

static FILE *infile;
static enum iq_format in_fmt = IQ_FMT_CS16;
static iq_upconv_t *upconv;
static uint8_t readbuf[READ_CHUNK];
static size_t pending;			/* bytes in readbuf not yet accepted */
static int input_done;

/* output buffers handed to the library callback, in a fixed cyclic order */
static struct {
	pthread_mutex_t lock;
	pthread_cond_t cond;
	int8_t **buf;
	enum buf_state *state;
	int num;
	int fill_pos;		/* next buffer to fill (producer) */
	int send_pos;		/* next buffer to send (callback) */
	int in_flight;		/* buffer given to the library in the last callback, or -1 */
	int started;		/* prefill done, the callback may send */
	int finished;		/* no more buffers will be filled */
	int drain_count;
	int drained;		/* everything sent, the device can be stopped */
	uint64_t underflows;
} ring;

static int8_t *zero_buf;

void usage(void)
{
	fprintf(stderr,
		"fl2k_iq, an I/Q transmitter for FL2K VGA dongles\n\n"
		"Usage:\n"
		"\t-c carrier frequency in Hz\n"
		"\t[-d device index (default: 0)]\n"
		"\t[-s samplerate in Hz (default: 4 * carrier, 10 to 100 MS/s)]\n"
		"\t[-i input sample rate in Hz (default: 48000)]\n"
		"\t[-F input format: cs8, cs16 or cf32 (default: cs16)]\n"
		"\t[-g gain in dB, 0 dB maps input full scale to the DAC full scale (default: 0)]\n"
		"\t[-p sample rate correction in ppm (default: 0)]\n"
		"\t[-b number of output buffers (default: 4)]\n"
		"\t[-w write the samples to a file instead of the device ('-' for stdout)]\n"
		"\t[--dither (add TPDF dither before quantisation)]\n"
		"\t[--swap-iq (swap I and Q, mirrors the spectrum)]\n"
		"\tfilename (use '-' to read from stdin)\n\n"
		"Input: interleaved I/Q, little-endian; full scale is +-127 (cs8),\n"
		"+-32767 (cs16) or +-1.0 (cf32). The carrier is the centre (0 Hz) of\n"
		"the input spectrum; the output goes to the red channel.\n\n"
	);
	exit(1);
}

#ifdef _WIN32
BOOL WINAPI
sighandler(int signum)
{
	if (CTRL_C_EVENT == signum) {
		fprintf(stderr, "Signal caught, exiting!\n");
		do_exit = 1;
		return TRUE;
	}
	return FALSE;
}
#else
static void sighandler(int signum)
{
	(void)signum;
	fprintf(stderr, "Signal caught, exiting!\n");
	do_exit = 1;
}
#endif

/* Gives the upconverter more input. Returns 0 if there is nothing left. */
static int feed(void)
{
	size_t size = iq_format_sample_size(in_fmt);
	size_t used;

	if (input_done)
		return 0;

	if (pending < size) {
		size_t r = fread(readbuf + pending, 1, sizeof(readbuf) - pending, infile);

		if (r == 0) {
			if (ferror(infile))
				fprintf(stderr, "Error reading the input.\n");
			else if (pending)
				fprintf(stderr, "Warning: %zu trailing byte(s) of an incomplete sample ignored.\n", pending);
			/* deliver the samples still inside the interpolation filter */
			iq_upconv_flush(upconv);
			input_done = 1;
			return 1;
		}
		pending += r;
	}

	used = iq_upconv_push(upconv, in_fmt, readbuf, pending);
	memmove(readbuf, readbuf + used, pending - used);
	pending -= used;
	return 1;
}

/* Fills buf with up to len output samples; less only at the end of the input. */
static size_t fill(int8_t *buf, size_t len)
{
	size_t got = 0;

	while (got < len && !do_exit) {
		got += iq_upconv_pull(upconv, buf + got, len - got);
		if (got < len && !feed())
			break;
	}

	return got;
}

static double dbfs(double amplitude)
{
	return amplitude > 0 ? 20.0 * log10(amplitude) : -INFINITY;
}

/* Status line on stderr: input level, clipping, underflows. */
static void print_stats(double out_rate, int final)
{
	static struct iq_upconv_stats last;
	static int have_last;
	struct iq_upconv_stats st;
	double power, clip_pct;
	uint64_t dn, dout;

	iq_upconv_get_stats(upconv, &st);
	dn = st.in_samples - (have_last ? last.in_samples : 0);
	dout = st.out_samples - (have_last ? last.out_samples : 0);

	if (final) {
		power = st.in_power;
		clip_pct = st.out_samples ? 100.0 * st.clipped / st.out_samples : 0.0;
	} else {
		double sum = st.in_power * st.in_samples - (have_last ? last.in_power * last.in_samples : 0.0);
		power = dn ? sum / dn : 0.0;
		clip_pct = dout ? 100.0 * (st.clipped - (have_last ? last.clipped : 0)) / dout : 0.0;
	}

	fprintf(stderr, "%s%.1f s: input %.1f dBFS rms, peak %.1f dBFS, clipped %.4f %%, underflows %llu\n",
		final ? "Total " : "", st.out_samples / out_rate, dbfs(sqrt(power)), dbfs(st.in_peak),
		clip_pct, (unsigned long long)ring.underflows);

	last = st;
	have_last = 1;
}

/* Library callback: never blocks, hands over the next filled buffer or silence. */
void fl2k_callback(fl2k_data_info_t *data_info)
{
	if (data_info->device_error) {
		fprintf(stderr, "Device error, exiting.\n");
		do_exit = 1;
	}

	pthread_mutex_lock(&ring.lock);

	/* the library has copied the previous buffer, it can be refilled */
	if (ring.in_flight >= 0) {
		ring.state[ring.in_flight] = BUF_FREE;
		ring.in_flight = -1;
	}

	data_info->sampletype_signed = 1;
	if (ring.started && ring.state[ring.send_pos] == BUF_FILLED) {
		data_info->r_buf = (char *)ring.buf[ring.send_pos];
		ring.in_flight = ring.send_pos;
		ring.send_pos = (ring.send_pos + 1) % ring.num;
	} else {
		data_info->r_buf = (char *)zero_buf;
		if (ring.finished) {
			if (++ring.drain_count >= DRAIN_BUFFERS)
				ring.drained = 1;
		} else if (ring.started) {
			ring.underflows++;
		}
	}

	pthread_cond_signal(&ring.cond);
	pthread_mutex_unlock(&ring.lock);
}

static int ring_init(int num)
{
	int i;

	ring.num = num;
	ring.in_flight = -1;
	ring.buf = calloc(num, sizeof(*ring.buf));
	ring.state = calloc(num, sizeof(*ring.state));
	zero_buf = calloc(1, FL2K_BUF_LEN);
	if (!ring.buf || !ring.state || !zero_buf)
		return -1;

	for (i = 0; i < num; i++) {
		ring.buf[i] = malloc(FL2K_BUF_LEN);
		if (!ring.buf[i])
			return -1;
		ring.state[i] = BUF_FREE;
	}

	pthread_mutex_init(&ring.lock, NULL);
	pthread_cond_init(&ring.cond, NULL);
	return 0;
}

static void ring_free(void)
{
	int i;

	if (ring.buf)
		for (i = 0; i < ring.num; i++)
			free(ring.buf[i]);
	free(ring.buf);
	free(ring.state);
	free(zero_buf);
}

/* Device mode: the main thread reads and converts, the callback sends. */
static int run_device(double out_rate)
{
	time_t last_stats = time(NULL);
	uint64_t filled = 0;

	while (!do_exit) {
		size_t n;
		int8_t *buf;

		pthread_mutex_lock(&ring.lock);
		while (ring.state[ring.fill_pos] != BUF_FREE && !do_exit)
			pthread_cond_wait(&ring.cond, &ring.lock);
		buf = ring.buf[ring.fill_pos];
		pthread_mutex_unlock(&ring.lock);
		if (do_exit)
			break;

		n = fill(buf, FL2K_BUF_LEN);
		if (n > 0 && n < FL2K_BUF_LEN)
			memset(buf + n, 0, FL2K_BUF_LEN - n);

		pthread_mutex_lock(&ring.lock);
		if (n > 0) {
			ring.state[ring.fill_pos] = BUF_FILLED;
			ring.fill_pos = (ring.fill_pos + 1) % ring.num;
			filled++;
		}
		if (n < FL2K_BUF_LEN)
			ring.finished = 1;
		/* start sending once all but one buffer are full (or the input ended) */
		if (!ring.started && (filled >= (uint64_t)(ring.num - 1) || ring.finished))
			ring.started = 1;
		pthread_mutex_unlock(&ring.lock);

		if (ring.finished)
			break;

		if (time(NULL) != last_stats) {
			last_stats = time(NULL);
			print_stats(out_rate, 0);
		}
	}

	/* let the callback send the remaining buffers */
	pthread_mutex_lock(&ring.lock);
	while (!ring.drained && !do_exit)
		pthread_cond_wait(&ring.cond, &ring.lock);
	pthread_mutex_unlock(&ring.lock);

	return 0;
}

/* File mode: no device, the samples go to a file or stdout. */
static int run_file(FILE *out, double out_rate)
{
	time_t last_stats = time(NULL);
	int8_t *buf = malloc(FL2K_BUF_LEN);
	int ret = 0;

	if (!buf) {
		fprintf(stderr, "malloc error!\n");
		return -1;
	}

	while (!do_exit) {
		size_t n = fill(buf, FL2K_BUF_LEN);

		if (n == 0)
			break;
		if (fwrite(buf, 1, n, out) != n) {
			fprintf(stderr, "Error writing the output.\n");
			ret = -1;
			break;
		}
		if (n < FL2K_BUF_LEN)
			break;

		if (time(NULL) != last_stats) {
			last_stats = time(NULL);
			print_stats(out_rate, 0);
		}
	}

	free(buf);
	return ret;
}

int main(int argc, char **argv)
{
#ifndef _WIN32
	struct sigaction sigact, sigign;
#endif
	struct iq_upconv_cfg cfg;
	const char *error;
	char *filename = NULL, *outname = NULL;
	FILE *outfile = NULL;
	double carrier = -1.0, samp_rate = 0.0, in_rate = DEFAULT_IN_RATE;
	double gain_db = 0.0, ppm = 0.0, out_rate;
	int dev_index = 0, buffers = DEFAULT_BUFFERS, dither = 0, swap_iq = 0;
	int opt, option_index = 0, r, ret = 1;

	static struct option long_options[] =
	{
		{"dither",  no_argument, NULL, 'D'},
		{"swap-iq", no_argument, NULL, 'S'},
		{0, 0, 0, 0}
	};

	while ((opt = getopt_long(argc, argv, "b:c:d:F:g:i:p:s:w:", long_options, &option_index)) != -1) {
		switch (opt) {
		case 'b':
			buffers = atoi(optarg);
			break;
		case 'c':
			carrier = atof(optarg);
			break;
		case 'd':
			dev_index = atoi(optarg);
			break;
		case 'F':
			if (iq_format_from_name(optarg, &in_fmt) < 0) {
				fprintf(stderr, "Unknown input format '%s'.\n", optarg);
				usage();
			}
			break;
		case 'g':
			gain_db = atof(optarg);
			break;
		case 'i':
			in_rate = atof(optarg);
			break;
		case 'p':
			ppm = atof(optarg);
			break;
		case 's':
			samp_rate = atof(optarg);
			break;
		case 'w':
			outname = optarg;
			break;
		case 'D':
			dither = 1;
			break;
		case 'S':
			swap_iq = 1;
			break;
		default:
			usage();
			break;
		}
	}

	if (argc <= optind || carrier < 0)
		usage();
	filename = argv[optind];

	if (dev_index < 0 || buffers < 2) {
		fprintf(stderr, "Invalid device index or buffer count (at least 2).\n");
		return 1;
	}

	if (samp_rate <= 0) {
		samp_rate = 4.0 * carrier;
		if (samp_rate < MIN_DEFAULT_RATE)
			samp_rate = MIN_DEFAULT_RATE;
		if (samp_rate > MAX_DEFAULT_RATE)
			samp_rate = MAX_DEFAULT_RATE;
	}

	if (strcmp(filename, "-") == 0) { /* Read samples from stdin */
		infile = stdin;
#ifdef _WIN32
		_setmode(_fileno(stdin), _O_BINARY);
#endif
	} else {
		infile = fopen(filename, "rb");
		if (!infile) {
			fprintf(stderr, "Failed to open %s\n", filename);
			return -ENOENT;
		}
	}

	if (outname) {
		if (strcmp(outname, "-") == 0) {
			outfile = stdout;
#ifdef _WIN32
			_setmode(_fileno(stdout), _O_BINARY);
#endif
		} else {
			outfile = fopen(outname, "wb");
			if (!outfile) {
				fprintf(stderr, "Failed to open %s\n", outname);
				goto out;
			}
		}
	}

#ifndef _WIN32
	sigact.sa_handler = sighandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigign.sa_handler = SIG_IGN;
	sigemptyset(&sigign.sa_mask);
	sigign.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	sigaction(SIGQUIT, &sigact, NULL);
	sigaction(SIGPIPE, &sigign, NULL);
#else
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) sighandler, TRUE );
#endif

	if (!outfile) {
		if (ring_init(buffers) < 0) {
			fprintf(stderr, "malloc error!\n");
			goto out;
		}

		fl2k_open(&dev, (uint32_t)dev_index);
		if (NULL == dev) {
			fprintf(stderr, "Failed to open fl2k device #%d.\n", dev_index);
			goto out;
		}

		/* the callback sends silence until the upconverter is ready */
		r = fl2k_start_tx(dev, fl2k_callback, NULL, 0);
		if (r < 0) {
			fprintf(stderr, "Failed to start the transmission.\n");
			goto out;
		}

		r = fl2k_set_sample_rate(dev, (uint32_t)samp_rate);
		if (r < 0)
			fprintf(stderr, "WARNING: Failed to set sample rate. %d\n", r);

		/* the PLL cannot hit every rate: work with the actual one */
		samp_rate = fl2k_get_sample_rate(dev);
	}

	out_rate = samp_rate * (1.0 + ppm * 1e-6);

	cfg.in_rate = in_rate;
	cfg.out_rate = out_rate;
	cfg.carrier = carrier;
	cfg.gain_db = gain_db;
	cfg.dither = dither;
	cfg.swap_iq = swap_iq;

	error = iq_upconv_check(&cfg);
	if (!error && carrier > 0.45 * out_rate)
		error = "carrier too close to half the sample rate (use at most 0.45 * samplerate)";
	if (!error && carrier < in_rate / 2)
		error = "carrier below half the input sample rate (the spectrum would fold over 0 Hz)";
	if (error) {
		fprintf(stderr, "Error: %s.\n", error);
		goto out;
	}

	upconv = iq_upconv_create(&cfg);
	if (!upconv) {
		fprintf(stderr, "malloc error!\n");
		goto out;
	}

	fprintf(stderr, "Samplerate:\t%3.6f MHz%s\n", out_rate / 1e6, outfile ? " (file output)" : "");
	fprintf(stderr, "Carrier:\t%3.6f MHz\n", carrier / 1e6);
	fprintf(stderr, "Input:\t\t%s at %.0f Hz, interpolation %d, gain %.1f dB%s%s\n",
		in_fmt == IQ_FMT_CS8 ? "cs8" : in_fmt == IQ_FMT_CS16 ? "cs16" : "cf32",
		in_rate, iq_upconv_interpolation(upconv), gain_db,
		dither ? ", dither" : "", swap_iq ? ", I/Q swapped" : "");
	fprintf(stderr, "Images:\t\t%3.6f MHz, %3.6f MHz (filter the output!)\n",
		(out_rate - carrier) / 1e6, (out_rate + carrier) / 1e6);

	if (outfile)
		ret = run_file(outfile, out_rate) < 0;
	else
		ret = run_device(out_rate) < 0;

	print_stats(out_rate, 1);

out:
	if (dev) {
		fl2k_stop_tx(dev);
		fl2k_close(dev);
	}

	iq_upconv_destroy(upconv);
	ring_free();

	if (outfile && outfile != stdout)
		fclose(outfile);

	if (infile && (infile != stdin))
		fclose(infile);

	return ret;
}
