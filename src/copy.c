#include "util.h"

void fill_int(int *dst, int n, int v) {
    if (n & 1) {
        *dst++ = v;
    }

    n >>= 1;

    if (n & 1) {
        *dst++ = v;
        *dst++ = v;
    }

    n >>= 1;

    while (n--) {
        dst[0] = v;
        dst[1] = v;
        dst[2] = v;
        dst[3] = v;
        dst += 4;
    }
}

void fill_dbl(double *dst, int n, double v) {
    if (n & 1) {
        *dst++ = v;
    }

    n >>= 1;

    if (n & 1) {
        *dst++ = v;
        *dst++ = v;
    }

    n >>= 1;

    while (n--) {
        dst[0] = v;
        dst[1] = v;
        dst[2] = v;
        dst[3] = v;
        dst += 4;
    }
}


void fill_float(float *dst, int n, float v) {
    if (n & 1) {
        *dst++ = v;
    }

    n >>= 1;

    if (n & 1) {
        *dst++ = v;
        *dst++ = v;
    }

    n >>= 1;

    while (n--) {
        dst[0] = v;
        dst[1] = v;
        dst[2] = v;
        dst[3] = v;
        dst += 4;
    }
}

void fill_char(char *dst, int n, char v) {
    if (n & 1) {
        *dst++ = v;
    }

    n >>= 1;

    if (n & 1) {
        *dst++ = v;
        *dst++ = v;
    }

    n >>= 1;

    while (n--) {
        dst[0] = v;
        dst[1] = v;
        dst[2] = v;
        dst[3] = v;
        dst += 4;
    }
}

void acopy_int(const int *src, int *dst, int n) {
    if (n & 1) {
        *dst++ = *src++;
    }

    n >>= 1;

    if (n & 1) {
        *dst++ = *src++;
        *dst++ = *src++;
    }

    n >>= 1;

    while (n--) {
        dst[0] = src[0];
        dst[1] = src[1];
        dst[2] = src[2];
        dst[3] = src[3];
        dst += 4;
        src += 4;
    }
}


void acopy_dbl(const double *src, double *dst, int n) {
    if (n & 1) {
        *dst++ = *src++;
    }

    n >>= 1;

    if (n & 1) {
        *dst++ = *src++;
        *dst++ = *src++;
    }

    n >>= 1;

    while (n--) {
        dst[0] = src[0];
        dst[1] = src[1];
        dst[2] = src[2];
        dst[3] = src[3];
        dst += 4;
        src += 4;
    }
}
