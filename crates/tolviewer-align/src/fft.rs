//! A small self-contained complex FFT, used by the MAFFT-style engine.
//!
//! Radix-2 decimation-in-time Cooley-Tukey, iterative, in place. Nothing here
//! is general-purpose numerics: it exists so that
//! [`crate::mafft`] can correlate two residue-vector signals without pulling in
//! an FFT dependency.

use std::f64::consts::PI;

/// Smallest power of two that is at least `n` (and at least 1).
pub fn next_pow2(n: usize) -> usize {
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

/// In-place forward FFT. `re` and `im` must be the same length and that length
/// must be a power of two; otherwise this is a no-op.
pub fn fft(re: &mut [f64], im: &mut [f64]) {
    transform(re, im, false);
}

/// In-place inverse FFT, including the `1/n` scaling.
pub fn ifft(re: &mut [f64], im: &mut [f64]) {
    transform(re, im, true);
    let n = re.len() as f64;
    if n > 0.0 {
        for (r, i) in re.iter_mut().zip(im.iter_mut()) {
            *r /= n;
            *i /= n;
        }
    }
}

fn transform(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();
    if n != im.len() || n < 2 || !n.is_power_of_two() {
        return;
    }

    // Bit-reversal permutation.
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }

    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2usize;
    while len <= n {
        let ang = sign * 2.0 * PI / len as f64;
        let (wr, wi) = (ang.cos(), ang.sin());
        let mut i = 0usize;
        while i < n {
            let (mut cr, mut ci) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let (ur, ui) = (re[i + k], im[i + k]);
                let (vr0, vi0) = (re[i + k + len / 2], im[i + k + len / 2]);
                let vr = vr0 * cr - vi0 * ci;
                let vi = vr0 * ci + vi0 * cr;
                re[i + k] = ur + vr;
                im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr;
                im[i + k + len / 2] = ui - vi;
                let ncr = cr * wr - ci * wi;
                ci = cr * wi + ci * wr;
                cr = ncr;
            }
            i += len;
        }
        len <<= 1;
    }
}

/// Circular cross-correlation of two real signals, zero-padded to `size`.
///
/// The result `c[k]` is `sum_i a[i] * b[i + k]` with indices taken modulo
/// `size`, so a peak at `k` means `a` lines up with `b` shifted left by `k`
/// residues, and a peak at `size - k` means the opposite shift.
pub fn cross_correlation(a: &[f64], b: &[f64], size: usize) -> Vec<f64> {
    let n = next_pow2(size.max(a.len() + b.len()).max(2));
    let mut ar = vec![0.0f64; n];
    let mut ai = vec![0.0f64; n];
    let mut br = vec![0.0f64; n];
    let mut bi = vec![0.0f64; n];
    ar[..a.len()].copy_from_slice(a);
    br[..b.len()].copy_from_slice(b);
    fft(&mut ar, &mut ai);
    fft(&mut br, &mut bi);
    // conj(A) * B
    for k in 0..n {
        let (x, y) = (ar[k], -ai[k]);
        let (u, v) = (br[k], bi[k]);
        ar[k] = x * u - y * v;
        ai[k] = x * v + y * u;
    }
    ifft(&mut ar, &mut ai);
    ar
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dft(re: &[f64], im: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = re.len();
        let mut or = vec![0.0; n];
        let mut oi = vec![0.0; n];
        for (k, (ork, oik)) in or.iter_mut().zip(oi.iter_mut()).enumerate() {
            for t in 0..n {
                let ang = -2.0 * PI * (k * t) as f64 / n as f64;
                *ork += re[t] * ang.cos() - im[t] * ang.sin();
                *oik += re[t] * ang.sin() + im[t] * ang.cos();
            }
        }
        (or, oi)
    }

    #[test]
    fn next_pow2_rounds_up() {
        assert_eq!(next_pow2(1), 1);
        assert_eq!(next_pow2(2), 2);
        assert_eq!(next_pow2(3), 4);
        assert_eq!(next_pow2(1000), 1024);
    }

    #[test]
    fn fft_matches_the_naive_dft() {
        let mut re: Vec<f64> = (0..32).map(|i| ((i * 7 % 11) as f64) - 5.0).collect();
        let mut im: Vec<f64> = (0..32).map(|i| ((i * 3 % 5) as f64) - 2.0).collect();
        let (er, ei) = dft(&re, &im);
        fft(&mut re, &mut im);
        for i in 0..32 {
            assert!((re[i] - er[i]).abs() < 1e-8, "re[{i}] {} vs {}", re[i], er[i]);
            assert!((im[i] - ei[i]).abs() < 1e-8, "im[{i}] {} vs {}", im[i], ei[i]);
        }
    }

    #[test]
    fn ifft_inverts_fft() {
        let orig_re: Vec<f64> = (0..64).map(|i| (i as f64).sin()).collect();
        let orig_im: Vec<f64> = (0..64).map(|i| (i as f64).cos()).collect();
        let mut re = orig_re.clone();
        let mut im = orig_im.clone();
        fft(&mut re, &mut im);
        ifft(&mut re, &mut im);
        for i in 0..64 {
            assert!((re[i] - orig_re[i]).abs() < 1e-9);
            assert!((im[i] - orig_im[i]).abs() < 1e-9);
        }
    }

    #[test]
    fn cross_correlation_finds_a_known_shift() {
        // b is a copy of a shifted right by 5, so the peak sits at lag 5.
        let a: Vec<f64> = (0..40).map(|i| if i % 7 == 0 { 1.0 } else { 0.0 }).collect();
        let mut b = vec![0.0; 5];
        b.extend(a.iter().copied());
        let c = cross_correlation(&a, &b, a.len() + b.len());
        let best = c
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert_eq!(best, 5, "correlation peaked at {best}");
    }

    #[test]
    fn odd_sizes_are_padded_not_rejected() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        let c = cross_correlation(&a, &b, 6);
        assert!(c.len().is_power_of_two());
        assert!(c[0] > c[1]);
    }
}
