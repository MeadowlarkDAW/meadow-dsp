//! An implementation of Andrew Simper's SVF (state variable filter) model (f64 version).
//! https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf

use std::f64::consts::PI;

use super::f32::SvfCoeff as SvfCoeffF32;

pub const Q_BUTTERWORTH_ORD2: f64 = 0.70710678118654752440;
pub const Q_BUTTERWORTH_ORD4: [f64; 2] = [0.54119610014619698440, 1.3065629648763765279];
pub const Q_BUTTERWORTH_ORD6: [f64; 3] = [
    0.51763809020504152470,
    0.70710678118654752440,
    1.9318516525781365735,
];
pub const Q_BUTTERWORTH_ORD8: [f64; 4] = [
    0.50979557910415916894,
    0.60134488693504528054,
    0.89997622313641570464,
    2.5629154477415061788,
];

pub const ORD4_Q_SCALE: f64 = 0.35;
pub const ORD6_Q_SCALE: f64 = 0.2;
pub const ORD8_Q_SCALE: f64 = 0.14;

/// The coefficients for an SVF (state variable filter) model.
#[derive(Default, Clone, Copy)]
pub struct SvfCoeff {
    pub a1: f64,
    pub a2: f64,
    pub a3: f64,

    pub m0: f64,
    pub m1: f64,
    pub m2: f64,
}

impl SvfCoeff {
    pub fn lowpass_ord2(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> Self {
        let g = g(cutoff_hz, sample_rate_recip);
        let k = 1.0 / q;

        Self::from_g_and_k(g, k, 0.0, 0.0, 1.0)
    }

    pub fn lowpass_ord4(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> [Self; 2] {
        let g = g(cutoff_hz, sample_rate_recip);
        let q_norm = scale_q_norm_for_order(q_norm(q), ORD4_Q_SCALE);

        std::array::from_fn(|i| {
            let q = q_norm * Q_BUTTERWORTH_ORD4[i];
            let k = 1.0 / q;

            Self::from_g_and_k(g, k, 0.0, 0.0, 1.0)
        })
    }

    pub fn lowpass_ord6(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> [Self; 3] {
        let g = g(cutoff_hz, sample_rate_recip);
        let q_norm = scale_q_norm_for_order(q_norm(q), ORD6_Q_SCALE);

        std::array::from_fn(|i| {
            let q = q_norm * Q_BUTTERWORTH_ORD6[i];
            let k = 1.0 / q;

            Self::from_g_and_k(g, k, 0.0, 0.0, 1.0)
        })
    }

    pub fn lowpass_ord8(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> [Self; 4] {
        let g = g(cutoff_hz, sample_rate_recip);
        let q_norm = scale_q_norm_for_order(q_norm(q), ORD8_Q_SCALE);

        std::array::from_fn(|i| {
            let q = q_norm * Q_BUTTERWORTH_ORD8[i];
            let k = 1.0 / q;

            Self::from_g_and_k(g, k, 0.0, 0.0, 1.0)
        })
    }

    pub fn highpass_ord2(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> Self {
        let g = g(cutoff_hz, sample_rate_recip);
        let k = 1.0 / q;

        Self::from_g_and_k(g, k, 1.0, -k, -1.0)
    }

    pub fn highpass_ord4(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> [Self; 2] {
        let g = g(cutoff_hz, sample_rate_recip);
        let q_norm = scale_q_norm_for_order(q_norm(q), ORD4_Q_SCALE);

        std::array::from_fn(|i| {
            let q = q_norm * Q_BUTTERWORTH_ORD4[i];
            let k = 1.0 / q;

            Self::from_g_and_k(g, k, 1.0, -k, -1.0)
        })
    }

    pub fn highpass_ord6(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> [Self; 3] {
        let g = g(cutoff_hz, sample_rate_recip);
        let q_norm = scale_q_norm_for_order(q_norm(q), ORD6_Q_SCALE);

        std::array::from_fn(|i| {
            let q = q_norm * Q_BUTTERWORTH_ORD6[i];
            let k = 1.0 / q;

            Self::from_g_and_k(g, k, 1.0, -k, -1.0)
        })
    }

    pub fn highpass_ord8(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> [Self; 4] {
        let g = g(cutoff_hz, sample_rate_recip);
        let q_norm = scale_q_norm_for_order(q_norm(q), ORD8_Q_SCALE);

        std::array::from_fn(|i| {
            let q = q_norm * Q_BUTTERWORTH_ORD8[i];
            let k = 1.0 / q;

            Self::from_g_and_k(g, k, 1.0, -k, -1.0)
        })
    }

    pub fn notch(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> Self {
        let g = g(cutoff_hz, sample_rate_recip);
        let k = 1.0 / q;

        Self::from_g_and_k(g, k, 1.0, -k, 0.0)
    }

    pub fn bell(cutoff_hz: f64, q: f64, gain_db: f64, sample_rate_recip: f64) -> Self {
        let a = gain_db_to_a(gain_db);

        let g = g(cutoff_hz, sample_rate_recip);
        let k = 1.0 / (q * a);

        Self::from_g_and_k(g, k, 1.0, k * (a * a - 1.0), 0.0)
    }

    pub fn low_shelf(cutoff_hz: f64, q: f64, gain_db: f64, sample_rate_recip: f64) -> Self {
        let a = gain_db_to_a(gain_db);

        let g = (PI * cutoff_hz * sample_rate_recip).tan() / a.sqrt();
        let k = 1.0 / q;

        Self::from_g_and_k(g, k, 1.0, k * (a - 1.0), a * a - 1.0)
    }

    pub fn high_shelf(cutoff_hz: f64, q: f64, gain_db: f64, sample_rate_recip: f64) -> Self {
        let a = gain_db_to_a(gain_db);

        let g = (PI * cutoff_hz * sample_rate_recip).tan() / a.sqrt();
        let k = 1.0 / q;

        Self::from_g_and_k(g, k, a * a, k * (1.0 - a) * a, 1.0 - a * a)
    }

    pub fn allpass(cutoff_hz: f64, q: f64, sample_rate_recip: f64) -> Self {
        let g = g(cutoff_hz, sample_rate_recip);
        let k = 1.0 / q;

        Self::from_g_and_k(g, k, 1.0, -2.0 * k, 0.0)
    }

    pub fn from_g_and_k(g: f64, k: f64, m0: f64, m1: f64, m2: f64) -> Self {
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;

        Self {
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
        }
    }

    pub fn to_f32(self) -> SvfCoeffF32 {
        SvfCoeffF32 {
            a1: self.a1 as f32,
            a2: self.a2 as f32,
            a3: self.a3 as f32,
            m0: self.m0 as f32,
            m1: self.m1 as f32,
            m2: self.m2 as f32,
        }
    }
}

/// The state of an SVF (state variable filter) model.
#[derive(Default, Clone, Copy)]
pub struct SvfState {
    pub ic1eq: f64,
    pub ic2eq: f64,
}

impl SvfState {
    #[inline(always)]
    pub fn tick(&mut self, input: f64, coeff: &SvfCoeff) -> f64 {
        let v3 = input - self.ic2eq;
        let v1 = coeff.a1 * self.ic1eq + coeff.a2 * v3;
        let v2 = self.ic2eq + coeff.a2 * self.ic1eq + coeff.a3 * v3;
        self.ic1eq = 2.0 * v1 - self.ic1eq;
        self.ic2eq = 2.0 * v2 - self.ic2eq;

        coeff.m0 * input + coeff.m1 * v1 + coeff.m2 * v2
    }
}

fn g(cutoff_hz: f64, sample_rate_recip: f64) -> f64 {
    (PI * cutoff_hz * sample_rate_recip).tan()
}

fn q_norm(q: f64) -> f64 {
    q * (1.0 / Q_BUTTERWORTH_ORD2)
}

fn gain_db_to_a(gain_db: f64) -> f64 {
    10.0f64.powf(gain_db * (1.0 / 40.0))
}

fn scale_q_norm_for_order(q_norm: f64, scale: f64) -> f64 {
    if q_norm > 1.0 {
        1.0 + ((q_norm - 1.0) * scale)
    } else {
        q_norm
    }
}

#[cfg(feature = "portable-simd")]
pub mod simd {
    use std::{
        array,
        simd::{f64x2, f64x4},
    };

    use super::{SvfCoeff, SvfState};

    /// The coefficients of two SVF (state variable filter) models packed
    /// into an SIMD vector.
    pub struct SvfCoeffx2 {
        pub a1: f64x2,
        pub a2: f64x2,
        pub a3: f64x2,

        pub m0: f64x2,
        pub m1: f64x2,
        pub m2: f64x2,
    }

    impl SvfCoeffx2 {
        pub const fn splat(coeffs: SvfCoeff) -> Self {
            Self {
                a1: f64x2::splat(coeffs.a1),
                a2: f64x2::splat(coeffs.a2),
                a3: f64x2::splat(coeffs.a3),
                m0: f64x2::splat(coeffs.m0),
                m1: f64x2::splat(coeffs.m1),
                m2: f64x2::splat(coeffs.m2),
            }
        }

        pub fn load(coeffs: &[SvfCoeff; 2]) -> Self {
            Self {
                a1: f64x2::from_array(array::from_fn(|i| coeffs[i].a1)),
                a2: f64x2::from_array(array::from_fn(|i| coeffs[i].a2)),
                a3: f64x2::from_array(array::from_fn(|i| coeffs[i].a3)),
                m0: f64x2::from_array(array::from_fn(|i| coeffs[i].m0)),
                m1: f64x2::from_array(array::from_fn(|i| coeffs[i].m1)),
                m2: f64x2::from_array(array::from_fn(|i| coeffs[i].m2)),
            }
        }
    }

    /// The coefficients of four SVF (state variable filter) models packed
    /// into an SIMD vector.
    pub struct SvfCoeffx4 {
        pub a1: f64x4,
        pub a2: f64x4,
        pub a3: f64x4,

        pub m0: f64x4,
        pub m1: f64x4,
        pub m2: f64x4,
    }

    impl SvfCoeffx4 {
        pub const fn splat(coeffs: SvfCoeff) -> Self {
            Self {
                a1: f64x4::splat(coeffs.a1),
                a2: f64x4::splat(coeffs.a2),
                a3: f64x4::splat(coeffs.a3),
                m0: f64x4::splat(coeffs.m0),
                m1: f64x4::splat(coeffs.m1),
                m2: f64x4::splat(coeffs.m2),
            }
        }

        pub fn load(coeffs: &[SvfCoeff; 4]) -> Self {
            Self {
                a1: f64x4::from_array(array::from_fn(|i| coeffs[i].a1)),
                a2: f64x4::from_array(array::from_fn(|i| coeffs[i].a2)),
                a3: f64x4::from_array(array::from_fn(|i| coeffs[i].a3)),
                m0: f64x4::from_array(array::from_fn(|i| coeffs[i].m0)),
                m1: f64x4::from_array(array::from_fn(|i| coeffs[i].m1)),
                m2: f64x4::from_array(array::from_fn(|i| coeffs[i].m2)),
            }
        }
    }

    /// The state of two SVF (state variable filter) models packed into an
    /// SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct SvfStatex2 {
        pub ic1eq: f64x2,
        pub ic2eq: f64x2,
    }

    impl SvfStatex2 {
        pub const fn splat(state: SvfState) -> Self {
            Self {
                ic1eq: f64x2::splat(state.ic1eq),
                ic2eq: f64x2::splat(state.ic2eq),
            }
        }

        pub fn load(&mut self, states: &[SvfState; 2]) -> Self {
            Self {
                ic1eq: f64x2::from_array(array::from_fn(|i| states[i].ic1eq)),
                ic2eq: f64x2::from_array(array::from_fn(|i| states[i].ic2eq)),
            }
        }

        pub fn store(&self, states: &mut [SvfState; 2]) {
            let ic1eq = self.ic1eq.to_array();
            let ic2eq = self.ic2eq.to_array();

            for (i, s) in states.iter_mut().enumerate() {
                s.ic1eq = ic1eq[i];
                s.ic2eq = ic2eq[i];
            }
        }

        #[inline(always)]
        pub fn tick(&mut self, input: f64x2, coeff: &SvfCoeffx2) -> f64x2 {
            const V_2: f64x2 = f64x2::from_array([2.0; 2]);

            let v3 = input - self.ic2eq;
            let v1 = coeff.a1 * self.ic1eq + coeff.a2 * v3;
            let v2 = self.ic2eq + coeff.a2 * self.ic1eq + coeff.a3 * v3;
            self.ic1eq = V_2 * v1 - self.ic1eq;
            self.ic2eq = V_2 * v2 - self.ic2eq;

            coeff.m0 * input + coeff.m1 * v1 + coeff.m2 * v2
        }
    }

    /// The state of four SVF (state variable filter) models packed into an
    /// SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct SvfStatex4 {
        pub ic1eq: f64x4,
        pub ic2eq: f64x4,
    }

    impl SvfStatex4 {
        pub const fn splat(state: SvfState) -> Self {
            Self {
                ic1eq: f64x4::splat(state.ic1eq),
                ic2eq: f64x4::splat(state.ic2eq),
            }
        }

        pub fn load(&mut self, states: &[SvfState; 4]) -> Self {
            Self {
                ic1eq: f64x4::from_array(array::from_fn(|i| states[i].ic1eq)),
                ic2eq: f64x4::from_array(array::from_fn(|i| states[i].ic2eq)),
            }
        }

        pub fn store(&self, states: &mut [SvfState; 4]) {
            let ic1eq = self.ic1eq.to_array();
            let ic2eq = self.ic2eq.to_array();

            for (i, s) in states.iter_mut().enumerate() {
                s.ic1eq = ic1eq[i];
                s.ic2eq = ic2eq[i];
            }
        }

        #[inline(always)]
        pub fn tick(&mut self, input: f64x4, coeff: &SvfCoeffx4) -> f64x4 {
            const V_2: f64x4 = f64x4::from_array([2.0; 4]);

            let v3 = input - self.ic2eq;
            let v1 = coeff.a1 * self.ic1eq + coeff.a2 * v3;
            let v2 = self.ic2eq + coeff.a2 * self.ic1eq + coeff.a3 * v3;
            self.ic1eq = V_2 * v1 - self.ic1eq;
            self.ic2eq = V_2 * v2 - self.ic2eq;

            coeff.m0 * input + coeff.m1 * v1 + coeff.m2 * v2
        }
    }
}
