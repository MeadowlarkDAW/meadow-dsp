use std::f64::consts::PI;

use super::f32::OnePoleIirCoeff as OnePoleIirCoeffF32;

/// The coefficients for a single-pole IIR filter.
#[derive(Default, Clone, Copy)]
pub struct OnePoleIirCoeff {
    pub a0: f64,
    pub b1: f64,

    pub m0: f64,
    pub m1: f64,
}

impl OnePoleIirCoeff {
    pub const NO_OP: Self = Self {
        a0: 0.0,
        b1: 0.0,
        m0: 1.0,
        m1: 0.0,
    };

    pub fn lowpass(cutoff_hz: f64, sample_rate_recip: f64) -> Self {
        let b1 = ((-2.0 * PI) * cutoff_hz * sample_rate_recip).exp();
        let a0 = 1.0 - b1;

        Self {
            a0: a0 as f64,
            b1: b1 as f64,
            m0: 0.0,
            m1: 1.0,
        }
    }

    pub fn highpass(cutoff_hz: f64, sample_rate_recip: f64) -> Self {
        let b1 = ((-2.0 * PI) * cutoff_hz * sample_rate_recip).exp();
        let a0 = 1.0 - b1;

        Self {
            a0: a0 as f64,
            b1: b1 as f64,
            m0: 1.0,
            m1: -1.0,
        }
    }

    pub fn to_f32(self) -> OnePoleIirCoeffF32 {
        OnePoleIirCoeffF32 {
            a0: self.a0 as f32,
            b1: self.b1 as f32,
            m0: self.m0 as f32,
            m1: self.m1 as f32,
        }
    }
}

/// The state of a single-pole IIR filter.
#[derive(Default, Clone, Copy)]
pub struct OnePoleIirState {
    z1: f64,
}

impl OnePoleIirState {
    #[inline(always)]
    pub fn tick(&mut self, input: f64, coeff: &OnePoleIirCoeff) -> f64 {
        self.z1 = (coeff.a0 * input) + (coeff.b1 * self.z1);
        coeff.m0 * input + coeff.m1 * self.z1
    }

    #[inline(always)]
    pub fn reset(&mut self) {
        self.z1 = 0.0;
    }
}

#[cfg(feature = "portable-simd")]
pub mod simd {
    use std::{
        array,
        simd::{f64x2, f64x4},
    };

    use super::{OnePoleIirCoeff, OnePoleIirState};

    /// The coefficients of two one-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleIirCoeffx2 {
        pub a0: f64x2,
        pub b1: f64x2,

        pub m0: f64x2,
        pub m1: f64x2,
    }

    impl OnePoleIirCoeffx2 {
        pub const fn splat(coeffs: OnePoleIirCoeff) -> Self {
            Self {
                a0: f64x2::splat(coeffs.a0),
                b1: f64x2::splat(coeffs.b1),
                m0: f64x2::splat(coeffs.m0),
                m1: f64x2::splat(coeffs.m1),
            }
        }

        pub fn load(coeffs: &[OnePoleIirCoeff; 2]) -> Self {
            Self {
                a0: f64x2::from_array(array::from_fn(|i| coeffs[i].a0)),
                b1: f64x2::from_array(array::from_fn(|i| coeffs[i].b1)),
                m0: f64x2::from_array(array::from_fn(|i| coeffs[i].m0)),
                m1: f64x2::from_array(array::from_fn(|i| coeffs[i].m1)),
            }
        }
    }

    /// The coefficients of four one-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleIirCoeffx4 {
        pub a0: f64x4,
        pub b1: f64x4,

        pub m0: f64x4,
        pub m1: f64x4,
    }

    impl OnePoleIirCoeffx4 {
        pub const fn splat(coeffs: OnePoleIirCoeff) -> Self {
            Self {
                a0: f64x4::splat(coeffs.a0),
                b1: f64x4::splat(coeffs.b1),
                m0: f64x4::splat(coeffs.m0),
                m1: f64x4::splat(coeffs.m1),
            }
        }

        pub fn load(coeffs: &[OnePoleIirCoeff; 4]) -> Self {
            Self {
                a0: f64x4::from_array(array::from_fn(|i| coeffs[i].a0)),
                b1: f64x4::from_array(array::from_fn(|i| coeffs[i].b1)),
                m0: f64x4::from_array(array::from_fn(|i| coeffs[i].m0)),
                m1: f64x4::from_array(array::from_fn(|i| coeffs[i].m1)),
            }
        }
    }

    /// The state of two single-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleIirStatex2 {
        z1: f64x2,
    }

    impl OnePoleIirStatex2 {
        pub const fn splat(state: OnePoleIirState) -> Self {
            Self {
                z1: f64x2::splat(state.z1),
            }
        }

        pub fn load(states: &[OnePoleIirState; 2]) -> Self {
            Self {
                z1: f64x2::from_array(array::from_fn(|i| states[i].z1)),
            }
        }

        #[inline(always)]
        pub fn tick(&mut self, input: f64x2, coeff: &OnePoleIirCoeffx2) -> f64x2 {
            self.z1 = (coeff.a0 * input) + (coeff.b1 * self.z1);
            coeff.m0 * input + coeff.m1 * self.z1
        }

        #[inline(always)]
        pub fn reset(&mut self) {
            self.z1 = f64x2::splat(0.0);
        }
    }

    /// The state of four single-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleIirStatex4 {
        z1: f64x4,
    }

    impl OnePoleIirStatex4 {
        pub const fn splat(state: OnePoleIirState) -> Self {
            Self {
                z1: f64x4::splat(state.z1),
            }
        }

        pub fn load(states: &[OnePoleIirState; 4]) -> Self {
            Self {
                z1: f64x4::from_array(array::from_fn(|i| states[i].z1)),
            }
        }

        #[inline(always)]
        pub fn tick(&mut self, input: f64x4, coeff: &OnePoleIirCoeffx4) -> f64x4 {
            self.z1 = (coeff.a0 * input) + (coeff.b1 * self.z1);
            coeff.m0 * input + coeff.m1 * self.z1
        }

        #[inline(always)]
        pub fn reset(&mut self) {
            self.z1 = f64x4::splat(0.0);
        }
    }
}
