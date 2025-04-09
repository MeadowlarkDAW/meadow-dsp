use std::f32::consts::PI;

/// The coefficients for a single-pole IIR filter.
#[derive(Default, Clone, Copy)]
pub struct OnePoleCoeff {
    pub a0: f32,
    pub b1: f32,

    pub m0: f32,
    pub m1: f32,
}

impl OnePoleCoeff {
    pub fn lowpass(cutoff_hz: f32, sample_rate_recip: f32) -> Self {
        let b1 = ((-2.0 * PI) * cutoff_hz * sample_rate_recip).exp();
        let a0 = 1.0 - b1;

        Self {
            a0: a0 as f32,
            b1: b1 as f32,
            m0: 0.0,
            m1: 1.0,
        }
    }

    pub fn highpass(cutoff_hz: f32, sample_rate_recip: f32) -> Self {
        let b1 = ((-2.0 * PI) * cutoff_hz * sample_rate_recip).exp();
        let a0 = 1.0 - b1;

        Self {
            a0: a0 as f32,
            b1: b1 as f32,
            m0: 1.0,
            m1: -1.0,
        }
    }
}

/// The state of a single-pole IIR filter.
#[derive(Default, Clone, Copy)]
pub struct OnePoleState {
    z1: f32,
}

impl OnePoleState {
    #[inline(always)]
    pub fn tick(&mut self, input: f32, coeff: &OnePoleCoeff) -> f32 {
        self.z1 = (coeff.a0 * input) + (coeff.b1 * self.z1);
        coeff.m0 * input + coeff.m1 * self.z1
    }
}

#[cfg(feature = "portable-simd")]
pub mod simd {
    use std::{
        array,
        simd::{f32x4, f32x8},
    };

    use super::{OnePoleCoeff, OnePoleState};

    /// The coefficients of four one-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleCoeffx4 {
        pub a0: f32x4,
        pub b1: f32x4,

        pub m0: f32x4,
        pub m1: f32x4,
    }

    impl OnePoleCoeffx4 {
        pub const fn splat(coeffs: OnePoleCoeff) -> Self {
            Self {
                a0: f32x4::splat(coeffs.a0),
                b1: f32x4::splat(coeffs.b1),
                m0: f32x4::splat(coeffs.m0),
                m1: f32x4::splat(coeffs.m1),
            }
        }

        pub fn load(coeffs: &[OnePoleCoeff; 4]) -> Self {
            Self {
                a0: f32x4::from_array(array::from_fn(|i| coeffs[i].a0)),
                b1: f32x4::from_array(array::from_fn(|i| coeffs[i].b1)),
                m0: f32x4::from_array(array::from_fn(|i| coeffs[i].m0)),
                m1: f32x4::from_array(array::from_fn(|i| coeffs[i].m1)),
            }
        }
    }

    /// The coefficients of eight one-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleCoeffx8 {
        pub a0: f32x8,
        pub b1: f32x8,

        pub m0: f32x8,
        pub m1: f32x8,
    }

    impl OnePoleCoeffx8 {
        pub const fn splat(coeffs: OnePoleCoeff) -> Self {
            Self {
                a0: f32x8::splat(coeffs.a0),
                b1: f32x8::splat(coeffs.b1),
                m0: f32x8::splat(coeffs.m0),
                m1: f32x8::splat(coeffs.m1),
            }
        }

        pub fn load(coeffs: &[OnePoleCoeff; 8]) -> Self {
            Self {
                a0: f32x8::from_array(array::from_fn(|i| coeffs[i].a0)),
                b1: f32x8::from_array(array::from_fn(|i| coeffs[i].b1)),
                m0: f32x8::from_array(array::from_fn(|i| coeffs[i].m0)),
                m1: f32x8::from_array(array::from_fn(|i| coeffs[i].m1)),
            }
        }
    }

    /// The state of four single-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleStatex4 {
        z1: f32x4,
    }

    impl OnePoleStatex4 {
        pub const fn splat(state: OnePoleState) -> Self {
            Self {
                z1: f32x4::splat(state.z1),
            }
        }

        pub fn load(states: &[OnePoleState; 4]) -> Self {
            Self {
                z1: f32x4::from_array(array::from_fn(|i| states[i].z1)),
            }
        }

        #[inline(always)]
        pub fn tick(&mut self, input: f32x4, coeff: &OnePoleCoeffx4) -> f32x4 {
            self.z1 = (coeff.a0 * input) + (coeff.b1 * self.z1);
            coeff.m0 * input + coeff.m1 * self.z1
        }
    }

    /// The state of eight single-pole IIR filters packed into an SIMD vector.
    #[derive(Default, Clone, Copy)]
    pub struct OnePoleStatex8 {
        z1: f32x8,
    }

    impl OnePoleStatex8 {
        pub const fn splat(state: OnePoleState) -> Self {
            Self {
                z1: f32x8::splat(state.z1),
            }
        }

        pub fn load(states: &[OnePoleState; 8]) -> Self {
            Self {
                z1: f32x8::from_array(array::from_fn(|i| states[i].z1)),
            }
        }

        #[inline(always)]
        pub fn tick(&mut self, input: f32x8, coeff: &OnePoleCoeffx8) -> f32x8 {
            self.z1 = (coeff.a0 * input) + (coeff.b1 * self.z1);
            coeff.m0 * input + coeff.m1 * self.z1
        }
    }
}
