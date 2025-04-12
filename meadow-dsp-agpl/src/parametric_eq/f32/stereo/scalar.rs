use crate::parametric_eq::f32::{coeff::MeadowEqDspCoeff, state::MeadowEqDspState, EqParams};

/// The DSP for a fully-featured parametric EQ. This version has two channels,
/// does not make use of SIMD optimizations (although the left and right channels
/// may be auto-vectorized together), and has zero latency. Both channels share
/// the same parameters.
///
/// TODO: Get rid of `NUM_BANDS_PLUS_8` const generic once const generic expressions
/// are stabilized. (please rust compiler team)
pub struct MeadowEqDspStereoLinked<const NUM_BANDS: usize, const NUM_BANDS_PLUS_8: usize> {
    coeff: MeadowEqDspCoeff<NUM_BANDS, NUM_BANDS_PLUS_8>,

    left_state: MeadowEqDspState<NUM_BANDS, NUM_BANDS_PLUS_8>,
    right_state: MeadowEqDspState<NUM_BANDS, NUM_BANDS_PLUS_8>,
}

impl<const NUM_BANDS: usize, const NUM_BANDS_PLUS_8: usize>
    MeadowEqDspStereoLinked<NUM_BANDS, NUM_BANDS_PLUS_8>
{
    pub const LATENCY: u32 = 0;

    pub fn new(sample_rate: f64) -> Self {
        Self {
            coeff: MeadowEqDspCoeff::new(sample_rate),
            left_state: MeadowEqDspState::new(),
            right_state: MeadowEqDspState::new(),
        }
    }

    pub fn params(&self) -> &EqParams<NUM_BANDS> {
        &self.coeff.params()
    }

    pub fn set_params(&mut self, params: &EqParams<NUM_BANDS>) {
        self.coeff.set_params(params);
    }

    pub fn needs_param_flush(&self) -> bool {
        self.coeff.needs_param_flush()
    }

    pub fn flush_param_changes(&mut self) {
        if let Some(info) = self.coeff.flush_param_changes() {
            self.left_state.sync(&info);
            self.right_state.sync(&info);
        }
    }

    pub fn process(&mut self, buf_l: &mut [f32], buf_r: &mut [f32]) {
        if self.needs_param_flush() {
            self.flush_param_changes();
        }

        let (one_pole_coeffs, svf_coeffs) = self.coeff.coeffs();

        let (l_one_pole_states, l_svf_states) = self.left_state.states_mut();
        let (r_one_pole_states, r_svf_states) = self.right_state.states_mut();

        if !one_pole_coeffs.is_empty() {
            // Hint to compiler to optimize loop;
            assert_eq!(one_pole_coeffs.len(), l_one_pole_states.len());
            assert_eq!(one_pole_coeffs.len(), r_one_pole_states.len());

            if one_pole_coeffs.len() == 1 {
                for (out_l, out_r) in buf_l.iter_mut().zip(buf_r.iter_mut()) {
                    *out_l = l_one_pole_states[0].tick(*out_l, &one_pole_coeffs[0]);
                    *out_r = r_one_pole_states[0].tick(*out_r, &one_pole_coeffs[0]);
                }
            } else {
                for (out_l, out_r) in buf_l.iter_mut().zip(buf_r.iter_mut()) {
                    let mut l = *out_l;
                    let mut r = *out_r;

                    l = l_one_pole_states[0].tick(l, &one_pole_coeffs[0]);
                    r = r_one_pole_states[0].tick(r, &one_pole_coeffs[0]);

                    l = l_one_pole_states[1].tick(l, &one_pole_coeffs[1]);
                    r = r_one_pole_states[1].tick(r, &one_pole_coeffs[1]);

                    *out_l = l;
                    *out_r = r;
                }
            }
        }

        if !svf_coeffs.is_empty() {
            // Hint to compiler to optimize loop;
            assert_eq!(svf_coeffs.len(), l_svf_states.len());
            assert_eq!(svf_coeffs.len(), r_svf_states.len());

            for (out_l, out_r) in buf_l.iter_mut().zip(buf_r.iter_mut()) {
                let mut l = *out_l;
                let mut r = *out_r;

                for (i, coeff) in svf_coeffs.iter().enumerate() {
                    l = l_svf_states[i].tick(l, coeff);
                    r = r_svf_states[i].tick(r, coeff);
                }

                *out_l = l;
                *out_r = r;
            }
        }
    }
}
