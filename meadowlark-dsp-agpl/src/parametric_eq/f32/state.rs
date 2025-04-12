use arrayvec::ArrayVec;
use meadowlark_dsp_mit::filter::{one_pole_iir::f32::OnePoleIirState, svf::f32::SvfState};

use super::{
    coeff::{StateSyncInfo, MAX_ONE_POLE_FILTERS},
    FilterOrder,
};

/// The struct that manages the filter states for a fully-featured
/// parametric equalizer. (For a single channel).
///
/// TODO: Get rid of `NUM_BANDS_PLUS_8` const generic once const generic expressions
/// are stabilized. (please rust compiler team)
pub struct MeadowEqDspState<const NUM_BANDS: usize, const NUM_BANDS_PLUS_8: usize> {
    lp_band: MultiOrderBand,
    hp_band: MultiOrderBand,

    bands: [SecondOrderBand; NUM_BANDS],

    one_pole_states: ArrayVec<OnePoleIirState, MAX_ONE_POLE_FILTERS>,
    svf_states: ArrayVec<SvfState, NUM_BANDS_PLUS_8>,
}

impl<const NUM_BANDS: usize, const NUM_BANDS_PLUS_8: usize>
    MeadowEqDspState<NUM_BANDS, NUM_BANDS_PLUS_8>
{
    pub fn new() -> Self {
        Self {
            lp_band: MultiOrderBand::default(),
            hp_band: MultiOrderBand::default(),
            bands: [SecondOrderBand::default(); NUM_BANDS],
            one_pole_states: ArrayVec::new(),
            svf_states: ArrayVec::new(),
        }
    }

    pub fn sync(&mut self, info: &StateSyncInfo<NUM_BANDS>) {
        let mut one_pole_iir_i = 0;
        let mut svf_i = 0;

        if self.lp_band.enabled {
            self.lp_band.sync_states(
                &mut self.one_pole_states,
                &mut self.svf_states,
                &mut one_pole_iir_i,
                &mut svf_i,
            );
        } else {
            self.lp_band.reset();
        }
        self.lp_band.enabled = info.lp_band_enabled;
        self.lp_band.order = info.lp_band_order;

        if self.hp_band.enabled {
            self.hp_band.sync_states(
                &mut self.one_pole_states,
                &mut self.svf_states,
                &mut one_pole_iir_i,
                &mut svf_i,
            );
        } else {
            self.hp_band.reset();
        }
        self.hp_band.enabled = info.hp_band_enabled;
        self.hp_band.order = info.hp_band_order;

        for i in 0..NUM_BANDS {
            if self.bands[i].enabled {
                self.bands[i].svf_state = self.svf_states[svf_i];
                svf_i += 1;
            } else {
                self.bands[i].reset();
            }

            self.bands[i].enabled = info.bands_enabled[i];
        }

        if self.lp_band.enabled {
            self.lp_band
                .add_states(&mut self.one_pole_states, &mut self.svf_states);
        }

        if self.hp_band.enabled {
            self.hp_band
                .add_states(&mut self.one_pole_states, &mut self.svf_states);
        }

        for i in 0..NUM_BANDS {
            if self.bands[i].enabled {
                self.svf_states[i] = self.bands[i].svf_state;
                svf_i += 1;
            }
        }
    }

    pub fn states_mut(
        &mut self,
    ) -> (
        &mut ArrayVec<OnePoleIirState, MAX_ONE_POLE_FILTERS>,
        &mut ArrayVec<SvfState, NUM_BANDS_PLUS_8>,
    ) {
        (&mut self.one_pole_states, &mut self.svf_states)
    }
}

#[derive(Default, Clone, Copy)]
struct SecondOrderBand {
    enabled: bool,
    svf_state: SvfState,
}

impl SecondOrderBand {
    fn reset(&mut self) {
        self.svf_state.reset();
    }
}

#[derive(Default, Clone, Copy)]
struct MultiOrderBand {
    enabled: bool,
    order: FilterOrder,

    one_pole_iir_state: OnePoleIirState,
    svf_states: [SvfState; 4],
}

impl MultiOrderBand {
    fn sync_states<const NUM_BANDS_PLUS_8: usize>(
        &mut self,
        one_pole_states: &mut ArrayVec<OnePoleIirState, MAX_ONE_POLE_FILTERS>,
        svf_states: &mut ArrayVec<SvfState, NUM_BANDS_PLUS_8>,
        one_pole_iir_i: &mut usize,
        svf_i: &mut usize,
    ) {
        match self.order {
            FilterOrder::X1 => {
                self.one_pole_iir_state = one_pole_states[*one_pole_iir_i];
                *one_pole_iir_i += 1;
            }
            FilterOrder::X2 => {
                self.svf_states[0] = svf_states[*svf_i];
                *svf_i += 1;
            }
            FilterOrder::X4 => {
                self.svf_states[0] = svf_states[*svf_i];
                self.svf_states[1] = svf_states[*svf_i + 1];
                *svf_i += 2;
            }
            FilterOrder::X6 => {
                self.svf_states[0] = svf_states[*svf_i];
                self.svf_states[1] = svf_states[*svf_i + 1];
                self.svf_states[2] = svf_states[*svf_i + 2];
                *svf_i += 3;
            }
            FilterOrder::X8 => {
                self.svf_states[0] = svf_states[*svf_i];
                self.svf_states[1] = svf_states[*svf_i + 1];
                self.svf_states[2] = svf_states[*svf_i + 2];
                self.svf_states[3] = svf_states[*svf_i + 3];
                *svf_i += 4;
            }
        }
    }

    fn add_states<const NUM_BANDS_PLUS_8: usize>(
        &self,
        one_pole_states: &mut ArrayVec<OnePoleIirState, MAX_ONE_POLE_FILTERS>,
        svf_states: &mut ArrayVec<SvfState, NUM_BANDS_PLUS_8>,
    ) {
        match self.order {
            FilterOrder::X1 => {
                one_pole_states.push(self.one_pole_iir_state);
            }
            FilterOrder::X2 => {
                svf_states.push(self.svf_states[0]);
            }
            FilterOrder::X4 => {
                svf_states.push(self.svf_states[0]);
                svf_states.push(self.svf_states[1]);
            }
            FilterOrder::X6 => {
                svf_states.push(self.svf_states[0]);
                svf_states.push(self.svf_states[1]);
                svf_states.push(self.svf_states[2]);
            }
            FilterOrder::X8 => {
                svf_states.push(self.svf_states[0]);
                svf_states.push(self.svf_states[1]);
                svf_states.push(self.svf_states[2]);
                svf_states.push(self.svf_states[3]);
            }
        }
    }

    fn reset(&mut self) {
        self.one_pole_iir_state.reset();
        self.svf_states = [SvfState::default(); 4];
    }
}
