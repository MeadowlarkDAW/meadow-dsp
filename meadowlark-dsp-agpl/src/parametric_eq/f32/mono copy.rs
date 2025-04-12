use meadowlark_dsp_mit::filter::{
    one_pole_iir::{
        f32::{OnePoleCoeff, OnePoleState},
        f64::OnePoleCoeff as OnePoleCoeffF64,
    },
    svf::{
        f32::{SvfCoeff, SvfState},
        f64::SvfCoeff as SvfCoeffF64,
    },
};

use super::{BandParams, BandType, EqParams, FilterOrder, LpOrHpBandParams};

pub const DEFAULT_Q: f32 = meadowlark_dsp_mit::filter::svf::f64::Q_BUTTERWORTH_ORD2 as f32;

const MAX_NUM_PACKED_ONE_POLE_FILTERS: usize = 2;

/// The struct that manages the filter coefficients for a fully-featured
/// parametric equalizer. (For a single channel).
pub struct MeadowEqDspCoeff<const NUM_BANDS: usize> {
    params: EqParams<NUM_BANDS>,

    lp_band: MultiOrderBand,
    hp_band: MultiOrderBand,

    bands: [SecondOrderBand; NUM_BANDS],

    packed_one_pole_filters: Vec<PackedOnePoleIIR>,
    packed_svf_filters: Vec<PackedSvf>,

    needs_param_flush: bool,
    lp_band_needs_recalc: bool,
    hp_band_needs_recalc: bool,
    bands_needing_recalc: [bool; NUM_BANDS],

    sample_rate_recip: f64,
}

impl<const NUM_BANDS: usize> MeadowEqDsp<NUM_BANDS> {
    pub fn new(sample_rate: f64) -> Self {
        let sample_rate_recip = sample_rate.recip();

        let max_num_packed_svf_filters = 4 + 4 + NUM_BANDS;

        Self {
            params: EqParams::default(),
            lp_band: MultiOrderBand::default(),
            hp_band: MultiOrderBand::default(),
            bands: [SecondOrderBand::default(); NUM_BANDS],
            packed_one_pole_filters: Vec::with_capacity(MAX_NUM_PACKED_ONE_POLE_FILTERS),
            packed_svf_filters: Vec::with_capacity(max_num_packed_svf_filters),
            needs_param_flush: false,
            lp_band_needs_recalc: false,
            hp_band_needs_recalc: false,
            bands_needing_recalc: [false; NUM_BANDS],
            sample_rate_recip,
        }
    }

    pub fn params(&self) -> &EqParams<NUM_BANDS> {
        &self.params
    }

    pub fn set_params(&mut self, params: &EqParams<NUM_BANDS>) {
        if self.params.lp_band != params.lp_band {
            self.params.lp_band = params.lp_band;
            self.lp_band_needs_recalc = true;
            self.needs_param_flush = true;
        }
        if self.params.hp_band != params.hp_band {
            self.params.hp_band = params.hp_band;
            self.hp_band_needs_recalc = true;
            self.needs_param_flush = true;
        }

        for i in 0..NUM_BANDS {
            if self.params.bands[i] != params.bands[i] {
                self.params.bands[i] = params.bands[i];
                self.bands_needing_recalc[i] = true;
                self.needs_param_flush = true;
            }
        }
    }

    pub fn needs_param_flush(&self) -> bool {
        self.needs_param_flush
    }

    pub fn flush_param_changes(&mut self) {
        if !self.needs_param_flush {
            return;
        }
        self.needs_param_flush = false;

        // -----------------------------------------------------------------------------------

        let mut one_pole_filter_i = 0;
        let mut svf_filter_i = 0;

        if self.lp_band.enabled {
            self.lp_band.sync_filter_states(
                &mut one_pole_filter_i,
                &mut svf_filter_i,
                &self.packed_one_pole_filters,
                &self.packed_svf_filters,
            );
        }
        if self.hp_band.enabled {
            self.hp_band.sync_filter_states(
                &mut one_pole_filter_i,
                &mut svf_filter_i,
                &self.packed_one_pole_filters,
                &self.packed_svf_filters,
            );
        }
        for band_i in 0..NUM_BANDS {
            if self.bands[band_i].enabled {
                self.bands[band_i].packed_svf.state = self.packed_svf_filters[svf_filter_i].state;
                svf_filter_i += 1;
            }
        }

        // -----------------------------------------------------------------------------------

        if self.lp_band_needs_recalc {
            self.lp_band_needs_recalc = false;
            self.lp_band
                .sync_params(&self.params.lp_band, self.sample_rate_recip, true);
        }
        if self.hp_band_needs_recalc {
            self.hp_band_needs_recalc = false;
            self.lp_band
                .sync_params(&self.params.lp_band, self.sample_rate_recip, false);
        }

        for band_i in 0..NUM_BANDS {
            if !self.bands_needing_recalc[band_i] {
                continue;
            }
            self.bands_needing_recalc[band_i] = false;

            self.bands[band_i].sync_params(&self.params.bands[band_i], self.sample_rate_recip);
        }

        // -----------------------------------------------------------------------------------

        self.packed_one_pole_filters.clear();
        self.packed_svf_filters.clear();

        if self.lp_band.enabled {
            self.lp_band.add_filter_states(
                &mut self.packed_one_pole_filters,
                &mut self.packed_svf_filters,
            );
        }
        if self.hp_band.enabled {
            self.hp_band.add_filter_states(
                &mut self.packed_one_pole_filters,
                &mut self.packed_svf_filters,
            );
        }
        for band_i in 0..NUM_BANDS {
            if self.bands[band_i].enabled {
                self.packed_svf_filters.push(self.bands[band_i].packed_svf);
            }
        }
    }

    pub fn sync_params_from(&mut self, other: &mut Self) {
        if other.needs_param_flush {
            other.flush_param_changes();
        }
        self.needs_param_flush = false;

        if !(!self.lp_band.enabled && !other.lp_band.enabled) {
            self.lp_band.sync_params_from(&other.lp_band);
        }
        if !(!self.hp_band.enabled && !other.hp_band.enabled) {
            self.hp_band.sync_params_from(&other.hp_band);
        }

        for i in 0..NUM_BANDS {
            self.bands[i].sync_params_from(&other.bands[i]);
        }

        self.packed_one_pole_filters.clear();
        self.packed_svf_filters.clear();

        self.packed_one_pole_filters
            .extend_from_slice(other.packed_one_pole_filters.as_slice());
        self.packed_svf_filters
            .extend_from_slice(other.packed_svf_filters.as_slice());
    }

    pub fn is_empty(&self) -> bool {
        self.packed_one_pole_filters.is_empty() && self.packed_svf_filters.is_empty()
    }

    // Process a single channel of audio. This version contains no SIMD optimizations, and
    // it has zero latency.
    pub fn process_scalar(&mut self, buffer: &mut [f32]) {
        if self.needs_param_flush {
            self.flush_param_changes();
        }

        if self.is_empty() {
            return;
        }

        for out_s in buffer.iter_mut() {
            let mut s = *out_s;

            for filter in self.packed_one_pole_filters.iter_mut() {
                s = filter.state.tick(s, &filter.coeff);
            }
            for filter in self.packed_svf_filters.iter_mut() {
                s = filter.state.tick(s, &filter.coeff);
            }

            *out_s = s;
        }
    }
}

#[derive(Default, Clone, Copy)]
struct PackedOnePoleIIR {
    coeff: OnePoleCoeff,
    state: OnePoleState,
}

#[derive(Default, Clone, Copy)]
struct PackedSvf {
    coeff: SvfCoeff,
    state: SvfState,
}

#[derive(Default, Clone, Copy)]
struct SecondOrderBand {
    enabled: bool,
    needs_param_sync: bool,

    packed_svf: PackedSvf,
}

impl SecondOrderBand {
    fn sync_params(&mut self, params: &BandParams, sample_rate_recip: f64) {
        self.needs_param_sync = false;

        if !params.enabled {
            self.enabled = false;
            self.packed_svf.coeff = SvfCoeff::NO_OP;
            self.packed_svf.state.reset();
            return;
        }

        self.enabled = true;

        self.packed_svf.coeff = match params.band_type {
            BandType::Bell => SvfCoeffF64::bell(
                params.cutoff_hz as f64,
                params.q as f64,
                params.gain_db as f64,
                sample_rate_recip,
            )
            .to_f32(),
            BandType::LowShelf => SvfCoeffF64::low_shelf(
                params.cutoff_hz as f64,
                params.q as f64,
                params.gain_db as f64,
                sample_rate_recip,
            )
            .to_f32(),
            BandType::HighShelf => SvfCoeffF64::high_shelf(
                params.cutoff_hz as f64,
                params.q as f64,
                params.gain_db as f64,
                sample_rate_recip,
            )
            .to_f32(),
            BandType::Notch => {
                SvfCoeffF64::notch(params.cutoff_hz as f64, params.q as f64, sample_rate_recip)
                    .to_f32()
            }
            BandType::Allpass => {
                SvfCoeffF64::allpass(params.cutoff_hz as f64, params.q as f64, sample_rate_recip)
                    .to_f32()
            }
        };
    }

    fn sync_params_from(&mut self, other: &Self) {
        self.enabled = other.enabled;
        self.needs_param_sync = false;

        if other.enabled {
            self.packed_svf.coeff = other.packed_svf.coeff;
        } else {
            self.packed_svf.state.reset();
        }
    }
}

#[derive(Default)]
struct MultiOrderBand {
    enabled: bool,
    needs_param_sync: bool,
    order: FilterOrder,

    packed_one_pole_iir: PackedOnePoleIIR,
    packed_svfs: [PackedSvf; 4],
}

impl MultiOrderBand {
    fn sync_params(&mut self, params: &LpOrHpBandParams, sample_rate_recip: f64, is_lowpass: bool) {
        self.needs_param_sync = false;

        if !params.enabled {
            self.enabled = false;
            self.packed_one_pole_iir = PackedOnePoleIIR::default();
            self.packed_svfs = [PackedSvf::default(); 4];
            return;
        }

        self.enabled = true;

        let order_changed = self.order != params.order;
        self.order = params.order;

        match params.order {
            FilterOrder::X1 => {
                self.packed_one_pole_iir.coeff = if is_lowpass {
                    OnePoleCoeffF64::lowpass(params.cutoff_hz as f64, sample_rate_recip).to_f32()
                } else {
                    OnePoleCoeffF64::highpass(params.cutoff_hz as f64, sample_rate_recip).to_f32()
                };

                if order_changed {
                    for f in self.packed_svfs.iter_mut() {
                        f.state.reset();
                    }
                }
            }
            FilterOrder::X2 => {
                self.packed_svfs[0].coeff = if is_lowpass {
                    SvfCoeffF64::lowpass_ord2(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                    .to_f32()
                } else {
                    SvfCoeffF64::highpass_ord2(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                    .to_f32()
                };

                if order_changed {
                    self.packed_one_pole_iir.state.reset();
                    for f in self.packed_svfs.iter_mut().skip(1) {
                        f.state.reset();
                    }
                }
            }
            FilterOrder::X4 => {
                let coeffs = if is_lowpass {
                    SvfCoeffF64::lowpass_ord4(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                } else {
                    SvfCoeffF64::highpass_ord4(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                };

                self.packed_svfs[0].coeff = coeffs[0].to_f32();
                self.packed_svfs[1].coeff = coeffs[1].to_f32();

                if order_changed {
                    self.packed_one_pole_iir.state.reset();
                    for f in self.packed_svfs.iter_mut().skip(2) {
                        f.state.reset();
                    }
                }
            }
            FilterOrder::X6 => {
                let coeffs = if is_lowpass {
                    SvfCoeffF64::lowpass_ord4(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                } else {
                    SvfCoeffF64::highpass_ord4(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                };

                self.packed_svfs[0].coeff = coeffs[0].to_f32();
                self.packed_svfs[1].coeff = coeffs[1].to_f32();
                self.packed_svfs[2].coeff = coeffs[2].to_f32();

                if order_changed {
                    self.packed_one_pole_iir.state.reset();
                    for f in self.packed_svfs.iter_mut().skip(3) {
                        f.state.reset();
                    }
                }
            }
            FilterOrder::X8 => {
                let coeffs = if is_lowpass {
                    SvfCoeffF64::lowpass_ord4(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                } else {
                    SvfCoeffF64::highpass_ord4(
                        params.cutoff_hz as f64,
                        params.q as f64,
                        sample_rate_recip,
                    )
                };

                self.packed_svfs[0].coeff = coeffs[0].to_f32();
                self.packed_svfs[1].coeff = coeffs[1].to_f32();
                self.packed_svfs[2].coeff = coeffs[2].to_f32();
                self.packed_svfs[3].coeff = coeffs[3].to_f32();

                if order_changed {
                    self.packed_one_pole_iir.state.reset();
                }
            }
        }
    }

    pub fn sync_params_from(&mut self, other: &Self) {
        self.enabled = other.enabled;
        self.needs_param_sync = false;

        let order_changed = self.order != other.order;
        self.order = other.order;

        if other.enabled {
            match other.order {
                FilterOrder::X1 => {
                    self.packed_one_pole_iir.coeff = other.packed_one_pole_iir.coeff;

                    if order_changed {
                        for f in self.packed_svfs.iter_mut() {
                            f.state.reset();
                        }
                    }
                }
                FilterOrder::X2 => {
                    self.packed_svfs[0].coeff = other.packed_svfs[0].coeff;

                    if order_changed {
                        self.packed_one_pole_iir.state.reset();
                        for f in self.packed_svfs.iter_mut().skip(1) {
                            f.state.reset();
                        }
                    }
                }
                FilterOrder::X4 => {
                    self.packed_svfs[0].coeff = other.packed_svfs[0].coeff;
                    self.packed_svfs[1].coeff = other.packed_svfs[1].coeff;

                    if order_changed {
                        self.packed_one_pole_iir.state.reset();
                        for f in self.packed_svfs.iter_mut().skip(2) {
                            f.state.reset();
                        }
                    }
                }
                FilterOrder::X6 => {
                    self.packed_svfs[0].coeff = other.packed_svfs[0].coeff;
                    self.packed_svfs[1].coeff = other.packed_svfs[1].coeff;
                    self.packed_svfs[2].coeff = other.packed_svfs[2].coeff;

                    if order_changed {
                        self.packed_one_pole_iir.state.reset();
                        for f in self.packed_svfs.iter_mut().skip(3) {
                            f.state.reset();
                        }
                    }
                }
                FilterOrder::X8 => {
                    self.packed_svfs[0].coeff = other.packed_svfs[0].coeff;
                    self.packed_svfs[1].coeff = other.packed_svfs[1].coeff;
                    self.packed_svfs[2].coeff = other.packed_svfs[2].coeff;
                    self.packed_svfs[3].coeff = other.packed_svfs[3].coeff;

                    if order_changed {
                        self.packed_one_pole_iir.state.reset();
                    }
                }
            }
        } else {
            self.packed_one_pole_iir = other.packed_one_pole_iir;
            self.packed_svfs = other.packed_svfs;
        }
    }

    fn add_filter_states(
        &self,
        packed_one_pole_filters: &mut Vec<PackedOnePoleIIR>,
        packed_svf_filters: &mut Vec<PackedSvf>,
    ) {
        match self.order {
            FilterOrder::X1 => {
                packed_one_pole_filters.push(self.packed_one_pole_iir);
            }
            FilterOrder::X2 => packed_svf_filters.push(self.packed_svfs[0]),
            FilterOrder::X4 => {
                packed_svf_filters.push(self.packed_svfs[0]);
                packed_svf_filters.push(self.packed_svfs[1]);
            }
            FilterOrder::X6 => {
                packed_svf_filters.push(self.packed_svfs[0]);
                packed_svf_filters.push(self.packed_svfs[1]);
                packed_svf_filters.push(self.packed_svfs[2]);
            }
            FilterOrder::X8 => {
                packed_svf_filters.push(self.packed_svfs[0]);
                packed_svf_filters.push(self.packed_svfs[1]);
                packed_svf_filters.push(self.packed_svfs[2]);
                packed_svf_filters.push(self.packed_svfs[3]);
            }
        }
    }

    fn sync_filter_states(
        &mut self,
        one_pole_filter_i: &mut usize,
        svf_filter_i: &mut usize,
        packed_one_pole_filters: &Vec<PackedOnePoleIIR>,
        packed_svf_filters: &Vec<PackedSvf>,
    ) {
        match self.order {
            FilterOrder::X1 => {
                self.packed_one_pole_iir.state = packed_one_pole_filters[*one_pole_filter_i].state;
                *one_pole_filter_i += 1;
            }
            FilterOrder::X2 => {
                self.packed_svfs[0].state = packed_svf_filters[*svf_filter_i].state;
                *svf_filter_i += 1;
            }
            FilterOrder::X4 => {
                self.packed_svfs[0].state = packed_svf_filters[*svf_filter_i].state;
                self.packed_svfs[1].state = packed_svf_filters[*svf_filter_i + 1].state;
                *svf_filter_i += 2;
            }
            FilterOrder::X6 => {
                self.packed_svfs[0].state = packed_svf_filters[*svf_filter_i].state;
                self.packed_svfs[1].state = packed_svf_filters[*svf_filter_i + 1].state;
                self.packed_svfs[2].state = packed_svf_filters[*svf_filter_i + 2].state;
                *svf_filter_i += 3;
            }
            FilterOrder::X8 => {
                self.packed_svfs[0].state = packed_svf_filters[*svf_filter_i].state;
                self.packed_svfs[1].state = packed_svf_filters[*svf_filter_i + 1].state;
                self.packed_svfs[2].state = packed_svf_filters[*svf_filter_i + 2].state;
                self.packed_svfs[3].state = packed_svf_filters[*svf_filter_i + 3].state;
                *svf_filter_i += 4;
            }
        }
    }
}
