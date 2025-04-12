use arrayvec::ArrayVec;
use meadowlark_dsp_mit::filter::{
    one_pole_iir::{f32::OnePoleIirCoeff, f64::OnePoleIirCoeff as OnePoleIirCoeffF64},
    svf::{f32::SvfCoeff, f64::SvfCoeff as SvfCoeffF64},
};

use super::{BandParams, BandType, EqParams, FilterOrder, LpOrHpBandParams};

pub const MAX_ONE_POLE_FILTERS: usize = 2;

/// The struct that manages the filter coefficients for a fully-featured
/// parametric equalizer. (For a single channel).
///
/// TODO: Get rid of `NUM_BANDS_PLUS_8` const generic once const generic expressions
/// are stabilized. (please rust compiler team)
pub struct MeadowEqDspCoeff<const NUM_BANDS: usize, const NUM_BANDS_PLUS_8: usize> {
    params: EqParams<NUM_BANDS>,

    lp_band: MultiOrderBand,
    hp_band: MultiOrderBand,

    bands: [SecondOrderBand; NUM_BANDS],

    one_pole_coeffs: ArrayVec<OnePoleIirCoeff, MAX_ONE_POLE_FILTERS>,
    svf_coeffs: ArrayVec<SvfCoeff, NUM_BANDS_PLUS_8>,

    needs_param_flush: bool,
    num_filters_changed: bool,
    lp_band_needs_param_sync: bool,
    hp_band_needs_param_sync: bool,
    bands_needing_param_sync: [bool; NUM_BANDS],

    sample_rate_recip: f64,
}

impl<const NUM_BANDS: usize, const NUM_BANDS_PLUS_8: usize>
    MeadowEqDspCoeff<NUM_BANDS, NUM_BANDS_PLUS_8>
{
    pub fn new(sample_rate: f64) -> Self {
        let sample_rate_recip = sample_rate.recip();

        Self {
            params: EqParams::default(),
            lp_band: MultiOrderBand::default(),
            hp_band: MultiOrderBand::default(),
            bands: [SecondOrderBand::default(); NUM_BANDS],
            one_pole_coeffs: ArrayVec::new(),
            svf_coeffs: ArrayVec::new(),
            needs_param_flush: false,
            num_filters_changed: false,
            lp_band_needs_param_sync: false,
            hp_band_needs_param_sync: false,
            bands_needing_param_sync: [false; NUM_BANDS],
            sample_rate_recip,
        }
    }

    pub fn params(&self) -> &EqParams<NUM_BANDS> {
        &self.params
    }

    pub fn set_params(&mut self, params: &EqParams<NUM_BANDS>) {
        if self.params.lp_band != params.lp_band {
            if self.params.lp_band.enabled != params.lp_band.enabled
                || self.params.lp_band.order != params.lp_band.order
            {
                self.num_filters_changed = true;
            }

            self.params.lp_band = params.lp_band;
            self.lp_band_needs_param_sync = true;
            self.needs_param_flush = true;
        }
        if self.params.hp_band != params.hp_band {
            if self.params.hp_band.enabled != params.hp_band.enabled
                || self.params.hp_band.order != params.hp_band.order
            {
                self.num_filters_changed = true;
            }

            self.params.hp_band = params.hp_band;
            self.hp_band_needs_param_sync = true;
            self.needs_param_flush = true;
        }

        for i in 0..NUM_BANDS {
            if self.params.bands[i] != params.bands[i] {
                if self.params.bands[i].enabled != params.bands[i].enabled {
                    self.num_filters_changed = true;
                }

                self.params.bands[i] = params.bands[i];
                self.bands_needing_param_sync[i] = true;
                self.needs_param_flush = true;
            }
        }
    }

    pub fn needs_param_flush(&self) -> bool {
        self.needs_param_flush
    }

    pub fn flush_param_changes(&mut self) -> Option<StateSyncInfo<NUM_BANDS>> {
        if !self.needs_param_flush {
            return None;
        }
        self.needs_param_flush = false;

        if self.num_filters_changed {
            self.num_filters_changed = false;
            self.one_pole_coeffs.clear();
            self.svf_coeffs.clear();
        }

        if self.lp_band_needs_param_sync {
            self.lp_band_needs_param_sync = false;

            self.lp_band.sync_params(
                &self.params.lp_band,
                self.sample_rate_recip,
                true,
                &mut self.one_pole_coeffs,
                &mut self.svf_coeffs,
            );
        }

        if self.hp_band_needs_param_sync {
            self.hp_band_needs_param_sync = false;

            self.hp_band.sync_params(
                &self.params.hp_band,
                self.sample_rate_recip,
                false,
                &mut self.one_pole_coeffs,
                &mut self.svf_coeffs,
            );
        }

        for band_i in 0..NUM_BANDS {
            if self.bands_needing_param_sync[band_i] {
                self.bands_needing_param_sync[band_i] = false;

                self.bands[band_i].sync_params(
                    &self.params.bands[band_i],
                    self.sample_rate_recip,
                    &mut self.svf_coeffs,
                );
            }
        }

        if self.num_filters_changed {
            Some(StateSyncInfo {
                lp_band_enabled: self.params.lp_band.enabled,
                lp_band_order: self.params.lp_band.order,
                hp_band_enabled: self.params.hp_band.enabled,
                hp_band_order: self.params.hp_band.order,
                bands_enabled: std::array::from_fn(|i| self.params.bands[i].enabled),
            })
        } else {
            None
        }
    }

    pub fn coeffs(
        &self,
    ) -> (
        &ArrayVec<OnePoleIirCoeff, MAX_ONE_POLE_FILTERS>,
        &ArrayVec<SvfCoeff, NUM_BANDS_PLUS_8>,
    ) {
        (&self.one_pole_coeffs, &self.svf_coeffs)
    }
}

#[derive(Default, Clone, Copy)]
struct SecondOrderBand {
    svf_filter_i: Option<usize>,
}

impl SecondOrderBand {
    fn sync_params<const NUM_BANDS_PLUS_8: usize>(
        &mut self,
        params: &BandParams,
        sample_rate_recip: f64,
        svf_filter_coeff: &mut ArrayVec<SvfCoeff, NUM_BANDS_PLUS_8>,
    ) {
        if !params.enabled {
            self.svf_filter_i = None;
            return;
        }

        let coeffs = match params.band_type {
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

        if let Some(i) = self.svf_filter_i {
            svf_filter_coeff[i] = coeffs;
        } else {
            self.svf_filter_i = Some(svf_filter_coeff.len());
            svf_filter_coeff.push(coeffs);
        }
    }
}

#[derive(Default)]
struct MultiOrderBand {
    order: FilterOrder,

    one_pole_iir_i: Option<usize>,
    svf_filter_i: Option<usize>,
}

impl MultiOrderBand {
    fn sync_params<const NUM_BANDS_PLUS_8: usize>(
        &mut self,
        params: &LpOrHpBandParams,
        sample_rate_recip: f64,
        is_lowpass: bool,
        one_pole_coeffs: &mut ArrayVec<OnePoleIirCoeff, MAX_ONE_POLE_FILTERS>,
        svf_coeffs: &mut ArrayVec<SvfCoeff, NUM_BANDS_PLUS_8>,
    ) {
        if !params.enabled {
            self.one_pole_iir_i = None;
            self.svf_filter_i = None;
            return;
        }

        self.order = params.order;

        match params.order {
            FilterOrder::X1 => {
                let coeffs = if is_lowpass {
                    OnePoleIirCoeffF64::lowpass(params.cutoff_hz as f64, sample_rate_recip).to_f32()
                } else {
                    OnePoleIirCoeffF64::highpass(params.cutoff_hz as f64, sample_rate_recip)
                        .to_f32()
                };

                if let Some(i) = self.one_pole_iir_i {
                    one_pole_coeffs[i] = coeffs;
                } else {
                    self.one_pole_iir_i = Some(one_pole_coeffs.len());
                    one_pole_coeffs.push(coeffs);
                }
            }
            FilterOrder::X2 => {
                let coeffs = if is_lowpass {
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

                if let Some(i) = self.svf_filter_i {
                    svf_coeffs[i] = coeffs;
                } else {
                    self.svf_filter_i = Some(svf_coeffs.len());
                    svf_coeffs.push(coeffs);
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

                if let Some(i) = self.svf_filter_i {
                    svf_coeffs[i] = coeffs[0].to_f32();
                    svf_coeffs[i + 1] = coeffs[1].to_f32();
                } else {
                    self.svf_filter_i = Some(svf_coeffs.len());
                    svf_coeffs.push(coeffs[0].to_f32());
                    svf_coeffs.push(coeffs[1].to_f32());
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

                if let Some(i) = self.svf_filter_i {
                    svf_coeffs[i] = coeffs[0].to_f32();
                    svf_coeffs[i + 1] = coeffs[1].to_f32();
                    svf_coeffs[i + 2] = coeffs[2].to_f32();
                } else {
                    self.svf_filter_i = Some(svf_coeffs.len());
                    svf_coeffs.push(coeffs[0].to_f32());
                    svf_coeffs.push(coeffs[1].to_f32());
                    svf_coeffs.push(coeffs[2].to_f32());
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

                if let Some(i) = self.svf_filter_i {
                    svf_coeffs[i] = coeffs[0].to_f32();
                    svf_coeffs[i + 1] = coeffs[1].to_f32();
                    svf_coeffs[i + 2] = coeffs[2].to_f32();
                    svf_coeffs[i + 3] = coeffs[3].to_f32();
                } else {
                    self.svf_filter_i = Some(svf_coeffs.len());
                    svf_coeffs.push(coeffs[0].to_f32());
                    svf_coeffs.push(coeffs[1].to_f32());
                    svf_coeffs.push(coeffs[2].to_f32());
                    svf_coeffs.push(coeffs[3].to_f32());
                }
            }
        }
    }
}

pub struct StateSyncInfo<const NUM_BANDS: usize> {
    pub lp_band_enabled: bool,
    pub lp_band_order: FilterOrder,
    pub hp_band_enabled: bool,
    pub hp_band_order: FilterOrder,

    pub bands_enabled: [bool; NUM_BANDS],
}

impl<const NUM_BANDS: usize> Default for StateSyncInfo<NUM_BANDS> {
    fn default() -> Self {
        Self {
            lp_band_enabled: false,
            lp_band_order: FilterOrder::X1,
            hp_band_enabled: false,
            hp_band_order: FilterOrder::X1,
            bands_enabled: [false; NUM_BANDS],
        }
    }
}
