pub mod coeff;
pub mod state;
pub mod stereo;

pub const DEFAULT_Q: f32 = meadow_dsp_mit::filter::svf::f64::Q_BUTTERWORTH_ORD2 as f32;

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FilterOrder {
    #[default]
    X1 = 0,
    X2,
    X4,
    X6,
    X8,
}

impl FilterOrder {
    pub fn from_u32(v: u32) -> Self {
        match v {
            0 => Self::X1,
            1 => Self::X2,
            2 => Self::X4,
            3 => Self::X6,
            _ => Self::X8,
        }
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BandType {
    #[default]
    Bell = 0,
    LowShelf,
    HighShelf,
    Notch,
    Allpass,
}

impl BandType {
    pub fn from_u32(v: u32) -> Self {
        match v {
            0 => Self::Bell,
            1 => Self::LowShelf,
            2 => Self::HighShelf,
            3 => Self::Notch,
            _ => Self::Allpass,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandParams {
    pub enabled: bool,
    pub band_type: BandType,
    pub cutoff_hz: f32,
    pub q: f32,
    pub gain_db: f32,
}

impl Default for BandParams {
    fn default() -> Self {
        Self {
            enabled: false,
            band_type: BandType::Notch,
            cutoff_hz: 1000.0,
            q: DEFAULT_Q as f32,
            gain_db: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LpOrHpBandParams {
    pub enabled: bool,
    pub cutoff_hz: f32,
    pub q: f32,
    pub order: FilterOrder,
}

impl Default for LpOrHpBandParams {
    fn default() -> Self {
        Self {
            enabled: false,
            cutoff_hz: 21_480.0,
            q: DEFAULT_Q as f32,
            order: FilterOrder::X2,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EqParams<const NUM_BANDS: usize> {
    pub lp_band: LpOrHpBandParams,
    pub hp_band: LpOrHpBandParams,

    pub bands: [BandParams; NUM_BANDS],
}

impl<const NUM_BANDS: usize> Default for EqParams<NUM_BANDS> {
    fn default() -> Self {
        Self {
            lp_band: LpOrHpBandParams {
                cutoff_hz: 20.0,
                ..Default::default()
            },
            hp_band: LpOrHpBandParams::default(),
            bands: [BandParams::default(); NUM_BANDS],
        }
    }
}
