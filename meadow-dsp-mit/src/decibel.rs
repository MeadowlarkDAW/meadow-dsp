pub mod f32 {
    /// Returns the raw amplitude from the given decibel value.
    #[inline]
    pub fn db_to_amp(db: f32) -> f32 {
        if db == f32::NEG_INFINITY {
            0.0
        } else {
            10.0f32.powf(0.05 * db)
        }
    }

    /// Returns the decibel value from the given raw amplitude.
    #[inline]
    pub fn amp_to_db(amp: f32) -> f32 {
        if amp == 0.0 {
            f32::NEG_INFINITY
        } else {
            20.0 * amp.log10()
        }
    }

    /// Returns the raw amplitude from the given decibel value.
    ///
    /// If `db == f32::NEG_INFINITY || db <= db_epsilon`, then `0.0` (silence) will be
    /// returned.
    #[inline]
    pub fn db_to_amp_clamped(db: f32, db_epsilon: f32) -> f32 {
        if db == f32::NEG_INFINITY || db <= db_epsilon {
            0.0
        } else {
            db_to_amp(db)
        }
    }

    /// Returns the decibel value from the given raw amplitude.
    ///
    /// If `amp <= amp_epsilon`, then `f32::NEG_INFINITY` (silence) will be returned.
    #[inline]
    pub fn amp_to_db_clamped(amp: f32, amp_epsilon: f32) -> f32 {
        if amp <= amp_epsilon {
            f32::NEG_INFINITY
        } else {
            amp_to_db(amp)
        }
    }

    /// Map the linear volume (where `0.0` means mute and `1.0` means unity
    /// gain) to the corresponding raw amplitude value (not decibels) for use in
    /// DSP. Values above `1.0` are allowed.
    ///
    /// If the resulting amplitude is `<= amp_epsilon`, then `0.0` (silence) will be
    /// returned.
    #[inline]
    pub fn linear_volume_to_amp_clamped(linear_volume: f32, amp_epsilon: f32) -> f32 {
        let v = linear_volume * linear_volume;
        if v <= amp_epsilon {
            0.0
        } else {
            v
        }
    }

    /// Map the raw amplitude (where `0.0` means mute and `1.0` means unity
    /// gain) to the corresponding linear volume.
    ///
    /// If the amplitude is `<= amp_epsilon`, then `0.0` (silence) will be
    /// returned.
    #[inline]
    pub fn amp_to_linear_volume_clamped(amp: f32, amp_epsilon: f32) -> f32 {
        if amp <= amp_epsilon {
            0.0
        } else {
            amp.sqrt()
        }
    }
}

pub mod f64 {
    /// Returns the raw amplitude from the given decibel value.
    #[inline]
    pub fn db_to_amp(db: f64) -> f64 {
        if db == f64::NEG_INFINITY {
            0.0
        } else {
            10.0f64.powf(0.05 * db)
        }
    }

    /// Returns the decibel value from the given raw amplitude.
    #[inline]
    pub fn amp_to_db(amp: f64) -> f64 {
        if amp == 0.0 {
            f64::NEG_INFINITY
        } else {
            20.0 * amp.log10()
        }
    }

    /// Returns the raw amplitude from the given decibel value.
    ///
    /// If `db == f64::NEG_INFINITY || db <= db_epsilon`, then `0.0` (silence) will be
    /// returned.
    #[inline]
    pub fn db_to_amp_clamped(db: f64, db_epsilon: f64) -> f64 {
        if db == f64::NEG_INFINITY || db <= db_epsilon {
            0.0
        } else {
            db_to_amp(db)
        }
    }

    /// Returns the decibel value from the given raw amplitude.
    ///
    /// If `amp <= amp_epsilon`, then `f64::NEG_INFINITY` (silence) will be returned.
    #[inline]
    pub fn amp_to_db_clamped(amp: f64, amp_epsilon: f64) -> f64 {
        if amp <= amp_epsilon {
            f64::NEG_INFINITY
        } else {
            amp_to_db(amp)
        }
    }

    /// Map the linear volume (where `0.0` means mute and `1.0` means unity
    /// gain) to the corresponding raw amplitude value (not decibels) for use in
    /// DSP. Values above `1.0` are allowed.
    ///
    /// If the resulting amplitude is `<= amp_epsilon`, then `0.0` (silence) will be
    /// returned.
    #[inline]
    pub fn linear_volume_to_amp_clamped(linear_volume: f64, amp_epsilon: f64) -> f64 {
        let v = linear_volume * linear_volume;
        if v <= amp_epsilon {
            0.0
        } else {
            v
        }
    }

    /// Map the raw amplitude (where `0.0` means mute and `1.0` means unity
    /// gain) to the corresponding linear volume.
    ///
    /// If the amplitude is `<= amp_epsilon`, then `0.0` (silence) will be
    /// returned.
    #[inline]
    pub fn amp_to_linear_volume_clamped(amp: f64, amp_epsilon: f64) -> f64 {
        if amp <= amp_epsilon {
            0.0
        } else {
            amp.sqrt()
        }
    }
}
