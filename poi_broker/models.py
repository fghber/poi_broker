from . import db
from flask_login import UserMixin, current_user
from functools import wraps
from flask import flash, redirect, url_for
from datetime import datetime, timezone


"""
BEGIN;
CREATE INDEX IF NOT EXISTS idx_featuretable_date_alert_mjd
ON featuretable(date_alert_mjd);

CREATE INDEX IF NOT EXISTS idx_featuretable_alert_id
ON featuretable(alert_id);

CREATE INDEX IF NOT EXISTS idx_featuretable_ztf_object_id
ON featuretable(ztf_object_id);

CREATE INDEX IF NOT EXISTS idx_featuretable_locus_id
ON featuretable(locus_id);

CREATE INDEX IF NOT EXISTS idx_featuretable_ant_passband
ON featuretable(ant_passband);

CREATE INDEX IF NOT EXISTS idx_featuretable_locus_ra
ON featuretable(locus_ra);

CREATE INDEX IF NOT EXISTS idx_featuretable_locus_dec
ON featuretable(locus_dec);

CREATE INDEX IF NOT EXISTS idx_classification_alert_id
ON classification(alert_id);
COMMIT;

ANALYZE;

CREATE INDEX idx_featuretable_ant_mag_corrected
ON featuretable(ant_mag_corrected)
"""
class Ztf(db.Model):
    __tablename__ = 'featuretable'
    #id = db.Column(db.Integer, primary_key=True)
    date_log = db.Column(db.String)
    date_alert_mjd = db.Column(db.Float, primary_key=True)
    alert_id = db.Column(db.String, primary_key=True)
    ztf_object_id = db.Column(db.String)
    num_alerts = db.Column(db.Integer)
    num_mag_values = db.Column(db.Integer)
    brightest_alert_id_ztf = db.Column(db.String)
    brightest_alert_magnitude_ztf = db.Column(db.Float)
    brightest_alert_observation_time_ztf = db.Column(db.Float)
    locus_id = db.Column(db.String, primary_key=True)
    locus_ra = db.Column(db.Float)
    locus_dec = db.Column(db.Float)
    ant_mag_corrected = db.Column(db.Float)
    ant_passband = db.Column(db.String)
    feature_amplitude_magn_r = db.Column(db.Float)
    feature_anderson_darling_normal_magn_r = db.Column(db.Float)
    feature_beyond_1_std_magn_r = db.Column(db.Float)
    feature_beyond_2_std_magn_r = db.Column(db.Float)
    feature_cusum_magn_r = db.Column(db.Float)
    feature_eta_e_magn_r = db.Column(db.Float)
    feature_inter_percentile_range_2_magn_r = db.Column(db.Float)
    feature_inter_percentile_range_10_magn_r = db.Column(db.Float)
    feature_inter_percentile_range_25_magn_r = db.Column(db.Float)
    feature_kurtosis_magn_r = db.Column(db.Float)
    feature_linear_fit_slope_magn_r = db.Column(db.Float)
    feature_linear_fit_slope_sigma_magn_r = db.Column(db.Float)
    feature_linear_fit_reduced_chi2_magn_r = db.Column(db.Float)
    feature_linear_trend_magn_r = db.Column(db.Float)
    feature_linear_trend_sigma_magn_r = db.Column(db.Float)
    feature_magnitude_percentage_ratio_40_5_magn_r = db.Column(db.Float)
    feature_magnitude_percentage_ratio_20_5_magn_r = db.Column(db.Float)
    feature_maximum_slope_magn_r = db.Column(db.Float)
    feature_mean_magn_r = db.Column(db.Float)
    feature_median_absolute_deviation_magn_r = db.Column(db.Float)
    feature_percent_amplitude_magn_r = db.Column(db.Float)
    feature_percent_difference_magnitude_percentile_5_magn_r = db.Column(db.Float)
    feature_percent_difference_magnitude_percentile_10_magn_r = db.Column(db.Float)
    feature_median_buffer_range_percentage_10_magn_r = db.Column(db.Float)
    feature_median_buffer_range_percentage_20_magn_r = db.Column(db.Float)
    feature_period_0_magn_r = db.Column(db.Float)
    feature_period_s_to_n_0_magn_r = db.Column(db.Float)
    feature_period_1_magn_r = db.Column(db.Float)
    feature_period_s_to_n_1_magn_r = db.Column(db.Float)
    feature_period_2_magn_r = db.Column(db.Float)
    feature_period_s_to_n_2_magn_r = db.Column(db.Float)
    feature_period_3_magn_r = db.Column(db.Float)
    feature_period_s_to_n_3_magn_r = db.Column(db.Float)
    feature_period_4_magn_r = db.Column(db.Float)
    feature_period_s_to_n_4_magn_r = db.Column(db.Float)
    feature_periodogram_amplitude_magn_r = db.Column(db.Float)
    feature_periodogram_beyond_2_std_magn_r = db.Column(db.Float)
    feature_periodogram_beyond_3_std_magn_r = db.Column(db.Float)
    feature_periodogram_standard_deviation_magn_r = db.Column(db.Float)
    feature_chi2_magn_r = db.Column(db.Float)
    feature_skew_magn_r = db.Column(db.Float)
    feature_standard_deviation_magn_r = db.Column(db.Float)
    feature_stetson_k_magn_r = db.Column(db.Float)
    feature_weighted_mean_magn_r = db.Column(db.Float)
    feature_anderson_darling_normal_flux_r = db.Column(db.Float)
    feature_cusum_flux_r = db.Column(db.Float)
    feature_eta_e_flux_r = db.Column(db.Float)
    feature_excess_variance_flux_r = db.Column(db.Float)
    feature_kurtosis_flux_r = db.Column(db.Float)
    feature_mean_variance_flux_r = db.Column(db.Float)
    feature_chi2_flux_r = db.Column(db.Float)
    feature_skew_flux_r = db.Column(db.Float)
    feature_stetson_k_flux_r = db.Column(db.Float)
    feature_amplitude_magn_g = db.Column(db.Float)
    feature_anderson_darling_normal_magn_g = db.Column(db.Float)
    feature_beyond_1_std_magn_g = db.Column(db.Float)
    feature_beyond_2_std_magn_g = db.Column(db.Float)
    feature_cusum_magn_g = db.Column(db.Float)
    feature_eta_e_magn_g = db.Column(db.Float)
    feature_inter_percentile_range_2_magn_g = db.Column(db.Float)
    feature_inter_percentile_range_10_magn_g = db.Column(db.Float)
    feature_inter_percentile_range_25_magn_g = db.Column(db.Float)
    feature_kurtosis_magn_g = db.Column(db.Float)
    feature_linear_fit_slope_magn_g = db.Column(db.Float)
    feature_linear_fit_slope_sigma_magn_g = db.Column(db.Float)
    feature_linear_fit_reduced_chi2_magn_g = db.Column(db.Float)
    feature_linear_trend_magn_g = db.Column(db.Float)
    feature_linear_trend_sigma_magn_g = db.Column(db.Float)
    feature_magnitude_percentage_ratio_40_5_magn_g = db.Column(db.Float)
    feature_magnitude_percentage_ratio_20_5_magn_g = db.Column(db.Float)
    feature_maximum_slope_magn_g = db.Column(db.Float)
    feature_mean_magn_g = db.Column(db.Float)
    feature_median_absolute_deviation_magn_g = db.Column(db.Float)
    feature_percent_amplitude_magn_g = db.Column(db.Float)
    feature_percent_difference_magnitude_percentile_5_magn_g = db.Column(db.Float)
    feature_percent_difference_magnitude_percentile_10_magn_g = db.Column(db.Float)
    feature_median_buffer_range_percentage_10_magn_g = db.Column(db.Float)
    feature_median_buffer_range_percentage_20_magn_g = db.Column(db.Float)
    feature_period_0_magn_g = db.Column(db.Float)
    feature_period_s_to_n_0_magn_g = db.Column(db.Float)
    feature_period_1_magn_g = db.Column(db.Float)
    feature_period_s_to_n_1_magn_g = db.Column(db.Float)
    feature_period_2_magn_g = db.Column(db.Float)
    feature_period_s_to_n_2_magn_g = db.Column(db.Float)
    feature_period_3_magn_g = db.Column(db.Float)
    feature_period_s_to_n_3_magn_g = db.Column(db.Float)
    feature_period_4_magn_g = db.Column(db.Float)
    feature_period_s_to_n_4_magn_g = db.Column(db.Float)
    feature_periodogram_amplitude_magn_g = db.Column(db.Float)
    feature_periodogram_beyond_2_std_magn_g = db.Column(db.Float)
    feature_periodogram_beyond_3_std_magn_g = db.Column(db.Float)
    feature_periodogram_standard_deviation_magn_g = db.Column(db.Float)
    feature_chi2_magn_g = db.Column(db.Float)
    feature_skew_magn_g = db.Column(db.Float)
    feature_standard_deviation_magn_g = db.Column(db.Float)
    feature_stetson_k_magn_g = db.Column(db.Float)
    feature_weighted_mean_magn_g = db.Column(db.Float)
    feature_anderson_darling_normal_flux_g = db.Column(db.Float)
    feature_cusum_flux_g = db.Column(db.Float)
    feature_eta_e_flux_g = db.Column(db.Float)
    feature_excess_variance_flux_g = db.Column(db.Float)
    feature_kurtosis_flux_g = db.Column(db.Float)
    feature_mean_variance_flux_g = db.Column(db.Float)
    feature_chi2_flux_g = db.Column(db.Float)
    feature_skew_flux_g = db.Column(db.Float)
    feature_stetson_k_flux_g = db.Column(db.Float)
    anomaly_score = db.Column(db.Float)
    anomaly_mask = db.Column(db.String)
    anomaly_type = db.Column(db.String)
    is_corrected = db.Column(db.Boolean)
    g_r_max = db.Column(db.Float)
    g_r_mean = db.Column(db.Float)
    best_period = db.Column(db.Float)
    best_period_significance = db.Column(db.Float)
    best_period_g = db.Column(db.Float)
    best_period_r = db.Column(db.Float)
    best_period_i = db.Column(db.Float)
    power_rate_0p5 = db.Column(db.Float)
    power_rate_0p333 = db.Column(db.Float)
    power_rate_0p25 = db.Column(db.Float)
    power_rate_2 = db.Column(db.Float)
    power_rate_3 = db.Column(db.Float)
    power_rate_4 = db.Column(db.Float)
    mhps_ratio_g = db.Column(db.Float)
    mhps_ratio_R = db.Column(db.Float)
    mhps_ratio_i = db.Column(db.Float)
    mhps_low_g = db.Column(db.Float)
    mhps_low_R = db.Column(db.Float)
    mhps_low_i = db.Column(db.Float)
    mhps_high_g = db.Column(db.Float)
    mhps_high_R = db.Column(db.Float)
    mhps_high_i = db.Column(db.Float)
    mhps_non_zero_g = db.Column(db.Float)
    mhps_non_zero_R = db.Column(db.Float)
    mhps_non_zero_i = db.Column(db.Float)
    mhps_pn_flag_g = db.Column(db.Float)
    mhps_pn_flag_R = db.Column(db.Float)
    mhps_pn_flag_i = db.Column(db.Float)
    drw_omega_g = db.Column(db.Float)
    drw_omega_R = db.Column(db.Float)
    drw_omega_i = db.Column(db.Float)
    drw_tau_g = db.Column(db.Float)
    drw_tau_R = db.Column(db.Float)
    drw_tau_i = db.Column(db.Float)
    amplitude_g = db.Column(db.Float)
    amplitude_R = db.Column(db.Float)
    amplitude_i = db.Column(db.Float)
    anderson_darling_g = db.Column(db.Float)
    anderson_darling_R = db.Column(db.Float)
    anderson_darling_i = db.Column(db.Float)
    autocorrelation_length_g = db.Column(db.Float)
    autocorrelation_length_R = db.Column(db.Float)
    autocorrelation_length_i = db.Column(db.Float)
    beyond1std_g = db.Column(db.Float)
    beyond1std_R = db.Column(db.Float)
    beyond1std_i = db.Column(db.Float)
    con_g = db.Column(db.Float)
    con_R = db.Column(db.Float)
    con_i = db.Column(db.Float)
    eta_e_g = db.Column(db.Float)
    eta_e_R = db.Column(db.Float)
    eta_e_i = db.Column(db.Float)
    gskew_g = db.Column(db.Float)
    gskew_R = db.Column(db.Float)
    gskew_i = db.Column(db.Float)
    maxslope_g = db.Column(db.Float)
    maxslope_R = db.Column(db.Float)
    maxslope_i = db.Column(db.Float)
    meanmag_g = db.Column(db.Float)
    meanmag_R = db.Column(db.Float)
    meanmag_i = db.Column(db.Float)
    meanvariance_g = db.Column(db.Float)
    meanvariance_R = db.Column(db.Float)
    meanvariance_i = db.Column(db.Float)
    medianabsdev_g = db.Column(db.Float)
    medianabsdev_R = db.Column(db.Float)
    medianabsdev_i = db.Column(db.Float)
    medianbrp_g = db.Column(db.Float)
    medianbrp_R = db.Column(db.Float)
    medianbrp_i = db.Column(db.Float)
    pairslopetrend_g = db.Column(db.Float)
    pairslopetrend_R = db.Column(db.Float)
    pairslopetrend_i = db.Column(db.Float)
    percent_amplitude_g = db.Column(db.Float)
    percent_amplitude_R = db.Column(db.Float)
    percent_amplitude_i = db.Column(db.Float)
    q31_g = db.Column(db.Float)
    q31_R = db.Column(db.Float)
    q31_i = db.Column(db.Float)
    rcs_g = db.Column(db.Float)
    rcs_R = db.Column(db.Float)
    rcs_i = db.Column(db.Float)
    skew_g = db.Column(db.Float)
    skew_R = db.Column(db.Float)
    skew_i = db.Column(db.Float)
    smallkurtosis_g = db.Column(db.Float)
    smallkurtosis_R = db.Column(db.Float)
    smallkurtosis_i = db.Column(db.Float)
    stdmag_g = db.Column(db.Float)
    stdmag_R = db.Column(db.Float)
    stdmag_i = db.Column(db.Float)
    stetsonk_g = db.Column(db.Float)
    stetsonk_R = db.Column(db.Float)
    stetsonk_i = db.Column(db.Float)
    p_chi_g = db.Column(db.Float)
    p_chi_R = db.Column(db.Float)
    p_chi_i = db.Column(db.Float)
    ex_var_g = db.Column(db.Float)
    ex_var_R = db.Column(db.Float)
    ex_var_i = db.Column(db.Float)
    iar_phi_g = db.Column(db.Float)
    iar_phi_R = db.Column(db.Float)
    iar_phi_i = db.Column(db.Float)
    regression_slope_g = db.Column(db.Float)
    regression_slope_R = db.Column(db.Float)
    regression_slope_i = db.Column(db.Float)
    delta_mag_fid_g = db.Column(db.Float)
    delta_mag_fid_R = db.Column(db.Float)
    delta_mag_fid_i = db.Column(db.Float)
    harmonics_mag_1_g = db.Column(db.Float)
    harmonics_mag_1_R = db.Column(db.Float)
    harmonics_mag_1_i = db.Column(db.Float)
    harmonics_mag_2_g = db.Column(db.Float)
    harmonics_mag_2_R = db.Column(db.Float)
    harmonics_mag_2_i = db.Column(db.Float)
    harmonics_mag_3_g = db.Column(db.Float)
    harmonics_mag_3_R = db.Column(db.Float)
    harmonics_mag_3_i = db.Column(db.Float)
    harmonics_mag_4_g = db.Column(db.Float)
    harmonics_mag_4_R = db.Column(db.Float)
    harmonics_mag_4_i = db.Column(db.Float)
    harmonics_mag_5_g = db.Column(db.Float)
    harmonics_mag_5_R = db.Column(db.Float)
    harmonics_mag_5_i = db.Column(db.Float)
    harmonics_mag_6_g = db.Column(db.Float)
    harmonics_mag_6_R = db.Column(db.Float)
    harmonics_mag_6_i = db.Column(db.Float)
    harmonics_mag_7_g = db.Column(db.Float)
    harmonics_mag_7_R = db.Column(db.Float)
    harmonics_mag_7_i = db.Column(db.Float)
    harmonics_phi_1_g = db.Column(db.Float)
    harmonics_phi_1_R = db.Column(db.Float)
    harmonics_phi_1_i = db.Column(db.Float)
    harmonics_phi_2_g = db.Column(db.Float)
    harmonics_phi_2_R = db.Column(db.Float)
    harmonics_phi_2_i = db.Column(db.Float)
    harmonics_phi_3_g = db.Column(db.Float)
    harmonics_phi_3_R = db.Column(db.Float)
    harmonics_phi_3_i = db.Column(db.Float)
    harmonics_phi_4_g = db.Column(db.Float)
    harmonics_phi_4_R = db.Column(db.Float)
    harmonics_phi_4_i = db.Column(db.Float)
    harmonics_phi_5_g = db.Column(db.Float)
    harmonics_phi_5_R = db.Column(db.Float)
    harmonics_phi_5_i = db.Column(db.Float)
    harmonics_phi_6_g = db.Column(db.Float)
    harmonics_phi_6_R = db.Column(db.Float)
    harmonics_phi_6_i = db.Column(db.Float)
    harmonics_phi_7_g = db.Column(db.Float)
    harmonics_phi_7_R = db.Column(db.Float)
    harmonics_phi_7_i = db.Column(db.Float)
    harmonics_mse_g = db.Column(db.Float)
    harmonics_mse_R = db.Column(db.Float)
    harmonics_mse_i = db.Column(db.Float)
    harmonics_chi_per_degree_g = db.Column(db.Float)
    harmonics_chi_per_degree_R = db.Column(db.Float)
    harmonics_chi_per_degree_i = db.Column(db.Float)

    # @property
    # def ra(self):
    #     ra = shape.to_shape(self.location).x
    #     if ra < 0:
    #         ra = ra + 360
    #     return ra
    # @property
    # def dec(self):
    #     return shape.to_shape(self.location).y

    # Relationship to Classification table
    classification = db.relationship("Classification", foreign_keys="[Classification.alert_id]", primaryjoin="Ztf.alert_id==Classification.alert_id", uselist=False, viewonly=True)

    def __str__(self):
        return self.ztf_object_id

class Crossmatches(db.Model):
    __tablename__ = 'crossmatches'
    id = db.Column(db.Integer, primary_key=True)
    locus_id = db.Column(db.String)
    catalog = db.Column(db.String)
    object = db.Column(db.String)
    ra_cat = db.Column(db.Float)
    dec_cat = db.Column(db.Float)
    separation = db.Column(db.Float)


class Classification(db.Model):
    __tablename__ = 'classification'

    alert_id = db.Column(db.String, primary_key=True)
    p_cvnova = db.Column(db.Float)
    p_e = db.Column(db.Float)
    p_lpv = db.Column(db.Float)
    p_puls = db.Column(db.Float)
    p_periodic_other = db.Column(db.Float)
    p_quas = db.Column(db.Float)
    p_sn = db.Column(db.Float)
    p_yso = db.Column(db.Float)
    prob_class = db.Column(db.String)

# class Gaiadr3_variability(db.Model):
#     __tablename__ = 'gaiadr3_variability'
#     source_id = db.Column(db.Integer, primary_key=True) #ot a real PK
#     ra = db.Column(db.Float)
#     dec = db.Column(db.Float)
#     phot_g_mean_mag = db.Column(db.Float)
#     phot_rp_mean_mag = db.Column(db.Float)
#     in_vari_classification_result = db.Column(db.Integer)
#     in_vari_rrlyrae = db.Column(db.Integer)
#     in_vari_cepheid = db.Column(db.Integer)
#     in_vari_planetary_transit = db.Column(db.Integer)
#     in_vari_short_timescale = db.Column(db.Integer)
#     in_vari_long_period_variable = db.Column(db.Integer)
#     in_vari_eclipsing_binary = db.Column(db.Integer)
#     in_vari_rotation_modulation = db.Column(db.Integer)
#     in_vari_ms_oscillator = db.Column(db.Integer)
#     in_vari_agn = db.Column(db.Integer)
#     in_vari_microlensing = db.Column(db.Integer)
#     in_vari_compact_companion = db.Column(db.Integer)

"""

TODO: Add table initialization scripts here or in a separate file, or use Flask-Migrate for migrations.
CREATE TABLE user (
    id INTEGER PRIMARY KEY,
    email VARCHAR(100) NOT NULL UNIQUE,
    password VARCHAR(100) NOT NULL,
    name VARCHAR(1000),
    role VARCHAR(20) NOT NULL DEFAULT 'user',
    email_verified INTEGER DEFAULT 0,
    email_verification_token VARCHAR(128),
    reset_token VARCHAR(128),
    reset_token_expires INTEGER,
    created_at INTEGER DEFAULT (strftime('%s', 'now'))
);

CREATE INDEX ix_user_email_verification_token ON user(email_verification_token);
CREATE INDEX ix_user_reset_token ON user(reset_token);

-- Add new columns to user table
ALTER TABLE user ADD COLUMN email_verified INTEGER DEFAULT 0;
ALTER TABLE user ADD COLUMN email_verification_token VARCHAR(128);

-- Index email_verification_token for fast lookups
CREATE INDEX IF NOT EXISTS ix_user_email_verification_token 
ON user(email_verification_token);
"""
class User(UserMixin, db.Model):
    __bind_key__ = 'users' # route this model to the users DB
    __tablename__ = 'user'

    id = db.Column(db.Integer, primary_key=True) # primary keys are required by SQLAlchemy
    email = db.Column(db.String(100), unique=True)
    password = db.Column(db.String(100))
    name = db.Column(db.String(1000))
    role = db.Column(db.String(20), default='user')    # Email verification
    email_verified = db.Column(db.Boolean, default=False, nullable=False)
    email_verification_token = db.Column(db.String(128), index=True, nullable=True)    # Used by the reset-password flow
    reset_token = db.Column(db.String(128), index=True, nullable=True)
    reset_token_expires = db.Column(db.Integer, nullable=True) # epoch seconds
    # created_at = db.Column(db.DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    # last_password_changed = db.Column(db.DateTime(timezone=True), nullable=True)

    def __repr__(self):
        return f"<User {self.email}>"
   
    def has_role(self, role):
        return self.role == role
    
    def is_admin(self):
        return self.role == 'admin'

# Custom decorator for role-based access
def role_required(role):
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            if not current_user.is_authenticated or not current_user.has_role(role):
                flash('Access denied. Insufficient permissions.')
                return redirect(url_for('main.start'))
            return f(*args, **kwargs)
        return decorated_function
    return decorator


"""
CREATE TABLE favorite_group (
  id INTEGER PRIMARY KEY,
  user_id INTEGER NOT NULL,
  name TEXT NOT NULL,
  created_at TEXT NOT NULL,
  UNIQUE(user_id, name),
  FOREIGN KEY(user_id) REFERENCES user(id) ON DELETE CASCADE
);
CREATE INDEX ix_favorite_group_user_id ON favorite_group(user_id);
"""
class FavoriteGroup(db.Model):
    __bind_key__ = 'users'
    __tablename__ = 'favorite_group'

    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', ondelete='CASCADE'), nullable=False, index=True)
    name = db.Column(db.String(128), nullable=False)
    created_at = db.Column(db.DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    __table_args__ = (db.UniqueConstraint('user_id', 'name', name='uix_user_group_name'),)

    def __repr__(self):
        return f"<FavoriteGroup {self.name}>"

"""
CREATE TABLE favorite (
  id INTEGER PRIMARY KEY,
  user_id INTEGER NOT NULL,
  locus_id TEXT NOT NULL,
  group_id INTEGER,
  created_at TEXT NOT NULL,
  UNIQUE(user_id, locus_id),
  FOREIGN KEY(user_id) REFERENCES user(id) ON DELETE CASCADE,
  FOREIGN KEY(group_id) REFERENCES favorite_group(id) ON DELETE SET NULL
);
CREATE INDEX ix_favorite_user_id ON favorite(user_id);
CREATE INDEX ix_favorite_locus_id ON favorite(locus_id);
CREATE INDEX ix_favorite_group_id ON favorite(group_id);
"""
class Favorite(db.Model):
    __bind_key__ = 'users'
    __tablename__ = 'favorite'

    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', ondelete='CASCADE'), nullable=False, index=True)
    locus_id = db.Column(db.String(128), nullable=False, index=True)
    group_id = db.Column(db.Integer, db.ForeignKey('favorite_group.id', ondelete='SET NULL'), nullable=True, index=True)
    created_at = db.Column(db.DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    __table_args__ = (db.UniqueConstraint('user_id', 'locus_id', name='uix_user_locus'),)

    def __repr__(self):
        return f"<Favorite {self.locus_id}>"

# DDL for watchlist table (run once against users.db):
# CREATE TABLE watchlist (
#     id INTEGER PRIMARY KEY,
#     user_id INTEGER NOT NULL,
#     name TEXT NOT NULL,
#     rules_json TEXT NOT NULL,
#     sql_where TEXT NOT NULL,
#     created_at INTEGER NOT NULL DEFAULT (strftime('%s', 'now')),
#     UNIQUE(user_id, name),
#     FOREIGN KEY(user_id) REFERENCES user(id) ON DELETE CASCADE
# );
# CREATE INDEX ix_watchlist_user_id ON watchlist(user_id);
class Watchlist(db.Model):
    __bind_key__ = 'users'
    __tablename__ = 'watchlist'

    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', ondelete='CASCADE'), nullable=False, index=True)
    name = db.Column(db.String(128), nullable=False)
    rules_json = db.Column(db.Text, nullable=False)
    sql_where = db.Column(db.Text, nullable=False)
    created_at = db.Column(db.Integer, nullable=False)

    __table_args__ = (db.UniqueConstraint('user_id', 'name', name='uix_watchlist_user_name'),)

    def __repr__(self):
        return f"<Watchlist {self.name}>"

"""
DDL for filter_bookmark (run once against users.db):
CREATE TABLE filter_bookmark (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    user_id INTEGER NOT NULL,
    name VARCHAR(200) NOT NULL,
    query_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    FOREIGN KEY (user_id) REFERENCES user(id) ON DELETE CASCADE
);
CREATE INDEX ix_filter_bookmark_user_id ON filter_bookmark(user_id);
"""
class FilterBookmark(db.Model):
    __bind_key__ = 'users'
    __tablename__ = 'filter_bookmark'

    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', ondelete='CASCADE'), nullable=False, index=True)
    name = db.Column(db.String(200), nullable=False)
    query_json = db.Column(db.Text, nullable=False)
    created_at = db.Column(db.DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    def __repr__(self):
        return f"<FilterBookmark {self.name!r}>"