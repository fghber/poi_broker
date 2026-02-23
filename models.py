from . import db
from flask_login import UserMixin, current_user
from functools import wraps
from flask import flash, redirect, url_for

class Ztf(db.Model):
    __tablename__ = 'featuretable'
    #id = db.Column(db.Integer, primary_key=True)
    date_alert_mjd = db.Column(db.Float, primary_key=True)
    alert_id = db.Column(db.String, primary_key=True)
    ztf_object_id = db.Column(db.String)
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

    # @property
    # def ra(self):
    #     ra = shape.to_shape(self.location).x
    #     if ra < 0:
    #         ra = ra + 360
    #     return ra
    # @property
    # def dec(self):
    #     return shape.to_shape(self.location).y

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

class Gaiadr3_variability(db.Model):
    __tablename__ = 'gaiadr3_variability'
    source_id = db.Column(db.Integer, primary_key=True) #ot a real PK
    ra = db.Column(db.Float)
    dec = db.Column(db.Float)
    phot_g_mean_mag = db.Column(db.Float)
    phot_rp_mean_mag = db.Column(db.Float)
    in_vari_classification_result = db.Column(db.Integer)
    in_vari_rrlyrae = db.Column(db.Integer)
    in_vari_cepheid = db.Column(db.Integer)
    in_vari_planetary_transit = db.Column(db.Integer)
    in_vari_short_timescale = db.Column(db.Integer)
    in_vari_long_period_variable = db.Column(db.Integer)
    in_vari_eclipsing_binary = db.Column(db.Integer)
    in_vari_rotation_modulation = db.Column(db.Integer)
    in_vari_ms_oscillator = db.Column(db.Integer)
    in_vari_agn = db.Column(db.Integer)
    in_vari_microlensing = db.Column(db.Integer)
    in_vari_compact_companion = db.Column(db.Integer)


class User(UserMixin, db.Model):
    __bind_key__ = 'users' # route this model to the users DB
    __tablename__ = 'user'

    id = db.Column(db.Integer, primary_key=True) # primary keys are required by SQLAlchemy
    email = db.Column(db.String(100), unique=True)
    password = db.Column(db.String(100))
    name = db.Column(db.String(1000))
    role = db.Column(db.String(20), default='user')
    # Used by the reset-password flow
    reset_token = db.Column(db.String(128), index=True, nullable=True)
    reset_token_expires = db.Column(db.DateTime(timezone=True), nullable=True)

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
