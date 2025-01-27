import math
import typing

KM_IN_A_PC = 3.086e13
KM_IN_AU = 1.496e8
DEG_IN_A_CIRC = 360
RAD_IN_A_CIRC = 2*math.pi
MIN_IN_A_H = 60
SEC_IN_A_MIN = 60
SEC_IN_A_H = 3600
DEC_IN_A_H = 15
H_IN_A_CIRC = 24
SEC_IN_A_YEAR = 3.154e7
RR_LYRAE_M = 0.75
SPEED_OF_LIGHT_KMS = 3.0e5
HUBBLE_CONSTANT_KMS = 70
PLANCK_H = 6.63e-34
J_IN_A_eV = 1.602e-19
SOLAR_RADIUS_KM = 696340


def to_sf(sf: int, res: float) -> str:
    """To N significan figures"""
    return f"%.{sf}g" % res


def to_dp(dp: int, res: float) -> str:
    """To N decimal places"""
    return f"%.{dp}f" % res


def format_result(
        res: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None,
    ) -> str:
    if dp is not None:
        return to_dp(dp=dp, res=res)
    elif sf is not None:
        return to_sf(sf=sf, res=res)
    else:
        return f"{res}"
    

def pc_to_au(pc: float, sf: typing.Optional[int] = None, dp: typing.Optional[int] = None) -> str:
    """Convert parsecs to astronomical units"""

    return format_result(res=( pc * KM_IN_A_PC ) / KM_IN_AU, sf=sf, dp=dp)


def to_radians(
        degs: typing.Optional[float] = 0.0, 
        arcmins: typing.Optional[float] = 0.0, 
        arcsecs: typing.Optional[float] = 0.0,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    '''Convert degrees, arcminutes, and arcsedons to radians'''
    
    degs_to_radians = (degs / DEG_IN_A_CIRC) * RAD_IN_A_CIRC
    arcmins_to_radians = (( arcmins / MIN_IN_A_H ) / DEG_IN_A_CIRC ) * RAD_IN_A_CIRC
    arcsecs_to_radians = (( arcsecs / SEC_IN_A_H ) / DEG_IN_A_CIRC ) * RAD_IN_A_CIRC

    res = degs_to_radians + arcmins_to_radians + arcsecs_to_radians

    return format_result(res=res, sf=sf, dp=dp)


def from_radians(
        radians: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    '''Convert radians to degrees, arcminutes, and arcseconds'''
    
    degrees_f = ( radians / RAD_IN_A_CIRC ) * DEG_IN_A_CIRC
    degrees_f, degrees = math.modf(degrees_f)
    arcmins_f = degrees_f * MIN_IN_A_H
    arcmins_f, arcmins = math.modf(arcmins_f)
    arcsecs = arcmins_f * SEC_IN_A_MIN

    return (
        format_result(res=degrees, sf=sf, dp=dp), 
        format_result(res=arcmins, sf=sf, dp=dp), 
        format_result(res=arcsecs, sf=sf, dp=dp)
    )


def ra_from_time_to_deg(
        h: typing.Optional[float] = 0.0,
        m: typing.Optional[float] = 0.0,
        s: typing.Optional[float] = 0.0,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None 
    ) -> str:
    '''Convert right ascension coordinates 
    in time to decimal degrees'''

    total_h = h + ( m / MIN_IN_A_H ) + ( s / SEC_IN_A_H )
    total_deg = total_h * DEC_IN_A_H

    return format_result(res=total_deg, sf=sf, dp=dp)


def angular_to_linear(
        theta_radians: float,
        distance_away_pc: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    '''Convert angular distance or size to linear'''
    res = math.sin(theta_radians) * distance_away_pc
    return format_result(res=res, sf=sf, dp=dp)


def magnitude_from_components(
        ra_comp = float,
        dec_comp = float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    """Compute magnitude using Pythagoras"""
    res = math.sqrt(ra_comp**2 + dec_comp**2)
    return format_result(res=res, sf=sf, dp=dp)


def apparent_to_absolute_magnitude(
        apparent_magnitude: float,
        distance_away_in_pc: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None,
    ) -> str:
    """Compute absolute magnitude from apparent magnitude via distance"""
    res = apparent_magnitude - 5 * math.log10(distance_away_in_pc) + 5
    return format_result(res=res, sf=sf, dp=dp)


def distance_from_magnitudes(
        apparent_magnitude: float,
        absolute_magnitude: float,
        interstellar_extinction: typing.Optional[float] = 0.0,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None,
    ) -> str:
    """Compute distance away from apparent and absolute magnitudes"""
    res = 10**(
        (apparent_magnitude - absolute_magnitude + 5 - interstellar_extinction) / 5
    )
    return format_result(res=res, sf=sf, dp=dp)


def magnitude_difference_from_fluxlum_ratio(
        fluxlum_ratio: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    """Calculate magnitude difference from relative fluxlum"""
    res = 2.5 * math.log10(fluxlum_ratio)
    return format_result(res=res, sf=sf, dp=dp)


def magnitude_from_fluxlum_ratio(
        magnitude_b: float,
        b_a_flux_ratio: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    """Calculate apparent magnitude from relative 
    flux and another magnitude
    """
    apparent_magnitude_difference = float(
        apparent_magnitude_difference_from_flux_ratio(
            flux_ratio=b_a_flux_ratio,
            sf=sf,
            dp=dp
        )
    )
    res = apparent_magnitude_difference + magnitude_b
    return format_result(res=res, sf=sf, dp=dp)


def flux_ratio_from_apparent_magnitude_difference(
        magnitude_a_b_difference: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    """"Calculate B/A flux ratio from A/B apparent magnitude difference"""
    res = 10**(magnitude_a_b_difference/2.5)
    return format_result(res=res, sf=sf, dp=dp)


def flux_ratio_from_apparent_magnitudes(
        magnitude_a: float,
        magnitude_b: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
    ) -> str:
    """""Calculate B/A flux ratio from two apparent magnitudes"""
    return flux_ratio_from_apparent_magnitude_difference(
        magnitude_a_b_difference=magnitude_a-magnitude_b,
        sf=sf,
        dp=dp
    )


def distance_from_parallax(
        parallax_angle_arcsec: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate distance in parsecs from parallax angle in arcseconds"""
    res = 1.0 / parallax_angle_arcsec
    return format_result(res=res, sf=sf, dp=dp)


def parallax_from_distance(
        distance_away: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate parallax in arcseconds from distance in parsecs"""
    res = 1.0 / distance_away
    return format_result(res=res, sf=sf, dp=dp)


def distance_from_rr_lyrae(
        apparent_magnitude: float,
        interstellar_extinction: typing.Optional[float] = 0.0,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate distance to an RR lyrae from its apparent magnitude"""
    res = 10**(
        (apparent_magnitude - RR_LYRAE_M + 5 - interstellar_extinction) / 5
    )
    return format_result(res=res, sf=sf, dp=dp)


def distance_from_cepheid(
        apparent_magnitude: float,
        period_days: float,
        interstellar_extinction: typing.Optional[float] = 0.0,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate distance to a Cepheid variable from its apparent
    magnitude and period"""
    absolute_magnitude = float(format_result(
        res=-2.43 * math.log10(period_days) - 1.62,
        sf=sf,
        dp=dp
    ))
    res = distance_from_magnitudes(
        apparent_magnitude=apparent_magnitude,
        absolute_magnitude=absolute_magnitude,
        interstellar_extinction=interstellar_extinction,
        sf=sf,
        dp=dp
    )
    return format_result(res=float(res), sf=sf, dp=dp)


def distance_from_redshift(
        redshift: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate distance from redshift and speed of recession"""
    speed_of_recession = float(format_result(
        res=redshift * SPEED_OF_LIGHT_KMS,
        sf=sf,
        dp=dp
    ))
    distance_away = float(format_result(
        res=speed_of_recession / HUBBLE_CONSTANT_KMS,
        sf=sf,
        dp=dp
    ))
    return (
        format_result(res=speed_of_recession, sf=sf, dp=dp),
        format_result(res=distance_away, sf=sf, dp=dp)
    )


def energy_from_wavelength(
        wavelength: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate a photon's energy in eV from its wavelength in meters"""
    res = PLANCK_H * SPEED_OF_LIGHT_KMS*1e3 / wavelength
    return format_result(res=res, sf=sf, dp=dp)


def wiens_law(
        input: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Use Wien's law to compute Temeprature or peak Wavelength"""
    res = 2.9e-3/input
    return format_result(res=res, sf=sf, dp=dp)


def radius_from_lt_relationship(
        luminosity: float,
        temperature: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Use LTR relationship to determine radius, in Solar units"""
    res = math.sqrt(luminosity) * (1/temperature)**2
    return format_result(res=res, sf=sf, dp=dp)


def luminosity_from_rt_relationship(
        radius: float,
        temperature: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Use LTR relationship to determine luminosity, in Solar units"""
    res = radius**2 * temperature**4
    return format_result(res=res, sf=sf, dp=dp)


def temperature_from_lr_relationship(
        luminosity: float,
        radius: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Use LTR relationship to determine temperature, in Solar units"""
    res = (luminosity*radius**2)**(1/4)
    return format_result(res=res, sf=sf, dp=dp)


def hydrogen_transition_wavelength(
        final_state: int,
        initial_state: int,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """"""
    initial_state_energy = float(format_result(
        res=-13.6/(initial_state**2),
        sf=sf,
        dp=dp
    ))
    final_state_energy = float(format_result(
        res=-13.6/(final_state**2),
        sf=sf,
        dp=dp
    ))
    res = ((PLANCK_H/J_IN_A_eV) * SPEED_OF_LIGHT_KMS * 1e3) \
        / (final_state_energy - initial_state_energy)
    return format_result(res=res, sf=sf, dp=dp)


def radial_velocity_from_doppler_shift(
        observed_wavelength: float,
        emitted_wavelength: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Compute radial velocity from doppler shift, in meters per second"""
    res = (SPEED_OF_LIGHT_KMS*1e3)*(observed_wavelength-emitted_wavelength) \
        / emitted_wavelength
    return format_result(res=res, sf=sf, dp=dp)


def doppler_shift_from_radial_velocity(
        radial_velocity: float,
        emitted_wavelength: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Compute doppler shift from radial veolocity"""
    res = radial_velocity*emitted_wavelength / (SPEED_OF_LIGHT_KMS*1e3)
    return format_result(res=res, sf=sf, dp=dp)


def emc2(
        mass: float, 
        sf: typing.Optional[int] = None, 
        dp: typing.Optional[int] = None
        ) -> str:
    """Compute energy from mass"""
    res = mass*(SPEED_OF_LIGHT_KMS*1e3)**2
    return format_result(res=res, sf=sf, dp=dp)

