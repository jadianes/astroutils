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
    

def pc_to_au(pc: float) -> float:
    """Convert parsecs to astronomical units"""

    return ( pc * KM_IN_A_PC ) / KM_IN_AU


def to_radians(
        sf: int,
        degs: typing.Optional[float] = 0.0, 
        arcmins: typing.Optional[float] = 0.0, 
        arcsecs: typing.Optional[float] = 0.0,
    ) -> str:
    '''Convert degrees, arcminutes, and arcsedons to radians'''
    
    degs_to_radians = (degs / DEG_IN_A_CIRC) * RAD_IN_A_CIRC
    arcmins_to_radians = (( arcmins / MIN_IN_A_H ) / DEG_IN_A_CIRC ) * RAD_IN_A_CIRC
    arcsecs_to_radians = (( arcsecs / SEC_IN_A_H ) / DEG_IN_A_CIRC ) * RAD_IN_A_CIRC

    res = degs_to_radians + arcmins_to_radians + arcsecs_to_radians

    return to_sf(sf=sf, res=res)


def from_radians(
        sf: int,
        radians: float
    ) -> str:
    '''Convert radians to degrees, arcminutes, and arcseconds'''
    
    degrees_f = ( radians / RAD_IN_A_CIRC ) * DEG_IN_A_CIRC
    degrees_f, degrees = math.modf(degrees_f)
    arcmins_f = degrees_f * MIN_IN_A_H
    arcmins_f, arcmins = math.modf(arcmins_f)
    arcsecs = arcmins_f * SEC_IN_A_MIN

    return (
        to_sf(sf=sf, res=degrees), 
        to_sf(sf=sf, res=arcmins), 
        to_sf(sf=sf, res=arcsecs)
    )


def ra_from_time_to_deg(
        sf: float, 
        h: typing.Optional[float] = 0.0,
        m: typing.Optional[float] = 0.0,
        s: typing.Optional[float] = 0.0
    ) -> str:
    '''Convert right ascension coordinates 
    in time to decimal degrees'''

    total_h = h + ( m / MIN_IN_A_H ) + ( s / SEC_IN_A_H )
    total_deg = total_h * DEC_IN_A_H

    return to_sf(sf=sf, res=total_deg)


def angular_to_linear(
        sf: int,
        theta: float,
        distance_away: float
    ) -> str:
    '''Convert angular distance or size to linear'''
    res = math.sin(theta) * distance_away
    return to_sf(sf=sf, res=res)


def magnitude_from_components(
        sf = int, 
        ra_comp = float,
        dec_comp = float
    ) -> str:
    """Compute magnitude using Pythagoras"""
    res = math.sqrt(ra_comp**2 + dec_comp**2)
    return to_sf(sf=sf, res=res)


def apparent_to_absolute_magnitude(
        apparent_magnitude: float,
        distance_away_in_pc: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None,
    ) -> str:
    """Compute absolute magnitude from apparent magnitude via distance"""
    res = apparent_magnitude - 5 * math.log10(distance_away_in_pc) + 5
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
        parallax_angle: float,
        sf: typing.Optional[int] = None,
        dp: typing.Optional[int] = None
) -> str:
    """Calculate distance in parsecs from parallax angle in arcseconds"""
    res = 1.0 / parallax_angle
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