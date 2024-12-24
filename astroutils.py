import math
import typing

def to_sf(sf: int, res: float) -> str:
    return f"%.{sf}f" % res

def to_radians(
        sf: int,
        degs: typing.Optional[float] = 0.0, 
        arcmins: typing.Optional[float] = 0.0, 
        arcsecs: typing.Optional[float] = 0.0,
        ) -> str:
    '''Convert degrees, arcminutes, and arcsedons to radians'''
    
    degs_to_radians = degs * (2*math.pi) / 360
    arcmins_to_radians = ( arcmins / 60 ) * (2*math.pi) / 360
    arcsecs_to_radians = ( arcsecs / 3600 ) * (2*math.pi) / 360

    res = degs_to_radians + arcmins_to_radians + arcsecs_to_radians

    return to_sf(sf=sf, res=res)


def from_radians(
        sf: int,
        radians: float
) -> str:
    '''Convert radians to degrees, arcminutes, and arcseconds'''
    degrees_f = (radians / (2*math.pi)) * 360
    print(degrees_f)
    degrees_f, degrees = math.modf(degrees_f)
    arcmins_f = degrees_f * 60
    arcmins_f, arcmins = math.modf(arcmins_f)
    arcsecs = arcmins_f * 60

    return (
        to_sf(sf=sf, res=degrees), 
        to_sf(sf=sf, res=arcmins), 
        to_sf(sf=sf, res=arcsecs)
    )
