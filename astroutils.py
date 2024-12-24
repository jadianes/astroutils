import math
import typing

def to_radians(
        sf: int,
        degs: typing.Optional[float] = 0.0, 
        arcmins: typing.Optional[float] = 0.0, 
        arcsecs: typing.Optional[float] = 0.0,
        ) -> float:
    '''Convert degrees, arcminutes, and arcsedons to radians'''
    
    degs_to_radians = degs * (2*math.pi) / 360
    arcmins_to_radians = ( arcmins / 60 ) * (2*math.pi) / 360
    arcsecs_to_radians = ( arcsecs / 3600 ) * (2*math.pi) / 360

    res = degs_to_radians + arcmins_to_radians + arcsecs_to_radians

    return f"%.{sf}f" % res
    
