from sched_object import schedule
from astropy.coordinates import EarthLocation
import astropy.units as u
import sys

# Birr
birr_loc = EarthLocation(lat=53.095*u.deg, lon=-7.922*u.deg)

# Onsala
onsala_loc = EarthLocation(lat=57.399*u.deg, lon=11.930*u.deg)

# Midpoint
mid_lon = (birr_loc.lon + onsala_loc.lon) / 2
mid_lat = (birr_loc.lat + onsala_loc.lat) / 2
mid_loc = EarthLocation(lat=mid_lat, lon=mid_lon)

sch = schedule(55, 2, 0.5, 4/3)
sch.set_times('Tuesday', '20:00', 8, buffer=True)
sch.find_targets_file('object\dyson_targets.csv')
print(sch.targets)