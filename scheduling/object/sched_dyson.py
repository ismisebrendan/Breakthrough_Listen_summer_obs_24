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

sch = schedule(mid_lat, mid_lon, 0.5, 4/3)
print('---Creating schedule---')
sch.set_times('Tuesday', '20:00', 8, True)
sch.print_info()
print('Finding appropriate targets')
sch.find_targets_file('dyson_targets.csv')
sch.make_realta(output=False)
print(sch.sched_realta)
print(sch.targets)