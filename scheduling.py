from astropy.io import ascii
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, Angle, SkyCoord
from astropy.time import Time
from astropy.table import Table
from funcs import *
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import datetime as dt
from psrqpy import QueryATNF    

##############
#            #
#   Set up   #
#            #
##############

###########
# Timings #
###########

# Find the dates of the next Tuesday (including possibly today) and Wednesday
next_Tuesday_days = (1 - dt.datetime.today().weekday()) % 7
next_Tuesday = dt.datetime.today() + dt.timedelta(days=next_Tuesday_days)
next_Wednesday = dt.datetime.today() + dt.timedelta(days=next_Tuesday_days + 1)

starting = next_Tuesday.strftime('%Y-%m-%d') + ' 20:02'
ending = next_Wednesday.strftime('%Y-%m-%d') + ' 03:59'
starting_time = Time(starting, scale='utc')
ending_time = Time(ending, scale='utc')

print('--- DEFAULT PARAMETERS ---')
print('Observe from 20:02 UTC (allowing two minutes grace) ' + next_Tuesday.strftime('%A %d/%m') +' to 03:59 UTC ' + next_Wednesday.strftime('%A %d/%m'))
print('First observe a test pulsar for 30 minutes')
print('Observe exoplanets for 80 minutes each for the rest of the time')

#############
# Locations #
#############

# Birr
birr_loc = EarthLocation(lat=53.095*u.deg, lon=-7.922*u.deg)

# Onsala
onsala_loc = EarthLocation(lat=57.399*u.deg, lon=11.930*u.deg)

# Midpoint
mid_lon = (birr_loc.lon + onsala_loc.lon) / 2
mid_lat = (birr_loc.lat + onsala_loc.lat) / 2
mid_loc = EarthLocation(lat=mid_lat, lon=mid_lon)
LST_start_mid = get_lst(mid_loc, starting)
LST_end_mid = get_lst(mid_loc, ending)

print(f'Starting at {starting} UTC')
print('Starting LST at midpoint:', LST_start_mid)
print(f'Ending at {ending} UTC')
print('Ending LST at midpoint:', LST_end_mid)

##################
# Pointing times #
##################

pointing_time_psr_cal = 0.5 + 1/60
pointing_time_planet = 4/3 + 1/60
planets = Table(names=['hostname', 'ra', 'dec', 'sy_dist'], dtype=[str, float, float, float])

# List of targets to ignore (uusally because they had been observed previously)
try:
    ignore = ascii.read('observed.csv')
    ignore.remove_column('date')
except:
    ignore = Table()

####################
#                  #
#   Observations   #
#                  #
####################

time_LST = LST_start_mid.value
time_offset = 0
target_list = []
target_type = []

# How long from the start we will observe for
obs_time = (ending_time.mjd - starting_time.mjd) * 24 # hours

################################
# Calibration pulsar initially #
################################

# Query a pulsar from the ATNF
mid_ra = Angle((time_LST + (pointing_time_psr_cal - 1/60)/2) * 15, u.deg)
zenith = SkyCoord(ra=mid_ra.to_string(u.hourangle), dec=mid_lat.to_string(), unit='deg')

c = [mid_ra.to_string(u.hourangle), mid_lat.to_string(), 20.]

pulsars = QueryATNF(params=['NAME', 'RAJ', 'DECJ', 'R_Lum'], circular_boundary=c, condition='R_Lum > 0').table
psr_coords = SkyCoord(ra=np.array(pulsars['RAJ']), dec=np.array(pulsars['DECJ']), unit=(u.hourangle, u.deg))

sep = zenith.separation(psr_coords)
ind = np.argmin(sep)
target_list.append(pulsars['NAME'][ind])
target_type.append('pulsar')

# Remove all other pulsars from pulsars and psr_coords
pulsars = pulsars[ind]
psr_coords = psr_coords[ind]

# Move the time on
time_offset += pointing_time_psr_cal
time_LST += pointing_time_psr_cal

#########################
# Go for obs_time hours #
#########################

while time_offset <= obs_time:
    # Find the RA that is directly overhead in the middle of the observation window
    # Choose to observe the target that is closest to the zenith at the middle of the observation window

    # Find zenith in the middle of observing
    mid_lst = (time_LST + (pointing_time_planet - 1/60)/2) * 15
    zenith = SkyCoord(ra=mid_lst, dec=mid_lat, unit='deg')

    # Query for list of planets to target near the zenith
    potential_targets = NasaExoplanetArchive.query_region(table='pscomppars', coordinates=zenith, radius=20*u.deg)
    potential_targets.keep_columns(('hostname', 'disc_facility', 'ra', 'dec', 'sy_dist'))

    # Only keep TESS targets
    ind = potential_targets['disc_facility'] == 'Transiting Exoplanet Survey Satellite (TESS)'
    potential_targets = potential_targets[ind]

    # Delete the specified targets to be ignored
    for i in range(len(ignore)):
        ind = potential_targets['hostname'] != ignore['name'][i]
        potential_targets = potential_targets[ind]
    
    # Make sure same target is not observed multiple times in one session
    for i in range(len(target_list)):
        ind = potential_targets['hostname'] != target_list[i]
        potential_targets = potential_targets[ind]

    # The coordinates of the planets
    planets_coords = SkyCoord(ra=potential_targets['ra'], dec=potential_targets['dec'], unit='deg')

    sep = zenith.separation(planets_coords)
    ind = np.argmin(sep)
    sep = sep[ind]

    # Record this planet
    target_list.append(potential_targets['hostname'][ind])
    target_type.append('planet')
    planets.add_row((potential_targets['hostname'][ind], potential_targets['ra'][ind], potential_targets['dec'][ind], potential_targets['sy_dist'][ind]))
    
    # Move on
    time_offset += pointing_time_planet
    time_LST += pointing_time_planet

#############
#           #
#   Plots   #
#           #
#############

###########################
# Optimum viewing windows #
###########################

fig, ax1 = plt.subplots(figsize=(8,6))
ax2 = ax1.twinx()

# Targeted planets
for i in range(len(planets)):
    # Plot from optimum start time of observation
    ax1.broken_barh([(planets['ra'][i] - (pointing_time_planet - 1/60) * 7.5, (pointing_time_planet - 1/60) * 15)], (planets['sy_dist'][i], 10), color='red')

# Plot the calibration pulsar
ax2.broken_barh([(Angle(psr_coords.ra.value, u.deg).value - (pointing_time_psr_cal - 1/60) * 7.5, (pointing_time_psr_cal - 1/60) * 15)], (pulsars['R_LUM'], 10), color='green')

ax1.set_title('Optimum observation windows')
ax1.set_xticks(np.arange(0, 361, 30))
ax1.set_xticklabels(np.arange(0, 25, 2))
ax1.grid(True)
ax1.set_xlabel('RA [hours]')
ax1.set_ylabel('Distance of exoplanets [pc]', color='red')
ax1.tick_params(axis='y', labelcolor='red')
ax2.set_ylabel(r'Luminosity of pulsars @ 400 MHz [mJy kpc$^2$]', color='green')
ax2.tick_params(axis='y', labelcolor='green')
ax1.set_ylim([0, max(planets['sy_dist']) + 50])
ax2.set_ylim([0, pulsars['R_LUM'] + 50])
plt.savefig('optimum.png')

########################
# Actual viewing times #
########################

fig, ax1 = plt.subplots(figsize=(8,6))

current_LST = LST_start_mid.value * 360/24
end_LST = LST_end_mid.value * 360/24

# Calibration pulsar
ax1.broken_barh([(current_LST, pointing_time_psr_cal * 15)], (0, 10), color='green')
ax1.text(current_LST, 5, pulsars['NAME'])
ax1.plot(psr_coords.ra.value, 0, 'bo')

current_LST += pointing_time_psr_cal * 15

# Plot planets
for i in range(len(planets)-1):
    ax1.broken_barh([(current_LST, (pointing_time_planet * 15))], ((i+1)*10, 10), color='red')
    ax1.text((current_LST) % 360, (i+1)*10+5, planets['hostname'][i])
    ax1.plot(planets['ra'][i], (i+1)*10, 'bo')

    current_LST += pointing_time_planet * 15
    current_LST %= 360

# Plot last planet (with shorter observing time) in different colour
ax1.broken_barh([(current_LST, end_LST - current_LST)], ((i+2)*10, 10), color='blue')
ax1.text((current_LST) % 360, (i+2)*10+5, planets['hostname'][-1])
ax1.plot(planets['ra'][-1], (i+2)*10, 'bo')

ax1.set_xticks(np.arange(0, 361, 30))
ax1.set_xticklabels(np.arange(0, 25, 2))
ax1.grid(True)
ax1.set_xlabel('LST [hours]')
ax1.set_title('Scheduled time, LST against Distance of stars, blue dots - when the star reaches its zenith')
plt.savefig('viewing.png')

###############
#             #
#   Outputs   #
#             #
###############

################
# iLiSA format #
################

# Take same frequency range for all observations
freq_range = '100e6:190e6'

# Table for output
sched_iLiSA = Table(names=('Name', 'Time', 'RA', 'DEC', 'freqrng', 'dur'), dtype=(str, str, float, float, str, str))

# Starting time
time = starting_time.mjd

# Record the calibraton pulsar
sched_iLiSA.add_row((pulsars['NAME'], Time(time, format='mjd').iso[11:16], psr_coords.ra.value*u.deg.to(u.rad), psr_coords.dec.value*u.deg.to(u.rad), freq_range, str(int(np.round((pointing_time_psr_cal - 1/60) * 60))) + 'm'))

# Wait for the 31 minutes
time += pointing_time_psr_cal/24

# Add the planets
for i in range(len(target_list) - 1):
    if i == len(target_list) - 2:
        dur = np.round((ending_time.mjd - time)*u.day.to(u.min))
        sched_iLiSA.add_row((planets['hostname'][-1], Time(time, format='mjd').iso[11:16], planets['ra'][-1]*u.deg.to(u.rad), planets['dec'][-1]*u.deg.to(u.rad), freq_range, str(dur) + 'm'))
    else:
        sched_iLiSA.add_row((planets['hostname'][i], Time(time, format='mjd').iso[11:16], planets['ra'][i]*u.deg.to(u.rad), planets['dec'][i]*u.deg.to(u.rad), freq_range, str(int(np.round((pointing_time_planet - 1/60) * 60))) + 'm'))

    # Wait the 1 hr 21 minutes
    time += pointing_time_planet/24

##################
# I-LOFAR format #
##################

sched_realta = Table(names=('start', '-', 'stop', ':', 'name', 'coords'), dtype=(str, str, str, str, str, str))

time = starting_time.mjd

end_time = time + (pointing_time_psr_cal - 1/60)/24

sched_realta.add_row((Time(time, format='mjd').iso, '-', Time(end_time, format='mjd').iso, ':', pulsars['NAME'], f'[{psr_coords.ra.value*u.deg.to(u.rad)}, {psr_coords.dec.value*u.deg.to(u.rad)}, \'J2000\']'))

time += pointing_time_psr_cal/24

for i in range(0, len(target_list) - 1):
    end_time = time + (pointing_time_planet - 1/60)/24
    # Make sure the end time doesn't overshoot
    if end_time > ending_time.mjd:
        end_time = ending_time.mjd

    sched_realta.add_row((Time(time, format='mjd').iso, '-', Time(end_time, format='mjd').iso, ':', planets['hostname'][i], f'[{planets[i][1]*u.deg.to(u.rad)}, {planets[i][2]*u.deg.to(u.rad)}, \'J2000\']'))

    # Wait the 1 hr 21 minutes
    time += pointing_time_planet/24     

################
# Output files #
################

ascii.write(sched_iLiSA, 'sched_iLiSA_temp.txt', overwrite=True)
ascii.write(sched_realta, 'sched_realta_temp.txt', overwrite=True)

# Need to remove the quotation marks
iLiSA = open('sched_iLiSA_temp.txt', 'r')
iLiSA_out = open('sched_iLiSA_out.txt', 'w')

rmv_quotes = iLiSA.read()
rmv_quotes = rmv_quotes.replace('"', '')

iLiSA_out.write(rmv_quotes)

realta = open('sched_realta_temp.txt', 'r')
realta_out = open('sched_realta_out.txt', 'w')

rmv_quotes = realta.read()
rmv_quotes = rmv_quotes.replace('"', '')

realta_out.write(rmv_quotes)

# Close the files
iLiSA.close()
iLiSA_out.close()

realta.close()
realta_out.close()