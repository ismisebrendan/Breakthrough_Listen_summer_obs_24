from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, Angle
from astropy.time import Time
from astropy.table import Table

# Birr
birr_loc = EarthLocation(lat=53.095*u.deg, lon=-7.922*u.deg)
observing_time_start_birr = Time('2024-06-25 07:00:00', scale='utc', location=birr_loc)
observing_time_end_birr = Time('2024-06-26 14:00:00', scale='utc', location=birr_loc)
LST_start_birr = observing_time_start_birr.sidereal_time('mean')
LST_end_birr = observing_time_end_birr.sidereal_time('mean')

# Onsala
onsala_loc = EarthLocation(lat=57.399*u.deg, lon=11.930*u.deg)
observing_time_start_onsala = Time('2024-06-25 07:00:00', scale='utc', location=onsala_loc)
observing_time_end_onsala = Time('2024-06-26 14:00:00', scale='utc', location=onsala_loc)
LST_start_onsala = observing_time_start_onsala.sidereal_time('mean')
LST_end_onsala = observing_time_end_onsala.sidereal_time('mean')

# Midpoint
mid_lon = (birr_loc.lon + onsala_loc.lon)/2
mid_loc = EarthLocation(lat=57.399*u.deg, lon=mid_lon)
observing_time_start_mid = Time('2024-06-25 07:00:00', scale='utc', location=mid_loc)
observing_time_end_mid = Time('2024-06-26 14:00:00', scale='utc', location=mid_loc)
LST_start_mid = observing_time_start_mid.sidereal_time('mean')
LST_end_mid = observing_time_end_mid.sidereal_time('mean')

t = Table(names=('date', 'utc', 'lst_birr', 'lst_onsala', 'lst_mid'))

times = np.array([])
for i in range(7, 24):
    a = Time('2024-06-25 '+str(i)+':00:00', scale='utc', location=birr_loc).sidereal_time('mean').value
    b = Time('2024-06-25 '+str(i)+':00:00', scale='utc', location=onsala_loc).sidereal_time('mean').value
    c = Time('2024-06-25 '+str(i)+':00:00', scale='utc', location=mid_loc).sidereal_time('mean').value
    t.add_row((25, i, a, b, c))

for i in range(0, 15):
    a = Time('2024-06-26 '+str(i)+':00:00', scale='utc', location=birr_loc).sidereal_time('mean').value
    b = Time('2024-06-26 '+str(i)+':00:00', scale='utc', location=onsala_loc).sidereal_time('mean').value
    c = Time('2024-06-25 '+str(i)+':00:00', scale='utc', location=mid_loc).sidereal_time('mean').value
    t.add_row((26, i, a, b, c))

print('Starting LST in Birr:', LST_start_birr)
print('Starting LST in Onsala:', LST_start_onsala)
print('Starting LST at midpoint:', LST_start_mid)

print('\n')

print('Ending LST in Birr:', LST_end_birr)
print('Ending LST in Onsala:', LST_end_onsala)
print('Ending LST at midpoint:', LST_end_mid)

# import data and clean it up

data = ascii.read('PS_2024.06.17_03.31.30.csv')
data_unmodified = ascii.read('PS_2024.06.17_03.31.30.csv')
data_coords = ascii.read('PS_2024.06.17_03.31.30.csv')
data.remove_columns(['pl_name', 'hd_name', 'hip_name', 'default_flag', 'sy_mnum', 'cb_flag', 'sy_pnum', 'disc_year', 'disc_refname', 'disc_pubdate', 'disc_locale', 'disc_facility', 'disc_telescope', 'disc_instrument', 'rv_flag', 'pul_flag', 'ptv_flag', 'tran_flag', 'ast_flag', 'obm_flag', 'micro_flag', 'etv_flag', 'ima_flag', 'dkin_flag', 'pl_controv_flag', 'gaia_id', 'pl_refname', 'pl_orbperlim', 'pl_orbsmaxlim', 'soltype', 'pl_radelim', 'pl_radjlim', 'pl_masselim', 'pl_massjlim', 'rowupdate', 'pl_pubdate', 'releasedate', 'pl_nnotes', 'st_nphot', 'st_nrvc', 'st_nspec', 'pl_nespec', 'pl_ntranspec', 'pl_msinie', 'pl_msinieerr1', 'pl_msinieerr2', 'pl_msinielim', 'pl_msinij', 'pl_msinijerr1', 'pl_msinijerr2', 'pl_msinijlim', 'pl_cmasse', 'pl_cmasseerr1', 'pl_cmasseerr2', 'pl_cmasselim', 'pl_cmassj', 'pl_cmassjerr1', 'pl_cmassjerr2', 'pl_cmassjlim', 'pl_bmasse', 'pl_bmasseerr1', 'pl_bmasseerr2', 'pl_bmasselim', 'pl_bmassj', 'pl_bmassjerr1', 'pl_bmassjerr2', 'pl_bmassjlim', 'pl_bmassprov', 'pl_dens', 'pl_denserr1', 'pl_denserr2', 'pl_denslim', 'pl_orbeccenlim', 'pl_insollim', 'pl_eqtlim', 'pl_orbincllim', 'pl_tranmid', 'pl_tranmiderr1', 'pl_tranmiderr2', 'pl_tranmidlim', 'pl_tsystemref', 'ttv_flag', 'pl_imppar', 'pl_impparerr1', 'pl_impparerr2', 'pl_impparlim', 'pl_trandep', 'pl_trandeperr1', 'pl_trandeperr2', 'pl_trandeplim', 'pl_trandur', 'pl_trandurerr1', 'pl_trandurerr2', 'pl_trandurlim', 'pl_ratdor', 'pl_ratdorerr1', 'pl_ratdorerr2', 'pl_ratdorlim', 'pl_ratror', 'pl_ratrorerr1', 'pl_ratrorerr2', 'pl_ratrorlim', 'pl_occdep', 'pl_occdeperr1', 'pl_occdeperr2', 'pl_occdeplim', 'pl_orbtper', 'pl_orbtpererr1', 'pl_orbtpererr2', 'pl_orbtperlim', 'pl_orblper', 'pl_orblpererr1', 'pl_orblpererr2', 'pl_orblperlim', 'pl_rvamp', 'pl_rvamperr1', 'pl_rvamperr2', 'pl_rvamplim', 'pl_projobliq', 'pl_projobliqerr1', 'pl_projobliqerr2', 'pl_projobliqlim', 'pl_trueobliq', 'pl_trueobliqerr1', 'pl_trueobliqerr2', 'pl_trueobliqlim', 'st_refname', 'st_spectype', 'st_tefflim', 'st_radlim', 'st_masslim', 'st_met', 'st_meterr1', 'st_meterr2', 'st_metlim', 'st_metratio', 'st_logg', 'st_loggerr1', 'st_loggerr2', 'st_logglim', 'st_age', 'st_ageerr1', 'st_ageerr2', 'st_agelim', 'st_dens', 'st_denserr1', 'st_denserr2', 'st_denslim', 'st_vsin', 'st_vsinerr1', 'st_vsinerr2', 'st_vsinlim', 'st_rotp', 'st_rotperr1', 'st_rotperr2', 'st_rotplim', 'st_radv', 'st_radverr1', 'st_radverr2', 'st_radvlim', 'sy_refname', 'st_lumlim', 'sy_pm', 'sy_pmerr1', 'sy_pmerr2', 'sy_pmra', 'sy_pmraerr1', 'sy_pmraerr2', 'sy_pmdec', 'sy_pmdecerr1', 'sy_pmdecerr2', 'sy_dist', 'sy_disterr1', 'sy_disterr2', 'sy_plx', 'sy_plxerr1', 'sy_plxerr2', 'sy_bmag', 'sy_bmagerr1', 'sy_bmagerr2', 'sy_vmag', 'sy_vmagerr1', 'sy_vmagerr2', 'sy_jmag', 'sy_jmagerr1', 'sy_jmagerr2', 'sy_hmag', 'sy_hmagerr1', 'sy_hmagerr2', 'sy_kmag', 'sy_kmagerr1', 'sy_kmagerr2', 'sy_umag', 'sy_umagerr1', 'sy_umagerr2', 'sy_gmag', 'sy_gmagerr1', 'sy_gmagerr2', 'sy_rmag', 'sy_rmagerr1', 'sy_rmagerr2', 'sy_imag', 'sy_imagerr1', 'sy_imagerr2', 'sy_zmag', 'sy_zmagerr1', 'sy_zmagerr2', 'sy_w1mag', 'sy_w1magerr1', 'sy_w1magerr2', 'sy_w2mag', 'sy_w2magerr1', 'sy_w2magerr2', 'sy_w3mag', 'sy_w3magerr1', 'sy_w3magerr2', 'sy_w4mag', 'sy_w4magerr1', 'sy_w4magerr2', 'sy_gaiamag', 'sy_gaiamagerr1', 'sy_gaiamagerr2', 'sy_icmag', 'sy_icmagerr1', 'sy_icmagerr2', 'sy_tmag', 'sy_tmagerr1', 'sy_tmagerr2', 'sy_kepmag', 'sy_kepmagerr1', 'sy_kepmagerr2'])


data_coords.keep_columns(['pl_name', 'hostname', 'pl_letter', 'ra', 'rastr', 'dec', 'pl_eqt', 'sy_dist'])
data_coords.sort(keys='ra')

# start and stop observing 30 minutes before/after zenith
start_obs = data_coords['ra'] - 7.5
stop_obs = data_coords['ra'] + 7.5

data_coords['start_obs'] = start_obs
data_coords['stop_obs'] = stop_obs

coords_pl = coord.SkyCoord(ra=data_coords['ra']*u.deg, dec=data_coords['dec']*u.deg)

# pulsars
pulsars = ascii.read('pulsar.csv')
pulsars
coords_psr = coord.SkyCoord(ra=pulsars['ra'], dec=pulsars['dec'], unit=(u.hourangle, u.deg))

start_psr = Angle(pulsars['ra'], u.hourangle).degree - 7.5
stop_psr = Angle(pulsars['ra'], u.hourangle).degree + 7.5

# plot exoplanet positions
ra = coord.Angle(data_coords['ra']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(data_coords['dec']*u.degree)

ra_psr = coords_psr.ra
ra_psr = ra_psr.wrap_at(180*u.degree)
dec_psr = coords_psr.dec

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")
ax.scatter(ra.radian, dec.radian, label='exoplanets')
ax.scatter(ra_psr.radian, dec_psr.radian, label='pulsars', color='red')
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
ax.grid(True)
ax.legend()
plt.show()


plt.rcParams.update({'font.size': 10})

fig, ax = plt.subplots(figsize=(8,6))

for i in range(len(start_obs)):
    ax.broken_barh([(start_obs[i], 15)], (i-0.4, 0.8), label=data_coords['pl_name'])
    ax.text(start_obs[i]+0.5, i-0.2, data_coords['pl_name'][i])

for i in range(len(start_psr)):
    ax.broken_barh([(start_psr[i], 15)], (25+i-0.4, 0.8), label=pulsars['name'], color='red')
    ax.text(start_psr[i]+0.5, 25+i-0.2, pulsars['name'][i])

ax.vlines(t['lst_mid']*15, 0, 80)
 

#ax.set_xticks(np.arange(0, 361, 30))
#ax.set_xticklabels(np.arange(0, 25, 2))
ax.grid(True)

print(t['lst_mid']*15)

plt.show()