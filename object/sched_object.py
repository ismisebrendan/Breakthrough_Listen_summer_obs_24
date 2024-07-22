import numpy as np
import datetime as dt
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, SkyCoord
import astropy.units as u
from psrqpy import QueryATNF
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.table import Table
from astropy.io import ascii


class schedule():
    """
    The schedule for the dual site observations.
    
    Parameters
    ----------
    latitude : float
        The latitude of the observatory.
    longitude : float
        The longitude of the observatory.
    calibraton : float
        The number of hours to observe the calibration pulsar for.
    target : float
        The number of hours to observe a target for.
    ignore : str, optional
        The path to a .csv that has a list of targets to not observe.
    """

    def __init__(self, latitude, longitude, calibration, target, ignore=None):
        self.location = EarthLocation(lat=latitude, lon=longitude)
        self.point_cal = calibration + 1/60
        self.point_target = target + 1/60
        self.start_obs = None
        self.end_obs = None
        self.start_LST = None
        self.end_LST = None
        self.obs_time = None
        self.planets = Table(names=['hostname', 'ra', 'dec', 'sy_dist'], dtype=[str, float, float, float])
        self.targets = Table(names=['hostname', 'ra', 'dec'], dtype=[str, float, float])
        self.pulsars = None
        self.psr_coords = None
        self.sched_iLiSA = Table(names=('Name', 'Time', 'RA', 'DEC', 'freqrng', 'dur'), dtype=(str, str, float, float, str, str))
        self.sched_realta = Table(names=('start', '-', 'stop', ':', 'name', 'coords'), dtype=(str, str, str, str, str, str))
        try:
            self.ignore = ascii.read(ignore)
            self.ignore.remove_column('date')
        except:
            self.ignore = Table()
    
    def print_info(self):
        """
        Print the information about the schedule.
        """
        print(f'Observing from {self.location.lat}, {self.location.lon}')
        print(f'\nStarting observation at {self.start_obs}')
        print(f'Ending observation at {self.end_obs}')
        print(f'\nObserving calibration pulsar for {self.point_cal - 1/60} hours')
        print(f'Observing targets for {self.point_target - 1/60} hours')
    
    def set_times(self, day, time, duration, buffer=False):
        """
        Get the starting and ending dates and times in UTC and LST.

        Parameters
        ----------
        day : string which is a day
            The day on which the observation is to start.
        time : string in the form HH:MM
            The starting time of the observation.
        duration : float
            The duration of the observations.
        buffer : bool, default=False
            Whether or not to add in a two minute buffer at the start and a 1 minute buffer at the end.
        
        Returns
        -------
        start_obs
            The starting time of the observation in UTC format
        end_obs
            The ending time of the observation in UTC format
        start_LST
            The starting time of the observation in LST format
        end_LST
            The ending time of the observation in LST format
        """
        # Split the time into hours and minutes
        start_hr, start_min = time.split(':')

        # Convert the days from strings to ints
        week_days = ['monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday']
        day_int = week_days.index(day.lower())

        # Find the date of the next occurrence of the day
        next_days = (day_int - dt.datetime.today().weekday()) % 7

        # Set the start and end time and date
        start_time = dt.datetime.today() + dt.timedelta(days=next_days)
        start_time = start_time.replace(hour=int(start_hr), minute=int(start_min), second=0, microsecond=0)
        end_time = start_time + dt.timedelta(hours=duration)

        # Add in the buffer if needed
        if buffer == True:
            start_time = start_time + dt.timedelta(minutes=2)
            end_time = end_time - dt.timedelta(minutes=1)

        # Convert to the correct format
        self.start_obs = Time(start_time.strftime(f'%Y-%m-%d %H:%M'), scale='utc', location=self.location)
        self.end_obs = Time(end_time.strftime('%Y-%m-%d %H:%M'), scale='utc', location=self.location)
        self.start_LST = self.start_obs.sidereal_time('mean')
        self.end_LST = self.end_obs.sidereal_time('mean')

        self.obs_time = duration

    def find_targets_query(self):
        """
        Find the list of targets that are to be observed by querying the ATNF and the NASA exoplanet catalogues.
        """
        # Keep track of the current time and the elapsed time
        time_LST = self.start_LST.value
        time_offset = 0

        # List of targets
        target_list = []

        # Observe the calibration pulsar from the ATNF database that is closest to zenith halfway through the observation window
        mid_ra = Angle((time_LST + (self.point_cal - 1/60)/2) * 15, u.deg)
        zenith = SkyCoord(ra=mid_ra.to_string(u.hourangle), dec=self.location.lat, unit='deg')
        c = [mid_ra.to_string(u.hourangle), self.location.lat.to_string(), 20.]
        pulsars = QueryATNF(params=['NAME', 'RAJ', 'DECJ', 'R_Lum'], circular_boundary=c, condition='R_Lum > 0').table
        psr_coords = SkyCoord(ra=np.array(pulsars['RAJ']), dec=np.array(pulsars['DECJ']), unit=(u.hourangle, u.deg))

        sep = zenith.separation(psr_coords)
        ind = np.argmin(sep)
        target_list.append(pulsars['NAME'][ind])

        # Remove all other pulsars from pulsars and psr_coords
        self.pulsars = pulsars[ind]
        self.psr_coords = psr_coords[ind]

        # Move the time on
        time_LST += self.point_cal
        time_offset += self.point_cal

        # Query the NASA exoplanet archive for a list of potential targets with dec +- 20 degrees from zenith
        potential_targets = NasaExoplanetArchive.query_criteria(table="pscomppars", where=f"disc_facility like '%TESS%' and dec between {self.location.lat.value - 20} and {self.location.lat.value + 20}")
        potential_targets.keep_columns(('hostname', 'disc_facility', 'ra', 'dec', 'sy_dist'))
        potential_targets_coords = SkyCoord(ra=potential_targets['ra'], dec=potential_targets['dec'], unit='deg')

        # The rest of the observation window
        while time_offset <= self.obs_time:
            # Find zenith in the middle of observing
            mid_lst = (time_LST + (self.point_target - 1/60)/2) * 15
            zenith = SkyCoord(ra=mid_lst, dec=self.location.lat.value, unit='deg')

            # Delete the specified targets to be ignored
            for i in range(len(self.ignore)):
                ind = potential_targets['hostname'] != self.ignore['name'][i]
                potential_targets = potential_targets[ind]
            
            # Make sure same target is not observed multiple times in one session
            for i in range(len(target_list)):
                ind = potential_targets['hostname'] != target_list[i]
                potential_targets = potential_targets[ind]
                potential_targets_coords = potential_targets_coords[ind]
            
            sep = zenith.separation(potential_targets_coords)
            ind = np.argmin(sep)
            sep = sep[ind]

            # Record the planet
            target_list.append(potential_targets['hostname'][ind])
            self.planets.add_row((potential_targets['hostname'][ind], potential_targets['ra'][ind], potential_targets['dec'][ind], potential_targets['sy_dist'][ind]))
            
            # Move on
            time_offset += self.point_target
            time_LST += self.point_target
    
    def find_targets_file(self, file):
        """
        Find the list of targets that are to be observed from an input file.

        Parameters
        ----------
        file : str
            The path to a file containing the potential targets.
        """
        # Keep track of the current time and the elapsed time
        time_LST = self.start_LST.value
        time_offset = 0

        # List of targets
        target_list = []

        # Observe the calibration pulsar from the ATNF database that is closest to zenith halfway through the observation window
        mid_ra = Angle((time_LST + (self.point_cal - 1/60)/2) * 15, u.deg)
        zenith = SkyCoord(ra=mid_ra.to_string(u.hourangle), dec=self.location.lat, unit='deg')
        c = [mid_ra.to_string(u.hourangle), self.location.lat.to_string(), 20.]
        pulsars = QueryATNF(params=['NAME', 'RAJ', 'DECJ', 'R_Lum'], circular_boundary=c, condition='R_Lum > 0').table
        psr_coords = SkyCoord(ra=np.array(pulsars['RAJ']), dec=np.array(pulsars['DECJ']), unit=(u.hourangle, u.deg))

        sep = zenith.separation(psr_coords)
        ind = np.argmin(sep)
        target_list.append(pulsars['NAME'][ind])

        # Remove all other pulsars from pulsars and psr_coords
        self.pulsars = pulsars[ind]
        self.psr_coords = psr_coords[ind]

        # Move the time on
        time_LST += self.point_cal
        time_offset += self.point_cal

        # The potential targets
        potential_targets = ascii.read(file)
        potential_targets_coords = SkyCoord(ra=potential_targets['RA'], dec=potential_targets['DEC'], unit='deg')

        # The rest of the observation window
        while time_offset <= self.obs_time:
            # Find zenith in the middle of observing
            mid_lst = (time_LST + (self.point_target - 1/60)/2) * 15
            zenith = SkyCoord(ra=mid_lst, dec=self.location.lat.value, unit='deg')
           
            # Make sure same target is not observed multiple times in one session
            for i in range(len(target_list)):
                ind = potential_targets['Gaia ID'] != target_list[i]
                potential_targets = potential_targets[ind]
                potential_targets_coords = potential_targets_coords[ind]

            sep = zenith.separation(potential_targets_coords)
            ind = np.argmin(sep)
            sep = sep[ind]

            # Record the planet
            target_list.append(potential_targets['Gaia ID'][ind])
            self.targets.add_row((potential_targets['Gaia ID'][ind], potential_targets['RA'][ind], potential_targets['DEC'][ind]))
            
            # Move on
            time_offset += self.point_target
            time_LST += self.point_target

    def make_iLiSA(self):
        """
        Produce the schedule in iLiSA format.
        """
        freq_range = '110e6:190e6'

        # Starting time
        time = self.start_obs.mjd

        # Record the calibraton pulsar
        self.sched_iLiSA.add_row((self.pulsars['NAME'], Time(time, format='mjd').iso[11:16], self.psr_coords.ra.value*u.deg.to(u.rad), self.psr_coords.dec.value*u.deg.to(u.rad), freq_range, str(int(np.round((self.point_cal - 1/60) * 60))) + 'm'))
        time += self.point_cal/24

        # Add the planets
        for i in range(len(self.planets)):
            if i == len(self.planets) - 1:
                dur = np.round((self.end_obs.mjd - time)*u.day.to(u.min))
                self.sched_iLiSA.add_row((self.planets['hostname'][-1], Time(time, format='mjd').iso[11:16], self.planets['ra'][-1]*u.deg.to(u.rad), self.planets['dec'][-1]*u.deg.to(u.rad), freq_range, str(dur) + 'm'))
            else:
                self.sched_iLiSA.add_row((self.planets['hostname'][i], Time(time, format='mjd').iso[11:16], self.planets['ra'][i]*u.deg.to(u.rad), self.planets['dec'][i]*u.deg.to(u.rad), freq_range, str(int(np.round((self.point_target - 1/60) * 60))) + 'm'))

            time += self.point_target/24
        
        # Output the file
        ascii.write(self.sched_iLiSA, 'sched_iLiSA.txt', overwrite=True)
    
    def make_realta(self):
        """
        Produce the schedule in realta format.
        """
        time = self.start_obs.mjd

        end_time = time + (self.point_cal - 1/60)/24

        self.sched_realta.add_row((Time(time, format='mjd').iso[11:16], '-', Time(end_time, format='mjd').iso[11:16], ':', self.pulsars['NAME'], f'[{self.psr_coords.ra.value*u.deg.to(u.rad)}, {self.psr_coords.dec.value*u.deg.to(u.rad)}, \'J2000\']'))

        time += self.point_cal/24

        for i in range(0, len(self.planets)):
            end_time = time + (self.point_target - 1/60)/24
            # Make sure the end time doesn't overshoot
            if end_time > self.end_obs.mjd:
                end_time = self.end_obs.mjd

            self.sched_realta.add_row((Time(time, format='mjd').iso[11:16], '-', Time(end_time, format='mjd').iso[11:16], ':', self.planets['hostname'][i], f'[{self.planets[i][1]*u.deg.to(u.rad)}, {self.planets[i][2]*u.deg.to(u.rad)}, \'J2000\']'))

            # Wait the 1 hr 21 minutes
            time += self.point_target/24

        # Output the file
        ascii.write(self.sched_realta, 'sched_realta.txt', overwrite=True)