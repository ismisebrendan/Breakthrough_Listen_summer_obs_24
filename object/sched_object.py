import datetime as dt
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, SkyCoord

class schedule():
    """
    The schedule for the dual site observations.

    Parameters
    ----------
    latitude : float
        The latitude of the observatory.
    longitude : float
        The longitude of the observatory.
    """

    def __init__(self, latitude, longitude):
        self.location = EarthLocation(lat=latitude, lon=longitude)
        self.start_obs = None
        self.end_obs = None
        self.start_LST = None
        self.end_LST = None
        self.pointing_time_psr_cal = None
        self.pointing_time_planet = None
    
    def print_info(self):
        print(f'Starting at {self.start_obs}')
        print(f'Ending at {self.end_obs}')
    
    def get_times(self, day, time, duration, buffer=False):
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
        week_days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
        day_int = week_days.index(day)

        # Find the date of the next occurrence of the day
        next_days = (day_int - dt.datetime.today().weekday()) % 7

        # Set the start and end time and date
        start_time = dt.datetime.today() + dt.timedelta(days=next_days)
        start_time.replace(hour=int(start_hr), minute=int(start_min), second=0)
        end_time = start_time + dt.timedelta(hours=duration)

        self.start_obs = Time(start_time.strftime('%Y-%m-%d %H:%M'), scale='utc', location=self.location)
        self.end_obs = Time(end_time.strftime('%Y-%m-%d %H:%M'), scale='utc', location=self.location)
        self.start_LST = self.start_obs.sidereal_time('mean')
        self.end_LST = self.end_obs.sidereal_time('mean')