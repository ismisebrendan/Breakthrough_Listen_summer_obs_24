# Import libraries

from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, Angle
from astropy.time import Time
from astropy.table import Table, unique
from datetime import datetime

# Functions

def check_date(message):
    '''
    Takes an inputted date and confirms that it is a valid future date.
    
    Parameters
    ----------
    message : str
        The message for the prompt, should accept a date in the form YYYY-MM-DD.

    Returns
    -------
    date : str
        The inputted valid date in YYYY-MM-DD format.
    '''
    # Get current date
    current_date = datetime.today()

    # Initially don't have a date
    have_date = False

    # Loop until given a valid date as input
    while have_date == False:
        # Get input and split into different parts
        date = input(message)
        try: 
            start_year = int(date[0:4])
            start_month = int(date[5:7])
            start_day = int(date[8:10])

            # Check if leap year
            if start_year % 4 == 0:
                if start_year % 400 != 0 and start_year % 100 == 0:
                    leap_year = False
                else:
                    leap_year = True
            else:
                leap_year = False

            # Check that the date is valid
            try:
                if (date[4] != '-') or (date[7] != '-') or (start_year < current_date.year) or (start_month > 12) or (start_month < current_date.month and start_year == current_date.year) or (start_day < current_date.day and start_month == current_date.month and start_year == current_date.year) or (start_month == 2 and leap_year == False and start_day > 28) or (start_month == 2 and leap_year == True and start_day > 29) or (start_day > 31) or ((start_month == 4 or start_month == 6 or start_month == 9 or start_month == 11) and start_day > 30) or (start_month <= 0) or (start_day <= 0):
                    print('Please provide a valid date in the future')
                else:
                    have_date = True
            except:
                print('Date not recognised please try again')
        
        except:
            print('Date not recognised please try again')

    return date

def check_time(message):
    '''
    Takes an inputted timestamp and confirms that it is a valid timestamp.
    
    Parameters
    ----------
    message : str
        The message for the prompt, should accept a timestamp in the form HH:MM:SS.

    Returns
    -------
    date : str
        The inputted valid timestamp in HH:MM:SS format.
    '''
    # Initially don't have a time
    have_time = False

    # Loop until given a valid time as input
    while have_time == False:
        # Get input and split into different parts
        time = input(message)
        try: 
            start_hour = int(time[0:2])
            start_minute = int(time[3:5])
            start_second = int(time[6:8])

            # Check that the date is valid
            try:
                if (time[2] != ':') or (time[5] != ':') or (start_hour > 23) or (start_minute > 59) or (start_second > 59):
                    print('Please provide a valid time')
                else:
                    have_time = True
            except:
                print('Time not recognised please try again')
        
        except:
            print('Time not recognised please try again')

    return time

def get_lst(loc, utc, name=None):
    '''
    Gets the LST at a point on the Earth's surface at a given time and prints it if given a name for the location, otherwise it returns the LST.

    Parameters
    ----------
    loc : EarthLocation
        The location in astropy EarthLocation format.
    utc
        The time and date in UTC that the LST is to be found at.
    name : str
        The name of the location that the LST is to be found at, default None.
    
    Returns
    -------
    lst
        The local siderial time at that location and time.
    '''
    time = Time(utc, scale='utc', location=loc)
    lst = time.sidereal_time('mean')
    if name != None: 
        print(f'The local siderial time at {name} at UTC: {utc} is: ' + str(lst))
    else:
        return lst
    
def check_file(message, columns=None):
    '''
    Checks an inputted path to a CSV and attempts to keep the specified columns, if none are specified, all columns are kept.

    Parameters
    ----------
    message : str
        The propmt for the CSV file to be dealt with.
    columns : array likee
        The list of columns to be kept in the CSV, None.
    Returns
    -------
    table
        An astropy table from the CSV.
    '''
    # Initially don't have a valid file
    have_file = False
    
    file = input(message)

    # If the file is valid we can quit, otherwise it tries again
    try:
        table = ascii.read(file)
        if columns != None:
            table.keep_columns(columns)
        have_file = True
    except:
        print('File name not recognised, please try again')
        check_file
    
    return table
    '''
    Checks an inputted path to a CSV and attempts to keep the specified columns, if none are specified, all columns are kept.

    Parameters
    ----------
    message : str
        The propmt for the CSV file to be dealt with.
    columns : array likee
        The list of columns to be kept in the CSV, None.
    Returns
    -------
    table
        An astropy table from the CSV.
    '''
    # Initially don't have a valid file
    have_file = False
    
    file = input(message)

    # If the file is valid we can quit, otherwise it tries again
    try:
        table = ascii.read(file)
        if columns != None:
            table.keep_columns(columns)
        have_file = True
    except:
        print('File name not recognised, please try again')
        check_file
    
    return table

##################################
#                                #
#   Set up times and locations   #
#                                #
##################################

valid_dates = False
valid_times = False

# Check that dates and times are in the right order
while valid_dates == False:
    start_date = check_date('Please enter starting date of observation (YYYY-MM-DD): ')
    end_date = check_date('Please enter ending date of observation (YYYY-MM-DD): ')

    if (int(start_date[0:4]) > int(end_date[0:4])) or (int(start_date[0:4]) == int(end_date[0:4]) and int(start_date[5:7]) > int(end_date[5:7])) or (start_date[0:4] == end_date[0:4] and start_date[5:7] == end_date[5:7] and int(start_date[8:10]) > int(end_date[8:10])):
        print('End date must be the same or later than the start date')
    else:
        valid_dates = True

while valid_times == False:
    start_time = check_time('Please enter starting time of observation (HH:MM:SS): ')
    end_time = check_time('Please enter ending time of observation (HH:MM:SS): ')

    if start_date == end_date:
        if (start_time == end_time) or (int(start_time[0:2]) > int(end_time[0:2])) or (int(start_time[0:2]) == int(end_time[0:2]) and int(start_time[3:5]) > int(end_time[3:5])) or (int(start_time[0:2]) == int(end_time[0:2]) and int(start_time[3:5]) == int(end_time[3:5]) and int(start_time[6:8]) >= int(end_time[6:8])):
            print('End time must be later than the start time')
        else:
            valid_times = True
    else:
        valid_times = True

# Start and end in right format
starting = start_date + ' ' + start_time
ending = end_date + ' ' + end_time

# Birr
birr_loc = EarthLocation(lat=53.095*u.deg, lon=-7.922*u.deg)
get_lst(birr_loc, starting, 'Birr')

# Onsala
onsala_loc = EarthLocation(lat=57.399*u.deg, lon=11.930*u.deg)

# Midpoint
mid_lon = (birr_loc.lon + onsala_loc.lon)/2
mid_loc = EarthLocation(lat=57.399*u.deg, lon=mid_lon)
LST_start_mid = get_lst(mid_loc, starting)
LST_end_mid = get_lst(mid_loc, ending)

print('Starting LST at midpoint:', LST_start_mid)
get_lst(birr_loc, starting, 'Birr')
get_lst(onsala_loc, starting, 'Onsala')

print('\n')

print('Ending LST at midpoint:', LST_end_mid)
get_lst(birr_loc, ending, 'Birr')
get_lst(onsala_loc, ending, 'Onsala')

