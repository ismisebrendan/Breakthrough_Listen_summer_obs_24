# Functions
from datetime import datetime
from astropy.time import Time

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
        The message for the prompt, should accept a timestamp in the form HH:MM.

    Returns
    -------
    date : str
        The inputted valid timestamp in the inputted format.
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

            # Check that the date is valid
            try:
                if (time[2] != ':') or (start_hour > 23) or (start_minute > 59):
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
        if columns != None:
            print(f'Please make sure the file has the columns {columns}')
        check_file
    
    return table

def want(objects):
    '''
    See if the user wants to observe a class of objects in the observation window.

    Parameters
    ----------
    objects : str
        The class of objects asked about.
    
    Returns
    -------
    bool
        Whether or not to observe the objects.
    '''
    obs = input(f'Do you want to observe {objects} in this observation? [Y/N]')

    if obs.upper() == 'Y':
        return True
    elif obs.upper() == 'N':
        return False
    else:
        want(objects)