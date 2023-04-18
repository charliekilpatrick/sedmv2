from datetime import tzinfo, timedelta, datetime
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time
from astropy.table import unique, Column
import csv, requests, sys, numpy as np

target_table_names = ('name', 'ra', 'dec', 'priority', 'date', 'mag', 'type')
target_table_row   = [['X' * 40], [0.], [0.], [0.0],
    [Time('2019-01-01T00:00:00')], [0.0], ['X' * 40]]
max_length = 10000

# Get a blank target table
def blank_target_table(length=0):
    if length==0:
        tab = Table(target_table_row, names=target_table_names)
        return(tab[:0].copy())
    else:
        data = [el*length for el in target_table_row]
        tab = Table(data, names=target_table_names)
        return(tab)

def get_exptime(mag, filt, min_exp=30., s_n=10.0):

    zeropoints = {
        # Zero point estimates for detection of a source at 3-sigma in 1s
        # TODO: Update these once we have actual SEDMv2 data
        'g': 17.0,
        'r': 17.0,
        'i': 17.0,
        'z': 17.0,
    }

    if filt not in zeropoints.keys():
        raise Exception(f'ERROR: unrecognized filter {filt}')

    zeropoint = zeropoints[filt]

    term1 = (s_n / 3.)**2
    term2 = 0.4*(mag - zeropoint)
    exp_time = term1*10**term2

    if exp_time < min_exp:
        exp_time = min_exp

    exp_time = int(np.round(exp_time))

    return exp_time

def download_targets(url):

        # Use requests to get list of targets
        try:
            data = requests.get(url, timeout=30)
        except:
            error = 'ERROR: could not get a response from YSE PZ.  Exiting...'
            print(error)
            sys.exit()

        # Format into a table with the same names as standard target file
        table = ascii.read(data.text)

        for key in table.keys():
            if 'name' in key:
                table.rename_column(key, 'name')
                break

        # Only want to compare to BVgri magnitudes
        match = [(f.startswith('B') or f.startswith('V')
            or f.startswith('g') or f.startswith('r') or f.startswith('i'))
            for f in table['filter']]
        table = table[match]

        # We can get duplicte targets from the queries, weed these by iterating
        # through targets by name and checking the name
        newtable = table[:0].copy()
        for name in list(set(table['name'])):
            # Find the row with the most recent
            subtable = table[table['name']==name]
            idx = np.argmax([Time(t).mjd for t in subtable['obs_date']])

            newtable.add_row(subtable[idx])

        # Start by generating a blank table
        targets = blank_target_table()

        # Now
        now = Time(datetime.now())

        # Add targets one by one from
        for row in newtable:
            # Get rough approximation of target priority
            pri = 1.0
            delta_t_days = Time(datetime.now()).mjd - Time(row['obs_date']).mjd
            effective_mag = row['Recent mag'] + 0.03 * delta_t_days

            # Prioritize interesting classes of transients
            if row['spec_class'] in ['SN IIn','SN Ib','SN Ic','SN IIb','LBV',
                'SN Ibn','SN Ia-91T-like','SN Ia-91bg-like','SN Icn']:
                pri += 1.0

            if not row['spec_class']: continue

            # Prioritize bright things
            delta_t_days = now.mjd - Time(row['obs_date']).mjd
            effective_mag = row['Recent mag'] + 0.03 * delta_t_days
            pri += 17.0-effective_mag

            # Prioritize follow up requested
            if row['status_id'] == 2:
                pri += 40.0
            if row['status_id'] == 4:
                pri += 20.0

            try:
                if (Time(datetime.now()).mjd - Time(row['disc_date']).mjd) < 50:
                    pri += 2.0
            except AttributeError:
                continue

            pri = round(pri, 5)

            time    = Time(row['obs_date'])
            add_row = [row['name'], row['ra'],row['dec'], pri, now,
                effective_mag,row['spec_class']]
            targets.add_row(add_row)

        # Rearrange priority for descending order
        max_priority = np.max(targets['priority'])+1
        targets['priority'] = max_priority - targets['priority']

        return(targets)
