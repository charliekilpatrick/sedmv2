import urllib
import os
import time
from astropy.time import Time,TimeDelta
from astropy import units as u
from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import numpy as np
import requests
import sys
import datetime
import util
from skyportal import skyportal

def main():

    if 'FRITZ_TOKEN' not in os.environ.keys():
        raise Exception('ERROR: add FRITZ_TOKEN to environmental variables')
    if 'SEDMV2_TARGETS' not in os.environ.keys():
        raise Exception('ERROR: add SEDMV2_TARGETS to environmental variables')

    api_token = os.environ['FRITZ_TOKEN']
    target_url = os.environ['SEDMV2_TARGETS']

    # Targets
    targets = util.download_targets(target_url)
    mask = (targets['type']!='SN Ia') & (targets['mag']<19.0)

    allocation_id = '1026'
    group_id = '1423'

    now=Time(datetime.datetime.utcnow())
    cadence = TimeDelta(4*86400*u.s)

    fritz = skyportal(os.environ['FRITZ_TOKEN'], 'sedmv2')
    print(fritz.targets)

    targets = targets[mask]
    for target in targets:
        obj_id = target['name']
        if obj_id.startswith('20'): obj_id = 'SN'+obj_id
        check = fritz.check_if_source_exists(obj_id, target['ra'], target['dec'])

        if not check:
            print(f'Adding {obj_id}')
            fritz.add_new_source(obj_id, target['ra'], target['dec'])

    data = fritz.get_observations(fritz.inst_id, now-cadence,now)

    curr_observations = []
    for d in data:
        if d['requester']['username']!='ckilpatrick': continue
        if 'expired' in d['status'].lower(): continue
        if 'deleted' in d['status'].lower(): continue

        curr_observations.append(d)

    for t in targets:
        obj_id = t['name']
        if obj_id.startswith('20'): obj_id = 'SN'+obj_id

        for filt in ['g','r','i','z']:

            needs_observation=True
            for c in curr_observations:
                if c['obj_id']==obj_id and c['payload']['observation_choice']==filt:
                    print(f'{obj_id} already has {filt} observation.  Skip...')
                    needs_observation = False
                    break

            if not needs_observation: continue

            exptime = util.get_exptime(t['mag'])

            status, obsdata = fritz.post_followup_request(allocation_id, obj_id,
                filt, exptime, now, now+TimeDelta(86400*u.s),
                priority=2, maximum_airmass=2, maximum_fwhm=1.5,
                minimum_lunar_distance=15.0, too='N', exposure_count=1,
                observation_type='transient')

            idnum = obsdata['id']
            if status:
                print(f'{obj_id} observation in {filt} successful.  obsid = {idnum}.')
            else:
                print(f'{obj_id} observation in {filt} failed.')


if __name__ == "__main__":
    main()
