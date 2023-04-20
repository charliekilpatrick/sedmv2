import os
from astropy.time import Time,TimeDelta
from astropy import units as u
import sys
import datetime

# Package dependencies
import util
from skyportal import skyportal
import option

def main():

    args = option.parse_arguments()

    api_token = option.parse_token(args)
    targets = util.parse_targets(args)

    # Handling arguments for each target
    #targets, kwargs = util.handle_arguments()

    # Targets - only observe bright sources and ignore SN Ia
    if 'type' in targets.keys() and args.mask_types:
        for typ in args.mask_types.split(','):
            mask = targets['type']!=typ
            targets = targets[mask]

    if args.max_mag:
        mask = targets['mag']<args.max_mag
        targets = targets[mask]

    if args.min_mag:
        mask = targets['mag']>args.min_mag
        targets = targets[mask]

    # Allocation and group for SEDMv2
    allocation_id = args.allocation_id
    group_id = args.group_id

    now=Time(datetime.datetime.utcnow())
    cadence = TimeDelta(args.cadence*86400*u.s)

    fritz = skyportal(api_token, 'sedmv2')

    for target in targets:
        obj_id = target['name']
        if obj_id.startswith('20'): obj_id = 'SN'+obj_id
        check = fritz.check_if_source_exists(obj_id,target['ra'],target['dec'])

        if not check:
            print(f'Adding {obj_id}')
            fritz.add_new_source(obj_id, target['ra'], target['dec'])

    data = fritz.get_observations(fritz.inst_id, now-cadence, now)

    curr_observations = []
    for d in data:
        if d['requester']['username']!=args.user: continue
        if 'expired' in d['status'].lower(): continue
        if 'deleted' in d['status'].lower(): continue

        curr_observations.append(d)

    if args.clean_requests:
        for d in curr_observations:
            fritz.delete_followup_request(d)
        curr_observations = []


    for t in targets:
        obj_id = t['name']
        if obj_id.startswith('20'): obj_id = 'SN'+obj_id

        for filt in args.filters.split(','):

            needs_observation=True
            for c in curr_observations:
                if c['obj_id']==obj_id and c['payload']['observation_choice']==filt:
                    print(f'{obj_id} already has {filt} observation.  Skip...')
                    needs_observation = False
                    break

            if not needs_observation: continue

            exptime = util.get_exptime(t['mag'], filt, min_exp=args.min_exptime)
            if exptime is None:
                continue
            elif exptime > args.max_exptime:
                continue

            status, obsdata = fritz.post_followup_request(allocation_id, obj_id,
                filt, exptime, now, now+TimeDelta(86400*u.s),
                priority=1, maximum_airmass=2, maximum_fwhm=1.5,
                minimum_lunar_distance=15.0, too='N', exposure_count=1,
                observation_type='transient')

            idnum = obsdata['id']
            if status:
                print(f'{obj_id} observation in {filt} successful.  obsid = {idnum}.')
            else:
                print(f'{obj_id} observation in {filt} failed.')


if __name__ == "__main__":
    main()
