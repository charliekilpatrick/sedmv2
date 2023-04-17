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

class skyportal(object):

    def __init__(self, token, instrument='sedmv2'):

        self.auth = {'Authorization': f'token {token}'}
        self.host = 'https://fritz.science'

        # Get target list and instrument list
        self.targets = self.get_objlist()
        self.instruments = self.get_instruments()

        self.inst_id = None
        self.inst_name = None
        if self.instruments is not None:
            for row in self.instruments:
                if row['instrument'].lower()==instrument.lower():
                    self.inst_id = row['id']
                    self.inst_name = row['instrument']
                    break

        if self.inst_id is None or self.inst_name is None:
            raise Exception('ERROR: could not get instrument data')

    def api(self, method, endpoint, data=None, return_type='json'):

        kwargs = {'headers': self.auth}
        if method.lower()=='get':
            kwargs['params']=data
        else:
            kwargs['json']=data

        url = urllib.parse.urljoin(self.host, f'/api/{endpoint}')

        print(f'{method} {url}')
        r = requests.request(method, url, **kwargs)

        if r.status_code!=200:
            return(False, None)
        else:
            if return_type=='json':
                return(True, r.json()['data'])
            else:
                return(True, r.text)

    def get_objlist(self):

        data = {'numPerPage':500, 'group_ids':['1423']}
        status, data = self.api('GET','sources', data=data)
        if status:

            outdata = Table([['X'*100],[0.],[0.],['X'*100],['X'*100]],
                names=('name','ra','dec','group_id','group_names')).copy()[:0]

            for source in data['sources']:
                group_id=','.join([str(g['id']) for g in source['groups']])
                group_name=','.join([str(g['name']) for g in source['groups']])
                outdata.add_row([source['id'],source['ra'],source['dec'],
                    group_id,group_name])

            return(outdata)

        else:
            return(None)

    def get_instruments(self):

        status, data = self.api('GET', 'instrument')
        if status:

            outdata = Table([['X'*100],['X'*100],['X'*100],[0]],
                names=('instrument','filters','telescope','id')).copy()[:0]

            for source in data:
                outdata.add_row([source['name'], ','.join(source['filters']),
                    source['telescope']['nickname'], source['id']])

            return(outdata)

        else:
            return(None)

    def get_observations(self, inst_id, start_date, end_date):

        t0 = Time(start_date)
        t1 = Time(end_date)

        followup_requests = []
        if t1-t0 > 4*86400*u.s:
            n_requests = int((t1-t0)/(4*86400*u.s))+1
            for i in np.arange(n_requests):
                tmin = t0 + TimeDelta(i*4*86400*u.s)
                tmax = t0 + TimeDelta((i+1)*4*86400*u.s)

                data = {'startDate': tmin.datetime.strftime('%Y-%m-%dT%H:%M:%S'),
                        'endDate': tmax.datetime.strftime('%Y-%m-%dT%H:%M:%S'),}

                status, outdata = self.api('GET', f'followup_request', data=data)
                if status:
                    followup_requests.extend(outdata['followup_requests'])

        else:
            data = {'startDate': t0.datetime.strftime('%Y-%m-%dT%H:%M:%S'),
                    'endDate': t1.datetime.strftime('%Y-%m-%dT%H:%M:%S'),}

            status, outdata = self.api('GET', f'followup_request', data=data)
            if status:
                followup_requests.extend(outdata['followup_requests'])

        return(followup_requests)

    def check_if_source_exists(self, obj_id, ra, dec):

        coords = SkyCoord(self.targets['ra'], self.targets['dec'], unit='deg')
        coord = SkyCoord(ra, dec, unit='deg')

        mask = (self.targets['name']==obj_id) & (coord.separation(coords)<2.0*u.arcsec)

        if len(self.targets[mask])==0:
            return(False)
        else:
            return(True)

    def add_new_source(self, obj_id, ra, dec, group_id='1423'):

        data = {
            'id':obj_id,
            'ra':ra,
            'dec':dec,
            'group_ids':[group_id],
        }

        status, data = self.api('POST', 'sources', data=data)

        return(status, data)


    def post_followup_request(self, allocation_id, obj_id,
        filt, exptime, start_date, end_date,
        priority=2, maximum_airmass=2, maximum_fwhm=1.5,
        minimum_lunar_distance=15.0, too='N', exposure_count=1,
        observation_type='transient'):


        payload = {
            'too':too,
            'exposure_time': exptime,
            'exposure_counts': exposure_count,
            'start_date': Time(start_date).datetime.strftime('%Y-%m-%dT%H:%M:%S'),
            'end_date': Time(end_date).datetime.strftime('%Y-%m-%dT%H:%M:%S'),
            'priority': priority,
            'maximum_airmass': maximum_airmass,
            'minimum_lunar_distance': minimum_lunar_distance,
            'observation_type': observation_type,
            'observation_choice': filt,
        }

        data = {'allocation_id': allocation_id,
                'obj_id': obj_id,
                'payload': payload}

        status, data = self.api('POST', 'followup_request', data=data)

        return(status, data)
