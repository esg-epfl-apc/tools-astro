# Licensed under a 3-clause BSD style license - see LICENSE.rst


# put all imports organized as shown below
# 1. standard library imports
import socket
# 2. third party imports
import urllib

from astropy.coordinates import SkyCoord
import requests
from astropy.table import Table, vstack
from astropy.io import fits
from astropy import config as _config
import astropy.utils.data as aud
import time

# 3. local imports - use relative imports
# commonly required local imports shown below as example
# all Query classes should inherit from BaseQuery.
from astroquery.exceptions import NoResultsWarning

# from ..query import BaseQuery
# has common functions required by most modules
# prepend_docstr is a way to copy docstrings between methods
# async_to_sync generates the relevant query tools from _async methods
# import configurable items declared in __init__.py
# from . import conf

import pyvo as vo
from numpy import cos, pi
# export all the public classes and methods
__all__ = ['DESILegacySurvey', 'DESILegacySurveyClass']

# declare global variables and constants if any


# Now begin your main class
# should be decorated with the async_to_sync imported previously
from astropy.io.fits import HDUList


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astroquery.desi`.
    """
    server = _config.ConfigItem(
        ['https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/',
         ],
        'base url')

    timeout = _config.ConfigItem(
        30,
        'Time limit for connecting to template_module server.')


class FileContainer:
    """
    A File Object container, meant to offer lazy access to downloaded FITS
    files.
    """

    def __init__(self, target, **kwargs):
        kwargs.setdefault('cache', True)
        self._target = target
        self._timeout = kwargs.get('remote_timeout', aud.conf.remote_timeout)
        self._readable_object = get_readable_fileobj(target, **kwargs)

    def get_fits(self):
        """
        Assuming the contained file is a FITS file, read it
        and return the file parsed as FITS HDUList
        """
        filedata = self.get_string()

        if len(filedata) == 0:
            raise TypeError("The file retrieved was empty.")

        self._fits = fits.HDUList.fromstring(filedata)

        return self._fits

    def get_string(self):
        """
        Download the file as a string
        """
        if not hasattr(self, '_string'):
            try:
                with self._readable_object as f:
                    data = f.read()
                    self._string = data
            except urllib.error.URLError as e:
                if isinstance(e.reason, socket.timeout):
                    raise TimeoutError("Query timed out, time elapsed {t}s".
                                       format(t=self._timeout))
                else:
                    raise e

        return self._string

    def __repr__(self):
        if hasattr(self, '_fits'):
            return f"Downloaded FITS file: {self._fits!r}"
        else:
            return f"Downloaded object from URL {self._target} with ID {id(self._readable_object)}"


def get_readable_fileobj(*args, **kwargs):
    """
    Overload astropy's get_readable_fileobj so that we can safely monkeypatch
    it in astroquery without affecting astropy core functionality
    """
    return aud.get_readable_fileobj(*args, **kwargs)


class DESILegacySurveyClass():

    """
    Not all the methods below are necessary but these cover most of the common
    cases, new methods may be added if necessary, follow the guidelines at
    <http://astroquery.readthedocs.io/en/latest/api.html>
    """
    # use the Configuration Items imported from __init__.py to set the URL,
    # TIMEOUT, etc.
    conf = Conf()
    URL = conf.server
    TIMEOUT = conf.timeout

    # all query methods are implemented with an "async" method that handles
    # making the actual HTTP request and returns the raw HTTP response, which
    # should be parsed by a separate _parse_result method.   The query_object
    # method is created by async_to_sync automatically.  It would look like
    # this:
    """
    def query_object(object_name, get_query_payload=False)
        response = self.query_object_async(object_name,
                                           get_query_payload=get_query_payload)
        if get_query_payload:
            return response
        result = self._parse_result(response, verbose=verbose)
        return result
    """

    # For services that can query coordinates, use the query_region method.
    # The pattern is similar to the query_object method. The query_region
    # method also has a 'radius' keyword for specifying the radius around
    # the coordinates in which to search. If the region is a box, then
    # the keywords 'width' and 'height' should be used instead. The coordinates
    # may be accepted as an `astropy.coordinates` object or as a string, which
    # may be further parsed.

    # similarly we write a query_region_async method that makes the
    # actual HTTP request and returns the HTTP response

    def query_region_async(self, coordinates, radius,
                           get_query_payload=False, cache=True, data_release=9, use_tap=True):
        """
        Queries a region around the specified coordinates.

        Parameters
        ----------
        coordinates : str or `astropy.coordinates`.
            coordinates around which to query
        radius : str or `astropy.units.Quantity`.
            the radius of the cone search
        get_query_payload : bool, optional
            Just return the dict of HTTP request parameters.
        verbose : bool, optional
            Display VOTable warnings or not.

        Returns
        -------
        response : `requests.Response`
            The HTTP response returned from the service.
            All async methods should return the raw HTTP response.
        """

        if use_tap:
            # TAP query
            # Download tractor catalogue
            url = 'https://datalab.noirlab.edu/tap'
            tap_service = vo.dal.TAPService(url)
            qstr = "SELECT all * FROM ls_dr" + str(data_release) + ".tractor WHERE dec>" + str(coordinates.dec.deg - radius.deg) + " and dec<" + str(
                coordinates.dec.deg + radius.deg) + " and ra>" + str(coordinates.ra.deg - radius.deg / cos(coordinates.dec.deg * pi / 180.)) + " and ra<" + str(
                coordinates.ra.deg + radius.deg / cos(coordinates.dec.deg * pi / 180))

            tap_result = tap_service.run_sync(qstr)
            tap_result = tap_result.to_table()
            # filter out duplicated lines from the table
            mask = tap_result['type'] != 'D'
            filtered_table = tap_result[mask]

            return filtered_table
        else:
            # call the brick list
            table_north = self.query_brick_list(data_release=data_release, emisphere="north")
            table_south = self.query_brick_list(data_release=data_release, emisphere="south")
            # needed columns: ra1, ra2, dec1, dec2 (corners of the bricks), and also brickname
            # radius not used for the moment, but it will be in the future
            # must find the brick within ra1 and ra2
            dec = coordinates.dec.deg
            ra = coordinates.ra.deg

            responses = []

            # north table extraction
            brick_name = table_north['brickname']
            ra1 = table_north['ra1']
            dec1 = table_north['dec1']
            ra2 = table_north['ra2']
            dec2 = table_north['dec2']

            corners1 = SkyCoord(ra1, dec1, unit="deg")
            corners2 = SkyCoord(ra1, dec2, unit="deg")
            corners3 = SkyCoord(ra2, dec1, unit="deg")
            corners4 = SkyCoord(ra2, dec2, unit="deg")

            sep1 = coordinates.separation(corners1)
            sep2 = coordinates.separation(corners2)
            sep3 = coordinates.separation(corners3)
            sep4 = coordinates.separation(corners4)

            t0 = time.time()
            print("Beginning processing bricks northern emishpere")
            for i in range(len(table_north)):
                if ((ra1[i] < ra < ra2[i]) and (dec1[i] < dec < dec2[i])) \
                        or (sep1[i] < radius) or (sep2[i] < radius) or (sep3[i] < radius) or (sep4[i] < radius):
                    # row_north_list.append(table_north[i])
                    brickname = brick_name[i]
                    raIntPart = "{0:03}".format(int(ra1[i]))
                    URL = f"{self.URL}/dr{data_release}/north/tractor/{raIntPart}/tractor-{brickname}.fits"

                    response = requests.get(URL)
                    if response is not None and response.status_code == 200:
                        responses.append(response)
            print("Completion processing bricks northern emishpere, total time: ", time.time() - t0)

            # south table extraction
            brick_name = table_south['brickname']
            ra1 = table_south['ra1']
            dec1 = table_south['dec1']
            ra2 = table_south['ra2']
            dec2 = table_south['dec2']

            corners1 = SkyCoord(ra1, dec1, unit="deg")
            corners2 = SkyCoord(ra1, dec2, unit="deg")
            corners3 = SkyCoord(ra2, dec1, unit="deg")
            corners4 = SkyCoord(ra2, dec2, unit="deg")

            sep1 = coordinates.separation(corners1)
            sep2 = coordinates.separation(corners2)
            sep3 = coordinates.separation(corners3)
            sep4 = coordinates.separation(corners4)

            t0 = time.time()
            print("Beginning processing bricks southern emisphere")
            for i in range(len(table_south)):
                if ((ra1[i] < ra < ra2[i]) and (dec1[i] < dec < dec2[i])) \
                        or (sep1[i] < radius) or (sep2[i] < radius) or (sep3[i] < radius) or (sep4[i] < radius):
                    # row_south_list.append(table_south[i])
                    brickname = brick_name[i]
                    raIntPart = "{0:03}".format(int(ra1[i]))
                    URL = f"{self.URL}/dr{data_release}/south/tractor/{raIntPart}/tractor-{brickname}.fits"

                    response = requests.get(URL)
                    if response is not None and response.status_code == 200:
                        responses.append(response)
            print("Completion processing bricks southern emisphere, total time: ", time.time() - t0)
            print("-----------------------------------------------------")

            return responses

    def query_region(self, coordinates, radius, get_query_payload=False, cache=True, data_release=9, use_tap=True):
        response = self.query_region_async(coordinates, radius, get_query_payload=get_query_payload, cache=cache, data_release=data_release, use_tap=use_tap)
        if get_query_payload:
            return response
        result = self._parse_result(response)
        return result

    def get_images_async(self, position, survey, coordinates=None, data_release=9,
                         projection=None, pixels=None, scaling=None,
                         sampler=None, resolver=None, deedger=None, lut=None,
                         grid=None, gridlabels=None, radius=None, height=None,
                         width=None, cache=True, show_progress=True, image_band='g'):
        """
        Returns
        -------
        A list of context-managers that yield readable file-like objects
        """

        image_size_arcsec = radius.arcsec
        pixsize = 2 * image_size_arcsec / pixels

        image_url = 'https://www.legacysurvey.org/viewer/fits-cutout?ra=' + str(position.ra.deg) + '&dec=' + str(position.dec.deg) + '&size=' + str(
            pixels) + '&layer=ls-dr' + str(data_release) + '&pixscale=' + str(pixsize) + '&bands=' + image_band

        print("image_url: ", image_url)

        file_container = FileContainer(image_url, encoding='binary', show_progress=show_progress)

        try:
            fits_file = file_container.get_fits()
        except (requests.exceptions.HTTPError, urllib.error.HTTPError) as e:
            # TODO not sure this is the most suitable exception
            raise NoResultsWarning(f"{str(e)} - Problem retrieving the file at the url: {str(image_url)}")

        return [fits_file]

    def get_images(self, position, survey, coordinates=None, data_release=9,
                   projection=None, pixels=None, scaling=None,
                   sampler=None, resolver=None, deedger=None, lut=None,
                   grid=None, gridlabels=None, radius=None, height=None,
                   width=None, cache=True, show_progress=True, image_band='g'):
        response = self.get_images_async(position, survey, coordinates=coordinates, data_release=data_release,
                                         projection=projection, pixels=pixels, scaling=scaling,
                                         sampler=sampler, resolver=resolver, deedger=deedger, lut=lut,
                                         grid=grid, gridlabels=gridlabels, radius=radius, height=height,
                                         width=width, cache=cache, show_progress=show_progress, image_band=image_band)
        # if get_query_payload:
        return response
        # result = self._parse_result(response)
        # return result

    # as we mentioned earlier use various python regular expressions, etc
    # to create the dict of HTTP request parameters by parsing the user
    # entered values. For cleaner code keep this as a separate private method:

    def _args_to_payload(self, *args, **kwargs):
        request_payload = dict()
        # code to parse input and construct the dict
        # goes here. Then return the dict to the caller

        return request_payload

    # the methods above call the private _parse_result method.
    # This should parse the raw HTTP response and return it as
    # an `astropy.table.Table`. Below is the skeleton:

    def _parse_result(self, responses):
        tables_list = []
        output_table = Table()

        # handles the cases of query_region with TAP which returns a Table
        # or get_images which should return a list of HDUList, and here in DESILegacySurvey
        # a list of one HDUList object is returned
        if isinstance(responses, Table) or \
                (isinstance(responses, list) and isinstance(responses[0], HDUList)):
            return responses

        # try to parse the result into an astropy.Table, else
        # return the raw result with an informative error message.
        try:
            if not isinstance(responses, list):
                responses = [responses]

            # do something with regex to get the result into
            # astropy.Table form. return the Table.
            # data = io.BytesIO(response.content)
            # table = Table.read(data)

            for r in responses:
                if r.status_code == 200:
                    # TODO figure out on how to avoid writing in a file
                    with open('/tmp/file_content', 'wb') as fin:
                        fin.write(r.content)

                    table = Table.read('/tmp/file_content', hdu=1)
                    tables_list.append(table)

            if len(tables_list) > 0:
                output_table = vstack(tables_list)

        except ValueError:
            # catch common errors here, but never use bare excepts
            # return raw result/ handle in some way
            pass

        return output_table


# the default tool for users to interact with is an instance of the Class
DESILegacySurvey = DESILegacySurveyClass()

# once your class is done, tests should be written
# See ./tests for examples on this

# Next you should write the docs in astroquery/docs/module_name
# using Sphinx.
