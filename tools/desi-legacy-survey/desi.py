# Licensed under a 3-clause BSD style license - see LICENSE.rst


# put all imports organized as shown below
# 1. standard library imports

# 2. third party imports
import urllib.error

import astropy.units as u
import astropy.coordinates as coord
import astropy.io.votable as votable
from astropy.coordinates import SkyCoord
import requests
from astropy.table import Table, vstack
from astropy.io import fits
from astropy import config as _config
import astropy.utils.data as aud
import io
import time
import os
from astropy.config import paths

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
from numpy import pi, cos
# export all the public classes and methods
__all__ = ['DESILegacySurvey', 'DESILegacySurveyClass']

# declare global variables and constants if any


# Now begin your main class
# should be decorated with the async_to_sync imported previously
from astropy.io.fits import HDUList

def suppress_vo_warnings():
    """
    Suppresses all warnings of the class
    `astropy.io.votable.exceptions.VOWarning`.
    """
    warnings.filterwarnings("ignore", category=votable.exceptions.VOWarning)


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


class BaseQuery:
    """
    This is the base class for all the query classes in astroquery. It
    is implemented as an abstract class and must not be directly instantiated.
    """

    def __init__(self):
        S = self._session = requests.Session()
        self._session.hooks['response'].append(self._response_hook)
        S.headers['User-Agent'] = (
            'astroquery/{vers} {olduseragent}'
            .format(vers="",
                    olduseragent=S.headers['User-Agent']))

        self.cache_location = os.path.join(
            paths.get_cache_dir(), 'astroquery',
            self.__class__.__name__.split("Class")[0])
        if not os.path.exists(self.cache_location):
            os.makedirs(self.cache_location)
        self._cache_active = True

    def __call__(self, *args, **kwargs):
        """ init a fresh copy of self """
        return self.__class__(*args, **kwargs)

    def _response_hook(self, response, *args, **kwargs):
        loglevel = log.getEffectiveLevel()

        if loglevel >= 10:
            # Log request at DEBUG severity
            request_hdrs = '\n'.join(f'{k}: {v}' for k, v in response.request.headers.items())
            request_log = textwrap.indent(
                f"-----------------------------------------\n"
                f"{response.request.method} {response.request.url}\n"
                f"{request_hdrs}\n"
                f"\n"
                f"{response.request.body}\n"
                f"-----------------------------------------", '\t')
            log.debug(f"HTTP request\n{request_log}")
        if loglevel >= 5:
            # Log response at super-DEBUG severity
            response_hdrs = '\n'.join(f'{k}: {v}' for k, v in response.headers.items())
            if kwargs.get('stream'):
                response_log = textwrap.indent(
                    f"-----------------------------------------\n"
                    f"{response.status_code} {response.reason} {response.url}\n"
                    f"{response_hdrs}\n"
                    "Streaming Data\n"
                    f"-----------------------------------------", '\t')
            else:
                response_log = textwrap.indent(
                    f"-----------------------------------------\n"
                    f"{response.status_code} {response.reason} {response.url}\n"
                    f"{response_hdrs}\n"
                    f"\n"
                    f"{response.text}\n"
                    f"-----------------------------------------", '\t')
            log.log(5, f"HTTP response\n{response_log}")

    def _request(self, method, url,
                 params=None, data=None, headers=None,
                 files=None, save=False, savedir='', timeout=None, cache=True,
                 stream=False, auth=None, continuation=True, verify=True,
                 allow_redirects=True,
                 json=None, return_response_on_save=False):
        """
        A generic HTTP request method, similar to `requests.Session.request`
        but with added caching-related tools

        This is a low-level method not generally intended for use by astroquery
        end-users.  However, it should _always_ be used by astroquery
        developers; direct uses of `urllib` or `requests` are almost never
        correct.

        Parameters
        ----------
        method : str
            'GET' or 'POST'
        url : str
        params : None or dict
        data : None or dict
        json : None or dict
        headers : None or dict
        auth : None or dict
        files : None or dict
            See `requests.request`
        save : bool
            Whether to save the file to a local directory.  Caching will happen
            independent of this parameter if `BaseQuery.cache_location` is set,
            but the save location can be overridden if ``save==True``
        savedir : str
            The location to save the local file if you want to save it
            somewhere other than `BaseQuery.cache_location`
        timeout : int
        cache : bool
        verify : bool
            Verify the server's TLS certificate?
            (see http://docs.python-requests.org/en/master/_modules/requests/sessions/?highlight=verify)
        continuation : bool
            If the file is partly downloaded to the target location, this
            parameter will try to continue the download where it left off.
            See `_download_file`.
        stream : bool
        return_response_on_save : bool
            If ``save``, also return the server response. The default is to only
            return the local file path.

        Returns
        -------
        response : `requests.Response`
            The response from the server if ``save`` is False
        local_filepath : list
            a list of strings containing the downloaded local paths if ``save``
            is True and ``return_response_on_save`` is False.
        (local_filepath, response) : tuple(list, `requests.Response`)
            a tuple containing a list of strings containing the downloaded local paths,
            and the server response object, if ``save`` is True and ``return_response_on_save``
            is True.
        """
        req_kwargs = dict(
            params=params,
            data=data,
            headers=headers,
            files=files,
            timeout=timeout,
            json=json
        )
        if save:
            local_filename = url.split('/')[-1]
            if os.name == 'nt':
                # Windows doesn't allow special characters in filenames like
                # ":" so replace them with an underscore
                local_filename = local_filename.replace(':', '_')
            local_filepath = os.path.join(savedir or self.cache_location or '.', local_filename)

            response = self._download_file(url, local_filepath, cache=cache,
                                           continuation=continuation, method=method,
                                           allow_redirects=allow_redirects,
                                           auth=auth, **req_kwargs)
            if return_response_on_save:
                return local_filepath, response
            else:
                return local_filepath
        else:
            query = AstroQuery(method, url, **req_kwargs)
            if ((self.cache_location is None) or (not self._cache_active) or (not cache)):
                with suspend_cache(self):
                    response = query.request(self._session, stream=stream,
                                             auth=auth, verify=verify,
                                             allow_redirects=allow_redirects,
                                             json=json)
            else:
                response = query.from_cache(self.cache_location)
                if not response:
                    response = query.request(self._session,
                                             self.cache_location,
                                             stream=stream,
                                             auth=auth,
                                             allow_redirects=allow_redirects,
                                             verify=verify,
                                             json=json)
                    to_cache(response, query.request_file(self.cache_location))
            self._last_query = query
            return response

    def _download_file(self, url, local_filepath, timeout=None, auth=None,
                       continuation=True, cache=False, method="GET",
                       head_safe=False, **kwargs):
        """
        Download a file.  Resembles `astropy.utils.data.download_file` but uses
        the local ``_session``

        Parameters
        ----------
        url : string
        local_filepath : string
        timeout : int
        auth : dict or None
        continuation : bool
            If the file has already been partially downloaded *and* the server
            supports HTTP "range" requests, the download will be continued
            where it left off.
        cache : bool
        method : "GET" or "POST"
        head_safe : bool
        """

        if head_safe:
            response = self._session.request("HEAD", url,
                                             timeout=timeout, stream=True,
                                             auth=auth, **kwargs)
        else:
            response = self._session.request(method, url,
                                             timeout=timeout, stream=True,
                                             auth=auth, **kwargs)

        response.raise_for_status()
        if 'content-length' in response.headers:
            length = int(response.headers['content-length'])
            if length == 0:
                log.warn('URL {0} has length=0'.format(url))
        else:
            length = None

        if ((os.path.exists(local_filepath)
             and ('Accept-Ranges' in response.headers)
             and continuation)):
            open_mode = 'ab'

            existing_file_length = os.stat(local_filepath).st_size
            if length is not None and existing_file_length >= length:
                # all done!
                log.info("Found cached file {0} with expected size {1}."
                         .format(local_filepath, existing_file_length))
                return
            elif existing_file_length == 0:
                open_mode = 'wb'
            else:
                log.info("Continuing download of file {0}, with {1} bytes to "
                         "go ({2}%)".format(local_filepath,
                                            length - existing_file_length,
                                            (length-existing_file_length)/length*100))

                # bytes are indexed from 0:
                # https://en.wikipedia.org/wiki/List_of_HTTP_header_fields#range-request-header
                end = "{0}".format(length-1) if length is not None else ""
                self._session.headers['Range'] = "bytes={0}-{1}".format(existing_file_length,
                                                                        end)

                response = self._session.request(method, url,
                                                 timeout=timeout, stream=True,
                                                 auth=auth, **kwargs)
                response.raise_for_status()
                del self._session.headers['Range']

        elif cache and os.path.exists(local_filepath):
            if length is not None:
                statinfo = os.stat(local_filepath)
                if statinfo.st_size != length:
                    log.warning("Found cached file {0} with size {1} that is "
                                "different from expected size {2}"
                                .format(local_filepath,
                                        statinfo.st_size,
                                        length))
                    open_mode = 'wb'
                else:
                    log.info("Found cached file {0} with expected size {1}."
                             .format(local_filepath, statinfo.st_size))
                    response.close()
                    return
            else:
                log.info("Found cached file {0}.".format(local_filepath))
                response.close()
                return
        else:
            open_mode = 'wb'
            if head_safe:
                response = self._session.request(method, url,
                                                 timeout=timeout, stream=True,
                                                 auth=auth, **kwargs)
                response.raise_for_status()

        blocksize = astropy.utils.data.conf.download_block_size

        log.debug(f"Downloading URL {url} to {local_filepath} with size {length} "
                  f"by blocks of {blocksize}")

        bytes_read = 0

        # Only show progress bar if logging level is INFO or lower.
        if log.getEffectiveLevel() <= 20:
            progress_stream = None  # Astropy default
        else:
            progress_stream = io.StringIO()

        with ProgressBarOrSpinner(
                length, ('Downloading URL {0} to {1} ...'
                         .format(url, local_filepath)),
                file=progress_stream) as pb:
            with open(local_filepath, open_mode) as f:
                for block in response.iter_content(blocksize):
                    f.write(block)
                    bytes_read += blocksize
                    if length is not None:
                        pb.update(bytes_read if bytes_read <= length else
                                  length)
                    else:
                        pb.update(bytes_read)

        response.close()
        return response


class FileContainer:
    """
    A File Object container, meant to offer lazy access to downloaded FITS
    files.
    """

    def __init__(self, target, **kwargs):
        kwargs.setdefault('cache', True)
        self._target = target
        self._timeout = kwargs.get('remote_timeout', aud.conf.remote_timeout)
        if (os.path.splitext(target)[1] == '.fits' and not
                ('encoding' in kwargs and kwargs['encoding'] == 'binary')):
            warnings.warn("FITS files must be read as binaries; error is "
                          "likely.", InputWarning)
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

    def save_fits(self, savepath, link_cache='hard'):
        """
        Save a FITS file to savepath

        Parameters
        ----------
        savepath : str
            The full path to a FITS filename, e.g. "file.fits", or
            "/path/to/file.fits".
        link_cache : 'hard', 'sym', or False
            Try to create a hard or symbolic link to the astropy cached file?
            If the system is unable to create a hardlink, the file will be
            copied to the target location.
        """
        self.get_fits()
        target_key = str(self._target)

        # There has been some internal refactoring in astropy.utils.data
        # so we do this check. Update when minimum required astropy changes.
        if ASTROPY_LT_4_0:
            if not aud.is_url_in_cache(target_key):
                raise IOError("Cached file not found / does not exist.")
            target = aud.download_file(target_key, cache=True)
        else:
            target = aud.download_file(target_key, cache=True, sources=[])

        if link_cache == 'hard':
            try:
                os.link(target, savepath)
            except (IOError, OSError, AttributeError):
                shutil.copy(target, savepath)
        elif link_cache == 'sym':
            try:
                os.symlink(target, savepath)
            except AttributeError:
                raise OSError('Creating symlinks is not possible on this OS.')
        else:
            shutil.copy(target, savepath)

    def get_string(self):
        """
        Download the file as a string
        """
        if not hasattr(self, '_string'):
            try:
                with self._readable_object as f:
                    data = f.read()
                    self._string = data
            except URLError as e:
                if isinstance(e.reason, socket.timeout):
                    raise TimeoutError("Query timed out, time elapsed {t}s".
                                       format(t=self._timeout))
                else:
                    raise e

        return self._string

    def get_stringio(self):
        """
        Return the file as an io.StringIO object
        """
        s = self.get_string()
        # TODO: replace with six.BytesIO
        try:
            return six.BytesIO(s)
        except TypeError:
            return six.StringIO(s)

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


def parse_votable(content):
    """
    Parse a votable in string format
    """
    tables = votable.parse(six.BytesIO(content), pedantic=False)
    return tables

class DESILegacySurveyClass(BaseQuery):

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

    def query_object_async(self, object_name, get_query_payload=False,
                           cache=True, data_release=9):
        """
        This method is for services that can parse object names. Otherwise
        use :meth:`astroquery.template_module.TemplateClass.query_region`.
        Put a brief description of what the class does here.

        Parameters
        ----------
        object_name : str
            name of the identifier to query.
        get_query_payload : bool, optional
            This should default to False. When set to `True` the method
            should return the HTTP request parameters as a dict.
        verbose : bool, optional
           This should default to `False`, when set to `True` it displays
           VOTable warnings.
        any_other_param : <param_type>
            similarly list other parameters the method takes

        Returns
        -------
        response : `requests.Response`
            The HTTP response returned from the service.
            All async methods should return the raw HTTP response.

        Examples
        --------
        While this section is optional you may put in some examples that
        show how to use the method. The examples are written similar to
        standard doctests in python.

        """
        # the async method should typically have the following steps:
        # 1. First construct the dictionary of the HTTP request params.
        # 2. If get_query_payload is `True` then simply return this dict.
        # 3. Else make the actual HTTP request and return the corresponding
        #    HTTP response
        # All HTTP requests are made via the `BaseQuery._request` method. This
        # use a generic HTTP request method internally, similar to
        # `requests.Session.request` of the Python Requests library, but
        # with added caching-related tools.

        # See below for an example:

        # first initialize the dictionary of HTTP request parameters
        request_payload = dict()

        # Now fill up the dictionary. Here the dictionary key should match
        # the exact parameter name as expected by the remote server. The
        # corresponding dict value should also be in the same format as
        # expected by the server. Additional parsing of the user passed
        # value may be required to get it in the right units or format.
        # All this parsing may be done in a separate private `_args_to_payload`
        # method for cleaner code.

        # request_payload['object_name'] = object_name
        # similarly fill up the rest of the dict ...

        if get_query_payload:
            return request_payload
        # BaseQuery classes come with a _request method that includes a
        # built-in caching system

        # TODO: implement here http query as needed
        # e.g. I suspect we get files like this one: https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/north/tractor/000/tractor-0001m002.fits
        # to confirm with AG, AN
        # if so:

        URL = f"{self.URL}/dr{data_release}/north/tractor/000/tractor-0001m002.fits"

        response = self._request('GET', URL, params={},
                                 timeout=self.TIMEOUT, cache=cache)
        return response

    def query_brick_list_async(self, data_release=9, get_query_payload=False, emisphere="north",
                           cache=True):
        """

        """
        request_payload = dict()

        if get_query_payload:
            return request_payload
        URL = f"{self.URL}/dr{data_release}/{emisphere}/survey-bricks-dr{data_release}-{emisphere}.fits.gz"
        # TODO make it work with the original request
        # response = self._request('GET', URL, params={},
        #                          timeout=self.TIMEOUT, cache=cache)

        response = requests.get(URL)

        print("completed fits file request")

        return response

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

    def _parse_result(self, responses, verbose=False):
        tables_list = []
        output_table = Table()

        # handles the cases of query_region with TAP which returns a Table
        # or get_images which should return a list of HDUList, and here in DESILegacySurvey
        # a list of one HDUList object is returned
        if isinstance(responses, Table) or \
                (isinstance(responses, list) and isinstance(responses[0], HDUList)):
            return responses

        # if verbose is False then suppress any VOTable related warnings
        if not verbose:
            commons.suppress_vo_warnings()
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