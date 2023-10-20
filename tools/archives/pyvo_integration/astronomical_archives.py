import errno
import functools
import json
import os
import sys
import signal
import urllib
from urllib import request

from astropy.coordinates import SkyCoord

import pyvo
from pyvo import DALAccessError, DALQueryError, DALServiceError
from pyvo import registry


MAX_ALLOWED_ENTRIES = 100
MAX_REGISTRIES_TO_SEARCH = 100


class TimeoutException(Exception):
    pass


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutException(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator


class Service:
    # https://pyvo.readthedocs.io/en/latest/api/pyvo.registry.Servicetype.html

    services = {
        'TAP': 'tap',
        'SIA': 'sia',
        'SIA2': 'sia2',
        'SPECTRUM': 'spectrum',
        'SCS': 'scs',
        'LINE': 'line'
    }

    supported_services = {
        'TAP': 'tap'
    }

    def __init__(self):
        pass

    @staticmethod
    def is_service_supported(service_type) -> bool:
        is_supported = True

        if service_type not in Service.services.keys():
            is_supported = False
        elif service_type not in Service.supported_services.keys():
            is_supported = False

        return is_supported


class Waveband:
    # https://pyvo.readthedocs.io/en/latest/api/pyvo.registry.Waveband.html
    # https://www.ivoa.net/rdf/messenger/2020-08-26/messenger.html

    wavebands = {
        'Extreme UV': 'EUV',
        'Gamma ray': 'Gamma-ray',
        'Infrared': 'Infrared',
        'Millimeter': 'Millimeter',
        'Neutrino': 'Neutrino',
        'Optical': 'Optical',
        'Photon': 'Photon',
        'Radio': 'Radio',
        'Ultra violet': 'UV',
        'X-ray': 'X-ray'
    }

    def __init__(self):
        pass

    @staticmethod
    def is_waveband_supported(waveband) -> bool:
        is_supported = True

        if waveband not in Waveband.wavebands.keys():
            is_supported = False

        return is_supported


class TapArchive:
    # https://www.ivoa.net/documents/ObsCore/20170509/REC-ObsCore-v1.1-20170509

    service_type = Service.services['TAP']

    def __init__(self,
                 id=1,
                 title="Unknown title",
                 name="Unknown name",
                 access_url=""):

        self.id = id,
        self.title = title,
        self.name = name,
        self.access_url = access_url
        self.initialized = False
        self.archive_service = None
        self.tables = None

    @timeout(10)
    def get_resources(self,
                      query,
                      number_of_results,
                      url_field='access_url'):

        resource_list_hydrated = []

        error_message = None

        if self.initialized:

            try:
                raw_resource_list = self.archive_service.search(query)

                for i, resource in enumerate(raw_resource_list):
                    if i < number_of_results:
                        resource_list_hydrated.append(
                            self._get_resource_object(resource))
                    else:
                        break

            except DALQueryError:
                if self.has_obscore_table():
                    error_message = "Error in query -> " + query
                    Logger.create_action_log(
                        Logger.ACTION_ERROR,
                        Logger.ACTION_TYPE_DOWNLOAD,
                        error_message)
                else:
                    error_message = "No obscore table in the archive"
                    Logger.create_action_log(
                        Logger.ACTION_ERROR,
                        Logger.ACTION_TYPE_DOWNLOAD,
                        error_message)

            except DALServiceError:
                error_message = "Error communicating with the service"
                Logger.create_action_log(
                    Logger.ACTION_ERROR,
                    Logger.ACTION_TYPE_DOWNLOAD,
                    error_message)

            except Exception:
                error_message = "Unknow error while querying the service"
                Logger.create_action_log(
                    Logger.ACTION_ERROR,
                    Logger.ACTION_TYPE_DOWNLOAD,
                    error_message)

        return resource_list_hydrated, error_message

    def _get_resource_object(self, resource):
        resource_hydrated = {}

        for key, value in resource.items():
            resource_hydrated[key] = value

        return resource_hydrated

    def initialize(self):
        error_message = None

        try:
            self._get_service()

            if self.archive_service:
                self._set_archive_tables()
                self.initialized = True

        except DALAccessError:
            error_message = \
                "A connection to the service could not be established"
            Logger.create_action_log(
                Logger.ACTION_ERROR,
                Logger.ACTION_TYPE_ARCHIVE_CONNECTION,
                error_message)

        except Exception:
            error_message = "Unknow error while initializing TAP service"
            Logger.create_action_log(
                Logger.ACTION_ERROR,
                Logger.ACTION_TYPE_ARCHIVE_CONNECTION,
                error_message)

        return self.initialized, error_message

    def _get_service(self):
        if self.access_url:
            self.archive_service = pyvo.dal.TAPService(self.access_url)

    def _set_archive_tables(self):

        self.tables = []

        for table in self.archive_service.tables:
            archive_table = {
                'name': table.name,
                'type': table.type,
                'fields': None
            }

            fields = []

            for table_field in table.columns:
                field = {
                    'name': table_field.name,
                    'description': table_field.description,
                    'unit': table_field.unit,
                    'datatype': table_field.datatype.content
                }

                fields.append(field)

            archive_table['fields'] = fields

            self.tables.append(archive_table)

    def _is_query_valid(self, query) -> bool:
        is_valid = True

        attribute_from = 'from'
        attribute_where = 'where'

        idx_from = query.index(attribute_from)
        idx_where = query.index(attribute_where)

        table_name = ''

        for idx in range(idx_from + len('from') + 1, idx_where):
            table_name = table_name + query[idx]

        if not next(
                (item for item in self.tables if
                 item["name"] == table_name),
                False):

            is_valid = False

        return is_valid

    def has_obscore_table(self) -> bool:
        has_obscore_table = self._has_table("ivoa.obscore")

        return has_obscore_table

    def _has_table(self, table_name) -> bool:
        _has_table = False

        _has_table = next(
            (item for item in self.tables if item["name"] == table_name),
            False)

        return _has_table

    def get_archive_name(self, archive_type):
        try:
            if archive_type == 'registry':
                name = str(self.title).strip("',()")
            else:
                name = self.access_url
        except Exception:
            name = 'Unknown archive title'

        return name


class ConeService(TapArchive):

    def _get_service(self):
        if self.access_url:
            self.archive_service = pyvo.dal.SCSService(self.access_url)

    def get_resources_from_service_list(self, service_list, target, radius):

        resource_list_hydrated = []

        for service in service_list:
            resources = service.search(target, radius)
            for i in range(resources.__len__()):
                resource_url = resources.getrecord(i).getdataurl()
                if resource_url:
                    resource_list_hydrated.append(resource_url)

        return resource_list_hydrated


class RegistrySearchParameters:

    def __init__(self, keyword=None, waveband=None, service_type=None):
        self.keyword = keyword
        self.waveband = waveband
        self.service_type = service_type

    def get_parameters(self):

        parameters = {
            'keywords': '',
            'waveband': '',
            'service_type': ''
        }

        if self.keyword:
            parameters['keywords'] = self.keyword

        if Waveband.is_waveband_supported(self.waveband):
            parameters['waveband'] = \
                Waveband.wavebands[self.waveband]

        if Service.is_service_supported(self.service_type):
            parameters['service_type'] = \
                Service.services[self.service_type]
        else:
            parameters['service_type'] = Service.services['TAP']

        return parameters


class Registry:

    def __init__(self):
        pass

    @staticmethod
    def search_registries(rsp: RegistrySearchParameters,
                          number_of_registries):

        parameters = rsp.get_parameters()

        keywords = parameters['keywords']
        waveband = parameters['waveband']
        service_type = parameters['service_type']

        if not waveband:
            registry_list = registry.search(
                keywords=keywords,
                servicetype=service_type)
        else:
            registry_list = registry.search(
                keywords=keywords,
                waveband=waveband,
                servicetype=service_type)

        if registry_list:
            registry_list = Registry._get_registries_from_list(
                registry_list,
                number_of_registries)

        return registry_list

    @staticmethod
    def _get_registries_from_list(registry_list, number_of_registries):

        archive_list = []

        for i, ivoa_registry in enumerate(registry_list):
            if i < number_of_registries:
                archive = TapArchive(ivoa_registry.standard_id,
                                     ivoa_registry.res_title,
                                     ivoa_registry.short_name,
                                     ivoa_registry.access_url)

                archive_list.append(archive)

        return archive_list


class ConeServiceRegistry:

    def __init__(self):
        pass

    @staticmethod
    def search_services(keyword, number_of_registries):

        service_list = []

        service_list = registry.search(servicetype="scs", keywords=keyword)

        if service_list:
            service_list = service_list[:number_of_registries]

        return service_list


class TapQuery:

    def __init__(self, query):
        self.raw_query = query

    def get_query(self):
        return urllib.parse.unquote(self.raw_query).replace("+", " ")


class BaseADQLQuery:

    def __init__(self):
        pass

    def _get_order_by_clause(self, order_type):
        order_by_clause = 'ORDER BY ' + order_type

        return order_by_clause

    def _get_where_clause(self, parameters):
        where_clause = ''
        is_first_statement = True

        for key, value in parameters.items():

            if value != '':
                statement = str(key) + ' = ' + '\'' + str(value) + '\' '

                if is_first_statement:
                    is_first_statement = False
                    where_clause += 'WHERE '
                else:
                    statement = 'AND ' + statement

                where_clause += statement

        return where_clause


class ToolRunner:

    def __init__(self,
                 run_parameters,
                 output,
                 output_csv,
                 output_html,
                 output_basic_html,
                 output_error):

        self._raw_parameters_path = run_parameters
        self._json_parameters = json.load(open(run_parameters, "r"))
        self._archive_type = ''
        self._query_type = ''
        self._archives = []
        self._adql_query = ''
        self._services_access_url = ''
        self._url_field = 'access_url'
        self._number_of_files = ''
        self._is_initialised = False

        self._csv_file = False
        self._image_file = False
        self._html_file = False
        self._basic_html_file = False

        self._output = output
        self._output_csv = output_csv
        self._output_html = output_html
        self._output_basic_html = output_basic_html
        self._output_error = output_error

        self._set_run_main_parameters()

        self._is_initialised, error_message = self._set_archive()

        if self._is_initialised and error_message is None:
            self._set_query()
            self._set_output()

    def _set_run_main_parameters(self):

        qs = "query_section"
        qsl = "query_selection"

        self._archive_type = \
            self._json_parameters['archive_selection']['archive_type']
        self._query_type = \
            self._json_parameters[qs][qsl]['query_type']

    def _set_archive(self):

        error_message = None

        if self._archive_type == 'archive':
            self._service_access_url =\
                self._json_parameters['archive_selection']['archive']

            self._archives.append(
                TapArchive(access_url=self._service_access_url))

        else:
            keyword = \
                self._json_parameters['archive_selection']['keyword']
            waveband = \
                self._json_parameters['archive_selection']['wavebands']
            service_type = \
                self._json_parameters['archive_selection']['service_type']

            rsp = RegistrySearchParameters(
                keyword=keyword,
                waveband=waveband,
                service_type=service_type)

            archive_list = Registry.search_registries(
                rsp,
                MAX_REGISTRIES_TO_SEARCH)

            if len(archive_list) >= 1:
                self._archives = archive_list
            else:
                error_message = "no archive matching search parameters"
                Logger.create_action_log(
                    Logger.ACTION_ERROR,
                    Logger.ACTION_TYPE_ARCHIVE_CONNECTION,
                    error_message)

        if error_message is None:

            self._archives[:] = \
                [archive for archive in self._archives if
                 archive.initialize()[0]]

            if len(self._archives) >= 1:
                return True, None
            else:
                return False, \
                    "no archive matching search" \
                    " parameters could be initialized"

        else:
            return False, error_message

    def _set_cone_service(self):

        qs = 'query_section'
        qsl = 'query_selection'
        csts = 'cone_search_target_selection'

        error_message = None
        is_service_initialised = True

        keyword = self._json_parameters[qs][qsl][csts]['keyword']

        service_list = ConeServiceRegistry.search_services(
            keyword,
            MAX_REGISTRIES_TO_SEARCH)

        if len(service_list) >= 1:
            self._services = service_list
        else:
            is_service_initialised = False
            error_message = "no services matching search parameters"
            Logger.create_action_log(
                Logger.ACTION_ERROR,
                Logger.ACTION_TYPE_ARCHIVE_CONNECTION,
                error_message)

        return is_service_initialised, error_message

    def _set_query(self):

        qs = 'query_section'
        qsl = 'query_selection'
        csts = 'cone_search_target_selection'
        cs = 'cone_section'
        ts = 'target_selection'
        con = 'cone_object_name'

        if self._query_type == 'obscore_query':

            dataproduct_type = \
                self._json_parameters[qs][qsl]['dataproduct_type']
            obs_collection = \
                self._json_parameters[qs][qsl]['obs_collection']
            obs_title = \
                self._json_parameters[qs][qsl]['obs_title']
            obs_id = \
                self._json_parameters[qs][qsl]['obs_id']
            facility_name = \
                self._json_parameters[qs][qsl]['facility_name']
            instrument_name = \
                self._json_parameters[qs][qsl]['instrument_name']
            em_min = \
                self._json_parameters[qs][qsl]['em_min']
            em_max = \
                self._json_parameters[qs][qsl]['em_max']
            target_name = \
                self._json_parameters[qs][qsl]['target_name']
            obs_publisher_id = \
                self._json_parameters[qs][qsl]['obs_publisher_id']
            s_fov = \
                self._json_parameters[qs][qsl]['s_fov']
            calibration_level = \
                self._json_parameters[qs][qsl]['calibration_level']
            t_min = \
                self._json_parameters[qs][qsl]['t_min']
            t_max = \
                self._json_parameters[qs][qsl]['t_max']
            order_by = \
                self._json_parameters[qs][qsl]['order_by']

            if self._json_parameters[qs][qsl][cs][csts][ts] == 'coordinates':
                ra = self._json_parameters[qs][qsl][cs][csts]['ra']
                dec = self._json_parameters[qs][qsl][cs][csts]['dec']
            else:
                obs_target = self._json_parameters[qs][qsl][cs][csts][con]

                if obs_target != 'none' and obs_target is not None:
                    target = CelestialObject(obs_target)
                    target_coordinates = target.get_coordinates_in_degrees()

                    ra = target_coordinates['ra']
                    dec = target_coordinates['dec']
                else:
                    ra = None
                    dec = None

            radius = self._json_parameters[qs][qsl][cs]['radius']

            if (ra != '' and ra is not None)\
                    and (dec != '' and dec is not None)\
                    and (radius != '' and radius is not None):
                cone_condition = \
                    ADQLConeSearchQuery.get_search_circle_condition(ra,
                                                                    dec,
                                                                    radius)
            else:
                cone_condition = None

            obscore_query_object = ADQLObscoreQuery(dataproduct_type,
                                                    obs_collection,
                                                    obs_title,
                                                    obs_id,
                                                    facility_name,
                                                    instrument_name,
                                                    em_min,
                                                    em_max,
                                                    target_name,
                                                    obs_publisher_id,
                                                    s_fov,
                                                    calibration_level,
                                                    t_min,
                                                    t_max,
                                                    cone_condition,
                                                    order_by)

            self._adql_query = obscore_query_object.get_query()

        elif self._query_type == 'raw_query':

            wc = 'where_clause'

            tap_table = \
                self._json_parameters[qs][qsl]['table']

            where_field = \
                self._json_parameters[qs][qsl][wc]['where_field']
            where_condition = \
                self._json_parameters[qs][qsl][wc]['where_condition']

            self._url_field = \
                self._json_parameters[qs][qsl]['url_field']

            self._adql_query = \
                ADQLTapQuery().get_query(
                    tap_table,
                    where_field,
                    where_condition)
        else:
            self._adql_query = ADQLObscoreQuery.base_query

    def _set_cone_query(self):

        qs = 'query_section'
        qsl = 'query_selection'
        csts = 'cone_search_target_selection'
        ts = 'target_selection'
        con = 'cone_object_name'

        search_radius = self._json_parameters[qs][qsl]['radius']
        time = None

        if self._json_parameters[qs][qsl][csts][ts] == 'coordinates':
            ra = self._json_parameters[qs][qsl][csts]['ra']
            dec = self._json_parameters[qs][qsl][csts]['dec']
            time = self._json_parameters[qs][qsl][csts]['time']
        else:
            target = CelestialObject(self._json_parameters[qs][qsl][csts][con])

            target_coordinates = target.get_coordinates_in_degrees()

            ra = target_coordinates['ra']
            dec = target_coordinates['dec']

        cone_query_object = ADQLConeSearchQuery(ra, dec, search_radius, time)

        self._adql_query = cone_query_object.get_query()

    def _set_output(self):
        self._number_of_files = \
            int(
                self._json_parameters['output_section']['number_of_files']
            )

        if self._number_of_files < 1:
            self._number_of_files = 1
        elif self._number_of_files > 100:
            self._number_of_files = MAX_ALLOWED_ENTRIES

        output_selection = \
            self._json_parameters['output_section']['output_selection']

        if output_selection is not None:
            if 'c' in output_selection:
                self._csv_file = True
            if 'i' in output_selection:
                self._image_file = True
            if 'h' in output_selection:
                self._html_file = True
            if 'b' in output_selection:
                self._basic_html_file = True

    def _validate_json_parameters(self, json_parameters):
        self._json_parameters = json.load(open(json_parameters, "r"))

    def run(self):
        if self._is_initialised:
            error_message = None
            file_url = []

            archive_name = self._archives[0].get_archive_name(
                self._archive_type)

            for archive in self._archives:
                try:
                    _file_url, error_message = archive.get_resources(
                        self._adql_query,
                        self._number_of_files,
                        self._url_field)

                    file_url.extend(_file_url)
                except TimeoutException:
                    error_message = \
                        "Archive is taking too long to respond (timeout)"
                    Logger.create_action_log(
                        Logger.ACTION_ERROR,
                        Logger.ACTION_TYPE_ARCHIVE_CONNECTION,
                        error_message)

                if len(file_url) >= int(self._number_of_files):
                    file_url = file_url[:int(self._number_of_files)]
                    break

            if file_url:

                if self._csv_file:
                    FileHandler.write_urls_to_output(
                        file_url,
                        self._output_csv,
                        self._url_field)

                if self._image_file:

                    try:
                        fits_file = FileHandler.download_file_from_url(
                            file_url[0][self._url_field])

                        FileHandler.write_file_to_output(
                            fits_file,
                            self._output, "wb")

                        log_message = "from url " +\
                                      file_url[0][self._url_field]

                        Logger.create_action_log(
                            Logger.ACTION_SUCCESS,
                            Logger.ACTION_TYPE_DOWNLOAD,
                            log_message)

                    except Exception:
                        error_message = "from url " + \
                                        file_url[0][self._url_field]

                        Logger.create_action_log(
                            Logger.ACTION_ERROR,
                            Logger.ACTION_TYPE_DOWNLOAD,
                            error_message)

                    for i, url in enumerate(file_url[1:], start=1):
                        try:
                            fits_file = \
                                FileHandler.download_file_from_url(
                                    url[self._url_field])

                            FileHandler.write_file_to_subdir(
                                fits_file,
                                FileHandler.get_file_name_from_url(
                                    url[self._url_field]))

                            log_message = "from url " + \
                                          url[self._url_field]

                            Logger.create_action_log(
                                Logger.ACTION_SUCCESS,
                                Logger.ACTION_TYPE_DOWNLOAD,
                                log_message)

                        except Exception:
                            error_message = "from url " + \
                                            url[self._url_field]

                            Logger.create_action_log(
                                Logger.ACTION_ERROR,
                                Logger.ACTION_TYPE_DOWNLOAD,
                                error_message)

                if self._html_file:
                    html_file = OutputHandler.generate_html_output(
                        file_url,
                        archive_name,
                        self._adql_query)

                    FileHandler.write_file_to_output(html_file,
                                                     self._output_html)

                if self._basic_html_file:
                    html_file = \
                        OutputHandler.generate_basic_html_output(
                            file_url,
                            archive_name,
                            self._adql_query)

                    FileHandler.write_file_to_output(
                        html_file,
                        self._output_basic_html)

                summary_file = Logger.create_log_file(archive_name,
                                                      self._adql_query)
                summary_file += "\n Tool run executed with success"

                FileHandler.write_file_to_output(summary_file,
                                                 self._output_error)

            else:

                summary_file = Logger.create_log_file(archive_name,
                                                      self._adql_query)

                if error_message is None:
                    summary_file += \
                        "\n No resources matching parameters found"
                else:
                    summary_file += error_message

                FileHandler.write_file_to_output(summary_file,
                                                 self._output_error)
        else:
            summary_file = Logger.create_log_file("Archive",
                                                  self._adql_query)

            summary_file += "Unable to initialize archives"

            FileHandler.write_file_to_output(summary_file,
                                             self._output_error)


class ADQLObscoreQuery(BaseADQLQuery):
    order_by_field = {
        'size': 'access_estsize',
        'collection': 'obs_collection',
        'object': 'target_name'
    }

    base_query = 'SELECT TOP ' + \
                 str(MAX_ALLOWED_ENTRIES) + \
                 ' * FROM ivoa.obscore '

    def __init__(self,
                 dataproduct_type,
                 obs_collection,
                 obs_title,
                 obs_id,
                 facility_name,
                 instrument_name,
                 em_min,
                 em_max,
                 target_name,
                 obs_publisher_id,
                 s_fov,
                 calibration_level,
                 t_min,
                 t_max,
                 cone_condition,
                 order_by):

        super().__init__()

        if calibration_level == 'none':
            calibration_level = ''

        if order_by == 'none':
            order_by = ''

        if t_min == 'None' or t_min is None:
            t_min = ''

        if t_max == 'None' or t_max is None:
            t_max = ''

        if em_min == 'None' or em_min is None:
            em_min = ''

        if em_max == 'None' or em_max is None:
            em_max = ''

        if dataproduct_type == 'none' or dataproduct_type is None:
            dataproduct_type = ''

        if cone_condition is not None:
            self.cone_condition = cone_condition
        else:
            self.cone_condition = None

        self.parameters = {
            'dataproduct_type': dataproduct_type,
            'obs_collection': obs_collection,
            'obs_title': obs_title,
            'obs_id': obs_id,
            'facility_name': facility_name,
            'instrument_name': instrument_name,
            'em_min': em_min,
            'em_max': em_max,
            'target_name': target_name,
            'obs_publisher_id': obs_publisher_id,
            's_fov': s_fov,
            'calib_level': calibration_level,
            't_min': t_min,
            't_max': t_max,
        }

        self.order_by = order_by

    def get_query(self):
        return ADQLObscoreQuery.base_query + \
            self.get_where_statement() + \
            self.get_order_by_statement()

    def get_order_by_statement(self):
        if self.order_by != '':
            return self._get_order_by_clause(self.order_by)
        else:
            return ''

    def _get_order_by_clause(self, order_type):

        obscore_order_type = ADQLObscoreQuery.order_by_field[order_type]

        return super()._get_order_by_clause(obscore_order_type)

    def get_where_statement(self):
        where_clause = self._get_where_clause(self.parameters)

        if where_clause == '' and self.cone_condition is not None:
            where_clause = 'WHERE ' + self.get_cone_condition()
        elif where_clause != '' and self.cone_condition is not None:
            where_clause += 'AND ' + self.get_cone_condition()

        return where_clause

    def _get_where_clause(self, parameters):
        return super()._get_where_clause(parameters)

    def get_cone_condition(self):
        return self.cone_condition


class ADQLTapQuery(BaseADQLQuery):
    base_query = 'SELECT TOP '+str(MAX_ALLOWED_ENTRIES)+' * FROM '

    def __init__(self):
        super().__init__()

    def get_order_by_clause(self, order_type):
        return super()._get_order_by_clause(order_type)

    def get_query(self, table, where_field, where_condition):
        if where_field != '' and where_condition != '':
            return ADQLTapQuery.base_query + \
                str(table) + \
                ' WHERE ' + \
                str(where_field) + ' = ' + '\'' + \
                str(where_condition) + '\''
        else:
            return ADQLTapQuery.base_query + str(table)


class ADQLConeSearchQuery:

    base_query = "SELECT TOP 100 * FROM ivoa.obscore"

    def __init__(self, ra, dec, radius, time=None):

        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.time = time

        self._query = ADQLObscoreQuery.base_query

        if self.ra and self.dec and self.radius:
            self._query += " WHERE "
            self._query += self._get_search_circle(ra, dec, radius)

            if self.time:
                self._query += self._get_search_time()

    def _get_search_circle(self, ra, dec, radius):
        return "(CONTAINS" \
               "(POINT('ICRS', s_ra, s_dec), " \
               "CIRCLE('ICRS', "+str(ra)+", "+str(dec)+", "+str(radius)+")" \
               ") = 1)"

    def _get_search_time(self):
        return " AND t_min <= "+self.time+" AND t_max >= "+self.time

    def get_query(self):
        return self._query

    @staticmethod
    def get_search_circle_condition(ra, dec, radius):
        return "(CONTAINS" \
               "(POINT('ICRS', s_ra, s_dec)," \
               "CIRCLE('ICRS', "+str(ra)+", "+str(dec)+", "+str(radius)+")" \
               ") = 1) "


class CelestialObject:

    def __init__(self, name):
        self.name = name
        self.coordinates = None

        self.coordinates = SkyCoord.from_name(self.name)

    def get_coordinates_in_degrees(self):

        coordinates = {
            'ra': '',
            'dec': ''
        }

        ra_dec = self.coordinates.ravel()

        coordinates['ra'] = ra_dec.ra.degree[0]
        coordinates['dec'] = ra_dec.dec.degree[0]

        return coordinates


class HTMLReport:
    _html_report_base_header = ''
    _html_report_base_body = ''
    _html_report_base_footer = ''
    _html_report_base_script = ''

    def __init__(self):
        pass


class OutputHandler:

    def __init__(self):
        pass

    @staticmethod
    def generate_html_output(urls_data, archive_name, adql_query):
        return OutputHandler.html_header + \
            OutputHandler.generate_html_content(
                urls_data,
                archive_name,
                adql_query,
                div_attr='class="title"',
                table_attr='class="fl-table"')

    @staticmethod
    def generate_basic_html_output(urls_data,
                                   archive_name,
                                   adql_query, ):
        return OutputHandler.generate_html_content(urls_data,
                                                   archive_name,
                                                   adql_query)

    @staticmethod
    def generate_html_content(urls_data, archive_name, adql_query,
                              div_attr="", table_attr="border='1'"):
        html_file = \
            f"""
                    <div {div_attr}>
                        <h2>Resources Preview archive:
                            <span>
                                {archive_name}
                            </span>
                        </h2>
                        <span>ADQL query : {adql_query}</span>
                    </div>"""

        html_file += f'<table {table_attr}><thead><tr>'

        for key in Utils.collect_resource_keys(urls_data):
            html_file += '<th>' + str(key) + '</th>'

        html_file += '</thead></tr><tbody>'

        for resource in urls_data:
            html_file += '<tr>'

            for key, value in resource.items():
                html_file += f'<td>{value}</td>'

            html_file += '<td>'
            for preview_key in \
                    ['preview', 'preview_url', 'postcard_url']:
                if preview_key in resource:
                    html_file += (
                        '<details><summary>Preview</summary>'
                        f'<img src="{resource[preview_key]}"/>'
                        '</details>'
                    )
            html_file += '</td>'
            html_file += '</tr>'

        html_file += '</tbody></table>'
        return html_file

    html_header = """ <head><style>

                    details {
                        padding: 10px;
                    }

                    .table-wrapper {
                        margin: 10px 70px 70px;
                        box-shadow: 0px 35px 50px rgba( 0, 0, 0, 0.2 );
                    }

                    .fl-table {
                        border-radius: 5px;
                        font-size: 12px;
                        font-weight: normal;
                        border: none;
                        border-collapse: collapse;
                        width: 100%;
                        max-width: 100%;
                        white-space: nowrap;
                        background-color: white;
                    }

                    .fl-table td, .fl-table th {
                        text-align: center;
                        padding: 8px;
                    }

                    .fl-table td {
                        border: 1px solid #999999;
                        font-size: 15px;
                    }

                    .fl-table thead th {
                        color: #ffffff;
                        background: #4FC3A1;
                        border: 1px solid #999999;
                    }


                    .fl-table thead th:nth-child(odd) {
                        color: #ffffff;
                        background: #324960;
                    }

                    .fl-table tr:nth-child(even) {
                        background: #F8F8F8;
                    }

                    .title h2 {
                      text-align: center;
                      font-size: 22px;
                      font-weight: 700; color:#202020;
                      text-transform: uppercase;
                      word-spacing: 1px; letter-spacing:2px;
                      margin-bottom: 50px;
                    }

                    .title h2 span {
                      padding-top: 40px;
                      text-transform: none;
                      font-size:.80em;
                      font-weight: bold;
                      font-family: "Playfair Display","Bookman",serif;
                      color:#999;
                      letter-spacing:-0.005em;
                      word-spacing:1px;
                      letter-spacing:none;
                    }

                    .title h1:before {
                      background-color: #dfdfdf;
                    }

                </style></head>"""


class FileHandler:

    def __init__(self):
        pass

    @staticmethod
    def download_file_from_url(file_url):
        with request.urlopen(file_url) as response:
            fits_file = response.read()

        return fits_file

    @staticmethod
    def write_file_to_output(file, output, write_type="w"):
        with open(output, write_type) as file_output:
            file_output.write(file)

    @staticmethod
    def write_urls_to_output(urls: [], output, access_url="access_url"):
        with open(output, "w") as file_output:
            for url in urls:
                try:
                    file_output.write(url[access_url] + ',')
                except Exception:
                    error_message = "url field not found for url"
                    Logger.create_action_log(
                        Logger.ACTION_ERROR,
                        Logger.ACTION_TYPE_WRITE_URL,
                        error_message)

    @staticmethod
    def write_file_to_subdir(file, index):
        dir = os.getcwd()

        dir += '/fits'

        upload_dir = os.path.join(dir, str(index) + '.fits')

        with open(upload_dir, "wb") as file_output:
            file_output.write(file)

    @staticmethod
    def get_file_name_from_url(url, index=None):
        url_parts = url.split('/')

        file_name = ''

        try:
            if (url_parts[-1]) != '':
                file_name = url_parts[-1]
            elif len(url_parts) > 1:
                file_name = url_parts[-2]
        except Exception:
            file_name = 'archive file '

        return file_name


class Utils:

    def __init__(self):
        pass

    @staticmethod
    def collect_resource_keys(urls_data: list) -> list:
        """
        Collect all the keys from the resources,
        keeping the order in the order of key appearance in the resources
        """

        resource_keys = []
        for resource in urls_data:
            for key in resource.keys():
                if key not in resource_keys:
                    resource_keys.append(key)
        return resource_keys


class Logger:
    _logs = []

    ACTION_SUCCESS = 1
    ACTION_ERROR = 2

    ACTION_TYPE = 1
    INFO_TYPE = 2

    ACTION_TYPE_DOWNLOAD = 1
    ACTION_TYPE_ARCHIVE_CONNECTION = 2
    ACTION_TYPE_WRITE_URL = 3
    ACTION_TYPE_WRITE_FILE = 4

    def __init__(self):
        pass

    @staticmethod
    def create_action_log(outcome, action, message) -> bool:

        is_log_created = False
        log = ""

        if action == Logger.ACTION_TYPE_DOWNLOAD:
            if outcome == Logger.ACTION_SUCCESS:
                log += "Success downloading file : " + message
            else:
                log += "Error downloading file : " + message

            is_log_created = True
        elif action == Logger.ACTION_TYPE_ARCHIVE_CONNECTION:
            if outcome == Logger.ACTION_SUCCESS:
                log += "Success connecting to archive : " + message
            else:
                log += "Error connecting to archive : " + message

            is_log_created = True
        elif action == Logger.ACTION_TYPE_WRITE_URL:
            if outcome == Logger.ACTION_SUCCESS:
                log += "Success writing url to file : " + message
            else:
                log += "Error writing to file : " + message

            is_log_created = True

        if is_log_created:
            Logger._insert_log(Logger.ACTION_TYPE, log)

        return is_log_created

    @staticmethod
    def create_info_log(message):
        pass

    @staticmethod
    def _insert_log(type, log):
        Logger._logs.append(log)

    @staticmethod
    def create_log_file(archive_name, query):
        log_file = ""

        log_file += "Run summary for archive : " + archive_name + "\n"
        log_file += "With query : " + query + "\n"

        for log in Logger._logs:
            log_file += log + "\n"

        return log_file


if __name__ == "__main__":
    output = sys.argv[1]
    output_csv = sys.argv[2]
    output_html = sys.argv[3]
    output_basic_html = sys.argv[4]
    output_error = sys.argv[5]

    inputs = sys.argv[6]

    tool_runner = ToolRunner(inputs,
                             output,
                             output_csv,
                             output_html,
                             output_basic_html,
                             output_error)

    tool_runner.run()
