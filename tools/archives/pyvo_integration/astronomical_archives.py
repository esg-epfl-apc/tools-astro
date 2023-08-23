import sys
import os

import pyvo
from pyvo import registry
from pyvo import DALQueryError, DALServiceError, DALAccessError

import urllib
from urllib import parse, request


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

    def __init__(self, id, title, name, access_url):
        self.id = id,
        self.title = title,
        self.name = name,
        self.access_url = access_url
        self.initialized = False
        self.archive_service = None
        self.tables = None

    def get_resources(self, query, number_of_results, url_field='access_url'):
        resource_list = []
        resource_list_hydrated = []

        error_message = None

        if self.initialized:

            try:
                raw_resource_list = self.archive_service.search(query)

                for i, resource in enumerate(raw_resource_list):
                    if i < number_of_results:
                        # resource_list.append(resource[url_field])
                        resource_list_hydrated.append(self._get_resource_object(resource))
            except DALQueryError as dqe:
                if self.has_obscore_table():
                    error_message = "Error in query -> " + query
                    Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
                else:
                    error_message = "No obscore table in the archive"
                    Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
            except DALServiceError as dse:
                error_message = "Error communicating with the service"
                Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
            except Exception as e:
                error_message = "Unknow error while querying the service"
                Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)

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
        except DALAccessError as dae:
            error_message = "A connection to the service could not be established"
            Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_ARCHIVE_CONNECTION, error_message)
        except Exception as e:
            error_message = "Unknow error while initializing TAP service"
            Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_ARCHIVE_CONNECTION, error_message)

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

        if not next((item for item in self.tables if item["name"] == table_name), False):
            is_valid = False

        return is_valid

    def has_obscore_table(self) -> bool:
        has_obscore_table = False

        has_obscore_table = self._has_table("ivoa.obscore")

        return has_obscore_table

    def _has_table(self, table_name) -> bool:
        _has_table = False

        _has_table = next((item for item in self.tables if item["name"] == table_name), False)

        return _has_table

    def get_archive_name(self, archive_type):
        name = ''

        try:
            if archive_type == 'registry':
                name = str(self.title).strip("',()")
            else:
                name = self.access_url
        except Exception as e:
            name = 'Unknown archive title'

        return name


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
            parameters['waveband'] = Waveband.wavebands[self.waveband]

        if Service.is_service_supported(self.service_type):
            parameters['service_type'] = Service.services[self.service_type]
        else:
            parameters['service_type'] = Service.services['TAP']

        return parameters


class Registry:

    def __init__(self):
        pass

    @staticmethod
    def search_registries(rsp: RegistrySearchParameters, number_of_registries):

        parameters = rsp.get_parameters()

        keywords = parameters['keywords']
        waveband = parameters['waveband']
        service_type = parameters['service_type']

        registry_list = []

        if not waveband:
            registry_list = registry.search(keywords=keywords, servicetype=service_type)
        else:
            registry_list = registry.search(keywords=keywords, waveband=waveband, servicetype=service_type)

        if registry_list:
            registry_list = Registry._get_registries_from_list(registry_list, 1)

        return registry_list

    @staticmethod
    def _get_registries_from_list(registry_list, number_of_registries):

        archive_list = []

        for i, registry in enumerate(registry_list):
            if i < number_of_registries:
                archive = TapArchive(registry.standard_id, registry.res_title, registry.short_name, registry.access_url)
                archive_list.append(archive)

        return archive_list


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
            statement = ''

            if value != '':
                statement = str(key) + ' = ' + '\'' + str(value) + '\' '

                if is_first_statement:
                    is_first_statement = False
                    where_clause += 'WHERE '
                else:
                    statement = 'AND ' + statement

                where_clause += statement

        return where_clause


class ADQLObscoreQuery(BaseADQLQuery):
    order_by_field = {
        'size': 'access_estsize',
        'collection': 'obs_collection',
        'object': 'target_name'
    }

    base_query = 'SELECT TOP 100 * FROM ivoa.obscore '

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
                 order_by):

        super().__init__()

        if calibration_level == 'none':
            calibration_level = ''

        if order_by == 'none':
            order_by = ''

        if dataproduct_type == 'none' or dataproduct_type is None:
            dataproduct_type = ''

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
            'calibration_level': calibration_level,
            't_min': t_min,
            't_max': t_max
        }

        self.order_by = order_by

    def get_query(self):
        return ADQLObscoreQuery.base_query + self.get_where_statement() + self.get_order_by_statement()

    def get_order_by_statement(self):
        if self.order_by != '':
            return self._get_order_by_clause(self.order_by)
        else:
            return ''

    def _get_order_by_clause(self, order_type):

        obscore_order_type = ADQLObscoreQuery.order_by_field[order_type]

        return super()._get_order_by_clause(obscore_order_type)

    def get_where_statement(self):
        return self._get_where_clause(self.parameters)

    def _get_where_clause(self, parameters):
        return super()._get_where_clause(parameters)


class ADQLTapQuery(BaseADQLQuery):
    base_query = 'SELECT TOP 100 * FROM '

    def __init__(self):
        super().__init__()

    def get_order_by_clause(self, order_type):
        return super().get_order_by_clause(order_type)

    def get_query(self, table, where_field, where_condition):
        if where_field != '' and where_condition != '':
            return ADQLTapQuery.base_query + str(table) + ' WHERE ' + str(where_field) + ' = ' + '\'' + str(
                where_condition) + '\''
        else:
            return ADQLTapQuery.base_query + str(table)


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
    def generateBasicUrlFile(urls_data):
        FileHandler.write_urls_to_output(urls_data)

    @staticmethod
    def generateCSVOutput(urls_data):
        csv_file = ''

        for i, url in enumerate(urls_data):
            csv_file += url['access_url']
            csv_file += ','
            csv_file += url['est_size']
            csv_file += ','
            csv_file += url['obs_collection']
            csv_file += ','
            csv_file += url['target_name']
            csv_file += '\n'

        return csv_file

    @staticmethod
    def generateHTMLOutput(urls_data, archive_name, adql_query):
        html_file = ''

        css = """ <head><style>

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

            .five h2 {
              text-align: center;
              font-size: 22px;
              font-weight: 700; color:#202020;
              text-transform: uppercase;
              word-spacing: 1px; letter-spacing:2px;
              margin-bottom: 50px;
            }

            .five h2 span {
              padding-top: 40px;
              text-transform: none;
              font-size:.80em;
              font-weight: bold;
              font-family: "Playfair Display","Bookman",serif;
              color:#999; letter-spacing:-0.005em; word-spacing:1px;
              letter-spacing:none;
            }

            .five h1:before {
              background-color: #dfdfdf;
            }

        </style></head>"""

        html_file += css

        html_file += """<div class="five">
                      <h2>Resources Preview from archive :
                        <br />
                        <br />
                        <span> """ + str(archive_name) + """</span>
                        <br />
                        <br />
                        <span> With ADQL query : """ + str(adql_query) + """</span>
                      </h2>
                    </div>"""

        for resource in urls_data:

            html_file += '<table class="fl-table"><thead><tr>'

            for key in resource.keys():
                html_file += '<th>' + str(key) + '</th>'

            html_file += '</thead></tr>'
            html_file += '<tbody><tr>'

            for key, value in resource.items():
                html_file += '<td>' + str(value) + '</td>'

            html_file += '</tr></tbody>'
            html_file += '</table>'

            if 'preview' in resource:
                html_file += '<details><summary>Preview</summary><img src="' + str(
                    resource['preview']) + '"/></details>'

            if 'preview_url' in resource:
                html_file += '<details><summary>Preview</summary><img src="' + str(
                    resource['preview_url']) + '"/></details>'

            if 'postcard_url' in resource:
                html_file += '<details><summary>Preview</summary><img src="' + str(
                    resource['postcard_url']) + '"/></details>'

            html_file += '<br />'
            html_file += '<br />'

        return html_file

    @staticmethod
    def generateBasicHTMLOutput(urls_data, archive_name, adql_query):
        html_file = ''

        html_file += """<div>
                          <h2>Resources Preview from archive :
                            <br />
                            <br />
                            <span> """ + str(archive_name) + """</span>
                            <br />
                            <br />
                            <span> With ADQL query : """ + str(adql_query) + """</span>
                          </h2>
                        </div>"""

        for resource in urls_data:

            html_file += '<table border="1"><thead><tr>'

            for key in resource.keys():
                html_file += '<th>' + str(key) + '</th>'

            html_file += '</thead></tr>'
            html_file += '<tbody><tr>'

            for key, value in resource.items():
                html_file += '<td>' + str(value) + '</td>'

            html_file += '</tr></tbody>'
            html_file += '</table>'

            if 'preview' in resource:
                html_file += '<details><summary>Preview</summary><img src="' + str(
                    resource['preview']) + '"/></details>'

            if 'preview_url' in resource:
                html_file += '<details><summary>Preview</summary><img src="' + str(
                    resource['preview_url']) + '"/></details>'

            if 'postcard_url' in resource:
                html_file += '<details><summary>Preview</summary><img src="' + str(
                    resource['postcard_url']) + '"/></details>'

            html_file += '<br />'
            html_file += '<br />'

        return html_file


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
                file_output.write(url[access_url] + ',')

    @staticmethod
    def write_multiple_outputs(output_id):

        out_files = {}

        for i in range(1, 10):
            out_files[i] = open(
                os.path.join(database_tmp_dir, "primary_%s_%s_visible_interval_%s" % (output_id, i, i)), "w+"
            )
            out_files[i].write("aaaaa")

        for file_out in out_files.values():
            file_out.close()

    @staticmethod
    def write_collection(output):
        dir = os.getcwd()

        dir += '/fits'

        upload_dir = os.path.join(dir, 'aaaaa.fits')

        with open(output, "w") as file_output:
            file_output.write(upload_dir)

    @staticmethod
    def write_collection1(index):
        dir = os.getcwd()

        dir += '/fits'

        upload_dir = os.path.join(dir, index + '.fits')

        with open(upload_dir, "w") as file_output:
            file_output.write(upload_dir)

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

        file_name = 'archive file '

        try:

            if (url_parts[-1]) != '':
                file_name = url_parts[-1]
            elif len(url_parts) > 1:
                file_name = url_parts[-2]
        except Exception:
            pass

        return file_name


class Logger:
    _logs = []

    ACTION_SUCCESS = 1
    ACTION_ERROR = 2

    ACTION_TYPE = 1
    INFO_TYPE = 2

    ACTION_TYPE_DOWNLOAD = 1
    ACTION_TYPE_ARCHIVE_CONNECTION = 2

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
    number_of_files = sys.argv[6]
    archive_type = sys.argv[7]

    archive = ''
    archive_name = 'No archive name'

    file_url = []
    error_message = None

    outputHandler = OutputHandler()

    if number_of_files is not None and number_of_files != '':
        if int(number_of_files) < 1:
            number_of_files = 1
        elif int(number_of_files) > 100:
            number_of_files = 100
    else:
        number_of_files = 1

    if archive_type == 'registry':

        keyword = sys.argv[8]
        waveband = sys.argv[9]
        service_type = sys.argv[10]
        query_type = sys.argv[11]

        if query_type == 'obscore_query':

            dataproduct_type = sys.argv[12]
            obs_collection = sys.argv[13]
            facility_name = sys.argv[14]
            instrument_name = sys.argv[15]
            em_min = sys.argv[16]
            em_max = sys.argv[17]
            target_name = sys.argv[18]
            obs_publisher_id = sys.argv[19]
            s_fov = sys.argv[20]
            calibration_level = sys.argv[21]
            order_by = sys.argv[22]
            obs_title = sys.argv[23]
            obs_id = sys.argv[24]
            t_min = sys.argv[25]
            t_max = sys.argv[26]

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
                                                    order_by)

            adql_query = obscore_query_object.get_query()

        elif query_type == 'raw_query':
            tap_table = sys.argv[12]
            where_field = sys.argv[13]
            where_condition = sys.argv[14]
            url_field = sys.argv[15]

            adql_query = ADQLTapQuery().get_query(tap_table, where_field, where_condition)

        else:
            adql_query = ADQLObscoreQuery.base_query

        rsp = RegistrySearchParameters(keyword=keyword, waveband=waveband, service_type=service_type)

        archive_list = Registry.search_registries(rsp, 1)

        if archive_list:
            archive = archive_list[0]

            archive_name = archive.get_archive_name(archive_type)

            is_initialisation_success, error_message = archive.initialize()

            if is_initialisation_success:
                if query_type == 'raw_query':
                    if url_field:
                        file_url, error_message = archive.get_resources(adql_query, int(number_of_files), url_field)
                    else:
                        error_message = "no url field specified"
                        Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
                else:
                    file_url, error_message = archive.get_resources(adql_query, int(number_of_files))

        else:
            error_message = "no archive matching search parameters"
            Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_ARCHIVE_CONNECTION, error_message)

    elif archive_type == 'archive':

        service_url = sys.argv[8]

        query_type = sys.argv[9]

        if query_type == 'obscore_query':

            dataproduct_type = sys.argv[10]
            obs_collection = sys.argv[11]
            facility_name = sys.argv[12]
            instrument_name = sys.argv[13]
            em_min = sys.argv[14]
            em_max = sys.argv[15]
            target_name = sys.argv[16]
            obs_publisher_id = sys.argv[17]
            s_fov = sys.argv[18]
            calibration_level = sys.argv[19]
            order_by = sys.argv[20]
            obs_title = sys.argv[21]
            obs_id = sys.argv[22]
            t_min = sys.argv[23]
            t_max = sys.argv[24]

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
                                                    order_by)

            adql_query = obscore_query_object.get_query()

        elif query_type == 'raw_query':

            tap_table = sys.argv[10]
            where_field = sys.argv[11]
            where_condition = sys.argv[12]
            url_field = sys.argv[13]

            adql_query = ADQLTapQuery().get_query(tap_table, where_field, where_condition)

        else:
            adql_query = ADQLObscoreQuery.base_query

        archive = TapArchive(1, 'name', 'title', service_url)

        archive_name = archive.get_archive_name(archive_type)

        is_initialisation_success, error_message = archive.initialize()

        if is_initialisation_success:
            if query_type == 'raw_query':
                if url_field:
                    file_url, error_message = archive.get_resources(adql_query, int(number_of_files), url_field)
                else:
                    error_message = "no url field specified"
                    Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
            else:
                file_url, error_message = archive.get_resources(adql_query, int(number_of_files))

    if file_url and output_csv != 'XXXX':

        if query_type == 'raw_query':
            if file_url:
                FileHandler.write_urls_to_output(file_url, output_csv, url_field)
            else:
                error_message = "no url field specified"
                Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
        else:
            FileHandler.write_urls_to_output(file_url, output_csv)

    if file_url and output != 'XXXX':

        if query_type == 'raw_query':
            access_url = url_field
        else:
            access_url = 'access_url'

        if access_url:

            try:
                fits_file = FileHandler.download_file_from_url(file_url[0][access_url])
                FileHandler.write_file_to_output(fits_file, output, "wb")

                log_message = "from url " + file_url[0][access_url]
                Logger.create_action_log(Logger.ACTION_SUCCESS, Logger.ACTION_TYPE_DOWNLOAD, log_message)
            except Exception as e:
                error_message = "from url " + file_url[0][access_url]
                Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)

            for i, url in enumerate(file_url[1:], start=1):
                try:
                    fits_file = FileHandler.download_file_from_url(url[access_url])
                    FileHandler.write_file_to_subdir(fits_file, FileHandler.get_file_name_from_url(url[access_url]))

                    log_message = "from url " + url[access_url]
                    Logger.create_action_log(Logger.ACTION_SUCCESS, Logger.ACTION_TYPE_DOWNLOAD, log_message)
                except Exception as e:
                    error_message = "from url " + url[access_url]
                    Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)
        else:
            error_message = "no url field specified"
            Logger.create_action_log(Logger.ACTION_ERROR, Logger.ACTION_TYPE_DOWNLOAD, error_message)

    if file_url and (output_html != 'XXXX' or output_basic_html != 'XXXX'):

        if output_html:
            html_file = OutputHandler.generateHTMLOutput(file_url, archive_name, adql_query)
            FileHandler.write_file_to_output(html_file, output_html)

        if output_basic_html:
            html_file = OutputHandler.generateBasicHTMLOutput(file_url, archive_name, adql_query)
            FileHandler.write_file_to_output(html_file, output_basic_html)

    if file_url is None or error_message:

        error_file = Logger.create_log_file(archive_name, adql_query)

        if error_message is None:
            error_file += "\n No resources matching parameters found"

        FileHandler.write_file_to_output(error_file, output_error)

    else:

        error_file = Logger.create_log_file(archive_name, adql_query)

        error_file += "\n Tool run executed with success"

        FileHandler.write_file_to_output(error_file, output_error)