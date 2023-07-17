import sys

import pyvo
from pyvo import registry

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

    def get_resources(self, query, number_of_results):
        resource_list = []

        if self.initialized:

            try:
                raw_resource_list = self.archive_service.search(query)

                for i, resource in enumerate(raw_resource_list):
                    if i < number_of_results:
                        resource_list.append(resource['access_url'])
            except Exception as e:
                pass

        return resource_list

    def initialize(self):
        self._get_service()

        if self.archive_service:
            self._set_archive_tables()
            self.initialized = True

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

        if table_name not in self.tables.values():
            is_valid = False

        return is_valid


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


class FileHandler:

    def __init__(self):
        pass

    @staticmethod
    def download_file_from_url(file_url):
        with request.urlopen(file_url) as response:
            fits_file = response.read()

        return fits_file

    @staticmethod
    def create_file_from_url_list(file_urls):
        url_file = ''

        for url in file_urls:
            url_file + url + ' '

        return url_file

    @staticmethod
    def write_file_to_output(file, output):
        with open(output, "wb") as file_output:
            file_output.write(file)


if __name__ == "__main__":

    output = sys.argv[1]

    keyword = sys.argv[2]
    waveband = sys.argv[3]
    service_type = sys.argv[4]
    adql_query = sys.argv[5]

    file_type = sys.argv[6]
    number_of_files = sys.argv[7]
    tabular_output = sys.argv[8]

    file_url = []

    rsp = RegistrySearchParameters(keyword=keyword, waveband=waveband, service_type='TAP')

    archive_list = Registry.search_registries(rsp, 1)

    if archive_list:
        archive_list[0].initialize()

        tap_query = TapQuery(adql_query)

        query = tap_query.get_query()

        file_url = archive_list[0].get_resources(query, 1)

        if file_url:
            fits_file = FileHandler.download_file_from_url(file_url[0])
            FileHandler.write_file_to_output(fits_file, output)

            # if file_type == 'urls':
            #     url_file = FileHandler.create_file_from_url_list(file_url)
            #     FileHandler.write_file_to_output(url_file, tabular_output)
