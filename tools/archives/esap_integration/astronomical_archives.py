import sys
import json
import re

from urllib import request, parse


class ArchiveRequest:

    def __init__(self, base_url, request_parameters):
        self._base_url = base_url
        self._request_parameters = dict(request_parameters)
        self._request = ''

    def create_request(self) -> str:

        for key in list(self._request_parameters.keys()):
            if self._request_parameters[key] == '':
                del self._request_parameters[key]

        querystring = parse.urlencode(self._request_parameters)

        request_url = self._base_url + '?' + querystring

        self._request = request_url

        return request_url


class Archive:

    def __init__(self):
        self._name = "archive"
        self.base_url = ""

    def query_archive1(self, base_url, url_params) -> None:
        archive_querystring = parse.urlencode(url_params)

        with request.urlopen(base_url + '?' + archive_querystring) as response:
            file = response.read()

    def query_archive_request(self, request) -> None:

        with request.urlopen(request) as response:
            file = response.read()


class AstronArchive(Archive):

    name = "Astron Archive"
    archive_type_name = "astron"

    archive_query_parameters = {
        'level': 'processed',
        'target': '',
        'ra': '',
        'dec': '',
        'fov': '',
        'collection': '',
        'page': '',
        'archive_uri': 'astron_vo'
    }

    def __init__(self):
        super().__init__()
        self._base_url = "https://sdc-dev.astron.nl/esap-api/query/query/"

    def query_archive(self, target, ra, dec, fov, collection, number_of_results) -> []:

        query_results = []

        self.archive_query_parameters['target'] = target
        self.archive_query_parameters['ra'] = ra
        self.archive_query_parameters['dec'] = dec
        self.archive_query_parameters['fov'] = fov
        self.archive_query_parameters['collection'] = collection

        astron_archive_query_url = self.__create_archive_query()

        with request.urlopen(astron_archive_query_url) as json_response:
            raw_astron_archive_response = json_response.read()

        astron_response = AstronArchiveResponse(raw_astron_archive_response)

        try:
            if astron_response.is_valid():
                query_results = astron_response.get_response_results(number_of_results)
        except InvalidAstronArchiveResponse:
            pass
        except NoResultsInAstronArchiveResponse:
            pass

        return query_results

    def __create_archive_query(self):
        astron_archive_request = ArchiveRequest(self._base_url, self.archive_query_parameters)
        query = astron_archive_request.create_request()

        return query


class AstronArchiveResponse:

    _allowed_dataproduct_type = ['image']

    def __init__(self, raw_response):
        self._raw_response = raw_response
        self._json_response = json.loads(raw_response.decode('utf-8'))

    def is_valid(self) -> bool:
        if "count" not in self._json_response or "results" not in self._json_response:
            raise InvalidAstronArchiveResponse
        elif self._json_response['count'] == 0:
            raise NoResultsInAstronArchiveResponse

        return True

    def get_response_results(self, number_of_results):

        results = []

        for i, ressource in enumerate(self._json_response['results']):
            if i < number_of_results and ressource['dataproduct_type'] in self._allowed_dataproduct_type:
                results.append(ressource['url'])
            else:
                break

        return results


class IVOAArchiveResponse:

    def __init__(self, raw_response):
        self._raw_response = raw_response
        self._json_response = json.loads(raw_response.decode('utf-8'))

    def is_valid(self) -> bool:
        if "count" not in self._json_response or "results" not in self._json_response:
            raise InvalidAstronArchiveResponse
        elif self._json_response['count'] == 0:
            raise NoResultsInAstronArchiveResponse

        return True

    def get_response_registries_results(self, number_of_results):

        results = []

        for i, ressource in enumerate(self._json_response['results']):
            if i < number_of_results:
                results.append(ressource['access_url'])
            else:
                break

        return results

    def get_response_ressource_results(self, number_of_results):

        results = []

        for i, ressource in enumerate(self._json_response['results']):
            if i < number_of_results:
                results.append(ressource['url'])
            else:
                break

        return results


class IVOAArchive(Archive):

    name = "IVOA Registry Archive"

    archive_type_name = "ivoa"

    _number_of_registry_to_query = 1

    _allowed_service_types = ['tap', 'scs', 'ssa', 'sia']

    registries_query_parameters = {
        'keyword': 'hess',
        'waveband': '',
        'service_type': 'tap',
        'dataset_uri': 'vo_reg'
    }

    registry_query_parameters = {
        'adql_query': "SELECT TOP 1 * from ivoa.obscore where dataproduct_type = 'image'",
        'access_url': '',
        'service_type': 'tap',
        'dataset_uri': 'vo_reg'
    }

    def __init__(self):
        super().__init__()
        self._base_registries_list_url = "https://sdc-dev.astron.nl/esap-api/query/get-services/"
        self._base_registry_ressource_list_url = "https://sdc-dev.astron.nl/esap-api/query/query/"

    def query_archive(self, keyword, waveband, service_type, adql_query, number_of_results) -> []:

        query_results = []

        self.registries_query_parameters['keyword'] = keyword
        self.registries_query_parameters['waveband'] = waveband
        self.registries_query_parameters['service_type'] = service_type
        self.registry_query_parameters['service_type'] = service_type
        self.registry_query_parameters['adql_query'] = adql_query

        self.check_parameters()

        ivoa_registries_query_url = self.__create_registries_query()

        with request.urlopen(ivoa_registries_query_url) as json_response:
            raw_ivoa_registries_response = json_response.read()

        ivoa_response = IVOAArchiveResponse(raw_ivoa_registries_response)

        registries_results = []

        try:
            if ivoa_response.is_valid():
                registries_results = ivoa_response.get_response_registries_results(self._number_of_registry_to_query)
        except InvalidAstronArchiveResponse:
            pass
        except NoResultsInAstronArchiveResponse:
            pass

        self.registry_query_parameters['access_url'] = registries_results[0]

        ivoa_ressource_query_url = self.__create_ressources_query()

        with request.urlopen(ivoa_ressource_query_url) as json_response:
            raw_ivoa_ressources_response = json_response.read()

        ivoa_response = IVOAArchiveResponse(raw_ivoa_ressources_response)

        try:
            if ivoa_response.is_valid():
                ressources_results = ivoa_response.get_response_ressource_results(number_of_results)
                query_results = ressources_results
        except InvalidIVOAArchiveResponse:
            pass
        except NoResultsInIVOAArchiveResponse:
            pass

        return query_results

    def __create_registries_query(self):
        ivoa_archive_request = ArchiveRequest(self._base_registries_list_url, self.registries_query_parameters)
        query = ivoa_archive_request.create_request()

        return query

    def __create_ressources_query(self):
        ivoa_archive_request = ArchiveRequest(self._base_registry_ressource_list_url, self.registry_query_parameters)
        query = ivoa_archive_request.create_request()

        return query

    def check_parameters(self):

        if self.registries_query_parameters['waveband'] == 'all':
            self.registries_query_parameters['waveband'] = ''

        if self.registries_query_parameters['service_type'] not in self._allowed_service_types:
            self.registries_query_parameters['service_type'] = 'tap'
            self.registry_query_parameters['service_type'] = 'tap'

        self.check_adql_query()

    def check_adql_query(self):
        if self.registry_query_parameters['adql_query'] == '':
            self.registry_query_parameters['adql_query'] = "SELECT TOP 1 * from ivoa.obscore where dataproduct_type = 'image'"


class ADQLQueryValidator:

    def __init__(self, query):
        self.query = query

    def is_valid(self):
        pattern = re.compile("^(?=.*SELECT.*FROM)(?!.*(?:WHERE|GROUP BY|LIMIT|ORDER BY)).*$")
        return pattern.match(self.query), pattern


class InvalidAstronArchiveResponse(Exception):

    def __init__(self, query, message="Invalid response from astron archive"):
        self.message = message
        super().__init__(self.message)


class NoResultsInAstronArchiveResponse(Exception):

    def __init__(self, query, message="No results from astron archive"):
        self.message = message
        super().__init__(self.message)


class InvalidIVOAArchiveResponse(Exception):

    def __init__(self, query, message="Invalid response from ivoa archive"):
        self.message = message
        super().__init__(self.message)


class NoResultsInIVOAArchiveResponse(Exception):

    def __init__(self, query, message="No results from ivoa archive"):
        self.message = message
        super().__init__(self.message)


class InvalidADQLQuery(Exception):

    def __init__(self, query, message="Syntax Error in the AQDL query : "):
        self.query = query
        self.message = message
        super().__init__(self.message + self.query)


class NoMatchingResultsFound(Exception):

    def __init__(self, query, message="No results found for query : "):
        self.query = query
        self.message = message
        super().__init__(self.message + self.query)


class QueryError(Exception):

    def __init__(self, query, message="Error while running query : "):
        self.query = query
        self.message = message
        super().__init__(self.message + self.query)


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
    archive_type = sys.argv[2]

    file_url = []

    if archive_type == AstronArchive.archive_type_name:

        target = sys.argv[3]
        ra = sys.argv[4]
        dec = sys.argv[5]
        search_radius = sys.argv[6]
        astron_collection = sys.argv[7]

        file_type = sys.argv[8]
        number_of_files = sys.argv[9]
        tabular_output = sys.argv[10]

        astron_archive = AstronArchive()

        file_url = astron_archive.query_archive(target, ra, dec, search_radius, astron_collection, 1)

    elif archive_type == IVOAArchive.archive_type_name:

        keyword = sys.argv[3]
        waveband = sys.argv[4]
        service_type = sys.argv[5]
        adql_query = sys.argv[6]

        file_type = sys.argv[7]
        number_of_files = sys.argv[8]
        tabular_output = sys.argv[9]

        ivoa_archive = IVOAArchive()

        file_url = ivoa_archive.query_archive(keyword, waveband, service_type, adql_query, 1)

    if file_url:
        fits_file = FileHandler.download_file_from_url(file_url[0])
        FileHandler.write_file_to_output(fits_file, output)

        # if file_type == 'urls':
        #     url_file = FileHandler.create_file_from_url_list(file_url)
        #     FileHandler.write_file_to_output(url_file, tabular_output)
