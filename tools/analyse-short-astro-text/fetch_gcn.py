import requests
import re
gcn_number = 40274


def fetch_gcn(gcn_number):
    """
    Fetches the GCN page for the given GCN number and returns the GCN text.
    
    input : gcn_number (int): The GCN number to fetch.
    output : response_text (str): The content of the GCN body.
    If an error occurs, it returns None.
    """

    # URL of the GCN page
    url = 'https://gcn.nasa.gov/circulars/{}.json'.format(gcn_number)

    response = requests.get(url)

    # Parse JSON
    data = response.json()

    # Get the "body" field
    body = data.get("body")
    
    # cleans non-ASCII characters
    cleaned_text = re.sub(r'[^\x00-\x7F]+', '', body)  # remove non-ASCII
    print(cleaned_text)

    return cleaned_text


if __name__ == "__main__":
    fetch_gcn(gcn_number)

