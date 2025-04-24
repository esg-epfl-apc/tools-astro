import requests
from bs4 import BeautifulSoup
import os
import re
atel_number = 16672


def fetch_atel(atel_number):
    """
    Fetches the ATel page for the given ATel number and returns the AteL text.
    It assumes that the paragraph is the first one after the paragraph that 
    contains the string "Tweet".
    input : atel_number (int): The ATel number to fetch.
    output : response_text (str): The HTML content of the ATel text.
    If an error occurs, it returns None.
    """
    
    # URL of the ATel page
    url = 'https://www.astronomerstelegram.org/?read={}'.format(atel_number)
    
    # To fake the User-Agent header
    # This is to avoid being blocked by the server for not having a User-Agent
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36'
    }

    # This is mainly for testing purposes
    # Check if the file already exists
    # If it does, read the content from the file
    # If it doesn't, fetch the page and save it to a file
    # The file name is based on the ATel number
    # For example, if the ATel number is 16672, the file name will be 'atel_16672.html'

    fname = 'atel_{}.html'.format(atel_number)

    if not os.path.isfile(fname):
        # Send a GET request to the URL
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            print("Page fetched successfully.")    
            with open(fname, 'w', encoding='utf-8') as f:
                f.write(response.text)
            response_text = response.text
        else:
            print(f"Failed to retrieve the page. Status code: {response.status_code}")
            return None
    elif os.path.isfile(fname):
        print("Page already fetched.")
        with open(fname, 'r', encoding='utf-8') as f:
            response_text = f.read()
    else:
        print("Page not found.")
        return None

    soup = BeautifulSoup(response_text, 'html.parser')

    # print(soup.prettify())

    tds = soup.body.find_all("p")
    twitter_index = -1
    for i, td in enumerate(tds):
        if 'Tweet' in td.get_text(strip=True):
            twitter_index = i

    para = tds[twitter_index + 1]

    cleaned_text = re.sub(r'[^\x00-\x7F]+', '', para.text)  # remove non-ASCII
    print(cleaned_text)
    
    return cleaned_text


if __name__ == "__main__":
    fetch_atel(atel_number)