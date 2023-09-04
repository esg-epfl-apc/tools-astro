# 

[![Galaxy Tool Linting and Tests for push and PR](https://github.com/esg-epfl-apc/tools-astro/actions/workflows/lint-and-test.yml/badge.svg?branch=main)](https://github.com/esg-epfl-apc/tools-astro/actions/workflows/lint-and-test.yml/badge.svg?branch=main)

[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](h[ttps://usegalaxy.eu/root?tool_id=interactive_tool_qgis](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fastroteam%2Fastropy_fits2bitmap%2Fastropy_fits2bitmap%2F0.1.0%2Bgalaxy0&version=0.1.0%20galaxy0))


## Astropy tools

[astropy](https://github.com/astropy/astropy) is a python package comprising a collection of libraries serving as a basis for much of the modern research in astrophysics.

Here we add some of the basic operations implemented in **astropy** as galaxy tools:

* *fitsinfo*: dumps information about fits file.
* *fits2bitmap*: converts fits image to a png which can be viewed in galaxy directly.

These tools provide scientifically meaningful outputs, but serve primarily as proof-of-concept for introducing astrophysical tools into galaxy.

## TODO: GammaPy simulation

This tool collection allows to produce synthetic gamma-ray observations with Cherenkov Telescope Array data. 
The tools rely on [GammaPy](https://github.com/gammapy/gammapy) and uses relatively little input data to yield scientifically useful outputs suitable for building further data analysis pipelines designed to assess telescope performance.

## TODO: ODA tool
* TODO: tool for oda api request
* QUESTION: how to explain queries to archives, i.e. implicit state dependencies?
* QUESTION: tool for each product or common tool? templates?

we can either generate a tool for each product, or 

uses type interpretation to match to oda-api python classes types

## ODA-annotated nb
* QUESTION: location of the nb file, git repo?
* TODO:
  tool for annotated nb
  annotates semantically inputs (tag:parameters) and outputs (tag:outputs)
  uses type interpretation to match to oda-api python classes types
  

TODO: remote data, gammapy responses, astro archives
TODO: complex operations to be expressed as tools? publish package every time?
TODO: understand how to deal with containers
TODO: add gammapy
TODO: supporting generic notebooks


https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/get-data/slides.html#4

## Big data issue

when we compute on data, compute should have access to it
*deferred data* approach, oakridge group

separate tool class for fetching data?

*filesources* where we can search for folders

*datasources* are directories

puslar arc is doing job staging and data staging


https://www.pyfilesystem.org/

may need to have admin training features on datasources

join WP4 for big data case, see chat
