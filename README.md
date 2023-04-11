# 

[![Galaxy Tool Linting and Tests for push and PR](https://github.com/esg-epfl-apc/tools-astro/actions/workflows/lint-and-test.yml/badge.svg?branch=main)](https://github.com/esg-epfl-apc/tools-astro/actions/workflows/lint-and-test.yml/badge.svg?branch=main)

## may try workflow synthesis

## ODA tool
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
