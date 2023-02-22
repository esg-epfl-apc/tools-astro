# Comment from Bjorn on Data Sources

As discussed today, here are a few "data_source" links. That might be relevant.

As recap, data can be imported in various ways into Galaxy.

* Download via URL / Paste Content / Upload from computer
* You can always create a tool that does in the background somekind of wget/curl to get the data
* if you have a pure API endpoint this is probably your way to go.
* If data is exposed with a standard protocoll, SFTP, S3, WebDAV ... if your data store can be queried via something like pyfilesystem, Galaxy has a nice integration for that
    Via a so called "data source"-tool a special Galaxy tool 
    
A data source tool can only be used if you have a website that is exposing the data somehow via a GUI.
You need code at website where you want to get data from. So a collaboration is needed.

What happens. You search in Galaxy for your data-source-tool and click on it.
You will be redirected to the data-store website and Galaxy will pass a few parameters to the external website.
The user is now on this external website, selects and filters datasets. The website will show
A button "Send data to Galaxy" this button can be only shown when coming from a Galaxy server, since Galaxy passes along some parameters.

If the user clicks this button, the user is redirected to Galaxy and Galaxy runs a job to get the data from this external website.

A few information here:
https://galaxyproject.org/admin/internals/data-sources/

The best resource to try it out is problably this repo:
https://github.com/hexylena/galaxy-data_source-examples

This repo includes a toy server (the external server that exposes data) and a
very simple Galaxy data-source-tool.

I hope that can get you started.

