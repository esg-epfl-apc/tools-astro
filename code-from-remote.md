

    Another related question, all what is needed for the tool to run should be either inside the tool definition or in requirements, right? E.g. the code is in a gitlab repository (or maybe just stored on zenodo), how do I make it part of the tool?

The best practice recommendation is to package your code at first. E.g. with the conda package manager. Then you simply define the dependency as a requirement of the tool.

    patrick-austin
    There's required files (https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-required-files). We're using it for ~100 line Python scripts for stuff like plotting and more complex logic for generating input files. Not sure how suitable it is for larger amounts of code though.

Those "required files" is an annotation feature. The idea here is, iff you ship a script next to your tool.xml then you annotate it, so Galaxy can better copy around those files. However, if you can package all the code in a package and "depend" on it.
Hope that makes sense.
Volodymyr

    bgruening

        Another related question, all what is needed for the tool to run should be either inside the tool definition or in requirements, right? E.g. the code is in a gitlab repository (or maybe just stored on zenodo), how do I make it part of the tool?

    The best practice recommendation is to package your code at first. E.g. with the conda package manager. Then you simply define the dependency as a requirement of the tool.

I see. It will end up with moving lot's of code to conda, and potentially frequently. But we can do it, need to try.
