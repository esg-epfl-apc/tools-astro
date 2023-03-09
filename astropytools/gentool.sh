wget -c https://fits.gsfc.nasa.gov/samples/WFPC2u5780205r_c0fx.fits 
# conda config --add channels anaconda

# TODO: try to build more recent astropy?

planemo tool_init --force \
                    --id 'fitsinfo' \
                    --name 'Convert to FASTA (seqtk)' \
                    --requirement astropy@5.1 \
                    --example_command 'fitsinfo WFPC2u5780205r_c0fx.fits > fitsinfo.out' \
                    --example_input WFPC2u5780205r_c0fx.fits \
                    --example_output fitsinfo.out \
                    --test_case \
                    --cite_url 'https://github.com/astropy/astropy' \
                    --help_from_command 'fitsinfo --help'