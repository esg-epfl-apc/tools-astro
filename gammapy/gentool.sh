# conda config --add channels anaconda

# TODO: try to build more recent astropy?

planemo tool_init --force \
                    --id 'gammapy' \
                    --name 'gammapy' \
                    --requirement gammapy@1.0 \
                    --example_command 'gammapy analysis config --filename myconfig.yaml --overwrite' \
                    --example_output myconfig.yaml \
                    --test_case \
                    --cite_url 'https://github.com/gammapy/gammapy' \
                    --help_from_command 'gammapy analysis config --help'
                                        # --example_input WFPC2u5780205r_c0fx.fits \
