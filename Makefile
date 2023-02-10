init:
	planemo tool_init --force \
                    --id 'odarequest' \
                    --name 'request ODA tool' \
                    --requirement oda-api@1.1.33 \
                    --example_command 'oda-api get -i spi_acs -p spi_acs_lc -a time_bin=`cat inp` > odaapi.stdout 2>&1' \
                    --example_input inp \
                    --example_output odaapi.stdout \
					--test_case \
					--cite_url 'https://github.com/oda-hub/oda_api' \
					--help_from_command 'oda-api'


lint:
	planemo l

test:
	planemo t

serve:
	planemo s