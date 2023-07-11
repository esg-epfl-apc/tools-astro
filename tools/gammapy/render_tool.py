import re

tool_xml = open("gammapysim.xml").read()

tool_code = open("gammapysim/gammapysim.py").read()
tool_code = re.sub('\$', '$$', tool_code, re.M)

with open("gammapysim.xml", "w") as f:
    f.write(re.sub(
            r'<configfile name="gammapysim.py">.*?</configfile>', 
            '<configfile name="gammapysim.py"><![CDATA[\n' + tool_code + '\n]]><', 
            tool_xml, 
            re.M | re.S
        )
    )