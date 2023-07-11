import re

tool_xml = open("gammapysim.xml").read()

with open("gammapysim.xml", "w") as f:
    f.write(re.sub(
            r'<configfile name="gammapysim.py">.*?<', 
            '<configfile name="gammapysim.py"><![CDATA[\n' + open("gammapysim/gammapysim.py").read() + '\n]]><', 
            tool_xml, 
            re.M | re.S
        )
    )