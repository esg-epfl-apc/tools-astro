import re

tool_xml = open("gammapysim.xml").read()

tool_code = open("gammapysim/gammapysim.py").read()
tool_code = re.sub('\$', '\\$', tool_code, re.M)

tool_xml_patched = re.sub(
        r'configfile name="gammapysim_py"', 
        # r'<configfile name="gammapysim_py">.*?', 
        # '<configfile name="gammapysim_py"><![CDATA[\n' + tool_code + '\n]]></configfile>', 
        '',
        tool_xml, 
        re.MULTILINE | re.DOTALL, 
    )

with open("gammapysim.xml.backup", "w") as f:
    f.write(tool_xml)

with open("gammapysim.xml", "w") as f:
    f.write(tool_xml_patched)