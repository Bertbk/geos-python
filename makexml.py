#! /usr/bin/python3

# exec(open("./makexml.py").read())

import geos_json
import writexml


data = geos_json.getData("data.json")

writexml.write_on_xml(data)
writexml.writebase(data)
