#!/usr/bin/env python3

'''Same as the .kt file but written in python'''

import requests
import xml.etree.ElementTree as ET

def main():
    payload = {"acc": "O22637", "output": "xml"} # TODO: what values can be passed here to get the correct protein?
    pfamXML = requests.get("http://pfam.xfam.org/protein/protein", params=payload).text

    root = ET.fromstring(pfamXML)
    print(root.tag)
    print(root.attrib)

    for child in root:
        print(child.attrib)

    print(root


if __name__=="__main__":
    main()