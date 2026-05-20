"""
This class creates an instance of the Trigno base. Put your key and license here.
"""
import clr
clr.AddReference("./resources/DelsysAPI")
clr.AddReference("System.Collections")

from Aero import AeroPy

key = "MIIBKjCB4wYHKoZIzj0CATCB1wIBATAsBgcqhkjOPQEBAiEA/////wAAAAEAAAAAAAAAAAAAAAD///////////////8wWwQg/////wAAAAEAAAAAAAAAAAAAAAD///////////////wEIFrGNdiqOpPns+u9VXaYhrxlHQawzFOw9jvOPD4n0mBLAxUAxJ02CIbnBJNqZnjhE50mt4GffpAEIQNrF9Hy4SxCR/i85uVjpEDydwN9gS3rM6D0oTlF2JjClgIhAP////8AAAAA//////////+85vqtpxeehPO5ysL8YyVRAgEBA0IABPlXtISQDgqrA7ZUGVzb7HYSuXGC2jy1ZWy360TC5vLNP3IjlunUeWCaYNcG7j7Zz7AQKfKjHdNrNUxfhRumJO8="
license =   "<License>"\
            "<Id>ae7183c1-0747-45d9-ac72-0253d2df2c04</Id>"\
            "<Type>Standard</Type>"\
            "<Quantity>10</Quantity>"\
            "<LicenseAttributes>"\
            "<Attribute name=\"Software\">VS2012</Attribute>"\
            "</LicenseAttributes>"\
            "<ProductFeatures>"\
            "<Feature name=\"Sales\">True</Feature>"\
            "<Feature name=\"Billing\">False</Feature>"\
            "</ProductFeatures>"\
            "<Customer>"\
            "<Name>Vrije Universiteit Brussel</Name>"\
            "<Email>ma.diaz@vub.be</Email>"\
            "</Customer>"\
            "<Expiration>Sun, 04 Apr 2032 23:00:00 GMT</Expiration>"\
            "<Signature>MEUCIHnqnM+RGHfxzfCblxG0V8qRmvClSWRs/Mt6hFIM+786AiEAp+S3QZl4gpQRULngJ5pzCXl7McChjlSsYz/bcOLP0gc=</Signature>"\
            "</License>"

class TrignoBase():
    def __init__(self):
        self.BaseInstance = AeroPy()