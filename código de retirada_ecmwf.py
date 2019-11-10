# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 16:33:50 2019

@author: leope
"""

#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1994-08-17/to/2004-08-17",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "165.128/166.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "area": "-2.661778/-39.734110/-4.161778/-38.234110",
    "type": "an",
    "target": "dados_ecmwf.nc",
})
