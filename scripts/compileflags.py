#!/usr/bin/python2


import re
import os
import os.path
import pickle
import sys
import json

compilation_database_pattern = re.compile(
    '(?<=\s)-[DIOUWfgs][^=\s]+(?:=\\"[^"]+\\"|=[^"]\S+)?')


def load_db(filename):
    compilation_database = {}
    with open(filename) as compilation_database_file:
        compilation_database_entries = json.load(compilation_database_file)

    entry = 0
    for compilation_entry in compilation_database_entries:
        entry = entry + 1
        compilation_database[compilation_entry["file"]] = [p.strip(
        ) for p in compilation_database_pattern.findall(compilation_entry["command"])]
    return compilation_database

scriptpath = os.path.dirname(os.path.abspath(sys.argv[1]))
cache_file = "%s/cached_options.txt" % (scriptpath)

db = None
if os.access(cache_file, os.R_OK) == 0:
    db = load_db(sys.argv[1])
    f = open(cache_file, "wb")
    pickle.dump(db, f)
    f.close()
else:
    f = open(cache_file)
    db = pickle.load(f)
    f.close()

if db and sys.argv[2] in db:
    for option in db[sys.argv[2]]:
        print option
