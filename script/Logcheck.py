#-------------------------------------------------------------------------------
# Name:        Logcheck
# Purpose:     Command-line utility that scans the per-dataset pipeline log files
#              (averaging / curation / processing / annotation) and prints their
#              start/finish status, optionally filtered by dataset, keyword, or
#              timestamp, or refreshed (cleared) entirely.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import ddatoolbox
import argparse
import glob, copy, pickle, csv, itertools, sys, datetime
from pathlib import Path
import os, re, shutil
from os.path import split
OVER_SIZE_LIMIT = 200_000_000
csv.field_size_limit(OVER_SIZE_LIMIT)


def listlogs(filelist, keyword, latest, newer):
    """Parse pipeline log files and print their "Started"/"Finished" entries.

    Reads each tab-separated log file in ``filelist``, keeps rows whose 3rd column
    is a valid ``YYYY.MM.DD.HH.MM.SS`` timestamp and whose 1st column contains
    ``"Started"`` or ``"Finished"``, then (optionally) reduces the list to only the
    single most-recent entry per ``(step, state, method)`` combination, and/or to
    only entries newer than a given timestamp. The remaining entries are sorted by
    timestamp and printed to stdout, with ``"Finished"`` states highlighted in red
    ANSI escape codes; when ``keyword`` is non-empty, only entries whose step,
    state, or raw row contains it are printed.

    Args:
        filelist (list[str]): Paths to log files to read. Each line is expected to
            be tab-separated, e.g. ``["Started"|"Finished", state, timestamp, ...]``
            or, for ``logstartstop.txt`` files,
            ``[state, step, timestamp, ...]``.
        keyword (str): Substring filter applied to the step name, state, or raw log
            row. If empty, all entries are printed.
        latest (bool): If ``True``, collapse the list to only the latest entry for
            each distinct ``(step, state, method)`` combination.
        newer (str): If non-empty, only entries whose timestamp string sorts after
            this value (e.g. ``"2025.12.31.17.49.36"``) are kept.

    Returns:
        None

    Note:
        This function references the module-level ``dataset`` variable (set by the
        surrounding ``for dataset in datasettemp:`` loop at module scope) rather
        than a parameter, so it must be called from within that loop.
    """
    #
    # Read record files
    #
    hitlist = []
    for logfilename in filelist:
        #print(logfilename)

        with open( logfilename,  encoding='UTF-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader):
                if len(row) >= 3:
                    datetimedata = row[2]
                    r'[0-9]{4}\.[0-9]{2}\.[0-9]{2}\.[0-9]{2}\.[0-9]{2}\.[0-9]{2}'
                    if re.match(r'[0-9]{4}\.[0-9]{2}\.[0-9]{2}\.[0-9]{2}\.[0-9]{2}\.[0-9]{2}', datetimedata):
                        #print(datetimedata)
                        if "Finished" in row[0]:
                            #print(logfilename, row)
                            hitlist.append([logfilename]+row)
                        if "Started" in row[0]:
                            hitlist.append([logfilename]+row)

    if latest == True:
        hitlistlatest = []
        hitlisttemp = {}
        hitlisttimetemp = {}



        for hititem in hitlist:
            #print(hititem)
            timestamp = hititem[3]

            if "logstartstop" in hititem[0]:
                step = hititem[2]
                state = hititem[1]
                if step == "annotation":
                    method = hititem[-1]
                else:
                    method = ""

            else:
                step = hititem[0].split("_")[-2]
                state = hititem[1]
                if step == "annotation":
                    method = hititem[-1]
                else:
                    method = ""

            #print(step, state, method, timestamp)
            itemid = (step, state, method)
            if not itemid in hitlisttemp:
                hitlisttemp[itemid] = hititem
                hitlisttimetemp[itemid] = timestamp
            if timestamp > hitlisttimetemp[itemid]:
                hitlisttemp[itemid] = hititem
                hitlisttimetemp[itemid] = timestamp
        for keyi , hititem in hitlisttemp.items():
            hitlistlatest.append(hititem)

        hitlist = hitlistlatest

    if not newer == "":
        hitlistlatest = []
        for hititem in hitlist:
            timestamp = hititem[3]
            #print(timestamp, newer, timestamp > newer, hititem)
            if timestamp > newer:
                hitlistlatest.append(hititem)
        hitlist = hitlistlatest




    hitlist_sorted = sorted(hitlist, key = lambda x:x[3])



    if keyword == "":
        for hititem in hitlist_sorted:

            if "logstartstop" in hititem[0]:
                step = hititem[2]
                state = hititem[1]
                if step == "annotation":
                    method = hititem[-1]
                else:
                    method = ""

            else:
                step = hititem[0].split("_")[-2]
                state = hititem[1]
                if step == "annotation":
                    method = hititem[-1]
                else:
                    method = ""



            #step = hititem[0].split("_")[-2]
            #state = hititem[1]
            RED = '\033[31m'
            END = '\033[0m'
            if "Finished" in state:
                state = RED+state+END
            print(dataset, step, state, hititem[3:])
    else:
        for hititem in hitlist_sorted:


            if "logstartstop" in hititem[0]:
                step = hititem[2]
                state = hititem[1]
                if step == "annotation":
                    method = hititem[-1]
                else:
                    method = ""

            else:
                step = hititem[0].split("_")[-2]
                state = hititem[1]
                if step == "annotation":
                    method = hititem[-1]
                else:
                    method = ""



            #step = hititem[0].split("_")[-2]
            #state = hititem[1]
            RED = '\033[31m'
            END = '\033[0m'
            if "Finished" in state:
                state = RED+state+END
            if keyword in step or keyword in state or keyword in hititem:
                print(dataset, step, state, hititem[3:])



parser = argparse.ArgumentParser(description='Progress checker for Shin-MassBank pipeline')    # Create the argument parser

parser.add_argument("target", nargs='*', default=[], help='List of target datasets. If nothing, all dataset in data foloder is processed.')    # positional argument
parser.add_argument("-c", "--check", action='store_true', help='Check all log files')    # required option
parser.add_argument("-r", "--refresh", action='store_true', help='Refresg all log files')    # required option
parser.add_argument("-f", "--latestfile", action='store_true', help='Only check latest log files')    # required option
parser.add_argument("-k", "--keyword", default="", help='Only includgin this keyword')    # required option
parser.add_argument("-l", "--latest", action='store_true', help='Only show latest log entries')    # required option
parser.add_argument("-n", "--newer", default="", help='earlier than 2025.12.31.17.49.36')    # required option
args = parser.parse_args()    # Parse command-line arguments
#print(args.target)
#print(args.check)
#print(args.refresh)
#print(args.latest)

targetdataset = args.target

keyword  = str(args.keyword)
data  = str(args.keyword)
latest = args.latest
newer = args.newer

logtypes = ["annotation", "averaging", "curation", "processing"]


if args.latestfile == True:
    # Mode: -f/--latestfile -- restrict to the two most recent log files per
    # log type (annotation/averaging/curation/processing) for each dataset.
    print("latest", targetdataset)
    #print("refresch", targetdataset)

    datasettemp = []
    if len(targetdataset) == 0:
        inputpath = "../data/*"
        files = glob.glob(inputpath)

        for file in files:
            #print(file)
            datasettemp.append(split(file)[-1])
    else:
        datasettemp = targetdataset



    for dataset in datasettemp:

        print(dataset, "**********************")
        logpath = "../data/"+dataset+"/log/*"
        #
        # Records to be analyzed
        #
        filelist = []
        files = glob.glob(logpath)
        for file in files:
            if not ".txt" in file:
                continue
            filelist.append(file)
        latestfilelist = set()

        for logtype in logtypes:
            filelisttemp = []
            for logfile in filelist:
                if logtype in logfile:
                    filelisttemp.append(logfile)
            #for logfile in sorted(filelisttemp):
            #    print(logfile)
            for logfile in sorted(filelisttemp)[-2:]:
                latestfilelist.add(logfile)
        finalfilelist = []
        for logfile in filelist:
            if logfile in latestfilelist:
                finalfilelist.append(logfile)

        listlogs(finalfilelist, keyword, latest, newer)


elif args.refresh == True:
    # Mode: -r/--refresh -- delete and recreate each dataset's log/ folder.
    #print("refresch", targetdataset)

    datasettemp = []
    if len(targetdataset) == 0:
        inputpath = "../data/*"
        files = glob.glob(inputpath)

        for file in files:
            #print(file)
            datasettemp.append(split(file)[-1])
    else:
        datasettemp = targetdataset


    for dataset in datasettemp:
        logpath = "../data/"+dataset+"/log"
        print("removed", logpath)

        if os.path.isdir(logpath):
            shutil.rmtree(logpath)
        os.mkdir(logpath)

else:
    # Default mode: list every log entry for each dataset (no filtering by
    # recency beyond what -l/--latest and -n/--newer already apply inside
    # listlogs()).
    datasettemp = []
    if len(targetdataset) == 0:
        inputpath = "../data/*"
        files = glob.glob(inputpath)

        for file in files:
            #print(file)
            datasettemp.append(split(file)[-1])
    else:
        datasettemp = targetdataset
    #print(datasettemp)

    hitlist = []
    for dataset in datasettemp:
        print(dataset, "**********************")
        inputpath = "../data/"+dataset+"/log/*"
        files = glob.glob(inputpath)


        #
        # List record files under /step4/
        #
        filelist = []
        for file in files:
            if not ".txt" in file:
                continue
            filelist.append(file)
        listlogs(filelist, keyword, latest, newer)

