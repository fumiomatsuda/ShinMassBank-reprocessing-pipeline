#-------------------------------------------------------------------------------
# Name:        DDAaveraging
# Purpose:     Average raw MS1/MS2 spectra sharing the same unit mass (or the same
#              target-metabolite m/z window) into consensus MassBank records,
#              sequentially or in parallel across CPU cores.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
# 231004 Show number of spectra
#-------------------------------------------------------------------------------

# History
# 240410
# 240907 New folder structure

import ddatoolbox, joblib
import numpy as np
import ctypes as c
import os, shutil, glob, csv, math, pickle, itertools,re, datetime, time, sys, subprocess, configparser



def averaging_parallel(args):
    """Unpack a job-argument list and run one unit-mass averaging job.

    Thin wrapper around ``ddatoolbox.ddatoolbox.averaging()`` so that it can be
    submitted to ``joblib.Parallel``/``joblib.delayed``, which requires a
    single-argument callable.

    Args:
        args (list): A 7-element list, in order:
            ``[msmsdata, msdata, targetunitmass, paths, thresholds, filenames, report]``.

            - msmsdata (dict): MS/MS spectral data for the target unit mass.
            - msdata (dict): MS1 spectral data for the target unit mass.
            - targetunitmass (int or float): Target mass (unit mass or m/z center).
            - paths (dict): Dictionary of relevant folder/file paths.
            - thresholds (dict): Threshold parameters for clique construction and
              spectral matching.
            - filenames (dict): Output and log file names.
            - report (str): Reporting mode, ``"on"`` or ``"off"``.

    Returns:
        bool: Result of ``ddatoolbox.ddatoolbox.averaging()``; ``False`` indicates
            the averaging process failed for this unit mass.
    """
    [msmsdata, msdata, targetunitmass, paths, thresholds, filenames, report] = args
    output = ddatoolbox.ddatoolbox.averaging(msmsdata, msdata, targetunitmass, paths, thresholds, filenames, report)
    return(output)


def averaging(dataset, batchsize, numcpu, cleanstart=True):
    """Average raw spectra for a dataset into consensus MassBank records.

    Sets up the averaging environment (via ``ddatoolbox.ddatoolbox.setup_averaging``),
    builds the list of m/z boundaries to process (either evenly spaced integer
    unit-mass windows, or windows centered on a fixed list of target metabolites,
    depending on ``settings["mode"]``), then loads and averages the raw spectra in
    batches of ``batchsize`` boundaries at a time. Averaging of each unit-mass batch
    is run sequentially when ``numcpu == 1``, or in parallel across ``numcpu`` worker
    processes (via ``joblib.Parallel``) otherwise. Finally, summary statistics for
    the whole run are written out.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).
        batchsize (int): Number of mass windows loaded and processed in a single
            batch, to limit peak memory usage.
        numcpu (int): Number of CPU cores to use. If ``1``, runs sequentially;
            otherwise dispatches unit-mass jobs to a ``joblib`` process pool.
        cleanstart (bool, optional): Whether to reset previous averaging output
            before starting. Defaults to ``True``.

    Returns:
        None

    Note:
        Two boundary-generation modes are supported, selected by
        ``settings["mode"]`` (from the dataset's ``RecordTemplate.txt``):

        1. ``"unitmass"``: generates boundaries based on integer mass units,
           spanning ``settings["unitmass_start"]`` to ``settings["unitmass_end"]``.
        2. Any other value (checked-metabolites mode): generates boundaries
           centered on each ``settings["checkmetabolites"]`` target compound's
           calculated exact m/z.
    """
    #
    # Generation of information and folders
    #
    paths, thresholds, filenames, usedfiles, settings = ddatoolbox.ddatoolbox.setup_averaging(dataset, batchsize, numcpu, cleanstart = cleanstart)
    #
    # Generation of list of boundaries
    #
    listofboundaries = []
    #
    # unit mass mode
    #
    mode = settings["mode"]
    if mode == "unitmass":
        unitmass_start = int(settings["unitmass_start"])
        unitmass_end = int(settings["unitmass_end"])
        for targetunitmass in range(unitmass_start , unitmass_end+1):
            targetmz = targetunitmass * 1.0007
            lowerboudary = (targetunitmass -0.4) * 1.0007
            upperboudary = (targetunitmass +0.6) * 1.0007
            listofboundaries.append((lowerboudary, upperboudary))
    #
    # target metbolites for check purpose
    #
    else:

        targetmetabolites = settings["checkmetabolites"]
        for (metname, metmf) in targetmetabolites:
            print(metname, metmf)
            if settings["polarity"] == "NEGATIVE":
                metmfH = ddatoolbox.mfsearch.deductProton(metmf)
            else:
                metmfH = ddatoolbox.mfsearch.addProton(metmf)
            metexactmass = ddatoolbox.mfsearch.calcmw(metmfH, polarity = settings["polarity"])
            threshold = metexactmass * thresholds["ms1_mz_raw_thres"]/1000000
            lowerboudary = metexactmass - threshold
            upperboudary = metexactmass + threshold
            listofboundaries.append((lowerboudary, upperboudary))
    #
    # Generation of averaged spectra
    #
    unfinished = []
    Totalspectracounter = 0
    #
    # Splitting listofboundaries by batchsize (batche data reading)
    #
    splited = [listofboundaries[i: i+batchsize] for i in range(0, len(listofboundaries), batchsize)]
    files = glob.glob(paths["path"])

    print("loading files", len(files))
    ddatoolbox.ddatoolbox.loadhivefiles_setup(files, listofboundaries, filenames, paths)



    for listofboundaries in splited:
        #
        # Data loading. How supperting hive format.
        #
        files = glob.glob(paths["path"])
        print("loading files", len(files))

        unitmsmsdata, unitmsdata = ddatoolbox.ddatoolbox.loadspecfromtext(files, listofboundaries, filenames, paths)
        #
        # Logging
        #
        ddatoolbox.ddatoolbox.writeLog(filenames["logfilename"], ["Loaded", listofboundaries])
        for targetunitmass in unitmsmsdata.keys():
            msmsdata = unitmsmsdata[targetunitmass]
            ddatoolbox.ddatoolbox.writeLog(filenames["logfilename"], ["Loaded", targetunitmass, str(len(msmsdata))])

        print(dataset, "Averaging", numcpu)
        if numcpu == 1:
            for targetunitmass in unitmsmsdata.keys():
                msmsdata = unitmsmsdata[targetunitmass]
                msdata = unitmsdata[targetunitmass]
                Totalspectracounter += len(msmsdata)

                output = ddatoolbox.ddatoolbox.averaging(msmsdata, msdata, targetunitmass, paths, thresholds, filenames, report = "on")
                if output == False:
                    unfinished.append(targetunitmass)
        else:
            joblist = []
            report = "off"
            for targetunitmass in unitmsmsdata.keys():
                msmsdata = unitmsmsdata[targetunitmass]
                msdata = unitmsdata[targetunitmass]
                Totalspectracounter += len(msmsdata)

                joblist.append([msmsdata, msdata, targetunitmass, paths, thresholds, filenames, report])

            results = joblib.Parallel(n_jobs=numcpu, verbose=10)(joblib.delayed(averaging_parallel)(job) for job in joblist)
            for result in results:
                output = result
                if output == False:
                    unfinished.append(targetunitmass)
        print(Totalspectracounter, unfinished)
    #
    # Statistics
    #
    ddatoolbox.ddatoolbox.averaging_statistics(paths, settings, filenames, Totalspectracounter)


def hiveconvert(dataset):
    """Convert raw vendor MS files for a dataset to Hive (``.hmd``) format.

    Note:
        This function currently returns immediately after printing the configured
        Hive executable path; the conversion logic below the early ``return()`` is
        unreachable and kept only as a reference / work-in-progress implementation.
        When enabled, it would scan the dataset folder for raw vendor files
        (``.raw``, ``.lcd``, ``.wiff``, ``.d``, ``.mzML``, ``.mzXML``), run the Hive
        conversion tool on each one, and move the resulting ``.hmd`` files into an
        ``hmd/`` subfolder.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).

    Returns:
        None
    """
    #hivepath = config["SOFTWAREPATH"]["hive"]
    hivepath = ddatoolbox.paths.getsettings("SOFTWAREPATH", "hive")
    path = '../data/'+dataset+'/'
    file_path_lists = glob.glob("{}/**".format(path), recursive=True)
    hmdpath = path+"hmd/"
    if os.path.isdir(hmdpath):
        shutil.rmtree(hmdpath)
    os.mkdir(hmdpath)

    print(file_path_lists)

    filetypes = [".raw", ".lcd", ".wiff", ".d", ".mzML", ".mzXML"]

    for file in file_path_lists:
        for filetype in filetypes:
            if not filetype in file:
                continue
            #
            # If please update hive folder upon Hive updata
            #
            subprocess.run([hivepath, file])

    file_path_lists = glob.glob("{}/**".format(path), recursive=True)

    for file in file_path_lists:
        if not ".hmd" in file:
            continue

        print(file, hmdpath)
        shutil.move(file, hmdpath)



if __name__ == '__main__':
    # Usage: python DDAAveraging.py <dataset> [batchsize numcpu]
    # Falls back to a hard-coded default (HSA001, batchsize=20, numcpu=1) when
    # called with no arguments.
    args = sys.argv
    if len(args) == 2:
        dataset = args[1]
        batchsize = 10
        numcpu = 1
    elif len(args) == 4:
        dataset = args[1]
        batchsize = int(args[2])
        numcpu = int(args[3])
    else:
        dataset = "HSA001"
        batchsize = 20
        numcpu = 1
    cleanstart = False
    averaging(dataset, batchsize, numcpu, cleanstart = cleanstart)