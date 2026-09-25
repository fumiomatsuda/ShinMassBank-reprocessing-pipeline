#-------------------------------------------------------------------------------
# Name:        DDAAnnotation
# Purpose:     Command-line entry point that runs one compound-annotation tool
#              (MS-FINDER, SIRIUS 5/6, CFM-ID, MetFrag, GNPS/MassBank similarity
#              search) over one or more MassBank datasets.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import ddatoolbox
import glob, copy, pickle, csv, itertools, sys, datetime
from pathlib import Path
import os, re, shutil, datetime
import subprocess


def annotation(datasets, software, containerlist = [], siriusid = "", reset = False):
    """Run one compound-annotation tool over a list of datasets.

    Prepares the ``annotation/<software>`` output folder structure for each
    dataset and dispatches to the matching ``ddatoolbox.compoundsearch`` routine
    (MS-FINDER, SIRIUS 6, CFM-ID, MetFrag against MS-FINDER's DB or PubChemLite,
    or GNPS/MassBank spectral-similarity search). SIRIUS logins are performed once,
    up front, before the per-dataset loop.

    Args:
        datasets (list[str]): Dataset identifiers to process (each must have a
            ``../data/<dataset>/RecordTemplate.txt``).
        software (str): Annotation tool to run. Must be one of ``"MS-FINDER"``,
            ``"SIRIUS5"``, ``"SIRIUS6"``, ``"CFM-ID"``, ``"MetFragMSFINDER"``,
            ``"MetFragPubChemLite"``, ``"SimilarityGNPS"``, or
            ``"SimilarityMassBank"``. Any other value causes the function to return
            immediately without doing anything.
        containerlist (list[str], optional): Docker container names used to
            parallelize CFM-ID execution. Defaults to ``[]``.
        siriusid (str, optional): Login credential/user ID for SIRIUS 5/6.
            Defaults to ``""``.
        reset (bool, optional): When ``True``, deletes and recreates the
            dataset's ``annotation/<software>`` output directories before running.
            Defaults to ``False``.

    Returns:
        None
    """
    if not software in ["MS-FINDER", "SIRIUS6", "CFM-ID", "MetFragMSFINDER", "MetFragPubChemLite", "SimilarityGNPS", "SimilarityMassBank"]:
        return()

    softwarepath = "../software/"
    softwarepath = Path(softwarepath).resolve()
    sirius6path = ddatoolbox.paths.getsettings("SOFTWAREPATH", "SIRIUS6")
    sirius6path = Path(sirius6path).resolve()
    path = {}

    path["software"] = softwarepath
    path["sirius6"] = sirius6path
    if "SIRIUS5" in software:
        ddatoolbox.compoundsearch.login_sirius(5, siriusid, path)

    if "SIRIUS6" in software:
        ddatoolbox.compoundsearch.login_sirius(6, siriusid, path)

    for dataset in datasets:
        settings = ddatoolbox.massbankio.read_RecordTemplate("../data/"+dataset+"/RecordTemplate.txt")

        accession_stem = settings["accession_stem"]

        dt_now = datetime.datetime.now()
        logfilename = "../data/"+dataset+"/log/MSBNK-SMB_Measured-"+dataset+"_annotation_"+ dt_now.strftime('%Y.%m.%d') + ".txt"

        startlogfilename = "../data/"+dataset+"/log/logstartstop.txt"

        #inputpath = "../data/"+dataset+"/step2/*"
        datasetpath = "../data/"+dataset+"/"
        datasetpath = Path(datasetpath).resolve()


        # annotation output folders, per dataset
        annotationpath = datasetpath /"annotation"
        outputpath = annotationpath / software
        softwareinputpath = outputpath / "file"
        softwareoutputpath = outputpath / "result"

        if not os.path.isdir(outputpath):
            os.mkdir(outputpath)
        if not os.path.isdir(softwareinputpath):
            os.mkdir(softwareinputpath)
        if not os.path.isdir(softwareoutputpath):
            os.mkdir(softwareoutputpath)

        if reset == True:
            print(dataset, outputpath, "Reset annotation results")
            if os.path.isdir(outputpath):
                shutil.rmtree(outputpath)
            if os.path.isdir(softwareinputpath):
                shutil.rmtree(softwareinputpath)
            if os.path.isdir(softwareoutputpath):
                shutil.rmtree(softwareoutputpath)

            os.mkdir(outputpath)
            os.mkdir(softwareinputpath)
            os.mkdir(softwareoutputpath)


        path = {}
        path["dataset"] = datasetpath
        path["output"] = outputpath
        path["softwareinput"] = softwareinputpath
        path["softwareoutput"] = softwareoutputpath
        path["software"] = softwarepath
        path["sirius6"] = sirius6path

        path["logfilename"] = logfilename
        path["startlogfilename"] = startlogfilename
                #
        if "MS-FINDER" == software:
            ddatoolbox.compoundsearch.search_msfinder(dataset, path)


        if "SIRIUS6" == software:
            ddatoolbox.compoundsearch.setup_sirius(6, dataset, path)
            ddatoolbox.compoundsearch.exec_sirius(6, dataset, path)
            ddatoolbox.compoundsearch.output_sirius6(dataset, path)

        if "CFM-ID" == software:
            #containername = "zealous_wu"
            numberofcpu = len(containerlist)
            ddatoolbox.compoundsearch.search_cfmid(dataset, path, containerlist, ncpu = numberofcpu)

        #
        if "MetFragMSFINDER" == software:
            ddatoolbox.compoundsearch.search_metfrag(dataset, path, database = "msfinder")

        if "MetFragPubChemLite" == software:
            ddatoolbox.compoundsearch.search_metfrag(dataset, path, database = "PubChemLite")


        if "SimilarityGNPS" == software:
            ddatoolbox.compoundsearch.search_similarity(dataset, path, database = "GNPS")
        if "SimilarityMassBank" == software:
            ddatoolbox.compoundsearch.search_similarity(dataset, path, database = "MassBank")

        accession_stem = settings["accession_stem"]
        dt_now = datetime.datetime.now()

        #ddatoolbox.ddatoolbox.ddatoolbox.ddatoolbox.writeLog(logfilename, ["Started", dataset], mode = "w")
        #ddatoolbox.ddatoolbox.writeLog(logfilename, ["finished", accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), software])
if __name__ == '__main__':

    # Usage: python DDAAnnotation.py <software> <dataset1> [dataset2 ...]
    # Falls back to a hard-coded default (MS-FINDER on HSA001) when called with no
    # arguments.
    args = sys.argv
    if len(args) >= 2:
        software = args[1]
        datasets = args[2:]
        containerlist = []
        siriusid = ""
    else:
        datasets = ["HSA001"]
        software = "MS-FINDER"
        containerlist = []
        siriusid = ""
        
    annotation(datasets, software, containerlist = containerlist , siriusid = "", reset = False)