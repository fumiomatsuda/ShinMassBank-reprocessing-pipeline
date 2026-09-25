#-------------------------------------------------------------------------------
# Name:        paths.py
# Purpose:     Central resolver for dataset, record, and resource folder/file
#              paths used throughout the pipeline, plus access to the shared
#              ``settings.ini`` configuration.
#
# Author:      Fumio Matsuda@ShinMassBank PJ
#
# Created:     23/07/2023
# Copyright:   (c) hinMassBank PJ 2023
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from pathlib import Path
import glob, datetime, os, shutil, configparser
from . import massbankio

#
# Read settings.ini
#
config = configparser.ConfigParser()
config.sections()
config.read('settings.ini')

def getsettings(key1, key2):
    """Look up a value from the shared ``settings.ini`` configuration file.

    Args:
        key1 (str): Section name in ``settings.ini`` (e.g. ``"SOFTWAREPATH"``).
        key2 (str): Option name within that section (e.g. ``"SIRIUS6"``).

    Returns:
        str: The configured value for ``config[key1][key2]``.

    Raises:
        KeyError: If ``key1`` or ``key2`` is not present in ``settings.ini``.
    """
    return(config[key1][key2])


def get_alldataset_names():
    """List every dataset folder name found under the data directory.

    Globs ``../data/*`` and extracts each match's final path component (folder
    name) as the dataset identifier.

    Returns:
        list[str]: Dataset folder names.

    Note:
        Splits paths on a literal backslash (``\\``), so this assumes a Windows-
        style path separator from ``glob.glob``.
    """
    datasettemp = []
    inputpath = "../data/*"
    files = glob.glob(inputpath)

    for file in files:
        #print(file)
        datasettemp.append(file.split('\\')[-1])
    return(datasettemp)

def generate_paths_data(dataset):
    """Build the full set of folder/file paths and settings for one dataset.

    Resolves every path the pipeline needs for a given dataset: the per-step
    record folders (``step1``-``step4``) and their glob patterns, the annotation
    folder, the Hive (``hmd``) folder, per-purpose log file paths (curation,
    annotation, processing, zipf), the dataset's parsed ``RecordTemplate.txt``
    settings, its MS1/MS2/record accession-stem prefixes, and pre-computed lists
    of the ``.txt``/``.hmd`` record files present in each step folder.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).

    Returns:
        dict: Dictionary of dataset paths and settings, including (non-exhaustive):

            - ``"root"``, ``"step1path"``-``"step4path"``, ``"annotation"``,
              ``"finalpath"``, ``"msp"``, ``"hmd"``, ``"logpath"`` (str): Folder paths.
            - ``"step1input"``-``"step4input"``, ``"hmdinput"`` (str): Glob patterns
              for the corresponding folder's files.
            - ``"step1inputlist"``-``"step4inputlist"``, ``"hmdinputlist"``
              (list[str]): Pre-filtered lists of MS2 record files (or ``.hmd``
              files) actually present in each folder.
            - ``"mztab"`` (str): Path to the dataset's ``.mztab`` export file.
            - ``"settings"`` (str): The dataset's parsed ``RecordTemplate.txt``
              settings (as returned by ``massbankio.read_RecordTemplate``),
              stringified via ``str()`` along with every other value in this
              dict (see the loop below). Callers that need it back as a dict
              must ``eval()`` it themselves — see e.g. ``DDAFinalize.py``,
              which does exactly that.
            - ``"logcuration"``, ``"logannotation"``, ``"logprocessing"``,
              ``"logzipf"`` (str): Today's per-purpose log file paths.
            - ``"accession_stem"``, ``"accession_stem_ms2"``,
              ``"accession_stem_ms1"``, ``"accession_stem_record"`` (str):
              MassBank accession-ID prefixes for this dataset.
    """
    paths = {}

    datapath = "../data/"
    dataroot = Path(datapath).resolve()

    mbhpath = "../MassBank Human/"
    mbhroot = Path(mbhpath).resolve()


    paths["root"] = dataroot/dataset

    paths["step1path"] = dataroot/dataset/"step1"
    paths["step2path"] = dataroot/dataset/"step2"
    paths["step3path"] = dataroot/dataset/"step3"
    paths["step4path"] = dataroot/dataset/"step4"

    paths["annotation"] = dataroot/dataset/"annotation"

    paths["finalpath"] = mbhroot/dataset

    paths["msp"] = mbhroot/"msp"


    mztabname = "ShinMassBank_"+dataset+".mztab"
    #paths["mztab"] = dataroot/dataset/mztabname
    paths["mztab"] = mbhroot/"mzTab-M"/mztabname

    paths["step1input"] = dataroot/dataset/"step1"/"*"
    paths["step2input"] = dataroot/dataset/"step2"/"*"
    paths["step3input"] = dataroot/dataset/"step3"/"*"
    paths["step4input"] = dataroot/dataset/"step4"/"*"

    paths["hmdinput"] = dataroot/dataset/"hmd"/"*"

    paths["logpath"] = dataroot/dataset/"log"
    #paths["logfile"] =     data    $datasetid$    log    logfilename

    paths["hmd"] = dataroot/dataset/"hmd"

    accession_stem_ms2 = dataset+"P"

    settings = massbankio.read_RecordTemplate(dataroot/dataset/"RecordTemplate.txt")
    paths["settings"] = settings

    accession_stem = settings["accession_stem"]
    dt_now = datetime.datetime.now()
    for purpose in ["curation", "annotation", "processing", "zipf"]:
        logfilename = accession_stem+"_"+purpose+"_"+ dt_now.strftime('%Y.%m.%d') + ".txt"
        logkey = "log"+purpose
        paths[logkey] = dataroot/dataset/"log"/logfilename

    for key, path in paths.items():
        paths[key] = str(path)

    paths["accession_stem"] = accession_stem
    paths["accession_stem_ms2"] = accession_stem+"P"
    paths["accession_stem_ms1"] = accession_stem+"S"
    paths["accession_stem_record"] = accession_stem+"R"


    for step in ["step1", "step2", "step3", "step4"]:
        #
        # Listing-up all record file
        #
        files = glob.glob(paths[step+"input"] )
        filelist = []


        for file in files:

            if not ".txt" in file:
                continue
            if not paths["accession_stem_ms2"] in file:
                continue
            filelist.append(file)

        paths[step+"inputlist"] = filelist

        files = glob.glob(paths[step+"input"] )
        filelist = []

    filelist = []
    files = glob.glob(paths["hmdinput"] )
    for file in files:
        if not ".hmd" in file:
            continue

        filelist.append(file)

    paths["hmdinputlist"] = filelist



    return(paths)


def generate_paths_record(execdate, buildid="Build260301"):
    """Build folder/file paths for the MassBank record-curation build output.

    Resolves the curation-level output folders (Level 3a-d, Level 4), the manual-
    curation checklist CSV, and the log file for a given build/execution date.

    Args:
        execdate (str): Execution date string used to build the log filename
            (e.g. ``"2026.08.14"``).
        buildid (str, optional): Build identifier, used as the output subfolder
            name under ``../records/``. Defaults to ``"Build260301"``.

    Returns:
        dict: Dictionary of paths, including:

            - ``"recordspath"`` (str): Root records folder (``../records/``).
            - ``"level3apath"``-``"level3dpath"``, ``"level4path"`` (str):
              Per-curation-level output folders.
            - ``"recordlogpath"``, ``"recordlogfile"`` (str): Log folder and
              today's log file path.
            - ``"records"`` (str): Build-specific output folder.
            - ``"manualcuration"`` (str): Path to ``fordatacheck.csv``.
            - ``"recordmassbankpath"``, ``"recordmanualcurationpath"`` (str):
              Build-specific MassBank output and manual-curation subfolders.
    """
    paths = {}
    logfilename = "MassBank"+ execdate + ".txt"

    recordpath = "../records/"
    root = Path(recordpath).resolve()


    paths["recordspath"] = root
    paths["level3apath"] = root/"level3a"
    paths["level3bpath"] = root/"level3b"
    paths["level3cpath"] = root/"level3c"
    paths["level3dpath"] = root/"level3d"

    paths["level4path"] = root/"level4"

    paths["recordlogpath"] = root/"log"
    paths["recordlogfile"] = root/"log"/logfilename

    paths["records"] = root/buildid
    paths["manualcuration"] = root/buildid/"fordatacheck.csv"
    paths["recordmassbankpath"] = root/buildid/"Massbank"
    paths["recordmanualcurationpath"] =  root/buildid/"fordatacheck"
    for key, path in paths.items():
        paths[key] = str(path)

    return(paths)


def generate_paths_resource():
    """Build paths to shared reference-resource files used across all datasets.

    Resolves the non-natural-product list, the cached PubChem lookup pickle, the
    Zipf/FDR resource root, and the InChIKey-to-ClassyFire lookup table location.

    Returns:
        dict: Dictionary of resource paths:


            - ``"pubchemchecked"`` (str): Path to the cached PubChem lookup pickle
              (see ``resources.load_pubchemchecked``/``save_pubchemchecked``).
            - ``"zipfroot"`` (str): Root folder for Zipf/FDR pipeline outputs
              (see ``DDAFdr.py``'s ``pathsetup``/``zipfsetup``).
            - ``"inchikey2classyfire"`` (str): Path to the ClassyFire lookup table
              (see ``resources.load_inchikey2classyfire``).
    """
    paths = {}
    recordpath = "../script/res/"
    root = Path(recordpath).resolve()
    

    paths["pubchemchecked"] = root/"pubchem"/"pubchem_checked.pkl"
    paths["zipfroot"] = root/"zipf"

    softwarepath = "../software/"
    softwareroot = Path(softwarepath).resolve()
    paths["inchikey2classyfire"] = softwareroot/"InchikeyClassyfireDB-VS5.txt"

    for key, path in paths.items():
        paths[key] = str(path)

    return(paths)



def setup_annotation_folder(dataset, software, reset = False):
    """Create (or reset) the ``annotation/<software>`` output folder structure.

    Ensures ``annotation/<software>/``, its ``file/`` (tool input) and ``result/``
    (tool output) subfolders exist for the dataset. When ``reset`` is ``True``, any
    existing folders are deleted and recreated from scratch first.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).
        software (str): Name of the annotation tool (used as the output subfolder
            name, e.g. ``"MS-FINDER"``, ``"FDR_2ndrank"``).
        reset (bool, optional): When ``True``, deletes and recreates the output,
            input, and result folders before returning. Defaults to ``False``.

    Returns:
        str: Path to the tool's summary output file (``annotation/<software>/output.txt``).
    """
    datapath = "../data/"
    dataroot = Path(datapath).resolve()
    annotationpath = dataroot/dataset/"annotation"
    outputpath = annotationpath / software
    softwareinputpath = outputpath / "file"
    softwareoutputpath = outputpath / "result"
    outputfile = outputpath / "output.txt"

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
    return(outputfile)
