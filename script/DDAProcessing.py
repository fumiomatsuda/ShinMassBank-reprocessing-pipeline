#-------------------------------------------------------------------------------
# Name:        DDAProcessing
# Purpose:     Check MS1 data for isotope (13C, 34S), charge-state, and adduct
#              ("ion form") signals around each MS2 record's precursor m/z, apply
#              m/z calibration, run molecular-formula and MassBank/GNPS spectral
#              similarity searches, and write the enriched (step2) records. Also
#              builds an m/z-accuracy calibration-curve summary CSV per dataset.
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#
# History
# 240907 New folder structure


import ddatoolbox
import glob, copy, pickle, csv, itertools, sys, datetime, re, ast, math
from statistics import stdev

# Significance threshold (alpha) used when de-noising a record's spectrum via
# Poisson-based correlation filtering before MassBank/GNPS similarity search.
production_correlation_alpha = 0.01


def generation_of_paths(dataset):
    """Build the folder paths and parsed record-template settings for one dataset.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).

    Returns:
        dict: Dictionary with keys:

            - ``"inputpath"`` (str): Glob pattern for ``step1`` (raw, averaged)
              record files.
            - ``"outputpath"`` (str): Output folder (``step2``) for processed
              records.
            - ``"logpath"`` (str): Dataset log folder.
            - ``"record_template"`` (str): Path to the shared
              ``recordformat.txt`` template.
            - ``"settings"`` (dict): The dataset's parsed
              ``RecordTemplate.txt`` settings, as returned by
              ``ddatoolbox.massbankio.read_RecordTemplate``.
    """
    paths = {}
    paths["inputpath"] = "../data/"+dataset+"/step1/*"
    paths["outputpath"] = "../data/"+dataset+"/step2/"
    paths["logpath"] = "../data/"+dataset+"/log"
    paths["record_template"] = "./recordformat.txt"
    #paths["record_templatefile"] = "./"+dataset+"/RecordTemplate.txt"
    paths["settings"] = ddatoolbox.massbankio.read_RecordTemplate("../data/"+dataset+"/RecordTemplate.txt")
    return(paths)

def generation_of_recordlist(dataset, paths):
    """List every MS2 (step1) record file for a dataset, paired with its MS1 sibling.

    Args:
        dataset (str): Dataset identifier (unused directly here; kept for
            interface consistency with sibling ``generation_of_*`` helpers).
        paths (dict): Paths dict, as returned by :func:`generation_of_paths`;
            must contain ``"inputpath"`` (glob pattern) and
            ``"settings"`` (parsed ``RecordTemplate.txt``, for
            ``"accession_stem"``).

    Returns:
        list[tuple]: List of ``(ms2_file, ms1_file)`` path pairs. ``ms1_file``
            is derived from ``ms2_file`` by swapping the MS2 accession-stem
            suffix (``"P"``) for the MS1 suffix (``"S"``); it is not verified
            to actually exist.
    """
    accession_stem = paths["settings"]["accession_stem"]

    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"
    filelist = []

    files = glob.glob(paths["inputpath"])

    for file in files:
        if not ".txt" in file:
            continue

        if not accession_stem_ms2 in file:
            continue
        # Derive the paired MS1 (precursor-scan) record filename from the MS2
        # filename by swapping the accession-stem suffix ("P" -> "S").
        ms1file = file.replace(accession_stem_ms2, accession_stem_ms1)
        filelist.append((file, ms1file))

    return(filelist)


def check_S(spec, ion, threshold):
    """Search MS1 for sulfur-isotope (34S), 13C-isotope, and charge-state satellites of an ion.

    Around ``ion``'s m/z, looks for six kinds of correlated satellite peaks in
    ``spec`` (the corrected, averaged MS1 spectrum): a 34S isotope peak
    (``+1.996`` Da, less intense than ``ion`` and at least 3% of its
    intensity -- used later to bound the maximum sulfur count for formula
    search), 13C isotope peaks consistent with a +1/+2/+3 charge state
    (``+1.003``, ``+0.502``, ``+0.334`` Da respectively, each less intense
    than ``ion``), and +1/+2 13C isotopologue peaks of ``ion`` itself (more
    intense than ``ion``, i.e. ``ion`` itself is the isotope and the peak
    found here is its more abundant, lower-mass sibling). A candidate only
    counts as a hit if its intensity-profile-vs-scan correlation with
    ``ion`` (via ``ddatoolbox.massbankio.calccorrel``) exceeds 0.6; among
    candidates passing that bar, the most intense one for each category is
    kept.

    Args:
        spec (dict): Corrected MS1 spectrum, mapping m/z to a dict with keys
            ``"intensity"``, ``"number"``, and ``"intlist"`` (a
            stringified list, ``eval``-parsed here).
        ion (dict): Query ion (typically the precursor found in MS1 via
            :func:`searchion`), with keys ``"mz"``, ``"intensity"``, and
            ``"intlist"``.
        threshold (float): m/z matching tolerance (Da) for each satellite search.

    Returns:
        tuple: Six dicts, each either populated with keys ``"mzdiff"``,
            ``"intensity"``, ``"mz"``, ``"number"``, ``"ratio"``, ``"cor"``,
            or ``{}`` if no qualifying peak was found:

            - mostintenseion: 34S isotope candidate.
            - singlediffmostintenseion: +1 charge-state (13C, singly charged) candidate.
            - doublediffmostintenseion: +2 charge-state candidate.
            - triplediffmostintenseion: +3 charge-state candidate.
            - singleisotopemostintenseion: ``ion``'s own +1 13C isotopologue partner.
            - doubleisotopemostintenseion: ``ion``'s own +2 13C isotopologue partner.

        Returns ``({}, {}, {}, {}, {}, {})`` immediately if ``ion["mz"] == 0``
        (no precursor was found in MS1).
    """
    if ion["mz"] == 0:
        return({}, {},{}, {},{}, {})

    mz = ion["mz"]
    intensity = ion["intensity"]
    number = ion["number"]
    intlist = ion["intlist"]

    Sdiff = 1.99579584
    Cdiff = 1.003354
    triplediff = Cdiff/3
    doublediff = Cdiff/2
    singlediff = Cdiff


    #
    # Finding most intense ion within a threshold level
    #
    mostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}
    triplediffmostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}
    doublediffmostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}
    singlediffmostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}

    doubleisotopemostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}
    singleisotopemostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}


    for tmz, datahash in spec.items():
        tintensity = datahash["intensity"]
        tnumber = datahash["number"]
        tintlist = eval(datahash["intlist"])

        #
        # Diff mz
        #
        mzdiff = abs(tmz - mz - Sdiff)
        ratio = tintensity / intensity
        if mzdiff <= threshold and tintensity < intensity and ratio > 0.03:
            correl = ddatoolbox.massbankio.calccorrel(intlist, tintlist)
            if correl > 0.6:
                mostintenseion["mzdiff"] = mzdiff
                mostintenseion["intensity"] = tintensity
                mostintenseion["mz"] = tmz
                mostintenseion["number"] = tnumber
                mostintenseion["ratio"] = ratio
                mostintenseion["cor"] = correl

        mzdiff = abs(tmz - mz - triplediff)
        if mzdiff <= threshold and tintensity < intensity:
            correl = ddatoolbox.massbankio.calccorrel(intlist, tintlist)
            if correl > 0.6:
                triplediffmostintenseion["mzdiff"] = mzdiff
                triplediffmostintenseion["intensity"] = tintensity
                triplediffmostintenseion["mz"] = tmz
                triplediffmostintenseion["number"] = tnumber
                triplediffmostintenseion["ratio"] = ratio
                triplediffmostintenseion["cor"] = correl


        mzdiff = abs(tmz - mz - doublediff)
        if mzdiff <= threshold and tintensity < intensity:
            correl = ddatoolbox.massbankio.calccorrel(intlist, tintlist)
            if correl > 0.6:
                doublediffmostintenseion["mzdiff"] = mzdiff
                doublediffmostintenseion["intensity"] = tintensity
                doublediffmostintenseion["mz"] = tmz
                doublediffmostintenseion["number"] = tnumber
                doublediffmostintenseion["ratio"] = ratio
                doublediffmostintenseion["cor"] = correl

        mzdiff = abs(tmz - mz - singlediff)
        if mzdiff <= threshold and tintensity < intensity:
            correl = ddatoolbox.massbankio.calccorrel(intlist, tintlist)
            if correl > 0.6:
                singlediffmostintenseion["mzdiff"] = mzdiff
                singlediffmostintenseion["intensity"] = tintensity
                singlediffmostintenseion["mz"] = tmz
                singlediffmostintenseion["number"] = tnumber
                singlediffmostintenseion["ratio"] = ratio
                singlediffmostintenseion["cor"] = correl


        mzdiff = abs(tmz - mz + Cdiff)
        if mzdiff <= threshold and tintensity > intensity:
            correl = ddatoolbox.massbankio.calccorrel(intlist, tintlist)
            if correl > 0.6:
                singleisotopemostintenseion["mzdiff"] = mzdiff
                singleisotopemostintenseion["intensity"] = tintensity
                singleisotopemostintenseion["mz"] = tmz
                singleisotopemostintenseion["number"] = tnumber
                singleisotopemostintenseion["ratio"] = ratio
                singleisotopemostintenseion["cor"] = correl

        mzdiff = abs(tmz - mz + 2*Cdiff)
        if mzdiff <= threshold and tintensity > intensity:
            correl = ddatoolbox.massbankio.calccorrel(intlist, tintlist)
            if correl > 0.6:
                doubleisotopemostintenseion["mzdiff"] = mzdiff
                doubleisotopemostintenseion["intensity"] = tintensity
                doubleisotopemostintenseion["mz"] = tmz
                doubleisotopemostintenseion["number"] = tnumber
                doubleisotopemostintenseion["ratio"] = ratio
                doubleisotopemostintenseion["cor"] = correl


    if mostintenseion["mz"] == 0:
        mostintenseion = {}

    if singlediffmostintenseion["mz"] == 0:
        singlediffmostintenseion = {}

    if doublediffmostintenseion["mz"] == 0:
        doublediffmostintenseion = {}

    if triplediffmostintenseion["mz"] == 0:
        triplediffmostintenseion = {}

    if singleisotopemostintenseion["mz"] == 0:
        singleisotopemostintenseion = {}

    if doubleisotopemostintenseion["mz"] == 0:
        doubleisotopemostintenseion = {}

    return(mostintenseion, singlediffmostintenseion, doublediffmostintenseion, triplediffmostintenseion, singleisotopemostintenseion, doubleisotopemostintenseion)


def searchion(spec, mz, threshold):
    """Locate a target m/z (typically the MS2 precursor) within an MS1 spectrum.

    Scans every peak in ``spec`` within ``threshold`` Da of ``mz`` and tracks
    two candidates: the most intense matching peak (``biggest``) and the
    closest-in-m/z matching peak (``closest``). Returns ``biggest`` when any
    match was found (preferring intensity over exact m/z proximity), falling
    back to ``closest`` (which stays at its untouched default, ``intensity``
    ``0``, ``mzdiff`` ``1.0``, if nothing matched at all) otherwise.

    Args:
        spec (dict): MS1 spectrum, mapping m/z to a dict with keys
            ``"intensity"``, ``"number"``, and ``"intlist"`` (a stringified
            list, ``eval``-parsed here).
        mz (float): Target m/z value to search for (typically the MS2
            record's corrected precursor m/z).
        threshold (float): m/z matching tolerance (Da); peaks farther than
            this from ``mz`` are ignored.

    Returns:
        dict: The matched ion, with keys ``"mz"`` (float), ``"intensity"``,
            ``"mzdiff"`` (``mz`` in MS1 minus the query ``mz``), ``"number"``,
            and ``"intlist"`` (list, not stringified). If no peak was within
            ``threshold``, returns the ``closest`` placeholder with
            ``"intensity"`` ``0`` and ``"mz"`` ``0``.
    """
    #annotation_hash = self.get_annotation_from_items()
    #
    # Finding most intense ion within a threshold level
    #


    biggest = {"mz":0, "intensity":0, "mzsd":0, "intsd":0, "number":0, "mzdiff":0, "intlist":""}
    closest = {"mz":0, "intensity":0, "mzsd":0, "intsd":0, "number":0, "mzdiff":1.0, "intlist":""}


    for tmz, datahash in spec.items():
        #
        # Diff mz
        #
        tintensity = datahash["intensity"]
        tnumber = datahash["number"]
        tlist = eval(datahash["intlist"])
        mzdiff = abs(tmz - mz)

        if mzdiff < threshold:

            if biggest["intensity"] < tintensity:
                biggest["mzdiff"] = mzdiff
                biggest["intensity"] = tintensity
                biggest["mz"] = float(tmz)

                biggest["number"] = tnumber
                biggest["intlist"] = tlist
            if closest["mzdiff"] > mzdiff:
                closest["mzdiff"] = mzdiff
                closest["intensity"] = tintensity
                closest["mz"] = float(tmz)

                closest["number"] = tnumber
                closest["intlist"] = tlist



    if biggest["intensity"] == 0:
        returndata = closest
    else:
        returndata = biggest

    return(returndata)


def check_ionform_in_MS1(mb1, mb2, precursorinms1, thresholds):
    """Check MS1 for alternative adduct ("ion form") satellites of the precursor.

    Always includes the default protonation/deprotonation adduct (``"M+H"``
    for positive mode, ``"M-H"`` for negative), then delegates to
    ``mb1.check_ionform_new`` to search MS1 for peaks matching every other
    configured ion form (e.g. ``M+Na``, ``M+NH4``) relative to
    ``precursorinms1``. Each additional ion form found is appended as a
    ``PROPOSED_ION_FORM`` comment on ``mb2`` (the MS2 record being built).

    Args:
        mb1 (ddatoolbox.massbankio.MassBankRecord): The MS1 record (source of
            the ``PROPOSED_ION_FORM`` comment's format string).
        mb2 (ddatoolbox.massbankio.MassBankRecord): The MS2 record being
            annotated; receives the new ``PROPOSED_ION_FORM`` comment lines.
        precursorinms1 (dict): The precursor ion as found in MS1, as returned
            by :func:`searchion`.
        thresholds (dict): m/z tolerance thresholds, as consumed by
            ``mb1.check_ionform_new``.

    Returns:
        tuple: A tuple of two elements.

            - ionformdata (dict): The loaded ion-form reference table (mass
              differences and multipliers), as returned by
              ``ddatoolbox.massbankio.read_ionformlist``.
            - ionform_list (list[str]): Ion forms detected for this
              precursor (always includes the default ``"M+H"``/``"M-H"``,
              plus any additional forms found in MS1).
    """
    #
    # Checking ion form
    #
    polarity = mb1.get_polarity()


    if polarity == "POSITIVE":
        ionform_list = ["M+H"]
    else:
        ionform_list = ["M-H"]

    thiskey = ("COMMENT:","PROPOSED_ION_FORM")

    ionformdata = ddatoolbox.massbankio.read_ionformlist(filename = "ion_forms.txt", polarity = polarity)

    proposedionforms = mb1.check_ionform_new(precursorinms1, thresholds, ionformdata)
    #ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, "PROPOSED_ION_FORM", proposedionforms])

    for ionform, datahash in proposedionforms.items():

        if1, if2 = ionform.split("_")
        ionform_list.append(if1)
        #print("PROPOSED_ION_FORM", accession, ionform)

        #texttemp = "{} due to occurence of {} at m/z {:.7f} in MS1, int {:.3f}, cor with prec {:.3f}".format(if1, if2, datahash['mz'], datahash['intensity'], datahash['cor'])
        formattext = mb1.record[thiskey]["format"]
        #print(formattext, if1, if2, datahash['mz'], datahash['intensity'], datahash['cor'])
        texttemp = formattext.format(if1, if2, datahash['mz'], datahash['intensity'], datahash['cor'])
        mb2.add2("COMMENT:","PROPOSED_ION_FORM", texttemp)


    #print("PROPOSED_ION_FORM", ionform_list)
    #ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, "PROPOSED_ION_FORM", ionform_list])
    return(ionformdata, ionform_list)




def processing(datasets):
    """Annotate MS1 signals onto every dataset's MS2 records and build a calibration report.

    Runs two phases per dataset:

    1. **Per-record MS1 annotation and search** (``step1`` -> ``step2``): for
       every MS2 record with enough raw spectra, loads its paired MS1 record,
       corrects both records' m/z values (:func:`ddatoolbox.mfsearch.mzcorrection`),
       locates the precursor in MS1 (:func:`searchion`), checks for sulfur-isotope/
       charge-state/13C-isotope satellites (:func:`check_S`) and alternative ion
       forms (:func:`check_ionform_in_MS1`), searches candidate molecular formulas
       for the (possibly adduct-shifted) precursor mass
       (:func:`ddatoolbox.mfsearch.mfsearch`), searches the local MassBank/GNPS
       spectral database (:func:`ddatoolbox.massbankserver.searchMassBank`), and
       writes the enriched record to ``step2``.
    2. **Calibration-curve summary** (``step1``/``step2`` -> ``<dataset>_mzcomparison.csv``):
       for every ``step2`` record with an existing MassBank/MF hit, looks up the
       matching MassBank record's formula, computes the theoretical precursor and
       fragment m/z values, and tabulates the raw/averaged/calibrated m/z deviation
       (in ppm) for both the precursor ion and every well-matched, non-isotope
       fragment. Datasets whose hit rate for a record falls at or below 50% (and
       isn't "all but one" hit) are excluded from the summary for that record.

    Args:
        datasets (list[str]): Dataset identifiers to process.

    Returns:
        None

    See Also:
        :func:`generation_of_paths`, :func:`generation_of_recordlist`: build the
        per-dataset paths/file list consumed here.
    """
    #
    # Folder including Hive files
    #
    minimunnumber = 10
    for dataset in datasets:
        execdate = datetime.datetime.now().strftime('%Y.%m.%d')
        logfilename = "../data/"+dataset+"/log/MSBNK-SMB_Measured-"+dataset+"_processing_"+ execdate + ".txt"

        paths = generation_of_paths(dataset)
        filelist = generation_of_recordlist(dataset, paths)

        accession_stem = paths["settings"]["accession_stem"]
        dt_now = datetime.datetime.now()

        #ddatoolbox.ddatoolbox.writeLog(logfilename, ["Started", dataset], mode = "w")
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Started", accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S')], mode = "w")


        # MassBank data
        mbdict = ""
        #
        totalfilenumber = len(filelist)
        progresscounter = 0



        for ms2filename, ms1filename in filelist:
            progresscounter += 1

            #if not "P0150001" in ms2filename:

            #if not "0170001" in ms2filename:
            #    continue


            print("{:s} {:d}/{:d} {:s}".format(dataset, progresscounter, totalfilenumber, ms2filename))
            #
            # Loading MassBank recored MS2
            #
            mb2 = ddatoolbox.massbankio.MassBankRecord()
            mb2.read_MassBankRecord(ms2filename)
            #
            # Number of raw MS2 specs have to be greater than minimun number
            #
            number = mb2.get_number()
            if number < minimunnumber:
                continue
            #
            # Setting corrections and mz_averaged_thres
            #
            mb2.add2("COMMENT:","MS1_M/Z_CORRECTION", paths["settings"]["ms1_mz_correction"])
            mb2.add2("COMMENT:","MS2_M/Z_CORRECTION", paths["settings"]["ms2_mz_correction"])
            mb2.add2("COMMENT:","MS1_M/Z_AVERAGED_THRES", paths["settings"]["ms1_mz_averaged_thres"])
            mb2.add2("COMMENT:","MS2_M/Z_AVERAGED_THRES", paths["settings"]["ms2_mz_averaged_thres"])
            #
            # Getting information
            #
            mz, mzstdev = mb2.get_precmz()
            number = mb2.get_number()
            accession = mb2.get_accession()
            thresholds = mb2.get_mzthresholds_averaged()
            polarity = mb2.get_polarity()

            #
            # I found a mistake here
            #precmz_corrected = ddatoolbox.mfsearch.mzcorrection(mz, paths["settings"]["ms2_mz_correction"])
            precmz_corrected = ddatoolbox.mfsearch.mzcorrection(mz, paths["settings"]["ms1_mz_correction"])
            #
            # Loading MassBank recored MS1
            #
            mb1 = ddatoolbox.massbankio.MassBankRecord()
            mb1.read_MassBankRecord(ms1filename)
            #
            # Correction of ms1 m/z
            #
            annotationlist = mb1.get_annotation_from_items()
            annoaationlist_corrected = {}
            spec_corrected = {}
            spec_corrected_extracted = {}

            for mzraw, datahash in annotationlist.items():
                mz_corrected = ddatoolbox.mfsearch.mzcorrection(mzraw, paths["settings"]["ms1_mz_correction"])
                annoaationlist_corrected[mz_corrected] = datahash
                annoaationlist_corrected[mz_corrected]["mzraw"] = mzraw
                inttemp = annoaationlist_corrected[mz_corrected]["intensity"]
                spec_corrected[mz_corrected] = inttemp

                if precmz_corrected -3 < mz_corrected < precmz_corrected+3:
                    spec_corrected_extracted[mz_corrected] = datahash

            mb1.add_annotation_from_items(annoaationlist_corrected, ["int", "annotation", "m/z_nocal", "n", "intlist"])
            mb1.add_spectra(spec_corrected)

            precms1threshold = ddatoolbox.massbankio.calcThresholdAveraged(precmz_corrected, "ms1", thresholds) * 2
            #
            #
            # Checking occurence of prec ion in MS1
            #precursorinms1, outputtext = search_prec_in_MS1(spec_corrected_extracted, precmz_corrected, precms1threshold)
            #
            precursorinms1 = searchion(spec_corrected_extracted, precmz_corrected, precms1threshold)
            #
            # To avoid zeto divition in calculating relative intensity ratio
            #
            if precursorinms1["intensity"] == 0:
                precursorinms1["intensity"] = 1.0
            #
            # Addition to MassBank record
            #
            if not precursorinms1["mz"] == 0:
                ppmtemp = (precursorinms1["mz"]-precmz_corrected)/precursorinms1['mz']*1000000
                thiskey = ("COMMENT:","PRECURSOR_IN_MS1")
                formattext = mb1.record[thiskey]["format"]
                outputtext = formattext.format(precursorinms1['mz'], precursorinms1['intensity'], ppmtemp)
            else:
                outputtext = "ND"

            mb2.add2("COMMENT:","PRECURSOR_IN_MS1", outputtext)
            #print("COMMENT:","PRECURSOR_IN_MS1", outputtext)


            # Checking occurence of sulfar atom
            #Smax = check_S_in_MS1(mb1, mb2, precursorinms1, thresholds)
            #
            # Checking occurence of sulfar atom
            sulfar, singlecs, doublecs, triplecs, singleit, doubleit = check_S(spec_corrected_extracted, precursorinms1, precms1threshold)
            #
            # Setting maximal S number for MF search
            #
            if sulfar != {}:
                Sratio = sulfar["ratio"]
                Smax = int(Sratio / 0.04 + 0.3)

            else:
                Smax = 0

            #
            # Addition to MassBank record
            #
            if sulfar != {}:
                thiskey = ("COMMENT:","PROPOSED_MAX_NUM_S")
                formattext = mb1.record[thiskey]["format"]
                outputtext = formattext.format(Smax, sulfar['mz'], Sratio, sulfar['cor'])

            else:
                outputtext = "{}".format(Smax)
            mb2.add2("COMMENT:","PROPOSED_MAX_NUM_S", outputtext)
            #print("COMMENT:","PROPOSED_MAX_NUM_S", outputtext)
            #
            #
            # charge state
            if polarity == "POSITIVE":
                symbol = "+"
                symbolp = 1
            else:
                symbol = "-"
                symbolp = -1

            thiskey = ("COMMENT:","PROPOSED_CHARGE_STATE")
            formattext = mb1.record[thiskey]["format"]

            if triplecs != {}:
                chargestate = triplecs
                texttemp = formattext.format(3*symbolp, chargestate['mz'], chargestate['intensity'], chargestate['cor'])
            elif doublecs != {}:
                chargestate = doublecs
                texttemp = formattext.format(2*symbolp, chargestate['mz'], chargestate['intensity'], chargestate['cor'])
            elif singlecs != {}:
                chargestate = singlecs
                texttemp = formattext.format(symbolp, chargestate['mz'], chargestate['intensity'], chargestate['cor'])
            else:
                texttemp = "ND"
            mb2.add2("COMMENT:","PROPOSED_CHARGE_STATE",texttemp)
            #print("COMMENT:","PROPOSED_CHARGE_STATE",texttemp)


            thiskey = ("COMMENT:","PROPOSED_ISOTOPE_STATE")
            formattext = mb1.record[thiskey]["format"]


            if doubleit != {}:
                chargestate = doubleit
                texttemp = formattext.format(2, chargestate['mz'], chargestate['intensity'], chargestate['cor'])
            elif singleit != {}:
                chargestate = singleit
                texttemp = formattext.format(1, chargestate['mz'], chargestate['intensity'], chargestate['cor'])
            else:
                texttemp = "ND"
            mb2.add2("COMMENT:","PROPOSED_ISOTOPE_STATE",texttemp)
            #print("COMMENT:","PROPOSED_ISOTOPE_STATE",texttemp)

            # Checking ion from, isotope, and charge states in MS1
            ionformdata, ionform_list = check_ionform_in_MS1(mb1, mb2,precursorinms1, thresholds)
            #
            # Correction of precursor ion m/z in MS2
            #
            mz_corrected = ddatoolbox.mfsearch.mzcorrection(mz, thresholds["ms1_mz_correction"])
            mb2.add2("COMMENT:","PRECURSOR_M/Z_RAW", mz)
            mb2.add2("MS$FOCUSED_ION:","PRECURSOR_M/Z", mz_corrected)

            #
            # Correction of fragment ion m/z  in MS2
            #
            annotationlist = mb2.get_annotation_from_items()
            annoaationlist_corrected = {}
            spec_corrected = {}

            for mzraw, datahash in annotationlist.items():
                mz_corrected = ddatoolbox.mfsearch.mzcorrection(mzraw, thresholds["ms2_mz_correction"])
                annoaationlist_corrected[mz_corrected] = datahash
                annoaationlist_corrected[mz_corrected]["mzraw"] = mzraw
                spec_corrected[mz_corrected] = annoaationlist_corrected[mz_corrected]["intensity"]
            mb2.add_annotation_from_items(annoaationlist_corrected, ["int", "annotation", "m/z_nocal", "m/z_std(ppm)", "int_std", "cor", "n"])

            mb2.add_spectra(spec_corrected)
            mb2.clear_anything("COMMENT:", "PRECURSOR_M/Z_LIST")
            #
            # Setting Up MassBank Server (Only one time)
            #
            polarity = mb2.get_polarity()
            if mbdict == "":

                #if polarity == "NEGATIVE":
                #    path = "./res/spectra/MassBanknegative/"
                #else:
                #    path = "./res/spectra/MassBankpositive/"

                if polarity == "NEGATIVE":
                    path = "./res/spectra/GNPSnegative/"
                else:
                    path = "./res/spectra/GNPSpositive/"

                mbdict = ddatoolbox.massbankserver.setupMassBank(path)


            #
            # Molecular formula searching of precursor ion
            #
            threshold_temp = thresholds["ms1_mz_averaged_thres"]
            mz, mzstdev = mb2.get_precmz()
            hitmfs = []

            for ion_form in ionformdata.keys():
                diff_temp  = ionformdata[ion_form]["mass"]
                diff_mult  = ionformdata[ion_form]['mult']

                if ion_form in ionform_list:
                    mz_Na = (mz - diff_temp)/diff_mult

                    result, result_list = ddatoolbox.mfsearch.mfsearch(mz_Na, threshold = threshold_temp, ihdtype = "integer", Climit = 99, Hlimit = 999, Nlimit = 5, Olimit = 99, Plimit = 2, Slimit = Smax, polarity = "None")
                    difflist = ["{:6f}".format(ddatoolbox.mfsearch.calcmw(mf) - mz_Na)  for mf in result_list]

                    logtemp = {}
                    logtemp["mz"] = mz_Na
                    logtemp["threshold"] = threshold_temp
                    logtemp["mz"] = mz_Na
                    hitlist = []
                    for (formula_temp, diff_temp) in zip(result_list, difflist):
                        hitlist.append((formula_temp, diff_temp))

                    logtemp["hits"] = hitlist

                    ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, "MFSEARCH_PREC", logtemp])

                    for hitmf, hitdiff in zip(result_list, difflist):
                        hitdiffppm = float(hitdiff)/mz_Na * 1000000

                        if polarity == "POSITIVE":
                            ionform = "[{:s}]+".format(ion_form)
                        else:
                            ionform = "[{:s}]1".format(ion_form)

                        texttemp = "{} of {} with delta m/z {:4f} ppm".format(ddatoolbox.mfsearch.cleanmf(hitmf), ion_form, hitdiffppm)
                        mb2.add2("COMMENT:","PROPOSED_MF_FROM_PRECMZ", texttemp)
                        name = "SevenGoldenRules"
                        formula = ddatoolbox.mfsearch.cleanmf(hitmf)
                        score = abs(hitdiffppm)

                        formattext = mb2.record[("COMMENT:", "PROPOSED_MF")]["format"]
                        texttemp = formattext.format(name, score, formula, ionform, float(hitdiffppm) * -1)
                        mb2.add2("COMMENT:", "PROPOSED_MF",texttemp)




                        hitmfs.append(ddatoolbox.mfsearch.cleanmf(hitmf))
            #
            # MassBank searcing
            #
            cleanedspectra = mb2.wash_spectra(alpha = production_correlation_alpha, intthres = 0)

            results = ddatoolbox.massbankserver.searchMassBank(mbdict, cleanedspectra, scorethreshold = 0.7, mzthreshold = 0.1, digit = 2, path = path)

            hittedMassBankrecordid = {}
            for sid, hittedspectra in results.items():
                hittedformula = ddatoolbox.mfsearch.cleanmf(hittedspectra["FORMULA"])
                hittedthormz = ddatoolbox.mfsearch.calcmw(hittedformula, polarity)
                recordtitle = hittedspectra["RECORD_TITLE"]
                score = hittedspectra["score"]
                inchikey = hittedspectra["INCHIKEY"]
                if inchikey == "":
                    continue
                if hittedformula in hitmfs:

                    #
                    # log file
                    #
                    logtemp = {}
                    logtemp["score"] = score
                    logtemp["id"] = sid
                    logtemp["recordtitle"] = recordtitle
                    logtemp["formula"] = hittedformula
                    logtemp["scorethreshold"] = 0.8
                    logtemp["mzthreshold"] = 0.1
                    logtemp["hittedthormz"] = hittedthormz

                    logtemp["inchikey"] = inchikey

                    ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, "MassBank", logtemp])



                    #ddatoolbox.ddatoolbox.writeLog(logfilename, [accession, "Cosine product", "MassBank", 0.7, 0.1, sid, score, hittedformula, recordtitle, inchikey])

                    if not inchikey in hittedMassBankrecordid:
                        hittedMassBankrecordid[inchikey]={
                        "recordtitle": recordtitle,
                        "score":score,
                        "sid": sid,
                        "mf":hittedformula,
                        "counter": 1
                        }
                    else:
                        hittedMassBankrecordid[inchikey]["counter"] += 1
            #
            # Sort result
            #
            counter = 0
            for inchikey in sorted(hittedMassBankrecordid.keys(), key=lambda x: hittedMassBankrecordid[x]["counter"] + hittedMassBankrecordid[x]["score"], reverse=True):
                recordtitle=hittedMassBankrecordid[inchikey]["recordtitle"]
                score=hittedMassBankrecordid[inchikey]["score"]
                sid=hittedMassBankrecordid[inchikey]["sid"]
                mf=hittedMassBankrecordid[inchikey]["mf"]
                hitcounter= hittedMassBankrecordid[inchikey]["counter"]
                print(dataset, "MassBank", inchikey, hitcounter, score, mf, sid, recordtitle)
                texttemp ="{} Total {:d} hit with max score {:.3f} of {} {} {}".format(inchikey, hitcounter, score, mf, sid, recordtitle)

                mb2.add2("COMMENT:","PROPOSED_FROM_SPEC_SIMILARITY", texttemp)

                #text = "{:s} {:s} {:.3f} {:s} {:s}".format("CosineProduct_GNPS", inchikey, hitcounter, mf, recordtitle)
                #mb2.add2("COMMENT:", "PROPOSED_COMPOUND", text)

                formattext = mb2.record[("COMMENT:", "PROPOSED_COMPOUND")]["format"]
                texttemp = formattext.format("SIMILARITY", "GNPS", hitcounter, inchikey, mf, recordtitle)
                mb2.add2("COMMENT:", "PROPOSED_COMPOUND",texttemp)



                counter += 1
                if counter > 10:
                    break

            #
            # Addition of list of original spectra
            #
            originals = ast.literal_eval(mb2.get_anything("COMMENT:", "COUNTER"))
            #print(originals.keys())
            #keytemp = list(originals['temp_record'].keys())[0]
            originalspecs = originals['intentisy']

            mb2.add2("COMMENT:", "RAW_SPEC", originalspecs)



            # clear RETENTION_TIME_LIST COUNTER
            mb2.clear_anything("COMMENT:", "RETENTION_TIME_LIST")
            mb2.clear_anything("COMMENT:", "COUNTER")
            mb2.clear_anything("COMMENT:", "MS1_M/Z_RAW_THRES")
            mb2.clear_anything("COMMENT:", "MS2_M/Z_RAW_THRES")
            #mb2.clear_anything("COMMENT:", "MS1_M/Z_AVERAGED_THRES")
            #mb2.clear_anything("COMMENT:", "MS2_M/Z_AVERAGED_THRES")
            mb2.clear_anything("COMMENT:", "MS1_M/Z_CORRECTION")
            mb2.clear_anything("COMMENT:", "MS2_M/Z_CORRECTION")
            mb2.clear_anything("COMMENT:", "DOTPRODUCT_THRES")
            mb2.clear_anything("COMMENT:", "RETENTION_TIME_THRES")



            #
            mb2.add_date(date = execdate, target = "Processed")
            mb2.show(path = paths["outputpath"])


        dt_now = datetime.datetime.now()
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Finished", accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S')])

        #
        # Annotation check
        # From here on: build the m/z-accuracy calibration-curve summary.
        #
        #

        dt_now = datetime.datetime.now()
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Calibration Started", accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S')])


        outputdata = []
        outputdata.append(["id", "MassBank hit", "MS", "type", "MF", "m/z measured", "m/z theol", "num averaging",\
        "diff prec raw corrected", "diff prec ave corrected", "diff frag raw corrected", "diff frag ave corrected", "diff prec raw", "diff prec ave", "diff frag raw", "diff frag ave"])

        inputpath1 = "../data/"+dataset+"/step1/*"
        inputpath2 = "../data/"+dataset+"/step2/*"
        outputfile = "../data/"+dataset+"/"+dataset+"_mzcomparison.csv"

        record_template = "./recordformat.txt"



        accession_stem = paths["settings"]["accession_stem"]

        accession_stem_ms2 = accession_stem+"P"
        accession_stem_ms1 = accession_stem+"S"

        filelist = []

        prec_ave = []
        prec_raw =[]
        frag_ave = []
        frag_raw = []


        counterall = 0
        counterhit = 0
        new_dict = ""

        #
        # Collect the step2 records to be analyzed.
        #
        files = glob.glob(inputpath2)
        for file in files:
            if not ".txt" in file:
                continue
            if not accession_stem_ms2 in file:
                continue
            step1file = file.replace("step2", "step1")
            filelist.append((file, step1file))

        #
        # Load each target record individually.
        #
        for ms2filename, step1file in filelist:

            mb2 = ddatoolbox.massbankio.MassBankRecord()
            mb2.read_MassBankRecord(ms2filename)

            mb2_step1 = ddatoolbox.massbankio.MassBankRecord()
            mb2_step1.read_MassBankRecord(step1file)

            temp_record = eval(mb2_step1.get_anything("COMMENT:", "COUNTER"))['temp_record']
            temp_recordkey = list(temp_record.keys())
            rawspecidlist = temp_record[temp_recordkey[0]]
            rawfilelist = []
            for rawspecid in rawspecidlist:
                filenametemp = re.sub(r'_[0-9]+$', '', rawspecid)
                rawfilelist.append(filenametemp)



            ionform_list = []


            #
            #raw mz
            #
            mz, mzstdev = mb2_step1.get_precmz()
            number = mb2.get_number()

            massbank_hits = mb2.get_MassBank()
            mf_hits = mb2.get_MF()

            accession = mb2.get_accession()
            polarity = mb2.get_polarity()
            #print(dataset, accession, mz, mzstdev)
            #thresholds = mb2.get_mzthresholds_averaged()


            correction_ms1 = paths["settings"]["ms1_mz_correction"]
            correction_ms2 = paths["settings"]["ms2_mz_correction"]
            mzthreshold_ms2 = paths["settings"]["ms2_mz_averaged_thres"]
            mzthreshold_ms1 = paths["settings"]["ms1_mz_averaged_thres"]

            thresholds ={}
            thresholds["ms2_mz_averaged_thres"] = float(mzthreshold_ms2)
            thresholds["ms1_mz_averaged_thres"] = float(mzthreshold_ms1)
            thresholds["ms2_mz_correction"] = correction_ms2
            thresholds["ms1_mz_correction"] = correction_ms1
            thresholds["unit"] = "ppm"


            #
            # Setting Up MassBank Server (Only one time)
            #
            if new_dict == "":

                #if polarity == "NEGATIVE":
                #    path = "./res/spectra/GNPSnegative/"
                #else:
                #    path = "./res/spectra/GNPSpositive/"


                new_dict = ddatoolbox.massbankserver.setupMassBank(path)

            #
            # Ionform data
            #
            ionformdata = ddatoolbox.massbankio.read_ionformlist(polarity = polarity)


            if number < minimunnumber:
                print("skipped", accession, number )
                break

            #print(accession, mf_hits)

            applied = []
            applied_deltamz = 10000000
            for massbank_hit, mf_hit in itertools.product(massbank_hits, mf_hits):
                print(massbank_hit, mf_hit)


                metname = massbank_hit

                recordid = massbank_hit.split(" ")[10]
                ionform = mf_hit.split(" ")[2]

                mfmf = mf_hit.split(" ")[0]

                mbmf = new_dict[recordid]["FORMULA"]
                diffmf = ddatoolbox.mfsearch.diffformula(mbmf, mfmf)

                if diffmf == "":
                    deltamz = abs(ddatoolbox.mfsearch.calcmw(mfmf, polarity = polarity)-mz)
                    if applied_deltamz > deltamz:

                        applied_deltamz = deltamz
                        applied = [massbank_hit, mf_hit]
            if applied == []:
                continue

            [massbank_hit, mf_hit] = applied
            metname = massbank_hit

            recordid = massbank_hit.split(" ")[10]
            ionform = mf_hit.split(" ")[2]

            mfmf = mf_hit.split(" ")[0]
            mbmf = new_dict[recordid]["FORMULA"]


            outputthisfragment = []
            counterall += 1

            print(dataset, accession, number, mz, mfmf, ddatoolbox.mfsearch.calcmw(mfmf, polarity = polarity)-mz, mbmf, massbank_hit)

            mz_corrected = ddatoolbox.mfsearch.mzcorrection(mz, correction_ms1)

            #exactprecmz = mfsearch.calcmw(mfsearch.ionform(mfmf, ionform), polarity = polarity)

            exactprecmz = ddatoolbox.mfsearch.calcmw(mfmf) + ionformdata[ionform]["mass"]

            print(mfmf,  ionform, ddatoolbox.mfsearch.ionform(mfmf, ionform))

            #
            # compariton of precmz
            #
            temp = [accession, metname, "p","averaged" , mfmf, str(mz_corrected),str(exactprecmz),str(number), "",str((mz_corrected - exactprecmz)/exactprecmz*1000000), "","",\
             "",str((mz - exactprecmz)/exactprecmz*1000000), "",""]
            print(temp)


            prec_ave.append(mz_corrected - exactprecmz)
            outputthisfragment.append(temp)
            mzlist =eval(mb2_step1.get_anything("COMMENT:", "PRECURSOR_M/Z_LIST"))




            for precmz, filenametemp in zip(mzlist, rawfilelist):
                if precmz == '': continue
                mz_corrected = ddatoolbox.mfsearch.mzcorrection(precmz, correction_ms1)

                temp = [accession, metname, "p",filenametemp, mfmf, str(mz_corrected),str(exactprecmz),str(number), str((mz_corrected - exactprecmz)/exactprecmz*1000000), "","",\
                 "",str((precmz - exactprecmz)/exactprecmz*1000000), "","",""]
                prec_raw.append(mz_corrected - exactprecmz)
                outputthisfragment.append(temp)


            result = ddatoolbox.mfsearch.elementnumber(mfmf)
            Climit = result["C"]
            Hlimit = result["H"]
            Nlimit = result["N"]
            Olimit = result["O"]
            Plimit = result["P"]
            Slimit = result["S"]

            intensity_threshold = 10

            temp_hash = mb2_step1.get_annotation_from_items()
            temp_fragmentmz = sorted(temp_hash.keys())

            allcounter = 0
            hitcounter = 0

            for fragmentmz in temp_fragmentmz:
                threshold_temp = ddatoolbox.massbankio.calcThresholdAveraged(fragmentmz, "ms2", thresholds)

                fragmentintensity = int(temp_hash[fragmentmz]["intensity"])
                fragmentmzstdev = float(temp_hash[fragmentmz]["mzsd"])
                fragmentnumber = float(temp_hash[fragmentmz]["number"])
                correlation = float(temp_hash[fragmentmz]["correlation"])

                if float(temp_hash[fragmentmz]["intensity"]) < intensity_threshold:
                    print("  ", "{:6f} {:3d} {:8f} {:8f}".format(fragmentmz,fragmentintensity,  fragmentmzstdev, correlation), "less than 1%")
                    continue
                if float(fragmentmz) > mz + threshold_temp * 2:
                    print("  ", "{:6f} {:3d} {:8f} {:8f}".format(fragmentmz,fragmentintensity,  fragmentmzstdev, correlation), "above prec")
                    continue
                if correlation < 0.6 :
                    print(" #", "{:6f} {:3d} {:8f} {:8f}".format(fragmentmz,fragmentintensity,  fragmentmzstdev, correlation), "Poor correl")
                    continue
                #
                # Ignore isotope-peak candidates (13C satellites of another fragment).
                #

                if ddatoolbox.mfsearch.check13C(fragmentmz, temp_hash, threshold = threshold_temp * 2):
                    print("  ", "{:6f} {:3d} {:8f} {:8f}".format(fragmentmz,fragmentintensity,  fragmentmzstdev, correlation),  "  ignored")

                    continue

                #
                # Counter of peaks
                #
                allcounter += 1

                #threshold_temp = massbankio.calcThreshold(fragmentmz, "ms2", thresholds)
                fragmentmz_corrected = ddatoolbox.mfsearch.mzcorrection(fragmentmz, correction_ms2)

                fragment_result, fragment_result_list = ddatoolbox.mfsearch.mfsearch(fragmentmz_corrected, threshold = float(mzthreshold_ms2), ihdtype = "non-integer", \
                Climit = Climit , Hlimit =Hlimit, Nlimit = Nlimit, Olimit = Olimit, Plimit = Plimit, Slimit = Slimit, polarity = polarity)

                difflist = ["{:8f}".format(fragment_result[mf]["delta"])  for mf in fragment_result_list]

                if len(fragment_result_list) > 0:
                    print("* ", "{:6f} {:3d} {:8f} {:8f}".format(fragmentmz,fragmentintensity,  fragmentmzstdev, correlation), fragment_result_list, difflist)
                    #
                    # Counter of peaks
                    #

                    hitcounter += 1
                    #
                    # compariton of precmz
                    #

                    for mftemp in fragment_result_list[0:1]:
                        exactprecmz = ddatoolbox.mfsearch.calcmw(mftemp, polarity = polarity)

                        mz_corrected =ddatoolbox.mfsearch.mzcorrection(fragmentmz, correction_ms2)
                        temp = [accession, metname, "f","averaged" , mfmf, str(mz_corrected),str(exactprecmz),str(number), "","", "",str((mz_corrected - exactprecmz)/exactprecmz*1000000),\
                         "","", "",str((fragmentmz - exactprecmz)/exactprecmz*1000000)]


                        frag_ave.append(mz_corrected - exactprecmz)

                        outputthisfragment.append(temp)

                        mzlist = eval(temp_hash[fragmentmz]["mzlist"])
                        #for precmz in mzlist:
                        for precmz, filenametemp in zip(mzlist, rawfilelist):

                            if precmz == '': continue
                            exactprecmz = ddatoolbox.mfsearch.calcmw(mftemp, polarity = polarity)
                            mz_corrected =ddatoolbox.mfsearch.mzcorrection(precmz, correction_ms2)

                            temp = [accession, metname, "f",filenametemp , mfmf, str(mz_corrected),str(exactprecmz),str(number),"", "",str((mz_corrected - exactprecmz)/exactprecmz*1000000), \
                            "","", "",str((precmz - exactprecmz)/exactprecmz*1000000), ""]

                            frag_raw.append(mz_corrected - exactprecmz)

                            outputthisfragment.append(temp)
                else:

                    print("--", fragmentmz, "{:3d} {:6f} {:6f}".format(fragmentintensity,  fragmentmzstdev, correlation), fragment_result_list, difflist)


                difflist = ["{:6f}".format(ddatoolbox.mfsearch.calcmw(mf, polarity = polarity) - float(fragmentmz))  for mf in fragment_result_list]
            print(hitcounter, allcounter)
            if allcounter == 0:
                continue
            if  hitcounter/ allcounter > 0.5 or hitcounter + 1 == allcounter:
                print("hit!")
                counterhit += 1
                outputdata.extend(outputthisfragment)

        rowtemp = []
        for row in outputdata:
            texttemp = '\t'.join(row)
            texttemp += "\n"
            rowtemp.append(texttemp)

        with open(outputfile, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, dialect='excel')
            writer.writerows(outputdata)

        #print(sum(prec_ave)/len(prec_ave),sum(prec_raw)/len(prec_raw), sum(frag_ave)/len(frag_ave),sum(frag_raw)/len(frag_raw))
        #print(stdev(prec_ave),stdev(prec_raw),stdev(frag_ave),stdev(frag_raw))

        print(counterhit, counterall)

        dt_now = datetime.datetime.now()
        ddatoolbox.ddatoolbox.writeLog(logfilename, ["Calibration Finished", accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), counterhit, counterall])

if __name__ == '__main__':
    args = sys.argv
    if len(args) >= 2:
        datasets = args[1:]
    else:
        datasets = ["HSA001"]
    processing(datasets)


