#-------------------------------------------------------------------------------
# Name:        ddatoolbox
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     07/10/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     MIT lisence
#-------------------------------------------------------------------------------
"""
ddatoolbox.py

Core module for Data-Dependent Acquisition (DDA) spectral processing.

This module provides functions for:
- Spectral similarity calculation (dot product)
- Decoding of compressed spectral data (Base64/Base32)
- Construction of MS/MS spectral cliques
- MS1 spectrum aggregation and averaging
- Ion clustering and peak grouping

The implemented algorithms are designed to:
- Handle large-scale LC-MS/MS datasets
- Identify structurally related spectra via similarity networks
- Improve signal quality through statistical aggregation

Dependencies
------------
- numpy
- numba
- scipy (optional)

Notes
-----
- Spectra are represented as dictionaries: {m/z: intensity}
- Intensities are often normalized before similarity calculations
- Threshold parameters control matching tolerance and filtering behavior

Author
------
Fumio Matsuda

Created
-------
2022

License
-------
MIT License

History
-------
240907 New folder structure

"""
import xml.etree.ElementTree as ET
import base64, struct,itertools, glob, re, shutil, time
import numpy as np
import zlib
from statistics import stdev
import os, csv, math, pickle, datetime
from scipy.stats import spearmanr
from . import massbankio

from numba import jit


#OVER_SIZE_LIMIT = 200_000_000
#csv.field_size_limit(OVER_SIZE_LIMIT)


#from hivetools.access import HiveFile, ValueRange
# hivetools 0.1.5
#from hivetools import HiveFile, ValueRange
# hivetools 0.1.8
from hive_data_access import HiveFile, ValueRange

def dot_product(mshash1, mshash2, thresholds):
    """Compute the normalized dot-product similarity between two MS/MS spectra.

    Short-circuits to ``0.0`` if the two spectra's precursor m/z or retention
    time differ by more than the configured thresholds. Otherwise, greedily
    matches each peak of ``mshash1`` (in descending intensity order) to the
    closest not-yet-used peak of ``mshash2`` within the m/z threshold, and sums
    the normalized intensity products of matched peaks.

    Args:
        mshash1 (dict): First spectrum hash, with keys ``"spectrum"`` (dict of
            m/z to intensity), ``"norm"`` (float, L2 norm used for
            normalization), ``"precursor_mz"`` (float), and ``"retentionTime"``
            (float or str).
        mshash2 (dict): Second spectrum hash, same structure as ``mshash1``.
        thresholds (dict): Dictionary containing ``"retentiontime_thres"``
            (float) and mass tolerances consumed via
            ``massbankio.calcThreshold``.

    Returns:
        float: Dot-product similarity score (0-1). Returns ``0.0`` if the
            precursor m/z or retention time mismatch exceeds the configured
            thresholds.
    """

    retentiontime_thres = thresholds["retentiontime_thres"]
    dotproduct_thres = thresholds["dotproduct_thres"]

    norm1 = mshash1['norm']
    norm2 = mshash2['norm']
    prec1 = mshash1['precursor_mz']
    prec2 = mshash2['precursor_mz']
    rt1 = float(mshash1['retentionTime'])
    rt2 = float(mshash2['retentionTime'])
    #
    # Check precursor m/z gap
    #
    threshold_temp = massbankio.calcThreshold(max(prec1, prec2), "ms1", thresholds) * 2
    if abs(prec1 - prec2) > threshold_temp:
        return(0.0)
    #
    # Check retention time gap
    #
    if abs(rt1 - rt2) > retentiontime_thres:
        return(0.0)

    mzlist1 = sorted(mshash1['spectrum'].keys(), key = lambda x: mshash1['spectrum'][x], reverse = True)
    mzlist2 = sorted(mshash2['spectrum'].keys(), key = lambda x: mshash2['spectrum'][x], reverse = True)

    dotproduct = 0.0
    used2 = {}
    for precmz1 in mzlist1:
        threshold_temp = massbankio.calcThreshold(precmz1, "ms2", thresholds) * 2
        for precmz2 in mzlist2:
            if precmz2 in used2:
                continue
            if abs(float(precmz1) - float(precmz2)) > threshold_temp:
                continue
            dotproduct = dotproduct + (mshash1['spectrum'][precmz1] / norm1) * (mshash2['spectrum'][precmz2] / norm2)
            used2[precmz2] = 1
            break
    return(dotproduct)

def decodeBase64(base64text, compressionType):
    """Decode a Base64-encoded (optionally zlib-compressed) 64-bit-float peak list.

    Adapted from
    https://groups.google.com/g/spctools-discuss/c/qK_QThoEzeQ. Decodes the
    Base64 payload into alternating m/z/intensity 64-bit float pairs, then drops
    peaks whose intensity is below 0.5% (5/1000) of the base peak.

    Args:
        base64text (str): Base64-encoded peak list (as found in mzML/mzXML files).
        compressionType (str): ``"zlib"`` to zlib-decompress the decoded bytes
            before unpacking; any other value skips decompression.

    Returns:
        tuple: A tuple of two elements.

            - norm (float): L2 norm of the retained (above-threshold) intensities.
            - spectra (dict): Mapping of m/z (float) to intensity (float) for the
              retained peaks.
    """

    decoded = base64.b64decode(base64text)
    if compressionType == "zlib":
        decoded = zlib.decompress(decoded)
    tmp_size = len(decoded)/8
    unpack_format1 = ">%dQ" % tmp_size

    idx = 0
    mz_list = []
    intensity_list = []

    for tmp in struct.unpack(unpack_format1,decoded):
        tmp_i = struct.pack("Q",tmp)
        tmp_f = struct.unpack("d",tmp_i)[0]
        if( idx % 2 == 0 ):
            mz_list.append( float(tmp_f) )
        else:
            intensity_list.append( float(tmp_f) )
        idx += 1
    basepeakintensity = max(intensity_list)
    mz_list_selected = []
    intensity_list_selected = []
    for mz, intensity in zip(mz_list, intensity_list):
        if intensity/ basepeakintensity * 1000 < 5:
            continue
        mz_list_selected.append(mz)
        intensity_list_selected.append(intensity)


    norm = np.linalg.norm(intensity_list_selected, ord=2)
    return(norm, dict(zip(mz_list_selected, intensity_list_selected)))


def decodeBase32(base64text, compressionType):
    """Decode a Base64-encoded (optionally zlib-compressed) 32-bit-float peak list.

    Adapted from
    https://groups.google.com/g/spctools-discuss/c/qK_QThoEzeQ. Same as
    :func:`decodeBase64` but for the 32-bit-float peak encoding (despite the
    function name, the wire encoding is still Base64; "32" refers to the float
    width). Decodes alternating m/z/intensity 32-bit float pairs, then drops
    peaks whose intensity is below 0.5% (5/1000) of the base peak.

    Args:
        base64text (str): Base64-encoded peak list (as found in mzML/mzXML files).
        compressionType (str): ``"zlib"`` to zlib-decompress the decoded bytes
            before unpacking; any other value skips decompression.

    Returns:
        tuple: A tuple of two elements.

            - norm (float): L2 norm of the retained (above-threshold) intensities.
            - spectra (dict): Mapping of m/z (float) to intensity (float) for the
              retained peaks.
    """
    decoded = base64.b64decode(base64text)
    if compressionType == "zlib":
        decoded = zlib.decompress(decoded)
    tmp_size = len(decoded)/4
    unpack_format1 = ">%dL" % tmp_size

    idx = 0
    mz_list = []
    intensity_list = []

    for tmp in struct.unpack(unpack_format1,decoded):

        tmp_i = struct.pack("I",tmp)
        tmp_f = struct.unpack("f",tmp_i)[0]
        if( idx % 2 == 0 ):
            mz_list.append( float(tmp_f) )
        else:
            intensity_list.append( float(tmp_f) )
        idx += 1
    basepeakintensity = max(intensity_list)
    mz_list_selected = []
    intensity_list_selected = []
    for mz, intensity in zip(mz_list, intensity_list):
        if intensity/ basepeakintensity * 1000 < 5:
            continue
        mz_list_selected.append(mz)
        intensity_list_selected.append(intensity)


    norm = np.linalg.norm(intensity_list_selected, ord=2)
    return(norm, dict(zip(mz_list_selected, intensity_list_selected)))

def writeLog(filename, contentslist, mode='a'):
    """Write one tab-separated log entry to a file.

    Args:
        filename (str): Path to the log file.
        contentslist (list): Elements to be stringified and joined with tabs
            into a single log line.
        mode (str, optional): File open mode, ``'a'`` to append or ``'w'`` to
            overwrite. Defaults to ``'a'``.

    Returns:
        None
    """

    f = open(filename, mode,  encoding='UTF-8')
    dt_now = datetime.datetime.now()

    #strlist = [dt_now.strftime('%Y.%m.%d')]+[str(x) for x in contentslist]
    strlist = [str(x) for x in contentslist]
    strtemp = "\t".join(strlist)
    #print(strtemp)
    f.write(strtemp+"\n")
    f.close()

def DDAclique_large(msmsdata, thresholds, filenames, report="off"):
    """Group MS/MS spectra into similarity cliques using a fast m/z-indexed heuristic.

    Designed for large datasets (>1000 spectra): rather than an all-pairs
    comparison, sorts all spectra by precursor m/z once, then processes them in
    descending base-peak-intensity order. For each not-yet-used spectrum, walks
    outward (both directions) through its m/z-sorted neighbors within the
    precursor matching threshold, adding any neighbor whose dot-product
    similarity (:func:`dot_product`) exceeds ``dotproduct_thres`` to the clique
    and marking it used, until the m/z window is exhausted.

    Args:
        msmsdata (dict): Mapping of spectrum ID to spectrum hash (see
            :func:`dot_product` for the expected shape).
        thresholds (dict): Threshold parameters (m/z tolerance, dot product,
            retention time), as consumed by :func:`dot_product`.
        filenames (dict): File paths (unused directly in this function; kept for
            interface consistency with sibling clique-construction functions).
        report (str, optional): ``"on"`` to print progress every 100 spectra
            processed. Defaults to ``"off"``.

    Returns:
        list[list]: List of cliques, each a list of spectrum IDs whose spectra
            were judged similar enough to average together.
    """
    allspecs = []
    msmslist = msmsdata.keys()
    msmsnumber = len(msmslist)


    retentiontime_thres = thresholds["retentiontime_thres"]
    dotproduct_thres = thresholds["dotproduct_thres"]


    #
    # Retrieve all peaks
    #
    for msms1 in msmslist:
        msms1_prec = float(msmsdata[msms1]['precursor_mz'])
        intensity = float(msmsdata[msms1]['basePeakIntensity'])

        allspecs.append([msms1, msms1_prec, intensity])
    #
    # Sorted [id, filename_scan, m/z, intensity, order in m/z, order in intensity (Reversed),]
    #
    allpeaks_temp = sorted(allspecs, key = lambda x:x[1])
    allpeaks_temp = [allpeaks_temp[x]+[x] for x in range(len(allpeaks_temp))]
    allpeaks_temp = sorted(allpeaks_temp, key = lambda x:x[2], reverse = True)
    allpeaks_temp = [allpeaks_temp[x]+[x] for x in range(len(allpeaks_temp))]

    allpeaks_mzsorted = sorted(allpeaks_temp, key = lambda x:x[1])
    allpeaks_intsorted = sorted(allpeaks_temp, key = lambda x:x[2], reverse = True)

    numberoftotalpeaks = len(allpeaks_mzsorted)

    cliques = []
    usedpeak = {}
    for i, keypeak in enumerate(allpeaks_intsorted):
        if report == "on" and i % 100 == 0:
            print(i, "/", msmsnumber)
        #
        # If this peak has already used in other clique
        #
        keyid = keypeak[0]
        if keyid in usedpeak:
            continue
        usedpeak[keyid] = i

        keymz = keypeak[1]
        threshold_temp = massbankio.calcThreshold(keymz, "ms2", thresholds)
        lowerboundary = keymz - threshold_temp
        upperboundary = keymz + threshold_temp
        #
        # Position in allpeaks_mzsorted
        #
        posinmz = keypeak[3]
        #
        # Preparing clique seed
        #
        clique_temp = [keypeak[0]]
        #
        # Searching to downward
        #
        step = 1

        while True:
            if (posinmz - step) < 0:
                break
            checkpeak = allpeaks_mzsorted[posinmz - step]

            chackmz = checkpeak[1]
            if chackmz < lowerboundary:
                break

            step +=1
            msms1 = checkpeak[0]
            if msms1 in usedpeak:
                continue
            dotproduct = dot_product(msmsdata[keyid], msmsdata[msms1], thresholds)
            if dotproduct < dotproduct_thres:
                continue
            clique_temp.append(msms1)
            usedpeak[msms1] = i



            #for msms2 in clique_temp:
            #    dotproduct = dot_product(msmsdata[msms1], msmsdata[msms2], thresholds)
            #    if dotproduct < 0.9:
            #        break
            #else:
            #    clique_temp.append(msms1)
            #    usedpeak[msms1] = i
        #
        # Searching to upward
        #
        step = 1

        while True:
            if (posinmz + step) >= numberoftotalpeaks:
                break
            checkpeak = allpeaks_mzsorted[posinmz + step]


            chackmz = checkpeak[1]
            if chackmz > upperboundary:
                break
            step +=1
            msms1 = checkpeak[0]
            if msms1 in usedpeak:
                continue

            dotproduct = dot_product(msmsdata[keyid], msmsdata[msms1], thresholds)
            if dotproduct < dotproduct_thres:
                continue
            clique_temp.append(msms1)
            usedpeak[msms1] = i

            #for msms2 in clique_temp:
            #    dotproduct = dot_product(msmsdata[msms1], msmsdata[msms2], thresholds)
            #    if dotproduct < 0.9:
            #        break
            #else:
            #    clique_temp.append(msms1)
            #    usedpeak[msms1] = 1

        cliques.append(clique_temp)
    return(cliques)


def DDAclique_small(msmsdata, thresholds, filenames, report = "off"):
    """Group MS/MS spectra into similarity cliques using an all-cliques-comparison heuristic.

    Intended for smaller datasets (< 1000 spectra), where the m/z-indexed
    approach of :func:`DDAclique_large` isn't needed. Processes spectra in
    descending base-peak-intensity order; for each spectrum, compares it against
    every existing clique's members (in precursor-m/z order, stopping early on
    the first mismatch or low dot product), and joins the best-scoring
    compatible clique (ranked by clique size plus average dot product), or
    starts a new clique if none qualify.

    Args:
        msmsdata (dict): Mapping of spectrum ID to spectrum hash (see
            :func:`dot_product` for the expected shape).
        thresholds (dict): Threshold parameters (m/z tolerance, dot product,
            retention time), as consumed by :func:`dot_product`.
        filenames (dict): File paths (unused directly in this function; kept for
            interface consistency with sibling clique-construction functions).
        report (str, optional): ``"on"`` to print progress every 100 spectra
            processed. Defaults to ``"off"``.

    Returns:
        list[list]: List of cliques, each a list of spectrum IDs whose spectra
            were judged similar enough to average together.
    """
    #
    #
    mzmzlist = msmsdata.keys()
    cliques = []
    msmslist = sorted(mzmzlist, key=lambda x: msmsdata[x]['basePeakIntensity'], reverse = True)
    msmsnumber = len(msmslist)
    counter = 0

    retentiontime_thres = thresholds["retentiontime_thres"]
    dotproduct_thres = thresholds["dotproduct_thres"]


    for msms1 in msmslist:
        if report == "on" and counter % 100 == 0:
            print(counter, "/", msmsnumber)
        counter += 1
        hitclique = []
        hitclique_averagescore = []
        prec_msms1 = msmsdata[msms1]['precursor_mz']




        for clique in cliques:
            averagescore_temp = []
            for msms2 in clique:
                prec_msms2 = msmsdata[msms2]['precursor_mz']

                threshold_temp = massbankio.calcThreshold(prec_msms1, "ms1", thresholds) * 2

                if abs(prec_msms1 - prec_msms2) > threshold_temp:
                    break

                dotproduct = dot_product(msmsdata[msms1], msmsdata[msms2], thresholds)
                averagescore_temp.append(dotproduct)
                if not dotproduct > dotproduct_thres:
                    break
            else:
                hitclique.append(clique)
                hitclique_averagescore.append(np.mean(averagescore_temp))
        if hitclique == []:
            cliques.append([msms1])
            continue
        hitclique_order_sorted = sorted(range(len(hitclique_averagescore)), key=lambda x: len(hitclique[x]) + hitclique_averagescore[x], reverse = True)
        clique = hitclique[hitclique_order_sorted[0]]
        clique.append(msms1)

    return(cliques)

def constructonofDDAclique(msmsdata, thresholds, filenames, report = "off"):
    """Dispatch to :func:`DDAclique_small`/:func:`DDAclique_large` and filter by minimum size.

    Chooses the appropriate clique-construction algorithm based on dataset size
    (``DDAclique_small`` when ``len(msmsdata) < 1000``, ``DDAclique_large``
    otherwise), then sorts the resulting cliques by descending size and drops
    any smaller than ``thresholds["size_minimum"]``.

    Args:
        msmsdata (dict): Mapping of spectrum ID to spectrum hash (see
            :func:`dot_product` for the expected shape).
        thresholds (dict): Threshold parameters, including
            ``"size_minimum"`` (int, minimum clique size to keep), plus the m/z
            tolerance/dot-product/retention-time thresholds consumed by
            :func:`dot_product`.
        filenames (dict): File paths, including ``"logfilename"`` (currently
            read but not used by this function directly).
        report (str, optional): ``"on"`` to print progress during clique
            construction. Defaults to ``"off"``.

    Returns:
        list[list]: Cliques with at least ``thresholds["size_minimum"]``
            members, sorted by descending size.
    """

    minimumsize = thresholds["size_minimum"]
    logfile = filenames["logfilename"]

    specnumber = len(msmsdata)
    if specnumber < 1000:
        cliques = DDAclique_small(msmsdata, thresholds, filenames, report)
    else:
        cliques = DDAclique_large(msmsdata, thresholds, filenames, report)

    cliquelist_sorted = sorted(cliques, key=lambda x: len(x), reverse = True)

    #
    # Selection of cluques by size (n)
    #
    cliquelist_selected = []
    for clique in cliquelist_sorted :
        if len(clique) < minimumsize:
            continue
        cliquelist_selected.append(clique)


    return(cliquelist_selected)


@jit(nopython=True, cache=True)
def findIonclique(numberofspectra, numberoftotalpeaks, allmz, allint,allid, allposition, allspec, allpeaks_mzsorted, allpeaks_intsorted, mzthreshold):
    """Group individual peaks (across many spectra) into m/z-consistent ion cliques.

    Numba-accelerated core of the MS1 ion-clique construction pipeline. Scans
    all peaks (across every spectrum being averaged) in descending intensity
    order; for each not-yet-used peak, walks outward through its m/z-sorted
    neighbors within ``mzthreshold`` (ppm), collecting candidate peaks. Within
    that m/z window, keeps only the single most intense peak per originating
    spectrum (so each clique has at most one peak per spectrum), marks all
    selected peaks used, and records the clique.

    Args:
        numberofspectra (int): Total number of spectra contributing peaks.
        numberoftotalpeaks (int): Total number of individual peaks across all
            spectra.
        allmz (ndarray): m/z value of each peak, indexed by peak ID.
        allint (ndarray): Intensity of each peak, indexed by peak ID.
        allid (ndarray): (Unused inside the function body; kept for interface
            compatibility with the caller's array layout.)
        allposition (ndarray): Index of each peak's position within
            ``allpeaks_mzsorted``, indexed by peak ID.
        allspec (ndarray): Originating spectrum index of each peak, indexed by
            peak ID.
        allpeaks_mzsorted (ndarray): Peak IDs sorted by ascending m/z.
        allpeaks_intsorted (ndarray): Peak IDs sorted by descending intensity.
        mzthreshold (float): m/z matching tolerance, in ppm.

    Returns:
        list[ndarray]: List of ion cliques; each clique is an array of peak IDs
            (at most one per originating spectrum) judged to represent the same
            ion.

    Note:
        Numba-jitted (``nopython=True``); array dtypes must stay as
        ``np.int64``/compatible numeric types for Numba's type inference to
        succeed.
    """
    cliques = []
    usedpeak = np.zeros(numberoftotalpeaks, dtype=np.int64) # np.int64 is the reliable choice under Numba

    threshold_numberofspecinoneclique = int(numberofspectra * 0.7)

    for i in range(numberoftotalpeaks):
        keyid = allpeaks_intsorted[i]
        #
        # If this peak has already used in other clique
        #
        if usedpeak[keyid] == 1:
            continue

        keymz = allmz[keyid]
        threshold_temp = keymz * mzthreshold/1000000
        lowerboundary = keymz - threshold_temp
        upperboundary = keymz + threshold_temp
        #
        # Position in allpeaks_mzsorted
        #
        posinmz = allposition[keyid]
        #
        # Preparing clique seed
        #
        clique_temp = [keyid]

        #
        # Searching to downward
        #
        step = 1

        while True:
            if (posinmz - step) < 0:
                break
            checkpeakid = allpeaks_mzsorted[posinmz - step]

            chackmz = allmz[checkpeakid] # 1: mz

            if chackmz < lowerboundary:
                break
            clique_temp.append(checkpeakid)
            #clique_temp = np.append(clique_temp,checkpeakid)
            step +=1
        #
        # Searching to upward
        #
        step = 1

        while True:
            if (posinmz + step) >= numberoftotalpeaks:
                break
            checkpeakid = allpeaks_mzsorted[posinmz + step]

            chackmz = allmz[checkpeakid] # 1: mz

            if chackmz > upperboundary:
                break
            clique_temp.append(checkpeakid)
            #clique_temp = np.append(clique_temp,checkpeakid)
            step +=1

        #if len(clique_temp) < threshold_numberofspecinoneclique:
        #    continue

        usedspecint = np.zeros(numberofspectra)

        usedspec = np.zeros(numberofspectra, dtype=np.int64) # np.int64 is the reliable choice under Numba

        for keyid in clique_temp:

            if usedpeak[keyid] == 1:
                continue


            specnum = int(allspec[keyid])
            thisint = allint[keyid]
            if usedspecint[specnum] < thisint:
                usedspecint[specnum] = thisint
                usedspec[specnum] = keyid
        #print(usedspec)

        clique = usedspec[np.nonzero(usedspec)]
        #print(clique)
        for keyid in clique:
            usedpeak[keyid] = 1

        #for keyid in usedspec:
        #    if not keyid == 0:
        #        clique.append(keyid)
        #        usedpeak[keyid] = 1
        #        continue

        cliques.append(clique)

    return(cliques)

def constructionofMS1Ionclique(msmslist, msmsdata, msdata, thresholds, filenames):
    """Average the MS1 (precursor) spectra of a DDA clique into one consensus MS1 spectrum.

    Flattens every raw MS1 peak across ``msmslist``'s member spectra into a
    single peak pool, groups them into ion cliques via :func:`findIonclique`
    (using the raw MS1 m/z threshold), keeps only ion cliques present in more
    than 70% of the member spectra, and for each surviving ion clique computes
    the mean (or, for >30 contributing peaks, median) m/z and the mean
    base-peak-normalized intensity. The resulting averaged spectrum is finally
    rescaled to a 0-999 relative-intensity range against its own base peak.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique, used to determine which MS1 spectra to average.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash, used for
            base-peak-intensity normalization.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash (with a ``"spectrum"`` sub-dict of m/z to intensity).
        thresholds (dict): Threshold parameters; must include
            ``"ms1_mz_raw_thres"`` (ppm).
        filenames (dict): File paths, including ``"logfilename"`` (currently
            read but not used by this function directly).

    Returns:
        tuple: A tuple of four elements.

            - averaged_spectra_hash (dict): Averaged MS1 spectrum, mapping mean
              m/z to relative intensity (0-999 scale).
            - stdev_hash (dict): Per-averaged-peak statistics, mapping mean m/z
              to a dict with ``"intensity"``, ``"number"`` (contributing peak
              count), ``"intlist"``/``"mzlist"`` (per-source-spectrum raw
              values, ``""`` where absent), and ``"annotation"`` (placeholder
              ``"?"``).
            - totalpeak (int): Total number of raw MS1 peaks considered.
            - usedpeak (int): Number of raw peaks retained in the averaged
              spectrum (``totalpeak`` minus peaks in cliques below the 70%
              presence threshold).

    Note:
        Returns ``({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)``
        if no ion cliques survive the presence-fraction filter, or if the
        resulting spectrum's base-peak intensity is ``0``.
    """

    logfile = filenames["logfilename"]
    ms1threshold = thresholds["ms1_mz_raw_thres"]


    allpeaks = {}
    allpeaks_list = []

    allmz = []
    allint = []
    allposition = []
    allid = []
    allspec = []
    #
    # Retrieve all peaks
    #
    counter = 0
    speccounter = 0
    for msms1 in msmslist:
        msms1_precs = list(msmsdata[msms1]['spectrum'].keys())
        for msms1_prec in msms1_precs:
            intensity = msmsdata[msms1]['spectrum'][float(msms1_prec)]

            allpeaks[counter] = [msms1, float(msms1_prec), intensity]

            allpeaks_list.append(counter)
            allmz.append(float(msms1_prec))
            allint.append(float(intensity))
            allposition.append(0)
            allid.append(msms1)
            allspec.append(speccounter)
            counter +=1
        speccounter += 1
    #
    # Sorted [id, filename_scan, m/z, intensity, order in m/z, order in intensity (Reversed),]
    #
    allpeaks_intsorted = sorted(allpeaks_list, key = lambda x:allint[x], reverse = True)
    allpeaks_mzsorted = sorted(allpeaks_list, key = lambda x:allmz[x])

    for i, targetid in enumerate(allpeaks_mzsorted):
        allposition[targetid] = i

    numberoftotalpeaks = len(allpeaks_mzsorted)

    allmz = np.array(allmz)
    allint = np.array(allint)
    allid = np.array(allid)
    allposition = np.array(allposition)
    allpeaks_mzsorted = np.array(allpeaks_mzsorted)
    allpeaks_intsorted = np.array(allpeaks_intsorted)
    allspec = np.array(allspec)
    numberofspectra = len(msmslist)

    cliques = findIonclique(numberofspectra, numberoftotalpeaks, allmz, allint,allid, allposition, allspec, allpeaks_mzsorted, allpeaks_intsorted, ms1threshold)

    numberofunusedpeaks = 0

    msmslist_number = len(msmslist)
    temp_hash = {}
    stdev_hash = {}

    for msmsclique in cliques:
        msmsclique_number = len(msmsclique)

        if msmsclique_number/msmslist_number <= 0.7:
            numberofunusedpeaks += msmsclique_number
            continue
        fragmentmzlist = []
        fragmentintensitylist = []

        intensity_hash_temp = {}
        mz_hash_temp = {}

        for keyid in msmsclique:

            nodeid = allid[keyid] # 0: id in msmsdata
            fragmentmz = allmz[keyid] # 1: mz
            intensity = allint[keyid] # 1: intensity

            fragmentmzlist.append(fragmentmz)

            basePeakIntensity = float(msmsdata[nodeid]['basePeakIntensity'])

            fragmentintensitylist.append(float(intensity)/basePeakIntensity)

            intensity_hash_temp[nodeid] = intensity
            mz_hash_temp[nodeid] = fragmentmz
        #
        # list of intensity
        #
        intensity_alllist = []
        mz_alllist = []
        for msmsid in msmslist:
            if msmsid in intensity_hash_temp:
                intensity_alllist.append(int(intensity_hash_temp[msmsid]))
                mz_alllist.append(float(mz_hash_temp[msmsid]))
            else:
                intensity_alllist.append("")
                mz_alllist.append("")



        if len(fragmentmzlist) > 30:
            mean_mz = np.median(fragmentmzlist)
        else:
            mean_mz = np.mean(fragmentmzlist)


        mean_intensity = sum(fragmentintensitylist)/len(fragmentintensitylist) * 999
        temp_hash[mean_mz] = mean_intensity
        stdev_hash[mean_mz] = {"intensity": mean_intensity, "number": len(fragmentmzlist), \
        "intlist": intensity_alllist, "mzlist": mz_alllist, "annotation": "?"}


    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)

    #
    # Searching Base Peak ion for MS1
    #
    basepeakintensity = max([x for x in temp_hash.values()])

    if basepeakintensity == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)


    temp_hash2 = {}
    for mz, intensity in temp_hash.items():
        relint= intensity/basepeakintensity * 999

        temp_hash2[mz] = relint
        stdev_hash[mz]["intensity"] = relint

    temp_hash_for_output = {}
    temp_stdevhash_for_output = {}
    mzlist = sorted(temp_hash2.keys())
    for mz in mzlist:
        temp_hash_for_output[mz] = temp_hash2[mz]
        temp_stdevhash_for_output[mz] = stdev_hash[mz]

    return(temp_hash_for_output, temp_stdevhash_for_output, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)



def constructionofMS1Ionclique_old(msmslist, msmsdata, msdata, thresholds, filenames):
    """Legacy (pure-Python, non-Numba) implementation of :func:`constructionofMS1Ionclique`.

    .. deprecated::
        Superseded by :func:`constructionofMS1Ionclique`, which uses the
        Numba-jitted :func:`findIonclique` for ion-clique grouping instead of a
        pure-Python m/z-window walk, making it substantially faster on large
        datasets. Kept for reference.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash.
        thresholds (dict): Threshold parameters; must include
            ``"ms1_mz_raw_thres"`` (ppm).
        filenames (dict): File paths, including ``"logfilename"``.

    Returns:
        tuple: Same 4-element shape as :func:`constructionofMS1Ionclique`'s
            return value: ``(averaged_spectra_hash, stdev_hash, totalpeak, usedpeak)``.

    See Also:
        :func:`constructionofMS1Ionclique`: current implementation.
    """

    logfile = filenames["logfilename"]

    allpeaks = {}
    allpeaks_list = []
    #
    # Retrieve all peaks
    #
    counter = 0
    for msms1 in msmslist:
        msms1_precs = list(msmsdata[msms1]['spectrum'].keys())
        for msms1_prec in msms1_precs:
            intensity = msmsdata[msms1]['spectrum'][float(msms1_prec)]

            allpeaks[counter] = [msms1, float(msms1_prec), intensity]
            allpeaks_list.append(counter)
            counter +=1
    #
    # Sorted [id, filename_scan, m/z, intensity, order in m/z, order in intensity (Reversed),]
    #
    allpeaks_intsorted = sorted(allpeaks_list, key = lambda x:allpeaks[x][2], reverse = True)
    #
    # Maybe sort() is faster and memory efficient than sorted()
    #allpeaks_mzsorted = sorted(allpeaks_list, key = lambda x:allpeaks[x][1])
    allpeaks_list.sort(key = lambda x:allpeaks[x][1])
    allpeaks_mzsorted = allpeaks_list

    for i, targetid in enumerate(allpeaks_mzsorted):
        allpeaks[targetid].append(i)

    numberoftotalpeaks = len(allpeaks_mzsorted)


    cliques = []
    usedpeak = set()
    for i, keyid in enumerate(allpeaks_intsorted):
        #
        # If this peak has already used in other clique
        #
        if keyid in usedpeak:
            continue

        keymz = allpeaks[keyid][1] # 1: mz
        threshold_temp = massbankio.calcThreshold(keymz, "ms1", thresholds)
        lowerboundary = keymz - threshold_temp
        upperboundary = keymz + threshold_temp
        #
        # Position in allpeaks_mzsorted
        #
        posinmz = allpeaks[keyid][3] # 3: pos in m/z list
        #
        # Preparing clique seed
        #
        clique_temp = [keyid]
        clique = []
        #
        # Searching to downward
        #
        step = 1

        while True:
            if (posinmz - step) < 0:
                break
            checkpeakid = allpeaks_mzsorted[posinmz - step]

            chackmz = allpeaks[checkpeakid][1] # 1: mz

            if chackmz < lowerboundary:
                break
            clique_temp.append(checkpeakid)
            step +=1
        #
        # Searching to upward
        #
        step = 1

        while True:
            if (posinmz + step) >= numberoftotalpeaks:
                break
            checkpeakid = allpeaks_mzsorted[posinmz + step]
            chackmz = allpeaks[checkpeakid][1] # 1: mz
            if chackmz > upperboundary:
                break
            clique_temp.append(checkpeakid)
            step +=1
        #
        # Sorting clique_temp with reved int order
        #
        #clique_temp_sorted = sorted(clique_temp, key = lambda x:allpeaks[x][2], reverse = True) # 2: intensity
        clique_temp.sort(key = lambda x:allpeaks[x][2], reverse = True) # 2: intensity
        #
        # Check duplication and
        # for keyid in clique_temp_sorted:
        usedfilenames = set()
        for keyid in clique_temp:
            if keyid in usedpeak:
                continue
            keyfilename = allpeaks[keyid][0] # 0: id in msmsdata
            if keyfilename in usedfilenames:
                continue
            usedpeak.add(keyid)
            usedfilenames.add(keyfilename)
            clique.append(keyid)

        cliques.append(clique)

    print(len(cliques), cliques[0])


    numberofunusedpeaks = 0

    msmslist_number = len(msmslist)
    temp_hash = {}
    stdev_hash = {}

    for msmsclique in cliques:
        msmsclique_number = len(msmsclique)

        if msmsclique_number/msmslist_number <= 0.7:
            numberofunusedpeaks += msmsclique_number
            continue
        fragmentmzlist = []
        fragmentintensitylist = []

        intensity_hash_temp = {}
        mz_hash_temp = {}

        for keyid in msmsclique:

            nodeid = allpeaks[keyid][0] # 0: id in msmsdata
            fragmentmz = allpeaks[keyid][1] # 1: mz
            intensity = allpeaks[keyid][2] # 1: intensity

            fragmentmzlist.append(fragmentmz)

            intensity = msmsdata[nodeid]['spectrum'][float(fragmentmz)]
            basePeakIntensity = float(msmsdata[nodeid]['basePeakIntensity'])

            fragmentintensitylist.append(float(intensity)/basePeakIntensity)

            intensity_hash_temp[nodeid] = intensity
            mz_hash_temp[nodeid] = fragmentmz
        #
        # list of intensity
        #
        intensity_alllist = []
        mz_alllist = []
        for msmsid in msmslist:
            if msmsid in intensity_hash_temp:
                intensity_alllist.append(int(intensity_hash_temp[msmsid]))
                mz_alllist.append(float(mz_hash_temp[msmsid]))
            else:
                intensity_alllist.append("")
                mz_alllist.append("")



        if len(fragmentmzlist) > 30:
            mean_mz = np.median(fragmentmzlist)
        else:
            mean_mz = np.mean(fragmentmzlist)


        mean_intensity = sum(fragmentintensitylist)/len(fragmentintensitylist) * 999
        temp_hash[mean_mz] = mean_intensity
        stdev_hash[mean_mz] = {"intensity": mean_intensity, "number": len(fragmentmzlist), \
        "intlist": intensity_alllist, "mzlist": mz_alllist, "annotation": "?"}


    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)

    #
    # Searching Base Peak ion for MS1
    #
    basepeakintensity = max([x for x in temp_hash.values()])

    #for mz, intensity in temp_hash.items():
    #    if intensity > basepeakintensity:
    #        basepeakintensity = intensity



    if basepeakintensity == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)


    temp_hash2 = {}
    for mz, intensity in temp_hash.items():
        relint= intensity/basepeakintensity * 999

        temp_hash2[mz] = relint
        stdev_hash[mz]["intensity"] = relint

    temp_hash_for_output = {}
    temp_stdevhash_for_output = {}
    mzlist = sorted(temp_hash2.keys())
    for mz in mzlist:
        temp_hash_for_output[mz] = temp_hash2[mz]
        temp_stdevhash_for_output[mz] = stdev_hash[mz]

    return(temp_hash_for_output, temp_stdevhash_for_output, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)


def constructionofIonclique(msmslist, msmsdata, msdata, thresholds,  msx, filenames):
    """Dispatch to the appropriate ion-clique averaging routine for MS1 or MS2 data.

    Routes to :func:`constructionofMS1Ionclique` for ``msx == "ms1"`` (precursor
    spectra) or to :func:`Ionclique_large` for ``msx == "ms2"`` (fragment
    spectra), unifying both under one call signature for callers such as
    :func:`generationMassBankrecord`.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash.
        thresholds (dict): Threshold parameters (``"ms1_mz_raw_thres"``/
            ``"ms2_mz_raw_thres"`` as appropriate for ``msx``).
        msx (str): ``"ms1"`` to average precursor spectra, ``"ms2"`` to average
            fragment spectra.
        filenames (dict): File paths, including ``"logfilename"``.

    Returns:
        tuple: A 4-element tuple ``(averaged_spectra_hash, stdev_hash, totalpeak,
            usedpeak)`` — same shape as :func:`constructionofMS1Ionclique`'s
            return value.
    """
    if msx == "ms1":
        averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = constructionofMS1Ionclique(msmslist, msmsdata, msdata, thresholds, filenames)
    elif msx == "ms2":
        averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = Ionclique_large(msmslist, msmsdata, msdata, thresholds, msx,filenames)

    #if msx == "ms1" or len(msmslist) > 0:
    #    averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = Ionclique_large(msmslist, msmsdata, msdata, thresholds, msx,filenames)
    #else:
    #    averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = Ionclique_small(msmslist, msmsdata, msdata, thresholds, msx,filenames)
    return(averaged_spectra_hash, stdev_hash, totalpeak, usedpeak)


def Ionclique_large(msmslist, msmsdata, msdata, thresholds, msx,filenames):
    """Average the MS2 (fragment) spectra of a DDA clique into one consensus MS2 spectrum.

    Flattens every raw fragment peak across ``msmslist``'s member spectra into a
    single peak pool, groups them into ion cliques via :func:`findIonclique`
    (using the raw ``ms1``/``ms2`` m/z threshold depending on ``msx``), keeps only
    ion cliques present in more than 70% of the member spectra, and for each
    surviving ion clique computes the mean (or, for >30 contributing peaks,
    median) m/z, mean base-peak-normalized intensity, and m/z/intensity standard
    deviations.

    For ``msx == "ms2"``, additionally identifies the most representative
    "base-peak" intensity profile to correlate every fragment against: the
    averaged precursor ion signal within the fragment spectrum if present near
    the mean precursor m/z, otherwise the fragment most correlated with the raw
    MS1 precursor-ion intensity profile across the clique's member spectra. Every
    averaged fragment's per-scan intensity list is then correlated against this
    reference profile and stored as ``"correlation"`` (consumed downstream by
    ``MassBankRecord.wash_spectra`` for statistical peak filtering).

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash.
        thresholds (dict): Threshold parameters; must include
            ``"ms1_mz_raw_thres"``/``"ms2_mz_raw_thres"`` (selected by ``msx``)
            and the thresholds consumed by ``massbankio.calcThreshold``.
        msx (str): ``"ms1"`` or ``"ms2"``, selecting the m/z threshold and
            whether the precursor-correlation logic above is applied.
        filenames (dict): File paths, including ``"logfilename"``.

    Returns:
        tuple: A tuple of four elements.

            - averaged_spectra_hash (dict): Averaged spectrum, mapping mean m/z
              to relative intensity (0-999 scale).
            - stdev_hash (dict): Per-averaged-peak statistics, mapping mean m/z
              to a dict with ``"intensity"``, ``"mzsd"``, ``"intsd"``,
              ``"number"``, ``"intlist"``/``"mzlist"``, ``"correlation"``
              (for ``msx == "ms2"``), and ``"annotation"`` (placeholder ``"?"``).
            - totalpeak (int): Total number of raw peaks considered.
            - usedpeak (int): Number of raw peaks retained in the averaged
              spectrum.

    Note:
        Returns ``({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)``
        if no ion cliques survive the presence-fraction filter, or if no base
        peak / no averaged peaks are found.
    """
    logfile = filenames["logfilename"]
    if msx == "ms1":
        mzthreshold = thresholds["ms1_mz_raw_thres"]
    else:
        mzthreshold = thresholds["ms2_mz_raw_thres"]

    allpeaks = {}
    allpeaks_list = []

    allmz = []
    allint = []
    allposition = []
    allid = []
    allspec = []
    #
    # Retrieve all peaks
    #
    counter = 0
    speccounter = 0
    for msms1 in msmslist:
        msms1_precs = list(msmsdata[msms1]['spectrum'].keys())
        for msms1_prec in msms1_precs:
            intensity = msmsdata[msms1]['spectrum'][float(msms1_prec)]

            allpeaks[counter] = [msms1, float(msms1_prec), intensity]

            allpeaks_list.append(counter)
            allmz.append(float(msms1_prec))
            allint.append(float(intensity))
            allposition.append(0)
            allid.append(msms1)
            allspec.append(speccounter)
            counter +=1
        speccounter += 1
    #
    # Sorted [id, filename_scan, m/z, intensity, order in m/z, order in intensity (Reversed),]
    #
    allpeaks_intsorted = sorted(allpeaks_list, key = lambda x:allint[x], reverse = True)
    allpeaks_mzsorted = sorted(allpeaks_list, key = lambda x:allmz[x])

    for i, targetid in enumerate(allpeaks_mzsorted):
        allposition[targetid] = i

    numberoftotalpeaks = len(allpeaks_mzsorted)

    allmz = np.array(allmz)
    allint = np.array(allint)
    allid = np.array(allid)
    allposition = np.array(allposition)
    allpeaks_mzsorted = np.array(allpeaks_mzsorted)
    allpeaks_intsorted = np.array(allpeaks_intsorted)
    allspec = np.array(allspec)
    numberofspectra = len(msmslist)

    cliques = findIonclique(numberofspectra, numberoftotalpeaks, allmz, allint,allid, allposition, allspec, allpeaks_mzsorted, allpeaks_intsorted, mzthreshold)


    numberofunusedpeaks = 0

    msmslist_number = len(msmslist)
    temp_hash = {}
    stdev_hash = {}

    intensity_hash = {}
    for msmsid in msmslist:
        intensity_hash[msmsid] = []

    for msmsclique in cliques:
        msmsclique_number = len(msmsclique)

        if msmsclique_number/msmslist_number <= 0.7:
            numberofunusedpeaks += msmsclique_number
            continue
        fragmentmzlist = []
        fragmentintensitylist = []
        precmzlist = []

        intensity_hash_temp = {}
        mz_hash_temp = {}

        for keyid in msmsclique:
            #[peakid, nodeid, fragmentmz, intensity, mzorder, intorder] = peak
            nodeid = allid[keyid] # 0: id in msmsdata
            fragmentmz = allmz[keyid] # 1: mz
            intensity = allint[keyid] # 1: intensity


            fragmentmzlist.append(fragmentmz)
            intensity = msmsdata[nodeid]['spectrum'][float(fragmentmz)]

            basePeakIntensity = float(msmsdata[nodeid]['basePeakIntensity'])
            fragmentintensitylist.append(float(intensity)/basePeakIntensity)

            intensity_hash_temp[nodeid] = intensity
            mz_hash_temp[nodeid] = fragmentmz
        #
        # list of intensity
        #
        intensity_alllist = []
        mz_alllist = []
        for msmsid in msmslist:
            if msmsid in intensity_hash_temp:
                intensity_alllist.append(int(intensity_hash_temp[msmsid]))
                mz_alllist.append(float(mz_hash_temp[msmsid]))
            else:
                intensity_alllist.append("")
                mz_alllist.append("")

        mean_mz = np.mean(fragmentmzlist)

        if len(fragmentmzlist) > 30:
            mean_mz = np.median(fragmentmzlist)


        mean_intensity = sum(fragmentintensitylist)/len(fragmentintensitylist) * 999
        temp_hash[mean_mz] = mean_intensity
        stdev_hash[mean_mz] = {"intensity": mean_intensity, "mzsd": stdev(fragmentmzlist), \
        "intsd": stdev(fragmentintensitylist)* 999, "number": len(fragmentmzlist), \
        "intlist": intensity_alllist, "mzlist": mz_alllist, "annotation":"?"}

    #
    # Searching Base Peak ion for MS1
    #
    basepeakintensity = 0
    basepeakmz = 0
    for mz, intensity in temp_hash.items():
        if intensity > basepeakintensity:
            basepeakintensity = intensity
            basepeakmz = mz

    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if basepeakmz == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if len(stdev_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)



    intlist_basepeak = stdev_hash[basepeakmz]["intlist"]

    #
    # For MS2 data
    #
    if msx == "ms2":
        #
        # Calculation of mean prec MZ
        #
        mean_precursormz, stdev_precursormz, list_precursormz = calcmz(msmslist, msmsdata)
        #
        # Calculation of a threshold level for m/z comparision
        #
        threshold_temp = massbankio.calcThreshold(mean_precursormz, msx, thresholds) *2
        #
        # Searching precursor ion in MS2
        #
        precmzinMS2 = 0
        precmzinMS2_diff = 1000000
        for mz, intensity in temp_hash.items():
            deltamz = abs(mean_precursormz-mz)
            if deltamz < precmzinMS2_diff:
                precmzinMS2_diff = deltamz
                precmzinMS2 = mz
        #
        # When precursor signal was existed in MS2, the intensity list of MS2 prec singal was used
        #
        if precmzinMS2_diff < threshold_temp:
            intlist_basepeak = stdev_hash[precmzinMS2]["intlist"]
        #
        # When precursor signal was existed in MS2,
        #
        else:
            #
            # Preparation of aintensity list of MS2 prec signal
            #
            msintlist =[]
            for msid in msmslist:
                ms1spec_hash = msdata[msid]["spectrum"]
                minimum = 1.0
                intmintemp = 0
                for mztemp, inttemp in ms1spec_hash.items():
                    difftemp = abs(mztemp - mean_precursormz)
                    if difftemp < threshold_temp:
                        if difftemp < minimum:
                            minimum = difftemp
                            intmintemp = float(inttemp)
                msintlist.append(intmintemp)
            if sum(msintlist) > 0:
                #
                # Searching MS2 fragment signal most correlated with MS2 prec signal
                #
                #
                # Searching Base Peak ion for MS1
                #
                mostcorrel = -5
                mostmz = 0
                for mz, intensity in temp_hash.items():
                    corcoeff_temp = massbankio.calccorrel(msintlist, stdev_hash[mz]["intlist"])
                    if corcoeff_temp > mostcorrel:
                        mostcorrel = corcoeff_temp
                        mostmz = mz
                intlist_basepeak = stdev_hash[mostmz]["intlist"]
            else:
                print("No prec in MS1")


    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if basepeakmz == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if len(stdev_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)

    temp_hash2 = {}
    for mz, intensity in temp_hash.items():
        #temp_hash2[mz] = int(intensity/basepeakintensity * 999)
        temp_hash2[mz] = intensity/basepeakintensity * 999
        stdev_hash[mz]["intsd"] = float(stdev_hash[mz]["intsd"]/basepeakintensity * 999.0)
        #stdev_hash[mz]["intensity"] = int(intensity/basepeakintensity * 999)
        stdev_hash[mz]["intensity"] = intensity/basepeakintensity * 999
        intlist_thisfragment = stdev_hash[mz]["intlist"]
        corcoeff = massbankio.calccorrel(intlist_basepeak, intlist_thisfragment)
        stdev_hash[mz]["correlation"] = corcoeff

    temp_hash_for_output = {}
    temp_stdevhash_for_output = {}
    mzlist = sorted(temp_hash2.keys())
    for mz in mzlist:
        temp_hash_for_output[mz] = temp_hash2[mz]
        temp_stdevhash_for_output[mz] = stdev_hash[mz]
    #print("Total:", numberoftotalpeaks, "Used:", numberoftotalpeaks-numberofunusedpeaks)

    return(temp_hash_for_output, temp_stdevhash_for_output, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)



def Ionclique_large_old(msmslist, msmsdata, msdata, thresholds, msx,filenames):
    """Legacy (pure-Python, non-Numba) implementation of :func:`Ionclique_large`.

    .. deprecated::
        Superseded by :func:`Ionclique_large`, which uses the Numba-jitted
        :func:`findIonclique` for ion-clique grouping instead of a pure-Python
        m/z-window walk. Kept for reference.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash.
        thresholds (dict): Threshold parameters (``"ms1_mz_raw_thres"``/
            ``"ms2_mz_raw_thres"`` as appropriate for ``msx``).
        msx (str): ``"ms1"`` or ``"ms2"``.
        filenames (dict): File paths, including ``"logfilename"``.

    Returns:
        tuple: Same 4-element shape as :func:`Ionclique_large`'s return value.

    See Also:
        :func:`Ionclique_large`: current implementation.
    """

    logfile = filenames["logfilename"]

    allpeaks = []

    #
    # Retrieve all peaks
    #
    counter = 0
    for msms1 in msmslist:
        msms1_precs = list(msmsdata[msms1]['spectrum'].keys())
        for msms1_prec in msms1_precs:
            intensity = msmsdata[msms1]['spectrum'][float(msms1_prec)]

            allpeaks.append([counter, msms1, float(msms1_prec), intensity])
            counter +=1
    #
    # Sorted [id, filename_scan, m/z, intensity, order in m/z, order in intensity (Reversed),]
    #
    allpeaks_temp = sorted(allpeaks, key = lambda x:x[2])
    allpeaks_temp = [allpeaks_temp[x]+[x] for x in range(len(allpeaks_temp))]
    allpeaks_temp = sorted(allpeaks_temp, key = lambda x:x[3], reverse = True)
    allpeaks_temp = [allpeaks_temp[x]+[x] for x in range(len(allpeaks_temp))]

    allpeaks_mzsorted = sorted(allpeaks_temp, key = lambda x:x[2])
    allpeaks_intsorted = sorted(allpeaks_temp, key = lambda x:x[3], reverse = True)

    numberoftotalpeaks = len(allpeaks_mzsorted)


    cliques = []
    usedpeak = {}
    for i, keypeak in enumerate(allpeaks_intsorted):
        #
        # If this peak has already used in other clique
        #
        keyid = keypeak[0]
        if keyid in usedpeak:
            continue

        keymz = keypeak[2]
        threshold_temp = massbankio.calcThreshold(keymz, msx, thresholds)
        lowerboundary = keymz - threshold_temp
        upperboundary = keymz + threshold_temp
        #
        # Position in allpeaks_mzsorted
        #
        posinmz = keypeak[4]
        #
        # Preparing clique seed
        #
        clique_temp = [keypeak]
        clique = []
        #
        # Searching to downward
        #
        step = 1

        while True:
            if (posinmz - step) < 0:
                break
            checkpeak = allpeaks_mzsorted[posinmz - step]
            chackmz = checkpeak[2]
            if chackmz < lowerboundary:
                break
            clique_temp.append(checkpeak)
            step +=1
        #
        # Searching to upward
        #
        step = 1

        while True:
            if (posinmz + step) >= numberoftotalpeaks:
                break
            checkpeak = allpeaks_mzsorted[posinmz + step]
            chackmz = checkpeak[2]
            if chackmz > upperboundary:
                break
            clique_temp.append(checkpeak)
            step +=1
        #
        # Sorting clique_temp with reved int order
        #
        clique_temp_sorted = sorted(clique_temp, key = lambda x:x[3], reverse = True)
        #
        #
        #
        usedfilenames = {}
        #
        # Check duplication and
        #
        for keypeak in clique_temp_sorted:
            keyid = keypeak[0]
            keyfilename = keypeak[1]
            if keyid in usedpeak:
                continue
            if keyfilename in usedfilenames:
                continue
            usedpeak[keyid] = i
            usedfilenames[keyfilename] = i
            clique.append(keypeak)

        cliques.append(clique)


    numberofunusedpeaks = 0

    msmslist_number = len(msmslist)
    temp_hash = {}
    stdev_hash = {}

    intensity_hash = {}
    for msmsid in msmslist:
        intensity_hash[msmsid] = []

    for msmsclique in cliques:
        msmsclique_number = len(msmsclique)

        if msmsclique_number/msmslist_number <= 0.7:
            numberofunusedpeaks += msmsclique_number
            continue
        fragmentmzlist = []
        fragmentintensitylist = []
        precmzlist = []

        intensity_hash_temp = {}
        mz_hash_temp = {}

        for peak in msmsclique:
            [peakid, nodeid, fragmentmz, intensity, mzorder, intorder] = peak
            fragmentmzlist.append(fragmentmz)
            intensity = msmsdata[nodeid]['spectrum'][float(fragmentmz)]

            basePeakIntensity = float(msmsdata[nodeid]['basePeakIntensity'])
            fragmentintensitylist.append(float(intensity)/basePeakIntensity)

            intensity_hash_temp[nodeid] = intensity
            mz_hash_temp[nodeid] = fragmentmz
        #
        # list of intensity
        #
        intensity_alllist = []
        mz_alllist = []
        for msmsid in msmslist:
            if msmsid in intensity_hash_temp:
                intensity_alllist.append(int(intensity_hash_temp[msmsid]))
                mz_alllist.append(float(mz_hash_temp[msmsid]))
            else:
                intensity_alllist.append("")
                mz_alllist.append("")

        mean_mz = np.mean(fragmentmzlist)

        if len(fragmentmzlist) > 30:
            mean_mz = np.median(fragmentmzlist)


        mean_intensity = sum(fragmentintensitylist)/len(fragmentintensitylist) * 999
        temp_hash[mean_mz] = mean_intensity
        stdev_hash[mean_mz] = {"intensity": mean_intensity, "mzsd": stdev(fragmentmzlist), \
        "intsd": stdev(fragmentintensitylist)* 999, "number": len(fragmentmzlist), \
        "intlist": intensity_alllist, "mzlist": mz_alllist, "annotation":"?"}

    #
    # Searching Base Peak ion for MS1
    #
    basepeakintensity = 0
    basepeakmz = 0
    for mz, intensity in temp_hash.items():
        if intensity > basepeakintensity:
            basepeakintensity = intensity
            basepeakmz = mz

    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if basepeakmz == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if len(stdev_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)



    intlist_basepeak = stdev_hash[basepeakmz]["intlist"]

    #
    # For MS2 data
    #
    if msx == "ms2":
        #
        # Calculation of mean prec MZ
        #
        mean_precursormz, stdev_precursormz, list_precursormz = calcmz(msmslist, msmsdata)
        #
        # Calculation of a threshold level for m/z comparision
        #
        threshold_temp = massbankio.calcThreshold(mean_precursormz, msx, thresholds) *2
        #
        # Searching precursor ion in MS2
        #
        precmzinMS2 = 0
        precmzinMS2_diff = 1000000
        for mz, intensity in temp_hash.items():
            deltamz = abs(mean_precursormz-mz)
            if deltamz < precmzinMS2_diff:
                precmzinMS2_diff = deltamz
                precmzinMS2 = mz
        #
        # When precursor signal was existed in MS2, the intensity list of MS2 prec singal was used
        #
        if precmzinMS2_diff < threshold_temp:
            intlist_basepeak = stdev_hash[precmzinMS2]["intlist"]
        #
        # When precursor signal was existed in MS2,
        #
        else:
            #
            # Preparation of aintensity list of MS2 prec signal
            #
            msintlist =[]
            for msid in msmslist:
                ms1spec_hash = msdata[msid]["spectrum"]
                minimum = 1.0
                intmintemp = 0
                for mztemp, inttemp in ms1spec_hash.items():
                    difftemp = abs(mztemp - mean_precursormz)
                    if difftemp < threshold_temp:
                        if difftemp < minimum:
                            minimum = difftemp
                            intmintemp = float(inttemp)
                msintlist.append(intmintemp)
            if sum(msintlist) > 0:
                #
                # Searching MS2 fragment signal most correlated with MS2 prec signal
                #
                #
                # Searching Base Peak ion for MS1
                #
                mostcorrel = -5
                mostmz = 0
                for mz, intensity in temp_hash.items():
                    corcoeff_temp = massbankio.calccorrel(msintlist, stdev_hash[mz]["intlist"])
                    if corcoeff_temp > mostcorrel:
                        mostcorrel = corcoeff_temp
                        mostmz = mz
                intlist_basepeak = stdev_hash[mostmz]["intlist"]
            else:
                print("No prec in MS1")


    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if basepeakmz == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if len(stdev_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)

    temp_hash2 = {}
    for mz, intensity in temp_hash.items():
        #temp_hash2[mz] = int(intensity/basepeakintensity * 999)
        temp_hash2[mz] = intensity/basepeakintensity * 999
        stdev_hash[mz]["intsd"] = float(stdev_hash[mz]["intsd"]/basepeakintensity * 999.0)
        #stdev_hash[mz]["intensity"] = int(intensity/basepeakintensity * 999)
        stdev_hash[mz]["intensity"] = intensity/basepeakintensity * 999
        intlist_thisfragment = stdev_hash[mz]["intlist"]
        corcoeff = massbankio.calccorrel(intlist_basepeak, intlist_thisfragment)
        stdev_hash[mz]["correlation"] = corcoeff

    temp_hash_for_output = {}
    temp_stdevhash_for_output = {}
    mzlist = sorted(temp_hash2.keys())
    for mz in mzlist:
        temp_hash_for_output[mz] = temp_hash2[mz]
        temp_stdevhash_for_output[mz] = stdev_hash[mz]
    #print("Total:", numberoftotalpeaks, "Used:", numberoftotalpeaks-numberofunusedpeaks)

    return(temp_hash_for_output, temp_stdevhash_for_output, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)



def Ionclique_small(msmslist, msmsdata, msdata, thresholds, msx, files):
    """Alternative small-scale ion-clique averaging using a direct m/z-window scan.

    Groups fragment (or precursor) peaks across ``msmslist``'s member spectra by
    direct m/z-window scanning (rather than the :func:`findIonclique`
    Numba-accelerated grouping used by :func:`Ionclique_large`/
    :func:`constructionofMS1Ionclique`) and averages each surviving group the
    same way (mean/median m/z, mean base-peak-normalized intensity). Intended
    for cliques with a small number of contributing peaks, where the simpler
    algorithm's lower constant overhead outweighs its worse asymptotic
    complexity.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash.
        thresholds (dict): Threshold parameters (``"ms1_mz_raw_thres"``/
            ``"ms2_mz_raw_thres"`` as appropriate for ``msx``).
        msx (str): ``"ms1"`` or ``"ms2"``.
        files (dict): File paths, including ``"logfilename"``.

    Returns:
        tuple: Same 4-element shape as :func:`Ionclique_large`'s return value:
            ``(averaged_spectra_hash, stdev_hash, totalpeak, usedpeak)``.

    Note:
        Not currently called by :func:`constructionofIonclique`, which always
        routes ``msx == "ms2"`` to :func:`Ionclique_large` regardless of clique
        size; kept available for callers that want the simpler algorithm
        explicitly.
    """
    logfile = files["logfilename"]
    #Dict of intensity
    allpeaks = {}
    #Dict of mz
    allpeaksmz = {}
    #Lower and Upper boundary
    minmz = 5000
    maxmz = 0
    #
    # Acuisition of all product ion signals
    # For each product ion spectra
    for msms1 in msmslist:
        #
        # list of mz of product ions
        #
        msms1_precs = list(msmsdata[msms1]['spectrum'].keys())
        #
        # For each mz of product ions
        #
        for msms1_prec in msms1_precs:
            if msms1_prec < minmz:
                minmz = msms1_prec
            if msms1_prec > maxmz:
                maxmz = msms1_prec
            #
            # Label Spec ID:mz of prec
            #
            label1 = msms1+":"+str(msms1_prec)
            #
            # Intensity of mz
            #
            intesity = msmsdata[msms1]['spectrum'][float(msms1_prec)]
            #
            # Storing data
            #
            allpeaks[label1] = intesity
            allpeaksmz[label1] = msms1_prec
    #
    # Sorted by mz values
    #
    allpeakmzsorted = sorted(allpeaks.keys(), key = lambda x:allpeaksmz[x])
    numberoftotalpeaks = len(allpeakmzsorted)
    #
    # Lower and Upper unit mass
    #
    minmz = int(minmz)
    maxmz = int(maxmz + 1)
    #
    #
    cliques = []
    #
    # For each unitmass
    for targetunitmass in range(minmz, maxmz):

        targetmz = targetunitmass * 1.0007
        #
        # Lower and upper baoundary of unis mass
        #
        lowerboundary = (targetunitmass -0.4) * 1.0007
        upperboundary = (targetunitmass +0.6) * 1.0007
        #
        # List of product ions within the boundaries
        #
        selectedpeakmz = []
        while True:
            #
            #
            #
            if len(allpeakmzsorted) == 0:
                break
            #
            # [101.11, 101.12, 101.13]
            # 101.11 is nowmz
            nowmz = allpeaksmz[allpeakmzsorted[0]]
            #
            #
            if lowerboundary < nowmz < upperboundary:
                #
                # 101.11 <= pop = [101.12, 101.13]
                # Move first mz to selectedpeakmz
                selectedpeakmz.append(allpeakmzsorted.pop(0))
            #
            # If value is out of range
            # upperboundary < nowmz
            #
            else:
                break
        #
        # Sort by peak intensity (Reversed )
        #
        selectedpeakmzsorted = sorted(selectedpeakmz, key = lambda x:allpeaks[x], reverse = True)
        #
        # for each produt ions
        #
        for label1 in selectedpeakmzsorted:
            (msms1, msms1_prec) = label1.split(":")
            msms1_prec = float(msms1_prec)
            hitclique = []
            hitclique_averagescore = []
            #
            # Distributed to suitable clique
            #
            for clique in cliques:
                averagescore_temp = []
                for label2 in clique:
                    (msms2, msms2_prec) = label2.split(":")

                    diff = abs(msms1_prec - float(msms2_prec))
                    averagescore_temp.append(diff)
                    #
                    # Ignored if any peak outside of the threshold level was existed
                    #
                    threshold_temp = massbankio.calcThreshold(max(msms1_prec, float(msms2_prec)), msx, thresholds) *2
                    if not diff < threshold_temp:
                        break
                    #
                    # Ignored if any peak derived from identical MSMS data was included in this clique
                    #
                    if msms1 == msms2:
                        break
                #
                # The mz can belong to the clique, this clique is stored in hitclisue
                #
                else:
                    hitclique.append(clique)
                    hitclique_averagescore.append(np.mean(averagescore_temp))
            #
            # If there is no hit clique, new clique was created
            #
            if hitclique == []:
                cliques.append([label1])
                continue
            #
            # mz is stored in the largest clique
            # If there are a few largest clique, mz is stored in a clique with the smallest averagescore_temp
            hitclique_order_sorted = sorted(range(len(hitclique_averagescore)), key=lambda x: len(hitclique[x]) - hitclique_averagescore[x], reverse = True)
            clique = hitclique[hitclique_order_sorted[0]]
            clique.append(label1)
    #
    #
    #
    numberofunusedpeaks = 0
    #
    # Number of all spectra
    #
    msmslist_number = len(msmslist)
    temp_hash = {}
    stdev_hash = {}
    for msmsclique in cliques:
        #
        # m/z observed more than 70% of product ion specs were employed
        #
        msmsclique_number = len(msmsclique)
        if msmsclique_number/msmslist_number <= 0.7:
            numberofunusedpeaks += msmsclique_number
            continue
        #
        # List of mz and intensity
        #
        fragmentmzlist = []
        fragmentintensitylist = []
        #
        # Dict id => mz or intensity
        #
        intensity_hash_temp = {}
        mz_hash_temp = {}

        for fragment in msmsclique:
            (nodeid, fragmentmz) = fragment.split(":")
            #print(nodeid, fragmentmz, msmsdata[nodeid]['spectrum'].keys())
            fragmentmzlist.append(float(fragmentmz))
            intensity = msmsdata[nodeid]['spectrum'][float(fragmentmz)]
            basePeakIntensity = float(msmsdata[nodeid]['basePeakIntensity'])
            fragmentintensitylist.append(float(intensity)/basePeakIntensity)

            intensity_hash_temp[nodeid] = int(intensity)
            mz_hash_temp[nodeid] = float(fragmentmz)
        #
        # Making a list of mz and intensity
        #
        intensity_alllist = []
        mz_alllist = []
        for msmsid in msmslist:
            if msmsid in intensity_hash_temp:
                intensity_alllist.append(int(intensity_hash_temp[msmsid]))
                mz_alllist.append(float(mz_hash_temp[msmsid]))
            else:
                intensity_alllist.append("")
                mz_alllist.append("")

        #
        # Calculation of mean or median mz
        #
        mean_mz = np.mean(fragmentmzlist)
        if len(fragmentmzlist) > 30:
            mean_mz = np.median(fragmentmzlist)


        mean_intensity = sum(fragmentintensitylist)/len(fragmentintensitylist) * 999
        temp_hash[mean_mz] = mean_intensity
        #stdev_hash[mean_mz] = {"meanintensity": mean_intensity, "mz": stdev(fragmentmzlist), "intensity": stdev(fragmentintensitylist)* 999, "number": len(fragmentmzlist)}
        stdev_hash[mean_mz] = {"intensity": mean_intensity, "mzsd": stdev(fragmentmzlist), \
        "intsd": stdev(fragmentintensitylist)* 999, "number": len(fragmentmzlist), \
        "intlist": intensity_alllist, "mzlist": mz_alllist, "annotation":"?"}



    #
    # Searching Base Peak ion for MS1
    #
    basepeakintensity = 0
    basepeakmz = 0
    for mz, intensity in temp_hash.items():
        if intensity > basepeakintensity:
            basepeakintensity = intensity
            basepeakmz = mz




    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if basepeakmz == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if len(stdev_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)


    intlist_basepeak = stdev_hash[basepeakmz]["intlist"]
    #
    # For MS2 data
    #
    if msx == "ms2":
        #
        # Calculation of mean prec MZ
        #
        mean_precursormz, stdev_precursormz, list_precursormz = calcmz(msmslist, msmsdata)
        #
        # Calculation of a threshold level for m/z comparision
        #
        threshold_temp = massbankio.calcThreshold(mean_precursormz, msx, thresholds) *2
        #
        # Searching precursor ion in MS2
        #
        precmzinMS2 = 0
        precmzinMS2_diff = 1000000
        for mz, intensity in temp_hash.items():
            deltamz = abs(mean_precursormz-mz)
            if deltamz < precmzinMS2_diff:
                precmzinMS2_diff = deltamz
                precmzinMS2 = mz
        #
        # When precursor signal was existed in MS2, the intensity list of MS2 prec singal was used
        #
        if precmzinMS2_diff < threshold_temp:
            intlist_basepeak = stdev_hash[precmzinMS2]["intlist"]
        #
        # When precursor signal was existed in MS2,
        #
        else:
            #
            # Preparation of aintensity list of MS2 prec signal
            #
            msintlist =[]
            for msid in msmslist:
                ms1spec_hash = msdata[msid]["spectrum"]
                minimum = 1.0
                intmintemp = 0
                for mztemp, inttemp in ms1spec_hash.items():
                    difftemp = abs(mztemp - mean_precursormz)
                    if difftemp < threshold_temp:
                        if difftemp < minimum:
                            minimum = difftemp
                            intmintemp = float(inttemp)
                msintlist.append(intmintemp)
            if sum(msintlist) > 0:
                #
                # Searching MS2 fragment signal most correlated with MS2 prec signal
                #
                #
                # Searching Base Peak ion for MS1
                #
                mostcorrel = -5
                mostmz = 0
                for mz, intensity in temp_hash.items():
                    corcoeff_temp = massbankio.calccorrel(msintlist, stdev_hash[mz]["intlist"])
                    if corcoeff_temp > mostcorrel:
                        mostcorrel = corcoeff_temp
                        mostmz = mz
                intlist_basepeak = stdev_hash[mostmz]["intlist"]
            else:
                print("No prec in MS1")


    if len(temp_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if basepeakmz == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)
    if len(stdev_hash) == 0:
        return({}, {}, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)



    temp_hash2 = {}
    for mz, intensity in temp_hash.items():
        #temp_hash2[mz] = int(intensity/basepeakintensity * 999)
        temp_hash2[mz] = intensity/basepeakintensity * 999
        stdev_hash[mz]["intsd"] = float(stdev_hash[mz]["intsd"]/basepeakintensity * 999.0)
        #stdev_hash[mz]["intensity"] = int(intensity/basepeakintensity * 999)
        stdev_hash[mz]["intensity"] = intensity/basepeakintensity * 999

        intlist_thisfragment = stdev_hash[mz]["intlist"]
        corcoeff = massbankio.calccorrel(intlist_basepeak, intlist_thisfragment)
        stdev_hash[mz]["correlation"] = corcoeff

    temp_hash_for_output = {}
    temp_stdevhash_for_output = {}
    mzlist = sorted(temp_hash2.keys())
    for mz in mzlist:
        temp_hash_for_output[mz] = temp_hash2[mz]
        temp_stdevhash_for_output[mz] = stdev_hash[mz]
    #print("Total:", numberoftotalpeaks, "Used:", numberoftotalpeaks-numberofunusedpeaks)

    return(temp_hash_for_output, temp_stdevhash_for_output, numberoftotalpeaks, numberoftotalpeaks-numberofunusedpeaks)




def generateSpecHash(spectrum, threshold = 5):
    """Convert a Hive-API spectrum object into a plain m/z-to-intensity dict.

    Args:
        spectrum: Hive-API spectrum object (as returned by
            ``qual_analysis.get_spectrum()``), providing
            ``get_spectrum_points()``.
        threshold (int, optional): Present for interface compatibility with
            :func:`generateSpecHash_old`; currently unused (no intensity
            filtering is applied). Defaults to ``5``.

    Returns:
        tuple: A tuple of three elements.

            - msmsdict (dict): Mapping of m/z (float) to intensity (float).
            - basePeakIntensity (float): Maximum intensity in the spectrum.
            - norm (float): L2 norm of the intensities.

        Returns ``({}, 0, [])`` if the spectrum has no data points.
    """

    intensity_list = []
    mz_list = []
    spectrum_point_collection = spectrum.get_spectrum_points()

    for spectrum_point in spectrum_point_collection:
        mz_list.append(spectrum_point.x)
        intensity_list.append(spectrum_point.intensity)
    #
    # If spectrum includes no data
    #
    if len(intensity_list) < 1:
        return({}, 0, [])

    basePeakIntensity = max(intensity_list)
    msmsdict = dict(zip(mz_list, intensity_list))
    norm = np.linalg.norm(intensity_list, ord=2)

    return(msmsdict, basePeakIntensity, norm)


def generateSpecHash_old(spectrum, threshold = 5):
    """Legacy implementation of :func:`generateSpecHash`.

    .. deprecated::
        Superseded by :func:`generateSpecHash`. Functionally near-identical
        (its intensity-threshold filtering is itself commented out, so no peaks
        are actually dropped by ``threshold`` here either); kept for reference.

    Args:
        spectrum: Hive-API spectrum object, providing ``get_spectrum_points()``.
        threshold (int, optional): Intended minimum intensity threshold; not
            currently applied (filtering logic is commented out).
            Defaults to ``5``.

    Returns:
        tuple: Same 3-element shape as :func:`generateSpecHash`'s return value:
            ``(msmsdict, basePeakIntensity, norm)``. Returns ``({}, 0, [])`` if
            the spectrum has no data points.

    See Also:
        :func:`generateSpecHash`: current implementation.
    """
    spectrum_point_collection = spectrum.get_spectrum_points()
    intensity_list = []
    intensity_list_selected = []
    mz_list_selected = []


    for spectrum_point in spectrum_point_collection:
        intensity_list.append(spectrum_point.intensity)
    #
    # If spectrum includes no data
    #
    if len(intensity_list) < 1:
        return({}, 0, [])

    basePeakIntensity = max(intensity_list)

    for spectrum_point in spectrum_point_collection:
        mz = spectrum_point.x
        intensity = spectrum_point.intensity

        #
        # Ignoring too small signals
        #
        #if intensity/ basePeakIntensity * 1000 < threshold:
        #    continue
        mz_list_selected.append(mz)
        intensity_list_selected.append(intensity)


    msmsdict = dict(zip(mz_list_selected, intensity_list_selected))
    norm = np.linalg.norm(intensity_list_selected, ord=2)
    return(msmsdict, basePeakIntensity, norm)


def loadhiveflles(files, mzlb, mzub, outputfiles):
    """Load MS2 spectra (and their paired MS1 precursor scans) from Hive files.

    Scans every ``.hmd`` file in ``files`` for MS2 spectra whose precursor m/z
    falls within ``[mzlb, mzub]`` (dropping any with fewer than 2 peaks),
    converts them via :func:`generateSpecHash`, and also fetches each MS2
    spectrum's corresponding MS1 (precursor) scan.

    Args:
        files (list[str]): Paths to candidate ``.hmd`` (Hive-format) files.
        mzlb (float): Lower boundary of precursor m/z to include.
        mzub (float): Upper boundary of precursor m/z to include.
        outputfiles (dict): File info dict; must include ``"logfilename"``
            (currently read but not used by this function directly).

    Returns:
        tuple: A tuple of two elements.

            - msmsdata (dict): Mapping of ``"<file_id>_<spectrum_index>"`` to a
              dict with ``"precursor_mz"``, ``"norm"``, ``"spectrum"``,
              ``"retentionTime"``, ``"basePeakIntensity"``, and
              ``"precursorIndex"``.
            - msdata (dict): Mapping of the same spectrum IDs to their
              corresponding MS1 spectrum data, each a dict with ``"spectrum"``,
              ``"retentionTime"``, ``"basePeakIntensity"``, and ``"scanindex"``.
    """
    logfile =outputfiles["logfilename"]


    msmsdata = {}
    msdata = {}


    lowerboundary = mzlb
    upperboundary = mzub
    usedfiles = []
    for file in files:

        if not ".hmd" in file:
            continue
        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".hmd", "")
        with HiveFile(file) as hive_file:
            ms2counter = 0
            for k in range(hive_file.get_measurement_count()):
                measurement = hive_file.get_measurement(k)

                for j in range(measurement.get_qual_analysis_count()):

                    qual_analysis = measurement.get_qual_analysis(j)
                    spectra_count = qual_analysis.get_spectra_count()
                    for i in range(spectra_count):
                        spectrum = qual_analysis.get_spectrum(i)
                        precursormz = spectrum.acquisition_condition.precursor_value
                        ms_level = str(spectrum.acquisition_condition.ms_level)
                        retentionTime = str(spectrum.retention_time)
                        precursorIndex = str(spectrum.precursor_index)

                        spectrum_number = str(i)

                        if not ms_level == "2":
                            continue
                        if not lowerboundary < float(precursormz) < upperboundary:
                            continue

                        msmsdict, basePeakIntensity, norm = generateSpecHash(spectrum, threshold = 0)

                        if len(msmsdict.keys()) < 2:
                            continue
                        # MS2 data
                        spectrum_number = str(i)
                        msmsspecid = data_number+"_"+spectrum_number
                        msmsdata[msmsspecid] = {'precursor_mz': float(precursormz), 'norm':norm,\
                        'spectrum': msmsdict, 'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity,\
                        'precursorIndex':precursorIndex}
                        ms2counter += 1

                        # MS1 data

                        spectrum = qual_analysis.get_spectrum(int(precursorIndex))

                        ms1retentionTime = str(spectrum.retention_time)
                        ms1dict, ms1basePeakIntensity, ms1norm = generateSpecHash(spectrum, threshold = 0)

                        msspecid = data_number+"_"+precursorIndex
                        msdata[msmsspecid] = {'spectrum': ms1dict, 'retentionTime': ms1retentionTime,\
                        'basePeakIntensity':ms1basePeakIntensity, 'scanindex':msspecid}


            writeLog(logfile, ["Loading file", lowerboundary, upperboundary, file, ms2counter])

    return(msmsdata,msdata)


def loadhivefilesbatch(files, listofboundaries, outputfiles):
    """Load MS2 spectra from Hive files, bucketed by multiple unit-mass boundary windows.

    Batch variant of :func:`loadhiveflles`: instead of a single ``[mzlb, mzub]``
    window, tests each MS2 spectrum's precursor m/z against every
    ``(lowerboundary, upperboundary)`` pair in ``listofboundaries`` and files it
    under the corresponding integer unit-mass bucket (a spectrum is only ever
    assigned to the first matching window). Also fetches each MS2 spectrum's
    corresponding MS1 (precursor) scan.

    Args:
        files (list[str]): Paths to candidate ``.hmd`` (Hive-format) files.
        listofboundaries (list[tuple[float, float]]): List of
            ``(lowerboundary, upperboundary)`` precursor m/z windows to bucket
            spectra into.
        outputfiles (dict): File info dict; must include ``"logfilename"``.

    Returns:
        tuple: A tuple of two elements.

            - unitmsmsdata (dict): ``unitmsmsdata[targetunitmass][spectrum_id]``
              -> MS2 spectrum dict (same shape as :func:`loadhiveflles`'s
              ``msmsdata`` values), where ``targetunitmass`` is the integer
              midpoint of the matched boundary window.
            - unitmsdata (dict): ``unitmsdata[targetunitmass][spectrum_id]`` ->
              corresponding MS1 spectrum dict (same shape as
              :func:`loadhiveflles`'s ``msdata`` values).
    """
    logfile =outputfiles["logfilename"]


    #msmsdata = {}
    #msdata = {}
    unitmsmsdata = {}
    unitmsdata = {}
    #lowerboundary = mzlb
    #upperboundary = mzub


    usedfiles = []
    for file in files:

        if not ".hmd" in file:
            continue
        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".hmd", "")
        with HiveFile(file) as hive_file:
            ms2counter = 0
            for k in range(hive_file.get_measurement_count()):
                measurement = hive_file.get_measurement(k)

                for j in range(measurement.get_qual_analysis_count()):

                    qual_analysis = measurement.get_qual_analysis(j)
                    spectra_count = qual_analysis.get_spectra_count()
                    for i in range(spectra_count):
                        spectrum = qual_analysis.get_spectrum(i)
                        precursormz = spectrum.acquisition_condition.precursor_value
                        ms_level = str(spectrum.acquisition_condition.ms_level)
                        retentionTime = str(spectrum.retention_time)
                        precursorIndex = str(spectrum.precursor_index)

                        spectrum_number = str(i)
                        if not ms_level == "2":
                            continue

                        for (lowerboundary, upperboundary) in listofboundaries:

                            if not lowerboundary < float(precursormz) <= upperboundary:
                                continue

                            targetunitmass = int((lowerboundary+upperboundary)/2)

                            msmsdict, basePeakIntensity, norm = generateSpecHash(spectrum, threshold = 0)

                            if len(msmsdict.keys()) < 2:
                                continue
                            # MS2 data
                            spectrum_number = str(i)
                            msmsspecid = data_number+"_"+spectrum_number
                            if not targetunitmass in unitmsmsdata:
                                unitmsmsdata[targetunitmass] = {}
                            if not targetunitmass in unitmsdata:
                                unitmsdata[targetunitmass] = {}



                            unitmsmsdata[targetunitmass][msmsspecid] = {'precursor_mz': float(precursormz), 'norm':norm,\
                            'spectrum': msmsdict, 'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity,\
                            'precursorIndex':precursorIndex}
                            ms2counter += 1

                            # MS1 data

                            spectrum = qual_analysis.get_spectrum(int(precursorIndex))

                            ms1retentionTime = str(spectrum.retention_time)
                            ms1dict, ms1basePeakIntensity, ms1norm = generateSpecHash(spectrum, threshold = 0)

                            msspecid = data_number+"_"+precursorIndex
                            unitmsdata[targetunitmass][msmsspecid] = {'spectrum': ms1dict, 'retentionTime': ms1retentionTime,\
                            'basePeakIntensity':ms1basePeakIntensity, 'scanindex':msspecid}
                            break


            writeLog(logfile, ["Loading file", lowerboundary, upperboundary, file, ms2counter])

    return(unitmsmsdata,unitmsdata)



def loadhivefiles_setup(files, listofboundaries, outputfiles, paths):
    """Index MS2 spectra from Hive files into per-unit-mass lookup text files.

    First-pass indexing step for :func:`loadspecfromtext`: for each ``.hmd``
    file, finds every MS2 spectrum matching one of ``listofboundaries``'
    windows, and appends a lightweight index row (spectrum ID, precursor scan
    index, precursor m/z, spectrum norm, retention time, base-peak intensity,
    measurement/qual-analysis indices) to
    ``paths["hmdtempms2"]<targetunitmass>.txt`` — one file per integer unit
    mass. Unlike :func:`loadhiveflles`/:func:`loadhivefilesbatch`, the full
    spectral peak data is not loaded or returned here.

    Args:
        files (list[str]): Paths to candidate ``.hmd`` (Hive-format) files.
        listofboundaries (list[tuple[float, float]]): List of
            ``(lowerboundary, upperboundary)`` precursor m/z windows to bucket
            spectra into.
        outputfiles (dict): File info dict; must include ``"logfilename"``.
        paths (dict): Path info dict; must include ``"hmdtempms2"`` (str),
            the filename prefix for the per-unit-mass index text files.

    Returns:
        tuple: ``(unitmsmsdata, unitmsdata)``, both empty dicts.

    Note:
        ``unitmsmsdata``/``unitmsdata`` are declared but never populated in this
        function; the actual index data is written only to the per-unit-mass
        text files under ``paths["hmdtempms2"]``, later read back by
        :func:`loadspecfromtext`.
    """
    logfile =outputfiles["logfilename"]
    hmdtempms2 = paths["hmdtempms2"]

    unitmassmin = 10000

    unitmassmax = 0

    #for (lowerboundary, upperboundary) in listofboundaries:
    #    targetunitmass = int((lowerboundary+upperboundary)/2)
    #    if unitmassmin > targetunitmass:
    #        unitmassmin = targetunitmass
    #    if unitmassmax < targetunitmass:
    #        unitmassmax = targetunitmass
    #print(unitmassmin, unitmassmax)

    unitmsmsdata = {}
    unitmsdata = {}

    usedfiles = []
    for file in files:

        if not ".hmd" in file:
            continue
        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".hmd", "")
        with HiveFile(file) as hive_file:
            print(file)
            outputtemp = {}
            ms2counter = 0
            for k in range(hive_file.get_measurement_count()):
                measurement = hive_file.get_measurement(k)

                for j in range(measurement.get_qual_analysis_count()):

                    qual_analysis = measurement.get_qual_analysis(j)
                    spectra_count = qual_analysis.get_spectra_count()
                    for i in range(spectra_count):
                        spectrum = qual_analysis.get_spectrum(i)
                        precursormz = spectrum.acquisition_condition.precursor_value
                        ms_level = str(spectrum.acquisition_condition.ms_level)
                        retentionTime = str(spectrum.retention_time)
                        precursorIndex = str(spectrum.precursor_index)

                        spectrum_number = str(i)
                        if not ms_level == "2":
                            continue

                        for (lowerboundary, upperboundary) in listofboundaries:
                            if not lowerboundary < float(precursormz) <= upperboundary:
                                continue

                            targetunitmass = int((lowerboundary+upperboundary)/2)


                            #targetunitmass = int(float(precursormz) / 1.0007 + 0.4)

                            #if not unitmassmin <= targetunitmass <= unitmassmax:
                            #    continue

                            msmsdict, basePeakIntensity, norm = generateSpecHash(spectrum, threshold = 0)

                            if len(msmsdict.keys()) < 2:
                                continue

                            spectrum_number = str(i)
                            msmsspecid = data_number+"_"+spectrum_number
                            if not targetunitmass in outputtemp:
                                outputtemp[targetunitmass] = []
                            templist = [msmsspecid, precursorIndex, float(precursormz), norm, retentionTime, basePeakIntensity, k, j]
                            outputtemp[targetunitmass].append(templist)


            for unitmasstemp, msmslisttemp in outputtemp.items():
                filenametemp = hmdtempms2+str(unitmasstemp)+".txt"
                if not os.path.isfile(filenametemp):
                    with open(filenametemp, mode="w",encoding='UTF-8', newline="") as f:
                        pass

                with open(filenametemp, mode="a",encoding='UTF-8', newline="") as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerows(msmslisttemp)
            writeLog(logfile, ["Loading file", lowerboundary, upperboundary, file, ms2counter])



    return(unitmsmsdata,unitmsdata)




def loadspecfromtext(files, listofboundaries, outputfiles, paths):
    """Load full MS2/MS1 spectra using the per-unit-mass index built by :func:`loadhivefiles_setup`.

    Second-pass loader: reads back the lightweight index rows written by
    :func:`loadhivefiles_setup` for every unit-mass window in
    ``listofboundaries``, then opens each ``.hmd`` file exactly once and fetches
    only the specific MS2 spectra (and their MS1 precursor scans, cached by
    precursor index to avoid re-fetching the same MS1 scan for multiple MS2
    children) referenced by the index. This two-pass design avoids scanning
    every spectrum in every file for every unit-mass window.

    Args:
        files (list[str]): Paths to the same ``.hmd`` (Hive-format) files that
            were previously indexed by :func:`loadhivefiles_setup`.
        listofboundaries (list[tuple[float, float]]): List of
            ``(lowerboundary, upperboundary)`` precursor m/z windows; must match
            those used when building the index.
        outputfiles (dict): File info dict; must include ``"logfilename"``.
        paths (dict): Path info dict; must include ``"hmdtempms2"`` (str), the
            filename prefix used for the per-unit-mass index text files.

    Returns:
        tuple: A tuple of two elements.

            - unitmsmsdata (dict): ``unitmsmsdata[targetunitmass][spectrum_id]``
              -> MS2 spectrum dict with ``"precursor_mz"``, ``"norm"``,
              ``"retentionTime"``, ``"basePeakIntensity"``, ``"precursorIndex"``,
              ``"msmsIndex"``, and ``"spectrum"`` (m/z-to-intensity dict).
            - unitmsdata (dict): ``unitmsdata[targetunitmass][spectrum_id]`` ->
              corresponding MS1 spectrum dict with ``"retentionTime"``,
              ``"scanindex"``, ``"spectrum"``, and ``"basePeakIntensity"``.

    See Also:
        :func:`loadhivefiles_setup`: builds the index this function reads.
    """
    logfile =outputfiles["logfilename"]

    hmdtempms2 = paths["hmdtempms2"]

    unitmsmsdata = {}
    unitmsdata = {}

    output_unitmsmsdata = {}

    #
    # List of hive files
    #
    hivefilenames = []
    for file in files:
        if not ".hmd" in file:
            continue
        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".hmd", "")
        hivefilenames.append(data_number)


    for (lowerboundary, upperboundary) in listofboundaries:
        targetunitmass = int((lowerboundary+upperboundary)/2)


        unitmsmsdata[targetunitmass] = {}
        unitmsdata[targetunitmass] = {}

        output_unitmsmsdata[targetunitmass] = {}
        for hivefilename in hivefilenames:
            output_unitmsmsdata[targetunitmass][hivefilename] = {}


        filenametemp = hmdtempms2+str(targetunitmass)+".txt"
        if not os.path.isfile(filenametemp):
            continue
        print(targetunitmass, filenametemp)
        with open(filenametemp, mode="r",encoding='UTF-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 8:
                    continue
                #templist = ["ms2", msmsspecid, data_number, spectrum_number,precursorIndex, float(precursormz), norm, retentionTime, basePeakIntensity]

                msmsspecid = row[0]
                msmsspecid_splitted = msmsspecid.split("_")
                msmsIndex = int(msmsspecid_splitted[-1])
                hivefilename = "_".join(msmsspecid_splitted[:-1])
                #print(msmsspecid, hivefilename, msmsIndex)

                precursormz = float(row[2])
                precursorIndex = int(row[1])

                norm = float(row[3])
                retentionTime = row[4]
                basePeakIntensity = float(row[5])
                k =row[6]
                j =row[7]

                output_unitmsmsdata[targetunitmass][hivefilename][msmsspecid] = {'precursorIndex':precursorIndex,'msmsIndex':msmsIndex, 'k':k, 'j':k}


                unitmsmsdata[targetunitmass][msmsspecid] = {'precursor_mz': float(precursormz), 'norm':norm,\
                'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity,\
                'precursorIndex':precursorIndex,'msmsIndex':msmsIndex}

                unitmsdata[targetunitmass][msmsspecid] = {'retentionTime': retentionTime, 'scanindex':hivefilename+"_"+str(precursorIndex)}


    ms2counter = 0
    for file in files:
        if not ".hmd" in file:
            continue
        data_number = file.split('\\')
        data_number = data_number[-1]
        hivefilename = data_number.replace(".hmd", "")
        print("Loading:", hivefilename)

        with HiveFile(file) as hive_file:
            #
            # record ms1
            #
            recordms1 = {}
            #
            #
            #
            measurement = hive_file.get_measurement(0)
            qual_analysis = measurement.get_qual_analysis(0)

            for targetunitmass, specdatas in output_unitmsmsdata.items():
                for msmsspecid, specdata in specdatas[hivefilename].items():
                    msmsIndex = specdata['msmsIndex']
                    precursorIndex = specdata['precursorIndex']
                    #
                    # ms2
                    #
                    #spectrum = hive_file.get_measurement(0).get_qual_analysis(0).get_spectrum(int(msmsIndex))
                    spectrum = qual_analysis.get_spectrum(int(msmsIndex))
                    ms2dict, ms2basePeakIntensity, ms2norm = generateSpecHash(spectrum, threshold = 0)
                    unitmsmsdata[targetunitmass][msmsspecid]['spectrum'] = ms2dict
                    #
                    # ms1
                    #
                    if precursorIndex in recordms1:
                        #print(precursorIndex)
                        unitmsdata[targetunitmass][msmsspecid]['spectrum'] = recordms1[precursorIndex]['spectrum']
                        unitmsdata[targetunitmass][msmsspecid]['basePeakIntensity'] = recordms1[precursorIndex]['basePeakIntensity']
                    else:
                        #print("*")
                        spectrum = qual_analysis.get_spectrum(int(precursorIndex))
                        #spectrum = hive_file.get_measurement(0).get_qual_analysis(0).get_spectrum(int(precursorIndex))
                        ms1dict, ms1basePeakIntensity, ms1norm = generateSpecHash(spectrum, threshold = 0)
                        unitmsdata[targetunitmass][msmsspecid]['spectrum'] = ms1dict
                        unitmsdata[targetunitmass][msmsspecid]['basePeakIntensity'] = ms1basePeakIntensity

                        recordms1[precursorIndex] = {'spectrum': ms1dict, 'basePeakIntensity': ms1basePeakIntensity}

                    ms2counter += 1

        writeLog(logfile, ["Loading file", file, ms2counter])

    return(unitmsmsdata,unitmsdata)


def checkhiveflles(path):
    """Count MS2 spectra per integer unit mass across a folder of Hive files.

    Diagnostic/inspection helper: scans every ``.hmd`` file matching ``path``,
    tallies how many MS2 spectra fall into each integer unit-mass bucket
    (``int(precursor_mz / 1.007 + 0.4)``), prints a per-file spectrum count and
    a sorted (ascending count) unit-mass histogram, and returns the histogram.

    Args:
        path (str): Glob pattern for the folder of ``.hmd`` files to scan.

    Returns:
        dict: Mapping of integer unit mass to the number of MS2 spectra found
            in that bucket.

    Note:
        Not currently called elsewhere in the pipeline; kept as a standalone
        diagnostic utility.
    """

    files = glob.glob(path)



    unitmass_hash = {}


    for file in files:
        numberofthisfile = 0
        if not ".hmd" in file:
            continue

        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".hmd", "")


        with HiveFile(file) as hive_file:
            for k in range(hive_file.get_measurement_count()):
                measurement = hive_file.get_measurement(k)

                for j in range(measurement.get_qual_analysis_count()):
                    qual_analysis = measurement.get_qual_analysis(j)
                    spectra_count = qual_analysis.get_spectra_count()
                    for i in range(spectra_count):
                        spectrum = qual_analysis.get_spectrum(i)
                        precursormz = spectrum.spectrum_condition.precursor_value
                        ms_level = str(spectrum.spectrum_condition.ms_level)

                        if not ms_level == "2":
                            continue
                        numberofthisfile += 1
                        unitmass = int(precursormz / 1.007+0.4)
                        if not unitmass in unitmass_hash:
                            unitmass_hash[unitmass] = 0
                        unitmass_hash[unitmass] += 1

        print(data_number, numberofthisfile)
    unitmasskeys = sorted(unitmass_hash.keys(), key = lambda x: unitmass_hash[x])
    for mass in unitmasskeys:
        print(mass, unitmass_hash[mass])


    return(unitmass_hash)



def loadShimadzumzXMLflles(path, targetmz, threshold = 0.02):
    """Load MS2 spectra near a target precursor m/z from Shimadzu-flavored mzXML files.

    .. deprecated::
        Obsolete; superseded by Hive-format loading
        (:func:`loadhiveflles`/:func:`loadspecfromtext`) for current pipeline
        use. Hard-codes the Shimadzu mzXML 3.2 schema namespace and only
        supports Base64-encoded (uncompressed) peak lists. Kept for reference.

    Args:
        path (str): Glob pattern for the folder of ``.mzXML`` files to scan.
        targetmz (float): Target precursor m/z to search for.
        threshold (float, optional): Half-width of the precursor m/z matching
            window (``targetmz ± threshold``). Defaults to ``0.02``.

    Returns:
        dict: Mapping of ``"<file_id>_<spectrum_number>"`` to a dict with
            ``"precursor_mz"``, ``"norm"``, ``"spectrum"``, ``"retentionTime"``,
            and ``"basePeakIntensity"``. Spectra with fewer than 2 peaks are
            excluded.

    See Also:
        :func:`loadmzXMLflles`: generalized (schema-detecting, multi-precision)
        successor.
    """
    msmsdata = {}
    files = glob.glob(path)


    for file in files:

        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".mzXML", "")

        print("loading files", file, data_number)


        tree = ET.parse(file)
        root = tree.getroot()
        attributes = root.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation']
        print(attributes)

        #for child in root:
        #    print(child.tag, child.attrib)
        for scan in root.iter('{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}scan'):

            if scan.attrib["msLevel"] != "2":
                continue
            if scan.attrib["peaksCount"] == "0":
                continue
            spectrum_number = scan.attrib["num"]
            basePeakIntensity = scan.attrib['basePeakIntensity']
            retentionTime = scan.attrib['retentionTime']

            for peaks in scan.findall('{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}precursorMz'):
                precursormz = peaks.text

                lowerboundary = targetmz - threshold
                upperboundary = targetmz + threshold

                if lowerboundary < float(precursormz) < upperboundary:
                    break
            else:
                continue


            for peaks in scan.findall('{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}peaks'):
                if peaks.text == None:
                    continue
                norm, msmsdict = decodeBase64(peaks.text, "base")
                #
                # Discurded if number of production ion is less than one
                #
                if len(msmsdict.keys()) >= 2:
                    msmsspecid = data_number+"_"+spectrum_number
                    msmsdata[msmsspecid] = {'precursor_mz': float(precursormz), 'norm':norm,\
                    'spectrum': msmsdict, 'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity}
    return(msmsdata)

def loadmzXMLflles(path, targetmz, threshold = 0.02):
    """Load MS2 spectra near a target precursor m/z from generic mzXML files.

    .. deprecated::
        Obsolete; superseded by Hive-format loading
        (:func:`loadhiveflles`/:func:`loadspecfromtext`) for current pipeline
        use. More general than :func:`loadShimadzumzXMLflles`: auto-detects the
        schema namespace from the file's ``schemaLocation`` attribute and
        supports both 64-bit (Base64) and 32-bit zlib-compressed peak encodings.
        Kept for reference.

    Args:
        path (str): Glob pattern for the folder of ``.mzXML`` files to scan.
        targetmz (float): Target precursor m/z to search for.
        threshold (float, optional): Half-width of the precursor m/z matching
            window (``targetmz ± threshold``). Defaults to ``0.02``.

    Returns:
        dict: Mapping of ``"<file_id>_<spectrum_number>"`` to a dict with
            ``"precursor_mz"``, ``"norm"``, ``"spectrum"``, ``"retentionTime"``,
            and ``"basePeakIntensity"``. Spectra with fewer than 2 peaks, or
            with an unsupported peak encoding, are excluded.
    """
    msmsdata = {}
    files = glob.glob(path)


    for file in files:

        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".mzXML", "")

        print("loading files", file, data_number)


        tree = ET.parse(file)
        root = tree.getroot()
        attributes = root.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation']
        sashimi = attributes.split(" ")[0]

        print(sashimi)
        for scan in root.iter('{'+sashimi+'}scan'):

            if scan.attrib["msLevel"] != "2":
                continue
            if scan.attrib["peaksCount"] == "0":
                continue
            spectrum_number = scan.attrib["num"]
            basePeakIntensity = scan.attrib['basePeakIntensity']
            retentionTime = scan.attrib['retentionTime']


            retentionTime = scan.attrib['retentionTime']

            for peaks in scan.findall('{'+sashimi+'}precursorMz'):
                precursormz = peaks.text



                lowerboundary = targetmz - threshold
                upperboundary = targetmz + threshold

                if lowerboundary < float(precursormz) < upperboundary:
                    break
            else:
                continue


            for peaks in scan.findall('{'+sashimi+'}peaks'):
                if peaks.text == None:
                    continue
                compressionType = "Base"
                if "compressionType" in peaks.attrib:
                    compressionType = peaks.attrib["compressionType"]


                precision = peaks.attrib["precision"]
                if precision == "64":
                    norm, msmsdict = decodeBase64(peaks.text, compressionType)
                elif precision == "32" and compressionType == "zlib":
                    norm, msmsdict = decodeBase32(peaks.text, compressionType)
                else:
                    continue
                #
                # Discurded if number of production ion is less than one
                #
                if len(msmsdict.keys()) >= 2:
                    msmsspecid = data_number+"_"+spectrum_number
                    msmsdata[msmsspecid] = {'precursor_mz': float(precursormz),
                    'norm':norm, 'spectrum': msmsdict, 'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity}
    return(msmsdata)

def calcrt(msmslist, msmsdata):
    """Compute the mean and standard deviation of retention time across a DDA clique.

    Parses each member spectrum's ``retentionTime`` string (stripping the
    ``"PT"``/``"S"`` ISO-8601-duration markers, e.g. ``"PT12.3S"`` -> ``12.3``)
    into a float before averaging.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash, with a
            ``"retentionTime"`` field.

    Returns:
        tuple: A tuple of three elements.

            - meanrt (float): Mean retention time across ``msmslist``.
            - stdevrt (float): Standard deviation of retention time.
            - listrt (list[float]): Individual retention time values, in the
              order iterated (sorted ``msmslist``).
    """

    rtlist = []
    for msmsid in sorted(msmslist):
        temprt =msmsdata[msmsid]['retentionTime']
        temprt = re.sub("PT", "", temprt)
        temprt = re.sub("S", "", temprt)
        rtlist.append(float(temprt))


    return(sum(rtlist)/len(rtlist), stdev(rtlist), rtlist)

def calcmz(msmslist, msmsdata):
    """Compute the mean (or median) and standard deviation of precursor m/z across a DDA clique.

    Uses the mean precursor m/z for cliques of 30 or fewer member spectra, or
    the median for larger cliques (more robust to outliers).

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash, with a
            ``"precursor_mz"`` field.

    Returns:
        tuple: A tuple of three elements.

            - meanmz (float): Mean (or median, for cliques > 30) precursor m/z.
            - stdevmz (float): Standard deviation of precursor m/z.
            - listmz (list[float]): Individual precursor m/z values, in
              ``msmslist`` order.
    """
    precursormzlist = []
    for msmsid in msmslist:
        precursormzlist.append(float(msmsdata[msmsid]['precursor_mz']))
    mean_precursormz = np.mean(precursormzlist)
    if len(precursormzlist) > 30:
        mean_precursormz = np.median(precursormzlist)

    return(mean_precursormz, stdev(precursormzlist), precursormzlist)


def calcpeaknum(msmslist, msmsdata):
    """Compute the mean and standard deviation of retention time across a DDA clique.

    Args:
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to one DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash, with a
            ``"retentionTime"`` field.

    Returns:
        tuple: A tuple of two elements ``(meannum, stdevnum)``: the mean and
            standard deviation of retention time across ``msmslist``.

    Note:
        Despite its name and docstring intent ("peak num"), this function's
        implementation is identical to :func:`calcrt` — it parses and averages
        each spectrum's ``retentionTime`` field, not its peak count. This looks
        like a copy-paste artifact from :func:`calcrt` rather than an actual
        peak-count calculation.
    """
    rtlist = []
    for msmsid in sorted(msmslist):
        temprt =msmsdata[msmsid]['retentionTime']
        temprt = re.sub("PT", "", temprt)
        temprt = re.sub("S", "", temprt)
        rtlist.append(float(temprt))


    return(sum(rtlist)/len(rtlist), stdev(rtlist))


def generationMassBankrecord(accession, msmslist, msmsdata, msdata, thresholds, filenames, report = "on"):
    """Build the paired MS2/MS1 MassBank records for one DDA (averaging) clique.

    Averages retention time and precursor m/z across the clique (via
    :func:`calcrt`/:func:`calcmz`), constructs the averaged MS2 fragment spectrum
    (via :func:`constructionofIonclique` with ``msx="ms2"``) and the averaged
    MS1 precursor spectrum (via :func:`constructionofIonclique` with
    ``msx="ms1"``), and populates two new ``MassBankRecord`` instances (one MS2,
    one MS1) from the dataset's record template with all the corresponding
    ``COMMENT:``/``AC$...``/``MS$...`` metadata fields (thresholds, corrections,
    retention time, precursor m/z, peak annotations, dates, record titles).
    Also updates the running ``counters`` dict with peak/spectrum statistics for
    this clique and logs a summary line for each generated record.

    Args:
        accession: An accession-ID generator object exposing a ``get()`` method
            returning ``(accession_ms2, accession_ms1)`` accession strings.
        msmslist (list): List of MS/MS (MS2) spectrum IDs belonging to this DDA
            clique.
        msmsdata (dict): Mapping of MS2 spectrum ID to spectrum hash.
        msdata (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash.
        thresholds (dict): Threshold parameters, including
            ``"ms1_mz_raw_thres"``, ``"ms2_mz_raw_thres"``,
            ``"ms1_mz_averaged_thres"``, ``"ms2_mz_averaged_thres"``,
            ``"ms1_mz_correction"``, ``"ms2_mz_correction"``,
            ``"retentiontime_thres"``, ``"dotproduct_thres"``, and ``"date"``.
        filenames (dict): File paths, including ``"logfilename"``,
            ``"record_templatefile"``, ``"record_templatefileMS1"``, and
            ``"numberoffiles"``.
        report (str, optional): ``"on"`` to print the averaged spectrum's peak
            statistics. Defaults to ``"on"``.

    Returns:
        tuple: A tuple of three elements.

            - mb (MassBankRecord): The generated MS2 record.
            - mb1 (MassBankRecord): The generated MS1 record.
            - counters (dict): Updated counter dict, as initialized by
              :func:`generationofcounters`.

    See Also:
        :func:`generationofcounters`: initializes the ``counters`` structure.
    """

    ms2_mz_raw_thres= thresholds["ms2_mz_raw_thres"]
    ms1_mz_raw_thres= thresholds["ms1_mz_raw_thres"]
    ms2_mz_averaged_thres= thresholds["ms2_mz_averaged_thres"]
    ms1_mz_averaged_thres= thresholds["ms1_mz_averaged_thres"]
    ms2_mz_correction= thresholds["ms2_mz_correction"]
    ms1_mz_correction= thresholds["ms1_mz_correction"]
    logfilename = filenames["logfilename"]
    record_templatefile = filenames["record_templatefile"]
    record_templatefileMS1 = filenames["record_templatefileMS1"]
    #record_template = filenames["template"]

    numberoffiles = filenames["numberoffiles"]

    counters = generationofcounters()

    accession_ms2, accession_ms1 = accession.get()
    counters["temp_record"][(accession_ms2, accession_ms1)]= msmslist

    #
    # Intensity of raw spectra
    #
    peakintensity = {}
    for msmsspecid in msmslist:
        basePeakIntensity = msmsdata[msmsspecid]['basePeakIntensity']
        totalintensity = 0
        for mz, peakmzintensity in msmsdata[msmsspecid]['spectrum'].items():
            totalintensity += peakmzintensity
        peakintensity[msmsspecid]= {"bp":basePeakIntensity, "sum":totalintensity}

        #unitmsmsdata[targetunitmass][msmsspecid] = {'precursor_mz': float(precursormz), 'norm':norm,\
        #'spectrum': msmsdict, 'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity,\
        #'precursorIndex':precursorIndex}
    counters["intentisy"]=peakintensity



    cliquesize = len(msmslist)
    counters["numberofclique"] += 1

    specinthisclique = str(cliquesize)
    #
    # Addtion to Tatal spec number
    #
    counters["spectra_counter"]["Used"] += cliquesize
    if not cliquesize in counters["cliquesize_counter"]:
        counters["cliquesize_counter"][cliquesize] = 0
    counters["cliquesize_counter"][cliquesize] += 1
    #
    #
    mean_rt, stdev_rt, list_rt = calcrt(msmslist, msmsdata)
    mean_precursormz, stdev_precursormz, list_precursormz = calcmz(msmslist, msmsdata)
    if report == "on":
        print(str(specinthisclique), "Retention time: {:.1f} {:.5f}".format(mean_rt, stdev_rt) ,\
         "Prec m/z: {:.5f} {:.5f}".format(mean_precursormz, stdev_precursormz))
    #
    # Generation of averaged spectra
    #
    averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = constructionofIonclique(msmslist, msmsdata, msdata, thresholds, "ms2", filenames)
    counters["peak_counter"]["Total"] += totalpeak
    counters["peak_counter"]["Used"] += usedpeak
    counters["numberofpeakms2"] += len(averaged_spectra_hash)

    if report == "on":
        for fragmentmz, mean_intensity in averaged_spectra_hash.items():
            print("\t {:.5f} {:>4f} {:.5f} {:.5f} {:.5f}  {:>4d}".format(fragmentmz, mean_intensity, stdev_hash[fragmentmz]["mzsd"]/fragmentmz*1000000,\
             stdev_hash[fragmentmz]["intsd"], stdev_hash[fragmentmz]["correlation"], stdev_hash[fragmentmz]["number"]))

    mb = massbankio.MassBankRecord()
    mb.read_template(record_templatefile)
    mb.add("ACCESSION:", accession_ms2)
    mb.add2("COMMENT:", "DATASET_FILE_NUM", numberoffiles)
    mb.add2("COMMENT:", "NUM_AVERAGED_SPEC", specinthisclique)
    mb.add2("COMMENT:", "RETENTION_TIME_STDEV", stdev_rt)
    mb.add2("COMMENT:", "RETENTION_TIME_LIST", list_rt)
    mb.add2("COMMENT:", "MS1_M/Z_RAW_THRES", ms1_mz_raw_thres)
    mb.add2("COMMENT:", "MS2_M/Z_RAW_THRES", ms2_mz_raw_thres)
    mb.add2("COMMENT:", "MS1_M/Z_AVERAGED_THRES", ms1_mz_averaged_thres)
    mb.add2("COMMENT:", "MS2_M/Z_AVERAGED_THRES", ms2_mz_averaged_thres)
    mb.add2("COMMENT:", "MS1_M/Z_CORRECTION", ms1_mz_correction)
    mb.add2("COMMENT:", "MS2_M/Z_CORRECTION", ms2_mz_correction)

    mb.add2("MS$DATA_PROCESSING:", "RECALIBRATE", "Data dependent recalibration using whole DDA dataset")
    thiskey = ("AC$MASS_SPECTROMETRY:", "MASS_ACCURACY")
    formattext = mb.record[thiskey]["format"]
    outputtext = formattext.format(ms1_mz_averaged_thres, ms2_mz_averaged_thres)
    mb.add2("AC$MASS_SPECTROMETRY:", "MASS_ACCURACY", outputtext)


    #mb.add2("AC$MASS_SPECTROMETRY:", "MASS_ACCURACY", "MS1 " + str(ms1_mz_averaged_thres) + " ppm, MS2 " + str(ms2_mz_averaged_thres) + " ppm.")

    mb.add2("COMMENT:", "RETENTION_TIME_THRES", thresholds["retentiontime_thres"])
    mb.add2("COMMENT:", "DOTPRODUCT_THRES", thresholds["dotproduct_thres"])
    mb.add2("AC$CHROMATOGRAPHY:", "RETENTION_TIME", mean_rt)
    mb.add2("COMMENT:", "PRECURSOR_M/Z_STDEV", stdev_precursormz)
    mb.add2("COMMENT:", "PRECURSOR_M/Z_LIST", list_precursormz)
    mb.add2("MS$FOCUSED_ION:", "PRECURSOR_M/Z", mean_precursormz)
    mb.add("PK$NUM_PEAK:", len(averaged_spectra_hash.keys()))
    mb.add_spectra(averaged_spectra_hash)
    #mb.add_annotation_from_items(stdev_hash, itemlist = ["int", "annotation",  "m/z_std", "m/z_std(ppm)","int_std","cor","n","mzlist"])
    # 251221 fix
    #
    mb.add_annotation_from_items(stdev_hash, itemlist = ["int", "annotation",  "m/z_std", "m/z_std(ppm)","int_std","cor","n","mzlist","intlist"])
    mb.add_date(date = thresholds["date"], target = "Averaged")
    mb.add_recordtitle()

    logtemp = {}
    logtemp["mean_precursormz"] = mean_precursormz
    logtemp["mean_rt"] = mean_rt
    logtemp["specinthisclique"] = specinthisclique
    logtemp["accession"] = accession_ms2
    logtemp["numofions"] = len(averaged_spectra_hash.keys())
    #logtemp["rawspectra"] = msmslist
    writeLog(filenames["logfilename"], ["MS2 Record generated", accession_ms2, logtemp])
    #
    #
    # Generation of averaged spectra
    #
    #averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = constructionofMS1Ionclique(msmslist, msdata, msdata, thresholds, filenames)
    averaged_spectra_hash, stdev_hash, totalpeak, usedpeak = constructionofIonclique(msmslist, msdata, msdata, thresholds, "ms1", filenames)
    counters["ms1peak_counter"]["Total"] += totalpeak
    counters["ms1peak_counter"]["Used"] += usedpeak
    counters["numberofpeakms1"] += len(averaged_spectra_hash)


    mb1 = massbankio.MassBankRecord()
    mb1.read_template(record_templatefile, msx = "ms1")
    mb1.add("ACCESSION:", accession_ms1)
    mb1.add2("COMMENT:", "DATASET_FILE_NUM", numberoffiles)
    mb1.add2("COMMENT:", "NUM_AVERAGED_SPEC", specinthisclique)
    mb1.add2("COMMENT:", "RETENTION_TIME_STDEV", stdev_rt)
    mb1.add2("COMMENT:", "MS1_M/Z_RAW_THRES", ms1_mz_raw_thres)
    mb1.add2("COMMENT:", "MS1_M/Z_AVERAGED_THRES", ms1_mz_averaged_thres)
    mb.add2("COMMENT:", "MS1_M/Z_CORRECTION", ms1_mz_correction)

    mb1.add2("AC$MASS_SPECTROMETRY:", "MASS_ACCURACY", str(ms1_mz_averaged_thres) + " ppm")
    mb1.add2("AC$CHROMATOGRAPHY:", "RETENTION_TIME", mean_rt)
    mb1.add("PK$NUM_PEAK:", len(averaged_spectra_hash.keys()))
    mb1.add_spectra(averaged_spectra_hash)
    #mb1.add_annotation_from_items(stdev_hash, itemlist = ["m/z(std)", "ppm(std)","int","int(std)","n","intlist"])

    mb1.add_annotation_from_items(stdev_hash, itemlist = ["int","annotation", "n","intlist"])# stdev data were removed 241216

    mb1.add_date(date = thresholds["date"], target = "Averaged")
    mb1.add_recordtitle_ms1(accession_ms2)



    mslist = []
    for msmsid in msmslist:
        mslist.append(msdata[msmsid]['scanindex'])
    countersms1 ={}
    countersms1["rawspec"] =mslist
    countersms1['numberofpeaks'] = counters['numberofpeakms1']
    countersms1['peak_counter'] = counters['ms1peak_counter']

    countersms2 ={}
    countersms2["rawspec"] =msmslist
    countersms2['numberofpeaks'] = counters['numberofpeakms2']
    countersms2['peak_counter'] = counters['peak_counter']


    logtemp = {}
    logtemp["mean_precursormz"] = ""
    logtemp["mean_rt"] = mean_rt
    logtemp["specinthisclique"] = specinthisclique
    logtemp["accession"] = accession_ms1
    logtemp["numofions"] = len(averaged_spectra_hash.keys())
    #logtemp["rawspectra"] = mslist
    writeLog(filenames["logfilename"], ["MS1 Record generated", accession_ms1, logtemp])

    mb.add2("COMMENT:", "COUNTER", counters)
    mb1.add2("COMMENT:", "COUNTER", countersms1)

    return(mb, mb1, counters)

def generationofcounters():
    """Initialize an empty counter dictionary for tracking averaging statistics.

    Builds the counter structure accumulated across calls to
    :func:`generationMassBankrecord` (and merged across batches via
    :func:`addtonofcounters`), tracking clique/spectrum/peak counts.

    Returns:
        dict: Dictionary with keys ``"temp_record"`` (dict, populated per
            clique), ``"spectra_counter"``, ``"peak_counter"``,
            ``"ms1peak_counter"`` (each ``{"Total": int, "Used": int}``),
            ``"numberofclique"`` (int), ``"numberofpeakms1"`` (int),
            ``"numberofpeakms2"`` (int), and ``"cliquesize_counter"`` (dict of
            clique size to occurrence count).
    """

    counters ={}
    counters["temp_record"] = {}
    counters["spectra_counter"] = {"Total": 0, "Used": 0}
    counters["peak_counter"] = {"Total": 0, "Used": 0}
    counters["ms1peak_counter"] = {"Total": 0, "Used": 0}
    counters["numberofclique"] = 0
    counters["numberofpeakms1"]= 0
    counters["numberofpeakms2"]= 0
    counters["cliquesize_counter"] = {}
    #counters["cliquesize_counter"] = {x:0 for x in range(50000)}
    return(counters)


def addtonofcounters(counters, counters_add):
    """Merge a second counter dict's totals into a running counter dict.

    Args:
        counters (dict): Running counter dict (as produced by
            :func:`generationofcounters`), updated and returned in place.
        counters_add (dict): Counter dict whose values are added into
            ``counters``.

    Returns:
        dict: The updated ``counters`` dict (same object as the input).
    """
    counters["temp_record"].update(counters_add["temp_record"])

    for words in ["Total", "Used"]:
        counters["spectra_counter"][words] += counters_add["spectra_counter"][words]
        counters["peak_counter"][words] += counters_add["peak_counter"][words]
        counters["ms1peak_counter"][words] += counters_add["ms1peak_counter"][words]
    counters["numberofclique"] += counters_add["numberofclique"]
    counters["numberofpeakms1"]+= counters_add["numberofpeakms1"]
    counters["numberofpeakms2"]+= counters_add["numberofpeakms2"]
    for x in counters_add["cliquesize_counter"].keys():
        if not x in counters["cliquesize_counter"]:
            counters["cliquesize_counter"][x] = 0
        counters["cliquesize_counter"][x] += counters_add["cliquesize_counter"][x]
    return(counters)


def averaging(msmsdata_temp, msdata_temp, targetunitmass, paths, thresholds, filenames, report):
    """Average one unit mass's worth of raw spectra into MassBank records.

    Discards spectra with a base-peak intensity of 200 or below, groups the
    remainder into similarity cliques via :func:`constructonofDDAclique`,
    generates one MS2/MS1 MassBank record pair per clique via
    :func:`generationMassBankrecord` (skipping cliques whose averaged spectrum
    ends up empty), writes each surviving record pair to
    ``paths["outputpath"]``, and merges each clique's per-clique counters into a
    running total via :func:`addtonofcounters`.

    Args:
        msmsdata_temp (dict): Mapping of MS2 spectrum ID to spectrum hash for
            this unit mass, before base-peak-intensity filtering.
        msdata_temp (dict): Mapping of MS2 spectrum ID to its corresponding MS1
            spectrum hash, before filtering.
        targetunitmass (int or float): Target unit mass being processed, used to
            seed the accession-number generator.
        paths (dict): Dictionary of relevant folder paths; must include
            ``"outputpath"`` for writing generated records.
        thresholds (dict): Threshold parameters, as consumed by
            :func:`generationMassBankrecord`.
        filenames (dict): File paths, including ``"accession_stem"`` and
            ``"logfilename"``, as consumed by :func:`generationMassBankrecord`.
        report (str): ``"on"`` to print progress (clique construction size,
            etc.); any other value suppresses it.

    Returns:
        bool: Always ``True``. ``output`` is initialized to ``False`` but is
            unconditionally set to ``True`` after the clique loop finishes, with
            no early-return path, so the return value does not currently signal
            whether any records were actually written.
    """
    class Accession_generator:
        """Sequential MassBank accession-ID generator for one unit-mass batch.

        Produces paired MS2 (``...P<number>``) / MS1 (``...S<number>``)
        accession strings with a shared, auto-incrementing numeric suffix.

        Attributes:
            x (int): Next accession number to assign.
            digit (int): Zero-padded width of the numeric suffix.
            stem (str): Accession prefix shared by all generated IDs.
        """
        def __init__(self,stem, startnumber, digit = 5):
            """Initialize the generator with a starting accession number.

            Args:
                stem (str): Accession prefix (e.g. dataset accession stem).
                startnumber (int): First numeric suffix to assign.
                digit (int, optional): Zero-padded width of the numeric suffix.
                    Defaults to ``5``.
            """
            self.x = startnumber
            self.digit = digit
            self.stem = stem
        def get(self):
            """Return the next paired MS2/MS1 accession IDs and advance the counter.

            Returns:
                tuple: A tuple of two elements ``(ms2, ms1)`` (str), the next
                    MS2 (``P``-prefixed) and MS1 (``S``-prefixed) accession IDs.
            """
            stem = self.stem
            digit = self.digit
            x = self.x
            self.x += 1
            text2 = "P{:0"+str(digit)+"d}"
            ms2 = stem+text2.format(x)

            text1 = "S{:0"+str(digit)+"d}"
            ms1 = stem+text1.format(x)

            return(ms2, ms1)

    output = False
    accession_stem = filenames["accession_stem"]
    accession_temp = Accession_generator(accession_stem, int(targetunitmass) * 1000, 7)

    #
    # Addtion to Tatal spec number
    #

    if report == "on":
        print("Construciton of DDA cliques", len(msmsdata_temp))
    #
    # Remove low basePeakIntensity data < 200
    #
    msmsdata = {}
    msdata = {}
    for msmsid, msmsspec in msmsdata_temp.items():
        basepeakintensity = float(msmsspec['basePeakIntensity'])
        if basepeakintensity > 200:
            msmsdata[msmsid] = msmsspec
            msdata[msmsid] = msdata_temp[msmsid]


    DDAcliques = constructonofDDAclique(msmsdata, thresholds, filenames, report)

    #
    # Generation of averaged spectra for each clique
    #
    for msmslist in DDAcliques:

        mb2, mb1, counters = generationMassBankrecord(accession_temp, msmslist, msmsdata, msdata, thresholds, filenames, report = report)

        spec = mb2.get_spectra()
        accession_ms2 = mb2.get_accession()
        if len(spec) == 0:
            writeLog(filenames["logfilename"], ["MS2 Record was cancelled due to no spec data", accession_ms2])
            continue
        spec = mb1.get_spectra()
        accession_ms1 = mb1.get_accession()
        if len(spec) == 0:
            writeLog(filenames["logfilename"], ["MS1 Record was cancelled due to no spec data", accession_ms1])
            continue

        mb2.show(path = paths["outputpath"])
        mb1.show(path = paths["outputpath"])

    output = True

    return(output)

def setup_averaging(dataset, batchsize, numcpu, cleanstart = True):
    """Set up folders, thresholds, and filenames for a dataset's spectral averaging run.

    Reads the dataset's ``RecordTemplate.txt`` settings, (re)creates the
    ``step1`` output folder (always) and, when ``cleanstart`` is ``True``, also
    the ``step2``/``step3``/``step4``/``annotation`` folders and the temporary
    Hive MS2 index folder (``hmd/ms2/``). Assembles the m/z/retention-time/
    dot-product thresholds and output file naming used throughout the averaging
    pipeline, and logs the run's configuration.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).
        batchsize (int): Batch size to log (not otherwise used by this function).
        numcpu (int): Number of CPU cores to log (not otherwise used by this
            function).
        cleanstart (bool, optional): When ``True``, deletes and recreates the
            ``step2``/``step3``/``step4``/``annotation`` output folders.
            Defaults to ``True``.

    Returns:
        tuple: A tuple of five elements.

            - paths (dict): ``{"outputpath": "<step1 dir>/", "path": "<hmd glob>",
              "hmdtempms2": "<hmd/ms2 dir>/"}``.
            - thresholds (dict): ``"ms2_mz_raw_thres"``, ``"ms1_mz_raw_thres"``,
              ``"ms2_mz_averaged_thres"``, ``"ms1_mz_averaged_thres"``
              (float, ppm thresholds), ``"ms2_mz_correction"``,
              ``"ms1_mz_correction"`` (calibration curves),
              ``"retentiontime_thres"``, ``"dotproduct_thres"`` (float),
              ``"size_minimum"`` (int, minimum clique size), ``"unit"``
              (always ``"ppm"``), and ``"date"`` (today's date string).
            - filenames (dict): ``"record_templatefile"``,
              ``"record_templatefileMS1"`` (both the dataset's
              ``RecordTemplate.txt`` path), ``"logfilename"`` (today's averaging
              log path), ``"accession_stem"``, and ``"numberoffiles"`` (int).
            - usedfiles (list[str]): File-ID stems of every ``.hmd`` file found
              under the dataset's ``hmd/`` folder.
            - settings (dict): Parsed ``RecordTemplate.txt`` settings, as
              returned by ``massbankio.read_RecordTemplate``.

    See Also:
        :func:`averaging`: consumes the paths/thresholds/filenames built here.
    """

    settings = massbankio.read_RecordTemplate("../data/"+dataset+"/RecordTemplate.txt")

    path = "../data/"+dataset+"/hmd/*"

    outputdir = "../data/"+dataset+"/step1"
    outputdir2 = "../data/"+dataset+"/step2"
    outputdir3 = "../data/"+dataset+"/step3"
    outputdir4 = "../data/"+dataset+"/step4"

    outputdirannotation = "../data/"+dataset+"/annotation"

    logpath = "../data/"+dataset+"/log"

    paths = {}

    paths["outputpath"] = outputdir + "/"
    paths["path"] = path
    #paths["logpath"] = logpath
    #
    # Preparation of folders

    if os.path.isdir(outputdir):
        shutil.rmtree(outputdir)
    os.mkdir(outputdir)

    if cleanstart == True:
        if os.path.isdir(outputdir2):
            shutil.rmtree(outputdir2)
        if os.path.isdir(outputdir3):
            shutil.rmtree(outputdir3)
        if os.path.isdir(outputdir4):
            shutil.rmtree(outputdir4)
        if os.path.isdir(outputdirannotation):
            shutil.rmtree(outputdirannotation)
        os.mkdir(outputdir2)
        os.mkdir(outputdir3)
        os.mkdir(outputdir4)
        os.mkdir(outputdirannotation)


    #hmdtempms1 = "../data/"+dataset+"/hmd/ms1/"
    hmdtempms2 = "../data/"+dataset+"/hmd/ms2/"
    #paths["hmdtempms1"] = hmdtempms1
    paths["hmdtempms2"] = hmdtempms2
    #if os.path.isdir(hmdtempms1):
    #    shutil.rmtree(hmdtempms1)
    if os.path.isdir(hmdtempms2):
        shutil.rmtree(hmdtempms2)
    #os.mkdir(hmdtempms1)
    os.mkdir(hmdtempms2)

    if not os.path.isdir(logpath):
        os.mkdir(logpath)

    thresholds ={}
    thresholds["ms2_mz_raw_thres"] = float(settings["ms2_mz_raw_thres"])
    thresholds["ms1_mz_raw_thres"] = float(settings["ms1_mz_raw_thres"])
    thresholds["ms2_mz_averaged_thres"] = float(settings["ms2_mz_averaged_thres"])
    thresholds["ms1_mz_averaged_thres"] = float(settings["ms1_mz_averaged_thres"])
    thresholds["ms2_mz_correction"] = settings["ms2_mz_correction"]
    thresholds["ms1_mz_correction"] = settings["ms1_mz_correction"]
    thresholds["retentiontime_thres"] = float(settings["retentiontime_thres"])
    thresholds["dotproduct_thres"] = float(settings["dotproduct_thres"])
    thresholds["size_minimum"] = int(settings["size_minimum"])
    thresholds["unit"] = "ppm"
    thresholds["date"] = datetime.datetime.now().strftime('%Y.%m.%d')


    accession_stem = settings["accession_stem"]
    dt_now = datetime.datetime.now()
    filenames = {}
    filenames["record_templatefile"] = "../data/"+dataset+"/RecordTemplate.txt"
    filenames["record_templatefileMS1"] = "../data/"+dataset+"/RecordTemplate.txt"
    filenames["logfilename"] = "../data/"+dataset+"/log/"+accession_stem+"_averaging_"+ dt_now.strftime('%Y.%m.%d') + ".txt"
    filenames["accession_stem"] = accession_stem

    mode = settings["mode"]
    accession_stem = settings["accession_stem"]

    files = glob.glob(paths["path"])

    usedfiles = []

    for file in files:
        if not ".hmd" in file:
            continue
        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".hmd", "")
        usedfiles.append(data_number)

    writeLog(filenames["logfilename"], ["Started"]+[accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), batchsize, numcpu], mode = "w")
    writeLog(filenames["logfilename"], ["Threshold"]+[accession_stem, thresholds["ms2_mz_raw_thres"],\
     thresholds["ms1_mz_raw_thres"], thresholds["size_minimum"], mode], mode = "w")

    writeLog(filenames["logfilename"], ["Target files", len(usedfiles), usedfiles])


    for threshold, value in thresholds.items():
        writeLog(filenames["logfilename"], ["Thresholds", threshold, value])

    filenames["numberoffiles"] =len(usedfiles)

    return(paths, thresholds, filenames, usedfiles, settings)



def averaging_statistics(paths, settings, filenames, Totalspectracounter):
    """Aggregate per-record counters across a dataset and log summary statistics.

    Re-reads every generated MS2 output record's ``COMMENT: COUNTER`` field
    (written by :func:`generationMassBankrecord`), merges them all via
    :func:`addtonofcounters`, and writes a series of summary log lines: total
    averaged spectra/peaks, total/used raw spectra and peaks, and a
    clique-size histogram with running (cumulative) totals.

    Args:
        paths (dict): Path info dict; must include ``"outputpath"`` (the step1
            output folder containing the generated MS2 records).
        settings (dict): Parsed ``RecordTemplate.txt`` settings, as returned by
            :func:`setup_averaging`; must include ``"accession_stem"``.
        filenames (dict): File paths, including ``"logfilename"``.
        Totalspectracounter (int): Total number of raw spectra loaded for this
            dataset (across all unit-mass batches), logged as-is.

    Returns:
        None

    See Also:
        :func:`setup_averaging`: builds ``paths``/``settings``/``filenames``.
    """
    #
    # Statistics
    #
    counters = generationofcounters()
    #
    files = glob.glob(paths["outputpath"]+"*")
    #
    for file in files:
        if not ".txt" in file:
            continue
        if not settings["accession_stem"]+"P" in file:
            continue
        mb2 = massbankio.MassBankRecord()
        mb2.read_MassBankRecord(file)
        counter = eval(mb2.get_anything("COMMENT:", "COUNTER"))
        counters = addtonofcounters(counters, counter)

    #
    # Log total and used peaks and spectra
    #
    templist = [counters["numberofclique"], counters["numberofpeakms2"], counters["numberofpeakms1"],Totalspectracounter, \
    counters["spectra_counter"]["Used"], counters["peak_counter"]["Total"], \
    counters["peak_counter"]["Used"], counters["ms1peak_counter"]["Total"], counters["ms1peak_counter"]["Used"]]

    writeLog(filenames["logfilename"], ["summary"]+templist)
    writeLog(filenames["logfilename"], ["summary", "Averaged spectra",counters["numberofclique"]])
    writeLog(filenames["logfilename"], ["summary", "Peaks in averaged MS2 spectra",counters["numberofpeakms2"]])
    writeLog(filenames["logfilename"], ["summary", "Peaks in averaged MS1 spectra",counters["numberofpeakms1"]])
    writeLog(filenames["logfilename"], ["summary", "Total raw spectra",Totalspectracounter])
    writeLog(filenames["logfilename"], ["summary", "Used raw spectra",counters["spectra_counter"]["Used"]])
    writeLog(filenames["logfilename"], ["summary", "Total peaks in raw MS2 spectra",counters["peak_counter"]["Total"]])
    writeLog(filenames["logfilename"], ["summary", "Used peaks in raw MS2 spectra",counters["peak_counter"]["Used"]])
    writeLog(filenames["logfilename"], ["summary", "Total peaks in raw MS1 spectra",counters["ms1peak_counter"]["Total"]])
    writeLog(filenames["logfilename"], ["summary", "Used peaks in raw MS1 spectra",counters["ms1peak_counter"]["Used"]])

    sum = 0
    sum_rawspec = 0
    cliquesize_counter = counters["cliquesize_counter"]
    for i in sorted(cliquesize_counter.keys()):
        if cliquesize_counter[i] == 0:
            continue
        sum += cliquesize_counter[i]
        sum_rawspec += cliquesize_counter[i] * i

        writeLog(filenames["logfilename"], [i, cliquesize_counter[i], sum, cliquesize_counter[i] * i, sum_rawspec])
    accession_stem = settings["accession_stem"]
    dt_now = datetime.datetime.now()
    writeLog(filenames["logfilename"], ["Finished", accession_stem, dt_now.strftime('%Y.%m.%d.%H.%M.%S'), counters["numberofclique"]])

def readsmalltable(path, keywords = set(), filetype = "tab"):
    """Read a tab- or comma-delimited table file into a list of row dicts.

    Uses the standard ``csv`` module (a leading UTF-8 BOM in the header row is
    stripped). For very large files where ``csv``'s field-size limits become a
    problem, use :func:`readtable` instead.

    Args:
        path (str or Path): Path of the target file.
        keywords (set or list, optional): Column names to extract; must match
            the file's header row exactly. Columns not in ``keywords`` are
            dropped from each row's dict. Defaults to ``set()`` (no columns kept).
        filetype (str, optional): ``"csv"`` for comma-delimited, or any other
            value for tab-delimited. Defaults to ``"tab"``.

    Returns:
        list[dict]: One dict per data row (excluding the header), mapping each
            requested column name to its string value for that row.
    """
    path = str(path)
    outputlist = []
    with open(path,  encoding='UTF-8') as f:
        if filetype == "csv":
            reader = csv.reader(f)
        else:
            reader = csv.reader(f, delimiter='\t')
        header = []
        for i, row in enumerate(reader):
            if i == 0:
                header = [x.lstrip('\ufeff') for x in list(row)]
                continue
            current_length = len(row)
            header_length = len(header)
            if current_length > header_length:
                current_length = header_length
            tempdict = {}
            for i in range(current_length):
                label = header[i]
                content = row[i]
                if label in keywords:
                    tempdict[label]=content
            outputlist.append(tempdict)
    return(outputlist)

def readtable(path, keywords = set(), filetype = "tab"):
    """Read a tab- or comma-delimited table file into a list of row dicts, via pandas.

    Uses ``pandas.read_csv`` (with quoting disabled) instead of the standard
    ``csv`` module, making it suitable for larger files or files with fields
    that would exceed ``csv``'s field-size limits. See :func:`readsmalltable`
    for the lighter-weight, ``csv``-based equivalent.

    Args:
        path (str or Path): Path of the target file.
        keywords (set or list, optional): Column names to extract; must match
            the file's header row exactly. Columns not found in a given row
            are filled with ``"N/A"``. Defaults to ``set()`` (no columns kept).
        filetype (str, optional): ``"csv"`` for comma-delimited, or any other
            value for tab-delimited. Defaults to ``"tab"``.

    Returns:
        list[dict]: One dict per data row, mapping each requested column name
            to its value for that row (``"N/A"`` if the column was absent).
    """

    import pandas as pd

    path = str(path)
    if filetype == "csv":
        df = pd.read_csv(path, quoting=csv.QUOTE_NONE)
    else:
        df = pd.read_csv(path, sep='\t', quoting=csv.QUOTE_NONE)
    outputlist = []
    for i, pdrow in df.iterrows():
        tempdict ={}
        for label in keywords:

            if label in pdrow:
                content = pdrow[label]
            else:
                content = "N/A"

            tempdict[label]=content
        outputlist.append(tempdict)
    return(outputlist)