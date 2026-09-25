#-------------------------------------------------------------------------------
# Name:        ddatoolbox
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     07/10/2022
# Updated:     03/11/2023
# Copyright:   (c) FumioMatsuda 2022
# Licence:     MIT lisence
#-------------------------------------------------------------------------------
"""
massbankserver.py

Module for building searchable spectral databases from MassBank and GNPS data.

This module provides:
- Parsing of MassBank text files
- Parsing of GNPS JSON spectral data
- Conversion into efficient binary search structures
- Fast m/z-based lookup using indexed peak tables

Main Outputs
------------
- peaktable.bin:
    Indexed m/z lookup table for fast search

- binary_format.bin:
    Serialized spectral peak data

- MassBankRecords.pkl:
    Metadata dictionary for all records

Key Features
------------
- Efficient large-scale spectral database construction
- Support for polarity and ionization filtering
- Normalization of spectral intensities
- Fast similarity search preparation

Binary Structures
-----------------
- SpectraPeaks:
    Stores (ID, m/z, intensity)

- PeakTable:
    Stores (m/z index, pointer position)

Dependencies
------------
- numpy
- json
- pickle
- ctypes

Notes
-----
- Designed for high-speed spectral matching workflows
- Compatible with downstream similarity search algorithms
- Uses normalized intensities (L2 normalization)

Author
------
Fumio Matsuda

Created
-------
2022

Updated
-------
# 240913 Errors in MassBankSearch was corrected.

License
-------
MIT License
"""


import xml.etree.ElementTree as ET
import base64, struct,itertools, glob, re, shutil, io
import numpy as np
import zlib, json
from statistics import stdev
import os, csv, math, pickle, sys
import ctypes as c


class SpectraPeaks(c.Structure):
    """Fixed-size binary record for one spectral peak, used in ``binary_format.bin``.

    A ``ctypes.Structure`` subclass whose in-memory/on-disk layout is exactly
    68 bytes (60 + 4 + 4) per record, allowing peaks to be seeked to and read
    directly from the binary spectral database by :func:`getpeaks`.

    Attributes:
        id (bytes): Record (spectrum) identifier, fixed-width 60-byte field.
        mz (float): Peak m/z value (32-bit float).
        intensity (float): Peak intensity (32-bit float).
    """
    _fields_ = [
        ("id", c.c_char*60), # 60-byte fixed-width character field
        ("mz", c.c_float), # float
        ("intensity", c.c_float) # float
    ]


class PeakTable(c.Structure):
    """Fixed-size binary index record for fast m/z-based lookup, used in ``peaktable.bin``.

    A ``ctypes.Structure`` subclass whose in-memory/on-disk layout is exactly
    8 bytes (4 + 4) per record. ``peaktable.bin`` holds one such record per
    integer-scaled m/z bin, mapping it to the starting byte offset of that bin's
    peaks within ``binary_format.bin`` (see :func:`getpeaks`).

    Attributes:
        mz (int): Integer-scaled m/z bin index (32-bit int).
        position (int): Index (record position, not byte offset) into the sorted
            ``binary_format.bin`` peak list where this bin's peaks begin
            (32-bit int).
    """
    _fields_ = [
        ("mz", c.c_int32), # int
        ("position", c.c_int32) # int
    ]




def readMassBankfile(file, ionization = "ESI", MSTYPE = "MS2"):
    """Parse a single MassBank flat-file record into a lightweight dict.

    Reads the whitespace-delimited MassBank record fields relevant to spectral
    search (title, accession, precursor m/z, formula, InChIKey, SMILES, ion mode,
    precursor type, and peak list), stopping early (returning the record
    accumulated so far) if the record's MS type, fragmentation mode (HCD is
    excluded), or instrument ionization does not match the requested filters.

    Args:
        file (str): Path to the MassBank record text file.
        ionization (str, optional): Required ionization method substring to match
            against ``AC$INSTRUMENT_TYPE:`` (``"APCI"`` is always additionally
            accepted). Defaults to ``"ESI"``.
        MSTYPE (str, optional): Required MS type to match against
            ``AC$MASS_SPECTROMETRY: MS_TYPE`` (e.g. ``"MS2"``). Defaults to ``"MS2"``.

    Returns:
        dict: Parsed record with keys including ``"RECORD_TITLE"``,
            ``"ACCESSION"``, ``"PRECURSOR_M/Z"``, ``"FORMULA"``, ``"INCHIKEY"``,
            ``"SMILES"``, ``"ION_MODE"``, ``"PRECURSOR_TYPE"``, and ``"Spectra"``
            (dict of m/z to intensity). Fields not yet reached when parsing
            stopped early are absent.

    See Also:
        :func:`readMassBankfiles`: batch caller of this function.
    """

    massbankrecord = {"RECORD_TITLE":"", "Spectra": {}, "PRECURSOR_M/Z":{}}
    massbankrecord_temp = {"RECORD_TITLE":"", "Spectra": {}, "PRECURSOR_M/Z":{}}

    mode = "before"
    with open(file, 'r', encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=' ')
        for i, row in enumerate(reader):
            #print(i, row)
            if row[0] == '//':
                continue

            if row[0] == 'RECORD_TITLE:':

                temptxt = " ".join(row[1:])
                temptxt1 = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
                record_title = "".join(temptxt1)

                massbankrecord_temp["RECORD_TITLE"] = record_title
            if row[0] == 'ACCESSION:':
                massbankrecord_temp["ACCESSION"] = row[1]
            #
            #
            if row[0] == 'MS$FOCUSED_ION:' and row[1] == 'PRECURSOR_M/Z':
                massbankrecord_temp["PRECURSOR_M/Z"] = row[2]
            #
            #
            if row[0] == 'CH$FORMULA:':
                temptxt = row[1]
                temptxt = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
                massbankrecord_temp["FORMULA"] = "".join(temptxt)
            #
            # Get INCHIKEY
            #
            if row[0] == 'CH$LINK:' and row[1] == 'INCHIKEY':
                temptxt = row[2]
                temptxt = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
                massbankrecord_temp["INCHIKEY"] = "".join(temptxt)

            #
            # Get SMILES
            #
            if row[0] == 'CH$SMILES:':
                temptxt = row[1]
                temptxt = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
                massbankrecord_temp["SMILES"] = "".join(temptxt)

            #
            # Porality
            #
            if row[0] == 'AC$MASS_SPECTROMETRY:' and row[1] == 'ION_MODE':
                massbankrecord_temp["ION_MODE"] = row[2]
                #print(file, massbankrecord["FORMULA"])
            #
            # Ignore if different polarity
            #
            if row[0] == 'MS$FOCUSED_ION:' and row[1] == 'PRECURSOR_TYPE':
                temptxt = row[2]
                temptxt = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
                massbankrecord_temp["PRECURSOR_TYPE"] = "".join(temptxt)
            #
            # Ignore if different polarity
            #
            if row[0] == 'AC$MASS_SPECTROMETRY:' and row[1] == 'MS_TYPE':

                if not row[2] == MSTYPE:
                    #print("w", row[2])
                    break
            #
            # Ignore if different polarity
            #
            if row[0] == 'AC$MASS_SPECTROMETRY:' and row[1] == 'FRAGMENTATION_MODE':

                if row[2] == "HCD":
                    #print(massbankrecord_temp["ACCESSION"], row[2])
                    break

            #
            # Ignore if ION_MODE is not ionization
            #
            if row[0] == 'AC$INSTRUMENT_TYPE:':
                if not ionization in row[1]:
                    if not "APCI" in row[1]:
                        #print("w", row[1])
                        break



            if mode == "data":
                mz = float(row[2])
                intensity = int(row[4])
                massbankrecord_temp["Spectra"][mz] = intensity


            if row[0] == 'PK$PEAK:':
                mode = "data"
        else:
            massbankrecord = massbankrecord_temp

    return(massbankrecord)

def progress_bar(current, total, bar_length=40):
    """Print (or update in place) a simple ASCII progress bar to stdout.

    Redraws the bar on the same terminal line (using a carriage return) each time
    it is called, so it should be called repeatedly with an increasing ``current``
    value from within a loop.

    Args:
        current (int): Number of items completed so far.
        total (int): Total number of items expected.
        bar_length (int, optional): Width of the bar in characters. Defaults to ``40``.

    Returns:
        None
    """
    fraction = current / total
    arrow = int(fraction * bar_length - 1) * '=' + '>'
    padding = int(bar_length - len(arrow)) * ' '

    # '\r' moves the cursor back to the start of the line; end='' avoids a newline
    print(f'\r[{arrow}{padding}] {int(fraction*100)}%', end='', file=sys.stdout)
    if current == total:
        print() # print a final newline once complete



def readMassBankfiles(path, targetpath, digit = 1, ionization = "ESI", polarity = "POSITIVE", MSTYPE = "MS2"):
    """Batch-parse MassBank record folders into a fast binary spectral search database.

    Recursively globs ``path`` for dataset subfolders, parses every record with
    :func:`readMassBankfile`, keeps only records matching ``polarity``, L2-
    normalizes each record's spectral intensities, and writes three output files
    to ``targetpath`` (recreated from scratch):

    - ``peaktable.bin``: A ``PeakTable`` index array, one entry per integer-scaled
      m/z bin, mapping to the starting index of that bin's peaks in
      ``binary_format.bin``.
    - ``binary_format.bin``: A ``SpectraPeaks`` array of all peaks
      (``id``, ``mz``, ``intensity``), sorted by ascending m/z.
    - ``MassBankRecords.pkl``: A pickled dict of full per-record metadata (as
      returned by :func:`readMassBankfile`), keyed by record ID.

    Args:
        path (str): Glob pattern for the folder(s) of MassBank record files to read.
        targetpath (str): Output folder for the three generated files (deleted and
            recreated if it already exists).
        digit (int, optional): Decimal resolution of the m/z index in
            ``peaktable.bin``. E.g. ``1`` indexes bins ``100.0, 100.1, 100.2, ...``;
            ``2`` indexes ``100.00, 100.01, 100.02, ...``. Larger values increase
            index file size but speed up search. Defaults to ``1``.
        ionization (str, optional): Required ionization method substring, passed
            through to :func:`readMassBankfile`. Defaults to ``"ESI"``.
        polarity (str, optional): Ion mode to keep, ``"POSITIVE"`` or ``"NEGATIVE"``.
            Records with a different (or missing) ``ION_MODE`` are skipped.
            Defaults to ``"POSITIVE"``.
        MSTYPE (str, optional): Required MS type, passed through to
            :func:`readMassBankfile`. Defaults to ``"MS2"``.

    Returns:
        None

    See Also:
        :func:`readMassBankfile`: per-record parser called by this function.
        :func:`searchMassBank`: consumes the database files produced here.

    Note:
        ⚠ Known issue: each ``PeakTable(mz, start)`` entry written to
        ``peaktable.bin`` is supposed to store its own bin index in ``mz``, but
        the code actually reuses a stale ``mz`` variable left over from an
        unrelated earlier loop (the last spectral peak of the last record
        processed), not the loop's own index variable. Every entry in
        ``peaktable.bin`` therefore ends up with the same, meaningless ``mz``
        field. This is harmless in practice because :func:`getpeaks` (the only
        reader of ``peaktable.bin``) never reads that field back, only
        ``position`` -- but the field should not be relied upon for anything.
    """



    data = []
    folders = glob.glob(path)
    recorddata = {}
    #print(folders)

    for folder in folders:
        if not os.path.isdir(folder):
            print(folder, "is not dir")
            continue
        print("Loading files in ", folder+"\\*")
        files = glob.glob(folder+"\\*")
        #print(files)

        total = len(files)
        current = 0

        for file in files:

            progress_bar(current, total, bar_length=40)
            current += 1

            data_id = file.split('\\')
            data_id = data_id[-1]
            data_id = data_id.replace(".txt", "")
            data_id = data_id.replace("MSBNK-", "")
            #print("loading file ", file, folder, data_id)
            massbankrecord = readMassBankfile(file, ionization, MSTYPE)
            #print(massbankrecord)
            if not "ION_MODE" in massbankrecord:
                continue
            if not massbankrecord["ION_MODE"] == polarity:
                continue
            recorddata[data_id] = massbankrecord


            spectra = {}
            squaresum = 0
            for mz, intensity in massbankrecord["Spectra"].items():
                squaresum += intensity**2
            for mz, intensity in massbankrecord["Spectra"].items():
                spectra[mz] = intensity/math.sqrt(squaresum)
            for mz, intensity in spectra.items():
                data.append([data_id, mz, intensity])

    data_sorted = sorted(data, key=lambda x: x[1])

    map_step1 = [0]*2000 * (10 ** digit)
    map_step2 = map_step1[:]
    for i in range(len(data_sorted)):
        peak = data_sorted[i]
        mz = int(round(peak[1], digit)* (10 ** digit)) # round at digit and multipyed into integer
        if mz > (len(map_step1)-10):
            continue
        map_step1[mz] = i+1

    pointer = 0

    for i, start in enumerate(map_step1):
        if pointer < start:
            pointer = start
        map_step2[i] = pointer

    table_data = []
    # Populate the binary record struct with values
    for i, start in enumerate(map_step2):
        table_data.append(PeakTable(mz, start))

    # Preparation of folders
    #
    if os.path.isdir(targetpath):
        shutil.rmtree(targetpath)
    os.mkdir(targetpath)

    # Write out the binary file
    with open(targetpath +"peaktable.bin", "bw") as f:
        for table_dataa in table_data:
            f.write(table_dataa)

    format_data = []
    # Populate the binary record struct with values
    for [spid, mz, intensity] in data_sorted:
        #print(spid, mz, intensity)
        format_data.append(SpectraPeaks(spid.encode(), mz, intensity))

    # Write out the binary file
    with open(targetpath +"binary_format.bin", "bw") as f:
        for format_dataa in format_data:
            f.write(format_dataa)

    with open(targetpath +"MassBankRecords.pkl", "wb") as tf:
        pickle.dump(recorddata,tf)


    return()




def readGNPSjsonrecord(jsondict, ionization = "ESI", MSTYPE = "2"):
    """Convert one GNPS JSON spectrum record into the common MassBank-like record dict.

    Filters out records whose ``Ion_Source`` does not match ``ionization``
    (``"APCI"`` is always additionally accepted) or whose ``ms_level`` does not
    match ``MSTYPE``, returning an empty-spectrum placeholder record in that case.
    Otherwise extracts title, instrument, accession, precursor m/z, formula,
    InChIKey, SMILES, ion mode, and adduct/precursor-type fields, and parses the
    GNPS ``peaks_json`` string into an m/z-to-intensity spectrum (dropping
    non-positive-intensity peaks).

    Args:
        jsondict (dict): One GNPS spectrum record, as parsed from the GNPS JSON
            export (one element of the top-level JSON array).
        ionization (str, optional): Required ionization source substring to match
            against ``jsondict["Ion_Source"]``. Defaults to ``"ESI"``.
        MSTYPE (str, optional): Required MS level to match against
            ``jsondict["ms_level"]``. Defaults to ``"2"``.

    Returns:
        dict: Parsed record with keys ``"RECORD_TITLE"``, ``"INSTRUMENT"``,
            ``"ACCESSION"``, ``"PRECURSOR_M/Z"``, ``"FORMULA"``, ``"INCHIKEY"``,
            ``"SMILES"``, ``"ION_MODE"``, ``"PRECURSOR_TYPE"``, and ``"Spectra"``
            (dict of m/z to intensity). When ``ionization``/``MSTYPE`` do not
            match, returns the empty placeholder
            ``{"RECORD_TITLE": "", "Spectra": {}, "PRECURSOR_M/Z": {}}``.

    See Also:
        :func:`readGNPSjsonfile`: batch caller of this function.
    """

    massbankrecord = {"RECORD_TITLE":"", "Spectra": {}, "PRECURSOR_M/Z":{}}

    Ion_Source = jsondict["Ion_Source"]

    if not ionization in Ion_Source:
        if not "APCI" in Ion_Source:
            return(massbankrecord)

    ms_level = jsondict["ms_level"]
    if not MSTYPE == ms_level:
        return(massbankrecord)

    #massbankrecord["RECORD_TITLE"] = jsondict["Compound_Name"]

    temptxt = jsondict["Compound_Name"]
    temptxt2 = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
    massbankrecord["RECORD_TITLE"] = "".join(temptxt2)

    massbankrecord["INSTRUMENT"] = jsondict['Instrument']


    massbankrecord["ACCESSION"] = jsondict['spectrum_id']
    massbankrecord["PRECURSOR_M/Z"] = float(jsondict["Precursor_MZ"])
    massbankrecord["FORMULA"] = jsondict["Formula_smiles"]
    #massbankrecord["INCHIKEY"] = jsondict["InChIKey_smiles"]

    temptxt = jsondict["InChIKey_smiles"]
    temptxt = re.findall(r'[A-Za-z0-9.,()/@ -=\[\]]', temptxt)
    massbankrecord["INCHIKEY"] = "".join(temptxt)


    massbankrecord["SMILES"] = jsondict["Smiles"]
    #print(massbankrecord["ACCESSION"], massbankrecord["RECORD_TITLE"], massbankrecord["INCHIKEY"], jsondict["Compound_Name"],  massbankrecord["PRECURSOR_M/Z"])

    # Porality
    massbankrecord["ION_MODE"] = jsondict["Ion_Mode"].upper()
    massbankrecord["PRECURSOR_TYPE"] = jsondict["Adduct"]


    peaks_json = eval(jsondict["peaks_json"])

    for (mz, intensity) in peaks_json:
        mz = float(mz)
        intensity = float(intensity)
        if intensity <= 0:
            continue
        massbankrecord["Spectra"][mz] = intensity


    return(massbankrecord)




def readGNPSjsonfile(filename, targetpath, digit = 1, ionization = "ESI", polarity = "POSITIVE", MSTYPE = "MS2"):
    """Batch-parse a GNPS spectral-library JSON export into a binary search database.

    Loads the GNPS JSON array from ``filename``, parses every entry with
    :func:`readGNPSjsonrecord`, keeps only records matching ``polarity``,
    L2-normalizes each record's spectral intensities, and writes the same three
    binary database files as :func:`readMassBankfiles` (``peaktable.bin``,
    ``binary_format.bin``, ``MassBankRecords.pkl``) to ``targetpath``.

    Args:
        filename (str): Path to the GNPS JSON export file (e.g.
            ``ALL_GNPS_NO_PROPOGATED.json``).
        targetpath (str): Output folder for the three generated files.
        digit (int, optional): Decimal resolution of the m/z index in
            ``peaktable.bin`` (see :func:`readMassBankfiles`). Defaults to ``1``.
        ionization (str, optional): Required ionization source substring, passed
            through to :func:`readGNPSjsonrecord`. Defaults to ``"ESI"``.
        polarity (str, optional): Ion mode to keep, ``"POSITIVE"`` or ``"NEGATIVE"``.
            Defaults to ``"POSITIVE"``.
        MSTYPE (str, optional): Required MS level, passed through to
            :func:`readGNPSjsonrecord`. Defaults to ``"MS2"``.

    Returns:
        None

    See Also:
        :func:`readGNPSjsonrecord`: per-record parser called by this function.
        :func:`searchMassBank`: consumes the database files produced here.

    Note:
        Internally calls :func:`readGNPSjsonrecord` with hard-coded
        ``ionization="ESI", MSTYPE="2"`` rather than forwarding this function's
        own ``ionization``/``MSTYPE`` arguments, so those two parameters
        currently have no effect on record-level filtering.
    """
    data = []
    recorddata = {}
    with open(filename) as f:
        di = json.load(f)

    for jsondict in di:

        data_id = jsondict['spectrum_id']
        print(data_id)

        massbankrecord = readGNPSjsonrecord(jsondict, ionization = "ESI", MSTYPE = "2")

        if not "ION_MODE" in massbankrecord:
            continue
        if not massbankrecord["ION_MODE"] == polarity:
            continue
        recorddata[data_id] = massbankrecord


        spectra = {}
        squaresum = 0
        for mz, intensity in massbankrecord["Spectra"].items():
            squaresum += intensity**2
        for mz, intensity in massbankrecord["Spectra"].items():
            spectra[mz] = intensity/math.sqrt(squaresum)
        for mz, intensity in spectra.items():
            data.append([data_id, mz, intensity])

    data_sorted = sorted(data, key=lambda x: x[1])

    map_step1 = [0]*2000 * (10 ** digit)
    map_step2 = map_step1[:]
    for i in range(len(data_sorted)):
        peak = data_sorted[i]
        mz = int(round(peak[1], digit)* (10 ** digit)) # round at digit and multipyed into integer
        if mz > (len(map_step1)-10):
            continue
        map_step1[mz] = i+1

    pointer = 0

    for i, start in enumerate(map_step1):
        if pointer < start:
            pointer = start
        map_step2[i] = pointer

    table_data = []
    # Populate the binary record struct with values
    for i, start in enumerate(map_step2):
        table_data.append(PeakTable(mz, start))

    # Preparation of folders
    #
    if os.path.isdir(targetpath):
        shutil.rmtree(targetpath)
    os.mkdir(targetpath)
    
    # Write out the binary file
    with open(targetpath +"peaktable.bin", "bw") as f:
        for table_dataa in table_data:
            f.write(table_dataa)

    format_data = []
    # Populate the binary record struct with values
    for [spid, mz, intensity] in data_sorted:
        #print(spid, mz, intensity)
        format_data.append(SpectraPeaks(spid.encode(), mz, intensity))

    # Write out the binary file
    with open(targetpath +"binary_format.bin", "bw") as f:
        for format_dataa in format_data:
            f.write(format_dataa)

    with open(targetpath +"MassBankRecords.pkl", "wb") as tf:
        pickle.dump(recorddata,tf)


    return()




def getpeaks(mz_query, mzthreshold, digit, path, fp, fb):
    """Fetch all indexed binary-database peaks within a threshold of a query m/z.

    Uses ``peaktable.bin`` (via the already-open file handle ``fp``) to look up
    the byte range in ``binary_format.bin`` (via ``fb``) covering
    ``[mz_query - mzthreshold, mz_query + mzthreshold]``, reads that whole chunk
    into memory in one call, and unpacks it into individual ``SpectraPeaks``
    records, filtering out any that fall outside the exact m/z window.

    Args:
        mz_query (float): Query m/z value to search around.
        mzthreshold (float): m/z tolerance (in the same unit as the index,
            typically Da) on either side of ``mz_query``.
        digit (int): Decimal resolution used when the index was built (must match
            the ``digit`` passed to :func:`readMassBankfiles`/:func:`readGNPSjsonfile`).
        path (str): Path to the folder containing ``peaktable.bin`` and
            ``binary_format.bin`` (unused directly here; kept for interface
            consistency since ``fp``/``fb`` are the actual open handles).
        fp (BinaryIO): Open file handle to ``peaktable.bin`` (read-binary mode).
        fb (BinaryIO): Open file handle to ``binary_format.bin`` (read-binary mode).

    Returns:
        list[list]: Hit peaks, each as ``[record_id (str), mz (float), intensity (float)]``.
            Returns ``[]`` immediately when ``mzthreshold == 0`` or when
            ``mz_query`` rounds to bin ``0``.

    See Also:
        :func:`searchMassBank`: caller of this function.
    """
    if mzthreshold == 0:
        return([])


    targetmz = int(round(mz_query, digit) * (10**digit)) # round at digit and multipyed into integer

    if targetmz == 0:
        return([])
    targetthreshold = int(round(mzthreshold, digit) * (10**digit))

    peaks = []
    #with open(path+"peaktable.bin", "rb") as f:
        # Create an instance of the defined binary struct
    format_data = PeakTable()
    # Read from the file into the struct instance
    #print(targetmz, targetthreshold)

    fp.seek( 8 * (targetmz - targetthreshold  -1), os.SEEK_SET)
    fp.readinto(format_data) # Load bytes into the pre-allocated struct instance

    #print(targetmz, f'mz: {format_data.mz}', f'position: {format_data.position}') # 100
    target_position_start = format_data.position

    fp.seek( 8 * (targetmz + targetthreshold ), os.SEEK_SET)
    fp.readinto(format_data) # Load bytes into the pre-allocated struct instance

    #print(targetmz, f'mz: {format_data.mz}', f'position: {format_data.position}') # 100
    target_position_end = format_data.position



    #print(target_position_start, target_position_end)
    hitted = {}
    #with open(path+"binary_format.bin", "rb") as f:
    # Create an instance of the defined binary struct
    format_data = SpectraPeaks()
    # Read from the file into the struct instance

    fb.seek(68 * target_position_start, os.SEEK_SET)
    chunk = fb.read(68 * (target_position_end - target_position_start))
    for position in range((target_position_end - target_position_start)):
        buffer = io.BytesIO(chunk[position*68: position*68 + 68])
        buffer.readinto(format_data)
        sid = format_data.id.decode()
        mz = format_data.mz
        if (mz_query - mz) > mzthreshold:
            continue
        if (mz - mz_query) > mzthreshold:
            break
        peaks.append([sid,mz, format_data.intensity])

    """
    for step in range(target_position_end - target_position_start):
        fb.seek( 68 * (target_position_start + step), os.SEEK_SET)

        fb.readinto(format_data) # Load bytes into the pre-allocated struct instance

        #print(f'id: {format_data.id}') # 100
        #print(f'type: {format_data.mz}') # b'DATA'
        #print(f'data: {format_data.intensity}') # sample binary data.
        sid = format_data.id.decode()
        #if sid in hitted:
        #    continue
        #hitted[sid] = 1
        mz = format_data.mz
        if (mz_query - mz) > mzthreshold:
            continue
        if (mz - mz_query) > mzthreshold:
            break

        peaks.append([sid, format_data.mz, format_data.intensity])
    """
    return(peaks)



def searchMassBank(new_dict, query, scorethreshold = 0.7, mzthreshold = 0.1, digit = 1, path = "./MassBankpositive/", unit = "Da"):
    """Search a query spectrum against the binary MassBank/GNPS spectral database.

    L2-normalizes ``query``, then for each query peak (processed in descending
    intensity order) looks up nearby database peaks via :func:`getpeaks`,
    accumulating a running dot-product per matched record while ensuring each
    database peak and each query peak is only counted once per record (greedy
    one-to-one matching within a record). The final cosine-product score for a
    record is the squared dot product (Stein's method); records above
    ``scorethreshold`` are kept and enriched with metadata from ``new_dict``.

    Args:
        new_dict (dict): Dict of MassBank/GNPS record metadata, as returned by
            :func:`setupMassBank`.
        query (dict): Query mass spectrum, mapping m/z (float) to intensity
            (int or float).
        scorethreshold (float, optional): Minimum cosine-product score (already
            squared) required to keep a hit. Defaults to ``0.7``.
        mzthreshold (float, optional): m/z matching tolerance, in the unit given
            by ``unit``. Defaults to ``0.1``.
        digit (int, optional): Decimal resolution used when the index was built
            (must match the ``digit`` passed to :func:`readMassBankfiles`/
            :func:`readGNPSjsonfile`). Defaults to ``1``.
        path (str, optional): Path to the folder containing ``peaktable.bin`` and
            ``binary_format.bin``. Defaults to ``"./MassBankpositive/"``.
        unit (str, optional): Unit of ``mzthreshold``, either ``"Da"`` (absolute)
            or ``"ppm"`` (relative to each query m/z). Defaults to ``"Da"``.

    Returns:
        dict: Mapping of record ID (str) to a dict with keys ``"score"`` (float,
            cosine-product score), ``"dice"`` (float, Dice coefficient),
            ``"RECORD_TITLE"``, ``"FORMULA"``, ``"ACCESSION"``, ``"hittedmz"``
            (list of matched query m/z values), ``"spectra"`` (the record's full
            spectrum), ``"forentoropy"`` (list of summed intensities per match,
            for downstream entropy-similarity calculations), ``"INCHIKEY"``
            (``""`` if absent), and ``"PRECURSOR_M/Z"`` (float, ``0`` if
            unavailable).

    See Also:
        :func:`readGNPSjsonfile`, :func:`readMassBankfiles`: build the database
        searched here. :func:`setupMassBank`: loads ``new_dict``.
    """
    spectra = {}
    results = {}
    output = {}
    squaresum = 0

    resultsforentoropy = {}
    for mz, intensity in query.items():
        if mz > 1950:
            continue
        squaresum += intensity**2
    for mz, intensity in query.items():
        if mz > 1950:
            continue

        spectra[mz] = intensity/math.sqrt(squaresum)

    spectra_keys_sorted = sorted(spectra.keys(), key = lambda x: spectra[x], reverse = True)


    hittothisMassBankpeak ={}
    hittedmz = {}


    with open(path+"peaktable.bin", "rb") as fp:
        with open(path+"binary_format.bin", "rb") as fb:
            for mz_query in spectra_keys_sorted:
                #print(mz_query)
                intensity_query = spectra[mz_query]
                if unit == "ppm":
                    mzthreshold_temp = mzthreshold /1000000 * mz_query
                else:
                    mzthreshold_temp = mzthreshold
                peaks = getpeaks(mz_query, mzthreshold_temp, digit, path, fp, fb)

                peaks_sorted = sorted(peaks, key = lambda x: x[2], reverse = True)

                #peaks_sorted = peaks

                hittothisquerypeak ={} #Check to avoid redundant hit

                for sid, mz, intensity in peaks_sorted:
                    if (sid, mz) in hittothisMassBankpeak:
                        continue
                    if sid in hittothisquerypeak:
                        continue
                    if not sid in results:
                        results[sid] = 0
                    if not sid in resultsforentoropy:
                        resultsforentoropy[sid] = []

                    results[sid] += intensity * intensity_query

                    resultsforentoropy[sid].append(intensity + intensity_query)

                    hittothisquerypeak[sid] = mz
                    hittothisMassBankpeak[(sid, mz)] = mz

                    if not sid in hittedmz:
                        hittedmz[sid] = []
                    hittedmz[sid].append(mz_query)


            #if "AC000114" in sid:
            #    print(sid, mz, mz_query, intensity, intensity_query)
    #
    # Sort result
    #
    results_selected = {}
    for sid, score in results.items():
        if score * score > scorethreshold:
            results_selected[sid] = score

    numoffrg_query = len(query)


    for sid in sorted(results_selected.keys(), key=lambda x: results[x], reverse=True):
        #
        #Stein style
        score = results[sid] * results[sid]
        #
        #
        numoffrg_hit = len(new_dict[sid]['Spectra'])
        numoffrg_common = len(hittedmz[sid])
        #print(numoffrg_query, numoffrg_hit, numoffrg_common)
        if score > scorethreshold:


            output[sid] = {}
            output[sid]["score"] = score
            output[sid]["dice"] = 2 * numoffrg_common/(numoffrg_query + numoffrg_hit)
            output[sid]["RECORD_TITLE"] = new_dict[sid]["RECORD_TITLE"]
            output[sid]["FORMULA"] = new_dict[sid]["FORMULA"]
            output[sid]["ACCESSION"] = new_dict[sid]["ACCESSION"]
            output[sid]["hittedmz"] = hittedmz[sid]
            output[sid]["spectra"] = new_dict[sid]["Spectra"]
            output[sid]["forentoropy"] =  resultsforentoropy[sid]

            if "INCHIKEY" in new_dict[sid]:
                output[sid]["INCHIKEY"] = new_dict[sid]["INCHIKEY"]
            else:
                output[sid]["INCHIKEY"] = ""
            if new_dict[sid]["PRECURSOR_M/Z"] == {}:
                output[sid]["PRECURSOR_M/Z"] = 0
            elif new_dict[sid]["PRECURSOR_M/Z"] == "NA":
                output[sid]["PRECURSOR_M/Z"] = 0
            else:
                output[sid]["PRECURSOR_M/Z"] = float(new_dict[sid]["PRECURSOR_M/Z"])
    return(output)




def getpeaks_old(mz_query, mzthreshold, digit, path):
    """Legacy (per-call file-open) implementation of :func:`getpeaks`.

    .. deprecated::
        Superseded by :func:`getpeaks`, which takes already-open file handles
        (``fp``, ``fb``) instead of re-opening ``peaktable.bin``/
        ``binary_format.bin`` on every call, and reads its byte range in one bulk
        read instead of per-record ``seek``/``readinto`` calls. Kept for
        reference; used only by the deprecated :func:`searchMassBank_old`.

    Args:
        mz_query (float): Query m/z value to search around.
        mzthreshold (float): m/z tolerance on either side of ``mz_query``.
        digit (int): Decimal resolution used when the index was built.
        path (str): Path to the folder containing ``peaktable.bin`` and
            ``binary_format.bin`` (opened fresh on every call).

    Returns:
        list[list]: Hit peaks, each as ``[record_id (str), mz (float), intensity (float)]``.
            Returns ``[]`` immediately when ``mzthreshold == 0`` or when
            ``mz_query`` rounds to bin ``0``.

    See Also:
        :func:`searchMassBank_old`: caller of this function.
    """
    if mzthreshold == 0:
        return([])


    targetmz = int(round(mz_query, digit) * (10**digit)) # round at digit and multipyed into integer

    if targetmz == 0:
        return([])
    targetthreshold = int(round(mzthreshold, digit) * (10**digit))

    peaks = []
    with open(path+"peaktable.bin", "rb") as f:
        # Create an instance of the defined binary struct
        format_data = PeakTable()
        # Read from the file into the struct instance
        #print(targetmz, targetthreshold)

        f.seek( 8 * (targetmz - targetthreshold  -1), os.SEEK_SET)
        f.readinto(format_data) # Load bytes into the pre-allocated struct instance

        #print(targetmz, f'mz: {format_data.mz}', f'position: {format_data.position}') # 100
        target_position_start = format_data.position

        f.seek( 8 * (targetmz + targetthreshold ), os.SEEK_SET)
        f.readinto(format_data) # Load bytes into the pre-allocated struct instance

        #print(targetmz, f'mz: {format_data.mz}', f'position: {format_data.position}') # 100
        target_position_end = format_data.position



    #print(target_position_start, target_position_end)
    hitted = {}
    with open(path+"binary_format.bin", "rb") as f:
        # Create an instance of the defined binary struct
        format_data = SpectraPeaks()
        # Read from the file into the struct instance

        for step in range(target_position_end - target_position_start):
            f.seek( 68 * (target_position_start + step), os.SEEK_SET)

            f.readinto(format_data) # Load bytes into the pre-allocated struct instance

            #print(f'id: {format_data.id}') # 100
            #print(f'type: {format_data.mz}') # b'DATA'
            #print(f'data: {format_data.intensity}') # sample binary data.
            sid = format_data.id.decode()
            #if sid in hitted:
            #    continue
            #hitted[sid] = 1
            mz = format_data.mz
            if (mz_query - mz) > mzthreshold:
                continue
            if (mz - mz_query) > mzthreshold:
                continue

            peaks.append([sid, format_data.mz, format_data.intensity])
    return(peaks)



def searchMassBank_old(new_dict, query, scorethreshold = 0.7, mzthreshold = 0.1, digit = 1, path = "./MassBankpositive/", unit = "Da"):
    """Legacy (per-peak file-open) implementation of :func:`searchMassBank`.

    .. deprecated::
        Superseded by :func:`searchMassBank`, which opens ``peaktable.bin``/
        ``binary_format.bin`` once and reuses the handles across all query peaks
        (via :func:`getpeaks`) instead of re-opening them for every peak (via
        :func:`getpeaks_old`), and additionally tracks per-hit intensity data for
        downstream entropy-similarity scoring. Kept for reference.

    Args:
        new_dict (dict): Dict of MassBank/GNPS record metadata, as returned by
            :func:`setupMassBank`.
        query (dict): Query mass spectrum, mapping m/z (float) to intensity
            (int or float).
        scorethreshold (float, optional): Minimum cosine-product score (already
            squared) required to keep a hit. Defaults to ``0.7``.
        mzthreshold (float, optional): m/z matching tolerance, in the unit given
            by ``unit``. Defaults to ``0.1``.
        digit (int, optional): Decimal resolution used when the index was built.
            Defaults to ``1``.
        path (str, optional): Path to the folder containing ``peaktable.bin`` and
            ``binary_format.bin``. Defaults to ``"./MassBankpositive/"``.
        unit (str, optional): Unit of ``mzthreshold``, either ``"Da"`` (absolute)
            or ``"ppm"`` (relative to each query m/z). Defaults to ``"Da"``.

    Returns:
        dict: Mapping of record ID (str) to a dict with keys ``"score"``,
            ``"dice"``, ``"RECORD_TITLE"``, ``"FORMULA"``, ``"ACCESSION"``,
            ``"hittedmz"``, ``"INCHIKEY"``, and ``"PRECURSOR_M/Z"`` (same meaning
            as in :func:`searchMassBank`, minus the ``"spectra"``/``"forentoropy"``
            keys).

    See Also:
        :func:`searchMassBank`: current implementation.

    Raises:
        TypeError: Always. This function calls ``getpeaks(mz_query,
            mzthreshold_temp, digit, path)`` with only 4 positional arguments,
            but :func:`getpeaks` requires 6 (it also needs open file handles
            ``fp``/``fb``, unlike :func:`getpeaks_old`, which this function was
            presumably meant to call instead). As written, ``searchMassBank_old``
            is dead/broken code that cannot run to completion; it is kept only
            for reference alongside :func:`getpeaks_old`.
    """
    spectra = {}
    results = {}
    output = {}
    squaresum = 0
    for mz, intensity in query.items():
        if mz > 1950:
            continue
        squaresum += intensity**2
    for mz, intensity in query.items():
        if mz > 1950:
            continue
        spectra[mz] = intensity/math.sqrt(squaresum)
    spectra_keys_sorted = sorted(spectra.keys(), key = lambda x: spectra[x], reverse = True)


    hittothisMassBankpeak ={}
    hittedmz = {}

    for mz_query in spectra_keys_sorted:
        #print(mz_query)
        intensity_query = spectra[mz_query]


        if unit == "ppm":
            mzthreshold_temp = mzthreshold /1000000 * mz_query
        else:

            mzthreshold_temp = mzthreshold


        peaks = getpeaks(mz_query, mzthreshold_temp, digit, path)


        peaks_sorted = sorted(peaks, key = lambda x: x[2], reverse = True)

        hittothisquerypeak ={} #Check to avoid redundant hit

        for sid, mz, intensity in peaks_sorted:
            if (sid, mz) in hittothisMassBankpeak:
                continue
            if sid in hittothisquerypeak:
                continue
            if not sid in results:
                results[sid] = 0
            results[sid] += intensity * intensity_query
            hittothisquerypeak[sid] = mz
            hittothisMassBankpeak[(sid, mz)] = mz

            if not sid in hittedmz:
                hittedmz[sid] = []
            hittedmz[sid].append(mz_query)


            #if "AC000114" in sid:
            #    print(sid, mz, mz_query, intensity, intensity_query)
    #
    # Sort result
    #
    results_selected = {}
    for sid, score in results.items():
        if score * score > scorethreshold:
            results_selected[sid] = score

    numoffrg_query = len(query)


    for sid in sorted(results_selected.keys(), key=lambda x: results[x], reverse=True):
        #
        #Stein style
        score = results[sid] * results[sid]
        #
        #
        numoffrg_hit = len(new_dict[sid]['Spectra'])
        numoffrg_common = len(hittedmz[sid])
        #print(numoffrg_query, numoffrg_hit, numoffrg_common)
        if score > scorethreshold:


            output[sid] = {}
            output[sid]["score"] = score
            output[sid]["dice"] = 2 * numoffrg_common/(numoffrg_query + numoffrg_hit)
            output[sid]["RECORD_TITLE"] = new_dict[sid]["RECORD_TITLE"]
            output[sid]["FORMULA"] = new_dict[sid]["FORMULA"]
            output[sid]["ACCESSION"] = new_dict[sid]["ACCESSION"]
            output[sid]["hittedmz"] = hittedmz[sid]
            if "INCHIKEY" in new_dict[sid]:
                output[sid]["INCHIKEY"] = new_dict[sid]["INCHIKEY"]
            else:
                output[sid]["INCHIKEY"] = ""
            if new_dict[sid]["PRECURSOR_M/Z"] == {}:
                output[sid]["PRECURSOR_M/Z"] = 0
            elif new_dict[sid]["PRECURSOR_M/Z"] == "NA":
                output[sid]["PRECURSOR_M/Z"] = 0
            else:
                output[sid]["PRECURSOR_M/Z"] = float(new_dict[sid]["PRECURSOR_M/Z"])
    return(output)



def getRecordData(recordid, path = "./MassBankpositive/"):
    """Look up a single record's metadata from the ``MassBankRecords.pkl`` database.

    Args:
        recordid (str): MassBank/GNPS record identifier to look up.
        path (str, optional): Path to the folder containing
            ``MassBankRecords.pkl`` (produced by :func:`readGNPSjsonfile` or
            :func:`readMassBankfiles`). Defaults to ``"./MassBankpositive/"``.

    Returns:
        dict: Metadata for ``recordid`` (same shape as one value returned by
            :func:`readMassBankfile`/:func:`readGNPSjsonrecord`).

    Raises:
        KeyError: If ``recordid`` is not present in the database.

    See Also:
        :func:`readGNPSjsonfile`, :func:`readMassBankfiles`: build the database.
        :func:`setupMassBank`: loads the whole database at once.
    """
    with open(path+"MassBankRecords.pkl", "rb") as tf:
        new_dict = pickle.load(tf)
    return(new_dict[recordid])


def setupMassBank(path = "./MassBankpositive/"):
    """Load the full ``MassBankRecords.pkl`` metadata database into memory.

    Args:
        path (str, optional): Path to the folder containing
            ``MassBankRecords.pkl`` (produced by :func:`readGNPSjsonfile` or
            :func:`readMassBankfiles`). Defaults to ``"./MassBankpositive/"``.

    Returns:
        dict: Mapping of record ID to record metadata, ready to be passed as
            ``new_dict`` to :func:`searchMassBank`.

    See Also:
        :func:`readGNPSjsonfile`, :func:`readMassBankfiles`: build the database.
        :func:`searchMassBank`: consumes the loaded database.
    """
    with open(path+"MassBankRecords.pkl", "rb") as tf:
        new_dict = pickle.load(tf)
    return(new_dict)
