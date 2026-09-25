#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""
massbankio.py

Module for reading, writing, and manipulating MassBank records.

This module defines the MassBankRecord class, which provides:
- Structured storage of MassBank-formatted metadata
- Parsing of MassBank text files
- Template-based record generation
- Spectral data handling and normalization
- Annotation generation and parsing

Key Features
------------
- Full compatibility with MassBank record format
- Support for iterative and unique fields
- Automated formatting and validation
- Flexible annotation system for processed spectra

Data Structure
--------------
MassBank records are internally stored as a dictionary:
    {(field, subfield): metadata_dict}

Spectral data:
    PK$PEAK → {m/z: relative intensity}

Annotations:
    PK$ANNOTATION → list of annotated peak entries

Dependencies
------------
- numpy
- scipy
- csv
- re

Notes
-----
- This module is designed for integration with DDA processing pipelines
- Supports both MS1 and MS2 record construction
- Includes utility functions for formatting MassBank-compliant outputs

Author
------
Fumio Matsuda

Created
-------
2022

License
-------
<Specify License>
"""
import csv, re, datetime, math, itertools, copy
import numpy as np
from scipy.stats import spearmanr, norm
from importlib.resources import read_text
from splash import Spectrum, SpectrumType, Splash

class MassBankRecord:
    """In-memory representation of one MassBank record, with read/write helpers.

    Backs a single MassBank flat-file record (MS1 or MS2) with a dictionary keyed
    by ``(field, subfield)`` tuples (e.g. ``("COMMENT:", "CONFIDENCE")``), built
    from the field definitions in the bundled ``recordformat.txt`` template. Each
    entry tracks whether the field is ``"unique"`` (single value) or ``"iterative"``
    (list of values, e.g. multiple ``COMMENT:`` lines), plus its required-ness,
    print format, and validation rule. Spectral peaks (``PK$PEAK:``) and peak
    annotations (``PK$ANNOTATION:``) are stored separately as dict/list data under
    their own record entries.

    Provides methods to populate a record from a dataset's ``RecordTemplate.txt``
    and/or an existing MassBank record file, to add/clear individual field values,
    to attach and manipulate spectral data, to run isotope/charge-state/ion-form
    checks against MS1 spectra, and to serialize the record back out to a
    MassBank-formatted text file.

    Attributes:
        record (dict): ``(fld1, fld2) -> {"order", "req", "uq", "fld1", "fld2",
            "data", "format", "check"}`` (plus ``"spectra"`` for the peak-table
            entry). The core record storage.
        MASS_SPECTROMETRY (list[str]): Known ``AC$MASS_SPECTROMETRY:`` subfield
            names, from the record-format template.
        LINK (list[str]): Known ``CH$LINK:`` subfield names.
        CHROMATOGRAPHY (list[str]): Known ``AC$CHROMATOGRAPHY:`` subfield names.
        COMMENT (list[str]): Known ``COMMENT:`` subfield names.
        FOCUSED_ION (list[str]): Known ``MS$FOCUSED_ION:`` subfield names.
        DATA_PROCESSING (list[str]): Known ``MS$DATA_PROCESSING:`` subfield names.
    """

    def __init__(self):
        """Initialize an empty record by loading the bundled field-format template.

        Reads ``recordformat.txt`` (packaged with ``ddatoolbox``) and builds
        :attr:`record` with one empty entry per defined ``(fld1, fld2)`` field,
        plus the per-category subfield-name lists (:attr:`MASS_SPECTROMETRY`,
        :attr:`LINK`, :attr:`CHROMATOGRAPHY`, :attr:`COMMENT`, :attr:`FOCUSED_ION`,
        :attr:`DATA_PROCESSING`) used elsewhere to validate/dispatch field writes.
        """

        self.record = {}
        self.MASS_SPECTROMETRY = []
        self.LINK = []
        self.CHROMATOGRAPHY = []
        self.COMMENT = []
        self.FOCUSED_ION = []
        self.DATA_PROCESSING = []

        reader = read_text("ddatoolbox", "recordformat.txt", encoding="utf-8").split('\n')
        for i, rowtex in enumerate(reader):
            row = rowtex.split('\t')
            if len(row) < 6:
                continue
            #print(i, row)
            fld1 = row[0]
            fld2 = row[1]
            fld1.replace(" ", "")
            fld2.replace(" ", "")
            req  = row[2]
            uq   = row[3]
            format = row[4]
            format_check = row[5]
            thiskey= (fld1, fld2)
            if uq == "iterative":
                self.record[thiskey] = {"order": i, "req": req, "uq": uq, "fld1": fld1, "fld2": fld2, "data":[], "format":format, "check":format_check}
            elif uq == "unique":
                self.record[thiskey] = {"order": i, "req": req, "uq": uq, "fld1": fld1, "fld2": fld2, "data":"", "format":format, "check":format_check}
            elif uq == "text" and fld1 == "PK$PEAK:":
                self.record[thiskey] = {"order": i, "req": req, "uq": uq, "fld1": fld1, "fld2": fld2, "data":{}, "spectra":{}}
            elif uq == "text" and fld1 == "PK$ANNOTATION:":
                self.record[thiskey] = {"order": i, "req": req, "uq": uq, "fld1": fld1, "fld2": fld2, "data":[]}
            if fld1 == "AC$MASS_SPECTROMETRY:":
                self.MASS_SPECTROMETRY.append(fld2)
            if fld1 == "AC$CHROMATOGRAPHY:":
                self.CHROMATOGRAPHY.append(fld2)
            if fld1 == "COMMENT:":
                if fld2 != "*":
                    self.COMMENT.append(fld2)
            if fld1 == "MS$FOCUSED_ION:":
                self.FOCUSED_ION.append(fld2)
            if fld1 == "CH$LINK:" and fld2 != "*":
                self.LINK.append(fld2)
            if fld1 == "MS$DATA_PROCESSING:" and fld2 != "*":
                self.DATA_PROCESSING.append(fld2)

    def read_template(self, deffile, msx = "ms2"):
        """Populate dataset-wide constant fields from a ``RecordTemplate.txt`` file.

        Reads the tab-separated ``(field, subfield, value)`` rows of a dataset's
        template file and writes each value into the matching record field via
        :meth:`add2`. Rows with ``fld1 == "PARAMETER:"`` are skipped (they
        configure the processing pipeline, not the MassBank record itself). A
        historical naming mistake is patched on the fly: ``SOLVENT A``/``SOLVENT B``
        subfields are merged into a single ``SOLVENT`` field prefixed with
        ``"A "``/``"B "``. When ``msx == "ms1"``, MS2-only fields
        (``FRAGMENTATION_MODE``, ``COLLISION_ENERGY``, ``COLLISION_GAS`` under
        ``AC$MASS_SPECTROMETRY:``) are skipped.

        Args:
            deffile (str): Path to the dataset's ``RecordTemplate.txt`` file.
            msx (str, optional): ``"ms2"`` or ``"ms1"``, controlling which fields
                are applicable. Defaults to ``"ms2"``.

        Returns:
            None
        """
        with open(deffile, 'r', encoding='UTF-8' ) as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader):
                if len(row) < 2:
                    continue
                #print(i, row)
                fld1 = row[0]
                fld2 = row[1]
                fld1.replace(" ", "")
                fld2.replace(" ", "")
                data  = row[2]

                #
                # Patch for mistake
                #
                if fld2 == "SOLVENT A":
                    fld2 = "SOLVENT"
                    data = "A "+data
                if fld2 == "SOLVENT B":
                    fld2 = "SOLVENT"
                    data = "B "+data

                thiskey= (fld1, fld2)
                #
                # Ignore "PARAMETER:"
                #
                if fld1 == "PARAMETER:":
                    continue
                #
                # Ignore following fields for MS1 data
                #
                if msx == "ms1":
                    if thiskey in [("AC$MASS_SPECTROMETRY:","FRAGMENTATION_MODE"), ("AC$MASS_SPECTROMETRY:","COLLISION_ENERGY"), ("AC$MASS_SPECTROMETRY:","COLLISION_GAS")]:
                        continue
                self.add2(fld1, fld2, data)

    def read_MassBankRecord(self, file):
        """Parse a MassBank-formatted text file into this record's field storage.

        Reads ``file`` line by line, dispatching each row to the appropriate
        record field via :meth:`add`/:meth:`add2` based on its ``fld1`` prefix
        (``AC$MASS_SPECTROMETRY:``, ``AC$CHROMATOGRAPHY:``, ``MS$FOCUSED_ION:``,
        ``CH$LINK:``, ``COMMENT:``, ``MS$DATA_PROCESSING:``). Special handling is
        applied to ``PK$PEAK:`` (switches to spectra-reading mode, storing each
        subsequent ``m/z, intensity, relative-intensity`` row) and
        ``PK$ANNOTATION:`` (switches to annotation-reading mode, appending each
        subsequent row as a raw text line). Reading stops at the ``//`` record
        terminator.

        Args:
            file (str): Path to the MassBank record text file to read.

        Returns:
            None
        """
        with open(file, 'r', encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=' ')
            mode = "none"

            for i, row in enumerate(reader):
                if row[0] == "//":
                    break
                #print(row)
                fld1 = row[0]
                fld2 = row[1]
                thiskey= (fld1, fld2)

                if fld1 == 'AC$MASS_SPECTROMETRY:':
                    mode = "none"
                    for key2 in self.MASS_SPECTROMETRY:
                        data = " ".join(row[2:])
                        if key2 == fld2:

                            self.add2( fld1, fld2, data)
                            break
                    continue
                if fld1 == 'AC$CHROMATOGRAPHY:':
                    mode = "none"
                    for key2 in self.CHROMATOGRAPHY:
                        data = " ".join(row[2:])
                        if key2 == fld2:
                            self.add2(fld1, fld2, data)
                            break
                    continue
                if fld1 == 'MS$FOCUSED_ION:':
                    mode = "none"
                    for key2 in self.FOCUSED_ION:
                        data = " ".join(row[2:])
                        if key2 == fld2:
                            self.add2( fld1, fld2, data)
                            break
                    continue
                if fld1 == 'CH$LINK:':
                    mode = "none"
                    for key2 in self.LINK:
                        data = " ".join(row[2:])
                        if key2 == fld2:
                            self.add2( fld1, fld2, data)
                            break
                    continue
                if fld1 == 'COMMENT:':
                    mode = "none"
                    for key2 in self.COMMENT:
                        data = " ".join(row[2:])
                        if key2 == fld2:
                            self.add2( fld1, fld2, data)
                            break
                    else:
                        data = " ".join(row[1:])
                        self.add(fld1, data)
                    continue
                if fld1 == 'PK$PEAK:':
                    mode = "spectra"
                    continue
                if fld1 == 'PK$ANNOTATION:':
                    mode = "annotation"
                    data = " ".join(row[1:])
                    thiskey = ("PK$ANNOTATION:","*")
                    self.record[thiskey]["data"].append(data)
                    continue

                if fld1 == 'PK$NUM_PEAK:':
                    mode = "none"
                    data = " ".join(row[1:])
                    self.add(fld1, data)
                    continue

                if fld1 == 'MS$DATA_PROCESSING:':
                    mode = "none"
                    for key2 in self.DATA_PROCESSING:
                        data = " ".join(row[2:])
                        if key2 == fld2:
                            self.add2( fld1, fld2, data)
                            break
                    continue


                if mode == "spectra" and fld1 == "":
                    mz = row[2]
                    intensity = row[3]
                    relint = row[4]
                    thiskey = ("PK$PEAK:","*")
                    self.record[thiskey]["data"][mz] = relint
                    self.record[thiskey]["spectra"][mz] = intensity
                    continue
                if mode == "annotation" and fld1 == "":
                    data = " ".join(row[:])
                    thiskey = ("PK$ANNOTATION:","*")
                    self.record[thiskey]["data"].append(data)
                    continue

                data = " ".join(row[1:])
                self.add(fld1, data)




    def add2(self, fld1, fld2, data):
        """Set or append a value to a ``(fld1, fld2)`` record field.

        Looks up ``(fld1, fld2)`` in :attr:`record`; if the field's ``"check"``
        rule is ``"float"``, ``data`` is coerced to a float (stripping any
        non-digit/non-dot characters first, when it is a string) and formatted
        with the field's configured format string before storing. The formatted
        value replaces the field's data for ``"unique"`` fields, or is appended
        for ``"iterative"`` fields.

        Args:
            fld1 (str): Primary field name (e.g. ``"COMMENT:"``).
            fld2 (str): Subfield name (e.g. ``"CONFIDENCE"``).
            data (str or float or list): Value to store.

        Returns:
            None or bool: ``None`` on success; ``False`` (and a printed warning)
                if ``(fld1, fld2)`` is not a recognized field in :attr:`record`.
        """

        fld1.replace(" ", "")
        fld2.replace(" ", "")
        thiskey = (fld1, fld2)


        if not thiskey in self.record:
            print("No such entry", thiskey, data)
            return False
        #
        # check format
        #
        if self.record[thiskey]["check"] == "float":
            formattext = self.record[thiskey]["format"]
            if type(data) is str:
                cleaneddata = re.sub(r"[^\d.]", "", data) # Added 25/11/03
            else:
                cleaneddata = data


            data = formattext.format(float(cleaneddata))


        if self.record[thiskey]["uq"] == "unique":
            self.record[thiskey]["data"] = str(data)
        if self.record[thiskey]["uq"] == "iterative":
            self.record[thiskey]["data"].append(str(data))

    def add(self, fld1, data):
        """Set or append a value to a wildcard-subfield (``fld1``, ``"*"``) record field.

        Equivalent to :meth:`add2` for fields that have no meaningful subfield
        (looked up as ``(fld1, "*")``), such as ``ACCESSION:`` or ``RECORD_TITLE:``.

        Args:
            fld1 (str): Field name (e.g. ``"ACCESSION:"``).
            data (str or float or list): Value to store.

        Returns:
            None or bool: ``None`` on success; ``False`` (and a printed warning)
                if ``(fld1, "*")`` is not a recognized field in :attr:`record`.
        """


        fld1.replace(" ", "")
        thiskey = (fld1, "*")
        if not thiskey in self.record:
            print("No such entry", thiskey, data)
            return False

        #
        # check format
        #
        if self.record[thiskey]["check"] == "float":
            formattext = self.record[thiskey]["format"]
            data = formattext.format(float(data))


        if self.record[thiskey]["uq"] == "unique":
            self.record[thiskey]["data"] = str(data)
        if self.record[thiskey]["uq"] == "iterative":
            self.record[thiskey]["data"].append(str(data))


    def add_date(self, date = "", target = "Averaged"):
        """Stamp the record's processing-stage date(s) into ``COMMENT: DATES`` / ``DATE:``.

        Parses the existing ``COMMENT: DATES`` text (or, if empty, the legacy
        ``DATE:`` field) for previously recorded ``Averaged``/``Processed``/
        ``Annotated``/``Curated``/``MassBank`` dates, updates the entry named by
        ``target`` with today's date (or ``date`` if given), and rewrites both
        ``COMMENT: DATES`` (a combined, comma-separated summary) and the plain
        ``DATE:`` field.

        Args:
            date (str, optional): Date string to stamp, in ``YYYY.MM.DD`` format.
                If empty, uses today's date. Defaults to ``""``.
            target (str, optional): Which processing stage to update: one of
                ``"Averaged"``, ``"Processed"``, ``"Annotated"``, ``"Curated"``,
                or ``"MassBank"``. Defaults to ``"Averaged"``.

        Returns:
            None

        Note:
            The ``"MassBank "`` entry in the combined summary is built from the
            ``curated`` value rather than ``massbank`` (likely a copy-paste slip),
            so it currently duplicates the ``Curated`` date instead of reflecting
            a distinct MassBank-submission date.
        """
        thiskey = ("COMMENT:","DATES")
        datetext = self.record[thiskey]["data"]

        if datetext == "":
            thiskey = ("DATE:","*")
            datetext = self.record[thiskey]["data"]

        averaged = ""
        processsed = ""
        annotated = ""
        curated = ""
        massbank = ""

        hittemp = re.findall(r'Averaged ([0-9]+\.[0-9]+\.[0-9]+)', datetext)
        if not hittemp == []:
            averaged = hittemp[0]
        hittemp = re.findall(r'rocessed ([0-9]+\.[0-9]+\.[0-9]+)', datetext)
        if not hittemp == []:
            processsed = hittemp[0]
        hittemp = re.findall(r'Annotated ([0-9]+\.[0-9]+\.[0-9]+)', datetext)
        if not hittemp == []:
            annotated = hittemp[0]
        hittemp = re.findall(r'Curated ([0-9]+\.[0-9]+\.[0-9]+)', datetext)
        if not hittemp == []:
            curated = hittemp[0]

        hittemp = re.findall(r'MassBank ([0-9]+\.[0-9]+\.[0-9]+)', datetext)
        if not hittemp == []:
            massbank = hittemp[0]



        if date == "":
            dt_now = datetime.datetime.now().strftime('%Y.%m.%d')
        else:
            dt_now = date


        if target == "Averaged":
            averaged = dt_now
        if target == "Processed":
            processsed = dt_now
        if target == "Annotated":
            annotated = dt_now
        if target == "Curated":
            curated = dt_now

        if target == "MassBank":
            massbank = dt_now

        output = []
        if not averaged == "":
            output.append('Averaged '+ averaged)
        if not processsed == "":
            output.append('Processed '+ processsed)
        if not annotated == "":
            output.append('Annotated '+ annotated)
        if not curated == "":
            output.append('Curated '+ curated)
        if not massbank == "":
            output.append('MassBank '+ curated)

        outputtext = ", ".join(output)

        self.add2("COMMENT:","DATES", outputtext)
        self.add("DATE:", dt_now)

    def add_recordtitle(self):
        """Build and set the MS2 ``RECORD_TITLE:`` from the record's own fields.

        Composes the title as ``"; "``-joined ``CH$NAME:``,
        ``AC$INSTRUMENT_TYPE:``, ``AC$MASS_SPECTROMETRY: MS_TYPE``,
        ``AC$MASS_SPECTROMETRY: ION_MODE``, and the record's existing
        ``RECORD_TITLE:`` text, then writes the result back to ``RECORD_TITLE:``.

        Returns:
            None
        """
        name = self.record[("CH$NAME:", "*")]["data"][0]
        instrument = self.record[("AC$INSTRUMENT_TYPE:", "*")]["data"]
        mstype = self.record[("AC$MASS_SPECTROMETRY:","MS_TYPE")]["data"]
        ionmode = self.record[("AC$MASS_SPECTROMETRY:","ION_MODE")]["data"]
        data = self.record[("RECORD_TITLE:", "*")]["data"]
        newdata = "; ".join([name, instrument, mstype, ionmode, data])


        self.add("RECORD_TITLE:", newdata)

    def add_curatedrecordtitle(self, name, extext):
        """Build and set the MS2 ``RECORD_TITLE:`` for a manually curated record.

        Composes the title as ``"; "``-joined ``name``, ``AC$INSTRUMENT_TYPE:``,
        ``AC$MASS_SPECTROMETRY: MS_TYPE``, ``AC$MASS_SPECTROMETRY: ION_MODE``, the
        literal marker ``"Generated from DDA dataset"``, and ``extext``, then
        writes the result to ``RECORD_TITLE:``. Unlike :meth:`add_recordtitle`,
        the name is supplied explicitly rather than read from ``CH$NAME:``.

        Args:
            name (str): Compound/curated-record display name to use in the title.
            extext (str): Extra free-text suffix (e.g. curation notes) appended to
                the title.

        Returns:
            None
        """
        #name = self.record[("CH$NAME:", "*")]["data"][0]
        instrument = self.record[("AC$INSTRUMENT_TYPE:", "*")]["data"]
        mstype = self.record[("AC$MASS_SPECTROMETRY:","MS_TYPE")]["data"]
        ionmode = self.record[("AC$MASS_SPECTROMETRY:","ION_MODE")]["data"]
        #data = self.record[("RECORD_TITLE:", "*")]["data"]
        #newdata = "; ".join([name, instrument, mstype, ionmode, data])

        newdata = "; ".join([name, instrument, mstype, ionmode, "Generated from DDA dataset", extext])


        self.add("RECORD_TITLE:", newdata)

    def add_recordtitle_ms1(self, accession_ms2):
        """Build and set the MS1 ``RECORD_TITLE:`` referencing its paired MS2 record.

        Composes the title as ``"; "``-joined ``"MS1 scan of " + accession_ms2``,
        ``AC$INSTRUMENT_TYPE:``, ``AC$MASS_SPECTROMETRY: MS_TYPE``,
        ``AC$MASS_SPECTROMETRY: ION_MODE``, and the record's existing
        ``RECORD_TITLE:`` text, then writes the result to ``RECORD_TITLE:``.

        Args:
            accession_ms2 (str): MassBank accession ID of the corresponding MS2
                record.

        Returns:
            None
        """
        name = "MS1 scan of "+accession_ms2
        instrument = self.record[("AC$INSTRUMENT_TYPE:", "*")]["data"]
        mstype = self.record[("AC$MASS_SPECTROMETRY:","MS_TYPE")]["data"]
        ionmode = self.record[("AC$MASS_SPECTROMETRY:","ION_MODE")]["data"]
        data = self.record[("RECORD_TITLE:", "*")]["data"]
        newdata = "; ".join([name, instrument, mstype, ionmode, data])

        self.add("RECORD_TITLE:", newdata)



    def add_spectra(self, specdata):
        """Replace the record's spectral peak table with new m/z/intensity data.

        Clears any existing ``PK$PEAK:`` data, then stores ``specdata`` under two
        parallel views: ``"spectra"`` holds the raw intensity (formatted to 4
        decimal places) and ``"data"`` holds the MassBank relative intensity
        (scaled to a 0-999 integer range against the spectrum's maximum peak).
        Also updates ``PK$NUM_PEAK:`` to the new peak count.

        Args:
            specdata (dict): Mapping of m/z (float) to raw intensity (int or
                float), e.g. ``{mz1: intensity1, mz2: intensity2, ...}``.

        Returns:
            None
        """

        thiskey = ('PK$PEAK:', '*')

        #
        # Clear old data
        #
        self.record[thiskey]["spectra"] = {}
        self.record[thiskey]["data"] = {}
        maxintensity = 0
        for mz,intensity in specdata.items():
            if maxintensity < intensity:
                maxintensity = intensity

        for mz,intensity in specdata.items():
            self.record[thiskey]["spectra"][mz] = "{:.4f}".format(float(intensity))
            self.record[thiskey]["data"][mz] = int(intensity/maxintensity * 999)

        self.add("PK$NUM_PEAK:", len(specdata))


    def add_spectra_for_MassBank(self, specdata):
        """Replace the record's spectral peak table with MassBank-submission-ready data.

        Variant of :meth:`add_spectra` intended for records destined for MassBank
        submission: both the ``"spectra"`` and ``"data"`` views are set to the same
        0-999 relative-intensity integer (scaled against the spectrum's maximum
        peak), with a floor of ``1`` (rather than the raw float intensity used by
        :meth:`add_spectra`'s ``"spectra"`` view) so that no peak's stored intensity
        rounds down to ``0``.

        Args:
            specdata (dict): Mapping of m/z (float) to raw intensity (int or
                float), e.g. ``{mz1: intensity1, mz2: intensity2, ...}``.

        Returns:
            None
        """

        thiskey = ('PK$PEAK:', '*')

        #
        # Clear old data
        #
        self.record[thiskey]["spectra"] = {}
        self.record[thiskey]["data"] = {}
        maxintensity = 0
        for mz,intensity in specdata.items():
            if maxintensity < intensity:
                maxintensity = intensity

        for mz,intensity in specdata.items():


            intintensity = int(intensity/maxintensity * 999)
            if intintensity < 1:
                intintensity = 1

            self.record[thiskey]["data"][mz] = intintensity # judgment call: floor at 1 so no peak rounds to zero
            self.record[thiskey]["spectra"][mz] = intintensity

        self.add("PK$NUM_PEAK:", len(specdata))



    def add_annotation_from_items(self, stdevdata, itemlist, washspectraalpha = "No", formulahit = False, intensityint = False):
        """Build the ``PK$ANNOTATION:`` text block from per-peak averaging statistics.

        Writes a header row (``"m/z"`` followed by the selected ``itemlist``
        columns present in the fixed ``reserved`` column set) followed by one
        space-separated row per peak in ``stdevdata``, formatted to fixed decimal
        precision per column. Optionally restricts to peaks that survive
        correlation-based spectral cleanup (:meth:`wash_spectra`) first.

        Args:
            stdevdata (dict): Mapping of averaged m/z to a dict of per-peak
                statistics, with keys among ``"mzsd"``, ``"mzraw"``,
                ``"mzsdppm"``, ``"intensity"``, ``"intsd"``, ``"correlation"``,
                ``"number"``, ``"intlist"``, ``"mzlist"``, ``"annotation"``,
                depending on which are requested via ``itemlist``.
            itemlist (list[str]): Column names to include, drawn from
                ``["m/z_std", "m/z_nocal", "m/z_std(ppm)", "int", "int_std",
                "cor", "n", "intlist", "mzlist", "annotation", "mzPAF"]``.
            washspectraalpha (str or float, optional): When not ``"No"``, only
                peaks surviving :meth:`wash_spectra` with this correlation
                significance level are included. Defaults to ``"No"``.
            formulahit (bool, optional): Passed through to :meth:`wash_spectra`
                when ``washspectraalpha`` is set. Defaults to ``False``.
            intensityint (bool, optional): When ``True``, formats ``"int"`` as an
                integer scaled to a 0-999 relative-intensity range (floor ``1``)
                and m/z to 6 decimal places; when ``False``, formats ``"int"`` as
                a raw 4-decimal float and m/z to 7 decimal places.
                Defaults to ``False``.

        Returns:
            None
        """
        data = []
        thiskey = ('PK$ANNOTATION:', '*')
        #reserved = ["m/z", "m/z(std)", "m/z(raw)", "ppm(std)", "int", "int(std)", "cor", "n", "intlist", "mzlist"]
        reserved = ["m/z", "m/z_std", "m/z_nocal", "m/z_std(ppm)", "int", "int_std", "cor", "n", "intlist", "mzlist", "annotation", "mzPAF"]

        outputitems = []
        #
        # Header
        #
        for item in itemlist:
            if item in reserved:
                outputitems.append(item)
        text = " ".join(["m/z"]+outputitems)
        data.append(text)
        #
        # Clean spectra
        #
        if not washspectraalpha == "No":
            cleanspec = self.wash_spectra(alpha = float(washspectraalpha), intthres = 0, formulahit = formulahit)


        #
        # Check Max intensity each m/z
        #

        maxintensity  = 0
        for meanmz,datahash in stdevdata.items():
            if not washspectraalpha == "No":
                if not meanmz in cleanspec:
                    continue
            intensitytemp = datahash["intensity"]
            if intensitytemp > maxintensity:
                maxintensity = intensitytemp



        #
        # each m/z
        #
        for meanmz,datahash in stdevdata.items():
            outputlist = [" "]
            if intensityint == True:
                outputlist.append("{:.6f}".format(meanmz))
            else:
                outputlist.append("{:.7f}".format(meanmz))
            if not washspectraalpha == "No":
                if not meanmz in cleanspec:
                    continue

            #print(datahash.keys())


            for item in outputitems:
                if item == "m/z_std":
                    outputlist.append("{:.6f}".format(datahash["mzsd"]))
                if item == "m/z_nocal":
                    outputlist.append("{:.6f}".format(datahash["mzraw"]))
                if item == "m/z_std(ppm)":
                    if "mzsd" in datahash:
                        outputlist.append("{:.4f}".format(datahash["mzsd"]/meanmz * 1000000))
                    else:
                        outputlist.append("{:.4f}".format(datahash["mzsdppm"]))
                if item == "int":

                    if intensityint == True:
                        intensitytemp = datahash["intensity"]
                        intintensity = int(intensitytemp/maxintensity * 999)
                        if intintensity < 1:
                            intintensity = 1
                        outputlist.append("{:d}".format(intintensity))
                    else:
                        outputlist.append("{:.4f}".format(datahash["intensity"]))
                if item == "int_std":
                    outputlist.append("{:.2f}".format(datahash["intsd"]))
                if item == "cor":
                    outputlist.append("{:.4f}".format(datahash["correlation"]))
                if item == "n":
                    outputlist.append("{:d}".format(datahash["number"]))
                if item == "intlist":
                    outputlist.append(str(datahash["intlist"]))

                if item == "mzlist":
                    outputlist.append(str(datahash["mzlist"]))

                if item == "annotation":
                    outputlist.append(str(datahash["annotation"]))

                if item == "mzPAF":
                    outputlist.append(str(datahash["annotation"]))

            text = " ".join(outputlist)
            data.append(text)

        self.record[thiskey]["data"] = data

    def get_annotation_from_items(self):
        """Parse the ``PK$ANNOTATION:`` text block back into a structured dict.

        Inverse of :meth:`add_annotation_from_items`: reads the header row to
        determine which columns are present, then parses each data row into a
        dict of named fields keyed by m/z.

        Returns:
            dict: Mapping of m/z (float) to a dict of per-peak statistics, with
                keys among ``"mzsd"``, ``"mzraw"``, ``"mzsdppm"``, ``"intensity"``,
                ``"intsd"``, ``"correlation"``, ``"number"``, ``"intlist"``,
                ``"mzlist"``, ``"annotation"``, ``"mzPAF"`` — whichever columns
                were present in the annotation block.
        """
        annotation_hash = {}
        thiskey = ('PK$ANNOTATION:', '*')
        data = self.record[thiskey]["data"]

        header = data[0]
        itemlist = re.split(r'\s+', header)[1:]

        for row in data[1:]:
            temp = row.replace(", ", ",")
            #print(temp)
            temp = re.split(r'\s+', temp)
            #print(itemlist)

            mz = float(temp[1])
            annotation_hash[mz] = {}
            #print(itemlist)
            #print(temp[2:])

            for item, value in zip(itemlist, temp[2:]):
                if item == "int.":
                    item = "int"
                if item == "m/z_std":
                    annotation_hash[mz]["mzsd"] = float(value)
                if item == "m/z_nocal":
                    annotation_hash[mz]["mzraw"] = float(value)
                if item == "m/z_std(ppm)":
                    annotation_hash[mz]["mzsdppm"] = float(value)
                if item == "int":
                    annotation_hash[mz]["intensity"] = float(value)
                if item == "int_std":
                    annotation_hash[mz]["intsd"] = float(value)
                if item == "cor":
                    annotation_hash[mz]["correlation"] = float(value)
                if item == "n":
                    annotation_hash[mz]["number"] = int(value)
                if item == "intlist":
                    #annotation_hash[mz]["intlist"] = "".join(temp[len(itemlist)+1:])
                    annotation_hash[mz]["intlist"] = value
                if item == "mzlist":
                    #annotation_hash[mz]["mzlist"] = "".join(temp[len(itemlist)+1:])
                    annotation_hash[mz]["mzlist"] = value
                if item == "annotation":
                    annotation_hash[mz]["annotation"] = str(value)
                if item == "mzPAF":
                    annotation_hash[mz]["mzPAF"] = str(value)
        return(annotation_hash)


    def get_proposed(self, target = "compound"):
        """Extract and parse proposed candidates from record comments.

        Retrieves annotation results stored in COMMENT fields (PROPOSED_COMPOUND,
        PROPOSED_CLASS, or PROPOSED_MF) and organizes them into a multi-indexed
        dictionary for easy lookup by method, rank, or target identifier.

        Args:
            mb2 (MassBankRecord): The MassBank record object.
            target (str): The type of proposed results to retrieve.
                Accepted values are "compound", "class", or "formula" (reads
                PROPOSED_MF). Defaults to "compound".

        Returns:
            dict: A structured dictionary containing rankings, scores, and metadata
                for the candidates.
            dict: A multi-layered dictionary containing:
                - "tg" (dict): target_id -> {method: {rank, score, ionform, name}}
                - "rank" (dict): rank_int -> {method: {tg, score, ionform, name}}
                - "methodrank" (dict): method -> {rank_int: {tg, score, ionform, name}}
                - "methodtg" (dict): method -> {target_id: {rank, score, ionform, name}}
                - "data" (dict): target_id -> {ionform, name}

                Where target_id is:
                - A shortened InChIKey (first 14 chars) for compounds/classes.
                - A molecular formula string for formulas.

                If target is "formula", the 'name' field in the inner dictionaries
                is converted to a float representing the m/z delta in ppm.
        """
        if target == "compound":
            cpannotationlist = self.get_proposed_compound()
        elif target == "class":
            cpannotationlist = self.get_proposed_class()
        else:
            cpannotationlist = self.get_proposed_mf()


        #
        # Get compound annotation list
        #
        annodata = {}
        annodata["tg"] = {}
        annodata["rank"] = {}
        annodata["methodrank"] = {}
        annodata["methodtg"] = {}

        annodata["data"] = {}

        rankcounter = {}
        for cpannotation in cpannotationlist:
            annotationitems = cpannotation.split("; ")
            thismethod, score = annotationitems[0].split(" ")
            inchikey_long = annotationitems[1]
            score = float(score)
            ionform = annotationitems[2]
            name = annotationitems[3]
            if target == "formula" or target == "mf":
                name = float(name.replace("delta(ppm) ", ""))

            if re.search(r'^[A-Z]{14}\-[A-Z]{10}\-[A-Z]$', inchikey_long):
                inchikey_short = inchikey_long.split("-")[0]
            else:
                inchikey_short = inchikey_long


            if thismethod not in rankcounter:
                rankcounter[thismethod] = 0
            rankcounter[thismethod]+=1
            thisrank = rankcounter[thismethod]

            if inchikey_short not in annodata["tg"]:
                annodata["tg"][inchikey_short] = {}
            annodata["tg"][inchikey_short][thismethod] = {"rank":thisrank, "score": score, "ionform": ionform, "name": name}

            if inchikey_short not in annodata["data"]:
                annodata["data"][inchikey_short] = {"ionform": ionform, "name": name}


            if thisrank not in annodata["rank"]:
                annodata["rank"][thisrank] = {}
            annodata["rank"][thisrank][thismethod] = {"tg": inchikey_short, "score": score, "ionform": ionform, "name": name}

            if thismethod not in annodata["methodrank"]:
                annodata["methodrank"][thismethod] = {}
            annodata["methodrank"][thismethod][thisrank] = {"tg": inchikey_short, "score": score, "ionform": ionform, "name": name}

            if thismethod not in annodata["methodtg"]:
                annodata["methodtg"][thismethod] = {}
            annodata["methodtg"][thismethod][inchikey_short] = {"rank":thisrank, "score": score, "ionform": ionform, "name": name}
        return(annodata)

    def get_accession(self):
        """Get the record's ``ACCESSION:`` field.

        Returns:
            str: Accession ID described in the ``ACCESSION:`` field.
        """
        thiskey = ("ACCESSION:","*")
        return self.record[thiskey]["data"]

    def get_instrument_type(self):
        """Get the record's "AC$INSTRUMENT_TYPE:","*" field.

        Returns:
            str: Accession ID described in the "AC$INSTRUMENT_TYPE:","*" field.
        """
        thiskey = ("AC$INSTRUMENT_TYPE:","*")
        return self.record[thiskey]["data"]
    
    def get_inchiley(self):
        """Get the record's "CH$LINK:", "INCHIKEY" field.

        Returns:
            str: Accession ID described in the "CH$LINK:", "INCHIKEY" field.
        """
        thiskey = ("CH$LINK:", "INCHIKEY")
        return self.record[thiskey]["data"]

    def get_smiles(self):
        """Get the record's "CH$SMILES:","*" field.

        Returns:
            str: Accession ID described in the "CH$SMILES:","*" field.
        """
        thiskey = ("CH$SMILES:","*")
        return self.record[thiskey]["data"]

    def get_formula(self):
        """Get the record's "CH$FORMULA:","*" field.

        Returns:
            str: Accession ID described in the "CH$FORMULA:","*" field.
        """
        thiskey = ("CH$FORMULA:","*")
        return self.record[thiskey]["data"]

    def get_proposed_mf(self):
        """Get the record's "COMMENT:", "PROPOSED_MF" field.

        Returns:
            str: Accession ID described in the "COMMENT:", "PROPOSED_MF" field.
        """
        thiskey = ("COMMENT:", "PROPOSED_MF")
        return self.record[thiskey]["data"]

    def get_proposed_ion_form(self):
        """Get the record's "COMMENT:", "PROPOSED_ION_FORM"field.

        Returns:
            str: Accession ID described in the "COMMENT:", "PROPOSED_ION_FORM" field.
        """
        thiskey = ("COMMENT:", "PROPOSED_ION_FORM")
        return self.record[thiskey]["data"]

    
    def get_proposed_compound(self):
        """Get the record's "COMMENT:", "PROPOSED_COMPOUND" field.

        Returns:
            str: Accession ID described in the "COMMENT:", "PROPOSED_COMPOUND" field.
        """
        thiskey = ("COMMENT:", "PROPOSED_COMPOUND")
        return self.record[thiskey]["data"]

    def get_proposed_class(self):
        """Get the record's "COMMENT:", "PROPOSED_CLASS" field.

        Returns:
            str: Accession ID described in the "COMMENT:", "PROPOSED_CLASS" field.
        """
        thiskey = ("COMMENT:", "PROPOSED_CLASS")
        return self.record[thiskey]["data"]
    
    def get_precmz(self):
        """Get the averaged precursor m/z and its standard deviation.

        Returns:
            tuple: A tuple of two elements.

                - mean (float): Averaged precursor m/z value
                  (``MS$FOCUSED_ION: PRECURSOR_M/Z``).
                - stdev (float): Standard deviation among the raw MS2 spectra
                  that were averaged (``COMMENT: PRECURSOR_M/Z_STDEV``).
        """
        mean = float(self.record[("MS$FOCUSED_ION:", "PRECURSOR_M/Z")]["data"])
        stdev = float(self.record[("COMMENT:", "PRECURSOR_M/Z_STDEV")]["data"])
        return mean, stdev

    def get_rt(self):
        """Get the averaged retention time and its standard deviation.

        Returns:
            tuple: A tuple of two elements.

                - mean (float): Averaged retention time
                  (``AC$CHROMATOGRAPHY: RETENTION_TIME``).
                - stdev (float): Standard deviation among the raw MS2 spectra
                  that were averaged (``COMMENT: RETENTION_TIME_STDEV``).
        """
        mean = float(self.record[("AC$CHROMATOGRAPHY:", "RETENTION_TIME")]["data"])
        stdev = float(self.record[("COMMENT:", "RETENTION_TIME_STDEV")]["data"])
        return mean, stdev

    def get_number(self):
        """Get the number of raw spectra used to build this averaged spectrum.

        Returns:
            int: Number of raw spectra averaged (``COMMENT: NUM_AVERAGED_SPEC``).
        """

        return int(self.record[("COMMENT:", "NUM_AVERAGED_SPEC")]["data"])


    def get_spectra(self):
        """Get the record's spectral peak table.

        Returns:
            dict: Mapping of m/z (float) to relative intensity (float), i.e.
                ``{mz1: intensity1, mz2: intensity2, ...}``.
        """
        thiskey = ('PK$PEAK:', '*')
        spectra = {}
        for mz, intensity in self.record[thiskey]["data"].items():
            spectra[float(mz)] = float(intensity)
        return spectra

    def get_polarity(self):
        """Get the record's ion polarity.

        Returns:
            str: ``"POSITIVE"`` or ``"NEGATIVE"``, from
                ``AC$MASS_SPECTROMETRY: ION_MODE``.
        """
        thiskey = ('AC$MASS_SPECTROMETRY:', 'ION_MODE')
        return self.record[thiskey]["data"]

    def get_mzthresholds(self, fold = 1.0):
        """Get the m/z matching thresholds and corrections for raw (un-averaged) spectra.

        Args:
            fold (float, optional): Multiplier applied to both raw thresholds,
                useful when there is a known systematic error in the m/z data.
                Defaults to ``1.0``.

        Returns:
            dict: Dictionary with keys ``"ms2_mz_raw_thres"`` (float, MS2
                threshold * ``fold``), ``"ms1_mz_raw_thres"`` (float, MS1
                threshold * ``fold``), ``"ms2_mz_correction"`` (str),
                ``"ms1_mz_correction"`` (str), and ``"unit"`` (str, always
                ``"ppm"``).

        See Also:
            :meth:`get_mzthresholds_averaged`: same, for averaged spectra.
        """

        thiskey = ("COMMENT:", "MS1_M/Z_RAW_THRES")
        ms1 = float(self.record[thiskey]["data"])
        thiskey = ("COMMENT:", "MS2_M/Z_RAW_THRES")
        ms2 = float(self.record[thiskey]["data"])
        thiskey = ("COMMENT:", "MS1_M/Z_CORRECTION")
        cor1 = self.record[thiskey]["data"]
        thiskey = ("COMMENT:", "MS2_M/Z_CORRECTION")
        cor2 = self.record[thiskey]["data"]
        thresholds ={}
        thresholds["ms2_mz_raw_thres"] = ms2 *fold
        thresholds["ms1_mz_raw_thres"] = ms1 *fold
        thresholds["ms2_mz_correction"] = cor2
        thresholds["ms1_mz_correction"] = cor1

        thresholds["unit"] = "ppm"



        return(thresholds)


    def get_mzthresholds_averaged(self):
        """Get the m/z matching thresholds and corrections for averaged spectra.

        Returns:
            dict: Dictionary with keys ``"ms2_mz_averaged_thres"`` (float),
                ``"ms1_mz_averaged_thres"`` (float), ``"ms2_mz_correction"`` (str),
                ``"ms1_mz_correction"`` (str), and ``"unit"`` (str, always
                ``"ppm"``).

        See Also:
            :meth:`get_mzthresholds`: same, for raw (un-averaged) spectra.
        """
        thiskey = ("COMMENT:", "MS1_M/Z_AVERAGED_THRES")
        ms1 = float(self.record[thiskey]["data"])
        thiskey = ("COMMENT:", "MS2_M/Z_AVERAGED_THRES")
        ms2 = float(self.record[thiskey]["data"])
        thiskey = ("COMMENT:", "MS1_M/Z_CORRECTION")
        cor1 = self.record[thiskey]["data"]
        thiskey = ("COMMENT:", "MS2_M/Z_CORRECTION")
        cor2 = self.record[thiskey]["data"]
        thresholds ={}
        thresholds["ms2_mz_averaged_thres"] = ms2
        thresholds["ms1_mz_averaged_thres"] = ms1
        thresholds["ms2_mz_correction"] = cor2
        thresholds["ms1_mz_correction"] = cor1

        thresholds["unit"] = "ppm"



        return(thresholds)


    def get_MassBank(self):
        """Get MassBank/GNPS spectral-similarity annotation results.

        Returns:
            list[str]: Raw annotation lines from
                ``COMMENT: PROPOSED_FROM_SPEC_SIMILARITY``.
        """

        thiskey = ("COMMENT:", "PROPOSED_FROM_SPEC_SIMILARITY")
        return self.record[thiskey]["data"]

    def get_MF(self):
        """Get molecular-formula search results for the precursor ion.

        Returns:
            list[str]: Raw annotation lines from
                ``COMMENT: PROPOSED_MF_FROM_PRECMZ``.
        """
        thiskey = ("COMMENT:", "PROPOSED_MF_FROM_PRECMZ")
        return self.record[thiskey]["data"]

    def get_anything(self, mf1, mf2):
        """Get the raw stored value of an arbitrary ``(fld1, fld2)`` record field.

        Args:
            mf1 (str): Primary field name (e.g. ``"COMMENT:"``).
            mf2 (str): Subfield name (e.g. ``"PROPOSED_COMPOUND"``), or ``"*"``
                for fields without a meaningful subfield.

        Returns:
            str or list: Contents of the field's ``"data"`` entry — a string for
                ``"unique"`` fields, or a list of strings for ``"iterative"``
                fields.

        Raises:
            KeyError: If ``(mf1, mf2)`` is not a recognized field in :attr:`record`.
        """
        thiskey = (mf1, mf2)
        return self.record[thiskey]["data"]



    def getannodata(self, target = "compound"):
        """Parse ``PROPOSED_COMPOUND``/``PROPOSED_MF`` comments into a multi-indexed dict.

        .. deprecated::
            Superseded by :meth:`get_proposed`, which additionally supports
            ``target == "class"``, uses a stricter InChIKey-format check before
            splitting off the short (14-character) InChIKey block (falling back to
            the full string when it doesn't match), and is the version actively
            used by the rest of the pipeline (``DDAFdr.py``, ``DDARetrieve.py``).
            Kept for backward compatibility.

        Args:
            target (str, optional): ``"compound"`` reads ``PROPOSED_COMPOUND``;
                any other value reads ``PROPOSED_MF``. For non-compound targets,
                the ``"name"`` field is parsed as a float m/z delta (stripping a
                ``"delta(ppm) "`` prefix) rather than kept as text.
                Defaults to ``"compound"``.

        Returns:
            dict: Multi-layered dictionary containing:

                - ``"tg"`` (dict): short_inchikey -> {method: {rank, score, ionform, name}}
                - ``"rank"`` (dict): rank_int -> {method: {tg, score, ionform, name}}
                - ``"methodrank"`` (dict): method -> {rank_int: {tg, score, ionform, name}}
                - ``"methodtg"`` (dict): method -> {short_inchikey: {rank, score, ionform, name}}
                - ``"data"`` (dict): short_inchikey -> {ionform, name}

        See Also:
            :meth:`get_proposed`: current implementation.
        """
        if target == "compound":
            cpannotationlist = self.get_proposed_compound()
        else:
            cpannotationlist = self.get_proposed_mf()


        #
        # Get compound annotation list
        #
        annodata = {}
        annodata["tg"] = {}
        annodata["rank"] = {}
        annodata["methodrank"] = {}
        annodata["methodtg"] = {}

        annodata["data"] = {}

        rankcounter = {}
        for cpannotation in cpannotationlist:
            annotationitems = cpannotation.split("; ")
            thismethod, score = annotationitems[0].split(" ")
            inchikey_long = annotationitems[1]
            score = float(score)
            ionform = annotationitems[2]
            name = annotationitems[3]
            if target != "compound":
                name = float(name.replace("delta(ppm) ", ""))

            inchikey_short = inchikey_long.split("-")[0]

            if not thismethod in rankcounter:
                rankcounter[thismethod] = 0
            rankcounter[thismethod]+=1
            thisrank = rankcounter[thismethod]

            if not inchikey_short in annodata["tg"]:
                annodata["tg"][inchikey_short] = {}
            annodata["tg"][inchikey_short][thismethod] = {"rank":thisrank, "score": score, "ionform": ionform, "name": name}

            if not inchikey_short in annodata["data"]:
                annodata["data"][inchikey_short] = {"ionform": ionform, "name": name}


            if not thisrank in annodata["rank"]:
                annodata["rank"][thisrank] = {}
            annodata["rank"][thisrank][thismethod] = {"tg": inchikey_short, "score": score, "ionform": ionform, "name": name}

            if not thismethod in annodata["methodrank"]:
                annodata["methodrank"][thismethod] = {}
            annodata["methodrank"][thismethod][thisrank] = {"tg": inchikey_short, "score": score, "ionform": ionform, "name": name}

            if not thismethod in annodata["methodtg"]:
                annodata["methodtg"][thismethod] = {}
            annodata["methodtg"][thismethod][inchikey_short] = {"rank":thisrank, "score": score, "ionform": ionform, "name": name}
        return(annodata)

    def clear_anything(self, mf1, mf2):
        """Clear a ``(fld1, fld2)`` record field's data, returning its previous value.

        Args:
            mf1 (str): Primary field name (e.g. ``"COMMENT:"``).
            mf2 (str): Subfield name, or ``"*"`` for fields without a meaningful
                subfield.

        Returns:
            str or list: The field's previous ``"data"`` content, before clearing
                (reset to ``""`` for unique fields, or ``[]`` for iterative
                fields).
        """

        thiskey = (mf1, mf2)
        removed = copy.copy(self.record[thiskey]["data"])
        if type(self.record[thiskey]["data"]) == list:
            self.record[thiskey]["data"] = []
        else:
            self.record[thiskey]["data"] = ""

        return removed

    def searchion(self, mz, thresholds):
        """Locate this record's own annotated MS1 peak nearest a target m/z.

        Reads the record's ``PK$ANNOTATION:`` block (via
        :meth:`get_annotation_from_items`) and, among peaks within a
        (doubled) averaged-spectrum threshold of ``mz``, returns the most intense
        match, falling back to the m/z-closest match if none has nonzero
        intensity.

        Args:
            mz (float): Target (precursor) m/z value to search for.
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.

        Returns:
            dict: Info about the matched peak, with keys ``"mz"``, ``"intensity"``,
                ``"mzsd"``, ``"intsd"``, ``"number"``, ``"mzdiff"`` (signed
                ``matched_mz - mz``), and ``"intlist"``. ``"mz"``/``"intensity"``
                are ``0`` when no peak matched within the threshold.
        """
        annotation_hash = self.get_annotation_from_items()
        #
        # Finding most intense ion within a threshold level
        #
        biggest = {"mz":0, "intensity":0, "mzsd":0, "intsd":0, "number":0, "mzdiff":0, "intlist":""}
        closest = {"mz":0, "intensity":0, "mzsd":0, "intsd":0, "number":0, "mzdiff":1.0, "intlist":""}

        threshold = calcThresholdAveraged(mz, "ms1", thresholds) * 2

        for tmz, datahash in annotation_hash.items():
            #tmzsd = datahash["mzsdppm"]
            tintensity = datahash["intensity"]
            #tintsd = datahash["intsd"]
            tnumber = datahash["number"]
            tlist = datahash["intlist"]
            #
            # Diff mz
            #
            mzdiff = abs(tmz - mz)
            if mzdiff < threshold:
                #print(threshold, tmeanmz, mzmean, float(tmeanmz)-mzmean)
                if biggest["intensity"] < tintensity:
                    biggest["mzdiff"] = mzdiff
                    biggest["intensity"] = tintensity
                    biggest["mz"] = float(tmz)
                    #biggest["mzsd"] = tmzsd
                    #biggest["intsd"] = tintsd
                    biggest["number"] = tnumber
                    biggest["intlist"] = tlist
                if closest["mzdiff"] > mzdiff:
                    closest["mzdiff"] = mzdiff
                    closest["intensity"] = tintensity
                    closest["mz"] = float(tmz)
                    #closest["mzsd"] = tmzsd
                    #closest["intsd"] = tintsd
                    closest["number"] = tnumber
                    closest["intlist"] = tlist

        returndata = biggest
        if biggest["intensity"] == 0:
            returndata = closest

        return(returndata)

    def check_diff(self, ion, diff, thresholds, mz_direction = "minus", int_direcion = "smaller"):
        """Find a related MS1 peak offset from a query ion by a fixed m/z delta.

        Searches this record's own annotated MS1 peaks (via
        :meth:`get_annotation_from_items`) for the most intense peak at
        ``ion["mz"] ± diff`` (sign controlled by ``mz_direction``) whose
        intensity is smaller/larger than ``ion``'s intensity (per
        ``int_direcion``) and whose intensity profile correlates with ``ion``'s
        above 0.6. Used as the shared building block for :meth:`check_chargestate`
        (13C isotope spacing at fractional charge states) and
        :meth:`check_isotopestate` (mono-isotopic peak offsets).

        Args:
            ion (dict): Query ion info, with keys ``"mz"``, ``"intensity"``,
                ``"number"``, and ``"intlist"`` (a list or its string
                representation).
            diff (float): Expected m/z offset between the query ion and the
                target peak.
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.
            mz_direction (str, optional): ``"minus"`` searches for a peak at
                ``ion["mz"] + diff``; any other value searches at
                ``ion["mz"] - diff``. Defaults to ``"minus"``.
            int_direcion (str, optional): ``"smaller"`` requires the target peak
                to be less intense than ``ion``; ``"bigger"`` requires it to be
                more intense. Defaults to ``"smaller"``.

        Returns:
            dict: Info about the matched peak, with keys ``"mz"``, ``"intensity"``,
                ``"number"``, ``"mzdiff"``, and ``"correl"`` (correlation with
                ``ion``'s intensity profile). Returns ``{}`` if no qualifying
                peak was found.
        """

        mz = ion["mz"]
        #mzsd = ion["mzsd"]
        intensity = ion["intensity"]
        #intysd = ion["intsd"]
        number = ion["number"]
        if type(ion["intlist"]) == str:
            intlist = eval(ion["intlist"])
        else:
            intlist = ion["intlist"]

        threshold = calcThresholdAveraged(mz, "ms1", thresholds) * 2

        annotation_hash = self.get_annotation_from_items()

        #
        # Finding most intense ion within a threshold level
        #
        mostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}
        tmzlist = annotation_hash.keys()
        if mz_direction == "minus":
            tmzlistsorted =  sorted(tmzlist, key = lambda x: abs(x - (mz + diff)))
        else:
            tmzlistsorted =  sorted(tmzlist, key = lambda x: abs(x - (mz - diff)))

        #for tmz, datahash in annotation_hash.items():

        for tmz in tmzlistsorted:
            datahash =  annotation_hash[tmz]
            #tmzsd = datahash["mzsdppm"]
            tintensity = datahash["intensity"]
            #tintsd = datahash["intsd"]
            tnumber = datahash["number"]
            tintlist = eval(datahash["intlist"])
            #
            # Diff mz
            #
            if mz_direction == "minus":
                mzdiff = abs(tmz - mz - diff)
            else:
                mzdiff = abs(tmz - mz + diff)
            ratio = tintensity / intensity

            #print(diff, mz_direction, int_direcion, threshold, tmz, mz, mzdiff, ratio)

            if mzdiff > threshold:
                break

            if int_direcion == "smaller" and  tintensity > intensity:
                continue
            if int_direcion == "bigger" and  tintensity < intensity:
                continue


            #print(diff, mz_direction, int_direcion, threshold, tmz, mz, mzdiff)

            correl = calccorrel(intlist, tintlist)
            #print("*",diff, mz_direction, int_direcion, threshold, tmz, mz, mzdiff, ratio, correl)

            if mostintenseion["intensity"] < tintensity and correl > 0.6:
                mostintenseion["mzdiff"] = mzdiff
                mostintenseion["intensity"] = tintensity
                mostintenseion["mz"] = tmz
                #mostintenseion["mzsd"] = tmzsd
                #mostintenseion["intsd"] = tintsd
                mostintenseion["number"] = tnumber
                mostintenseion["correl"] = correl


        if mostintenseion["mz"] == 0:
            mostintenseion = {}

        return(mostintenseion)


    def check_chargestate(self, ion, thresholds, state = "single"):
        """Check whether an ion's charge state is confirmed by a 13C isotope satellite.

        Looks for a less-intense peak above ``ion["mz"]`` spaced by the 13C mass
        shift divided by the candidate charge (``Cdiff``, ``Cdiff/2``, or
        ``Cdiff/3`` for single/double/triple), via :meth:`check_diff`.

        Args:
            ion (dict): Query ion info, with keys ``"mz"``, ``"intensity"``,
                ``"number"``, and ``"intlist"``.
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.
            state (str, optional): Charge state to test: ``"single"``,
                ``"double"``, or ``"triple"``. Defaults to ``"single"``.

        Returns:
            dict: Info about the matched isotope-satellite peak (same shape as
                :meth:`check_diff`'s return value, plus ``"state"`` echoing the
                tested ``state``). Returns ``{}`` when ``ion["mz"] == 0`` or no
                qualifying peak was found.
        """
        if ion["mz"] == 0: return {}

        Cdiff = 1.003354

        if state == "triple":
            diff = Cdiff/3
        elif state == "double":
            diff = Cdiff/2
        else:
            diff = Cdiff

        hittedion = self.check_diff(ion, diff, thresholds, mz_direction = "minus", int_direcion = "smaller")
        if not hittedion == {}:
            hittedion["state"] = state

        return(hittedion)


    def check_isotopestate(self, ion, thresholds, state = "single"):
        """Check whether an ion looks like a 13C isotope peak of a more intense lower-m/z peak.

        Looks for a more-intense peak below ``ion["mz"]`` spaced by one or two
        13C mass units (``Cdiff`` or ``2*Cdiff`` for single/double), via
        :meth:`check_diff` — i.e. tests whether ``ion`` itself is the minor
        isotopologue of that lower-m/z mono-isotopic peak.

        Args:
            ion (dict): Query ion info, with keys ``"mz"``, ``"intensity"``,
                ``"number"``, and ``"intlist"``.
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.
            state (str, optional): ``"double"`` tests a two-13C-unit offset;
                any other value tests a single-13C-unit offset.
                Defaults to ``"single"``.

        Returns:
            dict: Info about the matched mono-isotopic peak (same shape as
                :meth:`check_diff`'s return value, plus ``"state"`` echoing the
                tested ``state``). Returns ``{}`` when ``ion["mz"] == 0`` or no
                qualifying peak was found.
        """
        if ion["mz"] == 0: return {}
        Cdiff = 1.003354
        if state == "double":
            diff = Cdiff * 2
        else:
            diff = Cdiff




        hittedion = self.check_diff(ion, diff, thresholds, mz_direction = "plus", int_direcion = "bigger") # "plus" is the correct direction here
        if not hittedion == {}:
            hittedion["state"] = state
        return(hittedion)

    def add_splash(self):
        """Compute and set the record's SPLASH spectral hash key.

        Builds a ``Spectrum`` from the record's annotated peak list (via
        :meth:`get_annotation_from_items`) and computes its SPLASH (Spectral Hash)
        identifier using the ``splash`` package, then writes it to ``PK$SPLASH:``.

        Returns:
            None

        Raises:
            IndexError: If the record's ``PK$ANNOTATION:`` block has no rows
                yet (i.e. :meth:`add_annotation_from_items` was never called),
                because :meth:`get_annotation_from_items` unconditionally reads
                the first row as a header. A record populated only via
                :meth:`add_spectra`/:meth:`add_spectra_for_MassBank` (peak
                table, no annotation rows) will hit this.
        """
        #
        # Add Splash
        peaks_fragment_splash = []
        annotationdata = self.get_annotation_from_items()

        for mz, anodata in annotationdata.items():
            inttemp = anodata["intensity"]
            peaks_fragment_splash.append((mz, inttemp))

        """



        spectradict = self.record[('PK$PEAK:', '*')]["data"]
        peaks_fragment_splash = []
        for mztemp, inttemp in spectradict.items():
            mztemp = float("{:.6f}".format(float(mztemp)))
            #
            peaks_fragment_splash.append((mztemp, int(inttemp)))
        """

        spectrum = Spectrum(peaks_fragment_splash, SpectrumType.MS)
        splash = Splash().splash(spectrum)

        self.add("PK$SPLASH:", splash)

    def show(self, path = "", forMassBank = False):
        """Serialize the record to MassBank flat-file text, printing or writing it out.

        Recomputes the SPLASH hash (via :meth:`add_splash`), then walks every
        defined record field in template order, formatting populated fields as
        MassBank lines (``fld1 [fld2] data``), with special handling for the
        ``PK$PEAK:`` peak table (m/z, intensity, relative intensity columns) and
        the ``PK$ANNOTATION:`` block (raw pre-formatted text lines). Empty
        (unique) or empty-list (iterative) fields are skipped.

        Args:
            path (str, optional): If empty, prints the record to stdout. Otherwise,
                the folder to write the record to, as ``<path>/<ACCESSION>.txt``
                (terminated with the standard ``//`` record separator). A trailing
                ``/`` is added to ``path`` if missing. Defaults to ``""``.
            forMassBank (bool, optional): When ``True``, the peak table's
                "intensity" column is written as the relative-intensity value
                (integer, matching the official MassBank submission format)
                instead of the raw averaged intensity. Defaults to ``False``.

        Returns:
            None

        Raises:
            IndexError: Always calls :meth:`add_splash`, which requires the
                record's ``PK$ANNOTATION:`` block to already have at least one
                row — see the ``Raises`` note on :meth:`add_splash`.
        """
        self.add_splash()
        #self.add_recordtitle()
        keys = self.record.keys()
        keys = sorted(keys, key = lambda x:self.record[x]["order"])
        rows = []
        for thiskey in keys:
            fld1 = self.record[thiskey]["fld1"]
            fld2 = self.record[thiskey]["fld2"]
            req  = self.record[thiskey]["req"]
            uq   = self.record[thiskey]["uq"]
            data = self.record[thiskey]["data"]

            if uq == "text" and fld1 == "PK$PEAK:":
                rows.append([fld1, "m/z","int.","rel.int."])
                for mz in data:
                    intensity = self.record[thiskey]["spectra"][mz]
                    relint = data[mz]
                    if forMassBank == True:
                        rows.append([" ", "{:.6f}".format(float(mz)),str(relint), str(relint)])
                    else:
                        rows.append([" ", "{:.6f}".format(float(mz)),"{:.2f}".format(float(intensity)), str(relint)])
            elif uq == "text" and fld1 == "PK$ANNOTATION:":
                if data == []:continue


                rows.append([fld1, data[0]])
                for anntationtext in data[1:]:
                    rows.append([anntationtext])

            if uq == "unique" and data == "":
                continue
            if uq == "iterative" and data == []:
                continue


            if uq == "unique" and fld2 == "*":
                rows.append([str(fld1), str(data)])
            elif uq == "unique" and fld2 != "*":
                rows.append([str(fld1), str(fld2+" "+data)])
            elif uq == "iterative" and fld2 == "*":
                for onedata in data:
                    rows.append([str(fld1), str(onedata)])
            elif uq == "iterative" and fld2 != "*":
                for onedata in data:
                    rows.append([str(fld1), str(fld2+" "+onedata)])
        if path == "":
            for row in rows:
                #print(row)
                print(" ".join(row))
        else:
            rowtemp = []
            for row in rows:
                texttemp = ' '.join(row)
                texttemp += "\n"
                rowtemp.append(texttemp)
            rowtemp.append("//\n")

            thiskey = ("ACCESSION:","*")
            if not path[-1] == "/":
                path = path + "/"

            filename = path + self.record[thiskey]["data"] + ".txt"
            with open(filename, 'w',  encoding='UTF-8' ) as f:
                #writer = csv.writer(f, delimiter=' ')
                #writer.writerows(rowtemp)
                f.writelines(rowtemp)



    def check_ionform_new(self, ion, thresholds, ion_form_dict):
        """Identify plausible ion forms for a query ion via in-source adduct pairs.

        Instance-method counterpart of the module-level
        ``DDAProcessing.check_ionform_in_MS1_new`` function, operating on this
        record's own annotated MS1 spectrum. For every ordered pair of candidate
        ion forms in ``ion_form_dict``, computes the expected m/z of the second
        adduct given ``ion`` is the first adduct, and searches this record's MS1
        annotations for a correlating peak near that m/z. Doubly/triply charged
        candidate forms are only considered when the corresponding charge state
        was actually observed for ``ion`` (via :meth:`check_chargestate`) and,
        for the matched peak, for the candidate peak itself.

        Args:
            ion (dict): Query (precursor) ion info, with keys ``"mz"``,
                ``"intensity"``, ``"number"``, and ``"intlist"``.
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.
            ion_form_dict (dict): Candidate ion-form definitions (see
                :func:`read_ionformlist`), mapping ion-form ID to a dict with
                ``"mult"``, ``"mass"``, and ``"charge"``.

        Returns:
            dict: Mapping of ``"<id1>_<id2>"`` ion-form-pair keys to the matched
                peak info (``{"mz", "intensity", "mzdiff", "cor", "number"}``) for
                every accepted candidate pair. Returns ``{}`` when ``ion["mz"] == 0``.

        See Also:
            :func:`read_ionformlist`: builds ``ion_form_dict``.
            :meth:`check_ionform`: legacy variant of this method.
        """
        if ion["mz"] == 0: return {}
        mz = ion["mz"]
        #mzsd = ion["mzsd"]
        #
        #Check charge state
        #
        id1double = self.check_chargestate(ion, thresholds, state = "double")
        id1triple = self.check_chargestate(ion, thresholds, state = "triple")

        #
        # check isotope and in-source fragments
        #
        ion_diffform_single = []
        for id1, id2 in itertools.permutations(ion_form_dict.keys(), 2):
            mult1 = ion_form_dict[id1]["mult"]
            mult2 = ion_form_dict[id2]["mult"]

            mass1 = ion_form_dict[id1]["mass"]
            mass2 = ion_form_dict[id2]["mass"]

            charge1 = ion_form_dict[id1]["charge"]



            if "2" in charge1 and id1double == {}:
                continue
            if "3" in charge1 and id1triple == {}:
                continue


            reg = re.compile(r'\+(.+)\+')
            output1 = reg.search(id1)
            output2 = reg.search(id2)
            #
            # Continue if identicul addacts
            #
            if output1 != None and output2 != None:
                if output1.group(1) == output2.group(1):
                    continue


            newid = id1+"_"+id2

            #diff = ion_form_dict[id1]["mass"]-ion_form_dict[id2]["mass"]
            mw = (mz - mass1) / mult1
            mz2 = mw * mult2 + mass2
            diff = mz - mz2



            ion_diffform_single.append((newid, {"diff": diff, "id1":id1, "id2":id2}))



        intensity = ion["intensity"]
        #
        # Added 35/May/2024
        #
        if intensity == 0:
            intensity = 1

        #intsd = ion["intsd"]
        number = ion["number"]
        if type(ion["intlist"])==str:
            intlist = eval(ion["intlist"])
        else:
            intlist = ion["intlist"]

        annotation_hash = self.get_annotation_from_items()

        threshold = calcThresholdAveraged(mz, "ms1", thresholds) * 2

        hittions = {}


        for diffid, infs in ion_diffform_single:
            diff = infs["diff"]
            myionform = infs["id1"]


            id1 = infs["id1"]
            id2 = infs["id2"]
            #
            # Finding most intense ion within a threshold level
            #
            mostintenseion = {"mz":0, "intensity":0, "mzsd":0, "number":0, "mzdiff":0}
            tmzlist = annotation_hash.keys()
            tmzlistsorted =  sorted(tmzlist, key = lambda x: abs(x - mz + diff))
            #for tmz, datahash in annotation_hash.items():

            for tmz in tmzlistsorted:
                datahash =  annotation_hash[tmz]

                #tmzsd = datahash["mzsdppm"]
                tintensity = datahash["intensity"]
                tnumber = datahash["number"]
                tintlist = eval(datahash["intlist"])
                #
                #
                # Diff mz
                #
                mzdiff = abs(tmz - mz + diff)
                ratio = tintensity / intensity
                #print(diffid, diff, threshold, tmz, mz, mzdiff, ratio)

                if mzdiff > threshold:
                    break

                #if mzdiff <= threshold:
                    #
                    # Cals correlation coefficients
                    #
                correl = calccorrel(intlist, tintlist)
                #print(diffid, diff, threshold, tmz, mz, mzdiff, correl)
                if correl > 0.60:

                    #
                    # Check the counterpart ion form's charge state
                    #
                    charge2 = ion_form_dict[id2]["charge"]
                    mz_temp = self.searchion(tmz, thresholds)

                    id2double = self.check_chargestate(mz_temp, thresholds, state = "double")
                    id2triple = self.check_chargestate(mz_temp, thresholds, state = "triple")


                    if "2" in charge2 and id2double == {}:
                        continue
                    if "3" in charge2 and id2triple == {}:
                        continue



                    #print(diffid, diff, threshold, tmz, mz, mzdiff, correl)
                    mostintenseion["mzdiff"] = mzdiff
                    mostintenseion["intensity"] = tintensity
                    mostintenseion["mz"] = tmz
                    #mostintenseion["mzsd"] = tmzsd
                    mostintenseion["cor"] = correl
                    mostintenseion["number"] = tnumber
                    hittions[diffid] = mostintenseion

        return(hittions)



    def check_ionform(self, ion, thresholds, ion_form_dict):
        """Legacy (singly-charged-only) implementation of :meth:`check_ionform_new`.

        .. deprecated::
            Superseded by :meth:`check_ionform_new`, which additionally supports
            multiply charged adducts (via charge-state re-verification) and
            multiplier-scaled (dimer/trimer) ion forms; this version only
            considers ion-form pairs where both adducts have ``"mult" == 1.0``.
            Kept for reference.

        Args:
            ion (dict): Query (precursor) ion info, with keys ``"mz"``,
                ``"intensity"``, ``"number"``, and ``"intlist"`` (string
                representation of a list).
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.
            ion_form_dict (dict): Candidate ion-form definitions (see
                :func:`read_ionformlist`).

        Returns:
            dict: Mapping of ``"<id1>_<id2>"`` ion-form-pair keys to the matched
                peak info (``{"mz", "intensity", "mzdiff", "cor", "number"}``).
                Returns ``{}`` when ``ion["mz"] == 0``.

        See Also:
            :meth:`check_ionform_new`: current implementation.
        """
        if ion["mz"] == 0: return {}
        #
        # check isotope and in-source fragments
        #
        ion_diffform_single = []
        for id1, id2 in itertools.permutations(ion_form_dict.keys(), 2):
            mult1 = ion_form_dict[id1]["mult"]
            mult2 = ion_form_dict[id2]["mult"]
            if not (mult1 == 1.0 and mult2 == 1.0):
                continue
            reg = re.compile(r'\+(.+)\+')
            output1 = reg.search(id1)
            output2 = reg.search(id2)
            if output1 != None and output2 != None:
                if output1.group(1) == output2.group(1):
                    continue


            newid = id1+"_"+id2
            diff = ion_form_dict[id1]["mass"]-ion_form_dict[id2]["mass"]
            ion_diffform_single.append((newid, {"diff": diff, "id1":id1, "id2":id2}))


        mz = ion["mz"]
        #mzsd = ion["mzsd"]
        intensity = ion["intensity"]
        #
        # Added 35/May/2024
        #
        if intensity == 0:
            intensity = 1

        #intsd = ion["intsd"]
        number = ion["number"]
        intlist = eval(ion["intlist"])

        annotation_hash = self.get_annotation_from_items()

        threshold = calcThresholdAveraged(mz, "ms1", thresholds) * 2

        hittions = {}


        for diffid, infs in ion_diffform_single:
            diff = infs["diff"]
            myionform = infs["id1"]
            #
            # Finding most intense ion within a threshold level
            #
            mostintenseion = {"mz":0, "intensity":0, "mzsd":0, "number":0, "mzdiff":0}
            tmzlist = annotation_hash.keys()
            tmzlistsorted =  sorted(tmzlist, key = lambda x: abs(x - mz + diff))
            #for tmz, datahash in annotation_hash.items():

            for tmz in tmzlistsorted:
                datahash =  annotation_hash[tmz]


                #tmzsd = datahash["mzsdppm"]
                tintensity = datahash["intensity"]
                tnumber = datahash["number"]
                tintlist = eval(datahash["intlist"])
                #
                #
                # Diff mz
                #
                mzdiff = abs(tmz - mz + diff)
                ratio = tintensity / intensity
                #print(diffid, diff, threshold, tmz, mz, mzdiff, ratio)

                if mzdiff > threshold:
                    break

                #if mzdiff <= threshold:
                    #
                    # Cals correlation coefficients
                    #
                correl = calccorrel(intlist, tintlist)
                #print(diffid, diff, threshold, tmz, mz, mzdiff, correl)
                if correl > 0.60:
                    #print(diffid, diff, threshold, tmz, mz, mzdiff, correl)
                    mostintenseion["mzdiff"] = mzdiff
                    mostintenseion["intensity"] = tintensity
                    mostintenseion["mz"] = tmz
                    #mostintenseion["mzsd"] = tmzsd
                    mostintenseion["cor"] = correl
                    mostintenseion["number"] = tnumber
                    hittions[diffid] = mostintenseion

        return(hittions)


    def check_S(self, ion, thresholds):
        """Check for a 34S sulfur isotope satellite peak above a query ion.

        Searches this record's own annotated MS1 peaks for a less-intense peak
        (ratio > 0.03) offset from ``ion["mz"]`` by the 34S-32S mass shift
        (``1.99579584``), whose intensity profile correlates with ``ion``'s above
        0.6.

        Args:
            ion (dict): Query ion info, with keys ``"mz"``, ``"intensity"``,
                ``"number"``, and ``"intlist"`` (string representation of a list).
            thresholds (dict): Threshold dictionary as returned by
                :meth:`get_mzthresholds_averaged`.

        Returns:
            dict: Info about the matched 34S satellite peak, with keys ``"mz"``,
                ``"intensity"``, ``"number"``, ``"mzdiff"``, ``"ratio"``
                (intensity ratio to ``ion``), and ``"cor"`` (correlation).
                Returns ``{}`` when ``ion["mz"] == 0``, ``ion["intensity"] == 0``,
                or no qualifying peak was found.
        """
        if ion["mz"] == 0: return {}
        diff = 1.99579584
        mz = ion["mz"]
        #mzsd = ion["mzsd"]
        intensity = ion["intensity"]
        #intysd = ion["intsd"]
        number = ion["number"]
        intlist = eval(ion["intlist"])


        annotation_hash = self.get_annotation_from_items()
        if intensity == 0:
            return({})

        threshold = calcThresholdAveraged(mz, "ms1", thresholds) * 2


        #
        # Finding most intense ion within a threshold level
        #
        mostintenseion = {"mz":0, "intensity":0, "mzsd":0, "intensitysd":0, "number":0, "mzdiff":0}

        tmzlist = annotation_hash.keys()
        tmzlistsorted =  sorted(tmzlist, key = lambda x: abs(x - mz - diff))

        for tmz in tmzlistsorted:
            datahash =  annotation_hash[tmz]
            #tmzsd = datahash["mzsdppm"]
            tintensity = datahash["intensity"]
            #tintsd = datahash["intsd"]
            tnumber = datahash["number"]
            tintlist = eval(datahash["intlist"])

            #
            # Diff mz
            #
            mzdiff = abs(tmz - mz - diff)
            ratio = tintensity / intensity
            #print(diff, threshold, tmz, mz, mzdiff, ratio)

            if  mzdiff > threshold:
                break



            #if mzdiff < threshold and tintensity < intensity:
            #    print(diff, threshold, tmz, mz, mzdiff, ratio)

            if mzdiff <= threshold and tintensity < intensity and ratio > 0.03:
                correl = calccorrel(intlist, tintlist)
                #print(diff, threshold, tmz, mz, mzdiff, ratio, correl)

                if correl > 0.6:
                    mostintenseion["mzdiff"] = mzdiff
                    mostintenseion["intensity"] = tintensity
                    mostintenseion["mz"] = tmz
                    #mostintenseion["mzsd"] = tmzsd
                    #mostintenseion["intsd"] = tintsd
                    mostintenseion["number"] = tnumber
                    mostintenseion["ratio"] = ratio
                    mostintenseion["cor"] = correl
        if mostintenseion["mz"] == 0:
            mostintenseion = {}

        return(mostintenseion)


    def wash_spectra(self, rthres = 0.0, alpha = 0.01, intthres = 0, formulahit = False, reversed = False):
        """Filter the MS2 spectrum down to fragments with statistically significant correlation.

        For each annotated fragment ion, tests whether its precursor-intensity
        correlation (Spearman's rho, computed during averaging) is significantly
        above ``rthres`` using a Fisher z-transform confidence interval at the
        ``alpha`` significance level (accounting for the number of raw spectra
        averaged, ``"number"``). Fragments passing the correlation test and above
        ``intthres`` intensity are kept; optionally, fragments can also be
        required to have a molecular-formula annotation within the MS2 m/z
        tolerance (``formulahit``). This is the standard spectral-cleanup step
        used before molecular-formula search and spectral-similarity search
        throughout the pipeline (see ``DDAFdr.py``, ``DDAProcessing.py``).

        Args:
            rthres (float, optional): Baseline correlation threshold to test
                against (the significance test asks whether the true correlation
                exceeds this value). Defaults to ``0.0``.
            alpha (float, optional): Significance level for the one-sided
                correlation test. Defaults to ``0.01``.
            intthres (float, optional): Minimum intensity required to keep a
                fragment (in addition to passing the correlation test).
                Defaults to ``0``.
            formulahit (bool, optional): When ``True``, also requires the
                fragment to have a molecular-formula annotation
                (``mzPAF``-style ``"/<ppm>ppm"`` suffix) within the MS2 averaged
                m/z threshold. Defaults to ``False``.
            reversed (bool, optional): When ``True``, returns the *removed*
                fragments instead of the kept ones (useful for inspecting what
                was filtered out). Defaults to ``False``.

        Returns:
            dict: Mapping of m/z (float) to intensity, containing either the
                fragments that passed the filter (default) or the fragments that
                were removed (when ``reversed=True``).
        """
        annotationlist = self.get_annotation_from_items()
        cleanedspectra = {}
        removedspectra = {}
        zscore = norm.ppf(q= 1.0-alpha, loc=0, scale=1)

        mzthreshold_ms2 = float(self.get_anything("COMMENT:", "MS2_M/Z_AVERAGED_THRES"))


        for mz, datahash in annotationlist.items():
            #print(mz, datahash)
            spearman_r = float(datahash["correlation"])
            intensity = datahash["intensity"]
            count = int(datahash["number"])
            #print(zscore, count, spearman_r)

            stderr = 1.0 / math.sqrt(count - 3)
            delta = zscore * stderr
            #lower = math.tanh(math.atanh(rthres) - delta)
            upper = math.tanh(math.atanh(rthres) + delta)

            if formulahit == True:
                mzPAF = datahash["annotation"]
                hitppm = re.search(r'/([0-9.-]+)ppm', mzPAF)

                if not hitppm == None:
                    hitppm = float(hitppm[1])
                    if abs(hitppm) > mzthreshold_ms2:
                        continue
                else:
                    continue

            if upper <= spearman_r  and intensity > intthres:
                cleanedspectra[mz] = intensity
                #print()
            else:
                removedspectra[mz] = intensity

        if reversed == True:
            cleanedspectra = removedspectra

        return(cleanedspectra)



    def mostcorrelaed_spectra(self, top = 3):
        """Select the top-N fragment ions by precursor-intensity correlation.

        Alternative to :meth:`wash_spectra`'s threshold-based filtering: simply
        keeps the ``top`` fragments with the highest correlation coefficient
        (falling back to keeping all fragments, sorted by m/z, when there are
        fewer than ``top`` of them).

        Args:
            top (int, optional): Number of most-correlated fragment ions to keep.
                Defaults to ``3``.

        Returns:
            dict: Mapping of m/z (float) to intensity for the selected fragments.
        """
        annotationlist = self.get_annotation_from_items()
        cleanedspectra = {}


        fragment_sorted = sorted(annotationlist.keys(), key = lambda x: float(annotationlist[x]["correlation"]), reverse= True)
        if len(fragment_sorted) < top:
            fragment_sorted = sorted(fragment_sorted)
        else:
            fragment_sorted = sorted(fragment_sorted[:top])


        for mz in fragment_sorted:
            intensity = annotationlist[mz]["intensity"]
            cleanedspectra[mz] = intensity
        return(cleanedspectra)


def calccorrel(x, y):
    """Compute the Spearman rank correlation coefficient between two intensity lists.

    Drops any index where either list has a blank (``""``) entry before
    computing, and short-circuits to ``0`` when either cleaned list is empty or
    has zero variance (``scipy.stats.spearmanr`` is undefined in that case).

    Args:
        x (list[float]): First list of values (e.g. per-scan intensities of one peak).
        y (list[float]): Second list of values, same length as ``x``.

    Returns:
        float: Spearman correlation coefficient between ``x`` and ``y``
            (via ``scipy.stats.spearmanr``), or ``0`` if it cannot be computed.
    """

    #
    # Compute correlation coefficients
    #
    v1_temp = []
    vt_temp = []
    #
    # Remove blank data
    #
    for (v1, vt) in zip(x,y):
        if v1 == "":
            continue
        if vt == "":
            continue
        v1_temp.append(v1)
        vt_temp.append(vt)

    #cor = np.corrcoef(v1_temp, vt_temp)
    #correl =  cor[0][1]
    # 250514 Data check
    if len(v1_temp) == 0:
        correl = 0
    elif len(vt_temp) == 0:
        correl = 0
    elif np.var(v1_temp) == 0.0:
        correl = 0
    elif np.var(vt_temp) == 0.0:
        correl = 0
    else:
        correl, pvalue = spearmanr(v1_temp, vt_temp)
    return(correl)



def calcThreshold(mz, msx, thresholds):
    """Compute the m/z-comparison threshold (in Da) for raw (un-averaged) spectra.

    Converts the configured raw ppm threshold (``"ms1_mz_raw_thres"`` or
    ``"ms2_mz_raw_thres"``, from ``RecordTemplate.txt``) into an absolute m/z
    tolerance for the given ``mz``.

    Args:
        mz (float): m/z value of the target ion.
        msx (str): ``"ms1"`` or ``"ms2"``, selecting which threshold to use.
        thresholds (dict): Threshold dictionary, as generated by
            ``ddatoolbox.ddatoolbox.setup_averaging``.

    Returns:
        float: Threshold value in Da (``mz * ppm_threshold / 1e6``).

    See Also:
        :func:`calcThresholdAveraged`: same, for averaged spectra.
    """

    thresholdunit = thresholds["unit"]
    ms2threshold = thresholds["ms2_mz_raw_thres"]
    ms1threshold = thresholds["ms1_mz_raw_thres"]

    mz = float(mz)

    if msx == "ms1":
        thresholdtemp = ms1threshold
    else:
        thresholdtemp = ms2threshold
    thresholdoutput = mz * thresholdtemp/1000000

    #if thresholdunit == "ppm":
    #    thresholdoutput = mz * thresholdtemp/1000000
    #else:
    #    thresholdoutput = thresholdtemp
    return thresholdoutput

def calcThresholdAveraged(mz, msx, thresholds):
    """Compute the m/z-comparison threshold (in Da) for averaged spectra.

    Converts the configured averaged-spectrum ppm threshold
    (``"ms1_mz_averaged_thres"`` or ``"ms2_mz_averaged_thres"``, from
    ``RecordTemplate.txt``) into an absolute m/z tolerance for the given ``mz``.

    Args:
        mz (float): m/z value of the target ion.
        msx (str): ``"ms1"`` or ``"ms2"``, selecting which threshold to use.
        thresholds (dict): Threshold dictionary, as generated by
            ``ddatoolbox.ddatoolbox.setup_averaging`` or
            :meth:`MassBankRecord.get_mzthresholds_averaged`.

    Returns:
        float: Threshold value in Da (``mz * ppm_threshold / 1e6``).

    See Also:
        :func:`calcThreshold`: same, for raw (un-averaged) spectra.
    """

    thresholdunit = thresholds["unit"]
    ms2threshold = thresholds["ms2_mz_averaged_thres"]
    ms1threshold = thresholds["ms1_mz_averaged_thres"]

    mz = float(mz)

    if msx == "ms1":
        thresholdtemp = ms1threshold
    else:
        thresholdtemp = ms2threshold
    thresholdoutput = mz * thresholdtemp/1000000

    #if thresholdunit == "ppm":
    #    thresholdoutput = mz * thresholdtemp/1000000
    #else:
    #    thresholdoutput = thresholdtemp
    return thresholdoutput


def read_RecordTemplate_dataset(dataset):
    """Read a dataset's ``RecordTemplate.txt`` by dataset identifier.

    Convenience wrapper around :func:`read_RecordTemplate` that resolves the
    dataset's template path as ``../data/<dataset>/RecordTemplate.txt``.

    Args:
        dataset (str): Dataset identifier (folder name under ``../data/``).

    Returns:
        dict: Parsed template settings, as returned by :func:`read_RecordTemplate`.

    See Also:
        :func:`read_RecordTemplate`: does the actual parsing.
    """
    templatepath = "../data/"+dataset+"/RecordTemplate.txt"
    return(read_RecordTemplate(templatepath))

def read_RecordTemplate(file):
    """Parse a ``RecordTemplate.txt`` file into a flat settings dictionary.

    ``RecordTemplate.txt`` describes dataset-wide constant MassBank record
    fields and pipeline parameters. Rows are dispatched as follows: ``PARAMETER:``
    rows become top-level settings keyed by their subfield name (with
    ``"checkmetabolites"`` accumulated as a list of ``eval``-parsed tuples rather
    than overwritten); ``AC$MASS_SPECTROMETRY: ION_MODE`` sets ``"polarity"``;
    ``AC$INSTRUMENT_TYPE: *`` sets ``"INSTRUMENT_TYPE"``;
    ``AC$CHROMATOGRAPHY: SOLVENT A``/``SOLVENT B``/``SOLVENT`` are collected into
    the ``"SOLVENT"`` list; ``COMMENT: DATASET_ID``/``REPOSITORY`` set
    ``"dataset_id"``/``"repository"``. Parsing stops at the ``//`` record
    terminator.

    Args:
        file (str): Path to the ``RecordTemplate.txt`` file to read.

    Returns:
        dict: Parsed settings, including (non-exhaustive) ``"checkmetabolites"``
            (list), ``"SOLVENT"`` (list), ``"polarity"``, ``"INSTRUMENT_TYPE"``,
            ``"dataset_id"``, ``"repository"``, plus one key per
            ``PARAMETER:`` row (e.g. ``"accession_stem"``,
            ``"ms1_mz_raw_thres"``, ``"unitmass_start"``, ...), and an
            (unused/empty) ``"all"`` key.

    See Also:
        :func:`read_RecordTemplate_dataset`: convenience wrapper by dataset name.
    """
    output_hash = {"checkmetabolites":[]}
    output_hash["SOLVENT"] = []
    output_hash["all"] = {}
    with open(file, 'r', encoding='UTF-8' ) as f:
        reader = csv.reader(f, delimiter='\t')
        mode = "none"

        for i, row in enumerate(reader):
            if len(row) < 2:
                continue
            if row[0] == "//":
                break
            fld1 = row[0]
            fld2 = row[1]
            if fld1 == "PARAMETER:":
                if fld2 == "checkmetabolites":
                    output_hash["checkmetabolites"].append(eval(row[2]))

                else:
                    output_hash[fld2] = row[2]
            if fld1 == "AC$MASS_SPECTROMETRY:" and fld2 == "ION_MODE":
                output_hash["polarity"] = row[2]
            if fld1 == "AC$INSTRUMENT_TYPE:" and fld2 == "*":
                output_hash["INSTRUMENT_TYPE"] = row[2]

            if fld1 == "AC$CHROMATOGRAPHY:" and fld2 == "SOLVENT A":
                output_hash["SOLVENT"].append("A "+row[2])
            if fld1 == "AC$CHROMATOGRAPHY:" and fld2 == "SOLVENT B":
                output_hash["SOLVENT"].append("B "+row[2])
            if fld1 == "AC$CHROMATOGRAPHY:" and fld2 == "SOLVENT":
                output_hash["SOLVENT"].append(row[2])
            if fld1 == "COMMENT:" and fld2 == "DATASET_ID":
                output_hash["dataset_id"] = row[2]
            if fld1 == "COMMENT:" and fld2 == "REPOSITORY":
                output_hash["repository"] = row[2]

            if fld1 == "COMMENT:" and fld2 == "DATASET_AUTHORS":
                output_hash["dataset_authors"] = row[2]
                
            if fld1 == "COMMENT:" and fld2 == "DATASET_TITLE":
                output_hash["dataset_title"] = row[2]

            if fld1 == "COMMENT:" and fld2 == "DATASET_PUBLICATION":
                output_hash["dataset_publication"] = row[2]
            if fld1 == "COMMENT:" and fld2 == "DATASET_PUBLICATION_PMID":
                output_hash["dataset_publication_pmid"] = row[2]

            if fld1 == "SP$SAMPLE:" and fld2 == "*":
                output_hash["sample"] = row[2]



            if fld1 == "AC$CHROMATOGRAPHY:" and fld2 == "COLUMN_NAME":
                output_hash["column_name"] = row[2]
            
            output_hash["all"][(fld1, fld2)] = row[2]

    return(output_hash)

def read_ionformlist(filename = "ion_forms.txt", polarity = "POSITIVE"):
    """Load the reference table of candidate precursor ion forms (adducts).

    Reads the bundled ``ion_forms.txt`` table, keeping only rows matching
    ``polarity`` and excluding comment rows (``#``) and multi-molecule/
    multiply-charged adducts (``mult != 1``, e.g. ``[M+2H]2+`` or ``[2M]`` dimers,
    which SIRIUS does not support).

    Args:
        filename (str, optional): Name of the bundled ion-form reference file.
            Defaults to ``"ion_forms.txt"``.
        polarity (str, optional): Ion polarity to keep, ``"POSITIVE"`` or
            ``"NEGATIVE"``. Defaults to ``"POSITIVE"``.

    Returns:
        dict: Mapping of ion-form ID (str, e.g. ``"[M+H]+"``) to a dict with
            ``"charge"`` (str), ``"mult"`` (float, always ``1`` after filtering),
            and ``"mass"`` (float, adduct mass shift).
    """

    reader = read_text("ddatoolbox", filename, encoding="utf-8").split('\n')
    ion_form = {}
    for i, rowtex in enumerate(reader):
        row = rowtex.split('\t')
        if len(row) < 6:
            continue
        modeid = row[0]
        if "#" in modeid:
            continue
        if not row[1] == polarity:
            continue

        charge = row[3]

        mult = float(row[4])
        if mult == 0.33:
            mult = 1/3
        #
        # SIRIUS does not support adducts such as [M+2H]2+ and dimers, which are represented by [2M]
        #
        if not mult == 1:
            continue

        mass = float(row[5])

        ion_form[modeid] = {"charge":charge, "mult":mult, "mass":mass}

    return(ion_form)

