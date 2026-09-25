#-------------------------------------------------------------------------------
# Name:        DDAAnnotation
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
# This module finalizes DDA (Data-Dependent Acquisition) MassBank records
# for the ShinMassBank pipeline. It builds on top of the `ddatoolbox`
# package (paths, massbankio, mfsearch) to:
#   - clean up/finalize step4 MassBank records for distribution (finalize)
#   - export dataset-level MassBank records as an mzTab-M file (generate_mztab)
#   - fetch sample/experimental-group metadata from public repositories
#     (Metabolomics Workbench / MetaboLights) used to populate mzTab metadata
#   - export records as MSP spectral libraries (generage_msp)
#   - summarize per-dataset processing/QC statistics (statisticspreprocessing)
#
import ddatoolbox
import csv, itertools, sys, numpy, glob
import os, re, shutil, math, json
import urllib.request, urllib.parse
from numpy import mean, std

#
# mzTab-M 2.0.0-M specification: https://github.com/HUPO-PSI/mzTab
# PSI-MS controlled vocabulary:  https://raw.githubusercontent.com/hupo-psi/psi-ms-cv/master/psi-ms.obo
#
# Version of the PSI-MS controlled vocabulary declared in the mzTab-M metadata.
PSI_MS_VERSION = "4.2.2"
# Formula of the "Uncharacterized" placeholder compound in RecordTemplate.txt.
# Records that were not (re-)written with a real structure still carry it.
PLACEHOLDER_FORMULA = "C39H42O2"
# mzTab-M notation of a missing value.
NULL = "null"
# Reliability (Schymanski confidence level, 1-5) of each identification level
# used in the CONFIDENCE comment. L3-CP (class-proposed) only carries
# formula-level evidence, like L4, so it is grouped with L4 rather than L3-SP.
RELIABILITY_BY_LEVEL = {"L2": "2", "L3-SP": "3", "L4": "4", "L3-CP": "4", "L5": "5"}
# Database identifier of an identified small molecule without a PubChem CID:
# the database prefix followed by null (a bare "null" is only allowed for
# unidentified molecules).
UNKNOWN_CID = "CID:null"


def _text(value):
    """Returns `value` as a single-line, tab-free string (mzTab is tab-delimited)."""
    return re.sub(r"[\t\r\n]+", " ", str(value)).strip()


def _number(value):
    """Formats a number for mzTab; None/NaN/inf become "null"."""
    if value is None or math.isnan(value) or math.isinf(value):
        return NULL
    return format(value, ".12g")


def _cell(value):
    """Formats one table cell: numbers via `_number`, empty text as "null"."""
    if isinstance(value, str):
        return _text(value) or NULL
    return _number(value)


def _param(label, accession, name, value=""):
    """Builds an mzTab parameter "[label, accession, name, value]".

    Brackets are not allowed inside a field and fields containing a comma
    have to be double-quoted, as required by the mzTab-M specification.
    User-defined parameters are written with an empty label and accession.
    """
    fields = []
    for field in (label, accession, name, value):
        field = _text(field).replace("[", "(").replace("]", ")").replace('"', "'")
        if "," in field:
            field = '"{}"'.format(field)
        fields.append(field)
    return "[{}, {}, {}, {}]".format(*fields)


def _taxonid(text):
    """Extracts the NCBI taxonomy id from "9606", "NCBITaxon:9606" or a PURL ending in "NCBITaxon_9606"."""
    match = re.search(r"(\d+)\s*$", str(text))
    if match:
        return match.group(1)
    return ""


def _getfield(mb2, fld1, fld2):
    """Reads a MassBank record field; a missing field is returned as an empty string."""
    try:
        return mb2.get_anything(fld1, fld2)
    except KeyError:
        return ""


def _getproposed(mb2, target):
    """Reads PROPOSED_COMPOUND/PROPOSED_MF as {id: {method: {rank, score, ...}}}; {} if absent or unparsable."""
    try:
        return mb2.get_proposed(target=target)["tg"]
    except (KeyError, ValueError, IndexError):
        return {}


def _spectral_entropy_gnps(mb2):
    """GNPS spectral-library match score of a record's EVIDENCE_INVENTORY comment.

    Reads the "SL (Spectral library) Entropy:GNPS (score X)" entry written by
    DDACuration.py for L2 (spectral-library-supported) identifications.

    Returns:
        float or None: The score, or None if the comment is absent or has no
            such entry (e.g. any level other than L2).
    """
    inventory = _getfield(mb2, "COMMENT:", "EVIDENCE_INVENTORY")
    if isinstance(inventory, list):
        inventory = inventory[0] if len(inventory) > 0 else ""
    match = re.search(r"Entropy:GNPS \(score ([\d.eE+-]+)\)", inventory)
    if match:
        return float(match.group(1))
    return None


def _adduct(adduct_ion, polarity):
    """Normalizes an ion form to the mzTab-M adduct notation and gives its charge.

    mzTab-M requires "[M+H]+", "[M-2H]2-" etc. Some records carry a charge
    without its sign ("[M-H]1"); the sign is then taken from the polarity.

    Args:
        adduct_ion (str): Ion form of a record, e.g. "[M-H]1".
        polarity (str): "POSITIVE" or "NEGATIVE".

    Returns:
        tuple: (adduct in mzTab-M notation or "null", charge as an int).
    """
    default = 1 if polarity == "POSITIVE" else -1
    match = re.match(r"^(\[.+\])(\d*)([+-]*)$", adduct_ion.replace(" ", ""))
    if not match:
        return NULL, default
    body, digits, signs = match.groups()
    count = int(digits) if digits else max(len(signs), 1)
    sign = signs[0] if signs else ("+" if default > 0 else "-")
    return "{}{}{}".format(body, count if count > 1 else "", sign), count if sign == "+" else -count


def _theoretical_mz(formula, adduct_ion):
    """Theoretical m/z of `formula` with `adduct_ion`; None if it cannot be calculated."""
    if formula == NULL or adduct_ion == NULL:
        return None
    try:
        precmf = ddatoolbox.mfsearch.generateprecmf(formula, adduct_ion)
        return ddatoolbox.mfsearch.calcmw(precmf)
    except Exception:
        # Unknown element/adduct notation: theoretical m/z is optional in mzTab-M.
        return None


def _repository(recordtemplate):
    """Repository name of a dataset; inferred from the study id if RecordTemplate.txt lacks it."""
    repository = recordtemplate.get("repository", "")
    studyid = recordtemplate.get("dataset_id", "")
    if repository == "":
        if re.match(r"^MTBLS\d+$", studyid):
            repository = "MetaboLights"
        elif re.match(r"^ST\d+$", studyid):
            repository = "Metabolomics Workbench"
    return repository


def _neutral_mass(formula):
    """Exact neutral mass of `formula`; None if it cannot be calculated."""
    try:
        return ddatoolbox.mfsearch.calcmw(formula)
    except Exception:
        # Unknown element in the formula: the neutral mass is optional in mzTab-M.
        return None


def generate_mztab(dataset, keyword=""):
    """Builds an mzTab-M 2.0.0-M file summarizing a dataset's identification and quantification results.

    Collects sample/study metadata from the source repository (Metabolomics
    Workbench or MetaboLights; for other repositories, or if the repository
    cannot be read, from the list of raw data files of the dataset), reads
    every step4 (finalized) MassBank record of the dataset via
    `ddatoolbox.massbankio.MassBankRecord`, and assembles the mzTab MTD
    (metadata), SML (small molecule), SMF (small molecule feature) and SME
    (small molecule evidence) sections following the mzTab-M 2.0.0-M
    specification (https://github.com/HUPO-PSI/mzTab).

    One MassBank record corresponds to one SMF and one SML row. Records
    identified at L2/L3-SP (structure) or L4/L3-CP (formula) additionally
    get an SME row; L5 records ("exact mass of interest") are reported as
    unidentified features. The SML "reliability" is the Schymanski
    confidence level (1-5, MS:1002955) of the record. The abundance of a
    feature in an assay is the maximum base peak intensity of its MS2 scans
    in that raw data file ("null" if it was not detected in the file);
    abundance_study_variable/abundance_variation_study_variable are the mean
    and the standard deviation over the detected assays of each study
    variable. Records that were detected only in raw data files that are not
    part of the mzTab metadata are omitted.

    Args:
        dataset (str): Name of the dataset (e.g. "HSA007"), matching a
            folder under "../data/".
        keyword (str, optional): If given, only step4 record files whose
            filename contains this substring are processed. Defaults to ""
            (process all records).

    Returns:
        None: The mzTab file is written to `paths["mztab"]`; nothing is returned.
    """
    #
    # Get path informataion of the dataset
    #
    paths = ddatoolbox.paths.generate_paths_data(dataset)
    recordtemplate = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)

    # Get polarity of this dataset
    if recordtemplate["polarity"] == "POSITIVE":
        scan_polarity = _param("MS", "MS:1000130", "positive scan")
    else:
        scan_polarity = _param("MS", "MS:1000129", "negative scan")

    fragmentation_methods = {
        "CID": _param("MS", "MS:1000133", "collision-induced dissociation"),
        "HCD": _param("MS", "MS:1000422", "beam-type collision-induced dissociation"),
    }
    fragmentation = fragmentation_methods.get(
        recordtemplate["all"].get(("AC$MASS_SPECTROMETRY:", "FRAGMENTATION_MODE"), "").strip())
    #
    # get metadata information from data repository using APIs
    #
    #assay_dict[assaykey][assayvalue] = set()
    #assay_dict["sample_type"] = sample_type
    #assay_dict["species"] = (species, taxonid)
    #groups_dict[groupkey][groupvalue] = set()
    #id2file[sampleid] = filenamestem
    #file2id[filenamestem] = sampleid
    groups_dict, assay_dict, id2file, file2id = getmetadatafromrepository(dataset)

    # The repository sample ids are not necessarily consecutive (samples of
    # other datasets/polarities are dropped), whereas the mzTab-M indexes
    # sample[n], ms_run[n] and assay[n] are. Number the used raw files 1..n.
    run_index = {}      # raw file stem -> mzTab index (sample = ms_run = assay)
    sample2index = {}   # repository sample id -> mzTab index
    for index, (file, sampleid) in enumerate(sorted(file2id.items(), key=lambda x: x[1]), 1):
        run_index[file] = index
        sample2index[sampleid] = index
    numrun = len(run_index)
    if numrun == 0:
        print("No raw data files (hmd) found for {}: mzTab-M file is not generated.".format(dataset))
        return

    # Study variables (experimental groups) as (name, description, assay indexes).
    study_variables = []
    for groupkey, groupvalues in groups_dict.items():
        for groupvalue, sampleids in groupvalues.items():
            assays = sorted(sample2index[x] for x in sampleids if x in sample2index)
            if _text(groupvalue) == "" or len(assays) == 0:
                continue
            study_variables.append(("{}_{}".format(groupkey, groupvalue),
                                    "{}: {}".format(groupkey, groupvalue), assays))
    if len(study_variables) == 0:
        # mzTab-M requires at least one study variable.
        study_variables.append(("all samples", "All samples of dataset {}".format(dataset),
                                list(range(1, numrun + 1))))

    # Identification confidence measures. Each entry is (name, level of the
    # PROPOSED_* comment, method names in that comment, "rank" or "score");
    # the best (smallest) value among the listed methods is reported.
    qlike = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue")
    qlike2 = ddatoolbox.paths.getsettings("RESOURCES", "qlikevalue2")
    zbra = ddatoolbox.paths.getsettings("RESOURCES", "zbramethod")
    id_confidence_measures = [
        ("MS-FINDER formula rank", "formula", ["MS-FINDER"], "rank"),
        ("SIRIUS6 formula rank", "formula", ["SIRIUS6_Sirius"], "rank"),
        ("q-value formula", "formula", [qlike], "score"),
        ("MS-FINDER compound rank", "compound", ["MS-FINDER:MS-FINDER_DB"], "rank"),
        ("SIRIUS6 CSI compound rank", "compound", ["SIRIUS6_CSI:biodatabase"], "rank"),
        ("MetFrag MS-FINDER DB compound rank", "compound", ["MetFrag:MS-FINDER_DB"], "rank"),
        ("MetFrag PubChemLite compound rank", "compound", ["MetFrag:PubChemLite"], "rank"),
        ("q-value compound", "compound", ["{}:{}".format(qlike, zbra), "{}:{}".format(qlike2, zbra)], "score"),
        # Computed separately (see ENTROPY_INDEX below): the GNPS spectral
        # library match score of the EVIDENCE_INVENTORY comment, reported
        # only for L2 (spectral-library-supported) identifications.
        ("Spectral entropy GNPS", "compound", [], "score"),
    ]
    QVALUE_FORMULA = 3      # position of "q-value formula" in id_confidence_measures
    QVALUE_COMPOUND = 8     # position of "q-value compound" in id_confidence_measures
    ENTROPY_INDEX = 9       # position of "Spectral entropy GNPS" in id_confidence_measures

    #
    # Generation of MTD (metadata) section
    #
    MTD = {}
    MTD["mzTab-version"] = "2.0.0-M"
    MTD["mzTab-ID"] = "ShinMassBank_{}".format(dataset)
    studyid = recordtemplate.get("dataset_id", "")
    repository = _repository(recordtemplate)
    if "dataset_title" in recordtemplate:
        MTD["title"] = _text(recordtemplate["dataset_title"])
        MTD["description"] = "Shin-MassBank averaged DDA (MS2) data generated from {} {}: {}".format(
            repository, studyid, _text(recordtemplate["dataset_title"]))
    MTD["software[1]"] = _param("MS", "MS:1000799", "custom unreleased software tool", "ShinMassBank pipeline")

    if repository == "MetaboLights":
        MTD["uri[1]"] = "https://www.ebi.ac.uk/metabolights/{}".format(studyid)
    elif repository == "Metabolomics Workbench":
        MTD["uri[1]"] = "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={}".format(studyid)

    MTD["quantification_method"] = _param("MS", "MS:1002019", "label-free raw feature quantitation")
    #
    # Generation of sample list
    #
    species, taxonid = assay_dict["species"]
    taxonid = _taxonid(taxonid)
    sample_type = assay_dict["sample_type"]

    # Repository characteristics other than species/sample type become
    # user-defined sample parameters of the samples having that value.
    customs = {}
    for assaykey, assayvalues in assay_dict.items():
        if assaykey in ["species", "sample_type"]:
            continue
        for assayvalue, sampleids in assayvalues.items():
            if _text(assayvalue) == "":
                continue
            for sampleid in sampleids:
                if sampleid in sample2index:
                    customs.setdefault(sample2index[sampleid], []).append((assaykey, assayvalue))

    for file, index in run_index.items():
        sample = "sample[{}]".format(index)
        MTD[sample] = _text(file)
        if species != "":
            if taxonid != "":
                MTD[sample + "-species[1]"] = _param("NCBITaxon", "NCBITaxon:{}".format(taxonid), species)
            else:
                MTD[sample + "-species[1]"] = _param("", "", species)
        if sample_type != "":
            MTD[sample + "-tissue[1]"] = _param("", "", sample_type)
        for customcounter, (assaykey, assayvalue) in enumerate(customs.get(index, []), 1):
            MTD[sample + "-custom[{}]".format(customcounter)] = _param("", "", assaykey, assayvalue)
    #
    # Generation of ms_run and assay list
    #
    for file, index in run_index.items():
        ms_run = "ms_run[{}]".format(index)
        MTD[ms_run + "-location"] = "file:///data/hive/{}.hmd".format(urllib.parse.quote(file))
        MTD[ms_run + "-format"] = _param("", "", "hive")
        MTD[ms_run + "-id_format"] = _param("MS", "MS:1000776", "scan number only nativeID format")
        if fragmentation is not None:
            MTD[ms_run + "-fragmentation_method[1]"] = fragmentation
        MTD[ms_run + "-scan_polarity[1]"] = scan_polarity

        assay = "assay[{}]".format(index)
        MTD[assay] = _text(file)
        MTD[assay + "-sample_ref"] = "sample[{}]".format(index)
        MTD[assay + "-ms_run_ref"] = ms_run
    #
    # Generation of study_variable
    #
    for counter, (name, description, assays) in enumerate(study_variables, 1):
        study_variable = "study_variable[{}]".format(counter)
        MTD[study_variable] = _text(name)
        MTD[study_variable + "-assay_refs"] = " | ".join("assay[{}]".format(x) for x in assays)
        MTD[study_variable + "-average_function"] = _param("MS", "MS:1002962", "mean")
        MTD[study_variable + "-variation_function"] = _param("", "", "standard deviation")
        MTD[study_variable + "-description"] = _text(description)
    #
    # Database information: PubChem is the source of the CIDs, "no database"
    # is required for records without a database identifier.
    #
    MTD["database[1]"] = _param("", "", "PubChem")
    MTD["database[1]-prefix"] = "CID"
    MTD["database[1]-version"] = "Unknown"
    MTD["database[1]-uri"] = "https://pubchem.ncbi.nlm.nih.gov/"
    MTD["database[2]"] = _param("", "", "no database", NULL)
    MTD["database[2]-prefix"] = NULL
    MTD["database[2]-version"] = "Unknown"
    MTD["database[2]-uri"] = NULL
    #
    # ontologies used
    #
    MTD["cv[1]-label"] = "MS"
    MTD["cv[1]-full_name"] = "PSI-MS controlled vocabulary"
    MTD["cv[1]-version"] = PSI_MS_VERSION
    MTD["cv[1]-uri"] = "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"
    if taxonid != "":
        MTD["cv[2]-label"] = "NCBITaxon"
        MTD["cv[2]-full_name"] = "NCBI organismal classification"
        MTD["cv[2]-version"] = "Unknown"
        MTD["cv[2]-uri"] = "http://purl.obolibrary.org/obo/ncbitaxon.owl"
    #
    # other definitions
    #
    # The abundance of a feature in an assay is the (raw) base peak intensity
    # of the MS2 spectra of this feature in the raw data file.
    MTD["small_molecule-quantification_unit"] = _param("MS", "MS:1000131", "number of detector counts")
    MTD["small_molecule_feature-quantification_unit"] = _param("MS", "MS:1000131", "number of detector counts")
    MTD["small_molecule-identification_reliability"] = _param("MS", "MS:1002955", "hr-ms compound identification confidence level")
    for counter, (name, level, methods, kind) in enumerate(id_confidence_measures, 1):
        MTD["id_confidence_measure[{}]".format(counter)] = _param("", "", name)
    #
    # Column definitions of the SML, SMF and SME tables (opt_ columns last)
    #
    SMLitems = ["SML_ID",
    "SMF_ID_REFS",
    "database_identifier",
    "chemical_formula",
    "smiles",
    "inchi",
    "chemical_name",
    "uri",
    "theoretical_neutral_mass",
    "adduct_ions",
    "reliability",
    "best_id_confidence_measure",
    "best_id_confidence_value"]
    SMLitems += ["abundance_assay[{}]".format(x) for x in range(1, numrun + 1)]
    SMLitems += ["abundance_study_variable[{}]".format(x) for x in range(1, len(study_variables) + 1)]
    SMLitems += ["abundance_variation_study_variable[{}]".format(x) for x in range(1, len(study_variables) + 1)]
    # Only one opt_ column per table: jmzTab-m (the mzTab-M validator) fails
    # on tables with many abundance columns as soon as they have two or more
    # opt_ columns. The MassBank accession leads to all further information
    # (confidence level, InChIKey, compound class) of the record.
    SMLitems += ["opt_global_ShinMassBank_Record_ID"]

    SMFitems = ["SMF_ID",
    "SME_ID_REFS",
    "SME_ID_REF_ambiguity_code",
    "adduct_ion",
    "isotopomer",
    "exp_mass_to_charge",
    "charge",
    "retention_time_in_seconds",
    "retention_time_in_seconds_start",
    "retention_time_in_seconds_end"]
    SMFitems += ["abundance_assay[{}]".format(x) for x in range(1, numrun + 1)]
    SMFitems += ["opt_global_ShinMassBank_Record_ID"]

    SMEitems = ["SME_ID",
    "evidence_input_id",
    "database_identifier",
    "chemical_formula",
    "smiles",
    "inchi",
    "chemical_name",
    "uri",
    "derivatized_form",
    "adduct_ion",
    "exp_mass_to_charge",
    "charge",
    "theoretical_mass_to_charge",
    "spectra_ref",
    "identification_method",
    "ms_level"]
    SMEitems += ["id_confidence_measure[{}]".format(x) for x in range(1, len(id_confidence_measures) + 1)]
    SMEitems += ["rank", "opt_global_ShinMassBank_Record_ID"]

    # SML/SMF/SME hold one row (dict of column name -> value) per record.
    SML = []
    SMF = []
    SME = []
    skipped = 0
    # Get file path list in step 4
    filelist = sorted(paths["step4inputlist"])
    # Main loop: one iteration per finalized (step4) MS2 MassBank record.
    for ms2filename in filelist:
        if not keyword == "":
            if not keyword in ms2filename:
                continue
            else:
                print(keyword, ms2filename)
        #
        # Loading MassBank recored MS2
        #
        mb2 = ddatoolbox.massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        #
        # Collection of spec data
        #
        accession = mb2.get_accession() # Get the MassBank accession id
        exp_mass_to_charge, sd = mb2.get_precmz()
        rt, sd = mb2.get_rt()
        retention_time_in_seconds = int(rt*60 * 10)/10

        polarity = mb2.get_polarity()
        #
        # Gathering signal intensity data: the raw MS2 scans (RAW_SPEC keys are
        # "<raw file stem>_<scan number>") that were averaged into this record.
        #
        intrawdata = eval(mb2.get_anything("COMMENT:" ,"RAW_SPEC"))
        intensity = {}          # assay index -> maximum base peak intensity
        spectra_ref_list = []
        for filescan, temp in intrawdata.items():
            file, scan = filescan.rsplit("_", 1)
            #
            # if this datafiles are not in metabdata
            #
            if not file in run_index:
                continue
            index = run_index[file]
            spectra_ref_list.append("ms_run[{}]:scan={}".format(index, scan))
            intensity[index] = max(intensity.get(index, 0.0), float(temp['bp']))
        if len(spectra_ref_list) == 0:
            skipped += 1
            continue
        spectra_ref = " | ".join(spectra_ref_list)

        smf_id = len(SMF) + 1
        #
        # Identification level: "L2 (Probable structure); <InChIKey short>",
        # "L3-SP (...); <InChIKey short>; ...", "L4 (...); <formula>; <ion form>; ...",
        # "L3-CP (...); <class>; <formula>; <ion form>" or "L5 (...); prec m/z ...".
        #
        confidence = _getfield(mb2, "COMMENT:", "CONFIDENCE")
        if isinstance(confidence, list):
            confidence = confidence[0] if len(confidence) > 0 else ""
        confidenceitems = [x.strip() for x in confidence.split(";")]
        level = confidenceitems[0].split(" ")[0]
        # The reliability (Schymanski level, 1-5) is mandatory; a record without
        # a confidence level is at least an exact mass of interest (5). L3-CP
        # (class-proposed) only carries formula-level evidence, like L4, so it
        # is reported at reliability 4, not 3.
        reliability = RELIABILITY_BY_LEVEL.get(level, "5")

        chemical_formula = NULL
        adduct_ion = NULL
        smiles = NULL
        inchi = NULL
        inchikey = NULL
        chemical_name = NULL
        database_identifier = NULL
        uri = NULL
        inchikey_short = ""
        anno_compound = {}
        anno_formula = {}
        is_evidence = False

        if level in ["L2", "L3-SP"] and _getfield(mb2, "CH$FORMULA:", "*") not in ["", PLACEHOLDER_FORMULA]:
            #
            # Structure level identification: the compound information is in the CH$ fields
            #
            is_evidence = True
            chemical_formula = _getfield(mb2, "CH$FORMULA:", "*")
            smiles = _getfield(mb2, "CH$SMILES:", "*") or NULL
            inchi = _getfield(mb2, "CH$IUPAC:", "*") or NULL
            inchikey = _getfield(mb2, "CH$LINK:", "INCHIKEY") or NULL
            inchikey_short = confidenceitems[1] if len(confidenceitems) > 1 else inchikey.split("-")[0]
            names = [x for x in _getfield(mb2, "CH$NAME:", "*") if x not in ["", "nan"]]
            if len(names) > 0:
                # "|" separates the items of list columns in mzTab-M
                chemical_name = names[0].replace("|", "/")
            database_identifier = _getfield(mb2, "CH$LINK:", "PUBCHEM") or UNKNOWN_CID
            hmdb = _getfield(mb2, "CH$LINK:", "HMDB")
            if hmdb != "":
                uri = "https://hmdb.ca/metabolites/{}".format(hmdb)
            elif database_identifier != UNKNOWN_CID:
                uri = "https://pubchem.ncbi.nlm.nih.gov/compound/{}".format(database_identifier.replace("CID:", ""))
            adduct_ion = _getfield(mb2, "MS$FOCUSED_ION:", "ION_TYPE") or NULL
            anno_compound = _getproposed(mb2, "compound")
            anno_formula = _getproposed(mb2, "formula")

        elif level in ["L4", "L3-CP"]:
            #
            # Formula level identification: the record still has the placeholder
            # compound, the formula and the ion form are in the CONFIDENCE comment.
            #
            offset = 1 if level == "L4" else 2
            if len(confidenceitems) > offset + 1 and confidenceitems[offset] != "N/A":
                is_evidence = True
                database_identifier = UNKNOWN_CID
                chemical_formula = confidenceitems[offset]
                adduct_ion = _getfield(mb2, "MS$FOCUSED_ION:", "ION_TYPE") or confidenceitems[offset + 1]
                anno_formula = _getproposed(mb2, "formula")

        adduct_ion, charge = _adduct(adduct_ion, polarity)
        #
        # Confidence measures of this record: formula level ones and, for
        # structure level identifications, compound level ones.
        #
        confidence_values = {}
        for counter, (name, target, methods, kind) in enumerate(id_confidence_measures, 1):
            if target == "formula":
                found = anno_formula.get(chemical_formula, {})
            else:
                found = anno_compound.get(inchikey_short, {})
            values = [found[x][kind] for x in methods if x in found]
            if len(values) > 0:
                confidence_values[counter] = min(values)
        if reliability == "2":
            # L2 (spectral-library-supported) identification: the GNPS score
            # is the "SL (Spectral library) Entropy:GNPS (score X)" entry of
            # the EVIDENCE_INVENTORY comment; null for every other level.
            entropy = _spectral_entropy_gnps(mb2)
            if entropy is not None:
                confidence_values[ENTROPY_INDEX] = entropy

        #
        # Generation of SME field
        #
        sme_id = NULL
        theoretical_mass_to_charge = None
        if is_evidence:
            sme_id = len(SME) + 1
            theoretical_mass_to_charge = _theoretical_mz(chemical_formula, adduct_ion)
            sme = {}
            sme["SME_ID"] = sme_id
            sme["evidence_input_id"] = smf_id
            sme["database_identifier"] = database_identifier
            sme["chemical_formula"] = chemical_formula
            sme["smiles"] = smiles
            sme["inchi"] = inchi
            sme["chemical_name"] = chemical_name
            sme["uri"] = uri
            sme["derivatized_form"] = NULL
            sme["adduct_ion"] = adduct_ion
            sme["exp_mass_to_charge"] = exp_mass_to_charge
            sme["charge"] = charge
            sme["theoretical_mass_to_charge"] = theoretical_mass_to_charge
            sme["spectra_ref"] = spectra_ref
            sme["identification_method"] = _param("", "", "ShinMassBank reanalysis pipeline")
            sme["ms_level"] = _param("MS", "MS:1000511", "ms level", "2")
            for counter in range(1, len(id_confidence_measures) + 1):
                sme["id_confidence_measure[{}]".format(counter)] = confidence_values.get(counter)
            sme["rank"] = 1
            sme["opt_global_ShinMassBank_Record_ID"] = accession
            SME.append(sme)
        #
        # Generation of SMF field
        #
        smf = {}
        smf["SMF_ID"] = smf_id
        smf["SME_ID_REFS"] = sme_id
        smf["SME_ID_REF_ambiguity_code"] = NULL
        smf["adduct_ion"] = adduct_ion
        smf["isotopomer"] = NULL
        smf["exp_mass_to_charge"] = exp_mass_to_charge
        smf["charge"] = charge
        smf["retention_time_in_seconds"] = retention_time_in_seconds
        smf["retention_time_in_seconds_start"] = NULL
        smf["retention_time_in_seconds_end"] = NULL
        for index in range(1, numrun + 1):
            smf["abundance_assay[{}]".format(index)] = intensity.get(index)
        smf["opt_global_ShinMassBank_Record_ID"] = accession
        SMF.append(smf)
        #
        # Generation of SML field
        #
        sml = {}
        sml["SML_ID"] = smf_id
        sml["SMF_ID_REFS"] = smf_id
        sml["database_identifier"] = database_identifier
        sml["chemical_formula"] = chemical_formula
        sml["smiles"] = smiles
        sml["inchi"] = inchi
        sml["chemical_name"] = chemical_name
        sml["uri"] = uri
        if chemical_formula != NULL:
            sml["theoretical_neutral_mass"] = _neutral_mass(chemical_formula)
        sml["adduct_ions"] = adduct_ion
        sml["reliability"] = reliability
        # The best confidence measure is the q-value of the identification.
        if is_evidence and inchikey_short != "" and QVALUE_COMPOUND in confidence_values:
            best = QVALUE_COMPOUND
        elif is_evidence and QVALUE_FORMULA in confidence_values:
            best = QVALUE_FORMULA
        else:
            best = None
        if best is not None:
            sml["best_id_confidence_measure"] = _param("", "", id_confidence_measures[best - 1][0])
            sml["best_id_confidence_value"] = confidence_values[best]
        for index in range(1, numrun + 1):
            sml["abundance_assay[{}]".format(index)] = intensity.get(index)

        # Aggregate the intensities of the detected assays into per-study-variable
        # (group) mean and standard deviation.
        for counter, (name, description, assays) in enumerate(study_variables, 1):
            values = [intensity[x] for x in assays if x in intensity]
            if len(values) > 0:
                sml["abundance_study_variable[{}]".format(counter)] = sum(values)/len(values)
                sml["abundance_variation_study_variable[{}]".format(counter)] = float(std(values))
        sml["opt_global_ShinMassBank_Record_ID"] = accession
        SML.append(sml)

    #
    # Serialize MTD, then SML, SMF, SME sections (each preceded by its own
    # header row and separated by a blank line) into a tab-separated mzTab file.
    # Missing values are written as "null".
    #
    def table(headertag, rowtag, columns, rows):
        lines = [headertag + "\t" + "\t".join(columns)]
        for row in rows:
            lines.append(rowtag + "\t" + "\t".join(_cell(row.get(x)) for x in columns))
        return lines

    output = []
    for feat, text in MTD.items():
        output.append("MTD\t{}\t{}".format(feat, text))
    output.append("")
    output += table("SMH", "SML", SMLitems, SML)
    output.append("")
    output += table("SFH", "SMF", SMFitems, SMF)
    if len(SME) > 0:
        output.append("")
        output += table("SEH", "SME", SMEitems, SME)

    with open(paths["mztab"], 'w',  encoding='UTF-8', newline='') as f:
        f.write("\n".join(output) + "\n")
    print("{}: {} features, {} evidences, {} records skipped (no scan in the registered raw data files)".format(
        dataset, len(SMF), len(SME), skipped))

def finalize(dataset, keyword = ""):
    """Produces the final, publication-ready MassBank records for a dataset.

    For every step4 MS2 record: washes (filters) the averaged MS2 spectrum,
    rewrites its annotation comment fields, recalibrates the corresponding
    step1 MS1 spectrum m/z values around the precursor m/z, stores the
    recalibrated MS1 spectrum as a comment, strips intermediate/working
    COMMENT fields that are no longer needed, and writes the cleaned record
    to `paths["finalpath"]`.

    Args:
        dataset (str): Name of the dataset (e.g. "HSA007"), matching a
            folder under "../data/".
        keyword (str, optional): If given, only step4 record files whose
            filename contains this substring are processed. Defaults to ""
            (process all records).

    Returns:
        None: Finalized MassBank record files are written directly to
        `paths["finalpath"]`; nothing is returned.
    """

    paths = ddatoolbox.paths.generate_paths_data(dataset)

    filelist =paths["step4inputlist"]
    finalpath = paths["finalpath"]
    outputpath = paths["finalpath"]+"/"

    # Wipe out any previous (partial) output for this dataset before
    # regenerating it from scratch.
    removepath = paths["root"] + "/" + dataset
    if os.path.isdir(removepath):
       shutil.rmtree(removepath)

    if not os.path.isdir(finalpath):
        os.mkdir(finalpath)


    settings = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)
    accession_stem = settings["accession_stem"]
    accession_stem_ms2 = accession_stem+"P"
    accession_stem_ms1 = accession_stem+"S"


    # Main loop: one iteration per finalized (step4) MS2 MassBank record.
    for ms2filename in filelist:
        if not keyword == "":
            if not keyword in ms2filename:
                continue
            else:
                print(keyword, ms2filename)
        # Loading MassBank recored MS2
        #
        mb2 = ddatoolbox.massbankio.MassBankRecord()
        mb2.read_MassBankRecord(ms2filename)
        accession = mb2.get_accession()

        # Remove low-quality/uncorrelated peaks from the averaged MS2 spectrum.
        washedspec = mb2.wash_spectra(alpha=0.01, intthres= 0)
        mb2.add_spectra(washedspec)

        annotation = mb2.get_annotation_from_items()
        # reserved = ["m/z", "m/z_std", "m/z_nocal", "m/z_std(ppm)", "int", "int_std", "cor", "n", "intlist", "mzlist", "annotation", "mzPAF"]
        #print(annotation)
        # Re-write the peak annotation comment, keeping only these fields.
        mb2.add_annotation_from_items(annotation, ["int", "cor", "n", "mzPAF"])

        precmz, sd = mb2.get_precmz()

        # The corresponding raw MS1 record lives under step1, with the MS2
        # accession stem swapped for the MS1 one.
        ms1filename = ms2filename.replace(accession_stem_ms2, accession_stem_ms1)
        ms1filename = ms1filename.replace("step4", "step1")
        mb1 = ddatoolbox.massbankio.MassBankRecord()
        mb1.read_MassBankRecord(ms1filename)
        #
        # Recalibration of ms1 m/z
        #
        annotationlist = mb1.get_annotation_from_items()

        spec_corrected = {}

        # Keep only MS1 peaks within +/-50 m/z of the precursor and with
        # intensity > 1.0, then apply the dataset's m/z correction curve.
        for mzraw, datahash in annotationlist.items():
            intensity = datahash["intensity"]
            if precmz - 50 < mzraw < precmz + 50 and intensity > 1.0:
                mz_corrected = ddatoolbox.mfsearch.mzcorrection(mzraw, settings["ms1_mz_correction"])
                mz_corrected = int(mz_corrected * 1000000) / 1000000
                spec_corrected[mz_corrected] = int(intensity)

        mb2.add2("COMMENT:", "MS1_SPEC",str(spec_corrected))
        #
        # Discard working/intermediate COMMENT fields that were only
        # needed during processing and should not appear in the final record.
        #
        mb2.clear_anything("COMMENT:", "MS1_M/Z_AVERAGED_THRES")
        mb2.clear_anything("COMMENT:", "MS2_M/Z_AVERAGED_THRES")
        mb2.clear_anything("COMMENT:", "DATASET_FILE_NUM")
        mb2.clear_anything("COMMENT:", "NUM_AVERAGED_SPEC")
        mb2.clear_anything("COMMENT:", "RAW_SPEC")
        mb2.clear_anything("COMMENT:", "PROPOSED_ISOTOPE_STATE")
        mb2.clear_anything("COMMENT:", "PROPOSED_MAX_NUM_S")
        mb2.clear_anything("COMMENT:", "PRECURSOR_M/Z_RAW")
        mb2.clear_anything("COMMENT:", "PRECURSOR_M/Z_STDEV")
        mb2.clear_anything("COMMENT:", "RETENTION_TIME_STDEV")
        mb2.show(path = outputpath)

        print('\r', "{}".format(accession), end = "")

def getmetadatafromrepository(dataset):
    """Dispatches metadata retrieval to the repository named in the dataset's RecordTemplate.

    Reads "dataset_id" and "repository" from the dataset's RecordTemplate.txt
    (via `ddatoolbox.massbankio.read_RecordTemplate_dataset`) and delegates to
    the matching repository-specific fetcher: Metabolomics Workbench
    (`getmwtabfromMetabolomicsWB`) or MetaboLights
    (`getmetadatafromMetabolights`). Also builds the set of raw data file
    stems ("hmd" files) actually present for this dataset, which the
    delegate functions use to filter out samples not included in it.

    Args:
        dataset (str): Name of the dataset (e.g. "HSA007"), matching a
            folder under "../data/".

    Returns:
        tuple: A 4-tuple of
            - groups_dict (dict): Maps experimental factors to sets of sample IDs.
            - assay_dict (dict): Maps assay metadata to sets of sample IDs,
              including "species" and "sample_type".
            - rawfile_dict (dict): Maps sample IDs to raw file stems.
            - file2id (dict): Maps raw file stems to their sample IDs.
            When the repository is neither "Metabolomics Workbench" nor
            "MetaboLights", cannot be read, or has none of the raw data
            files of the dataset registered, the metadata is generated from
            the raw data files of the dataset by `_default_metadata`.
    """
    recordtemplate = ddatoolbox.massbankio.read_RecordTemplate_dataset(dataset)
    hivefilelist = ddatoolbox.paths.generate_paths_data(dataset)["hmdinputlist"]
    # Build stem -> filename map of raw data files actually present in this
    # dataset, used by the repository-specific parsers to skip unused samples.
    datafiles = {}
    for file in hivefilelist:
        file = file.replace("\\", "/")
        filename = file.split("/")[-1]
        stem = filename.split(".")[0]
        datafiles[stem] = filename

    studyid = recordtemplate.get("dataset_id", "")
    repository = _repository(recordtemplate)
    groups_dict = {}
    assay_dict = {}
    rawfile_dict = {}
    file2id = {}
    try:
        if repository == "Metabolomics Workbench":
            groups_dict, assay_dict, rawfile_dict, file2id = getmwtabfromMetabolomicsWB(studyid, datafiles)
        elif repository == "MetaboLights":
            groups_dict, assay_dict, rawfile_dict, file2id = getmetadatafromMetabolights(studyid, datafiles)
    except (OSError, ValueError, KeyError) as e:
        # Network (urllib errors are OSErrors) or unexpected repository format.
        print("Could not read the metadata of {} {} ({}: {}).".format(repository, studyid, type(e).__name__, e))
        file2id = {}

    if len(file2id) == 0:
        # Unsupported/unknown repository or none of the raw data files is
        # registered in the repository: use the raw data files themselves.
        print("Sample metadata of {} is generated from the raw data files.".format(dataset))
        groups_dict, assay_dict, rawfile_dict, file2id = _default_metadata(datafiles, recordtemplate)

    return(groups_dict, assay_dict, rawfile_dict, file2id)


def _default_metadata(datafiles, recordtemplate):
    """Generates sample metadata without a repository: one sample per raw data file.

    Used for datasets whose repository is not supported (e.g. MB-POST) or
    unavailable. Every raw data file becomes a sample, all samples form one
    study variable, and species and sample type are taken from the SP$
    fields of the dataset's RecordTemplate.txt.

    Args:
        datafiles (dict): Raw data file stems of the dataset (see
            `getmetadatafromrepository`).
        recordtemplate (dict): The dataset's parsed RecordTemplate.txt.

    Returns:
        tuple: The same 4-tuple as `getmetadatafromrepository`.
    """
    rawfile_dict = {}
    file2id = {}
    for sampleid, file in enumerate(sorted(datafiles), 1):
        rawfile_dict[sampleid] = file
        file2id[file] = sampleid

    fields = recordtemplate["all"]
    species = fields.get(("SP$SCIENTIFIC_NAME:", "*"), "").strip()
    taxonid = ""
    match = re.search(r"NCBI-TAXONOMY\s+(\d+)", fields.get(("SP$LINK:", "*"), ""))
    if match:
        taxonid = match.group(1)

    groups_dict = {"dataset": {recordtemplate.get("dataset_id", "all samples"): set(rawfile_dict)}}
    assay_dict = {}
    assay_dict["sample_type"] = fields.get(("SP$SAMPLE:", "*"), "").strip()
    assay_dict["species"] = (species, taxonid)
    return(groups_dict, assay_dict, rawfile_dict, file2id)


def getmetadatafromMetabolights(studyid, datafiles):
    """Fetches and parses study metadata from MetaboLights.

    Downloads the sample sheet ("s_<studyid>.txt") for the given study from
    the MetaboLights web service, evaluates it as a Python literal (the
    endpoint returns a JSON-like structure using Python's True/False), and
    extracts, per sample row, the experimental factors ("Factor Value[...]"
    columns) and characteristics ("Characteristics[...]" columns, including
    species/organism part), restricted to samples whose source file stem is
    present in `datafiles`.

    Args:
        studyid (str): The MetaboLights study identifier (e.g. 'MTBLS...').
        datafiles (dict): Dict of available raw data file stems in the dataset
            (as produced by `getmetadatafromrepository`), used to skip
            samples not included in this dataset.

    Returns:
        tuple: A 4-tuple of
            - groups_dict (dict): Maps experimental factors ("Factor Value")
              to {value: set(sample_id)}.
            - assay_dict (dict): Maps characteristics ("Characteristics") to
              {value: set(sample_id)}, plus "species" (tuple of species name
              and taxon id) and "sample_type".
            - rawfile_dict (dict): Maps sample IDs to raw file stems.
            - file2id (dict): Maps raw file stems to their sample IDs.
    """

    url = "https://www.ebi.ac.uk/metabolights/ws/studies/{}/s_{}.txt".format(studyid, studyid)
    # Send a GET request to the MetaboLights API

    groups_dict = {}
    assay_dict = {}
    rawfile_dict = {}
    file2id = {}

    species = ""
    sample_type = ""
    taxonid = ""

    with urllib.request.urlopen(url, timeout=60) as response:
        # Read and decode the HTTP response body
        lines = response.read().decode('utf-8').split('\n')

    # The response body is a single line of JSON text.
    lines = json.loads(lines[0])
    for datafile in lines["data"]["rows"]:
        sampleid = datafile["index"]
        sampleid = int(sampleid) + 1
        filename = datafile["Source Name"]
        species = datafile.get("Characteristics[Organism]", "")
        sample_type = datafile.get("Characteristics[Organism part]", "")
        taxonid = datafile.get("Term Accession Number", "")
        filenamestem = filename.split(".")[0]

        #print(sampleid, filename, species, sample_type, taxonid)
        #
        # Skipping datafiles not in the dataset (pos/neg, analysis methods)
        #

        if not filenamestem in datafiles:
            continue
        rawfile_dict[sampleid] = filenamestem
        file2id[filenamestem] = sampleid

        # "Factor Value[...]" columns become experimental groups;
        # "Characteristics[...]" columns become assay metadata.
        for keyofdatafile in datafile.keys():
            if "Factor Value" in keyofdatafile:
                group = keyofdatafile.replace("Factor Value[","")
                group = group.replace("]","")
                if group == "":
                    group = "No label"+str(sampleid)

                if not group in groups_dict:
                    groups_dict[group] = {}
                groupvalue = datafile[keyofdatafile]
                if not groupvalue in groups_dict[group]:
                    groups_dict[group][groupvalue] = set()
  
                groups_dict[group][groupvalue].add(sampleid)
                #print(group,groupvalue,sampleid)
            if "Characteristics" in keyofdatafile:
                assay = keyofdatafile.replace("Characteristics[","")
                assay = assay.replace("]","")
                if not assay in assay_dict:
                    assay_dict[assay] = {}
                assayvalue = datafile[keyofdatafile]
                if not assayvalue in assay_dict[assay]:
                    assay_dict[assay][assayvalue] = set()
  
                assay_dict[assay][assayvalue].add(sampleid)
                #print(assay,assayvalue,sampleid)

    assay_dict["sample_type"] = sample_type
    assay_dict["species"] = (species, taxonid)




    """
    {"index": 0,
    "Source Name": "CS61791_pos_IDA",
    "Characteristics[Organism]": "Homo sapiens",
    "Term Source REF": "NCBITAXON",
    "Term Accession Number": "http://purl.obolibrary.org/obo/NCBITaxon_9606",
    "Characteristics[Organism part]": "blood serum",
    "Term Source REF.1": "BTO",
    "Term Accession Number.1": "http://purl.obolibrary.org/obo/BTO_0000133",
    "Characteristics[Variant]": "",
    "Term Source REF.2": "",
    "Term Accession Number.2": "",
    "Characteristics[Sample type]": "",
    "Term Source REF.3": "",
    "Term Accession Number.3": "",
    "Protocol REF": "Sample collection",
    "Sample Name": "CS61791_pos_IDA",
    "Factor Value[Sample group]": "Control",
    "Term Source REF.4": "",
    "Term Accession Number.4": "",
    "Factor Value[Gender]": "",
    "Term Source REF.5": "",
    "Term Accession Number.5": ""},


    for assaykey in assay_dict:
        for assayvalue in assay_dict[assaykey]:
            print(assaykey, assayvalue)
            #print(assay_dict[assaykey][assayvalue])
    for groupkey in groups_dict:
        for groupvalue in groups_dict[groupkey]:
            print(groupkey, groupvalue)
            print(groups_dict[groupkey][groupvalue])
    for fileid in rawfile_dict:
        print(studyid, fileid, rawfile_dict[fileid], datafiles[rawfile_dict[fileid]])
    """
         
    return(groups_dict, assay_dict, rawfile_dict,file2id)


def getmwtabfromMetabolomicsWB(studyid, datafiles):
    """Fetches and parses study metadata from the Metabolomics Workbench.

    Retrieves metadata in mwTab text format for a specific study ID from the
    Metabolomics Workbench REST API. It parses the metadata to extract
    experimental groups (factors), assay details, and maps sample IDs to
    the provided data files.

    Args:
        studyid (str): The Metabolomics Workbench study identifier (e.g., 'ST002338').
        datafiles (dict): A dictionary of available data file stems in the dataset.

    Returns:
        tuple: A tuple containing:
            - groups_dict (dict): Maps experimental factors to sets of sample IDs.
            - assay_dict (dict): Maps assay metadata to sets of sample IDs,
              including species and sample type.
            - rawfile_dict (dict): Maps sample IDs to raw file stems found in datafiles.
            - file2id (dict): Maps raw file stems to their corresponding sample IDs.
    """


    url = "https://www.metabolomicsworkbench.org/rest/study/study_id/{}/mwtab/txt".format(studyid)
    # Send a GET request to the Metabolomics Workbench API

    groups_dict = {}
    assay_dict = {}
    rawfile_dict = {}
    file2id = {}

    species = ""
    sample_type = ""
    taxonid = ""

    with urllib.request.urlopen(url, timeout=60) as response:
        # Read and decode the HTTP response body
        lines = response.read().decode('utf-8').split('\n')
    sampleid = 1
    # For study ST002338 the raw file names recorded in mwTab use the
    # opposite polarity suffix from the actual files on disk, so
    # replaceword1/2 are used below to swap "neg"<->"pos" when matching.
    replaceword1 = ""
    replaceword2 = ""
    # mwTab is a tab-separated, section-based text format; each row starting
    # with an empty first column marks the start of a new section, and rows
    # starting with "SUBJECT_SAMPLE_FACTORS" carry per-sample metadata.
    for line in lines:

        line = re.sub(' +\t', '\t', line)
        line = re.sub('\r', '', line)
        rows = line.split('\t')
        #print(line)
        if rows[0] == "":
            sampleid = 1
            if studyid == "ST002338":
                replaceword1 = "neg"
                replaceword2 = "pos"                
            #break

        if rows[0] == "SU:SUBJECT_SPECIES":
            species = rows[1]
        if rows[0] == "SU:TAXONOMY_ID":
            taxonid = rows[1]
        if rows[0] == "CO:SAMPLE_TYPE":
            sample_type = rows[1]

        if rows[0] == "SUBJECT_SAMPLE_FACTORS":
            sampleid += 1
            groups = rows[3].split(" | ")
            for group in groups:
                groupkey, groupvalue = group.split(":", 1)
                if not groupkey in groups_dict:
                    groups_dict[groupkey] = {}
                if not groupvalue in groups_dict[groupkey]:
                    groups_dict[groupkey][groupvalue] = set()
                groups_dict[groupkey][groupvalue].add(sampleid)
            print(rows)
            
            if rows[4] == "":
                filenamestem = rows[2]
                if filenamestem in datafiles:
                    rawfile_dict[sampleid] = filenamestem
                    file2id[filenamestem] = sampleid
                continue

            assays = rows[4].split("; ")
            for assay in assays:
                assaykey, assayvalue = assay.split("=", 1)

                if assaykey == "RAW_FILE_NAME":
                    filenames = assayvalue.split(" ")
                    for filename in filenames:
                        
                        filename = filename.replace(replaceword1 , replaceword2)
                        print(filename)

                        filenamestem = filename.split(".")[0]

                        
                        if filenamestem in datafiles:
                            rawfile_dict[sampleid] = filenamestem
                            file2id[filenamestem] = sampleid
                            break
                    continue
                if not assaykey in assay_dict:
                    assay_dict[assaykey] = {}
                if not assayvalue in assay_dict[assaykey]:
                    assay_dict[assaykey][assayvalue] = set()
                assay_dict[assaykey][assayvalue].add(sampleid)

            if not sampleid in rawfile_dict:
                filenamestem = rows[2]
                if filenamestem in datafiles:
                    rawfile_dict[sampleid] = filenamestem
                    file2id[filenamestem] = sampleid
    assay_dict["sample_type"] = sample_type
    assay_dict["species"] = (species, taxonid)


    for assaykey in assay_dict:
        for assayvalue in assay_dict[assaykey]:
            print(assaykey, assayvalue)
            #print(assay_dict[assaykey][assayvalue])
    for groupkey in groups_dict:
        for groupvalue in groups_dict[groupkey]:
            print(groupkey, groupvalue)
            print(groups_dict[groupkey][groupvalue])
    for fileid in rawfile_dict:
        print(studyid, fileid, rawfile_dict[fileid], datafiles[rawfile_dict[fileid]])

         
    return(groups_dict, assay_dict, rawfile_dict,file2id)



if __name__ == '__main__':
    # Simple manual entry point: `datasets` from CLI args is currently
    # unused below, and `dataset`/`keyword` are hardcoded for ad-hoc runs.
    args = sys.argv
    if len(args) >= 2:
        datasets = args[1:]
    else:
        dataset = "HSA001"


    keyword = "205016"
    keyword = ""

    #finalize(dataset, keyword=keyword)
    generate_mztab(dataset, keyword=keyword)


